!> Contains the routines for the radiation processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_rad_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control, ONLY: acc_stream
  USE mo_kind,        ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC ::                        &
    & calc_radiation_surface_net,  &
    & calc_par,                    &
    & calc_alb_lwtr, calc_alb_lice, calc_glacier_albedo, Has_minimal_radiation, calc_snow_age,                 &
    & calc_soil_albedo, calc_pond_albedo, calc_snow_albedo, calc_sky_view_fractions, Merge_albedos_of_vegtile, &
    & get_surface_albedo_simple, get_litter_carbon, get_soil_carbon
#ifndef __NO_QUINCY__
  PUBLIC :: calc_canopy_radiation_budget
#endif

  INTERFACE Has_minimal_radiation
    MODULE PROCEDURE Has_minimal_radiation_1
    MODULE PROCEDURE Has_minimal_radiation_2
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_process'

CONTAINS

  SUBROUTINE calc_radiation_surface_net( &
    & swvis_down,    & ! in
    & swnir_down,    & ! in
    & alb_vis,       & ! in
    & alb_nir,       & ! in
    & lw_down,       & ! in
    & t,             & ! in
    & rad_net,       & ! output
    & swvis_net,     & ! output, optional
    & swnir_net,     & ! output, optional
    & sw_net,        & ! output, optional
    & lw_net         & ! output, optional
    & )

    USE mo_phy_schemes,   ONLY: lwnet_from_lwdown

    REAL(wp),  INTENT(in), DIMENSION(:) :: &
      & swvis_down,          &
      & swnir_down,          &
      & alb_vis,             &
      & alb_nir,             &
      & lw_down,             &
      & t

    REAL(wp),  INTENT(out) :: &
      & rad_net(:)

    REAL(wp),  INTENT(out), OPTIONAL, DIMENSION(:) :: &
      & swvis_net,                      & ! net shortwave radiation in visible range
      & swnir_net,                      & ! net shortwave radiation in near infrared range
      & sw_net,                         &
      & lw_net

    REAL(wp) :: zswvis_net, zswnir_net, zsw_net, zlw_net

    INTEGER :: nc, ic

    nc = SIZE(swvis_down)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(zswvis_net, zswnir_net, zsw_net, zlw_net)
    DO ic=1,nc

      ! Compute net SW radiation from downward SW and albedo
      zswvis_net = swvis_down(ic) * (1._wp - alb_vis(ic))
      zswnir_net = swnir_down(ic) * (1._wp - alb_nir(ic))
      zsw_net    = zswvis_net + zswnir_net

      ! Compute LW net radiation from incoming and the thermal radiation from the surface
      zlw_net = lwnet_from_lwdown(lw_down(ic), t(ic))

      ! Compute net radiation
      rad_net(ic) = zsw_net + zlw_net

      IF (PRESENT(swvis_net)) swvis_net(ic) = zswvis_net
      IF (PRESENT(swnir_net)) swnir_net(ic) = zswnir_net
      IF (PRESENT(sw_net))    sw_net(ic)    = zsw_net
      IF (PRESENT(lw_net))    lw_net(ic)    = zlw_net

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_radiation_surface_net

  SUBROUTINE calc_par( &
    & ncanopy,                  & ! in
    & icanopy,                  & ! in
    & use_alb_veg_simple,       & ! in
    & par,                      & ! in
    & alb_vis_soil,             & ! in
    & lai,                      & ! in
    & cos_zenith_angle,         & ! in
    & fract_par_direct,         & ! in
    & canopy_bound_lai,         & ! in
    & canopy_bound_lai_delta,   & ! in
    & B4_layer_above,           & ! inout
    & par_down_mol,             & ! out
    & soil_reflectivity_par,    & ! out
    & faPAR_cl,                 & ! out
    & lai_cl,                   & ! out
    & apar_per_lai_cl  )          ! out

    !$ACC ROUTINE SEQ

    USE mo_rad_constants,  ONLY: &
      & Epar,                    & ! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1
      & SoilReflectivityParMin,  & ! Minimum soil reflectivity in PAR region
      & FcMax,                   &
      & FcMin,                   &
      & ZenithMinPar,            &
      & LaiMin,                  &
      & LaiLimit

    INTEGER,   INTENT(in)    :: ncanopy, &
                                icanopy
    LOGICAL,   INTENT(in)    :: use_alb_veg_simple
    REAL(wp),  INTENT(in)    :: canopy_bound_lai,      &
                                canopy_bound_lai_delta

    REAL(wp),  INTENT(in)    :: par,                   &
                                alb_vis_soil,          &
                                lai,                   &
                                cos_zenith_angle,      &
                                fract_par_direct
    REAL(wp),  INTENT(inout) :: B4_layer_above
    REAL(wp),  INTENT(out)   :: par_down_mol,          &
                                soil_reflectivity_par, &
                                faPAR_cl,              &
                                lai_cl,                &
                                apar_per_lai_cl

    REAL(wp) :: ZH, K0, ZP1, Q0, Q1, B0, B1,    & ! These variables vary per grid cell
                EKL, EHL, B4, FC, FSHD,         & ! These variables vary per grid cell, FC =fractional vegetation cover
                ZP0, EKL0, EHL0, X0, X1, X2, F    ! These variables do not vary per grid cell

    REAL(wp), PARAMETER :: OMEGA = 0.12_wp   ! single leaf scattering albedo

    ! ---------------------------
    ! Go

    ! Expresse radiation in  mol(photons)/(m^2 s)
    par_down_mol = par / Epar

    ! Compute soil reflectivity of PAR
    IF (.NOT. use_alb_veg_simple) THEN
      ! soil reflectivity is set to soil albedo of the visible range
      soil_reflectivity_par = alb_vis_soil
!        soil_reflectivity_par = alb_background ! USE ONLY FOR COMPARISON WITH JSBACH3
    ELSE
      ! soil reflectivity is derived from background albedo of the whole solar spectrum (compare eq. (122) in Knorr)
      soil_reflectivity_par = MAX(0.92_wp * alb_vis_soil - 0.015_wp, SoilReflectivityParMin)
!       soil_reflectivity_par = MAX(0.92_wp * alb_background - 0.015_wp, SoilReflectivityParMin) ! USE ONLY FOR COMPARISON WITH JS3
    END IF

    ! Initialize ABSORBED PAR PER LEAF AREA,  IF cos_zenith_angle<1e-3
    ! Spread LAI equally around the three layers
    faPAR_cl = 0._wp
    apar_per_lai_cl = 0._wp
    lai_cl   = lai * canopy_bound_lai_delta
    lai_cl   = MAX(lai_cl, LaiMin) ! ensure minimum LAI

    ! IF (cos_zenith_angle < ZenithMinPar) THEN ! do not calculate PAR for zenith angles below ZenithMinPar (0.001)
    !   RETURN
    ! ELSE
    !   IF (icanopy > 1) THEN  ! R: necessary for elemental subroutine as B4 from layer above comes in ...
    !     B4 = B4_layer_above                 !    with the second level use level above...
    !   ENDIF
    ! END IF

    ! Compute "fractional vegetation cover" per grid cell
    ! (not the same as fractional cover in JSBACH!)
    IF (lai < LaiLimit) THEN
      FC = lai/LaiLimit * FcMax
    ELSE
      FC = FcMax
    END IF
    FC = MAX(FC, FcMin)

    ! COMPUTE ABSORBED PAR PER LEAF AREA (!) AND lai PARTITIONING

     !--------------------------------------------------------------------
     ! The Absorbed Par per laef Area which is used later for the Net Assimilation,
     ! is calculated via the two stream approximation of Sellers (1985):
     !  muq * d(Rdown)/dl + (1-(1-b)*omega)*Rdown - omega*b*Rup   = omega*muq*k*(1-b0)*R
     ! -muq * d(Rup)/dl   + (1-(1-b)*omega)*Rup   - omega*b*Rdown = omega*muq*k*b0*R
     !  with
     !   Rdown - downwards diffusive flux
     !   Rup   - upwards diffusive Flux
     !   R     - direct flux, with R(0) = dPAR * RPAR with RPAR incoming irridiance in PAR
     !           and R = R0 * EXP(-kl) with R0=R(0) - exponential extinction law of
     !           Monsi and Racki (1953)
     !           or similar Lambert-Beer's law
     !   b     - forward scatter fraction of diffusive flux
     !   b0    - forward scatter fraction of direct flux
     !   k     - extinction coefficient for the direct flux
     !   muq = int(1/K(mu))dmu|_0^1 - the integral of 1/K over the downward hemisphere
     !   omega - the single leaf albedo
     !
     !  The general solutions are (kl,hl=k*l,h*l):
     !  Rup   = q2*R0*EXP(-kl) + p1*B1*EXP(hl) + p2*B2*EXP(-hl)
     !  Rdown =-q1*R0*EXP(-kl) +    B1*EXP(hl) +    B2*EXP(-hl)
     !   with
     !    h  = sqrt( (1-(1-b)*omega)^2/muq^2 - omega^2*b^2/muq^2 )
     !    p1 = ( (1-(1-b)*omega) + muq * h )/omega/b
     !    p2 = ( (1-(1-b)*omega) - muq * h )/omega/b
     !   -q1 = (omega^2*b*muq*k* (1-b0) + omega*    b0  *muq*k*(1-(1-b)*omega - muq*k))/
     !         ((1-(1-b)*omega)^2-muq^2*k^2-omega^2*b^2)
     !    q2 = (omega^2*b*muq*k*    b0  + omega* (1-b0) *muq*k*(1-(1-b)*omega + muq*k))/
     !         ((1-(1-b)*omega)^2-muq^2*k^2-omega^2*b^2)
     !    B1/B2 from boundary conditions
     !-------------------------------------------------------------------
     !  Make two assumptions:
     !  1) the distribution of leaf angles is isotropic
     !  2) the leaf reflectivity and transmissivity are equal (the sum = omega)
     !  => b=0.5, b0=0.5, k=0.5/mu with mu=cos(theta) solar zenith angle => muq=1
     !
     !  => k  = 1/2/mu
     !     h  = sqrt( 1 - omega )
     !     p1 = ( 1-omega/2 + h )/omega/2
     !     p2 = ( 1-omega/2 - h )/omega/2
     !   ! p2 = 1 / p1 !
     !     q1 = ( (1 + 2*mu)*omega/2 )/(1-4*mu^2*(1-omega)) = ( k*(k + 1)*omega/2 )/
     !                                                        (k^2-1-omega)
     !     q2 = ( (1 - 2*mu)*omega/2 )/(1-4*mu^2*(1-omega)) = ( k*(k - 1)*omega/2 )/
     !                                                        (k^2-1-omega)
     !
     ! Determine B1 and B2 from the boundary conditions:
     !  1) Rdown(0) equals the incoming diffuse radiation
     !     Rdown(0) = (1-dPAR)*RPAR
     !  => Rdown(0) + R(0) = (1-dPAR)*RPAR + dPAR*RPAR = RPAR as total incoming PAR
     !  2) the reflection at the lower boundary of the canopy is given by the soil
     !     reflectance
     !     Rup(LAI) = soil_reflectivity_parPAR * (R(LAI) + Rdown(LAI))
     !  Here: FfaPAR_clL gets soil_reflectivity_parPAR as Variable soil_reflectivity_par, LAI is the total canopy LAI, lai
     !
     !  => B1 = + ( eta*R0 - (Rd+q1*R0) * gamma2 )/(gamma1 - gamma2)
     !     B2 = - ( eta*R0 - (Rd+q1*R0) * gamma1 )/(gamma1 - gamma2)
     !   with
     !     eta    = soil_reflectivity_parPAR * (1-q1)-q2) * EXP(-k*LAI)
     !     gamma1 = ( p1 - soil_reflectivity_parPAR) * EXP( + h*LAI)
     !     gamma2 = ( p2 - soil_reflectivity_parPAR) * EXP( - h*LAI)
     !     Rd     = Rdown(0) = (1-dPAR)*RPAR
     !------------------------------------------------------------------
     ! THAT IS THE COMPLETE SOLUTION OF THE TWO STREAM APPROXIMATION UNDER THE BOUNDARY
     ! CONDITIONS AND ASSUMPTIONS MENTIONED ABOVE !!!!!!!!!!!!!!!!!!
     !
     ! Therefore, the absorbed Radiation inside the canopy is:
     !   faPAR_cl = -d/dl ( R(l) + Rdown - Rup)
     !        = (1-q1-q2)*k*R0*EXP(-kl) - (1-p1)*h*B1*EXP(hl) + (1-p2)*h*B2*EXP(-hl)
     ! But the absorbed PAR per canopy layer in discrete steps is:
     !   faPAR_cl = 1/(D-D(iGP-1)) * int(-d/dl(R(l) + Rdown - Rup))dl|_D^D(iGP-1)
     !            = (R(D(iGP-1)) + Rdown(D(iGP-1)) - Rup(D(iGP-1)) - R(D) + Rdown(D)
     !              - Rup(D) / ((D-D(iGP-1))
     !  and  R(l)+Rdown-Rup = (1-q1-q2)*  R0*EXP(-kl) + (1-p1)*  B1*EXP(hl) + (1-p2) *
     !                        B2*EXP(-hl)
     !------------------------------------------------------------------
     ! The clumping of the vegetation is taken into account, defining LAIc = LAI / fc
     ! as an effective LAI.
     ! Taken this into account, l = l/fc but the solutions stay as they are because
     ! of the differentiations are take d/dl according to the NEW l=l/fc
     ! Only faPAR_cl has to be multiplied with fc at the END because D-D(iGP-1) is still
     ! the old l
     !------------------------------------------------------------------

     IF (cos_zenith_angle .GE. ZenithMinPar) THEN ! do not calculate PAR for zenith angles below ZenithMinPar (0.001)
        !----------------------------------
        !  h = sqrt( 1 - omega )
        !----------------------------------
        ZH = SQRT (1._wp - OMEGA)
        !----------------------------------
        !  p1 = ( 1-omega/2 + h )/omega/2
        !----------------------------------
        ZP1 = (1._wp - OMEGA / 2._wp + ZH) &
              / OMEGA * 2._wp
        !----------------------------------------
        ! p2 = ( 1-omega/2 - h )/omega/2 = 1 / p1
        !----------------------------------------
        ZP0 = 1._wp / ZP1
        !----------------------------------
        ! k = 0.5/mu
        !----------------------------------
        K0 = 0.5_wp / cos_zenith_angle
        IF (K0 == ZH)  K0 = K0 + 1.E-12_wp
        IF (K0 == -ZH) K0 = K0 + 1.E-12_wp
        !--------------------------------------------------------------
        ! denominator of q1 and q2
        !--------------------------------------------------------------
        X0 = (1._wp - 4._wp * cos_zenith_angle**2 * ZH**2)
        !--------------------------------------------------------------
        ! q1 = ( (1 + 2*mu)*omega/2 )/(1-4*mu^2*(1-omega))
        !    = ( k*(k + 1)*omega/2 )/(k^2-1-omega)
        !--------------------------------------------------------------
        Q1 = ((1._wp + 2._wp * cos_zenith_angle) &
                   * OMEGA / 2._wp) / X0
        !--------------------------------------------------------------
        ! q2 = ( (1 - 2*mu)*omega/2 )/(1-4*mu^2*(1-omega))
        !    = ( k*(k - 1)*omega/2 )/(k^2-1-omega)
        !--------------------------------------------------------------
        Q0 = ((1._wp - 2._wp * cos_zenith_angle) * OMEGA / 2._wp) / X0

        FSHD = MAX (FC, FcMin)
        !-----------------------------------------------------------
        ! EXP(-k*LAI/fc)
        !-----------------------------------------------------------
        EKL = EXP(-K0 / FSHD * lai)
        !------------------------------------------------------------
        ! EXP(-h*LAI/fc)
        !-----------------------------------------------------------
        EHL = EXP(-ZH / FSHD * lai)
        !-----------------------------------------------------------
        ! gamma1 = ( p1 - soil_reflectivity_parPAR) * EXP( + h*LAI)
        !-----------------------------------------------------------
        X1 = (ZP1 - soil_reflectivity_par) / EHL
        !-----------------------------------------------------------
        ! gamma2 = ( p2 - soil_reflectivity_parPAR) * EXP( - h*LAI)
        !-----------------------------------------------------------
        X0 = (ZP0 - soil_reflectivity_par) * EHL
        !-----------------------------------------------------------
        ! eta = (soil_reflectivity_parPAR * (1-q1)-q2) * EXP(-k*LAI)
        !-----------------------------------------------------------
        X2 = (soil_reflectivity_par * (1._wp - Q1) - Q0) * EKL
        !------------------------------------------------------------
        ! F = 1 - dPAR + dPAR * q1
        ! => F*RPAR = Rd + q1*R0
        ! i.e. calculation takes RPAR=1
        !-----------------------------------------------------------
        F = 1._wp - fract_par_direct + Q1 * fract_par_direct
        !-------------------------------------------------------------
        ! B1*RPAR = B1, B2*RPAR = B2, B4*RPAR = R(0) + Rdown(0) - Rup(0)
        !  B1 = + ( eta*R0 - (Rd+q1*R0) * gamma2 )/(gamma1 - gamma2)
        !------------------------------------------------------------
        B1 = (X2 * fract_par_direct - F * X0) / (X1 - X0)
        !-----------------------------------------------------------
        !  B2 = - ( eta*R0 - (Rd+q1*R0) * gamma1 )/(gamma1 - gamma2)
        !     = eta*R0/(gamma2-gamma1) + (Rd+q1*R0) / (1-gamma2/gamma1)
        !  Note: the second form is used to avoid compiler-dependent ambiguity in the
        !        result of gamma1/(gamma1-gamma2) if gamma1 = infinity
        !-----------------------------------------------------------
        B0 = X2 * fract_par_direct / (X0 - X1) + F / (1.0_wp - X0/X1)
        !-----------------------------------------------------------
        !  R(l)+Rdown-Rup = (1-q1-q2)* R0*EXP(-kl) + (1-p1)* B1*EXP(hl) + (1-p2)* B2*EXP(-hl)
        !----------------------------------------------------------------------------
        B4 =   (1._wp - Q0 - Q1) * fract_par_direct &
          &  + (1._wp - ZP1) * B1 + (1._wp - ZP0) * B0

        IF (B4_layer_above < 7777776._wp) THEN  ! R: necessary for elemental subroutine as B4 from layer above comes in ...
          B4 = B4_layer_above                  !    with the second level use level above...
        ENDIF

        IF (icanopy < ncanopy) THEN
           !----------------------------------------------------------------------------
           ! p2=1/p1
           !----------------------------------------------------------------------------
           ZP0 = 1._wp / ZP1
          !----------------------------------------------------------------------------
           ! EXP(-k*l/fc)
           ! with l=LAI*canopy_bound_lai, i.e. l is element of [0,LAI], i.e. 0,LAI/3,2LAI/3,LAI
           !----------------------------------------------------------------------------
           EKL0 =  EXP(-K0 / FSHD * canopy_bound_lai * lai)
           !----------------------------------------------------------------------------
           ! EXP(-h*l/fc)
           !----------------------------------------------------------------------------
           EHL0 = EXP(-ZH / FSHD * canopy_bound_lai * lai)
           !----------------------------------------------------------------------------
           ! R(D)+ Rdown(D)- Rup(D)=
           !      (1-q1-q2)*R0*EXP(-kl)+(1-p1)*B1*EXP(hl)+(1-p2)*B2*EXP(-hl)
           !  i.e. X0*RPAR = above
           !----------------------------------------------------------------------------
           X0 = (1._wp - Q0 - Q1) * EKL0 * fract_par_direct &
                + (1._wp - ZP1) * B1 / EHL0 &
                + (1._wp - ZP0) * B0 * EHL0
           !----------------------------------------------------------------
           ! faPAR_cl = (R(D(iGP-1))+Rdown(D(iGP-1))-Rup(D(iGP-1)) - R(D)+
           ! Rdown(D)-Rup(D) / ((D-D(iGP-1))
           ! Here faPAR_cl only the nominator; the division is made outside FfaPAR_clL
           !---------------------------------------------------------------
           faPAR_cl = B4 - X0
           !---------------------------------------------------------------
           ! Partition evenly LAI, lai_cl=LAI/3
           !          lai_cl = (canopy_bound_lai - canopy_bound_lai(icanopy-1)) * lai
           ! R(D(iGP-1))+Rdown(D(iGP-1))-Rup(D(iGP-1) in next step
           !---------------------------------------------------------------
           B4_layer_above = X0
        ELSE
           !---------------------------------------------------------------
           ! Now the same for the last layer
           !---------------------------------------------------------------

           !-------------------------------------------------------------
           ! R(D(NL))+ Rdown(D(NL))- Rup(D(NL))=
           !      (1-q1-q2)*R0*EXP(-kLAI)+(1-p1)*B1*EXP(hLAI)+(1-p2)*B2*EXP(-hLAI)
           !  i.e. X0*RPAR = above for the lower boundary
           !------------------------------------------------------------
           X0 = (1._wp - Q0 - Q1) * EKL * fract_par_direct &
                + (1._wp - ZP1) * B1 / EHL &
                + (1._wp - ZP0) * B0 * EHL          ! these variables do not change with canopy level
           !------------------------------------------------------------------
           ! faPAR_cl(NL) = (R(D(NL-1))+Rdown(D(NL-1))-Rup(D(NL-1)) - R(D(NL))+
           ! Rdown(D(NL))-Rup(D(NL)) / ((D(NL)-D(NL-1))
           ! Here faPAR_cl only the nominator; the division is made outside FfaPAR_clL
           !------------------------------------------------------------------
           faPAR_cl = B4_layer_above - X0
           !B4_layer_above = X0
        ENDIF ! icanopy
        !------------------------------------------------------------
        ! Multiplication of faPAR_cl with fc because D-D(iGP-1)
        ! division is made outside FfaPAR_clL and is still the old l ???
        !------------------------------------------------------------
        faPAR_cl = faPAR_cl * FSHD
     ENDIF  ! cos_zenith_angle

    ! Compute absorbed PAR per leaf area in canopy layer [units: (absorbed photons) / (m^2(leaf area) s)] from
    ! par and fraction of absorbed PAR
    apar_per_lai_cl = par_down_mol * faPAR_cl / (MAX(lai_cl, 1.e-10_wp))

  END SUBROUTINE calc_par

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Compute non-dimensional snow age
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_snow_age( &
    & dtime,                & ! in
    & temperature_surface,  & ! in
    & snow,                 & ! in
    & fract_snow_soil,      & ! in
    & snow_age              & ! inout
    & )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants,  ONLY: tmelt     ! tmelt =273.15_wp

    REAL(wp), INTENT(in)    :: dtime,               & !< model time step [s]
      &                        temperature_surface, & !< surface temperature [K]
      &                        snow,                & !< snowfall [kg m-2 s-1]
      &                        fract_snow_soil        !< Fraction of snow covered ground [-]

    REAL(wp), INTENT(inout) :: snow_age               !< non-dimensional age of snow

    REAL(wp)                :: hlp1                   !< temporary variable

    ! ---------------------------
    ! Update snow age where ever there is snow on the ground
    ! Snow age increases with time and temperatures closer to the melting point
    ! It decreases with new snowfall or resets to zero in case all snow is melted

    IF (fract_snow_soil > EPSILON(1.0_wp)) THEN
      hlp1     = MIN(1._wp, EXP(5000._wp * ((1._wp / tmelt) - (1._wp / temperature_surface))))
      snow_age = MAX(0.0_wp, (snow_age +                          &
        & ((hlp1**10._wp + hlp1 + 0.3_wp) * dtime * 1.e-06_wp)) * &
        & (1._wp - (snow * dtime * 0.5_wp)))
    ELSE
      snow_age = 0._wp
    END IF

  END SUBROUTINE calc_snow_age

  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_alb_lwtr( &
    & l_day,                               &
    & alb_lwtr                             &
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants, ONLY: AlbedoLakeWater

    ! Arguments
    LOGICAL,  INTENT(in)    :: &
      & l_day
    REAL(wp), INTENT(inout) :: &
      & alb_lwtr

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! @todo: use method from mo_surface_ocean:update_albedo_ocean from ECHAM6
    alb_lwtr = AlbedoLakeWater

  END SUBROUTINE calc_alb_lwtr
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_alb_lice( &
    & l_day,         &
    & t,             &
    & weq_snow_lice, &
    & alb_lice       &
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants,          ONLY: AlbedoLakeIceMin, AlbedoLakeIceMax, AlbedoLakeSnowMin, AlbedoLakeSnowMax
    USE mo_jsb_physical_constants, ONLY: tmelt

    ! Arguments
    LOGICAL,  INTENT(in)    :: &
      & l_day
    REAL(wp), INTENT(in)    :: &
      & t,                     &
      & weq_snow_lice
    REAL(wp), INTENT(inout) :: &
      & alb_lice

    ! Local variables
    REAL(wp) :: t_upper_limit    ! Upper temperature limit for cold snow albedo
    REAL(wp) :: alb_min, alb_max ! Minimum and maximum snow albedo

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    t_upper_limit = tmelt - 1._wp

    IF (weq_snow_lice > 0.01_wp) THEN
      alb_min = AlbedoLakeSnowMin
      alb_max = AlbedoLakeSnowMax
    ELSE
      alb_min = AlbedoLakeIceMin
      alb_max = AlbedoLakeIceMax
    END IF

    ! Temperature-dependent snow albedo
    IF (t >= tmelt) THEN
      alb_lice = alb_min
    ELSE IF (t < t_upper_limit) THEN
      alb_lice = alb_max
    ELSE
      alb_lice = alb_min + &
        &      ( ((alb_max - alb_min) / (tmelt - t_upper_limit)) &  ! Change of snow albedo per deg C
        &        * (tmelt - t))
    END IF

  END SUBROUTINE calc_alb_lice
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_glacier_albedo( &
    & l_day,   &
    & t,       &
    & alb_vis, &
    & alb_nir  &
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants, ONLY: AlbedoGlacierVisMin,  & ! Albedo of glacier in the visible range at the melting point
                                AlbedoGlacierVisMax,  & ! Albedo of glacier in the visible range at hard frost
                                TempAlbedoGlacierMax, & ! Maximum glacier albedo at this temperature below melting point of H2O
                                AlbedoGlacierNirMin,  & ! Albedo of glacier in the NIR range at at the melting point
                                AlbedoGlacierNirMax     ! Albedo of glacier in the NIR range at hard frost
    USE mo_jsb_physical_constants, ONLY: tmelt

    ! Arguments
    LOGICAL,  INTENT(in)    :: &
      & l_day
    REAL(wp), INTENT(in)    :: &
      & t
    REAL(wp), INTENT(inout) :: &
      & alb_vis,               & ! Albedo of canopy in the VIS range without snow cover.
      & alb_nir                  ! Albedo of canopy in the NIR range without snow cover.

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! calculate new glacier albedo (visible and NIR) only at grid points with solar radiation (i.e. during daytime)
    IF (t >= tmelt) THEN
       alb_vis = AlbedoGlacierVisMin
       alb_nir = AlbedoGlacierNirMin
    ELSE IF (t < tmelt - TempAlbedoGlacierMax) THEN
       alb_vis = AlbedoGlacierVisMax
       alb_nir = AlbedoGlacierNirMax
    ELSE
       alb_vis = AlbedoGlacierVisMin +                                  &
            (tmelt - t) * (AlbedoGlacierVisMax - AlbedoGlacierVisMin) / &
            TempAlbedoGlacierMax
       alb_nir = AlbedoGlacierNirMin +                                  &
            (tmelt - t) * (AlbedoGlacierNirMax - AlbedoGlacierNirMin) / &
             TempAlbedoGlacierMax
    END IF

  END SUBROUTINE calc_glacier_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
  ! Calculate the albedo of the soil surface
  ! Note, this subroutine should not be called for any soil or vegetation parent tile but only for tiles of the lowest level

  ! In JSBACH3 in the subroutine update_soil_albedo calculated the albedo of soils. The scheme accumulates all carbon from
  ! boxC_litter_green_ag and boxC_slow (canopy) of all tiles and builds a whole grid cell average. (This is done in subroutine
  ! update_albedo. Note, the glacier tile implementation of the scheme is arbitary as JSBACH3 could never handle a glacier tile.)
  ! Then: From boxC_slow an albedo reduction may be calculated for the soil (linear or log).
  ! From c_litter_green_ag/c_ag_sum_1 the albedo of litter may be calculted (linear).
  ! Then soil and litter albedos may be summed up according to the litter view fraction.
  ! However, the JSBACH3 scheme is not useable for the new JSBACH4 concept of child tiles. To keep the JSBACH3 concept in JSBACH4
  ! the albedo calculation could be done only on the highest parent tile (the grid cell). The childs would not have a meaningful
  ! albedo. Thus in JSBACH4 the albedo calculations are done only on the lowest child tiles and are then accumulated on the parent
  ! tiles. This represents a "well separated" assumption in contrast of the "well mixed" assumtion in JSBACH3. In the result it
  ! does not matter for linear relations between carbon (in and on the soil) and soil albedo. (Exception: the litter for any
  ! vegetation tile in JSBACH4 would cover more area than available for this tile and the litter of any other vegetation tile
  ! covers less than its area available.) Only for the logarithmic case (albedo_options%UseSOC .eq. 'log') the soil albedo result
  ! of JSBACH4 differs from JSBACH3. In this case this "well separated" result is as reasonable as the "well mixed" result of JS3.

#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_soil_albedo( &
    & l_day,                     &
    & use_alb_mineralsoil_const, &
    & use_alb_soil_organic_int,  &
    & use_alb_soil_litter,       &
    & litter_carbon,             &
    & soil_carbon,               &
    & SpecificLeafArea_C,        &
    & AlbedoLitterVIS,           &
    & AlbedoLitterNIR,           &
    & alb_background,            &
    & alb_vis_mineralsoil,       & ! InOut
    & alb_nir_mineralsoil,       & ! InOut
    & alb_vis_soil,              & ! InOut
    & alb_nir_soil               & ! InOut
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants, ONLY: AlbedoSoilConstVis, AlbedoSoilConstNir ! Albedo constant alternatively usable for
                                                                       ! alb_vis_mineralsoil

    ! Arguments
    LOGICAL,          INTENT(in)  :: l_day

    LOGICAL,          INTENT(in)  :: use_alb_mineralsoil_const
    INTEGER,          INTENT(in)  :: use_alb_soil_organic_int
    LOGICAL,          INTENT(in)  :: use_alb_soil_litter
    REAL(wp),         INTENT(in)  :: &
      & litter_carbon,               &  ! Amount of litter carbon [mol(C)/m2(canopy)], used with use_alb_soil_litter
      & soil_carbon                     ! Amount of soil carbon [mol(C)/m2(canopy)], used with use_alb_soil_organic_C
    REAL(wp),         INTENT(in)  :: &
      & SpecificLeafArea_C, &
      & AlbedoLitterVIS,    &
      & AlbedoLitterNIR
    REAL(wp),         INTENT(in)  :: &
      & alb_background
    REAL(wp),         INTENT(inout) :: &
      & alb_vis_mineralsoil,   & ! Albedo of soil without carbon or litter.
                                 ! "out" because it should be changed if use_alb_mineralsoil_const =.TRUE.
      & alb_nir_mineralsoil,   & ! Albedo of soil without carbon or litter.
      & alb_vis_soil,          & ! Albedo of soil (may incl. litter layer and/or carbon within soil)
      & alb_nir_soil             ! Albedo of soil (may incl. litter layer and/or carbon within soil)

    ! Local variables
    REAL(wp), PARAMETER ::                                   &
                           Scalfact_SOC_to_alb_vis = 3.e-4 , & ! Factor to scale soil organic content (SOC)
                                                                      ! to a albedo reduction for the vis range
                           Scalfact_SOC_to_alb_nir = 3.e-4 , & ! See above, for the nir range
                           Soil_C_min              = 500.      !

    REAL(wp)             ::                         &
                            alb_vis_soil_nolitter,  & ! Albedo of soil for the vis range (incl. the carbon within soil but
                                                      ! without the litter layer above)
                            alb_nir_soil_nolitter,  & ! See above, for the nir range
                            litter_area,            &
                            litter_view_fract,      &
                            alb_vis_litter,         &
                            alb_nir_litter
    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! use constant mineral soil albedos?
    IF (use_alb_mineralsoil_const) THEN ! all background soil albedo values are set to the same value
       alb_vis_mineralsoil = AlbedoSoilConstVis
       alb_nir_mineralsoil = AlbedoSoilConstNir
    ELSE
       alb_vis_mineralsoil = alb_background ! R: This is Rainers invention to mitigate the "error" in
                                            !    the JSBACH3 version. But still it should not be used.
                                            !    In JSBACH3 it was filled with albedo_soil_vis_base,
                                            !    which included soil and litter carbon. alb_background
                                            !    is the albedo without snow or canopy, BUT it still
       alb_nir_mineralsoil = alb_background !    includes carbon. However, this all is done to
                                            !    implement the possibility to use mineral soil albedos
                                            !    which would then be consitent with the soil albedo parametrisations!
    END IF

    ! albedo of soil with carbon content but without litter
    ! R: In JSBACH3 albedo_desert_dry_vis is either set to (1) AlbedoSoilConstVis or to albedo_soil_vis_base, which is (dependent
    !    on the jsbach.nc file) either (2) the mineral soil albedo or (3) the soil albedo with soil C and litter. Thus reading in
    !    the "wrong" kind of albedo from the jsbach.nc file wasn't so unlikely.
    !    However, the variable was also used inconsistent in the scheme and senseless for glacier.

    IF(use_alb_soil_organic_int .eq. 1) THEN ! linear
      alb_vis_soil_nolitter = alb_vis_mineralsoil - Scalfact_SOC_to_alb_vis * MIN(soil_carbon,Soil_C_min)
      alb_nir_soil_nolitter = alb_nir_mineralsoil - Scalfact_SOC_to_alb_nir * MIN(soil_carbon,Soil_C_min)
    ELSE IF(use_alb_soil_organic_int .eq. 2) THEN ! log
      alb_vis_soil_nolitter = MAX(6._wp * alb_vis_mineralsoil / (6._wp + log(1._wp + soil_carbon)),0.1_wp)
      alb_nir_soil_nolitter = MAX(6._wp * alb_nir_mineralsoil / (6._wp + log(1._wp + soil_carbon)),0.2_wp)
    ELSE
      alb_vis_soil_nolitter = alb_vis_mineralsoil ! R: This is Rainers approximation, because in JSBACH3
                                                  !    the result was set to albedo_desert_dry_vis,
                                                  !    which could include different kinds of albedos
                                                  !    (see comment above) but in each case not the
                                                  !    correct one.

      alb_nir_soil_nolitter = alb_nir_mineralsoil
    END IF ! use_alb_soil_organic_C

    ! calculate the dependence of albedo on soil moisture (to be implemented later); Note also not existent in JSBACH3

    ! albedo of litter on soil
    IF(use_alb_soil_litter) THEN
      litter_area = litter_carbon * specificLeafArea_C
      alb_vis_litter = litter_area * AlbedoLitterVIS
      alb_nir_litter = litter_area * AlbedoLitterNIR

      IF (litter_area <= EPSILON(1._wp)) THEN
        alb_vis_litter = 0.1_wp
        alb_nir_litter = 0.2_wp
      END IF

      ! calculate the litter_view_fraction, i.e. what fraction of the box is covered by litter
      ! this determines how much soil is seen through litter
      litter_view_fract = 1.0_wp - EXP(-0.5_wp * litter_area) ! R: because this graph is not linear (logarithmical) and in JSBACH4
                                                              !    the litter_areas are not summed up the result in JS4 will be
                                                              !    smaller here as it is calcul. for each single tile but the sum
                                                              !    of tiles will give a much higher "litter_view_fract" or - to
                                                              !    name it better - a higher contribution of the litter part to
                                                              !    alb_soil_vis However, therefore it meight be necessary to
                                                              !    rescale this funktion ...

      alb_vis_soil = (1.0_wp - litter_view_fract) * alb_vis_soil_nolitter + litter_view_fract * alb_vis_litter
      alb_nir_soil = (1.0_wp - litter_view_fract) * alb_nir_soil_nolitter + litter_view_fract * alb_nir_litter
    ELSE
       alb_vis_soil = alb_vis_soil_nolitter
       alb_nir_soil = alb_nir_soil_nolitter
    END IF

  END SUBROUTINE calc_soil_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_snow_albedo( &
    & l_day,             & ! in
    & fract_snow_soil,   & ! in
    & snow_age,          & ! in
    & cos_zenith_angle,  & ! in
    & albedo_age_weight, & ! in
    & t_srf,             & ! in
    & alb_vis_snow,      & ! inout
    & alb_nir_snow       & ! inout
    & )

    !$ACC ROUTINE SEQ

    ! Declarations
    USE mo_rad_constants, ONLY: &
      & ZenithAngleFactor,      & ! Factor in solar zenith angle dependence of snow albedo
                                  ! (the increase of snow albedo is the higher this factor)
      & AlbedoSnowVisMax_age,   &
      & AlbedoSnowNirMax_age,   &
      & AlbedoSnowVisAge,       & ! Maximal rel. reduction of snow albedo by aging in the visible range
      & AlbedoSnowNirAge,       & ! Maximal rel. reduction of snow albedo by aging in the NIR range
      & AlbedoSnowAngle,        & ! Maximal rel. reduction of snow absorption by large solar zenith angle
      & TempAlbedoSnowMax,      & ! Maximum snow albedo at this temperature below melting point of H2O
      & AlbedoSnowVisMin_temp,  &
      & AlbedoSnowNirMin_temp,  &
      & AlbedoSnowVisMax_temp,  &
      & AlbedoSnowNirMax_temp

    USE mo_jsb_physical_constants,  ONLY: tmelt     ! tmelt =273.15_wp

    ! Arguments
    LOGICAL,  INTENT(in)  :: l_day
    REAL(wp), INTENT(in)  :: fract_snow_soil,   & !< Fraction of snow covered ground
                             snow_age,          & !< non-dimensional age of snow
                             cos_zenith_angle,  & !< cosinus of zenith angle of incoming radiation
                             t_srf

    REAL(wp), INTENT(in)  ::                    & !
                             albedo_age_weight    ! 0 < albedo_age_weight < 1: snow albedo is calcualted by linearly weighting
                                                  !                            the snow albedo resulting from both age and
                                                  !                            temperature schemes
    REAL(wp), INTENT(inout) ::              &
                             alb_vis_snow,  & ! Albedo of canopy in the VIS range without snow cover.
                             alb_nir_snow     ! Albedo of canopy in the NIR range without snow cover.

    ! Local variables
    REAL(wp)              ::                    &
                             alb_vis_snow_age,  &
                             alb_nir_snow_age,  &
                             alb_vis_snow_temp, &
                             alb_nir_snow_temp, &
                             snow_age_factor,   &
                             snow_albedo_angle_factor

    ! Initialize local variables
    alb_vis_snow_age  = 0.0_wp
    alb_nir_snow_age  = 0.0_wp
    alb_vis_snow_temp = 0.0_wp
    alb_nir_snow_temp = 0.0_wp
    alb_vis_snow      = 0.0_wp
    alb_nir_snow      = 0.0_wp

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! calculate albedo of snow dependent on age (snow age scheme, from BATS model of NCAR)
    IF (fract_snow_soil > EPSILON(1._wp)) THEN
       snow_age_factor = snow_age / (1._wp + snow_age)
       IF (cos_zenith_angle < 0.5_wp) THEN
         snow_albedo_angle_factor = ((1.0_wp + ZenithAngleFactor)/ &
                                    (1.0_wp + 2.0_wp * ZenithAngleFactor * cos_zenith_angle) - 1.0_wp) / &
                                  ZenithAngleFactor
       ELSE
          snow_albedo_angle_factor = 0._wp
       END IF
       alb_vis_snow_age = AlbedoSnowVisMax_age * (1._wp - AlbedoSnowVisAge * snow_age_factor)
       alb_vis_snow_age = alb_vis_snow_age + AlbedoSnowAngle * snow_albedo_angle_factor * &
                                                   (1._wp - alb_vis_snow_age)
       alb_nir_snow_age = AlbedoSnowNirMax_age * (1._wp - AlbedoSnowNirAge * snow_age_factor)
       alb_nir_snow_age = alb_nir_snow_age + AlbedoSnowAngle * snow_albedo_angle_factor * &
                                                   (1._wp - alb_nir_snow_age)
    ELSE
       alb_vis_snow_age = 0._wp
       alb_nir_snow_age = 0._wp
    END IF


    ! calculate albedo of snow dependent on temperature (snow temperature scheme, from echam5 model of MPI-M)
    IF (t_srf >= tmelt) THEN
          alb_vis_snow_temp = AlbedoSnowVisMin_temp
          alb_nir_snow_temp = AlbedoSnowNirMin_temp
    ELSE IF (t_srf < tmelt - TempAlbedoSnowMax) THEN
          alb_vis_snow_temp = AlbedoSnowVisMax_temp
          alb_nir_snow_temp = AlbedoSnowNirMax_temp
    ELSE
          alb_vis_snow_temp = AlbedoSnowVisMin_temp +     &
               (tmelt - t_srf) * (AlbedoSnowVisMax_temp - &
               AlbedoSnowVisMin_temp) / TempAlbedoSnowMax
          alb_nir_snow_temp = AlbedoSnowNirMin_temp +     &
               (tmelt - t_srf) * (AlbedoSnowNirMax_temp - &
               AlbedoSnowNirMin_temp) / TempAlbedoSnowMax
    END IF

    ! weight age vs. temperature of snow albedo to receive alb_vis_snow and alb_nir_snow
    alb_vis_snow = albedo_age_weight * alb_vis_snow_age + &
                              (1._wp - albedo_age_weight) * alb_vis_snow_temp
    alb_nir_snow = albedo_age_weight * alb_nir_snow_age + &
                              (1._wp - albedo_age_weight) * alb_nir_snow_temp

  END SUBROUTINE calc_snow_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_pond_albedo( &
    & l_day,             & ! in
    & cos_zenith_angle,  & ! in
    & t,                 & ! in
    & wtr_pond,          & ! in
    & ice_pond,          & ! in
    & alb_pond_vis,      & ! inout
    & alb_pond_nir       & ! inout
    & )

    !$ACC ROUTINE SEQ

    ! Declarations
    USE mo_rad_constants,           ONLY: AlbedoLakeIceMin, AlbedoLakeIceMax

    USE mo_jsb_physical_constants,  ONLY: tmelt     ! tmelt =273.15_wp

    ! Arguments
    LOGICAL,  INTENT(in)  :: l_day
    REAL(wp), INTENT(in)  :: cos_zenith_angle,  & ! cosinus of zenith angle of incoming radiation
                             t,                 &
                             wtr_pond,          &
                             ice_pond

    REAL(wp), INTENT(inout) :: alb_pond_vis,    & ! Albedo of ponds in the VIS range without snow cover
                               alb_pond_nir       ! Albedo of ponds in the NIR range without snow cover

    ! Local variables

    REAL(wp) :: t_upper_limit    ! Upper temperature limit for cold snow albedo

    REAL(wp)              :: alb_min,           &
                             alb_max,           &
                             alb_wtr_pond_vis,  &
                             alb_wtr_pond_nir,  &
                             alb_ice_pond_vis,  &
                             alb_ice_pond_nir,  &
                             fact_ice

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    !!!! Albedo formulation for unfrozen ponds is taken from CLM
    alb_wtr_pond_vis = 0.05_wp /(MAX(0.001_wp,cos_zenith_angle) + 0.15_wp)
    alb_wtr_pond_nir = alb_pond_vis


    ! Temperature-dependent ice albedo (taken from lake ice fomulation)
    t_upper_limit     = tmelt - 1._wp
    IF (t >= tmelt) THEN
      alb_ice_pond_vis = AlbedoLakeIceMin
    ELSE IF (t < t_upper_limit) THEN
      alb_ice_pond_vis = AlbedoLakeIceMax
    ELSE
      alb_ice_pond_vis = AlbedoLakeIceMin + &
        &      ( ((AlbedoLakeIceMax - AlbedoLakeIceMin) / (tmelt - t_upper_limit)) &  ! Change of snow albedo per deg C
        &        * (tmelt - t))
    END IF
    ! Albedo nir from CLM
    alb_ice_pond_nir = 0.4_wp

    ! Agregation of albedo of ice and wtr under the assumption that ice is predominantly
    ! near the surfce hence has a stronger influence on albedo
    IF (ice_pond > 0._wp) THEN
      fact_ice     = MAX(0._wp, MIN(1._wp, 1._wp - (wtr_pond/(wtr_pond+ice_pond))**2._wp))
      alb_pond_vis = fact_ice * alb_ice_pond_vis + (1.0_wp - fact_ice) * alb_wtr_pond_vis
      alb_pond_nir = fact_ice * alb_ice_pond_nir + (1.0_wp - fact_ice) * alb_wtr_pond_nir
    ELSE
      alb_pond_vis = alb_wtr_pond_vis
      alb_pond_nir = alb_wtr_pond_nir
    END IF

  END SUBROUTINE calc_pond_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE calc_sky_view_fractions( &
    & l_day,              & ! in
    & lai,                & ! in
    & fract_fpc_max,      & ! in
    & cos_zenith_angle,   & ! in
    & StemArea,           & ! in
    & sky_view_fract,     & ! inout
    & sky_view_fract_stem & ! inout
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants, ONLY: SkyViewFactor ! Constant in calculating the sky view fractions

    ! Arguments
    LOGICAL,  INTENT(in)  :: l_day
    REAL(wp), INTENT(in)  :: lai,              & ! leaf area index
                             fract_fpc_max,    & ! maximum foliage projected cover
                             cos_zenith_angle    ! cosinus of zenith angle of incoming radiationa
    REAL(wp), INTENT(in)  :: StemArea
    REAL(wp), INTENT(inout) :: &
      & sky_view_fract,   & ! Fraction of "visible" bare ground below vegetation without accounting for stems.
                            ! (Falling with rising LAI.)
      & sky_view_fract_stem   ! Fraction of soil below vegetation which is covered by stems but not by LAI.

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! Fraction of "visible" bare ground below vegetation regarding only the LAI coverage
    sky_view_fract = 1._wp - fract_fpc_max * (1._wp - EXP(- SkyViewFactor * lai))

    ! Fraction of soil below vegetation covered only by stems
    ! R: later: skip "sky_view_fract -" => changes meaning of sky_view_fract_stem but later shortens calculations
    sky_view_fract_stem = &
      & sky_view_fract - (1._wp - fract_fpc_max * (1._wp - EXP(- SkyViewFactor * &
                                                     (lai + StemArea) * (2.0_wp - cos_zenith_angle))))

  END SUBROUTINE calc_sky_view_fractions
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE Merge_albedos_of_vegtile( &
    & l_day,               & ! in
    & l_forest,            & ! in
    & sky_view_fract,      & ! in
    & sky_view_fract_stem, & ! in
    & AlbedoCanopySnow,    & ! in
    & fract_snow_soil,     & ! in
    & fract_snow_can,      & ! in
    & fract_pond,          & ! in
    & alb_vis_snow,        & ! in
    & alb_nir_snow,        & ! in
    & alb_vis_soil,        & ! in
    & alb_nir_soil,        & ! in
    & alb_vis_can,         & ! in
    & alb_nir_can,         & ! in
    & alb_pond_vis,        & ! in
    & alb_pond_nir,        & ! in
    & alb_vis,             & ! inout
    & alb_nir              & ! inout
    & )

    !$ACC ROUTINE SEQ

    ! Arguments
    LOGICAL,  INTENT(in) :: l_day ! flag to indicate if it is day or night
    LOGICAL,  INTENT(in) :: l_forest

    REAL(wp), INTENT(in) :: &
      sky_view_fract,       &
      sky_view_fract_stem
    REAL(wp), INTENT(in) :: &
      AlbedoCanopySnow
    REAL(wp), INTENT(in) :: &
      fract_snow_soil,      & ! Snow covered area fraction of soil
      fract_snow_can,       & ! Fraction of snow covered canopy (forest)
      fract_pond,           & ! Fraction of ponds
      alb_vis_snow,         & ! Note, this snow albedo is used only for snow on soil.
                              ! The albedo for snow on canopy (AlbedoCanopySnow) is not the same.
      alb_nir_snow,         & ! See comment above
      alb_vis_soil,         &
      alb_nir_soil,         &
      alb_vis_can,          & ! Albedo of canopy in the VIS range without snow cover.
      alb_nir_can,          & ! Albedo of canopy in the NIR range without snow cover.
      alb_pond_vis,         & ! Albedo of ponds in the VIS range without snow cover.
      alb_pond_nir            ! Albedo of ponds in the NIR range without snow cover.

    REAL(wp), INTENT(inout) :: &
      & alb_vis, &
      & alb_nir

    ! Local variables
    REAL(wp) :: &
      & alb_vis_nostems, & ! albedo vis of surface without considering the snow cover
                           ! when it is assumed that no stems (but only the canopy) cover the soil
      & alb_nir_nostems, &
      & alb_vis_grnd,    &
      & alb_nir_grnd

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    alb_vis_grnd    = fract_pond           * alb_pond_vis &
                    + (1._wp - fract_pond) * alb_vis_soil
    alb_nir_grnd    = fract_pond           * alb_pond_nir &
                    + (1._wp - fract_pond) * alb_nir_soil

    alb_vis_nostems = (1.0_wp - sky_view_fract) * alb_vis_can + sky_view_fract * alb_vis_grnd
    alb_nir_nostems = (1.0_wp - sky_view_fract) * alb_nir_can + sky_view_fract * alb_nir_grnd

    ! albedo of forests = weighted mean of albedo of ground below canopy and albedo of canopy
    IF (l_forest) THEN
      alb_vis = MAX( &
          ! Albedo of the (snow covered and uncoverd) visible soil:
        & (sky_view_fract - sky_view_fract_stem)                                                  &
        &   * (fract_snow_soil * alb_vis_snow        + (1.0_wp - fract_snow_soil) * alb_vis_grnd) &
        & +                                                                                       &
          ! Albedo of the visible canopy (snow covered and uncoverd) but only stems:
          ! alb_vis_grnd represents here the stem albedo (consistent with the input datasets
          ! where the albedo soil represents the ground albedo for LAI=0 which includes stems)
        & sky_view_fract_stem                                                                     &
        &   * (fract_snow_can * AlbedoCanopySnow     + (1.0_wp - fract_snow_can) * alb_vis_grnd)  &
        & +                                                                                       &
          ! Albedo of the visible canopy (snow covered and uncoverd) without stems:
        & (1.0_wp - sky_view_fract)                                                               &
        &   * (fract_snow_can * AlbedoCanopySnow     + (1.0_wp - fract_snow_can) * alb_vis_can)   &
          !
          ! OR use the albedo without regarding stems if then the albedo is higher:
        & , alb_vis_nostems)

      alb_nir =  MAX( &
          ! Albedo of the (snow covered and uncoverd) visible soil:
        & (sky_view_fract - sky_view_fract_stem)                                                   &
        &   * (fract_snow_soil * alb_nir_snow        +  (1.0_wp - fract_snow_soil) * alb_nir_grnd) &
        & +                                                                                        &
          ! Albedo of the visible canopy (snow covered and uncoverd) but only stems:
        & sky_view_fract_stem                                                                      &
        &   * (fract_snow_can * AlbedoCanopySnow     +  (1.0_wp - fract_snow_can) * alb_nir_grnd)  &
        & +                                                                                        &
          ! Albedo of the visible canopy (snow covered and uncoverd) without stems:
        & (1.0_wp - sky_view_fract)                                                                &
        &   * (fract_snow_can * AlbedoCanopySnow     +  (1.0_wp - fract_snow_can) * alb_nir_can)   &
          !
          ! OR use the albedo without regarding stems if then the albedo is higher:
        & , alb_nir_nostems)
    ELSE
      alb_vis = MAX( &
        ! Albedo of the (snow covered and uncoverd) visible soil:
        & fract_snow_soil * alb_vis_snow         + (1.0_wp - fract_snow_soil) * alb_vis_nostems &
        &      , alb_vis_nostems)
       alb_nir = MAX( &
        ! Albedo of the (snow covered and uncoverd) visible soil:
        & fract_snow_soil * alb_nir_snow         + (1.0_wp - fract_snow_soil) * alb_nir_nostems &
        &      , alb_nir_nostems)
    END IF ! l_forest

  END SUBROUTINE Merge_albedos_of_vegtile
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
  ! This is the echam5 scheme
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  REAL(wp) FUNCTION get_surface_albedo_simple( &
    & t_srf,             &
    & fract_forest,      &
    & grnd_snow_fract,   & ! R: fract_snow_soil umbenannt, da grnd_snow_fract nicht unbedingt dasselbe sein muss...
    & fract_snow_can,    &
    & fract_pond,        &
    & lai,               &
    & alb_background     &
    & )

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    ! Use declarations
    USE mo_jsb_physical_constants,  ONLY: tmelt
    USE mo_rad_constants,           ONLY: AlbedoCanopySnow_simple, AlbedoSnowMin_simple, AlbedoSnowMax_simple, &
      &                                   AlbedoLakeWater, AlbedoLakeIceMin

    ! Arguments
    REAL(wp), INTENT(in) ::                  &
                            t_srf,           &
                            fract_forest,    &
                            grnd_snow_fract, &
                            fract_snow_can,  &
                            fract_pond,      &
                            lai,             &
                            alb_background

    ! Local variables
    REAL(wp) ::                                   &
                            sky_view_fract,       &  ! R: how much of light reaches the ground through the canopy
                            forest_view_fract,    &  ! R: how much of light is catched only by forest
                            albedo_snow_grnd,     &  ! Albedo of snow on ground
                            min_temp_snow_albedo, &  ! Temperature below which albedo of snow is at maximum
                            albedo_grnd,          &  ! Albedo on ground
                            albedo_can,           &  ! Albedo on canopy
                            albedo_pond              ! Albedo of ponds


    ! ---------------------------
    ! Go

    ! Calculate albedo_grnd incl. snowcover
    min_temp_snow_albedo = tmelt - 5._wp

    IF (t_srf >= tmelt) THEN
      albedo_snow_grnd = AlbedoSnowMin_simple
      albedo_pond      = AlbedoLakeWater
    ELSE IF (t_srf < min_temp_snow_albedo) THEN
      albedo_snow_grnd = AlbedoSnowMax_simple
      albedo_pond      = AlbedoLakeIceMin
    ELSE
      albedo_snow_grnd = AlbedoSnowMin_simple &
                         + (tmelt - t_srf) * (AlbedoSnowMax_simple - AlbedoSnowMin_simple) / (tmelt - min_temp_snow_albedo)
      albedo_pond      = AlbedoLakeIceMin &
                         - (tmelt - t_srf) * (AlbedoLakeIceMin - AlbedoLakeWater) / (tmelt - min_temp_snow_albedo)
    END IF

    ! Calclulate albedo of ground incl. snowcover
    albedo_grnd = grnd_snow_fract                                   * albedo_snow_grnd &
                + (1._wp - grnd_snow_fract) * fract_pond            * albedo_pond      &
                + (1._wp - grnd_snow_fract) * (1.0_wp - fract_pond) * alb_background

    ! Calculate albedo_can incl. snowcover
    albedo_can = fract_snow_can * AlbedoCanopySnow_simple + (1._wp - fract_snow_can) * alb_background

    ! Calculate albedo of surface
    sky_view_fract = EXP(-1.0_wp * MAX(lai, 2._wp)) ! The value 1.0 is the SkyViewFactor.
                                                    ! The value 2.0 represents the stem_area.
                                                    ! Keep these values constant:
                                                    ! The function get_surface_albedo_simple_echam5 is not thought as an actual
                                                    ! alternative albedo scheme but to provide comparability to old ECHAM versions

    forest_view_fract = fract_forest * (1._wp - sky_view_fract)


!   R: falls man Glaciermaske verwenden will, siehe auch radiation interface...
!     If (glaciermaske) THEN
!        get_surface_albedo_simple = albedo_snow_grnd
!     ELSE
    get_surface_albedo_simple = MAX(alb_background, &
                                    (1._wp - forest_view_fract) * albedo_grnd &
                                    + forest_view_fract * albedo_can)

  END FUNCTION get_surface_albedo_simple
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
  FUNCTION Has_minimal_radiation_2(swvis, swnir) RESULT(l_day)

    !$ACC ROUTINE SEQ

    USE mo_rad_constants, ONLY: MinRadiation ! Minimum amount of radiation for albedo calculations

    REAL(wp), INTENT(in) :: &
      & swvis, &
      & swnir
    LOGICAL :: l_day

    l_day = swvis + swnir > MinRadiation

  END FUNCTION Has_minimal_radiation_2
  !
  FUNCTION Has_minimal_radiation_1(sw) RESULT(l_day)

    !$ACC ROUTINE SEQ

    USE mo_rad_constants, ONLY: MinRadiation ! Minimum amount of radiation for albedo calculations

    REAL(wp), INTENT(in) :: &
      & sw
    LOGICAL :: l_day

    l_day = sw > MinRadiation

  END FUNCTION Has_minimal_radiation_1

  ! ================================================================================================================== !
  !>
  !> Calculate the sums of below ground C pools - relevant for soil albedo schemes
  !!
  FUNCTION get_soil_carbon(c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1,  &
    &                      c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2,  &
    &                      c_humus_1, c_humus_2)     RESULT(soil_carbon)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(in) :: &
      & c_acid_bg1,         &
      & c_water_bg1,        &
      & c_ethanol_bg1,      &
      & c_nonsoluble_bg1,   &
      & c_acid_bg2,         &
      & c_water_bg2,        &
      & c_ethanol_bg2,      &
      & c_nonsoluble_bg2,   &
      & c_humus_1,          &
      & c_humus_2

    REAL(wp) :: soil_carbon  ! Amount of soil carbon relevant for soil albedo

    soil_carbon =  c_acid_bg1 + c_water_bg1 + c_ethanol_bg1 + c_nonsoluble_bg1  &
      &          + c_acid_bg2 + c_water_bg2 + c_ethanol_bg2 + c_nonsoluble_bg2  &
      &          + c_humus_1  + c_humus_2

  END FUNCTION get_soil_carbon

  ! ================================================================================================================== !
  !>
  !> Calculate the sums of above ground green litter and green humus C pools - relevant for soil albedo schemes
  !!
  FUNCTION get_litter_carbon(c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1,  &
    &                        c_humus_1)     RESULT(ag_green_litter_carbon)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(in) :: &
      & c_acid_ag1,         &
      & c_water_ag1,        &
      & c_ethanol_ag1,      &
      & c_nonsoluble_ag1,   &
      & c_humus_1

    REAL(wp) :: ag_green_litter_carbon  ! Amount of litter carbon relevant for soil albedo

    ag_green_litter_carbon =  c_acid_ag1 + c_water_ag1 + c_ethanol_ag1 + c_nonsoluble_ag1  &
      &                     + c_humus_1

  END FUNCTION get_litter_carbon


  ! ======================================================================================================= !
  !>
  !> quincy model
  !>
  ! ======================================================================================================= !

#ifndef __NO_QUINCY__
  !-----------------------------------------------------------------------------------------------------
  !> calculate shortwave radiation budget of the canopy, given soil / snow albedo and incident radiation
  !!
  !! @par info on routine
  !!  This routine has been adapted from ADFs implementation of the radiation calculations
  !!  by Spitters, C J T. 1986. \n Separating the Diffuse and Direct Component of Global Radiation
  !!  and Its Implications for Modeling Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis.
  !!  Agricultural and Forest Meteorology 38 (1-3): 231-42. \n
  !!  Equations numbers refer to the Spitters paper \n
  !!  - The routine has been extended to diagnose albedo in the visible and NIR range. The implicit assumption
  !!  is made that albedo in the PAR and VIS range are identical \n
  !!  - The routine has been adjusted to a two-stream approximation scheme by adding the upward flux
  !!  from soil reflection \n
  !!  - The routine has been adjusted to account for canopy clumping (Campbell & Norman 1998) \n
  !!  - The routine has been adjusted to account for light absorption and reflection on the non-
  !!    green stem area
  !!  - The output has been modified in that arad_shaded|sunlit_cl is still light on canopy
  !!    but arad_can is the energy adsorbed by leaves+stem to properly calculate energy fluxes,
  !!    albedo and fapar
  !!
  !! @par Spitters Eq. 12 is wrong
  !!  accordingly to Wang2003, "radf4" was replaced by "-kbl" when calculating "idrdra"
  !!
  !! @par canopy clumping factor
  !!  to disable clumping set lctlib(:)%omega_clumping = 1.0_wp
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_canopy_radiation_budget( nc     , &   ! in
                              ncanopy             , &   ! in
                              lctlib_sigma        , &   ! in
                              lctlib_omega_clumping     , &   ! in
                              lctlib_crown_shape_factor , &   ! in
                              vis_nir_id          , &   ! in
                              sw_srf_down         , &   ! in
                              fract_par_diffuse   , &   ! in
                              cos_zenith_angle    , &   ! in
                              lai                 , &   ! in
                              sai                 , &   ! in
                              lai_cl              , &   ! in
                              cumm_lai_cl         , &   ! in
                              alb_surface         , &   ! in
                              appfd_sunlit_cl     , &   ! out
                              appfd_shaded_cl     , &   ! out
                              fleaf_sunlit_cl     , &   ! out
                              fapar               , &   ! out
                              alb_can             , &   ! out
                              alb_lnd               )   ! out

    USE mo_rad_parameters,          ONLY: min_cos_zenith_angle_rad, kbl0_vis, kdf0_vis, kbl0_nir, kdf0_nir, rho2sbeta, &
      &                                   k_csf, albedo_stem_vis, albedo_stem_nir, kbl_stem, rad2ppfd
    USE mo_jsb_math_constants,      ONLY: eps8
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER                                , INTENT(in)  :: nc                 !< dimensions
    INTEGER                                , INTENT(in)  :: ncanopy            !< number of canopy layers
    REAL(wp)                               , INTENT(in)  :: lctlib_sigma       !< lctlib parameter
    REAL(wp)                               , INTENT(in)  :: lctlib_omega_clumping     !< lctlib parameter
    REAL(wp)                               , INTENT(in)  :: lctlib_crown_shape_factor !< lctlib parameter
    CHARACTER(len=3)                       , INTENT(in)  :: vis_nir_id         !< vis nir identifier for this function only
    REAL(wp), DIMENSION(:)                 , INTENT(in)  :: sw_srf_down        !< shortwave downward radiation flux (W m-2)
    REAL(wp), DIMENSION(:)                 , INTENT(in)  :: fract_par_diffuse
      !< fraction of diffuse radiation (--). Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
    REAL(wp), DIMENSION(:)                 , INTENT(in)  :: cos_zenith_angle   !< solar angle (radians)
    REAL(wp), DIMENSION(:)                 , INTENT(in)  :: lai                !< leaf area index (m-2 m-2)
    REAL(wp), DIMENSION(:)                 , INTENT(in)  :: sai                !< stem area index (m-2 m-2)
    REAL(wp), DIMENSION(:,:)               , INTENT(in)  :: lai_cl             !< leaf area index of canopy layer (m-2 m-2)
    REAL(wp), DIMENSION(:,:)               , INTENT(in)  :: cumm_lai_cl        !< cummulative leaf area index above centre of canopy layer
    REAL(wp), DIMENSION(:)                 , INTENT(in)  :: alb_surface        !< merged soil and snow albedo
    REAL(wp), DIMENSION(:,:), OPTIONAL     , INTENT(out) :: fleaf_sunlit_cl    !< fraction of leaves that is sunlit (--)
    REAL(wp), DIMENSION(:,:), OPTIONAL     , INTENT(out) :: appfd_sunlit_cl    !< absorbed radiation on sunlit leaves (mu mol m-2 s-1)
    REAL(wp), DIMENSION(:,:), OPTIONAL     , INTENT(out) :: appfd_shaded_cl    !< absorbed radiation of shaded leaves (mu mol m-2 s-1)
    REAL(wp), DIMENSION(:),   OPTIONAL     , INTENT(out) :: fapar              !< fraction of absorbed radiation by the canopy+stem area (--)
    REAL(wp), DIMENSION(:),   OPTIONAL     , INTENT(out) :: alb_can            !< albedo of the canopy (including green and stem area)
    REAL(wp), DIMENSION(:),   OPTIONAL     , INTENT(out) :: alb_lnd            !< albedo of the tile (including soil albedo)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                 :: tmp_cos_zenith_angle            ! local variable of solar angle (radians)
    REAL(wp), DIMENSION(nc)                 :: kdf, kbl                        ! extinction coefficients for radiation over black canopy
    REAL(wp), DIMENSION(nc)                 :: irad_dir_down, &                ! downward direct and diffuse radiation
                                               irad_dif_down
    REAL(wp), DIMENSION(nc)                 :: irad_dir_down_bc, &             ! light intensity at the bottom of the canopy
                                               irad_dif_down_bc
    REAL(wp), DIMENSION(nc)                 :: irad_dif_up                     ! upward diffuse shortwave radiation
    REAL(wp), DIMENSION(nc)                 :: arad_can                        ! absorbed radiation by the canopy (W m-2)
    REAL(wp), DIMENSION(nc)                 :: arad_soil                       ! absorbed radiation by the soil (W m-2)
    REAL(wp), DIMENSION(nc)                 :: arad_can_up                     ! absorbed radiation by the canopy from the upward stream (W m-2)
    REAL(wp), DIMENSION(nc,ncanopy)         :: arad_can_shaded_up_cl           ! idem for canopy+stem layers
    REAL(wp), DIMENSION(nc,ncanopy)         :: arad_can_shaded_cl              ! absorbed shaded radiation for canopy+stem layers
    REAL(wp), DIMENSION(nc,ncanopy)         :: arad_can_sunlit_cl              ! absorbed sunlit radiation for canopy+stem layers
    REAL(wp), DIMENSION(nc,ncanopy)         :: arad_shaded_cl                  ! absorbed shaded radiation for canopy layers
    REAL(wp), DIMENSION(nc,ncanopy)         :: arad_sunlit_cl                  ! absorbed sunlit radiation for canopy layers
    REAL(wp), DIMENSION(nc,ncanopy)         :: fsunlit_cl                      ! fraction of sunlit canopy layer
    REAL(wp), DIMENSION(nc,ncanopy)         :: sai_cl                          ! stem area index per canopy layer
    REAL(wp), DIMENSION(nc,ncanopy)         :: cumm_sai_cl                     ! cummulative stem area index per canopy layer
    REAL(wp), DIMENSION(nc)                 :: plai                            ! leaf+stem area index
    REAL(wp), DIMENSION(nc,ncanopy)         :: plai_cl                         ! leaf+stem area index per canopy layer
    REAL(wp), DIMENSION(nc,ncanopy)         :: cumm_plai_cl                    ! cummulative leaf+stem area index per canopy layer
    REAL(wp), DIMENSION(nc)                 :: radf1, radf2, radf3             ! dummy variables
    REAL(wp), DIMENSION(nc)                 :: radf4, radf5                    ! dummy variables
    REAL(wp), DIMENSION(nc)                 :: radf1_leaf, radf2_leaf          ! dummy variables
    REAL(wp), DIMENSION(nc)                 :: radf3_leaf, radf4_leaf          ! dummy variables
    REAL(wp)                                :: alb_stem                        ! stem albedo (--)
    REAL(wp)                                :: temp1, temp2, radf3b            ! dummy variables
    REAL(wp)                                :: kbl0, kdf0                      ! hlp for local rad_constants variables, differentiating vis/nir
    REAL(wp)                                :: idfa                            ! diffuse radiation penetrating foliage at lc
    REAL(wp)                                :: idra                            ! direct radiation penetrating sunlit foliage at lc
    REAL(wp)                                :: idrdra                          ! direct radiation remaining direct at lc
    REAL(wp)                                :: omega                           ! clumping factor
    REAL(wp)                                :: rho_leaf                        ! fraction reflected by green canopy (--)
    REAL(wp)                                :: rho                             ! fraction reflected by canopy+stems (--)
    REAL(wp)                                :: kdf_leaf, kbl_leaf              ! extinction coefficients for radiation over black canopy
    INTEGER                                 :: ic                              ! iteration over grid points
    INTEGER                                 :: icanopy                         ! iteration over canopy layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_canopy_radiation_budget'

    !$ACC DATA CREATE(arad_can, arad_can_up, arad_soil) &
    !$ACC   CREATE(arad_sunlit_cl(1:nc,1:ncanopy), arad_shaded_cl(1:nc,1:ncanopy), arad_can_sunlit_cl(1:nc,1:ncanopy)) &
    !$ACC   CREATE(arad_can_shaded_cl(1:nc,1:ncanopy),arad_can_shaded_up_cl(1:nc,1:ncanopy),sai_cl(1:nc,1:ncanopy)) &
    !$ACC   CREATE(cumm_sai_cl(1:nc,1:ncanopy),fsunlit_cl(1:nc,1:ncanopy), cumm_plai_cl(1:nc,1:ncanopy)) &
    !$ACC   CREATE(plai_cl(1:nc,1:ncanopy), plai, tmp_cos_zenith_angle, kdf, kbl) &
    !$ACC   CREATE(irad_dir_down, irad_dif_down, irad_dir_down_bc, irad_dif_down_bc, irad_dif_up) &
    !$ACC   CREATE(radf1, radf2, radf3, radf4, radf5, radf1_leaf, radf2_leaf, radf3_leaf, radf4_leaf) ASYNC(acc_stream)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.9 set local variables according to vis/nir function call and initialise output
    !!
    SELECT CASE(vis_nir_id)
      CASE("vis")
        kbl0         = kbl0_vis                 ! extinction coefficient over black leaves
        kdf0         = kdf0_vis                 ! extinction coefficient for diffuse radiation
        alb_stem     = albedo_stem_vis
      CASE("nir")
        kbl0         = kbl0_nir                 ! extinction coefficient over black leaves
        kdf0         = kdf0_nir                 ! extinction coefficient for diffuse radiation
        alb_stem     = albedo_stem_nir
    END SELECT

    ! init out vars
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO icanopy = 1,ncanopy
      DO ic = 1,nc
        IF (PRESENT(fleaf_sunlit_cl)) fleaf_sunlit_cl(ic,icanopy) = 0.0_wp
        IF (PRESENT(appfd_sunlit_cl)) appfd_sunlit_cl(ic,icanopy) = 0.0_wp
        IF (PRESENT(appfd_shaded_cl)) appfd_shaded_cl(ic,icanopy) = 0.0_wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1,nc
      IF (PRESENT(fapar)) fapar(ic) = 0.0_wp
      IF (PRESENT(alb_can)) alb_can(ic) = 0.0_wp
      IF (PRESENT(alb_lnd)) alb_lnd(ic) = 0.0_wp
    END DO
    !$ACC END PARALLEL LOOP

    ! init local vars
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1,nc
      arad_can(ic) = 0.0_wp
      arad_can_up(ic) = 0.0_wp
      arad_soil(ic) = 0.0_wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO icanopy = 1,ncanopy
      DO ic = 1,nc
        arad_sunlit_cl(ic,icanopy) = 0.0_wp
        arad_shaded_cl(ic,icanopy) = 0.0_wp
        arad_can_sunlit_cl(ic,icanopy) = 0.0_wp
        arad_can_shaded_cl(ic,icanopy) = 0.0_wp
        arad_can_shaded_up_cl(ic,icanopy) = 0.0_wp
        sai_cl(ic,icanopy) = 0.0_wp
        cumm_sai_cl(ic,icanopy) = 0.0_wp
        fsunlit_cl(ic,icanopy) = 0.0_wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !> 1.0 calculations
    !!
    ! calculate stem and plant area index per canopy layer
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1,nc
      IF (lai(ic) < eps8) THEN ! case of vegetation with no leaves, condense SAI to first layer
        sai_cl(ic,1)      = sai(ic)
        cumm_sai_cl(ic,1) = sai(ic) / 2.0_wp ! sai mid-point evaluation
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO icanopy = 1,ncanopy
      DO ic = 1,nc
        IF (lai(ic) >= eps8) THEN ! case of vegetation with leaves
          ! assign SAI that is given by current leaf area
          sai_cl(ic,icanopy)      = lai_cl(ic,icanopy) * sai(ic) / lai(ic)
          cumm_sai_cl(ic,icanopy) = cumm_lai_cl(ic,icanopy) * sai(ic) / lai(ic)
        END IF

        plai_cl(ic,icanopy)      = lai_cl(ic,icanopy) + sai_cl(ic,icanopy)
        cumm_plai_cl(ic,icanopy) = cumm_lai_cl(ic,icanopy) + cumm_sai_cl(ic,icanopy)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1,nc
      plai(ic) = lai(ic) + sai(ic)
    END DO
    !$ACC END PARALLEL LOOP

    !> enforce tmp_cos_zenith_angle to be equal/larger 'min_cos_zenith_angle_rad'
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      tmp_cos_zenith_angle(ic) = MAX(cos_zenith_angle(ic), min_cos_zenith_angle_rad)
    END DO
    !$ACC END PARALLEL LOOP

    !> Eq.1 constants for albedo calculations
    radf3b = (1.0_wp - lctlib_sigma)
    temp1  = SQRT(radf3b)
    temp2  = 2.0_wp * (1.0_wp - temp1) / (1.0_wp + temp1)

    !> @par IF only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(omega, rho_leaf, rho, kdf_leaf, kbl_leaf)
    DO ic = 1, nc
      IF (sw_srf_down(ic) > eps8) THEN

        !> Downwelling direct radiation (Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation)
        irad_dir_down(ic) = (1.0_wp - fract_par_diffuse(ic)) * sw_srf_down(ic)
        !> Downwelling diffuse radiation
        irad_dif_down(ic) = fract_par_diffuse(ic) * sw_srf_down(ic)

        !> @par IF lai > eps8
        IF (plai(ic) > eps8) THEN

          !> Eq.1  Fraction of radiation reflected by green closed canopy (fraction)
          rho_leaf = temp2 / (1.0_wp + rho2sbeta * tmp_cos_zenith_angle(ic))
          rho = (rho_leaf * lai(ic) + alb_stem * sai(ic)) / plai(ic)

          !> clumping factor (eq 15.35) Campbell & Norman (1998) as a function of sun-angle
          omega = lctlib_omega_clumping / (lctlib_omega_clumping + (1.0_wp - lctlib_omega_clumping) &
            &        * EXP(-k_csf * ACOS(tmp_cos_zenith_angle(ic)) ** lctlib_crown_shape_factor))

          !> Extinction coefficient for radiation over black leaves in canopy
          kdf_leaf = kdf0 * temp1 * omega ! temp1 == SQRT(1.0_wp - lctlib_sigma)
          kbl_leaf = kbl0 * kdf_leaf / (0.8_wp * temp1 * tmp_cos_zenith_angle(ic))
          kdf(ic)      = ( kdf_leaf * lai(ic) + kbl_stem * sai(ic) ) / plai(ic)
          kbl(ic)      = ( kbl_leaf * lai(ic) + kbl_stem * sai(ic) ) / plai(ic)

          !> Eq.10 & 11 Some factors for radiation calculations that change with time (for canopy and canopy+stem)
          radf1_leaf(ic) = (1.0_wp - rho_leaf) * irad_dif_down(ic) * kdf_leaf            ! eq 10
          radf2_leaf(ic) = (1.0_wp - rho_leaf) * temp1 * irad_dir_down(ic) * kbl_leaf    ! eq 11
          radf1(ic)      = (1.0_wp - rho) * irad_dif_down(ic) * kdf(ic)                      ! eq 10
          radf2(ic)      = (1.0_wp - rho) * temp1 * irad_dir_down(ic) * kbl(ic)              ! eq 11
          !> Eq.12 absorbed direct radiation
          radf3_leaf(ic) = radf3b * irad_dir_down(ic) * kbl_leaf                    ! eq 12
          radf3(ic)      = radf3b * irad_dir_down(ic) * kbl(ic)                         ! eq 12
          radf4(ic)      = -temp1 * kbl(ic)

          !> light intensity at bottom of canopy
          irad_dif_down_bc(ic) = (1.0_wp - rho) * irad_dif_down(ic) * EXP(-kdf(ic)  * plai(ic))
          irad_dir_down_bc(ic) = (1.0_wp - rho) * irad_dir_down(ic) * EXP(radf4(ic) * plai(ic))

          !> diffuse upward shortwave flux from soil reflection and absorbed radiation of the soil or snow
          irad_dif_up(ic) = alb_surface(ic) * (irad_dif_down_bc(ic) + irad_dir_down_bc(ic))
          arad_soil(ic)   = (1.0_wp - alb_surface(ic)) * (irad_dif_down_bc(ic) + irad_dir_down_bc(ic))

          !> factor for the absorption calculation of the upward stream
          radf4_leaf(ic) = (1.0_wp - rho_leaf) * irad_dif_up(ic) * kdf_leaf            ! eq 10
          radf5(ic)      = (1.0_wp - rho) * irad_dif_up(ic) * kdf(ic)                      ! eq 10

          !> upward shortwave diffuse flux reaching the top of the canopy
          irad_dif_up(ic) = (1.0_wp - rho) * irad_dif_up(ic) * EXP(-kdf(ic)*plai(ic))
        ELSE
           !> case of no PLAI: absorped radiation and upward flux from soil only
           arad_soil(ic)   = (1.0_wp - alb_surface(ic) ) * (irad_dif_down(ic) + irad_dir_down(ic))
           irad_dif_up(ic) = alb_surface(ic)             * (irad_dif_down(ic) + irad_dir_down(ic))
        END IF
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !> @par loop over canopy layers. All fluxes are evaluated at the centre of the canopy layer
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream) &
    !$ACC   PRIVATE(idfa, idra, idrdra)
    DO icanopy = 1, ncanopy
      DO ic = 1, nc
        !> @par IF leaf area index of canopy layer > eps8: Integrate absorbed radiation over actually existing canopy layers
        IF (sw_srf_down(ic) > eps8 .AND. plai_cl(ic,icanopy) > eps8) THEN
          !> @par IF radf2 <= eps8
          ! case of only diffuse radiation
          IF (radf2(ic) <= eps8) THEN
            !> Eq.10 Mean radiation penetrating foliage at cumm_plai_cl
            arad_shaded_cl(ic,icanopy)     = radf1_leaf(ic) * EXP(-kdf(ic) * cumm_plai_cl(ic,icanopy)) + &    ! eq 10
                                            radf4_leaf(ic) * EXP(-kdf(ic) * (MAX(plai(ic)-cumm_plai_cl(ic,icanopy),0.0_wp)))
            arad_can_shaded_cl(ic,icanopy) = radf1(ic) * EXP(-kdf(ic) * cumm_plai_cl(ic,icanopy)) + &         ! eq 10
                                            radf5(ic) * EXP(-kdf(ic) * (MAX(plai(ic)-cumm_plai_cl(ic,icanopy),0.0_wp)))
          !> @par ELSE Case of direct and diffuse radiation
          ELSE
            !> Sunlit foliage at cumm_plai_cl(icanopy) (fraction)
            !! kbl: Extinction coefficient for radiation over black leaves in canopy
            fsunlit_cl(ic,icanopy) = EXP(-kbl(ic) * cumm_plai_cl(ic,icanopy))

            !> first calculate absorption of radiation for ecosystem level (leaves+stem)
            !> Eq.10 Diffuse radiation penetrating foliage at cumm_plai_cl(icanopy)
            idfa = radf1(ic) * EXP(-kdf(ic) * cumm_plai_cl(ic,icanopy)) &
              &   + radf5(ic) * EXP(-kdf(ic) * (MAX(plai(ic)-cumm_plai_cl(ic,icanopy), 0.0_wp)))
            !> Direct radiation penetrating sunlit foliage at cumm_plai_cl(icanopy)
            idra = radf2(ic) * EXP(radf4(ic) * cumm_plai_cl(ic,icanopy))
            !> Eq.12 Mean direct radiation remaining direct at cumm_plai_cl(icanopy)
            !! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
            idrdra = radf3(ic) * EXP(-kbl(ic) * cumm_plai_cl(ic,icanopy))
            !> Radiation penetrating shaded foliage at cumm_plai_cl(icanopy) (umol/m2/s)
            !! Eq.13 Sum of diffuse radiation and diffuse component of direct flux
            arad_can_shaded_cl(ic,icanopy) = idfa + MAX((idra - idrdra), 0.0_wp)
            !> Radiation penetrating sunlit foliage at lc (umol/m2/s)
            !! sum of diffuse and direct radiation
            arad_can_sunlit_cl(ic,icanopy) = arad_can_shaded_cl(ic,icanopy) + radf3(ic)

            !> now calculate radiation absorption of green foliage only
            !> Eq.10 Diffuse radiation penetrating foliage at cumm_plai_cl(icanopy)
            idfa = radf1_leaf(ic) * EXP(-kdf(ic) * cumm_plai_cl(ic,icanopy)) &
              &    + radf4_leaf(ic) * EXP(-kdf(ic) * (MAX(plai(ic) - cumm_plai_cl(ic,icanopy), 0.0_wp)))
            !> Direct radiation penetrating sunlit foliage at cumm_plai_cl(icanopy)
            idra = radf2_leaf(ic) * EXP(radf4(ic) * cumm_plai_cl(ic,icanopy))
            !> Eq.12 Mean direct radiation remaining direct at cumm_plai_cl(icanopy)
            !! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
            idrdra = radf3_leaf(ic) * EXP(-kbl(ic) * cumm_plai_cl(ic,icanopy))
            !> Radiation penetrating shaded foliage at cumm_plai_cl(icanopy) (umol/m2/s)
            !! Eq.13 Sum of diffuse radiation and diffuse component of direct flux
            arad_shaded_cl(ic,icanopy) = idfa + MAX((idra - idrdra), 0.0_wp)
            !> Radiation penetrating sunlit foliage at lc (umol/m2/s)
            !! sum of diffuse and direct radiation
            arad_sunlit_cl(ic,icanopy) = arad_shaded_cl(ic,icanopy) + radf3_leaf(ic)
          END IF

          !> for albedo calculation: radiation absorbed from the upward stream
          arad_can_shaded_up_cl(ic,icanopy) = radf5(ic) * EXP(-kdf(ic) * (MAX(plai(ic) - cumm_plai_cl(ic,icanopy), 0.0_wp)))
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO icanopy = 1, ncanopy
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        IF (sw_srf_down(ic) > eps8 .AND. plai_cl(ic,icanopy) > eps8) THEN
          !> total absorbed radiation (canopy+stem)
          arad_can(ic) = arad_can(ic) + ( arad_can_shaded_cl(ic,icanopy) * (1.0_wp - fsunlit_cl(ic,icanopy)) &
            &                         + arad_can_sunlit_cl(ic,icanopy) * fsunlit_cl(ic,icanopy) )            &
            &                         * plai_cl(ic,icanopy)
          !> radiation penetrating foliage from below for albedo calculation (only needed for canopy+stem)
          arad_can_up(ic) = arad_can_up(ic) + arad_can_shaded_up_cl(ic,icanopy) * plai_cl(ic,icanopy)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !> 2.0 Diagnostics

    !> 2.1 Photosynthetically Active Photon Flux Density (ppfd) and fraction of sunlit leaves per canopy layer
    !>
    IF (PRESENT(appfd_sunlit_cl) .AND. PRESENT(appfd_shaded_cl) .AND. PRESENT(fleaf_sunlit_cl)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO icanopy = 1, ncanopy
        DO ic = 1, nc
          appfd_sunlit_cl(ic,icanopy) = arad_sunlit_cl(ic,icanopy) * rad2ppfd
          appfd_shaded_cl(ic,icanopy) = arad_shaded_cl(ic,icanopy) * rad2ppfd
          fleaf_sunlit_cl(ic,icanopy) = fsunlit_cl(ic,icanopy)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    !> 2.2 Fraction of absorbed radiation (faPAR)
    !>
    IF (PRESENT(fapar)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF (sw_srf_down(ic) > eps8) THEN
          fapar(ic) = arad_can(ic) / sw_srf_down(ic)
        ELSE
          fapar(ic) = 0.0_wp
        END IF
      END DO
      !$ACC END PARALLEL
    END IF

    !> 2.3 Canopy and surface albedo
    !>
    !> @par IF only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    IF (PRESENT(alb_lnd) .AND. PRESENT(alb_can)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF (sw_srf_down(ic) > eps8) THEN
          !> @par IF lai > eps8
          IF (plai(ic) > eps8) THEN
            !> albedo of the canopy, as seen from the top of the canopy
            !! inferred from the difference between top and bottom of the
            !! canopy light intensities minus the amount of absorbed downward radiation
            alb_can(ic) = 1.0_wp - (irad_dif_down_bc(ic) + irad_dir_down_bc(ic) + (arad_can(ic) - arad_can_up(ic))) &
              &                   / (irad_dir_down(ic) + irad_dif_down(ic))
            !> total upward shortwave radiation flux at the top of the canopy:
            !! soil-based flux + reflectance of the canopy
            irad_dif_up(ic) = irad_dif_up(ic) &
              &              + alb_can(ic) * (irad_dir_down(ic) + irad_dif_down(ic) &
              &                               - (irad_dif_down_bc(ic) + irad_dir_down_bc(ic)))
          END IF
          !> albedo and net shortwave radiation of the tile
          alb_lnd(ic) = irad_dif_up(ic) / (irad_dir_down(ic) + irad_dif_down(ic))
        END IF
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC END DATA
  END SUBROUTINE calc_canopy_radiation_budget

#endif
#endif
END MODULE mo_rad_process
