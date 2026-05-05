! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_art_oem_vprm

  !------------------------------------------------------------------------------
  !
  ! Description:
  !   This module contains subroutines for the reading in of gridded vegetation 
  !   indices (LSWI and EVI) and computation of the biospheric fluxes for the 
  !   Online Emission Module (OEM).
  !
  !! Modifications: 
  !! 2024: Arash Hamzehloo, Empa
  !! - VPRM was refactored & ported to GPUs.
  !==============================================================================
    ! ICON
    USE mo_kind,                   ONLY: wp, i8
    USE mo_exception,              ONLY: message, message_text
    USE mo_var_list,               ONLY: t_var_list_ptr
    USE mo_model_domain,           ONLY: t_patch, p_patch
    USE mo_nonhydro_state,         ONLY: p_nh_state
    USE mo_parallel_config,        ONLY: nproma, idx_1d, blk_no,     &
                                     &   idx_no
    USE mo_math_constants,         ONLY: rad2deg
  
    USE mo_time_config,            ONLY: time_config !configure_time
    USE mtime,                     ONLY: MAX_DATETIME_STR_LEN, &
                                     &   datetimeToString,     &
                                     &   julianday,            &
                                     &   newJulianday,         &
                                     &   getJulianDayFromDatetime,  &
                                     &   getDatetimeFromJulianDay,  &
                                     &   datetime,             &
                                     &   no_of_ms_in_a_day,    &
                                     &   timedelta,            &
                                     &   newTimedelta,         &
                                     &   OPERATOR(+),          &
                                     &   OPERATOR(==)
  
    USE mo_nwp_phy_state,          ONLY: prm_diag
  
    ! ART
    USE mo_art_atmo_data,          ONLY: t_art_atmo
    USE mo_art_data,               ONLY: p_art_data
    USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
  
    ! OEM
    USE mo_art_oem_types,          ONLY: p_art_oem_data,          &
                                     &   t_art_oem_data,          &
                                     &   t_art_oem_config,        &
                                     &   t_art_oem_ensemble
  
    USE mo_oem_config,             ONLY: vprm_par,          & 
                                     &   vprm_lambda,       & 
                                     &   vprm_alpha,        & 
                                     &   vprm_beta,         & 
                                     &   vprm_tmin,         & 
                                     &   vprm_tmax,         & 
                                     &   vprm_topt,         & 
                                     &   vprm_tlow,         & 
                                     &   lcut_area,         & 
                                     &   lon_cut_start,     &
                                     &   lon_cut_end,       &
                                     &   lat_cut_start,     & 
                                     &   lat_cut_end
  USE netcdf,                      ONLY: nf90_noerr, nf90_open, nf90_strerror, &
                                     &   nf90_inq_dimid, nf90_inq_varid,       &
                                     &   nf90_get_var, nf90_inquire_dimension, &
                                     &   nf90_nowrite, nf90_close

!---------------------------------------------------------------------------

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: art_oem_compute_biosphere_fluxes, &
      &       art_oem_extract_dos
  


    ! Constant variable
    INTEGER,  PARAMETER :: tp_param_hourofday = 24
    INTEGER,  PARAMETER :: tp_param_dayofweek = 7
    INTEGER,  PARAMETER :: tp_param_monthofyear = 12
    INTEGER(KIND=2), PARAMETER :: tp_param_hour = 8784

    INTEGER :: global_iteration_v = 0

    INTEGER, ALLOCATABLE :: is_vm(:), ie_vm(:)
  
    TYPE(julianday), POINTER :: jd, jdref
  
    CHARACTER(LEN=3), DIMENSION(7) :: day_of_week = &
      & (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)
  
  
  !==============================================================================
  ! Module procedures
  !==============================================================================
  
  
  CONTAINS
  
  
  SUBROUTINE art_oem_compute_biosphere_fluxes(p_tracer_now,p_patch,dtime,mtime_current,ierror,yerrmsg)
  
  !-----------------------------------------------------------------------------
  ! Description: This subroutine computes the VPRM flux field and adds them to the OEM-tracer
  !-----------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    REAL(wp), INTENT(inout) :: &
     &  p_tracer_now(:,:,:,:)     !< tracer mixing ratio
  
    TYPE(t_patch), INTENT(IN) :: &
     &  p_patch                   !< patch on which computation is performed
  
    REAL(wp), INTENT(IN) :: &
     &  dtime                     !< time step    
  
    TYPE(datetime), POINTER ::  &
     &  mtime_current             !< current datetime
  
    INTEGER, INTENT(OUT)             :: ierror
    CHARACTER(LEN= *), INTENT(OUT)   :: yerrmsg
  
  
    !---------------------------------------------------------------------------
    ! Local variables
    INTEGER :: dos, is, ie, i_startblk, i_endblk, &
      &        jb, jc, jg, k, l, nr, nt, nblks_c, &
      &        table_nr, trcr_idx, i, ens_count, nt_old, & 
      &        nens, ens_int_idx, nc, nlev 
  
    REAL(KIND=wp) :: z_mass, t_degc, veg_frac, a1, a2, a3, &
      &              t_scale, w_scale, p_scale, gee, evi_thresh, pabs, lambda, newflux, new_ens_flux
  
    REAL(KIND=wp), PARAMETER :: eps_div = 1.e-12_wp
  
    CHARACTER(LEN=2) :: numstring
  
    TYPE(datetime) :: datetime_next
  
    TYPE(timedelta), POINTER :: mtime_td
  
  
    CHARACTER(*), PARAMETER :: routine = "art_oem_compute_biosphere_fluxes"

    INTEGER,  PARAMETER :: veg_class = 8

    REAL(KIND=wp), DIMENSION(nproma, 8) :: temp_flux
  
  
  !- End of header
  !==============================================================================
 
      ierror= 0

      jg   = p_patch%id
      nblks_c = p_art_data(jg)%atmo%nblks
      CALL art_oem_extract_dos(time_config%tc_current_date,dos)

      i_startblk = p_art_data(jg)%atmo%i_startblk
      i_endblk   = p_art_data(jg)%atmo%i_endblk

      nlev = p_art_data(jg)%atmo%nlev

      global_iteration_v = global_iteration_v + 1

      IF (global_iteration_v==1) THEN 
        ALLOCATE(is_vm(i_endblk))
        ALLOCATE(ie_vm(i_endblk))

        DO jb = i_startblk, i_endblk

          CALL art_get_indices_c(jg, jb, is, ie)
          is_vm(jb) = is
          ie_vm(jb) = ie

        END DO

        !$ACC ENTER DATA COPYIN(is_vm, ie_vm)

      ENDIF

      !$ACC DATA COPYIN(dos) CREATE(temp_flux, newflux, new_ens_flux)
      IF (p_art_oem_data%configure%vprm_tracer>0) THEN
        nt_old = 0
        ens_count = 1
        

        DO jb = i_startblk, i_endblk

          is = is_vm(jb) 
          ie = ie_vm(jb)

          DO nt = 1, p_art_oem_data%configure%vprm_tracer
            trcr_idx = p_art_oem_data%configure%vprm_idx(nt)
            !------------------------------------------------------------------------------
            ! Section 1: Compute the respiration fluxes
            !------------------------------------------------------------------------------ 
            !$ACC PARALLEL DEFAULT(PRESENT) 
            IF (p_art_oem_data%configure%vprm_flux_type(nt) == 'resp' ) THEN
              !$ACC LOOP SEQ 
              DO l = 1, veg_class
                !$ACC LOOP GANG VECTOR PRIVATE(t_degc, veg_frac) 
                DO jc = is, ie
                  ! Get 2-meter temperature in degC
                  t_degc = p_art_data(jg)%atmo%t_2m(jc,jb) - 273.15_wp
  
                  veg_frac = p_art_oem_data%data_fields%vprm_lu_class_fraction(jc,jb,l)

                  temp_flux(jc, l) = (vprm_alpha(l) * t_degc + vprm_beta(l)) * veg_frac * 44.01_wp * 1.e-9_wp

                  ! Check if the flux is negative or NaN (the latter is important for class 8 - Urban in the case of GPP)
                  IF (.NOT. (temp_flux(jc, l)>0)) THEN
                    temp_flux(jc, l) = 0._wp
                  ENDIF
                END DO
              ENDDO
            END IF ! p_art_oem_data%configure%vprm_flux_type(nt) == 'resp'

            !------------------------------------------------------------------------------
            ! Section 2: Compute the gpp fluxes
            !------------------------------------------------------------------------------

            IF (p_art_oem_data%configure%vprm_flux_type(nt) == 'gpp') THEN
              p_scale = 0.0_wp
              w_scale = 0.0_wp
              ! Radiation
              !$ACC LOOP SEQ
              DO l = 1, veg_class 

                !$ACC LOOP GANG VECTOR PRIVATE(w_scale, p_scale, evi_thresh,veg_frac, pabs, a1, a2, a3, t_degc, t_scale, w_scale, p_scale )
                DO jc = is, ie

                  t_degc = p_art_data(jg)%atmo%t_2m(jc,jb) - 273.15_wp

                  a1 = t_degc - vprm_tmin(l)
                  a2 = t_degc - vprm_tmax(l)
                  a3 = t_degc - vprm_topt(l)
                  IF (a1 < 0._wp) CYCLE ! No photosynthesis

                  ! Temperature sensitivity on photosynthesis
                  t_scale = (a1 * a2 / (a1 * a2 - a3**2 + eps_div))

                  w_scale = (1._wp + p_art_oem_data%data_fields%lswi(jc,jb,dos)) / (1._wp + p_art_oem_data%data_fields%lswi_max(jc,jb))
      
                  ! Effect of leaf phenology
                  IF (l == p_art_oem_data%data_fields%i_vprm_lc_evergreen) THEN          ! Evergreen
                    p_scale = 1._wp
                    
                  ELSE                                                ! Other vegetation types
                    evi_thresh = p_art_oem_data%data_fields%evi_min(jc,jb) + &
                      &          0.55_wp * (p_art_oem_data%data_fields%evi_max(jc,jb) - p_art_oem_data%data_fields%evi_min(jc,jb))
                    IF (p_art_oem_data%data_fields%evi(jc,jb,dos) >= evi_thresh) THEN  ! Full canopy period
                      p_scale = 1._wp
                    ELSE
                      p_scale = (1._wp + p_art_oem_data%data_fields%lswi(jc,jb,dos)) / 2._wp  ! Bad-burst to full canopy period                  
                    END IF
                  END IF
                  ! VPRM vegetation class fraction
                  veg_frac = p_art_oem_data%data_fields%vprm_lu_class_fraction(jc,jb,l)
                  ! Shortwave downward photosynthetically active flux at the surface [W/m2]
                  pabs = (prm_diag(jg)%swflxsfc(jc,jb) + prm_diag(jg)%swflx_up_sfc(jc,jb))/ 0.505_wp
                  temp_flux(jc, l) = -1.0_wp * (vprm_lambda(l) * t_scale * p_scale * w_scale *  &
                    &          p_art_oem_data%data_fields%evi(jc,jb,dos) * 1._wp / (1._wp + pabs/vprm_par(l)) * pabs) &
                    &          * veg_frac * 44.01_wp * 1.e-9_wp
                  ! Check if the flux is negative or NaN (the latter is important for class 8 - Urban in the case of GPP)
                  IF (.NOT. (temp_flux(jc, l)>0)) THEN
                    temp_flux(jc, l) = 0._wp
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            !$ACC END PARALLEL
            lambda = 1._wp

            !------------------------------------------------------------------------------
            ! Section 3: Compute the emission fields for ensemble members
            !------------------------------------------------------------------------------

            IF ( ANY( p_art_oem_data%ensemble%ens_name==p_art_oem_data%configure%vprm_name(nt) ) .OR. &
              &  ANY( p_art_oem_data%ensemble%vprm_bg_ens)) THEN
                 
              !$ACC PARALLEL DEFAULT(PRESENT)
              !$ACC LOOP GANG COLLAPSE(2) PRIVATE(ens_int_idx, jc)
              DO nens = 1,SIZE(p_art_oem_data%data_fields%lambda_mat, dim=4) 
                DO table_nr=1,p_art_oem_data%configure%ens_tracer            
                  IF ((p_art_oem_data%ensemble%ens_table(1,table_nr)==nens .AND. &
                    &  p_art_oem_data%ensemble%ens_name(table_nr)==p_art_oem_data%configure%vprm_name(nt)) .OR.            &
                    & (p_art_oem_data%ensemble%ens_table(1,table_nr)==nens .AND. &
                    &  p_art_oem_data%ensemble%vprm_bg_ens(table_nr) .AND. (                &
                    &  p_art_oem_data%configure%vprm_name(nt) == p_art_oem_data%ensemble%vprm_bg_ens_name(2,table_nr) .OR. &
                    &  p_art_oem_data%configure%vprm_name(nt) == p_art_oem_data%ensemble%vprm_bg_ens_name(1,table_nr))))   &
                  & THEN   

                    ens_int_idx = p_art_oem_data%ensemble%ens_table(2,table_nr) 
                    !$ACC LOOP VECTOR PRIVATE(nr, lambda, l, nc, new_ens_flux)     
                    DO jc = is, ie           
 
                      new_ens_flux = 0.0_wp
                      nr = p_art_oem_data%data_fields%reg_map(jc,jb)

                        ! Scale newfluxes with their corresponding category lambdas
                        !$ACC LOOP SEQ 
                        DO l=1, veg_class
                          IF (p_art_oem_data%configure%vprm_flux_type(nt) == 'resp') THEN
                            nc = p_art_oem_data%data_fields%lambda_categories_ids(l)
                          ELSE
                            nc = p_art_oem_data%data_fields%lambda_categories_ids(l + veg_class)
                          END IF
                          lambda = p_art_oem_data%data_fields%lambda_mat(ens_count,nc,nr,nens)
                            
                          new_ens_flux = new_ens_flux + temp_flux(jc, l) * lambda 

                        ENDDO ! l, number of vprm classes

                      ! Update concentration tendencies with GPP or RE ensembles, in the case of GPP it should be substructed 
                      IF (p_art_oem_data%ensemble%vprm_bg_ens(table_nr) .AND. &
                        &  p_art_oem_data%configure%vprm_name(nt) == p_art_oem_data%ensemble%vprm_bg_ens_name(2,table_nr)) THEN 
                          p_tracer_now(jc,nlev,jb,ens_int_idx) = p_tracer_now(jc,nlev,jb,ens_int_idx) - &
                            &                             dtime * new_ens_flux  /p_nh_state(jg)%diag%airmass_now(jc,nlev,jb)
                      ELSE
                          p_tracer_now(jc,nlev,jb,ens_int_idx) = p_tracer_now(jc,nlev,jb,ens_int_idx) + &
                            &                             dtime * new_ens_flux  /p_nh_state(jg)%diag%airmass_now(jc,nlev,jb)
                      ENDIF
                    ENDDO
                  ENDIF ! ens_table(1,table_nr)==nens
                ENDDO ! table_nr=1,p_art_oem_data%configure%ens_tracer 
              ENDDO ! nens = 1,SIZE(lambda_mat, dim=4)
              !$ACC END PARALLEL
            ENDIF !lens(nt)==.TRUE.

            !------------------------------------------------------------------------------
            ! Section 4: Compute the emission fields for non-ensemble members
            !------------------------------------------------------------------------------
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR PRIVATE(k, newflux) 
            DO jc = is, ie
              newflux = 0.0_wp 
              !$ACC LOOP SEQ 
              DO k = 1, veg_class
                newflux = newflux + temp_flux(jc, k)
              ENDDO 
              newflux = MAX(newflux, 0.0_wp)

              IF ( ANY( p_art_oem_data%configure%emis_diag_out==trcr_idx ) ) THEN  

                !$ACC LOOP SEQ
                DO k = 1, veg_class
                  p_tracer_now(jc,k,jb,trcr_idx) = temp_flux(jc, k)
                ENDDO

              ELSE   

                p_tracer_now(jc,nlev,jb,trcr_idx) = p_tracer_now(jc,nlev,jb,trcr_idx) + &
                  &                             dtime * newflux /p_nh_state(jg)%diag%airmass_now(jc,nlev,jb)

              END IF

            ENDDO ! jc = is, nproma
            !$ACC END PARALLEL
          ENDDO ! nt = 1, p_art_oem_data%configure%vprm_tracer
        ENDDO ! jb = i_startblk, i_endblk
      ENDIF ! p_art_oem_data%configure%vprm_tracer>0
      !$ACC END DATA
      CALL message(routine, "VPRM executed")

  
  
  !------------------------------------------------------------------------------
  ! End of the Subroutine
  !------------------------------------------------------------------------------
  
    END SUBROUTINE art_oem_compute_biosphere_fluxes
  
  !==============================================================================
  !==============================================================================
  !+ Extract which day of the simulation we are at [1, 2, 3, ...]
  !------------------------------------------------------------------------------
  
    SUBROUTINE art_oem_extract_dos(date, day_of_year)
  
      IMPLICIT NONE
  
      ! Parameters
      TYPE(datetime), INTENT(IN) :: date
      INTEGER, INTENT(OUT) :: day_of_year
  
      ! Local variables
      REAL(KIND=wp) :: y1, y2
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: time_string
      INTEGER :: errno, hour_of_year
      TYPE(datetime) :: refdate

      TYPE(timedelta), POINTER :: mtime_td
  
  
      CALL datetimeToString(date, time_string)
      !format: 2017-01-12T05:57:00.000
  
      ! Get current date
      jd => newJulianday(0_i8, 0_i8)
      mtime_td => newTimedelta("PT12H") ! Add 12 hours to the current date and to the reference in order to take into account Julian day starting from 12 and not from 00
      CALL getJulianDayFromDatetime(date + mtime_td,jd,errno) !(dt, jd, errno)
    
      refdate = time_config%tc_startdate
      ! Julian day of refdate:
      jdref => newJulianday(0_i8, 0_i8)
      mtime_td => newTimedelta("PT12H")
      CALL getJulianDayFromDatetime(refdate + mtime_td,jdref,errno) !(dt, jd, errno)
  
      y1 = INT(jd%day,wp)
      y2 = INT(jdref%day,wp)
      day_of_year = (y1 - y2) + 1
      ! print*,'!!!DAY12!!!',y1,y2
      ! y1 = REAL(jd%day,wp) + REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)
      ! y2 = REAL(jdref%day,wp) + REAL(jdref%ms,wp)/REAL(no_of_ms_in_a_day,wp)
      ! hour_of_year = int(24._wp*(y1-y2)) !24 * (nactday-1) + hour_of_day ! 1 Jan, 00 UTC -> 1
      ! day_of_year = hour_of_year/24_i8 + 1
  
    END SUBROUTINE art_oem_extract_dos
  !==============================================================================
  !==============================================================================
  
  END MODULE mo_art_oem_vprm
