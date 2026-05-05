!> Contains the routines for the hydro processes
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

!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"

MODULE mo_hydro_process
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_impl_constants, ONLY: WB_IGNORE, WB_LOGGING, WB_ERROR
  USE mo_jsb_control, ONLY: acc_stream

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_surface_hydrology_land, calc_surface_hydrology_glacier, calc_soil_hydrology,             &
    & get_soilhyd_properties, calc_wskin_fractions_lice, calc_wet_fractions_veg, calc_wet_fractions_bare, &
    & get_canopy_conductance, get_water_stress_factor, calc_orographic_features

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  !TODO: This parameter is only used in arno_scheme and should be defined there.
  REAL(wp), PARAMETER ::  zwdmin = 0.05_wp
  !$ACC DECLARE COPYIN(zwdmin)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'

CONTAINS

  ! ===============================================================================================================================
  !>
  !> #### Compute surface hydrology on glacier-free land
  !>
  !> Surface and pond storages are updated from precipitation, evapotranspiration and sublimation
  !> fluxes, and depending on the surface conditions and scale assumptions water infiltration into
  !> the soil and surface runoff are calculated.

#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE calc_surface_hydrology_land (                            &
    & lstart, dtime, ltpe_closed, hydro_scale, l_dynsnow,             &
    & l_infil_subzero, l_latflow_to_streamflow,                       &
    & snow_depth_max, steepness, flowlag, t_soil_sl1, t_snow_mean,    &
    & wind_10m, t_air, skinres_canopy_max, skinres_max,               &
    & weq_pond_max, fract_snow, fract_skin, fract_pond,               &
    & fract_pond_max, evapotrans, evapopot,                           &
    & transpiration, rain, snow, wpi_rootzone, wpi_rootzone_max,      &
    & wtr_soil, hyd_cond_sat_sl1, ice_impedance,                      &
    & wtr_skin, weq_snow_soil, weq_snow_can, snow_soil_dens,          &
    & wtr_pond, ice_pond, wtr_latflow_res_srf, le_pc_remain,          &
    & q_snocpymlt, snow_accum, snowmelt_soil, pond_freeze, pond_melt, &
    & evapotrans_soil, evapo_skin, evapo_snow, evapo_pond,            &
    & wtr_pond_net_flx, water_to_soil, evapo_deficit, infilt,         &
    & runoff, runoff_horton, wtr_latflow_srf)

    USE mo_jsb_math_constants,     ONLY: pi, seconds_per_day
    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf, dens_snow_min, dens_snow
    USE mo_hydro_constants,        ONLY: InterceptionEfficiency, Semi_Distributed_, Uniform_, &
      &                                  dens_snow_max, crhosmint, crhosmaxt, csnow_tmin, crhosmax_tmin
    USE mo_sse_constants,          ONLY: snow_depth_min

    LOGICAL, INTENT(in) ::       &
      & lstart,                  & !< T: Start of experiment
      & ltpe_closed,             & !< T: Terraplanet setup with closed water balance
      & l_dynsnow,               & !< T: Compute snow density dynamically
      & l_infil_subzero,         & !< T: Allow infiltration at temperatures below 0 degC
      & l_latflow_to_streamflow    !< T: Outflow of intermediary surface runoff storage goes to streamflow;
                                   !< True, unless HydroTiles are used.
    INTEGER, INTENT(in) ::       &
      & hydro_scale                !< Hydrology scale (Semi_distributed: with ARNO scheme; Uniform: e.g.
                                   !< for site level or HydroTiles)
    REAL(wp), INTENT(in) ::      &
      & dtime,                   & !< Time step length [s]
      & snow_depth_max,          & !< Maximum snow depth [m water equivalent]; -1. for no limit
      & steepness(:),            & !< Parameter representing subgrid scale slopes (for runoff calculation) []
      & flowlag(:),              & !< Lag factor accounting for retention in surface runoff []
      & t_soil_sl1(:),           & !< Temperature of the uppermost soil layer [K]
      & t_snow_mean(:),          & !< Level weighted mean snow temperature [K]
      & wind_10m(:),             & !< Wind speed at 10m height [m/s]
      & t_air(:),                & !< Lowest layer atmosphere temperature [K]
      & skinres_canopy_max(:),   & !< Capacity of the canopy skin reservoir [m water equivalent]
      & skinres_max(:),          & !< Total capacity of the skin reservoirs, i.e. soil and canopy
                                   !< [m water equivalent]
      & weq_pond_max(:),         & !< Maximum pond water storage [m water equivalent]
      & fract_snow(:),           & !< Snow cover fraction (not incl. canopy) []
      & fract_skin(:),           & !< Wet skin fraction (not incl. ponds) []
      & fract_pond(:),           & !< Actual pond fraction []
      & fract_pond_max(:),       & !< Maximum pond fraction []
      & evapotrans(:),           & !< Evapotranspiration (including sublimation) [kg m-2 s-1]
      & evapopot(:),             & !< Potential evaporation/sublimation (if there was enough water/ice)
                                   !< [kg m-2 s-1]
      & transpiration(:),        & !< Transpiration [kg m-2 s-1]
      & rain(:),                 & !< Liquid precipitation [kg m-2 s-1]
      & snow(:),                 & !< Solid precipitation [kg m-2 s-1]
      & wpi_rootzone(:),         & !< Water and ice in the root zone [m]
      & wpi_rootzone_max(:),     & !< Maximum amount of water or ice in the root zone [m]
      & wtr_soil(:),             & !< (Liquid) water content of the soil [m]
      & hyd_cond_sat_sl1(:),     & !< Saturated hydraulic conductivity of the uppermost soil layer [m/s]
      & ice_impedance(:)           !< Impedance of infiltration due to frozen soil moisture []
    REAL(wp), INTENT(inout) ::   &
      & wtr_skin(:),             & !< Water content of the skin reservoir (canopy and soil) [m]
      & weq_snow_soil(:),        & !< Amount of snow on the ground [m water equivalent]
      & weq_snow_can(:),         & !< Amount of snow on canopy [m water equivalent]
      & snow_soil_dens(:),       & !< Snow density (not incl. snow on canopy) [kg m-3]
      & wtr_pond(:),             & !< Water content of the pond reservoir [m]
      & ice_pond(:),             & !< Ice content of the pond reservoir [m water equivalent]
      & wtr_latflow_res_srf(:),  & !< Intermediary storage of surface runoff [m]
      & le_pc_remain(:)            !< Latent energy available for phase change [J m-2]
    REAL(wp), INTENT(out) ::     &
      & q_snocpymlt(:),          & !< Heating due to snow melt on canopy [W m-2]
      & snow_accum(:),           & !< Snow budget change within time step [m water equivalent]
      & snowmelt_soil(:),        & !< Snow/ice melt at land points (excluding canopy) [kg m-2 s-1]
      & pond_freeze(:),          & !< Amount of pond water freezing; on return [kg m-2 s-1]
      & pond_melt(:),            & !< Amount of pond ice melting; on return [kg m-2 s-1]
      & evapotrans_soil(:),      & !< Evapotranspiration from soil w/o snow, pond and skin reservoirs
                                   !< [kg m-2 s-1]
      & evapo_skin(:),           & !< Evaporation/sublimation from skin reservoir [kg m-2 s-1]
      & evapo_snow(:),           & !< Evaporation/sublimation from snow [kg m-2 s-1]
      & evapo_pond(:),           & !< Evaporation/sublimation from pond reservoir [kg m-2 s-1]
      & wtr_pond_net_flx(:),     & !< Net flux into pond reservoir; on return [kg m-2 s-1]
      & water_to_soil(:),        & !< Water available for infiltration into the soil [m /(time step)]
      & evapo_deficit(:),        & !< Evaporation from different storage than intended [m]
      & infilt(:),               & !< Infiltration into the soil [m /(time step)]
      & runoff(:),               & !< Surface runoff [m /(time step)]
      & runoff_horton(:),        & !< Hortonian surface runoff [m /(time step)]
      & wtr_latflow_srf(:)         !< Outflow from intermediary storage [m /(time step)]
    !
    !  local variables
    !
    REAL(wp) ::               &
      & rain_in_m,            & !< Rainfall within time step [m]
      & new_snow,             & !< Snowfall within time step [m water equivalent]
      & weq_snow_soil_old,    & !< Amount of snow on the ground before updating it [m water equivalent]
      & evapotrans_in_m,      & !< Evapotranspiration [m /(time step)]
      & evapotrans_soil_in_m, & !< Evapotranspiration from soil (w/o snow, skin and ponds) [m/(time step)]
      & evapo_soil_in_m,      & !< Soil evaporation within time step [m]
      & evapo_skin_in_m,      & !< Evaporation from skin reservoir within time step [m]
      & evapo_pond_in_m,      & !< Evaporation from pond reservoir within time step [m]
      & transpiration_in_m,   & !< Transpiration within time step [m]
      & evapo_snow_in_m,      & !< Sublimation of snow within time step [m]
      & snowmelt_can,         & !< Snow melt on canopy within time step [m water equivalent]
      & new_snow_can,         & !< Snowfall on canopy within time step [m]
      & new_snow_soil,        & !< Snowfall to soil within time step [m]
      & exp_t,                & !< Exponent for unloading of snow due to temperature
      & exp_w,                & !< Exponent for unloading of snow due to wind
      & snow_blown,           & !< Snow blown from canopy to the ground within time step [m water equivalent]
      & canopy_pre_snow,      & !< Snow depth on the canopy prior to its update [m water equivalent]
      & evapotrans_no_snow,   & !< Evapotranspiration without snow evaporation within the time step [m]
      & evapo_snow_pot_in_m,  & !< Potential sublimation of snow within time step - if there was enough snow
                                !< [m water equivalent]
      & evapo_skin_pot_in_m,  & !< Potential evaporation from skin reservoir within time step - if there
                                !< was enough water/ice in the skin reservoir [m water equivalent]
      & evapo_pond_pot_in_m,  & !< Potential evaporation/sublimation from pond reservoir within time step
                                !< - if there was enough water/ice in the pond reservoir [m water equivalent]
      & rain_to_skinres,      & !< Rainfall going into the skin reservoir within the time step [m]
      & wtr_skin_pre_rain,    & !< Amount of water in the skin reservoir prior to the update [m]
      & wtr_pond_inflow,      & !< Pond inflow within time step [m]
      & fract_pond_infil,     & !< Approximation of pond fraction for infiltration (with inflow included) [] ???
      & wtr_pond_infilt,      & !< Pond water infiltration into the soil within the time step [m]
      & wtr_pond_overflow,    & !< Water exceeding the maximum allowed pond water content (liquid + ice) [m]
                                !< within the time step [m water equivalent]
      ! TODO: remove pond_change_pot (never used)
      & pond_change_pot,      & !< Amount of water that could potentially be frozen or thawed [m]
      & t_snow_rel,           & !< Relative snow temperature
      & tau_snow_days,        & !< Relaxation/ageing constant
      & rho_snow_max            !< Temperature dependent maximum value of snow density

    INTEGER :: &
      & ic,                   & !< Looping index for grid cells
      & nc                      !< Number of grid cells
    !
    !  Parameters - compare JSBACH3 documentation, chapter 2, section 2.3.1.1 "Interception of snow by the canopy"
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp     !< Parameter for unloading of snow from canopy [K];
                                                   !< comp. JSBACH3 documentation
    REAL(wp), PARAMETER :: zc2 = 1.87E5_wp         !< Parameter for unloading of snow from canopy [Ks];
                                                   !< comp. JSBACH3 documentation
    REAL(wp), PARAMETER :: zc3 = 1.56E5_wp         !< Parameter for unloading of snow from canopy [m];
                                                   !< comp. JSBACH3 documentation
    REAL(wp), PARAMETER :: k_sat_clay = 8.50E-8_wp !< Saturated hydraulic conductivity for clay (Terra model) [m s-1]

    nc = SIZE(snowmelt_soil)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(rain_in_m, new_snow, weq_snow_soil_old) &
    !$ACC   PRIVATE(evapotrans_in_m, evapotrans_soil_in_m, evapo_soil_in_m, evapo_skin_in_m) &
    !$ACC   PRIVATE(transpiration_in_m, evapo_snow_in_m, snowmelt_can) &
    !$ACC   PRIVATE(new_snow_can, new_snow_soil, exp_t, exp_w, snow_blown, canopy_pre_snow) &
    !$ACC   PRIVATE(evapotrans_no_snow, evapo_snow_pot_in_m, evapo_skin_pot_in_m) &
    !$ACC   PRIVATE(evapo_pond_pot_in_m, rain_to_skinres, wtr_skin_pre_rain)

    DO ic = 1, nc
      snow_accum(ic)          = 0._wp
      weq_snow_soil_old       = weq_snow_soil(ic)
      snowmelt_soil(ic)       = 0._wp
      evapotrans_soil(ic)     = 0._wp
      evapo_skin(ic)          = 0._wp
      evapo_snow(ic)          = 0._wp
      pond_freeze(ic)         = 0._wp
      pond_melt(ic)           = 0._wp
      water_to_soil(ic)       = 0._wp
      evapo_deficit(ic)       = 0._wp
      wtr_pond_net_flx(ic)    = 0._wp
      wtr_latflow_srf(ic)     = 0._wp

      !-----------------------------------------------------------------
      ! Convert water fluxes to m water equivalent within the time step
      !-----------------------------------------------------------------
      rain_in_m           = rain(ic)          * dtime/rhoh2o
      new_snow            = snow(ic)          * dtime/rhoh2o   !in_m       ! snow_fall
      evapotrans_in_m     = evapotrans(ic)    * dtime/rhoh2o
      transpiration_in_m  = transpiration(ic) * dtime/rhoh2o
      evapo_snow_pot_in_m = fract_snow(ic)    * evapopot(ic) * dtime/rhoh2o
      evapotrans_no_snow  = evapotrans_in_m - evapo_snow_pot_in_m
      evapo_skin_pot_in_m = (1._wp - fract_snow(ic)) * fract_skin(ic) * evapopot(ic) * dtime/rhoh2o
      evapo_pond_pot_in_m = (1._wp - fract_snow(ic)) * fract_pond(ic) * evapopot(ic) * dtime/rhoh2o

      !---------------------------------------------
      !  Budgets of snow (canopy, ground)
      !---------------------------------------------
      !> We first treat the snow. Snowfall is intercepted by the canopy as long as the canopy skin
      !> reservoir is not completely filled. Excess snow is falling on the ground. Sublimation,
      !> melting and wind blow are reducing the snow on the canopy.

      ! Amount of snow intercepted by the canopy
      new_snow_can = MIN(new_snow * InterceptionEfficiency, MAX(skinres_canopy_max(ic) - weq_snow_can(ic), 0._wp))

      ! The remaining snow falls on the ground
      new_snow_soil = new_snow - new_snow_can

      ! Update snow on the canopy
      ! Note: evaporation happens from canopy, as long as there is enough snow
      canopy_pre_snow  = weq_snow_can(ic)
      weq_snow_can(ic) = MIN(MAX(0._wp, canopy_pre_snow + new_snow_can + evapo_snow_pot_in_m), skinres_canopy_max(ic))
      evapo_snow_in_m  = evapo_snow_pot_in_m - (weq_snow_can(ic) - canopy_pre_snow - new_snow_can)

      ! Unloading of snow from the canopy due to melting (compare JSBACH3 documentation, section 2.3.1.1)
      exp_t            = MAX(0._wp, t_air(ic) - zc1) / zc2 * dtime
      snowmelt_can     = weq_snow_can(ic) * (1._wp-EXP(-exp_t))
      weq_snow_can(ic) = weq_snow_can(ic) - snowmelt_can

      ! Unloading of snow from the canopy due to wind (compare JSBACH3 documentation, section 2.3.1.1)
      exp_w            = wind_10m(ic) / zc3 * dtime
      snow_blown       = weq_snow_can(ic) * (1._wp-EXP(-exp_w))
      weq_snow_can(ic) = weq_snow_can(ic)  - snow_blown
      new_snow_soil    = new_snow_soil + snow_blown

      ! Heating due to snow melt on canopy
      !-----------------------------------------------
      ! Note: This energy sink is not part of the surface energy balance and therefore not
      !     subtracted from the latent energy available for phase change. Instead, it is given
      !     directly to the atmosphere model.
      q_snocpymlt(ic) = snowmelt_can * rhoh2o * alf / dtime

      !> The amount of snow on the ground (soil) depends on the amount of new snow (including
      !> snow falling from the canopy), sublimation and melting.

      !  Snowfall and sublimation
      !-----------------------------------------------
      weq_snow_soil(ic) = weq_snow_soil(ic) + new_snow_soil + evapo_snow_in_m

      ! Correction if there was too much snow evaporation
      IF (weq_snow_soil(ic) < 0._wp) THEN
        evapotrans_no_snow = evapotrans_no_snow + weq_snow_soil(ic)
        evapo_deficit(ic)  = weq_snow_soil(ic)
        weq_snow_soil(ic)  = 0._wp
      ELSE
        evapo_deficit(ic)  = 0._wp
      END IF

      !  Snow melt
      !------------------------
      IF (.NOT. lstart) THEN      ! TODO: remove this condition
        IF (le_pc_remain(ic) > 0._wp .AND. weq_snow_soil(ic) > 0._wp) THEN
          snowmelt_soil(ic) = MIN(le_pc_remain(ic) / (rhoh2o * alf), weq_snow_soil(ic))
          le_pc_remain(ic)  = le_pc_remain(ic) - snowmelt_soil(ic) * rhoh2o * alf
          weq_snow_soil(ic) = weq_snow_soil(ic) - snowmelt_soil(ic) ! Reduce snow depth according to melting
        END IF
      END IF

      !  Snow budget and meltwater
      !-----------------------------------------------------

      !> Meltwater fills the soil skin reservoir. With [[t_hydro_config:l_ponds]]=true, excess water
      !> can fill surface water ponds. The remaining water will infiltrate the soil
      !> or go into the runoff.

      ! Add melt water from canopy to skin reservoir
      wtr_skin(ic)      = wtr_skin(ic) + snowmelt_can
      ! Excess water available for ponding, runoff, or infiltration
      water_to_soil(ic) = snowmelt_soil(ic) + MAX(0._wp, wtr_skin(ic) - skinres_max(ic))
      wtr_skin(ic)      = MIN(skinres_max(ic), wtr_skin(ic))

      ! Snow budget change of this time step
      !> Snowfall increases the amount of snow, while it is reduced by sublimation and melting.
      snow_accum(ic)    = new_snow + evapo_snow_pot_in_m - snowmelt_soil(ic) - snowmelt_can

      !> If [[t_hydro_config:snow_depth_max]] is set to a positive value, the snow amount is limited
      !> to this maximum value to avoid grid cells with infinitely growing snow depths in cooler climates.
      IF (snow_depth_max >= 0._wp) THEN      ! Namelist key: -1 for no limitation
        water_to_soil(ic) = water_to_soil(ic) + MAX(0._wp, weq_snow_soil(ic) - snow_depth_max)
        weq_snow_soil(ic) = MIN(snow_depth_max, weq_snow_soil(ic))
      END IF

      !  Freezing and melting of surface pond reservoir
      !------------------------
      ! 1. Melt pond ice using the available latent energy
      IF (le_pc_remain(ic) > 0._wp .AND. ice_pond(ic) > 0._wp) THEN
        pond_melt(ic)        = MIN(le_pc_remain(ic) / (rhoh2o * alf), ice_pond(ic))
        le_pc_remain(ic)     = le_pc_remain(ic) - pond_melt(ic) * rhoh2o * alf
        ice_pond(ic)         = ice_pond(ic) - pond_melt(ic)
        wtr_pond(ic)         = wtr_pond(ic) + pond_melt(ic)
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + pond_melt(ic)

      ! 2. Freeze pond water to compensate the latent energy deficit
      ELSE IF (le_pc_remain(ic) < 0._wp .AND. wtr_pond(ic) > 0._wp) THEN
        pond_freeze(ic)      = MIN(ABS(le_pc_remain(ic) / (rhoh2o * alf)), wtr_pond(ic))
        le_pc_remain(ic)     = le_pc_remain(ic) + pond_freeze(ic) * rhoh2o * alf
        wtr_pond(ic)         = wtr_pond(ic) - pond_freeze(ic)
        ice_pond(ic)         = ice_pond(ic) + pond_freeze(ic)
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) - pond_freeze(ic)
      END IF

      !-----------------------------------------------------------
      !   Budget of water in skin and pond reservoirs
      !-----------------------------------------------------------
      !> Similar to the snow, also a certain fraction of rainwater is intercepted by the canopy,
      !> while the rest reaches the ground. If the pond scheme is active ([[t_hydro_config:l_ponds]]=true)
      !> water exceeding the skin reservoir capacities of canopy and soil will flow into ponds,
      !> until the maximum storage capacity is reached. Evaporation and sublimation are reducing
      !> pond storages.

      ! Interception of rainwater in the canopy
      rain_to_skinres = MIN(rain_in_m * InterceptionEfficiency, MAX(skinres_max(ic) - wtr_skin(ic), 0._wp))
      ! Remaining rainwater reaches the soil
      ! (Note: the MAX in following line accounts for possible precision error (small negative values))
      water_to_soil(ic) = water_to_soil(ic) + MAX(rain_in_m - rain_to_skinres, 0._wp)
      wtr_skin_pre_rain = wtr_skin(ic)
      wtr_skin(ic)      = MIN(skinres_max(ic), MAX(0._wp, wtr_skin_pre_rain + rain_to_skinres + evapo_skin_pot_in_m))

      ! Note: At this point we are ignoring the amount of dew that exceeds the skin reservoir
      !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
      evapo_skin_in_m = wtr_skin(ic) - (wtr_skin_pre_rain + rain_to_skinres)

      ! Limit pond evaporation to the available pond water and ice
      evapo_pond_in_m = MAX(-(wtr_pond(ic) + ice_pond(ic)), evapo_pond_pot_in_m)
      ! Only evaporate from ponds if demand cannot be fulfilled from skin reservoir and transpiration
      evapo_pond_in_m = MAX(evapo_pond_in_m, &
        & MIN(0._wp, evapotrans_no_snow - evapo_skin_in_m - transpiration_in_m))

      ! Update evaporation required from the soil
      evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m - evapo_pond_in_m
      IF (wtr_soil(ic) + evapotrans_soil_in_m < 0._wp) THEN
        ! If the soil does not contain enough moisture, evaporate from ponds instead if possible
        evapo_pond_in_m  = MAX(-(wtr_pond(ic) + ice_pond(ic)), &
          &                    evapo_pond_in_m + (wtr_soil(ic) + evapotrans_soil_in_m))
        evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m - evapo_pond_in_m
      END IF
      evapo_soil_in_m = evapotrans_soil_in_m - transpiration_in_m

      ! Add dew to liquid pond storage, sublimate pond ice or evaporate pond water
      !   Note that pond ice sublimation has priority over pond water evaporation because
      !   pond ice is expected to form as a top layer and would inhibit evaporation of
      !   the water below.
      IF (evapo_pond_in_m > 0._wp) THEN
        ! Adding dew to pond water storage
        wtr_pond(ic) = wtr_pond(ic) + evapo_pond_in_m
      ELSE IF (-evapo_pond_in_m <= ice_pond(ic)) THEN
        ! Sublimating pond ice (if existing)
        ice_pond(ic) = ice_pond(ic) + evapo_pond_in_m
      ELSE IF (-evapo_pond_in_m < wtr_pond(ic) + ice_pond(ic)) THEN
        ! Sublimating pond ice (if existing) and/or evaporating pond water
        wtr_pond(ic) = wtr_pond(ic) + (evapo_pond_in_m + ice_pond(ic))
        ice_pond(ic) = 0._wp
      ELSE
        ! Evapotranspiration cannot be satisfied with pond content, remainder goes into deficit
        evapo_deficit(ic) = evapo_deficit(ic) + (evapo_pond_in_m + wtr_pond(ic) + ice_pond(ic))
        wtr_pond(ic) = 0._wp
        ice_pond(ic) = 0._wp
      END IF
      wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + evapo_pond_in_m

      IF (evapo_skin_pot_in_m < 0._wp) THEN
        evapo_deficit(ic) = evapo_deficit(ic) + evapo_skin_pot_in_m - evapo_skin_in_m
      END IF

      ! Positive values of evaporation and transpiration (dew) are added to water_to_soil.
      ! Negative fluxes change soil moisture later in calc_soil_hydrology.
      IF (evapo_soil_in_m > 0._wp)    water_to_soil(ic) = water_to_soil(ic) + evapo_soil_in_m
      IF (transpiration_in_m > 0._wp) water_to_soil(ic) = water_to_soil(ic) + transpiration_in_m

      !-----------------------------------------------------------
      !   Compute snow density
      !-----------------------------------------------------------

      IF (l_dynsnow) THEN
        !> Dynamic snow density calculations follow eq. (3) by E. Heise et al. (2006). The density
        !> depends on a medium snow temperature and snow age. This parametrization is also implemented
        !> in the multi layer snow model by DWD for operational usage.
        IF (weq_snow_soil_old * rhoh2o / snow_soil_dens(ic) > snow_depth_min) THEN
          t_snow_rel         = (MIN(tmelt, t_snow_mean(ic)) - csnow_tmin) / (tmelt - csnow_tmin)
          tau_snow_days      = MIN(MAX(0.05_wp, crhosmint + (crhosmaxt - crhosmint) * t_snow_rel), crhosmaxt)
          rho_snow_max       = crhosmax_tmin  + MAX(-0.25_wp, t_snow_rel) * (dens_snow_max - crhosmax_tmin)
          snow_soil_dens(ic) = MAX(snow_soil_dens(ic), rho_snow_max + &
            &                  (snow_soil_dens(ic) - rho_snow_max) * EXP(-tau_snow_days * dtime/seconds_per_day))
        ELSE
          snow_soil_dens(ic) = dens_snow_min  ! fresh snow
        END IF
        ! Calculate weighted mean density from old and fresh snow if any compaction happened already
        IF (weq_snow_soil(ic) > weq_snow_soil_old .AND. snow_soil_dens(ic) > dens_snow_min) THEN
          snow_soil_dens(ic) = (snow_soil_dens(ic) * weq_snow_soil_old                    &
            &                  + dens_snow_min * (weq_snow_soil(ic) - weq_snow_soil_old)) &
            &                 / weq_snow_soil(ic)
        END IF
      ELSE
        !> If snow density dynamics are switched off ([[t_sse_config:l_dynsnow]]=false), snow density is fixed
        !> to a constant value.
        snow_soil_dens(ic) = dens_snow
      END IF

      !---------------------
      !  Pond water, infiltration and surface runoff
      !---------------------
      !> So far, the snow amount and skin reservoir surface water storages have been updated. We now need to
      !> handle the excess water. Depending on surface conditions and soil properties it fills surface
      !> depressions forming ponds, it infiltrates into the soil or feeds into rivers via surface runoff.

      ! Note, that infiltration, runoff, and pond storage are further modified in the soil hydrology routine.

      IF (ltpe_closed) THEN
        ! Terra planet setup without runoff and ponds
        runoff(ic) = 0._wp
        infilt(ic) = water_to_soil(ic)
        runoff_horton(ic) = 0._wp
      ELSE
        !> The flow scheme for ponds is based on the [WEED scheme](https://doi.org/10.5194/tc-15-1097-2021)
        !> and is modified to work with the semi_distributed (ARNO) scale and the uniform (e.g. point scale or
        !> very high resolution) scale assumptions.

        ! Calculate pond water inflow and infiltration from ponds
        wtr_pond_inflow = 0._wp
        IF (hydro_scale == Semi_Distributed_) THEN
          ! The lateral inflow into ponds originates from the whole grid cell and flows towards the
          ! depressions. Thus, most of the available water is expected to end up in ponds even if not
          ! the whole cell is flooded.
          wtr_pond_inflow = water_to_soil(ic) * fract_pond(ic)**(1.0_wp/3.0_wp)

          ! As ponds form in depressions, no lateral outflow is assumed but pond water rather infiltrates
          ! into the soil.
          ! We use the assumptions
          ! - Ponds usually correspond to clay rich soils, thus saturated hydraulic conductivity for
          !   clay is used to limit infiltration.
          ! - Infiltration is only allowed for ponds without ice (used as proxy for frozen ground).
          ! - Inflow might not be distributed equally over the whole time step, thus we assume
          !   only half of it is available for infiltration all the time.
          IF (ice_pond(ic) > 0._wp) THEN
            wtr_pond_infilt = 0._wp
          ELSE
            wtr_pond_infilt = MIN(wtr_pond(ic) + 0.5_wp * wtr_pond_inflow, &
              &                   k_sat_clay * dtime * fract_pond(ic))
          END IF

        ELSE IF (hydro_scale == Uniform_) THEN
          IF (weq_pond_max(ic) > 1.0e-10_wp) THEN
            ! If any potential pond fraction exists for the actual grid cell (or tile) the available water
            ! first goes into the pond storage.
            wtr_pond_inflow  = water_to_soil(ic) * fract_pond_max(ic)**(1.0_wp/3.0_wp)

            ! Approximation of the actual pond fraction
            fract_pond_infil = MIN(1._wp, ((wtr_pond(ic) + wtr_pond_inflow) / weq_pond_max(ic))**0.5_wp) &
              &              * fract_pond_max(ic)
            ! For the uniform scale no assumption is made on the subgrid scale distribution of soil textures.
            ! Infiltration is only limited by ice_impedance, and not cut in case of pond ice.
            wtr_pond_infilt = MIN(wtr_pond(ic) + wtr_pond_inflow, &
              &                   ice_impedance(ic) * hyd_cond_sat_sl1(ic) * dtime * fract_pond_infil)
          ELSE
            wtr_pond_inflow  = 0._wp
            wtr_pond_infilt  = 0._wp
          END IF
        END IF

        ! Update pond storages and handle overflow
        wtr_pond(ic) = wtr_pond(ic) + wtr_pond_inflow - wtr_pond_infilt
        IF (wtr_pond(ic) + ice_pond(ic) > weq_pond_max(ic)) THEN
          ! If the water/ice exceeds the pond storage capacity, first liquid water goes to the overflow.
          wtr_pond_overflow = MAX(0._wp, wtr_pond(ic) + ice_pond(ic) - weq_pond_max(ic))
          wtr_pond(ic)      = weq_pond_max(ic) - ice_pond(ic)
          IF (wtr_pond(ic) < 0.0_wp) THEN
            ice_pond(ic) = ice_pond(ic) + wtr_pond(ic)
            wtr_pond(ic) = 0._wp
          END IF
        ELSE
          ! No overflow as pond storage capacity is sufficient.
          wtr_pond_overflow = 0._wp
        END IF
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + wtr_pond_inflow  &
          &                  - wtr_pond_infilt - wtr_pond_overflow

        ! Update amount of water available for runoff and infiltration
        IF (hydro_scale == Semi_Distributed_) THEN
          ! Water infiltrating from ponds and pond overflow is added to water_to_soil and is available
          ! for runoff and infiltration.
          water_to_soil(ic) = water_to_soil(ic) - wtr_pond_inflow + wtr_pond_infilt &
            &               + wtr_pond_overflow
        ELSE IF (hydro_scale == Uniform_) THEN
          ! Here we assume that the outflow from ponds does not uniformly inundate the surrounding surfaces
          ! and can infiltrate into the soil, but rather that it runs off as surface runoff into the stream
          ! or the connected downstream tile.
          water_to_soil(ic) = water_to_soil(ic) - wtr_pond_inflow + wtr_pond_infilt
        END IF

        !---------------------
        ! Infiltration and surface runoff
        !---------------------
        !> When considering larger scales ([[t_hydro_config:hydro_scale]]=semi_distributed) runoff is
        !> computed via the ARNO Scheme, which implicitly considers subgrid scale orographic heterogeneity.
        !> For very high resolution, site level, or HydroTile simulations, surface conditions are assumed
        !> to be uniform ([[t_hydro_config:hydro_scale]]=uniform). With this approach two kinds of runoff
        !> are calculated: Hortonien runoff depending on terrain steepness and the runoff depending on
        !> soil moisture.

        IF (hydro_scale == Semi_Distributed_) THEN
          !---------------------------------------
          ! Compute fluxes based on ARNO scheme (semi-distributed scale approach)
          !---------------------------------------
          ! The Arno-scheme is designed for large scales and includes assumptions for subgrid scale
          ! variablility. This may need reconsideration with very high resolution simulations.

          ! Infiltration is based on the root zone, exceeding water goes to the runoff. With the ARNO
          ! scheme, individual runoff components cannot be distinguished and the horton runoff
          ! diagnostic flux is set to zero.
          CALL arno_scheme(l_infil_subzero, t_soil_sl1(ic),        &
            &              water_to_soil(ic), wpi_rootzone(ic),    &
            &              steepness(ic), wpi_rootzone_max(ic),    &
            &              runoff(ic), infilt(ic))
          runoff_horton(ic) = 0._wp

        ELSE IF (hydro_scale == Uniform_) THEN
          !---------------------------------------
          ! Compute fluxes based on uniform scale approach, e.g. for HydroTiles or site level setup
          !---------------------------------------
          ! The runoff has two components:
          !   1) Orographic runoff which depends only on steepness and hydraulic conductivity
          !      (Horton runoff)
          !   2) Excess water which exceeds the water holding capacity of the root zone

          ! This approach uses the saturated hydraulic conductivity as upper limit for infiltration
          ! assuming any precipitation event quickly fills up the thin top layer soil within our
          ! standard model time step length.
          infilt(ic) = MIN(water_to_soil(ic), &
            & ice_impedance(ic) * hyd_cond_sat_sl1(ic) * (1.0_wp - SIN(steepness(ic)*pi/2._wp)) * dtime)
          runoff(ic) = MAX(water_to_soil(ic) - infilt(ic), 0._wp) + wtr_pond_overflow
          runoff_horton(ic) = runoff(ic)

          ! Handling of quasi lateral fluxes using a simple reservoir to account for the lag between
          ! runoff generation and water reaching the river / downstream tile.
          CALL surfhyd_lat(flowlag(ic), runoff(ic), &
            &              wtr_latflow_res_srf(ic), wtr_latflow_srf(ic))

          IF (l_latflow_to_streamflow) THEN
            runoff(ic)          = wtr_latflow_srf(ic)
            wtr_latflow_srf(ic) = 0._wp
          END IF
        END IF

      END IF   ! not in terraplanet setup

      ! Transform fluxes from [m (time step)] to [kg m-2 s-1]
      snowmelt_soil(ic)    = snowmelt_soil(ic)    * rhoh2o / dtime
      pond_freeze(ic)      = pond_freeze(ic)      * rhoh2o / dtime
      pond_melt(ic)        = pond_melt(ic)        * rhoh2o / dtime
      evapotrans_soil(ic)  = evapotrans_soil_in_m * rhoh2o / dtime
      evapo_skin(ic)       = evapo_skin_in_m      * rhoh2o / dtime
      evapo_snow(ic)       = evapo_snow_in_m      * rhoh2o / dtime
      evapo_pond(ic)       = evapo_pond_in_m      * rhoh2o / dtime
      wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) * rhoh2o / dtime
      wtr_latflow_srf(ic)  = wtr_latflow_srf(ic)  * rhoh2o / dtime

    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

  END SUBROUTINE calc_surface_hydrology_land

  ! ===============================================================================================================================
  !>
  !> #### Compute surface hydrology on glaciers
  !>
  !> On glaciers, we assume an unlimited amount of ice/snow. This amount increases with snow fall
  !> and decreases by sublimation and snow melt. Rain fall and snow melt directly go into runoff,
  !> there is no glacier infiltration.

#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_surface_hydrology_glacier ( &
    & lstart, dtime,                          &
    & fract_snow,                             &
    & evapotrans, evapopot, rain, snow,       &
    & weq_glac, le_pc_remain, q_snocpymlt,    &
    & snowmelt, runoff_glac,                  &
    & pme_glacier)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: rhoh2o, alf

    LOGICAL,  INTENT(in) :: &
      & lstart                 !< True: beginning of the experiment
    REAL(wp), INTENT(in) :: &
      & dtime, &               !< Time stp length [s]
      & fract_snow, &          !< Snow cover fraction []
      & evapotrans, &          !< Evapotranspiration incl. sublimation [kg m-2 s-1]
      & evapopot, &            !< Potential evaporation/sublimation [kg m-2 s-1]
      & rain, &                !< Liquid precipitation [kg m-2 s-1]
      & snow                   !< Solid precipitation [kg m-2 s-1]
    REAL(wp), INTENT(inout) :: &
      & weq_glac, &            !< Glacier depth (snow and ice) [m water equivalent]
      & le_pc_remain           !< Latent energy available for snow melt [J m-2]
    REAL(wp), INTENT(out) :: &
      & q_snocpymlt, &         !< Heating due to snow melt on canopy [W m-2]
      & snowmelt, &            !< Snow/ice melt at glacier points [kg m-2 s-1]
      & runoff_glac, &         !< Glacier runoff (rain+snow/ice melt, no calving) [kg m-2 s-1]
      & pme_glacier            !< Precipitation minus sublimation on glacier [kg m-2 s-1]

    !
    !  local variables
    !
    REAL(wp) ::              &
      & rain_in_m,           & !< Amount of rainfall within time step [m]
      & new_snow,            & !< Amount of snowfall within time step [m water equivalent]
      & evapotrans_in_m,     & !< Amount of sublimation within time step [m water equivalent]
      & evapo_snow_pot_in_m, & !< Potential snow sublimation [m water equivalent /(time step)]
      & pme_glacier_in_m,    & !< P-E on glacier [m /(time step)]
      & runoff_glac_in_m,    & !< Glacier runoff [m /(time step)]
      & snowmelt_in_m          !< Snow/ice melt [m /(time step)]

    !----------------------------------------------------------------------------------------------

    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m          = rain       * dtime/rhoh2o
    new_snow           = snow       * dtime/rhoh2o
    evapotrans_in_m    = evapotrans * dtime/rhoh2o
    evapo_snow_pot_in_m = fract_snow * evapopot * dtime/rhoh2o

    !  Snowfall and sublimation on glaciers
    !---------------------------------------

    weq_glac         = weq_glac  + new_snow + evapo_snow_pot_in_m   ! Glacier depth [m water equivalent]
    pme_glacier_in_m = rain_in_m + new_snow + evapotrans_in_m       ! P-E on glaciers [m water equivalent]
    runoff_glac_in_m = rain_in_m                                    ! No infiltration on glaciers


    !  Snow and glacier melt
    !------------------------

    snowmelt_in_m = 0._wp
    IF (.NOT. lstart) THEN      ! TODO: condition not necessary anymore, remove
      IF (le_pc_remain > 0._wp) THEN
        snowmelt_in_m     = le_pc_remain / (rhoh2o * alf)     ! There is an unlimited amount of snow
        weq_glac          = weq_glac - snowmelt_in_m          ! Reduce glacier depth according to melting
        runoff_glac_in_m  = runoff_glac_in_m + snowmelt_in_m  ! Add melt water to the runoff
      END IF
    END IF
    le_pc_remain = 0._wp   ! All latent energy used for snow melt, no remaining residual.

    q_snocpymlt = 0._wp    ! No canopy and thus now snow melt on canopy.


    ! Unit conversion [(m water equivalent)/(time step)] -> [kg m-2 s-1]
    !--------------------------------------------------------------------

    snowmelt    = snowmelt_in_m    * rhoh2o / dtime
    runoff_glac = runoff_glac_in_m * rhoh2o / dtime
    pme_glacier = pme_glacier_in_m * rhoh2o / dtime

  END SUBROUTINE calc_surface_hydrology_glacier

  ! ===============================================================================================================================
  !>
  !> #### Compute soil hydrology on non-glacier land
  !>
  !> In this subroutine we calculate the vertical movement of water in the soil column. Main calculations
  !> actually happens in subroutine [[soilhyd]]. We here prepare the subroutine call and perform
  !> some processing afterwards.
  !>
  SUBROUTINE calc_soil_hydrology(                                                        &
    ! in
    & nc, l_fract, l_pf_soil, lat, lon, nsoil, dtime,                                    &
    & ltpe_closed, ltpe_open, enforce_water_budget, l_latflow_to_streamflow,             &
    & soilhydmodel, interpol_mean, hydro_scale,                                          &
    & model_scheme,                                                                      &
    & w_soil_wilt_fract,                                                                 &
    & soil_depth_sl, root_depth_sl,                                                      &
    & hyd_cond_sat_sl, matric_pot_sl, bclapp_sl, pore_size_index_sl, vol_porosity_sl,    &
    & vol_field_cap_sl, vol_p_wilt_sl, vol_wres_sl, wtr_soil_pot_scool_sl,  slope,       &
    & flowlag, fract_pond_max, evapotrans_soil, transpiration, ice_pond, weq_pond_max,   &
    ! inout
    & infilt, runoff, wtr_soil_sl, ice_soil_sl, wtr_pond, wtr_pond_net_flx,              &
    & tpe_overflow, evapo_deficit, wtr_latflow_res_srf, wtr_latflow_res_sl,              &
    ! out
    & wtr_wsat_sl, wtr_field_cap_sl, wtr_p_wilt_sl, wtr_wres_sl, runoff_dunne, drainage, &
    & drainage_sl, wtr_transp_down_sl, wtr_soilhyd_res, drainage_lowest_soil_layer,      &
    & wtr_latflow_sl,                                                                    &
    ! optional: only used with QUINCY
    & ftranspiration_sl                                                                  &
    & )

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_constants,        ONLY: Uniform_

    !TODO: sequence of arguments is NOT consistent with the list of declaration
    INTEGER, INTENT(in) ::       &
      & nc,                      & !< Vector length
      & nsoil,                   & !< Number of soil layers (vertical grid dimension)
      & enforce_water_budget,    & !< Handling of water balance errors
      & soilhydmodel,            & !< Model scheme for soil hydraulic properties
      & interpol_mean,           & !< Interpolation scheme for hydraulic properties
      & hydro_scale,             & !< Effective area represented in calculations:
                                   !< Semi_Distributed: ARNO scheme; Uniform: e.g. site level
      & model_scheme               !< Land surface model
    REAL(wp), INTENT(in) ::      &
      & dtime,                   & !< Time step length [s]
      & w_soil_wilt_fract          !< Soil moisture at wilting point [m]
    LOGICAL, INTENT(in) ::       &
      & ltpe_closed,             & !< Terraplanet setup: with closed water balance
      & ltpe_open,               & !< Terraplanet setup: lower soil layer kept wet
      & l_latflow_to_streamflow    !< T: drainage goes directly into streamflow
    LOGICAL, INTENT(in), DIMENSION(:) :: &
      & l_fract,                 & !< Tile has a non-zero grid cell fraction
      & l_pf_soil                  !< Permafrost within the soil above bedrock
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & soil_depth_sl,           & !< Soil depth until bedrock within each layer [m]
      & root_depth_sl,           & !< Root depth within each soil layer [m]
      & hyd_cond_sat_sl,         & !< Saturated hydraulic conductivity on soil layers [m s-1]
      & matric_pot_sl,           & !< Saturated matric potential on soil layers [m]
      & bclapp_sl,               & !< Clapp & Hornberger exponent b for soil layer []
      & pore_size_index_sl,      & !< Pore size index of soil layer []
      & vol_porosity_sl,         & !< Volumetric porosity of soil layer []
      & vol_field_cap_sl,        & !< Volumetric field capacity of soil layer []
      & vol_p_wilt_sl,           & !< Volumetric wilting point of soil layers []
      & vol_wres_sl,             & !< Volumetric residual soil water content on soil layers []
      & wtr_soil_pot_scool_sl      !< Potential amount of supercooled water [m]
    REAL(wp), INTENT(in), DIMENSION(:) :: &
      & lat,                     & !< Latitude [degree]
      & lon,                     & !< Longitude [degree]
      & slope,                   & !< Slope of grid cell or tile []
      & flowlag,                 & !< Lag factor accounting for retention in surface runoff []
      & fract_pond_max,          & !< Maximum tile fraction with ponds []
      & evapotrans_soil,         & !< Evapotranspiration from soil [kg m-2 s-1]
      & transpiration,           & !< Transpiration [kg m-2 s-1]
      & ice_pond,                & !< Amount of pond ice [m water equivalent]
      & weq_pond_max               !< Maximum amount of pond water or ice [m water equivalent]
    REAL(wp), INTENT(inout), DIMENSION(:) :: &
      & tpe_overflow,            & !< Terraplanet: Overflow reservoir for soil water [m]
      & infilt,                  & !< Infiltration [kg m-2 s-1]
      & runoff,                  & !< Surface runoff [kg m-2 s-1]
      & evapo_deficit,           & !< Evaporation deficit due to inconsistencies [m]
      & wtr_pond,                & !< Amount of pond water [m]
      & wtr_pond_net_flx,        & !< Net inflow into surface water ponds [kg m-2 s-1]
      & wtr_latflow_res_srf        !< Water in intermediate storage for surface runoff [m]
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & wtr_soil_sl,             & !< Amount of water in soil layer [m]
      & ice_soil_sl,             & !< Amount of ice in soil layer [m]
      & wtr_latflow_res_sl         !< Water in "lateral" subsurface drainage reservoirs [m]
    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & runoff_dunne,            & !< Dunne component of surface runoff (saturation excess) [kg m-2 s-1]
      & drainage,                & !< Drainage [kg m-2 s-1]
      & wtr_soilhyd_res,         & !< Residual of vertical soil water transport scheme [m]
      & drainage_lowest_soil_layer !< Bottom layer drainage [kg m-2 s-1]
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & wtr_wsat_sl,             & !< Saturation capacity of the soil layer (reduced by ice) [m]
      & wtr_field_cap_sl,        & !< Field capacity of the soil layer (reduced by ice) [m]
      & wtr_p_wilt_sl,           & !< Water content at wilting point in soil layer (reduced by ice) [m]
      & wtr_wres_sl,             & !< Residual water content in soil layer (reduced by ice) [m]
      & wtr_transp_down_sl,      & !< Vertical water transport to the below soil layer [kg m-2 s-1]
      & drainage_sl,             & !< Drainage [kg m-2 s-1]
      & wtr_latflow_sl             !< Outflow from the intermediary reservoir representing
                                   !< lateral flow of subsurface drainage [kg m-2 s-1]
    REAL(wp), OPTIONAL, INTENT(in), DIMENSION(:,:) :: &
      & ftranspiration_sl          !< Transpiration fraction []

    ! Local variables
    REAL(wp) ::                  &
      & evapo_soil_in_m(nc),     & !< Evaporation from the soil (without transpiration) [m /(time step)]
      & transpiration_in_m(nc),  & !< Transpiration within time step [m /(time step)]
      & drain_bot(nc),           & !< Drainage towards bedrock [m /(time step)]
      & wtr_pond_corr(nc),       & !< Correction term for pond fluxes with saturated soils [m /(time step)]
      & hyd_cond_bot(nc)           !< Hydraulic conductivity of the bottom soil layer [m/s]
    INTEGER :: ic, is

    INTEGER, PARAMETER :: ilog = 0 !< Switch for debugging output

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(evapo_soil_in_m, transpiration_in_m, wtr_pond_corr, drain_bot, hyd_cond_bot)


    ! Preparation: Calculation of soil saturation, field capacity and wilting point
    !--------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        wtr_wsat_sl(ic,is) = MAX(vol_porosity_sl (ic,is) * soil_depth_sl(ic,is) - ice_soil_sl(ic,is), 0._wp)
        IF (vol_porosity_sl(ic,is) > 0._wp) THEN
          wtr_field_cap_sl(ic,is) = wtr_wsat_sl(ic,is) * (vol_field_cap_sl(ic,is) / vol_porosity_sl(ic,is))
          wtr_p_wilt_sl(ic,is)    = wtr_wsat_sl(ic,is) * (vol_p_wilt_sl(ic,is)    / vol_porosity_sl(ic,is))
          wtr_wres_sl(ic,is)      = wtr_wsat_sl(ic,is) * (vol_wres_sl(ic,is)      / vol_porosity_sl(ic,is))
        ELSE
          wtr_field_cap_sl(ic,is) = 0._wp
          wtr_p_wilt_sl(ic,is)    = 0._wp
          wtr_wres_sl(ic,is)      = 0._wp
        END IF
        wtr_latflow_sl(ic,is)     = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL


    ! Preparation: Initialize diagnostic runoff and drainage variables
    !--------------

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      runoff_dunne(ic) = 0._wp
      drain_bot(ic)    = 0._wp
      drainage(ic)     = 0._wp
      drainage_lowest_soil_layer(ic) = 0._wp      ! output variable, equals drain_bot
      ! Note: any dew is removed from these fluxes, as dew is already applied
      !       during the surface hydrology computation
      evapo_soil_in_m(ic)    = MIN(0._wp, evapotrans_soil(ic) - transpiration(ic)) * dtime / rhoh2o
      transpiration_in_m(ic) = MIN(0._wp, transpiration(ic))  * dtime / rhoh2o
    END DO
    !$ACC END PARALLEL LOOP


    ! Calculation of vertical soil water movement for multiple soil layers
    !---------------------------------------------
    CALL soilhyd( &
      &    nc, l_fract, l_pf_soil, lat, lon, nsoil, dtime, enforce_water_budget,       &
      &    slope, soilhydmodel, interpol_mean, hydro_scale,                            &
      &    model_scheme,                                                               &
      &    w_soil_wilt_fract,                                                          &
      &    soil_depth_sl, wtr_wres_sl, wtr_p_wilt_sl, wtr_field_cap_sl, wtr_wsat_sl,   &
      &    hyd_cond_sat_sl, vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl,          &
      &    vol_wres_sl, bclapp_sl, matric_pot_sl, pore_size_index_sl,                  &
      &    wtr_soil_pot_scool_sl,                                                      &  ! in
      &    ice_soil_sl, wtr_soil_sl,                                                   &  ! inout
      &    transpiration_in_m, evapo_soil_in_m, root_depth_sl,                         &  ! in
      &    infilt, runoff_dunne,                                                       &  ! inout
      &    drain_bot, drainage_sl, wtr_transp_down_sl,                                 &  ! out
      &    tpe_overflow, evapo_deficit,                                                &  ! inout
      &    wtr_soilhyd_res, hyd_cond_bot,                                              &  ! out
      &    ltpe_closed, ltpe_open, ftranspiration_sl)                                     ! in

    IF (hydro_scale == Uniform_) THEN
      ! Handling of quasi lateral fluxes (within the grid cell) and respective reservoirs in case
      ! of uniform scale calculations (e.g. with HydroTiles).
      CALL calc_soilhyd_lateral( &
        &    l_pf_soil, nc, nsoil, dtime,                    &
        &    flowlag, hyd_cond_bot,                          &
        &    soil_depth_sl, vol_field_cap_sl,                &
        &    ice_soil_sl, wtr_soil_sl,                       &
        &    wtr_latflow_res_sl, drainage_sl, drain_bot,     &
        &    wtr_latflow_sl)

      IF (l_latflow_to_streamflow) THEN
        ! Subsurface lateral flow goes to the streamflow, i.e. it is handled as drainage.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is = 1, nsoil
          DO ic = 1, nc
            drainage_sl(ic,is)    = wtr_latflow_sl(ic,is)
            wtr_latflow_sl(ic,is) = 0._wp
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF


    ! Postprocessing: Aggregate or distribute the different runoff flux components
    !-----------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        drainage(ic) = drainage(ic) + drainage_sl(ic, is)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      drainage(ic)  = drainage(ic) + drain_bot(ic)

      ! Handle infiltration overflow (i.e. Dunne runoff)
      IF (fract_pond_max(ic) > EPSILON(1._wp) .AND. runoff_dunne(ic) > 0._wp) THEN
        ! If the pond scheme is active, Dunne runoff goes into the pond storages.

        ! Add Dunne type runoff to pond storage if capacity is available
        wtr_pond_corr(ic)    = MIN(MAX(0._wp, weq_pond_max(ic) - wtr_pond(ic) - ice_pond(ic)), runoff_dunne(ic))
        wtr_pond(ic)         = wtr_pond(ic) + wtr_pond_corr(ic)
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + wtr_pond_corr(ic) * rhoh2o / dtime
        ! Check for remaining Dunne type runoff and add to surface runoff
        runoff_dunne(ic)     = MAX(0._wp, runoff_dunne(ic) - wtr_pond_corr(ic))
        runoff(ic)           = runoff(ic) + runoff_dunne(ic)

      ELSE
        ! Without ponds, all Dunne runoff goes directly into surface runoff
        runoff(ic)    = runoff(ic) + runoff_dunne(ic)
      END IF

      IF (hydro_scale == Uniform_) THEN
        ! In case of site level or HydroTile setups we do not add Dunne type runoff to runoff
        ! but to the lateral surface flow reservoir.
        ! (The outflow from this reservoir is calculated later in update_surface_hydrology.)
        wtr_latflow_res_srf(ic) = wtr_latflow_res_srf(ic) + runoff_dunne(ic)
        runoff(ic)              = runoff(ic)              - runoff_dunne(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP


    ! Unit conversion: [m /(time step)] -> [kg m-2 s-1]
    !------------------

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      infilt(ic)        = infilt       (ic) * rhoh2o / dtime
      runoff(ic)        = runoff       (ic) * rhoh2o / dtime
      runoff_dunne(ic)  = runoff_dunne (ic) * rhoh2o / dtime
      drainage(ic)      = drainage     (ic) * rhoh2o / dtime
      evapo_deficit(ic) = evapo_deficit(ic) * rhoh2o / dtime
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil
      DO ic = 1, nc
        wtr_transp_down_sl(ic,is) = wtr_transp_down_sl(ic,is) * rhoh2o / dtime
        drainage_sl(ic,is)        = drainage_sl(ic,is)     * rhoh2o / dtime
        wtr_latflow_sl(ic,is)     = wtr_latflow_sl(ic,is)  * rhoh2o / dtime
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      drainage_lowest_soil_layer(ic) = drain_bot(ic) * rhoh2o / dtime
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE calc_soil_hydrology

  ! ===============================================================================================================================
  !>
  !> #### Calculation of subgrid scale orographic features
  !>
  !> In this routine we derive a steepness factor defining the subgrid scale steepness of the orography
  !> from the subgrid scale orographic standard deviation following the ARNO Scheme (compare E. Todini,
  !> The ARNO rainfall-runoff model 1996 [DOI](https://doi.org/10.1016/S0022-1694(96)80016-3)).
  !>
  SUBROUTINE calc_orographic_features(nc, nlat, oro_stddev, steepness)

    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max

    INTEGER,  INTENT(in) :: &
      & nc, &            !< Vector length
      & nlat             !< Effective number of latitudes (for ICON grid)

    REAL(wp), INTENT(in) :: &
      & oro_stddev(:)    !< Standard deviation of orography [m]
    REAL(wp), INTENT(out) :: &
      & steepness(:)     !< Parameter defining the subgrid slope distribution []

    INTEGER :: ic        !< Grid cell index

    REAL(wp) :: &
      & sigma_0, &       !< Minimum value (100 m); below b = 0.01
      & sigma_max        !< Resolution-dependent maximum

    sigma_0   = oro_var_min
    sigma_max = oro_var_max * 64._wp / REAL(nlat, wp) !Todo: site level

    !----------------------------------------------------------
    ! The steepness parameter (b) is defined as
    !    b = (oro_stddev - sigma_0) / (oro_stddev + sigma_max)
    !----------------------------------------------------------

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      steepness(ic) = MAX(0._wp, oro_stddev(ic) - sigma_0) / (oro_stddev(ic) + sigma_max)
      ! Limit parameter to realistic bounds
      steepness(ic) = MAX(MIN(steepness(ic), 0.5_wp), 0.01_wp)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_orographic_features

  ! =========================================================================================================
  !>
  !> #### Calculation of surface runoff based on the ARNO Scheme
  !>
  !> Only part of the water reaching the surface (rain, snowmelt, ...) infiltrates into the soil. The
  !> remaining part leaves the grid cell as surface runoff. Surface runoff and infiltration flux are
  !> calculated following the ARNO Scheme (compare E. Todini, The ARNO rainfall-runoff model 1996
  !> [DOI](https://doi.org/10.1016/S0022-1694(96)80016-3)). Further information on the implementation
  !> details is given in [JSBACH3 documentation section 2.3.2.2 ](
  !> https://pure.mpg.de/rest/items/item_3279802_23/component/file_3316522/content#subsubsection.2.3.2.2)
  !>
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE arno_scheme(l_infil_subzero, t_soil_sl1, water_to_soil, wpi_rootzone, &
    &                    steepness, wpi_rootzone_max,                              &
    &                    runoff, infilt)

  !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt

    LOGICAL, INTENT(in)  ::        &
      & l_infil_subzero              !< Allow surface water infiltration at temperatures below 0 degree C

    REAL(wp), INTENT(in) ::        &
      & t_soil_sl1,                & !< Temperature of uppermost soil layer [K]
      & water_to_soil,             & !< Amount of water reaching the ground [m /(time step)]
      & wpi_rootzone,              & !< Liquid water + ice within the root zone [m]
      & steepness,                 & !< Parameter defining the subgrid slope distribution []
      & wpi_rootzone_max             !< Maximum water holding capacity (liquid + ice) of the root zone [m]

    REAL(wp), INTENT(out) ::       &
      & runoff,                    & !< Surface runoff [m /(time step)]
      & infilt                       !< Infiltration of water into the soil [m /(time step)]

    ! Local variables
    REAL(wp) ::                    &
      & ws_min_drain,              & !< Minimum amount of root zone soil water for drainage [m]
      & ws_rel,                    & !< Relative root zone soil moisture (water + ice)
      & zb1, zbm, zconw1, zvol       !< Parameters

    ! todo: parameters should be defined as parameters
    zb1    = 1._wp + steepness
    zbm    = 1._wp / zb1
    zconw1 = wpi_rootzone_max * zb1

    ! Surface runoff and infiltration
    ! -----------------------------------
    !   f(w) = 1 - (1 - w/w_max)**b

    ! For very dry soils (relative soil moisture below zwdmin) we assume that all water infiltrates
    ! the soil and there is no runoff.
    ! TODO: rename ws_min_drain to ws_min_runoff
    ws_min_drain = zwdmin * wpi_rootzone_max

    ! Temperature dependence of infiltration
    ! Note: It is recommend to allow infiltration also at temperatures below 0.C (l_infil_subzero=true)
    !       with the multi layer snow scheme. Otherwise all snow melt is going to surface runoff
    !       if tsurf = 0.C and T_soil < 0.C leading to a severe dry bias of the soil.
    IF (.NOT. l_infil_subzero .AND. (t_soil_sl1 < tmelt)) THEN
      ! No infiltration as uppermost soil layer is frozen -> all water goes into runoff
      runoff = water_to_soil
      infilt = 0._wp

    ELSE ! Infiltration is possible

      IF (water_to_soil > 0._wp .AND. wpi_rootzone > ws_min_drain) THEN
        ! There is enough soil moisture and surface runoff is possible.

        ws_rel = MIN(1._wp, wpi_rootzone / wpi_rootzone_max)         ! Relative root zone soil moisture
        zvol   = (1._wp - ws_rel)**zbm - water_to_soil / zconw1      ! Factor accounting for subgrid
                                                                     !  scale inhomogeneity and steepness
        runoff = water_to_soil - (wpi_rootzone_max - wpi_rootzone)   ! runoff > 0: water exceeding soil capacity
        IF (zvol > 0._wp) runoff = runoff + wpi_rootzone_max * zvol**zb1
        runoff = MAX(MIN(runoff, water_to_soil), 0._wp)              ! No runoff < 0, limited by available water
        infilt = water_to_soil - runoff                              ! Remaining water is infiltrated.

      ELSE
        ! No runoff as soil is too dry. All water infiltrates.
        runoff = 0._wp
        infilt = water_to_soil
      END IF
    END IF

  END SUBROUTINE arno_scheme

  ! =========================================================================================================
  !>
  !> #### Calculation of vertical soil water movement for multiple soil layers
  !>
  !> This is the core routine for vertical water transport between the soil layers. The routine starts
  !> with sanity checks for the soil moisture.
  !> Before actually calculating the vertical transport, we check for potential water deficits in case
  !> evapotranspiration exceeds soil moisture and reduce soil evaporation and transpiration used in the
  !> vertical transport scheme accordingly ([[diagnose_evapotrans]]). Hydraulic conductivity and diffusivity
  !> of the soil layers and at the layer interfaces are updated in [[get_soilhyd_properties]].
  !> Only then the vertical movement of mobile water is calculated in subroutine [[soilhyd]].
  !> Finally, several sanity checks are performed and corrections are applied if needed.
  !>
  SUBROUTINE soilhyd(                                                &
    & nc, l_fract, l_pf_soil, lat, lon, nsoil, dtime, enforce_water_budget,  &
    & slope, soilhydmodel, interpol_mean, hydro_scale,                       &
    & model_scheme,                                                          &
    & w_soil_wilt_fract,   &
    & dsoil, wtr_wres_sl, wtr_p_wilt_sl,                                     &
    & wtr_field_cap_sl, wtr_wsat_sl, hyd_cond_sat_sl,                        &
    & vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl,                      &
    & vol_wres_sl,                                                           &
    & bclapp_sl, matric_pot_sl, pore_size_index_sl,                          &
    & wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl,                       &
    & transpiration_in_m, evapo_soil_in_m, root_depth_sl,                    &
    & infilt, runoff_dunne, drain_bot, drainage_sl, wtr_transp_down_sl,      &
    & tpe_overflow, evapo_deficit, wtr_soilhyd_res, hyd_cond_bot,            &
    & ltpe_closed, ltpe_open, ftranspiration_sl                              &
    & )

    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_jsb_physical_constants, ONLY: rhoh2o, rhoi
    USE mo_jsb_impl_constants,     ONLY: WB_LOGGING, WB_ERROR
    USE mo_hydro_constants,        ONLY: Semi_Distributed_, Uniform_, &
      & drain_min, drain_max, drain_exp, k_brock

    !TODO: sequence of arguments is NOT consistent with this list of argument declarations
    INTEGER, INTENT(IN) ::          &
      & nc                            !< Vector length
    LOGICAL, INTENT(IN) ::          &
      & l_fract(:),                 & !< True: grid cell fraction of tile is > 0
      & l_pf_soil(:)                  !< True: permafrost within the soil (above bedrock)
    REAL(wp), INTENT(IN) ::         &
      & lat(:),                     & !< Grid cell latitudes [degree]
      & lon(:)                        !< Grid cell longitudes [degree]
    INTEGER, INTENT(IN) ::          &
      & nsoil,                      & !< Number of soil layers (vertical dimension)
      & enforce_water_budget,       & !< Action in case of water balance issues
      & soilhydmodel,               & !< Model scheme for soil hydraulic properties
      & interpol_mean,              & !< Vertical interpolation scheme for soil hydraulic properties
      & hydro_scale,                & !< Hydrology scale scheme
      & model_scheme                  !< ICON-Land model scheme
    REAL(wp), INTENT(in) ::         &
      & w_soil_wilt_fract,          & !< Relative soil moisture at wilting point []
      & dtime,                      & !< Model time step legth [s]
      & slope(:),                   & !< Slope of gridcell / tile []
      & dsoil(:,:),                 & !< Soil depth until bedrock within each layer [m]
      & wtr_wres_sl(:,:),           & !< Residual water content of soil layer (reduced by ice) [m]
      ! Note: wtr_p_wilt_sl is currently not used. We still keep it, as we are considering updating
      !   the water stress function to use the actual instead of the global uniform values for the
      !   wilting point.
      & wtr_p_wilt_sl(:,:),         & !< Wilting point of the soil layer (reduced by ice) [m]
      & wtr_field_cap_sl(:,:),      & !< Field capacity of the soil layer (reduced by ice) [m]
      & wtr_wsat_sl(:,:),           & !< Saturation capacity of the soil layer (reduced by ice) [m]
      & hyd_cond_sat_sl(:,:),       & !< Hydraulic conductivity of saturated soil [m/s]
      & vol_porosity_sl(:,:),       & !< Volumetric soil porosity [m/m]
      & vol_field_cap_sl(:,:),      & !< Volumetric field capacity [m/m]
      & vol_p_wilt_sl(:,:),         & !< Volumetric permanent wilting point [m/m]
      & vol_wres_sl(:,:),           & !< Volumetric residual water content [m/m]
      & bclapp_sl(:,:),             & !< Exponent B in Clapp and Hornberger
      & matric_pot_sl(:,:),         & !< Saturated soil matric potential [m]
      & pore_size_index_sl(:,:),    & !< Soil pore size distribution index used in Van Genuchten
      & wtr_soil_pot_scool_sl(:,:), & !< Potentially supercooled water [m]
      & root_depth_sl(:,:),         & !< Thicknesses of soil layers until rooting depth [m]
      & transpiration_in_m(:),      & !< Amount of transpiration within timestep [m]
      & evapo_soil_in_m(:)            !< Amount of soil evaporation within timestep [m]
    REAL(wp), INTENT(INOUT) ::      &
      & infilt(:),                  & !< Amount of infiltration within timestep [m]
      & ice_soil_sl(:,:),           & !< Amount of ice in the soil [m]
      & wtr_soil_sl(:,:),           & !< Soil moisture of each layer [m]
      & runoff_dunne(:),            & !< Amount of saturation overflow (Dunne runoff) [m]
      & tpe_overflow(:),            & !< Terraplanet: soil water overflow reservoir [m]
      & evapo_deficit(:)              !< Water evaporated from unintended sources [m]
    REAL(wp), INTENT(OUT) ::        &
      & drain_bot(:),               & !< Drainage towards bedrock [m]
      & drainage_sl(:,:),           & !< Subsurface drainage on soil layers [m]
      & wtr_transp_down_sl(:,:),    & !< Vertical water transport into the next deeper soil layer [m]
      & wtr_soilhyd_res(:),         & !< Residual error of the vertical transport scheme
      & hyd_cond_bot(:)               !< hydraulic conductivity of the bottom soil layer [m/s]
    LOGICAL, INTENT(IN) ::          &
      & ltpe_closed,                & !< Terraplanet setup with closed water balance
      & ltpe_open                     !< Terraplanet setup with additional moisture in deep soil layers
    REAL(wp), OPTIONAL, INTENT(IN) :: &
      & ftranspiration_sl(:,:)        !< Fraction oftranspiration [] ???

    ! Local variables
    INTEGER ::                      &
      & ic,                         & !< Index for grid cells
      & is,                         & !< Index for soil layers
      & last_soil_layer(nc)           !< Index for deepest soil layer above the bedrock

    CHARACTER(len=4096) :: message_text_long !< long string for soil hydrology checks

    REAL(wp) ::                     &
      & weq_soil(nc),               & !< Total column soil moisture (water + ice) [m water equivalent]
      & ws_inter_sl(nc,nsoil),      & !< Soil moisture with infiltration added to first layer [m]
      & wtr_soil_sl_t0(nc,nsoil),   & !< Soil moisture state prior to vertical transport [m]
      & ws_vol_sl(nc,nsoil),        & !< Volumetric soil moisture of soil layers []
      & ice_impedance_sl(nc,nsoil), & !< Ice impedance factor on soil layers []
      & hyd_cond_sl(nc,nsoil),      & !< Hydraulic conductivity of soil layers [m/s]
      & diffus_sl(nc,nsoil),        & !< Diffusivity at mid-depth of soil layers [m2/s]
      & hyd_cond_li(nc,nsoil),      & !< Hydraulic conductivity at upper layer interface [m/s]
      & diffus_li(nc,nsoil),        & !< Diffusivity at upper layer interface [m2/s]
      & evapo_soil_sl1(nc),         & !< Bare soil evaporation extracted from top layer (mobile part) [m]
      & transpiration_sl(nc,nsoil), & !< Transpiration extracted from the soil layers [m]
      & deficit_evapotrans(nc),     & !< Water that could not be extracted from the soil [m]
      & deficit_trans(nc),          & !< Transpiration that could not be extracted from the soil [m]
      & deficit_sevap(nc),          & !< Soil evaporation that could not be extracted from the soil [m]
      & remoist(nc),                & !< Terraplanet: water recycled for soil moisture
      & wpi_wsat_sl(nc,nsoil),      & !< Maximum storage in each soil layer (porosity; not reduced by ice) [m]
      & wpi_field_cap_sl(nc,nsoil), & !< Field capacity of the soil layer (not reduced by ice) [m]
      & wpi_p_wilt_sl(nc,nsoil),    & !< Wilting point (not reduced with respect to soil ice fraction) [m]
      & wpi_wres_sl(nc,nsoil),      & !< Residual water content (not reduced with respect to soil ice fraction) [m]
      & wtr_transp_residual(nc),    & !< Numerical error caused by the vertical transport scheme [m]
      & wtr_transp_corr(nc),        & !< Uncorrected part of the transport residual [m]
      & wtr_residual_sl(nc,nsoil),  & !< Immobile soil water (supercooled or below residual water content) [m]
      & drain_slow_sl(nc,nsoil),    & !< Slow drainage component (below field capacity) [m]
      & drain_fast_sl(nc,nsoil),    & !< Fast drainage component (above field capacity) [m]
      & storage_change_sl(nc,nsoil),& !< Soil storage change term used in vertical water transport (m)
      & hlp1(nc,nsoil),             & !< Helper array
      & hlp4(nc,nsoil),             & !< Helper array
      & hlp2(nc),                   & !< Helper array
      & hlp3(nc),                   & !< Helper array
      & hlp5,                       & !< Helper array
      & wtr_flux_top,               & !< Water flux entering a soil layer from above [m]
      & drain_excess,               & !< Additional drainage due to soil oversaturation [m]
      & ice2weq,                    & !< Conversion factor: [m ice] to [m water equivalent]
      & weq2ice                       !< Conversion factor: [m water equivalent] to [m ice]

    ! Temporary index lists for error detection
    INTEGER :: ie             !< Error cell index
    INTEGER :: n_errors       !< Number of error cells
    INTEGER :: error_idx(nc)  !< Cell indices of error cells
    LOGICAL :: error_flag(nc) !< Flag that given cell index is already in list (to suppress double-counting)

    CHARACTER(len=*), PARAMETER :: routine = modname//':soilhyd'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(weq_soil, ws_inter_sl, wtr_soil_sl_t0, ws_vol_sl, ice_impedance_sl, hyd_cond_sl, diffus_sl) &
    !$ACC   CREATE(hyd_cond_li, diffus_li, evapo_soil_sl1, transpiration_sl, drain_slow_sl, drain_fast_sl) &
    !$ACC   CREATE(deficit_evapotrans, deficit_trans, deficit_sevap, remoist, storage_change_sl) &
    !$ACC   CREATE(wtr_residual_sl, wpi_wsat_sl, wpi_field_cap_sl, wpi_p_wilt_sl, wpi_wres_sl) &
    !$ACC   CREATE(wtr_transp_residual, wtr_transp_corr, last_soil_layer, hlp1, hlp2, hlp3, hlp4)

    ! Conversion factors: [m ice] to [m water equivalent] and vice versa
    ice2weq = rhoi / rhoh2o
    weq2ice = rhoh2o / rhoi

    !-------------------
    ! Sanity checks
    !-------------------
    ! Check for negative soil moisture values and handle according to "enforce_water_budget".
    ! Note: This is not possible on GPUs
#ifndef _OPENACC
    IF (enforce_water_budget /= WB_IGNORE) THEN
      ! Assure, all soil moisture values are positive.
      n_errors = 0
      error_flag(:) = .FALSE.
      DO is = 1, nsoil
        DO ic = 1, nc
          ! Check for negative soil moisture (in cells with a relevant grid cell fraction
          ! that had not yet been counted)
          IF (wtr_soil_sl(ic,is) < 0._wp .AND. l_fract(ic) .AND. .NOT. error_flag(ic)) THEN
            n_errors = n_errors + 1    ! Count the soil moisture errors
            error_idx(n_errors) = ic   ! Keep the respective cell indices
            error_flag(ic) = .TRUE.    ! Set error flag
          END IF
        END DO
      END DO

      ! Finish or write at least an error messages in case errors have been found.
      DO ie = 1, n_errors
        ic = error_idx(ie)
        WRITE (message_text,*) 'negative soil moisture (routine start) at ', &
          & lat(ic), 'N and ', lon(ic), 'E:', NEW_LINE('a'), &
          & 'Soil moisture:    ', wtr_soil_sl(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish ('soilhyd', message_text)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message ('soilhyd', message_text, all_print=.TRUE.)
        END IF
      END DO
    END IF
#endif

    !-------------------
    ! Preparations
    !-------------------

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil
      DO ic = 1, nc

        ! Calculate soil water and/or ice content at saturation, at field capacity, at wilting point
        ! and at residual point. These are volumetric quantities depending on the filling of the
        ! soil pore volume. (The amount of liquid water [m] fitting into the specific volume is
        ! greater than the amount of ice [m water equivalent].)
        wpi_wsat_sl     (ic,is) = vol_porosity_sl (ic,is) * dsoil(ic,is)
        wpi_field_cap_sl(ic,is) = vol_field_cap_sl(ic,is) * dsoil(ic,is)
        wpi_p_wilt_sl   (ic,is) = vol_p_wilt_sl(ic,is)    * dsoil(ic,is)
        wpi_wres_sl     (ic,is) = vol_wres_sl(ic,is)      * dsoil(ic,is)

        ! Supercooled soil water and water below the residual water content are excluded
        ! from vertical water movement.
        wtr_residual_sl(ic,is) = MIN(wtr_soil_sl(ic,is), &
          &                      MAX(wtr_soil_pot_scool_sl(ic,is), wtr_wres_sl(ic,is)))

        ! Initialize subsurface drainage, storage change term and vertical water transport
        drainage_sl(ic,is)        = 0._wp
        drain_slow_sl(ic,is)      = 0._wp
        drain_fast_sl(ic,is)      = 0._wp
        storage_change_sl(ic,is)  = 0._wp
        wtr_transp_down_sl(ic,is) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ! Initialize different runoff and drainage components
      runoff_dunne(ic) = 0._wp
      drain_bot(ic)    = 0._wp
      ! Initialize variable that stores the number of active layers
      last_soil_layer(ic) = 0
    END DO
    !$ACC END PARALLEL LOOP

    ! Find out the lowest "active" soil layer, i.e. with a soil fraction above bedrock
    !----------------------------------------
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (dsoil(ic,is) > 0._wp) last_soil_layer(ic) = is
      END DO
    END DO
    !$ACC END PARALLEL
#ifndef _OPENACC
    IF (ANY(last_soil_layer(:) < 1)) CALL finish('soilhyd', 'Problem with no. of active soil layers (=0)')
#endif

    ! Liquid soil moisture including infiltrating surface water
    ! Water infiltration from the surface is added to the uppermost soil layer.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ws_inter_sl(ic,1) = wtr_soil_sl(ic,1) + infilt(ic)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 2, nsoil
      DO ic = 1, nc
        ws_inter_sl(ic,is) = wtr_soil_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Find out potential water deficits in case evapotranspiration exceeds soil moisture and reduce
    ! soil evaporation and transpiration fluxes intended for the vertical transport scheme accordingly.
    ! Actual evapo- and transpiration ('transpiration_in_m' and 'evapo_soil_in_m') remain unchanged.
    ! Attention: the diagnostic for transpiration needs to be consistent with the water_stress
    ! computation. Besides the assumptions on inaccessible water amounts (below wilting point
    ! or supercooled), this also comprises the treatment of soil ice.
    CALL diagnose_evapotrans(nc, l_fract, lat, lon, nsoil, enforce_water_budget,        &
      &                      model_scheme,                                              &
      &                      w_soil_wilt_fract, root_depth_sl, dsoil, wtr_field_cap_sl, &
      &                      transpiration_in_m, evapo_soil_in_m, wtr_soil_sl,          &
      &                      evapo_deficit, evapo_soil_sl1, transpiration_sl,           &
      &                      deficit_sevap, deficit_trans, ftranspiration_sl)

    ! Determine hydraulic conductivity and diffusivity on layers and at layer interfaces
    !   Note that the indexing uses the upper interface -> interface Ii = I(i,i-1)
    !   For numerical reasons, soil loops always cover all soil layers, even if they are not active (bedrock)
    CALL get_soilhyd_properties(                                            &
          & soilhydmodel, interpol_mean,                                    &
          & nc, nsoil, dsoil(:,:),                                          &
          & ws_inter_sl(:,:), ice_soil_sl(:,:), wtr_soil_pot_scool_sl(:,:), &
          & wpi_wsat_sl(:,:), wpi_wres_sl(:,:), hyd_cond_sat_sl(:,:),       &
          & matric_pot_sl(:,:), bclapp_sl(:,:), pore_size_index_sl(:,:),    &
          & dt=dtime, last_soil_layer=last_soil_layer(:),                   &
          & ice_impedance=ice_impedance_sl(:,:),                            &
          & K=hyd_cond_sl(:,:),       D=diffus_sl(:,:),                     &
          & K_inter=hyd_cond_li(:,:), D_inter=diffus_li(:,:)                &
          )

    ! ----------------------
    !  Subsurface drainage
    ! ----------------------
    IF (hydro_scale == Semi_Distributed_) THEN

      ! Calculate drainage from each soil layer
      ! Note: the semi_distributed case only computes soil layer drainage but no bottom drainage
      !   because it is no part of the original ARNO scheme

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (wpi_wsat_sl(ic,is)      > wpi_field_cap_sl(ic,is) .AND. &
            & wpi_field_cap_sl(ic,is) > wtr_residual_sl(ic,is)) THEN

            ! Slow drainage resulting from water between residual and field capacity fractions
            drain_slow_sl(ic,is) = (drain_min * dtime)                                  &
              & * MIN(1._wp, MAX(0._wp, (wtr_soil_sl(ic,is) - wtr_residual_sl(ic,is)))) &
              & / (wpi_field_cap_sl(ic,is) - wtr_residual_sl(ic,is))

            ! Fast drainage resulting from water above the field capacity
            drain_fast_sl(ic,is) = (drain_max - drain_min) * dtime &
              & * (MIN(1._wp, MAX(0._wp, (                         &
              &    wtr_soil_sl(ic,is) - wpi_field_cap_sl(ic,is)))) &
              & / (wpi_wsat_sl(ic,is) - wpi_field_cap_sl(ic,is)))**drain_exp
            drain_fast_sl(ic,is) = MAX(0._wp, MIN(drain_fast_sl(ic,is), &
              & wtr_soil_sl(ic,is) - wpi_field_cap_sl(ic,is)))

            ! Add drainage components and apply ice impedance
            drainage_sl(ic,is) = ice_impedance_sl(ic,is) * (drain_slow_sl(ic,is) + drain_fast_sl(ic,is))
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

    ELSE IF (hydro_scale == Uniform_) THEN

      ! Compute drainage from active soil layers and bottom drainage from lowest layer

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (is < last_soil_layer(ic)) THEN
            ! Calculate drainage for all soil layers above lowest hydrologically active layer
            drainage_sl(ic,is) = MIN(hyd_cond_sl(ic,is) * SIN(slope(ic)*pi/2._wp) * dtime, &
              &                      MAX(0._wp, wtr_soil_sl(ic,is)-wtr_residual_sl(ic,is)))

          ELSE IF (is == last_soil_layer(ic)) THEN
            ! Prescribe drainage from hydraulic conductivity of the lowest layer [m s-1] --> [m]
            ! or assumed hydraulic conductivity of (fractured) bedrock
            IF (l_pf_soil(ic)) THEN
              drain_bot(ic) = 0._wp    ! No bottom layer drainage in case of permafrost
            ELSE
              drain_bot(ic) = MIN(k_brock, hyd_cond_sl(ic,last_soil_layer(ic))) * dtime
            END IF

          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)

    ! ----------------------------------------------
    !  Initialize water transport conservation test
    ! ----------------------------------------------
    ! First part: We store the actual soil water content including the amount of infiltration and
    ! evaporation of this time step. We will re-calculate the amount of soil water below, after the
    ! vertical transport calculation to find out the residual (resulting from computational
    ! precision errors).

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      ! Sum up surface fluxes of this time step (for water conservation test)
      wtr_transp_residual(ic) = infilt(ic) + evapo_soil_sl1(ic)
      ! Hydrological conductivity of the lowest soil layer
      hyd_cond_bot(ic) = hyd_cond_sl(ic,last_soil_layer(ic))
    END DO

    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        ! Store soil moisture state for later diagnosis of vertical water fluxes
        wtr_soil_sl_t0(ic,is) = wtr_soil_sl(ic,is)

        ! Sum up surface fluxes calculated above, this layer's water content and transpiration from
        ! the layer to get the total actual soil moisture (for water conservation test).
        wtr_transp_residual(ic) = wtr_transp_residual(ic) + &
            & wtr_soil_sl(ic,is) + transpiration_sl(ic,is) - drainage_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL

    ! ----------------------------------------------------
    !  Vertical movement of mobile water through the soil
    ! ----------------------------------------------------

    ! Prepare the arguments needed to call subroutine calc_vertical_transport:
    !   Note that normalized storage values (volumetric soil moisture and storage change) are
    !   needed for the transport and immobile water is temporarily removed from the soil.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (dsoil(ic,is) > 0._wp) THEN
          ws_vol_sl(ic,is)         = (wtr_soil_sl(ic,is) - wtr_residual_sl(ic,is))  / dsoil(ic,is)
          storage_change_sl(ic,is) = (transpiration_sl(ic,is) - drainage_sl(ic,is)) / dsoil(ic,is)
        END IF
        hlp1(ic,is) = 1._wp
        hlp4(ic,is) = (storage_change_sl(ic,is)) / dtime ! Normalized storage change term
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      hlp2(ic) = (infilt(ic) + evapo_soil_sl1(ic)) / dtime ! Upper boundary condition
      hlp3(ic) = drain_bot(ic) / dtime                     ! Lower boundary condition
    END DO
    !$ACC END PARALLEL LOOP

    hlp5 = 1._wp   !TODO: Check if this is options%alpha and if we should use it here.

    CALL calc_vertical_transport( &
                          ! in
                        & dtime, hlp5,                            & ! Time step, alpha
                        & last_soil_layer(1:nc),                  & ! Number of active soil layers above bedrock
                        & dsoil(1:nc,1:nsoil),                    & ! Soil layer thickness
                        & hlp2(1:nc),                             & ! Upper boundary: Infiltration and evaporation flux
                        & hlp3(1:nc),                             & ! Lower boundary: Bottom layer drainage flux
                        & hlp4(1:nc,1:nsoil),                     & ! Normalized storage change term (Transpiration)
                        & hlp1(1:nc,1:nsoil),                     & ! Should be 1 for water transport
                        & diffus_li(1:nc,1:nsoil),                & ! Diffusivity
                        & hyd_cond_li(1:nc,1:nsoil),              & ! Hydraulic conductivity
                          ! inout
                        & ws_vol_sl(1:nc,1:nsoil)                 & ! Normalized soil moisture
                        & )

    ! Converting back to absolute water content and add residual water again
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (dsoil(ic,is) > 0._wp) THEN
          wtr_soil_sl(ic,is) = ws_vol_sl(ic,is) * dsoil(ic,is) + wtr_residual_sl(ic,is)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! ----------------------------------
    !  Final testing
    ! ----------------

    ! Assure there is no negative soil moisture

    ! In case the soil moisture of the lowest soil layer becomes negative, we reduce the bottom
    ! drainage accordingly. In rare cases - if the bottom layer is very thin - also above layers
    ! might have a negative soil moisture. Generally, this only happens right after initialization.
    !
    ! But before applying the bottom drainage correction, we diagnose the problem and print a
    ! corresponding message - depending  on namelist parameter enforce_water_budget.
#ifndef _OPENACC
    IF (enforce_water_budget /= WB_IGNORE) THEN
      n_errors = 0
      error_flag(:) = .FALSE.
      DO is = nsoil, 1, -1
        DO ic = 1, nc
          ! Check for negative soil moisture in lowest and above soil layers (in cells with a
          ! relevant land fraction that had not been counted before) and count affected cells.
          IF (wtr_soil_sl(ic, last_soil_layer(ic)) < 0._wp .AND. l_fract(ic) .AND. .NOT. error_flag(ic)) THEN
            IF (is <= last_soil_layer(ic) .AND. wtr_soil_sl(ic,is) < 0._wp) THEN
              n_errors = n_errors + 1
              error_idx(n_errors) = ic
              error_flag(ic) = .TRUE.
            END IF
          END IF
        END DO
      END DO

      ! Write error messages in case cells with negative soil moisture have been found.
      DO ie = 1, n_errors
        ic = error_idx(ie)

        WRITE (message_text_long,*) 'Soil moisture correction needed at ',   &
          &  lat(ic), 'N and ', lon(ic), 'E:' ,      NEW_LINE('a'),          &
          &  'Soil moisture:   ', wtr_soil_sl(ic,:), NEW_LINE('a'),          &
          &  'Bottom drainage: ', drain_bot(ic)
        CALL message('soilhyd', message_text_long, all_print=.TRUE.)
      END DO
    END IF
#endif

    ! Apply botton layer drainage correction
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = nsoil, 1, -1
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (wtr_soil_sl(ic, last_soil_layer(ic)) < 0._wp .AND. l_fract(ic)) THEN
          IF (is <= last_soil_layer(ic) .AND. wtr_soil_sl(ic,is) < 0._wp) THEN
            drain_bot(ic) = MAX((drain_bot(ic) + wtr_soil_sl(ic,is)), 0._wp)
            wtr_soil_sl(ic,is) = 0._wp
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)

    ! Water transport conservation test
    ! ----------------------------------
    ! Second part: We subtract the current total soil moisture from the state calculated above,
    ! prior to the vertical transport. The difference is the vertical transport residual. It will
    ! be diagnosed corrected below.
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        wtr_transp_residual(ic) = wtr_transp_residual(ic) - wtr_soil_sl(ic,is)
      END DO
    END DO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      wtr_transp_residual(ic) = wtr_transp_residual(ic) - drain_bot(ic)
      wtr_soilhyd_res(ic) = wtr_transp_residual(ic)
    END DO

    ! ------------------------------------
    !  Handling of water transport errors
    ! ------------------------------------

#ifndef _OPENACC
    ! 1. Complain if the overall transport error is too large
    IF (enforce_water_budget /= WB_IGNORE) THEN
      n_errors = 0
      DO ic = 1, nc
        IF (ABS(wtr_transp_residual(ic)) > 1.0e-10_wp .AND. l_fract(ic)) THEN
          n_errors = n_errors + 1
          error_idx(n_errors) = ic
        END IF
      END DO

      DO ie = 1, n_errors
        ic = error_idx(ie)
        WRITE (message_text_long,*) 'Soil water transport residual too large at ', &
          & lat(ic), 'N and ', lon(ic), 'E:',                NEW_LINE('a'), &
          & 'Transport residual: ', wtr_transp_residual(ic), NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(ic),              NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(ic),      NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_bot(ic),           NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(ic,:),  NEW_LINE('a'), &
          & 'Layer drainage:     ', drainage_sl(ic,:),       NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish ('soilhyd', message_text_long)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message ('soilhyd', message_text_long, all_print=.TRUE.)
        END IF
      END DO
    END IF
#endif

    ! 2. Modify lower boundary condition (drainage) as much as possible
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      IF (ABS(wtr_transp_residual(ic)) > 0._wp) THEN
        wtr_transp_corr(ic)     = MIN(wtr_transp_residual(ic) + drain_bot(ic), 0._wp)
        drain_bot(ic)           = MAX(drain_bot(ic) + wtr_transp_residual(ic), 0._wp)
        wtr_transp_residual(ic) = wtr_transp_corr(ic)
      END IF
    END DO

    !$ACC END PARALLEL

    ! 3. Modify soil layers if necessary
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = nsoil, 1, -1
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (ABS(wtr_transp_residual(ic)) > 0._wp) THEN
          ! Make sure to correct the residual and not to put negative soil moisture onto the residual
          IF (wtr_soil_sl(ic,is) >= 0._wp .OR. wtr_soil_sl(ic,is) + wtr_transp_residual(ic) > 0._wp) THEN
            wtr_transp_corr(ic)     = MIN(wtr_transp_residual(ic) + wtr_soil_sl(ic,is), 0._wp)
            wtr_soil_sl(ic,is)       = MAX(wtr_soil_sl(ic,is) + wtr_transp_residual(ic), 0._wp)
            wtr_transp_residual(ic) = wtr_transp_corr(ic)
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    ! 4. Complain if correction is not sufficient
    IF (enforce_water_budget /= WB_IGNORE) THEN
      n_errors = 0
      DO ic = 1, nc
        IF (ABS(wtr_transp_residual(ic)) > 1.0e-10_wp .AND. l_fract(ic)) THEN
          n_errors = n_errors + 1
          error_idx(n_errors) = ic
        END IF
      END DO

      DO ie = 1, n_errors
        ic = error_idx(ie)
        WRITE (message_text_long,*) 'Cannot correct soil water transport residual at ', &
          & lat(ic),'N and ',lon(ic),'E:',                   NEW_LINE('a'), &
          & 'Transport residual: ', wtr_transp_residual(ic), NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(ic),              NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(ic),      NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_bot(ic),           NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(ic,:),  NEW_LINE('a'), &
          & 'Layer drainage:     ', drainage_sl(ic,:),       NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish ('soilhyd', message_text_long)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message ('soilhyd', message_text_long, all_print=.TRUE.)
        END IF
      END DO
    END IF
#endif

    ! Diagnose vertical water transport
    !-----------------------------------
    ! by comparing the current soil moisture state and fluxes with the state saved at the beginning
    ! of the routine, before the vertical transport.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR PRIVATE(wtr_flux_top)
      DO ic = 1, nc
        IF (is == 1) THEN
          ! We start from the top layer: vertical transport is fed from the surface fluxes.
          wtr_flux_top = infilt(ic) + evapo_soil_sl1(ic)
        ELSE IF (is <= last_soil_layer(ic)) THEN
          ! In the layers below we get the inflow from the downward transport of the above layer
          ! (just calculated in the lines below).
          wtr_flux_top = wtr_transp_down_sl(ic,is-1)
        END IF
        IF (is <= last_soil_layer(ic)) THEN
          ! Downward transport is calculated from the difference to the state before the transport
          ! calculations, the flux from above and the loss from transpiration.
          wtr_transp_down_sl(ic,is) = wtr_soil_sl_t0(ic,is) - wtr_soil_sl(ic,is) + wtr_flux_top + transpiration_sl(ic,is)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL


    !-----------------------------------------------------------------
    ! Final corrections in case of soil water deficit or excess water
    !-----------------------------------------------------------------

    ! We now account for the water deficit diagnosed in 'diagnose_evapotrans', which
    ! also includes possible water deficits from the surface hydrology calculations.

    ! We extract the bare soil evaporation deficit from the top layer (using mobile and immobile soil water,
    ! because this has to be consistent with the computation of the upper soil layer relative humidity
    ! in sse_processes:relative_humidity_soil).
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      IF (deficit_sevap(ic) < 0._wp) THEN                              ! Remaining soil evaporation demand
        wtr_soil_sl(ic,1) = wtr_soil_sl(ic,1) + deficit_sevap(ic)      ! Take water from uppermost soil layer.
        IF (wtr_soil_sl(ic,1) < 0._wp) THEN                            ! If there is not enough liquid soil water,
          deficit_sevap(ic) = wtr_soil_sl(ic,1)                        ! this is the amount of water missing.
          wtr_soil_sl(ic,1) = 0._wp                                    ! Soil moisture cannot be negative.
        ELSE                                                           ! No deficit remains.
          deficit_sevap(ic) = 0._wp
        END IF
      END IF
      evapo_deficit(ic) = evapo_deficit(ic) + deficit_sevap(ic)        ! Remaining demand is added to deficit
    END DO

    ! Correct evapotranspiration deficit residual if existent
    ! @todo this is an unfortunate necessity as sometimes the evaporative demand is not fully
    !       backed by the water available in the land surface reservoirs.
    !       Changing the soil ice content here obviously violates the energy balance. A correction
    !       term is needed for this, if no better solution can be found to avoid the evaporation deficit.
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      weq_soil(ic) = 0._wp   ! Initialization for diagnostic below
      deficit_evapotrans(ic) = deficit_sevap(ic) + deficit_trans(ic)
    END DO

    ! Diagnose total amount of soil water and ice
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        weq_soil(ic) = weq_soil(ic) + wtr_soil_sl(ic,is) + ice_soil_sl(ic,is) * ice2weq
      END DO
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    IF (ANY(weq_soil(:) + deficit_evapotrans(:) < -EPSILON(1.0_wp) .AND. l_fract(:))) THEN
      IF (enforce_water_budget == WB_ERROR .OR. enforce_water_budget == WB_LOGGING) THEN
        WRITE (message_text,*) 'ET deficit cannot be compensated with the soil moisture storage'
        CALL message ('soilhyd', message_text, all_print=.TRUE.)
      END IF
    END IF
#endif

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (deficit_evapotrans(ic) < 0._wp .AND. weq_soil(ic) > 0._wp) THEN
          ! reduce soil moisture and ice relative to the water content of each layer
          wtr_soil_sl(ic,is) = MAX(0._wp, wtr_soil_sl(ic,is) + &
            & deficit_evapotrans(ic) * wtr_soil_sl(ic,is) / weq_soil(ic))
          ice_soil_sl(ic,is) = MAX(0._wp, ice_soil_sl(ic,is) + &
            & (deficit_evapotrans(ic) * (ice_soil_sl(ic,is) * ice2weq) / weq_soil(ic)) * weq2ice )
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    ! Additional corrections for excess soil moisture

    IF (hydro_scale == Uniform_) THEN
      ! Initial estimates of layered drainage are too low if soils were dry at the beginning of time step.
      ! Hence, add water to layered drainage according to k_sat - k_act
      ! In case of last soil layer only k_sat since no layered drainage was calculated before
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          ! If saturation moisture is exceeded, add exess water to layered drainage (according to ksat)
          IF (wtr_soil_sl(ic,is) > wtr_wsat_sl(ic,is)) THEN
            IF (is < last_soil_layer(ic)) THEN
              drainage_sl(ic,is) = drainage_sl(ic,is)                                                   &
                & + MIN(wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is), SIN(slope(ic)*pi/2._wp) * dtime        &
                & * MAX(0._wp,(hyd_cond_sat_sl(ic,is) * ice_impedance_sl(ic,is) - hyd_cond_sl(ic,is))))
              wtr_soil_sl(ic,is) = wtr_soil_sl(ic,is)                                                   &
                & - MIN(wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is), SIN(slope(ic)*pi/2._wp) * dtime        &
                & * MAX(0._wp,(hyd_cond_sat_sl(ic,is) * ice_impedance_sl(ic,is) - hyd_cond_sl(ic,is))))
            ELSE
              drainage_sl(ic,is) = drainage_sl(ic,is)                                                   &
                & + MIN(wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is), SIN(slope(ic)*pi/2._wp) * dtime        &
                & * hyd_cond_sat_sl(ic,is) * ice_impedance_sl(ic,is))
              wtr_soil_sl(ic,is) = wtr_soil_sl(ic,is)                                                   &
                & - MIN(wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is), SIN(slope(ic)*pi/2._wp) * dtime        &
                & * hyd_cond_sat_sl(ic,is) * ice_impedance_sl(ic,is))
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    ! In case of over-saturation, increase layer drainage to maximum and, if needed,
    ! pile water upwards
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = nsoil, 2, -1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        ! If saturation moisture is exceeded, then ...
        IF (wtr_soil_sl(ic,is) > wtr_wsat_sl(ic,is)) THEN
          ! first compute additional drainage as if soil was already saturated from the start
          drain_excess = MIN(wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is), MAX(0._wp, &
            &  ((drain_max - drain_min) * dtime - drain_fast_sl(ic,is)) * ice_impedance_sl(ic,is)))
          wtr_soil_sl(ic,is) = wtr_soil_sl(ic,is) - drain_excess
          drainage_sl(ic,is) = drainage_sl(ic,is) + drain_excess
          IF (wtr_soil_sl(ic,is) > wtr_wsat_sl(ic,is)) THEN
            ! and if the soil is still oversaturated, pile moisture upwards, as it should not
            ! have percolated in the first place
            wtr_soil_sl(ic,is-1) = wtr_soil_sl(ic,is-1) + (wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is))
            wtr_soil_sl(ic,is)   = wtr_wsat_sl(ic,is)
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ! Top layer saturation excess is added to the surface runoff (Dunne runoff)
      IF (wtr_soil_sl(ic,1) > wtr_wsat_sl(ic,1)) THEN
        runoff_dunne(ic)  = wtr_soil_sl(ic,1) - wtr_wsat_sl(ic,1)
        wtr_soil_sl(ic,1) = wtr_wsat_sl(ic,1)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !------------------------------------------
    ! Modifications for the Terraplanet setups without global water conservation (TPE open)
    ! or with global water conservation (TPE closed)

    IF (ltpe_open) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          IF (is == nsoil-1 .OR. is == nsoil) THEN
            IF (dsoil(ic,is) > 0._wp) THEN
              ! Keep the two lowest soil layers always very wet
              wtr_soil_sl(ic,is) = MAX(wtr_soil_sl(ic,is), wtr_wsat_sl(ic,is) * 0.9_wp)
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    ELSE IF (ltpe_closed) THEN
      ! no drainage for TPE closed case --> redirect overflow from drainage to overflow pool
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          tpe_overflow(ic) = tpe_overflow(ic) + drainage_sl(ic,is)
        END DO
      END DO

      !$ACC LOOP SEQ
      DO is = nsoil, 1, -1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          ! refill soil moisture from overflow pool
          remoist(ic)        = MIN(tpe_overflow(ic), wtr_field_cap_sl(ic,is) - wtr_soil_sl(ic,is))
          wtr_soil_sl(ic,is) = wtr_soil_sl(ic,is) + remoist(ic)
          tpe_overflow(ic)   = tpe_overflow(ic) - remoist(ic)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

#ifndef _OPENACC
    !------------------------------------------
    !  Check value range for soil moisture
    ! -------------------------------------
    ! In case of errors, handle according to hydro namelist key "enforce_water_balance".


    ! Make sure there are no negative soil moisture values on a tile with non-zero fraction.
    n_errors = 0
    error_flag(:) = .FALSE.

    DO is = 1, nsoil
      DO ic = 1, nc
        IF (wtr_soil_sl(ic,is) < -EPSILON(1.0_wp) .AND. l_fract(ic)) THEN
          IF (.NOT. error_flag(ic)) THEN
            n_errors = n_errors + 1
            error_idx(n_errors) = ic
            error_flag(ic) = .TRUE.
          END IF
        ELSE IF (wtr_soil_sl(ic,is) < 0._wp .AND. l_fract(ic)) THEN
          wtr_soil_sl(ic,is) = 0._wp
        END IF
      END DO
    END DO

    IF (enforce_water_budget /= WB_IGNORE) THEN
      DO ie = 1, n_errors
        ic = error_idx(ie)

        WRITE (message_text_long,*) 'negative soil moisture (routine end) at ', &
          & lat(ic),'N and ',lon(ic),'E:',                NEW_LINE('a'), &
          & 'Infiltration:     ', infilt(ic),             NEW_LINE('a'), &
          & 'Soil evaporation: ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Bottom drainage:  ', drain_bot(ic),          NEW_LINE('a'), &
          & 'Transpiration:    ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Layer drainage:   ', drainage_sl(ic,:),      NEW_LINE('a'), &
          & 'Soil moisture:    ', wtr_soil_sl(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish ('soilhyd', message_text_long)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message ('soilhyd', message_text_long, all_print=.TRUE.)
        END IF
      END DO
    END IF

    ! Make sure soil moisture never exceeds saturation capacity on a tile with non-zero fraction.
    n_errors = 0
    error_flag(:) = .FALSE.

    DO is = 1, nsoil
      DO ic = 1, nc
        IF (wtr_soil_sl(ic,is) > wtr_wsat_sl(ic,is) + EPSILON(1.0_wp) .AND. l_fract(ic) .AND. .NOT. error_flag(ic)) THEN
          n_errors = n_errors + 1
          error_idx(n_errors) = ic
        END IF
      END DO
    END DO

    IF (enforce_water_budget /= WB_IGNORE) THEN
      DO ie = 1, n_errors
        ic = error_idx(ie)

        WRITE (message_text_long,*) 'soil moisture exceeds saturation capacity(routine end) at ', &
          &  lat(ic),'N and ',lon(ic),'E:',                 NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(ic),             NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_bot(ic),          NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Layer drainage:     ', drainage_sl(ic,:),      NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(ic,:),      NEW_LINE('a'), &
          & 'Saturation capcity: ', wtr_wsat_sl(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish ('soilhyd', message_text_long)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message ('soilhyd', message_text_long, all_print=.TRUE.)
        END IF
      END DO
    END IF
#else
    ! On GPUs negative soil moisture values are just set to zero, without message or recording.
    ! Soil water exceeding saturation capacity is completely ignored.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is=1,nsoil
      DO ic=1,nc
        IF (wtr_soil_sl(ic,is) < 0._wp) THEN
          wtr_soil_sl(ic,is) = 0._wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP
#endif

    !$ACC END DATA

  END SUBROUTINE soilhyd

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate lateral water fluxes
  !>
  !> This routine is part of the uniform scale soil hydrology, which is used with HydroTiles or in
  !> site level configurations.
  !>
  !> We regard here the lateral flow of water within each soil layer of grid cell.
  !
  SUBROUTINE calc_soilhyd_lateral(                                          &
    & l_pf_soil, nc, nsoil, delta_time,                                     &
    & flowlag, hyd_cond_bot, soil_depth_sl, vol_field_cap_sl,               &
    & ice_soil_sl, wtr_soil_sl, wtr_latflow_res_sl, drainage_sl, drain_bot, &
    & wtr_latflow_sl)

    USE mo_hydro_constants,    ONLY: k_brock

    LOGICAL, INTENT(IN) ::          &
      & l_pf_soil(:)                  !< True: Permafrost within the soil above bedrock
    INTEGER,  INTENT(IN) ::         &
      & nc,                         & !< Number of cells
      & nsoil                         !< Number of soil layers
    REAL(wp), INTENT(IN) ::         &
      & delta_time,                 & !< Time step length [s]
      & flowlag(:),                 & !< Lag factor determining the outflow []
      & hyd_cond_bot(:),            & !< Hydraulic conductivity of the bottom soil layer [m/s]
      & soil_depth_sl(:,:),         & !< Thicknesses of soil layers until bedrock [m]
      & vol_field_cap_sl(:,:),      & !< Volumetric field capacity []
      & ice_soil_sl(:,:)              !< Ice content [m]
    REAL(wp), INTENT(INOUT) ::      &
      & wtr_soil_sl(:,:),           & !< Soil moisture of each layer [m]
      & wtr_latflow_res_sl(:,:),    & !< Water content in lateral flow reservoir [m]
      & drainage_sl(:,:),           & !< Drainage on layers [m /(time step)]
      & drain_bot(:)                  !< Bottom layer drainage [m /(time step)]
    REAL(wp), INTENT(OUT) ::        &
      & wtr_latflow_sl(:,:)           !< Lateral flux to downstream unit / river [m /(time step)]

    ! Local variables
    INTEGER ::                      &
      & ic,                         & !< Grid ell index
      & is                            !< Soil layer index
    REAL(wp) ::                     &
      & drain_bot_act,              & !< Contribution of layer to bottom layer drainage [m /(time step)]
      & drain_bot_pot(nc),          & !< Potential bottom layer drainage from column [m /(time step)]
      & weq_sl(nc,nsoil),           & !< Total water content of soil layer [m]
      & field_cap_sl(nc,nsoil)        !< Field capacity of soil layer [m]

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(drain_bot_pot, weq_sl, field_cap_sl)

    ! Initialization
    !----------------
    drain_bot_act = 0._wp

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil
      DO ic = 1, nc
        wtr_latflow_sl(ic,is) = 0._wp
        weq_sl(ic,is)         = ice_soil_sl(ic,is) + wtr_soil_sl(ic,is)
        field_cap_sl(ic,is)   = vol_field_cap_sl(ic,is) * soil_depth_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Determine the potential amount of bottom layer drainage, i.e. the potential loss from the
    ! intermediate lateral subsurface drainage reservoir, if this reservoir was sufficiently filled.

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      IF (l_pf_soil(ic)) THEN
        ! No bottom layer drainage in case of permafrost
        drain_bot_pot(ic) = 0._wp
      ELSE
        ! The drainage is either limited by the conductivity of the bedrock or of the lowest soil layer.
        drain_bot_pot(ic) = MIN(k_brock, hyd_cond_bot(ic)) * delta_time
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !> The lateral flow reservoir is filled with drainage from the above soil layer. Until field capacity
    !> is reached, the water will remoist the soil layer. Additional water is assumed to directly drain
    !> to the bedrock, bottom layer drainage is however limited by the conductivity of the bedrock
    !> interface. The remaining water feeds the lateral flux to a downstream tile or to the rivers.

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is=1,nsoil
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        ! Adding layered drainage to storage reservoir and setting drainage to 0.
        wtr_latflow_res_sl(ic,is)     = wtr_latflow_res_sl(ic,is) + drainage_sl(ic,is)
        drainage_sl(ic,is)            = 0._wp
        ! Allowing to remoist the soil in case soil moisture and ice are below field capacity.
        IF (weq_sl(ic,is) < field_cap_sl(ic,is)) THEN
          wtr_soil_sl(ic,is)          = wtr_soil_sl(ic,is)        &
            &                         + MIN(field_cap_sl(ic,is) - weq_sl(ic,is), wtr_latflow_res_sl(ic,is))
          wtr_latflow_res_sl(ic,is)   = wtr_latflow_res_sl(ic,is) &
            &                         - MIN(field_cap_sl(ic,is) - weq_sl(ic,is), wtr_latflow_res_sl(ic,is))
        END IF
        ! Calculating additional bottom layer drainage assuming water can directly move to the bedrock-boundary.
        IF (drain_bot_pot(ic) > 0._wp .AND.  wtr_latflow_res_sl(ic,is) > 0._wp) THEN
          drain_bot_act               = MIN(drain_bot_pot(ic), wtr_latflow_res_sl(ic,is))
          drain_bot_pot(ic)           = drain_bot_pot(ic) - drain_bot_act
          wtr_latflow_res_sl(ic,is)   = wtr_latflow_res_sl(ic,is) - drain_bot_act
          drain_bot(ic)               = drain_bot(ic) + drain_bot_act
        END IF
        ! Calculating the lateral flux to downstream unit / river
        IF (wtr_latflow_res_sl(ic,is) > 0._wp) THEN
          wtr_latflow_sl(ic,is)     = wtr_latflow_res_sl(ic,is) * flowlag(ic)
          ! Removing outflow from reservoir
          wtr_latflow_res_sl(ic,is) = wtr_latflow_res_sl(ic,is) - wtr_latflow_sl(ic,is)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE calc_soilhyd_lateral

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate lateral surface water fluxes
  !>
  !> This routine is part of the uniform scale soil hydrology, which is used with HydroTiles or in
  !> site level configurations.
  !>
  !> We regard the lateral surface water flow within a grid cell: Surface runoff feeds the
  !> the lateral flow reservoir of surface water. The flow to the downstream tile or to the
  !> rivers depends on the assumed retention time for surface runoff (compare namelist parameter
  !> [[t_hydro_config:ret_macro_srf]].

#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE surfhyd_lat(flowlag, runoff, wtr_latflow_res_srf, wtr_latflow_srf)

  !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN)            :: flowlag               !< Factor accounting for retention in surface flow []
    REAL(wp), INTENT(INOUT)         :: runoff                !< Surface runoff [m /(time step)]
    REAL(wp), INTENT(INOUT)         :: wtr_latflow_res_srf   !< Water content in lateral flow reservoir [m]
    REAL(wp), INTENT(OUT)           :: wtr_latflow_srf       !< Lateral flux to downstream unit / river [m /(time step)]

    ! Adding surface runoff to storage reservoir and setting runoff to 0.
    wtr_latflow_res_srf = wtr_latflow_res_srf + runoff
    runoff              = 0._wp
    ! Calculating the lateral flux to downstream unit / river
    wtr_latflow_srf     = wtr_latflow_res_srf * flowlag
    wtr_latflow_res_srf = wtr_latflow_res_srf - wtr_latflow_srf

  END SUBROUTINE surfhyd_lat

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Diagnose evaporation and transpiration fluxes
  !>
  !> This routine tries to distribute the evaporative demand of soil evaporation and transpiration
  !> between the soil layers of the root zone - without actually affecting soil moisture or the
  !> fluxes.
  !> In case there is not enough available root zone water we diagnose the water deficit and
  !> accordingly reduce soil evaporation and transpiration fluxes that are later used in the
  !> vertical transport routine [[calc_vertical_transport]], where the soil moisture is updated.
  !> The actual fluxes ('transpiration_in_m' and 'evapo_soil_in_m') remain unchanged. The water
  !> deficit will be accounted for later at the end of routine [[soilhyd]].
  !>
  !> Note, that only negative evapotranspiration fluxes (actual evapotranspiration) are considered
  !> here as positive fluxes (dew) are considered in the 'calc_surface_hydrology' routines.
  !>
  !> The routine performs several sanity checks, and - depending on
  !> [[t_hydro_config:enforce_water_budget]] - might stop model execution in case of negative soil
  !> moisture or other water balance issues.
  !>
  SUBROUTINE diagnose_evapotrans( &
    & nc,                         &
    & l_fract,                    &
    & lat,                        &
    & lon,                        &
    & nsoil,                      &
    & enforce_water_budget,       &
    & model_scheme,               &
    & w_soil_wilt_fract,          &
    & droot,                      &
    & dsoil,                      &
    & wtr_field_cap_sl,           &
    & transpiration_in_m,         &
    & evapo_soil_in_m,            &
    & wtr_soil_sl,                &
    & evapo_deficit,              &
    & evapo_soil_sl1,             &
    & transpiration_sl,           &
    & deficit_sevap,              &
    & deficit_trans,              &
    & ftranspiration_sl           &
    & )

    USE mo_jsb_model_class,    ONLY: MODEL_QUINCY

    INTEGER,  INTENT(in) ::        &
      & nc                           !< Vector length
    LOGICAL,  INTENT(in) ::        &
      & l_fract(:)                   !< Tile fraction is grater than 0.
    REAL(wp), INTENT(in) ::        &
      & lat(:),                    & !< Grid cell latitudes
      & lon(:)                       !< Grid cell longitudes
    INTEGER,  INTENT(in) ::        &
      & nsoil,                     & !< Number of soil layers (vertical grid dimension)
      & enforce_water_budget,      & !< Switch defining consequence of water balance issues
      & model_scheme                 !< Model scheme (e.g. jsbach / quincy)
    REAL(wp), INTENT(in) ::        &
      & w_soil_wilt_fract,         & !< Relative soil moisture at wilting point []
      & droot(:,:),                & !< Root depth within layer [m]
      & dsoil(:,:),                & !< Soil depth (until bedrock) within layer [m]
      & wtr_field_cap_sl(:,:),     & !< Field capacity of the layer (reduced by ice) [m]
      & transpiration_in_m(:),     & !< Transpiration [m /(time step)]
      & evapo_soil_in_m(:),        & !< Bare soil evaporation [m /(time step)]
      & wtr_soil_sl(:,:)             !< Water content of the soil layer [m]

    REAL(wp), INTENT(inout) ::     &
      & evapo_deficit(:)             !< Evaporation from unintended sources [m]

    REAL(wp), INTENT(out) ::       &
      & evapo_soil_sl1(:),         & !< Evaporation from the top soil layer [m/(time step)]
      & transpiration_sl(:,:),     & !< Transpiration from soil layers [m/(time step)]
      & deficit_trans(:),          & !< Unaccounted transpiration [m/(time step)]
      & deficit_sevap(:)             !< Unaccounted soil evaporation [m/(time step)]

    REAL(wp), OPTIONAL, INTENT(in) :: &
      & ftranspiration_sl(:,:)       !< Fraction of transpiration per soil layer [m/(time step)] (fraction of what???)

    ! Local variables

    INTEGER ::                     &
      & ic,                        & !< Index for grid cells
      & is,                        & !< Index for soil layers
      & ie,                        & !< Index for errors
      & n_errors                     !< Number of errors

    CHARACTER(len=4096) ::         &
      & message_text_long            !< Extra long string for messages

    REAL(wp) ::                    &
      & dummy_wtr_soil(nc, nsoil), & !< Dummy soil moisture storage [m]
      & deficit,                   & !< Water deficit within the actual soil layer [m]
      & remaining,                 & !< Water remaining in the soil layer [m]
      & rootfract,                 & !< Fraction of the soil layer within the root zone []
      & fixed,                     & !< Amount of water below the wilting point, not available for plants [m]
      & root_depth(nc),            & !< Total depth of the root zone [m]
      & transpiration_sum(nc)        !< Total transpiration needed for error check [m]

    INTEGER :: error_idx(nc)         !< Index of grid cell with water balance issue
    LOGICAL :: error_flag(nc)        !< Flag indicating water balance issues

    CHARACTER(len=*), PARAMETER :: routine = modname//':diagnose_evapotrans'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(dummy_wtr_soil, root_depth)

#ifndef _OPENACC
    IF (enforce_water_budget /= WB_IGNORE) THEN ! We do not ignore water balance issues

      ! Check for grid cells with negative soil moisture
      n_errors = 0
      error_flag(:) = .FALSE.
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (wtr_soil_sl(ic,is) < 0._wp .AND. l_fract(ic) .AND. .NOT. error_flag(ic)) THEN
            n_errors = n_errors + 1
            error_idx(n_errors) = ic
            error_flag(ic) = .TRUE.
          END IF
        END DO
      END DO

      ! Write an error message in case negative soil moisture was found.
      DO ie = 1, n_errors
        ic = error_idx(ie)
        WRITE (message_text,*) 'negative soil moisture (routine start) at ', &
          & lat(ic),'N and ',lon(ic),'E:', NEW_LINE('a'), &
          & 'Soil moisture: ', wtr_soil_sl(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish (TRIM(routine), message_text)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message (TRIM(routine), message_text, all_print=.TRUE.)
        END IF
      END DO
    END IF
#endif

    !  Preparations
    ! --------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      root_depth(ic) = 0._wp
    END DO

    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        root_depth(ic) = root_depth(ic) + droot(ic,is)
      END DO
    END DO

    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        ! Use dummy soil moisture storage for diagnostics, because the real soil moisture storage
        !   is only to be changed during the vertical transport (routine soilhyd).
        dummy_wtr_soil(ic,is)   = wtr_soil_sl(ic,is)
        transpiration_sl(ic,is) = 0._wp
      END DO
    END DO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc

      evapo_soil_sl1(ic) = 0._wp     ! Water extracted from upper soil layer by soil evaporation
      deficit_sevap(ic)  = 0._wp     ! Water missing for soil evaporation due to limited soil water
      deficit_trans(ic)  = 0._wp     ! Water missing for transpiration due to limited soil water

      ! 1. Bare soil evaporation (mobile water only)
      ! -------------------------

      IF (evapo_soil_in_m(ic) < 0._wp) THEN                               ! Actual evaporation flux (not dew)
        dummy_wtr_soil(ic,1) = dummy_wtr_soil(ic,1) + evapo_soil_in_m(ic) ! Take water from uppermost soil layer.
        IF (dummy_wtr_soil(ic,1) < 0._wp) THEN                            ! If there is not enough mobile soil water
          deficit_sevap(ic) = dummy_wtr_soil(ic,1)                        !   this is the amount of water missing.
          dummy_wtr_soil(ic,1) = 0._wp
        END IF
        evapo_soil_sl1(ic) = evapo_soil_in_m(ic) - deficit_sevap(ic)      ! Reduce soil evaporation accordingly.
      END IF

      !  2. Transpiration
      ! ------------------

      ! Transpiration should not occur from sealed grid cells (with rooting_depth zero), as no water
      ! is reachable for transpiration. ==> add flux to deficit term.
      ! TODO: what does rooting depth=zero for bare soil tile?
      IF(root_depth(ic) <= 0._wp .AND. transpiration_in_m(ic) < 0._wp) THEN
        deficit_trans(ic) = transpiration_in_m(ic)
      END IF
    END DO
    !$ACC END PARALLEL

    IF (model_scheme == MODEL_QUINCY .AND. PRESENT(ftranspiration_sl)) THEN
      ! With QUINCY, we do not need the redistribution and deficit calculation.
      ! Slight excess transpiration below the wilting point is tolerated, in which case the transpiration function
      ! will then be zero for this layer in the next time step, so that there is no excess transpiration possible
      ! ftranspiration_sl is calculated in q_assimilation/mo_q_assimi_process
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          transpiration_sl(ic,is) = transpiration_in_m(ic) * ftranspiration_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

    ELSE ! jsbach

      ! Get the water needed for transpiration equally from all over the root zone.

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(deficit, rootfract, fixed, remaining)
        DO ic = 1, nc

          IF (root_depth(ic) > 0._wp .AND. transpiration_in_m(ic) < 0._wp) THEN
            ! Grid cell not sealed and actual transpiration (not dew) takes place

            ! Water needed for transpiration from the current layer
            transpiration_sl(ic,is) = transpiration_in_m(ic) * droot(ic,is) / root_depth(ic)

            IF (droot(ic,is) > 0._wp) THEN  ! Roots are reaching into the soil layer.

              rootfract = droot(ic,is) / dsoil(ic,is)    ! Fraction of the soil layer that is within the root zone
              fixed = wtr_field_cap_sl(ic,is) * w_soil_wilt_fract * rootfract     ! Amount of unavailable root zone water

              IF (dummy_wtr_soil(ic,is) * rootfract >= fixed) THEN  ! There is plant available water in the root zone
                ! Remaining water in the root zone after transpiration
                remaining = dummy_wtr_soil(ic,is) * rootfract + transpiration_sl(ic,is)
                IF (remaining < fixed) THEN
                  ! In case transpiration exceeds the amount of available water we diagnose the deficit and
                  ! update soil water: the fixed amount + the potentially available soil water below the root zone.
                  deficit = remaining - fixed
                  dummy_wtr_soil(ic,is) = fixed + dummy_wtr_soil(ic,is) * (1._wp-rootfract)
                ELSE
                  ! If there is enough water in the root zone there is no deficit and we get the soil layer
                  ! water content from the remaining water in the root zone and the soil layer water below
                  ! the root zone.
                  deficit = 0._wp
                  dummy_wtr_soil(ic,is) = remaining + dummy_wtr_soil(ic,is) * (1._wp-rootfract)
                END IF
              ELSE  ! There is no plant available water
                deficit = transpiration_sl(ic,is)
              END IF
              ! The transpiration flux from the soil layer is reduced in case of a water deficit.
              transpiration_sl(ic,is) = transpiration_sl(ic,is) - deficit

              ! Sum up deficits from the different layers
              deficit_trans(ic) = deficit_trans(ic) + deficit

            END IF   ! Roots reaching into the soil layer
          END IF   ! Transpiration actually happening
        END DO
      END DO
      !$ACC END PARALLEL
    END IF ! jsbach

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      evapo_deficit(ic) = evapo_deficit(ic) + deficit_trans(ic)
    END DO

    ! If there is a transpiration deficit diagnosed for some soil layers, but there is still
    ! available water in other layers of the root zone, we reduce the deficit by also transpiring
    ! water from those layers.
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(remaining)
      DO ic = 1, nc
        IF (deficit_trans(ic) < 0._wp .AND. dummy_wtr_soil(ic,is) > wtr_field_cap_sl(ic,is) * w_soil_wilt_fract) THEN
          ! Add available water to the transpiration to reduce the overall deficit.
          remaining          = dummy_wtr_soil(ic,is) - &
            & MAX(dummy_wtr_soil(ic,is) + deficit_trans(ic), wtr_field_cap_sl(ic,is) * w_soil_wilt_fract)
          deficit_trans(ic)      = deficit_trans(ic)      + remaining
          dummy_wtr_soil(ic,is)   = dummy_wtr_soil(ic,is)   - remaining
          transpiration_sl(ic,is) = transpiration_sl(ic,is) - remaining
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

    ! As writing is not possible on GPUs this code only runs on CPUs
#ifndef _OPENACC
    IF (enforce_water_budget /= WB_IGNORE) THEN
      n_errors = 0
      error_flag(:) = .FALSE.
      transpiration_sum(:) = 0._wp

      ! Make sure there is no negative soil moisture diagnosed after applying ET
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (dummy_wtr_soil(ic,is) < 0._wp .AND. l_fract(ic) .AND. .NOT. error_flag(ic)) THEN
            n_errors = n_errors + 1
            error_idx(n_errors) = ic
            error_flag(ic) = .TRUE.
          END IF
          transpiration_sum(ic) = transpiration_sum(ic) + transpiration_sl(ic,is)
        END DO
      END DO

      DO ie = 1, n_errors
        ic= error_idx(ie)
        WRITE (message_text_long,*) 'negative soil moisture (routine end) at ', &
          & lat(ic),'N and ',lon(ic),'E:',                        NEW_LINE('a'), &
          & 'Soil evaporation:         ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Soil evaporation deficit: ', deficit_sevap(ic),      NEW_LINE('a'), &
          & 'Transpiration:            ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Transpiration deficit:    ', deficit_trans(ic),      NEW_LINE('a'), &
          & 'Soil moisture:            ', dummy_wtr_soil(ic,:)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish (TRIM(routine), message_text_long)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message (TRIM(routine), message_text_long, all_print=.TRUE.)
        END IF
      END DO

      ! Check that diagnosed ET and expected ET are still identical
      n_errors = 0
      DO ic = 1, nc
        IF ( ABS((evapo_soil_sl1(ic) + deficit_sevap(ic) + transpiration_sum(ic) + deficit_trans(ic)) &
          &    - (evapo_soil_in_m(ic) + transpiration_in_m(ic)) ) > 1.0e-13_wp &
          &  .AND. l_fract(ic)) THEN
            n_errors = n_errors + 1
            error_idx(n_errors) = ic
        END IF
      END DO

      DO ie = 1, n_errors
        ic = error_idx(ie)
        WRITE (message_text_long,*) 'ET fluxes mismatch at ', &
          & lat(ic),'N and ',lon(ic),'E:',                                  NEW_LINE('a'), &
          & 'Initial mobile soil water:          ', wtr_soil_sl(ic,:),      NEW_LINE('a'), &
          & 'Expected soil evaporation:          ', evapo_soil_in_m(ic),    NEW_LINE('a'), &
          & 'Expected transpiration:             ', transpiration_in_m(ic), NEW_LINE('a'), &
          & 'Diagnosed soil evaporation :        ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Diagnosed transpiration:            ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Diagnosed soil moisture:            ', dummy_wtr_soil(ic,:),   NEW_LINE('a'), &
          & 'Remaining soil evaporation deficit: ', deficit_sevap(ic),      NEW_LINE('a'), &
          & 'Remaining transpiration deficit:    ', deficit_trans(ic)
        IF (enforce_water_budget == WB_ERROR) THEN
          CALL finish (TRIM(routine), message_text_long)
        ELSE IF (enforce_water_budget == WB_LOGGING) THEN
          CALL message (TRIM(routine), message_text_long, all_print=.TRUE.)
        END IF
      END DO
    END IF
#endif

  END SUBROUTINE diagnose_evapotrans

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate snow fraction on lake ice
  !>
  !> We derive the snow fraction on lake ice from snow depth using a hyperbolic tangent function.
  !>
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_wskin_fractions_lice( &
    & weq_snow_lice,                    & ! in
    & fract_snow_lice                   & ! out
    & )

    !$ACC ROUTINE SEQ

    REAL(wp),  INTENT(in)    :: weq_snow_lice    !< Snow depth on lake ice [m water equivalent]
    REAL(wp),  INTENT(out)   :: fract_snow_lice  !< Snow fraction on lake ice []

    fract_snow_lice = TANH(weq_snow_lice * 100._wp)

  END SUBROUTINE calc_wskin_fractions_lice

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate wet and frozen surface fractions on tiles with vegetation
  !>
  !> We calculate the wet or snow covered surface and canopy fractions. Depending on
  !> [[t_hydro_config:pond_dynamics]] different pond shapes are assumed.
  !>
  SUBROUTINE calc_wet_fractions_veg(   &
    & dtime,                           & ! in
    & use_tmx,                         & ! in
    & skinres_max,                     & ! in
    & weq_pond_max,                    & ! in
    & fract_pond_max,                  & ! in
    & pond_dynamics_scheme,            & ! in
    & oro_stddev,                      & ! in
    & t_srf_old,                       & ! in
    & press_srf,                       & ! in
    & heat_tcoef,                      & ! in
    & q_air,                           & ! in
    & wtr_skin,                        & ! in
    & weq_pond,                        & ! in
    & weq_snow_soil,                   & ! in
    & weq_snow_can,                    & ! in
    & fract_snow_can,                  & ! out
    & fract_skin,                      & ! out
    & fract_pond,                      & ! out
    & fract_wet,                       & ! out
    & fract_snow_soil                  & ! out
    & )

    USE mo_phy_schemes,            ONLY: qsat_water, qsat_mixed
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_hydro_constants,        ONLY: Quad_, Tanh_, oro_crit

    REAL(wp), INTENT(in) :: &
      & dtime

    INTEGER, INTENT(in)  :: pond_dynamics_scheme

    LOGICAL, INTENT(in) :: use_tmx

    REAL(wp), INTENT(in), DIMENSION(:) :: &
      & oro_stddev,         & !< Standard deviation of the subgrid scale orography [m]
      & t_srf_old,          & !< Surface temperature prior to update [K]
      & press_srf,          & !< Surface pressure [Pa]
      & heat_tcoef,         & !< Heat transfer coefficient []
      & q_air,              & !< Specific humidity of air []
      & skinres_max,        & !< Maximum capacity of the skin reservoir [m]
      & weq_pond_max,       & !< Maximum capacity of the pond reservoir [m]
      & fract_pond_max,     & !< Maximum pond fraction []
      & wtr_skin,           & !< Amount of water in the skin reservoir [m]
      & weq_pond,           & !< Amount of water/ice in ponds reservoir [m water equivalent]
      & weq_snow_soil,      & !< Amount of snow on the ground [m water equivalent]
      & weq_snow_can          !< Amount of snow on the canopy [m water equivalent]

    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & fract_snow_can,     & !< Snow fraction on the canopy []
      & fract_skin,         & !< Wet skin fraction (not including ponds) []
      & fract_pond,         & !< Pond fraction []
      & fract_wet,          & !< Wet surface fraction (including ponds) []
      & fract_snow_soil       !< Snow fraction on the ground []

    ! Local variables
    REAL(wp) ::             &
      & qsat_srf_old,       & !< Surface saturated humidity from beginning of time step []
      & evapo_pot             !< Potential evaporation

    INTEGER ::              &
      & nc,                 & !< Vector length
      & ic                    !< Grid cell index

    nc = SIZE(press_srf)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) PRIVATE(qsat_srf_old, evapo_pot)
    DO ic=1,nc

      ! Fractional filling of the skin reservoirs for snow on canopy and for liquid water on soil and canopy
      IF (skinres_max(ic) > EPSILON(1._wp)) THEN
        fract_snow_can(ic) = MIN(1._wp, weq_snow_can(ic) / skinres_max(ic))
        fract_skin(ic)     = MIN(1._wp, wtr_skin(ic)     / skinres_max(ic))
      ELSE
        fract_snow_can(ic) = 0._wp
        fract_skin(ic)     = 0._wp
      END IF

      ! Snow fraction on soil
      ! In case the soil is snow free but there is snow on the canopy, we add the canopy snow fraction to the
      ! soil snow fraction (???)
      fract_snow_soil(ic) = Get_snow_fract_noforest(weq_snow_soil(ic), oro_stddev(ic)) ! snow on soil below forest (???)
      fract_snow_soil(ic) = MERGE(fract_snow_can(ic), fract_snow_soil(ic), &
                               fract_snow_soil(ic) < EPSILON(1._wp) .AND. fract_snow_can(ic) > EPSILON(1._wp))

      ! Update pond fraction
      IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        & pond_dynamics_scheme == Quad_) THEN
        fract_pond(ic) = MIN(1._wp, (weq_pond(ic) / weq_pond_max(ic))**0.5_wp) * fract_pond_max(ic)
      ELSE IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        &      pond_dynamics_scheme == Tanh_) THEN
        fract_pond(ic) = MIN(1._wp, (TANH(pi * (weq_pond(ic) / weq_pond_max(ic))))**(oro_stddev(ic)/oro_crit)) &
          &              * fract_pond_max(ic)
      ELSE
        fract_pond(ic) = 0._wp
      END IF

      ! In case potential evaporation within the time step would lead to a loss of snow or water exceeding the
      ! available amounts, we reduce the wet and frozen fractions accordingly. (Water budgets remain unchanged.)

      ! Potential evaporation using old values of air and surface humidity
      IF (use_tmx) THEN
        qsat_srf_old = qsat_water(t_srf_old(ic), press_srf(ic))
      ELSE
        qsat_srf_old = qsat_mixed(t_srf_old(ic), press_srf(ic))
      END IF
      evapo_pot    = -1._wp * heat_tcoef(ic) * (q_air(ic) - qsat_srf_old)  ! Positive upwards

      IF (fract_snow_soil(ic) > 0._wp) THEN
        ! @todo Shouldn't one take rhoice here insteady of rhoh2o?
        fract_snow_soil(ic) = fract_snow_soil(ic) / MAX(1._wp, fract_snow_soil(ic) * evapo_pot * dtime             &
          &                                                    / (rhoh2o * (weq_snow_soil(ic) + weq_snow_can(ic))) &
          &                                            )
      END IF
      IF (fract_skin(ic) > 0._wp ) THEN
        fract_skin(ic) = fract_skin(ic) / MAX(1._wp, (1._wp - fract_snow_soil(ic)) * evapo_pot * dtime      &
          &                                          / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin(ic)))  &
          &                                )
      END IF
      IF (fract_pond(ic) > 0._wp ) THEN
        fract_pond(ic) = fract_pond(ic) / MAX(1._wp, fract_pond(ic) * evapo_pot * dtime                &
        &                                          / (rhoh2o * MAX(EPSILON(1._wp), weq_pond(ic)))    &
        &                                  )
      END IF

      ! Combine the different wet surface fractions assuming that skin and pond locations
      ! are evenly distributed within the tile
      fract_wet(ic) = fract_pond(ic) + fract_skin(ic) * (1._wp - fract_pond(ic))

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_wet_fractions_veg

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate the wet and frozen surface fractions on tiles without vegetation
  !>
  !> We calculate the wet or snow covered surface fractions. Depending on
  !> [[t_hydro_config:pond_dynamics]] different pond shapes are assumed.
  !>
  ! TODO: merge subroutine with calc_wet_fractions_veg to reduce code duplications
  SUBROUTINE calc_wet_fractions_bare( &
    & dtime,                            & ! in
    & use_tmx,                          & ! in
    & skinres_max,                      & ! in
    & weq_pond_max,                     & ! in
    & fract_pond_max,                   & ! in
    & pond_dynamics_scheme,             & ! in
    & oro_stddev,                       & ! in
    & t_srf_old,                        & ! in
    & press_srf,                        & ! in
    & heat_tcoef,                       & ! in
    & q_air,                            & ! in
    & wtr_skin,                         & ! in
    & weq_pond,                         & ! in
    & weq_snow_soil,                    & ! in
    & fract_skin,                       & ! out
    & fract_pond,                       & ! out
    & fract_wet,                        & ! out
    & fract_snow_soil                   & ! out
    & )

    USE mo_phy_schemes,            ONLY: qsat_water, qsat_mixed
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_hydro_constants,        ONLY: Quad_, Tanh_, oro_crit

    REAL(wp), INTENT(in) :: &
      & dtime                 ! Time step length [s]

    INTEGER, INTENT(in) ::  &
      & pond_dynamics_scheme  ! Scheme for pond dynamics with assumptions on the pond shape

    LOGICAL, INTENT(in) ::  &
      & use_tmx               ! Use tmx scheme for turbulent vertical mixing of the atmosphere

    REAL(wp), INTENT(in), DIMENSION(:) :: &
      & oro_stddev,         & ! Standard deviation of the subgrid scale topography [m]
      & t_srf_old,          & ! Surface temperature prior to being updated [K]
      & press_srf,          & ! Surface pressure [Pa]
      & heat_tcoef,         & ! Heat transfere coefficient []
      & q_air,              & ! Specific moisture of the lowest atmosphere level []
      & skinres_max,        & ! Maximum capacity of the skin reservoir [m]
      & weq_pond_max,       & ! Maximum amount of water/ice in ponds [m water equivalent]
      & fract_pond_max,     & ! Maximum surface fraction covered by ponds []
      & wtr_skin,           & ! Amount of water in the skin reservoir (without ponds) [m]
      & weq_pond,           & ! Amount of water/ice in ponds [m water equivalent]
      & weq_snow_soil         ! Amount of snow on the ground [m water equivalent]

    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & fract_skin,         & ! Wet skin fraction (not including ponds) []
      & fract_pond,         & ! Pond fraction []
      & fract_wet,          & ! Wet surface fraction (skin and ponds) []
      & fract_snow_soil       ! Snow fraction []

    REAL(wp) ::             &
      & qsat_srf_old,       & ! Specific humidity at saturation []
      & evapo_pot             ! Potential evaporation [kg m-2 s-1]

    INTEGER ::              &
      & nc,                 & ! Number of grid cells in the vector
      & ic                    ! Index for grid cells

    nc = SIZE(press_srf)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) PRIVATE(qsat_srf_old, evapo_pot)
    DO ic=1,nc

      ! The wet skin fraction is simply calculated from the fractional filling of the skin reservoir.
      IF (skinres_max(ic) > EPSILON(1._wp)) THEN
        fract_skin(ic) = MIN(1._wp, wtr_skin(ic) / skinres_max(ic))
      ELSE
        fract_skin(ic) = 0._wp
      END IF

      ! Snow fraction on soil
      fract_snow_soil(ic) = Get_snow_fract_noforest(weq_snow_soil(ic), oro_stddev(ic))

      ! Update pond fraction
      IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        & pond_dynamics_scheme == Quad_) THEN
        fract_pond(ic) = MIN(1._wp, (weq_pond(ic) / weq_pond_max(ic))**0.5_wp) * fract_pond_max(ic)
      ELSE IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        &      pond_dynamics_scheme == Tanh_) THEN
        fract_pond(ic) = MIN(1._wp, (TANH(pi * (weq_pond(ic) / weq_pond_max(ic))))**(oro_stddev(ic)/oro_crit)) &
          &              * fract_pond_max(ic)
      ELSE
        fract_pond(ic) = 0._wp
      END IF

      ! In case potential evaporation within the time step would lead to a loss of snow or water exceeding the
      ! available amounts, we reduce the wet and frozen fractions accordingly. (Water budgets remain unchanged.)

      ! Potential evaporation using old values of air and surface humidity
      IF (use_tmx) THEN
        qsat_srf_old = qsat_water(t_srf_old(ic), press_srf(ic))
      ELSE
        qsat_srf_old = qsat_mixed(t_srf_old(ic), press_srf(ic))
      END IF
      evapo_pot    = heat_tcoef(ic) * (qsat_srf_old - q_air(ic)) ! Positive upwards

      IF (fract_snow_soil(ic) > 0._wp .AND. weq_snow_soil(ic) > EPSILON(1._wp)) THEN
        fract_snow_soil(ic) = fract_snow_soil(ic) / MAX(1._wp, fract_snow_soil(ic) * evapo_pot * dtime            &
          &                                                    / (rhoh2o * MAX(EPSILON(1._wp), weq_snow_soil(ic))) )
      END IF
      IF (fract_skin(ic) > 0._wp .AND. wtr_skin(ic) > EPSILON(1._wp)) THEN
        fract_skin(ic) = fract_skin(ic) / MAX(1._wp, (1._wp - fract_snow_soil(ic)) * evapo_pot * dtime      &
          &                                                 / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin(ic))) )
      END IF
      IF (fract_pond(ic) > 0._wp ) THEN
        fract_pond(ic) = fract_pond(ic) / MAX(1._wp, fract_pond(ic) * evapo_pot * dtime                &
          &                                          / (rhoh2o * MAX(EPSILON(1._wp), weq_pond(ic)))    &
          &                                  )
      END IF

      ! Combine the different wet surface fractions assuming that skin and pond locations
      ! are evenly distributed within the tile
      fract_wet(ic) = fract_pond(ic) + fract_skin(ic) * (1._wp - fract_pond(ic))

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_wet_fractions_bare

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate the snow fraction
  !>
  !> We calculate the snow fraction following [Roesch et al. 2001](https://link.springer.com/article/10.1007/s003820100153)
  !>
#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION Get_snow_fract_noforest(snow, orodev)

    !$ACC ROUTINE SEQ

  USE mo_hydro_constants, ONLY: wsn2fract_const, wsn2fract_eps, wsn2fract_sigfac

    REAL(wp), INTENT(in) :: &
      & snow,   & !< Amount of snow on the ground [m water equivalent]
      & orodev    !< Standard deviation of the subgrid scale orography [m]

    Get_snow_fract_noforest = wsn2fract_const * TANH(snow * 100._wp) &
                               & * SQRT(snow * 1000._wp / (snow * 1000._wp + wsn2fract_eps + wsn2fract_sigfac * orodev))

  END FUNCTION Get_snow_fract_noforest

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate unstressed canopy conductance
  !>
  !> We here calculate canopy conductance as if there was no water stress.
  !>
#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION get_canopy_cond_unstressed_simple(lai, par) RESULT(conductance)

    USE mo_hydro_constants, ONLY: k => conductance_k, a => conductance_a, b => conductance_b, c => conductance_c

    REAL(wp), INTENT(in) :: &
      & lai,   & !< Leaf area index []
      & par      !< Photosynthetically active radiation [w/m2]

    ! Local variables
    !
    REAL(wp) :: &
      & d, zpar  !< Helper variables

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_canopy_cond_unstressed_simple'

    !$ACC ROUTINE SEQ

    zpar = MAX(1.E-10_wp, par)

    IF (lai > EPSILON(1._wp)) THEN
      d = (a + b*c) / (c * zpar)
      conductance = ( LOG((d * EXP(k*lai) + 1._wp) / (d + 1._wp)) * b / (d * zpar) - &
                    & LOG((d + EXP(-k*lai)) / (d + 1._wp))                           &
                    ) / (k * c)
    ELSE
      conductance = EPSILON(1._wp)
    END IF

  END FUNCTION get_canopy_cond_unstressed_simple

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate canopy conductance under water stress
  !>
  !> Calculation of the stomatal conductance accounting for water stress.
  !>
#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION get_canopy_cond_stressed_simple(cond_unstressed, water_stress, air_is_saturated) RESULT(conductance)

    REAL(wp), INTENT(in) :: &
      & cond_unstressed,    & !< Canopy conductance without water stress [m/s]
      & water_stress          !< Water stress factor []
    LOGICAL,  INTENT(in) :: &
      & air_is_saturated      !< True: air in lowest atmosphere level is water saturated

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_canopy_cond_stressed_simple'

    !$ACC ROUTINE SEQ

    IF (air_is_saturated) THEN
      conductance = EPSILON(1._wp)
    ELSE
      conductance = cond_unstressed * water_stress
    END IF

  END FUNCTION get_canopy_cond_stressed_simple

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate the water stress factor
  !>
  !> We calculate the water stress factor from the actual amount of available water in the root
  !> zone and the maximum possible amount.
  !>
#ifndef _OPENACC
  ELEMENTAL &
#endif
  FUNCTION get_water_stress_factor ( &
    & w_soil, w_soil_max, w_soil_crit_fract, w_soil_wilt_fract) RESULT(water_stress_factor)

    !$ACC ROUTINE SEQ

    ! TODO
    ! TBD: Maybe change later so that it uses geographically dependent critical values (like wilting point)

    REAL(wp), INTENT(in) ::  &
      & w_soil,              & !< Plant available water in the root zone [m]
      & w_soil_max,          & !< Maximum possible available water in root zone [m]
      & w_soil_crit_fract,   & !< Relative root zone available water content at critical point,
                               !< i.e. below which plants sturt to suffer from water stress []
      & w_soil_wilt_fract      !< Relative root zone available water content at wilting point []
    REAL(wp) ::              &
      & water_stress_factor    !< RESULT: water stress factor (0: extreme stress, 1: no stress) []

    ! Local variables
    REAL(wp) ::              &
      & w_crit,              & !< Root zone available water content at critical point [m]
      & w_wilt                 !< Root zone available water content at wilting point [m]

    w_crit = w_soil_max * w_soil_crit_fract
    w_wilt = w_soil_max * w_soil_wilt_fract

    IF (w_crit - w_wilt > 0._wp) THEN
      water_stress_factor = MAX (0._wp, &
        &                        MIN (1._wp, (w_soil - w_wilt) / (w_crit - w_wilt) ))
    ELSE
      water_stress_factor = 0._wp
    END IF

  END FUNCTION get_water_stress_factor

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate vertical transport
  !>
  !> This is the core routine of vertical transport of liquid water through the soil layers.
  !>
  ! TODO: Add link documentation by Christian Reick
  !
  SUBROUTINE calc_vertical_transport(dt, alpha, n_act, dzf, top_bound, bot_bound, S, C, P, Q, X)

    USE mo_util, ONLY: tdma_solver_vec ! OpenACC-enabled solver for tridiagonal matrix (Thomas algorithm)

    INTEGER, INTENT(in) ::  &
      & n_act(:)              !< Number of layers
    REAL(wp), INTENT(in) :: &
      & dt,                 & !< Time step
      & alpha,              & !< Alpha
      & top_bound(:),       & !< Upper boundary fluxes
      & bot_bound(:),       & !< Lower boundary fluxes
      & dzf(:,:),           & !< Layer thickness
      & C(:,:),             & !< Heat capacity if routine is used for heat transport,
                              !< set to 1 for water transport
      & S(:,:),             & !< Prescribed water/heat change flux per layer
      & Q(:,:),             & !< Diffusivity at (upper) layer boundaries
      & P(:,:)                !< Conductivity at (upper) layer boundaries

    REAL(wp), INTENT(inout) :: &
      & X(:,:)                !< Normalized soil state for every layer

    ! Local variables
    INTEGER :: &
      & j,     & !< Looping index
      & ic,    & !< Looping index for grid cells
      & nc,    & !< Vector length
      & nsoil    !< Number of soil layers (vertical dimension)
    REAL(wp) :: &
      & dzh(SIZE(S,1),SIZE(S,2)),     & !< Distance between mid layer depth of current and above layer
      & matrix(SIZE(S,1),SIZE(S,2),4)   !< Matrix containing sub-diagonal (1), diagonal (2),
                                        !< super-diagonal (3) of triangular matrix and rhs (4) of
                                        !< linear equation system

    nc    = SIZE(S,1)
    nsoil = SIZE(S,2)

    !$ACC DATA CREATE(dzh, matrix) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO j=1,nsoil
      DO ic=1,nc
        dzh(ic,j) = 0.0
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO j=2,nsoil
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        IF (j <= n_act(ic)) THEN
          dzh(ic,j) = (dzf(ic,j-1) + dzf(ic,j)) / 2._wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    ! Note that P & Q refer to the inter face above a given lvl not below.
    ! Hence for j == 1, there is no P(j) and Q(j)
    ! Make sure that this is consistent with the definition of input params.

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO j = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        IF(j == 1) THEN
          ! Top soil layer
          matrix(ic,j,1) =  0.0_wp
          matrix(ic,j,2) =  P(ic,j+1) / dzh(ic,j+1) + dzf(ic,j) * C(ic,j) / (alpha * dt)
          matrix(ic,j,3) = -P(ic,j+1) / dzh(ic,j+1)
          matrix(ic,j,4) = -Q(ic,j+1) + dzf(ic,j) * S(ic,j) + dzf(ic,j) * C(ic,j) * X(ic,j) / (alpha * dt) + top_bound(ic)
        ELSE IF (j == n_act(ic)) THEN
          ! Bottom soil layers
          matrix(ic,j,1) = -P(ic,j) / dzh(ic,j)
          matrix(ic,j,2) =  P(ic,j) / dzh(ic,j) + dzf(ic,j) * C(ic,j) / (alpha * dt)
          matrix(ic,j,3) =  0.0_wp
          matrix(ic,j,4) =  Q(ic,j)          + dzf(ic,j) * S(ic,j) + dzf(ic,j) * C(ic,j) * X(ic,j) / (alpha * dt) - bot_bound(ic)
        ELSE IF (j < n_act(ic)) THEN
          ! Middle soil layers
          matrix(ic,j,1) = -P(ic,j) / dzh(ic,j)
          matrix(ic,j,2) =  P(ic,j) / dzh(ic,j) +  P(ic,j+1) / dzh(ic,j+1) + dzf(ic,j) * C(ic,j) / (alpha * dt)
          matrix(ic,j,3) = -P(ic,j+1) / dzh(ic,j+1)
          matrix(ic,j,4) =  Q(ic,j) - Q(ic,j+1) + dzf(ic,j) * S(ic,j) + dzf(ic,j) * C(ic,j) * X(ic,j) / (alpha * dt)
        ELSE
          matrix(ic,j,1) = 0._wp
          matrix(ic,j,2) = 1._wp
          matrix(ic,j,3) = 0._wp
          matrix(ic,j,4) = 1._wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    DO ic=1,nc
      IF( matrix(ic,1,2) .EQ. 0.0) THEN
        WRITE (message_text,*) 'something went terribly wrong ... rewrite equations!', C(ic,1), dzf(ic,1)
        CALL finish('calc_vertical_transport', message_text)
      END IF
    END DO
#endif

    CALL tdma_solver_vec(matrix(:,:,1), matrix(:,:,2), matrix(:,:,3), matrix(:,:,4), &
      &                  1, nsoil, 1, nc, X(:,:), opt_acc_queue=acc_stream)

  !$ACC END DATA
  END SUBROUTINE calc_vertical_transport

  !----------------------------------------------------------------------------------------------
  !>
  !> #### Calculate soil hydraulic properties
  !>
  !> In this routine we update hydraulic conductivity and diffusivity on soil layers and at layer
  !> interface. The properties depend on ice and organic matter content of the soils and thus
  !> need to be updated each time step. Different parametrization schemes are supported (compare
  !> [[t_hydro_config:soilhydmodel]]).
  !>
  SUBROUTINE get_soilhyd_properties(                &
            & soilhydmodel, interpol_mean,          &
            & nc, nsoil, dsoil,                     &
            & wtr, ice, wscools, wsat, wres,        &
            & k_sat, mpot_sat, bclapp, ps_index,    &
            & dt, last_soil_layer,                  &
            & ice_impedance, K, D, K_inter,         &
            & D_inter, mpot_act                     &
            )

    USE mo_hydro_constants,    ONLY: BrooksCorey_, Campbell_, VanGenuchten_, &
      &                              Upstream_, Arithmetic_, &
      &                              matric_pot_min

    INTEGER,  INTENT(in)            :: soilhydmodel       !< Model to determine conductivites (K) & diffusivities (D)
    INTEGER,  INTENT(in)            :: interpol_mean      !< Method to derive values at layer interface (avg scheme)
    INTEGER,  INTENT(in)            :: nc                 !< Vector length
    INTEGER,  INTENT(in)            :: nsoil              !< Number of below ground layers
    REAL(wp), INTENT(in)            :: dsoil(:,:)         !< Soil depth (until bedrock) per layer [m]
    REAL(wp), INTENT(in)            :: wtr(:,:)           !< Soil water content [m]
    REAL(wp), INTENT(in)            :: ice(:,:)           !< Soil ice content [m]
    REAL(wp), INTENT(in)            :: wscools(:,:)       !< Amount of potentially supercooled water [m]
    REAL(wp), INTENT(in)            :: wsat(:,:)          !< Soil depth (until bedrock) per layer [m]
    REAL(wp), INTENT(in)            :: wres(:,:)          !< Residual water content [m]
    REAL(wp), INTENT(in)            :: k_sat(:,:)         !< Hydraulic conductivity at saturation [m/s]
    REAL(wp), INTENT(in)            :: mpot_sat(:,:)      !< Matric potential for saturated soils [m]
    REAL(wp), INTENT(in)            :: bclapp(:,:)        !< Clapp & Hornberger exponent b []
    REAL(wp), INTENT(in)            :: ps_index(:,:)      !< Pore size index []
    REAL(wp), INTENT(in),  OPTIONAL :: dt                 !< Time step length [s]
    INTEGER,  INTENT(in),  OPTIONAL :: last_soil_layer(:) !< Index of deepest soil layers (above bedrock)
    REAL(wp), INTENT(out), OPTIONAL :: ice_impedance(:,:) !< Impedance factor to account for ice blocking flowpaths []
    REAL(wp), INTENT(out), OPTIONAL :: K(:,:)             !< Hydraulic conductivity of the soil layer [m/s]
    REAL(wp), INTENT(out), OPTIONAL :: D(:,:)             !< Diffusivity of the soil layer [m2/s]
    REAL(wp), INTENT(out), OPTIONAL :: K_inter(:,:)       !< Hydraulic conductivity at (upper) layer interface [m/s]
    REAL(wp), INTENT(out), OPTIONAL :: D_inter(:,:)       !< Diffusivity at (upper) layer interface [m2/s]
    REAL(wp), INTENT(out), OPTIONAL :: mpot_act(:,:)      !< Matric potential at actual soil state [m]

    ! Local variables
    INTEGER  ::                  &
      & ic,                      & !< Grid cell index
      & is                         !< Soil layer index
    REAL(wp) ::                  &
      & ck,                      & !< Parameter for matric potential formula (Zhang 2007)
      & local_dt                   !< Time step, copy needed for numerical reasons
    REAL(wp) ::                  &
      & wtr_vol(nc,nsoil),       & !< Volumetric soil water content []
      & ice_vol(nc,nsoil),       & !< Volumetric soil ice content []
      & scool_vol(nc,nsoil),     & !< Volumetric amount of supercooled water []
      & nvgn(nc,nsoil),          & !< Parameter with van Genuchten scheme
      & mvgn(nc,nsoil),          & !< Parameter with van Genuchten scheme
      & ws_rel_vol(nc,nsoil),    & !< Relative volumetric soil moisture []
      & ws_range_vol(nc,nsoil),  & !< Reference maximum soil moisture [m]
      & wsat_vol(nc,nsoil),      & !< Volumetric maximum soil moisture []
      & wres_vol(nc,nsoil),      & !< Volumetric minimum soil moisture []
      & ws_free_vol(nc,nsoil),   & !< Volumetric mobile soil moisture []
      & ice_imp(nc,nsoil)          !< Ice impedance factor []

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(wtr_vol, ice_vol, scool_vol, nvgn, mvgn, ws_rel_vol, ws_range_vol) &
    !$ACC   CREATE(wsat_vol, wres_vol, ws_free_vol, ice_imp)

    ck = 8.0_wp ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc

        ! Set defaults
        ws_rel_vol(ic,is) = 0._wp
        ice_imp(ic,is)    = 1._wp
        nvgn(ic,is)       = ps_index(ic,is) + 1._wp
        mvgn(ic,is)       = ps_index(ic,is) / nvgn(ic,is)

        ! Convert to volumetric quantities
        IF (dsoil(ic,is) > 0._wp) THEN
          wtr_vol(ic,is)   =  wtr(ic,is)                      / dsoil(ic,is)
          ice_vol(ic,is)   =  ice(ic,is)                      / dsoil(ic,is)
          scool_vol(ic,is) =  MIN(wscools(ic,is), wtr(ic,is)) / dsoil(ic,is)
          wsat_vol(ic,is)  =  wsat(ic,is)                     / dsoil(ic,is)
          wres_vol(ic,is)  =  wres(ic,is)                     / dsoil(ic,is)
        ELSE
          wtr_vol(ic,is)   =  0._wp
          ice_vol(ic,is)   =  0._wp
          scool_vol(ic,is) =  0._wp
          wsat_vol(ic,is)  =  0._wp
          wres_vol(ic,is)  =  0._wp
        END IF
      END DO
    END DO
    !$ACC END LOOP

    IF (PRESENT(K)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          K(ic,is)       = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(D)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          D(ic,is)       = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(K_inter)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          K_inter(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(D_inter)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          D_inter(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(mpot_act)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          mpot_act(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    !$ACC END PARALLEL

    ! Compute impedance factor due to soil ice (e.g. Hansson et al., 2004)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        IF (wsat_vol(ic,is) > 1.0e-10_wp) THEN
          ice_imp(ic,is)  = 10._wp**(-6._wp * MIN(1._wp, (ice_vol(ic,is) + scool_vol(ic,is)) / wsat_vol(ic,is)))
        ELSE
          ice_imp(ic,is)  = 1._wp
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    IF (PRESENT(ice_impedance)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          ice_impedance(ic,is) = ice_imp(ic,is)
        END DO
      END DO
      !$ACC END LOOP
    END IF

    ! Compute relative soil moisture depending on soil hydrological model - default is "Van Genuchten"
    IF (soilhydmodel == Campbell_) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_free_vol(ic,is)   = MAX(0._wp, wtr_vol(ic,is) - scool_vol(ic,is))
          ws_range_vol(ic,is)  = wsat_vol(ic,is)
        END DO
      END DO
      !$ACC END LOOP
    ELSE  ! Brooks & Corey; Van Genuchten
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_free_vol(ic,is)   = MAX(0._wp, wtr_vol(ic,is) - MAX(scool_vol(ic,is), wres_vol(ic,is)))
          ws_range_vol(ic,is)  = wsat_vol(ic,is) - wres_vol(ic,is)
        END DO
      END DO
      !$ACC END LOOP
    END IF
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (ws_range_vol(ic,is) > 0._wp) THEN
          ws_rel_vol(ic,is)    = ws_free_vol(ic,is) / ws_range_vol(ic,is)
        END IF
      END DO
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    ! Making sure ws_rel is in the physical range
    ! This is especially important when infiltration is added to the first layer
    IF (soilhydmodel == VanGenuchten_) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_rel_vol(ic,is) = MIN(0.9999_wp, MAX(0._wp, ws_rel_vol(ic,is))) ! doesn't work for wrel == 1;
                                                                            ! (1-ws_rel(ic,is)**(1/mvgn(ic,is)))**(-mvgn)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ELSE  ! Brooks & Corey; Van Genuchten
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_rel_vol(ic,is) = MIN(1._wp, MAX(0._wp, ws_rel_vol(ic,is)))
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Compute hydrological conductivity (K) for every layer
    IF (PRESENT(K)) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)

      IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == Campbell_) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO ic = 1, nc
            K(ic,is) =  k_sat(ic,is)                                                   &
              & * ws_rel_vol(ic,is)**(2._wp * bclapp(ic,is) + 3._wp)
          END DO
        END DO
        !$ACC END LOOP

      ELSE IF (soilhydmodel == VanGenuchten_) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO ic = 1, nc
            IF (mvgn(ic,is) > 0._wp) THEN
              K(ic,is) = k_sat(ic,is) * ws_rel_vol(ic,is)**(0.5_wp) &
                & *(1._wp-(1._wp-ws_rel_vol(ic,is)**(1._wp/mvgn(ic,is)))**mvgn(ic,is))**2._wp
            ELSE
              K(ic,is) = 0._wp
            END IF
          END DO
        END DO
        !$ACC END LOOP
      END IF

      ! Reduce transport velocity in presence of ice
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          K(ic,is) = K(ic,is) * ice_imp(ic,is)
        END DO
      END DO
      !$ACC END LOOP

      !$ACC END PARALLEL

      ! Make sure percolation from top layer is limited to water content !!!
      IF (PRESENT(dt)) THEN
        local_dt = dt ! Optional scalar variables are not treated correctly by NVHPC OpenACC
        IF (local_dt > 0._wp) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1, nc
            K(ic,1) = MIN(K(ic,1), MAX(0._wp,(wtr(ic,1)-wres(ic,1)) / local_dt))
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      END IF
    END IF

    ! Compute hydrological diffusivity (D) for every layer
    IF (PRESENT(D)) THEN

      IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == Campbell_) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is = 1, nsoil
          DO ic = 1, nc
            IF (ws_range_vol(ic,is) > 0._wp) THEN
              D(ic,is) = -k_sat(ic,is) * mpot_sat(ic,is) * bclapp(ic,is) / ws_range_vol(ic,is) &
                & * ws_rel_vol(ic,is)**(1._wp * bclapp(ic,is) + 2._wp)
            ELSE
              D(ic,is) = 0._wp
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (soilhydmodel == VanGenuchten_) THEN
#ifndef __OPENACC__
        IF (.NOT. PRESENT(K)) THEN
          WRITE (message_text,*) 'VanGenuchten diffusivity computation requires hydrological conductivity'
          CALL finish ('get_soilhyd_properties', message_text)
        END IF
#endif
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is = 1, nsoil
          DO ic = 1, nc
            IF (ws_range_vol(ic,is) > 0._wp .AND. ws_rel_vol(ic,is) > 0._wp) THEN
              D(ic,is) = -K(ic,is)* mpot_sat(ic,is) * bclapp(ic,is) / ws_range_vol(ic,is) &
                & * ws_rel_vol(ic,is)**(-1._wp/mvgn(ic,is)) * (1._wp-ws_rel_vol(ic,is)**(1._wp/mvgn(ic,is)))**(-mvgn(ic,is))
            ELSE
              D(ic,is) = 0._wp
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Reduce transport velocity in presence of ice
      ! However, ice seems to attract water --> is this the right approach for diffusivity?
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          D(ic,is) = D(ic,is) * ice_imp(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Determining K & D on layer interface
    ! Possibility
    ! 1) Upstream mean [K_lev; Dmax(x_lev, x_lev-1)}
    ! 2) Adjusted arithmetic mean A (of k, D)
    ! Note that in the vertical transport scheme the notation refers to the upper interface of a level.
    ! Hence the loop goes from last_soil_layer to 2

    IF (PRESENT(K_inter)) THEN

#ifndef __OPENACC__
      IF (.NOT. PRESENT(K) .OR. .NOT. PRESENT(last_soil_layer)) THEN
        WRITE (message_text,*) 'K at layer interface cannot be computed without '//&
          & 'information about actual soil depth and computing K on layers'
        CALL finish ('get_soilhyd_properties', message_text)
      END IF
#endif

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = nsoil, 2, -1
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          IF (is <= last_soil_layer(ic)) THEN
            IF (interpol_mean == Upstream_) THEN
              K_inter(ic,is)  =   K(ic,is-1)
            ELSE IF (interpol_mean == Arithmetic_) THEN
              IF (K(ic,is-1) > K(ic,is)) THEN ! MOVEMENT OF WETTING FRONTS LIMITED BY SATURATION OF SUBJACENT LAYER. ELSE ...
                K_inter(ic,is)  =   (K(ic,is)*dsoil(ic,is)+K(ic,is-1)*dsoil(ic,is-1))/(dsoil(ic,is)+dsoil(ic,is-1))
              ELSE ! ... PREVENT FLUX FROM COMPLETLY DRYING OUT UPPER LAYER
                K_inter(ic,is)  = K(ic,is-1)
              END IF
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END IF

    IF (PRESENT(D_inter)) THEN

#ifndef __OPENACC__
      IF (.NOT. PRESENT(D) .OR. .NOT. PRESENT(last_soil_layer)) THEN
        WRITE (message_text,*) 'D at layer interface cannot be computed without '//&
          & 'information about actual soil depth and computing D on layers'
        CALL finish ('get_soilhyd_properties', message_text)
      END IF
#endif
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = nsoil, 2, -1
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          IF (is <= last_soil_layer(ic)) THEN
            IF (interpol_mean == Upstream_) THEN
              D_inter(ic,is)  = MAX(D(ic,is),D(ic,is-1))
            ELSE IF (interpol_mean == Arithmetic_) THEN
              D_inter(ic,is)  = (D(ic,is)*dsoil(ic,is)+D(ic,is-1)*dsoil(ic,is-1))/(dsoil(ic,is)+dsoil(ic,is-1))
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END IF

    ! Determine saturated soil matric potential for actual soil moisture state
    ! (i.e., matric potential for saturated soils)
    ! Formula from Stuurop 2021 (https://doi.org/10.1016/j.coldregions.2021.103456),
    ! adjusted from Zhang 2007 (https://doi.org/10.1175/JHM605.1)
    IF (PRESENT(mpot_act)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (ws_rel_vol(ic,is) > 0._wp) THEN
            IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == Campbell_) THEN
              mpot_act(ic,is) = mpot_sat(ic,is) * (ws_rel_vol(ic,is)**(-bclapp(ic,is)))
            ELSE
              mpot_act(ic,is) = mpot_sat(ic,is) * (ws_rel_vol(ic,is)**(-1._wp/mvgn(ic,is)) - 1._wp)**(1._wp/nvgn(ic,is))
            END IF
            mpot_act(ic,is) = MAX(matric_pot_min, mpot_act(ic,is) * ((1._wp + ck * ice_vol(ic,is))**2._wp))
          ELSE
            mpot_act(ic,is) = matric_pot_min
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC END DATA

  END SUBROUTINE get_soilhyd_properties

#endif
END MODULE mo_hydro_process
