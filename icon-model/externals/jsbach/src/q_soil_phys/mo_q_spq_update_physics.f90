!> QUINCY soil-physics calculation
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### calculate the task update_spq_physics, i.e., snow accumulation, soil hydrology, surface energy balance
!>
MODULE mo_q_spq_update_physics
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, message_text

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_spq_physics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_spq_update_physics'

CONTAINS

  ! ======================================================================================================= !
  !>update soil moisture and theta for the plant only model - NOTE update/improve docu
  !>
  ! two cases considered:
  !   a) dynamic water uptake calculation
  !   b) prescribed values from namelist
  SUBROUTINE update_spq_physics(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_control,           ONLY: jsbach_runs_standalone
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, Q_ASSIMI_, RAD_, SPQ_, VEG_, TURB_, SEB_, SSE_, HYDRO_
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_grid, Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SPQ_)
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(TURB_)
    dsl4jsb_Use_memory(SEB_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model          !< the model
    TYPE(t_jsb_grid),       POINTER       :: grid
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_w   !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_snow_spq !< Vertical grid snow
    INTEGER                               :: nsoil_w        !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: nsnow          !< number of snow layers
    INTEGER                               :: isoil
    INTEGER                               :: iblk, ics, ice, nc
    INTEGER                               :: ic
    REAL(wp)                              :: dtime          !< timestep length
    REAL(wp)                              :: alpha          !< implicitness factor
    REAL(wp), ALLOCATABLE                 :: drag_srf_local(:)    !< surface drag, provided by atmosphere or TURB_
    REAL(wp), ALLOCATABLE                 :: t_acoef_local(:)     !< variable provided by atmosphere or TURB_
    REAL(wp), ALLOCATABLE                 :: t_bcoef_local(:)     !< variable provided by atmosphere or TURB_
    REAL(wp), ALLOCATABLE                 :: q_acoef_local(:)     !< variable provided by atmosphere or TURB_
    REAL(wp), ALLOCATABLE                 :: q_bcoef_local(:)     !< variable provided by atmosphere or TURB_
    LOGICAL                               :: jsb_standalone       !< model runs standalone?
    REAL(wp), POINTER                     :: lat(:), lon(:)       !< helpers for debugging
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_spq_physics'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SPQ_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk                :: t_air
    dsl4jsb_Real2D_onChunk                :: q_air
    dsl4jsb_Real2D_onChunk                :: press_srf
    dsl4jsb_Real2D_onChunk                :: wind_10m
    dsl4jsb_Real2D_onChunk                :: rain
    dsl4jsb_Real2D_onChunk                :: snow
    dsl4jsb_Real2D_onChunk                :: t_acoef
    dsl4jsb_Real2D_onChunk                :: t_bcoef
    dsl4jsb_Real2D_onChunk                :: q_acoef
    dsl4jsb_Real2D_onChunk                :: q_bcoef
    dsl4jsb_Real2D_onChunk                :: drag_srf
    ! TURB_ 2D
    dsl4jsb_Real2D_onChunk                :: fact_q_air
    dsl4jsb_Real2D_onChunk                :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk                :: fact_qsat_trans_srf
    dsl4jsb_Real2D_onChunk                :: t_acoef_stdalone
    dsl4jsb_Real2D_onChunk                :: t_bcoef_stdalone
    dsl4jsb_Real2D_onChunk                :: q_acoef_stdalone
    dsl4jsb_Real2D_onChunk                :: q_bcoef_stdalone
    dsl4jsb_Real2D_onChunk                :: drag_srf_stdalone
    ! SEB_ 2D
    dsl4jsb_Real2D_onChunk                :: t
    dsl4jsb_Real2D_onChunk                :: t_old
    dsl4jsb_Real2D_onChunk                :: qsat_star
    dsl4jsb_Real2D_onChunk                :: q_sat_srf
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk                :: s_star
    dsl4jsb_Real2D_onChunk                :: evapotranspiration
    dsl4jsb_Real2D_onChunk                :: temp_srf_eff_4
    dsl4jsb_Real2D_onChunk                :: zril_old
    dsl4jsb_Real2D_onChunk                :: interception
    dsl4jsb_Real2D_onChunk                :: evapopot
    dsl4jsb_Real2D_onChunk                :: evaporation
    dsl4jsb_Real2D_onChunk                :: evaporation_snow
    dsl4jsb_Real2D_onChunk                :: srf_runoff
    dsl4jsb_Real2D_onChunk                :: drainage
    dsl4jsb_Real2D_onChunk                :: gw_runoff
    dsl4jsb_Real2D_onChunk                :: drainage_fraction
    dsl4jsb_Real2D_onChunk                :: ground_heat_flx_old
    dsl4jsb_Real2D_onChunk                :: ground_heat_flx
    dsl4jsb_Real2D_onChunk                :: latent_heat_flx
    dsl4jsb_Real2D_onChunk                :: sensible_heat_flx
    dsl4jsb_Real2D_onChunk                :: snow_height
    dsl4jsb_Real2D_onChunk                :: snow_soil_heat_flux
    dsl4jsb_Real2D_onChunk                :: snow_melt_to_soil
    dsl4jsb_Real2D_onChunk                :: snow_srf_heat_flux
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk                :: heat_capa_sl
    dsl4jsb_Real3D_onChunk                :: therm_cond_sl
    dsl4jsb_Real3D_onChunk                :: saxtonA
    dsl4jsb_Real3D_onChunk                :: saxtonB
    dsl4jsb_Real3D_onChunk                :: saxtonC
    dsl4jsb_Real3D_onChunk                :: kdiff_sat_sl
    dsl4jsb_Real3D_onChunk                :: gw_runoff_sl
    dsl4jsb_Real3D_onChunk                :: drainage_sl
    dsl4jsb_Real3D_onChunk                :: w_soil_freeze_flux
    dsl4jsb_Real3D_onChunk                :: w_soil_melt_flux
    dsl4jsb_Real3D_onChunk                :: snow_present_snl
    dsl4jsb_Real3D_onChunk                :: w_snow_snl
    dsl4jsb_Real3D_onChunk                :: t_snow_snl
    dsl4jsb_Real3D_onChunk                :: snow_lay_thickness_snl
    dsl4jsb_Real3D_onChunk                :: w_liquid_snl
    ! SSE_ 3D
    dsl4jsb_Real3D_onChunk                :: t_soil_sl
    dsl4jsb_Real3D_onChunk                :: bulk_dens_sl
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk                :: wtr_soil_root_pot
    ! Q_ASSIMI_ 3D
    dsl4jsb_Real3D_onChunk                :: ftranspiration_sl
    ! RAD_
    dsl4jsb_Real2D_onChunk                :: rad_srf_net
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk                :: height
    dsl4jsb_Real2D_onChunk                :: lai
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk                :: root_fraction_sl
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onChunk                :: fract_wet
    dsl4jsb_Real2D_onChunk                :: fract_snow
    dsl4jsb_Real2D_onChunk                :: fract_snow_soil
    dsl4jsb_Real2D_onChunk                :: fract_snow_can
    dsl4jsb_Real2D_onChunk                :: num_sl_above_bedrock
    dsl4jsb_Real2D_onChunk                :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk                :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk                :: wtr_skin
    dsl4jsb_Real2D_onChunk                :: wtr_rootzone
    dsl4jsb_Real2D_onChunk                :: wtr_rootzone_rel
    dsl4jsb_Real2D_onChunk                :: wtr_plant_avail_rel
    dsl4jsb_Real2D_onChunk                :: transpiration
    dsl4jsb_Real2D_onChunk                :: root_depth
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onChunk                :: soil_depth_sl
    dsl4jsb_Real3D_onChunk                :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onChunk                :: soil_lay_width_sl
    dsl4jsb_Real3D_onChunk                :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onChunk                :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk                :: ice_soil_sl
    dsl4jsb_Real3D_onChunk                :: wtr_soil_pwp_sl
    dsl4jsb_Real3D_onChunk                :: wtr_soil_fc_sl
    dsl4jsb_Real3D_onChunk                :: wtr_soil_sat_sl
    dsl4jsb_Real3D_onChunk                :: wtr_soil_pot_sl
    dsl4jsb_Real3D_onChunk                :: wtr_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    alpha   = options%alpha
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SPQ_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model          => Get_model(tile%owner_model_id)
    grid           => Get_grid(model%grid_id)
    vgrid_soil_w   => Get_vgrid('soil_depth_water')
    vgrid_snow_spq => Get_vgrid('snow_layer_spq')
    nsoil_w        =  vgrid_soil_w%n_levels
    nsnow          =  vgrid_snow_spq%n_levels
    lat            => grid%lat(ics:ice, iblk)
    lon            => grid%lon(ics:ice, iblk)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SPQ_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2l_, t_acoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2l_, t_bcoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2l_, q_acoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2l_, q_bcoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2l_, drag_srf)             ! in
    ! ---------------------------
    ! TURB_ 2D
    dsl4jsb_Get_var2D_onChunk(TURB_, fact_q_air)          ! inout
    dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_srf)       ! inout
    dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_trans_srf) ! in
    dsl4jsb_Get_var2D_onChunk(TURB_, t_acoef_stdalone)    ! in
    dsl4jsb_Get_var2D_onChunk(TURB_, t_bcoef_stdalone)    ! in
    dsl4jsb_Get_var2D_onChunk(TURB_, q_acoef_stdalone)    ! in
    dsl4jsb_Get_var2D_onChunk(TURB_, q_bcoef_stdalone)    ! in
    dsl4jsb_Get_var2D_onChunk(TURB_, drag_srf_stdalone)   ! in
    ! ---------------------------
    ! SEB_ 2D
    dsl4jsb_Get_var2D_onChunk(SEB_, t)                    ! inout
    dsl4jsb_Get_var2D_onChunk(SEB_, t_old)                ! inout
    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_star)            ! inout
    dsl4jsb_Get_var2D_onChunk(SEB_, q_sat_srf)            ! out
    ! ---------------------------
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, s_star)               ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, evapotranspiration)   ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, temp_srf_eff_4)       ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, zril_old)             ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, evapopot)             ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, interception)         ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, evaporation)          ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, evaporation_snow)     ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, srf_runoff)           ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, drainage)             ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, gw_runoff)            ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, drainage_fraction)    ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, ground_heat_flx_old)  ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, ground_heat_flx)      ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, latent_heat_flx)      ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, sensible_heat_flx)    ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, snow_height)          ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, snow_soil_heat_flux)  ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, snow_melt_to_soil)    ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, snow_srf_heat_flux)   ! inout
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_, heat_capa_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, therm_cond_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, saxtonA)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, saxtonB)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, saxtonC)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, kdiff_sat_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, gw_runoff_sl)             ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, drainage_sl)              ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_freeze_flux)       ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_melt_flux)         ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_snow_snl)               ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, t_snow_snl)               ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, snow_lay_thickness_snl)   ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, snow_present_snl)         ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_liquid_snl)             ! inout
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, bulk_dens_sl)             ! in
    ! ---------------------------
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, wtr_soil_root_pot)   ! out
    ! Q_ASSIMI_ 3D
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, ftranspiration_sl)   ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)            ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                    ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)       ! in
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_wet)              ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)             ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)        ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_can)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_soil_root_pwp)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_soil_root_fc)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_skin)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_rootzone_rel)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_plant_avail_rel)    ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_rootzone)           ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, transpiration)          ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, root_depth)             ! in
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)                  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_lay_depth_ubound_sl)       ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_lay_width_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)        ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_w_lat_loss_sl)             ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)                    ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pwp_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_fc_sl)                 ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sat_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)                ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)                    ! inout
    ! ----------------------------------------------------------------------------------------------------- !

    jsb_standalone = jsbach_runs_standalone()

    !> 0.9 use variables from A2L_ or TURB_ - differ between standalone and coupled simulations
    !>
    ALLOCATE(drag_srf_local(nc))
    ALLOCATE(t_acoef_local(nc))
    ALLOCATE(t_bcoef_local(nc))
    ALLOCATE(q_acoef_local(nc))
    ALLOCATE(q_bcoef_local(nc))
    !$ACC ENTER DATA CREATE(drag_srf_local, t_acoef_local, t_bcoef_local, q_acoef_local, q_bcoef_local)
    IF (jsb_standalone) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        drag_srf_local(ic) = drag_srf_stdalone(ic)
        t_acoef_local(ic)  = t_acoef_stdalone(ic)
        t_bcoef_local(ic)  = t_bcoef_stdalone(ic)
        q_acoef_local(ic)  = q_acoef_stdalone(ic)
        q_bcoef_local(ic)  = q_bcoef_stdalone(ic)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        drag_srf_local(ic) = drag_srf(ic)
        t_acoef_local(ic)  = t_acoef(ic)
        t_bcoef_local(ic)  = t_bcoef(ic)
        q_acoef_local(ic)  = q_acoef(ic)
        q_bcoef_local(ic)  = q_bcoef(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !> 1.0 Compute snow accumulation (in the case of snowfall), used later in calc_surface_energy_balance()
    !!
    IF (dsl4jsb_Config(SPQ_)%flag_snow) THEN

      !$ACC UPDATE HOST(snow(:), snow_height(:), snow_present_snl(:,:), w_snow_snl(:,:), snow_lay_thickness_snl(:,:)) ASYNC(1)
      !$ACC WAIT(1)
      CALL calc_snow_accumulation(nc, &                                   ! in
                                  nsoil_w, &                              ! in
                                  nsnow, &                                ! in
                                  dtime, &                                ! in
                                  snow(:), &                              ! in
                                  snow_height(:), &                       ! inout
                                  snow_present_snl(:,:), &                ! inout
                                  w_snow_snl(:,:), &                      ! inout
                                  snow_lay_thickness_snl(:,:))            ! inout
      !$ACC UPDATE DEVICE(snow(:), snow_height(:), snow_present_snl(:,:), w_snow_snl(:,:), snow_lay_thickness_snl(:,:)) ASYNC(1)
#ifdef _OPENACC
      snow_height(:) = -1.0E+22; snow_lay_thickness_snl(:,:) = -2.0E+33
#endif
    END IF

    !> 2.0 Calculate surface energy budget
    !!
    !$ACC UPDATE HOST(rad_srf_net(:), t_air(:), t_acoef_local(:), t_bcoef_local(:), q_air(:)) &
    !$ACC   HOST(q_acoef_local(:), q_bcoef_local(:)) &
    !$ACC   HOST(press_srf(:), wind_10m(:), height(:), lai(:), wtr_skin(:), drag_srf_local(:)) &
    !$ACC   HOST(soil_depth_sl(:,:), heat_capa_sl(:,:), bulk_dens_sl(:,:), fact_q_air(:), fact_qsat_srf(:)) &
    !$ACC   HOST(fact_qsat_trans_srf(:)) &
    !$ACC   HOST(qsat_star(:), s_star(:), evapotranspiration(:), ground_heat_flx_old(:)) &
    !$ACC   HOST(zril_old(:), t(:), t_old(:), temp_srf_eff_4(:), evapopot(:), evaporation(:)) &
    !$ACC   HOST(evaporation_snow(:), interception(:), transpiration(:), sensible_heat_flx(:), latent_heat_flx(:)) &
    !$ACC   HOST(ground_heat_flx(:), snow_height(:), snow_soil_heat_flux(:), snow_srf_heat_flux(:)) &
    !$ACC   HOST(snow_melt_to_soil(:)) &
    !$ACC   HOST(therm_cond_sl(:,:), t_soil_sl(:,:), wtr_soil_sl(:,:), wtr_soil_fc_sl(:,:)) &
    !$ACC   HOST(wtr_soil_pwp_sl(:,:), ice_soil_sl(:,:), w_soil_freeze_flux(:,:), w_soil_melt_flux(:,:), t_snow_snl(:,:)) &
    !$ACC   HOST(w_snow_snl(:,:), snow_present_snl(:,:), snow_lay_thickness_snl(:,:), q_sat_srf(:)) &
    !$ACC   HOST(fract_wet(:), fract_snow(:), fract_snow_soil(:), fract_snow_can(:)) ASYNC(1)
    !$ACC WAIT(1)
    CALL calc_surface_energy_balance(nc, nsoil_w, nsnow, dtime, &
                                     lat, lon, &
                                     alpha, &
                                     jsb_standalone, &
                                     rad_srf_net(:), &
                                     t_air(:), t_acoef_local(:), t_bcoef_local(:), &
                                     q_air(:), q_acoef_local(:), q_bcoef_local(:), &
                                     press_srf(:), wind_10m(:), &
                                     height(:), lai(:), wtr_skin(:), drag_srf_local(:), &
                                     soil_depth_sl(:,:), heat_capa_sl(:,:), bulk_dens_sl(:,:), &                       ! in
                                     fact_q_air(:), fact_qsat_srf(:), &                                                ! inout
                                     fact_qsat_trans_srf(:), &                                                         ! inout
                                     qsat_star(:), s_star(:), &                                                        ! inout
                                     evapotranspiration(:), &                                                          ! inout
                                     ground_heat_flx_old(:), zril_old(:), &                                            ! inout
                                     t(:), t_old(:), temp_srf_eff_4(:), &                                              ! inout
                                     evapopot(:), evaporation(:), evaporation_snow, &                                  ! inout
                                     interception(:), transpiration(:), &                                              ! inout
                                     sensible_heat_flx(:), latent_heat_flx(:), ground_heat_flx(:), &                   ! inout
                                     snow_height(:), snow_soil_heat_flux(:), snow_srf_heat_flux (:), &                 ! inout
                                     snow_melt_to_soil(:), &                                                           ! inout
                                     therm_cond_sl(:,:), t_soil_sl(:,:), wtr_soil_sl(:,:), &                           ! inout
                                     wtr_soil_fc_sl(:,:), wtr_soil_pwp_sl(:,:), &                                      ! inout
                                     ice_soil_sl(:,:), &                                                               ! inout
                                     w_soil_freeze_flux(:,:), &                                                        ! inout
                                     w_soil_melt_flux(:,:), &                                                          ! inout
                                     t_snow_snl(:,:), &                                                                ! inout
                                     w_snow_snl(:,:), &                                                                ! inout
                                     snow_present_snl, &                                                               ! inout
                                     snow_lay_thickness_snl(:,:), &                                                    ! inout
                                     q_sat_srf(:), &                                                                   ! out
                                     fract_wet(:), &                                                                   ! out
                                     fract_snow(:), &                                                                  ! out
                                     fract_snow_soil(:), &                                                             ! out
                                     fract_snow_can(:))                                                                ! out
    !$ACC UPDATE DEVICE(rad_srf_net(:), t_air(:), t_acoef_local(:), t_bcoef_local(:), q_air(:)) &
    !$ACC   DEVICE(q_acoef_local(:), q_bcoef_local(:)) &
    !$ACC   DEVICE(press_srf(:), wind_10m(:), height(:), lai(:), wtr_skin(:), drag_srf_local(:)) &
    !$ACC   DEVICE(soil_depth_sl(:,:), heat_capa_sl(:,:), bulk_dens_sl(:,:), fact_q_air(:), fact_qsat_srf(:)) &
    !$ACC   DEVICE(fact_qsat_trans_srf(:)) &
    !$ACC   DEVICE(qsat_star(:), s_star(:), evapotranspiration(:), ground_heat_flx_old(:)) &
    !$ACC   DEVICE(zril_old(:), t(:), t_old(:), temp_srf_eff_4(:), evapopot(:), evaporation(:)) &
    !$ACC   DEVICE(evaporation_snow(:), interception(:), transpiration(:), sensible_heat_flx(:), latent_heat_flx(:)) &
    !$ACC   DEVICE(ground_heat_flx(:), snow_height(:), snow_soil_heat_flux(:), snow_srf_heat_flux(:)) &
    !$ACC   DEVICE(snow_melt_to_soil(:)) &
    !$ACC   DEVICE(therm_cond_sl(:,:), t_soil_sl(:,:), wtr_soil_sl(:,:), wtr_soil_fc_sl(:,:)) &
    !$ACC   DEVICE(wtr_soil_pwp_sl(:,:), ice_soil_sl(:,:), w_soil_freeze_flux(:,:), w_soil_melt_flux(:,:), t_snow_snl(:,:)) &
    !$ACC   DEVICE(w_snow_snl(:,:), snow_present_snl(:,:), snow_lay_thickness_snl(:,:), q_sat_srf(:)) &
    !$ACC   DEVICE(fract_wet(:), fract_snow(:), fract_snow_soil(:), fract_snow_can(:)) ASYNC(1)
#ifdef _OPENACC
    lai(:) = -1.0E+11; wtr_soil_fc_sl(:,:) = -2.0E+36
#endif

    !> 3.0 Calculate soil hydrology processes
    !!
    !$ACC UPDATE HOST(num_sl_above_bedrock(:), soil_lay_depth_ubound_sl(:,:), rain(:), snow(:), interception(:)) &
    !$ACC   HOST(evaporation(:), evaporation_snow(:), transpiration(:), w_soil_root_fc(:), w_soil_root_pwp(:)) &
    !$ACC   HOST(lai(:), root_depth(:), snow_melt_to_soil(:), soil_depth_sl(:,:), wtr_soil_sat_sl(:,:)) &
    !$ACC   HOST(wtr_soil_fc_sl(:,:), wtr_soil_pwp_sl(:,:), saxtonA(:,:), saxtonB(:,:), saxtonC(:,:)) &
    !$ACC   HOST(kdiff_sat_sl(:,:), root_fraction_sl(:,:), w_soil_freeze_flux(:,:), w_soil_melt_flux(:,:)) &
    !$ACC   HOST(wtr_skin(:), wtr_rootzone(:), wtr_rootzone_rel(:), wtr_plant_avail_rel(:), wtr_soil_root_pot(:)) &
    !$ACC   HOST(srf_runoff(:), drainage(:), gw_runoff(:), drainage_fraction(:), wtr_soil_sl(:,:), wtr_soil_pot_sl(:,:)) &
    !$ACC   HOST(drainage_sl(:,:), frac_wtr_transp_down_sl(:,:), gw_runoff_sl(:,:), t_soil_sl(:,:), ice_soil_sl(:,:)) &
    !$ACC   HOST(frac_w_lat_loss_sl(:,:), ftranspiration_sl(:,:)) ASYNC(1)
    !$ACC WAIT(1)
    CALL calc_soil_hydrology(nc, nsoil_w, dtime, &
                             num_sl_above_bedrock(:), &
                             soil_lay_depth_ubound_sl(:,:), &                                      ! (upper bound > lower bound !)
                             dsl4jsb_Config(VEG_)%flag_dynamic_roots, &
                             rain(:),snow(:),interception(:),evaporation(:),evaporation_snow(:),transpiration(:), &
                             w_soil_root_fc(:),w_soil_root_pwp(:),lai(:),root_depth(:), &
                             snow_melt_to_soil(:), soil_depth_sl(:,:), ftranspiration_sl(:,:), &
                             wtr_soil_sat_sl(:,:),wtr_soil_fc_sl(:,:),wtr_soil_pwp_sl(:,:), &
                             saxtonA(:,:),saxtonB(:,:),saxtonC(:,:),kdiff_sat_sl(:,:),root_fraction_sl(:,:),&
                             w_soil_freeze_flux(:,:), w_soil_melt_flux(:,:), &                                              ! in
                             wtr_skin(:), wtr_rootzone(:), wtr_rootzone_rel(:), wtr_plant_avail_rel(:), &                   ! inout
                             wtr_soil_root_pot(:), srf_runoff(:), drainage(:), gw_runoff(:), drainage_fraction(:), &        ! inout
                             wtr_soil_sl(:,:), wtr_soil_pot_sl(:,:), drainage_sl(:,:), frac_wtr_transp_down_sl(:,:), &      ! inout
                             gw_runoff_sl(:,:), t_soil_sl(:,:), ice_soil_sl(:,:), &                                         ! inout
                             frac_w_lat_loss_sl(:,:))                                                                       ! out
    !$ACC UPDATE DEVICE(num_sl_above_bedrock(:), soil_lay_depth_ubound_sl(:,:), rain(:), snow(:), interception(:)) &
    !$ACC   DEVICE(evaporation(:), evaporation_snow(:), transpiration(:), w_soil_root_fc(:), w_soil_root_pwp(:)) &
    !$ACC   DEVICE(lai(:), root_depth(:), snow_melt_to_soil(:), soil_depth_sl(:,:), wtr_soil_sat_sl(:,:)) &
    !$ACC   DEVICE(wtr_soil_fc_sl(:,:), wtr_soil_pwp_sl(:,:), saxtonA(:,:), saxtonB(:,:), saxtonC(:,:)) &
    !$ACC   DEVICE(kdiff_sat_sl(:,:), root_fraction_sl(:,:), w_soil_freeze_flux(:,:), w_soil_melt_flux(:,:)) &
    !$ACC   DEVICE(wtr_skin(:), wtr_rootzone(:), wtr_rootzone_rel(:), wtr_plant_avail_rel(:), wtr_soil_root_pot(:)) &
    !$ACC   DEVICE(srf_runoff(:), drainage(:), gw_runoff(:), drainage_fraction(:), wtr_soil_sl(:,:), wtr_soil_pot_sl(:,:)) &
    !$ACC   DEVICE(drainage_sl(:,:), frac_wtr_transp_down_sl(:,:), gw_runoff_sl(:,:), t_soil_sl(:,:), ice_soil_sl(:,:)) &
    !$ACC   DEVICE(frac_w_lat_loss_sl(:,:), ftranspiration_sl(:,:)) ASYNC(1)
#ifdef _OPENACC
    wtr_rootzone(:) = -1.0E+21; frac_w_lat_loss_sl(:,:) = -2.0E+14
#endif

    !$ACC WAIT(1)

    !$ACC EXIT DATA DELETE(drag_srf_local, t_acoef_local, t_bcoef_local, q_acoef_local, q_bcoef_local)
    DEALLOCATE(drag_srf_local, t_acoef_local, t_bcoef_local, q_acoef_local, q_bcoef_local)

  END SUBROUTINE update_spq_physics

  !> calc snow accumulation
  !!
  SUBROUTINE calc_snow_accumulation(nc, &
                                    nsoil_w, &
                                    nsnow, &
                                    dtime, &
                                    snow, &
                                    snow_height, &
                                    snow_present_snl, &
                                    w_snow_snl, &
                                    snow_lay_thickness_snl)

    USE mo_jsb_math_constants,        ONLY: eps8
    USE mo_spq_constants,             ONLY: w_density, snow_dens, albedo_snow, snow_height_min, w_snow_max

    IMPLICIT NONE
    ! --------------------------------------------
    ! 0.1 In
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimension chunk
                                                                       nsoil_w, &             !< dimension vertical soil layers
                                                                       nsnow                  !< dimension vertical snow layers
    REAL(wp),                                         INTENT(in)    :: dtime                  !< timestep length
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: snow                   !< Snowfall           [kg m-2 s-1]
    ! 0.2 Inout
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: snow_height            !< height (thickness) of all snow layers           [m]
    REAL(wp), DIMENSION(nc,nsnow),                    INTENT(inout) :: snow_present_snl, &    !< check if there is snow in every layer           [-]
                                                                       w_snow_snl, &          !< snow (volumetric) water content of snow layer   [m]
                                                                       snow_lay_thickness_snl !< snow depth of every layer                       [m]
    ! 0.3 Local
    INTEGER                                                         :: isnow
    REAL(wp), DIMENSION(nc)                                         :: hlp1,hlp2
    REAL(wp), DIMENSION(nc,nsnow)                                   :: w_snow_max_snl         !< maximum amount of water in snow layer at snow_dens [m]
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_snow_accumulation'
    !------------------------------------------------------------------------------------------------------ !

    !> 0.9 init local var
    !>
    hlp2(:)   = 0.0_wp

    !>1.0 Snow inputs
    ! some prerequisites
    w_snow_max_snl(:,:) = w_snow_max
    hlp1(:) = snow(:) * dtime / 1000._wp  ! recalculate snow to [m]

    ! rain below a threshold temperature could be added, but this should be done with reading in inputs, ask Soenke Zaehle
    ! not really needed if the forcing input is correctly put together

    !>1.1 Add new snow to snow layers
    DO isnow=1,nsnow-1
      WHERE (hlp1(:) > eps8)
        WHERE (w_snow_snl(:,isnow)<w_snow_max_snl(:,isnow)) ! find first layer where snow can be deposited
          hlp2(:)             = MIN(w_snow_max_snl(:,isnow) - w_snow_snl(:,isnow),hlp1(:))
          w_snow_snl(:,isnow) = w_snow_snl(:,isnow) + hlp2(:)
          hlp1(:)             = hlp1(:) - hlp2(:)
        END WHERE
      END WHERE
    END DO
    ! calculate seperately for top layer, since it is not limited in height
    WHERE (hlp1(:) > eps8)
      w_snow_snl(:,nsnow) = w_snow_snl(:,nsnow) + hlp1(:)
    END WHERE

    ! 1.2. Calculate snow height
    snow_height = 0.0_wp
    snow_lay_thickness_snl = 0.0_wp
    DO isnow=1,nsnow
      WHERE (w_snow_snl(:,isnow) > 0.0_wp)
        snow_present_snl(:,isnow)    = 1.0_wp ! there is snow in this layer
        snow_lay_thickness_snl(:,isnow) = w_snow_snl(:,isnow) * w_density / snow_dens
        snow_height(:)                  = snow_height(:) + snow_lay_thickness_snl(:,isnow)
      ELSEWHERE
        snow_present_snl(:,isnow) = 0.0_wp ! there isn't snow in this layer
      END WHERE
    END DO
  END SUBROUTINE calc_snow_accumulation


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_spq_physics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to provide an update of the soil temperature given the net radiation, surface drag, \n
  !!   surface wetness, and canopy conductance.
  !!
  !!
  !! Calculates surface energy fluxes and deep soil temperatures
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_surface_energy_balance(nc, nsoil_w, nsnow, dtime, &
                                         lat, lon, &
                                         alpha, &
                                         jsb_standalone, &
                                         rad_srf_net, &
                                         t_air, t_acoef, t_bcoef, &
                                         q_air, q_acoef, q_bcoef, &
                                         press_srf, wind_10m, &
                                         height, lai, wtr_skin, drag_srf, &
                                         soil_depth_sl, heat_capa_sl, bulk_dens_sl, &                                 ! in
                                         fact_q_air, fact_qsat_srf, &                                                 ! inout
                                         fact_qsat_trans_srf, &                                                       ! inout
                                         qsat_star, s_star, &                                                         ! inout
                                         evapotranspiration, &                                                        ! inout
                                         ground_heat_flx_old, zril_old, &                                             ! inout
                                         t, t_old, temp_srf_eff_4, &                                                  ! inout
                                         evapopot, evaporation, evaporation_snow, &                                   ! inout
                                         interception, transpiration, &                                               ! inout
                                         sensible_heat_flx, latent_heat_flx, ground_heat_flx, &                       ! inout
                                         snow_height, snow_soil_heat_flux, snow_srf_heat_flux, &                      ! inout
                                         snow_melt_to_soil, &                                                         ! inout
                                         therm_cond_sl,t_soil_sl,wtr_soil_sl, &                                       ! inout
                                         wtr_soil_fc_sl, wtr_soil_pwp_sl, &                                           ! inout
                                         ice_soil_sl, &                                                               ! inout
                                         w_soil_freeze_flux, w_soil_melt_flux, &                                      ! inout
                                         t_snow_snl, &                                                                ! inout
                                         w_snow_snl, &                                                                ! inout
                                         snow_present_snl, &                                                          ! inout
                                         snow_lay_thickness_snl, &                                                    ! inout
                                         q_sat_srf, &                                                                 ! out
                                         fract_wet, &                                                                 ! out
                                         fract_snow, &                                                                ! out
                                         fract_snow_soil, &                                                           ! out
                                         fract_snow_can)                                                              ! out

    USE mo_jsb_math_constants,        ONLY: eps8
    USE mo_atmland_constants,         ONLY: min_wind, max_wind
    USE mo_spq_constants,             ONLY: w_skin_max, soil_therm_cond, soil_frozen_therm_cond, latent_heat_fusion, &
                                            w_density, water2ice_density_ratio, fact_water_supercooled, w_snow_min, snow_dens, &
                                            snow_heat_capa, snow_therm_cond, t_soil_max_increase, t_soil_max_decrease
    USE mo_jsb_physical_constants,    ONLY: r_gas_dryair, LatentHeatSublimation, LatentHeatVaporization, &
                                            stbo, zemiss_def, grav, cpd, cpvd1, Tzero
    USE mo_phy_schemes,               ONLY: qsat_mixed, heat_transfer_coef, q_effective, surface_dry_static_energy


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_w, &
                                                                       nsnow
    REAL(wp),                                         INTENT(in)    :: dtime                  !< timestep length
    REAL(wp),                                         INTENT(in)    :: lat(:)                 !< timestep length
    REAL(wp),                                         INTENT(in)    :: lon(:)                 !< timestep length
    REAL(wp),                                         INTENT(in)    :: alpha                  !< implicitness factor
    LOGICAL,                                          INTENT(in)    :: jsb_standalone         !< standalone run?
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: rad_srf_net, &         !< surface net radiation [W/m2]
                                                                       t_air, &               !< air temperature [K]
                                                                       t_acoef, &             !< Richtmyer-Morton Coefficient Ta
                                                                       t_bcoef, &             !< Richtmyer-Morton Coefficient Tb
                                                                       q_air, &               !< air specific humidity [g/g]
                                                                       q_acoef, &             !< Richtmyer-Morton Coefficient Qa
                                                                       q_bcoef, &             !< Richtmyer-Morton Coefficient Qb
                                                                       press_srf, &           !< surface pressure [Pa]
                                                                       wind_10m, &            !< wind speed [m/2]
                                                                       height, &              !< height
                                                                       lai, &                 !< LAI
                                                                       wtr_skin, &            !< skin water reservoir on leaves [m]
                                                                       drag_srf               !< surface drag coefficient (JSBACH-style!)
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(in)    :: soil_depth_sl, &       !< soil depth per layer [m]
                                                                       heat_capa_sl, &        !< heat capacity of the soil layer [J/kg]
                                                                       bulk_dens_sl           !< bulk density of the soil layer [kg/m3]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: fact_q_air, &          !<
                                                                       fact_qsat_srf, &       !<
                                                                       fact_qsat_trans_srf, & !<
                                                                       qsat_star, &           !<
                                                                       s_star, &              !<
                                                                       evapotranspiration, &  !<
                                                                       ground_heat_flx_old, & !< ground heat flux from previous time step [W/m2]
                                                                       zril_old, &            !< previous timestep's Reynolds number
                                                                       t, &                   !< current surface temperature [K]
                                                                       t_old, &               !< previous timestep's surface temperature
                                                                       temp_srf_eff_4, &      !< effective surface temperature ** 4.0
                                                                       evapopot, &            !< potential evaporation rate [kg / m2 /s]
                                                                       evaporation, &         !< evaporation rate [kg / m2 /s]
                                                                       evaporation_snow, &
                                                                       interception, &        !< interception loss rate [kg / m2 /s]
                                                                       transpiration, &       !< transpiration rate [kg / m2 /s]
                                                                       sensible_heat_flx, &   !< sensible heat flux [W/m2]
                                                                       latent_heat_flx, &     !< latent heat flux [W/m2]
                                                                       ground_heat_flx, &     !< ground heat flux [W/m2]
                                                                       snow_height, &         !< snow depth as metric value [m]
                                                                       snow_soil_heat_flux, & !< snow soil heat flux [W/m2]
                                                                       snow_srf_heat_flux, &  !< snow heat flux from top snow layer to layer below [W / m2]
                                                                       snow_melt_to_soil      !< liquid water flux to soil [m]
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(inout) :: therm_cond_sl, &       !< thermal conductivity of the soil layer [J / m / s ]
                                                                       t_soil_sl, &           !< soil temperature [K]
                                                                       wtr_soil_sl, &         !< soil water content [m]
                                                                       wtr_soil_fc_sl, &      !< soil water field capacity [m]
                                                                       wtr_soil_pwp_sl, &     !< soil water pwp [m]
                                                                       ice_soil_sl, &         !< ice in water equivalent [m]
                                                                       w_soil_freeze_flux, &  !< water freezing flux in water equ.      [m]
                                                                       w_soil_melt_flux       !< ice melting flux in water equ.         [m]
    REAL(wp), DIMENSION(nc,nsnow),                    INTENT(inout) :: t_snow_snl, &          !< snow layer temperature [K]
                                                                       w_snow_snl, &          !< snow liquid water content equivalent [m]
                                                                       snow_present_snl, &    !< indicator of presence of snow  [-]
                                                                       snow_lay_thickness_snl !< snow layer depth                      [m]
    REAL(wp), DIMENSION(:),                           INTENT(out)   :: q_sat_srf              !< surface saturated specific humidity
    REAL(wp), DIMENSION(:),                           INTENT(out)   :: fract_wet              !< wet (skin and pond reservoir) fraction of tile (soil and canopy)
    REAL(wp), DIMENSION(:),                           INTENT(out)   :: fract_snow             !< snow fraction of tile
    REAL(wp), DIMENSION(:),                           INTENT(out)   :: fract_snow_soil        !< snow fraction on soil
    REAL(wp), DIMENSION(:),                           INTENT(out)   :: fract_snow_can         !< snow fraction on canopy
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                       ::     nsnow_top                            !< integer to determine top snow layer
    REAL(wp), DIMENSION(nc)                       ::     q_sat, &
                                                         heat_transfer_coefficient , &
                                                         w_snow_total , &
                                                         hlp1         , &
                                                         hlp2
    REAL(wp), DIMENSION(nc,nsoil_w)               :: ground_heat_flx_sl,&             !< Ground heat flux                                [J/m2/s]
                                                     t_soil_old_sl, &                 !< soil layer temperature from previous time step  [m]
                                                     latent_heat_freeze_flx, &        !< latent heat flux freezing                       [J/m2]
                                                     heat_excess, &                   !< Heat excess/deficit after freezing/melting      [J/m2]
                                                     latent_heat_thaw_flx, &          !< latent heat flux melting
                                                     hlp1_2d, &
                                                     hlp2_2d
    REAL(wp), DIMENSION(nc)                       :: dQdT,s_old,t2s_conv
    REAL(wp), DIMENSION(nc)                       :: q_air_eff,q_sat_srf_eff
    REAL(wp), DIMENSION(nc)                       :: evaporation_bare
    REAL(wp), DIMENSION(nc)                       :: srf_heat_capa,srf_depth,srf_dens,srf_heat_flx_old !< Local variables to define surface properties
                                                                                                       !  for implicit surface heat  exchange

    INTEGER                                       :: ichunk, ic     ! loop over chunk
    INTEGER                                       :: isoil, isnow   ! loop over layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_surface_energy_balance'
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.1 save old soil temperature, will be needed later
    !>
    t_soil_old_sl(:,:) = t_soil_sl(:,:)

    !> 0.8 init out var
    !>
    q_sat_srf(:)        = 0.0_wp
    fract_wet(:)        = 0.0_wp
    fract_snow(:)       = 0.0_wp
    fract_snow_soil(:)  = 0.0_wp
    fract_snow_can(:)   = 0.0_wp

    !> 0.9 init local var
    !>
    hlp1_2d   = 0.0_wp
    hlp2_2d   = 0.0_wp

    !> 1.0 calculate surface characteristic given old surface state
    !>
    DO ic = 1, nc
      q_sat(ic)      = qsat_mixed(t_air(ic), press_srf(ic))
      q_sat_srf(ic)  = qsat_mixed(t_old(ic), press_srf(ic))
    END DO

    !>   1.1 Calculate heat transfer coefficient
    !>
    DO ic = 1, nc
      heat_transfer_coefficient(ic) = heat_transfer_coef(drag_srf(ic), dtime, alpha)
    END DO

    !>   1.2 Calculate potential evaporation rate
    !>
    evapopot(:)      = -1._wp * heat_transfer_coefficient(:) * (q_air(:) - q_sat_srf(:))

    !>   1.3 calculate snow and water-saturated fraction of tile
    !>     following JSBACH4 hydrology
    !>
    !>   NOTE: this function corresponds to `mo_hydro_process:calc_wet_fractions_veg` and is functionally equivalent
    !>
    hlp1(:) = 0.0_wp
    CALL calc_wskin_fractions_veg(dtime, &
                                  w_skin_max * lai(:), &
                                  0.01_wp, &              ! IQ: oro_stddev from HYDRO_ init; QS: assume almost flat terrain
                                  evapopot(:), &
                                  wtr_skin(:), &
                                  w_snow_snl(:,1), &      ! jsb4: weq_snow_soil calc in HYDRO_ calc_surface_hydrology_land
                                  0.0_wp, &               ! jsb4: weq_snow_can  calc in HYDRO_ calc_surface_hydrology_land
                                  fract_snow_can(:), &    ! out
                                  fract_wet(:), &         ! out
                                  fract_snow_soil(:))     ! out
    ! for consistency with jsbach: 'fract_snow(ic) = fract_snow_soil(ic)' in HYDRO_ interface
    fract_snow(:) = fract_snow_soil(:)

    !> 2.0 JSBACH4: update_surface_energy_land
    !>

    !>   2.1 Compute sensitivity of surface saturation specific humidity to temperature
    !>
    DO ic = 1, nc
      dQdT(ic) = (qsat_mixed(t_old(ic) + 0.001_wp, press_srf(ic)) - q_sat_srf(ic)) * 1000._wp
    END DO

    !>   2.2 Old dry static energy
    !>
    DO ic = 1, nc
      s_old(ic) = surface_dry_static_energy( &
        & t_old(ic), &
        & q_effective(q_sat_srf(ic), q_air(ic), fact_qsat_srf(ic), fact_q_air(ic)), &
        & cpd, jsb_standalone)
    END DO
    t2s_conv(:) = s_old(:) / t_old(:)

    !>   2.3 solve energy balance calculation using implicit formulation (JSBACH4)
    !>
    ! calculate surface values depending on snow/bare coverage for implicit surface exchange scheme
    DO ichunk=1,nc
      nsnow_top = MAX(NINT(SUM(snow_present_snl(ichunk,:))), 1)
      IF (w_snow_snl(ichunk,1) > w_snow_min) THEN
        srf_heat_capa(ichunk) = snow_heat_capa * fract_snow_soil(ichunk) + (1 - fract_snow_soil(ichunk)) * heat_capa_sl(ichunk,1)
        srf_dens(ichunk)  = snow_dens * fract_snow_soil(ichunk) + (1 - fract_snow_soil(ichunk)) * bulk_dens_sl(ichunk,1)
        srf_depth(ichunk) = snow_lay_thickness_snl(ichunk,nsnow_top) &
          & * fract_snow_soil(ichunk) + (1 - fract_snow_soil(ichunk)) * soil_depth_sl(ichunk,1)
        srf_heat_flx_old(ichunk) = snow_srf_heat_flux(ichunk) &
          & * fract_snow_soil(ichunk) + (1 - fract_snow_soil(ichunk)) * ground_heat_flx(ichunk)
      ELSE
        srf_heat_capa(ichunk) = heat_capa_sl(ichunk,1)
        srf_dens(ichunk)  = bulk_dens_sl(ichunk,1)
        srf_depth(ichunk) = soil_depth_sl(ichunk,1)
        srf_heat_flx_old(ichunk) = ground_heat_flx(ichunk)
      ENDIF
    ENDDO

    hlp1 = calc_temp_srf_new_implicit(dtime, &
                                      alpha, &
                                      t2s_conv(:), &
                                      t_acoef(:), &
                                      t_bcoef(:), &
                                      q_acoef(:), &
                                      q_bcoef(:), &
                                      s_old(:), &
                                      q_sat_srf(:), &
                                      dQdT(:), &
                                      rad_srf_net(:), &
                                      srf_heat_flx_old(:), &
                                      heat_transfer_coefficient(:), &
                                      fact_q_air(:), &
                                      fact_qsat_srf(:), &
                                      fract_snow_soil(:), &
                                      srf_heat_capa(:) * srf_dens(:) * srf_depth(:))

    ! time-filtering as in JSBACH3
    s_star(:)    = 1._wp / alpha * hlp1(:) + (1._wp - 1._wp / alpha) * s_old(:)
    t(:)         = s_star(:) / t2s_conv(:)
    qsat_star(:) = q_sat_srf(:) + dQdT(:) * (s_star(:) - s_old(:)) / t2s_conv(:)

    ! debug output
    ! write lon/lat of grid cell with surface temperature below 100 K
    ! DO ic = 1, nc
    !   IF (t(ic) < 100.0_wp) THEN
    !     WRITE (message_text,*) 't below 100 K at ', &
    !       & lat(ic), 'N and ', lon(ic), 'E:', NEW_LINE('a'), &
    !       & '  surface temperature:   ', t(ic)
    !       CALL message (routine, message_text, all_print = .TRUE.)
    !   END IF
    ! END DO

    ! temporary bugfix after implementing TURB_ with quincy (and using shared/mo_phy_schemes:update_drag -> calc_qsat)
    ! limit surface temperature to [50,400] Kelvin to avoid potential segmentation fault in shared/mo_phy_schemes:calc_qsat
    !   consistently with "seb_check_temperature_range()"
    t(:) = MAX(50.0_wp, t(:))
    t(:) = MIN(400.0_wp, t(:))

    !> 3.0 calculate water fluxes
    !>   This is from JSBACH4 hydrology/mo_hydro_interface:update_evaporation
    !>
    DO ic = 1, nc
      q_air_eff(ic)       = q_effective(0.0_wp, q_air(ic), 1.0_wp, 0.0_wp)
      q_sat_srf_eff(ic)   = q_effective(qsat_star(ic), q_air(ic), fact_qsat_srf(ic), fact_q_air(ic))
    END DO
    evapotranspiration(:) = heat_transfer_coefficient(:) * (q_air_eff(:) - q_sat_srf_eff(:))
    evapopot(:)           = heat_transfer_coefficient(:) * (q_air(:) - qsat_star(:))

    !> 3.1 convert to positive fluxes for water module
    !! needed for QUINCY, but not JSBACH4
    transpiration(:)    = -1._wp * fact_qsat_trans_srf(:) * evapopot(:)
    interception(:)     = -1._wp * (1._wp - fract_snow_soil(:)) * fract_wet(:) * evapopot(:)
    evaporation_snow(:) = -1._wp * fract_snow_soil(:) * evapopot(:)
    ! interception evaporation must be less than the total skin water
    interception(:) = MIN(interception(:) * dtime / 1000._wp, wtr_skin(:)) * 1000._wp / dtime
    ! snow evaparation must be less than the total snow water
    evaporation_snow(:) = MIN(evaporation_snow(:) * dtime / 1000._wp, SUM(w_snow_snl(:,:),DIM=2)) * 1000._wp / dtime

    ! no dew into stomates or onto canopy
    DO ic = 1,nc
      IF (evapopot(ic) > 0.0_wp) THEN
         transpiration(ic)    = 0.0_wp
         interception(ic)     = 0.0_wp
         evaporation_snow(ic) = 0.0_wp
      END IF
    END DO
    evaporation_bare(:) = -1._wp * evapotranspiration(:) - transpiration(:) - interception(:) - evaporation_snow(:)
    ! BPLARR @TODO, as long as there is no snow
    evaporation = evaporation_bare  ! there should be an output stream for evaporation_snow (i.e. sublimation), but
                                    ! but this cannot be added here (!), because this goes into soil hydrology calculations

    !> 4.0 calculate sensible and latent heat flux
    !>  from JSBACH4 srf_energy_bal/mo_seb_land:update_surface_fluxes_land
    !>

    ! Compute sensible heat flux
    sensible_heat_flx(:) = heat_transfer_coefficient(:) * (t_acoef(:) * s_old(:) + t_bcoef(:) - s_star(:))

    ! Compute latent heat flux
    latent_heat_flx(:) = LatentHeatVaporization * evapotranspiration(:) &
      &                  + (LatentHeatVaporization - LatentHeatSublimation) * evaporation_snow(:)

    !> 5.0 update surface temperature and remember previous value
    !>

    ! call calc_snow_dynamics to calculate:
    !  snow top layer temperature,
    !  heat transfer through snow, and
    !  soil surface temperature
    CALL calc_snow_dynamics(nc, &                   !in
                             nsoil_w,&              !in
                             nsnow,  &              !in
                             dtime,  &              !in
                             t_air, &               !in
                             t_old(:), &            !in
                             t(:), &                !in
                             evaporation_snow(:), & !in
                             heat_capa_sl(:,:), &   !in
                             bulk_dens_sl(:,:), &   !in
                             soil_depth_sl(:,:), &  !in
                             snow_height(:), &          !inout
                             snow_soil_heat_flux(:), &  !inout
                             snow_srf_heat_flux(:), &   !inout
                             snow_melt_to_soil(:), &    !inout
                             t_soil_sl(:,:), &          !inout
                             t_snow_snl(:,:), &         !inout
                             w_snow_snl(:,:), &         !inout
                             snow_present_snl(:,:), &   !inout
                             snow_lay_thickness_snl)    !inout

    ! use surface temperature as calculated in calc_snow_dynamics()
    ! if there is no snow: overwrite t_soil_sl(:,1) with t(:)
    ! modify only soil layer: 1
    DO ic = 1,nc
      IF (w_snow_snl(ic,1) < w_snow_min) THEN
        t_soil_sl(ic,1)      = t(ic)
      END IF
    END DO

    !>   5.1 calculate ground heat flux and update deep surface temperatures
    !>
    ground_heat_flx_sl(:,:) = 0.0_wp
    DO ic = 1,nc
      DO isoil = 1,nsoil_w
        IF (t_soil_sl(ic,isoil) <= Tzero) THEN
          therm_cond_sl(ic,isoil) = soil_frozen_therm_cond
        ELSE
          therm_cond_sl(ic,isoil) = soil_therm_cond
        END IF
      END DO
    END DO

    DO ic = 1,nc
      DO isoil = 1,nsoil_w-1
        IF (soil_depth_sl(ic,isoil) > eps8) THEN
          ground_heat_flx_sl(ic,isoil) = therm_cond_sl(ic,isoil) &
            &                           * (t_soil_sl(ic,isoil + 1) - t_soil_sl(ic,isoil)) &
            &                           / ((soil_depth_sl(ic,isoil) + soil_depth_sl(ic,isoil + 1)) / 2.0_wp)
        END IF
      END DO
    END DO

    ground_heat_flx(:)     = ground_heat_flx_sl(:,1)
    ground_heat_flx_old(:) = ground_heat_flx(:)
    ! loop soil layers 2:max
    DO ic = 1,nc
      DO isoil = 2,nsoil_w
        IF(soil_depth_sl(ic,isoil) > eps8) THEN
          t_soil_sl(ic,isoil) = t_soil_sl(ic,isoil) &
            &                  + (ground_heat_flx_sl(ic,isoil) - ground_heat_flx_sl(ic,isoil - 1)) * dtime &
            &                  / (heat_capa_sl(ic,isoil) * bulk_dens_sl(ic,isoil) * soil_depth_sl(ic,isoil))
          ! limit temperature changes to a certain maximum
          ! TODO https://gitlab.dkrz.de/jsbach/jsbach/-/issues/192
          ! increase by more than t_soil_max_increase
          IF ((t_soil_sl(ic,isoil) - t_soil_old_sl(ic,isoil)) > t_soil_max_increase) THEN
            t_soil_sl(ic,isoil) = t_soil_old_sl(ic,isoil) + t_soil_max_increase
          ! decrease by more than t_soil_max_decrease (a negative value)
          ELSE IF ((t_soil_sl(ic,isoil) - t_soil_old_sl(ic,isoil)) < t_soil_max_decrease) THEN
            t_soil_sl(ic,isoil) = t_soil_old_sl(ic,isoil) + t_soil_max_decrease
          END IF
        END IF
      END DO
    END DO

    !>   5.2 Check for freezing temperatures & compute freezing dynamics
    !>
    w_soil_freeze_flux(:,:) = 0.0_wp
    w_soil_melt_flux(:,:)   = 0.0_wp
    DO ic = 1,nc
      DO isoil = 1,nsoil_w
        ! freezing
        IF (t_soil_sl(ic,isoil) <= Tzero) THEN
          ! if there is water available for freezing
          ! calculate heat deficit past freezing point
          IF ((wtr_soil_sl(ic,isoil) - fact_water_supercooled * wtr_soil_pwp_sl(ic,isoil)) > eps8) THEN
            hlp1_2d(ic,isoil) = heat_capa_sl(ic,isoil) * bulk_dens_sl(ic,isoil) * soil_depth_sl(ic,isoil) &
              &                 * (t_soil_sl(ic,isoil) - Tzero) ! heat flux deficit
            ! Maximum latent heat release from freezing per soil layer
            hlp2_2d(ic,isoil) = latent_heat_fusion * w_density &
              &                 * (wtr_soil_sl(ic,isoil)- fact_water_supercooled * wtr_soil_pwp_sl(ic,isoil))
            ! if heat deficit is less than heat release of fusion of available layer soil water
            IF (-hlp1_2d(ic,isoil) <= hlp2_2d(ic,isoil)) THEN
              latent_heat_freeze_flx(ic,isoil) = hlp1_2d(ic,isoil)
              w_soil_freeze_flux(ic,isoil)     = MIN(wtr_soil_sl(ic,isoil) - fact_water_supercooled * wtr_soil_pwp_sl(ic,isoil), &
                &                                - latent_heat_freeze_flx(ic,isoil) / (w_density * latent_heat_fusion))
              w_soil_freeze_flux(ic,isoil)     = MIN(w_soil_freeze_flux(ic,isoil), &
                &                                wtr_soil_fc_sl(ic,isoil) - ice_soil_sl(ic,isoil) &
                &                                * water2ice_density_ratio - fact_water_supercooled * wtr_soil_pwp_sl(ic,isoil))
              IF (t_soil_old_sl(ic,isoil) >= Tzero) THEN
                t_soil_sl(ic,isoil) = Tzero
              ! if temperature is already below freezing point, stays the same
              ELSE
                t_soil_sl(ic,isoil) = t_soil_sl(ic,isoil)
              END IF
            ! if heat deficit is more than heat release of fusion of all soil water
            ELSE
              w_soil_freeze_flux(ic,isoil)     = MIN(wtr_soil_sl(ic,isoil) - fact_water_supercooled * wtr_soil_pwp_sl(ic,isoil), &
                &                                wtr_soil_fc_sl(ic,isoil) - ice_soil_sl(ic,isoil) &
                &                                * water2ice_density_ratio - fact_water_supercooled &
                &                                * wtr_soil_pwp_sl(ic,isoil))
              latent_heat_freeze_flx(ic,isoil) = latent_heat_fusion * w_density * w_soil_freeze_flux(ic,isoil)
              IF (t_soil_old_sl(ic,isoil) >= Tzero) THEN
                heat_excess(ic,isoil) = hlp1_2d(ic,isoil) + hlp2_2d(ic,isoil)
                ! Soil T cooles past freezing point
                t_soil_sl(ic,isoil)   = Tzero + heat_excess(ic,isoil) / (heat_capa_sl(ic,isoil) * bulk_dens_sl(ic,isoil) &
                  &                     * soil_depth_sl(ic,isoil))
              ELSE
                IF (t_soil_sl(ic,isoil) < t_soil_old_sl(ic,isoil)) THEN
                  t_soil_sl(ic,isoil) = t_soil_sl(ic,isoil) + latent_heat_freeze_flx(ic,isoil) &
                    &                   / (heat_capa_sl(ic,isoil) * bulk_dens_sl(ic,isoil) * soil_depth_sl(ic,isoil))
                END IF
              END IF
            END IF
          END IF
        ! thawing (t_soil_sl(ic,isoil) > Tzero)
        ELSE
          ! thaw when soil T below freezing point and ice exists
          ! same as above the other way around
          IF (ice_soil_sl(ic,isoil) > eps8) THEN
            ! heat flux excess
            hlp1_2d(ic,isoil) = heat_capa_sl(ic,isoil) * bulk_dens_sl(ic,isoil) * soil_depth_sl(ic,isoil) &
              &                 * (t_soil_sl(ic,isoil) - Tzero)
            hlp2_2d(ic,isoil) = latent_heat_fusion * w_density * ice_soil_sl(ic,isoil)
            ! Excess heat is less than latent heat of layer ice melting, -> temperature remains at Tzero
            IF (hlp1_2d(ic,isoil) <= hlp2_2d(ic,isoil)) THEN
              latent_heat_thaw_flx(ic,isoil) = hlp1_2d(ic,isoil)
              w_soil_melt_flux(ic,isoil)     = MIN(latent_heat_thaw_flx(ic,isoil) &
                &                              / (w_density * latent_heat_fusion), ice_soil_sl(ic,isoil)) ! unit of f_thw is [m] water not [m] ice(!)
              t_soil_sl(ic,isoil)            = Tzero
            ! Excess heat is more than latent heat release of layer ice melting -> temperature continues to increase after ice melting
            ELSE
              heat_excess(ic,isoil)      = hlp1_2d(ic,isoil) - latent_heat_fusion * w_density * ice_soil_sl(ic,isoil)
              w_soil_melt_flux(ic,isoil) = ice_soil_sl(ic,isoil)
              t_soil_sl(ic,isoil)        = Tzero + heat_excess(ic,isoil) / (heat_capa_sl(ic,isoil) * bulk_dens_sl(ic,isoil) &
                &                          * soil_depth_sl(ic,isoil))
            END IF
          END IF
        END IF   ! freezing/thawing
      END DO     ! soil layer
    END DO       ! nc

    ! diagnose effective surface temperature ** 4.0 (needed for land2atm)
    DO ic = 1,nc
      temp_srf_eff_4(ic) = t_old(ic) ** 3.0_wp * (4.0_wp * t(ic) - 3.0_wp * t_old(ic))
    END DO

    ! remember old surface temperature for next timestep
    ! this depends whethere there is snow (in snow layer 1, i.e., the one at the soil surface) or not
    DO ic = 1,nc
      IF (w_snow_snl(ic,1) > w_snow_min) THEN
        nsnow_top = MAX(NINT(SUM(snow_present_snl(ic,:))), 1)
        t_old(ic) = t_snow_snl(ic,nsnow_top)
      ELSE
        t_old(ic) = t_soil_sl(ic,1)
      END IF
    END DO
  END SUBROUTINE calc_surface_energy_balance

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_spq_physics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to provide an update of the soil water given rates of rain and snow fall \n
  !!  as well as interception loss, evaporation and transpiration.
  !!
  !! Calculates srf_runoff, drainage and an updated state of soil water variables (theta, potential, etc.)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_soil_hydrology(nc, nsoil_w, dtime, &
                                 num_sl_above_bedrock, &
                                 soil_layer_w_ubounds, &
                                 flag_dynamic_roots, &
                                 rain,snow,interception,evaporation,evaporation_snow,transpiration, &
                                 w_soil_root_fc,w_soil_root_pwp,lai,root_depth, &
                                 snow_melt_to_soil, &
                                 soil_depth_sl, ftranspiration_sl, &
                                 wtr_soil_sat_sl,wtr_soil_fc_sl,wtr_soil_pwp_sl, &
                                 saxtonA,saxtonB,saxtonC,kdiff_sat_sl,root_fraction_sl, &
                                 w_soil_freeze_flux, w_soil_melt_flux, &                           ! in
                                 wtr_skin, wtr_rootzone, wtr_rootzone_rel, wtr_plant_avail_rel, &  ! inout
                                 wtr_soil_root_pot, srf_runoff, drainage, gw_runoff, &             ! inout
                                 drainage_fraction, wtr_soil_sl, wtr_soil_pot_sl, drainage_sl, &   ! inout
                                 frac_wtr_transp_down_sl, &                                        ! inout
                                 gw_runoff_sl, t_soil_sl, ice_soil_sl, &                           ! inout
                                 frac_w_lat_loss_sl)                                               ! out

    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_jsb_physical_constants,      ONLY: grav, Tzero
    USE mo_spq_constants,               ONLY: drain_frac, fdrain_srf_runoff, &
                                              InterceptionEfficiency, w_skin_max, &
                                              wpot_min, water2ice_density_ratio, &
                                              frac_wtr_vertical_transport_up_max, frac_wtr_vertical_transport_down_max, &
                                              frac_w_lat_loss_max

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                   !< dimensions
                                                                       nsoil_w                 !<
    REAL(wp),                                         INTENT(in)    :: dtime                   !< timestep length
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: num_sl_above_bedrock    !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(in)    :: soil_layer_w_ubounds    !< upper bound of soil layers water (upper bound > lower bound !)
    LOGICAL,                                          INTENT(in)    :: flag_dynamic_roots      !< VEG_ config: use dynmic root scheme if TRUE
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: rain, &                 !< rainfall [kg / m2 / s ]
                                                                       snow, &                 !< snowfall [kg / m2 / s ]
                                                                       interception, &         !< interception [kg / m2 / s ]
                                                                       evaporation, &          !< evaporation [kg / m2 / s ]
                                                                       evaporation_snow, &     !< evaporation snow [kg / m2 / s ]
                                                                       transpiration, &        !< transpiration [kg / m2 / s ]
                                                                       w_soil_root_fc, &       !< root moisture at field capacity [m]
                                                                       w_soil_root_pwp, &      !< root moisture at wilting point [m]
                                                                       lai, &                  !< leaf area index
                                                                       root_depth, &           !< rooting depth [m]
                                                                       snow_melt_to_soil       !< snow melt water flux to soil in liquid water equivalent [m]
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(in)    :: soil_depth_sl, &        !< soil depth per layer [m]
                                                                       ftranspiration_sl, &    !< fraction of transpiration []
                                                                       wtr_soil_sat_sl, &      !< soil moisture per layer at saturation [m]
                                                                       wtr_soil_fc_sl, &       !< soil moisture per layer at field capacity [m]
                                                                       wtr_soil_pwp_sl, &      !< soil moisture per layer at wilting point [m]
                                                                       kdiff_sat_sl, &         !< saturated water diffusivity [m/s]
                                                                       saxtonA, &              !< pedotransfer parameter
                                                                       saxtonB, &              !< pedotransfer parameter
                                                                       saxtonC, &              !< pedotransfer parameter
                                                                       root_fraction_sl, &     !< fraction of root mass per soil layer
                                                                       w_soil_freeze_flux, &   !< water flux from liquid to ice [m]
                                                                       w_soil_melt_flux        !< water flux from ice to liquid [m]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: wtr_skin, &             !< water in skin reservoir [m]
                                                                       wtr_rootzone, &         !< water in root zone [m]
                                                                       wtr_rootzone_rel, &     !< fractional water availability in root zone [-]
                                                                       wtr_plant_avail_rel, &  !< fractional plant available water availability in root zone [-]
                                                                       wtr_soil_root_pot, &    !< soil water potential in root zone [MPa]
                                                                       srf_runoff, &           !< surface runoff [kg / m2 / s]
                                                                       drainage, &             !< deep drainage [kg / m2 / s]
                                                                       gw_runoff, &            !< sum of lateral runoff across all soil layers [kg / m2 /s]
                                                                       drainage_fraction       !< fraction of soil column lost to drainage [1/s]
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(inout) :: wtr_soil_sl, &                     !< soil moisture per layer [m]
                                                                       wtr_soil_pot_sl, &                 !< soil water potential per layer [MPa]
                                                                       drainage_sl, &                     !< drainage per layer [kg / m2 / s]
                                                                       frac_wtr_transp_down_sl, &         !< fractional water loss [1000/s]
                                                                       gw_runoff_sl, &                    !< groundwater runoff per layer [kg / m2 /s]
                                                                       t_soil_sl, &                       !< soil layer temperature
                                                                       ice_soil_sl                        !< ice content of soil layer [m]
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(out)   :: frac_w_lat_loss_sl                 !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                         :: &
                                                       infiltration , &
                                                       water_to_skin, &
                                                       water_to_soil, &
                                                       w_soil_root_old, &
                                                       hlp1, hlp2, hlp3
    REAL(wp), DIMENSION(nc)                         :: fdrain_srf_runoff_mod        ! modified fdrain_srf_runoff
    REAL(wp),DIMENSION(nc,nsoil_w)                  :: &
                                                       fevaporation   , &           ! fraction of evaporation
                                                       w_soil_sl_old  , &
                                                       k_diff         , &
                                                       hlp1_2d
    INTEGER                                         :: isoil, ichunk, ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_hydrology'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.8 safe old water content for conservation tests
    !!
    w_soil_root_old(:) = wtr_rootzone(:)
    w_soil_sl_old(:,:) = wtr_soil_sl(:,:)

    !> 0.9 init local var
    !>
    hlp1                     = 0.0_wp
    hlp2                     = 0.0_wp
    hlp3                     = 0.0_wp
    hlp1_2d                  = 0.0_wp

    !> 1.1 Update skin reservoir and calculate interception loss
    !!
    WHERE((w_skin_max * lai(:)) > eps8)
      ! m / timestep
      water_to_skin(:) = MIN(InterceptionEfficiency * rain(:) / 1000.0_wp * dtime, &
        &                 MAX(interception(:) * dtime / 1000.0_wp + w_skin_max * lai(:) - wtr_skin(:), 0._wp))
    ELSEWHERE
      water_to_skin(:) = 0.0_wp
    ENDWHERE
    water_to_soil(:) = MAX(rain(:) + snow_melt_to_soil(:) - water_to_skin(:) / dtime * 1000.0_wp, 0._wp)
    wtr_skin(:)      = MAX(wtr_skin(:) + water_to_skin(:) - interception(:) * dtime / 1000.0_wp, 0.0_wp)

    !>  1.2. Freeze thaw water fluxes
    !!
    wtr_soil_sl(:,:) = wtr_soil_sl(:,:) + w_soil_melt_flux(:,:) - w_soil_freeze_flux(:,:)
    ice_soil_sl(:,:) = ice_soil_sl(:,:) - w_soil_melt_flux(:,:) + w_soil_freeze_flux(:,:)

    !> 2.0 update soil moisture
    !!
    fevaporation(:,:) = 0.0_wp
    fevaporation(:,1) = 1.0_wp

    !> 2.1 fraction of transpiration per layer is now determined in assimi
    !!

    !> 2.2 determine water balance of the top layer
    !!
    !! infiltration rate (limited to attaining field capacity), surface runoff and deep drainage
    !! a fraction of srf_runoff generation is leaked below the first layer
    !! a fraction of infiltration is leaked below the first layer (preferential flow)
    ! first guess

    infiltration(:) = MAX(MIN((wtr_soil_fc_sl(:,1) - (wtr_soil_sl(:,1) + ice_soil_sl(:,1) * water2ice_density_ratio)) &
      &               * 1000._wp / dtime                                                                       &
      &               + (evaporation(:) * fevaporation(:,1) + transpiration(:) * ftranspiration_sl(:,1)),      &
      &                water_to_soil(:)), 0.0_wp)
    srf_runoff(:)   = MAX(water_to_soil(:) - infiltration(:), 0.0_wp)

    DO ichunk = 1,nc
       ! If soil is frozen, then stop draingage to deeper layers
       IF (t_soil_sl(ichunk,1) <= Tzero) THEN
         hlp1(ichunk) = 0.0_wp ! fraction of infiltration that is draining
         fdrain_srf_runoff_mod(ichunk) = 0.0_wp ! fraction of surface runoff that is draining
         drainage_sl(ichunk,1) = 0.0_wp ! drainage
       ELSE
         hlp1(ichunk) = drain_frac ** soil_depth_sl(ichunk,1) ! fraction of infiltration that is draining
         fdrain_srf_runoff_mod(ichunk) = fdrain_srf_runoff    ! fraction of surface runoff that is draining
         drainage_sl(ichunk,1) = hlp1(ichunk) * infiltration(ichunk) + fdrain_srf_runoff_mod(ichunk) * srf_runoff(ichunk)
       END IF
       ! update the first soil layer water content and surface runoff
       wtr_soil_sl(ichunk,1) = wtr_soil_sl(ichunk,1)                                    &
         &                   + ((1._wp - hlp1(ichunk)) * infiltration(ichunk)       &
         &                   - evaporation(ichunk) * fevaporation(ichunk,1)         &
         &                   - transpiration(ichunk) * ftranspiration_sl(ichunk,1)) &
         &                   / 1000._wp * dtime
       srf_runoff(ichunk)  = (1.0_wp - fdrain_srf_runoff_mod(ichunk)) * srf_runoff(ichunk)
    END DO

    gw_runoff_sl(:,1) = 0.0_wp
    ! Check if water content (including ice volume) is above field capacity
    WHERE(wtr_soil_sl(:,1) + ice_soil_sl(:,1) * water2ice_density_ratio > wtr_soil_fc_sl(:,1))
      gw_runoff_sl(:,1) = (wtr_soil_sl(:,1) + ice_soil_sl(:,1) &
        &                  * water2ice_density_ratio - wtr_soil_fc_sl(:,1)) * 1000._wp / dtime
      wtr_soil_sl(:,1) = wtr_soil_fc_sl(:,1) - ice_soil_sl(:,1) * water2ice_density_ratio
    END WHERE

    !> 2.3 update deeper soil layers given infiltration, transpiration and evaporation
    !> 2.3.1 Compute infiltration, transpiration, evaporation and drainage
    DO isoil = 2,nsoil_w-1
      ! turn off drainage if layer below is frozen/saturated
      hlp1(:) = drain_frac ** soil_depth_sl(:,isoil) * drainage_sl(:,isoil-1)
      hlp2(:) = SUM(wtr_soil_fc_sl(:,isoil:nsoil_w),DIM=2) - (SUM(wtr_soil_sl(:,isoil:nsoil_w),DIM=2) + &
        &       SUM(water2ice_density_ratio * ice_soil_sl(:,isoil:nsoil_w),DIM=2))
      WHERE (t_soil_sl(:,isoil) <= Tzero)
        ! Available pore space in layer below
        drainage_sl(:,isoil)  = 0.0_wp
        gw_runoff_sl(:,isoil) = 0.0_wp
      ELSEWHERE (t_soil_sl(:,nsoil_w) <= Tzero.AND.hlp1(:)>hlp2(:))
          drainage_sl(:,isoil) =  MAX(hlp2(:),0.0_wp)
      ELSEWHERE
        ! drainage from this layer due to preferential flow of infiltering water flux
        drainage_sl(:,isoil)  = drain_frac ** soil_depth_sl(:,isoil) * drainage_sl(:,isoil-1)
      END WHERE
      ! groundwater runoff
      gw_runoff_sl(:,isoil) = 0.0_wp
      ! update water balance of lower soil layers
      wtr_soil_sl(:,isoil)    = wtr_soil_sl(:,isoil) + &
                              (drainage_sl(:,isoil-1) - drainage_sl(:,isoil) - &
                              evaporation(:) * fevaporation(:,isoil) - &
                              transpiration(:) * ftranspiration_sl(:,isoil)) / &
                              1000._wp * dtime

      ! Compute flows to avoid wtr_soil_sl being above field capacity
      ! Check if water content (including ice volume) is above field capacity
      WHERE((wtr_soil_sl(:,isoil) + ice_soil_sl(:,isoil) * water2ice_density_ratio) > wtr_soil_fc_sl(:,isoil))
        ! Calculate excess water in current soil layer (over field capacity)
        hlp1(:) = wtr_soil_sl(:,isoil) + ice_soil_sl(:,isoil) * water2ice_density_ratio - wtr_soil_fc_sl(:,isoil)
        ! Calculate  remaining free water uptake capacity in all lower layers
        hlp2(:) = SUM(wtr_soil_fc_sl(:,isoil:nsoil_w),DIM=2) - (SUM(wtr_soil_sl(:,isoil:nsoil_w),DIM=2) + &
          &       SUM(water2ice_density_ratio * ice_soil_sl(:,isoil:nsoil_w),DIM=2))
        ! Calculate remaining free water uptake capacity in upper layer
        hlp3(:) = wtr_soil_fc_sl(:,isoil-1) - (wtr_soil_sl(:,isoil-1) + ice_soil_sl(:,isoil-1) * water2ice_density_ratio)
        ! Aggregated water capacity of lower layers has been reached either from liquid or frozen water (impermeable layer)
        WHERE (hlp2(:) <= 0.0_wp)
          WHERE (hlp3(:) < hlp1(:))  ! Water excess over saturation of current layer exceeds free capacity in upper layer
            ! Upper layer takes up water until it reaches field capacity, rest goes to lateral groundwater
            wtr_soil_sl(:,isoil-1)  = MIN(wtr_soil_sl(:,isoil-1) + hlp1(:),wtr_soil_fc_sl(:,isoil-1) - ice_soil_sl(:,isoil-1) * &
              &                           water2ice_density_ratio)!
            drainage_sl(:,isoil-1)  = drainage_sl(:,isoil-1) - hlp3(:) * 1000._wp / dtime
            gw_runoff_sl(:,isoil)   = gw_runoff_sl(:,isoil) + (hlp1(:) - hlp3(:)) * 1000._wp / dtime
          ELSEWHERE ! Upper layer has capacity to take up water excess from current isoil layer
            wtr_soil_sl(:,isoil-1)  = MIN(wtr_soil_sl(:,isoil-1) + hlp1(:),wtr_soil_fc_sl(:,isoil-1) - ice_soil_sl(:,isoil-1) * &
              &                           water2ice_density_ratio)
            drainage_sl(:,isoil-1)  = drainage_sl(:,isoil-1) - hlp1(:) * 1000._wp / dtime
          END WHERE
        ! Aggregated water content of all layers beneath isn't at field capacity
        ! There is capacity in the lower layers to take up excess water  from layer isoil
        ELSEWHERE (hlp2(:) >= hlp1(:))
          ! Excess water is passed on through drainage
          drainage_sl(:,isoil)  = drainage_sl(:,isoil) + hlp1(:) * 1000._wp / dtime
        ELSEWHERE (hlp2(:) < hlp1(:))
          ! Excess water from current layer exceeds aggregated available capacity in lower layers
          ! -> drainage of water up to maximum capacity of aggregated lower layers, then upward inundation up to maximum capacity of layer isoil-1,
          !    then lateral groundwater runoff
          drainage_sl(:,isoil)   = drainage_sl(:,isoil) + hlp2(:) * 1000._wp / dtime
          hlp3(:)                = MIN(hlp1(:) - hlp2(:), wtr_soil_fc_sl(:,isoil-1) - (wtr_soil_sl(:,isoil-1) + &
            &                          ice_soil_sl(:,isoil-1) * water2ice_density_ratio))
          wtr_soil_sl(:,isoil-1) = MAX(MIN(wtr_soil_sl(:,isoil-1) + hlp3, &
                                           wtr_soil_fc_sl(:,isoil-1) - ice_soil_sl(:,isoil-1) * water2ice_density_ratio), 0.0_wp)
          drainage_sl(:,isoil-1) = drainage_sl(:,isoil) - hlp3(:) * 1000._wp / dtime
          gw_runoff_sl(:,isoil)  = gw_runoff_sl(:,isoil) + (hlp1(:) - hlp2(:) - hlp3(:)) * 1000._wp / dtime
        END WHERE
        ! water content in layer isoil should now be down to field capacity
        wtr_soil_sl(:,isoil)      = wtr_soil_fc_sl(:,isoil) - ice_soil_sl(:,isoil) * water2ice_density_ratio
      END WHERE
    ENDDO

    !> 2.3.2 Update deepest soil layer
    ! Compute drainage for deepest layer
    ! Boundary condition for wetland to due frozen bottom soil layers
    WHERE (t_soil_sl(:,nsoil_w) <= Tzero)
      drainage_sl(:,nsoil_w)  = 0.0_wp
    ! Otherwise regular drainage takes place
    ELSEWHERE
      drainage_sl(:,nsoil_w)  = MAX(drain_frac ** soil_depth_sl(:,nsoil_w) * drainage_sl(:,nsoil_w-1),0.0_wp)
    END WHERE
    ! Update water content of deepest layer
    wtr_soil_sl(:,nsoil_w)    = wtr_soil_sl(:,nsoil_w) + &
                                (drainage_sl(:,nsoil_w-1) - drainage_sl(:,nsoil_w) - &
                                evaporation(:) * fevaporation(:,nsoil_w) - &
                                transpiration(:) * ftranspiration_sl(:,nsoil_w)) /  &
                                1000._wp * dtime

    ! conditions in the case that water content in lowest layer is above field capcity
    WHERE (wtr_soil_sl(:,nsoil_w) + ice_soil_sl(:,nsoil_w) * water2ice_density_ratio > wtr_soil_fc_sl(:,nsoil_w))
      hlp1(:) = (wtr_soil_sl(:,nsoil_w) + ice_soil_sl(:,nsoil_w) * water2ice_density_ratio) - wtr_soil_fc_sl(:,nsoil_w)
      hlp2(:) = wtr_soil_fc_sl(:,nsoil_w-1) - (wtr_soil_sl(:,nsoil_w-1) + ice_soil_sl(:,nsoil_w-1) * water2ice_density_ratio)
      WHERE (t_soil_sl(:,nsoil_w) <= Tzero)
        WHERE (hlp1(:) > hlp2(:))
          wtr_soil_sl(:,nsoil_w-1) = wtr_soil_sl(:,nsoil_w-1) + hlp2(:)
          drainage_sl(:,nsoil_w-1) = drainage_sl(:,nsoil_w-1) - hlp2(:) * 1000._wp / dtime
        ELSEWHERE
          wtr_soil_sl(:,nsoil_w-1) = wtr_soil_sl(:,nsoil_w-1) + hlp1(:)
          drainage_sl(:,nsoil_w-1) = drainage_sl(:,nsoil_w-1) - hlp1(:) * 1000._wp / dtime
        END WHERE
      ELSEWHERE
        drainage_sl(:,nsoil_w) = drainage_sl(:,nsoil_w) + hlp1(:) * 1000._wp / dtime
      END WHERE
      wtr_soil_sl(:,nsoil_w) = wtr_soil_fc_sl(:,nsoil_w) - ice_soil_sl(:,nsoil_w) * water2ice_density_ratio
    END WHERE

    !> 2.4. Diffusion fluxes
    !!
    !! lowest value of wtr_soil_sl given minimum number allowed as outcome of theta**C
    hlp1_2d(:,:) = EXP(LOG(eps8)/saxtonC(:,:)) * wtr_soil_sat_sl(:,:)
    WHERE(wtr_soil_sat_sl(:,:) > eps8 .AND. (wtr_soil_sl(:,:)) > hlp1_2d(:,:))
       k_diff(:,:) = kdiff_sat_sl(:,:) * ((wtr_soil_sl(:,:)) / wtr_soil_sat_sl(:,:)) ** saxtonC(:,:)
    ELSEWHERE
       k_diff(:,:) = 0.0_wp
    END WHERE
    !! update soil water potential for diffusion calculation
    WHERE(soil_depth_sl(:,:) > eps8)
       ! Campbell applies only water content below field capacity (water content from ice also considered here)
       hlp1_2d(:,:) = MIN(wtr_soil_fc_sl(:,:), wtr_soil_sl(:,:) + ice_soil_sl(:,:) * water2ice_density_ratio) / soil_depth_sl(:,:)
       ! for numerical reasons, constrain to minimum water potential
       hlp1_2d(:,:) = MAX(EXP(LOG(wpot_min/saxtonA(:,:)) /saxtonB(:,:)), hlp1_2d(:,:))
    ELSEWHERE
       hlp1_2d(:,:) = EXP(LOG(wpot_min/saxtonA(:,:)) /saxtonB(:,:))
    ENDWHERE
    wtr_soil_pot_sl(:,:) = saxtonA(:,:) * (hlp1_2d(:,:) ** saxtonB(:,:))
    ! compute diffusive fluxes and limit them by field capacity
    DO ic = 1,nc
      DO is = 2,INT(num_sl_above_bedrock(ic))
        ! weighted mean of diffusivity (m/s)
        hlp1(ic) = (k_diff(ic,is-1) * soil_depth_sl(ic,is-1) + k_diff(ic,is) * soil_depth_sl(ic,is)) / &
                  (soil_depth_sl(ic,is-1) + soil_depth_sl(ic,is))
        ! flux resulting from the pressure gradient between the centres of the soil layer
        ! note that matrix potential is converted to pressure head MPa -> J/N = m
        hlp1(ic) = hlp1(ic) * dtime * &
                  (wtr_soil_pot_sl(ic,is-1) - wtr_soil_pot_sl(ic,is)) * &
                  1000.0_wp / grav / &
                  (0.5_wp * (soil_depth_sl(ic,is-1) + soil_depth_sl(ic,is)))

        ! limit diffusion rate to the rate that produces equal fractional soil moisture in the two layers
        hlp3(ic) =  soil_depth_sl(ic,is-1) / (soil_depth_sl(ic,is))
        hlp2(ic) = ((wtr_soil_sl(ic,is-1) + ice_soil_sl(ic,is-1) * water2ice_density_ratio) - &
          &        (wtr_soil_sl(ic,is) + ice_soil_sl(ic,is)* water2ice_density_ratio) * hlp3(ic)) / (1._wp + hlp3(ic))
        ! ...
        IF ((hlp1(ic) < 0.0_wp  .AND. hlp1(ic) < hlp2(ic)) .OR. (hlp1(ic) >= 0.0_wp .AND. hlp1(ic) > hlp2(ic))) THEN
          hlp1(ic) = hlp2(ic)
        ENDIF
        ! avoid soil water content exceeding soil water capacity
        IF (hlp1(ic) < 0.0_wp) THEN ! a negative hlp1 means an upward diffusion flux
          ! Limit upward diffusion to available pore space in top layer in MIN() function
          ! MAX() function: hlp1 needs to be large enough to reduce water content of bottom layer to wtr_soil_fc_sl
          hlp1(ic) = -MAX( &
            &       MIN(wtr_soil_fc_sl(ic,is-1) - (wtr_soil_sl(ic,is-1) + ice_soil_sl(ic,is-1) &
            &       * water2ice_density_ratio), -hlp1(ic)), &
            &       wtr_soil_sl(ic,is) + ice_soil_sl(ic,is) * water2ice_density_ratio - wtr_soil_fc_sl(ic,is))
        ELSE ! downward diffusion
          ! Limit downward diffusion to available pore space in bottom layer in MIN() function
          ! MAX() function: hlp1 needs to be large enough to reduce water content of top layer to wtr_soil_fc_sl
          hlp1(ic) = MAX(MIN(hlp1(ic),wtr_soil_fc_sl(ic,is) &
            &        - (wtr_soil_sl(ic,is) + ice_soil_sl(ic,is) * water2ice_density_ratio)), &
            &        wtr_soil_sl(ic,is-1) + ice_soil_sl(ic,is-1) * water2ice_density_ratio - wtr_soil_fc_sl(ic,is-1))
        ENDIF
        ! update soil layers with diffusion
        wtr_soil_sl(ic,is-1)   = wtr_soil_sl(ic,is-1) - hlp1(ic)
        wtr_soil_sl(ic,is)     = wtr_soil_sl(ic,is) + hlp1(ic)
        drainage_sl(ic,is-1) = drainage_sl(ic,is-1) + hlp1(ic) * 1000._wp / dtime
      ENDDO
    ENDDO

    !> 2.5. diagnose deep drainage and the total lateral runoff
    !>
    !>   The threshold 0.75 is to prevent the entire layer from being washed out,
    !>   mainly for the top layers but applied in the whole soil profile, just in case
    !>
    drainage(:) = drainage_sl(:,nsoil_w)
    gw_runoff(:) = SUM(gw_runoff_sl(:,1:nsoil_w), DIM=2)
    ! Layer one: water balance of surface layer
    WHERE(wtr_soil_sl(:,1) > eps8)
      frac_wtr_transp_down_sl(:,1) = MIN( MAX((drainage_sl(:,1) / 1000._wp * dtime) / &
        &                                   (w_soil_sl_old(:,1)), -0.75_wp), 0.75_wp) * 1000._wp / dtime
    ELSEWHERE
      frac_wtr_transp_down_sl(:,1) = 0.0_wp
    ENDWHERE

    ! layer 1 to nsoil_w
    !   using the soil water content of the previous timestep here (w_soil_sl_old) because
    !   the frac_w_lat_loss_sl is applied to the material (elements / nutrients) concentration
    !   of the previous timestep
    DO isoil = 1,nsoil_w
      WHERE(wtr_soil_sl(:,isoil) > eps8)
        frac_wtr_transp_down_sl(:,isoil) = MIN( MAX((drainage_sl(:,isoil) / 1000._wp * dtime) &
          &                                     / (w_soil_sl_old(:,isoil)), frac_wtr_vertical_transport_up_max), &
          &                                frac_wtr_vertical_transport_down_max) &
          &                                * 1000._wp / dtime
        frac_w_lat_loss_sl(:,isoil) = MIN( MAX((gw_runoff_sl(:,isoil) / 1000._wp * dtime) / &
          &                           (w_soil_sl_old(:,isoil)), 0.0_wp), frac_w_lat_loss_max) * 1000._wp / dtime
      ELSEWHERE
        frac_wtr_transp_down_sl(:,isoil) = 0.0_wp
        frac_w_lat_loss_sl(:,isoil)      = 0.0_wp
      ENDWHERE
    ENDDO


    !> 3.0 final diagnostics
    !!

    !> 3.1 calculate rooting zone water content and theta for compatability
    !!       vegetation water stress calculation.
    !!       If a soil layer is only partially filled with roots, the plant-accessible
    !!       water is limited to the AWC of the root zone. All of the layer's water
    !!       (fc-w) deficit is substracted from the root accessible AWC (fc-pwp)*root_depth
    !!
    !! note: upper bound is the larger value compared to lower bound !
    !!
    wtr_rootzone(:) = wtr_soil_sl(:,1)
    DO isoil = 2,nsoil_w
       WHERE(root_depth(:) > soil_layer_w_ubounds(:,isoil-1))
          WHERE(root_depth(:) >= soil_layer_w_ubounds(:,isoil))
            wtr_rootzone(:) = wtr_rootzone(:) + wtr_soil_sl(:,isoil)
          ELSEWHERE
            ! fraction of the layer that has roots
            hlp1(:) = (root_depth(:) - soil_layer_w_ubounds(:,isoil-1)) / &
                      (soil_layer_w_ubounds(:,isoil) - soil_layer_w_ubounds(:,isoil-1))

            wtr_rootzone(:) = wtr_rootzone(:) + wtr_soil_fc_sl(:,isoil) * hlp1(:) - &
                              MIN((wtr_soil_fc_sl(:,isoil) - wtr_soil_sl(:,isoil)), wtr_soil_pwp_sl(:,isoil) * hlp1(:))
          ENDWHERE
       ENDWHERE
    ENDDO

    !> 3.2 calculate plant fractional water availability
    !!
    wtr_rootzone_rel(:) = wtr_rootzone(:) / w_soil_root_fc(:)
    wtr_rootzone_rel(:) = MAX(wtr_rootzone_rel(:),0.0_wp)
    wtr_plant_avail_rel(:) = ( wtr_rootzone(:) - w_soil_root_pwp(:) ) / &
                             ( w_soil_root_fc(:) - w_soil_root_pwp(:) )
    wtr_plant_avail_rel(:) = MAX(wtr_plant_avail_rel(:),0.0_wp)

    !> 3.3 calculate root zone and layer soil water potential according to Campbell et al 1974
    !!
    hlp1(:) = MIN(w_soil_root_fc(:),wtr_rootzone(:)) / root_depth(:)
    ! for numerical reasons, constrain to minimum water potential
    hlp1(:) = MAX(exp(log(wpot_min/saxtonA(:,1))/saxtonB(:,1)),hlp1(:))
    wtr_soil_root_pot(:) = saxtonA(:,1) * (hlp1(:) ** saxtonB(:,1))

    WHERE(soil_depth_sl(:,:) > eps8)
       ! Campbell applies only water content below field capacity
       hlp1_2d(:,:) = MIN(wtr_soil_fc_sl(:,:),wtr_soil_sl(:,:)) / soil_depth_sl(:,:)
       ! for numerical reasons, constrain to minimum water potential
       hlp1_2d(:,:) = MAX(exp(log(wpot_min/saxtonA(:,:))/saxtonB(:,:)),hlp1_2d(:,:))
    ELSEWHERE
       hlp1_2d(:,:) = exp(log(wpot_min/saxtonA(:,:))/saxtonB(:,:))
    ENDWHERE
    wtr_soil_pot_sl(:,:) = saxtonA(:,:) * (hlp1_2d(:,:) ** saxtonB(:,:))

    !> 3.4 Diagnose "soil root water potential"
    !>
    !>   i.e., the soil water potential that the plants detect for growth and scenescence
    !>
    !>   this is calculated assuming plants are more sensitive to the soil moisture
    !>   where they have more root mass
    !>
    !>   the soil root water potential is thus weighted by the fraction of roots per layer
    !>
    wtr_soil_root_pot(:) = SUM(wtr_soil_pot_sl(:,:) * root_fraction_sl(:,:), DIM=2)

    !> 3.5 Diagnostic fraction of water mass lost from the soil colum per second
    !!
    WHERE(w_soil_root_old(:) > eps8)
       drainage_fraction(:) = drainage(:) / 1000._wp / w_soil_root_old(:)
    ELSEWHERE
       drainage_fraction(:) = 0.0_wp
    ENDWHERE

  END SUBROUTINE calc_soil_hydrology


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_surface_energy_balance
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate snow and wet fraction of a tile
  !!  mimicks the JSBACH4 routine calc_wskin_fractions_veg in
  !!   hydroloogy/mo_hyrology_interface:calc_wskin_fractions_veg
  !!
  !! Calculates scaling factors fract_snow_soil and fract_wet
  !-----------------------------------------------------------------------------------------------------
 ELEMENTAL SUBROUTINE calc_wskin_fractions_veg( &
    & dtime,                                     & ! in
    & wsr_max,                                   & ! in
    & oro_stddev,                                & ! in
    & evapo_pot,                                 & ! in
    & wtr_skin,                                   & ! in
    & w_snow_soil,                               & ! in
    & w_snow_can,                                & ! in
    & fract_snow_can,                            & ! out
    & fract_wet    ,                             & ! out
    & fract_snow_soil                            & ! out
    & )

    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_spq_constants,          ONLY: w_density, w_snow_min, snow_dens, w_snow_max


    REAL(wp), INTENT(in) :: &
      & dtime
    REAL(wp), INTENT(in) :: &
      & oro_stddev,         &
      & evapo_pot,          &
      & wsr_max,            &
      & wtr_skin,           &
      & w_snow_soil,        &
      & w_snow_can
    REAL(wp), INTENT(out) :: &
      & fract_snow_can,      &
      & fract_wet    ,     &
      & fract_snow_soil


    IF (w_snow_soil > w_snow_min) THEN
      fract_snow_soil = (w_snow_soil - w_snow_min) / (0.05_wp * snow_dens / w_density - w_snow_min)
    ELSE
      fract_snow_soil = 0.0_wp
    END IF

    ! Snow fraction on canopy
    ! Wet skin reservoir fraction on soil and canopy (between 0 and 1)
    IF (wsr_max > eps8) THEN
      fract_snow_can = MIN(1._wp, w_snow_can / wsr_max)
      fract_wet      = MIN(1._wp, wtr_skin   / wsr_max)
    ELSE
      fract_snow_can = 0._wp
      fract_wet      = 0._wp
    END IF
    fract_snow_can = 0._wp
!    fract_snow_soil = 0._wp


    ! Snow fraction on soil
    fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, &
                               fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))


    ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy; same for skin reservoir
!    IF (fract_snow_soil > 0._wp) THEN
!      ! @TODO Shouldn't one take rhoice here insteady of rhoh2o?
!      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evapo_pot * dtime &
!        &                                            / (rhoh2o * (w_snow_soil + w_snow_can)) &
!        &                                    )
!    END IF
    IF (fract_wet > eps8 ) THEN
      fract_wet = fract_wet    / MAX(1._wp, (1._wp - fract_snow_soil) * evapo_pot * dtime &
        &                                       / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin)) &
        &                                    )
    END IF

  END SUBROUTINE calc_wskin_fractions_veg

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_surface_energy_balance
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate new surface dry static energy using JSBACHs implicit formulation
  !!  mimicks the routine JSBACH4 surface_temp_implicit
  !!   srf_energy_bal/mo_seb_land:update_surface_energy_land:surface_temp_implicit
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_temp_srf_new_implicit(dtime, alpha, t2s_conv,t_Acoef,t_Bcoef,q_Acoef,q_Bcoef, &
                                                  s_old,qsat_srf_old,dQdT, &
                                                  rad_srf_net,ground_heat_flx,heat_transfer_coefficient, &
                                                  cair,csat,snow_cover_fract,heat_capa_srf) RESULT (s_new)

    USE mo_jsb_physical_constants,    ONLY: LatentHeatVaporization, LatentHeatSublimation, stbo, zemiss_def

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 In
    REAL(wp), INTENT(in) :: dtime, &
                            alpha, &
                            t2s_conv, &
                            t_Acoef, &
                            t_Bcoef, &
                            q_Acoef, &
                            q_Bcoef, &
                            s_old, &
                            qsat_srf_old, &
                            dQdT, &
                            rad_srf_net, &
                            ground_heat_flx, &
                            heat_transfer_coefficient, &
                            cair, &
                            csat, &
                            snow_cover_fract, &
                            heat_capa_srf
     ! out
     REAL(wp)            :: s_new
     ! local
     REAL(wp) :: pdt, zicp, zca, zcs, zcolin, zcohfl, zcoind

     pdt = alpha * dtime

     zicp = 1._wp / t2s_conv

     zca = LatentHeatSublimation * snow_cover_fract + LatentHeatVaporization * (cair - snow_cover_fract)
     zcs = LatentHeatSublimation * snow_cover_fract + LatentHeatVaporization * (csat - snow_cover_fract)

     zcolin = heat_capa_srf * zicp + &
              pdt * (zicp * 4.0_wp * stbo * zemiss_def * ((zicp * s_old) ** 3.0_wp)        &
       &              - heat_transfer_coefficient * ( zca * q_Acoef - zcs) * zicp * dQdT)

     zcohfl = -pdt * heat_transfer_coefficient * (t_Acoef - 1.0_wp)

     zcoind = pdt * (rad_srf_net + heat_transfer_coefficient * t_Bcoef + heat_transfer_coefficient &
       &              * ((zca * q_Acoef - zcs) * qsat_srf_old + zca * q_Bcoef) + ground_heat_flx)

     s_new  = (zcolin * s_old + zcoind) / (zcolin + zcohfl)

  END FUNCTION calc_temp_srf_new_implicit

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_spq_physics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate root depth and raction per soil layer
  !!
  !-----------------------------------------------------------------------------------------------------
    SUBROUTINE calc_snow_dynamics(nc, &
                              nsoil_w, &
                              nsnow,  &
                              dtime,  &
                              t_air, &
                              t_old, &
                              t, &
                              evaporation_snow, &
                              heat_capa_sl, &
                              bulk_dens_sl, &
                              soil_depth_sl, &
                              snow_height, &
                              snow_soil_heat_flux, &
                              snow_srf_heat_flux, &
                              snow_melt_to_soil, &
                              t_soil_sl, &
                              t_snow_snl, &
                              w_snow_snl, &
                              snow_present_snl, &
                              snow_lay_thickness_snl)

    USE mo_jsb_math_constants,        ONLY: eps8
    USE mo_spq_constants,             ONLY: latent_heat_fusion,w_density,w_snow_min, snow_dens, &
                                            snow_heat_capa,snow_therm_cond,w_snow_max
    USE mo_jsb_physical_constants,    ONLY: Tzero

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_w, &
                                                                       nsnow
    REAL(wp),                                         INTENT(in)    :: dtime                  !< timestep length
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: t_air,     &           !< air temperature
                                                                       t_old, &               !< old surface temperature [K]
                                                                       t, &                   !< new surface temperature [K]
                                                                       evaporation_snow       !< evaporation snow    [kg / m2]
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(in)    :: heat_capa_sl, &        !< soil heat capacity  [J/kg]
                                                                       bulk_dens_sl, &        !< soil density        [kg m-3]
                                                                       soil_depth_sl          !< soil layer width    [m]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: snow_height, &         !< snow depth          [m]
                                                                       snow_soil_heat_flux, & !< snow soil heat flux [W/m2]
                                                                       snow_srf_heat_flux, &  !< Heat flux from top snow layer
                                                                                              !  to lower snow layer/ soil [W/m2]
                                                                       snow_melt_to_soil      !< snow soil heat flux [m/s]
    REAL(wp), DIMENSION(nc,nsoil_w),                  INTENT(inout) :: t_soil_sl              !< soil temperature     [K]
    REAL(wp), DIMENSION(nc,nsnow),                    INTENT(inout) :: t_snow_snl, &          !< snow layer temperature [K]
                                                                       w_snow_snl, &          !< snow water [m]
                                                                       snow_present_snl, &    !< flag if there is snow [-]
                                                                       snow_lay_thickness_snl !< snow thickness at every layer
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                       ::     w_snow_total, &                      !< total snow water in liquid equivalent[m]
                                                         hlp1, &
                                                         hlp2
    REAL(wp), DIMENSION(nc,nsnow)                 ::     snow_heat_flux, &                !< Snow layer heat flux                             [J/m2/s]
                                                         snow_melt_snl, &                 !< Snow melt in water equivalent                    [m]
                                                         hlp1_2d, &
                                                         hlp2_2d
    INTEGER                                       :: nsnow_top
    INTEGER                                       :: ichunk, isnow   ! looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_snow_dynamics'

    hlp1_2d = 0.0_wp
    hlp2_2d = 0.0_wp

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! -----------------------------------------------------------------------------------------------------------

    !> 0.8 init local var
    !>
    hlp1                = 0.0_wp
    hlp2                = 0.0_wp
    hlp1_2d             = 0.0_wp
    hlp2_2d             = 0.0_wp
    snow_heat_flux      = 0.0_wp
    snow_soil_heat_flux = 0.0_wp
    snow_melt_to_soil   = 0.0_wp

    !> 0.9
    !> make sure snow-layer temperature is zero
    !>   if snow is assumed to be not present in the layer (i.e., less than minimum threshold w_snow_min)
    !>
    WHERE (w_snow_snl(:,:) < w_snow_min)
      t_snow_snl(:,:) = Tzero
    ENDWHERE

    DO ichunk=1,nc
      ! 1.0 Set snow surface temperature as calculated in calc_surface_energy()
      nsnow_top = MAX(NINT(SUM(snow_present_snl(ichunk,:))), 1)
      IF (w_snow_snl(ichunk,nsnow_top) > w_snow_min) THEN
        t_snow_snl(ichunk,nsnow_top) = t(ichunk)           ! temperature of top layer equal that of air
      ELSE ! to melt the rest of the snow even when energy exchanges are assumed to be too small for vertical conductivity/implicit scheme!
        t_snow_snl(ichunk,nsnow_top) = t_air(ichunk)
      END IF
      ! 2.0 Calculate heat transfer through snow layers
      IF (w_snow_snl(ichunk,1) > w_snow_min) THEN
        DO isnow = nsnow_top, 2, -1  ! Loop from the top soil layer
          IF (w_snow_snl(ichunk,isnow) >= w_snow_min) THEN
            snow_heat_flux(ichunk,isnow) = snow_therm_cond * (t_snow_snl(ichunk,isnow) - t_snow_snl(ichunk,isnow - 1)) / &
              &                            ((snow_lay_thickness_snl(ichunk,isnow) / &
              &                            2._wp  + snow_lay_thickness_snl(ichunk,isnow - 1) / 2._wp))
          END IF
        END DO
        ! 2.1 Update snow temperatures
        !
        DO isnow = nsnow_top-1, 2, -1
          t_snow_snl(ichunk,isnow) = t_snow_snl(ichunk,isnow) + (snow_heat_flux(ichunk,isnow + 1) - &
                                     snow_heat_flux(ichunk,isnow)) * dtime / &
                                     (snow_heat_capa * snow_dens *  snow_lay_thickness_snl(ichunk,isnow))
        END DO

        ! 3.0 Calculate heat transfer between soil surface layer and first snow layer
        IF (w_snow_snl(ichunk,1) > w_snow_min) THEN
          snow_soil_heat_flux(ichunk) = snow_therm_cond * (t_snow_snl(ichunk,1) - t_soil_sl(ichunk,1)) / &
            &                           (snow_lay_thickness_snl(ichunk,1))
          IF ((snow_heat_flux(ichunk,2) - snow_soil_heat_flux(ichunk)) /= 0.0_wp) THEN
            t_snow_snl(ichunk,1) = t_snow_snl(ichunk,1) + (snow_heat_flux(ichunk,2) - snow_soil_heat_flux(ichunk)) * dtime / &
                                   (snow_heat_capa * snow_dens * snow_lay_thickness_snl(ichunk,1))
            t_soil_sl(ichunk,1)  = t_soil_sl(ichunk,1)  + snow_soil_heat_flux(ichunk) * dtime / &
                                  (heat_capa_sl(ichunk,1) * bulk_dens_sl(ichunk,1) * soil_depth_sl(ichunk,1))
          END IF
        END IF
      END IF
      ! 4.0 Snow latent heat
      ! 4.1 Consider snow evaporation that is calculated in calc_surface_energy()
      hlp1(ichunk) = evaporation_snow(ichunk) * dtime / 1000._wp
      DO isnow = nsnow_top, 1, -1
        IF (hlp1(ichunk) > 0.0_wp) THEN
          hlp2(ichunk) = MIN(hlp1(ichunk), w_snow_snl(ichunk,isnow))
          w_snow_snl(ichunk,isnow) = w_snow_snl(ichunk,isnow) - hlp2(ichunk)
          hlp1(ichunk) = hlp1(ichunk) - hlp2(ichunk)
        END IF
      END DO
      ! 4.2 Snow melting
      snow_melt_snl(ichunk,:) = 0.0_wp
      ! snow melting
      WHERE(w_snow_snl(ichunk,:) > 0.0_wp)
        WHERE(t_snow_snl(ichunk,:) > Tzero)
        ! Energy excess
          hlp1_2d(ichunk,:) = snow_heat_capa * snow_dens * snow_lay_thickness_snl(ichunk,:) * (t_snow_snl(ichunk,:) - Tzero)
          hlp2_2d(ichunk,:) = latent_heat_fusion * w_snow_snl(ichunk,:) * w_density
          WHERE (hlp1_2d(ichunk,:) <= hlp2_2d(ichunk,:))
            snow_melt_snl(ichunk,:)            = hlp1_2d(ichunk,:) / (latent_heat_fusion * w_density)
            w_snow_snl(ichunk,:)               = w_snow_snl(ichunk,:) - snow_melt_snl(ichunk,:)
            t_snow_snl(ichunk,:)               = Tzero
          ELSEWHERE
            snow_melt_snl(ichunk,:)            = w_snow_snl(ichunk,:)
            w_snow_snl(ichunk,:)               = 0.0_wp
            t_snow_snl(ichunk,:)               = Tzero
          END WHERE
        END WHERE
      END WHERE
      snow_melt_to_soil(ichunk) = SUM(snow_melt_snl(ichunk,:) * 1000.0_wp / dtime)
    END DO

    ! 5.0 Redistribute snow in the case lower layers have melted
    w_snow_total(:) = SUM(w_snow_snl(:,:),DIM=2)
    hlp1(:)         = w_snow_total(:)
    w_snow_snl(:,:) = 0.0_wp
    DO isnow=1,nsnow-1
      WHERE (hlp1(:) > 0.0_wp)
        hlp2(:)             = MIN(w_snow_max,hlp1(:))
        w_snow_snl(:,isnow) = hlp2(:)
        hlp1(:)             = hlp1(:) - hlp2(:)
      END WHERE
    END DO
    ! top snow layer isn't limited in height
    WHERE (hlp1(:) > 0.0_wp)
      w_snow_snl(:,nsnow) = w_snow_snl(:,nsnow) + hlp1(:)
    END WHERE

    ! 6.0 Calculate snow depth again based on changes in snow water, this snow depth should be given to output
    snow_height = 0.0_wp
    snow_lay_thickness_snl = 0.0_wp
    DO isnow=1,nsnow
      WHERE (w_snow_snl(:,isnow) > 0.0_wp .AND. w_snow_snl(:,1) > w_snow_min)
        snow_present_snl(:,isnow) = 1.0_wp
        snow_lay_thickness_snl(:,isnow) = w_snow_snl(:,isnow) * w_density / snow_dens
        snow_height(:)                  = snow_height(:) + snow_lay_thickness_snl(:,isnow)
        WHERE (w_snow_snl(:,isnow) > w_snow_min)
          snow_present_snl(:,isnow) = 1.0_wp
        END WHERE
      ELSEWHERE
        snow_present_snl(:,isnow) = 0.0_wp
      END WHERE
    END DO
    ! compute snow surface (top snow layer) heat flx for implicit scheme
    DO ichunk=1,nc
      nsnow_top = MAX(NINT(SUM(snow_present_snl(ichunk,:))), 1)
      snow_srf_heat_flux(ichunk) = snow_heat_flux(ichunk,nsnow_top)
    ENDDO
  END SUBROUTINE

#endif
END MODULE mo_q_spq_update_physics
