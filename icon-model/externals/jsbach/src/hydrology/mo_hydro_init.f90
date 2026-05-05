!> Initialization of the hydrology process
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
!>#### Initialize soil properties and initial conditions of the hydrology process
!>
!> The module provides the subroutines used to initialize the HYDRO process.
!> Main subroutine is [[hydro_init]], which is calling the other initialization routines.
!>
!> Depending on the HYDRO configuration, i.e. the namelist parameters defined for the
!> hydrology in namelist 'jsb_hydro_nml', initial and boundary condition data is read
!> from specified input data files in subroutine [[hydro_read_init_vars]]. Reading
!> the data happens only once, by the top tile. Based on this data boundary conditions
!> such as the soil properties are set for all tiles in subroutine [[hydro_init_bc]]
!> while initial conditions as e.g. the initial soil moisture are set in subroutine
!> [[hydro_init_ic]]. After the last tile has been initialized, the initialization
!> prozess is finalized by subroutine [[hydro_finalize_init_vars]], where the memory
!> that was only required in the initialization phase is released again.
!>
MODULE mo_hydro_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, warning, finish
  USE mo_jsb_control,        ONLY: debug_on, jsbach_runs_standalone

  USE mo_jsb_model_class,    ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, soil_init_from_texture

  dsl4jsb_Use_processes HYDRO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)

#ifdef __QUINCY_STANDALONE__
  dsl4jsb_Use_processes SSE_
  dsl4jsb_Use_config(SSE_)
#else
  dsl4jsb_Use_processes PHENO_
  dsl4jsb_Use_memory(PHENO_)
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hydro_init, hydro_sanitize_state

  ! -------------------------------------------------------------------------------------------------- !
  ! Lookup table for deriving soil properties from soil textures

  ! TODO move the below parameters to constants / parameter module

  ! Parameters for sand, loam (here called silt), clay, peat (here called oc) are used from the TERRA
  ! land model in ICON (see src/lnd_phy_schemes/sfc_terra_data.f90 in ICON)

  ! Pore volume (called cporv in TERRA)
  REAL(wp), PARAMETER ::  porosity_sand = 0.364_wp      !< pore volume of sand [vol. fraction]
  REAL(wp), PARAMETER ::  porosity_silt = 0.455_wp      !< pore volume of loam [vol. fraction]
  REAL(wp), PARAMETER ::  porosity_clay = 0.507_wp      !< pore volume of clay [vol. fraction]
  REAL(wp), PARAMETER ::  porosity_oc   = 0.863_wp      !< pore volume of peat [vol. fraction]

  ! Field capacity (called cfcap in TERRA)
  REAL(wp), PARAMETER ::  field_cap_sand = 0.196_wp     !< field capacity of sand [vol. fraction]
  REAL(wp), PARAMETER ::  field_cap_silt = 0.340_wp     !< field capacity of loam [vol. fraction]
  REAL(wp), PARAMETER ::  field_cap_clay = 0.463_wp     !< field capacity of clay [vol. fraction]
  REAL(wp), PARAMETER ::  field_cap_oc   = 0.763_wp     !< field capacity of peat [vol. fraction]

  ! Hydraulic conductivity at saturation (called ckw0 in TERRA)
  REAL(wp), PARAMETER ::  hyd_cond_sand = 4.79e-5_wp    !< hydraulic conductivity at saturation of sand [m/s]
  REAL(wp), PARAMETER ::  hyd_cond_silt = 5.31e-6_wp    !< hydraulic conductivity at saturation of loam [m/s]
  REAL(wp), PARAMETER ::  hyd_cond_clay = 8.50e-8_wp    !< hydraulic conductivity at saturation of clay [m/s]
  REAL(wp), PARAMETER ::  hyd_cond_oc   = 5.80e-8_wp    !< hydraulic conductivity at saturation of peat [m/s]

  ! Wilting point (called cpwp in TERRA)
  REAL(wp), PARAMETER ::  wilt_sand = 0.042_wp          !< wilting point of sand [vol. fraction]
  REAL(wp), PARAMETER ::  wilt_silt = 0.11_wp           !< wilting point of loam [vol. fraction]
  REAL(wp), PARAMETER ::  wilt_clay = 0.257_wp          !< wilting point of clay [vol. fraction]
  REAL(wp), PARAMETER ::  wilt_oc   = 0.265_wp          !< wilting point of peat [vol. fraction]

  ! Pore size index
  REAL(wp), PARAMETER :: pore_size_index_sand = 0.35_wp !< pore size index of sand []
  REAL(wp), PARAMETER :: pore_size_index_silt = 0.2_wp  !< pore size index of loam []
  REAL(wp), PARAMETER :: pore_size_index_clay = 0.13_wp !< pore size index of clay []
  REAL(wp), PARAMETER :: pore_size_index_oc   = 0.65_wp !< pore size index of loam []

  ! Clapp & Hornberger exponent b (see Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_sand = 3.39_wp   !< Clapp & Hornberger exponent b for sand (type 1 in Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_silt = 4.98_wp   !< Clapp & Hornberger exponent b for loam (type 5 in Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_clay = 10.38_wp  !< Clapp & Hornberger exponent b for clay (type 10 in Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_oc   = 4._wp     !< Clapp & Hornberger exponent b for peat (type 12 in Beringer et al. 2001)

  ! Soil matric potential (see Beringer et al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_sand = -0.04729_wp !< soil matric potential of sand [m] (type 1 in Beringer at al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_silt = -0.45425_wp !< soil matric potential of loam [m] (type 5 in Beringer at al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_clay = -0.633_wp   !< soil matric potential of clay [m] (type 10 in Beringer at al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_oc = -0.12_wp      !< soil matric potential of peat [m] (type 12 in Beringer at al. 2001)

  REAL(wp), PARAMETER ::  hyd_cond_sat_profile = 0.432332_wp
                          !< factor to compensate profile of hyd_cond_sat with depth in TERRA

  !> Residual soil moisture (fraction of volume; Maidment, Handbook of Hydrology, 1993)
  REAL(wp), PARAMETER ::  wres_sand = 0.020_wp          !< residual soil moisture (type sand) [vol. fraction]
  REAL(wp), PARAMETER ::  wres_silt = 0.015_wp          !< residual soil moisture (type silt loam) [vol. fraction]
  REAL(wp), PARAMETER ::  wres_clay = 0.090_wp          !< residual soil moisture (type clay) [vol. fraction]
  REAL(wp), PARAMETER ::  wres_oc   = 0.150_wp          !< residual soil moisture [vol. fraction] (Letts et al., 2000)

  !
  ! -------------------------------------------------------------------------------------------------- !
  !
  !> Type to hold variables read from input files
  TYPE t_hydro_init_vars
    REAL(wp), POINTER ::                  &
      & elevation        (:,:) => NULL(), & !< geometric height of surface above sea level [m]
      & oro_stddev       (:,:) => NULL(), & !< standard deviation of subgrid-scale orography [m]
      & soil_depth       (:,:) => NULL(), & !< soil depth until bedrock [m]
      & vol_porosity     (:,:) => NULL(), & !< volumetric soil porosity []
      & hyd_cond_sat     (:,:) => NULL(), & !< saturated hydraulic conductivity [m s-1]
      & matric_pot       (:,:) => NULL(), & !< saturated matric potential [m]
      & bclapp           (:,:) => NULL(), & !< Clapp and Hornberger exponent b []
      & pore_size_index  (:,:) => NULL(), & !< pore size distribution index []
      & vol_field_cap    (:,:) => NULL(), & !< volumetric soil field capacity []
      & vol_p_wilt       (:,:) => NULL(), & !< volumetric wilting point []
      & vol_wres         (:,:) => NULL(), & !< volumetric residual water content []
      & root_depth       (:,:) => NULL(), & !< root depth [m]
      & weq_snow_soil    (:,:) => NULL(), & !< initial snow on surface [m water equivalent]
      & fract_pond_max   (:,:) => NULL(), & !< maximum pond fraction []
      & depth_pond_max   (:,:) => NULL(), & !< maximum pond depth [m]
      & wtr_soil_sl    (:,:,:) => NULL(), & !< initial soil moisture within soil layers [m water equivalent]
      & fract_org_sl   (:,:,:) => NULL(), & !< fraction of soil organic matter within soil layers []
      & fr_sand          (:,:) => NULL(), & !< fraction of sand for soil index []
      & fr_silt          (:,:) => NULL(), & !< fraction of silt for soil index []
      & fr_clay          (:,:) => NULL(), & !< fraction of clay for soil index []
      & fr_sand_deep     (:,:) => NULL(), & !< fraction of sand for deep soil index []
      & fr_silt_deep     (:,:) => NULL(), & !< fraction of silt for deep soil index []
      & fr_clay_deep     (:,:) => NULL(), & !< fraction of clay for deep soil index []
      & ifs_smi_sl     (:,:,:) => NULL(), & !< volumetric soil water (soil moisture index) from IFS []
      & ifs_weq_snow     (:,:) => NULL(), & !< snow on surface from IFS [m water equivalent]
      & soil_sat_water_content(:,:) => NULL(), &  !< saturated water content (awc)
      & soil_theta_prescribe(:,:)   => NULL()     !< soil theta necessary for calc intial root-zone soil water
  END TYPE t_hydro_init_vars

  TYPE(t_hydro_init_vars) :: hydro_init_vars                !< Module variable holding input fields

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_init'  !< Name of this module

CONTAINS

  ! ====================================================================================================== !
  !>
  !>#### Intialize the hydrology process
  !>
  !> Hydrological initial and boundary condition data is read only once - by the top tile - in subroutine
  !> [[hydro_read_init_vars]] and is made avilable to the other tiles in subroutines [[hydro_init_bc]]
  !> and [[hydro_init_ic]].
  !>
  SUBROUTINE hydro_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile               !< current tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init' !< name of this routine

#ifdef __QUINCY_STANDALONE__
    model  => Get_model(tile%owner_model_id)
    IF (model%config%use_soil_phys_jsbach) THEN
      CALL qs_hydro_read_init_vars(tile)
      CALL qs_hydro_init_ic_bc(tile)
      CALL hydro_finalize_init_vars()
    END IF
#else
    ! Read initial and boundary condition data by the top tile, only.
    ! Note: this requires that the init routines in mo_jsb_model_init:jsbach_init are called starting
    !       from the top tile!
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN  ! This is the top tile.
      CALL hydro_read_init_vars(tile)
    END IF

    ! Set soil properties (on all tiles)
    CALL hydro_init_bc(tile)

    ! Set soil initial conditions (on all tiles)
    CALL hydro_init_ic(tile)

    ! Finalize hydrology initialization after the last tile has been initialized
    IF (tile%Is_last_process_tile(HYDRO_)) THEN
      CALL hydro_finalize_init_vars()
    END IF
#endif
  END SUBROUTINE hydro_init

  ! ====================================================================================================== !
  !>
  !>#### Finalize hydrology initialization
  !>
  SUBROUTINE hydro_finalize_init_vars

    !> Deallocation of variables that had only been needed during the initialization phase

    DEALLOCATE( &
      & hydro_init_vars%oro_stddev      ,       &
      & hydro_init_vars%soil_depth      ,       &
      & hydro_init_vars%vol_porosity    ,       &
      & hydro_init_vars%hyd_cond_sat    ,       &
      & hydro_init_vars%matric_pot      ,       &
      & hydro_init_vars%bclapp          ,       &
      & hydro_init_vars%pore_size_index ,       &
      & hydro_init_vars%vol_field_cap   ,       &
      & hydro_init_vars%vol_p_wilt      ,       &
      & hydro_init_vars%vol_wres        ,       &
      & hydro_init_vars%root_depth      ,       &
      & hydro_init_vars%fract_org_sl    ,       &
      & hydro_init_vars%fract_pond_max  ,       &
      & hydro_init_vars%depth_pond_max  ,       &
      & hydro_init_vars%weq_snow_soil)

    IF (ASSOCIATED(hydro_init_vars%elevation))    DEALLOCATE(hydro_init_vars%elevation)
    IF (ASSOCIATED(hydro_init_vars%fr_sand))      DEALLOCATE(hydro_init_vars%fr_sand)
    IF (ASSOCIATED(hydro_init_vars%fr_silt))      DEALLOCATE(hydro_init_vars%fr_silt)
    IF (ASSOCIATED(hydro_init_vars%fr_clay))      DEALLOCATE(hydro_init_vars%fr_clay)
    IF (ASSOCIATED(hydro_init_vars%fr_sand_deep)) DEALLOCATE(hydro_init_vars%fr_sand_deep)
    IF (ASSOCIATED(hydro_init_vars%fr_silt_deep)) DEALLOCATE(hydro_init_vars%fr_silt_deep)
    IF (ASSOCIATED(hydro_init_vars%fr_clay_deep)) DEALLOCATE(hydro_init_vars%fr_clay_deep)
    IF (ASSOCIATED(hydro_init_vars%wtr_soil_sl))  DEALLOCATE(hydro_init_vars%wtr_soil_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_smi_sl))   DEALLOCATE(hydro_init_vars%ifs_smi_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_weq_snow)) DEALLOCATE(hydro_init_vars%ifs_weq_snow)
    IF (ASSOCIATED(hydro_init_vars%soil_sat_water_content)) DEALLOCATE(hydro_init_vars%soil_sat_water_content)
    IF (ASSOCIATED(hydro_init_vars%soil_theta_prescribe))   DEALLOCATE(hydro_init_vars%soil_theta_prescribe)

  END SUBROUTINE hydro_finalize_init_vars

#ifdef __QUINCY_STANDALONE__
  ! ============================================================================================== !
  !>
  !>#### Read boundary and initial conditions of the hydrology for quincy standalone
  !>
  SUBROUTINE qs_hydro_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w
    INTEGER                    :: nproma, nblks, nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':qs_hydro_read_init_vars'


    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Reading/setting hydrology init vars')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_w => Get_vgrid('soil_depth_water')
    nsoil  = soil_w%n_levels

    ALLOCATE( &
      & hydro_init_vars%oro_stddev      (nproma, nblks),       &
      & hydro_init_vars%soil_depth      (nproma, nblks),       &
      & hydro_init_vars%vol_porosity    (nproma, nblks),       &
      & hydro_init_vars%hyd_cond_sat    (nproma, nblks),       &
      & hydro_init_vars%matric_pot      (nproma, nblks),       &
      & hydro_init_vars%bclapp          (nproma, nblks),       &
      & hydro_init_vars%pore_size_index (nproma, nblks),       &
      & hydro_init_vars%vol_field_cap   (nproma, nblks),       &
      & hydro_init_vars%vol_p_wilt      (nproma, nblks),       &
      & hydro_init_vars%vol_wres        (nproma, nblks),       &
      & hydro_init_vars%root_depth      (nproma, nblks),       &
      & hydro_init_vars%fract_pond_max  (nproma, nblks),       &
      & hydro_init_vars%depth_pond_max  (nproma, nblks),       &
      & hydro_init_vars%weq_snow_soil   (nproma, nblks),       &
      & hydro_init_vars%fr_sand         (nproma, nblks),       &
      & hydro_init_vars%fr_silt         (nproma, nblks),       &
      & hydro_init_vars%fr_clay         (nproma, nblks),       &
      & hydro_init_vars%fr_sand_deep    (nproma, nblks),       &
      & hydro_init_vars%fr_silt_deep    (nproma, nblks),       &
      & hydro_init_vars%fr_clay_deep    (nproma, nblks),       &
      & hydro_init_vars%elevation       (nproma, nblks),       &
      & hydro_init_vars%soil_sat_water_content(nproma, nblks), &
      & hydro_init_vars%soil_theta_prescribe(nproma, nblks) &
      & )

      hydro_init_vars%oro_stddev     (:,:)   = missval
      hydro_init_vars%soil_depth     (:,:)   = missval
      hydro_init_vars%vol_porosity   (:,:)   = missval
      hydro_init_vars%hyd_cond_sat   (:,:)   = missval
      hydro_init_vars%matric_pot     (:,:)   = missval
      hydro_init_vars%bclapp         (:,:)   = missval
      hydro_init_vars%pore_size_index(:,:)   = missval
      hydro_init_vars%vol_field_cap  (:,:)   = missval
      hydro_init_vars%vol_p_wilt     (:,:)   = missval
      hydro_init_vars%vol_wres       (:,:)   = missval
      hydro_init_vars%root_depth     (:,:)   = missval
      hydro_init_vars%fract_pond_max (:,:)   = missval
      hydro_init_vars%depth_pond_max (:,:)   = missval
      hydro_init_vars%weq_snow_soil  (:,:)   = missval
      hydro_init_vars%fr_sand        (:,:)   = missval
      hydro_init_vars%fr_silt        (:,:)   = missval
      hydro_init_vars%fr_clay        (:,:)   = missval
      hydro_init_vars%fr_sand_deep   (:,:)   = missval
      hydro_init_vars%fr_silt_deep   (:,:)   = missval
      hydro_init_vars%fr_clay_deep   (:,:)   = missval
      hydro_init_vars%elevation      (:,:)   = missval
      hydro_init_vars%elevation      (:,:)   = missval
      hydro_init_vars%soil_sat_water_content(:,:) = missval
      hydro_init_vars%soil_theta_prescribe(:,:)   = missval

    ALLOCATE(hydro_init_vars%fract_org_sl(nproma, nsoil, nblks))
    hydro_init_vars%fract_org_sl(:,:,:) = missval
    ALLOCATE(hydro_init_vars%wtr_soil_sl(nproma, nsoil, nblks))
    hydro_init_vars%wtr_soil_sl(:,:,:) = missval

    ! ----------------------------------------------------------------------------------------------------- !

    !> total soil depth until bedrock
    !>
    ! quincy standalone soil_depth for one site (read from namelist)
    hydro_init_vars%soil_depth     (:,:)   = dsl4jsb_Config(HYDRO_)%qs_soil_depth
    ! for all points of the domain soil_depth(:,:) is at least 0.1
    hydro_init_vars%soil_depth     (:,:)   = MAX(0.1_wp, MERGE(hydro_init_vars%soil_depth(:,:), 0._wp, &
      &                                      hydro_init_vars%soil_depth(:,:) >= 0._wp))

    !> sand / silt / clay
    !>
    hydro_init_vars%fr_sand        (:,:)   = dsl4jsb_Config(SSE_)%qs_soil_sand
    hydro_init_vars%fr_silt        (:,:)   = dsl4jsb_Config(SSE_)%qs_soil_silt
    hydro_init_vars%fr_clay        (:,:)   = dsl4jsb_Config(SSE_)%qs_soil_clay
    hydro_init_vars%fr_sand_deep   (:,:)   = dsl4jsb_Config(SSE_)%qs_soil_sand
    hydro_init_vars%fr_silt_deep   (:,:)   = dsl4jsb_Config(SSE_)%qs_soil_silt
    hydro_init_vars%fr_clay_deep   (:,:)   = dsl4jsb_Config(SSE_)%qs_soil_clay

    !> saturated water content (awc)
    !>
    hydro_init_vars%soil_sat_water_content(:,:) = dsl4jsb_Config(HYDRO_)%qs_soil_awc_prescribe

    !> soil_theta_prescribe - needed to calc initial water amount in root zone
    hydro_init_vars%soil_theta_prescribe(:,:) = dsl4jsb_Config(HYDRO_)%qs_soil_theta_prescribe

    !> root depth - is calculated below
    !>
    ! hydro_init_vars%root_depth(:,:) = x

    !> Water content of soil layers [m] - is calculated below
    !>
    ! hydro_init_vars%wtr_soil_sl(:,:,:)  = x

    !> other variables - set useful default values
    !>
    hydro_init_vars%oro_stddev     (:,:)   = 30.0_wp  !< Standard deviation of the orography [m] - estimated from icon-land input file
    hydro_init_vars%weq_snow_soil  (:,:)   = 0.0_wp   !< initial snow on surface [m water equivalent]
    hydro_init_vars%elevation      (:,:)   = 60.0_wp  !< Topographic height [m] - estimated from icon-land input file
    ! hydro_init_vars%vol_porosity   (:,:)   = x   !< is calculated below
    ! hydro_init_vars%hyd_cond_sat   (:,:)   = x   !< is calculated below
    ! hydro_init_vars%matric_pot     (:,:)   = x   !< is calculated below
    ! hydro_init_vars%bclapp         (:,:)   = x   !< is calculated below
    ! hydro_init_vars%pore_size_index(:,:)   = x   !< is calculated below
    ! hydro_init_vars%vol_field_cap  (:,:)   = x   !< is calculated below
    ! hydro_init_vars%vol_p_wilt     (:,:)   = x   !< is calculated below
    ! hydro_init_vars%vol_wres       (:,:)   = x   !< is calculated below

    !> set to arbitrary value - not relevant for quincy standalone
    !>
    hydro_init_vars%fract_pond_max (:,:)   = 0.0_wp   !< maximum pond fraction []
    hydro_init_vars%depth_pond_max (:,:)   = 0.0_wp   !< maximum pond depth [m]

    !> Soil organic carbon fractions for each soil layer []
    !>
    !>   NOTE, needs improvment / proper parameterization
    !>   IQ uses fract_org_sl data from bc_land_soil.nc (rawdata = "GLOBCOVER2009, HWSD, GLOBE, Lake Database")
    !>
    hydro_init_vars%fract_org_sl(:,:,:) = 0.0_wp      ! all soil layers
    hydro_init_vars%fract_org_sl(:,1,:) = 0.1_wp      ! 1st soil layer
    hydro_init_vars%fract_org_sl(:,2,:) = 0.01_wp     ! 2nd soil layer
    hydro_init_vars%fract_org_sl(:,3,:) = 0.01_wp     ! 3rd soil layer

  END SUBROUTINE qs_hydro_read_init_vars

  ! ====================================================================================================== !
  !
  !>#### Set initial and boundary conditions for the hydrology process - quincy standalone
  !>
  SUBROUTINE qs_hydro_init_ic_bc(tile)
    USE mo_util,                  ONLY: soil_depth_to_layers_2d
    USE mo_hydro_util,            ONLY: get_amount_in_rootzone
    USE mo_hydro_process,         ONLY: calc_orographic_features
    USE mo_hydro_constants,       ONLY: &
      &                           k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
      &                           k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
      &                           k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt
    USE mo_jsb_math_constants,    ONLY: eps8
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER              :: model
    TYPE(t_jsb_grid),  POINTER              :: hgrid
    TYPE(t_jsb_vgrid), POINTER              :: soil_w
    INTEGER                                 :: nsoil, nblks, nproma
    INTEGER                                 :: ic, iblk, is
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: soil_awc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: hlp2
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: hlp1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_pwp_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_fc_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_sat_sl
    CHARACTER(len=*), PARAMETER :: routine = modname//':qs_hydro_init_ic_bc'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onDomain    :: soil_depth
    dsl4jsb_Real2D_onDomain    :: num_sl_above_bedrock
    dsl4jsb_Real2D_onDomain    :: root_depth
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pwp
    dsl4jsb_Real2D_onDomain    :: w_soil_root_fc
    dsl4jsb_Real2D_onDomain    :: wtr_skin
    dsl4jsb_Real2D_onDomain    :: wtr_rootzone
    dsl4jsb_Real2D_onDomain    :: wtr_rootzone_rel
    dsl4jsb_Real2D_onDomain    :: wtr_plant_avail_rel
    dsl4jsb_Real2D_onDomain    :: oro_stddev
    dsl4jsb_Real2D_onDomain    :: weq_snow_soil
    dsl4jsb_Real2D_onDomain    :: elevation
    dsl4jsb_Real2D_onDomain    :: fract_pond_max
    dsl4jsb_Real2D_onDomain    :: weq_pond_max
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_depth_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain    :: wtr_soil_sl
    dsl4jsb_Real3D_onDomain    :: wtr_soil_pwp_sl
    dsl4jsb_Real3D_onDomain    :: wtr_soil_fc_sl
    dsl4jsb_Real3D_onDomain    :: wtr_soil_sat_sl
    dsl4jsb_Real3D_onDomain    :: fract_org_sl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN
    IF (debug_on()) CALL message(routine, 'Setting  QUINCY hydrology intial and boundary conditions for tile '//TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    model   => Get_model(tile%owner_model_id)
    hgrid   => Get_grid(model%grid_id)
    soil_w  => Get_vgrid('soil_depth_water')
    nsoil   =  soil_w%n_levels
    nblks   =  hgrid%nblks
    nproma  =  hgrid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(soil_awc(nproma, nblks))
    ALLOCATE(hlp2(nproma, nblks))
    ALLOCATE(hlp1(nproma, nsoil, nblks))
    ALLOCATE(theta_pwp_sl(nproma, nsoil, nblks))
    ALLOCATE(theta_fc_sl(nproma, nsoil, nblks))
    ALLOCATE(theta_sat_sl(nproma, nsoil, nblks))
    ! ----------------------------------------------------------------------------------------------------- !
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   soil_depth)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   num_sl_above_bedrock)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   root_depth)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   w_soil_root_pwp)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   w_soil_root_fc)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   wtr_skin)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   wtr_rootzone)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   wtr_rootzone_rel)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   wtr_plant_avail_rel)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   oro_stddev)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   weq_snow_soil)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   elevation)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   fract_pond_max)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   weq_pond_max)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_width_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_lbound_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_ubound_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_center_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   wtr_soil_pwp_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   wtr_soil_fc_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   wtr_soil_sat_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   fract_org_sl)
    ! ----------------------------------------------------------------------------------------------------- !

    !> total soil depth until bedrock (from namelist)
    !>
    soil_depth(:,:) = hydro_init_vars%soil_depth

    !> calc actual soil layer depths based on the fixed layer thicknesses from namelist and the bedrock depth
    !>
    ! NOTE: soil_depth_to_layers_2d() returns a 3D variable
    soil_depth_sl(:,:,:) = soil_depth_to_layers_2d(soil_depth(:,:), &   ! Total soil depth until bedrock (from textures), 2D variable
      &                                            soil_w%dz(:))        ! Soil layer thicknesses from namelist
    ! pass these values to a HYDRO_ variable with an improved name
    soil_lay_width_sl(:,:,:) = soil_depth_sl(:,:,:)

    !>  calc more metrics of the soil layers
    !>
    ! i)  lower & upper boundary
    ! ii) depth at the center of each layer
    DO iblk = 1,nblks
      DO is = 1,nsoil
        ! lower & upper boundary
        IF (is == 1) THEN
          soil_lay_depth_lbound_sl(:, is, iblk) = 0.0_wp
          soil_lay_depth_ubound_sl(:, is, iblk) = soil_lay_width_sl(:, is, iblk)
        ELSE
          soil_lay_depth_lbound_sl(:, is, iblk) = soil_lay_depth_lbound_sl(:, is-1, iblk) &
            &                                     + soil_lay_width_sl(:, is-1, iblk)
          soil_lay_depth_ubound_sl(:, is, iblk) = soil_lay_depth_ubound_sl(:, is-1, iblk) &
            &                                     + soil_lay_width_sl(:, is, iblk)
        ENDIF
        ! soil-layer center
        soil_lay_depth_center_sl(:, is, iblk) = (soil_lay_depth_lbound_sl(:, is, iblk) &
          &                                     + soil_lay_depth_ubound_sl(:, is, iblk)) &
          &                                     * 0.5_wp
      END DO
    END DO

    !> get number of soil layers above bedrock for each gridcell
    !>
    !>    use this variable in quincy for looping over soil layers,
    !>    excluding soil layers with a width smaller/equal eps8
    !>
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,nsoil
          IF (soil_lay_width_sl(ic, is, iblk) > eps8) THEN
            num_sl_above_bedrock(ic, iblk) = REAL(is, wp)
          ELSE
            EXIT  ! stop looping over further soil layers with thickness < eps8
          END IF
        END DO
      END DO
    END DO

    !> calc depth of roots
    !>
    !>   given by prescribed AWC (soil water in the rooting zone [m])
    !>
    !>   the init of 'root_fraction_sl(:,:,:)' below depends on this init of the root_depth(:)
    !>

    ! calc fractional water holding capacity at permanent wilting point
    ! NOTE that this assumes that the soil hydraulic properties are identical between the layers
    !
    hlp1(:,:,:)         = k_pwp_s    * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_pwp_c  * SPREAD(hydro_init_vars%fr_clay(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_pwp_sc * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) &
      &                   * SPREAD(hydro_init_vars%fr_clay(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_pwp_a
    theta_pwp_sl(:,:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:,:)
    ! calc fractional water holding capacity at field capacity
    !
    hlp1(:,:,:)         = k_fc_s    * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_fc_c  * SPREAD(hydro_init_vars%fr_clay(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_fc_sc * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) &
      &                   * SPREAD(hydro_init_vars%fr_clay(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_fc_a
    theta_fc_sl(:,:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:,:) + k_fc_ct * hlp1(:,:,:) ** 2._wp
    ! calc fractional water holding capacity at saturation
    !
    hlp1(:,:,:)         = k_sat_s    * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_sat_c  * SPREAD(hydro_init_vars%fr_clay(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_sat_sc * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) &
      &                   * SPREAD(hydro_init_vars%fr_clay(:,:), DIM = 2, ncopies = nsoil) &
      &                   + k_sat_a
    theta_sat_sl(:,:,:) = theta_fc_sl(:,:,:) + k_sat_at + (1.0_wp + k_sat_bt) * hlp1(:,:,:) - &
                          k_sat_ct * SPREAD(hydro_init_vars%fr_sand(:,:), DIM = 2, ncopies = nsoil) + k_sat_dt
    ! calc root depth and soil properties in root zone
    !
    soil_awc(:,:)         = 0._wp
    w_soil_root_pwp(:,:)  = 0._wp
    w_soil_root_fc(:,:)   = 0._wp
    root_depth(:,:)       = 0._wp
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          hlp1(ic, is, iblk) = soil_lay_width_sl(ic, is, iblk) &
            &                     * (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
          soil_awc(ic, iblk)    = soil_awc(ic, iblk) + hlp1(ic, is, iblk)
          IF (soil_awc(ic, iblk) < (hydro_init_vars%soil_sat_water_content(ic, iblk) / 1000._wp)) THEN
            root_depth(ic, iblk)      = root_depth(ic, iblk) + soil_lay_width_sl(ic, is, iblk)
            w_soil_root_pwp(ic, iblk) = w_soil_root_pwp(ic, iblk) + theta_pwp_sl(ic, is, iblk) &
              &                        * soil_lay_width_sl(ic, is, iblk)
            w_soil_root_fc(ic, iblk)  = w_soil_root_fc(ic, iblk) + theta_fc_sl(ic, is, iblk) &
              &                        * soil_lay_width_sl(ic, is, iblk)
          ELSE
            hlp2(ic, iblk)            = MAX(0._wp, hydro_init_vars%soil_sat_water_content(ic, iblk) &
              &                        / 1000._wp - soil_awc(ic, iblk) + hlp1(ic, is, iblk))
            root_depth(ic, iblk)      = root_depth(ic, iblk) + hlp2(ic, iblk) &
              &                        / (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
            w_soil_root_pwp(ic, iblk) = w_soil_root_pwp(ic, iblk) + theta_pwp_sl(ic, is, iblk) * hlp2(ic, iblk) &
              &                        / (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
            w_soil_root_fc(ic, iblk)  = w_soil_root_fc(ic, iblk) + theta_fc_sl(ic, is, iblk) * hlp2(ic, iblk) &
              &                        / (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
          END IF
        END DO
      END DO
    END DO

    ! Constrain root depth max to soil depth
    root_depth(:,:) = MIN(root_depth(:,:), soil_depth(:,:))

    !> Root depth within the soil layers:
    !> calculated from the total root depth and the layer thicknesses
    !>
    dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl) = soil_depth_to_layers_2d( &
      & root_depth(:,:), &      ! Total rooting depth
      & soil_w%dz(:))           ! Soil layer thicknesses from vertical grid

    !> water in the surface reservoir [m]
    !>
    wtr_skin(:,:) = 0.0_wp

    !> adjust initial water amount to prescribed theta
    !>
    wtr_rootzone_rel(:,:)    = hydro_init_vars%soil_theta_prescribe(:,:)
    wtr_plant_avail_rel(:,:) = hydro_init_vars%soil_theta_prescribe(:,:)
    wtr_rootzone(:,:)       = w_soil_root_pwp(:,:) &
      &                       + hydro_init_vars%soil_theta_prescribe(:,:) * (w_soil_root_fc(:,:) - w_soil_root_pwp(:,:))

    !> soil water and ice pwp, fc and sat per layer depth
    !>
    wtr_soil_pwp_sl(:,:,:)  = 0.0_wp
    wtr_soil_fc_sl(:,:,:)   = 0.0_wp
    wtr_soil_sat_sl(:,:,:)  = 0.0_wp
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          IF (soil_depth_sl(ic, is, iblk) > 0.0_wp) THEN
            wtr_soil_pwp_sl(ic, is, iblk) = theta_pwp_sl(ic, is, iblk) * soil_depth_sl(ic, is, iblk)
            wtr_soil_fc_sl(ic, is, iblk)  = theta_fc_sl(ic, is, iblk)  * soil_depth_sl(ic, is, iblk)
            wtr_soil_sat_sl(ic, is, iblk) = theta_sat_sl(ic, is, iblk) * soil_depth_sl(ic, is, iblk)
          END IF
        END DO
      END DO
    END DO

    !> adjust initial water amount to prescribed theta
    !> avoid water stress (i.e., dry soil) of plants at the first timestep
    !>
    wtr_soil_sl(:,:,:) = wtr_soil_pwp_sl(:,:,:) &
      &                  + (SPREAD(hydro_init_vars%soil_theta_prescribe(:,:), DIM = 2, ncopies = nsoil) &
      &                  * (wtr_soil_fc_sl(:,:,:) - wtr_soil_pwp_sl(:,:,:)))

    !> init various HYDRO_ variables
    !>
    oro_stddev(:,:)     = hydro_init_vars%oro_stddev(:,:)
    weq_snow_soil(:,:)  = hydro_init_vars%weq_snow_soil(:,:)
    elevation(:,:)      = hydro_init_vars%elevation(:,:)
    fract_pond_max(:,:) = hydro_init_vars%fract_pond_max(:,:)
    weq_pond_max(:,:)   = hydro_init_vars%depth_pond_max(:,:) * (-1._wp)
    fract_org_sl(:,:,:) = hydro_init_vars%fract_org_sl(:,:,:)

    !> Steepness parameter for surface runoff
    !> (calculated from standard deviation of the orography; compare subroutine [[arno_scheme]])
    DO iblk = 1, nblks
      CALL calc_orographic_features(                     &
        ! in
        & nproma,                                        &
        & 180, & !hgrid%nlat_g,                                  & ! (Effective) number of latitudes for steepness parameter
        & oro_stddev(:,iblk),                            &
        ! out
        & dsl4jsb_var_ptr(HYDRO_, steepness  ) (:,iblk)  &
        & )
    END DO

    !> Calculate the soil properties from soil textures
    !>
    ! Calculate volumetric soil porosity of mineral soil
    CALL soil_init_from_texture( &
      & porosity_sand, porosity_silt, porosity_clay, porosity_oc,                                                &  ! in
      & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &  ! in
      & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &  ! in
      & hydro_init_vars%vol_porosity(:,:))                                                                          ! out

    ! Calculate saturated hydraulic conductivity of mineral soil
    CALL soil_init_from_texture( &
      & hyd_cond_sand, hyd_cond_silt, hyd_cond_clay, hyd_cond_oc,                                                &
      & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
      & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
      & hydro_init_vars%hyd_cond_sat(:,:))
    hydro_init_vars%hyd_cond_sat(:,:) = hydro_init_vars%hyd_cond_sat(:,:) * hyd_cond_sat_profile

    ! Calculate volumetric field capacity [m/m] of mineral soil
    CALL soil_init_from_texture( &
      & field_cap_sand, field_cap_silt, field_cap_clay, field_cap_oc,                                            &
      & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
      & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
      & hydro_init_vars%vol_field_cap(:,:))

    ! Calculate plant wilting point of mineral soil
    CALL soil_init_from_texture( &
      & wilt_sand, wilt_silt, wilt_clay, wilt_oc,                                                                &
      & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
      & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
      & hydro_init_vars%vol_p_wilt(:,:))

    ! Calculate residual soil water content of mineral soil
    CALL soil_init_from_texture( &
       & wres_sand, wres_silt, wres_clay, wres_oc,                                                                &
       & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
       & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
       & hydro_init_vars%vol_wres(:,:))

    ! Calculate matric potential of mineral soil
    CALL soil_init_from_texture( &
       & matric_pot_sand, matric_pot_silt, matric_pot_clay, matric_pot_oc,                                        &
       & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
       & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
       & hydro_init_vars%matric_pot(:,:))

    ! Calculate pore_size_index of mineral soil
    CALL soil_init_from_texture( &
      & pore_size_index_sand, pore_size_index_silt, pore_size_index_clay, pore_size_index_oc,                    &
      & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
      & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
      & hydro_init_vars%pore_size_index(:,:))

    ! Calculate exponent b in Clapp & Hornberger of mineral soil
    CALL soil_init_from_texture( &
      & bclapp_sand, bclapp_silt, bclapp_clay, bclapp_oc,                                                        &
      & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
      & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
      & hydro_init_vars%bclapp(:,:))

    ! init memory variables
    dsl4jsb_var2D_onDomain(HYDRO_, vol_porosity)     = hydro_init_vars%vol_porosity(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, hyd_cond_sat)     = hydro_init_vars%hyd_cond_sat(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, vol_field_cap)    = hydro_init_vars%vol_field_cap(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, vol_p_wilt)       = hydro_init_vars%vol_p_wilt(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, vol_wres)         = hydro_init_vars%vol_wres(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, matric_pot)       = hydro_init_vars%matric_pot(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, pore_size_index)  = hydro_init_vars%pore_size_index(:,:)
    dsl4jsb_var2D_onDomain(HYDRO_, bclapp)           = hydro_init_vars%bclapp(:,:)

    !> Residual water content:
    !> calculated from volumetric residual soil moisture and soil layer thicknesses
    ! Note: At this point a potential ice fraction is ignored.
    DO is = 1,nsoil
      dsl4jsb_var_ptr(HYDRO_, wtr_soil_res_sl) (:,is,:) = hydro_init_vars%vol_wres(:,:) * &
        & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl) (:,is,:)
    END DO

    !> Set maximum root zone soil moisture (water plus ice)
    !>
    ! Note: This corresponds to the water content at field capacity.
    dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = &
      & hydro_init_vars%vol_field_cap * hydro_init_vars%root_depth
    ! Recalculate maximum root zone soil moisture in case a maximum was set (default is no limit)
    IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
      dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MIN( &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
    END IF

    !> Update the soil properties according to boundary and initial conditions just set
    !>
    ! output variables:
    !     hyd_cond_sat_sl,    vol_porosity_sl,  bclapp_sl,     matric_pot_sl,
    !     pore_size_index_sl, vol_field_cap_sl, vol_p_wilt_sl, vol_wres_sl
    CALL init_soil_properties(tile)

    ! Make sure the soil water and ice content is not larger than the pore volume given by the volumetric soil
    ! porosity. Soil water and ice content need to be reduced if this limit is exceeded.
    IF (ANY(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)     + dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) > &
      &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl))) THEN
      ! First limit soil water to maximum soil water capacity
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = MAX(0._wp, &
        & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl),     &
        &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)))
      ! Then limit soil ice to the remaining capacity (without soil water)
      dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) = MAX(0._wp, &
        & MIN(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl),     &
        &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) &
        &     - dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)))
    END IF

    !> Define the amount of soil water in saturated soils, at field capacity, at the permanent wilting point
    !> and the residual water content according to the soil properties updated above.
    ! TODO: At the time hydro_init is called ice_soil_sl hasn't been initialized, yet, i.e. it should be 0 ?
    dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl) = MAX(0._wp,                                          &
      & dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
      & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
    dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) =  MAX(0._wp,                                          &
      & dsl4jsb_var3D_onDomain(HYDRO_, vol_field_cap_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
      & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
    dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl) = MAX(0._wp,                                          &
      & dsl4jsb_var3D_onDomain(HYDRO_, vol_p_wilt_sl)    * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
      & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
    dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl) = MAX(0._wp,                                          &
      & dsl4jsb_var3D_onDomain(HYDRO_, vol_wres_sl)      * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
      & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))

    ! Soils with minimum field capacity do not have organic fractions (numerical reasons)
    WHERE (dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) <= 1.0E-10_wp)
      dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
    ENDWHERE

    !> Re-compute maximum root zone moisture based on the updated soil parameters
    !>
    ! Note: All plant related computations use wpi_rootzone_max reduced by the ice and supercooled water
    !       content, which represents the unfrozen part of the soil.
    !       This assumes that plants retract/extents their roots immediately if the soil freezes/thaws.
    DO iblk = 1,SIZE(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max), 2)
        CALL get_amount_in_rootzone( &
        & dsl4jsb_var_ptr(HYDRO_, wtr_soil_fc_sl)   (:,:,iblk), &   ! in
        & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)    (:,:,iblk), &   ! in
        & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)    (:,:,iblk), &   ! in
        & dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max) (:,  iblk)  )   ! inout
    END DO
    ! Set maximum root zone moisture to 0.2 for glacier and ocean areas (where wpi_rootzone_max == 0)
    dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MERGE(       &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), 0.2_wp, &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) > 0._wp )
    ! Re-calculate the maximum rootzone moisture in case an upper limit was defined by namelist.
    IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
      dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MIN( &
      & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
    END IF

    ! Deallocate local allocatable variables
    DEALLOCATE(soil_awc)
    DEALLOCATE(hlp2)
    DEALLOCATE(hlp1)
    DEALLOCATE(theta_pwp_sl)
    DEALLOCATE(theta_fc_sl)
    DEALLOCATE(theta_sat_sl)

    IF (debug_on()) CALL message(TRIM(routine), 'Done with tile '//TRIM(tile%name)//' ...')

  END SUBROUTINE qs_hydro_init_ic_bc
#else
  ! ====================================================================================================== !
  !
  !>#### Read boundary and initial conditions of the hydrology
  !>
  !> In this subroutine we read boundary and initial condition data from the respective input files
  !> and, if necessary, derive further boundary conditions from this data.
  !> All input file names are configurable via namelist: Soil properties are read from
  !> [[t_hydro_config:bc_filename]], orographic data from [[t_hydro_config:bc_sso_filename]],
  !> and initial conditions from [[t_hydro_config:ic_filename]] or [[t_jsb_model_config:ifs_filename]],
  !> depending on parameter [[t_jsb_model_config:init_from_ifs]] (compare namelists 'jsb_hydro_nml'
  !> and 'jsb_model_nml').
  !>
  !> Soil parameters and initial conditions need to be set on all tiles on which the HYDRO_ process
  !> is active. Currently, the data is identical for all tiles. In this routine we therefore read
  !> the data just once - by the top tile - and store it in data structure
  !> [[mo_hydro_init:hydro_init_vars]]. The remaining tiles then use this data structure to initialize
  !> their memories in subroutines [[hydro_init_ic]] and [[hydro_init_bc]].
  !>
  SUBROUTINE hydro_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile

    dsl4jsb_Def_config(HYDRO_)                  !< Configurable parameters of the HYDRO process

    REAL(wp), POINTER :: &
      & ptr_1D(:    ),   &                      !< Pointer to 1d variable
      & ptr_2D(:,  :),   &                      !< Pointer to 2d variable
      & ptr_3D(:,:,:)                           !< Pointer to 3d variable
    LOGICAL,  POINTER   :: is_in_domain(:,:)    !< True: the cell belongs to domain; False: halo cell

    TYPE(t_jsb_model), POINTER :: model             !< Current instance of the model
    TYPE(t_jsb_grid),  POINTER :: grid              !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w            !< Vertical grid used for hydrology
    TYPE(t_input_file)         :: input_file        !< File containing input data to be read
    INTEGER                    :: nproma, nblks     !< Dimensions of the horizontal grid
    INTEGER                    :: nsoil             !< Number of vertical soil layers of the vgrid 'soil_depth_water'
    INTEGER                    :: nsoil_input_file  !< Number of vertical soil layers as provided by the ic input file
    INTEGER                    :: isoil             !< Looping index for soil layers
    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_read_init_vars'

    ! Assertion: do nothing if hydrology process is not active on this tile (should not happen)
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Reading/setting hydrology init vars')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> - Allocate memory for the variables of the inital data structure [[mo_hydro_init:hydro_init_vars]]
    !>
    ALLOCATE( &
      & hydro_init_vars%oro_stddev      (nproma, nblks),       &
      & hydro_init_vars%soil_depth      (nproma, nblks),       &
      & hydro_init_vars%vol_porosity    (nproma, nblks),       &
      & hydro_init_vars%hyd_cond_sat    (nproma, nblks),       &
      & hydro_init_vars%matric_pot      (nproma, nblks),       &
      & hydro_init_vars%bclapp          (nproma, nblks),       &
      & hydro_init_vars%pore_size_index (nproma, nblks),       &
      & hydro_init_vars%vol_field_cap   (nproma, nblks),       &
      & hydro_init_vars%vol_p_wilt      (nproma, nblks),       &
      & hydro_init_vars%vol_wres        (nproma, nblks),       &
      & hydro_init_vars%root_depth      (nproma, nblks),       &
      & hydro_init_vars%fract_pond_max  (nproma, nblks),       &
      & hydro_init_vars%depth_pond_max  (nproma, nblks),       &
      & hydro_init_vars%weq_snow_soil   (nproma, nblks)        & ! TODO only needed for init_from_ifs == .FALSE.
      & )

      ! Initialize input variables with missing value
      hydro_init_vars%oro_stddev     (:,:)   = missval
      hydro_init_vars%soil_depth     (:,:)   = missval
      hydro_init_vars%vol_porosity   (:,:)   = missval
      hydro_init_vars%hyd_cond_sat   (:,:)   = missval
      hydro_init_vars%matric_pot     (:,:)   = missval
      hydro_init_vars%bclapp         (:,:)   = missval
      hydro_init_vars%pore_size_index(:,:)   = missval
      hydro_init_vars%vol_field_cap  (:,:)   = missval
      hydro_init_vars%vol_p_wilt     (:,:)   = missval
      hydro_init_vars%vol_wres       (:,:)   = missval
      hydro_init_vars%root_depth     (:,:)   = missval
      hydro_init_vars%fract_pond_max (:,:)   = missval
      hydro_init_vars%depth_pond_max (:,:)   = missval
      hydro_init_vars%weq_snow_soil  (:,:)   = missval     ! TODO only needed for init_from_ifs == .FALSE.

    ! Further memory allocation - depending on the setup

    ! TODO only necessary for l_socmap == .TRUE.
    ALLOCATE(hydro_init_vars%fract_org_sl(nproma, nsoil, nblks))
    hydro_init_vars%fract_org_sl(:,:,:) = missval

    IF (jsbach_runs_standalone()) THEN
       ALLOCATE(hydro_init_vars%elevation(nproma,nblks))
       hydro_init_vars%elevation(:,:) = missval
    END IF
    IF (model%config%init_from_ifs) THEN
      ALLOCATE(hydro_init_vars%ifs_weq_snow(nproma,nblks))
      hydro_init_vars%ifs_weq_snow(:,:) = missval
    END IF
    IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
      ALLOCATE(hydro_init_vars%ifs_smi_sl(nproma,ifs_nsoil,nblks))
      hydro_init_vars%ifs_smi_sl(:,:,:) = missval
    ELSE
      ALLOCATE(hydro_init_vars%wtr_soil_sl(nproma, nsoil, nblks))
      hydro_init_vars%wtr_soil_sl(:,:,:) = missval
    END IF
    IF (dsl4jsb_Config(HYDRO_)%l_soil_texture) THEN
      ALLOCATE( &
        & hydro_init_vars%fr_sand     (nproma, nblks), &
        & hydro_init_vars%fr_silt     (nproma, nblks), &
        & hydro_init_vars%fr_clay     (nproma, nblks), &
        & hydro_init_vars%fr_sand_deep(nproma, nblks), &
        & hydro_init_vars%fr_silt_deep(nproma, nblks), &
        & hydro_init_vars%fr_clay_deep(nproma, nblks)  &
        & )
      hydro_init_vars%fr_sand     (:,:) = missval
      hydro_init_vars%fr_silt     (:,:) = missval
      hydro_init_vars%fr_clay     (:,:) = missval
      hydro_init_vars%fr_sand_deep(:,:) = missval
      hydro_init_vars%fr_silt_deep(:,:) = missval
      hydro_init_vars%fr_clay_deep(:,:) = missval
    END IF ! l_soil_texture

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> - Read orographic properties from input file [[t_hydro_config:bc_sso_filename]]
    !>
    ! Open input file - bc_sso
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_sso_filename), model%grid_id)

    !>    - Orography and its standard deviation
    IF (jsbach_runs_standalone()) THEN
       ptr_2D => input_file%Read_2d(  &
         & variable_name='elevation', &
         & fill_array = hydro_init_vars%elevation)
       hydro_init_vars%elevation = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
    END IF

    ptr_2D => input_file%Read_2d(  &
      & variable_name='orostd',    &
      & fill_array = hydro_init_vars%oro_stddev)
    hydro_init_vars%oro_stddev = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    !>     - Maximum fraction and depth of local surface depressions needed with the ponds scheme
    IF (dsl4jsb_Config(HYDRO_)%l_ponds .AND. .NOT. tile%is_lake .AND. .NOT. tile%is_glacier) THEN
      IF (input_file%Has_var('surf_depr_fract') .AND. input_file%Has_var('surf_depr_depth')) THEN
        ! Use depression area fraction as maximum pond fraction
        ptr_2D => input_file%Read_2d(                                 &
          & variable_name='surf_depr_fract',                           &
          & fill_array=hydro_init_vars%fract_pond_max)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        ! Use depression depth as maximum pond depth (Attention: Values are negative!)
        ptr_2D => input_file%Read_2d(                                 &
          & variable_name='surf_depr_depth',                           &
          & fill_array=hydro_init_vars%depth_pond_max)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) <= 0._wp)
      ELSE
        CALL finish(TRIM(routine), '*** Error: surf_depr_fract and/or surf_depr_depth not found in bc_sso file. '// &
          & 'Pond scheme cannot run without these fields')
      END IF
    END IF

    ! Close bc_sso input file
    CALL input_file%Close()

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> - Read soil properties or soil textures from input file [[t_hydro_config:bc_filename]]
    !>
    !>    Depending on namelist parameter [[t_hydro_config:l_soil_texture]] (namelist 'jsb_hydro_nml')
    !>    soil properties are either read directly or they are calculated from soil textures.

    IF (tile%contains_soil) THEN ! for tiles containing soil

      ! Open ic input file to get the number of soil layers
      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%ic_filename))
      ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of the bottom of each soil layer
      nsoil_input_file = SIZE(ptr_1D)                          ! number of soil layers provided by the ic input file
      ! Close input file - ic
      CALL input_file%Close()

      ! Open input file - bc
      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename), model%grid_id)

      !>
      !>    A: Soil properties read in all configurations
      !>

      !>     - Soil depth until bedrock
      ptr_2D => input_file%Read_2d(   &
        & variable_name='soil_depth', &
        & fill_array = hydro_init_vars%soil_depth)
      hydro_init_vars%soil_depth = MAX(0.1_wp, MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp))

      IF (.NOT. dsl4jsb_Config(HYDRO_)%l_soil_texture) THEN
        !>
        !>    B: Soil properties only read if soil textures are ignored
        !>
        !>     - Volumetric soil porosity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & fill_array = hydro_init_vars%vol_porosity, &
            & variable_name='soil_porosity_mineral')
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & fill_array = hydro_init_vars%vol_porosity, &
            & variable_name='soil_porosity')
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

        !>     - Saturated hydraulic conductivity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='hyd_cond_sat_mineral',      &
            & fill_array = hydro_init_vars%hyd_cond_sat)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='hyd_cond_sat',              &
            & fill_array = hydro_init_vars%hyd_cond_sat)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 4.e-6_wp, ptr_2D(:,:) >= 0._wp)

        !>     - Volumetric field capacity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='soil_field_cap_mineral',    &
            & fill_array = hydro_init_vars%vol_field_cap)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='soil_field_cap',            &
            & fill_array = hydro_init_vars%vol_field_cap)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Volumetric permanent wilting point
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='wilting_point_mineral',     &
            & fill_array=hydro_init_vars%vol_p_wilt)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='wilting_point',             &
            & fill_array=hydro_init_vars%vol_p_wilt)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Volumetric residual water content (if missing deduce from porosity)
        ! TODO Remove deduction of residual water content. This was only temporarily implemented
        !      while bc_land files without this quantity were still around and used
        IF (input_file%Has_var('residual_water')) THEN
          IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
            ptr_2D => input_file%Read_2d(                &
              & variable_name='residual_water_mineral',  &
              & fill_array=hydro_init_vars%vol_wres)
          ELSE
            ptr_2D => input_file%Read_2d(                &
              & variable_name='residual_water',          &
              & fill_array=hydro_init_vars%vol_wres)
          END IF
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        ELSE
          WRITE (message_text,*) 'BC File does not contain data on residual water content. Will use approximation instead!'
          CALL message (TRIM(routine), message_text)
          hydro_init_vars%vol_wres = hydro_init_vars%vol_porosity * 0.2_wp
        END IF

        !>     - Soil matric potential
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='moisture_pot_mineral',      &
            & fill_array = hydro_init_vars%matric_pot)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='moisture_pot',              &
            & fill_array = hydro_init_vars%matric_pot)
        END IF
        ! @todo There is still old bc data around where the matric potential is given with a positiv sign while
        !       it should be negativ from a soil physics perspective and is expected as such from the soil hydrology
        !       routines. Until this bc data is not used anymore, the following check is required to ensure the
        !       sign of the soil matric potential is (converted) correct.
        ! Exclude values outside of the domain
        WHERE(.NOT. is_in_domain(:,:))
          ptr_2D(:,:) = 0._wp
        ENDWHERE
        ! Exclude (standard) missing values
        WHERE(ABS(ptr_2D(:,:)) > 1.0e+20_wp)
          ptr_2D(:,:) = 0._wp
        ENDWHERE
        ! Check if all values are negative or positive (negative values are needed)
        !   and convert positive -> negative values
        ! TODO please revise if still is correct and necessary with current ic/bc files
        IF (ALL(ptr_2D(:,:) < 1.0e-10_wp)) THEN
          ! Correct: all values are negative
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), -0.2_wp, ptr_2D(:,:) <= 0._wp)
        ELSE IF (ALL(ptr_2D(:,:) > -1.0e-10_wp)) THEN
          ! False: all values are positive --> convert to negative values
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp) * (-1._wp)
        ELSE
          ! Totally wrong: dataset contains positive and negativ values
          WRITE (message_text,*) 'Found positive and negative values in soil moisture matric potential: ', &
            & MINVAL(ptr_2D(:,:)),' <-> ',MAXVAL(ptr_2D(:,:))
          CALL finish (TRIM(routine), message_text)
        END IF

        !>     - Exponent B in Clapp and Hornberger
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='bclapp_mineral',            &
            & fill_array = hydro_init_vars%bclapp)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='bclapp',                    &
            & fill_array = hydro_init_vars%bclapp)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 6._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Pore size distribution index
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='pore_size_index_mineral',   &
            & fill_array = hydro_init_vars%pore_size_index)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='pore_size_index',           &
            & fill_array = hydro_init_vars%pore_size_index)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.25_wp, ptr_2D(:,:) >= 0.05_wp)

      ELSE   ! l_soil_texture
        !>
        !>    C: Soil textur maps - read if soil properties are to be calculated,
        !>       and then calculate the soil properties from these textures.
        !>
        !>     - Fraction of sand in upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_SAND',                   &
          & fill_array=hydro_init_vars%fr_sand)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)  ! ptr_2D points to hydro_init_vars%fr_sand, hence
                                                                       !   hydro_init_vars%fr_sand is modified here

        !>     - Fraction of silt in upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_SILT',                   &
          & fill_array=hydro_init_vars%fr_silt)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Fraction of clay in upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_CLAY',                   &
          & fill_array=hydro_init_vars%fr_clay)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Fraction of sand in deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_SAND',               &
          & fill_array=hydro_init_vars%fr_sand_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Fraction of silt in deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_SILT',               &
          & fill_array=hydro_init_vars%fr_silt_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !>     - Fraction of clay in deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_CLAY',               &
          & fill_array=hydro_init_vars%fr_clay_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        !
        ! Calculate the soil properties from soil textures
        !
        ! Calculate volumetric soil porosity of mineral soil
        CALL soil_init_from_texture( &
          & porosity_sand, porosity_silt, porosity_clay, porosity_oc,                                                &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%vol_porosity(:,:))

        ! Calculate saturated hydraulic conductivity of mineral soil
        CALL soil_init_from_texture( &
          & hyd_cond_sand, hyd_cond_silt, hyd_cond_clay, hyd_cond_oc,                                                &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%hyd_cond_sat(:,:))
        hydro_init_vars%hyd_cond_sat(:,:) = hydro_init_vars%hyd_cond_sat(:,:) * hyd_cond_sat_profile

        ! Calculate volumetric field capacity [m/m] of mineral soil
        CALL soil_init_from_texture( &
          & field_cap_sand, field_cap_silt, field_cap_clay, field_cap_oc,                                            &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%vol_field_cap(:,:))

        ! Calculate plant wilting point of mineral soil
        CALL soil_init_from_texture( &
          & wilt_sand, wilt_silt, wilt_clay, wilt_oc,                                                                &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%vol_p_wilt(:,:))

        ! Calculate residual soil water content of mineral soil
        CALL soil_init_from_texture( &
           & wres_sand, wres_silt, wres_clay, wres_oc,                                                                &
           & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           & hydro_init_vars%vol_wres(:,:))

        ! Calculate matric potential of mineral soil
        CALL soil_init_from_texture( &
           & matric_pot_sand, matric_pot_silt, matric_pot_clay, matric_pot_oc,                                        &
           & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           & hydro_init_vars%matric_pot(:,:))

        ! Calculate pore_size_index of mineral soil
        CALL soil_init_from_texture( &
          & pore_size_index_sand, pore_size_index_silt, pore_size_index_clay, pore_size_index_oc,                    &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%pore_size_index(:,:))

        ! Calculate exponent b in Clapp & Hornberger of mineral soil
        CALL soil_init_from_texture( &
          & bclapp_sand, bclapp_silt, bclapp_clay, bclapp_oc,                                                        &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%bclapp(:,:))

      END IF ! l_soil_texture

      !>
      !>    D: Further soil data
      !>
      !>     - Soil organic carbon fractions for each soil layer
      IF (dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        DO isoil = 1,nsoil_input_file
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='fract_org_sl',    &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%fract_org_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO
        ! apply the value from the lowest soil layer (of the input file) to all soil layers below
        IF (nsoil > nsoil_input_file) THEN
          DO isoil = nsoil_input_file+1, nsoil
            hydro_init_vars%fract_org_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
          END DO
        END IF
      END IF

      ! TODO Move up reading root depth to section of variables read in all configurations to
      !      get documentation (and code) clearer.
      !>     - Root depth
      ! TODO root depth should become a function of PFTs eventually, instead of being a soil property
      ptr_2D => input_file%Read_2d(   &
        & variable_name='root_depth', &
        & fill_array = hydro_init_vars%root_depth)
      hydro_init_vars%root_depth = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
      ! Make sure there is root depth equivalent to 0.2 m weq in all cells
      !   that contain field capacity data but no root depth
      WHERE (hydro_init_vars%vol_field_cap(:,:) > 1.0e-10_wp)
        hydro_init_vars%root_depth(:,:) = MAX(hydro_init_vars%root_depth(:,:), &
          & 0.2_wp / hydro_init_vars%vol_field_cap(:,:))
      ELSEWHERE
        hydro_init_vars%root_depth(:,:) = 0._wp
      ENDWHERE
      ! Constrain root depth because in the bc file it may be larger then soil depth
      hydro_init_vars%root_depth = MIN(hydro_init_vars%root_depth, hydro_init_vars%soil_depth)

      ! Close input file - bc
      CALL input_file%Close()

      !>
      !> -  Read initial conditons for the hydrology
      !>
      !>    Depending on the configuration, the hydrology is either initialized with data gained
      !>    from previous simulations (read from file [[t_hydro_config:ic_filename]]), or from
      !>    IFS data (read from file [[t_jsb_model_config:ifs_filename]]).
      !>
      ! Note: The initial conditions are also read in restarted simulations, but the restart
      !       file will be read later and the variables read here will be overwritten.

      IF (.NOT. model%config%init_from_ifs .OR. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        input_file = jsb_netcdf_open_input(dsl4jsb_Config(HYDRO_)%ic_filename, model%grid_id)
      END IF

      !>     - Initial snow depth
      IF (model%config%init_from_ifs) THEN
        ptr_3D => model%config%ifs_input_file%Read_2d_time( &
          & variable_name='W_SNOW', start_time_step=1, end_time_step=1)
        hydro_init_vars%ifs_weq_snow(:,:) = ptr_3D(:,:,1)

      ELSE
        IF (input_file%Has_var('snow')) THEN
          ptr_2D => input_file%Read_2d( &
            & variable_name='snow',     &
            & fill_array = hydro_init_vars%weq_snow_soil)
          hydro_init_vars%weq_snow_soil = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
        ELSE
          hydro_init_vars%weq_snow_soil = 0._wp
          CALL message(TRIM(routine), 'Initializing snow to zero')
        END IF

      END IF

      !>     - Initial water content
      IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN

        DO isoil=1,ifs_nsoil
          ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
            & variable_name='SMIL'//int2string(isoil))
          hydro_init_vars%ifs_smi_sl(:,isoil,:) = ptr_2D(:,:)
        END DO

      ELSE

        DO isoil = 1,nsoil_input_file
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='layer_moist',     &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%wtr_soil_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO
        ! apply the value from the lowest soil layer (of the input file) to all soil layers below
        IF (nsoil > nsoil_input_file) THEN
          DO isoil = nsoil_input_file+1, nsoil
            hydro_init_vars%wtr_soil_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
          END DO
        END IF

      END IF ! init_from_ifs

      IF (.NOT. model%config%init_from_ifs .OR. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        CALL input_file%Close()
      END IF

    END IF ! tile contains_soil

  END SUBROUTINE hydro_read_init_vars

  ! ====================================================================================================== !
  !
  !>#### Set boundary conditions for the hydrology process on each tile
  !>
  !> The subroutine makes use of the data structure [[mo_hydro_init:hydro_init_vars]] containing the
  !> input variables that had been read in subroutine [[hydro_read_init_vars]] and provides this
  !> data to all tiles on which the HYDRO process is active and the data is needed.
  !
  SUBROUTINE hydro_init_bc(tile)

    USE mo_util,                ONLY: soil_depth_to_layers_2d
    USE mo_hydro_process,       ONLY: calc_orographic_features
    USE mo_jsb_math_constants,  ONLY: eps8

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Current tile

    dsl4jsb_Def_config(HYDRO_)                        !< Configuration of the HYDRO process
    dsl4jsb_Def_memory(HYDRO_)                        !< Memory of the HYDRO process

    TYPE(t_jsb_model), POINTER :: model               !< Current instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER :: grid                !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w              !< Vertical grid used for hydrology
    INTEGER                    :: nsoil               !< Number of soil layers

    INTEGER :: i                                      !< Looping index
    INTEGER :: nblks, nproma                          !< dimensions
    INTEGER :: iblk, ic, is                           !< Looping indices

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_bc'

    ! HYDRO_ 2D
    dsl4jsb_Real2D_onDomain    :: num_sl_above_bedrock
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_center_sl

    ! Do nothing if hydrology process is not active on this tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting hydrology boundary conditions for tile '//TRIM(tile%name))

    model  => Get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil  = soil_w%n_levels
    nblks  = grid%nblks
    nproma = grid%nproma

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)

#ifndef __NO_QUINCY__
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      ! HYDRO_ 2D
      dsl4jsb_Get_var2D_onDomain(HYDRO_,   num_sl_above_bedrock)
      ! HYDRO_ 3D
      dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_width_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_lbound_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_ubound_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_center_sl)
    END SELECT
#endif

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> - Set orographic parameters
    !
    !>     - Elevation (for standalone model)
    IF (jsbach_runs_standalone()) dsl4jsb_var2D_onDomain(HYDRO_,elevation) = hydro_init_vars%elevation

    !>     - Standard deviation of the orography
    dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev) = hydro_init_vars%oro_stddev

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_,elevation)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev))

    !>     - Steepness parameter for surface runoff
    !>       (calculated from standard deviation of the orography; compare subroutine [[arno_scheme]])
    DO i = 1, grid%Get_nblks()
      CALL calc_orographic_features(                    &
        ! in
        & grid%Get_nproma(),                            &
        & grid%nlat_g,                                  & ! (Effective) number of latitudes for steepness parameter
        & dsl4jsb_var_ptr(HYDRO_, oro_stddev ) (:,i),   &
        ! out
        & dsl4jsb_var_ptr(HYDRO_, steepness  ) (:,i)    &
        & )
    END DO

    !>     - Maximum depth and extent of surface water ponds
    IF (.NOT. tile%is_lake .AND. .NOT. tile%is_glacier) THEN

      IF (dsl4jsb_Config(HYDRO_)%l_ponds) THEN
        ! Surface depression depth is given as negative value in the data, but represents
        ! a storage on top of the soil and thus needs to be converted to positive values.
        dsl4jsb_var_ptr(HYDRO_, fract_pond_max) = hydro_init_vars%fract_pond_max(:,:)
        dsl4jsb_var_ptr(HYDRO_, weq_pond_max)   = hydro_init_vars%depth_pond_max(:,:) * (-1._wp)
      ELSE
        ! Note, that ponds are implicitely disabled by setting the maximum allowed pond fraction
        !   and depth to zero. Thus, we can avoid a large number of l_ponds conditions for surface
        !   and soil hydrology processes.
        dsl4jsb_var_ptr(HYDRO_, fract_pond_max) = 0._wp
        dsl4jsb_var_ptr(HYDRO_, weq_pond_max)   = 0._wp
      END IF
    END IF

    IF (tile%contains_soil) THEN

      ! ------------------------------------------------------------------------------------------------ !
      !>
      !> - Set soil depths
      !
      !>     - Total soil depth until bedrock
      dsl4jsb_var2D_onDomain(HYDRO_, soil_depth) = hydro_init_vars%soil_depth

      !>     - Soil depth within each soil layer:
      !>       calculated from total soil depth until bedrock and the fixed layer thicknesses of the
      !>       vertical grid.
      dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, soil_depth), & ! Total soil depth until bedrock
        & soil_w%dz(:))                                 ! Fixed soil layer thicknesses from vertical grid
      !>

#ifndef __NO_QUINCY__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        ! pass these values to a HYDRO_ variable with an improved name
        soil_lay_width_sl(:,:,:) = dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)

        !>  calc more metrics of the soil layers
        !>
        ! i)  lower & upper boundary
        ! ii) depth at the center of each layer
        DO iblk = 1,nblks
          DO is = 1,nsoil
            ! lower & upper boundary
            IF (is == 1) THEN
              soil_lay_depth_lbound_sl(:, is, iblk) = 0.0_wp
              soil_lay_depth_ubound_sl(:, is, iblk) = soil_lay_width_sl(:, is, iblk)
            ELSE
              soil_lay_depth_lbound_sl(:, is, iblk) = soil_lay_depth_lbound_sl(:, is-1, iblk) &
                &                                     + soil_lay_width_sl(:, is-1, iblk)
              soil_lay_depth_ubound_sl(:, is, iblk) = soil_lay_depth_ubound_sl(:, is-1, iblk) &
                &                                     + soil_lay_width_sl(:, is, iblk)
            ENDIF
            ! soil-layer center
            soil_lay_depth_center_sl(:, is, iblk) = (soil_lay_depth_lbound_sl(:, is, iblk) &
              &                                     + soil_lay_depth_ubound_sl(:, is, iblk)) &
              &                                     * 0.5_wp
          END DO
        END DO

        !> get number of soil layers above bedrock for each gridcell
        !>
        !>    use this variable in quincy for looping over soil layers,
        !>    excluding soil layers with a width smaller/equal eps8
        !>
        DO iblk = 1,nblks
          DO ic = 1,nproma
            DO is = 1,nsoil
              IF (soil_lay_width_sl(ic, is, iblk) > eps8) THEN
                num_sl_above_bedrock(ic, iblk) = REAL(is, wp)
              ELSE
                EXIT  ! stop looping over further soil layers with thickness < eps8
              END IF
            END DO
          END DO
        END DO
      END SELECT

      !$ACC UPDATE DEVICE(soil_lay_width_sl(:,:,:), soil_lay_depth_lbound_sl(:,:,:)) ASYNC(1)
      !$ACC UPDATE DEVICE(soil_lay_depth_ubound_sl(:,:,:), soil_lay_depth_center_sl(:,:,:), num_sl_above_bedrock(:,:)) ASYNC(1)
#endif

      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set hydrological properties
      !
      !>     - Volumetric soil porosity
      ! TODO: move down setting volumetric soil porosity to where it is used: just before the block
      !       where wtr_soil_sat_s is calculated. That would make the code clearer.
      dsl4jsb_var_ptr(HYDRO_, vol_porosity) = hydro_init_vars%vol_porosity

      !>     - Saturated hydraulic conductivity
      dsl4jsb_var_ptr(HYDRO_, hyd_cond_sat) = hydro_init_vars%hyd_cond_sat

      !>     - Matric potential
      dsl4jsb_var_ptr(HYDRO_, matric_pot) = hydro_init_vars%matric_pot

      !>     - Clapp and Hornberger exponent b
      dsl4jsb_var_ptr(HYDRO_, bclapp) = hydro_init_vars%bclapp

      !>     - Pore size index
      dsl4jsb_var_ptr(HYDRO_, pore_size_index) = hydro_init_vars%pore_size_index

      !>     - Volumetric field capacity
      dsl4jsb_var_ptr(HYDRO_, vol_field_cap) = hydro_init_vars%vol_field_cap

      !>     - Water content at field capacity:
      !>       calculated from volumetric field capacity and soil layer thicknesses
      !   Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,i,:) = hydro_init_vars%vol_field_cap(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Volumetric permanent wilting point
      ! TODO: Do we need this for bare soil tiles?
      dsl4jsb_var_ptr(HYDRO_, vol_p_wilt) = hydro_init_vars%vol_p_wilt

      !>     - Permanent wilting point:
      !>       calculated from volumetric permanent wilting point and soil layer thicknesses
      !   Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl) (:,i,:) = hydro_init_vars%vol_p_wilt(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Volumetric residual water content
      dsl4jsb_var_ptr(HYDRO_, vol_wres) = hydro_init_vars%vol_wres

      !>     - Residual water content:
      !>       calculated from volumetric residual soil moisture and soil layer thicknesses
      !   Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_res_sl) (:,i,:) = hydro_init_vars%vol_wres(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Water holding capacity (water content at saturation):
      !>       calculated from volumetric porosity and soil layer thicknesses
      !   Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_sat_sl) (:,i,:) = dsl4jsb_var2D_onDomain(HYDRO_,vol_porosity) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Soil organic carbon fraction
      IF (dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        dsl4jsb_var_ptr(HYDRO_, fract_org_sl) = hydro_init_vars%fract_org_sl
      END IF

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_sat_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_pwp_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_fc_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_res_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_field_cap)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_p_wilt)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_wres)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, pore_size_index)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, bclapp)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, matric_pot)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, hyd_cond_sat)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_porosity)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, fract_org_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, fract_pond_max)), &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, weq_pond_max)), &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, soil_depth))

    END IF ! contains_soil

    IF (tile%contains_vegetation) THEN ! tile contains vegetation
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set rooting depths
      !
      !>     - Total root depth
      dsl4jsb_var_ptr(HYDRO_, root_depth) = hydro_init_vars%root_depth

      !>     - Root depth within the soil layers:
      !>       calculated from the total root depth and the layer thicknesses
      dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, root_depth), & ! Total rooting depth
        & soil_w%dz(:))                                 ! Soil layer thicknesses from vertical grid

      ! -------------------------------------------------------------------------------------------------- !
      !
      !> - Set maximum root zone soil moisture (water plus ice)
      !  Note: This corresponds to the water content at field capacity.
      dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = &
        & hydro_init_vars%vol_field_cap * hydro_init_vars%root_depth

      ! Recalculate maximum root zone soil moisture in case a maximum was set (default is no limit)
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MIN( &
          & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, root_depth)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max))
    END IF

    !$ACC WAIT(1)

  END SUBROUTINE hydro_init_bc

  ! ====================================================================================================== !
  !
  !>#### Set initial conditions for the hydrology process
  !>
  !> In this suboutine the initial state of the hydrology is set (e.g. initial soil water or ice content).
  !>
  !> The routine makes use of the data structure [[mo_hydro_init:hydro_init_vars]] containing the
  !> initial data read in subroutine [[hydro_read_init_vars]] and provides this data to all tiles on which
  !> the HYDRO process is active and the data is needed.
  !>
  !
  SUBROUTINE hydro_init_ic(tile)

    USE mo_hydro_util,          ONLY: get_amount_in_rootzone
    USE mo_util,                ONLY: ifs2soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Current tile

    dsl4jsb_Def_config(HYDRO_)                        !< Configuration of the HYDRO process
    dsl4jsb_Def_memory(HYDRO_)                        !< Memory of the HYDRO process
    dsl4jsb_Def_memory(PHENO_)                        !< Memory of the PHENO process

    TYPE(t_jsb_model), POINTER :: model               !< Current instance of ICON-Land
    TYPE(t_jsb_vgrid), POINTER :: soil_w              !< Vertical grid used for hydrology

    INTEGER :: i                                      !< Looping index

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_ic'

    ! Do nothing if hydrology process is not active on this tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting hydrology initial conditions for tile '//TRIM(tile%name))

    model  => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation .AND. model%config%model_scheme == MODEL_JSBACH) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    IF (tile%contains_soil) THEN ! tile contains soil, i.e. non-glacier land

      ! -------------------------------------------------------------------------------------------------- !
      !> There are two options to set the initial snow amount and soil water content. Either the
      !> initial values are calculated from IFS (ECMWF Integrated Forecasting System) data
      !> ([[t_jsb_model_config:init_from_ifs]]=true; namelist 'jsb_model_nml'), or we use data from
      !> file [[t_hydro_config:ic_filename]] containing hydrological conditions based on previous
      !> model simulations. IFS data is generally available for specific initialization dates. As
      !> the IFS data is defined on different soil layers, the data needs to be interpolated to the
      !> ICON-Land soil layers.
      !>
      !> - Set amount of snow (in m water equivalent)

      IF (model%config%init_from_ifs) THEN
        dsl4jsb_var_ptr(HYDRO_, weq_snow_soil)(:,:) = hydro_init_vars%ifs_weq_snow(:,:)
        dsl4jsb_var2D_onDomain(HYDRO_, weq_snow) = dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)
      ELSE
        dsl4jsb_var_ptr(HYDRO_, weq_snow_soil) = hydro_init_vars%weq_snow_soil
        dsl4jsb_var2D_onDomain(HYDRO_, weq_snow) = dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)
      END IF

      !> - Set soil water content
      !>
      IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        !>    In case soil water is initialized from IFS data, a soil moisture index first needs to
        !>    be interpolated from the IFS soil layers to the jsbach soil layers.
        !>    Then the soil moisture index is converted to water content using field capacity and
        !>    wilting point.
        ! Note: wtr_soil_sl is used here to temporarily store the soil moisture index
        CALL ifs2soil(hydro_init_vars%ifs_smi_sl(:,:,:), ifs_soil_depth(:),        &
          &           dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:), soil_w%ubounds(:))

        ! Conversion of the soil moisture index (defined as (SM-PWP)/(FC-PWP)) to water content
        WHERE (dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) < -999._wp) ! Fill missing values
          dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) = &
            &  0.5_wp * (dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:) + dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,:,:))
        ELSE WHERE
          dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) = &
            & MIN(dsl4jsb_var_ptr(HYDRO_,wtr_soil_sat_sl)(:,:,:), &
            &     MAX(0._wp, &
            &         dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) &
            &         * (dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,:,:) - dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:)) &
            &         + dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:) &
            &        ) &
            &    )
        END WHERE

      ELSE
        ! In case soil water content is read from "ic_land_soil" no processing is needed.
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl) (:,:,:) = hydro_init_vars%wtr_soil_sl(:,:,:)

      END IF ! init_from_ifs

      ! Initial soil water overflow (used in terra planet experiments)
      dsl4jsb_var2D_onDomain(HYDRO_, tpe_overflow) = 0._wp

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_snow)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, tpe_overflow))

      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set the fraction of organic matter within the soil layers
      !>
      !>    If the soil organic matter fraction is to be explicitly taken into account
      !>    ([[t_hydro_config:l_organic]]=true, namelist 'jsb_hydro_nml') but the fractions map has not
      !>    been read from file ([[t_hydro_config:l_socmap]]=false, namelist 'jsb_hydro_nml'), the organic
      !>    fraction of the uppermost soil layer is initialized according to the forest fraction, while
      !>    it is zero in the deeper soil layers.
      ! Note: Currently fract_org_sl does not change with time

      IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN        ! Explicitly represent organic soil fractions

        SELECT CASE (model%config%model_scheme)
        CASE (MODEL_JSBACH)
          IF (.NOT. dsl4jsb_Config(HYDRO_)%l_socmap) THEN ! Not using organic soil fractions from a map

            ! Initialize with zero
            dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp

            ! A forest fraction can only be calculated on vegetation tiles
            IF (tile%is_vegetation) THEN          ! This is a vegetation tile
              IF (tile%lcts(1)%lib_id /= 0)  THEN ! lctlib is only defined for PFTs
                ! PFT tile
                ! Note: in current setup hydrology does not run on PFTs
                ! TODO: call finish unless it is approved we want to get the organic soil fraction of PFT tiles
                !       from the maximum vegetation fraction (and not from the forest fraction).
                IF (dsl4jsb_Lctlib_param(ForestFlag)) THEN
                  dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max)
                END IF
              ELSE
                ! non-PFT vegetation tile (typically 'veg')
                ! Define the organic fraction of the uppermost soil layer from the forest fraction.
                dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_forest)
              END IF
            END IF ! is_vegetation
          END IF ! l_socmap
        CASE (MODEL_QUINCY)
          ! Not using organic soil fractions from a map
          IF (.NOT. dsl4jsb_Config(HYDRO_)%l_socmap) THEN
            CALL finish(TRIM(routine), 'Error: init of HYDRO_ fract_org_sl from PHENO_ variables is not yet working with QUINCY')
          END IF
        END SELECT

        ! Soils with minimum field capacity do not have organic fractions (numerical reasons)
        WHERE (dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) <= 1.0E-10_wp)
          dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
        ENDWHERE

        !$ACC UPDATE DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)) ASYNC(1)

      END IF ! l_organic

      !> - Update the soil properties according to boundary and initial conditions just set
      !>
      ! output variables:
      !     hyd_cond_sat_sl,    vol_porosity_sl,  bclapp_sl,     matric_pot_sl,
      !     pore_size_index_sl, vol_field_cap_sl, vol_p_wilt_sl, vol_wres_sl
      CALL init_soil_properties(tile)

      ! Make sure the soil water and ice content is not larger than the pore volume given by the volumetric soil
      ! porosity. Soil water and ice content need to be reduced if this limit is exceeded.
      IF (ANY(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)     + dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) > &
        &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl))) THEN
        ! First limit soil water to maximum soil water capacity
        dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = MAX(0._wp, &
          & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl),     &
          &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)))
        ! Then limit soil ice to the remaining capacity (without soil water)
        dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) = MAX(0._wp, &
          & MIN(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl),     &
          &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) &
          &     - dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)))
      END IF

      !> - Define the amount of soil water in saturated soils, at field capacity, at the permanent wilting point
      !>   and the residual water content according to the soil properties updated above.
      ! TODO: At the time hydro_init is called ice_soil_sl hasn't been initialized, yet, i.e. it should be 0 ?
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) =  MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_field_cap_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_p_wilt_sl)    * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_wres_sl)      * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))

      !> water in the surface reservoir [m]
      !>
      dsl4jsb_var2D_onDomain(HYDRO_, wtr_skin) = 0.0_wp
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_skin)) ASYNC(1)

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl))
      !$ACC WAIT(1)

    END IF ! contains_soil

    IF (tile%contains_vegetation) THEN   ! tile contains vegetation
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Re-compute maximum root zone moisture based on the updated soil parameters
      !>
      ! Note: All plant related computations use wpi_rootzone_max reduced by the ice and supercooled water
      !       content, which represents the unfrozen part of the soil.
      !       This assumes that plants retract/extents their roots immediately if the soil freezes/thaws.
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max), 2)
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, wtr_soil_fc_sl)   (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max) (:,  i)  )
      END DO
      !$ACC UPDATE HOST(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max)) ASYNC(1)
      !$ACC WAIT(1)

      ! Set maximum root zone moisture to 0.2 for glacier and ocean areas (where wpi_rootzone_max == 0)
      dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MERGE(       &
          & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), 0.2_wp, &
          & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) > 0._wp )

      ! Re-calculate the maximum rootzone moisture in case an upper limit was defined by namelist.
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MIN( &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max)) ASYNC(1)

      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Initialize the relative root zone soil moisture, needed for the LoGro-P phenology
      !>
      !>    _Note_: It will be re-calculated every time step in subroutine
      !>    [[mo_hydro_interface:update_water_stress]].
      !>    It needs to be calculated here as the phenology is updated before the hydrology.

      ! Calculate current water and ice content of the root zone - only used to calculate relative soil
      ! moisture for phenology. It is re-calculated in mo_hydro_interface:update_water_stress.
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone), 2)
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, wtr_soil_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)   (:,  i)  )
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, ice_soil_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, ice_rootzone)   (:,  i)  )
      END DO
      !$ACC UPDATE ASYNC(1) &
      !$ACC   HOST(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)) &
      !$ACC   HOST(dsl4jsb_var_ptr(HYDRO_, ice_rootzone))
      !$ACC WAIT(1)

      ! Make sure the root zone water content is smaller than the maximum possible amount of
      ! liquid water in the root zone.
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) > &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) = MAX(0._wp, &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
      END WHERE
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)) ASYNC(1)

      ! Initialization of the relative root zone moisture and plant available root zone moisture
      ! Note: It is re-calculated every time step in mo_hydro_interface:update_water_stress.
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone) > 0._wp)
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel) = &
          &  dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) / &
          & (dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
      ELSE WHERE
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel) = 0._wp
      END WHERE
      dsl4jsb_var2D_onDomain(HYDRO_, wtr_plant_avail_rel) = dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel)

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone_rel)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_plant_avail_rel))
      !$ACC WAIT(1)

    END IF ! Tile with vegetation
  END SUBROUTINE hydro_init_ic
#endif

  ! ====================================================================================================== !
  !
  !>#### Initialize soil properties
  !>
  !> The soil properties had been set in [[hydro_init_bc]] and [[hydro_init_ic]]  as two
  !> dimensional maps. Here, they are defined for all soil levels.
  !>
  ! Note: Routine [[mo_hydro_interface:update_soil_properties]] does similar calculations at each time
  ! step. We also need the routine here - within the initialization phase - because the soil snow and
  ! energy (SSE) processes make use of the 3d soil properties and the SSE processes are updated before
  ! the soil properties are calculated the first time in [[mo_hydro_interface:update_soil_properties]].
  !
  SUBROUTINE init_soil_properties(tile)

    USE mo_hydro_constants, ONLY: &
      ! Parameters for organic soil component in top layer
      & vol_porosity_org_top, hyd_cond_sat_org_top, bclapp_org_top, matric_pot_org_top, &
      & pore_size_index_org_top, vol_field_cap_org_top, vol_p_wilt_org_top, vol_wres_org_top, &
      ! Parameters for organic soil component in layers below top layer
      & vol_porosity_org_below, hyd_cond_sat_org_below, bclapp_org_below, matric_pot_org_below, &
      & pore_size_index_org_below, vol_field_cap_org_below, vol_p_wilt_org_below, vol_wres_org_below

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_) !< Configuration of the HYDRO process
    dsl4jsb_Def_memory(HYDRO_) !< Memory of the HYDRO process

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & hyd_cond_sat,          & !< Hydraulic conductivity [m s-1]
      & vol_porosity,          & !< Volumetric soil porosity []
      & bclapp,                & !< Clapp and Hornberger exponent b []
      & matric_pot,            & !< Soil matric potential [m]
      & pore_size_index,       & !< Soil pore size index []
      & vol_field_cap,         & !< Volumetric soil field capacity []
      & vol_p_wilt,            & !< Volumetric soil wilting point []
      & vol_wres                 !< Volumetric residual soil water content []
    dsl4jsb_Real3D_onDomain :: &
      & fract_org_sl,          & !< Volumetric fraction of soil layer organic carbon []
      & hyd_cond_sat_sl,       & !< Saturated hydraulic conductivity on soil layers [m s-1]
      & vol_porosity_sl,       & !< Volumetric soil porosity on soil layers []
      & bclapp_sl,             & !< Clapp and Hornberger exponent b on soil layers []
      & matric_pot_sl,         & !< Matric potential on soil layers [m]
      & pore_size_index_sl,    & !< Pore size index on soil layers []
      & vol_field_cap_sl,      & !< Volumetric soil field capacity on soil layers []
      & vol_p_wilt_sl,         & !< Volumetric soil wilting point on soil layers []
      & vol_wres_sl              !< Volumetric residual soil water content on soil layers []

    ! Locally allocated vectors
    TYPE(t_jsb_model), POINTER :: model    !< This instence of ICON-Land
    TYPE(t_jsb_vgrid), POINTER :: soil_w   !< Vertical grid for the hydrology

    REAL(wp) ::                &
      & hyd_cond_sat_org,      & !< Hydraulic conductivity of organic matter [m s-1]
      & vol_porosity_org,      & !< Volumetric soil porosity of organic matter []
      & bclapp_org,            & !< Clapp and Hornberger exponent b for organic matter []
      & matric_pot_org,        & !< Soil matric potential of organic matter [m]
      & pore_size_index_org,   & !< Soil pore size index of organic matter []
      & vol_field_cap_org,     & !< Volumetric soil field capacity of organic matter []
      & vol_p_wilt_org,        & !< Volumetric soil wilting point of organic matter []
      & vol_wres_org             !< Volumetric residual soil water content of organic matter []

    INTEGER  ::                &
      & nsoil,                 & !< Number of soil layers
      & i                        !< Looping index

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    ! Do nothing if hydrology process is not active on this tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    ! Do nothing if this is a lake tile
    IF (tile%is_lake) RETURN

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, hyd_cond_sat)       ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_porosity)       ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, bclapp)             ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, matric_pot)         ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, pore_size_index)    ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_field_cap)      ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_p_wilt)         ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_wres)           ! in
    IF (dsl4jsb_Config(HYDRO_)%l_organic)   &
     &   dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)  ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, hyd_cond_sat_sl)    ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)    ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)          ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, matric_pot_sl)      ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, pore_size_index_sl) ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_p_wilt_sl)      ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_wres_sl)        ! out

    ! -------------------------------------------------------------------------------------------------- !
    !
    !> - Distribute 2d input parameters over all soil layers
    !
    hyd_cond_sat_sl   (:,:,:) = SPREAD(hyd_cond_sat   (:,:), DIM=2, ncopies=nsoil)
    vol_porosity_sl   (:,:,:) = SPREAD(vol_porosity   (:,:), DIM=2, ncopies=nsoil)
    bclapp_sl         (:,:,:) = SPREAD(bclapp         (:,:), DIM=2, ncopies=nsoil)
    matric_pot_sl     (:,:,:) = SPREAD(matric_pot     (:,:), DIM=2, ncopies=nsoil)
    pore_size_index_sl(:,:,:) = SPREAD(pore_size_index(:,:), DIM=2, ncopies=nsoil)
    vol_field_cap_sl  (:,:,:) = SPREAD(vol_field_cap  (:,:), DIM=2, ncopies=nsoil)
    vol_p_wilt_sl     (:,:,:) = SPREAD(vol_p_wilt     (:,:), DIM=2, ncopies=nsoil)
    vol_wres_sl       (:,:,:) = SPREAD(vol_wres       (:,:), DIM=2, ncopies=nsoil)

    ! -------------------------------------------------------------------------------------------------- !
    !
    !> - Update soil parameters with organic fractions: Use the parameters for upper soil organic matter
    !>   in the top soil layer and the parameters for deeper soils in the layers below.
    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
      DO i=1,nsoil
        IF (i == 1) THEN
          hyd_cond_sat_org    = hyd_cond_sat_org_top
          vol_porosity_org    = vol_porosity_org_top
          bclapp_org          = bclapp_org_top
          matric_pot_org      = matric_pot_org_top
          pore_size_index_org = pore_size_index_org_top
          vol_field_cap_org   = vol_field_cap_org_top
          vol_p_wilt_org      = vol_p_wilt_org_top
          vol_wres_org        = vol_wres_org_top
        ELSE
          hyd_cond_sat_org    = hyd_cond_sat_org_below
          vol_porosity_org    = vol_porosity_org_below
          bclapp_org          = bclapp_org_below
          matric_pot_org      = matric_pot_org_below
          pore_size_index_org = pore_size_index_org_below
          vol_field_cap_org   = vol_field_cap_org_below
          vol_p_wilt_org      = vol_p_wilt_org_below
          vol_wres_org        = vol_wres_org_below
        END IF

        !>     The actual soil properties are calculated as weighted means from the properties of mineral
        !> and organic soils.
        ! Attention: hyd_cond_sat_sl is calculated in a more complex way in update_soil_properties!
        hyd_cond_sat_sl   (:,i,:) = (1 - fract_org_sl(:,i,:)) * hyd_cond_sat_sl   (:,i,:) &
          & + fract_org_sl(:,i,:) * hyd_cond_sat_org
        vol_porosity_sl   (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_porosity_sl   (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_porosity_org
        bclapp_sl         (:,i,:) = (1 - fract_org_sl(:,i,:)) * bclapp_sl         (:,i,:) &
          & + fract_org_sl(:,i,:) * bclapp_org
        matric_pot_sl     (:,i,:) = (1 - fract_org_sl(:,i,:)) * matric_pot_sl     (:,i,:) &
          & + fract_org_sl(:,i,:) * matric_pot_org
        pore_size_index_sl(:,i,:) = (1 - fract_org_sl(:,i,:)) * pore_size_index_sl(:,i,:) &
          & + fract_org_sl(:,i,:) * pore_size_index_org
        vol_field_cap_sl  (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_field_cap_sl  (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_field_cap_org
        vol_p_wilt_sl     (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_p_wilt_sl     (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_p_wilt_org
        vol_wres_sl       (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_wres_sl       (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_wres_org
      END DO
    END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl) &
    !$ACC   DEVICE(matric_pot_sl, pore_size_index_sl, vol_field_cap_sl) &
    !$ACC   DEVICE(vol_p_wilt_sl, vol_wres_sl)
    !$ACC WAIT(1)

  END SUBROUTINE init_soil_properties

  ! ====================================================================================================== !
  !
  !>#### Sanitize HYDRO variables
  !>
  !> Loading a model state from a file with lower than double precision can lead to model aborts because
  !> the state violates model constraints due to rounding errors in soil water content and root depth.
  !> This routine sanitizes the current state by limiting soil water and ice to saturation capacity, and
  !> by recomputing the root depth.
  !>
  SUBROUTINE hydro_sanitize_state(tile)

    USE mo_util, ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_sanitize_state'

    REAL(wp), PARAMETER :: warn_level_wpi_soil_excess = 1.05_wp
                           !< If water content exceeds saturation by more than this factor prior
                           !< to sanitization, a warning is issued.
    REAL(wp), PARAMETER :: zfcmin = 1.e-10_wp   !TODO: Remove

    INTEGER :: maxloc_wpi_soil_excess(3)  !< Location of the maximum excess soil water
    TYPE(t_jsb_vgrid), POINTER :: soil_w  !< Vertical grid of the hydrology

    dsl4jsb_Def_memory(HYDRO_)            !< Memory of the HYDRO process

    dsl4jsb_Real2D_onDomain :: &
      & root_depth,            &          !< Rooting depth [m]
      & soil_depth                        !< Soil depth until bedrock [m]

    dsl4jsb_Real3D_onDomain :: &
      & root_depth_sl,         &          !< Rooting depth within the layer [m]
      & soil_depth_sl,         &          !< Soil depth within the layer (until bedrock) [m]
      & vol_porosity_sl,       &          !< Volumetric soil porosity []
      & wtr_soil_sat_sl,       &          !< Amount of water within the layer at saturation [m water equivalent]
      & wtr_soil_sl,           &          !< Actual amount of water within the layer [m water equivalent]
      & ice_soil_sl                       !< Actual amount of ice within the layer [m ice]

    ! Do nothing if the hydrology is not calculated on the current tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, root_depth)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, soil_depth)

    dsl4jsb_Get_var3D_onDomain(HYDRO_, root_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sat_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, ice_soil_sl)

    IF (tile%contains_soil) THEN

      !> A warning is issued in case the current soil water anywhere exceeds saturation capacity by more than
      !> a certain factor.
      IF (ANY(wtr_soil_sl(:,:,:) > warn_level_wpi_soil_excess * wtr_soil_sat_sl(:,:,:))) THEN
        maxloc_wpi_soil_excess(:) = MAXLOC(wtr_soil_sl(:,:,:) - wtr_soil_sat_sl(:,:,:))
        WRITE (message_text,*) 'liquid water content above saturation capacity at ', maxloc_wpi_soil_excess(:), &
          & ': wtr_soil_sl = ',                                                                                 &
          & wtr_soil_sl    (maxloc_wpi_soil_excess(1), maxloc_wpi_soil_excess(2), maxloc_wpi_soil_excess(3)),   &
          & ', wtr_soil_sat_sl = ',                                                                             &
          & wtr_soil_sat_sl(maxloc_wpi_soil_excess(1), maxloc_wpi_soil_excess(2), maxloc_wpi_soil_excess(3))
        CALL warning(routine, message_text)
      END IF

      ! Limit soil water and soil ice to saturation amounts.
      IF (ANY((wtr_soil_sl(:,:,:) + ice_soil_sl(:,:,:))                &
        &  > MAX(vol_porosity_sl (:,:,:) * soil_depth_sl(:,:,:), 0._wp))) THEN
        ! Limit liquid soil water to saturation amount
        wtr_soil_sl(:,:,:) = MIN(wtr_soil_sl(:,:,:), wtr_soil_sat_sl(:,:,:))

        ! Limit ice content such that soil water + soil ice <= saturation.
        ice_soil_sl(:,:,:) = MAX(0._wp, MIN(ice_soil_sl(:,:,:), &
          &                             MAX(vol_porosity_sl (:,:,:) * soil_depth_sl(:,:,:), 0._wp) - wtr_soil_sl(:,:,:)))

        !$ACC UPDATE DEVICE(wtr_soil_sl, ice_soil_sl) ASYNC(1)
      END IF
    END IF

    IF (tile%contains_vegetation) THEN
      ! Due to rounding, the root depth in a layer might differ from the soil depth even though the layer is contained
      ! in the root zone. The might not add up to the total root depth for the same reason. Fix by recomputing.
      soil_w => Get_vgrid('soil_depth_water')
      root_depth(:,:) = MIN(root_depth(:,:), soil_depth(:,:))
      root_depth_sl(:,:,:) = soil_depth_to_layers_2d(root_depth(:,:), soil_w%dz(:))

      !$ACC UPDATE DEVICE(root_depth, root_depth_sl) ASYNC(1)
    END IF
    !$ACC WAIT(1)

  END SUBROUTINE hydro_sanitize_state
  !
  ! ====================================================================================================== !

#endif
END MODULE mo_hydro_init
