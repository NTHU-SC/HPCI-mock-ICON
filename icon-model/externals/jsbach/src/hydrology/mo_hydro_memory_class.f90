!> Contains the memory class for the hydrology process.
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
!>#### Memory definition of the hydrology process
!>
!> The module defines the data structure [[t_hydro_memory]] that contains the hydrology
!> memory. Each tile of the HYDRO process has its own version of this memory.
!> Variables added to the structure are easily accessible in the process routines of the
!> HYDRO and other processes, e.g. with 'dsl4jsb_Get_var2D_onChunk(HYDRO_, water_stress)'
!> for using variable 'water_stress'.
!>
!> With the [[t_jsb_memory:Add_var]] calls, memory is allocated for
!> the global variable arrays on CPUs. It is also ensured that the correct variable
!> patches are available on different MPI processes, and the variable is made available
!> on GPUs, in case the model is compiled to be run on GPUs.
!> Besides, calling [[t_jsb_memory:Add_var]] enables writing selected variables to output
!> files and adds the necessary variables to the list of restart file variables.
!>
!> Subroutine [[Init_hydro_memory]] is a type-bound procedure of [[t_hydro_memory]]. It is
!> is called within the initialization phase for all tiles of the HYDRO process.
!>
MODULE mo_hydro_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_control,      ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,  ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE, LAKE_TYPE, GLACIER_TYPE
  USE mo_jsb_varlist,      ONLY: NONE, BASIC, MEDIUM, FULL

  ! Use of processes in this module
  dsl4jsb_Use_processes HYDRO_, SSE_, SEB_
#ifndef __QUINCY_STANDALONE__
  dsl4jsb_Use_processes ASSIMI_
#endif

  ! Use process configurations
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(SEB_)
#ifndef __QUINCY_STANDALONE__
  dsl4jsb_Use_config(ASSIMI_)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 130   !< Maximum number of 'Add_var' calls in this module

  !> Memory of the HYDRO process:
  !> _Note for the ford documentation_: The table lists the general jsbach memory variables
  !>   ([[t_jsb_memory]]) as well as its extension by the HYDRO process variables.
  TYPE, EXTENDS(t_jsb_memory) :: t_hydro_memory

    TYPE(t_jsb_var_real2d) :: &
      & fract_snow,           & !< Snow area fraction (not including snow on canopy)     []
      & weq_snow,             & !< Snow amount (including snow on canopy)                [m water equivalent]
      & evapo_snow,           & !< Evaporation from snow                                 [kg m-2 s-1]
      & snowmelt,             & !< Snow melt                                             [kg m-2 s-1]
      & le_pc_remain,         & !< Remaining latent energy after phase change            [J m-2]
      & evapopot,             & !< Potential evaporation                                 [kg m-2 s-1]
      & evapotrans,           & !< Evapotranspiration                                    [kg m-2 s-1]
      & fract_pond,           & !< Inudated area fraction (area of temporary ponds)      []
      & fract_pond_max,       & !< Maximum possible inudated tile fraction               []
      & weq_pond,             & !< Water content of pond reservoir                       [m water equivalent]
      & weq_pond_max,         & !< Maximum water content of pond reservoir               [m water equivalent]
      & wtr_pond,             & !< Liquid water content of pond reservoir                [m]
      & pond_melt,            & !< Melting of ice in pond reservoir                      [kg m-2 s-1]
      & pond_freeze,          & !< Freezing of liquid water in pond reservoir            [kg m-2 s-1]
      & wtr_pond_net_flx,     & !< Net water flux into pond reservoir                    [kg m-2 s-1]
      & ice_pond,             & !< Frozen water content of pond reservoir                [m water equivalent]
      & weq_fluxes,           & !< Sum of all land water fluxes; needed to test water conservation   [m3/s]
      & weq_land,             & !< Total amount of land water and ice                    [m3]
      & weq_balance_err,      & !< Land water balance error within a time step           [m3/(time step)]
      & weq_balance_err_count   !< Number of time steps with water balance error         []

    TYPE(t_jsb_var_real2d) :: &
      & infilt,               & !< Infiltration                                          [kg m-2 s-1]
      & runoff,               & !< Surface runoff                                        [kg m-2 s-1]
      & runoff_horton,        & !< Horton component of surface runoff (infiltration excess)   [kg m-2 s-1] TODO: Reference?
      & runoff_dunne,         & !< Dunne component of surface runoff (saturation excess)      [kg m-2 s-1] TODO: Reference?
      & drainage,             & !< Total (subsurface) drainage                           [kg m-2 s-1]
      & discharge,            & !< Discharge (local)                                     [m3 s-1]
      & discharge_ocean,      & !< Discharge to the ocean                                [m3 s-1]
      & internal_drain          !< Internal drainage                                     [kg m-2 s-1]

    TYPE(t_jsb_var_real2d) :: &
      & elevation,            & !< Topographic height                                    [m]
      & oro_stddev,           & !< Standard deviation of the orography                   [m]
      & steepness               !< Parameter defining subgrid scale slope                []

    ! Additional variables for land type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_soil,      & !< Snow fraction on soil or glacier                      []
      & weq_snow_soil,        & !< Amount of snow on soil                                [m water equivalent]
      & snow_soil_dens,       & !< Density of snow on soil                               [kg m-3]
      & evapotrans_lnd,       & !< Evapotranspiration of vegetation and bare soil        [kg m-2 s-1]
      & q_snocpymlt             !< Required energy flux for melting snow on canopy       [W m-2]

    ! Additional variables for soil
    TYPE(t_jsb_var_real3d) :: &
      & soil_depth_sl,        & !< Thickness of each soil layer until bedrock            [m]
      & fract_org_sl            !< Fractions of organic material in soil layers          []

    TYPE(t_jsb_var_real3d) :: &
      & soil_lay_width_sl, &        !< thickness of each soil layer that can be saturated with water, identical with soil_depth_sl (until bedrock) [m]
      & soil_lay_depth_center_sl, & !< depth at the center of each soil layer that can be saturated with water (until bedrock) [m]
      & soil_lay_depth_ubound_sl, & !< depth at the upper bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]
      & soil_lay_depth_lbound_sl    !< depth at the lower bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]

    TYPE(t_jsb_var_real2d) :: &
      & soil_depth,           & !< Soil depth until bedrock                              [m]
      & num_sl_above_bedrock, & !< number of soil layers (thickness > eps8) above bedrock[]
      & wpi_rootzone_max ,    & !< Maximum amount of water or ice in the root zone       [m]
      & vol_field_cap,        & !< Volumetric soil field capacity                        [m/m]
      & vol_p_wilt,           & !< Volumetric permanent wilting point                    [m/m]
      & vol_porosity,         & !< Volumetric soil porosity                              [m/m]
      & vol_wres,             & !< Volumetric residual water content                     [m/m]
      & pore_size_index,      & !< Soil pore size distribution index                     []
      & bclapp,               & !< Exponent B in Clapp and Hornberger                    []
      & matric_pot,           & !< Soil matric potential                                 [m]
      & hyd_cond_sat            !< Saturated hydraulic conductivity                      [m/s]

    TYPE(t_jsb_var_real2d) :: &
      & fract_wet,            & !< Wet tile fraction (skin and ponds; soil and canopy)   []
      & fract_skin,           & !< Wet skin fraction (w/o ponds; soil and canopy)        []
      & wtr_skin,             & !< Water content in skin reservoir (soil and canopy)     [m]
      & snow_accum,           & !< Change in snow amount on non-glacier/non-lake land
                                !< tiles                                   [m water equivalent / (time step)]
      & evapotrans_soil,      & !< Evapotranspiration from soil (w/o snow, skin or pond
                                !< evaporation)                                          [kg m-2 s-1]
      & evapo_pond,           & !< Evaporation from pond reservoir                       [kg m-2 s-1]
      & evapo_skin,           & !< Evaporation from skin reservoir                       [kg m-2 s-1]
      & evapo_deficit,        & !< Evaporation deficit flux due to inconsistent
                                !< treatment of snow evaporation           [m water equivalent / (time step)]
      & water_to_soil,        & !< Water available for infiltration into the soil        [m water equivalent]
      & wtr_soilhyd_res,      & !< Residual of vertical soil water transport scheme      [m3 / (time step)]
      & tpe_overflow,         & !< Water content of terraplanet reservoir for soil water
                                !< overflow                                              [m]
      & wtr_rootzone,         & !< Liquid water content in the root zone                 [m]
      & wtr_rootzone_rel,     & !< Water content of the root zone relative to the maximum
                                !< possible water content                                []
      & wtr_latflow_res_srf,  & !< Water content of intermediary reservoir representing
                                !< lateral flow of surface runoff (only with
                                !< [[t_hydro_config:hydro_scale=Uniform]])               [m]
      & wtr_latflow_srf         !< Outflow from the intermediary reservoir representing
                                !< lateral flow of surface runoff                        [kg m-2 s-1]

    TYPE(t_jsb_var_real3d) :: &
      & vol_field_cap_sl,     & !< Volumetric soil field capacity                        [m/m]
      & vol_p_wilt_sl,        & !< Volumetric permanent wilting point                    [m/m]
      & vol_porosity_sl,      & !< Volumetric soil porosity                              [m/m]
      & vol_wres_sl,          & !< Volumetric residual water content                     [m/m]
      & pore_size_index_sl,   & !< Soil pore size distribution index                     []
      & bclapp_sl,            & !< Exponent B in Clapp and Hornberger                    []
      & matric_pot_sl,        & !< Soil matric potential                                 [m]
      & wtr_soil_pot_sl,      & !< Soil water potential per layer                        [MPa]
      & hyd_cond_sat_sl,      & !< Saturated hydraulic conductivity                      [m/s]
      & wtr_soil_sl,          & !< Water content in soil layers                          [m]
      & ice_soil_sl,          & !< Ice content in soil layers                            [m]
      & vol_weq_soil_sl,      & !< Volumetric water content (liquid+ice) in soil
                                !< layers (Water content below bedrock depth is zero.)   [m3/m3]
      & wtr_freeze_sl,        & !< Freezing water in soil layers                         [kg m-2 s-1]
      & ice_melt_sl,          & !< Melting ice in soil layers                            [kg m-2 s-1]
      & drainage_sl,          & !< Subsurface drainage (horizontal) from soil layers     [kg m-2 s-1]
      & wtr_transp_down_sl,   & !< Vertical water transport into the below soil layer    [kg m-2 s-1]
                                !< negative values refer to upwards transport accordingly
      & wtr_soil_sat_sl,      & !< Soil water content at saturation (calculated from
                                !< soil porosity, reduced by ice)                        [m]
      & wtr_soil_fc_sl,       & !< Water content at field capacity (reduced by ice)      [m]
      & wtr_soil_pwp_sl,      & !< Water content at perm. wilting point (reduced by ice) [m]
      & wtr_soil_res_sl,      & !< Residual soil water content (reduced by ice)          [m]
      & wtr_soil_pot_scool_sl,& !< Potential amount of supercooled water                 [m]
      & wtr_latflow_res_sl,   & !< Water content of intermediary reservoir representing
                                !< lateral flow of subsurface runoff (only with
                                !< [[t_hydro_config:hydro_scale=Uniform]])               [m]
      & wtr_latflow_sl          !< Outflow from the intermediary reservoir representing
                                !< lateral flow of subsurface runoff                     [kg m-2 s-1]

    TYPE(t_jsb_var_real2d) :: &
      & wtr_plant_avail_rel, &  !< relative plant available water (difference between saturation and PWP in root zone) [unitless]
      & w_soil_root_fc,   &     !< Water content at field capacity in root zone of the soil           [m]
      & w_soil_root_pwp         !< Water content at permanent wilting point in root zone of the soil  [m]

    TYPE(t_jsb_var_real3d) :: &
      & frac_w_lat_loss_sl, &   !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep) [unitless]
      & frac_wtr_transp_down_sl !< fraction of material transferred (vertical) to below layer                                  [unitless]
                                !< negative values refer to upwards transport accordingly

    ! Additional variables for PFT lct_type
    TYPE(t_jsb_var_real3d) ::   &
      & root_depth_sl             !< Rooted depth per soil layer (until rooting depth)   [m]
    TYPE(t_jsb_var_real2d) ::   &
      & root_depth,             & !< Rooting depth                                       [m]
      & fract_snow_can,         & !< Snow fraction on canopy                             []
      & weq_snow_can,           & !< Snow amount on canopy                               [m water equivalent]
      & ice_rootzone,           & !< Ice content of the root zone (not water equivalent) [m]
      & wtr_rootzone_scool_pot, & !< Potential amount of supercooled water in the root   [m]
                                  !< zone
      & wtr_rootzone_scool_act, & !< Supercooled water content in root zone              [m]
      & wtr_rootzone_avail,     & !< Plant available water content in root zone          [m]
      & wtr_rootzone_avail_max, & !< Maximum plant available water content in root zone  [m]
      & water_stress,           & !< Water stress factor used in assimilation and
                                  !< canopy conductance (1: none, 0: infinite stress)    []
      & canopy_cond_unlimited,  & !< Canopy conductance without water limitation         [m/s]
      & canopy_cond_limited,    & !< Canopy conductance accounting for water stress      [m/s]
      & transpiration             !< Transpiration                                       [kg m-2 s-1]

    ! Additional variables for GLACIER lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_glac,      & !< Snow fraction on glacier                              []  TODO: remove
      & weq_glac,             & !< Glacier depth (snow is considered as glacier)         [m water equivalent]
      & runoff_glac             !< Runoff from glacier (rain + ice melt, no calving)     [m water equivalent]

    ! Additional variables if no separate glacier lct_type is used and glaciers are treated as part of SOIL lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_glac              !< Glacier fraction                                      []

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & evapo_wtr,            & !< Evaporation from lake water                           [kg m-2 s-1]
      & evapo_ice,            & !< Evaporation from lake ice                             [kg m-2 s-1]
      & fract_snow_lice,      & !< Snow fraction on lake ice                             []
      & weq_snow_lice           !< Water content of snow on lake ice                     [m water equivalent]

    ! Diagnostic global land means/sums e.g. for monitoring
    TYPE(t_jsb_var_real1d) :: &
      trans_gmean,            & !< Global mean transpiration                [kg m-2 s-1]
      evapotrans_gmean,       & !< Global land mean evapotranspiration      [kg m-2 s-1]
      weq_land_gsum,          & !< Global land water and ice content        [km3]
      discharge_ocean_gsum,   & !< Global water discharge to the oceans     [Sv]
      wtr_rootzone_rel_gmean, & !< Global mean relative root zone moisture  []
      fract_snow_gsum,        & !< Global snow area on non-glacier land     [Mio m2]
      weq_snow_gsum,          & !< Global snow amount on non-glacier land   [Gt]
      weq_balance_err_gsum      !< Global sum of water balance error        [m3/(time step)]

  CONTAINS
    PROCEDURE :: Init => Init_hydro_memory
  END TYPE t_hydro_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_memory_class'

CONTAINS

  ! ============================================================================================== !
  !>
  !>#### Initialization of the hydrology process memory
  !>
  !> In this subroutine we initialize the memory of the hydrology process.
  !> The subroutine is a type-bound procedure of [[t_hydro_memory]]. It is called within the
  !> initialization phase for all tiles needed by the hydrology process ([[mo_jsb_tile:init_tile]]:
  !> 'CALL this%mem(iproc)%p%Init(TRIM(process_name), ...)').
  !>
  SUBROUTINE Init_hydro_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    & TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    USE mo_jsb_physical_constants, ONLY: &
      & dens_snow,                       & ! Density of snow (l_dynsnow=.FALSE.)
      & dens_snow_min                      ! Density of fresh snow (l_dynsnow=.TRUE.)

    USE mo_hydro_constants,        ONLY: Semi_Distributed_, Uniform_

    CLASS(t_hydro_memory), INTENT(inout), TARGET :: mem   !< Memory of the hydrology process
    CHARACTER(len=*),      INTENT(in)    :: prefix        !< Prefix for variable names in output/restart files,
                                                          !< i.e. process name
    CHARACTER(len=*),      INTENT(in)    :: suffix        !< Suffix for variable names in output/restart files,
                                                          !< i.e. corresponding tile name
    INTEGER,               INTENT(in)    :: lct_ids(:)    !< Land cover type IDs
    INTEGER,               INTENT(in)    :: lib_id        !< ID of the tile's primary land cover type (e.g. LAND_TYPE)
    INTEGER,               INTENT(in)    :: model_id      !< ID of the current model instance

    dsl4jsb_Def_config(HYDRO_)              !< Configurable settings of the hydrology process
    dsl4jsb_Def_config(SSE_)                !< Configurable settings of the soil and snow energy process
    dsl4jsb_Def_config(SEB_)                !< Configurable settings of surface energy balance process
#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_config(ASSIMI_)             !< Configurable settings of the assimilation process
#endif

    TYPE(t_jsb_model), POINTER :: model         !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER :: hgrid         !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, &    !< Vertical grid for surface variables (vertical dimension = 1)
      &                           soil_w        !< Vertical grid for soil hydrology

    INTEGER :: table                            !< Grib code table ID
    TYPE(t_grib2) :: grib2_desc                 !< Grib2 related metadata

    REAL(wp) :: ini_snow_dens                   !< Initial value for snow density
    LOGICAL  :: l_ponds                         !< Abreviation: Local variable for [[t_hydro_config:l_ponds]]
    LOGICAL  :: l_uniform                       !< Abreviation: True if [[t_hydro_config:hydro_scale]]=Uniform
    LOGICAL  :: lrestart_local                  !< add var to restart file, model specific (quincy/jsbach)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hydro_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    soil_w       => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(SEB_)
#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Get_config(ASSIMI_)
#endif

    ! Initial snow density depends on whether or not snow density is to be calculated dynamically
    IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
      ini_snow_dens = dens_snow_min
    ELSE
      ini_snow_dens = dens_snow
    END IF

    ! Define logicals to be used below to exclude some variables from restart file
    IF (dsl4jsb_Config(HYDRO_)%hydro_scale == Uniform_) THEN
      l_uniform = .TRUE.
    ELSE
      l_uniform = .FALSE.
    END IF
    l_ponds = dsl4jsb_Config(HYDRO_)%l_ponds

    ! set local lrestart switch to TRUE for the MODEL_QUINCY
    lrestart_local = .FALSE.
    IF (model%config%model_scheme == MODEL_QUINCY) THEN
      lrestart_local = .TRUE.
    END IF

    ! Variables defined on all tiles
    ! -------------------------------

    !> Only a few selected variables are available in grib format. For these variables official WMO grib codes
    !> are defined. Note, that the code triplet '255, 255, 255' used for most jsbach output variables is rather
    !> a placeholder.

    !TODO: Change unit for unitless variables '-' -> ''
    ! Define grib codes for snow fraction if we are on the box tile, i.e. the name of the tile "owning" this
    ! memory is 'box'.
    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,202, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF

    CALL mem%Add_var('fract_snow', mem%fract_snow,                             &
      & hgrid, surface,                                                        &
      & t_cf('fract_snow', '-', 'Snow area fraction'),                         &
      & t_grib1(table, 255, grib_bits), grib2_desc,                            &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,212, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF
    CALL mem%Add_var( 'weq_snow', mem%weq_snow,                                &
      & hgrid, surface,                                                        &
      & t_cf('weq_snow', 'm (water equivalent)', 'Snow amount'),               &
      & t_grib1(table, 255, grib_bits), grib2_desc,                            &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapotrans', mem%evapotrans,                                      &
      & hgrid, surface,                                                                  &
      & t_cf('evapotrans', 'kg m-2 s-1', 'Evapotranspiration, including sublimation'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapopot', mem%evapopot,                                          &
      & hgrid, surface,                                                                  &
      & t_cf('evapopot', 'kg m-2 s-1', 'Potential evaporation (including sublimation)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapo_snow', mem%evapo_snow,                                  &
      & hgrid, surface,                                                              &
      & t_cf('evapo_snow', 'kg m-2 s-1', 'Evaporation from snow'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'snowmelt', mem%snowmelt,                                      &
      & hgrid, surface,                                                              &
      & t_cf('snowmelt', 'kg m-2 s-1', 'Snow melt'),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'infilt', mem%infilt,                                          &
      & hgrid, surface,                                                              &
      & t_cf('infilt', 'kg m-2 s-1', 'Water infiltrating into the soil'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE., output_level=BASIC,                                        &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'runoff', mem%runoff,                                          &
      & hgrid, surface,                                                              &
      & t_cf('runoff', 'kg m-2 s-1', 'Total surface runoff'),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & loutput=.TRUE., output_level=BASIC,                                          &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'runoff_horton', mem%runoff_horton,                            &
      & hgrid, surface,                                                              &
      & t_cf('runoff_horton', 'kg m-2 s-1', 'Horton component of surface runoff'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'runoff_dunne', mem%runoff_dunne,                              &
      & hgrid, surface,                                                              &
      & t_cf('runoff_dunne', 'kg m-2 s-1', 'Dunne component of surface runoff'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'drainage', mem%drainage,                                      &
      & hgrid, surface,                                                              &
      & t_cf('drainage', 'kg m-2 s-1', 'Total drainage'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & loutput=.TRUE., output_level=BASIC,                                          &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'wtr_transp_down_sl', mem%wtr_transp_down_sl,                  &
      & hgrid, soil_w,                                                               &
      & t_cf('wtr_transp_down_sl', 'kg m-2 s-1',                                     &
      &      'Downwards transport of water into next deeper soil layer'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'discharge', mem%discharge,                                    &
      & hgrid, surface,                                                              &
      & t_cf('discharge', 'm3 s-1', 'Local discharge'),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, output_level=BASIC,                                          &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'discharge_ocean', mem%discharge_ocean,                        &
      & hgrid, surface,                                                              &
      & t_cf('discharge_ocean', 'm3 s-1', 'Discharge to the ocean'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, output_level=BASIC,                                          &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'internal_drain', mem%internal_drain,                          &
      & hgrid, surface,                                                              &
      & t_cf('internal_drain', 'kg m-2 s-1', 'Internal drainage'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'weq_fluxes', mem%weq_fluxes,                                  &
      & hgrid, surface,                                                              &
      & t_cf('weq_fluxes', 'm3 s-1',                                                 &
      &      'All land water fluxes (for water balance check)'),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'weq_land', mem%weq_land,                                      &
      & hgrid, surface,                                                              &
      & t_cf('weq_land', 'm3', 'Total land water and ice content'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'weq_balance_err', mem%weq_balance_err,                        &
      & hgrid, surface,                                                              &
      & t_cf('weq_balance_err', 'm3/(time step)', 'Land water balance error'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & loutput=.TRUE., output_level=FULL,                                           &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    !TODO: Change cf names: t_cf('weq_balance_err_count', '',                                            &
    !      'Number of time steps with water balance error since (re)start')
    CALL mem%Add_var( 'weq_balance_err_count', mem%weq_balance_err_count,            &
      & hgrid, surface,                                                              &
      & t_cf('weq_balance_err_count', 'time steps',                                  &
      &      'Amount of time steps with water balance errors'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    IF (jsbach_runs_standalone()) THEN
      CALL mem%Add_var( 'elevation', mem%elevation,                                  &
        & hgrid, surface,                                                            &
        & t_cf('elevation', 'm', 'Mean orography'),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE.,                                                          &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )
    END IF

    CALL mem%Add_var( 'oro_stddev', mem%oro_stddev,                                  &
      & hgrid, surface,                                                              &
      & t_cf('oro_stddev', 'm', 'Orographic standard deviation'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'steepness', mem%steepness,                                    &
      & hgrid, surface,                                                              &
      & t_cf('steepness', '/', 'subgrid slope distribution shape parameter'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    IF (model%config%use_tmx) THEN
      CALL mem%Add_var( 'heating_snow_cpy_melt', mem%q_snocpymlt,                       &
        & hgrid, surface,                                                               &
        & t_cf('heating_snow_cpy_melt', 'W m-2', 'Heating due to snow melt on canopy'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & prefix, suffix,                                                               &
        & lrestart=.FALSE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    END IF

    ! Additional variables if the tile contains land, i.e. the tile or sub-tiles are land tiles
    ! ------------------------------------------------
    IF ( (     One_of(VEG_TYPE,     lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE,    lct_ids(:)) > 0 &
      &   .OR. One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'fract_snow_soil', mem%fract_snow_soil,                      &
        & hgrid, surface,                                                            &
        & t_cf('fract_snow_soil', '-', 'Snow fraction on soil'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & output_level=BASIC,                                                        &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'weq_snow_soil', mem%weq_snow_soil,                          &
        & hgrid, surface,                                                            &
        & t_cf('weq_snow_soil', 'm (water equivalent)', 'Snow amount on soil'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & output_level=BASIC,                                                        &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_soil_dens', mem%snow_soil_dens,                        &
        & hgrid, surface,                                                            &
        & t_cf('snow_soil_dens', 'kg m-3', 'Density of snow on soil'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=ini_snow_dens, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'evapotrans_lnd', mem%evapotrans_lnd,                                                     &
        & hgrid, surface,                                                                                         &
        & t_cf('surface_evapotranspiration_lnd', 'kg m-2 s-1', 'Evapotranspiration from land surface w/o lakes'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.FALSE.,                                                                                       &
        & loutput=.TRUE., output_level=BASIC,                                                                     &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'le_pc_remain', mem%le_pc_remain,                            &
        & hgrid, surface,                                                            &
        & t_cf('le_pc_remain', 'J m-2',                                              &
        &      'Remaining latent energy after near-surface water phase change'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE., output_level=MEDIUM,                                     &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (.NOT. model%config%use_tmx) THEN
        CALL mem%Add_var( 'q_snocpymlt', mem%q_snocpymlt,                                 &
          & hgrid, surface,                                                               &
          & t_cf('heating_snow_cpy_melt', 'W m-2', 'Heating due to snow melt on canopy'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
          & prefix, suffix,                                                               &
          & lrestart=.FALSE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,201, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'fract_wet', mem%fract_wet,                              &
        & hgrid, surface,                                                        &
        & t_cf('fract_wet', '-', 'Wet surface fraction'),                        &
        & t_grib1(table, 255, grib_bits), grib2_desc,                            &
        & prefix, suffix,                                                        &
        & output_level=BASIC,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_skin', mem%fract_skin,                            &
        & hgrid, surface,                                                        &
        & t_cf('fract_skin', '-', 'Wet skin reservoir fraction'),                &
        & t_grib1(table, 255, grib_bits), grib2_desc,                            &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE.,                                                      &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables if the tile contains soil (i.e. non-glacier land)
    ! ------------------------------------------------
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'soil_depth', mem%soil_depth,                          &
        & hgrid, surface,                                                      &
        & t_cf('soil_depth', 'm', 'Depth of soil until bedrock'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

#ifndef __NO_QUINCY__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        CALL mem%Add_var('num_sl_above_bedrock', mem%num_sl_above_bedrock, &
          & hgrid, surface, &
          & t_cf('num_sl_above_bedrock', '[unitless]', 'number of soil layers above bedrock'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = BASIC, &
          & loutput = .TRUE., &
          & lrestart=.TRUE., &
          & initval_r=0.0_wp)
      END SELECT
#endif

      CALL mem%Add_var( 'soil_depth_sl', mem%soil_depth_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('soil_depth_sl', 'm', 'Thickness of each soil layer'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
        CALL mem%Add_var( 'fract_org_sl', mem%fract_org_sl,                          &
          & hgrid, soil_w,                                                           &
          & t_cf('fract_org_sl', '-', 'Fraction of organic material in soil layer'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
          & prefix, suffix,                                                          &
          & lrestart=.TRUE.,                                                         &
          & initval_r=0.0_wp )
      END IF

#ifndef __NO_QUINCY__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        CALL mem%Add_var('soil_lay_width_sl', mem%soil_lay_width_sl, &
          & hgrid, soil_w, &
          & t_cf('soil_lay_width_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = BASIC, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp, isteptype=TSTEP_CONSTANT )

        CALL mem%Add_var('soil_lay_depth_center_sl', mem%soil_lay_depth_center_sl, &
          & hgrid, soil_w, &
          & t_cf('soil_lay_depth_center_sl', 'm', 'Depth at the center of each soil layer, until bedrock'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('soil_lay_depth_ubound_sl', mem%soil_lay_depth_ubound_sl, &
          & hgrid, soil_w, &
          & t_cf('soil_lay_depth_ubound_sl', 'm', 'Depth at the upper boundary of each soil layer, until bedrock'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('soil_lay_depth_lbound_sl', mem%soil_lay_depth_lbound_sl, &
          & hgrid, soil_w, &
          & t_cf('soil_lay_depth_lbound_sl', 'm', 'Depth at the lower boundary of each soil layer, until bedrock'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)
      END SELECT
#endif

      CALL mem%Add_var( 'wpi_rootzone_max', mem%wpi_rootzone_max,                    &
        & hgrid, surface,                                                            &
        & t_cf('wpi_rootzone_max', 'm', 'Maximum amount of water/ice in root zone'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE.,                                                          &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_field_cap', mem%vol_field_cap,                    &
        & hgrid, surface,                                                      &
        & t_cf('vol_field_capacity', 'm/m', 'Volumetric soil field capacity'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_field_cap_sl', mem%vol_field_cap_sl,              &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_field_capy_sl', 'm/m', 'Volumetric soil field capacity'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_p_wilt', mem%vol_p_wilt,                          &
        & hgrid, surface,                                                      &
        & t_cf('vol_p_wilt', 'm/m', 'Volumetric permanent wilting point'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_p_wilt_sl', mem%vol_p_wilt_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_p_wilt_sl', 'm/m', 'Volumetric permanent wilting point'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_porosity', mem%vol_porosity,                      &
        & hgrid, surface,                                                      &
        & t_cf('vol_porosity', 'm/m', 'Volumetric soil porosity'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_porosity_sl', mem%vol_porosity_sl,                &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_porosity_sl', 'm/m', 'Volumetric soil porosity'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_wres', mem%vol_wres,                              &
        & hgrid, surface,                                                      &
        & t_cf('vol_wres', 'm/m', 'Volumetric residual water content'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_wres_sl', mem%vol_wres_sl,                        &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_wres_sl', 'm/m', 'Volumetric residual water content'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'pore_size_index', mem%pore_size_index,                &
        & hgrid, surface,                                                      &
        & t_cf('pore_size_index', '', 'Soil pore size distribution index'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'pore_size_index_sl', mem%pore_size_index_sl,          &
        & hgrid, soil_w,                                                       &
        & t_cf('pore_size_index_sl', '', 'Soil pore size distribution index'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'bclapp', mem%bclapp,                                  &
        & hgrid, surface,                                                      &
        & t_cf('bclapp', '', 'Clapp and Hornberger (1978) exponent b'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'bclapp_sl', mem%bclapp_sl,                            &
        & hgrid, soil_w,                                                       &
        & t_cf('bclapp_sl', '', 'Clapp and Hornberger (1978) exponent b'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'matric_pot', mem%matric_pot,                          &
        & hgrid, surface,                                                      &
        & t_cf('matric_pot', 'm', 'Soil saturated matric potential'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'matric_pot_sl', mem%matric_pot_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('matric_pot_sl', 'm', 'Soil saturated matric potential'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart = lrestart_local,                                           &
        & initval_r=0.0_wp )

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        CALL mem%Add_var( 'wtr_soil_pot_sl', mem%wtr_soil_pot_sl,                &
          & hgrid, soil_w,                                                       &
          & t_cf('wtr_soil_pot_sl', 'MPa', 'Soil water potential per layer'),    &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & output_level=BASIC,                                                  &
          & initval_r=0.0_wp)
      END SELECT

      CALL mem%Add_var( 'hyd_cond_sat', mem%hyd_cond_sat,                           &
        & hgrid, surface,                                                           &
        & t_cf('hyd_cond_sat', 'm s-1', 'Hydraulic conductivity at saturation'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE.,                                                         &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'hyd_cond_sat_sl', mem%hyd_cond_sat_sl,                     &
        & hgrid, soil_w,                                                            &
        & t_cf('hyd_cond_sat_sl', 'm s-1', 'Hydraulic conductivity at saturation'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE.,                                                         &
        & initval_r=0.0_wp )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,211, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'wtr_skin', mem%wtr_skin,                                   &
        & hgrid, surface,                                                           &
        & t_cf('wtr_skin', 'm', 'Water content of soil and canopy skin reservoir'), &
        & t_grib1(table, 255, grib_bits), grib2_desc,                               &
        & prefix, suffix,                                                           &
        & output_level=BASIC,                                                       &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_accum', mem%snow_accum,                               &
        & hgrid, surface,                                                           &
        & t_cf('snow_accum', 'm (water equivalent)',                                &
        &      'Snow budget change within the time step'),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'evapotrans_soil', mem%evapotrans_soil,                &
        & hgrid, surface,                                                      &
        & t_cf('evapotrans_soil', 'kg m-2 s-1',                                &
        &      'Evapotranspiration from soil (w/o snow and skin res.)'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'evapo_skin', mem%evapo_skin,                          &
        & hgrid, surface,                                                      &
        & t_cf('evapo_skin', 'kg m-2 s-1', 'Evaporation from skin reservoir'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'evapo_deficit', mem%evapo_deficit,                    &
        & hgrid, surface,                                                      &
        & t_cf('evapo_deficit', 'm', 'Water which evaporated from other sources&
              &/soillayers then intented'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'water_to_soil', mem%water_to_soil,                    &
        & hgrid, surface,                                                      &
        & t_cf('water_to_soil', 'm', 'Amount of water entering the soil'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Overflow pool is only used in the terra planet setup (TPE closed).
      CALL mem%Add_var( 'tpe_overflow', mem%tpe_overflow,                      &
        & hgrid, surface,                                                      &
        & t_cf('tpe_overflow', 'm', 'Terra planet soil water overflow pool'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soilhyd_res', mem%wtr_soilhyd_res,                        &
        & hgrid, surface,                                                              &
        & t_cf('wtr_soilhyd_res', 'm3/(time step)', 'Vertical transport residual'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,213, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'wtr_rootzone', mem%wtr_rootzone,                      &
        & hgrid, surface,                                                      &
        & t_cf('wtr_rootzone', 'm', 'Liquid water content of the root zone'),  &
        & t_grib1(table, 255, grib_bits), grib2_desc,                          &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'ice_rootzone', mem%ice_rootzone,                                &
        & hgrid, surface,                                                                &
        & t_cf('ice_rootzone', 'm', 'Ice content (column) in the root zone'),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soil_sl', mem%wtr_soil_sl,                        &
        & hgrid, soil_w,                                                       &
        & t_cf('wtr_soil_sl', 'm', 'Water content of soil layers'),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'ice_soil_sl', mem%ice_soil_sl,                        &
        & hgrid, soil_w,                                                       &
        & t_cf('ice_soil_sl', 'm', 'Ice content of soil layers'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'vol_weq_soil_sl', mem%vol_weq_soil_sl,                       &
        & hgrid, soil_w,                                                              &
        & t_cf('vol_weq_soil_sl', 'm', 'Volumetric water equivalent of soil layers'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix,                                                             &
        & loutput=.TRUE., lrestart=.FALSE.,                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_freeze_sl', mem%wtr_freeze_sl,                              &
        & hgrid, soil_w,                                                                 &
        & t_cf('wtr_freeze_sl', 'kg m-2 s-1', 'Freezing water flux in soil layers'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'ice_melt_sl', mem%ice_melt_sl,                                  &
        & hgrid, soil_w,                                                                 &
        & t_cf('ice_melt_sl', 'kg m-2 s-1', 'Melting ice flux in soil layers'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'drainage_sl', mem%drainage_sl,                                  &
        & hgrid, soil_w,                                                                 &
        & t_cf('drainage_sl', 'kg m-2 s-1', 'Subsurface drainage on soil layers'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soil_sat_sl', mem%wtr_soil_sat_sl,                          &
        & hgrid, soil_w,                                                                 &
        & t_cf('wtr_soil_sat_sl', 'm', 'Water content in soil layers at saturation'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soil_fc_sl', mem%wtr_soil_fc_sl,                            &
        & hgrid, soil_w,                                                                 &
        & t_cf('wtr_soil_fc_sl', 'm', 'Water content in soil layers at field capacity'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soil_pwp_sl', mem%wtr_soil_pwp_sl,                          &
        & hgrid, soil_w,                                                                 &
        & t_cf('wtr_soil_pwp_sl', 'm',                                                   &
        &      'Water content in soil layers at permanent wilting point'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soil_res_sl', mem%wtr_soil_res_sl,                          &
        & hgrid, soil_w,                                                                 &
        & t_cf('wtr_soil_res_sl', 'm',                                                   &
        &      'Residual water content in soil layers'),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_soil_pot_scool_sl', mem%wtr_soil_pot_scool_sl,                   &
        & hgrid, soil_w,                                                                      &
        & t_cf('wtr_soil_pot_scool_sl', 'm', 'Potentially supercooled water on soil layers'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix,                                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_latflow_res_sl', mem%wtr_latflow_res_sl,                         &
        & hgrid, soil_w,                                                                      &
        & t_cf('wtr_latflow_res_sl', 'm', 'Water content of latflow res. on soil layers'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix,                                                                     &
        & lrestart=l_uniform, lrestart_cont=.TRUE.,                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_latflow_sl', mem%wtr_latflow_sl,                                 &
        & hgrid, soil_w,                                                                      &
        & t_cf('wtr_latflow_sl', 'kg m-2 s-1', 'Outflow from latflow res. on soil layers'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, lrestart=.FALSE.,                                                   &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_latflow_res_srf', mem%wtr_latflow_res_srf,             &
        & hgrid, surface,                                                           &
        & t_cf('wtr_latflow_res_srf', 'm', 'Water content of latflow res. srf.'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=l_uniform, lrestart_cont=.TRUE.,                                 &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_latflow_srf', mem%wtr_latflow_srf,                     &
        & hgrid, surface,                                                           &
        & t_cf('wtr_latflow_srf', 'kg m-2 s-1', 'Outflow from latflow res. srf'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix, lrestart=.FALSE.,                                         &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('wtr_rootzone_rel', mem%wtr_rootzone_rel,                         &
        & hgrid, surface,                                                                &
        & t_cf('wtr_rootzone_rel', '-', 'Relative root zone moisture'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'fract_pond', mem%fract_pond,                                  &
        & hgrid, surface,                                                              &
        & t_cf('fract_pond', '-', 'Inundated surface fraction'),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'weq_pond', mem%weq_pond,                                      &
        & hgrid, surface,                                                              &
        & t_cf('weq_pond', 'm (water equivalent)', 'Surface pond reservoir'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_pond', mem%wtr_pond,                                      &
        & hgrid, surface,                                                              &
        & t_cf('wtr_pond', 'm (water equivalent)', 'Liquid water in pond reservoir'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'ice_pond', mem%ice_pond,                                      &
        & hgrid, surface,                                                              &
        & t_cf('ice_pond', 'm (water equivalent)', 'Frozen water in pond reservoir'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'pond_melt', mem%pond_melt,                                    &
        & hgrid, surface,                                                              &
        & t_cf('pond_melt', 'kg m-2 s-1', 'Melting of pond ice'),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE.,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'pond_freeze', mem%pond_freeze,                                &
        & hgrid, surface,                                                              &
        & t_cf('pond_freeze', 'kg m-2 s-1', 'Freezing of pond water'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE.,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_pond_net_flx', mem%wtr_pond_net_flx,                      &
        & hgrid, surface,                                                              &
        & t_cf('wtr_pond_net_flx', 'kg m-2 s-1',                                       &
        &      'Net water flux into surface water ponds'),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE.,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_pond_max', mem%fract_pond_max,                               &
        & hgrid, surface,                                                                   &
        & t_cf('fract_pond_max', '-', 'Maximum surface fraction available for inundation'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & lrestart=.FALSE.,                                                                 &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE., isteptype=TSTEP_CONSTANT  )

      CALL mem%Add_var( 'weq_pond_max', mem%weq_pond_max,                                   &
        & hgrid, surface,                                                                   &
        & t_cf('weq_pond_max', 'm (water equivalent)', 'Maximum surface pond reservoir'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & lrestart=.FALSE.,                                                                 &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE., isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'evapo_pond', mem%evapo_pond,                          &
        & hgrid, surface,                                                      &
        & t_cf('evapo_pond', 'kg m-2 s-1',                                     &
        &      'Evaporation from ponds (w/o snow and skin res.'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    END IF

    ! Additional variables if the tile contains vegetation
    ! --------------------------
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'root_depth', mem%root_depth,                                    &
        & hgrid, surface,                                                                &
        & t_cf('root_depth', 'm', 'Rooting depth'),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'root_depth_sl', mem%root_depth_sl,                              &
        & hgrid, soil_w,                                                                 &
        & t_cf('root_depth_sl', 'm', 'Rooting depth within the soil layer'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'fract_snow_can', mem%fract_snow_can,                            &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_can', '-', 'Snow fraction on canopy'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=MEDIUM,                                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'weq_snow_can', mem%weq_snow_can,                                &
        & hgrid, surface,                                                                &
        & t_cf('weq_snow_can', 'm (water equivalent)', 'Snow amount on canopy'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=MEDIUM,                                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_rootzone_scool_pot', mem%wtr_rootzone_scool_pot,                        &
        & hgrid, surface,                                                                            &
        & t_cf('wtr_rootzone_scool_pot', 'm', 'Potentially supercooled water content in root zone'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
        & prefix, suffix,                                                                            &
        & lrestart=lrestart_local,                                                                   &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_rootzone_scool_act', mem%wtr_rootzone_scool_act,            &
        & hgrid, surface,                                                                &
        & t_cf('wtr_rootzone_scool_act', 'm', 'Supercooled water content in root zone'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_rootzone_avail', mem%wtr_rootzone_avail,                    &
        & hgrid, surface,                                                                &
        & t_cf('wtr_rootzone_avail', 'm', 'Plant available water content in root zone'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=FULL,                                                             &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_rootzone_avail_max', mem%wtr_rootzone_avail_max,                        &
        & hgrid, surface,                                                                            &
        & t_cf('wtr_rootzone_avail_max', 'm', 'Maximum plant available water content in root zone'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
        & prefix, suffix,                                                                            &
        & output_level=FULL,                                                                         &
        & lrestart=.FALSE.,                                                                          &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'water_stress', mem%water_stress,                                &
        & hgrid, surface,                                                                &
        & t_cf('water_stress', '-', 'Water stress factor'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & loutput=.TRUE., lrestart=.FALSE., initval_r=0.0_wp )

#ifndef __QUINCY_STANDALONE__
      IF (.NOT. model%Is_process_enabled(ASSIMI_)) THEN
#endif
        CALL mem%Add_var( 'canopy_cond_unlimited', mem%canopy_cond_unlimited,            &
          & hgrid, surface,                                                              &
          & t_cf('canopy_cond_unlimited', 'm/s',                                         &
          &      'Canopy conductance ignoring water limitation'),                        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & loutput=.TRUE., initval_r=0.0_wp )

        CALL mem%Add_var( 'canopy_cond_limited', mem%canopy_cond_limited,                &
          & hgrid, surface,                                                              &
          & t_cf('canopy_cond_limited', 'm/s',                                           &
          &      'Canopy conductance accounting for water limitation'),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & loutput=.TRUE., initval_r=0.0_wp )
#ifndef __QUINCY_STANDALONE__
      END IF
#endif

      CALL mem%Add_var( 'transpiration', mem%transpiration,                              &
        & hgrid, surface,                                                                &
        & t_cf('transpiration', 'kg m-2 s-1', 'Transpiration'),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & loutput=.TRUE., output_level=BASIC,                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('wtr_plant_avail_rel', mem%wtr_plant_avail_rel, &
        & hgrid, surface, &
        & t_cf('wtr_plant_avail_rel', '', 'relative plant available soil moisture in root zone'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = NONE, &
        & loutput = .FALSE., &
        & lrestart = lrestart_local, &
        & initval_r = 0.0_wp)

#ifndef __NO_QUINCY__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        CALL mem%Add_var('w_soil_root_fc', mem%w_soil_root_fc, &
          & hgrid, surface, &
          & t_cf('w_soil_root_fc', 'm', 'Water content at field capacity in root zone of the soil'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

        CALL mem%Add_var('w_soil_root_pwp', mem%w_soil_root_pwp, &
          & hgrid, surface, &
          & t_cf('w_soil_root_pwp', 'm', 'Water content at permanent wilting point in root zone of the soil'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

        CALL mem%Add_var('frac_w_lat_loss_sl', mem%frac_w_lat_loss_sl, &
          & hgrid, soil_w, &
          & t_cf('frac_w_lat_loss_sl', '', 'fraction of lateral (horizontal) water loss of w_soil_sl_old'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & loutput = .FALSE., &
          & lrestart = .FALSE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('frac_wtr_transp_down_sl', mem%frac_wtr_transp_down_sl, &
          & hgrid, soil_w, &
          & t_cf('frac_wtr_transp_down_sl', 'fraction s-1', 'fraction of mass transferred to below layer'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .TRUE., &
          & lrestart = .FALSE., &
          & initval_r = 0.0_wp)
      END SELECT
#endif

    END IF ! VEG_TYPE/LAND_TYPE

    ! Additional variables for tiles containing GLACIER lct
    ! --------------------------
    IF ( (     One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      ! TODO: Remove fract_snow_glac - never used
      CALL mem%Add_var( 'fract_snow_glac', mem%fract_snow_glac,                          &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_glac', '-', 'Snow fraction on glacier'),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'weq_glac', mem%weq_glac,                                        &
        & hgrid, surface,                                                                &
        & t_cf('weq_glac', 'm (water equivalent)', 'Glacier depth (including snow)'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & initval_r=20.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'runoff_glac', mem%runoff_glac,                                  &
        & hgrid, surface,                                                                &
        & t_cf('runoff_glac', 'kg m-2 s-1', 'Runoff from glacier (rain + melt)'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & initval_r=20.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for tiles containing LAKE lct
    ! --------------------------
    IF (      One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 'evapo_wtr', mem%evapo_wtr,                                &
        & hgrid, surface,                                                          &
        & t_cf('evapo_wtr', 'kg m-2 s-1', 'Evaporation from lake water'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & lrestart=.FALSE.,                                                        &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        CALL mem%Add_var( 'evapo_ice', mem%evapo_ice,                              &
          & hgrid, surface,                                                        &
          & t_cf('evapo_ice', 'kg m-2 s-1', 'Evaporation from lake ice'),          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & lrestart=.FALSE.,                                                      &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'fract_snow_lice', mem%fract_snow_lice,                  &
          & hgrid, surface,                                                        &
          & t_cf('fract_snow_lice', '-',                                           &
          &      'Snow fraction on lake ice (rel. to ice fraction)'),              &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & output_level=MEDIUM,                                                   &
          & initval_r=0.0_wp )

        CALL mem%Add_var( 'weq_snow_lice', mem%weq_snow_lice,                        &
          & hgrid, surface,                                                          &
          & t_cf('weq_snow_lice', 'm water equivalent', 'Snow amount on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
          & prefix, suffix,                                                          &
          & output_level=MEDIUM,                                                     &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

    END IF

    ! Diagnostic 0 dim. global land variables for experiment monitoring
    ! ------------------------------------------------------------------
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var('trans_gmean', mem%trans_gmean,                                   &
        & hgrid, surface,                                                                &
        & t_cf('trans_gmean', 'kg m-2 s-1', 'Global mean transpiration'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('evapotrans_gmean', mem%evapotrans_gmean,                         &
        & hgrid, surface,                                                                &
        & t_cf('evapotrans_gmean', 'kg m-2 s-1', 'Global land evapotranspiration'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('weq_land_gsum', mem%weq_land_gsum,                               &
        & hgrid, surface,                                                                &
        & t_cf('weq_land_gsum', 'km3', 'Global amount of land water and ice'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('discharge_ocean_gsum', mem%discharge_ocean_gsum,                 &
        & hgrid, surface,                                                                &
        & t_cf('discharge_ocean_gsum', 'Sv', 'Global water discharge to the oceans'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('wtr_rootzone_rel_gmean', mem%wtr_rootzone_rel_gmean,             &
        & hgrid, surface,                                                                &
        & t_cf('wtr_rootzone_rel_gmean', '-', 'Global mean relative rootzone moisture'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('fract_snow_gsum', mem%fract_snow_gsum,                           &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_gsum', 'Mio km2', 'Global snow area (including glaciers)'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('weq_snow_gsum', mem%weq_snow_gsum,                               &
        & hgrid, surface,                                                                &
        & t_cf('weq_snow_gsum', 'Gt', 'Global snow amount on non-glacier land'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('weq_balance_err_gsum', mem%weq_balance_err_gsum,                 &
        & hgrid, surface,                                                                &
        & t_cf('weq_balance_err_gsum', 'm3/(time step)', 'Global land water balance error'),&
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    END IF

    END SUBROUTINE Init_hydro_memory
#endif
END MODULE mo_hydro_memory_class
