!> QUINCY soil-physics process memory
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
!>#### definition and init of (memory) variables for the soil-physics-quincy process
!>
MODULE mo_spq_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE, LAKE_TYPE, BARE_TYPE, GLACIER_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_spq_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars      = 240


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for spq memory
  !!
  !! @par includes: \n
  !!    spq variables spq\%VAR
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_spq_memory


    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth      , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    temp_srf_eff_4  , & !< effective surface temperature ** 4.0 [K**4.0] (in jsbach: t_eff4)
                                    zril_old        , & !< previous timestep's Reynold's number
                                    elevation       , & !< currently 0.0 [m]
                                    s_star              !< surface dry static energy               [m2 s-2 ?]

    TYPE(t_jsb_var_real2d)          :: &
      & evapotranspiration                !< evapotranspiration [kg m-2 s-1]

    ! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d)          :: &
                                    interception,     & !< surface interception evaporation                                         [kg m-2 s-1]
                                    evapopot,         & !< potential evaporation                                                    [kg m-2 s-1]
                                    evaporation,      & !< surface evaporation                                                      [kg m-2 s-1]
                                    evaporation_snow, & !< snow surface evaporation                                                 [kg m-2 s-1] unit correct ?
                                    srf_runoff,       & !< surface runoff                                                           [kg m-2 s-1]
                                    drainage,         & !< drainage                                                                 [kg m-2 s-1]
                                    gw_runoff,        & !< lateral runoff                                                           [kg m-2 s-1]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage                [1 s-1]

    ! Additional variables for spq
    TYPE(t_jsb_var_real3d)          :: &
                                    gw_runoff_sl,     & !< lateral (horizontal) soilwater to groundwater runoff per layer           [kg m-2 yr-1]
                                    drainage_sl,      & !< vertical water flow (drainage) per soil layer across soil layers         [kg m-2 s-1]
                                    saxtonA,          & !< coefficient in moisture-tension relationship [unitless]
                                    saxtonB,          & !< coefficient in moisture-tension relationship [unitless]
                                    saxtonC,          & !< exponent of moisture-conductivity relationship [unitless]
                                    kdiff_sat_sl        !< saturated hydraulic conductivity                                         [m s-1]

    ! heat fluxes
    TYPE(t_jsb_var_real2d)          :: &
                                    sensible_heat_flx     , &  !< sensible heat flux [W m-2]
                                    latent_heat_flx       , &  !< latent heat flux [W m-2]
                                    ground_heat_flx       , &  !< ground heat flux [W m-2]
                                    ground_heat_flx_old        !< ground heat flux of previous timestep [W m-2]

    TYPE(t_jsb_var_real3d)          :: &
      & w_soil_freeze_flux, &                            !< soil water flux from liquid water to ice       [m]
      & w_soil_melt_flux                                 !< soil water flux from ice to liquid water       [m]

    ! soil physical properties (not sure where they should belong to, but they are needed in the long-term)
    TYPE(t_jsb_var_real3d)          :: &
                                    heat_capa_sl,         &  !< soil heat capacity [J kg-1]
                                    therm_cond_sl            !< soil thermal conductivity [W m-2 K-1]

    ! variables for snow
    TYPE(t_jsb_var_real2d)          :: snow_height, &           !< height (thickness) of all snow layers                         [m]
                                       snow_soil_heat_flux, &   !< heat flux between snow and soil                          [W m-2]
                                       snow_srf_heat_flux,  &   !< heat flux between atmoshere and snow                     [W m-2]
                                       snow_melt_to_soil        !< Snow melt water to soil flux                              [kg m-2 s-1]

    TYPE(t_jsb_var_real3d)          :: snow_present_snl       , &  !< variable to check if there is any snow for each layer     [unitless]
                                       t_snow_snl             , &  !< temperature of snow layer                                 [K]
                                       w_snow_snl             , &  !< water content of snow layer                               [m]
                                       snow_lay_thickness_snl , &  !< snow layer thickness                                      [m]
                                       w_snow_max_snl         , &  !< maximum water content of snow layer                       [m]
                                       w_liquid_snl                !< Liquid water in snow layers                               [m]

  CONTAINS
    PROCEDURE :: Init => Init_spq_memory
  END TYPE t_spq_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SPQ_ process
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_spq_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,         ONLY: BASIC , MEDIUM, FULL
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_CONSTANT
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_quincy_output_class, ONLY: unitless
    USE mo_jsb_physical_constants, ONLY: &
      & Tzero,                           &
      & dens_snow                          ! Density of snow (l_dynsnow=.FALSE.)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_spq_memory), INTENT(inout), TARGET  :: mem             !< spq memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !< process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !< tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: lib_id          !< id of primary lct in lctlib
    INTEGER,              INTENT(in)            :: model_id        !< model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_w                 ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_snow_spq               ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_spq_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model          => Get_model(model_id)
    table          =  tables(1)
    hgrid          => Get_grid(model%grid_id)
    surface        => Get_vgrid('surface')
    vgrid_soil_w   => Get_vgrid('soil_depth_water')
    vgrid_snow_spq => Get_vgrid('snow_layer_spq')

    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('soil_depth', mem%soil_depth, &
        & hgrid, surface, &
        & t_cf('soil_depth', 'm', 'Soil depth derived from textures (bedrock), for hydrology'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('temp_srf_eff_4', mem%temp_srf_eff_4, &
        & hgrid, surface, &
        & t_cf('temp_srf_eff_4', 'K**4.0', 'effective surface temperature ** 4.0'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 273.15_wp ** 4.0_wp)

      CALL mem%Add_var('zril_old', mem%zril_old, &
        & hgrid, surface, &
        & t_cf('zril_old', '', 'previous timesteps Reynolds number'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('elevation', mem%elevation, &
        & hgrid, surface, &
        & t_cf('elevation', 'm', 'Mean orography'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart=.FALSE., &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT)

      CALL mem%Add_var('s_star', mem%s_star, &
        & hgrid, surface, &
        & t_cf('s_star', 'm2 s-2 ?', 'surface dry static energy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=2.9E5_wp)

      CALL mem%Add_var('evapotranspiration', mem%evapotranspiration, &
        & hgrid, surface, &
        & t_cf('evapotranspiration', 'kg m-2 s-1 ?', 'evapotranspiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('interception', mem%interception, &
        & hgrid, surface, &
        & t_cf('interception', 'kg m-2 s-1', 'surface interception evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('evapopot', mem%evapopot, &
        & hgrid, surface, &
        & t_cf('evapopot', 'kg m-2 s-1', 'potential evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('evaporation', mem%evaporation, &
        & hgrid, surface, &
        & t_cf('evaporation', 'kg m-2 s-1', 'surface evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('evaporation_snow', mem%evaporation_snow, &
        & hgrid, surface, &
        & t_cf('evaporation_snow', 'kg m-2 s-1', 'snow surface evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('srf_runoff', mem%srf_runoff, &
        & hgrid, surface, &
        & t_cf('srf_runoff', 'kg m-2 s-1', 'surface runoff'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('drainage', mem%drainage, &
        & hgrid, surface, &
        & t_cf('drainage', 'kg m-2 s-1', 'drainage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('gw_runoff', mem%gw_runoff, &
        & hgrid, surface, &
        & t_cf('gw_runoff', 'kg m-2 s-1', 'sum of lateral runoff across all soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('drainage_fraction', mem%drainage_fraction, &
        & hgrid, surface, &
        & t_cf('drainage_fraction', '1 s-1', 'water mass fraction lost from the soil column by drainage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gw_runoff_sl', mem%gw_runoff_sl, &
        & hgrid, vgrid_soil_w, &
        & t_cf('gw_runoff_sl', 'kg m-2 s-1', 'lateral groundwater runoff per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('drainage_sl', mem%drainage_sl, &
        & hgrid, vgrid_soil_w, &
        & t_cf('drainage_sl', 'kg m-2 s-1', 'drainage of soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('saxtonA', mem%saxtonA, &
        & hgrid, vgrid_soil_w, &
        & t_cf('saxtonA', unitless, 'coefficient in moisture-tension relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonB', mem%saxtonB, &
        & hgrid, vgrid_soil_w, &
        & t_cf('saxtonB', unitless, 'coefficient in moisture-tension relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonC', mem%saxtonC, &
        & hgrid, vgrid_soil_w, &
        & t_cf('saxtonC', unitless, 'exponent of moisture-conductivity relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('kdiff_sat_sl', mem%kdiff_sat_sl, &
        & hgrid, vgrid_soil_w, &
        & t_cf('kdiff_sat_sl', 'm s-1', 'saturated hydraulic conductivity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all = .TRUE.)

      CALL mem%Add_var('sensible_heat_flx', mem%sensible_heat_flx, &
        & hgrid, surface, &
        & t_cf('sensible_heat_flx', 'W m-2', 'sensible heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('latent_heat_flx', mem%latent_heat_flx, &
        & hgrid, surface, &
        & t_cf('latent_heat_flx', 'W m-2', 'latent heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('ground_heat_flx', mem%ground_heat_flx, &
        & hgrid, surface, &
        & t_cf('ground_heat_flx', 'W m-2', 'ground heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('ground_heat_flx_old', mem%ground_heat_flx_old, &
        & hgrid, surface, &
        & t_cf('ground_heat_flx_old', 'W m-2', 'ground heat flux of previous timestep'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('w_soil_freeze_flux', mem%w_soil_freeze_flux, &
        & hgrid, vgrid_soil_w, &
        & t_cf('w_soil_freeze_flux', 'm', 'Water flux from liquid to ice from freezing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('w_soil_melt_flux', mem%w_soil_melt_flux, &
        & hgrid, vgrid_soil_w, &
        & t_cf('w_soil_melt_flux', 'm', 'Water transfer from ice to liquid from melting'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp, &
        & l_aggregate_all=.TRUE.)

      CALL mem%Add_var('heat_capa_sl', mem%heat_capa_sl, &
        & hgrid, vgrid_soil_w, &
        & t_cf('heat_capa_sl', 'J kg-1', 'soil heat capacity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('therm_cond_sl', mem%therm_cond_sl, &
        & hgrid, vgrid_soil_w, &
        & t_cf('therm_cond_sl', 'W m-2 K-1', 'soil thermal conductivity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('snow_height', mem%snow_height, &
          &  hgrid, surface, &
          &  t_cf('snow_height','m', 'snow depth'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .TRUE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &  ! initvalue
          & l_aggregate_all=.TRUE.)

      CALL mem%Add_var( 'snow_soil_heat_flux', mem%snow_soil_heat_flux, &
          &  hgrid, surface, &
          &  t_cf('snow_soil_heat_flux','W m-2', 'snow soil heat flux'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, & ! initvalue
          &  l_aggregate_all=.TRUE.)

      CALL mem%Add_var('snow_srf_heat_flux', mem%snow_srf_heat_flux, &
          &  hgrid, surface, &
          &  t_cf('snow_srf_heat_flux','W m-2', 'snow srf heat flux'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &  ! initvalue
          &  l_aggregate_all=.TRUE.)

      CALL mem%Add_var('snow_melt_to_soil', mem%snow_melt_to_soil, &
          &  hgrid, surface, &
          &  t_cf('snow_melt_to_soil','kg m-2 s-1', 'snow melt water to soil'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &  ! initvalue
          &  l_aggregate_all=.TRUE.)

      CALL mem%Add_var('snow_present_snl', mem%snow_present_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('snow_present_snl',unitless, 'check if there is snow'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('t_snow_snl', mem%t_snow_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('t_snow_snl','m', 'temperature snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

      CALL mem%Add_var('w_snow_snl', mem%w_snow_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_snow_snl','m', 'water content snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

      CALL mem%Add_var('snow_lay_thickness_snl', mem%snow_lay_thickness_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('snow_lay_thickness_snl','m', 'snow layer thickness'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

      CALL mem%Add_var('w_snow_max_snl', mem%w_snow_max_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_snow_max_snl','m', 'maximum water content snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

      CALL mem%Add_var('w_liquid_snl', mem%w_liquid_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_liquid_snl','m', 'w_liquid_snl'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp, &
          & l_aggregate_all = .TRUE.)

    END IF ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)

  END SUBROUTINE Init_spq_memory

#endif
END MODULE mo_spq_memory_class
