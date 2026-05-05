!> Contains structures and methods for soil and snow energy config
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
MODULE mo_sse_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string, int2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_config

  TYPE, EXTENDS(t_jsb_config) :: t_sse_config
    INTEGER  :: nsnow           !< Maximum number of snow layers
    LOGICAL  :: l_snow          !< T: Use multi-layer snow scheme;
                                !  Note: This affects thermal snow and soil properties:
                                !  l_dynsnow, l_heat_cap_dyn and l_heat_cond_dyn are only effective if l_snow=.true.
    LOGICAL  :: l_dynsnow       !< T: Calculate snow heat cond. and heat cap. dynamically depending on snow density
    LOGICAL  :: l_heat_cap_dyn  !< T: Dynamic calculation of soil heat capacity;
                                !  F: Static map or FAO data, depending on l_heat_cap_map
    LOGICAL  :: l_heat_cond_dyn !< T: Dynamic calculation of soil heat conductivity;
                                !  F: Static map or FAO data, depending on l_heat_cond_map
    !
    LOGICAL  :: l_heat_cap_map  !< T: Use soil heat capacity from input map, F: derive from FAO
    LOGICAL  :: l_heat_cond_map !< T: Use soil heat conductivity from input map, F: derive from heat cap. and FAO thermal diff.
    LOGICAL  :: l_heat_cap_pond !< T: Account for the heat capacity of ponds in the top soil layer heat capacity.
    LOGICAL  :: l_soil_texture  !< T: Deduce mineral soil thermal parameters from soil texture
    LOGICAL  :: l_freeze        !< T: Consider freezing and thawing of soil water
    LOGICAL  :: l_supercool     !< T: Allow for supercooled soil water
    LOGICAL  :: l_uniform_snow  !< T: Temporary bugfix to use zero or one snow fraction for soil temperature computation
    REAL(wp) :: w_soil_critical !< Critical water/ice content in upper soil layer for correction of new surface
                                !! temperature for freezing/melting [m water equivalent]
    ! quincy
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename  !< IQ (QUINCY) soil input data
    INTEGER                     :: quincy_nsoil_e           !< number of soil layers for soil energy calculations (equivalent 'quincy_nsoil_w' in HYDRO_)
#ifdef __QUINCY_STANDALONE__
    ! quincy standalone
    REAL(wp) :: qs_soil_sand         !< soil sand proportion - site specific from forcing info file
    REAL(wp) :: qs_soil_silt         !< soil silt proportion - site specific from forcing info file
    REAL(wp) :: qs_soil_clay         !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
    REAL(wp) :: qs_bulk_density      !< soil bulk density
#endif

  CONTAINS
    PROCEDURE :: Init => Init_sse_config
  END TYPE t_sse_config

  INTEGER, PARAMETER :: max_snow_layers = 20
  INTEGER, PARAMETER :: max_soil_layers = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_config_class'

CONTAINS

  SUBROUTINE Init_sse_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_sse_constants,      ONLY: snow_depth_min
    USE mo_jsb_math_constants, ONLY: eps8

    CLASS(t_sse_config), INTENT(inout) :: config

    LOGICAL  :: active, lrestart_cont
    LOGICAL  :: l_snow, l_dynsnow, l_heat_cap_dyn, l_heat_cond_dyn, l_heat_cap_map, l_heat_cond_map, l_heat_cap_pond
    LOGICAL  :: l_freeze, l_supercool, l_soil_texture, l_uniform_snow
    REAL(wp) :: w_soil_critical
    INTEGER  :: nsnow
    REAL(wp) :: dz_snow(max_snow_layers)
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    ! quincy
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename
    REAL(wp) :: soil_layer_profile_ubound_estimate(max_soil_layers) !< intial estimate of soil profile using the upper bound of each soil layer [m]
    REAL(wp) :: k_soil_profile  !< for soil layer calculations
    REAL(wp) :: min_layer_depth !< for soil layer calculations
    INTEGER  :: isoil           !< loop over soil layers
    INTEGER  :: quincy_nsoil_e  !< quincy specific number of soil layers
#ifdef __QUINCY_STANDALONE__
    ! quincy standalone
    REAL(wp) :: qs_soil_sand
    REAL(wp) :: qs_soil_silt
    REAL(wp) :: qs_soil_clay
    REAL(wp) :: qs_bulk_density
#endif

    NAMELIST /jsb_sse_nml/                      &
      & active,                                 &
      & lrestart_cont,                          &
      & nsnow, dz_snow, l_snow, l_dynsnow,      &
      & l_heat_cap_dyn, l_heat_cond_dyn,        &
      & l_heat_cap_map, l_heat_cond_map,        &
      & l_heat_cap_pond,                        &
      & l_soil_texture, l_uniform_snow,         &
      & l_freeze, l_supercool, w_soil_critical, &
      & ic_filename,                            &
      & bc_filename,                            &
      & bc_quincy_soil_filename,                &
#ifdef __QUINCY_STANDALONE__
      & qs_soil_sand, &
      & qs_soil_silt, &
      & qs_soil_clay, &
      & qs_bulk_density, &
#endif
      & quincy_nsoil_e

    INTEGER :: nml_handler, nml_unit, istat, i
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    REAL(wp), ALLOCATABLE :: dz_energy(:)               !< thickness of soil layers
    REAL(wp), ALLOCATABLE :: ubounds_soil_lay_energy(:) !< upper bound of soil layers (larger value compared to lower bound)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_e, vgrid_snow_e

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting soil and snow energy configuration')

    ! Set defaults
    active           = .TRUE.
    lrestart_cont    = .FALSE.  ! TRUE: Continue although SSE variables are missing in restart file
    nsnow            = 5
    dz_snow(:)       = 0._wp
    dz_snow(1:nsnow) = 0.05_wp
    l_snow           = .TRUE.
    l_dynsnow        = .TRUE.
    l_heat_cap_dyn   = .TRUE.
    l_heat_cond_dyn  = .TRUE.
    l_heat_cap_map   = .FALSE.
    l_heat_cond_map  = .FALSE.
    l_heat_cap_pond  = .FALSE.
    l_soil_texture   = .FALSE.
    l_freeze         = .TRUE.
    l_supercool      = .TRUE.
    l_uniform_snow   = .TRUE.
    w_soil_critical  = 5.85036E-3_wp
    ic_filename      = 'ic_land_soil.nc'
    bc_filename      = 'bc_land_soil.nc'
    bc_quincy_soil_filename = 'bc_quincy_soil.nc'
    quincy_nsoil_e          = 0           ! by default the values is not set via namelist, but used from 'ic_filename'
#ifdef __QUINCY_STANDALONE__
    qs_soil_sand            = 0.3_wp
    qs_soil_silt            = 0.4_wp
    qs_soil_clay            = 0.3_wp
    qs_bulk_density         = 1500.0_wp
#endif

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_sse_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_sse_nml)

    CALL close_nml(nml_handler)

    config%active        = active
    config%lrestart_cont = lrestart_cont
    config%ic_filename             = ic_filename
    config%bc_filename             = bc_filename
    config%bc_quincy_soil_filename = bc_quincy_soil_filename

    config%nsnow         = nsnow
    IF (l_snow .AND. nsnow > max_snow_layers) THEN
      CALL finish(TRIM(routine), 'Too many snow layers, maximum is '//TRIM(int2string(max_snow_layers)))
    END IF
    IF (l_snow .AND. nsnow < 3) THEN
      CALL finish(TRIM(routine), 'Number of snow layers must be larger than 2.')
    END IF
    config%l_snow        = l_snow
    IF (l_snow) THEN
      CALL message(TRIM(routine), 'Using '//TRIM(int2string(nsnow))//' snow layers.')
    ELSE
      CALL message(TRIM(routine), 'Not using multi-layer snow model')
    END IF

    IF (.NOT. l_snow .AND. l_dynsnow) THEN
      config%l_dynsnow = .FALSE.
      CALL message(TRIM(routine), 'l_dynsnow not in effect for l_snow=.false.')
    ELSE
      config%l_dynsnow       = l_dynsnow
    END IF
    IF (.NOT. l_snow .AND. l_heat_cap_dyn) THEN
      config%l_heat_cap_dyn = .FALSE.
      CALL message(TRIM(routine), 'l_heat_cap_dyn not in effect for l_snow=.false.')
    ELSE
      config%l_heat_cap_dyn  = l_heat_cap_dyn
    END IF
    IF (.NOT. l_snow .AND. l_heat_cond_dyn) THEN
      config%l_heat_cond_dyn = .FALSE.
      CALL message(TRIM(routine), 'l_heat_cond_dyn not in effect for l_snow=.false.')
    ELSE
      config%l_heat_cond_dyn = l_heat_cond_dyn
    END IF
    config%l_heat_cap_map  = l_heat_cap_map
    config%l_heat_cond_map = l_heat_cond_map
    config%l_soil_texture  = l_soil_texture
    IF (l_soil_texture) THEN
      CALL message(TRIM(routine), 'Using soil texture to derive heat capacity and heat conductivity')
    END IF
    config%l_heat_cap_pond = l_heat_cap_pond
    IF (l_heat_cap_pond) THEN
      CALL message(TRIM(routine), 'The heat capacity of ponds is accounted for in top soil layer heat capacity.')
    END IF
    config%l_freeze        = l_freeze
    config%l_supercool     = l_supercool
    config%l_uniform_snow  = l_uniform_snow
    IF (l_uniform_snow) THEN
      CALL message(TRIM(routine), 'Assuming uniform snow cover for soil temperature computation')
    END IF

    config%w_soil_critical     = w_soil_critical
    CALL message(TRIM(routine), 'Critical water/ice content in upper soil layer for correction of '// &
      &                         'surface temperature for freezing/melting: '//TRIM(real2string(w_soil_critical)))

    ! quincy
    config%quincy_nsoil_e  = quincy_nsoil_e
#ifdef __QUINCY_STANDALONE__
    ! quincy standalone
    config%qs_soil_sand    = qs_soil_sand
    config%qs_soil_silt    = qs_soil_silt
    config%qs_soil_clay    = 1.0_wp - qs_soil_sand - qs_soil_silt
    config%qs_bulk_density = qs_bulk_density

    ! check for an error in clay content of soil
    IF(config%qs_soil_clay < eps8) THEN
      WRITE(message_text,'(a)') 'Invalid proportion of qs_soil_clay (value < eps8)'
      CALL finish(routine, message_text)
    END IF
#endif

    IF (.NOT. active) RETURN

    !
    !> create vertical grid for soil energy calculations
    !
    ! init vgrid with values from 'ic_filename' (default MODEL_JSBACH)
    IF (config%quincy_nsoil_e == 0) THEN
#ifndef __QUINCY_STANDALONE__
      input_file = jsb_netcdf_open_input(ic_filename)

      ! @todo: At the moment, the soil layers for the energy calculations are the same as for the hydrology
      ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
      nsoil = SIZE(ptr_1D)
      ALLOCATE(depths(nsoil+1))
      ALLOCATE(mids(nsoil))
      depths(1) = 0._wp
      depths(2:nsoil+1) = ptr_1D(1:nsoil)
      mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
      DEALLOCATE(ptr_1D)

      CALL input_file%Close()

      vgrid_soil_e  => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
        & levels    = mids                 (1:nsoil  ),                              &
        & lbounds   = depths               (1:nsoil  ),                              &
        & ubounds   = depths               (2:nsoil+1),                              &
        & units='m')
      CALL register_vgrid(vgrid_soil_e)

      DEALLOCATE(depths, mids)

      CALL message(TRIM(routine), 'Init vgrid soil_depth_energy with values from input file: '//TRIM(ic_filename))
#endif
    ! init vgrid with both from namelist and predefined (default MODEL_QUINCY)
    ELSE
      ! finish if number of soil layers has an unexpected value
      IF (config%quincy_nsoil_e < 5 .OR. config%quincy_nsoil_e > 20) THEN
        WRITE(message_text,'(a)') 'Invalid number of soil layers quincy_nsoil_e defined in namelist (outside [5,20])'
        CALL finish(routine, message_text)
      END IF

      ! parameter for soil layer calculations
      k_soil_profile  = 0.25_wp
      min_layer_depth = 0.065_wp

      nsoil = config%quincy_nsoil_e
      ALLOCATE(ubounds_soil_lay_energy(nsoil))
      ALLOCATE(dz_energy(nsoil))

      IF (nsoil > 5) THEN
        ! initial calculation of soil profile (for soil energy calculations) calculating the upper bound of each layer
        DO isoil = 1,nsoil
          soil_layer_profile_ubound_estimate(isoil) = &
            &   min_layer_depth &
            &   * EXP(k_soil_profile * isoil * REAL(max_soil_layers, wp) / REAL(nsoil, wp)) &
            &   - min_layer_depth
        END DO
        ! final calculation of soil profile, calculating the layer thickness (dz(:)) and correcting for minimum soil-layer thickness
        dz_energy(1) = min_layer_depth
        DO isoil = 2,nsoil
          dz_energy(isoil) = MAX(min_layer_depth, &
            &   soil_layer_profile_ubound_estimate(isoil) - soil_layer_profile_ubound_estimate(isoil-1))
        END DO
      ELSE
        dz_energy(1:nsoil) = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
      END IF

      ! calc upper bound of soil layers (water and energy) from layer thickness
      ! note: upper bound is the larger value compared to lower bound of the same layer!
      ubounds_soil_lay_energy(1)        = dz_energy(1)
      DO isoil = 2,nsoil
        ubounds_soil_lay_energy(isoil)  = ubounds_soil_lay_energy(isoil-1) + dz_energy(isoil)
      END DO

      !< Create vertical soil-layer axis (vgrid)
      !!
      !! these are "infrastructure" values, ubounds and dz must be corrected/limited by site-specific bedrock depth "soil_depth" !
      !!
      !! site-specific soil-layer thickness values are stored in soil_depth_sl (and calc in spq_init) \n
      !! the function new_vgrid() does: (a) check for 'dz(:) <= zero', and (b) calc levels & lbounds from ubounds and dz \n
      !! levels is defined as the depth at the center of the layer: levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
      vgrid_soil_e => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
        &   longname='Soil (energy) layers from SSE_ config', &
        &   units='m', &
            ! levels=  , &      ! is calculated from ubounds and dz
            ! lbounds=  , &     ! is calculated from ubounds and dz
        &   ubounds=ubounds_soil_lay_energy(:), &
        &   dz=dz_energy(:))
      CALL register_vgrid(vgrid_soil_e)

      DEALLOCATE(ubounds_soil_lay_energy)
      DEALLOCATE(dz_energy)
    END IF

    ! Info saved to logfile
    WRITE(message_text, *) 'Soil levels in soil energy (upper) [m]: ', vgrid_soil_e%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (mid)   [m]: ', vgrid_soil_e%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (lower) [m]: ', vgrid_soil_e%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil level depths in soil energy   [m]: ', vgrid_soil_e%dz
    CALL message(TRIM(routine), message_text)

    !
    !> create vertical grid for snow energy calculations
    !
    IF (l_snow) THEN
      IF (dz_snow(1) < snow_depth_min) THEN
        CALL finish(TRIM(routine), 'Depth of first snow layer should be larger than snow_depth_min=' &
          &                         //real2string(snow_depth_min))
      END IF
      ALLOCATE(depths(nsnow+1))
      ALLOCATE(mids(nsnow))
      depths(1) = 0._wp
      DO i=1,nsnow
        depths(i+1) = depths(i) + dz_snow(i)
      END DO
      mids(1:nsnow) = (depths(1:nsnow) + depths(2:nsnow+1)) / 2._wp

      vgrid_snow_e  => new_vgrid('snow_depth_energy', ZAXIS_GENERIC, nsnow,     &
        & levels    = mids                 (1:nsnow  ),                         &
        & lbounds   = depths               (1:nsnow  ),                         &
        & ubounds   = depths               (2:nsnow+1),                         &
        & units='m')
      CALL register_vgrid(vgrid_snow_e)

      ! Info saved to logfile
      WRITE(message_text, *) 'Snow levels in soil energy (upper) [m]: ', vgrid_snow_e%lbounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow levels in soil energy (mid)   [m]: ', vgrid_snow_e%levels
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow levels in soil energy (lower) [m]: ', vgrid_snow_e%ubounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow level depths in soil energy   [m]: ', vgrid_snow_e%dz
      CALL message(TRIM(routine), message_text)
      DEALLOCATE(depths, mids)
    END IF

  END SUBROUTINE Init_sse_config

#endif
END MODULE mo_sse_config_class
