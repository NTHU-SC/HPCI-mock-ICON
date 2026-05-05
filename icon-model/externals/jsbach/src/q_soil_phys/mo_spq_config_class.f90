!> QUINCY soil-physics process config
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
!>#### define soil-physics-quincy config structure, read soil-physics-quincy namelist and init configuration parameters
!>
MODULE mo_spq_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_spq_config, max_soil_layers

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the spq process, derived from t_jsb_config
  !!
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_spq_config
    LOGICAL                     :: spq_deactivate_spq   !< only used with quincy-standalone: deactivate SPQ_ consistently with use_soil_phys_jsbach
    LOGICAL                     :: flag_snow            !< on/off snow accumulation
    REAL(wp)                    :: soil_depth           !< actual soil and rooting depth
    REAL(wp)                    :: wsr_capacity         !< Water holding capacity of ground skin reservoir [m water equivalent]
    REAL(wp)                    :: wsn_capacity         !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
    REAL(wp)                    :: soil_awc_prescribe, &
                                   soil_theta_prescribe
    REAL(wp)                    :: elevation
    REAL(wp)                    :: soil_sand            !< soil sand proportion - site specific from forcing info file
    REAL(wp)                    :: soil_silt            !< soil silt proportion - site specific from forcing info file
    REAL(wp)                    :: soil_clay            !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
    REAL(wp)                    :: bulk_density
    CHARACTER(len=filename_max) :: bc_sso_filename          !< elevation and oro_stddev
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename  !< IQ (QUINCY) soil input data
  CONTAINS
    PROCEDURE :: Init => Init_spq_config
  END TYPE t_spq_config

  INTEGER, PARAMETER :: max_snow_layers = 20  ! consistent with JSBACH4
  INTEGER, PARAMETER :: max_soil_layers = 20  ! consistent with JSBACH4

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_spq_config
  !!
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_spq_config(config)
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
    USE mo_spq_constants,      ONLY: snow_height_min
    USE mo_jsb_math_constants, ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_spq_config), INTENT(inout) :: config    !< config type for spq
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    LOGICAL                     :: spq_deactivate_spq
    LOGICAL                     :: lrestart_cont
    LOGICAL                     :: flag_snow
    REAL(wp)                    :: wsr_capacity    , &
                                   wsn_capacity    , &
                                   spq_soil_depth
    REAL(wp)                    :: spq_soil_awc_prescribe, &
                                   spq_soil_theta_prescribe
    REAL(wp)                    :: spq_elevation
    REAL(wp)                    :: spq_soil_sand       , &
                                   spq_soil_silt       , &
                                   spq_soil_clay
    REAL(wp)                    :: spq_bulk_density
    INTEGER                     :: isoil   ! looping
    INTEGER                     :: i       ! looping over snow layers, consistent with jsb4
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(len=filename_max) :: bc_sso_filename
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename
    INTEGER                     :: nsnow
    REAL(wp)                    :: dz_snow(max_snow_layers)   ! width of each snow layer
    REAL(wp), ALLOCATABLE       :: depths(:), mids(:)         ! snow-layer vgrid init
    ! TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb              ! soil layers
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_snow_spq             ! snow layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_spq_config'

    NAMELIST /lnd_spq_nml/   &
                          active                  , &
                          spq_deactivate_spq      , &
                          lrestart_cont           , &
                          ic_filename             , &
                          bc_filename             , &
                          bc_sso_filename         , &
                          bc_quincy_soil_filename , &
                          flag_snow               , &
                          wsr_capacity            , &
                          wsn_capacity            , &
                          spq_soil_depth          , &
                          spq_soil_awc_prescribe  , &
                          spq_soil_theta_prescribe, &
                          spq_elevation           , &
                          spq_soil_sand           , &
                          spq_soil_silt           , &
                          spq_soil_clay           , &
                          spq_bulk_density

    ! variables for reading model-options from namelist
    INTEGER :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    IF (config%model_config%model_scheme == MODEL_QUINCY) THEN
      CALL message(TRIM(routine), 'Starting spq configuration')
    ELSE
      CALL message(TRIM(routine), 'NOT starting spq configuration - not running MODEL_QUINCY')
      RETURN
    END IF
    ! ----------------------------------------------------------------------------------------------------- !

    ! Set defaults
#ifdef __QUINCY_STANDALONE__
    active                        = .TRUE.   ! activated by default with quincy standalone
#else
    active                        = .FALSE.  ! de-activated by default with quincy in icon-land
#endif
    spq_deactivate_spq            = .FALSE.  ! may need to be consistent with namelist option "use_soil_phys_jsbach" from model config
    lrestart_cont                 = .FALSE.  ! TRUE: Continue although SPQ variables are missing in restart file
    ic_filename                   = 'ic_land_spq.nc'
    bc_filename                   = 'bc_land_spq.nc'
    bc_sso_filename               = 'bc_land_sso.nc'
    bc_quincy_soil_filename       = 'bc_quincy_soil.nc'
    flag_snow                     = .TRUE.
    wsr_capacity                  = 2.E-4_wp
    wsn_capacity                  = 2.E-4_wp
    spq_soil_depth                = 9.5_wp
    spq_soil_awc_prescribe        = 300.0_wp
    spq_soil_theta_prescribe      = 1.0_wp
    spq_elevation                 = 0.0_wp
    spq_soil_sand                 = 0.3_wp
    spq_soil_silt                 = 0.4_wp
    spq_soil_clay                 = 0.3_wp
    spq_bulk_density              = 1500.0_wp

    ! read the namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('lnd_spq_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_spq_nml)

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active                    = active
    config%spq_deactivate_spq        = spq_deactivate_spq
    config%lrestart_cont             = lrestart_cont
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%bc_sso_filename           = bc_sso_filename
    config%bc_quincy_soil_filename   = bc_quincy_soil_filename
    config%flag_snow                 = flag_snow
    config%soil_awc_prescribe        = spq_soil_awc_prescribe
    config%soil_theta_prescribe      = spq_soil_theta_prescribe
    config%elevation                 = spq_elevation
    config%wsr_capacity              = wsr_capacity
    config%wsn_capacity              = wsn_capacity
    config%soil_depth                = spq_soil_depth
    config%soil_sand                 = spq_soil_sand
    config%soil_silt                 = spq_soil_silt
    !config%soil_clay                = not using the input value, as soil_clay must be calculated from soil_sand & soil_silt to guarantee: 1=sand+silt+clay
    config%soil_clay                 = 1.0_wp - spq_soil_sand - spq_soil_silt
    config%bulk_density              = spq_bulk_density

#ifdef __QUINCY_STANDALONE__
    ! de-activate SPQ_ process
    ! is used consistently with namelist option "use_soil_phys_jsbach"
    IF (config%spq_deactivate_spq) THEN
      config%active = .FALSE.
    ELSE
      config%active = .TRUE.
    END IF
#endif

    ! check for an error in clay content of soil
    IF(config%soil_clay < eps8) THEN
       WRITE(message_text,'(a)') 'invalid proportion of soil_clay (value < eps8)'
       CALL finish(routine, message_text)
    END IF

    !< Create vertical snow-layer axis (vgrid) - consistent with jsb4
    !!
    !! these are "infrastructure" values \n
    !! layer 1 is the one at the ground \n
    !! snow_lay_thickness_snl(:,:) defines how much each of the layers if filled with snow
    !!
    nsnow             = 5       ! by default 5 snow layer (in jsb4 it is a namelist value, which is actually not used)
    dz_snow(:)        = 0._wp   ! jsb4 default
    dz_snow(1:nsnow)  = 0.05_wp ! jsb4 default
    IF (dz_snow(1) < snow_height_min) THEN
      CALL finish(TRIM(routine), 'Depth of first snow layer should be larger than snow_height_min=' &
      &                            //real2string(snow_height_min))
    END IF
    ALLOCATE(depths(nsnow+1))
    ALLOCATE(mids(nsnow))
    depths(1) = 0._wp
    DO i=1,nsnow
      depths(i+1) = depths(i) + dz_snow(i)
    END DO
    mids(1:nsnow) = (depths(1:nsnow) + depths(2:nsnow+1)) / 2._wp

    vgrid_snow_spq  => new_vgrid('snow_layer_spq', ZAXIS_GENERIC, nsnow,     &
      & levels  = mids                 (1:nsnow  ),                         &
      & lbounds = depths               (1:nsnow  ),                         &
      & ubounds = depths               (2:nsnow+1),                         &
      & units='m')
    CALL register_vgrid(vgrid_snow_spq)
    WRITE(message_text, *) 'Snow layer SPQ_ (upper)  [m]: ', vgrid_snow_spq%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow layer SPQ_ (mid)    [m]: ', vgrid_snow_spq%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow layer SPQ_ (lower)  [m]: ', vgrid_snow_spq%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow layer thickness SPQ_ [m]: ', vgrid_snow_spq%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

  END SUBROUTINE Init_spq_config

#endif
END MODULE mo_spq_config_class
