!> QUINCY sylviculture process config
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
!>#### define sylviculture config structure, read sylviculture namelist and init configuration parameters
!>
MODULE mo_q_syl_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_math_constants,ONLY: eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_q_syl_config, CONST_HARV, REF_HARV, TRANS_HARV

  ! ======================================================================================================= !
  !> Q_SYL_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_syl_config
    CHARACTER(len=filename_max) :: harvest_filename_prefix !< IQ. Note: should NOT include an ending "_"

    LOGICAL   :: l_daily_harvest   !< determines if harvest is to be applied on a daily (TRUE) or annual (FALSE) basis
    INTEGER   :: harvest_scheme    !< one of: CONST_HARV=1, REF_HARV=2, TRANS_HARV=3
                                   !< - in case of REF_HARV and TRANS_HARV the data is read from an external file (IQ only)
    REAL(wp)  :: harvest_fraction  !< default fraction used if not read from a map (annual change, relative to whole gridcell)
      !< for stand replacing harvest events, in QS, this is ignored and will be set to 1.0 instead

    ! QS specific configurations
    LOGICAL   :: harvest_active_in_qs     !< QS: if harvest should be active in QS (required for QS run script environment)
    LOGICAL   :: flag_stand_harvest       !< QS: if stand replacing harvest is enabled
    INTEGER   :: stand_replacing_year     !< QS: site specific stand replacing harvest year
    LOGICAL   :: flag_stand_harvest_event !< QS: if a site specific stand replacing harvest happened in this timestep
  CONTAINS
    PROCEDURE :: Init => Init_q_syl_config
  END TYPE t_q_syl_config

  !> harvest schemes (possible values for harvest_scheme)
  ENUM, BIND(C)
    ENUMERATOR ::       &
      & CONST_HARV = 1, &  !< use global constant grid harvest fraction
      & REF_HARV,       &  !< initially read one map with grid harvest fractions to be used each year (IQ)
      & TRANS_HARV         !< read a map with grid harvest fractions in each simulation year (IQ)
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_syl_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read Q_SYL_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_q_syl_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_util_string,           ONLY: tolower
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_GENERIC
    USE mo_jsb_math_constants,    ONLY: eps4
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_syl_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    LOGICAL                     :: lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(len=filename_max) :: harvest_filename_prefix
    CHARACTER(9)                :: harvest_scheme           !< one of constant, reference, transient

    LOGICAL   :: l_daily_harvest          !< determines if harvest is to be applied on a daily (TRUE) or annual (FALSE) basis
    REAL(wp)  :: harvest_fraction         !< default fraction used if not read from a map (annual change, relative to whole gridcell)
      !< for stand replacing harvest events, in QS, this is ignored and will be set to 1.0 instead

    ! QS specific configurations
    LOGICAL   :: harvest_active_in_qs     !< QS: if harvest should be active in QS (required for QS run script environment)
    LOGICAL   :: flag_stand_harvest       !< QS: if stand replacing harvest is enabled
    INTEGER   :: stand_replacing_year     !< QS: site specific stand replacing harvest year

    NAMELIST /lnd_q_syl_nml/  &

      & active,                   &
      & lrestart_cont,            &
      & ic_filename,              &
      & bc_filename,              &
      & harvest_filename_prefix,  &
      & l_daily_harvest,          &
      & harvest_scheme,           &
      & harvest_fraction,         &
      & harvest_active_in_qs,     &
      & flag_stand_harvest,       &
      & stand_replacing_year


    INTEGER                     :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_syl_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting sh configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                        = .FALSE.
    lrestart_cont                 = .TRUE.          ! TRUE: Continue even if Q_SYL_ variables are missing in restart file
    bc_filename                   = 'bc_land_use_static.nc'
    ic_filename                   = 'ic_land_syl.nc'
    !
    harvest_filename_prefix       = 'bc_land_use_quincy'
    harvest_scheme                = 'none'  ! possible values: none, constant, reference or transient
    l_daily_harvest               = .TRUE.  ! determines if harvest is to be applied on a daily (TRUE) or annual (FALSE) basis
    ! for now arbitrary default fraction used if not read from a map (annual change, relative to whole gridcell)
    harvest_fraction              = 0.05_wp

    ! QS specific configurations
    harvest_active_in_qs          = .FALSE. ! (required for QS run script environment)
    flag_stand_harvest            = .FALSE.
    stand_replacing_year          = 1500
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('lnd_q_syl_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_syl_nml)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                    = active
    config%lrestart_cont             = lrestart_cont
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%harvest_filename_prefix   = harvest_filename_prefix
    config%l_daily_harvest           = l_daily_harvest
    config%harvest_fraction          = harvest_fraction
    config%harvest_active_in_qs      = harvest_active_in_qs
    config%flag_stand_harvest        = flag_stand_harvest
    config%stand_replacing_year      = stand_replacing_year

    ! QS only: Test if we have a stand replacing harvest event (default: FALSE) -- does not change in IQ!
    config%flag_stand_harvest_event  = .FALSE.

    CALL close_nml(nml_handler)

    ! Set harvest scheme
    SELECT CASE (tolower(TRIM(harvest_scheme)))
    CASE ("none")
      CALL message(TRIM(routine), 'Running with process active but no global harvest.')
      config%harvest_scheme = 0.0_wp
    CASE ("constant")
      CALL message(TRIM(routine), 'Running with constant global harvest.')
      config%harvest_scheme = CONST_HARV
    CASE ("reference")
      CALL message(TRIM(routine), 'Running with harvest from a reference file.')
      config%harvest_scheme = REF_HARV
    CASE ("transient")
      CALL message(TRIM(routine), 'Running with transient harvest.')
      config%harvest_scheme = TRANS_HARV
    CASE DEFAULT
      CALL finish(TRIM(routine), 'harvest_scheme == '//tolower(TRIM(harvest_scheme))//' not available.')
    END SELECT

  END SUBROUTINE Init_q_syl_config

#endif
END MODULE mo_q_syl_config_class
