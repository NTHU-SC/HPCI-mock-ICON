!> alcc config - defines vgrid required in alcc memory
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
!>#### Could contain namelist info (e.g. alcc scheme), but currently of no use beyond infrastructural requirements
!>
MODULE mo_alcc_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_string,        ONLY: tolower
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_kind,               ONLY: wp
  USE mo_jsb_config_class,   ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_alcc_config

  TYPE, EXTENDS(t_jsb_config)      :: t_alcc_config
     CHARACTER(len=filename_max-7) :: alcc_filename_prefix
       !< Prefix of alcc filename prefix with len filename_max-7 (to add yyyy.nc)
     CHARACTER(len=SHORT_NAME_LEN) :: scheme
     LOGICAL                       :: l_daily_alcc
       !< Determines if land use change is to be applied on a daily (TRUE) or annual (FALSE) basis.
       ! Note: This also determines the year of the land use maps file to read in: in case of "annual" the file for
       ! the current year is read, while in case of "daily" the next year's file is read as target land cover map.
   CONTAINS
     PROCEDURE :: Init => Init_alcc_config
  END type t_alcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize alcc process
  !
  ! -------------------------------------------------------------------------------------------------------
  SUBROUTINE Init_alcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_alcc_config), INTENT(inout) :: config !<Configuration type of process (t_alcc_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                       :: active, lrestart_cont
    CHARACTER(len=filename_max)   :: ic_filename, bc_filename
    CHARACTER(len=filename_max-7) :: alcc_filename_prefix
    CHARACTER(len=SHORT_NAME_LEN) :: scheme
    LOGICAL                       :: l_daily_alcc
    TYPE(t_jsb_vgrid), POINTER    :: pft_layer
    INTEGER                       :: i

    NAMELIST /jsb_alcc_nml/      &
         active,                 &
         lrestart_cont,          &
         ic_filename,            &
         bc_filename,            &
         alcc_filename_prefix,   &
         scheme,                 &
         l_daily_alcc

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_alcc_config'
    ! -------------------------------------------------------------------------------------------------- !

    CALL message(TRIM(routine), 'Starting alcc configuration')

    ! Set defaults
    active               = .FALSE.
    lrestart_cont        = .FALSE.       ! TRUE: Continue although ALCC variables are missing in restart file
    bc_filename          = 'bc_land_alcc.nc'
    ic_filename          = 'ic_land_alcc.nc'
    alcc_filename_prefix = 'bc_land_frac_'
    scheme               = 'maps'
    l_daily_alcc         = .TRUE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_alcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_alcc_nml)

    CALL close_nml(nml_handler)

    config%active               = active
    config%lrestart_cont        = lrestart_cont
    config%ic_filename          = ic_filename
    config%bc_filename          = bc_filename
    config%alcc_filename_prefix = alcc_filename_prefix
    config%scheme               = scheme
    config%l_daily_alcc         = l_daily_alcc

    IF (config%active) THEN
      SELECT CASE (TRIM(tolower(config%scheme)))
      CASE ('maps')
        CALL message(TRIM(routine), 'Reading land use data from maps.')
        config%scheme=tolower(config%scheme)

        IF (l_daily_alcc) THEN
          CALL message(TRIM(routine), '... land use change will be applied on a daily basis, reading next years land cover map' &
            &                       //' as target land cover to reach by the end of the year.')
        ELSE
          CALL message(TRIM(routine), '... land use change will be applied on an annual basis,' &
            &                       //' reading the land cover map of the current year to be applied in the first timestep.')
        END IF

      CASE DEFAULT
        CALL finish(TRIM(routine), 'Not implemented alcc scheme specified: '// TRIM(config%scheme))
      END SELECT
    END IF

  END SUBROUTINE Init_alcc_config

#endif
END MODULE mo_alcc_config_class
