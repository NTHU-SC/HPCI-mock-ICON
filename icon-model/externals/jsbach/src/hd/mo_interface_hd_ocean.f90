!> Interface between HD and the ocean, through a coupler
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
MODULE mo_interface_hd_ocean
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_coupling_utils,      ONLY: cpl_def_cell_field_mask, &
                                    cpl_def_field, &
                                    cpl_put_field
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_control,        ONLY: acc_stream

  dsl4jsb_Use_processes HYDRO_, HD_
  dsl4jsb_Use_config(HD_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(HD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsb_fdef_hd_fields
  PUBLIC :: interface_hd_ocean

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_hd_ocean'

  INTEGER :: field_id_river_runoff

CONTAINS

  SUBROUTINE jsb_fdef_hd_fields(comp_id, cell_point_id, grid_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid
    USE mo_jsb_grid,           ONLY: Get_grid

    INTEGER, INTENT(in) :: comp_id       ! yac component id
    INTEGER, INTENT(in) :: cell_point_id ! yac point is
    INTEGER, INTENT(in) :: grid_id       ! yac grid id
    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    INTEGER :: cell_mask_id
    INTEGER :: model_id
    INTEGER :: jb, jc
    INTEGER :: nproma, nblks
    LOGICAL,  ALLOCATABLE :: is_valid(:)
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    CHARACTER(LEN=*), PARAMETER :: routine = modname // ':jsb_fdef_hd_fields'

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    dsl4jsb_Real2D_onDomain :: hd_mask

    ! It is assumed that only one domain is active (jg=1) i.e. no nesting.
    ! Currently, this requirement is fulfilled since ocean coupling is only set up for one domain.
    model_id = 1

    model => Get_model(model_id)

    dsl4jsb_Get_config(HD_)
    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    CALL model%Get_top_tile(tile)
    dsl4jsb_Get_memory(HD_)

    grid  => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    ! Define the HD mask for YAC.
    ! It shall contain ocean coast points only for
    ! source point mapping (source_to_target_map)

    ALLOCATE(is_valid(nproma*nblks))

    is_valid(:) = .FALSE.

    dsl4jsb_Get_var2D_onDomain(HD_, hd_mask)

    DO jb = 1, nblks
      DO jc = 1, nproma
        IF ( hd_mask(jc, jb) .EQ. 0.0 ) THEN
          is_valid((jb-1)*nproma+jc) = .TRUE.
        END IF
      END DO
    END DO

    CALL cpl_def_cell_field_mask(routine, grid_id, is_valid, cell_mask_id)

    DEALLOCATE(is_valid)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'river_runoff', 1, field_id_river_runoff)

  END SUBROUTINE jsb_fdef_hd_fields

  SUBROUTINE interface_hd_ocean(tile, options)

    USE mo_fortran_tools, ONLY: init

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':interface_hd_ocean'

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain ::discharge_ocean

    ! Local variables

    LOGICAL :: write_coupler_restart
    INTEGER :: n

    n =options%nc ! only to avoid compiler warnings

    IF ( .NOT. is_coupled_run() ) RETURN

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HD_)
    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(routine, 'HD model works on root tile only!')

    grid => Get_grid(model%grid_id)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge_ocean) ! IN

    ! Send river discharge
    ! --------------------

    !$ACC UPDATE HOST(discharge_ocean)

    CALL cpl_put_field( &
      & routine, field_id_river_runoff, 'river_runoff', grid%ntotal, &
      & discharge_ocean, write_restart=write_coupler_restart)

    IF ( write_coupler_restart ) &
      & CALL message(routine, 'YAC says it is put for restart')

  END SUBROUTINE interface_hd_ocean
#endif

END MODULE mo_interface_hd_ocean
