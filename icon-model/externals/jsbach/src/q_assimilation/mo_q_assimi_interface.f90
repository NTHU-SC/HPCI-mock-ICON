!> QUINCY assimilation process interface
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
!>#### definition and init of tasks for the assimilation process (photosynthesis) incl. update and aggregate routines
!>
MODULE mo_q_assimi_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes Q_ASSIMI_
  dsl4jsb_Use_config(Q_ASSIMI_)
  dsl4jsb_Use_memory(Q_ASSIMI_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_q_assimi_tasks
  PUBLIC ::  update_canopy_fluxes
  PUBLIC ::  update_ftranspiration_per_sl   ! used only with quincy + jsbach soil physics processes
  PUBLIC ::  update_reset_q_assimi_fluxes, update_time_average_q_assimilation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_assimi_interface'

  ! ======================================================================================================= !
  !> Type definition: reset_q_assimi_fluxes task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_q_assimi_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_q_assimi_fluxes     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_q_assimi_fluxes  !< Aggregates computed task variables
  END TYPE tsk_reset_q_assimi_fluxes
  !> Constructor interface: reset_q_assimi_fluxes task
  !>
  INTERFACE tsk_reset_q_assimi_fluxes
    PROCEDURE Create_task_reset_q_assimi_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_q_assimi_fluxes

  ! ======================================================================================================= !
  !> Type definition: canopy_fluxes task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_fluxes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_fluxes !< Aggregates computed task variables
  END TYPE tsk_canopy_fluxes
  !> Constructor interface: update_canopy_fluxes task
  !>
  INTERFACE tsk_canopy_fluxes
    PROCEDURE Create_task_canopy_fluxes         !< Constructor function for task
  END INTERFACE tsk_canopy_fluxes

  ! ======================================================================================================= !
  !> Type definition: ftranspiration_per_sl task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_ftranspiration_per_sl
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_ftranspiration_per_sl    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_ftranspiration_per_sl !< Aggregates computed task variables
  END TYPE tsk_ftranspiration_per_sl
  !> Constructor interface: update_ftranspiration_per_sl task
  !>
  INTERFACE tsk_ftranspiration_per_sl
    PROCEDURE Create_task_ftranspiration_per_sl         !< Constructor function for task
  END INTERFACE tsk_ftranspiration_per_sl

  ! ======================================================================================================= !
  !> Type definition: time_average_assimilation task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_assimilation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_q_assimilation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_q_assimilation !< Aggregates computed task variables
  END TYPE tsk_time_average_assimilation
  !> Constructor interface: time_average_assimilation task
  !>
  INTERFACE tsk_time_average_assimilation
    PROCEDURE Create_task_time_average_assimilation         !< Constructor function for task
  END INTERFACE tsk_time_average_assimilation

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_q_assimi_fluxes task
  !>
  FUNCTION Create_task_reset_q_assimi_fluxes(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_q_assimi_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_q_assimi_fluxes', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_reset_q_assimi_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_canopy_fluxes task
  !>
  FUNCTION Create_task_canopy_fluxes(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_fluxes::return_ptr)
    CALL return_ptr%Construct(name='canopy_fluxes', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_canopy_fluxes

  ! ======================================================================================================= !
  !> Constructor: update_ftranspiration_per_sl task
  !>
  FUNCTION Create_task_ftranspiration_per_sl(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_ftranspiration_per_sl::return_ptr)
    CALL return_ptr%Construct(name='ftranspiration_per_sl', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_ftranspiration_per_sl

  ! ======================================================================================================= !
  !> Constructor: tsk_time_average_assimilation task
  !>
  FUNCTION Create_task_time_average_assimilation(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_assimilation::return_ptr)
    CALL return_ptr%Construct(name='tavrg_assimilation', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_time_average_assimilation

  ! ======================================================================================================= !
  !> Register tasks for assimi process
  !>
  SUBROUTINE Register_q_assimi_tasks(this, model_id)
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER :: model

    model => Get_model(model_id)

    CALL this%Register_task(tsk_reset_q_assimi_fluxes(model_id))
    CALL this%Register_task(tsk_canopy_fluxes(model_id))
    CALL this%Register_task(tsk_ftranspiration_per_sl(model_id))
    CALL this%Register_task(tsk_time_average_assimilation(model_id))
  END SUBROUTINE Register_q_assimi_tasks

  ! ======================================================================================================= !
  !> Implementation of "update": reset_q_assimi_fluxes task
  !>
  SUBROUTINE update_reset_q_assimi_fluxes(tile, options)
    USE mo_q_assimi_util, ONLY: reset_q_assimi_fluxes
    dsl4jsb_Use_processes Q_ASSIMI_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_q_assimi_fluxes'

    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    iblk = options%iblk

    IF (debug_on()) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL reset_q_assimi_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_reset_q_assimi_fluxes

  ! ======================================================================================================= !
  !> Implementation of "aggregate": reset_q_assimi_fluxes task
  !>
  SUBROUTINE aggregate_reset_q_assimi_fluxes(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_q_assimi_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)

    ! Q_ASSIMI_ 2D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C13                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C14                , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_reset_q_assimi_fluxes

  ! ======================================================================================================= !
  !> Implementation of "update": canopy_fluxes task
  !>
  SUBROUTINE update_canopy_fluxes(tile, options)
    USE mo_q_assimi_process, ONLY: real_update_canopy_fluxes => update_canopy_fluxes

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_fluxes'

    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    iblk = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_update_canopy_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_canopy_fluxes

  ! ======================================================================================================= !
  !>Implementation of "aggregate": canopy_fluxes task
  !>
  !>  plus: correction of ftranspiration_sl after aggregation
  !>
  SUBROUTINE aggregate_canopy_fluxes(tile, options)

    USE mo_jsb_math_constants,     ONLY: eps8

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model
    CLASS(t_jsb_aggregator),  POINTER :: weighted_by_fract
    REAL(wp)                          :: ftranspiration_sum(options%nc)           !< SUM(ftranspiration_sl) across soil layers
    INTEGER                           :: iblk, nc, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    nc      = options%nc
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)

    ! Q_ASSIMI_ 2D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C13      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C14      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, canopy_cond                 , weighted_by_fract)
    IF (model%config%use_soil_phys_jsbach) THEN
      ! only with jsbach soil physics processes enabled, otherwise calculated with SPQ_
      dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, wtr_soil_root_pot           , weighted_by_fract)
    END IF
    ! Q_ASSIMI_ 3D
    IF (.NOT. model%config%use_soil_phys_jsbach) THEN
      ! only with SPQ_ calculated in update_canopy_fluxes
      dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, ftranspiration_sl           , weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_canopy_fluxes

  ! ======================================================================================================= !
  !>Implementation of "update": ftranspiration_per_sl task
  !>
  !>  called only with QUINCY and jsbach soil physics processes (not with SPQ_)
  !>
  SUBROUTINE update_ftranspiration_per_sl(tile, options)
    USE mo_q_assimi_process, ONLY: real_update_ftranspiration_per_sl => update_ftranspiration_per_sl

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_ftranspiration_per_sl'

    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    iblk = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_update_ftranspiration_per_sl(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_ftranspiration_per_sl

  ! ======================================================================================================= !
  !>Implementation of "aggregate": ftranspiration_per_sl task
  !>
  !>  plus: correction of ftranspiration_sl after aggregation
  !>
  !>  called only with QUINCY and jsbach soil physics processes (not with SPQ_)
  !>
  SUBROUTINE aggregate_ftranspiration_per_sl(tile, options)
    USE mo_jsb_math_constants,     ONLY: eps8

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model
    CLASS(t_jsb_aggregator),  POINTER :: weighted_by_fract
    REAL(wp)                          :: ftranspiration_sum(options%nc)           !< SUM(ftranspiration_sl) across soil layers
    REAL(wp)                          :: ftranspiration_sum_mult(options%nc)      !< SUM(ftranspiration_sl) across soil layers
                                                                                  !! multiplied by 1.E6_wp to avoide precision issues at some grid cells
    INTEGER                           :: iblk, ics, ice
    INTEGER                           :: nc, ic, is
    INTEGER                           :: nsoil
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_ftranspiration_per_sl'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onChunk :: ftranspiration_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    nc      = options%nc
    ics     = options%ics
    ice     = options%ice
    nsoil   = options%nsoil_w
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, ftranspiration_sl)     ! inout

    ! Q_ASSIMI_ 3D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, ftranspiration_sl           , weighted_by_fract)

    ! correction to ensure ftranspiration_sl sums up to 1.0 across soil layers
    !
    ! necessary because the aggregation from PFT to upper tiles
    ! does not ensure that the sum across the soil profile is 1.0
    !
    ! TODO  improve Aggregate_weighted_by_fract_3d() routine
    !
    !$ACC DATA ASYNC(1) CREATE(ftranspiration_sum, ftranspiration_sum_mult)
    ! init ftranspiration_sum with zero
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ftranspiration_sum(ic)      = 0.0_wp
      ftranspiration_sum_mult(ic) = 0.0_wp
    END DO
    !$ACC END PARALLEL LOOP
    ! calc ftranspiration_sum across soil layers
    !
    ! use ftranspiration_sum_mult(:,:) to deal with very low values due to very large fractions of the bare tile
    !    and, hence, very low fractions of vegetated tiles of a grid cell to avoid water balance issues
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        ftranspiration_sum(ic)      = ftranspiration_sum(ic)      + ftranspiration_sl(ic, is)
        ftranspiration_sum_mult(ic) = ftranspiration_sum_mult(ic) + ftranspiration_sl(ic, is) * 1.E6_wp
      END DO
    END DO
    !$ACC END PARALLEL
    ! scale ftranspiration_sl values to sum up to 1.0 across soil layers
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (ftranspiration_sum(ic) > eps8) THEN
          ftranspiration_sl(ic, is) = ftranspiration_sl(ic, is) / ftranspiration_sum(ic)
        ELSE IF (ftranspiration_sum_mult(ic) > eps8) THEN
          ftranspiration_sl(ic, is) = ftranspiration_sl(ic, is) * 1.E6_wp / ftranspiration_sum_mult(ic)
        ELSE
        ! this may lead to a detectable mass balance error in the water balance
        ! the code cannot deal with fractions of a bare tile larger than 0.9999 (E-14)
          IF (is == 1) THEN
            ftranspiration_sl(ic, is) = 1.0_wp
          ELSE
            ftranspiration_sl(ic, is) = 0.0_wp
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_ftranspiration_per_sl

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for Q_ASSIMI_
  !>
  SUBROUTINE update_time_average_q_assimilation(tile, options)
    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_q_assimi_util,         ONLY: calculate_time_average_q_assimilation
    dsl4jsb_Use_processes Q_ASSIMI_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_q_assimilation'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) THEN
      CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    END IF

    CALL calculate_time_average_q_assimilation(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_time_average_q_assimilation

  ! ======================================================================================================= !
  !> Implementation of "aggregate": time_average_assimilation task
  !>
  SUBROUTINE aggregate_time_average_q_assimilation(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_q_assimilation'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_time_average_q_assimilation

#endif
END MODULE mo_q_assimi_interface
