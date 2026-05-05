!> QUINCY vegetation process interface
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
!>#### definition and init of tasks for the vegetation process incl. update and aggregate routines
!>
!> includes plant growth and turnover
!>
MODULE mo_q_veg_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_grid_class,              ONLY: t_jsb_grid
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process

  ! "Integrate" routines
  USE mo_q_veg_turnover,                ONLY: update_veg_turnover_real => update_veg_turnover
  USE mo_q_veg_dynamics,                ONLY: update_veg_dynamics_real => update_veg_dynamics
  USE mo_q_veg_growth,                  ONLY: update_veg_growth_real => update_veg_growth
  USE mo_q_veg_update_pools,            ONLY: update_veg_pools_real => update_veg_pools
  USE mo_q_veg_plant_characteristics,   ONLY: update_plant_characteristics_real => update_plant_characteristics
  USE mo_q_veg_products_decay,          ONLY: update_products_decay_real => update_products_decay
  USE mo_q_veg_update_pools_on_harvest, ONLY: update_pools_on_harvest_real => update_pools_on_harvest

  dsl4jsb_Use_processes VEG_, SB_, Q_ASSIMI_
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_memory(SB_)
  dsl4jsb_Use_memory(Q_ASSIMI_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_veg_tasks_quincy
  PUBLIC :: global_q_veg_diagnostics



  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_veg_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_veg_fluxes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_veg_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_turnover task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_turnover
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_turnover     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_turnover  !< Aggregates computed task variables
  END TYPE tsk_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_dynamics task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_dynamics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_dynamics     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_dynamics  !< Aggregates computed task variables
  END TYPE tsk_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_growth task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_growth
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_growth     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_growth  !< Aggregates computed task variables
  END TYPE tsk_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_pools task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_pools      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_pools   !< Aggregates computed task variables
  END TYPE tsk_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_plant_characteristics task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_plant_characteristics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_plant_characteristics      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_plant_characteristics   !< Aggregates computed task variables
  END TYPE tsk_plant_characteristics

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: time_average_vegetation task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_vegetation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_vegetation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_vegetation !< Aggregates computed task variables
  END TYPE tsk_time_average_vegetation

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_products_decay task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_products_decay
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_products_decay    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_products_decay !< Aggregates computed task variables
  END TYPE tsk_products_decay

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_pools_on_harvest task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_pools_on_harvest
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_pools_on_harvest    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_pools_on_harvest !< Aggregates computed task variables
  END TYPE tsk_pools_on_harvest

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_veg_fluxes
    PROCEDURE Create_task_reset_veg_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_turnover task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_turnover
    PROCEDURE Create_task_update_veg_turnover         !< Constructor function for task
  END INTERFACE tsk_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_dynamics task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_dynamics
    PROCEDURE Create_task_update_veg_dynamics         !< Constructor function for task
  END INTERFACE tsk_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_growth task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_growth
    PROCEDURE Create_task_update_veg_growth         !< Constructor function for task
  END INTERFACE tsk_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_pools task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_pools
    PROCEDURE Create_task_update_veg_pools         !< Constructor function for task
  END INTERFACE tsk_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_plant_characteristics task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_plant_characteristics
    PROCEDURE Create_task_update_plant_characteristics         !< Constructor function for task
  END INTERFACE tsk_plant_characteristics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_time_average_vegetation task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_time_average_vegetation
    PROCEDURE Create_task_update_time_average_vegetation         !< Constructor function for task
  END INTERFACE tsk_time_average_vegetation

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_products_decay task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_products_decay
    PROCEDURE Create_task_update_products_decay         !< Constructor function for task
  END INTERFACE tsk_products_decay

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_pools_on_harvest task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_pools_on_harvest
    PROCEDURE Create_task_update_pools_on_harvest         !< Constructor function for task
  END INTERFACE tsk_pools_on_harvest

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_veg_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: VEG_ process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_veg_tasks_quincy(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_veg_fluxes(model_id))
    CALL this%Register_task(tsk_veg_turnover(model_id))
    CALL this%Register_task(tsk_veg_dynamics(model_id))
    CALL this%Register_task(tsk_veg_growth(model_id))
    CALL this%Register_task(tsk_veg_pools(model_id))
    CALL this%Register_task(tsk_time_average_vegetation(model_id))
    CALL this%Register_task(tsk_products_decay(model_id))
    CALL this%Register_task(tsk_plant_characteristics(model_id))
    CALL this%Register_task(tsk_pools_on_harvest(model_id))

  END SUBROUTINE Register_veg_tasks_quincy

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_veg_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_veg_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_veg_fluxes', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_veg_fluxes


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_turnover task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_turnover(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_turnover::return_ptr)
    CALL return_ptr%Construct(name='veg_turnover', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_dynamics task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_dynamics(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_dynamics::return_ptr)
    CALL return_ptr%Construct(name='veg_dynamics', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_growth task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_growth(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_growth::return_ptr)
    CALL return_ptr%Construct(name='veg_growth', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_pools task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_pools(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_pools::return_ptr)
    CALL return_ptr%Construct(name='veg_pools', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_plant_characteristics task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_plant_characteristics(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_plant_characteristics::return_ptr)
    CALL return_ptr%Construct(name='plant_characteristics', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_plant_characteristics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: tsk_time_average_vegetation task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_time_average_vegetation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_vegetation::return_ptr)
    CALL return_ptr%Construct(name='tavrg_vegetation', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_time_average_vegetation

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_products_decay task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_products_decay(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_products_decay::return_ptr)
    CALL return_ptr%Construct(name='products_decay', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_products_decay

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_pools_on_harvest task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_pools_on_harvest(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_pools_on_harvest::return_ptr)
    CALL return_ptr%Construct(name='pools_on_harvest', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_pools_on_harvest

  !> Implementation of "update": reset_veg_fluxes task
  !!
  SUBROUTINE update_reset_veg_fluxes(tile, options)

    USE mo_veg_util, ONLY: reset_veg_fluxes

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_veg_fluxes'

    iblk    = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL reset_veg_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_reset_veg_fluxes

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrappers for update routines from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_veg_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_pools_real(tile, options)

  END SUBROUTINE update_veg_pools

  SUBROUTINE update_plant_characteristics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_plant_characteristics_real(tile, options)

  END SUBROUTINE update_plant_characteristics

  SUBROUTINE update_veg_growth(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_growth_real(tile, options)

  END SUBROUTINE update_veg_growth

  SUBROUTINE update_veg_dynamics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_dynamics_real(tile, options)

  END SUBROUTINE update_veg_dynamics

  SUBROUTINE update_veg_turnover(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_turnover_real(tile, options)

  END SUBROUTINE update_veg_turnover

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrappers for update routines from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_products_decay(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_products_decay_real(tile, options)

  END SUBROUTINE update_products_decay

  SUBROUTINE update_pools_on_harvest(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_pools_on_harvest_real(tile, options)

  END SUBROUTINE update_pools_on_harvest

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_veg_fluxes task
  !>
  SUBROUTINE aggregate_reset_veg_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_veg_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_veg_fluxes


  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_turnover task
  !>
  SUBROUTINE aggregate_veg_turnover(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_turnover'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_veg_turnover

  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_dynamics task
  !>
  SUBROUTINE aggregate_veg_dynamics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_dynamics'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_veg_dynamics

  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_growth task
  !>
  SUBROUTINE aggregate_veg_growth(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_growth'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_veg_growth

  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_pools task
  !>
  SUBROUTINE aggregate_veg_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_config(VEG_)

    ! C fluxes of vegetation
    dsl4jsb_Aggregate_onChunk(VEG_, npp                                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_biosphere_production            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_c                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_c              , weighted_by_fract)

    ! C isotopes needed for global analysis (not all fluxes required!)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_c13                             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_biosphere_production_c13        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_c14                             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_biosphere_production_c14        , weighted_by_fract)

    ! N fluxes
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation                          , weighted_by_fract)
    ! needed here only because otherwise it would not be aggregated properly in PLANT mode
    dsl4jsb_Aggregate_onChunk(VEG_, biological_n_fixation               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_nh4                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_no3                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_n                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_n              , weighted_by_fract)

    ! P fluxes
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_po4                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_p                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_p              , weighted_by_fract)

    ! Fluxes due to product pool decay if product pool is present
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_decay_c              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_decay_n              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_decay_p              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_decay_c13            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_decay_c14            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_decay_n15            , weighted_by_fract)
    END IF

    ! C fluxes from herbivory if herbivory is calculated
    IF (dsl4jsb_Config(VEG_)%flag_herbivory) THEN
      dsl4jsb_Aggregate_onChunk(VEG_, herbivory_leaf_resp               , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, herbivory_fruit_resp              , weighted_by_fract)
    END IF

    ! VEG_ bgcm pools needed to close global mass balance
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_c            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_n            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_p            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_c13          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_c14          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_n15          , weighted_by_fract)
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_c      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_n      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_p      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_c13    , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_c14    , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_n15    , weighted_by_fract)
    END IF
    ! special diagnostic output for specific pools
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_leaf_c             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_leaf_n             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_leaf_p             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_wood_c             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_wood_n             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_wood_p             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_fine_root_c        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_fine_root_n        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_fine_root_p        , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_veg_pools

  ! ======================================================================================================= !
  !>Implementation of "aggregate": plant_characteristics task
  !>
  SUBROUTINE aggregate_plant_characteristics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_config(VEG_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, diameter                            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, height                              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, dens_ind                            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, lai                                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fract_fpc                           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, rough_veg_star                      , weighted_by_fract)
    ! VEG_ 3D
    dsl4jsb_Aggregate_onChunk(VEG_, root_fraction_sl            , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_plant_characteristics

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for VEG_
  !>
  SUBROUTINE update_time_average_vegetation(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_veg_util,            ONLY: calculate_time_average_vegetation
    dsl4jsb_Use_processes VEG_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_vegetation'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL calculate_time_average_vegetation(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_vegetation


  ! ======================================================================================================= !
  !> Implementation of "aggregate": time_average_vegetation task
  !>
  SUBROUTINE aggregate_time_average_vegetation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_vegetation'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_vegetation

  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for products decay task
  !
  SUBROUTINE aggregate_products_decay(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_products_decay'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    !TODO: implement -- only required with IQ
    !dsl4jsb_Aggregate_onChunk(VEG_, XXX , weighted_by_fract)

  END SUBROUTINE aggregate_products_decay


  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for pools on harvest
  !
  SUBROUTINE aggregate_pools_on_harvest(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_products_decay'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    !TODO: implement -- only required with IQ
    !dsl4jsb_Aggregate_onChunk(VEG_, XXX , weighted_by_fract)

  END SUBROUTINE aggregate_pools_on_harvest

  ! ======================================================================================================= !
  !>
  !> calculations of vegetation diagnostic global sums
  !>
  SUBROUTINE global_q_veg_diagnostics(tile)

    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N
    USE mo_jsb_math_constants,      ONLY: one_day, one_year
    USE mo_sync,                    ONLY: global_sum_array
    USE mo_jsb_grid,                ONLY: Get_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER       :: total_veg_c_gsum(:)
    REAL(wp), POINTER       :: total_prod_c_gsum(:)
    REAL(wp), POINTER       :: lai_gmean(:)
    REAL(wp), POINTER       :: fract_fpc_gmean(:)
    REAL(wp), POINTER       :: gpp_gsum(:)
    REAL(wp), POINTER       :: npp_gsum(:)
    REAL(wp), POINTER       :: nuptake_gsum(:)
    REAL(wp)                :: global_land_area

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_q_veg_diagnostics'

    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_ASSIMI_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: veg_pool_total_c
    dsl4jsb_Real2D_onDomain :: veg_products_total_c
    dsl4jsb_Real2D_onDomain :: lai
    dsl4jsb_Real2D_onDomain :: fract_fpc
    dsl4jsb_Real2D_onDomain :: npp
    dsl4jsb_Real2D_onDomain :: uptake_nh4
    dsl4jsb_Real2D_onDomain :: uptake_no3
    dsl4jsb_Real2D_onDomain :: gross_assimilation

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_config(VEG_)

    dsl4jsb_Get_var2D_onDomain(VEG_,  veg_pool_total_c)                ! in
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      dsl4jsb_Get_var2D_onDomain(VEG_,  veg_products_total_c)          ! in
    END IF
    dsl4jsb_Get_var2D_onDomain(VEG_,      lai)                         ! in
    dsl4jsb_Get_var2D_onDomain(VEG_,      fract_fpc)                   ! in
    dsl4jsb_Get_var2D_onDomain(VEG_,      npp)                         ! in
    dsl4jsb_Get_var2D_onDomain(VEG_,      uptake_nh4)                  ! in
    dsl4jsb_Get_var2D_onDomain(VEG_,      uptake_no3)                  ! in
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, gross_assimilation)          ! in

    total_veg_c_gsum   => VEG__mem%total_veg_c_gsum%ptr(:)             ! out
    total_prod_c_gsum  => VEG__mem%total_prod_c_gsum%ptr(:)            ! out
    lai_gmean          => VEG__mem%lai_gmean%ptr(:)                    ! out
    fract_fpc_gmean    => VEG__mem%fract_fpc_gmean%ptr(:)              ! out
    gpp_gsum           => VEG__mem%gpp_gsum%ptr(:)                     ! out
    npp_gsum           => VEG__mem%npp_gsum%ptr(:)                     ! out
    nuptake_gsum       => VEG__mem%nuptake_gsum%ptr(:)                 ! out

    ! Domain Mask - to mask all halo cells for global sums (otherwise these cells are counted twice)
    ALLOCATE(in_domain(grid%nproma,grid%nblks))
    ALLOCATE(scaling(grid%nproma,grid%nblks))

    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ! Calculate global carbon inventories, if requested for output
    !  => Conversion from [mol(C)/m^2] to [PgC]
    !     1 mol C = molar_massC g C   => 1 mol C = molar_mass_C * e-15 Pg C
    scaling(:,:) = molar_mass_C * 1.e-15_wp * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (VEG__mem%total_veg_c_gsum%is_in_output)        &
          &  total_veg_c_gsum        = global_sum_array(veg_pool_total_c(:,:)        * scaling(:,:))
    IF (VEG__mem%total_prod_c_gsum%is_in_output)        &
          &  total_prod_c_gsum        = global_sum_array(veg_products_total_c(:,:)   * scaling(:,:))

    ! Calculate global land means variables, if requested for output
    global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))
    scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (VEG__mem%lai_gmean%is_in_output)     &
      &  lai_gmean     = global_sum_array(lai(:,:)    * scaling(:,:)) / global_land_area
    IF (VEG__mem%fract_fpc_gmean%is_in_output)  &
      &  fract_fpc_gmean  = global_sum_array(fract_fpc(:,:) * scaling(:,:)) / global_land_area

    ! Calculate global BGC fluxes, if requested for output
    !  => Conversion from [ummol/m^2/s] to [PgC/yr]
    !     1 mu mol C = molar_massC g C e-6  => 1 Pg = molar_mass_C * e-15 * e-6
    !     1 s times one_day (60*60*24) time one_year (365)
    scaling(:,:) = molar_mass_C * 1.e-21_wp * (one_day*one_year) * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (VEG__mem%gpp_gsum%is_in_output)        &
          &  gpp_gsum        = global_sum_array(gross_assimilation(:,:) * scaling(:,:))
    IF (VEG__mem%npp_gsum%is_in_output)        &
          &  npp_gsum        = global_sum_array(npp(:,:) * scaling(:,:))
    !  => Conversion from [ummol/m^2/s] to [TgN/yr]
    !     1 mu mol N = molar_massN g N e-6  => 1 Tg = molar_mass_N * e-12 * e-6
    !     1 s times one_day (60*60*24) time one_year (365)
    scaling(:,:) = molar_mass_N * 1.e-18_wp * (one_day*one_year) * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (VEG__mem%nuptake_gsum%is_in_output)        &
          &  nuptake_gsum        = global_sum_array((uptake_nh4(:,:)+uptake_no3(:,:)) * scaling(:,:))

    DEALLOCATE (scaling, in_domain)

  END SUBROUTINE global_q_veg_diagnostics


#endif
END MODULE mo_q_veg_interface
