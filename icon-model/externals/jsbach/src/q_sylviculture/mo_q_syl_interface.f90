!> QUINCY sylviculture process interface
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
!>#### definition and init of tasks for the quincy sylviculture process
!>
MODULE mo_q_syl_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_math_constants,  ONLY: eps8
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_time,            ONLY: is_newyear, is_newday, get_year, get_year_length
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes Q_SYL_, VEG_
  dsl4jsb_Use_config(Q_SYL_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_memory(Q_SYL_)

#ifndef __QUINCY_STANDALONE__
  dsl4jsb_Use_processes PPLCC_
  dsl4jsb_Use_memory(PPLCC_)
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_q_syl_tasks
  PUBLIC ::  update_sylviculture

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_syl_interface'

  ! ======================================================================================================= !
  !> Type definition: sylviculture task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sylviculture
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sylviculture
    PROCEDURE, NOPASS :: Aggregate => aggregate_sylviculture
  END TYPE tsk_sylviculture

  !> Constructor interface: update_sylviculture task
  !>
  INTERFACE tsk_sylviculture
    PROCEDURE Create_task_sylviculture
  END INTERFACE tsk_sylviculture

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: Q_SYL_
  !>
  SUBROUTINE Register_q_syl_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_sylviculture(model_id))

  END SUBROUTINE Register_q_syl_tasks

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_sylviculture task
  !>
  FUNCTION Create_task_sylviculture(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_sylviculture::return_ptr)
    CALL return_ptr%Construct(name='sylviculture', process_id=Q_SYL_, owner_model_id=model_id)
  END FUNCTION Create_task_sylviculture

  ! ======================================================================================================= !
  !> Implementation of "update": sylviculture task
  !>
  SUBROUTINE update_sylviculture(tile, options)
    USE mo_q_syl_config_class,   ONLY: CONST_HARV
    USE mo_q_syl_constants,      ONLY: fract_wood_to_pp_fuel, fract_wood_to_pp_paper, fract_wood_to_pp_fiberboard,   &
      &                                fract_wood_to_pp_oirw, fract_wood_to_pp_pv, fract_wood_to_pp_sawnwood
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: box_tile
    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(Q_SYL_)

#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_memory_tile(PPLCC_, box_tile)
    dsl4jsb_Real2D_onChunk :: tree_fract_box_tile

    dsl4jsb_Def_memory_tile(Q_SYL_, box_tile)
    dsl4jsb_Real2D_onChunk :: fract_forest_harvest_y_read
#endif

    dsl4jsb_Def_memory(Q_SYL_)
    dsl4jsb_Real2D_onChunk :: fract_forest_harvest   ! fraction harvested from this tile in this timestep (rel to grid cell)

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                         :: ic, iblk, ics, ice, nc !< loop counter / dimensions
    REAL(wp)                        :: dtime                  !< timestep length
    INTEGER                         :: number_of_days         !< number of days in this year
    INTEGER                         :: current_year           !< the current year
    REAL(wp)                        :: fract_sum              !< sum of allocation fractions to different pools
    REAL(wp), DIMENSION(options%nc) :: cover_fraction         !< cover fraction of given tile
    REAL(wp), DIMENSION(options%nc) :: fact_harvest_share     !< factor to determine the harvest share of this tile
                                                              !< depends on cover fraction and total area with PFTs with trees
    REAL(wp), DIMENSION(options%nc) :: fract_harvest_box      !< total area relative to grid cell harvest from all forested PFTS
    TYPE(t_jsb_model), POINTER      :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_sylviculture'

    ! ----------------------------------------------------------------------------------------------------- !
    ! Get pointers to process configs and memory
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(Q_SYL_)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_SYL_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    nc    = options%nc
    ics   = options%ics
    ice   = options%ice
    iblk  = options%iblk
    dtime = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_SYL_)
    dsl4jsb_Get_var2D_onChunk(Q_SYL_, fract_forest_harvest)          ! out
    ! ----------------------------------------------------------------------------------------------------- !
    ! Since harvest only is happening daily or annually (or in QS even just once at the stand replacing year)
    ! the default would be that nothing is harvested in this timestep
    fract_forest_harvest(:) = 0._wp

    ! Check if some harvest is to be done or if we can RETURN here
    IF (dsl4jsb_Config(Q_SYL_)%flag_stand_harvest) THEN
      IF (dsl4jsb_Config(Q_SYL_)%flag_stand_harvest_event) THEN
        fract_forest_harvest(:) =  1._wp
      ELSE
        ! If we run QS and use stand replacing harvest nothing needs to be done
        IF (debug_on() .AND. iblk==1) CALL message(routine, 'Returned.')
        RETURN
      END IF
    ELSE IF (dsl4jsb_Config(Q_SYL_)%l_daily_harvest) THEN
      ! In case of daily harvest we have to check if this is a new day
      IF (.NOT. is_newday(options%current_datetime, dtime)) THEN
        ! If not then we do not have a harvest event and can return
        IF (debug_on() .AND. iblk==1) CALL message(routine, 'Returned.')
        RETURN
      END IF
    ELSE
      ! In case of annual harvest we have to check if this is a new year
      IF (.NOT. is_newyear(options%current_datetime, dtime)) THEN
        ! If not then we do not have a harvest event and can return
        IF (debug_on() .AND. iblk==1) CALL message(routine, 'Returned.')
        RETURN
      END IF
    END IF

    ! If this is not a stand replacing harvest event in which the whole tile is to be harvested,
    ! the to be harvested area still needs to be determined
    IF (.NOT. dsl4jsb_Config(Q_SYL_)%flag_stand_harvest) THEN
      ! Init harvested box area with the global constant value from the config (used in QS and in IQ in case of constant harvest)
      fract_harvest_box(:) = dsl4jsb_Config(Q_SYL_)%harvest_fraction

      ! Calculate the factor with which to derive the share to be harvested from this tile
      CALL tile%Get_fraction(ics, ice, iblk, fract=cover_fraction)
      fact_harvest_share(:) = cover_fraction(:)

#ifndef __QUINCY_STANDALONE__
      CALL model%Get_top_tile(box_tile)

      IF (.NOT. dsl4jsb_Config(Q_SYL_)%harvest_scheme == CONST_HARV) THEN
        dsl4jsb_Get_memory_tile(Q_SYL_, box_tile)
        dsl4jsb_Get_var2D_onChunk_tile(Q_SYL_, fract_forest_harvest_y_read, box_tile)
        fract_harvest_box(:) = fract_forest_harvest_y_read(:)
      END IF

      dsl4jsb_Get_memory_tile(PPLCC_, box_tile)
      dsl4jsb_Get_var2d_onChunk_tile_name(PPLCC_, tree_fract, box_tile) ! in

      DO ic = 1,nc
        IF (fract_harvest_box(ic) > tree_fract_box_tile(ic)) THEN
          IF (debug_on()) THEN
            WRITE(message_text,*) 'The to be harvested fraction (', fract_harvest_box(ic), &
              &                    ') exceeds the current tree area in the grid-cell, harvesting available area only...'
            CALL message(TRIM(routine),message_text)
          END IF
          fract_harvest_box(ic) = tree_fract_box_tile(ic) - eps8
        END IF

        IF (tree_fract_box_tile(ic) > eps8) THEN
          fact_harvest_share(ic) = cover_fraction(ic) / tree_fract_box_tile(ic)
        ELSE
          fact_harvest_share(ic) = 0.0_wp
        END IF
      END DO
#endif

      IF (dsl4jsb_Config(Q_SYL_)%l_daily_harvest) THEN
        ! In case of daily harvest the harvest needs to be distributed over the number of days in this year
        current_year  = get_year(options%current_datetime)
        number_of_days = get_year_length(current_year)
        fact_harvest_share(:) = fact_harvest_share(:) / REAL(number_of_days, wp)
      END IF

      ! Now the grid-cell harvest fraction for this tile can be determined
      fract_forest_harvest(:) = fract_harvest_box(:) * fact_harvest_share(:)

      ! Assertions
      IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        ! Assert: sum of allocation fractions to different product pools should sum up to 1._
        !         !NOTE: i.e. not including the slash fraction!
        fract_sum = fract_wood_to_pp_fuel + fract_wood_to_pp_paper + fract_wood_to_pp_fiberboard   &
          &         + fract_wood_to_pp_oirw + fract_wood_to_pp_pv + fract_wood_to_pp_sawnwood
        IF (ABS(fract_sum - 1._wp) >= eps8) THEN
          WRITE(message_text,*) 'A sum of allocation fractions to different pools does not sum up to one, but was: ', &
            &                    fract_sum,'! Please check!'
          CALL finish(TRIM(routine),message_text)
        END IF
      END IF
      DO ic = 1,nc
        ! Assert: harvest fraction should not be negative
        IF (fract_forest_harvest(ic) < 0._wp) THEN
          WRITE(message_text,*) 'The harvested fraction (', fract_forest_harvest(ic), &
            &                    ') should not be negative. Please check!'
          CALL finish(TRIM(routine),message_text)
        END IF
      END DO

      ! Assert: harvest fraction is not allowed to exceed the available cover fraction
      DO ic = 1,nc
        IF (fract_forest_harvest(ic) > cover_fraction(ic)) THEN
          IF (fract_forest_harvest(ic) - eps8 < cover_fraction(ic)) THEN
            fract_forest_harvest(ic) = cover_fraction(ic)
          ELSE
            WRITE(message_text,*) 'The harvested fraction (', fract_forest_harvest(ic), &
              & ') is not allowed to exceed the actually available fraction in this tile (', &
              & cover_fraction(ic),' in ', TRIM(tile%name),'). Please check!'
            CALL finish(TRIM(routine),message_text)
          END IF
        END IF
      END DO
    END IF ! .NOT. dsl4jsb_Config(Q_SYL_)%flag_stand_harvest

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_sylviculture

  ! ======================================================================================================= !
  !> Implementation of "aggregate": sylviculture task
  !>
  SUBROUTINE aggregate_sylviculture(tile, options)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_SYL_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sylviculture'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_SYL_)
    dsl4jsb_Aggregate_onChunk(Q_SYL_, fract_forest_harvest, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_sylviculture

#endif
END MODULE mo_q_syl_interface
