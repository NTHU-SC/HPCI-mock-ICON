!> Interface to the alcc process (anthropogenic land cover change)
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
!>#### Contains the interfaces to the alcc process, determining area changes and moving to-be-conserved matter
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_alcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,       ONLY: debug_on
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_util,              ONLY: one_of

  USE mo_jsb_model_class,    ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
!  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

#ifndef __NO_QUINCY__
  USE mo_quincy_model_config, ONLY: QLAND, QPLANT
#endif

  ! Use of processes in this module
  dsl4jsb_Use_processes ALCC_, PPLCC_, FAGE_, NLCC_

  ! Use of process configurations
  dsl4jsb_Use_config(ALCC_)

  ! Use of process memories
  dsl4jsb_Use_memory(ALCC_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_alcc_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_interface'

  ! -------------------------------------------------------------------------------------------------------
  !> Type definition: task to determine ALCC cover fraction changes
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_alcc_delta_cover_fractions
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_alcc_delta_cover_fractions     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_alcc_delta_cover_fractions  !< Aggregates computed task variables
  END TYPE tsk_alcc_delta_cover_fractions

  !> Constructor interface for ALCC task deriving cover fraction changes
  INTERFACE tsk_alcc_delta_cover_fractions
    PROCEDURE Create_task_alcc_delta_cover_fractions                        !< Constructor function for task
  END INTERFACE tsk_alcc_delta_cover_fractions

  ! -------------------------------------------------------------------------------------------------------
  !> Type definition: alcc translocation task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_alcc_translocation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_alcc_translocation     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_alcc_translocation  !< Aggregates computed task variables
  END TYPE tsk_alcc_translocation

  !> Constructor interface for ALCC translocation task
  INTERFACE tsk_alcc_translocation
    PROCEDURE Create_task_alcc_translocation                        !< Constructor function for task
  END INTERFACE tsk_alcc_translocation

CONTAINS

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for ALCC task task deriving cover fraction changes
  !>
  FUNCTION Create_task_alcc_delta_cover_fractions(model_id) RESULT(return_ptr)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "alcc_delta_cover_fractions"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_alcc_delta_cover_fractions::return_ptr)
    CALL return_ptr%Construct(name='alcc_delta_cover_fractions', process_id=ALCC_, owner_model_id=model_id)
  END FUNCTION Create_task_alcc_delta_cover_fractions

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for ALCC translocation task
  !>
  FUNCTION Create_task_alcc_translocation(model_id) RESULT(return_ptr)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "alcc_translocation"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_alcc_translocation::return_ptr)
    CALL return_ptr%Construct(name='alcc_translocation', process_id=ALCC_, owner_model_id=model_id)
  END FUNCTION Create_task_alcc_translocation

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for ALCC process
  !>
  SUBROUTINE Register_alcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of alcc process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_alcc_delta_cover_fractions(model_id))
    CALL this%Register_task(tsk_alcc_translocation(model_id))

  END SUBROUTINE Register_alcc_tasks


  ! ====================================================================================================== !
  !>
  !> #### Update ALCC cover fraction changes
  !>
  !> The area fraction changes for each PFT tile are derived from prescribed anthropogenic land cover change (ALCC) forcing, read
  !> in [[mo_alcc_init:read_land_use_data]].
  !> Updating of cover fraction changes is executed at the beginning of each year.
  !> The calculation of land cover fractions and the used forcing depends on the configuration (annual or daily changes).
  !> Note: while this lcc process changes the PFT/age class tiles, it runs and is controlled on their parent tile (the VEG tile).
  !>
  SUBROUTINE update_alcc_delta_cover_fractions(tile, options)

    USE mo_util,              ONLY: one_of
    USE mo_jsb_time,          ONLY: is_newyear, get_year, get_year_length
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc, min_daily_cf_change, min_annual_cf_change

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_)                !< Memory of the anthropogenic land use change process
    dsl4jsb_Def_config(ALCC_)                !< Configuration of the anthropogenic land use change process
    dsl4jsb_Real3D_onChunk :: cf_target_year !< Target cover fractions derived from annually read land use data []
    dsl4jsb_Real3D_onChunk :: cf_delta       !< Change in cover fractions derived from cf_target_year (annual or daily change) []

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile !< Pointer to the tile for which fraction changes are determined
    TYPE(t_jsb_model), POINTER          :: model        !< Current instance of the model

    REAL(wp), DIMENSION(options%nc) :: current_fract    !< Current cover fraction

    INTEGER  :: &
      & iblk, &         !< Current block index
      & ics, &          !< Index of first cell of block
      & ice, &          !< Index of last cell of block
      & nc, &           !< Number of cells in current block
      & ic, &           !< Cell index
      & i_cf_tile, &    !< Index used for forest age classes
      & current_year, & !< Index of a PFT tile in the array with cover fractions changes
      & number_of_days  !< Number of tiles participating in the ALCC process

    REAL(wp) :: delta_threshold !< Threshold below which cover fraction changes are not conducted
    REAL(wp) :: dtime   !< Time step length

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc_delta_cover_fractions'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    current_year  = get_year(options%current_datetime)
    number_of_days = get_year_length(current_year)

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_target_year)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_delta)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(alcc_)) RETURN

    ! Assert that PPLCC process is active
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
    END IF

    ! Cover fraction changes only need to be calculated once at the beginning of a new year
    IF (.NOT. is_newyear(options%current_datetime,dtime)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    ! IF (dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
    !   delta_threshold = min_daily_cf_change
    ! ELSE
    !   delta_threshold = min_annual_cf_change
    ! END IF
    delta_threshold = EPSILON(1._wp)

    ! Calculate cover fraction changes
    current_tile => tile%Get_first_child_tile()
    i_cf_tile = 0

    DO WHILE (ASSOCIATED(current_tile))
      i_cf_tile = i_cf_tile + 1

      CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract)
      ! Init the boolean to false which tracks if a tile will change due to alcc on a grid-cell
      current_tile%l_fract_alcc_change(ics:ice,iblk) = .FALSE.

      ! Derive area changes
      DO ic = 1,nc
        IF (cf_target_year(ic, i_cf_tile) < 0.0_wp) THEN
          ! This should only happen upon init with -0.1 / i.e. if not read from a file
          cf_delta(ic, i_cf_tile) = 0.0_wp
        ELSE
          cf_delta(ic, i_cf_tile) = cf_target_year(ic, i_cf_tile) - current_fract(ic)
        END IF

        IF (dsl4jsb_Config(ALCC_)%l_daily_alcc .AND. ABS(cf_delta(ic, i_cf_tile)) >= 1.e-15_wp) THEN
          ! In case that area changes are daily the difference has to be distributed over the year
          ! such that the target is reached at the end of the year
          cf_delta(ic, i_cf_tile) = cf_delta(ic, i_cf_tile) / REAL(number_of_days, wp)
        END IF

!          IF (ABS(cf_delta(ic, i_cf_tile)) >= EPSILON(1._wp)) THEN
        IF (ABS(cf_delta(ic, i_cf_tile)) >= delta_threshold) THEN
          ! If there are area changes above the threshold set l_fract_alcc_change flag (just required for QUINCY)
          current_tile%l_fract_alcc_change(ics+ic-1,iblk) = .TRUE.
        ELSE
          ! else set to 0.0_wp
          cf_delta(ic, i_cf_tile) = 0.0_wp
        END IF

      END DO

      current_tile => current_tile%Get_next_sibling_tile()
    END DO ! current_tile

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_alcc_delta_cover_fractions


  ! ====================================================================================================== !
  !>
  !> #### Aggregation of variables from calculation of ALCC cover fraction change -- currently not used.
  !>
  SUBROUTINE aggregate_alcc_delta_cover_fractions(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_alcc_delta_cover_fractions'
    ! -------------------------------------------------------------------------------------------------- !

  END SUBROUTINE aggregate_alcc_delta_cover_fractions


  ! ====================================================================================================== !
  !
  !> #### Calculate daily/annual relocation of fractions and matter
  !>
  !> This task applies the actual fraction changes as derived in [[update_alcc_delta_cover_fractions]] for each PFT/ or age class,
  !> as well as the corresponding translocation of matter from tiles loosing area to tiles gaining area. If running with natural
  !> land cover changes (NLCC) dynamic pfts are handled separately [[handle_dynamic_pfts]].
  !>
  !> Note:
  !>
  !>   - While this lcc process changes the PFT/age class tiles, it runs and is controlled on their parent tile (the VEG tile).
  !>   - The ALCC process requires that the PPLCC process is active, too (translocated matter is assumed to be relative to the box).
  !>
  SUBROUTINE update_alcc_translocation(tile, options)

    USE mo_jsb_time,          ONLY: is_newday, is_newyear
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc, max_tolerated_cf_mismatch
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_)          !< Memory of the anthropogenic land use change process
    dsl4jsb_Def_config(ALCC_)          !< Configuration of the anthropogenic land use change process
    dsl4jsb_Real3D_onChunk :: cf_delta !< Derived change in cover fractions (annual or daily change) []

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile    !< Pointer to the tile for which relocations are determined
    CLASS(t_jsb_tile_abstract), POINTER :: age_class_tile  !< Pointer to for age class tiles
    TYPE(t_jsb_model), POINTER          :: model           !< Current instance of the model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations !< lcc relocations instance of the ALCC process

    REAL(wp), POINTER ::   &
      & initial_area(:,:), & !< Pointer to initial area in lcc relocations instance of the ALCC process
      & lost_area(:,:),    & !< Pointer to lost area in lcc relocations instance of the ALCC process
      & gained_area(:,:)     !< Pointer to gained area in lcc relocations instance of the ALCC process

    REAL(wp), DIMENSION(options%nc) :: &
      & cf_diff,      & !< Change in cover fractions
      & current_fract   !< Current cover fraction

    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice, &        !< Index of last cell of block
      & nc, &         !< Number of cells in current block
      & ic, &         !< Cell index
      & i_tile, &     !< Index of a tile in the lcc_relocations
      & i_ac_index, & !< Index used for forest age classes
      & i_cf_tile, &  !< Index of a PFT tile in the array with cover fractions changes
      & nr_of_tiles   !< Number of tiles participating in the ALCC process
    REAL(wp) :: dtime !< Time step length

    LOGICAL :: dynamic_tile !< Boolean indicating if a tile is a dynamic tile (participating in NLCC; only if NLCC is active)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc_translocation'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_delta)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(ALCC_)) RETURN

    ! Assert that PPLCC process is active
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
    END IF

    ! Redistributions are either only conducted at the start of each year or each day
    IF (.NOT. dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
      IF (.NOT. is_newyear(options%current_datetime,dtime)) RETURN
    ELSE
      IF (.NOT. is_newday(options%current_datetime,dtime)) RETURN
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    ! 1. Get lcc structure
    dsl4jsb_Get_lcc_relocations(ALCC_, lcc_relocations)
    nr_of_tiles = lcc_relocations%nr_of_tiles

    ! Init area vectors
    dsl4jsb_Get_lcc_area_matrix(initial_area)
    dsl4jsb_Get_lcc_area_matrix(lost_area)
    dsl4jsb_Get_lcc_area_matrix(gained_area)
    initial_area(:,:) = 0.0_wp
    lost_area(:,:) = 0.0_wp
    gained_area(:,:) = 0.0_wp

    ! Collect current area of all involved tiles
    CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)

    ! 2. Calculate area changes
    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    i_cf_tile = 0

    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      i_cf_tile = i_cf_tile + 1

      ! Find out if there are dynamic PFT tiles, whose fractions should not be replaced by the map's fractions.
      dynamic_tile = .FALSE.
      IF (model%Is_process_enabled(NLCC_) .AND. dsl4jsb_Lctlib_param_tile(dynamic_PFT, current_tile)) THEN
        dynamic_tile = .TRUE.
      END IF

      ! Area changes need to be determined differently for forest age classes than for the pft usecase
      IF ((TRIM(model%config%usecase) == 'jsbach_forest_age_classes') .AND. (current_tile%Has_children())) THEN

        ! Assertion: age-classes are subsequent parts of the lcc structure
        age_class_tile => current_tile%Get_first_child_tile()
        i_ac_index = i_tile
        DO WHILE (ASSOCIATED(age_class_tile))
          IF (.NOT. TRIM(age_class_tile%name) .EQ. TRIM(lcc_relocations%tile_names(i_ac_index))) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: age class (' //TRIM(age_class_tile%name) &
              & //') is not listed in the expected place in the alcc lcc-structure. Found ' &
              & // TRIM(lcc_relocations%tile_names(i_ac_index)) // ' instead. Please check!')
          END IF
          age_class_tile => age_class_tile%Get_next_sibling_tile()
          i_ac_index = i_ac_index + 1
        END DO

        ! Forest pfts are not part of the lcc structure (i.e. not involved tiles)
        ! and thus the initial area has not yet been automatically collected above
        CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract)

        ! The area loss or gain of a forest pft needs to be redistributed to its children (the age classes)
        CALL redistribute_cover_fraction_changes_of_forest_pft( &
            & current_tile, options, cf_delta(:, i_cf_tile), i_tile, gained_area, lost_area)
        i_tile = i_tile + current_tile%Get_no_of_children() - 1

        ! Finally, update the fraction on the forest pft tile itself
        current_fract(:) = current_fract(:) + cf_delta(:, i_cf_tile)
        CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))
      ELSE
        ! Pft without children (either pft usecase or non forest pfts in fage usecase!)
        ! in the pft usecase i_cf_tile equals i_tile, however, this is not the case for non forest pfts in fage usecase
        ! therefore cf_delta is are used with pft index (i_cf_tile), while the lcc arrays (initial_area, lost_area and gained_area)
        ! are used with the current tile index (i_tile)

        ! Assert: pft without children
        IF (current_tile%Has_children()) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: alcc only expects children of a lcc tile '  &
            & //'in case of the forest age class usecase, but usecase was '// TRIM(model%config%usecase)  &
            & // '. Tile with children: ' // trim(current_tile%name) // '. Please check.')
        END IF

        ! Assert: In mo_alcc_init read_land_use_data assumes that all child tiles of the veg tile
        ! are part of the lcc structure
        IF (one_of(current_tile%name, lcc_relocations%tile_names) <= 0) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: child tile of veg tile(' //TRIM(current_tile%name) &
          & //') is not part of the alcc lcc-structure. Current implementation assumes that all pfts are part of it.')
        END IF

        DO ic = 1,nc
          cf_diff(ic) = cf_delta(ic, i_cf_tile)

          ! Determine new cover fractions, while keeping within the boundaries valid for cover fractions [0,1]
          IF ((cf_diff(ic) < 0.0_wp) .AND. (ABS(cf_diff(ic)) >= initial_area(ic, i_tile))) THEN
            ! If area is lost it should not exceed the initial area
            current_fract(ic) = 0.0_wp
            cf_diff(ic) = -1.0_wp * initial_area(ic, i_tile)
          ELSE IF ((cf_diff(ic) > 0.0_wp) .AND. ((cf_diff(ic) + initial_area(ic, i_tile)) >= 1.0_wp)) THEN
            ! If area is gained the sum of initial area and gained area should not exceed 1 (total available area)
            current_fract(ic) = 1.0_wp
            cf_diff(ic) = 1.0_wp - initial_area(ic, i_tile)
          ELSE
            current_fract(ic) = initial_area(ic, i_tile) + cf_diff(ic)
          END IF

          ! Derive gained and lost area
          IF (cf_diff(ic) >= 0.0_wp) THEN
            gained_area(ic, i_tile) = cf_diff(ic)
          ELSE
            ! The lcc framework expects also lost area to be a positive value, therefore the sign is swapped here
            lost_area(ic, i_tile) = -1.0_wp * cf_diff(ic)
          END IF
        END DO

        IF (.NOT. dynamic_tile) THEN
          CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))
        END IF
      END IF ! If running with jsbach_forest_age_classes usecase and a forest pft or ELSE a pft without leaves as current tile

      current_tile => current_tile%Get_next_sibling_tile()
    END DO ! Current_tile

    ! If demanded: handle the dynamic PFT tiles
    IF (model%Is_process_enabled(NLCC_)) THEN
      CALL handle_dynamic_pfts(tile, options, initial_area, gained_area, lost_area)
    END IF

    ! Assertion: the total gained area needs to equal the total lost area to enable matter conservation
    cf_diff(:) = SUM(gained_area(:,:),DIM=2) - SUM(lost_area(:,:),DIM=2)
    ! Check for violations outside of the tolerated range (which are probably rather due to errors than precision)
    IF (ANY(ABS(cf_diff(:)) > max_tolerated_cf_mismatch)) THEN
      WRITE (message_text,*) 'Violation of assertion: gained area does not equal lost area! Please check input data. ' &
        & // 'Max difference of: ', MAXVAL(cf_diff(:))
      CALL finish(TRIM(routine), TRIM(message_text))
    END IF
    IF (model%config%check_child_cf) THEN
      ! Check and correct remaining mismatches
      CALL check_and_correct_cf_mismatches(tile, options, nr_of_tiles, gained_area, lost_area)
    END IF

    ! Now do the matter translocation
    IF (model%config%model_scheme == MODEL_JSBACH) THEN
      CALL alcc_matter_translocation_jsb(tile, options)
#ifndef __NO_QUINCY__
    ELSE IF (model%config%model_scheme == MODEL_QUINCY) THEN
      IF (model%config%qmodel_id == QLAND) THEN
        CALL alcc_matter_translocation_iq(tile, options)
      END IF
#endif
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_alcc_translocation


  ! ====================================================================================================== !
  !>
  !> #### Aggregation of variables from calculation of ALCC translocation
  !>
  SUBROUTINE aggregate_alcc_translocation(tile, options)

    dsl4jsb_Use_memory(FAGE_)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice           !< Index of last cell of block

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions
    dsl4jsb_Def_memory(FAGE_)                              !< Memory of the forest age classes process
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_alcc_translocation'
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (tile%Is_process_active(FAGE_)) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      dsl4jsb_Get_memory(FAGE_)
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_alcc_translocation


  ! ====================================================================================================== !
  !
  !> Do the matter translocation for jsbach.
  !
  SUBROUTINE alcc_matter_translocation_jsb(tile, options)

    USE mo_jsb_time,          ONLY: is_newday, is_newyear
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc
    USE mo_jsb_lcc,           ONLY: start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_) !< Memory of the anthropogenic land use change process
    dsl4jsb_Def_config(ALCC_) !< Configuration of the anthropogenic land use change process

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile    !< Pointer to the tile for which relocations are applied
    CLASS(t_jsb_tile_abstract), POINTER :: forest_pft_tile !< Pointer for forest pft tiles
    TYPE(t_jsb_model), POINTER          :: model           !< Current instance of the model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations !< lcc relocations instance of the ALCC process

    REAL(wp), POINTER ::   &
      & initial_area(:,:), & !< Pointer to initial area in lcc relocations instance of the ALCC process
      & lost_area(:,:),    & !< Pointer to lost area in lcc relocations instance of the ALCC process
      & gained_area(:,:)     !< Pointer to gained area in lcc relocations instance of the ALCC process

    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice, &        !< Index of last cell of block
      & nc, &         !< Number of cells in current block
      & ic, &         !< Cell index
      & i_tile        !< Index of a tile in the lcc_relocations
    REAL(wp) :: dtime !< Time step length

    LOGICAL  :: is_age_class !< Boolean indicating if a tile is a forest age class

    CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_matter_translocation_jsb'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_memory(ALCC_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(ALCC_)) RETURN

    !>
    !> Assert that PPLCC process is active
    !>
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
    END IF

    ! Redistributions are either only conducted at the start of each year or each day
    IF (.NOT. dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
      IF (.NOT. is_newyear(options%current_datetime,dtime)) RETURN
    ELSE
      IF (.NOT. is_newday(options%current_datetime,dtime)) RETURN
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    !>
    !> 1. Get lcc structure
    !>
    dsl4jsb_Get_lcc_relocations(ALCC_, lcc_relocations)

    ! Get area vectors
    dsl4jsb_Get_lcc_area_matrix(initial_area)
    dsl4jsb_Get_lcc_area_matrix(lost_area)
    dsl4jsb_Get_lcc_area_matrix(gained_area)

    !>
    !> 2. Collect to be transferred matter
    !>
    CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    !>
    !> 3. Transfer matter from active to passive vars
    !>
    current_tile => tile%Get_first_child_tile()
    is_age_class = .FALSE.
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))

      ! In case of forest age classes we have to decent one level further
      IF ((TRIM(model%config%usecase) == 'jsbach_forest_age_classes') .AND. (current_tile%Has_children())) THEN
        is_age_class = .TRUE.
        forest_pft_tile => current_tile
        current_tile => current_tile%Get_first_child_tile()
      END IF

      IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
        i_tile = i_tile + 1
        IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
          CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
        END IF
      END IF

      current_tile => current_tile%Get_next_sibling_tile()

      ! In case of forest age classes we have to go up one level and go to the next tile there
      IF (is_age_class .AND. .NOT. ASSOCIATED(current_tile)) THEN
        current_tile => forest_pft_tile
        current_tile => current_tile%Get_next_sibling_tile()
        is_age_class = .FALSE.
      END IF
    END DO

    !>
    !> 4. make the passive matter transfer
    !>
    CALL end_lcc_reloc(lcc_relocations, options, gained_area)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE alcc_matter_translocation_jsb


#ifndef __NO_QUINCY__
  ! ====================================================================================================== !
  !
  !> Do the matter translocation for quincy.
  !> Notes:
  !>
  !>   - Due to the complexity of the bgcms and because of the call to calc_litter_patritioning the lcc framework cannot be used
  !>     to translocate active to passive CQs. Instead the active translocation is done explicitly and the framework is only used
  !>     for the passive relocation!
  !>   - To reduce the number of bgcm variables and thereby the size of the memory two existing flux bgcms are used here in a
  !>     temporal manner: the harvest litter flux and the harvest flux to product pools. After usage they are reset to zero.
  !>     If these fluxes would be needed for diagnostical purposses additional bgcms would be required.
  !
  SUBROUTINE alcc_matter_translocation_iq(tile, options)
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid

    USE mo_jsb_lcc_class,         ONLY: t_jsb_lcc_proc
    USE mo_jsb_lcc,               ONLY: start_lcc_reloc, end_lcc_reloc

    USE mo_lnd_bgcm_idx,          ONLY: ix_sap_wood, ix_heart_wood, LAST_ELEM_ID
    USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,  ONLY: VEG_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, VEG_BGCM_PP_FUEL_ID, &
        &                               VEG_BGCM_HARVEST_LITTER_ID, VEG_BGCM_HARVEST_TO_PROD_ID

    USE mo_q_sb_litter_processes, ONLY: calc_litter_partitioning

    dsl4jsb_Use_processes VEG_, HYDRO_, SB_, Q_SYL_
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(Q_SYL_)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)                                !< Configuration of the soil beigeochemistry process
    dsl4jsb_Def_config(VEG_)                               !< Configuration of the vegetation process

    dsl4jsb_Def_memory_tile(Q_SYL_, box_tile)              !< Memory of the sylviculture process holding the slash fraction
    dsl4jsb_Real2D_onChunk :: fract_wood_to_slash

    dsl4jsb_Def_memory_tile(VEG_, current_tile)            !< Memory of the vegetation process on a tile
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: dens_ind
    dsl4jsb_Real2D_onChunk :: delta_dens_ind
    dsl4jsb_Real3D_onChunk :: root_fraction_sl

    dsl4jsb_Def_memory_tile(HYDRO_, current_tile)          !< Memory of the hydrology process on a tile loosing area
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: soil_depth_sl

    CLASS(t_jsb_tile_abstract), POINTER :: box_tile        !< Pointer to the box tile to get the harvest slash fraction
    CLASS(t_jsb_tile_abstract), POINTER :: current_tile    !< Pointer to the tile for which relocations are applied
    TYPE(t_jsb_model),          POINTER :: model           !< Current instance of the model
    TYPE(t_jsb_lcc_proc),       POINTER :: lcc_relocations !< lcc relocations instance of the ALCC process
    TYPE(t_lnd_bgcm_store),     POINTER :: bgcm_store      !< The bgcm store of the tile for which relocations are applied
    TYPE(t_lctlib_element),     POINTER :: lctlib          !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_soil_w    !< Vertical grid

    REAL(wp)                            :: zero_mt(LAST_ELEM_ID, options%nc) !< helper array to pass zero seedbed litter

    REAL(wp), POINTER ::   &
      & initial_area(:,:), & !< Pointer to initial area in lcc relocations instance of the ALCC process
      & lost_area(:,:),    & !< Pointer to lost area in lcc relocations instance of the ALCC process
      & gained_area(:,:)     !< Pointer to gained area in lcc relocations instance of the ALCC process

    REAL(wp) :: scaling !< scaling of new to old cover fraction required to calculate the density change in extending forest tiles
    REAL(wp) :: dtime   !< Time step length

    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice, &        !< Index of last cell of block
      & nc, &         !< Number of cells in current block
      & ic, &         !< Cell index
      & nsoil_w, &    !< Number of soil layers (water)
      & i_tile        !< Index of a tile in the lcc_relocations

    dsl4jsb_Def_mt2L2D :: veg_litter_flux_harvest_mt !< Litter flux used temporarily to move matter from veg pool to formation flux
    dsl4jsb_Def_mt2L2D :: veg_pool_mt                !< Explicily relocated to formation (and product pool flux, if enabled)
    dsl4jsb_Def_mt2L3D :: sb_formation_mt            !< Formation flux filled with litter from veg pool
    dsl4jsb_Def_mt1L2D :: veg_pp_flux_harvest_mt     !< Harvest flux to product pools (synchronised and distributed)
    dsl4jsb_Def_mt1L2D :: veg_pp_fuel_mt             !< Fuel product pool (filled on extending tiles after translocation)

    CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_matter_translocation_iq'
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels

    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !>
    !> 1. Get lcc structure
    !>
    dsl4jsb_Get_lcc_relocations(ALCC_, lcc_relocations)

    ! get area vectors
    dsl4jsb_Get_lcc_area_matrix(initial_area)
    dsl4jsb_Get_lcc_area_matrix(lost_area)
    dsl4jsb_Get_lcc_area_matrix(gained_area)

    !>
    !> 2. Matter translocation for "active" variables needs to be done explicitly for quincy (if bgcm matrices are involved)
    !>
    zero_mt(:,:) = 0.0_wp

    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1

      IF (.NOT. ANY(current_tile%name == lcc_relocations%tile_names)) THEN
        CALL finish(TRIM(routine), 'Violation of precondition: for quincy all pft tiles are expected to be involved with ALCC,' &
          & //' not found: '// trim(current_tile%name))
      END IF

      ! Get lctlib of this tile
      lctlib => model%lctlib(current_tile%lcts(1)%lib_id)

      ! Get required vegetation variables for this tile
      dsl4jsb_Get_memory_tile(VEG_, current_tile)
      dsl4jsb_Get_var2D_onChunk_tile(VEG_, dens_ind, current_tile)
      dsl4jsb_Get_var2D_onChunk_tile(VEG_, delta_dens_ind, current_tile)
      dsl4jsb_Get_var3D_onChunk_tile(VEG_, root_fraction_sl, current_tile)

      ! Relocation of "active" variables only happens on tiles with grid cells that are loosing area
      IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
        ! ---------------------------------------------------------------------------------------------------!
        ! For tiles with grid cells that loose area also hydrology variables are requried
        dsl4jsb_Get_memory_tile(HYDRO_, current_tile)
        dsl4jsb_Get_var2D_onChunk_tile(HYDRO_, num_sl_above_bedrock, current_tile)
        dsl4jsb_Get_var3D_onChunk_tile(HYDRO_, soil_depth_sl, current_tile)

        ! Get involved bgcm matrices for this tile
        bgcm_store => current_tile%bgcm_store
        dsl4jsb_Get_mt2L2D_tile(VEG_BGCM_HARVEST_LITTER_ID, veg_litter_flux_harvest_mt, current_tile)
        dsl4jsb_Get_mt2L2D_tile(VEG_BGCM_POOL_ID, veg_pool_mt, current_tile)
        dsl4jsb_Get_mt2L3D_tile(SB_BGCM_FORMATION_ID, sb_formation_mt, current_tile)

        IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
          dsl4jsb_Get_mt1L2D_tile(VEG_BGCM_HARVEST_TO_PROD_ID, veg_pp_flux_harvest_mt, current_tile)
        END IF

        ! ---------------------------------------------------------------------------------------------------!
        ! - Transfer cleared amount of each element of each compartment of the vegetation to litter flux (including seeds)
        !   Note: matter translocation tasks of land cover change processes are executed in-between the two pplcc tasks, in which
        !         all variables that contain conserved matter are in m2 relative to the box area and not the tile area.
        !         (Converted in the pre-lcc task of the pplcc process.)
        DO ic = 1,nc
          veg_litter_flux_harvest_mt(:,:,ic) = veg_pool_mt(:,:,ic) * lost_area(ic, i_tile)
          veg_pool_mt(:,:,ic) = veg_pool_mt(:,:,ic) - veg_litter_flux_harvest_mt(:,:,ic)
        END DO

        IF (dsl4jsb_Config(VEG_)%l_use_product_pools .AND. (dsl4jsb_Lctlib_param_tile(ForestFlag, current_tile))) THEN
          CALL model%Get_top_tile(box_tile)
          dsl4jsb_Get_memory_tile(Q_SYL_, box_tile)
          dsl4jsb_Get_var2D_onChunk_tile(Q_SYL_, fract_wood_to_slash, box_tile)

          DO ic = 1,nc
            ! sapwood: lctlib_frac_sapwood_branch goes to litter, rest to product pool -> but mind: slash fraction
            veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
              & + (veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
              &    * (1._wp - dsl4jsb_Lctlib_param_tile(frac_sapwood_branch, current_tile)) * (1._wp - fract_wood_to_slash(ic)))
            veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) = veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
              & - (veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
              &    * (1._wp - dsl4jsb_Lctlib_param_tile(frac_sapwood_branch, current_tile)) * (1._wp - fract_wood_to_slash(ic)))

            ! heartwood: goes to product pool -> but mind: slash fraction
            veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
              & + veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) * (1._wp - fract_wood_to_slash(ic))
            veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) = veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) &
              & - veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) * (1._wp - fract_wood_to_slash(ic))
          END DO
        END IF

        ! ---------------------------------------------------------------------------------------------------!
        ! After potential reduction due to product usage, the alcc litter flux is put into the sb formation flux
        ! which in turn is subject to passive relocation
        ! NOTE: if the matter relocation for quincy should be called in a loop (e.g. transitions), an additional
        !       alcc formation flux variable would be required which would need to be treated similarly to the product pool flux
        !       variable (veg_pp_flux_harvest_mt), i.e. it would need to start empty, be filled in calc_litter_partitioning,
        !       passively relocated and then added to the formation flux which might potentially not be empty due to previous
        !       iterations within the transition loop.
        CALL calc_litter_partitioning( &
            & nc, &                                         ! in
            & nsoil_w, &
            & num_sl_above_bedrock(:), &
            & dsl4jsb_Lctlib_param_tile(sla, current_tile), &
            & dsl4jsb_Lctlib_param_tile(growthform, current_tile), &
            & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
            & soil_depth_sl(:,:), &
            & root_fraction_sl(:,:), &
            & veg_litter_flux_harvest_mt(:,:,:), &          ! in
            & zero_mt(:,:), &
            & sb_formation_mt(:,:,:,:) )                    ! inout

        ! veg_litter_flux_harvest_mt is only temporarily required to carry the litter flux into the partitioning routine
        ! it will not be synchronised back and is set to zero here
        veg_litter_flux_harvest_mt(:,:,:) = 0.0_wp

        ! ---------------------------------------------------------------------------------------------------!
        ! Sync matrices back to variables
        dsl4jsb_Write_mt2L3D_to_vars_tile(SB_BGCM_FORMATION_ID, sb_formation_mt, current_tile)
        dsl4jsb_Write_mt2L2D_to_vars_tile(VEG_BGCM_POOL_ID, veg_pool_mt, current_tile)

        IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
          dsl4jsb_Write_mt1L2D_to_vars_tile(VEG_BGCM_HARVEST_TO_PROD_ID, veg_pp_flux_harvest_mt, current_tile)
        END IF
      END IF ! IF (ANY(lost_area(:, i_tile) > 0.0_wp))

      ! -----------------------------------------------------------------------------------------------------!
      ! The density on a tile only needs to be adapted on tiles which gain area
      ! -> adapt the change of individual density to account for the gained "empty" fraction (only for extending forest tiles)
      !    Note: while the to be conserved matter pools have been converted to grid area in the encompassing pplcc tasks,
      !          this is not the case for the density such that a scaling to the change in area is required
      IF (dsl4jsb_Lctlib_param_tile(ForestFlag, current_tile)) THEN
        DO ic = 1,nc
          IF (gained_area(ic, i_tile) > 0.0_wp) THEN
            scaling = initial_area(ic,i_tile) / (initial_area(ic,i_tile) + gained_area(ic, i_tile))
            delta_dens_ind(ic) = dens_ind(ic) * (scaling - 1._wp)
          END IF
        END DO
      END IF

      current_tile => current_tile%Get_next_sibling_tile()
    END DO

    !>
    !> 3. Collect to be transferred matter
    !>
    CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    !>
    !> 4. make the passive matter transfer
    !>
    CALL end_lcc_reloc(lcc_relocations, options, gained_area)

    !>
    !> 5. Finally - if product pools are used - add content of product pool flux to product pools on the gaining tiles
    !>
    ! Note: Alternatively to updating the product pool flux matrix from the associated bgcm variables within the current alcc task
    !       one could do the redistribution later after the post task of the pplcc process as it is done for the formation flux.
    !       This could be done via implementing an additional alcc task for Quincy or by generalising the updates_pools_on_harvest
    !       task to include this redistribution such that the slash and product pool use is more clustered.
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))
        i_tile = i_tile + 1

        IF (ANY(gained_area(:, i_tile) > 0.0_wp)) THEN
          bgcm_store => current_tile%bgcm_store

          dsl4jsb_Get_mt1L2D_tile(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt, current_tile)
          dsl4jsb_Get_mt1L2D_tile(VEG_BGCM_HARVEST_TO_PROD_ID, veg_pp_flux_harvest_mt, current_tile)

          dsl4jsb_Update_mt1L2D_tile(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt, current_tile)
          dsl4jsb_Update_mt1L2D_tile(VEG_BGCM_HARVEST_TO_PROD_ID, veg_pp_flux_harvest_mt, current_tile)

          veg_pp_fuel_mt(:,:) = veg_pp_fuel_mt(:,:) + veg_pp_flux_harvest_mt(:,:)
          veg_pp_flux_harvest_mt(:,:) = 0.0_wp

          dsl4jsb_Write_mt1L2D_to_vars_tile(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt, current_tile)
          dsl4jsb_Write_mt1L2D_to_vars_tile(VEG_BGCM_HARVEST_TO_PROD_ID, veg_pp_flux_harvest_mt, current_tile)
        END IF

        current_tile => current_tile%Get_next_sibling_tile()
      END DO
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE alcc_matter_translocation_iq
#endif

  ! ====================================================================================================== !
  !
  !> Redistribution of losses and gains of a forest pft to its age-classes if running with age classes.
  !
  SUBROUTINE redistribute_cover_fraction_changes_of_forest_pft(tile, options, cf_diff, i_first_child, gained_area, lost_area)
    USE mo_fage_interface,         ONLY: apply_fract_per_age_change_for_alcc_forest_gain, &
      &                                  distribute_forest_loss_from_alcc
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile          !< Forest pft for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options       !< Additional run-time parameters
    REAL(wp),                   INTENT(in)    :: cf_diff(:)    !< To be distributed changes of cover fractions
    INTEGER,                    INTENT(in)    :: i_first_child
       !< Index of first child tile in gained and lost area matrices
    REAL(wp),                   INTENT(inout) :: gained_area(:,:) !< Gained area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)   !< Lost area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':redistribute_cover_fraction_changes_of_forest_pft'

    REAL(wp), DIMENSION(options%nc) :: this_loss !< Area loss of forest PFT to be distributed to its children (age classes)
    REAL(wp), DIMENSION(options%nc) :: this_gain !< Area gains of forest PFT to be distributed to its children (age classes)
    ! -------------------------------------------------------------------------------------------------- !
    this_loss(:) = -1._wp * MIN(cf_diff(:), 0._wp)
    this_gain(:) = MAX(cf_diff(:), 0._wp)

    ! per definition, all gains are distributed to the first age class
    gained_area(:, i_first_child) = gained_area(:, i_first_child) + this_gain(:)
    CALL apply_fract_per_age_change_for_alcc_forest_gain(tile, options, this_gain(:))

    ! how losses are distributed depends on configurations of the FAGE process
    CALL distribute_forest_loss_from_alcc(tile, options, this_loss(:), i_first_child, lost_area(:, :))

  END SUBROUTINE redistribute_cover_fraction_changes_of_forest_pft


  ! ====================================================================================================== !
  !> #### Handle anthropogenic LCC if NLCC is active
  !
  !>  In case natural landcover change (nlcc) is active, the natural PFT fractions should not be
  !>  set to the map's fractions. The fractions of the dynamic PFT tiles need to be rescaled all
  !>  together, corresponding to the changes of the non-dynamic PFTs.
  !
  SUBROUTINE handle_dynamic_pfts(tile, options, initial_area, gained_area, lost_area)
    USE mo_nlcc_process,      ONLY: fract_small
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile              !< Parent tile for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< Additional run-time parameters
    REAL(wp),                   INTENT(inout) :: initial_area(:,:) !< Initial_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: gained_area(:,:)  !< Gained_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)    !< Lost_area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: ic, i_tile                                !< Looping indices
    TYPE(t_jsb_model), POINTER          :: model          !< Current instance of the model
    CLASS(t_jsb_tile_abstract), POINTER :: current_tile   !< Pointer on currently investigated tile7
    REAL(wp), DIMENSION(options%nc)     :: current_fract  !< Cover fraction if current tile
    REAL(wp), DIMENSION(options%nc)     :: area_dynamic   !< Cover fraction correction
    REAL(wp), DIMENSION(options%nc)     :: gained_dynamic !< Cover fraction correction
    REAL(wp), DIMENSION(options%nc)     :: lost_dynamic   !< Cover fraction correction
    CHARACTER(len=*), PARAMETER :: routine = modname//':handle_dynamic_pfts'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)

    ! Init area vectors
    lost_dynamic(:)   = 0.0_wp
    gained_dynamic(:) = 0.0_wp
    area_dynamic(:)   = 0.0_wp

    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1

      IF (dsl4jsb_Lctlib_param_tile(dynamic_PFT, current_tile)) THEN
        ! Sum up the total area fraction of the dynamic PFT tiles
        area_dynamic(:) = area_dynamic(:) + initial_area(:,i_tile)
      ELSE
        ! The area losses/gains of the non-dynamic PFTs define the joint gains/losses
        ! of the dynamic PFTs
        gained_dynamic(:) = gained_dynamic(:) + lost_area(:,i_tile)
        lost_dynamic(:)   = lost_dynamic(:)   + gained_area(:,i_tile)
      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1

      IF (dsl4jsb_Lctlib_param_tile(dynamic_PFT, current_tile)) THEN
        DO ic = 1, options%nc
          ! The dynamic PFTs are gaining or losing area proportional to their current area fractions
          IF (area_dynamic(ic) > fract_small) THEN
            gained_area(ic,i_tile) = MIN(gained_dynamic(ic) * initial_area(ic,i_tile) / area_dynamic(ic), &
              &                          1._wp - initial_area(ic,i_tile))
            lost_area(ic,i_tile)   = MIN(lost_dynamic(ic)   * initial_area(ic,i_tile) / area_dynamic(ic), &
              &                          initial_area(ic,i_tile))
          ELSE
            gained_area(ic,i_tile) = 0._wp
            lost_area(ic,i_tile)   = 0._wp
          END IF

          current_fract(ic) = initial_area(ic, i_tile) + gained_area(ic,i_tile) - lost_area(ic,i_tile)
        END DO
        ! Update the fractions of dynamic PFTs
        CALL current_tile%Set_fraction(options%ics, options%ice, options%iblk, fract=current_fract(:))

      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

  END SUBROUTINE handle_dynamic_pfts


  ! ====================================================================================================== !
  !> #### Check and if required correct gained and lost areas
  !>
  !> Verify that the total gained area fraction equals the lost area fractions and apply corrections in
  !> case of a mismatch.
  !
  SUBROUTINE check_and_correct_cf_mismatches(tile, options, nr_of_tiles, gained_area, lost_area)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile             !< Parent tile for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options          !< Additional run-time parameters
    INTEGER,                    INTENT(in)    :: nr_of_tiles      !< Number of tiles involved in the lcc calculation
    REAL(wp),                   INTENT(inout) :: gained_area(:,:) !< Gained_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)   !< Lost_area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL  :: gained_exceeds_loss !< Logical indicating if an area gain exceeds an area loss; used to assert carbon conservation
    INTEGER  :: ic, i_tile, i_sort  !< Looping indices
    REAL(wp) :: remaining_diff      !< Remaining difference between area losses and gains
    INTEGER  :: i_max_changing_tile !< Tile with the largest area change left
    CLASS(t_jsb_tile_abstract), POINTER            :: current_tile  !< Pointer on currently investigated tile
    REAL(wp), DIMENSION(options%nc)                :: current_fract !< Cover fraction if current tile
    REAL(wp), DIMENSION(options%nc)                :: cf_diff       !< Difference between gained and lost areas
    LOGICAL,  DIMENSION(nr_of_tiles)               :: masked_tiles  !< Mask of already considered tiles
    REAL(wp), DIMENSION(nr_of_tiles)               :: area_change   !< Area change to sort for (either losses or gains)
    REAL(wp), DIMENSION(options%nc, nr_of_tiles) :: cf_correction   !< Cover fraction correction
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_and_correct_cf_mismatches'
    ! -------------------------------------------------------------------------------------------------- !
    cf_diff(:) = SUM(gained_area(:,:),DIM=2) - SUM(lost_area(:,:),DIM=2)
    cf_correction(:,:) = 0.0_wp
    area_change(:) = 0.0_wp

    IF (ANY(ABS(cf_diff(:)) > EPSILON(1._wp))) THEN  !EPSILON(1._wp) ! 1.E-15_wp
      DO ic = 1, options%nc
        IF (ABS(cf_diff(ic)) > EPSILON(1._wp)) THEN
          masked_tiles(:) = .FALSE.
          remaining_diff = cf_diff(ic)
          DO i_sort = 1, nr_of_tiles
            IF (remaining_diff > EPSILON(1._wp)) THEN
              ! I.e. gain > loss
              area_change(:) = gained_area(ic,:)
              i_tile = MAXLOC(area_change(:), DIM=1, MASK=masked_tiles(:))
              cf_correction(ic,i_tile) = MIN(remaining_diff, gained_area(ic,i_tile))
              remaining_diff = MAX(0.0_wp, remaining_diff - gained_area(ic,i_tile))
              gained_area(ic,i_tile) = gained_area(ic,i_tile) - cf_correction(ic,i_tile)
              masked_tiles(i_tile) = .TRUE.
            ELSE IF (remaining_diff < -EPSILON(1._wp)) THEN
              ! I.e. loss > gain
              area_change(:) = lost_area(ic,:)
              i_tile = MAXLOC(area_change(:), DIM=1, MASK=masked_tiles(:))
              cf_correction(ic,i_tile) = - MIN(ABS(remaining_diff), lost_area(ic,i_tile))
              remaining_diff = - MAX(0.0_wp, ABS(remaining_diff) - lost_area(ic,i_tile))
              lost_area(ic,i_tile) = lost_area(ic,i_tile) - cf_correction(ic,i_tile)
              masked_tiles(i_tile) = .TRUE.
            ELSE
              EXIT
            END IF
          END DO
        END IF
      END DO ! ic = 1, nc

      ! If requried modify tiles and area changes accordingly
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))
        i_tile = i_tile + 1
        CALL current_tile%Get_fraction(options%ics, options%ice, options%iblk, fract=current_fract(:))
        current_fract(:) = current_fract(:) - cf_correction(:,i_tile)
        CALL current_tile%Set_fraction(options%ics, options%ice, options%iblk, fract=current_fract(:))
        current_tile => current_tile%Get_next_sibling_tile()
      END DO ! WHILE (ASSOCIATED(current_tile))
    END IF ! (ANY(ABS(cf_diff) > EPSILON(1._wp))) THEN

  END SUBROUTINE check_and_correct_cf_mismatches

#endif
END MODULE mo_alcc_interface
