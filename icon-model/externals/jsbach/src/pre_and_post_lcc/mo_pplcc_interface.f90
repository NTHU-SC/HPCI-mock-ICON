!> interface to the pplcc process
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
!>#### Contains the interfaces to the pplcc process which does the pre- and postprocessing for all
!>     land cover change processes
!>
!> Note:
!>
!>       This module takes advantage of the handling of the integrate and aggregate routines, which both
!>       call the same routines in this module (integrate = aggregate).
!>       Thereby, each tile is processed by the same routine for the tasks of this process.
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_pplcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish
  USE mo_util,            ONLY: int2string

  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  !USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_var_class,       ONLY: REAL2D, REAL3D

  USE mo_jsb_time,          ONLY: is_newday
  USE mo_jsb_cqt_class,     ONLY: Get_cqt_name, LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                             BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE, &
    &                             IQ_2L2D_POOL_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE, &
    &                             IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE

#ifndef __NO_QUINCY__
  USE mo_quincy_model_config, ONLY: QCANOPY
  USE mo_quincy_output_class, ONLY: unit_sb_pool, unit_sb_flux, unit_veg_pool_flux
#endif

  USE mo_carbon_interface,  ONLY: recalc_carbon_per_tile_vars

  ! Use of processes in this module
  dsl4jsb_Use_processes PPLCC_, PHENO_, CARBON_

  ! Use of process configurations
  !  dsl4jsb_Use_config(PPLCC_)

  ! Use of process memories
  dsl4jsb_Use_memory(PPLCC_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pplcc_tasks, global_pplcc_diagnostics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_interface'

  !> Type definition for pplcc pre task (preprocessing for lcc)
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_pplcc_pre_lcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_pplcc_pre_lcc  !< prepares lcc processes for leaves
    PROCEDURE, NOPASS :: Aggregate => update_pplcc_pre_lcc  !< prepares lcc processes for all other nodes
  END TYPE tsk_pplcc_pre_lcc

  !> Constructor interface for pplcc pre task (preprocessing for lcc)
  INTERFACE tsk_pplcc_pre_lcc
    PROCEDURE Create_task_pplcc_pre_lcc                      !< Constructor function for task
  END INTERFACE tsk_pplcc_pre_lcc

  !> Type definition for pplcc post task (postprocessing for lcc)
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_pplcc_post_lcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_pplcc_post_lcc  !< post-processing for lcc processes on leaves
    PROCEDURE, NOPASS :: Aggregate => update_pplcc_post_lcc  !< post-processing for lcc processes on other nodes
  END TYPE tsk_pplcc_post_lcc

  !> Constructor interface for pplcc post task (postprocessing for lcc)
  INTERFACE tsk_pplcc_post_lcc
    PROCEDURE Create_task_pplcc_post_lcc                        !< Constructor function for task
  END INTERFACE tsk_pplcc_post_lcc

CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for pplcc pre lcc task (preprocessing for lcc)
  !
  FUNCTION Create_task_pplcc_pre_lcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "pplcc_pre_lcc"
    ! -------------------------------------------------------------------------------------------------- !

    ALLOCATE(tsk_pplcc_pre_lcc::return_ptr)
    CALL return_ptr%Construct(name='pplcc_pre_lcc', process_id=PPLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_pplcc_pre_lcc

  ! ====================================================================================================== !
  !
  !> Constructor for pplcc post lcc task (postprocessing for lcc)
  !
  FUNCTION Create_task_pplcc_post_lcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "pplcc_post_lcc"
    ! -------------------------------------------------------------------------------------------------- !

    ALLOCATE(tsk_pplcc_post_lcc::return_ptr)
    CALL return_ptr%Construct(name='pplcc_post_lcc', process_id=PPLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_pplcc_post_lcc

  ! ====================================================================================================== !
  !
  !> Register tasks for pplcc process
  !
  SUBROUTINE Register_pplcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this      !< Instance of pplcc process class
    INTEGER,              INTENT(in)    :: model_id  !< Model id
    ! -------------------------------------------------------------------------------------------------- !

    CALL this%Register_task(tsk_pplcc_pre_lcc(model_id))
    CALL this%Register_task(tsk_pplcc_post_lcc(model_id))

  END SUBROUTINE Register_pplcc_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "pplcc_pre_lcc"
  !>
  !> Task "pplcc_pre_lcc" prepares lcc processes
  !> -- for now: convert all cq vars to per grid cell.
  !
  SUBROUTINE update_pplcc_pre_lcc(tile, options)

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model

    REAL(wp), DIMENSION(options%nc) :: conversion_factor

    INTEGER  :: nc   !< Number of cells in current block
    INTEGER  :: iblk, ics, ice, ic, cqt_id, i_cqt, i_cq_var
    REAL(wp) :: dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_pplcc_pre_lcc'

    CHARACTER(len=:), ALLOCATABLE :: this_unit, expected_string_in_unit

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(PHENO_)
!    dsl4jsb_Def_memory(PPLCC_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(PPLCC_)) RETURN

    ! If not newday, do nothing
    IF (.NOT. is_newday(options%current_datetime,dtime)) RETURN

    ! If no lcc process is running this task does not need to be executed
    IF (.NOT. model%Do_fractions_change()) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (tile%Has_conserved_quantities()) THEN

#ifndef __NO_QUINCY__
      IF (model%config%model_scheme == MODEL_QUINCY) THEN

        IF (model%config%qmodel_id == QCANOPY) THEN
          ! If running in CANOPY configuration nothing needs to be conserved because the states are prescribed
          RETURN
        END IF

        ! For the quincy model the variables need to be up to date with the content of the matrices before they are converted
        IF (ASSOCIATED(tile%bgcm_store)) THEN
          CALL tile%bgcm_store%Write_matrices_to_bgcm_vars(tile%name, ics, ice, iblk)
        END IF
      END IF
#endif

      ! Iterate over all cq types and their vars and convert
      DO i_cqt = 1,tile%nr_of_cqts

        ! conversion of m2 tile to m2 box is required for several types of conserved quantities
        CALL tile%Get_fraction(ics, ice, iblk, fract=conversion_factor(:))

        ! conversion depends on the cq type
        cqt_id = tile%conserved_quantities(i_cqt)%p%type_id
        SELECT CASE (cqt_id)
          CASE (LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE)

            ! in the current implementation of jsb4 carbon variables are only available on the pft tiles
            IF (.NOT. index(tile%name, 'pft') /= 0) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types only expected for pft tiles, ' &
                & //' but found on '// trim(tile%name))
            END IF

            ! in jsb4 the carbon variables are specified per canopy
            expected_string_in_unit = 'm-2(canopy)'

            ! therefore, PHENO variables are required to determine the conversion factor
            ! -> Assert that PHENO is active for this tile
            IF (.NOT. tile%Has_process_memory(PHENO_)) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types require pheno vars for pplcc conversion, ' &
                & //' but pheno memory is not available for '// trim(tile%name))
            END IF

            !- Set pointers to memory on the current tile
            dsl4jsb_Get_memory(PHENO_)
            ! and get required variables on this tile
            dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction)    ! in
            dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)           ! in

            conversion_factor(:) = conversion_factor(:) * veg_fract_correction(:) * fract_fpc_max(:)

#ifndef __NO_QUINCY__
          CASE (IQ_2L2D_POOL_CQ_TYPE, IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE)
            ! in the current implementation of quincy to be conserved variables are only available on the pft tiles
            IF (.NOT. index(tile%name, 'pft') /= 0) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: for quincy cq types are only expected for pft tiles, ' &
                & //' but found on '// trim(tile%name))
            END IF

            IF ((cqt_id == IQ_2L2D_POOL_CQ_TYPE) .OR. (cqt_id == IQ_1L2D_POOL_CQ_TYPE)) THEN
              ! in quincy the 2D pool type variables are in mol per tile
              expected_string_in_unit = 'mol m-2'
            ELSE IF(cqt_id == IQ_FLUX_CQ_TYPE) THEN
              expected_string_in_unit = unit_veg_pool_flux
            ELSE IF(cqt_id == IQ_SL_POOL_CQ_TYPE) THEN
              ! in quincy the 3D pool type variables are in mol per tile and layer
              expected_string_in_unit = unit_sb_pool
            ELSE IF(cqt_id == IQ_SL_FLUX_CQ_TYPE) THEN
              ! in quincy the 3D flux type variables are in mol per tile and layer per timestep-1
              expected_string_in_unit = unit_sb_flux
            END IF

            ! Tile fractions in quincy are currently only modified by the alcc process
            ! thus only conserved quantities on cells where a tile fraction is modified by alcc need to be converted
            DO ic = 1,nc
              IF (.NOT. tile%l_fract_alcc_change(ics+ic-1, iblk)) THEN
                conversion_factor(ic) = 1.0_wp
              END IF
            END DO
#endif

          ! CASE (WATER_CQ_TYPE)
          CASE DEFAULT
            !>
            !> Assertion: for each CQT a conversion needs to be specified in pplcc_pre_lcc
            !>
            CALL finish(TRIM(routine), 'Conversion factors unspecified for '//Get_cqt_name(cqt_id))
        END SELECT

        ! Convert all variables of this type
        DO i_cq_var = 1,tile%conserved_quantities(i_cqt)%p%no_of_vars
          ! Assert expected unit
          this_unit = tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%unit
          IF (.NOT. INDEX(TRIM(this_unit), expected_string_in_unit) > 0) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: variable does not have the expected unit, ' &
              & // 'unit is expected to contain ' // expected_string_in_unit // ' but '// trim(this_unit) // ' does not.')
          END IF

          IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), '..... and var '// &
              & TRIM(tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%full_name))

          !JN-TODO: make a function in a new process module? -> DSL?!
          ! Note: for bookkeeping I would have wanted to change the unit, but seems not to be intended to be changed
          SELECT CASE (tile%conserved_quantities(i_cqt)%p%var_type)
          CASE (REAL2D)
            tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr2d(ics:ice,iblk) &
              & = tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr2d(ics:ice,iblk) &
              &   * conversion_factor(:)
          CASE (REAL3D)
            DO ic = ics,ice
              tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr3d(ic,:,iblk) &
                & = tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr3d(ic,:,iblk) &
                &   * conversion_factor(ic-ics+1)
            END DO
          CASE DEFAULT
            CALL finish(TRIM(routine), 'Unexpected variable type for conserved quantity variable, ' &
                & // 'expected 2 or 3 dimensions but found ' // int2string(tile%conserved_quantities(i_cqt)%p%var_type)// &
                & ' for cqt '//TRIM(Get_cqt_name(tile%conserved_quantities(i_cqt)%p%type_id)) // ', please check!')
          END SELECT
        END DO
      END DO

#ifndef __NO_QUINCY__
      IF (model%config%model_scheme == MODEL_QUINCY) THEN
        ! The bgcm matrices need to be updated after conversion because they are required in the alcc interface
        IF (ASSOCIATED(tile%bgcm_store)) THEN
          CALL tile%bgcm_store%Store_bgcm_vars_in_matrices(tile%name, ics, ice, iblk)
        END IF
      END IF
#endif

    END IF ! IF (tile%Has_conserved_quantities())

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_pplcc_pre_lcc

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "pplcc_post_lcc"
  !>
  !> Task "pplcc_post_lcc" does a postprocessing after all lcc processes were executed
  !> -- for now: convert all variables of conserved quantities (cq) back to their initial reference area.
  !> -- and if the carbon process is active: recalculate the tile average (ta) variables.
  !
  SUBROUTINE update_pplcc_post_lcc(tile, options)

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model

    REAL(wp), DIMENSION(options%nc) :: conversion_factor

    INTEGER  :: iblk, ics, ice, ic, i_cqt, i_cq_var, cqt_id, nc
    REAL(wp) :: dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_pplcc_post_lcc'

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    dtime   = options%dtime
    nc      = options%nc

    model => Get_model(tile%owner_model_id)

    ! If not newday, do nothing
    IF (.NOT. is_newday(options%current_datetime,dtime)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! PPLCC fraction diagnostics are calculated on the box tile, although the PPLCC process is calculated
    ! on leafs (compare usecase definition). We only enter this routine on the box tile, because
    ! update_pplcc_post_lcc is not only defined as Integrate but also as Aggregate procedure and thus is
    ! called from the box tile on the way up the tile tree.
    IF (.NOT. ASSOCIATED(tile%parent)) CALL pplcc_fraction_diagnostics(tile, options)

    ! If process is not calculated on this tile, conversions are not necessary
    IF (.NOT. tile%Is_process_calculated(PPLCC_)) RETURN

    ! Also: if no lcc process is running this task does not further needs to be executed
    IF (.NOT. model%Do_fractions_change()) RETURN

    IF (tile%Has_conserved_quantities()) THEN

#ifndef __NO_QUINCY__
      IF (model%config%model_scheme == MODEL_QUINCY) THEN
        IF (model%config%qmodel_id == QCANOPY) THEN
          ! If running in CANOPY configuration nothing needs to be conserved because the states are prescribed
          RETURN
        END IF
        ! Since the lcc framework mainly working on the variables and only on additionally synchronised
        ! chuncks of selected matrices no sync from matrix to variables should be conducted before post lcc conversion
      END IF
#endif

      ! Iterate over all cq types and their vars and convert
      DO i_cqt = 1,tile%nr_of_cqts

        ! conversion of m2 tile to m2 box is required for several types of conserved quantities
        CALL tile%Get_fraction(ics, ice, iblk, fract=conversion_factor(:))

        ! conversion depends on the cq type
        cqt_id = tile%conserved_quantities(i_cqt)%p%type_id

        SELECT CASE (cqt_id)
          CASE (LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE)

            ! in the current implementation of jsb4 carbon variables are only available on the pft tiles
            IF (.NOT. index(tile%name, 'pft') /= 0) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types only expected for pft tiles, ' &
                  & //' but found on '// trim(tile%name))
            END IF

            ! in jsb4 the carbon variables are specified as mol per canopy, therefore adapt the conversion factor
            ! -> Assert that PHENO memory is available on this tile
            IF (.NOT. tile%Has_process_memory(PHENO_)) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types require pheno vars for conversion, ' &
                  & //' but pheno memory is not available for '// trim(tile%name))
            END IF

            !- Set pointers to memory on the current tile
            dsl4jsb_Get_memory(PHENO_)
            ! and get required variables on this tile
            dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction)    ! in
            dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)           ! in
            conversion_factor(:) = conversion_factor(:) * veg_fract_correction(:) * fract_fpc_max(:)

          CASE (IQ_2L2D_POOL_CQ_TYPE, IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE)
            ! in the current implementation of quincy to be conserved variables are only available on the pft tiles
            IF (.NOT. index(tile%name, 'pft') /= 0) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: for quincy cq types are only expected for pft tiles, ' &
                  & //' but found on '// trim(tile%name))
            END IF
            ! Tile fractions in quincy are currently only modified by the alcc process
            ! thus only conserved quantities on cells where a tile fraction is modified by alcc need to be converted
            DO ic = 1,nc
              IF (.NOT. tile%l_fract_alcc_change(ics+ic-1, iblk)) THEN
                conversion_factor(ic) = 1.0_wp
              END IF
            END DO

          ! CASE (WATER_CQ_TYPE)
          CASE DEFAULT
            !>
            !> Assertion: for each CQT a conversion needs to be specified in pplcc_post_lcc
            !>
            CALL finish(TRIM(routine), 'Conversion factors unspecified for '//Get_cqt_name(cqt_id))
        END SELECT

        ! Convert all variables of this type
        DO i_cq_var = 1,tile%conserved_quantities(i_cqt)%p%no_of_vars

          IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), '..... and var '// &
              & TRIM(tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%full_name))

          !JN-TODO: make a function in new process module?
          SELECT CASE (tile%conserved_quantities(i_cqt)%p%var_type)
          CASE (REAL2D)
            ! only do the conversion where the conversion factor is not zero
            WHERE (conversion_factor(:) > 0.0_wp)
              tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr2d(ics:ice,iblk) &
                & = tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr2d(ics:ice,iblk) &
                &   / conversion_factor(:)
            ELSEWHERE
              tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr2d(ics:ice,iblk) = 0.0_wp
            END WHERE
          CASE (REAL3D)
            DO ic = ics,ice
              IF (conversion_factor(ic-ics+1) > 0.0_wp) THEN
                ! only do the conversion where the conversion factor is not zero
                tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr3d(ic,:,iblk) &
                  & = tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr3d(ic,:,iblk) &
                  &   / conversion_factor(ic-ics+1)
              ELSE
                tile%conserved_quantities(i_cqt)%p%cq_vars(i_cq_var)%p%ptr3d(ic,:,iblk) = 0.0_wp
              END IF
            END DO
          CASE DEFAULT
            CALL finish(TRIM(routine), 'Unexpected variable type for conserved quantity variable, ' &
                & // 'expected 2 or 3 dimensions but found ' // int2string(tile%conserved_quantities(i_cqt)%p%var_type)// &
                & ' for cqt '//TRIM(Get_cqt_name(tile%conserved_quantities(i_cqt)%p%type_id)) // ', please check!')
          END SELECT
        END DO
      END DO

      IF (tile%Is_process_calculated(CARBON_)) THEN
        CALL recalc_carbon_per_tile_vars(tile, options)
      END IF

#ifndef __NO_QUINCY__
      IF (model%config%model_scheme == MODEL_QUINCY) THEN
        ! For the quincy model the matrices need to be updated after post lcc conversion
        IF (ASSOCIATED(tile%bgcm_store)) THEN
          CALL tile%bgcm_store%Store_bgcm_vars_in_matrices(tile%name, ics, ice, iblk)
        END IF
      END IF
#endif

    END IF ! IF (tile%Has_conserved_quantities())

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_pplcc_post_lcc

  !-----------------------------------------------------------------------------------------------------
  !> Post LCC diagnostics of land cover fractions, i.e. aggregation of different land cover types
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE pplcc_fraction_diagnostics(tile, options)

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options  !< Additional run-time parameters

    ! Local variables
    TYPE(t_jsb_model), POINTER           :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: ptr_tile, ptr_pft
    REAL(wp), DIMENSION(options%nc)      :: pft_fract, veg_fract, bare_fract
    INTEGER                              :: iblk, ics, ice
    INTEGER                              :: ilct, nlct
    INTEGER                              :: ipft, npft
    REAL(wp), DIMENSION(options%nc)      :: correction_factor

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(PPLCC_)
    dsl4jsb_Def_memory(PHENO_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':pplcc_fraction_diagnostics'

    ! Pointers to variables in memory
    !
    dsl4jsb_Real2D_onChunk :: fract_fpc_max

    dsl4jsb_Real2D_onChunk :: tree_fract
    dsl4jsb_Real2D_onChunk :: shrub_fract
    dsl4jsb_Real2D_onChunk :: grass_fract
    dsl4jsb_Real2D_onChunk :: crop_fract
    dsl4jsb_Real2D_onChunk :: pasture_fract
    dsl4jsb_Real2D_onChunk :: baresoil_fract
    dsl4jsb_Real2D_onChunk :: C3pft_fract
    dsl4jsb_Real2D_onChunk :: C4pft_fract


    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_memory(PPLCC_)
    model => Get_model(tile%owner_model_id)

    ! The fraction variables are defined as fractions relative to the grid cell fraction. Thus
    ! only calculations on the box tile make sense.
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    IF (model%config%model_scheme == MODEL_JSBACH) THEN
      dsl4jsb_Get_memory(PHENO_)
      ! Get required variables on this tile
      dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)        ! in
      correction_factor(:) = fract_fpc_max(:)
    ELSE
      correction_factor(:) = 1.0_wp ! Note: In quincy all state variables are defined for 100% of the cover fraction
    END IF

    dsl4jsb_Get_var2D_onChunk(PPLCC_, tree_fract)           ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, shrub_fract)          ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, grass_fract)          ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, crop_fract)           ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, pasture_fract)        ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, baresoil_fract)       ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, C3pft_fract)          ! out
    dsl4jsb_Get_var2D_onChunk(PPLCC_, C4pft_fract)          ! out

    ! Initialization
    tree_fract     = 0._wp
    shrub_fract    = 0._wp
    grass_fract    = 0._wp
    crop_fract     = 0._wp
    pasture_fract  = 0._wp
    baresoil_fract = 0._wp
    C3pft_fract    = 0._wp
    C4pft_fract    = 0._wp


    ! To get the PFT fractions we first need a pointer to the land tile,
    ! and from the land tile to the veg tile, which has PFT children

    ! Get pointer to the land tile
    nlct = tile%Get_no_of_children()
    DO ilct=1,nlct
      IF (ilct == 1) THEN
        ptr_tile => tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! Only the land tile needs to be considered
      IF (TRIM(ptr_tile%name) == 'land') EXIT
    END DO

    ! Get pointer to the veg tile
    nlct = ptr_tile%Get_no_of_children()
    DO ilct=1,nlct
      IF (ilct == 1) THEN
        ptr_tile => ptr_tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! Only the veg tile needs to be considered
      IF (TRIM(ptr_tile%name) == 'veg') EXIT
    END DO

    npft = ptr_tile%Get_no_of_children()
    DO ipft=1,npft

      ! Point to the correct pft
      IF (ipft == 1) THEN
        ptr_pft => ptr_tile%Get_first_child_tile()
      ELSE
        ptr_pft => ptr_pft%Get_next_sibling_tile()
      END IF

      ! Get fraction for this PFT
      CALL ptr_pft%Get_fraction(ics, ice, iblk, fract=pft_fract(:))

      IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%NaturalVegFlag) THEN
        IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%ForestFlag) THEN
          tree_fract(:) = tree_fract(:) + pft_fract(:) * correction_factor(:)
        ELSE IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%GrassFlag) THEN
          grass_fract(:) = grass_fract(:) + pft_fract(:) * correction_factor(:)
        ELSE   ! Current lctlib file does not include a shrub flag
          shrub_fract(:) = shrub_fract(:) + pft_fract(:) * correction_factor(:)
        END IF
      ELSE IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%CropFlag) THEN
        crop_fract(:) = crop_fract(:) + pft_fract(:) * correction_factor(:)
      ELSE IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%PastureFlag) THEN
        pasture_fract(:) = pasture_fract(:) + pft_fract(:) * correction_factor(:)
      END IF
      IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%C4flag) THEN
        C4pft_fract(:) = C4pft_fract(:) + pft_fract(:) * correction_factor(:)
      ELSE
        C3pft_fract(:) = C3pft_fract(:) + pft_fract(:) * correction_factor(:)
      END IF

    END DO ! ipft

    ! Fraction of the veg tile
    CALL ptr_tile%Get_fraction(ics, ice, iblk, fract=veg_fract(:))

    ! Bare land tile
    ! Note: In the current usecase the land tile does not have a bare soil child, thus
    !       bare_fract is set to zero.
    bare_fract(:)=0._wp

    ! Baresoil fraction
    baresoil_fract(:) = bare_fract(:) + veg_fract(:) * (1._wp-correction_factor(:))


  END SUBROUTINE pplcc_fraction_diagnostics

  !-----------------------------------------------------------------------------------------------------
  !> Global land mean pplcc output, i.e. global mean area of different land cover types
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_pplcc_diagnostics(tile)

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_memory(PPLCC_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_pplcc_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: tree_fract
    dsl4jsb_Real2D_onDomain :: shrub_fract
    dsl4jsb_Real2D_onDomain :: grass_fract
    dsl4jsb_Real2D_onDomain :: crop_fract
    dsl4jsb_Real2D_onDomain :: pasture_fract
    dsl4jsb_Real2D_onDomain :: baresoil_fract
    dsl4jsb_Real2D_onDomain :: C3pft_fract
    dsl4jsb_Real2D_onDomain :: C4pft_fract

    REAL(wp), POINTER       :: tree_area_gsum(:)
    REAL(wp), POINTER       :: shrub_area_gsum(:)
    REAL(wp), POINTER       :: grass_area_gsum(:)
    REAL(wp), POINTER       :: crop_area_gsum(:)
    REAL(wp), POINTER       :: pasture_area_gsum(:)
    REAL(wp), POINTER       :: baresoil_area_gsum(:)
    REAL(wp), POINTER       :: C3pft_area_gsum(:)
    REAL(wp), POINTER       :: C4pft_area_gsum(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)


    dsl4jsb_Get_memory(PPLCC_)
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  tree_fract)                 ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  shrub_fract)                ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  grass_fract)                ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  crop_fract)                 ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  pasture_fract)              ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  baresoil_fract)             ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  C3pft_fract)                ! in
    dsl4jsb_Get_var2D_onDomain(PPLCC_,  C4pft_fract)                ! in


    tree_area_gsum     => PPLCC__mem%tree_area_gsum%ptr(:)        ! out
    shrub_area_gsum    => PPLCC__mem%shrub_area_gsum%ptr(:)       ! out
    grass_area_gsum    => PPLCC__mem%grass_area_gsum%ptr(:)       ! out
    crop_area_gsum     => PPLCC__mem%crop_area_gsum%ptr(:)        ! out
    pasture_area_gsum  => PPLCC__mem%pasture_area_gsum%ptr(:)     ! out
    baresoil_area_gsum => PPLCC__mem%baresoil_area_gsum%ptr(:)    ! out
    C3pft_area_gsum    => PPLCC__mem%C3pft_area_gsum%ptr(:)       ! out
    C4pft_area_gsum    => PPLCC__mem%C4pft_area_gsum%ptr(:)       ! out


    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Domain Mask - to mask all halo cells for global sums (otherwise these
    ! cells are counted twice)
    ALLOCATE (in_domain(grid%nproma,grid%nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ALLOCATE (scaling(grid%nproma,grid%nblks))

    ! Calculate 1d global land variables, if requested for output
    ! Unit transformation from [m2] to [Mio km2]: 1.e-12
    scaling(:,:) = area(:,:) * in_domain(:,:) * 1.e-12_wp
    IF (PPLCC__mem%tree_area_gsum%is_in_output)     &
      &  tree_area_gsum       = global_sum_array(tree_fract(:,:)     * scaling(:,:))
    IF (PPLCC__mem%shrub_area_gsum%is_in_output)    &
      &  shrub_area_gsum      = global_sum_array(shrub_fract(:,:)    * scaling(:,:))
    IF (PPLCC__mem%grass_area_gsum%is_in_output)    &
      &  grass_area_gsum      = global_sum_array(grass_fract(:,:)    * scaling(:,:))
    IF (PPLCC__mem%crop_area_gsum%is_in_output)     &
      &  crop_area_gsum       = global_sum_array(crop_fract(:,:)     * scaling(:,:))
    IF (PPLCC__mem%pasture_area_gsum%is_in_output)  &
      &  pasture_area_gsum    = global_sum_array(pasture_fract(:,:)  * scaling(:,:))
    IF (PPLCC__mem%baresoil_area_gsum%is_in_output) &
      &  baresoil_area_gsum   = global_sum_array(baresoil_fract(:,:) * scaling(:,:))
    IF (PPLCC__mem%C3pft_area_gsum%is_in_output)    &
      &  C3pft_area_gsum      = global_sum_array(C3pft_fract(:,:)    * scaling(:,:))
    IF (PPLCC__mem%C4pft_area_gsum%is_in_output)    &
      &  C4pft_area_gsum      = global_sum_array(C4pft_fract(:,:)    * scaling(:,:))

    DEALLOCATE (scaling, in_domain)

  END SUBROUTINE global_pplcc_diagnostics

#endif
END MODULE mo_pplcc_interface
