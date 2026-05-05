!> QUINCY agriculture variables init
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
!>#### initialization of agriculture memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_agr_init
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,   ONLY: Q_AGR_

  dsl4jsb_Use_config(Q_AGR_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_agr_init
#ifndef __QUINCY_STANDALONE__
  PUBLIC read_fertiliser_input
#endif

  ! -------------------------------------------------------------------------------------------------- !
  !
  !> Type to hold variables read from input files
  TYPE t_q_agr_init_vars
    REAL(wp), POINTER ::  &
      & c3_crop_type(:,:), &    !< major c3 crop
      & c4_crop_type(:,:)       !< major c4 crop
  END TYPE t_q_agr_init_vars

  TYPE(t_q_agr_init_vars) :: q_agr_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_init'

CONTAINS

  ! ======================================================================================================= !
  !> Run agriculture init
  !>
  SUBROUTINE q_agr_init(tile)
    USE mo_q_agr_constants,       ONLY: qs_def_n_fertiliser, ix_wheat, ix_tecorn, ix_trcorn, ix_trsoy
    USE mo_q_assimi_constants,    ONLY: ic3phot
    dsl4jsb_Use_memory(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: box_tile
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid), POINTER           :: grid
    INTEGER                             :: iblk, ic, is         !< loop dimensions
    INTEGER                             :: nproma, nblks        !< dimensions
    INTEGER                             :: lctlib_ps_pathway
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_agr_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory_tile(Q_AGR_, box_tile)
    dsl4jsb_Real2D_onDomain :: n_fertiliser_c3
    dsl4jsb_Real2D_onDomain :: n_fertiliser_c4
    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onDomain :: crop_type_index
    dsl4jsb_Def_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of agr memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    grid   => Get_grid(model%grid_id)
    nblks  =  grid%nblks
    nproma =  grid%nproma
    CALL model%Get_top_tile(box_tile)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory_tile(Q_AGR_, box_tile)
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, n_fertiliser_c3, box_tile) ! out
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, n_fertiliser_c4, box_tile) ! out
    dsl4jsb_Get_memory(Q_AGR_)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, crop_type_index)                ! out
    dsl4jsb_Get_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !

#ifdef __QUINCY_STANDALONE__
    lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)

    !---- For QS: Determine crop type index from latitudes
    ! set defaults based on photosynthesis type and latitude
    IF (lctlib_ps_pathway == ic3phot) THEN ! C3 Crops
      DO ic = 1,nproma
        DO iblk = 1,nblks
          IF (grid%lat(ic,iblk) < 30.0_wp .AND. grid%lat(ic,iblk) > -30.0_wp ) THEN
            crop_type_index(ic,iblk) = ix_trsoy ! tropical C3 crop = soybean
          ELSE
            crop_type_index(ic,iblk) = ix_wheat ! temperate C3 crop = wheat
          END IF
        END DO
      END DO
    ELSE ! C4 crops
      DO ic = 1,nproma
        DO iblk = 1,nblks
          IF (grid%lat(ic,iblk) < 30.0_wp .AND. grid%lat(ic,iblk) > -30.0_wp ) THEN
            crop_type_index(ic,iblk) = ix_trcorn ! tropical C4 crop = corn
          ELSE
            crop_type_index(ic,iblk) = ix_tecorn ! temperate C4 crop = corn
          END IF
        END DO
      END DO
    END IF

    ! default fertiliser used in QS (overwritten in IQ with constant or transient data in read_fertiliser_input)
    n_fertiliser_c3(:,:) = qs_def_n_fertiliser
    n_fertiliser_c4(:,:) = qs_def_n_fertiliser
    ! ----------------------------------------------------------------------------------------------------- !
#else
    !---- For IQ: Read crop type from bc file (fertiliser is also read but in the time-loop not in the init)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      ! Read bc file on the box tile
      CALL q_agr_read_init_vars(tile)
    ELSE IF (tile%lcts(1)%lib_id /= 0) THEN
      IF (dsl4jsb_Lctlib_param(CropFlag)) THEN
        ! Determine crop_type_index for crops
        lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)

        IF (lctlib_ps_pathway == ic3phot) THEN
          ! C3 Crop
          crop_type_index(:,:) = q_agr_init_vars%c3_crop_type(:,:)
        ELSE
          ! C4 crop
          crop_type_index(:,:) = q_agr_init_vars%c4_crop_type(:,:)
        END IF
      END IF
    END IF

    IF (tile%Is_last_process_tile(Q_AGR_)) THEN
      CALL q_agr_finalise_init_vars()
    END IF

    !$ACC UPDATE DEVICE(crop_type_index) ASYNC(1)
    !$ACC WAIT(1)

#endif
  END SUBROUTINE q_agr_init

#ifndef __QUINCY_STANDALONE__
  ! ====================================================================================================== !
  !
  !> Reading of this years fertiliser input from input file
  !
  SUBROUTINE read_fertiliser_input(model_id, current_datetime)
    USE mo_io_units,                  ONLY: filename_max
    USE mo_jsb_time_iface,            ONLY: t_datetime
    USE mo_jsb_time,                  ONLY: is_time_experiment_start, get_year
    USE mo_jsb_io_netcdf,             ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_q_agr_constants,           ONLY: iq_background_n_fertiliser_rate
    USE mo_jsb_physical_constants,    ONLY: molar_mass_N

    dsl4jsb_Use_memory(Q_AGR_)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),            POINTER :: model
    TYPE(t_jsb_grid),             POINTER :: hgrid
    CLASS(t_jsb_tile_abstract),   POINTER :: tile

    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onDomain :: n_fertiliser_c3
    dsl4jsb_Real2D_onDomain :: n_fertiliser_c4

    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp),           POINTER :: ptr_2D(:,:)  ! tmp pointer
    INTEGER                     :: current_year
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=filename_max) :: filename_fertiliser_data
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_fertiliser_input'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_memory(Q_AGR_)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, n_fertiliser_c3)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, n_fertiliser_c4)

    IF (.NOT. dsl4jsb_Config(Q_AGR_)%l_read_fertiliser) THEN
      ! ... use background rate
      n_fertiliser_c3(:,:) = iq_background_n_fertiliser_rate
      n_fertiliser_c4(:,:) = iq_background_n_fertiliser_rate

    ELSE
      ! search a file ending on the current year
      current_year  = get_year(current_datetime)

      !>
      !> Assertion: routine currently expects filenames with 4 digits
      !>
      IF (( current_year > 9999) .OR. (current_year < 1000)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: this routine currently expects filenames with 4 digits.')
      END IF
      WRITE (filename_fertiliser_data,'(a,a,I4.4,a)') &
        & TRIM(dsl4jsb_Config(Q_AGR_)%fertiliser_filename_prefix), '_', current_year, ".nc"

      input_file = jsb_netcdf_open_input(TRIM(filename_fertiliser_data), model%grid_id)
      ptr_2D => input_file%Read_2d(variable_name='fertl_c3crops', fill_array = n_fertiliser_c3(:,:))
      ptr_2D => input_file%Read_2d(variable_name='fertl_c4crops', fill_array = n_fertiliser_c4(:,:))
      CALL input_file%Close()

      ! The actually applied fertiliser is the maximum of the default and the (converted) read fertiliser
      n_fertiliser_c3(:,:) = MAX(iq_background_n_fertiliser_rate, n_fertiliser_c3(:,:) / 10._wp / molar_mass_N)
      n_fertiliser_c4(:,:) = MAX(iq_background_n_fertiliser_rate, n_fertiliser_c4(:,:) / 10._wp / molar_mass_N)
    END IF

    !$ACC UPDATE DEVICE(n_fertiliser_c3, n_fertiliser_c4) ASYNC(1)
    !$ACC WAIT(1)

    IF (debug_on()) CALL message(TRIM(routine), 'Finishing routine')

  END SUBROUTINE read_fertiliser_input

  ! ====================================================================================================== !
  !>
  !> Finalise initialisation variables
  !>
  SUBROUTINE q_agr_finalise_init_vars

    !> Deallocation of variables that had only been needed during the initialisation phase for reading from input
    DEALLOCATE( &
      & q_agr_init_vars%c3_crop_type, &
      & q_agr_init_vars%c4_crop_type  &
      & )

  END SUBROUTINE q_agr_finalise_init_vars


  ! ====================================================================================================== !
  !
  !> Read variables from bc file to module type variable used to initialise fertiliser input
  !>
  SUBROUTINE q_agr_read_init_vars(tile)
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,              ONLY: missval
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    REAL(wp), POINTER  :: ptr_2D(:,:)
    TYPE(t_input_file) :: input_file
    INTEGER            :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':q_agr_read_init_vars'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    model => get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE( &
      & q_agr_init_vars%c3_crop_type(nproma, nblks), &
      & q_agr_init_vars%c4_crop_type(nproma, nblks)  &
      & )

    q_agr_init_vars%c3_crop_type(:,:) = missval
    q_agr_init_vars%c4_crop_type(:,:) = missval

    dsl4jsb_Get_config(Q_AGR_)

    IF (debug_on()) CALL message(TRIM(routine), 'Reading agriculture init vars from ' &
      &                          //TRIM(dsl4jsb_Config(Q_AGR_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(Q_AGR_)%bc_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading c3 crop type ...')
    ptr_2D => input_file%Read_2d(     &
      & variable_name='c3_crop_type', &
      & fill_array = q_agr_init_vars%c3_crop_type)

    IF (debug_on()) CALL message(TRIM(routine), 'reading c4 crop type ...')
    ptr_2D => input_file%Read_2d(     &
      & variable_name='c4_crop_type', &
      & fill_array = q_agr_init_vars%c4_crop_type)

  END SUBROUTINE q_agr_read_init_vars

#endif
#endif
END MODULE mo_q_agr_init
