!> Initialization of the the seb memory.
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
MODULE mo_seb_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid
  USE mo_jsb_grid,          ONLY: Get_grid

  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_memory(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: seb_init

  TYPE t_seb_init_vars
    REAL(wp), ALLOCATABLE :: skin_conductivity(:,:)
  END TYPE

  TYPE(t_seb_init_vars) :: seb_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_init'

CONTAINS

  !
  !> Intialize soil process (after memory has been set up)
  !
  SUBROUTINE seb_init(tile)
    USE mo_jsb_time,          ONLY: timestep_in_days
#ifdef __QUINCY_STANDALONE__
    USE mo_qs_forcing,          ONLY: forcing_options
    USE mo_jsb_math_constants,  ONLY: eps8
#else
    USE mo_pheno_parameters,  ONLY: pheno_param_jsbach
#endif
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init'

    dsl4jsb_Real2D_onDomain :: &
      & F_pseudo_soil_temp,    &
      & N_pseudo_soil_temp

    dsl4jsb_Def_memory(SEB_)

    model => get_model(tile%owner_model_id)

#ifdef __QUINCY_STANDALONE__
    CALL qs_seb_init_bc(tile)
#else
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL seb_read_init_vars(tile)
    END IF

    CALL seb_init_bc(tile)
    !CALL seb_init_ic(tile)

    IF (tile%Is_last_process_tile(SEB_)) THEN
      CALL seb_finalize_init_vars()
    END IF
#endif

    ! Initial values for the calculation of pseudo_soil_temp in SUBROUTINE calc_pseudo_soil_temp
    dsl4jsb_Get_memory(SEB_)


    IF (tile%contains_land .OR. model%config%use_tmx) THEN
      dsl4jsb_Get_var2D_onDomain(SEB_, N_pseudo_soil_temp) ! OUT
      dsl4jsb_Get_var2D_onDomain(SEB_, F_pseudo_soil_temp) ! OUT
#ifdef __QUINCY_STANDALONE__
      F_pseudo_soil_temp = EXP(- MAX(eps8, timestep_in_days(forcing_options%forcing_hour)) / 10._wp)
#else
      F_pseudo_soil_temp = EXP(- timestep_in_days(model_id=tile%owner_model_id)/pheno_param_jsbach%EG_SG%tau_pseudo_soil)
#endif
      N_pseudo_soil_temp = 1._wp / (1._wp - F_pseudo_soil_temp)

      !$ACC UPDATE DEVICE(F_pseudo_soil_temp, N_pseudo_soil_temp) ASYNC(1)
    END IF
  END SUBROUTINE seb_init

#ifdef __QUINCY_STANDALONE__
#else
  SUBROUTINE seb_read_init_vars(tile)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    TYPE(t_input_file) :: input_file
    INTEGER :: nproma, nblks

    REAL(wp), POINTER :: ptr(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_read_init_vars'

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
      ALLOCATE(seb_init_vars%skin_conductivity(nproma, nblks))

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SEB_)%bc_filename), model%grid_id)

      IF (input_file%Has_var('skin_conductivity')) THEN
        ptr => input_file%Read_2d('skin_conductivity', fill_array=seb_init_vars%skin_conductivity)
      ELSE
        CALL finish(TRIM(routine), '*** Error: BC file does not contain skin_conductivity variable.')
      END IF

      CALL input_file%Close()
    END IF

  END SUBROUTINE seb_read_init_vars
#endif

  SUBROUTINE seb_finalize_init_vars
    IF (ALLOCATED(seb_init_vars%skin_conductivity)) DEALLOCATE(seb_init_vars%skin_conductivity)
  END SUBROUTINE seb_finalize_init_vars

  ! ======================================================================================================= !
  !>
  !> SEB_ init bc - quincy standalone
  !>
  SUBROUTINE qs_seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':qs_seb_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
      IF (tile%is_land .OR. model%config%use_tmx) THEN
        dsl4jsb_var2D_onDomain(SEB_, skin_conductivity) = 40.0_wp   ! estimated from icon-land input file
      END IF
    END IF

  END SUBROUTINE qs_seb_init_bc

  ! ======================================================================================================= !
  !>
  !> SEB_ init bc
  !>
  SUBROUTINE seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
      IF (tile%is_land .OR. model%config%use_tmx) THEN
        dsl4jsb_var2D_onDomain(SEB_, skin_conductivity) = seb_init_vars%skin_conductivity
        !$ACC UPDATE ASYNC(1) DEVICE(dsl4jsb_var2D_onDomain(SEB_, skin_conductivity))
      END IF
    END IF

  END SUBROUTINE seb_init_bc

  SUBROUTINE seb_init_ic(tile)

    !USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! TYPE(t_jsb_model),       POINTER :: model

    ! dsl4jsb_Def_config(SEB_)
    !dsl4jsb_Def_memory(SEB_)

    !REAL(wp), POINTER :: &
    !  & ptr_2D(:,  :)       !< temporary pointer

    !TYPE(t_input_file) :: input_file

    !INTEGER  :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_ic'

    !model => get_model(tile%owner_model_id)

    ! Get seb config
    ! dsl4jsb_Get_config(SEB_)

    ! IF (.NOT. model%Is_process_enabled(SEB_)) RETURN

    ! Get seb memory of the tile
    !dsl4jsb_Get_memory(SEB_)

    ! Initialize parameters
    !
    ! IF (debug_on()) CALL message(TRIM(routine), 'Setting initial state of seb memory for tile '// &
      ! &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(SEB_)%ic_filename))

    ! Initial surface temperature
    ! Is initialized in sse_init_ic together with soil temperatures

  END SUBROUTINE seb_init_ic

#endif
END MODULE mo_seb_init
