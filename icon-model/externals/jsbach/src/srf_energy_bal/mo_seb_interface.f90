!> Contains the interfaces for the surface energy balance process.
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_seb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_jsb_control,     ONLY: debug_on, acc_stream
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SEB_, HYDRO_, SSE_, A2L_, RAD_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)

  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(RAD_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_seb_tasks, seb_check_temperature_range
#ifndef __QUINCY_STANDALONE__
  PUBLIC :: global_seb_diagnostics
#endif
#ifdef __QUINCY_STANDALONE__
  ! USEd by the mo_qs_model_interface
  PUBLIC :: update_surface_energy, update_asselin, update_surface_fluxes, update_energy_balance
#endif

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_energy
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_energy
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_energy
  END TYPE tsk_surface_energy

  INTERFACE tsk_surface_energy
    PROCEDURE Create_task_surface_energy
  END INTERFACE tsk_surface_energy

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_asselin_filter
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_asselin
    PROCEDURE, NOPASS :: Aggregate => aggregate_asselin
  END TYPE tsk_asselin_filter

  INTERFACE tsk_asselin_filter
    PROCEDURE Create_task_asselin_filter
  END INTERFACE tsk_asselin_filter

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_fluxes
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_fluxes
  END TYPE tsk_surface_fluxes

  INTERFACE tsk_surface_fluxes
    PROCEDURE Create_task_surface_fluxes
  END INTERFACE tsk_surface_fluxes

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_energy_balance
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_energy_balance
    PROCEDURE, NOPASS :: Aggregate => aggregate_energy_balance
  END TYPE tsk_energy_balance

  INTERFACE tsk_energy_balance
    PROCEDURE Create_task_energy_balance
  END INTERFACE tsk_energy_balance

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_energy"
  !!
  FUNCTION Create_task_surface_energy(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_energy::return_ptr)
    CALL return_ptr%Construct(name='surface_energy', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for asselin_filter task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "asselin_filter"
  !!
  FUNCTION Create_task_asselin_filter(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_asselin_filter::return_ptr)
    CALL return_ptr%Construct(name='asselin_filter', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_asselin_filter

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy fluxes
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_fluxes"
  !!
  FUNCTION Create_task_surface_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_fluxes::return_ptr)
    CALL return_ptr%Construct(name='surface_fluxes', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for energy balance residual computation
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "energy_balance"
  !!
  FUNCTION Create_task_energy_balance(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_energy_balance::return_ptr)
    CALL return_ptr%Construct(name='energy_balance', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_energy_balance

  ! =======================================================================================================
  !> Register tasks for surface energy process
  !!
  !! @param[in,out] this      Instance of surface energy process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_seb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_energy(model_id))
    CALL this%Register_task(tsk_asselin_filter(model_id))
    CALL this%Register_task(tsk_surface_fluxes(model_id))
    CALL this%Register_task(tsk_energy_balance(model_id))

  END SUBROUTINE Register_seb_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "surface energy"
  !! Task "update_surface_energy" calculates the new surface temperature from the surface energy balance (sensible and latent heat,
  !! net radiation, ground heat).
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_energy(tile, options)

    USE mo_seb_land, ONLY: update_surface_energy_land
    USE mo_seb_lake, ONLY: update_surface_energy_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%is_lake) THEN
      CALL update_surface_energy_lake(tile, options)
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type '//TRIM(tile%lcts(1)%name)//' on tile '//TRIM(tile%name))
    END IF

  END SUBROUTINE update_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_energy"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_surface_energy(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_energy'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SEB_, t,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_old,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_unfilt,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_srf,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_eff4,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, qsat_star,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, s_star,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, forc_hflx,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, heat_cap,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, le_phase_change, weighted_by_fract)

    IF (model%config%use_lakes) THEN
      dsl4jsb_Aggregate_onChunk(SEB_, t_lwtr,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, s_lwtr,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, qsat_lwtr,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, fract_lice, weighted_by_fract)
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Aggregate_onChunk(SEB_, t_lice,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(SEB_, s_lice,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(SEB_, qsat_lice,  weighted_by_fract)
      END IF
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_energy

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "asselin_filter"
  !! Task "asselin_filter" calculates applies the Asselin time filter to the new surface temperature, if applicable
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_asselin(tile, options)

    USE mo_seb_land, ONLY: update_asselin_land

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_asselin'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%contains_lake) THEN
      ! No Asselin filter for lake tile
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_asselin_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_asselin_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_asselin

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "asselin_filter"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_asselin(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_asselin'

    INTEGER :: iblk , ics, ice
    TYPE(t_jsb_model), POINTER :: model

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SEB_, t,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_filt, weighted_by_fract)
    IF (model%config%use_tmx) THEN
      dsl4jsb_Aggregate_onChunk(SEB_, t_rad4, weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_asselin

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "surface fluxes"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_fluxes(tile, options)

    USE mo_seb_land, ONLY: update_surface_fluxes_land
    USE mo_seb_lake, ONLY: update_surface_fluxes_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%contains_lake) THEN
      CALL update_surface_fluxes_lake(tile, options)
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_fluxes"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_surface_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_fluxes'

    INTEGER :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx,      weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_lnd,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_lnd,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_wtr,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_wtr,  weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_ice,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_ice,  weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "update" for task "energy balance"
  !>
  !> This diagnostic computes the balance of the net surface radiation fluxes at the land surface.
  !>
  SUBROUTINE update_energy_balance(tile, options)

    USE mo_jsb_grid,               ONLY: Get_grid
    USE mo_jsb_grid_class,         ONLY: t_jsb_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & rad_srf_net,        &               !< Surface net radiation                               [W m-2]
      & lw_srf_down,        &               !< Longwave downwards radiation                        [W m-2]
      & swvis_srf_down,     &               !< Shortwave visible downwards radiation               [W m-2]
      & swnir_srf_down,     &               !< Shortwave near infrared downwards radiation         [W m-2]
      & swpar_srf_down,     &               !< Shortwave downwards photosynthetic active radiation [W m-2]
      & lw_srf_net,         &               !< Longwave net radiation on land surface              [W m-2]
      & sw_srf_net,         &               !< Shortwave net radiation on land surface             [W m-2]
      & sensible_hflx,      &               !< Sensible heat flux                                  [W m-2]
      & latent_hflx,        &               !< Latent heat flux                                    [W m-2]
      & forc_hflx,          &               !< Ground heat flux                                    [W m-2]
      & le_phase_change,    &               !< Latent energy required/released due to phase change [J m-2]
      & t,                  &               !< Current time steps temperature                          [K]
      & t_old,              &               !< Previous time steps temperature                         [K]
      & net_srf_rad_balance                 !< Net surface radiation balance                       [W m-2]

    REAL(wp), POINTER :: &
      & tile_fract(:), lat(:), lon(:)

    ! Locally allocated vectors
    !

    REAL(wp), DIMENSION(options%nc) :: &
      & lw_srf_up,       &
      & sw_srf_up,       &
      & sw_srf_down,     &
      & delta_t_srf

    INTEGER  :: iblk, ics, ice, nc, ic

    INTEGER  :: n_cells, ie
    INTEGER  :: cell_idx(options%nc)

    REAL(wp) :: dtime
    LOGICAL  :: is_experiment_start

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_energy_balance'
    CHARACTER(len=5000)         :: message_text_long

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    is_experiment_start = is_time_experiment_start(options%current_datetime)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    tile_fract => tile%fract(ics:ice, iblk)
    lat        => grid%lat(ics:ice, iblk)
    lon        => grid%lon(ics:ice, iblk)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      lw_srf_down)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      swvis_srf_down)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      swnir_srf_down)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      swpar_srf_down)             ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      rad_srf_net)                ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      lw_srf_net)                 ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      sw_srf_net)                 ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      latent_hflx)                ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      sensible_hflx)              ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      forc_hflx)                  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      le_phase_change)            ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                          ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_old)                      ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      net_srf_rad_balance)        ! out

    n_cells = 0   ! initialize counter for cells exceeding energy balance threshold

    !$ACC DATA CREATE(lw_srf_up, sw_srf_up, sw_srf_down, delta_t_srf)  ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      ! Compute fluxes used in log message
      lw_srf_up(ic)   = lw_srf_net(ic) - lw_srf_down(ic)
      sw_srf_down(ic) = swvis_srf_down(ic) + swnir_srf_down(ic) + swpar_srf_down(ic)
      sw_srf_up(ic)   = sw_srf_net(ic) - sw_srf_down(ic)

      ! Compute surface/skin temperature change
      delta_t_srf(ic) = t(ic) - t_old(ic)

      IF (.NOT. is_experiment_start .AND. tile_fract(ic) > EPSILON(1.0_wp)) THEN
        ! Compute energy flux residual (without surface/skin temperature change)
        net_srf_rad_balance(ic) = rad_srf_net(ic) + sensible_hflx(ic) + latent_hflx(ic) &
          &                     + forc_hflx(ic) + le_phase_change(ic) / dtime
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! Note that contrary to the water balance, the net surface radiation balance just
    ! checks the residual of the fluxes while ignoring the heat storages. Therefore,
    ! a large residual can still be acceptable if the heat storage is large enough to
    ! compensate for it, and large values do not necessarily indicate a problem.
#ifndef _OPENACC
    DO ic=1,nc
      IF (dsl4jsb_Config(SEB_)%eb_threshold > 0._wp) THEN
        IF (ABS(net_srf_rad_balance(ic)) > dsl4jsb_Config(SEB_)%eb_threshold .OR. &
          & ((t(ic) > 400._wp .OR. t(ic) < 50._wp).AND. tile_fract(ic) > EPSILON(1.0_wp))) THEN
          n_cells = n_cells + 1       ! Count number of cells with energy balance issues of this time step
          cell_idx(n_cells) = ic      ! Index of the respective grid cell (on chunk)
        END IF
      END IF
    END DO

    DO ie = 1, n_cells
      ic = cell_idx(ie)

      WRITE (message_text_long,*) 'Net surface radiation balance [W m-2]', NEW_LINE('a'), &
        & 'on ',TRIM(tile%name),' tile at', lat(ic),'N and ',lon(ic),'E',  NEW_LINE('a'), &
        & '(ic: ',ic,' iblk: ',iblk, ' tile_fract: ',tile_fract(ic),'):',  NEW_LINE('a'), &
        & 'Net surface radiation balance: ', net_srf_rad_balance(ic),      NEW_LINE('a'), &
        & 'Longwave downwards:            ', lw_srf_down(ic),              NEW_LINE('a'), &
        & 'Longwave upwards:              ', lw_srf_up(ic),                NEW_LINE('a'), &
        & 'Shortwave downwards:           ', sw_srf_down(ic),              NEW_LINE('a'), &
        & 'Shortwave upwards:             ', sw_srf_up(ic),                NEW_LINE('a'), &
        & 'Net surface radiation:         ', rad_srf_net(ic),              NEW_LINE('a'), &
        & 'Sensible heat flux:            ', sensible_hflx(ic),            NEW_LINE('a'), &
        & 'Latent heat flux:              ', latent_hflx(ic),              NEW_LINE('a'), &
        & 'Ground heat flux:              ', forc_hflx(ic),                NEW_LINE('a'), &
        & 'Phase change flux:             ', le_phase_change(ic) / dtime,  NEW_LINE('a'), &
        & 'Temperature change [K]:        ', delta_t_srf(ic)
      CALL message (TRIM(routine), message_text_long, all_print=.TRUE.)
    END DO
#endif

    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_energy_balance

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "energy_balance"
  !!
  SUBROUTINE aggregate_energy_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_energy_balance'

    INTEGER :: iblk

    iblk = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Don't aggregate, but explicitely compute energy balance on each tile
    CALL update_energy_balance(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_energy_balance

  !-----------------------------------------------------------------------------------------------------
  !> Make sure temperature is within a realistic range and finish the simulation with a meaningfull
  !! message if not.
  !!
  !! The routine is called at the beginning and at the end of each time step.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE seb_check_temperature_range(model_id, no_omp_thread)

    USE mo_jsb_grid,            ONLY: Get_grid
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid
    USE mo_jsb_parallel,        ONLY: Get_omp_thread

    ! Argument
    INTEGER, INTENT(in) :: model_id, no_omp_thread

    ! Local variables
    TYPE(t_jsb_model), POINTER           :: model
    TYPE(t_jsb_grid),  POINTER           :: grid
    CLASS(t_jsb_tile_abstract), POINTER  :: tile

    dsl4jsb_Def_memory(SEB_)

    INTEGER                       :: iblk, ic
    REAL(wp), POINTER             :: lat(:,:), lon(:,:), tile_fract(:,:)
    CHARACTER(len=*),  PARAMETER  :: routine = modname//':seb_check_temperature_range'

    dsl4jsb_Real2D_onDomain :: t

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => Get_model(model_id)
    grid  => Get_grid(model%grid_id)
    lat   => grid%lat(:,:)
    lon   => grid%lon(:,:)

    !vg no_omp_thread = Get_omp_thread()
#ifndef __QUINCY_STANDALONE__
    CALL model%Get_top_tile(tile)
    DO WHILE (ASSOCIATED(tile))
      IF ((tile%Has_children() .AND. .NOT. tile%visited(no_omp_thread)) .OR. .NOT. tile%Is_process_active(SEB_)) THEN
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF
#endif

      ! We're on a leaf or on the way up.

      dsl4jsb_Get_memory(SEB_)
      dsl4jsb_Get_var2D_onDomain(SEB_, t)

      tile_fract => tile%fract(:,:)
      IF (ANY(t(:,:) < 50._wp .OR. t(:,:) > 400._wp)) THEN
        DO iblk = 1, grid%nblks
          DO ic = 1, grid%nproma
            IF ((t(ic,iblk) < 50._wp .OR. t(ic,iblk) > 400._wp) .AND. tile_fract(ic,iblk) > 0._wp) THEN
              WRITE (message_text,*) 'Temperature out of bound: ', t(ic,iblk), 'K',             NEW_LINE('a'), &
                & 'on ',TRIM(tile%name),' tile at', lat(ic,iblk), 'N and ', lon(ic,iblk), 'E',  NEW_LINE('a'), &
                & '(ic: ',ic,' iblk: ',iblk, ' tile_fract:', tile_fract(ic,iblk),'):',          NEW_LINE('a'), &
                & 'One thing to check: consistency of land-sea mask and forcing.'
              CALL finish(TRIM(routine), TRIM(message_text))
            END IF
          END DO
        END DO
      END IF
#ifndef __QUINCY_STANDALONE__
      CALL model%Goto_next_tile(tile)
    ENDDO
#endif
  END SUBROUTINE seb_check_temperature_range

  !-----------------------------------------------------------------------------------------------------
  !> Calculations of diagnostic global land mean output
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
#ifndef __QUINCY_STANDALONE__
  SUBROUTINE global_seb_diagnostics(tile)

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_memory(SEB_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_seb_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: t

    REAL(wp), POINTER       :: t_gmean(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)
    REAL(wp)               :: global_land_area
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_var2D_onDomain(SEB_,  t)        ! in

    t_gmean        => SEB__mem%t_gmean%ptr(:)   ! out

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea


    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')


    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    IF (SEB__mem%t_gmean%is_in_output) THEN

      ! Domain Mask - to mask all halo cells for global sums (otherwise these
      ! cells are counted twice)
      ALLOCATE (in_domain(grid%nproma,grid%nblks))
      WHERE (is_in_domain(:,:))
        in_domain = 1._wp
      ELSEWHERE
        in_domain = 0._wp
      END WHERE

      ALLOCATE (scaling(grid%nproma,grid%nblks))

      ! Calculate global land mean seb variables
      global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))
      scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)
      t_gmean      = global_sum_array(t(:,:) * scaling(:,:)) / global_land_area

      DEALLOCATE (scaling, in_domain)
    END IF

  END SUBROUTINE global_seb_diagnostics
#endif

  !
  ! ================================================================================================================================
  !
#endif
END MODULE mo_seb_interface
