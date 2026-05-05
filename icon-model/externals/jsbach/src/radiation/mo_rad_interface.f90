!> Contains the interfaces to the radiation processes
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
!NEC$ options "-finline-file=externals/jsbach/src/radiation/mo_rad_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_rad_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,         ONLY: debug_on, acc_stream
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish

  USE mo_jsb_control,         ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, RAD_, SEB_, HYDRO_

  ! Use process configurations
  dsl4jsb_Use_config(RAD_)
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
#ifdef __QUINCY_STANDALONE__
#else
  ! Use of processes in this module
  dsl4jsb_Use_processes PHENO_, ASSIMI_, CARBON_

  ! Use process configurations
  dsl4jsb_Use_config(ASSIMI_)

  ! Use of process memories
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(CARBON_)
#endif
#ifndef __NO_QUINCY__
  dsl4jsb_Use_processes SB_, VEG_
  dsl4jsb_Use_memory(SB_)
  dsl4jsb_Use_memory(VEG_)
#endif
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_rad_tasks

#ifdef __QUINCY_STANDALONE__
  ! USEd by the mo_qs_model_interface
  PUBLIC :: update_surface_radiation, update_radiation_par, update_albedo, update_time_average_radiation
#endif

  !> Type definition for radiation_surface task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_radiation
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_radiation
  END TYPE tsk_surface_radiation

  INTERFACE tsk_surface_radiation
    PROCEDURE Create_task_surface_radiation
  END INTERFACE tsk_surface_radiation

  !> Type definition for radiation_par task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_radiation_par
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_par
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_par
  END TYPE tsk_radiation_par

  INTERFACE tsk_radiation_par
    PROCEDURE Create_task_radiation_par
  END INTERFACE tsk_radiation_par


  !> Type definition for albedo task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_albedo
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_albedo
    PROCEDURE, NOPASS :: Aggregate => aggregate_albedo
  END TYPE tsk_albedo

  INTERFACE tsk_albedo
    PROCEDURE Create_task_albedo
  END INTERFACE tsk_albedo

  ! ======================================================================================================= !
  !> time_average_radiation task
  !>
  ! TYPE definition
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_radiation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_radiation !< Aggregates computed task variables
  END TYPE tsk_time_average_radiation
  ! Constructor interface
  INTERFACE tsk_time_average_radiation
    PROCEDURE Create_task_update_time_average_radiation                !< Constructor function for task
  END INTERFACE tsk_time_average_radiation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_radiation task
  !>
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_surface_radiation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_radiation::return_ptr)
    CALL return_ptr%Construct(name='surface_radiation', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for radiation_par task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_radiation_par(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_radiation_par::return_ptr)
    CALL return_ptr%Construct(name='radiation_par', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_radiation_par

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for albedo task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "albedo"
  !!
  FUNCTION Create_task_albedo(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_albedo::return_ptr)
    CALL return_ptr%Construct(name='albedo', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_albedo

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for time_average_radiation
  !>
  FUNCTION Create_task_update_time_average_radiation(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_radiation::return_ptr)
    CALL return_ptr%Construct(name='rad_time_average', process_id=RAD_, owner_model_id=model_id)
  END FUNCTION Create_task_update_time_average_radiation

  ! ======================================================================================================= !
  !> Register tasks for radiation process
  !>
  SUBROUTINE Register_rad_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_surface_radiation       (model_id))
    CALL this%Register_task(tsk_radiation_par           (model_id))
    CALL this%Register_task(tsk_albedo                  (model_id))
    CALL this%Register_task(tsk_time_average_radiation(model_id))

  END SUBROUTINE Register_rad_tasks

  ! ================================================================================================================================
  !>
  !> Implementation to update the surface radiation
  !! Task "surface_radiation" calculates the net radiations at the surface.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_radiation(tile, options)

    USE mo_rad_process,             ONLY: calc_radiation_surface_net
    USE mo_jsb_physical_constants,  ONLY: zemiss_def, stbo

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    ! Declare pointers for variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & swvis_srf_net,        &
      & swnir_srf_net,        &
      & sw_srf_net,           &
      & lw_srf_net,           &
      & rad_srf_net,          &
      & lw_srf_down,          &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
      & alb_vis_lnd,          &
      & alb_nir_lnd,          &
      & fract_lice,           &
      & rad_net_lwtr,         &
      & rad_net_lice,         &
      & albedo_lwtr,          &
      & albedo_lice,          &
      & t,                    &
      & t_old,                &
      & t_lwtr,               &
      & t_lice

    REAL(wp) ::                  &
      & lw_net_lwtr(options%nc), &
      & lw_net_lice(options%nc), &
      & sw_net_lwtr(options%nc), &
      & sw_net_lice(options%nc)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: iblk , ics, ice, nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_radiation'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(RAD_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)    ! in

    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, lw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)    ! out

    ! ---------------------------
    ! Go

    !> 1.0 lakes
    !!
    IF (tile%is_lake) THEN

      ! Water
      dsl4jsb_Get_var2D_onChunk(SEB_, t_lwtr)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lwtr)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lwtr) ! out

      !$ACC DATA CREATE(sw_net_lwtr, lw_net_lwtr, sw_net_lice, lw_net_lice) ASYNC(acc_stream)

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),           & ! in
        & swnir_srf_down(:),           & ! in
        & albedo_lwtr(:),              & ! in
        & albedo_lwtr(:),              & ! in
        & lw_srf_down(:),              & ! in
        & t_lwtr(:),                   & ! in
        & rad_net_lwtr(:),             & ! out
        & sw_net=sw_net_lwtr(:),       & ! out
        & lw_net=lw_net_lwtr(:)        & ! out
        & )

      ! Ice on lake
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(SEB_, t_lice)       ! in
        dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice)   ! in
        dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lice)  ! in
        dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lice) ! out

        CALL calc_radiation_surface_net( &
          & swvis_srf_down(:),           & ! in
          & swnir_srf_down(:),           & ! in
          & albedo_lice(:),              & ! in
          & albedo_lice(:),              & ! in
          & lw_srf_down(:),              & ! in
          & t_lice(:),                   & ! in
          & rad_net_lice(:),             & ! out
          & sw_net=sw_net_lice(:),       & ! out
          & lw_net=lw_net_lice(:)        & ! out
          & )

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          rad_srf_net(ic) = (1._wp - fract_lice(ic)) * rad_net_lwtr(ic) + fract_lice(ic) * rad_net_lice(ic)
          sw_srf_net (ic) = (1._wp - fract_lice(ic)) * sw_net_lwtr(ic)  + fract_lice(ic) * sw_net_lice(ic)
          lw_srf_net (ic) = (1._wp - fract_lice(ic)) * lw_net_lwtr(ic)  + fract_lice(ic) * lw_net_lice(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          rad_srf_net(ic) = rad_net_lwtr(ic)
          sw_srf_net (ic) = sw_net_lwtr (ic)
          lw_srf_net (ic) = lw_net_lwtr (ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC END DATA

    !> 2.0 all but lakes
    !!
    !!
    ELSE
      dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! in
#ifndef __NO_QUINCY__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        dsl4jsb_Get_var2D_onChunk(SEB_,    t_old)          ! in
      END SELECT
#endif
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_lnd)    ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_lnd)    ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      swvis_srf_net)  ! out
      dsl4jsb_Get_var2D_onChunk(RAD_,      swnir_srf_net)  ! out

      CALL calc_radiation_surface_net(  &
        & swvis_srf_down(:),            & ! in
        & swnir_srf_down(:),            & ! in
        & alb_vis_lnd(:),               & ! in
        & alb_nir_lnd(:),               & ! in
        & lw_srf_down(:),               & ! in
        & t(:),                         & ! in
        & rad_srf_net(:),               & ! out
        & swvis_net = swvis_srf_net(:), & ! out
        & swnir_net = swnir_srf_net(:), & ! out
        & sw_net    = sw_srf_net(:),    & ! out
        & lw_net    = lw_srf_net(:)     & ! out
        & )

#ifndef __NO_QUINCY__
      ! re-calculate rad_srf_net for QUINCY
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          rad_srf_net(ic) = sw_srf_net(ic) + lw_srf_down(ic) - stbo * zemiss_def * t_old(ic) ** 4._wp
        END DO
        !$ACC END PARALLEL LOOP
      END SELECT
#endif
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_radiation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_radiation(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(RAD_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_radiation'

    INTEGER :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(RAD_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(RAD_, rad_srf_net,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, sw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, swvis_srf_net,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, swnir_srf_net,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, lw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lwtr,   weighted_by_fract)

    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lice,   weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_radiation

  ! ================================================================================================================================
  !>
  !> Task to update the photosynthetically active radiation for jsbach and quincy
  !>
  SUBROUTINE update_radiation_par(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER  :: model                        !< model
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_radiation_par'

    model        => Get_model(tile%owner_model_id)


    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      CALL update_radiation_par_jsb(tile, options)
#ifndef __NO_QUINCY__
    CASE (MODEL_QUINCY)
      CALL update_radiation_par_quincy(tile, options)
#endif
    END SELECT

  END SUBROUTINE update_radiation_par


  ! ======================================================================================================= !
  !> calculates PAR (Photosynthetically Active Photon Flux Density)
  !>
#ifndef __NO_QUINCY__
  SUBROUTINE update_radiation_par_quincy(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_process_class,     ONLY: A2L_, RAD_, VEG_, HYDRO_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_control,           ONLY: debug_on
    USE mo_jsb_math_constants,    ONLY: eps1, seconds_per_day
    USE mo_exception,             ONLY: message
    USE mo_rad_process,           ONLY: calc_canopy_radiation_budget
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                  !< the model
    TYPE(t_lctlib_element), POINTER   :: lctlib                 !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_canopy_q_assimi  !< Vertical grid
    INTEGER                           :: ncanopy                !< number of canopy layers
    INTEGER                           :: icanopy                !< loop 1,ncanopy
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    INTEGER                           :: ic                     !< looping
    REAL(wp)                          :: dtime                  !< timestep
    REAL(wp)                          :: fapar(options%nc)      !< fraction of absorbed photosynthetically active radiation
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_radiation_par_quincy'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: fract_par_diffuse
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk :: local_time_day_seconds
    ! RAD_ 2D
    dsl4jsb_Real2D_onChunk :: alb_vis_soil
    dsl4jsb_Real2D_onChunk :: fapar_noon
    ! RAD_ 3D
    dsl4jsb_Real3D_onChunk :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk :: ppfd_shaded_cl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: sai
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk :: lai_cl
    dsl4jsb_Real3D_onChunk :: cumm_lai_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(RAD_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! run only at PFT tiles (i.e. leafs of vegetation tile) because this routine is tailored towards vegetation only
    ! in the QUINCY usecases RAD_ is defined to run at the leafs of box tile, which includes also other tiles than PFT
    IF (tile%lcts(1)%lib_id == 0) THEN
      IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Run only at PFT tiles. '//TRIM(tile%name)//' is not a PFT. ')
      RETURN
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    lctlib                => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, fract_par_diffuse)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)         ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)   ! in
    ! ---------------------------
    ! RAD_ 2D
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil)             ! in
    dsl4jsb_Get_var2D_onChunk(RAD_, fapar_noon)               ! out
    ! RAD_ 3D
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)           ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)           ! out
    ! ---------------------------
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, sai)                      ! in
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)          ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                   ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)              ! in
    ! ----------------------------------------------------------------------------------------------------- !

    !$ACC DATA CREATE(fapar) ASYNC(acc_stream)

    ! Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
    CALL calc_canopy_radiation_budget(nc, ncanopy, &
      &                                dsl4jsb_Lctlib_param(sigma_vis), &
      &                                dsl4jsb_Lctlib_param(omega_clumping), &
      &                                dsl4jsb_Lctlib_param(crown_shape_factor), &
      &                                "vis", &
      &                                swpar_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &
      &                                lai(:), sai(:), lai_cl(:,:), cumm_lai_cl(:,:), &
      &                                alb_vis_soil(:), &
      &                                appfd_sunlit_cl = ppfd_sunlit_cl(:,:), &
      &                                appfd_shaded_cl = ppfd_shaded_cl(:,:), &
      &                                fleaf_sunlit_cl = fleaf_sunlit_cl(:,:), &
      &                                fapar = fapar(:))

    ! at noon: set alb_noon to current albedo (noon equals half a day plus one timestep)
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1,nc
      IF (ABS(local_time_day_seconds(ic) - ((seconds_per_day / 2._wp) + dtime)) < eps1) THEN
        fapar_noon(ic) = fapar(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(acc_stream)

    !$ACC END DATA
  END SUBROUTINE update_radiation_par_quincy
#endif
  ! ================================================================================================================================
  !>
  !> Implementation to update the photosynthetically active radiation
  !! Task "radiation_par" calculates PAR and APAR.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_radiation_par_jsb(tile, options)

    ! Use declarations
    USE mo_rad_process,    ONLY: calc_par

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
#ifndef __QUINCY_STANDALONE__
    ! Local variables
    INTEGER                           :: icanopy

    INTEGER                           :: iblk , ics, ice, nc, ic
    REAL(wp)                          :: dtime

    REAL(wp)                          :: canopy_bound_lai_tmp, canopy_bound_lai_delta
    REAL(wp), POINTER :: canopy_bound_lai(:)

    REAL(wp), DIMENSION(options%nc)   :: B4_layer_above ! R: cannot be scalar as it will be used as inout argum. for elemen. proce.
                                                        !    and thus has to have the same dimensions as the other out arguments.

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ncanopy               ! number of canopy layers
    !INTEGER, POINTER       :: acc_counter           ! counter for each delta_time step since the last variable output
    LOGICAL :: use_alb_veg_simple    !

    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: par                   ! Downward PAR flux in [W/m^2], Photosynth. Active Radiation (400-700nm)
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    ! dsl4jsb_Real2D_onChunk :: par_down_mol_nacc     ! PAR accumulated [MOL PHOTONS/M^2 ]
    ! dsl4jsb_Real2D_onChunk :: par_down_mol_tavg     ! PAR [MOL PHOTONS/M^2 S] (time mean value)
    dsl4jsb_Real2D_onChunk :: soil_reflectivity_par ! Soil reflectivity within the PAR spectrum [-]
    dsl4jsb_Real2D_onChunk :: alb_vis_soil          ! Soil albedo in the visible range
    dsl4jsb_Real2D_onChunk :: lai                   ! Leaf area index [-]
    dsl4jsb_Real3D_onChunk :: lai_cl                ! Leaf area index of canopy layer

    dsl4jsb_Real2D_onChunk :: cos_zenith_angle      ! Angle of incoming radiation [-]
    dsl4jsb_Real2D_onChunk :: fract_par_direct      ! Fraction of direct (in contrast to diffuse) PAR [-]
    dsl4jsb_Real2D_onChunk :: fract_par_diffuse     ! Fraction of diffuse (in contrast to direct) PAR [-]
    ! dsl4jsb_Real2D_onChunk :: apar_nacc
    dsl4jsb_Real2D_onChunk :: faPAR

    dsl4jsb_Real3D_onChunk :: faPAR_cl              ! Fraction of absorbed PAR per leaf area for each canopy layer
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Absorbed PAR of canopy layer
    ! dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_nacc  ! Absorbed PAR accumulated [MOL PHOTONS/M^2 ]
    ! dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_tavg  ! Absorbed PAR [MOL PHOTONS/M^2 S] (time mean value)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_radiation_par_jsb'

    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    IF (.NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get pointers to process configs and memory
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(A2L_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(RAD_)) RETURN

    ncanopy                     =  dsl4jsb_Config(ASSIMI_)%ncanopy
    use_alb_veg_simple          =  dsl4jsb_Config(RAD_)%use_alb_veg_simple
    canopy_bound_lai(0:ncanopy) => dsl4jsb_Config(ASSIMI_)% canopy_bound_lai(0:ncanopy)

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   cos_zenith_angle)      ! in (received from the atmosphere)
    dsl4jsb_Get_var2D_onChunk(A2L_,   fract_par_diffuse)     ! in (received from the atmosphere)
    dsl4jsb_Get_var2D_onChunk(A2L_,   swpar_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                   ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   fract_par_direct)      ! out
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_soil)          ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   par)                   ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   par_down_mol)          ! out

    ! R: These variables could also be local as they are not used outside this procedure in the moment.
    !    However I want it in the std model output and thus it has to be in the memory
    ! dsl4jsb_Get_var2D_onChunk(RAD_,   par_down_mol_nacc)     ! inout
    ! dsl4jsb_Get_var2D_onChunk(RAD_,   par_down_mol_tavg)     ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   soil_reflectivity_par) ! out
    dsl4jsb_Get_var3D_onChunk(RAD_,   lai_cl)                ! out
    dsl4jsb_Get_var3D_onChunk(RAD_,   faPAR_cl)              ! out
    ! used outside for assimilation:
    dsl4jsb_Get_var3D_onChunk(RAD_,   apar_per_lai_cl)       ! out
    ! dsl4jsb_Get_var3D_onChunk(RAD_,   apar_per_lai_cl_nacc)  ! inout
    ! dsl4jsb_Get_var3D_onChunk(RAD_,   apar_per_lai_cl_tavg)  ! inout
    ! dsl4jsb_Get_var2D_onChunk(RAD_,   apar_nacc)             ! inout

    dsl4jsb_Get_var2D_onChunk(RAD_,   faPAR)                 ! out

    ! ---------------------------
    ! Go

    !$ACC DATA CREATE(B4_layer_above) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
    DO icanopy=1,ncanopy
      DO ic=1,nc
        apar_per_lai_cl(ic,icanopy) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      par_down_mol(ic)          = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
      soil_reflectivity_par(ic) = 0._wp ! for now also in the memory but is only for current time step/"delta_time"

      par(ic) = swpar_srf_down(ic)

      ! Compute direct fraction of PAR.
      ! R: fract_par_direct = in ECHAM jsswdifpar, which is zsw_par_fract_diffuse = fract_diffuse_par  =sw_par_fract_direct
      fract_par_direct(ic) = 1._wp - fract_par_diffuse(ic)

      B4_layer_above(ic) = 7777777._wp ! Used to transfer variable B4 from one canopy layer to the next. Startvalue only to indicate
                                      ! first layer
      fapar(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    ! Compute faPAR_cl (absorbed photsynthetic active radiation for each canopy layer) and
    ! lai_cl (leaf area index for each canopy layer)
    DO icanopy=1,ncanopy
        canopy_bound_lai_tmp = canopy_bound_lai(icanopy)
        canopy_bound_lai_delta = canopy_bound_lai(icanopy) - canopy_bound_lai(icanopy-1)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
       !NEC$ ivdep
        DO ic=1,nc
          CALL calc_par(                 &
            & ncanopy,                   & ! in
            & icanopy,                   & ! in
            & use_alb_veg_simple,        & ! in
            & par(ic),                    & ! in
            & alb_vis_soil(ic),           & ! in
            & lai(ic),                    & ! in
            & cos_zenith_angle(ic),       & ! in
            & fract_par_direct(ic),       & ! in
            & canopy_bound_lai_tmp,       & ! in
            & canopy_bound_lai_delta,    & ! in
            & B4_layer_above(ic),         & ! inout
            & par_down_mol(ic),           & ! out
            & soil_reflectivity_par(ic),  & ! out
            & faPAR_cl(ic,icanopy),       & ! out
            & lai_cl(ic,icanopy),         & ! out
            & apar_per_lai_cl(ic,icanopy) ) ! out, used outside for assimilation
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    ! R: Note, the following code is preliminary and should be replaced later by an laccu-procedure.
    !    For the correct integration over time as it is done now, this task must be called only once per timestep.
    ! Accumulate incoming par in mol(photons)/(m^2) over time
    ! par_down_mol_nacc(:) = par_down_mol_nacc(:) + par_down_mol(:) * dtime

    ! R: This whole workaround should get abitary soon...
    !acc_counter = acc_counter + 1
    !par_down_mol_tavg = par_down_mol_nacc / dtime / acc_counter
    ! apar_per_lai_cl_nacc(:,:)  = apar_per_lai_cl_nacc(:,:) + apar_per_lai_cl(:,:) * dtime

    ! R: This whole workaround (acc_counter) should get abitary soon...
    !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
    ! DO icanopy=1,ncanopy
      ! apar_nacc(:) = apar_nacc(:) + par_down_mol(:) * faPAR_cl(:,icanopy) * dtime
    ! END DO

    ! Create faPAR
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO icanopy=1,ncanopy
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        ! faPAR_cl is always the fraction relativ to incoming radiation on the top of the whole canopy
        faPAR(ic) = faPAR(ic) + faPAR_cl(ic,icanopy)
      END DO
    END DO
    !$ACC END PARALLEL


!        ! Compute absorbed PAR per leaf area in canopy layer [units: (absorbed photons) / (m^2(leaf area) s)] from
!        ! par and fraction of absorbed PAR
!        apar_per_lai_cl(:,icanopy) = par_down_mol(:) * faPAR_cl(:,icanopy) / (MAX(lai_cl(:,icanopy),1.e-10_wp))
!        ! R: Folgendes faellt weg:
!        ! mask(1:nidx,itile) <--theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles)
!        ! apar_per_lai_cl(:,icanopy) = MERGE(apar_per_lai_cl(:,icanopy),0._dp,mask(1:nidx,itile))
!        ! Da damit die Wuestenflaechen auf 0 gesetzt werden und JSBACH4 diese Prozedur dann
!        ! gar nicht erst aufruft fuer diese tile.
!        IF (waterLimitationFlag) THEN ! to make sure to accumulate only once per timestep
!           apar_per_lai_cl_nacc(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) + apar_per_lai_cl(:,icanopy) * dtime
!           ! R: This whole workaround (acc_counter) should get abitary soon...
!           !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
!           apar_nacc(:)                     = apar_nacc(:) + par_down_mol(:) * faPAR_cl(:,icanopy) * dtime
!        END IF

    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
#endif
  END SUBROUTINE update_radiation_par_jsb

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "radiation_par"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_radiation_par(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(RAD_)
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    TYPE(t_jsb_model),       POINTER :: model
    INTEGER                          :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_radiation_par'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Get_memory(RAD_)

    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      ! RAD_ 2D
      dsl4jsb_Aggregate_onChunk(RAD_, fapar_noon,            weighted_by_fract)
    END SELECT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_radiation_par

  ! ================================================================================================================================
  !>
  !> Implementation to update the albedo
  !! Task "albedo" calculates several albedos.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_albedo(tile, options)

    ! Use declarations
    USE mo_jsb_lctlib_class,    ONLY: t_lctlib_element
    USE mo_rad_process,         ONLY: &
      & get_surface_albedo_simple, calc_soil_albedo, calc_pond_albedo, calc_alb_lwtr, calc_alb_lice, &
      & calc_glacier_albedo, calc_sky_view_fractions, calc_snow_albedo, Merge_albedos_of_vegtile, &
      & Has_minimal_radiation, calc_snow_age, get_litter_carbon, get_soil_carbon
#ifndef __NO_QUINCY__
    USE mo_rad_process,         ONLY: calc_canopy_radiation_budget
#endif
    USE mo_jsb_time,            ONLY: is_time_experiment_start, is_time_ltrig_rad_m1
    USE mo_jsb_math_constants,  ONLY: eps8, eps1, seconds_per_day
    USE mo_jsb_control,         ONLY: jsbach_runs_standalone
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_vgrid
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_lctlib_element), POINTER :: lctlib                !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER :: vgrid_soil_w          !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER :: vgrid_canopy_q_assimi !< Vertical grid
    INTEGER                         :: nsoil_w      !< number of soil layers as used/defined by the SB_ process
    INTEGER                         :: ncanopy      !< number of canopy layers as used/defined by the Q_ASSIMI process
    INTEGER                         :: iblk, ic, is
    INTEGER                         :: ics, ice, nc
    INTEGER                         :: model_scheme
    LOGICAL, ALLOCATABLE            :: l_day(:)
    REAL(wp)                        :: dtime                     !< timestep
    REAL(wp)                        :: lctlib_specific_leaf_area !< use for jsbach & quincy, as helper variable
    REAL(wp), ALLOCATABLE           :: temperature_surface(:)    !< use for jsbach & quincy, as helper variable
    CHARACTER(len=*), PARAMETER     :: routine = modname//':update_albedo'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(CARBON_)
#endif
#ifndef __NO_QUINCY__
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SB_)
#endif
    ! Declare configuration variables
    CHARACTER(len=:), ALLOCATABLE :: use_alb_soil_organic_C
    INTEGER :: use_alb_soil_organic_int

    REAL(wp) ::            &
      & albedo_age_weight, &    ! 0 < albedo_age_weight < 1: snow albedo is calculated by linearly weighting
                                !                            the snow albedo resulting from both schemes
      & AlbedoCanopySnow

    LOGICAL ::                     &
      & use_alb_veg_simple,        &
      & use_alb_soil_scheme,       &
      & use_alb_soil_litter,       &
      & use_alb_mineralsoil_const

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  &
      & cos_zenith_angle,      &
      & t_unfilt,              &
      & snow,                  &
      & alb_vis,               &
      & alb_nir,               &
      & alb_vis_lnd,           &
      & alb_nir_lnd,           &
      & alb_vis_can,           &
      & alb_nir_can,           &
      & alb_vis_soil,          &
      & alb_nir_soil,          &
      & alb_vis_soil_and_snow, &
      & alb_nir_soil_and_snow, &
      & alb_vis_snow,          &
      & alb_nir_snow,          &
      & alb_vis_pond,          &
      & alb_nir_pond,          &
      & alb_background,        &
      & alb_vis_mineralsoil,   &
      & alb_nir_mineralsoil,   &
      & sky_view_fract,        &
      & sky_view_fract_stem,   &
      & snow_age,              &
      & lai,                   &
      & fract_fpc_max,         &
      & fract_forest,          &
      & fract_snow_can,        &
      & fract_snow_soil,       &
      & fract_pond,            &
      & wtr_pond,              &
      & ice_pond,              &
      & swvis_srf_down,        &
      & swnir_srf_down,        &
      & fract_par_diffuse,     &
      & sai
    ! Lakes
    dsl4jsb_Real2D_onChunk ::  &
      & albedo_lwtr,           &
      & fract_lice,            &
      & t_lice,                &
      & weq_snow_lice,         &
      & albedo_lice

    ! Carbon pools needed with alb_soil_scheme
    dsl4jsb_Real2D_onChunk ::  &
      c_acid_ag1,              &
      c_water_ag1,             &
      c_ethanol_ag1,           &
      c_nonsoluble_ag1,        &
      c_acid_bg1,              &
      c_water_bg1,             &
      c_ethanol_bg1,           &
      c_nonsoluble_bg1,        &
      c_acid_bg2,              &
      c_water_bg2,             &
      c_ethanol_bg2,           &
      c_nonsoluble_bg2,        &
      c_humus_1,               &
      c_humus_2

    REAL(wp), ALLOCATABLE ::   &
      & litter_carbon(:),      &  ! Amount of litter carbon [mol(C)/m2(canopy)], used with use_alb_soil_litter
      & soil_carbon(:)            ! Amount of soil carbon [mol(C)/m2(canopy)], used with use_alb_soil_organic_C

    ! quincy
    dsl4jsb_Real2D_onChunk :: &
      & sb_pool_total_ag_litter_c, &
      & t, &
      & num_sl_above_bedrock, &
      & albedo_noon, &
      & alb_vis_noon, &
      & alb_nir_noon, &
      & local_time_day_seconds
    dsl4jsb_Real3D_onChunk ::       &
      & sb_pool_total_bg_soil_c,    &
      & soil_depth_sl,              &
      & lai_cl,                     &
      & cumm_lai_cl

    REAL(wp) :: &
      & specificLeafArea_C_param, &
      & AlbedoLitterVIS_param,    &
      & AlbedoLitterNIR_param,    &
      & StemArea_param
    LOGICAL :: &
      & ForestFlag_param,         &
      & l_ponds_param

    ! Set local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set model
    model => Get_model(tile%owner_model_id)
    model_scheme = model%config%model_scheme

#ifndef __NO_QUINCY__
    SELECT CASE (model_scheme)
    CASE (MODEL_QUINCY)
      vgrid_soil_w  => Get_vgrid('soil_depth_water')
      nsoil_w       =  vgrid_soil_w%n_levels
    END SELECT
#endif

#ifndef __QUINCY_STANDALONE__
    ! Update albedo only at experiment start, if model runs standalone or one time step before the radiation time step
    IF (.NOT. is_time_experiment_start(options%current_datetime) .AND. &
      & .NOT. is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id, jsbach_runs_standalone()) .AND. &
      & .NOT. jsbach_runs_standalone()) THEN
      IF (debug_on() .AND. iblk == 1) THEN
        WRITE (message_text,*) 'Update albedo only at experiment start, ', &
          &                    'if model runs standalone or one time step before the radiation time step ', &
          &                    'RETURN'
        CALL message(TRIM(routine), message_text)
      END IF
      RETURN
    END IF
#endif

    ! Enable for debugging:
    ! IF (iblk == 1 .AND. .NOT. ASSOCIATED(tile%next_sibling_tile)) THEN
    !   IF (is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id, jsbach_runs_standalone()) .OR. &
    !     & is_time_experiment_start(options%current_datetime)) THEN
    !     CALL message('update_albedo ', 'compute surface albedo')
    !   END IF
    ! END IF

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    IF (tile%contains_vegetation) THEN
      SELECT CASE (model_scheme)
#ifndef __QUINCY_STANDALONE__
      CASE (MODEL_JSBACH)
        dsl4jsb_Get_memory(PHENO_)
#endif
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        dsl4jsb_Get_memory(VEG_)
#endif
      END SELECT
    END IF

    ! Set pointers to variables in config
    use_alb_veg_simple        = dsl4jsb_Config(RAD_) %use_alb_veg_simple
    use_alb_soil_scheme       = dsl4jsb_Config(RAD_) %use_alb_soil_scheme
    use_alb_mineralsoil_const = dsl4jsb_Config(RAD_) %use_alb_mineralsoil_const
    use_alb_soil_litter       = dsl4jsb_Config(RAD_) %use_alb_soil_litter
    use_alb_soil_organic_C    = dsl4jsb_Config(RAD_) %use_alb_soil_organic_C
    albedo_age_weight         = dsl4jsb_Config(RAD_) %albedo_age_weight
    AlbedoCanopySnow          = dsl4jsb_Config(RAD_) %AlbedoCanopySnow

    l_ponds_param             = dsl4jsb_Config(HYDRO_)%l_ponds

    IF(use_alb_soil_organic_C .eq. 'linear') THEN
       use_alb_soil_organic_int = 1
    ELSE IF(use_alb_soil_organic_C .eq. 'log') THEN
       use_alb_soil_organic_int = 2
    END IF


    ! Prepare carbon arrays needed with soil albedo scheme

    IF (use_alb_soil_scheme .AND. (tile%is_bare .OR. tile%is_vegetation)) THEN

      ALLOCATE(litter_carbon(nc))
      ALLOCATE(soil_carbon(nc))
      !$ACC ENTER DATA CREATE(litter_carbon, soil_carbon) ASYNC(acc_stream)

      SELECT CASE (model_scheme)
#ifndef __QUINCY_STANDALONE__
      CASE (MODEL_JSBACH)
        IF (.NOT. tile%Is_process_calculated(CARBON_)) THEN
          CALL finish(TRIM(routine), 'Carbon process must be active with use_alb_soil_scheme')
        END IF

        dsl4jsb_Get_memory(CARBON_)

        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! in

        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! in

        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)            ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)           ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)         ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)      ! in

        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! in
        dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)             ! in

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          litter_carbon(ic) =         &  ! above ground green litter and green humus C pools
            &  get_litter_carbon(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), &
            &                    c_humus_1(ic))
          soil_carbon(ic) =           &  ! below ground litter and both humus C pools
            &    get_soil_carbon(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), &
            &                    c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), &
            &                    c_humus_1(ic), c_humus_2(ic))
        END DO
        !$ACC END PARALLEL LOOP
#endif
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        dsl4jsb_Get_memory(SB_)
        dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)      ! in
        dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)             ! in
        dsl4jsb_Get_var3D_onChunk(SB_,    sb_pool_total_bg_soil_c)   ! in
        dsl4jsb_Get_var2D_onChunk(SB_,    sb_pool_total_ag_litter_c) ! in

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          litter_carbon(ic) = sb_pool_total_ag_litter_c(ic)
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          soil_carbon(ic) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
        !$ACC LOOP SEQ
        DO is = 1, nsoil_w
          !$ACC LOOP GANG VECTOR
          DO ic = 1, nc
            IF (is <= NINT(num_sl_above_bedrock(ic))) THEN
              soil_carbon(ic) = soil_carbon(ic) + (sb_pool_total_bg_soil_c(ic, is) * soil_depth_sl(ic, is))
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
#endif
      END SELECT
    END IF   ! perparations for alb_soil_scheme

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   cos_zenith_angle)      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   swvis_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   swnir_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   fract_par_diffuse)     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)                  ! in

    ALLOCATE(temperature_surface(nc))
    !$ACC ENTER DATA CREATE(temperature_surface) ASYNC(acc_stream)
    SELECT CASE (model_scheme)
#ifndef __QUINCY_STANDALONE__
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)              ! in
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1,nc
        temperature_surface(ic) = t_unfilt(ic)
      END DO
      !$ACC END PARALLEL LOOP
#endif
#ifndef __NO_QUINCY__
    CASE (MODEL_QUINCY)
      dsl4jsb_Get_var2D_onChunk(SEB_,   t)                   ! in
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1,nc
        temperature_surface(ic) = t(ic)
      END DO
      !$ACC END PARALLEL LOOP

      dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_noon)           ! out
      dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_noon)          ! out
      dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_noon)          ! out
#endif
    END SELECT

    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_background)        ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_snow)          ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_snow)          ! inout
    IF (.NOT. tile%is_lake) THEN
      dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_lnd)         ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_lnd)         ! inout
      IF (.NOT. tile%is_glacier) THEN
        dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_pond)      ! inout
        dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_pond)      ! inout
      END IF
    END IF

    ! ---------------------------
    ! Go

    !> 0.9  Set day-or-night-switch
    !!
    ALLOCATE(l_day(nc))
    !$ACC ENTER DATA CREATE(l_day) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      l_day(ic) = Has_minimal_radiation(swvis_srf_down(ic), swnir_srf_down(ic))
    END DO
    !$ACC END PARALLEL LOOP

    !> 1.0 snow and pont albedo
    !!
    !! Calculate snow albedos. Used later only for bare and vegetation tile types
    !! Note, for vegetation tiles this snow albedos are used only for snow on soil
    IF (tile%is_bare .OR. tile%is_vegetation) THEN

      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)  ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_pond)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)         ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,   snow_age)         ! in

      ! Update (non-dimensional) snow age
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_snow_age(          &
          & dtime,                   & ! in
          & temperature_surface(ic), & ! in
          & snow(ic),                & ! in
          & fract_snow_soil(ic),     & ! in
          & snow_age(ic)             & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_snow_albedo(    &
          & l_day(ic),            & ! in
          & fract_snow_soil(ic),  & ! in
          & snow_age(ic),         & ! in
          & cos_zenith_angle(ic), & ! in
          & albedo_age_weight,    & ! in
          & temperature_surface(ic), & ! in
          & alb_vis_snow(ic),     & ! inout
          & alb_nir_snow(ic)      & ! inout
          & )

        IF (l_ponds_param) THEN
          CALL calc_pond_albedo( &
            & l_day(ic),            & ! in
            & cos_zenith_angle(ic), & ! in
            & temperature_surface(ic), & ! in
            & wtr_pond(ic),         & ! in
            & ice_pond(ic),         & ! in
            & alb_vis_pond(ic),     & ! inout
            & alb_nir_pond(ic)      & ! inout
            & )
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END IF  ! IF tile%is_bare .OR. tile%is_vegetation

    !> 2.0 lakes & lakes with ice
    !!
    !!
    IF (tile%is_lake) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,   t_lice)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_lwtr)  ! inout

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_alb_lwtr( &
          & l_day(ic),       & ! in
          & albedo_lwtr(ic)  & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

      ! Ice on lake
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_lice)  ! in
        dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)     ! in
        dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_lice)    ! inout

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          CALL calc_alb_lice(    &
            & l_day(ic),         & ! in
            & t_lice(ic),        & ! in
            & weq_snow_lice(ic), & ! in
            & albedo_lice(ic)    & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          IF(l_day(ic)) THEN
            alb_vis(ic) = (1._wp - fract_lice(ic)) * albedo_lwtr(ic) + fract_lice(ic) * albedo_lice(ic)
            alb_nir(ic) = alb_vis(ic)  ! TODO
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          IF(l_day(ic)) THEN
            alb_vis(ic) = albedo_lwtr(ic)
            alb_nir(ic) = alb_vis(ic)  ! TODO
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    !> 3.0 glacier
    !!
    ELSE IF (tile%is_glacier) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_glacier_albedo(  &
          & l_day(ic),             & ! in
          & temperature_surface(ic), & ! in
          & alb_vis_lnd(ic),       & ! inout
          & alb_nir_lnd(ic)        & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

    !> 4.0 soil
    !! currently jsbach4 does not consider "bare tiles"
    !!
    ELSE IF (tile%is_bare) THEN
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_soil)          ! inout

      !!       needs the active CARBON_ process
      !!       by default: use_alb_soil_scheme = .FALSE.
      !!       also the soil albedo is calc below for vegetation tiles
      IF (use_alb_soil_scheme) THEN
        IF (tile%lcts(1)%lib_id == 0) CALL finish(TRIM(routine), 'lctlib not available on tile')

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        !NEC$ ivdep
        DO ic=1,nc
          CALL calc_soil_albedo(                        &
            & l_day(ic),                                & ! in
            & use_alb_mineralsoil_const,                & ! in
            & use_alb_soil_organic_int,                 & ! in
            & use_alb_soil_litter,                      & ! in
            & litter_carbon(ic),                        & ! in
            & soil_carbon(ic),                          & ! in
            ! @todo: on a bare soil tile, there is no vegetation and therefore also no lctlib entry!?
            & lctlib_specific_leaf_area,                & ! in
            & dsl4jsb_Lctlib_param(AlbedoLitterVIS),    & ! in
            & dsl4jsb_Lctlib_param(AlbedoLitterNIR),    & ! in
            & alb_background(ic),                       & ! in
            & alb_vis_mineralsoil(ic),                  & ! inout
            & alb_nir_mineralsoil(ic),                  & ! inout
            & alb_vis_soil(ic),                         & ! inout
            & alb_nir_soil(ic)                          & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from bc_land_phy.nc file
        ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
        ! the canopy above.
        IF (debug_on() .AND. iblk == 1) THEN
          CALL message(TRIM(routine), 'Soil albedo scheme not used for bare soil tiles but taken from bc_land_phy.nc file')
        END IF
      END IF

      ! Merge the albedo of soil with the snow on it for the overall bare soil tile albedo
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        alb_vis_lnd(ic) = fract_snow_soil(ic)                                       * alb_vis_snow(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * fract_pond(ic)            * alb_vis_pond(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * (1.0_wp - fract_pond(ic)) * alb_vis_soil(ic)

        alb_nir_lnd(ic) = fract_snow_soil(ic)                                       * alb_nir_snow(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * fract_pond(ic)            * alb_nir_pond(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * (1.0_wp - fract_pond(ic)) * alb_nir_soil(ic)
      END DO
      !$ACC END PARALLEL LOOP

    !> 5.0 vegetation / canopy
    !!
    ELSE IF (tile%is_vegetation) THEN

      ! TODO: correct indentation >>>
    SELECT CASE (model_scheme)
#ifndef __QUINCY_STANDALONE__
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)                   ! in
      dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max)         ! in
      IF (use_alb_veg_simple) THEN
        dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_forest)          ! in
      ENDIF
#endif
    END SELECT
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_can)           ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_can)           ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_soil)          ! inout

      !> 5.1  simple calculation of alb_vis_lnd / alb_nir_lnd
      !!
      !! TRUE for jsbach_lite, FALSE for jsbach_pfts and quincy use cases
      IF (use_alb_veg_simple) THEN
#ifndef __QUINCY_STANDALONE__
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          alb_vis_lnd(ic) = get_surface_albedo_simple( &
            & t_unfilt(ic),                            &
            & fract_forest(ic),                        &
            & fract_snow_soil(ic),                     &
            & fract_snow_can(ic),                      &
            & fract_pond(ic),                          &
            & lai(ic),                                 &
            & alb_background(ic)                       & ! Note, alb_background <-- bc_land_phy.nc
                        )
          alb_nir_lnd(ic) = alb_vis_lnd(ic)
        END DO
        !$ACC END PARALLEL LOOP
#endif
      !> 5.2  detailed calculation of alb_vis_lnd & alb_nir_lnd
      !!
      !! merging albedo of snow, soil and canopy
      ELSE

        !> 5.2.1 Get bare soil albedo for vegetation tiles
        !!       needs the active CARBON_ process
        !!       by default: use_alb_soil_scheme = .FALSE.
        !!       also the soil albedo is calc above with the same routine in case "IF (tile%is_bare)"
        IF (use_alb_soil_scheme) THEN
          ! The albedo for bare soil under the canopy of a vegetation tile is the same as for the bare soil tile:

          IF (tile%lcts(1)%lib_id == 0) CALL finish(TRIM(routine), 'lctlib not available on tile')
          SELECT CASE (model_scheme)
#ifndef __QUINCY_STANDALONE__
          CASE (MODEL_JSBACH)
            specificLeafArea_C_param = dsl4jsb_Lctlib_param(specificLeafArea_C)
#endif
#ifndef __NO_QUINCY__
          CASE (MODEL_QUINCY)
            specificLeafArea_C_param = dsl4jsb_Lctlib_param(sla)
#endif
          END SELECT
          AlbedoLitterVIS_param    = dsl4jsb_Lctlib_param(AlbedoLitterVIS)
          AlbedoLitterNIR_param    = dsl4jsb_Lctlib_param(AlbedoLitterNIR)

          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          !NEC$ ivdep
          DO ic=1,nc
            CALL calc_soil_albedo(                        &
              & l_day(ic),                                & ! in
              & use_alb_mineralsoil_const,                & ! in
              & use_alb_soil_organic_int,                 & ! in
              & use_alb_soil_litter,                      & ! in
              & litter_carbon(ic),                        & ! in
              & soil_carbon(ic),                          & ! in
              & specificLeafArea_C_param,                 & ! in
              & AlbedoLitterVIS_param,                    & ! in
              & AlbedoLitterNIR_param,                    & ! in
              & alb_background(ic),                       & ! in
              & alb_vis_mineralsoil(ic),                  & ! inout
              & alb_nir_mineralsoil(ic),                  & ! inout
              & alb_vis_soil(ic),                         & ! inout
              & alb_nir_soil(ic)                          & ! inout
              & )
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from the bc_land_phy.nc file
          ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
          ! the canopy above.

          IF (debug_on() .AND. iblk == 1) THEN
            CALL message(TRIM(routine), 'Soil albedo scheme not used for veg tiles but taken from bc_land_phy.nc file')
          END IF
        END IF

        SELECT CASE (model_scheme)
#ifndef __QUINCY_STANDALONE__
        ! TODO: adapt indent for jsbach case
        CASE (MODEL_JSBACH)
        ! Get the fractions of the vegetation tile that is bare soil from the sky point of view
        ! Needed later for merging the contributing snow-, bare soil- and canopy-albedos for this vegetation tile
        dsl4jsb_Get_var2D_onChunk(RAD_, sky_view_fract)       ! inout
        dsl4jsb_Get_var2D_onChunk(RAD_, sky_view_fract_stem)  ! inout

        ! double check if tile%is_vegetation is PFT
        IF (tile%lcts(1)%lib_id == 0) CALL finish(TRIM(routine), 'lctlib not available on tile')

        !> 5.2.2 calc sky-view fractions
        !!
        StemArea_param = dsl4jsb_Lctlib_param(StemArea)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          CALL calc_sky_view_fractions(       &
            & l_day(ic),                      & ! in
            & lai(ic),                        & ! in
            & fract_fpc_max(ic),              & ! in
            & cos_zenith_angle(ic),           & ! in
            & StemArea_param,                 & ! in (Area of stems and branches of woody plants.)
            & sky_view_fract(ic),             & ! inout
            & sky_view_fract_stem(ic)         & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP

        !> 5.2.3 merge albedos of snow, soil and canopy into land
        !!
        ForestFlag_param = dsl4jsb_Lctlib_param(ForestFlag)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          CALL Merge_albedos_of_vegtile(        &
            & l_day(ic),                        & ! in
            & ForestFlag_param,                 & ! in  ("true" for forest landcover types)
            & sky_view_fract(ic),               & ! in
            & sky_view_fract_stem(ic),          & ! in
            & AlbedoCanopySnow,                 & ! in
            & fract_snow_soil(ic),              & ! in
            & fract_snow_can(ic),               & ! in
            & fract_pond(ic),                   & ! in
            & alb_vis_snow(ic),                 & ! in
            & alb_nir_snow(ic),                 & ! in
            & alb_vis_soil(ic),                 & ! in
            & alb_nir_soil(ic),                 & ! in
            & alb_vis_can(ic),                  & ! in
            & alb_nir_can(ic),                  & ! in
            & alb_vis_pond(ic),                 & ! in
            & alb_nir_pond(ic),                 & ! in
            & alb_vis_lnd(ic),                  & ! inout
            & alb_nir_lnd(ic)                   & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP
#endif
#ifndef __NO_QUINCY__
        CASE (MODEL_QUINCY)
          vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
          ncanopy               =  vgrid_canopy_q_assimi%n_levels
          dsl4jsb_Get_var2D_onChunk(VEG_, sai)                   ! in
          dsl4jsb_Get_var2D_onChunk(VEG_, lai)                   ! in
          dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                ! in
          dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)           ! in
          dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_soil_and_snow) ! out
          dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil_and_snow) ! out
          ! for now, prescribed, should be based on JSBACH 4 routines, and include the effect of snow on the
          ! soil's albedo
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1, nc
            alb_vis_soil_and_snow(ic) = fract_snow_soil(ic) * alb_vis_snow(ic) &
              &                + (1.0_wp - fract_snow_soil(ic)) * alb_vis_soil(ic)
            alb_nir_soil_and_snow(ic) = fract_snow_soil(ic) * alb_nir_snow(ic) &
              &                + (1.0_wp - fract_snow_soil(ic)) * alb_nir_soil(ic)
          END DO
          !$ACC END PARALLEL LOOP

          !> 1.0 calculate VIS radiation budget (absorbed radiation on canopy layers and soil + albedo)
          ! Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
          CALL calc_canopy_radiation_budget(nc, ncanopy, &
            &                                dsl4jsb_Lctlib_param(sigma_vis), &
            &                                dsl4jsb_Lctlib_param(omega_clumping), &
            &                                dsl4jsb_Lctlib_param(crown_shape_factor), &
            &                                "vis", &
            &                                swvis_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &
            &                                lai(:), sai(:), lai_cl(:,:), cumm_lai_cl(:,:), &
            &                                alb_vis_soil_and_snow(:), &
            &                                alb_can = alb_vis_can(:), alb_lnd = alb_vis_lnd(:))

          ! correct surface albedo by canopy fraction covered with snow
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic=1,nc
            alb_vis_lnd(ic) = alb_vis_lnd(ic) * (1._wp - fract_snow_can(ic)) &
              &               + AlbedoCanopySnow * fract_snow_can(ic)
          END DO
          !$ACC END PARALLEL LOOP

          !> 2.0 calculate NIR radiation budget (absorbed radiation on canopy layers and soil + albedo)
          !!
          !! Following Weiss & Norman, Agri. For. Met. 1985, the fraction of VIS and NIR
          !! diffuse radiation are very similar
          ! Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
          CALL calc_canopy_radiation_budget(nc, ncanopy, &
            &                                dsl4jsb_Lctlib_param(sigma_nir), &
            &                                dsl4jsb_Lctlib_param(omega_clumping), &
            &                                dsl4jsb_Lctlib_param(crown_shape_factor), &
            &                                "nir", &
            &               swnir_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &
            &               lai(:), sai(:), lai_cl(:,:), cumm_lai_cl(:,:), &
            &               alb_nir_soil_and_snow(:), &
            &               alb_can = alb_nir_can(:), alb_lnd = alb_nir_lnd(:))

          ! correct canopy albedo by fraction covered with snow
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic=1,nc
            alb_nir_lnd(ic) = alb_nir_lnd(ic) * (1._wp - fract_snow_can(ic)) &
              &               + AlbedoCanopySnow * fract_snow_can(ic)
          END DO
          !$ACC END PARALLEL LOOP

#endif
        END SELECT  ! quincy / jsbach
      END IF      ! NOT use_alb_veg_simple
    END IF      ! IF (tile%is_vegetation)

    !> 6.0 pass the "land albedo" to the "tile albedo"
    !! this is needed because for lake tiles, the "tile albedo" is calculated directly
    !! but for vegetation, bare and glacier, the "land albedo" is calculated
    !!
#ifndef __NO_QUINCY__
    IF (model_scheme == MODEL_QUINCY) THEN
      dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)         ! in
    END IF
#endif
    IF (tile%is_bare .OR. tile%is_glacier .OR. tile%is_vegetation) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        ! for lake tile, alb_vis/nur are already set above to albedo_lwtr
        alb_vis(ic) = alb_vis_lnd(ic)
        alb_nir(ic) = alb_nir_lnd(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

#ifndef __NO_QUINCY__
    ! calculate albedo and albedo_noon for all tiles
    !
    IF (model_scheme == MODEL_QUINCY) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        ! at noon: set alb_noon to current albedo (noon equals half a day plus one timestep)
        !
        IF (ABS(local_time_day_seconds(ic) - ((seconds_per_day / 2._wp) + dtime)) < eps1) THEN
          alb_vis_noon(ic) = alb_vis(ic)
          alb_nir_noon(ic) = alb_nir(ic)
          ! calculate albedo for total shortwave radiation
          !
          IF ((swvis_srf_down(ic) + swnir_srf_down(ic)) > eps8) THEN
            albedo_noon(ic) = ( alb_vis(ic) * swvis_srf_down(ic) + alb_nir(ic) * swnir_srf_down(ic) ) &
              &               / (swvis_srf_down(ic) + swnir_srf_down(ic))
          ELSE
            ! albedo is not defined when radiation <= 0 and set to zero
            albedo_noon(ic) = 0.0_wp
          END IF
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END IF
#endif

    !$ACC EXIT DATA DELETE(l_day, temperature_surface) ASYNC(acc_stream)
    DEALLOCATE(l_day, temperature_surface)

    IF (use_alb_soil_scheme .AND. (tile%is_bare .OR. tile%is_vegetation)) THEN
      !$ACC EXIT DATA DELETE(litter_carbon, soil_carbon) ASYNC(acc_stream)
      DEALLOCATE(litter_carbon, soil_carbon)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_albedo

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "albedo"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_albedo(tile, options)

    USE mo_jsb_time,          ONLY: is_time_experiment_start, is_time_ltrig_rad_m1
    USE mo_jsb_control,       ONLY: jsbach_runs_standalone
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(RAD_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    TYPE(t_jsb_model),       POINTER :: model
    INTEGER                          :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_albedo'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Get_memory(RAD_)

    CALL weighted_by_fract%BeginAggregate()

    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)

      ! equivalent to jsbach
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_pond,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_pond,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, snow_age,            weighted_by_fract)
      IF (.NOT. tile%is_vegetation) THEN
        dsl4jsb_Aggregate_onChunk(RAD_, albedo_lwtr,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(RAD_, albedo_lice,         weighted_by_fract)
      ENDIF

      ! quincy specific
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_can,            weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_can,            weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, albedo_noon,            weighted_by_fract)
#ifndef __NO_QUINCY__
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_noon,           weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_noon,           weighted_by_fract)
#endif

    CASE (MODEL_JSBACH)
#ifndef __QUINCY_STANDALONE__
    ! Update albedo only at experiment start, if model runs standalone or one time step before the radiation time step
    IF (.NOT. is_time_experiment_start(options%current_datetime) .AND. &
      & .NOT. is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id, jsbach_runs_standalone()) .AND. &
      & .NOT. jsbach_runs_standalone()) RETURN
#endif

      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_pond,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_pond,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, snow_age,            weighted_by_fract)
      IF (.NOT. tile%is_vegetation) THEN
        dsl4jsb_Aggregate_onChunk(RAD_, albedo_lwtr,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(RAD_, albedo_lice,         weighted_by_fract)
      ENDIF

    END SELECT

    CALL weighted_by_fract%EndAggregate()


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_albedo

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for RAD_
  !>
#ifndef __NO_QUINCY__
  SUBROUTINE update_time_average_radiation(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8, eps1
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_lnd_time_averages        ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    USE mo_lnd_bgcm_idx
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID
    dsl4jsb_Use_processes A2L_, VEG_, RAD_, Q_PHENO_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_canopy_q_assimi  !< Vertical grid
    TYPE(t_jsb_model),      POINTER         :: model
    TYPE(t_lnd_bgcm_store), POINTER         :: bgcm_store             !< the bgcm store of this tile
    LOGICAL,  DIMENSION(options%nc)         :: l_growing_season       !< growing_season LOGICAL
    LOGICAL,  ALLOCATABLE, DIMENSION(:,:)   :: l_growing_season_cl    !< growing_season LOGICAL, at vgrid_canopy_q_assimi
    INTEGER                                 :: icanopy                !< looping
    INTEGER                                 :: ncanopy                !< number of canopy layers, from vgrid
    INTEGER                                 :: iblk, ic, ics, ice, nc !< dimensions, looping
    REAL(wp)                                :: dtime                  !< timestep
    CHARACTER(len=*), PARAMETER             :: routine = modname//':update_time_average_radiation'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: daytime_counter
    dsl4jsb_Real2D_onChunk      :: local_time_day_seconds
    ! Q_PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! 3D RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_dacc_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_dacc_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(RAD_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)                 RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(Q_PHENO_)

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, daytime_counter)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)         ! in
    ! Q_PHENO_
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)             ! in
    ! 3D RAD_
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_daytime_dacc_cl)    ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_daytime_dacc_cl)    ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_daytime_cl)         ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_daytime_cl)         ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)                 ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)                 ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tcnl_mavg_cl)       ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tcnl_mavg_cl)       ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! inout
    ! ----------------------------------------------------------------------------------------------------- !

    ALLOCATE(l_growing_season_cl(nc, ncanopy))
    !$ACC ENTER DATA CREATE(l_growing_season_cl, l_growing_season) ASYNC(acc_stream)

    !>0.9 transform REAL growing_season(:) into LOGICAL l_growing_season(:)/_cl(:,:)
    !>
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      IF (growing_season(ic) > test_false_true) THEN
        l_growing_season(ic)      = .TRUE.
      ELSE
        l_growing_season(ic)      = .FALSE.
      END IF
    END DO
    !$ACC LOOP SEQ
    DO icanopy = 1, ncanopy
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        l_growing_season_cl(ic,icanopy) = l_growing_season(ic)
      END DO
    END DO
    !$ACC END PARALLEL

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 daytime averages
    !>
    !>

    !>  1.1 accumulate values - if light is available (i.e. daytime)
    !>
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO icanopy = 1, ncanopy
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (swpar_srf_down(ic) > eps8) THEN
        ! add values
          ppfd_sunlit_daytime_dacc_cl(ic,icanopy)  = ppfd_sunlit_daytime_dacc_cl(ic,icanopy) &
            & + ppfd_sunlit_cl(ic,icanopy)
          ppfd_shaded_daytime_dacc_cl(ic,icanopy)  = ppfd_shaded_daytime_dacc_cl(ic,icanopy) &
            & + ppfd_shaded_cl(ic,icanopy)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !>  1.2 calc daytime averages - at the timestep after local midnight (1st timestep of the new day)
    !>
    !>    calculate the average of the previous day and pass this value to the according variable
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO icanopy = 1, ncanopy
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (ABS(local_time_day_seconds(ic) - dtime) < eps1) THEN
          ! check if at least one value had been assigned at the current day, avoiding division by zero
          IF (daytime_counter(ic) > eps8) THEN
            ppfd_sunlit_daytime_cl(ic,icanopy) = ppfd_sunlit_daytime_dacc_cl(ic,icanopy) / daytime_counter(ic)
            ppfd_shaded_daytime_cl(ic,icanopy) = ppfd_shaded_daytime_dacc_cl(ic,icanopy) / daytime_counter(ic)
          ELSE
            ppfd_sunlit_daytime_cl(ic,icanopy) = 0.0_wp
            ppfd_shaded_daytime_cl(ic,icanopy) = 0.0_wp
          END IF
          ! zero the accumulation variables after daily average has been calculated
          ppfd_sunlit_daytime_dacc_cl(ic,icanopy) = 0.0_wp
          ppfd_shaded_daytime_dacc_cl(ic,icanopy) = 0.0_wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    ! ----------------------------------------------------------------------------------------------------- !
    !>2.0 moving averages
    !>
    !>

    ! docu:
    ! calc_time_mavg(dtime, current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO icanopy = 1, ncanopy
      DO ic = 1, nc
        !>  2.1 tfrac (averages at the time-scale of within-leaf N allocation fractions)
        !>
        ppfd_sunlit_tfrac_mavg_cl(ic,icanopy) = calc_time_mavg(dtime, ppfd_sunlit_tfrac_mavg_cl(ic,icanopy), &
          & ppfd_sunlit_daytime_cl(ic,icanopy), mavg_period_tfrac, do_calc=l_growing_season_cl(ic,icanopy))
        ppfd_shaded_tfrac_mavg_cl(ic,icanopy) = calc_time_mavg(dtime, ppfd_shaded_tfrac_mavg_cl(ic,icanopy), &
          & ppfd_shaded_daytime_cl(ic,icanopy), mavg_period_tfrac, do_calc=l_growing_season_cl(ic,icanopy))
        !>  2.2 tcnl (averages at the time-scale of leaf N allocation fractions)
        !>
        ppfd_sunlit_tcnl_mavg_cl(ic,icanopy) = calc_time_mavg(dtime, ppfd_sunlit_tcnl_mavg_cl(ic,icanopy), &
          & ppfd_sunlit_cl(ic,icanopy), mavg_period_tcnl, do_calc=l_growing_season_cl(ic,icanopy))
        ppfd_shaded_tcnl_mavg_cl(ic,icanopy) = calc_time_mavg(dtime, ppfd_shaded_tcnl_mavg_cl(ic,icanopy), &
          & ppfd_shaded_cl(ic,icanopy), mavg_period_tcnl, do_calc=l_growing_season_cl(ic,icanopy))
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(acc_stream)

    !$ACC EXIT DATA DELETE(l_growing_season_cl, l_growing_season) ASYNC(acc_stream)
    DEALLOCATE(l_growing_season_cl)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_radiation

#else

  SUBROUTINE update_time_average_radiation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

  END SUBROUTINE update_time_average_radiation

#endif
  ! ======================================================================================================= !
  !> aggregate: time moving averages and daytime averages for RAD_
  !>
  SUBROUTINE aggregate_time_average_radiation(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_radiation'
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
    dsl4jsb_Get_memory(RAD_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_time_average_radiation

#endif
END MODULE mo_rad_interface
