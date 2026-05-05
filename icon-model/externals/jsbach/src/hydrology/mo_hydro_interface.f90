!> Contains the interfaces to the hydro processes
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
!NEC$ options "-finline-file=externals/jsbach/src/hydrology/mo_hydro_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_hydro_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on, acc_stream
  USE mo_jsb_time,           ONLY: is_time_experiment_start
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_impl_constants, ONLY: WB_IGNORE

  ! Use of processes in this module
#ifndef __QUINCY_STANDALONE__
  dsl4jsb_Use_processes PHENO_
#endif
  dsl4jsb_Use_processes SEB_, SSE_, TURB_, A2L_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(TURB_)
#ifndef __QUINCY_STANDALONE__
  dsl4jsb_Use_memory(PHENO_)
#endif

#ifndef __NO_QUINCY__
  dsl4jsb_Use_processes VEG_, Q_ASSIMI_
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(Q_ASSIMI_)
#endif

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks
#ifndef __QUINCY_STANDALONE__
  PUBLIC :: global_hydrology_diagnostics
#endif
#ifdef __QUINCY_STANDALONE__
  ! USEd by the mo_qs_model_interface
  PUBLIC :: update_surface_hydrology, update_soil_properties, update_soil_hydrology, update_evaporation, &
    &       update_canopy_cond_unstressed, update_water_stress, update_canopy_cond_stressed, &
    &       update_snow_and_wet_fraction, update_water_balance
  ! Currently done in SEB process
  ! PUBLIC :: update_snow_and_ice_hydrology
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

  !> Type definition for surface_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_hydrology
  END TYPE tsk_surface_hydrology

  !> Constructor interface for surface_hydrology
  INTERFACE tsk_surface_hydrology
    PROCEDURE Create_task_surface_hydrology
  END INTERFACE tsk_surface_hydrology

  !> Type definition for soil_properties
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_properties
  END TYPE tsk_soil_properties

  !> Constructor interface for soil_properties
  INTERFACE tsk_soil_properties
    PROCEDURE Create_task_soil_properties
  END INTERFACE tsk_soil_properties

  !> Type definition for soil_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_hydrology
  END TYPE tsk_soil_hydrology

  !> Constructor interface for soil_hydrology
  INTERFACE tsk_soil_hydrology
    PROCEDURE Create_task_soil_hydrology
  END INTERFACE tsk_soil_hydrology

  !> Type definition for evaporation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_evaporation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_evaporation
    PROCEDURE, NOPASS :: Aggregate => aggregate_evaporation
  END TYPE tsk_evaporation

  !> Constructor interface for evaporation task
  INTERFACE tsk_evaporation
    PROCEDURE Create_task_evaporation
  END INTERFACE tsk_evaporation

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for water_stress task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_stress
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_stress
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_stress
  END TYPE tsk_water_stress

  !> Constructor interface for water_stress task
  INTERFACE tsk_water_stress
    PROCEDURE Create_task_water_stress
  END INTERFACE tsk_water_stress

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for snow_and_ice_hydrology task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_ice_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_ice_hydrology     ! Currently done in SEB process
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_ice_hydrology
  END TYPE tsk_snow_and_ice_hydrology

  !> Constructor interface for snow_and_ice_hydrology task
  INTERFACE tsk_snow_and_ice_hydrology
    PROCEDURE Create_task_snow_and_ice_hydrology
  END INTERFACE tsk_snow_and_ice_hydrology

  !> Type definition for snow_and_wet_fraction task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_wet_fraction
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_wet_fraction    !< Advances task computation for one time step
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_wet_fraction !< Aggregates computed task variables
  END TYPE tsk_snow_and_wet_fraction

  !> Constructor interface for snow_and_wet_fraction task
  INTERFACE tsk_snow_and_wet_fraction
    PROCEDURE Create_task_snow_and_wet_fraction                       !< Constructor function for task
  END INTERFACE tsk_snow_and_wet_fraction

  !> Type definition for water_balance task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_balance
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_balance
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_balance
  END TYPE tsk_water_balance

  !> Constructor interface for water balance task
  INTERFACE tsk_water_balance
    PROCEDURE Create_task_water_balance
  END INTERFACE tsk_water_balance

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_hydrology task
  !>
  FUNCTION Create_task_surface_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "surface_hydrology"

    ALLOCATE(tsk_surface_hydrology::return_ptr)
    CALL return_ptr%Construct(name='surface_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_properties task
  !>
  FUNCTION Create_task_soil_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "soil_properties"

    ALLOCATE(tsk_soil_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_properties', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_properties

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_hydrology task
  !>
  FUNCTION Create_task_soil_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "soil_hydrology"

    ALLOCATE(tsk_soil_hydrology::return_ptr)
    CALL return_ptr%Construct(name='soil_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_unstressed task
  !>
  FUNCTION Create_task_canopy_cond_unstressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "canopy_cond_unstressed"

    ALLOCATE(tsk_canopy_cond_unstressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_unstressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_stress task
  !>
  FUNCTION Create_task_water_stress(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "water_stress"

    ALLOCATE(tsk_water_stress::return_ptr)
    CALL return_ptr%Construct(name='water_stress', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_stress

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_stressed task
  !>
  FUNCTION Create_task_canopy_cond_stressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "canopy_cond_stressed"

    ALLOCATE(tsk_canopy_cond_stressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_stressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for evaporation task
  !>
  FUNCTION Create_task_evaporation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "evaporation"

    ALLOCATE(tsk_evaporation::return_ptr)
    CALL return_ptr%Construct(name='evaporation', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_evaporation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow_and_ice_hydrology task
  !>
  FUNCTION Create_task_snow_and_ice_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "snow_and_ice_hydrology"

    ALLOCATE(tsk_snow_and_ice_hydrology::return_ptr)
    ! TODO: Rename process (snow_AND_ice_) or rename all function names and comments
    CALL return_ptr%Construct(name='snow_ice_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_ice_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow_and_wet_fraction task
  !>
  FUNCTION Create_task_snow_and_wet_fraction(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "snow_and_wet_fraction"

    ALLOCATE(tsk_snow_and_wet_fraction::return_ptr)
    CALL return_ptr%Construct(name='snow_and_wet_fraction', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_wet_fraction

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_balance task
  !>
  FUNCTION Create_task_water_balance(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "water_balance"

    ALLOCATE(tsk_water_balance::return_ptr)
    CALL return_ptr%Construct(name='water_balance', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_balance

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for hydrology process
  !>
  SUBROUTINE Register_hydro_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process         !< Process
    INTEGER,              INTENT(in)    :: model_id        !< Model ID

    CALL process%Register_task( tsk_surface_hydrology      (model_id))
    CALL process%Register_task( tsk_soil_properties        (model_id))
    CALL process%Register_task( tsk_soil_hydrology         (model_id))
    CALL process%Register_task( tsk_water_stress           (model_id))
    CALL process%Register_task( tsk_canopy_cond_unstressed (model_id))
    CALL process%Register_task( tsk_canopy_cond_stressed   (model_id))
    CALL process%Register_task( tsk_evaporation            (model_id))
    !CALL process%Register_task( tsk_snow_and_ice_hydrology (model_id))
    CALL process%Register_task( tsk_snow_and_wet_fraction  (model_id))
    CALL process%Register_task( tsk_water_balance          (model_id))

  END SUBROUTINE Register_hydro_tasks

  ! ======================================================================================================== !
  !>
  !>#### Update surface hydrology
  !>
  !> In this subroutine we handle water fluxes in all phases at the surface - on the land (soil, canopy
  !> and ponds), lake and glacier tiles - and calculate the amount of water available for infiltration into
  !> the soil.

  SUBROUTINE update_surface_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: get_soilhyd_properties, calc_surface_hydrology_land, &
      &                                  calc_surface_hydrology_glacier
    USE mo_jsb_physical_constants, ONLY: dens_snow, grav
    USE mo_hydro_util,             ONLY: get_amount_in_rootzone

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile          !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options       !< Options of the current block

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(HYDRO_)     !< Configurable hydrology parameters
    dsl4jsb_Def_config(SSE_)       !< Configurable parameters of the soil and snow energy process
    dsl4jsb_Def_memory(HYDRO_)     !< Memory of the hydrology process
    dsl4jsb_Def_memory(SSE_)       !< Memory of the soil and snow energy process
    dsl4jsb_Def_memory(SEB_)       !< Memory of process surface energy balance
    dsl4jsb_Def_memory(A2L_)       !< Memory of the atmosphere to land interface
#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_memory(PHENO_)     !< Memory of the phenology process
#endif
#ifndef __NO_QUINCY__
    dsl4jsb_Def_memory(VEG_)       !< Memory of the vegetation process
#endif

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: q_snocpymlt        !< Heat required to melt snow on canopy [W m-2]
    dsl4jsb_Real2D_onChunk :: t_air              !< Air temperature (lowest level or 2m) [K]
    dsl4jsb_Real2D_onChunk :: wind_10m           !< Near surface wind speed [m s-1]
    dsl4jsb_Real2D_onChunk :: rain               !< Rain fall rate [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: snow               !< Snow fall rate [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: le_phase_change    !< Latent energy used for/gained from phase change [J m-2]
    dsl4jsb_Real2D_onChunk :: le_pc_remain       !< Latent energy available for phase change [J m-2]
    dsl4jsb_Real2D_onChunk :: wtr_skin           !< Amount of water in skin reservoir (soil and canopy) [m]
    dsl4jsb_Real2D_onChunk :: fract_skin         !< Wet skin fraction (soil and canopy, excluding ponds) []
    dsl4jsb_Real2D_onChunk :: weq_snow_soil      !< Amount of snow on the ground [m water equivalent]
    dsl4jsb_Real2D_onChunk :: snow_soil_dens     !< Snow density [kg m-3]
    dsl4jsb_Real2D_onChunk :: steepness          !< Parameter representing subgrid scale slopes []
    dsl4jsb_Real2D_onChunk :: fract_pond         !< Actual pond fraction []
    dsl4jsb_Real2D_onChunk :: fract_pond_max     !< Maximum pond fraction []
    dsl4jsb_Real2D_onChunk :: weq_pond_max       !< Maximum Water/ice in pond reservoir [m water equivalent]
    dsl4jsb_Real2D_onChunk :: wtr_pond           !< Current amount of liquid water in ponds [m]
    dsl4jsb_Real2D_onChunk :: ice_pond           !< Current amount of ice in ponds [m water equivalent]
    dsl4jsb_Real2D_onChunk :: weq_pond           !< Current amount of water or ice in ponds [m water equiv.]
    dsl4jsb_Real2D_onChunk :: pond_freeze        !< Amount of pond water freezing [m water equivalent]
    dsl4jsb_Real2D_onChunk :: pond_melt          !< Amount of pond ice melting [m water equivalent]
    dsl4jsb_Real2D_onChunk :: wtr_pond_net_flx   !< Net water inflow into ponds [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: weq_snow           !< Amount of snow [m water equivalent]
    dsl4jsb_Real2D_onChunk :: snow_accum         !< Snow budget change within time step [m water equivalent]
    dsl4jsb_Real2D_onChunk :: fract_snow         !< Snow cover fraction
    dsl4jsb_Real2D_onChunk :: weq_snow_can       !< Amount of snow on canopy [m water equivalent]
    dsl4jsb_Real2D_onChunk :: evapotrans_soil    !< Evapotranspiration from soil [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: evapo_skin         !< Evaporation from skin reservoir [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: evapo_deficit      !< Evaporation deficit due to inconsistencies
                                                 !< [m water eq./(time step)]
    dsl4jsb_Real2D_onChunk :: lai                !< LAI, leaf area index referring to the canopy area []
    dsl4jsb_Real2D_onChunk :: fract_fpc_max      !< Maximum foliage projected cover fraction []
    dsl4jsb_Real2D_onChunk :: evapotrans         !< Evapotranspiration incl. sublimation [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: transpiration      !< Transpiration [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: evapopot           !< Potential evapotranspiration (i.e without water limitation)
                                                 !< [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: evapo_snow         !< Sublimation from snow [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: evapo_pond         !< Evaporation/sublimation from ponds [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: snowmelt           !< Snow/ice melt flux [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: infilt             !< Infiltration of water into the soil [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: runoff             !< Surface runoff [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: runoff_horton      !< Horton component of surface runoff [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: drainage           !< Drainage [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: runoff_glac        !< Runoff from glaciers [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: water_to_soil      !< Water available for infiltration into the soil [m]
    dsl4jsb_Real2D_onChunk :: weq_glac           !< Glacier depth [m water equivalent]

    dsl4jsb_Real3D_onChunk :: wtr_soil_sl           !< Amount of soil water [m]
    dsl4jsb_Real3D_onChunk :: ice_soil_sl           !< Amount of soil ice [m water equivalent]
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl !< Potential amount of supercooled water [m]
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl       !< Volumetric soil porosity [m/m]
    dsl4jsb_Real3D_onChunk :: vol_wres_sl           !< Volumetric residual water content [m/m]
    dsl4jsb_Real3D_onChunk :: hyd_cond_sat_sl       !< Hydraulic conductivity of soil layers [m s-1]
    dsl4jsb_Real3D_onChunk :: matric_pot_sl         !< Matric potential of soil layers at saturation [m]
    dsl4jsb_Real3D_onChunk :: bclapp_sl             !< Clapp and Hornberger exponent b of soil layers []
    dsl4jsb_Real3D_onChunk :: pore_size_index_sl    !< Pore size index of soil layers []
    dsl4jsb_Real3D_onChunk :: soil_depth_sl         !< Soil depth until bedrock within each layer [m]
    dsl4jsb_Real3D_onChunk :: t_soil_sl             !< Soil temperature of the layer [K]
    dsl4jsb_Real3D_onChunk :: t_snow                !< Snow temperature (on snow layers) [K]
    dsl4jsb_Real3D_onChunk :: snow_depth_sl         !< Snow depth within the snow layers [m]
    ! quincy
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_sl       !< Soil water potential per layer [MPa]

    dsl4jsb_Real2D_onChunk :: wtr_latflow_res_srf   !< Intermediate storage for surface runoff [m]
    dsl4jsb_Real2D_onChunk :: wtr_latflow_srf       !< Outflow from intermediate surface runoff storage
                                                    !< [kg m-2 s-1]

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & skinres_max,               & !< Capacity of the skin reservoir (soil and canopy) [m]
      & skinres_canopy_max,        & !< Capacity of the canopy skin reservoir [m]
      & wpi_rootzone,              & !< Amount of water and ice in the root zone [m]
      & wpi_rootzone_max,          & !< Maximum amount of water and ice in root zone [m]
      & wtr_soil,                  & !< Amount of liquid soil water [m]
      & t_snow_mean,               & !< Level weighted mean snow temperature [K]
      & total_snow_depth,          & !< Snow depth [m]
      & trans_tmp,                 & !< Transpiration [kg m-2 s-1]
      & flowlag_used,              & !< Lag factor to account for retention in surface runoff []
      & slope_used                   !< Effective slope on tile or steepness for ARNO scheme []

    REAL(wp), DIMENSION(options%nc, options%nsoil_w) :: &
      & wpi_soil_sl_tmp,           & !< Actual amount of water and ice in soil layers [m]
      & wpi_wsat_soil_sl_tmp,      & !< Maximum amount of water and ice in soil layers [m]
      & wpi_wres_soil_sl_tmp,      & !< Residual water/ice content of soil layer [m]
      & ice_impedance_sl,          & !< Ice impedance factor for soil layers []
      & mpot_act_sl                  !< Soil matric potential of soil layer [m]

    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice, &        !< Index of last cell of block
      & nc, &         !< Number of cells in current block
      & nsoil, &      !< Number of soil layers
      & nsnow, &      !< Number of snow layers
      & ic, &         !< Cell index
      & is, &         !< Soil layer index
      & isnow         !< Snow layer index

    REAL(wp) :: dtime !< Time step length

    ! Variables from process configuration
    REAL(wp) :: &
      & config_w_skin_max,       & !< Maximum amount of snow on canopy and liquid water in skin
                                   !< reservoirs (canopy and soil); snow on soil is not included
                                   !< [m water equivalent]
      & snow_depth_max_config,   & !< Maximum snow depth [m]
      & ret_macro_srf              !< Retention time for surface runoff, i.e. time until water
                                   !< reaches tributary [s]
    INTEGER  :: &
      & hydro_scale                !< Assumed scale for HYDRO process
    LOGICAL  :: &
      & is_experiment_start,     & !< True: This is the first time step of an experiment
      & l_dynsnow_config,        & !< True: Snow density calculated dynamically
      & l_infil_subzero_config,  & !< True: Infiltration is possible below zero degree Celsius
      & l_latflow_to_streamflow, & !< With "[[t_hydro_config:hydro_scale]]=Uniform"
                                   !< True: outflow of intermediary storage goes to runoff/drainage
                                   !< False: outflow of intermediary storage goes to downstream tile
      & ltpe_closed,             & !< Terraplanet setup with closed water balance
      & use_tmx                    !< Use turbulent mixing scheme tmx

    TYPE(t_jsb_model), POINTER :: model   !< Current instance of the model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    nsoil = options%nsoil_w
    nsnow = options%nsnow_e
    dtime = options%dtime
    is_experiment_start = is_time_experiment_start(options%current_datetime)

    ! Return, if hydrology does not need to be calculated on this tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    use_tmx = model%config%use_tmx

    ! Define local variables from HYDRO namelist parameters (needed on GPUs)
    dsl4jsb_Get_config(HYDRO_)
    config_w_skin_max      = dsl4jsb_Config(HYDRO_)%w_skin_max
    snow_depth_max_config  = dsl4jsb_Config(HYDRO_)%snow_depth_max
    l_infil_subzero_config = dsl4jsb_Config(HYDRO_)%l_infil_subzero
    hydro_scale            = dsl4jsb_Config(HYDRO_)%hydro_scale

    ! Define local variable from SSE namelist parameter (needed on GPUs)
    dsl4jsb_Get_config(SSE_)
    l_dynsnow_config = dsl4jsb_Config(SSE_)%l_dynsnow

    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    ! Only used with terraplanet setup
    ltpe_closed = .FALSE.
    IF (model%config%tpe_scheme == 'closed') ltpe_closed = .TRUE.

    IF (tile%is_lake) THEN
      !> On the lake tile we just compute runoff from P-E, which is needed for the water balance.
      !> There is no infiltration nor drainage.

      ! TODO: compute actual water balance (snow/ice melt) in lake model and use for runoff. The corresponding energy
      ! fluxes and snow/ice budgets are already computed in the lake model.
      dsl4jsb_Get_var2D_onChunk(A2L_,   rain)               ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,   snow)               ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)           ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, infilt)             ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)             ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)           ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)           ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_snow)         ! out
      IF (use_tmx) dsl4jsb_Get_var2D_onChunk(HYDRO_, q_snocpymlt)    ! out

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(acc_stream) GANG VECTOR
      DO ic=1,nc
        runoff  (ic) = rain(ic) + snow(ic) + evapotrans(ic)
        infilt(ic)   = 0._wp
        drainage(ic) = 0._wp
        weq_snow(ic) = 0._wp            ! Later updated in lake energy balance
        evapo_snow(ic) = evapopot(ic)   ! TODO: check why this is set to evapopot
      END DO
      !$ACC END PARALLEL LOOP
      IF (use_tmx) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(acc_stream) GANG VECTOR
        DO ic=1,nc
          q_snocpymlt(ic) = 0._wp       ! no canopy on lake ice --> no heat release from canopy snow
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Nothing more to do on lake tile
      RETURN
    END IF

    ! Set pointers to variables needed on non-lake tiles
    dsl4jsb_Get_memory(SSE_)
    IF (tile%contains_vegetation) THEN
      SELECT CASE (model%config%model_scheme)
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

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   wind_10m)       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   rain)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)           ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   le_phase_change)! in

    dsl4jsb_Get_var2D_onChunk(HYDRO_, le_pc_remain)   ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, steepness)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens) ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_snow)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, q_snocpymlt)    ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, infilt)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff_horton)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)       ! out

    IF (tile%is_glacier) THEN
      ! Variables only needed on glacier tile
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_glac)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   runoff_glac)      ! inout
    ELSE
      ! Variables only needed on tiles with soil
      dsl4jsb_Get_var3D_onChunk(SSE_,     t_soil_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(SSE_,     t_snow)           ! in
      dsl4jsb_Get_var3D_onChunk(SSE_,     snow_depth_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   wtr_soil_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   ice_soil_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   wtr_soil_pot_scool_sl) ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   vol_porosity_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   vol_wres_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   hyd_cond_sat_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   matric_pot_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   bclapp_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   pore_size_index_sl) ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   soil_depth_sl)    ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   fract_skin)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_skin)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   water_to_soil)    ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   snow_accum)       ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapotrans_soil)  ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_skin)       ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_deficit)    ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   fract_pond)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   fract_pond_max)   ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_pond_max)     ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_pond)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   ice_pond)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_latflow_res_srf) ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_pond)         ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_pond)       ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   pond_freeze)      ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   pond_melt)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_pond_net_flx) ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_latflow_srf)  ! out
#ifndef __NO_QUINCY__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        dsl4jsb_Get_var3D_onChunk(HYDRO_,   wtr_soil_pot_sl)  ! out  [MPa]
      END SELECT
#endif

      IF (tile%contains_vegetation) THEN
        ! Variables only needed on tiles with vegetation
        dsl4jsb_Get_var2D_onChunk(HYDRO_,   transpiration)    ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_snow_can)     ! in
        SELECT CASE (model%config%model_scheme)
#ifndef __QUINCY_STANDALONE__
        CASE (MODEL_JSBACH)
          dsl4jsb_Get_var2D_onChunk(PHENO_,   lai)              ! in
          dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max)    ! in
#endif
#ifndef __NO_QUINCY__
        CASE (MODEL_QUINCY)
          dsl4jsb_Get_var2D_onChunk(VEG_,     lai)              ! in
#endif
        END SELECT
      END IF
    END IF

    ! Set to true, until HydroTiles are implemented
    l_latflow_to_streamflow = .TRUE.

    ! TODO: Is this needed?
    IF (is_experiment_start) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        weq_snow(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Switch sign of phase change flux for the hydrology related processes. From the
    !   energy balance perspective, ice melting is an energy sink (therefore negative)
    !   and freezing a source (therefore positive). However, from the hydrology perspective
    !   melting occurs when enough energy is available (positive sign) and freezing if
    !   energy is required (therefore negative).
    !   Thus, we switch the sign now to make the hydrology equations more intuitive.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      le_pc_remain(ic) = -1._wp * le_phase_change(ic)
    END DO
    !$ACC END PARALLEL LOOP

    IF (tile%is_glacier) THEN
      !> On the glacier tile runoff is calculated from rain and snow melt. As on the lake tile,
      !> there is no infiltration or drainage.

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
#ifdef NFORT_BROKEN_INLINES
        !NEC$ noinline
#endif
        CALL calc_surface_hydrology_glacier( &
          & is_experiment_start,             & ! in
          & dtime,                           & ! in
          & fract_snow(ic),                  & ! in
          & evapotrans(ic),                  & ! in
          & evapopot(ic),                    & ! in
          & rain(ic),                        & ! in
          & snow(ic),                        & ! in
          & weq_glac(ic),                    & ! inout
          & le_pc_remain(ic),                & ! inout
          & q_snocpymlt(ic),                 & ! out
          & snowmelt(ic),                    & ! out
          & runoff_glac(ic),                 & ! out
          & pme_glacier = runoff(ic)         & ! out, P-E ... used as runoff for HD model
          & )
        infilt(ic)   = 0._wp
        drainage(ic) = 0._wp
        weq_snow_soil(ic) = 0._wp
        evapo_snow(ic) = evapopot(ic)
        snow_soil_dens(ic) = dens_snow
      END DO
      !$ACC END PARALLEL LOOP
    ! NOT tile%is_glacier
    ELSE
      !> Things are getting more complicated on the land tile: after updating the soil and
      !> snow properties, hydrological fluxes and the hydrological state at the land surface
      !> are calculated in [[calc_surface_hydrology_land]].

      !$ACC DATA CREATE(skinres_max, skinres_canopy_max, wpi_rootzone, wpi_rootzone_max) &
      !$ACC   CREATE(trans_tmp, wpi_soil_sl_tmp(1:nc,1:nsoil), wtr_soil(1:nc)) &
      !$ACC   CREATE(wpi_wsat_soil_sl_tmp(1:nc,1:nsoil), wpi_wres_soil_sl_tmp(1:nc,1:nsoil)) &
      !$ACC   CREATE(ice_impedance_sl(1:nc,1:nsoil), t_snow_mean(1:nc), total_snow_depth(1:nc)) &
      !$ACC   CREATE(mpot_act_sl(1:nc,1:nsoil)) &
      !$ACC   CREATE(slope_used(1:nc), flowlag_used(1:nc)) ASYNC(acc_stream)

      !TODO: Move up to the other config parameter
      ret_macro_srf = dsl4jsb_Config(HYDRO_)%ret_macro_srf
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        slope_used(ic)   = steepness(ic)
        flowlag_used(ic) = dtime / ret_macro_srf
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate absolute layer values from the volumetric variables
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          wpi_wsat_soil_sl_tmp(ic,is)      = vol_porosity_sl(ic,is)  * soil_depth_sl(ic,is)
          wpi_wres_soil_sl_tmp(ic,is)      = vol_wres_sl(ic,is)      * soil_depth_sl(ic,is)
          ! @todo There should be no water in cells with field_cap == 0. Still, this is currently
          !       necessary. WHY?
          IF (wpi_wsat_soil_sl_tmp(ic,is) > 0._wp) THEN
            wpi_soil_sl_tmp(ic,is) = wtr_soil_sl(ic,is) + ice_soil_sl(ic,is)
          ELSE
            wpi_soil_sl_tmp(ic,is) = 0._wp
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        wtr_soil(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          wtr_soil(ic) = wtr_soil(ic) + wtr_soil_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL

      ! Compute actual and maximum total root zone moisture (liquid + ice)
      ! Note: here wpi_rootzone_max is computed based on water saturation (i.e. soil porosity), not on field
      ! capacity. It is used to compute the bucket overflow in the ARNO scheme for infiltration and runoff,
      ! and therefore should consider the maximum amount of water that can fit into the root zone.
      CALL get_amount_in_rootzone(wpi_soil_sl_tmp(:,:),                                                  &
        &  dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl), dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl), &
        &  wpi_rootzone(:))
      CALL get_amount_in_rootzone(wpi_wsat_soil_sl_tmp(:,:),                                             &
        &  dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl), dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl), &
        &  wpi_rootzone_max(:))

      ! Compute the capacity of the canopy skin reservoir, i.e. the maximum amount of water/snow that
      ! can be held by the canopy.
      IF (tile%contains_vegetation) THEN
        SELECT CASE (model%config%model_scheme)
#ifndef __QUINCY_STANDALONE__
        CASE (MODEL_JSBACH)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1,nc
            skinres_canopy_max(ic) = config_w_skin_max * lai(ic) * fract_fpc_max(ic)
          END DO
          !$ACC END PARALLEL LOOP
#endif
#ifndef __NO_QUINCY__
        CASE (MODEL_QUINCY)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1,nc
            ! no scaling with 'fract_fpc(:)' because:
            ! max LAI fraction of the tile area is one (but not fract_fpc(:)) as vegetation spreads across the whole tile area
            skinres_canopy_max(ic) = config_w_skin_max * lai(ic)
          END DO
          !$ACC END PARALLEL LOOP
#endif
        END SELECT
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1,nc
          skinres_canopy_max(ic) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! The total skin reservoir capacity consists of the skin reservoirs of the ground and of the canopy.
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1,nc
        IF (tile%contains_soil) THEN
          skinres_max(ic) =  config_w_skin_max + skinres_canopy_max(ic)
        ELSE
          skinres_max(ic) = 0._wp
        END IF
        ! Note: Transpiration is not defined on tiles without vegetation, we thus need this local variable.
        IF (tile%is_vegetation) THEN
          trans_tmp(ic) = transpiration(ic)
        ELSE
          trans_tmp(ic) = 0._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate weighted mean snow temperature from the temperatures of each snow layer

      ! Note: For the snow density calculation in calc_surface_hydrology_land we need a temperature for
      !       all grid cells with weq_snow_soil > 0. There is however a mismatch between weq_snow_soil
      !       (calculated in the routine) and snow_depth_sl used below, which was calculated earlier in
      !       SSE:update_soil_and_snow_temperature based on weq_snow_soil of the previous time step.
      !       To also provide a temperature to snow in grid cells that did not have snow in the time
      !       step before, we here initialize t_snow_mean for all grid cells.
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        total_snow_depth(ic) = SUM(snow_depth_sl(ic, :))  ! Preparation for the below loop
        t_snow_mean(ic)      = t_snow(ic,1)               ! Surface temperature, in case uppermost layer has no snow.
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      ! Note: Layer nsnow is always the lowest snow layer. If only one snow layer is needed, this is
      !       layer nsnow. Snow layer 1 is only needed if snow layers nsnow to 2 are filled.
      DO isnow = nsnow, 1, -1
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          IF (snow_depth_sl(ic, isnow) > EPSILON(1.0_wp)) THEN
            IF (isnow == nsnow) THEN
              ! Initialize with weighted average of snow layer temperature for the lowest snow layer
              t_snow_mean(ic) = t_snow(ic, isnow) * (snow_depth_sl(ic, isnow) / total_snow_depth(ic))
            ELSE
              ! Add weighted average of snow layer temperature using layer depth fractions
              t_snow_mean(ic) = t_snow_mean(ic) + t_snow(ic, isnow) * (snow_depth_sl(ic, isnow) / total_snow_depth(ic))
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

      ! Compute soil hydrological properties to get top layer ice impedance and
      ! soil hydrological conductivity for point scale (i.e. uniform) infiltration
      CALL get_soilhyd_properties(                                                           &
        & dsl4jsb_Config(HYDRO_)%soilhydmodel,                                               &
        & dsl4jsb_Config(HYDRO_)%interpol_mean,                                              &
        & nc, nsoil, soil_depth_sl(:,:),                                                     &
        & wtr_soil_sl(:,:), ice_soil_sl(:,:), wtr_soil_pot_scool_sl(:,:),                    &
        & wpi_wsat_soil_sl_tmp(:,:), wpi_wres_soil_sl_tmp(:,:),                              &
        & hyd_cond_sat_sl(:,:), matric_pot_sl(:,:), bclapp_sl(:,:), pore_size_index_sl(:,:), &  ! in
        & ice_impedance=ice_impedance_sl(:,:)                                                &  ! out
        )

      ! Update hydrological fluxes and pools at the land surface
      CALL calc_surface_hydrology_land(   &
        ! in
        & is_experiment_start,            & ! First time step of the experiment?
        & dtime,                          & ! Time step length [s]
        & ltpe_closed,                    & ! With Terraplanet: closed water balance?
        & hydro_scale,                    & ! Hydrological scale scheme: uniform or semi_distributed
        & l_dynsnow_config,               & ! Dynamic snow density?
        & l_infil_subzero_config,         & ! Infiltration at temperature below 0 C?
        & l_latflow_to_streamflow,        & ! Pass runoff directly to stream flow?
        & snow_depth_max_config,          & ! Maximum allowed snow depth [m]
        & slope_used(:),                  & ! effective slope of tile / steepness []
        & flowlag_used(:),                & ! Factor to account for retention in surface runoff []
        & t_soil_sl(:,1),                 & ! Temperature of the uppermost soil layer [K]
        & t_snow_mean(:),                 & ! Mean Temperature of the snow layers [K]
        & wind_10m(:),                    & ! Wind speed at 10m height [m/s]
        & t_air(:),                       & ! Lowest layer atmosphere temperature [K]
        & skinres_canopy_max(:),          & ! Capacity of the canopy skin reservoir (also used to limit snow on canopy) [m]
        & skinres_max(:),                 & ! Total capacity of the skin reservoirs, i.e. soil and canopy [m]
        & weq_pond_max(:),                & ! Total capacity of pond reservoir (water + ice) [m]
        & fract_snow(:),                  & ! Snow cover fraction (not including canopy) []
        & fract_skin(:),                  & ! Wet skin fraction []
        & fract_pond(:),                  & ! Actual pond fraction []
        & fract_pond_max(:),              & ! Maximum pond fraction []
        & evapotrans(:),                  & ! Evapotranspiration (including sublimation) [kg m-2 s-1]
        & evapopot(:),                    & ! Potential evaporation (if there was enough water/ice) [kg m-2 s-1]
        & trans_tmp(:),                   & ! Transpiration [kg m-2 s-1]
        & rain(:),                        & ! Liquid precipitation [kg m-2 s-1]
        & snow(:),                        & ! Solid precipitation [kg m-2 s-1]
        & wpi_rootzone(:),                & ! Actual content of water plus ice in the root zone [m]
        & wpi_rootzone_max(:),            & ! Maximum possible content of water plus ice in the root zone [m]
        & wtr_soil(:),                    & ! Soil water content until bedrock (aggregated from soil layers) [m]
        & hyd_cond_sat_sl(:,1),           & ! Saturated hydraulic conductivity of the top soil layer [m s-1]
        & ice_impedance_sl(:,1),          & ! Ice impedance factor for the top soil layer []
        ! inout
        & wtr_skin(:),                    & ! Water content of the skin reservoir (canopy and bare soil) [m]
        & weq_snow_soil(:),               & ! Amount of snow on the ground [m water equivalent]
        & weq_snow_can(:),                & ! Amount of snow on the canopy [m water equivalent]
        & snow_soil_dens(:),              & ! Density of snow on the ground [kg m-3]
        & wtr_pond(:),                    & ! Water content of pond reservoir [m]
        & ice_pond(:),                    & ! Ice content of pond reservoir [m]
        & wtr_latflow_res_srf(:),         & ! Intermediate storage for surface runoff [m]
        & le_pc_remain(:),                & ! (Remaining) latent energy for phase change [J m-2]
        ! out
        & q_snocpymlt(:),                 & ! Heating due to snow melt on canopy [W m-2]
        & snow_accum(:),                  & ! Snow budget change within time step [m water equivalent]
        & snowmelt(:),                    & ! Snow/ice melt at land points (excluding canopy) [kg m-2 s-1]
        & pond_freeze(:),                 & ! Amount of water freezing within ponds [kg m-2 s-1]
        & pond_melt(:),                   & ! Amount of ice melting within ponds [kg m-2 s-1]
        & evapotrans_soil(:),             & ! Evapotranspiration from soil (excluding skin, snow and pond
                                            ! reservoirs) [kg m-2 s-1]
        & evapo_skin(:),                  & ! Evaporation from skin reservoir [kg m-2 s-1]
        & evapo_snow(:),                  & ! Evaporation/sublimation from snow [kg m-2 s-1]
        & evapo_pond(:),                  & ! Evaporation/sublimation from ponds [kg m-2 s-1]
        & wtr_pond_net_flx(:),            & ! Diagnostic net inflow into surface water ponds [kg m-2 s-1]
        & water_to_soil(:),               & ! Water available for infiltration into the soil [m]
        & evapo_deficit(:),               & ! Evaporation deficit from inconsistencies [m water eq./(time step)]
        & infilt(:),                      & ! Infiltration flux into the soil [m / (time step)]
        & runoff(:),                      & ! Surface runoff [m / (time step)]
        & runoff_horton(:),               & ! Horton component of surface runoff [m / (time step)]
        & wtr_latflow_srf(:)              & ! outflow from intermediate storage for surface runoff [kg m-2 s-1]
        & )

      ! Update the pond storage
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        weq_pond(ic) = wtr_pond(ic) + ice_pond(ic)
      END DO
      !$ACC END PARALLEL LOOP

#ifndef __NO_QUINCY__
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)

      ! Compute soil hydrological properties to get top layer ice impedance and
      ! soil hydrological conductivity for point scale infiltration
      !
      ! called here to get an updated value of mpot_act_sl(:,:)
      !
      CALL get_soilhyd_properties(                                                           &
        & dsl4jsb_Config(HYDRO_)%soilhydmodel,                                               &
        & dsl4jsb_Config(HYDRO_)%interpol_mean,                                              &
        & nc, nsoil, soil_depth_sl(:,:),                                                     &
        & wtr_soil_sl(:,:), ice_soil_sl(:,:), wtr_soil_pot_scool_sl(:,:),                    &
        & wpi_wsat_soil_sl_tmp(:,:), wpi_wres_soil_sl_tmp(:,:),                              &
        & hyd_cond_sat_sl(:,:), matric_pot_sl(:,:), bclapp_sl(:,:), pore_size_index_sl(:,:), &  ! in
        & mpot_act=mpot_act_sl(:,:)                                                          &  ! out
        )

      ! Update soil water potential based on "soil matric potential"
      !   includes unit conversion from: m -> J/N -> MPa
      ! TODO: replace the 1000 with a parameter (rhoh2o from mo_jsb_physical_constants?)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          wtr_soil_pot_sl(ic, is) = mpot_act_sl(ic, is) * grav / 1000.0_wp
        END DO
      END DO
      !$ACC END PARALLEL
    END SELECT
#endif

      !$ACC END DATA
    END IF ! NOT tile%is_glacier

    ! Calculate the total amount of snow from snow on ground and snow on canopy
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    IF (tile%is_glacier) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic=1,nc
          weq_snow(ic) = 0._wp
      END DO
      !$ACC END LOOP
    ELSE IF (tile%is_bare) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic=1,nc
        weq_snow(ic) = weq_snow_soil(ic)
      END DO
      !$ACC END LOOP
    ELSE IF (tile%contains_vegetation) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic=1,nc
        weq_snow(ic) = weq_snow_soil(ic) + weq_snow_can(ic)
      END DO
      !$ACC END LOOP
    END IF
    !$ACC END PARALLEL

  END SUBROUTINE update_surface_hydrology

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of variables from surface hydrology
  !>
  !> Surface hydrology variables are aggregated in this routine.
  !>
  ! -------------------------------------------------------------------------------------------------------
  SUBROUTINE aggregate_surface_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    dsl4jsb_Def_memory(HYDRO_)                             !< Memory of the hydrology process

    TYPE(t_jsb_model),       POINTER :: model              !< This instance of ICON-Land
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_hydrology'

    INTEGER :: &
      & iblk, &             !< Number of current block
      & ics, &              !< Index of first cell of current block
      & ice                 !< Index of last cell of current block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_skin,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_soil,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_can,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snowmelt,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, q_snocpymlt,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_glac,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_accum,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_soil_dens,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans_soil,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_skin,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_snow,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_pond,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_to_soil,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, infilt,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff_horton,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pond_freeze,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pond_melt,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond_net_flx,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_res_srf, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_srf,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, le_pc_remain,        weighted_by_fract)
#ifndef __NO_QUINCY__
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pot_sl,     weighted_by_fract)
    END SELECT
#endif

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_hydrology

  ! ============================================================================================== !
  !>
  !> Update routine for soil properties
  !>
  !> Soil properties are calculated for each soil layer, with [[t_hydro_config:l_organic]]=true
  !> taking into account organic soil fractions. As organic soil fractions might change with
  !> time - which is however not the case in current setups - the properties are updated on every
  !> time step.
  !>
  SUBROUTINE update_soil_properties(tile, options)

    USE mo_hydro_constants, ONLY: vol_porosity_org_top, hyd_cond_sat_org_top, bclapp_org_top, matric_pot_org_top, &
      & pore_size_index_org_top, vol_field_cap_org_top, vol_p_wilt_org_top, vol_wres_org_top, &
      & vol_porosity_org_below, hyd_cond_sat_org_below, bclapp_org_below, matric_pot_org_below, &
      & pore_size_index_org_below, vol_field_cap_org_below, vol_p_wilt_org_below, vol_wres_org_below, &
      & thresh_org, beta_perc


    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options      !< Options of current block

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)          !< Configurable parameters of the hydrology
    dsl4jsb_Def_memory(HYDRO_)          !< Memory of the hydrology

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & hyd_cond_sat,         &         !< Saturated hydraulic conductivity [m/s]
      & vol_porosity,         &         !< Volumetric soil porosity [m3/m3]
      & bclapp,               &         !< Clapp and Hornberger exponent b []
      & matric_pot,           &         !< Saturated soil matric potential [m]
      & pore_size_index,      &         !< Soil pore size index []
      & vol_field_cap,        &         !< Volumetric field capacity [m/m]
      & vol_p_wilt,           &         !< Volumetric moisture at permanent wilting point [m/m]
      & vol_wres                        !< Volumetric residual water content [m/m]
    dsl4jsb_Real3D_onChunk :: &
      & fract_org_sl,         &         !< Soil layer organic matter fraction []
      & hyd_cond_sat_sl,      &         !< Saturated hydraulic conductivity of soil layer [m/s]
      & vol_porosity_sl,      &         !< Volumetric soil layer porosity [m3/m3]
      & bclapp_sl,            &         !< Clapp and Hornberger exponent b in soil layer []
      & matric_pot_sl,        &         !< Saturated soil layer matric potential [m]
      & pore_size_index_sl,   &         !< Soil layer pore size index []
      & vol_field_cap_sl,     &         !< Volumetric field capacity of soil layer [m/m]
      & vol_p_wilt_sl,        &         !< Volumetric moisture at permanent wilting point in layer [m/m]
      & vol_wres_sl                     !< Volumetric residual water content in soil layer [m/m]

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model !< This instance of jsbach

    INTEGER  :: &
      & iblk,   &                       !< Index of the current block
      & ics,    &                       !< Index of first grid cell in block
      & ice,    &                       !< Index of last grid cell in block
      & nc,     &                       !< Number of grid cells in block
      & ic,     &                       !< Looping index for grid cells
      & nsoil,  &                       !< Number of soil layers (vertical dimension)
      & is                              !< Looping index for soil layers

    REAL(wp) :: N_perc                                 !< Variable used in percolation theory
    REAL(wp), ALLOCATABLE ::  fract_perc(:,:)          !< Fraction of the soil with connected organic pathways
    REAL(wp), ALLOCATABLE ::  fract_uncon(:,:)         !< Fraction of the soil with unconnected organic matter
    REAL(wp), ALLOCATABLE ::  hyd_cond_sat_uncon(:,:)  !< Saturated hydraulic conductivity of 'unconnected'
                                                       !< soil fraction
    REAL(wp), ALLOCATABLE ::  hyd_cond_sat_org(:,:)    !< Saturated hydraulic conductivity of organic matter [m/s]
    REAL(wp), ALLOCATABLE ::  vol_porosity_org(:,:)    !< Volumetric porosity of organic matter [m3/m3]
    REAL(wp), ALLOCATABLE ::  bclapp_org(:,:)          !< Clapp and Hornberger exponent B for organic matter []
    REAL(wp), ALLOCATABLE ::  matric_pot_org(:,:)      !< Soil matric potential for organic matter [m]
    REAL(wp), ALLOCATABLE ::  pore_size_index_org(:,:) !< Pore size index of organic matter []
    REAL(wp), ALLOCATABLE ::  vol_field_cap_org(:,:)   !< Volumetric field capacity of organic matter [m/m]
    REAL(wp), ALLOCATABLE ::  vol_p_wilt_org(:,:)      !< Volumetric soil moisture at permanent wilting point
                                                       !< in organic matter [m/m]
    REAL(wp), ALLOCATABLE ::  vol_wres_org(:,:)        !< Volum. residual water content in organic matter [m/m]

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    ! Return if hydrology calculations are not needed on this tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    ! Lake and glacier tiles have no soil.
    IF (tile%is_lake .OR. tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    nsoil = options%nsoil_w

    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(HYDRO_, hyd_cond_sat)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_porosity)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, bclapp)             ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, matric_pot)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, pore_size_index)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_field_cap)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_p_wilt)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_wres)           ! in
    IF (dsl4jsb_Config(HYDRO_)%l_organic)   &
      &  dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, hyd_cond_sat_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, matric_pot_sl)      ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, pore_size_index_sl) ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_p_wilt_sl)      ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_wres_sl)        ! out

    !> The soil properties that had been set as two dimensional maps in [[mo_hydro_init:hydro_init_bc]]
    !> are spread to all soil layers.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is=1,nsoil
      DO ic=1,nc
        hyd_cond_sat_sl   (ic,is) = hyd_cond_sat(ic)
        vol_porosity_sl   (ic,is) = vol_porosity(ic)
        bclapp_sl         (ic,is) = bclapp(ic)
        matric_pot_sl     (ic,is) = matric_pot(ic)
        pore_size_index_sl(ic,is) = pore_size_index(ic)
        vol_field_cap_sl  (ic,is) = vol_field_cap(ic)
        vol_p_wilt_sl     (ic,is) = vol_p_wilt(ic)
        vol_wres_sl       (ic,is) = vol_wres(ic)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN

      ! Allocation of corresponding maps for the properties of soil organic matter fractions
      ALLOCATE(hyd_cond_sat_org   (nc, nsoil))
      ALLOCATE(vol_porosity_org   (nc, nsoil))
      ALLOCATE(bclapp_org         (nc, nsoil))
      ALLOCATE(matric_pot_org     (nc, nsoil))
      ALLOCATE(pore_size_index_org(nc, nsoil))
      ALLOCATE(vol_field_cap_org  (nc, nsoil))
      ALLOCATE(vol_p_wilt_org     (nc, nsoil))
      ALLOCATE(vol_wres_org       (nc, nsoil))
      ! TODO: these are not needed as arrays
      ALLOCATE(fract_perc         (nc, nsoil))  ! Fraction of the soil with connected organic pathways
      ALLOCATE(fract_uncon        (nc, nsoil))  ! Fraction of the soil with unconnected organic matter
      ALLOCATE(hyd_cond_sat_uncon (nc, nsoil))  ! Saturated hydraulic conductivity of 'unconnected' fraction

      !$ACC DATA ASYNC(acc_stream) &
      !$ACC   CREATE(fract_perc, fract_uncon, hyd_cond_sat_uncon, hyd_cond_sat_org, vol_porosity_org) &
      !$ACC   CREATE(bclapp_org, matric_pot_org, pore_size_index_org, vol_field_cap_org, vol_p_wilt_org, vol_wres_org)

      !> In simulation with [[t_hydro_config:l_organic]]=true, these soil properties are updated using
      !> parameters for soil organic matter, that are available for the top layer and deeper soil layers.

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        hyd_cond_sat_org(ic,1) = hyd_cond_sat_org_top
        vol_porosity_org(ic,1) = vol_porosity_org_top
        bclapp_org(ic,1) = bclapp_org_top
        matric_pot_org(ic,1) = matric_pot_org_top
        pore_size_index_org(ic,1) = pore_size_index_org_top
        vol_field_cap_org(ic,1) = vol_field_cap_org_top
        vol_p_wilt_org(ic,1) = vol_p_wilt_org_top
        vol_wres_org(ic,1) = vol_wres_org_top
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is=2,nsoil
        DO ic=1,nc
          hyd_cond_sat_org(ic,is) = hyd_cond_sat_org_below
          vol_porosity_org(ic,is) = vol_porosity_org_below
          bclapp_org(ic,is) = bclapp_org_below
          matric_pot_org(ic,is) = matric_pot_org_below
          pore_size_index_org(ic,is) = pore_size_index_org_below
          vol_field_cap_org(ic,is) = vol_field_cap_org_below
          vol_p_wilt_org(ic,is) = vol_p_wilt_org_below
          vol_wres_org(ic,is) = vol_wres_org_below
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !> Saturated hydraulic conductivity is calculated from the parameters for mineral and organic soils
      !> following the percolation theory (see CLM45 Tech Note, section 7.4.1 Hydraulic Properties, p.161 f.).

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is=1,nsoil
        DO ic=1,nc
          IF (fract_org_sl(ic,is) > thresh_org) THEN
            ! Connected flow pathways consisting of organic material only exist
            N_perc = (1._wp - thresh_org)**(-beta_perc)
            fract_perc(ic,is)  = MIN(1._wp, N_perc * fract_org_sl(ic,is) * (fract_org_sl(ic,is)-thresh_org)**beta_perc)
            fract_uncon(ic,is) = 1._wp - fract_perc(ic,is)
          ELSE
            ! No connected organic matter pathways exist, and flow passes mineral and organic soil components in series
            fract_perc(ic,is)  = 0._wp
            fract_uncon(ic,is) = 1._wp
          END IF

          IF (fract_org_sl(ic,is) > 0._wp) THEN
            hyd_cond_sat_uncon(ic,is) = fract_uncon(ic,is)                                           &
              &                       * ( (1._wp - fract_org_sl(ic,is)) / hyd_cond_sat_sl(ic,is)     &
              &                            + (fract_org_sl(ic,is) - fract_perc(ic,is)) / hyd_cond_sat_org(ic,is) )**(-1._wp)
            hyd_cond_sat_sl(ic,is)    = fract_uncon(ic,is) * hyd_cond_sat_uncon(ic,is)  &
              &                       + fract_perc(ic,is) * hyd_cond_sat_org(ic,is)
          END IF
        END DO
      END DO
      !$ACC END LOOP

      !> The other soil properties are calculated from organic and mineral values according to the respective
      !> fractions.
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is=1,nsoil
        DO ic=1,nc
          vol_porosity_sl   (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_porosity_sl   (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_porosity_org(ic,is)
          bclapp_sl         (ic,is) = (1._wp - fract_org_sl(ic,is)) * bclapp_sl         (ic,is) &
            &                       + fract_org_sl(ic,is) * bclapp_org(ic,is)
          matric_pot_sl     (ic,is) = (1._wp - fract_org_sl(ic,is)) * matric_pot_sl     (ic,is) &
            &                       + fract_org_sl(ic,is) * matric_pot_org(ic,is)
          pore_size_index_sl(ic,is) = (1._wp - fract_org_sl(ic,is)) * pore_size_index_sl(ic,is) &
            &                       + fract_org_sl(ic,is) * pore_size_index_org(ic,is)
          vol_field_cap_sl  (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_field_cap_sl  (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_field_cap_org(ic,is)
          vol_p_wilt_sl     (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_p_wilt_sl     (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_p_wilt_org(ic,is)
          vol_wres_sl       (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_wres_sl       (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_wres_org(ic,is)
        END DO
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

      !$ACC END DATA

      ! Clean up
      DEALLOCATE(hyd_cond_sat_org)
      DEALLOCATE(vol_porosity_org)
      DEALLOCATE(bclapp_org)
      DEALLOCATE(matric_pot_org)
      DEALLOCATE(pore_size_index_org)
      DEALLOCATE(vol_field_cap_org)
      DEALLOCATE(vol_p_wilt_org)
      DEALLOCATE(vol_wres_org)

      DEALLOCATE(fract_perc)
      DEALLOCATE(fract_uncon)
      DEALLOCATE(hyd_cond_sat_uncon)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_properties

  ! ============================================================================================== !
  !>
  !> Aggregation of the soil properties
  !>
  !> This is the aggregation routine for variables of task "soil_properties". The variables are
  !> aggregate on the parent tiles as area weighted means of the child tile values.
  !>
  SUBROUTINE aggregate_soil_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Options of the current block

    dsl4jsb_Def_memory(HYDRO_)                              !< Memory of the hydrology process

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract   !< Aggregation method: area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_properties'

    INTEGER :: &
      & iblk, &                                             !< Current block index
      & ics, &                                              !< Index of first cell in block
      & ice                                                 !< Index of last cell in block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, hyd_cond_sat_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_porosity_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, bclapp_sl,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, matric_pot_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pore_size_index_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_field_cap_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_p_wilt_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_wres_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_org_sl,       weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_properties

  ! ============================================================================================== !
  !>
  !> Update routine for the soil hydrology
  !>
  !> The soil hydrology comprises hydrological processes in the soil. In this interface routine
  !> we prepare the data, but the actual calculations happen in process routine
  !> [[mo_hydro_process:calc_soil_hydrology]].
  !>
  SUBROUTINE update_soil_hydrology(tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o, rhoi
    USE mo_hydro_process,          ONLY: calc_soil_hydrology
    USE mo_hydro_constants,        ONLY: frac_wtr_vertical_transport_up_max, frac_wtr_vertical_transport_down_max, &
      &                                  frac_w_lat_loss_max
    USE mo_jsb_math_constants,     ONLY: eps8

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of the current block

    ! Local variables
    dsl4jsb_Def_config(HYDRO_)                   !< Configurable parameters of the hydrology process
    dsl4jsb_Def_memory(HYDRO_)                   !< Variables of the hydrology process
    dsl4jsb_Def_memory(SSE_)                     !< Variables of the "snow and snow energy" process
#ifndef __NO_QUINCY__
    dsl4jsb_Def_memory(Q_ASSIMI_)
#endif

    ! Pointers to variables in memory
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl        !< Amount of liquid water in soil layers [m]
    dsl4jsb_Real3D_onChunk :: ice_soil_sl        !< Amount of ice in soil layers [m water equivalent]
    dsl4jsb_Real3D_onChunk :: vol_weq_soil_sl    !< Volumetric water content (liquid+ice) in soil layers
                                                 !< [m3 water equiv. / m3]
    dsl4jsb_Real3D_onChunk :: frac_wtr_transp_down_sl !< Faction of water transported vertically to the below layer []
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl !< Fraction of soil water that is lost laterally (horizontally) []
    dsl4jsb_Real3D_onChunk :: drainage_sl        !< Drainage from the soil layers [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: wtr_transp_down_sl !< Downward flux of water to the soil layer below [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: soil_depth_sl      !< (Active) soil depth within the soil layers [m]

    dsl4jsb_Real2D_onChunk :: wtr_soilhyd_res    !< Residual of the vertical soil water transport scheme [m]
    dsl4jsb_Real2D_onChunk :: wtr_pond           !< Water content of the pond reservoir [m]
    dsl4jsb_Real2D_onChunk :: ice_pond           !< Amount of ice in the pond reservoir [m water equivalent]
    dsl4jsb_Real2D_onChunk :: weq_pond           !< Water/ice content of the pond reservoir [m water equiv.]

    dsl4jsb_Real2D_onChunk :: thaw_depth_max     !< Maximum thaw depth in current year [m]
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1 !< Maximum thaw depth of previous year [m]
    dsl4jsb_Real2D_onChunk :: steepness          !< Parameter representing subgrid scale slopes []
    dsl4jsb_Real2D_onChunk :: soil_depth         !< Soil depth until bedrock [m]

    ! Locally allocated vectors
    TYPE(t_jsb_model), POINTER :: model          !< This instance of jsbach
    TYPE(t_jsb_grid),  POINTER :: grid           !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w         !< Vertical grid used in hydrology

    LOGICAL  :: ltpe_open                        !< For terraplanet setup: Deep soils are kept wet
    LOGICAL  :: ltpe_closed                      !< For terraplanet setup: Water balance is closed
    LOGICAL  :: l_latflow_to_streamflow          !< With "[[t_hydro_config:hydro_scale]]=Uniform" True:
                                                 !< outflow of intermediary runoff/drainage storage goes to streamflow

    INTEGER  :: &
      & iblk,   &                                !< Index of current block
      & ics,    &                                !< Index of first cell in block
      & ice,    &                                !< Index of last cell in block
      & nc,     &                                !< Number of cells in block
      & ic,     &                                !< Looping index for grid cells
      & is,     &                                !< Looping index for soil layers
      & nsoil                                    !< Number of vertical layers (grid dimension)
    REAL(wp) :: dtime                            !< Time step length [s]

    INTEGER :: last_soil_layer(options%nc)       !< Deepest active soil layer index (above the bedrock boundary)

    REAL(wp), POINTER :: &
      & area(:), &                               !< Grid cell area [m2]
      & lat(:), &                                !< Grid cell center latitude [deg.]
      & lon(:), &                                !< Grid cell center longitude [deg.]
      & tile_fract(:)                            !< Grid cell fraction of the tile []

    REAL(wp) :: &
      & slope_used(options%nc), &                !< Effective slope on tile []
      & flowlag_used(options%nc), &              !< Lag factor for surface runoff (time step dependent) []
      & depth_thrsh(options%nc), &               !< Maximum active layer depth up to which permafrost
                                                 !< affects bottom layer drainage [m]
      & drainage_lowest_soil_layer(options%nc), & !< Bottom layer drainage [kg m-2 s-1]
      & ret_macro_blg, &                         !< Retention time for below ground runoff [s]
      & wtr_flux_in_m                            !< Variable used for water fluxes in [m / (time step)]

    LOGICAL :: &
      & l_fract(options%nc), &                   !< Tile has a non-zero grid cell fraction
      & l_pf_soil(options%nc)                    !< Permafrost affects bottom layer drainage
    REAL(wp), ALLOCATABLE :: wtr_soil_old_sl(:,:) !< wtr_soil_sl of previous timestep

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    nsoil = options%nsoil_w
    dtime = options%dtime

    ! Only proceed if the hydrology calculations are needed on the current tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    ! As lakes have no soil, the soil hydrology does not need to be calculated on lake tiles.
    IF (tile%is_lake) RETURN

    ! Runoff and drainage for glaciers already computed in update_surface_hydrology
    ! Glaciers have no water in soil
    IF (tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)
    soil_w => Get_vgrid('soil_depth_water')

    tile_fract => tile%fract(ics:ice, iblk)
    area       => grid%area (ics:ice, iblk)
    lat        => grid%lat  (ics:ice, iblk)
    lon        => grid%lon  (ics:ice, iblk)

    dsl4jsb_Get_config(HYDRO_)

    ! Variables used in terraplanet setup, i.e. in a land-only world without oceans
    ltpe_open = .FALSE.      ! Deep soils are artificially kept wet
    ltpe_closed = .FALSE.    ! Closed water balance: no runoff or drainage
    IF (model%config%tpe_scheme == 'open')   ltpe_open = .TRUE.
    IF (model%config%tpe_scheme == 'closed') ltpe_closed = .TRUE.

    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)      ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)      ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_weq_soil_sl)  ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, soil_depth)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, steepness)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_soilhyd_res)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)         ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)         ! out

    dsl4jsb_Get_var2D_onChunk(SSE_, thaw_depth_max)     ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, thaw_depth_max_ym1) ! in

    !$ACC DATA ASYNC(acc_stream) CREATE(l_fract) &
    !$ACC   CREATE(l_pf_soil, depth_thrsh, slope_used, flowlag_used) &
    !$ACC   CREATE(drainage_lowest_soil_layer, last_soil_layer)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ! Initialize variable that stores the number of active layers
      last_soil_layer(ic) = 0
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (soil_depth_sl(ic, is) > 0._wp) THEN
          last_soil_layer(ic) = is
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    IF (ANY(last_soil_layer(:) < 1)) CALL finish('soilhyd', 'Problem with no. of active soil layers (=0)')
#endif

    dsl4jsb_Get_var3D_onChunk(HYDRO_, drainage_sl)              ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_transp_down_sl)       ! out

    SELECT CASE (model%config%model_scheme)
#ifndef __NO_QUINCY__
    CASE (MODEL_QUINCY)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)  ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_w_lat_loss_sl)       ! out
      ! Save the soil-water state at the beginning of the time step
      !   for the calculation of fractional water transport in the water column
      ALLOCATE(wtr_soil_old_sl(nc, nsoil))
      !$ACC ENTER DATA CREATE(wtr_soil_old_sl) ASYNC(acc_stream)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          wtr_soil_old_sl(ic, is) = wtr_soil_sl(ic, is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
#endif
    END SELECT

    ! Find out the tile's area fraction. Some processes are only calculated in grid cells with a relevant
    ! tile fraction, e.g. not over complete ocean cells.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      l_fract(ic) = tile_fract(ic) > 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    ! In deep soils permafrost affects bottom layer drainage as long as the active layer is smaller than soil
    ! depth. In shallow soils, bottom layer drainage is also affected if the whole soil column is in the
    ! active zone, but permafrost remains (below bedrock) in less then 1m depth.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      depth_thrsh(ic) = MAX(soil_depth(ic), 1.0_wp)
      l_pf_soil(ic)   = (thaw_depth_max_ym1(ic) < depth_thrsh(ic) .AND. &
        &                thaw_depth_max(ic)     < depth_thrsh(ic))
    END DO
    !$ACC END PARALLEL LOOP

    ! Set to true, until HydroTiles are implemented
    l_latflow_to_streamflow = .TRUE.
    ret_macro_blg = dsl4jsb_Config(HYDRO_)%ret_macro_blg
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      slope_used(ic)   = steepness(ic)
      flowlag_used(ic) = dtime / ret_macro_blg
    END DO
    !$ACC END PARALLEL LOOP

    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      CALL calc_soil_hydrology( &
        ! in
        & nc, &
        & l_fract(:), &
        & l_pf_soil(:), &
        & lat(:), &
        & lon(:), &
        & nsoil, &
        & dtime, &
        & ltpe_closed, ltpe_open,                                 & ! Options for terraplanet setup
        & dsl4jsb_Config(HYDRO_)%enforce_water_budget,            & ! Options for water balance check
        & l_latflow_to_streamflow,                                & ! True: water from lateral flow reservoir
                                                                    ! directly goes to streamflow
        & dsl4jsb_Config(HYDRO_)%soilhydmodel,                    & ! Model for soil hydraulic properties
        & dsl4jsb_Config(HYDRO_)%interpol_mean,                   & ! Interpolation scheme for hydraulic
                                                                    ! conductivity
        & dsl4jsb_Config(HYDRO_)%hydro_scale,                     & ! Subgrid scale inhomogeneity option
        & model%config%model_scheme,                              &
        & dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract,               & ! Rel. root zone moisture at wilting point []
        & dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl     ),     & ! Soil depth until bedrock [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl     ),     & ! Root depth [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, hyd_cond_sat_sl   ),     & ! Saturated hydraulic conductivity [m/s]
        & dsl4jsb_var3D_onChunk (HYDRO_, matric_pot_sl     ),     & ! Saturated soil matric potential [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, bclapp_sl         ),     & ! Clapp and Hornberger exponent b []
        & dsl4jsb_var3D_onChunk (HYDRO_, pore_size_index_sl),     & ! Soil pore size index []
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_porosity_sl   ),     & ! Volumetric soil porosity []
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_field_cap_sl  ),     & ! Volumetric soil field capacity []
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_p_wilt_sl     ),     & ! Volumetric soil wilting point []
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_wres_sl       ),     & ! Volumetric residual soil water content [m/m]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pot_scool_sl ), & ! Potential amount of supercooled water [m]
        & slope_used(:),                                          & ! Effective slope on tile []
        & flowlag_used(:),                                        & ! Time lag factor for below ground drainage []
        & dsl4jsb_var2D_onChunk (HYDRO_, fract_pond_max    ),     & ! Area fraction with surface depressions []
        & dsl4jsb_var2D_onChunk (HYDRO_, evapotrans_soil   ),     & ! Evapotranspiration from soil (w/o snow and
                                                                    ! skin res.) [kg m-2 s-1]
        & dsl4jsb_var2D_onChunk (HYDRO_, transpiration     ),     & ! Transpiration [kg m-2 s-1]
        & dsl4jsb_var2D_onChunk (HYDRO_, ice_pond          ),     & ! Amount of ice in surface depressions
                                                                    ! [m water equivalent]
        & dsl4jsb_var2D_onChunk (HYDRO_, weq_pond_max      ),     & ! Maximum water/ice content in pond reservoir
                                                                    ! [m water equivalent]
        ! inout
        & dsl4jsb_var2D_onChunk (HYDRO_, infilt            ),     & ! Infiltration [kg m-2 s-1]
        & dsl4jsb_var2D_onChunk (HYDRO_, runoff            ),     & ! (Surface) runoff [kg m-2 s-1]
        & wtr_soil_sl,                                            & ! Amount of liquid soil water [m]
        & ice_soil_sl,                                            & ! Amount of soil ice [m water equivalent]
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond          ),     & ! Amount of water in surface depressions [m]
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond_net_flx  ),     & ! Net inflow into surface water ponds
                                                                    ! [kg m-2 s-1]
        & dsl4jsb_var2D_onChunk (HYDRO_, tpe_overflow      ),     & ! Terraplanet setup: Overflow reservoir [m]
        & dsl4jsb_var2D_onChunk (HYDRO_, evapo_deficit     ),     & ! Evaporation deficit due to inconsistencies
                                                                    ! [m water equivalent / (time step)]
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_latflow_res_srf),    & ! Water in lateral surface runoff reservoir [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_latflow_res_sl),     & ! Water in lateral subsurface drainage
                                                                    ! reservoirs [m]
        ! out
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_sat_sl   ),     & ! Soil water content at saturation [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_fc_sl    ),     & ! Soil field capacity [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pwp_sl   ),     & ! Soil water content at wilting point [m]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_res_sl   ),     & ! Residual soil water content [m]
        & dsl4jsb_var2D_onChunk (HYDRO_, runoff_dunne      ),     & ! Dunne component of surface runoff [kg m-2 s-1]
        & dsl4jsb_var2D_onChunk (HYDRO_, drainage          ),     & ! Subsurface drainage [kg m-2 s-1]
        & dsl4jsb_var3D_onChunk (HYDRO_, drainage_sl       ),     & ! Subsurface drainage on soil layers [kg m-2 s-1]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_transp_down_sl),     & ! Vertical flux to the below soil layer [kg m-2 s-1]
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_soilhyd_res   ),     & ! Residual of vertical soil water transport
                                                                    ! scheme [m]
        & drainage_lowest_soil_layer(:),                          & ! Bottom layer drainage [kg m-2 s-1]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_latflow_sl    )      & ! Outflow of wtr_latflow_res_sl [kg m-2 s-1]
        & )

#ifndef __NO_QUINCY__
    CASE (MODEL_QUINCY)
      dsl4jsb_Get_memory(Q_ASSIMI_)

      CALL calc_soil_hydrology( &
        ! in
        & nc, &
        & l_fract(:), &
        & l_pf_soil(:), &
        & lat(:), &
        & lon(:), &
        & nsoil, &
        & dtime, &
        & ltpe_closed, ltpe_open,                                 &
        & dsl4jsb_Config(HYDRO_)%enforce_water_budget,            & ! Water balance check setting
        & l_latflow_to_streamflow,                                &
        & dsl4jsb_Config(HYDRO_)%soilhydmodel,                    &
        & dsl4jsb_Config(HYDRO_)%interpol_mean,                   &
        & dsl4jsb_Config(HYDRO_)%hydro_scale,                     & ! choice of hydrological scale scheme
        & model%config%model_scheme,                              &
        & dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract,               &
        & dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl     ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl     ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, hyd_cond_sat_sl   ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, matric_pot_sl     ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, bclapp_sl         ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, pore_size_index_sl),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_porosity_sl   ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_field_cap_sl  ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_p_wilt_sl     ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, vol_wres_sl       ),     & ! volumetric residual soil water content [m m-1]
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pot_scool_sl ), &
        & slope_used(:),                                          &
        & flowlag_used(:),                                        &
        & dsl4jsb_var2D_onChunk (HYDRO_, fract_pond_max    ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, evapotrans_soil   ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, transpiration     ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, ice_pond          ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, weq_pond_max      ),     &
        ! inout
        & dsl4jsb_var2D_onChunk (HYDRO_, infilt            ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, runoff            ),     &
        & wtr_soil_sl,                                            &
        & ice_soil_sl,                                            &
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond          ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond_net_flx  ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, tpe_overflow      ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, evapo_deficit     ),     &
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_latflow_res_srf),    &
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_latflow_res_sl),     &
        ! out
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_sat_sl   ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_fc_sl    ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pwp_sl   ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_res_sl   ),     & ! residual soil water content
        & dsl4jsb_var2D_onChunk (HYDRO_, runoff_dunne      ),     & ! dunne component of surface runoff
        & dsl4jsb_var2D_onChunk (HYDRO_, drainage          ),     &
        & dsl4jsb_var3D_onChunk (HYDRO_, drainage_sl       ),     & ! subsurface drainage on each soil layer
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_transp_down_sl),     & ! vertical downward transport of water
        & dsl4jsb_var2D_onChunk (HYDRO_, wtr_soilhyd_res   ),     &
        & drainage_lowest_soil_layer(:),                          & ! Bottom layer drainage
        & dsl4jsb_var3D_onChunk (HYDRO_, wtr_latflow_sl    ),     &
        ! optional in
        & dsl4jsb_var3D_onChunk (Q_ASSIMI_, ftranspiration_sl)    & ! only used with QUINCY
        & )

      ! calculate vertical and horizontal fraction of water transport for the use in BGCM soil transport calculations
      !   (using the soil water content of the previous timestep here (wtr_soil_old_sl) because the frac_w_lat_loss_sl
      !    is used for BGCM soil transport calculations of the material (elements / nutrients) from the previous timestep)
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG VECTOR PRIVATE(wtr_flux_in_m)
        DO ic = 1, nc
          IF (wtr_soil_old_sl(ic, is) > eps8) THEN
            ! calculation of the two fractions
            ! vertical (across soil layers)
            wtr_flux_in_m = wtr_transp_down_sl(ic, is) * dtime / rhoh2o
            frac_wtr_transp_down_sl(ic, is) = &
              &  MIN( MAX( wtr_flux_in_m / wtr_soil_old_sl(ic, is), frac_wtr_vertical_transport_up_max), &
              &       frac_wtr_vertical_transport_down_max) * 1000._wp / dtime
            ! lateral (i.e., horizontal)
            wtr_flux_in_m = drainage_sl(ic, is) * dtime / rhoh2o
            IF (is == last_soil_layer(ic)) THEN
              wtr_flux_in_m = wtr_flux_in_m + drainage_lowest_soil_layer(ic) * dtime / rhoh2o
            END IF
            frac_w_lat_loss_sl(ic, is) = &
              &  MIN( MAX(wtr_flux_in_m / wtr_soil_old_sl(ic, is), 0.0_wp), frac_w_lat_loss_max) * 1000._wp / dtime
            ! ensure the sum of the two fractions would not increase beyond frac_w_lat_loss_max
            IF ((frac_wtr_transp_down_sl(ic, is) + frac_w_lat_loss_sl(ic, is)) > frac_w_lat_loss_max) THEN
              frac_wtr_transp_down_sl(ic, is) = frac_wtr_transp_down_sl(ic, is) * frac_w_lat_loss_max &
                &                               / (frac_wtr_transp_down_sl(ic, is) + frac_w_lat_loss_sl(ic, is))
              frac_w_lat_loss_sl(ic, is)      = frac_w_lat_loss_sl(ic, is) * frac_w_lat_loss_max &
                &                               / (frac_wtr_transp_down_sl(ic, is) + frac_w_lat_loss_sl(ic, is))
            END IF
          ELSE
            frac_wtr_transp_down_sl(ic, is) = 0.0_wp
            frac_w_lat_loss_sl(ic, is)      = 0.0_wp
          END IF
        END DO
      END DO
      !$ACC END PARALLEL
#endif
    END SELECT

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)

    ! Update weq_pond storage due to potential soil infiltration overflow
    !$ACC LOOP GANG VECTOR
    DO ic = 1, nc
      weq_pond(ic) = wtr_pond(ic) + ice_pond(ic)
    END DO

    ! Convert soil hydrology vertical transport residual to
    ! volume [m] -> [m3]
    !$ACC LOOP GANG VECTOR
    DO ic = 1, nc
      wtr_soilhyd_res(ic) = wtr_soilhyd_res(ic) * tile_fract(ic) * area(ic)
    END DO
    !ACC END LOOP

    !$ACC END PARALLEL

    ! Compute volumetric soil moisture (liquid + ice) in m3/m3 (= m/m) water equivalent
    ! Note that the soil ice content needs to be converted to water equivalent first.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil
      DO ic = 1, nc
        vol_weq_soil_sl(ic,is) = (wtr_soil_sl(ic,is) + ice_soil_sl(ic,is) * (rhoi / rhoh2o)) / soil_w%dz(is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    SELECT CASE (model%config%model_scheme)
#ifndef __NO_QUINCY__
    CASE (MODEL_QUINCY)
      !$ACC EXIT DATA DELETE(wtr_soil_old_sl) ASYNC(acc_stream)
      DEALLOCATE(wtr_soil_old_sl)
#endif
    END SELECT

    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_hydrology

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of the soil hydrology variables
  !>
  !> This is the aggregation routine for the variables of the "soil_hydrology" task. The variables are
  !> aggregate on the parent tiles as area weighted means of the child tile values.
  !>
  SUBROUTINE aggregate_soil_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< This tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Current block

    dsl4jsb_Def_memory(HYDRO_)                              !< Configurable parameters of the HYDRO process

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract   !< Aggregation method for area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_hydrology'

    INTEGER :: &
      & iblk,  &                                           !< Block index
      & ics,   &                                           !< Index of first cell in block
      & ice                                                !< Index of last cell in block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_soil_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_weq_soil_sl,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sat_sl,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_fc_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pwp_sl,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_res_sl,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond_net_flx,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, infilt,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff_dunne,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_deficit,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soilhyd_res,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_transp_down_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_res_srf, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_res_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_sl,      weighted_by_fract)
#ifndef __NO_QUINCY__
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      dsl4jsb_Aggregate_onChunk(HYDRO_, frac_wtr_transp_down_sl,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(HYDRO_, frac_w_lat_loss_sl,          weighted_by_fract)
    END SELECT
#endif
    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_hydrology

  ! ======================================================================================================== !
  !>
  !>#### Calculation of canopy conductance assuming no water stress
  !>
  !> The "update" routine of task "canopy_cond_unstressed" calculates the canopy conductance in case of
  !> unlimited water availability, which is a first step to calculate the actual canopy conductance in routine
  !> [[update_canopy_cond_stressed]].
  !>
  SUBROUTINE update_canopy_cond_unstressed(tile, options)

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Options of the current block

    dsl4jsb_Def_memory(HYDRO_)                              !< Memory of hydrology process
    dsl4jsb_Def_memory(A2L_)                                !< Memory of "atmosphere to land" process
#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_memory(PHENO_)                              !< Memory of phenology process
#endif
#ifndef __NO_QUINCY__
    dsl4jsb_Def_memory(VEG_)                                !< Memory of vegetation process
#endif

    TYPE(t_jsb_model), POINTER :: model   !!!vg check !!!

    dsl4jsb_Real2D_onChunk :: lai                           !< Leaf area index []
    dsl4jsb_Real2D_onChunk :: swpar_srf_down                !< Photosynthetically active radiation [W/m2]
    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited         !< Canopy conductivity with unlimited water supply

    INTEGER :: &
      & iblk,  &                                            !< Current block index
      & ics,   &                                            !< Index of first cell in block
      & ice,   &                                            !< Index of last cell in block
      & nc,    &                                            !< Number of cells in block
      & ic                                                  !< Looping index for the cells

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_unstressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! Do nothing unless we are on a vegetation tile and hydrology calculations need to be done on this tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)
    SELECT CASE (model%config%model_scheme)
#ifndef __QUINCY_STANDALONE__
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_memory(PHENO_)
      dsl4jsb_Get_var2D_onChunk(PHENO_,  lai)                   ! in
#endif
#ifndef __NO_QUINCY__
    CASE (MODEL_QUINCY)
      dsl4jsb_Get_memory(VEG_)
      dsl4jsb_Get_var2D_onChunk(VEG_,    lai)                   ! in
#endif
    END SELECT

    dsl4jsb_Get_var2D_onChunk(A2L_,      swpar_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_unlimited) ! out

    !> The main calculation of the stomatal canopy conduction for unstressed conditions happens in
    !> process routine [[mo_hydro_process:get_canopy_conductance]].
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      canopy_cond_unlimited(ic) = get_canopy_conductance( lai(ic), swpar_srf_down(ic) )
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE update_canopy_cond_unstressed

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of canopy conductance assuming no water stress
  !>
  !> This is the aggregation routine for variables of task "canopy_cond_unstressed". Canopy conductance
  !> of different vegetation tiles is aggregated to the parent tiles as area weighted mean.
  !>
  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile          !< This tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options       !< Options of the current block

    dsl4jsb_Def_memory(HYDRO_)                                 !< Memory with the hydrology variables

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract      !< Aggregation method for area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER :: &
      & iblk, &                                                !< Index of this block
      & ics, &                                                 !< Index of first grid cell in the block
      & ice                                                    !< Index of last grid cell in the block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, canopy_cond_unlimited, weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  ! ======================================================================================================== !
  !>
  !>#### Calculate root zone soil moisture and water stress
  !>
  !> Actual and relative root zone water and/or ice content as well as the water stress factor are updated
  !> in this routine.
  !>
  SUBROUTINE update_water_stress(tile, options)

    ! Used variables
    USE mo_hydro_process,      ONLY: get_water_stress_factor
    USE mo_hydro_util,         ONLY: get_amount_in_rootzone

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile        !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options     !< Parameters for the current block

    TYPE(t_jsb_model), POINTER :: model                      !< This instance of ICON-Land

    dsl4jsb_Def_config(HYDRO_)                !< Configurable parameters for the hydrology
    dsl4jsb_Def_config(SSE_)                  !< Configurable parameters of the "soil and snow energy" process

    dsl4jsb_Def_memory(HYDRO_)                !< Memory of the hydrology process

    dsl4jsb_Real3D_onChunk :: soil_depth_sl          !< Soil depth until bedrock within the soil layers [m]
    dsl4jsb_Real3D_onChunk :: root_depth_sl          !< Root depth within the soil layers [m]
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl       !< Volumetric soil field capacity []
    dsl4jsb_Real3D_onChunk :: vol_p_wilt_sl          !< Volumetric permanent wilting point []
    dsl4jsb_Real2D_onChunk :: wtr_rootzone           !< Liquid water content in the root zone [m]
    dsl4jsb_Real2D_onChunk :: ice_rootzone           !< Frozen water content in the root zone [m]
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl            !< Liquid water content in soil layers [m]
    dsl4jsb_Real3D_onChunk :: ice_soil_sl            !< Frozen water content in soil layers [m]
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl  !< Potentially supercooled water content in soil layers [m]
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_scool_pot !< Potentially supercooled water content in root zone [m]
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_scool_act !< Actual supercooled water content in root zone [m]
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_avail     !< Plant available liquid water content in root zone [m]
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_avail_max !< Maximum plant available liquid water in root zone [m]
    dsl4jsb_Real2D_onChunk :: wpi_rootzone_max       !< Maximum possible amount of water and/or ice in the
                                                     !< root zone [m]
    dsl4jsb_Real2D_onChunk :: water_stress           !< Water stress factor: 0 no stress, 1 maximum stress []
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_rel       !< Relative root zone soil moisture []
    dsl4jsb_Real2D_onChunk :: wtr_plant_avail_rel    !< Relative plant available soil moisture in root zone []
                                                     !< (used with quincy; not including water below the wilting
                                                     !< point.)
                                                     ! Note: with jsbach the wilting point in account for later.
                                                     ! TODO: do we need both variables?

    INTEGER :: &
      & iblk, &                                      !< Current block index
      & ics, &                                       !< First cell in block
      & ice, &                                       !< Last cell in block
      & nc, &                                        !< Number of cells in block
      & ic, &                                        !< Looping index for cells
      & nsoil, &                                     !< Number of vertical soil layers (grid dimension)
      & is                                           !< Looping index for soil layers

    REAL(wp) :: &
      & config_w_soil_crit_fract, &                  !< Below this critical relative soil moisture, plants
                                                     !< start to suffer from water stress.
      & config_w_soil_wilt_fract                     !< Below this relative soil moisture (wilting point)
                                                     !< transpiration is no longer possible.
    REAL(wp) :: &
      & wtr_soil_pot_scool_sl_min(options%nc, options%nsoil_w), & !< Actual amount of supercooled water [m]
      & wpi_fcap_soil_sl(options%nc, options%nsoil_w)             !< Volumetric amount of water and/or ice
                                                                  !< at field capacity [m] (???)
    REAL(wp) :: wpi_pwp_soil_sl(options%nc, options%nsoil_w)      !< Volumetric amount of water and/or ice
                                                                  !< at wilting point [m] (???)
    REAL(wp) :: wtr_rootzone_pwp(options%nc)                      !< Amount of liquid water at wilting point [m]
    REAL(wp) :: hlp1                                              !< helper variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_stress'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    nsoil = options%nsoil_w

    ! Do nothing unless this is a vegetation tile and hydrology calculations are needed on this tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Set pointers to the types containing configurable parameters of the respective processes
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var3D_onChunk(HYDRO_,    soil_depth_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    vol_field_cap_sl)       ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    vol_p_wilt_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_pot_scool_sl)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone)           ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    root_depth_sl)          ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wpi_rootzone_max)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    ice_rootzone)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_scool_pot) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_scool_act) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_avail)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_avail_max) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_rel)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_plant_avail_rel)    ! out

    ! Define plain local variables, needed on GPUs
    config_w_soil_crit_fract = dsl4jsb_Config(HYDRO_)%w_soil_crit_fract
    config_w_soil_wilt_fract = dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract

    !$ACC DATA CREATE(wtr_soil_pot_scool_sl_min, wpi_fcap_soil_sl) &
    !$ACC   CREATE(wpi_pwp_soil_sl, wtr_rootzone_pwp) ASYNC(acc_stream)

    !> Water stress is calculated from the ratio between the actual plant available water
    !> "wtr_rootzone_avail" and the maximum potentielly available water "wtr_rootzone_avail_max" between the
    !> wilting point and the "critical" soil moisture, where plants feel no water limitation.
    !>

    ! Calculate actual plant available water
    !
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:),  &
      &                         soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    ! Calculate root zone field capacity
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil
      DO ic = 1, nc
        wpi_fcap_soil_sl(ic,is) = vol_field_cap_sl(ic,is) * soil_depth_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    CALL get_amount_in_rootzone(wpi_fcap_soil_sl(:,:),  &
      &  soil_depth_sl(:,:), root_depth_sl(:,:), wpi_rootzone_max(:))

#ifndef __NO_QUINCY__
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      ! Compute minimum plant available water in the root zone water based on soil permanent wilting point
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          wpi_pwp_soil_sl(ic,is) = vol_p_wilt_sl(ic,is) * soil_depth_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      CALL get_amount_in_rootzone( &
        &   wpi_pwp_soil_sl(:,:), &   ! in
        &   soil_depth_sl(:,:), &     ! in
        &   root_depth_sl(:,:), &     ! in
        &   wtr_rootzone_pwp(:) &     ! inout
        &   )

    !> Available water is corrected by ice and supercooled water content below,
    !> helping to ensure a consistent definition of water at permanent wilting point.
    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:), &
        &                         soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
      !$ACC   PRIVATE(hlp1)
      DO ic = 1, nc
        ! assume ice and supercooled water are distributed equally between water below and above the PWP
        hlp1 = wtr_rootzone_pwp(ic) / wpi_rootzone_max(ic)
        wtr_rootzone_pwp(ic) = MAX(wtr_rootzone_pwp(ic) - ice_rootzone(ic) * hlp1 - wtr_rootzone_scool_pot(ic) * hlp1, 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE IF (dsl4jsb_Config(SSE_)%l_freeze .AND. .NOT. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:), &
        &                         soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
      !$ACC   PRIVATE(hlp1)
      DO ic = 1, nc
        hlp1 = wtr_rootzone_pwp(ic) / wpi_rootzone_max(ic)
        wtr_rootzone_pwp(ic) = MAX(wtr_rootzone_pwp(ic) - ice_rootzone(ic) * hlp1, 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    END IF
    END SELECT
#endif

    !> By plant available water we mean the liquid water content of the root zone, not including supercooled
    !> water. However, water below the wilting point is included - although it is not available to plants.

    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      ! Calculate the potential amount of supercooled water within the root zone. This is the amount of
      ! supercooled water that could remain liquid with current temperatures, if there was enough soil water
      ! (compare [[mo_sse_process:Get_liquid_max]]).
      CALL get_amount_in_rootzone(wtr_soil_pot_scool_sl(:,:),  &
        &                         soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone_scool_pot(:))

      ! Calculate the actual amount of supercooled water within the root zone.
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is=1,nsoil
        DO ic=1,nc
          wtr_soil_pot_scool_sl_min(ic,is) = MIN(wtr_soil_pot_scool_sl(ic,is), wtr_soil_sl(ic,is))
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      CALL get_amount_in_rootzone(wtr_soil_pot_scool_sl_min(:,:),  &
        &                         soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone_scool_act(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)

      ! Subtract the supercooled water from the plant available water, as supercooled water is not available
      ! for plants.
      DO ic=1,nc
        wtr_rootzone_avail(ic) = MAX(wtr_rootzone(ic) - wtr_rootzone_scool_act(ic), 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      ! In setups without soil water phase changes all (liquid) water of the root zone is available to plants
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        wtr_rootzone_avail(ic) = wtr_rootzone(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Calculate the soil field capacity
    !>
    !> As the actual root zone moisture, the potential root zone moisture is reduced by the amount of ice
    !> - otherwise plants suffer constant water stress in permafrost-regions. Besides, we also subtract the
    !> amount of potential supercooled water, as supercooled water is not accounted as plant available water.
    !> As for the actual available water, water below the wilting point is not subtracted from the
    !> potentially available water.

    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        wtr_rootzone_avail_max(ic) = MAX(wpi_rootzone_max(ic) - ice_rootzone(ic) - wtr_rootzone_scool_pot(ic), 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE IF (dsl4jsb_Config(SSE_)%l_freeze .AND. .NOT. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        wtr_rootzone_avail_max(ic) = MAX(wpi_rootzone_max(ic) - ice_rootzone(ic), 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        wtr_rootzone_avail_max(ic) = wpi_rootzone_max(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Calculate water stress factor
    !
    ! TODO: maybe use the distributed permanent wilting point/field cap from above and change get_water_stress_factor function
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      water_stress(ic) = &
        & get_water_stress_factor(wtr_rootzone_avail(ic), wtr_rootzone_avail_max(ic), &
            &                 config_w_soil_crit_fract, config_w_soil_wilt_fract)
    END DO
    !$ACC END PARALLEL LOOP

    ! Calculate relative soil moisture
    !>
    !> The logistic growth phenology (LoGro-P; compare [[mo_pheno_interface:update_phenology_logrop]]) has
    !> a different approach to account for water stress: there we use the relative soil moisture as measure
    !> for water stress directly, instead of using the water stress factor.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic= 1, nc
      IF (wtr_rootzone_avail_max(ic) > 0._wp) THEN
        wtr_rootzone_rel(ic) =  MIN(wtr_rootzone_avail(ic) / wtr_rootzone_avail_max(ic), 1._wp)
      ELSE
        ! wpi_rootzone_max is zero for glacier tiles
        wtr_rootzone_rel(ic) = 0._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

#ifndef __NO_QUINCY__
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      ! Calculate relative soil moisture (used in Quincy Vegetation)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF ((wtr_rootzone_avail_max(ic) - wtr_rootzone_pwp(ic)) > 0.0_wp) THEN
          wtr_plant_avail_rel(ic) = MIN(MAX(wtr_rootzone_avail(ic) - wtr_rootzone_pwp(ic), 0.0_wp) &
            &                       / (wtr_rootzone_avail_max(ic) - wtr_rootzone_pwp(ic)), 1.0_wp)
        ELSE
          ! wpi_rootzone_max is zero for glacier tiles
          wtr_plant_avail_rel(ic) = 0.0_wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END SELECT
#endif

    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_stress

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of water stress and root zone moisture
  !>
  !> In this routine we aggregate variables that are updated in the "water stress" task. This is not only
  !> water stress but also several variables for root zone water or ice content. The variables are aggregate
  !> as area weighted means of the child tile values.
  !>
  SUBROUTINE aggregate_water_stress(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Parameters of the current block

    dsl4jsb_Def_memory(HYDRO_)                              !< Memory of the hydrology process

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract   !< Method for area weighted aggregation

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_stress'

    INTEGER :: &
      & iblk, &             !< Index of current block
      & ics, &              !< First grid cell of the block
      & ice                 !< Last grid cell of the block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set reference to the memory
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wpi_rootzone_max,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_rootzone,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_rel,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_avail,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_avail_max, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_stress,           weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_stress

  ! ======================================================================================================== !
  !>
  !>#### Calculate canopy conductance under water stress
  !>
  !> In this routine we calculate canopy conductance accounting for water stress.
  !>
  SUBROUTINE update_canopy_cond_stressed(tile, options)

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance
    USE mo_phy_schemes,   ONLY: qsat_water, qsat_mixed

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Parameters of nproma block

    TYPE(t_jsb_model), POINTER :: model                     !< This ICON-Land instance

    dsl4jsb_Def_memory(SEB_)           !< Memory of process "surface energy balance"
    dsl4jsb_Def_memory(A2L_)           !< Memory of process "atmosphere to land"
    dsl4jsb_Def_memory(HYDRO_)         !< Memory of the hydrology process

    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited   !< Canopy conductivity without water stress [m/s]
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited     !< Canopy conductivity accounting for water stress [m/s]
    dsl4jsb_Real2D_onChunk :: water_stress            !< Water stress factor []
    dsl4jsb_Real2D_onChunk :: t                       !< Surface temperature [K]
    dsl4jsb_Real2D_onChunk :: q_air                   !< Specific humidity of the air []
    dsl4jsb_Real2D_onChunk :: press_srf               !< Surface pressure [Pa]

    INTEGER :: &
      & iblk, &                        !< Current block index
      & ics, &                         !< Index of first cell of the block
      & ice, &                         !< Index of last cell of the block
      & nc, &                          !< Number of cells in this block
      & ic                             !< Looping index for grid cells

    LOGICAL :: &
      & q_air_gt_qsat_tmp, &           !< True: spec. humidity exceeds saturation
      & use_tmx                        !< True: use TMX scheme

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! Only proceed if this is a vegetation tile and hydrology is calculated on this tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    use_tmx = model%config%use_tmx

    ! Get reference to process memories
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Set pointers to variables for current block
    dsl4jsb_Get_var2D_onChunk(A2L_,      q_air)                  ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)              ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)           ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_unlimited)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_limited)    ! out

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(q_air_gt_qsat_tmp)
    DO ic=1,nc
      ! Find out if atmospheric moisture is saturated
      IF (use_tmx) THEN
        q_air_gt_qsat_tmp = q_air(ic) > qsat_water(t(ic),press_srf(ic))
      ELSE
        q_air_gt_qsat_tmp = q_air(ic) > qsat_mixed(t(ic),press_srf(ic))
      END IF
      ! Compute the actual canopy (stomatal) conductance under water stress
      canopy_cond_limited(ic) = get_canopy_conductance(canopy_cond_unlimited(ic), & ! in, unstressed canopy conductance
                                                       water_stress(ic),          & ! in, water stress factor
                                                       q_air_gt_qsat_tmp          & ! in, atmosphere saturated?
                                                      )
    END DO
    !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_stressed

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of water limited canopy conductance
  !>
  !> In this routine we aggregate water limited canopy conductance: We calculate the variable
  !> on parent tiles as area weighted means from the child tile values.
  !>
  SUBROUTINE aggregate_canopy_cond_stressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Parameters for current block

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract   !< Aggregation method for area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: &
      & iblk, &    !< Current block ID
      & ics, &     !< Index of first cell in the block
      & ice        !< Index of last cell in the block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, canopy_cond_limited, weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  ! ======================================================================================================== !
  !>
  !>#### Update evaporation variables
  !>
  !> In this routine we calculate potential evaporation, evaporation and transpiration for the different
  !> tiles.
  !>
  !> The calculations depend on the turbulent mixing scheme: With tmx (explicit land-atmosphere coupling)
  !> we use exchange coefficients calculated in [[mo_turb_interface:update_exchange_coefficients]] for
  !> each leaf tile. The implicit coupling makes use of Richtmyer Morton coefficients provided from the
  !> atmosphere vertical diffusion (vdiff) turbulence scheme.
  !>
  SUBROUTINE update_evaporation(tile, options)

    USE mo_phy_schemes,            ONLY: q_effective, heat_transfer_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of the current block

    ! Local variables
    dsl4jsb_Def_config(SEB_)    !< Configurable variables of the "surface energy balance"
    dsl4jsb_Def_memory(A2L_)    !< Memory of "atmosphere to land" process
    dsl4jsb_Def_memory(HYDRO_)  !< Memory of "hydrology" process
    dsl4jsb_Def_memory(SEB_)    !< Memory of "surface energy balance" process
    dsl4jsb_Def_memory(TURB_)   !< Memory of "turbulence" process

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & q_acoef, &              !< Richtmyer-Morton A coefficient for specific humidity []
      & q_bcoef, &              !< Richtmyer-Morton B coefficient for specific humidity []
      & q_acoef_wtr, &          !< q_acoef for water surfaces []
      & q_bcoef_wtr, &          !< q_bcoef for water surfaces []
      & q_acoef_ice, &          !< q_acoef for ice surfaces []
      & q_bcoef_ice, &          !< q_bcoef for ice surfaces []
      & drag_srf, &             !< Surface drag [???]
      & drag_wtr, &             !< Drag for water surfaces [???]
      & drag_ice, &             !< Drag for ice surfaces [???]
      & ch, &                   !< Surface transfer coefficient for heat [???]
      & qsat_star, &            !< Surface saturation specific humidity []
      & qsat_lwtr, &            !< Saturation specific humidity over lake water []
      & qsat_lice, &            !< Saturation specific humidity over lake ice []
      & fact_qsat_srf, &        !< Factor for surface saturation specific energy []
      & fact_qsat_trans_srf, &  !< Factor to calculate transpiration from evaporation []
      & fact_q_air, &           !< Factor to calculate air moisture accounting for partially wet surfaces []
      & fract_lice, &           !< Lake ice fraction []
      & evapotrans_lnd, &       !< Evapotranspiration from land [kg m-2 s-1]
      & evapo_wtr, &            !< Evaporation from lake water [kg m-2 s-1]
      & evapo_ice, &            !< Evaporation from lake ice [kg m-2 s-1]
      & evapopot, &             !< Potential evaporation [kg m-2 s-1]
      & evapotrans, &           !< Potential evapotranspiration [kg m-2 s-1]
      & transpiration           !< Transpiration [kg m-2 s-1]

    REAL(wp) :: &
      & q_air, &                !< Specific humidity at lowest atmospheric level [kg kg-1]
      & q_air_eff, &            !< Effective specific humidity at lowest atmospheric level []
      & qsat_srf_eff, &         !< Effective surface saturation specific humidity []
      & heat_tcoef              !< Heat transfer coefficient (rho*C_h*|v|) [???]

    INTEGER :: &
      & iblk, &                 !< Current block ID
      & ics, &                  !< Index of first cell in block
      & ice, &                  !< Index of last cell in block
      & nc, &                   !< Number of cells in block
      & ic                      !< Looping index for grid cells
    REAL(wp) :: &
      & dtime, &                !< Time step length [s]
      & alpha                   !< Implicitness factor []

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_evaporation'

    ! Do nothing if hydrology calculations are not needed on this tile
    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    alpha   = options%alpha

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)

    ! Pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans) ! out

    ! --------------
    !  Lakes
    ! --------------
    ! As the HYDRO process is not calculated on the box tile, currently only the lake tile contains lakes.
    IF (tile%contains_lake) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,    qsat_lwtr)     ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  evapo_wtr)     ! out
      IF (model%config%use_tmx) THEN
        ! With tmx exchange coefficients are calculated on each leaf tile
        dsl4jsb_Get_var2D_onChunk(TURB_, ch)
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef)       ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef)       ! in
      ELSE
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef_wtr)   ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef_wtr)   ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,  drag_wtr)      ! in
      END IF

      ! Calculate (potential) evaporation from lakes
      !
      IF (model%config%use_tmx) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
        !$ACC   PRIVATE(q_air, heat_tcoef)
        DO ic=1,nc
          q_air = q_bcoef(ic)                                   ! Old moisture at lowest atmospheric level
          heat_tcoef = ch(ic)                                   ! Transfer coefficient for heat
                                                                ! TODO: distinguish wtr and ice?
          ! Evaporation from lake water corresponds to potential evaporation (no water limitation)
          evapo_wtr(ic) = heat_tcoef * (q_air - qsat_lwtr(ic))
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
        !$ACC   PRIVATE(q_air, heat_tcoef)
        DO ic=1,nc
          q_air = q_acoef_wtr(ic) * qsat_lwtr(ic) + q_bcoef_wtr(ic)      ! New moisture at lowest atmospheric
                                                                         ! level by back-substitution
          heat_tcoef = heat_transfer_coef(drag_wtr(ic), dtime, alpha)    ! Heat transfer coefficient
          ! Evaporation from lake water corresponds to potential evaporation (no water limitation)
          evapo_wtr (ic) = heat_tcoef * (q_air - qsat_lwtr(ic))
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN          ! Lake ice is activated in the model setup
        IF (.NOT. model%config%use_tmx) THEN
          dsl4jsb_Get_var2D_onChunk(A2L_,   q_acoef_ice)     ! in
          dsl4jsb_Get_var2D_onChunk(A2L_,   q_bcoef_ice)     ! in
          dsl4jsb_Get_var2D_onChunk(A2L_,   drag_ice)        ! in
        END IF
        dsl4jsb_Get_var2D_onChunk(SEB_,   qsat_lice)         ! in
        dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)        ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_ice)         ! out

        ! Calculation (potential) evaporation from lake ice
        !
        IF (model%config%use_tmx) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
          !$ACC   PRIVATE(q_air, heat_tcoef)
          DO ic=1,nc
            q_air = q_bcoef(ic)                                   ! Old moisture at lowest atmospheric level
            heat_tcoef = ch(ic)                                   ! Transfer coefficient; TODO: distinguish between wtr and ice?
            evapo_ice (ic) = heat_tcoef * (q_air - qsat_lice(ic)) ! Potential evaporation
            evapopot  (ic) = (1._wp - fract_lice(ic)) * evapo_wtr(ic) + fract_lice(ic) * evapo_ice(ic)
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
          !$ACC   PRIVATE(q_air, heat_tcoef)
          DO ic=1,nc
            q_air = q_acoef_ice(ic) * qsat_lice(ic) + q_bcoef_ice(ic) ! New moisture at lowest atmospheric level by back-substitution
            heat_tcoef = heat_transfer_coef(drag_ice(ic), dtime, alpha)    ! Transfer coefficient
            evapo_ice (ic) = heat_tcoef * (q_air - qsat_lice(ic))          ! Potential evaporation
            evapopot  (ic) = (1._wp - fract_lice(ic)) * evapo_wtr(ic) + fract_lice(ic) * evapo_ice(ic)
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      ELSE

        ! Without lake ice, lake potential evaporation equals evaporation from lake water.

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          evapopot(ic) = evapo_wtr(ic)
        END DO
        !$ACC END PARALLEL LOOP

      END IF

      ! There is no transpiration from lakes, thus evapotranspiration equals evaporation.

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        evapotrans(ic) = evapopot(ic)
      END DO
      !$ACC END PARALLEL LOOP

    ELSE IF (tile%contains_soil .OR. tile%contains_glacier) THEN
      ! --------------
      !  Land tiles
      ! --------------

      dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef)             ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef)             ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,  qsat_star)           ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_q_air)          ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_srf)       ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_trans_srf) ! in

      IF (model%config%use_tmx) THEN
        dsl4jsb_Get_var2D_onChunk(TURB_, ch)
      ELSE
        dsl4jsb_Get_var2D_onChunk(A2L_,  drag_srf)          ! in
      END IF

      dsl4jsb_Get_var2D_onChunk(HYDRO_,  evapotrans_lnd)    ! out
      IF (tile%contains_vegetation) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_, transpiration)    ! out
      END IF

      IF (model%config%use_tmx) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
        !$ACC   PRIVATE(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
        DO ic=1,nc
          q_air = q_bcoef(ic)  ! Old moisture at lowest atmospheric level
          heat_tcoef = ch(ic)  ! Transfer coefficient

          ! Compute effective air moisture and surface saturation humidity
          q_air_eff      = q_effective( 0._wp, q_air, 1._wp, 0._wp)
          qsat_srf_eff   = q_effective(qsat_star(ic), q_air, fact_qsat_srf(ic), fact_q_air(ic))
          evapotrans_lnd(ic) = heat_tcoef * (q_air_eff - qsat_srf_eff)  ! Evapotranspiration
          evapotrans(ic)     = evapotrans_lnd(ic)
          evapopot(ic)       = heat_tcoef * (q_air     - qsat_star(ic)) ! Potential evaporation
          IF (tile%contains_vegetation) THEN
            transpiration(ic) = fact_qsat_trans_srf(ic) * evapopot(ic)  ! Transpiration
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
        !$ACC   PRIVATE(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
        DO ic=1,nc
          q_air = q_acoef(ic) * qsat_star(ic) + q_bcoef(ic)              ! New moisture at lowest atmospheric
                                                                         ! level by back-substitution
          heat_tcoef = heat_transfer_coef(drag_srf(ic), dtime, alpha)    ! Transfer coefficient for heat

          ! Compute effective air moisture and surface saturation humidity
          q_air_eff      = q_effective( 0._wp, q_air, 1._wp, 0._wp)
          qsat_srf_eff   = q_effective(qsat_star(ic), q_air, fact_qsat_srf(ic), fact_q_air(ic))
          evapotrans_lnd(ic) = heat_tcoef * (q_air_eff - qsat_srf_eff)  ! Evapotranspiration
          evapotrans(ic)     = evapotrans_lnd(ic)
          evapopot(ic)       = heat_tcoef * (q_air     - qsat_star(ic)) ! Potential evaporation
          IF (tile%contains_vegetation) THEN
            transpiration(ic) = fact_qsat_trans_srf(ic) * evapopot(ic)  ! Transpiration
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_evaporation

  ! ======================================================================================================== !
  !>
  !>#### Aggregate evaporation variables
  !>
  !> In this routine we aggregate variables of task "evaporation": evaporation and transpiration
  !> variables are calculated on parent tiles as area weighted means from the child tile values.
  !>
  SUBROUTINE aggregate_evaporation(tile, options)

    TYPE(t_jsb_model), POINTER :: model                   !< This instance of jsbach
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options for the current block

    ! Local variables
    dsl4jsb_Def_config(SEB_)                   !< Configurable variables of process "surface energy balance"
    dsl4jsb_Def_memory(HYDRO_)                 !< Memory of the hydrology process

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract !< Aggregator to calculate area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_evaporation'

    INTEGER :: &
      & iblk, &                                !< Current block ID
      & ics, &                                 !< Index of first cell in the block
      & ice                                    !< Index of last cell in the block

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapopot,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, transpiration,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans_lnd,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_wtr,          weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_ice,        weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_evaporation

  ! ======================================================================================================== !
  !>
  !>#### Update snow and ice hydrology
  !>
  !> The routine is intended to update variables of the "snow_and_ice_hydrology" task. It is currently
  !> not used.
  !>
  !TODO: Remove routine
  SUBROUTINE update_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_ice_hydrology'

    IF (options%nc > 0) CONTINUE ! avoid compiler warnings about dummy arguments not being used

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%contains_lake) THEN
      !! Currently done in SEB process
    ELSE IF (tile%contains_land) THEN
!!$      CALL update_snow_and_ice_hydrology_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_snow_and_ice_hydrology

  ! ======================================================================================================== !
  !>
  !>#### Aggregate variables of task "snow_and_ice_hydrology"
  !>
  !> The routine is intended to aggregate variables of the "snow_and_ice_hydrology" task. It is currently
  !> not used.
  !>
  !TODO: Remove routine
  SUBROUTINE aggregate_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_ice_hydrology'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_ice_hydrology

  ! ======================================================================================================== !
  !>
  !>#### Update snow and wet fractions
  !>
  !> In this subroutine wet and frozen land surface fractions are calculated for the different tiles.
  !> For vegetation tiles, the main calculations are done in [[mo_hydro_process:calc_wet_fractions_veg]],
  !> while calculations for bare soil tiles happen in [[mo_hydro_process:calc_wet_fractions_bare]].
  !>
  SUBROUTINE update_snow_and_wet_fraction(tile, options)

    USE mo_hydro_process,   ONLY: calc_wskin_fractions_lice, calc_wet_fractions_veg, calc_wet_fractions_bare
    USE mo_phy_schemes,     ONLY: heat_transfer_coef

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile       !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options    !< Options for current block

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model                     !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER :: grid                      !< Horizontal grid
    REAL(wp), DIMENSION(options%nc) :: &
      & skinres_canopy_max, &       !< Maximum amount of water/snow on canopy
      & skinres_max, &              !< Maximum amount of surface water/ice
      & heat_tcoef                  !< Heat transfer coefficient

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_wet_fraction'

    INTEGER :: &
      & iblk, &                     !< Current block
      & ics, &                      !< Index of first cell in block
      & ice, &                      !< Index of last cell in block
      & nc, &                       !< Number of cells in block
      & ic, &                       !< Looping index for cells
      & pond_dynamics_scheme        !< Selected scheme for pond dynamics

    REAL(wp) :: &
      & dtime, &                    !< Time step length - exception in initialization time step [s]
      & alpha, &                    !< Implicitness factor [] (compare [[mo_jsb_interface:interface_full]])
      & config_w_skin_max           !< Maximum amount of water in skin reservoir [m water equivalent]

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(SEB_)        !< Configurable parameters of surface energy balance
    dsl4jsb_Def_memory(A2L_)        !< Memory of the atmosphere to land interface
    dsl4jsb_Def_config(HYDRO_)      !< Configurable parameters of the hydrology process
    dsl4jsb_Def_memory(SEB_)        !< Memory of surface energy balance process
    dsl4jsb_Def_memory(TURB_)       !< Memory of turbulence process
    dsl4jsb_Def_memory(HYDRO_)      !< Memory of hydrology process
#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_memory(PHENO_)      !< Memory of phenology process
#endif
#ifndef __NO_QUINCY__
    dsl4jsb_Def_memory(VEG_)
#endif

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: lai                 !< Leaf area index []
    dsl4jsb_Real2D_onChunk :: fract_snow          !< Snow fraction (not including snow on canopy) []
    dsl4jsb_Real2D_onChunk :: fract_skin          !< Wet skin fraction (soil and canopy; without ponds) []
    dsl4jsb_Real2D_onChunk :: fract_wet           !< Wet tile fraction (soil and canopy; incl. ponds) []
    dsl4jsb_Real2D_onChunk :: fract_pond          !< Inundated area fraction (area of temporary ponds) []
    dsl4jsb_Real2D_onChunk :: fract_pond_max      !< Maximum inundated area fraction []
    dsl4jsb_Real2D_onChunk :: fract_fpc_max       !< Maximum foliage projected cover []
    dsl4jsb_Real2D_onChunk :: fract_snow_can      !< Snow fraction on canopy []
    dsl4jsb_Real2D_onChunk :: fract_snow_soil     !< Snow fraction on soil []
    dsl4jsb_Real2D_onChunk :: fract_snow_lice     !< Snow fraction on lake ice []
    dsl4jsb_Real2D_onChunk :: weq_snow_can        !< Amount of snow on canopy [m water equivalent]
    dsl4jsb_Real2D_onChunk :: weq_snow_soil       !< Amount of snow on soil [m water equivalent]
    dsl4jsb_Real2D_onChunk :: weq_snow_lice       !< Amount of snow on lake ice [m water equivalent]
    dsl4jsb_Real2D_onChunk :: weq_pond            !< Amount of water/ice in pond reservoir [m water equiv.]
    dsl4jsb_Real2D_onChunk :: weq_pond_max        !< Maximum capacity of pond reservoir [m water equivalent]
    dsl4jsb_Real2D_onChunk :: wtr_skin            !< Liq. water in skin reservoir (soil and canopy) [m]
    dsl4jsb_Real2D_onChunk :: oro_stddev          !< Standard deviation of the topography [m]
    dsl4jsb_Real2D_onChunk :: drag_srf            !< Surface drag coefficient []
    dsl4jsb_Real2D_onChunk :: ch                  !< Surface heat transfer coefficient []
    dsl4jsb_Real2D_onChunk :: t                   !< Surface temperature [K]
    dsl4jsb_Real2D_onChunk :: press_srf           !< Surface pressure [Pa]
    dsl4jsb_Real2D_onChunk :: q_air               !< Specific humidity of the air []


    ! --------------------------------------------------------------------------------------------------- !

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    alpha   = options%alpha

    ! Do nothing if calculations are not needed on the current tile.
    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Set pointers to process configurations
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(HYDRO_)

    ! Use simple scalar for GPU
    config_w_skin_max = dsl4jsb_Config(HYDRO_)%w_skin_max       ! Maximum capacity of skin reservoir (soil + canopy)
    pond_dynamics_scheme = dsl4jsb_Config(HYDRO_)%pond_dynamics

    ! Set pointer to process memory
    dsl4jsb_Get_memory(HYDRO_)

    !> On the lake tile we only need to calculate the snow fraction on lake ice. (The lake ice fraction is
    !> calculated in [[mo_seb_lake:update_surface_energy_lake]].)
    IF (tile%is_lake) THEN
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_lice)    ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_lice)  ! out
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          CALL calc_wskin_fractions_lice( &
            & weq_snow_lice(ic),          & ! in
            & fract_snow_lice(ic)         & ! out
            & )
        END DO
        !$ACC END PARALLEL LOOP
      END IF
      ! Nothing else needs to be done on lake tiles
      RETURN
    END IF

    ! Set pointers to additional process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      SELECT CASE (model%config%model_scheme)
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

    ! Set pointers to specific variables needed from the process memories
    IF (model%config%use_tmx) THEN
      dsl4jsb_Get_memory(TURB_)
      dsl4jsb_Get_var2D_onChunk(TURB_, ch)            ! in
    ELSE
      dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)       ! in
    END IF
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)        ! in

    dsl4jsb_Get_var2D_onChunk(SEB_, t)                ! in

    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    oro_stddev)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_pond_max)    ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)          ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond)          ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond_max)      ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_skin)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_pond)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_wet)         ! out
    END IF
    IF (tile%contains_vegetation) THEN

      SELECT CASE (model%config%model_scheme)
#ifndef __QUINCY_STANDALONE__
      CASE (MODEL_JSBACH)
        dsl4jsb_Get_var2D_onChunk(PHENO_,   lai)              ! in
        dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max)    ! in
#endif
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        dsl4jsb_Get_var2D_onChunk(VEG_,     lai)              ! in
#endif
      END SELECT
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)    ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_soil)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)  ! out

#ifndef _OPENACC
    ! Test sanity of the surface temperature to prevent the unspecific "lookup table overflow" error message
    IF (ANY(t(:) < 50._wp .OR. t(:) > 400._wp)) THEN
      IF (ANY(t(:) < 50._wp)) THEN
        ic = MINLOC(t(:), DIM=1)
      ELSE
        ic = MaxLOC(t(:), DIM=1)
      END IF
      WRITE (message_text,*) 'Temperature out of bounds on tile ', tile%name, ' at ', '(', &
        &                    grid%lon(ic,iblk), ';', grid%lat(ic,iblk), '): t: ', t(ic)
      CALL finish(TRIM(routine), message_text)
    END IF
#endif

    !$ACC DATA CREATE(skinres_canopy_max, skinres_max, heat_tcoef) ASYNC(acc_stream)

    !
    ! Maximum capacity of the soil and canopy skin reservoir
    !
    ! tile: contains vegetation
    IF (tile%contains_vegetation) THEN
      SELECT CASE (model%config%model_scheme)
#ifndef __QUINCY_STANDALONE__
      CASE (MODEL_JSBACH)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1,nc
          skinres_canopy_max(ic) = config_w_skin_max * lai(ic) * fract_fpc_max(ic)
        END DO
        !$ACC END PARALLEL LOOP
#endif
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1,nc
          ! no scaling with 'fract_fpc(:)' because:
          ! max LAI fraction of the tile area is one (but not fract_fpc(:)) as vegetation spreads across the whole tile area
          skinres_canopy_max(ic) = config_w_skin_max * lai(ic)
        END DO
        !$ACC END PARALLEL LOOP
#endif
      END SELECT
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1,nc
        skinres_canopy_max(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! tile: contains soil and is glacier
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1,nc
      IF (tile%contains_soil) THEN
        skinres_max(ic) = config_w_skin_max + skinres_canopy_max(ic)
      ELSE
        skinres_max(ic) = 0._wp
      END IF

      !> Glaciers are entirely covered with snow. They do not have a wet fraction.
      !> Note: Although glaciers are snow covered, snow depth and the amount
      !> of snow on glaciers is zero.
      IF (tile%is_glacier) THEN
        fract_snow(ic)      = 1._wp
        fract_snow_soil(ic) = 1._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! Transfer coefficient
    IF (model%config%use_tmx) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        heat_tcoef(ic) = ch(ic)  ! TODO: distinguish between wtr and ice?
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        heat_tcoef(ic) = heat_transfer_coef(drag_srf(ic), dtime, alpha)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_can)       ! in
      CALL calc_wet_fractions_veg( &
        & dtime,                     & ! in
        & model%config%use_tmx,      & ! in
        & skinres_max(:),            & ! in
        & weq_pond_max(:),           & ! in
        & fract_pond_max(:),         & ! in
        & pond_dynamics_scheme,      & ! in
        & oro_stddev(:),             & ! in
        & t(:),                      & ! in, from the previous time step as long as this is called before the asselin filter
        & press_srf(:),              & ! in
        & heat_tcoef(:),             & ! in
        & q_air(:),                  & ! in
        & wtr_skin(:),               & ! in
        & weq_pond(:),               & ! in
        & weq_snow_soil(:),          & ! in
        & weq_snow_can(:),           & ! in
        & fract_snow_can(:),         & ! out
        & fract_skin(:),             & ! out
        & fract_pond(:),             & ! out
        & fract_wet(:),              & ! out
        & fract_snow_soil(:)         & ! out
        & )

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        fract_snow(ic) = fract_snow_soil(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (tile%is_bare) THEN
      CALL calc_wet_fractions_bare( &
        & dtime,                      & ! in
        & model%config%use_tmx,       & ! in
        & skinres_max,                & ! in
        & weq_pond_max(:),            & ! in
        & fract_pond_max(:),          & ! in
        & pond_dynamics_scheme,       & ! in
        & oro_stddev(:),              & ! in
        & t(:),                       & ! in
        & press_srf(:),               & ! in
        & heat_tcoef(:),              & ! in
        & q_air(:),                   & ! in
        & wtr_skin(:),                & ! in
        & weq_pond(:),                & ! in
        & weq_snow_soil(:),           & ! in
        & fract_skin(:),              & ! out
        & fract_pond(:),              & ! out
        & fract_wet(:),               & ! out
        & fract_snow_soil(:)          & ! out
        & )
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        fract_snow(ic) = fract_snow_soil(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC END DATA

  END SUBROUTINE update_snow_and_wet_fraction

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of snow and wet fractions
  !>
  !> Variables representing wet or frozen land surface fractions are aggregated here.
  !>
  SUBROUTINE aggregate_snow_and_wet_fraction(tile, options)

    TYPE(t_jsb_model), POINTER                :: model    !< This instance of ICON-Land
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of the current block

    dsl4jsb_Def_config(SEB_)                      !< Configurable parameters of the surface energy balance
    dsl4jsb_Def_memory(HYDRO_)                    !< Memory of the hydrology process

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract !< Aggregator to calculate area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_wet_fraction'

    INTEGER :: &
      & iblk, &                                   !< Number of current block
      & ics,  &                                   !< Index of first cell in the block
      & ice                                       !< Index of last cell in the block
    !------------------------------------------------------------------------------------------------------!

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_wet,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_skin,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_pond,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_soil,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_can,   weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_lice,  weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_wet_fraction

  ! ======================================================================================================== !
  !>
  !>#### Update and check the water balance
  !>
  !> Within this routine we check the water balance on each of the HYDRO process tiles.
  !>
  SUBROUTINE update_water_balance( tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o, rhoi
    USE mo_jsb_impl_constants,     ONLY: WB_LOGGING, WB_ERROR

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of the current block

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model                !< Current instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER    :: grid                 !< Horizontal grid

    dsl4jsb_Def_config(HYDRO_)                            !< Configurable parameters of the HYDRO process
    dsl4jsb_Def_memory(HYDRO_)                            !< Memory of the HYDRO process
    dsl4jsb_Def_memory(A2L_)                              !< Memory of the atmosphere to land interface

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & rain, &                      !< Rain (precipitation flux) [kg m-2 s-1]
      & snow, &                      !< Snow (precipitation flux) [kg m-2 s-1]
      & runoff, &                    !< (Surface) runoff [kg m-2 s-1]
      & drainage, &                  !< (Subsurface) drainage [kg m-2 s-1]
      & evapotrans, &                !< Evapotranspiration [kg m-2 s-1]
      & evapo_deficit, &             !< Evaporation deficit due to inconsistencies [m water eq. / (time step)]
      & wtr_skin, &                  !< (Liquid) water in the skin reservoir [m]
      & weq_snow, &                  !< Amount of snow [m water equivalent]
      & weq_pond, &                  !< Water/ice in pond reservoir [m water equivalent]
      & wtr_latflow_res_srf, &       !< Water content of intermediary reservoir representing lateral flow of
                                     !< surface runoff [m] (only with [[t_hydro_config:hydro_scale]]='Uniform')
      & weq_fluxes, &                !< Sum of all water fluxes [m3 s-1]
      & weq_land, &                  !< All liquid and frozen land water [m3]
      & weq_balance_err, &           !< Water imbalance within the time step [m3/(time step)]
      & weq_balance_err_count        !< Number of water balance errors since (re)start
    dsl4jsb_Real3D_onChunk :: &
      & wtr_soil_sl, &               !< Amount of liquid water in the soil layers [m]
      & ice_soil_sl, &               !< Amount of ice in the soil layers [m]
      & wtr_latflow_res_sl           !< Amount of water in lateral flow reservoir [m]

    REAL(wp), POINTER :: &
      & tile_fract(:), &             !< Fraction of the tile (rel. to the grid cell) []
      & area(:), &                   !< Grid cell area [m2]
      & lat(:), &                    !< Grid cell center latitude
      & lon(:)                       !< Grid cell center longitude

    ! Local variables
    REAL(wp) :: &
      & tile_area,    &              !< Area of the tile fraction [m2]
      & wb_threshold, &              !< Threshold for water balance errors
      & weq_land_new                 !< Total amount of land water, updated in this time step []

    REAL(wp), DIMENSION(options%nc) :: &
      & weq_soil, &                  !< Amount of water or ice in the soil [m water equivalent]
      & wtr_latflow_res              !< Amount of water in the reservoir representing lateral flow [m]

    INTEGER :: &
      & iblk, &                      !< Number of current nproma block
      & ics, &                       !< Index of first cell in current block
      & ice, &                       !< Index of last cell of current block
      & nc, &                        !< Number of cells in current block
      & ic, &                        !< Looping index for grid cells
      & is, &                        !< Looping index for soil layers
      & n_errors, &                  !< Counter for water balance errors
      & ie                           !< Looping index for errors
    INTEGER :: error_idx(options%nc) !< Vector listing grid cells with water balance issue

    REAL(wp) :: &
      & dtime, &                     !< Time step length [s]
      & ice2weq                      !< Factor to convert ice column height to m water equivalent []
    LOGICAL :: &
      & is_experiment_start, &       !< True: This is the first time step of an experiment
      & tile_contains_soil           !< True: This tile or child tiles have soil

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_balance'
    CHARACTER(len=5000)         :: message_text_long   !< Extra long messages

    !------------------------------------------------------------------------------------------------------!

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    is_experiment_start = is_time_experiment_start(options%current_datetime)

    ! Note: In other update routines we return here, if the process is not to be calculated on this tile
    !       (compare "IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN"). In this case however, the
    !       routine is also called from "aggregate_water_balance" - on the parent tile - to explicitly check
    !       the water balance on that tile.

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    tile_fract => tile%fract(ics:ice, iblk)
    area       => grid%area (ics:ice, iblk)
    lat        => grid%lat  (ics:ice, iblk)
    lon        => grid%lon  (ics:ice, iblk)

    tile_contains_soil = tile%contains_soil

    ! Scaling factor to convert ice column into water equivalent
    ice2weq = (rhoi / rhoh2o)

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      rain)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      snow)                ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    drainage)            ! in
    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_soil_sl)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_deficit)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_latflow_res_srf) ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_latflow_res_sl)  ! in
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_fluxes)            ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_land)              ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_balance_err)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_balance_err_count) ! inout

    n_errors = 0   ! initialize error counter

    !$ACC PARALLEL DEFAULT(PRESENT) CREATE(weq_soil, wtr_latflow_res) ASYNC(acc_stream)

    ! Calculation of the amount of liquid or frozen water at the surface or in the soil as well as in
    ! the reservoir representing lateral flow.
    IF (tile_contains_soil) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        weq_soil(ic)        = 0._wp
        wtr_latflow_res(ic) = wtr_latflow_res_srf(ic)
      END DO

      !$ACC LOOP SEQ
      DO is = 1, SIZE(wtr_soil_sl, DIM=2)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          weq_soil(ic)        = weq_soil(ic) + wtr_soil_sl(ic,is) + ice_soil_sl(ic,is) * ice2weq
          wtr_latflow_res(ic) = wtr_latflow_res(ic) + wtr_latflow_res_sl(ic,is)
        END DO
      END DO
    END IF

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tile_area, weq_land_new)
    DO ic=1,nc
      tile_area  = tile_fract(ic) * area(ic)

      ! Sum up all water fluxes (solid and liquid) of this time step
      weq_fluxes(ic) = rain(ic) + snow(ic) + evapotrans(ic) - runoff(ic) - drainage(ic)
      weq_fluxes(ic) = weq_fluxes(ic) * tile_area / rhoh2o               ! kg m-2 s-1 -> m3/s

      IF (tile_contains_soil) THEN
        weq_land_new = (wtr_skin(ic) + weq_snow(ic) + weq_pond(ic) &
          &            + weq_soil(ic) + wtr_latflow_res(ic))       &
          &          * tile_area          ! m -> m^3
      ELSE
        weq_land_new = 0._wp
      END IF

      !  Calculate the water balance error
      !> For each grid cell the current water budget should match the water budget of the previous time step
      !> plus the sum of the water fluxes within the time step. For tiles without soil, i.e. lake or glacier
      !> tiles, the water budget is zero, thus in and outgoing fluxes need to match.
      !  As weq_land, the water budget of the previous time step, is not available at the very first time
      !  step, no error is calculated right after mode initialization.
      IF (.NOT. is_experiment_start) THEN
        weq_balance_err(ic) = weq_land(ic) + weq_fluxes(ic) * dtime - weq_land_new
      END IF
      ! Define weq_land for the next time step
      weq_land(ic) = weq_land_new

#ifndef _OPENACC
      !> The threshold indicating a water balance error depends on model resolution and time step. It
      !> is defined such, that a violation in all grid cells would lead to a 10 cm sea level rise or fall
      !> after 1000 years.
      wb_threshold = 1.0e-11_wp * tile_area * dtime

      IF (ABS(weq_balance_err(ic)) > wb_threshold) THEN
        n_errors = n_errors + 1       ! Count number of water balance issues on this chunk
        error_idx(n_errors) = ic      ! Index of the respective grid cell (on chunk)
        weq_balance_err_count(ic) = weq_balance_err_count(ic) + 1.0_wp  ! Count water budget issues per cell
      END IF
#endif
    END DO
    !$ACC END PARALLEL

    !> In case of significant water balance problems detailed error messages can be printed to help debugging
    !> the code. Depending on parameter [[t_hydro_config:enforce_water_budget]] (namelist 'jsb_hydro_nml')
    !> the simulation is stopped after the first error message, detailed or short messages are printed,
    !> or the messages are suppressed.
#ifndef _OPENACC
    IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget /= WB_IGNORE) THEN
      DO ie = 1, n_errors
        ic = error_idx(ie)
        tile_area  = tile_fract(ic) * area(ic)

        WRITE (message_text_long,*) 'Water balance violation [m3 dt-1]',               NEW_LINE('a'), &
          & 'on ',TRIM(tile%name),' tile at', lat(ic),'N and ',lon(ic),'E',            NEW_LINE('a'), &
          & '(ic: ',ic,' iblk: ',iblk, ' tile_fract:',tile_fract(ic),'):',             NEW_LINE('a'), &
          & 'WB Error:           ', weq_balance_err(ic),                               NEW_LINE('a'), &
          & 'Rainfall:           ', rain(ic)       * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Snowfall:           ', snow(ic)       * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Evapotranspiration: ', evapotrans(ic) * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Runoff:             ', runoff(ic)     * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Drainage:           ', drainage(ic)   * tile_area / rhoh2o * dtime
        IF (tile%contains_soil) THEN
          WRITE (message_text_long,*) TRIM(message_text_long),                         NEW_LINE('a'), &
            & 'Skin reservoir:     ', wtr_skin(ic)           * tile_area,              NEW_LINE('a'), &
            & 'Snow reservoir:     ', weq_snow(ic)           * tile_area
          IF (.NOT. tile%is_lake) THEN
            WRITE (message_text_long,*) TRIM(message_text_long),                       NEW_LINE('a'), &
              & 'Pond reservoir:     ', weq_pond(ic)           * tile_area
          END IF
          WRITE (message_text_long,*) TRIM(message_text_long),                         NEW_LINE('a'), &
            & 'Soil water:         ', SUM(wtr_soil_sl(ic,:)) * tile_area,              NEW_LINE('a'), &
            & 'Soil ice weq:       ', SUM(ice_soil_sl(ic,:)) * ice2weq * tile_area,    NEW_LINE('a'), &
            & 'Lat. flow. res. srf: ', wtr_latflow_res_srf(ic) * tile_area,            NEW_LINE('a'), &
            & 'Lat. flow. res. blg: ', SUM(wtr_latflow_res_sl(ic,:)) * tile_area,      NEW_LINE('a'), &
            & 'Old reservoirs:     ', weq_land(ic),                                    NEW_LINE('a'), &
            & 'ET deficit:         ', evapo_deficit(ic) * tile_area / rhoh2o * dtime
        END IF
        IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget == WB_ERROR) THEN
          CALL finish (TRIM(routine), message_text_long)
        ELSE IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget == WB_LOGGING) THEN
          CALL message (TRIM(routine), message_text_long, all_print=.TRUE.)
        END IF
      END DO
    END IF
#endif

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_balance

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of the water balance
  !>
  !> Instead of aggregating the water balance from the child tiles, we explicitly compute it here also
  !> for the parent tile by calling [[update_water_balance]].
  !>
  SUBROUTINE aggregate_water_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of current block

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_balance'

    INTEGER :: iblk       ! Number of the current block

    iblk = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Explicitly compute the water balance on this (parent) tiles
    CALL update_water_balance(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_balance

#ifndef __QUINCY_STANDALONE__
  ! ======================================================================================================== !
  !>
  !>#### Global diagnostics of the hydrology process
  !>
  !> We calculate here global means or global sums of selected key variables of the hydrology process, that
  !> are intended to be used to monitor the performance of a running experiment. To have the global data
  !> available, the routine is called from [[mo_jsb_interface:jsbach_finish_timestep]], after the loop
  !> over the nproma blocks.
  !>
  SUBROUTINE global_hydrology_diagnostics(tile)

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)      !< Configurable parameters of the hydrology process
    dsl4jsb_Def_memory(HYDRO_)      !< Memory of the hydrology process

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_hydrology_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: transpiration    !< Transpiration [kg m-2 s-1]
    dsl4jsb_Real2D_onDomain :: evapotrans       !< Evapotranspiration [kg m-2 s-1]
    dsl4jsb_Real2D_onDomain :: weq_land         !< Total amount of liquid or frozen water of the land [m3]
    dsl4jsb_Real2D_onDomain :: discharge_ocean  !< Discharge to the ocean [m3 s-1]
    dsl4jsb_Real2D_onDomain :: wtr_rootzone_rel !< Relative root zone moisture []
    dsl4jsb_Real2D_onDomain :: fract_snow       !< Snow area fraction (not including snow on canopy) []
    dsl4jsb_Real2D_onDomain :: weq_snow         !< Snow amount on non-glacier land [m water equivalent]
    dsl4jsb_Real2D_onDomain :: weq_balance_err  !< Land water balance error [m3/(time step)]
    dsl4jsb_Real2D_onDomain :: weq_balance_err_count       !< Number of time steps with water balance error []

    LOGICAL, SAVE           :: print_wb_warning = .TRUE.   !< Print a warning in case of water balance issues

    REAL(wp), POINTER       :: trans_gmean(:)              !< Global mean transpiration on land [kg m-2 s-1]
    REAL(wp), POINTER       :: evapotrans_gmean(:)         !< Global mean evapotranspiration on land [kg m-2 s-1]
    REAL(wp), POINTER       :: weq_land_gsum(:)            !< Total amount of liquid or frozen land water [km3]
    REAL(wp), POINTER       :: discharge_ocean_gsum(:)     !< Global sum of discharge to the ocean [Sv]
    REAL(wp), POINTER       :: wtr_rootzone_rel_gmean(:)   !< Global mean relative root zone moisture []
    REAL(wp), POINTER       :: fract_snow_gsum(:)          !< Global snow area [Mio m2]
    REAL(wp), POINTER       :: weq_snow_gsum(:)            !< Global snow amount on non-glacier land [Gt]
    REAL(wp), POINTER       :: weq_balance_err_gsum(:)     !< Global land water balance error [m3/(time step)]

    TYPE(t_jsb_model), POINTER      :: model               !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER      :: grid                !< The ICON-Land horizontal grid

    REAL(wp), POINTER      :: area(:,:)                    !< Grid cell area [m2]
    REAL(wp), POINTER      :: notsea(:,:)                  !< Land mask (1: land incl. lakes; 0: ocean)
    LOGICAL,  POINTER      :: is_in_domain(:,:)            !< True: cell in domain; False: halo cell
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)              !< 1: cell in domain; 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)                !< Scaling factor: Land area of the cells
    REAL(wp)               :: global_land_area             !< Global land are [m3]


    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  transpiration)              ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  evapotrans)                 ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_land)                   ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  discharge_ocean)            ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  wtr_rootzone_rel)           ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  fract_snow)                 ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_snow)                   ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_balance_err)            ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_balance_err_count)      ! in

    trans_gmean          => HYDRO__mem%trans_gmean%ptr(:)           ! out
    evapotrans_gmean     => HYDRO__mem%evapotrans_gmean%ptr(:)      ! out
    weq_land_gsum        => HYDRO__mem%weq_land_gsum%ptr(:)         ! out
    discharge_ocean_gsum => HYDRO__mem%discharge_ocean_gsum%ptr(:)  ! out
    wtr_rootzone_rel_gmean => HYDRO__mem%wtr_rootzone_rel_gmean%ptr(:)  ! out
    fract_snow_gsum      => HYDRO__mem%fract_snow_gsum%ptr(:)       ! out
    weq_snow_gsum        => HYDRO__mem%weq_snow_gsum%ptr(:)         ! out
    weq_balance_err_gsum => HYDRO__mem%weq_balance_err_gsum%ptr(:)  ! out

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)     ! Fraction of the box tile

    dsl4jsb_Get_config(HYDRO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Domain Mask - to mask all halo cells for global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(grid%nproma,grid%nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ALLOCATE (scaling(grid%nproma,grid%nblks))


    ! Calculate global land variables, if requested for output
    ! -------------------------------

    IF (HYDRO__mem%trans_gmean%is_in_output .OR. HYDRO__mem%evapotrans_gmean%is_in_output &
      & .OR. HYDRO__mem%wtr_rootzone_rel_gmean%is_in_output) THEN

      ! Global land area
      global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))

      ! Weighting factor for each grid cell corresponding to its area fraction
      scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)

    END IF

    IF (HYDRO__mem%trans_gmean%is_in_output)           &
      &  trans_gmean          = global_sum_array(transpiration(:,:)   * scaling(:,:)) / global_land_area
    IF (HYDRO__mem%evapotrans_gmean%is_in_output)      &
      &  evapotrans_gmean     = global_sum_array(evapotrans(:,:)      * scaling(:,:)) / global_land_area
    ! Unit transformation from [m3] to [km3]: 1.e-9
    IF (HYDRO__mem%weq_land_gsum%is_in_output)         &
      &  weq_land_gsum        = global_sum_array(weq_land(:,:)  * notsea(:,:) * in_domain(:,:)) * 1.e-9_wp
    ! Unit transformation from [m3/s] to [Sv] (1 Sv = 1.e6 m3/s)
    IF (HYDRO__mem%discharge_ocean_gsum%is_in_output)  &
      &  discharge_ocean_gsum = global_sum_array(discharge_ocean(:,:) * in_domain(:,:)) * 1.e-6_wp
    IF (HYDRO__mem%wtr_rootzone_rel_gmean%is_in_output)  &
      &  wtr_rootzone_rel_gmean = global_sum_array(wtr_rootzone_rel(:,:) * scaling(:,:)) / global_land_area
    ! Unit transformation from [m2] to [Mio km2]: 1.e-12
    IF (HYDRO__mem%fract_snow_gsum%is_in_output)       &
      &  fract_snow_gsum      = global_sum_array(fract_snow(:,:)      * scaling(:,:)) * 1.e-12_wp
    ! Unit transformation from [m water equivalent](= [t]) to [Gt]: 1.e-9
    IF (HYDRO__mem%weq_snow_gsum%is_in_output)         &
      &  weq_snow_gsum        = global_sum_array(weq_snow(:,:)        * scaling(:,:)) * 1.e-9_wp
    ! No unit transformation [m3]: 1
    IF (HYDRO__mem%weq_balance_err_gsum%is_in_output)  &
      &  weq_balance_err_gsum = global_sum_array(weq_balance_err(:,:) * notsea(:,:) * in_domain(:,:))

    IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget == WB_IGNORE .AND. print_wb_warning) THEN
      IF (global_sum_array(weq_balance_err_count(:,:) * in_domain(:,:)) > 0.5_wp) THEN
        CALL message (TRIM(routine), 'Water balance issues detected. '// &
          & 'Consider rerun the simulation with enforce_water_budget = "logging" for more details.')
        ! Don't repeat this warning during this simulation period
        print_wb_warning = .FALSE.
      END IF
    END IF

    DEALLOCATE (scaling, in_domain)

  END SUBROUTINE global_hydrology_diagnostics
#endif

#endif
END MODULE mo_hydro_interface
