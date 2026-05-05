!> vegetation process config (QUINCY)
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
!>#### define vegetation config structure, read vegetation namelist and init configuration parameters
!>
MODULE mo_veg_config_class
#ifndef __NO_QUINCY__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters
  USE mo_jsb_config_class,    ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_config, Get_number_of_veg_compartments
  PUBLIC :: VEG_PART_LEAF_IDX , VEG_PART_FINE_ROOT_IDX, VEG_PART_COARSE_ROOT_IDX, VEG_PART_SAP_WOOD_IDX, &
    &       VEG_PART_HEART_WOOD_IDX, VEG_PART_LABILE_IDX, VEG_PART_RESERVE_IDX, VEG_PART_FRUIT_IDX

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the vegetation process, derived from t_jsb_config
  !!
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_veg_config
    INTEGER          :: pft_id                     !< ID of the PFT (at the moment [1,8] see lctlib file for details)
    CHARACTER(15)    :: bnf_scheme                 !< select scheme to simulate symbiotic BNF (biological nitrogen fixation)
                                                   !! identical with the asymbiotic N fixation in SB_
                                                   !! none:      symbiotic N fixation = 0
                                                   !! fixed:     symbiotic N fixation described from lctlib, but capped if plants become N saturated
                                                   !! dynamic:   symbiotic N fixation as a dynamic trade-off of carbon and nitrogen opportunity costs;
                                                   !!              based on Rastaetter et al. 2001, Meyerholt et al. 2016, Kern, 2021
                                                   !! unlimited: symbiotic N fixation set to satify plant growth demands at any point in time
    CHARACTER(15)    :: veg_dynamics_scheme        !< select scheme to calculate within-tile vegetation dynamics + start from bareground
                                                   !! none: constant mortality
                                                   !! "none population"
      ! Note: in previous versions of QS an additional "cohorts" scheme was available.
      !       SZ: "To enable a cohort like behaviour of the population mode in QS only (not possible in IQ),
      !       force the stand_replacing harvest to be triggered as a function of time since last disturbance
      !       (this will likely require to add a new memory variable to determine the tile age)."
      !       The initial planting (mo_veg_dynamics) may need to be adjusted to represent initial stand density.
      !       The function calc_veg_establishment needs to be disabled at all time-steps.
      !       And some parameters need to be adjusted (such as background_mort rates and stand density).

    CHARACTER(15)    :: biomass_alloc_scheme       !< select scheme for veg biomass allocation: fixed dynamic optimal
    CHARACTER(15)    :: leaf_stoichom_scheme       !< select scheme for C/N leaf stoichometry:  fixed dynamic optimal
    LOGICAL          :: read_veg_state             !< TRUE: call veg_read_states (mo_veg_init) after init/restart from file
    LOGICAL          :: flag_dynamic_roots         !< roots across soil layers: fixed after init or dynamic over time
                                                   !! for the models QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION fixed roots are used
                                                   !! by default; this is hard-coded in SPQ_ init
    LOGICAL          :: l_use_product_pools        !< enable product pools (required for crop/ and wood harvest processes); default: FALSE
    LOGICAL          :: flag_dynroots_h2o_n_limit  !< root growth across layers limited/affected by H2O and Nitrogen availability (needs flag_dynamic_roots="T")
    LOGICAL          :: flag_herbivory             !< whether herbivory (on grass PFT) is enabled
    LOGICAL          :: flag_log_negative_vegpool  !< debug option: test for negative veg_pool and message to LOG file for any negative value
    LOGICAL          :: flag_log_c_conservation    !< debug option: test for C conservation tests and message to LOG file if failing
    LOGICAL          :: flag_veg_interactive_n     !< TRUE: N cycle state affects C,P processes via stoichiometry or allocation changes
    LOGICAL          :: flag_veg_interactive_p     !< TRUE: P cycle state affects C,N processes via stoichiometry or allocation changes

    ! parameter calibration
    LOGICAL          :: flag_apply_nml_parameters  !< replace values of veg_constants parameters by the values of the below nml_* variables
    REAL(wp)         :: nml_fmaint_rate_base
    REAL(wp)         :: nml_fstore_sap_wood_max
    REAL(wp)         :: nml_k_f_demand
    REAL(wp)         :: nml_omega_nutrient_demand
    REAL(wp)         :: nml_background_mort_rate_tree
    REAL(wp)         :: nml_k_herbivory_grass
    REAL(wp)         :: nml_k_herbivory_pasture
    REAL(wp)         :: nml_sm2lm_grass
    REAL(wp)         :: nml_k1_root_alloc

   CONTAINS
     PROCEDURE :: Init => Init_veg_config
  END type t_veg_config

  !> INDs of vegetation compartments (since vegetation compartments are static, these are also used as IDs)
  ENUM, BIND(C)
    ENUMERATOR ::               &
      & VEG_PART_LEAF_IDX = 1,    &   !< leaves
      & VEG_PART_FINE_ROOT_IDX,   &   !< fine roots
      & VEG_PART_COARSE_ROOT_IDX, &   !< coarse roots
      & VEG_PART_SAP_WOOD_IDX,    &   !< sap wood
      & VEG_PART_HEART_WOOD_IDX,  &   !< heart wood
      & VEG_PART_LABILE_IDX,      &   !< labile
      & VEG_PART_RESERVE_IDX,     &   !< reserve
      & VEG_PART_FRUIT_IDX,       &   !< fruit
      & LAST_VEG_PART_IDX ! needs to be the last -- it is used to determine the number of compartments
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_config_class'

CONTAINS
  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_veg_config
  !!
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_veg_config), INTENT(inout) :: config    !< config type for veg
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namelist, identical to variable-name in namelist
    LOGICAL                     :: active
    LOGICAL                     :: lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename,   &
      &                            bc_filename
    INTEGER                     :: plant_functional_type_id
    CHARACTER(15)               :: veg_bnf_scheme
    CHARACTER(15)               :: veg_dynamics_scheme
    CHARACTER(15)               :: biomass_alloc_scheme
    CHARACTER(15)               :: leaf_stoichom_scheme
    LOGICAL                     :: read_veg_state
    LOGICAL                     :: flag_dynamic_roots
    LOGICAL                     :: flag_dynroots_h2o_n_limit
    LOGICAL                     :: flag_herbivory
    LOGICAL                     :: flag_log_negative_vegpool
    LOGICAL                     :: flag_log_c_conservation
    LOGICAL                     :: flag_veg_interactive_n
    LOGICAL                     :: flag_veg_interactive_p
    LOGICAL                     :: l_use_product_pools
    LOGICAL                     :: flag_apply_nml_parameters
    REAL(wp)                    :: nml_fmaint_rate_base
    REAL(wp)                    :: nml_fstore_sap_wood_max
    REAL(wp)                    :: nml_k_f_demand
    REAL(wp)                    :: nml_omega_nutrient_demand
    REAL(wp)                    :: nml_background_mort_rate_tree
    REAL(wp)                    :: nml_k_herbivory_grass
    REAL(wp)                    :: nml_k_herbivory_pasture
    REAL(wp)                    :: nml_sm2lm_grass
    REAL(wp)                    :: nml_k1_root_alloc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    NAMELIST /lnd_veg_nml/          &
      active,                       &
      lrestart_cont,                &
      ic_filename,                  &
      bc_filename,                  &
      plant_functional_type_id,     &
      veg_bnf_scheme,               &
      veg_dynamics_scheme,          &
      biomass_alloc_scheme,         &
      leaf_stoichom_scheme,         &
      read_veg_state,               &
      flag_dynamic_roots,           &
      flag_dynroots_h2o_n_limit,    &
      flag_herbivory,               &
      flag_log_negative_vegpool,    &
      flag_log_c_conservation,      &
      flag_veg_interactive_n,       &
      flag_veg_interactive_p,       &
      l_use_product_pools,          &
      &   flag_apply_nml_parameters, &
      &   nml_fmaint_rate_base, &
      &   nml_fstore_sap_wood_max, &
      &   nml_k_f_demand, &
      &   nml_omega_nutrient_demand, &
      &   nml_background_mort_rate_tree, &
      &   nml_k_herbivory_grass, &
      &   nml_k_herbivory_pasture, &
      &   nml_sm2lm_grass, &
      &   nml_k1_root_alloc

    INTEGER  :: nml_handler, nml_unit, istat      ! variables for reading model-options from namelist
    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    active                    = .TRUE.
    lrestart_cont             = .FALSE.  ! TRUE: Continue although VEG variables are missing in restart file
    bc_filename               = 'bc_land_phys.nc'
    ic_filename               = 'ic_land_veg.nc'
    plant_functional_type_id  = 1
    veg_bnf_scheme            = "dynamic"
    veg_dynamics_scheme       = "population"
    biomass_alloc_scheme      = "dynamic"
    leaf_stoichom_scheme      = "dynamic"
    read_veg_state            = .FALSE.
    flag_dynamic_roots        = .TRUE.
    flag_dynroots_h2o_n_limit = .FALSE.
    flag_herbivory            = .FALSE.
    flag_log_negative_vegpool = .FALSE.
    flag_log_c_conservation   = .FALSE.
    flag_veg_interactive_n    = .TRUE.
    flag_veg_interactive_p    = .FALSE.
    l_use_product_pools       = .FALSE.
    flag_apply_nml_parameters      = .FALSE.
    nml_fmaint_rate_base           = def_parameters
    nml_fstore_sap_wood_max        = def_parameters
    nml_k_f_demand                 = def_parameters
    nml_omega_nutrient_demand      = def_parameters
    nml_background_mort_rate_tree  = def_parameters
    nml_k_herbivory_grass          = def_parameters
    nml_k_herbivory_pasture        = def_parameters
    nml_sm2lm_grass                = def_parameters
    nml_k1_root_alloc              = def_parameters

    ! read the namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('lnd_veg_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_veg_nml)

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active                     = active
    config%lrestart_cont              = lrestart_cont
    config%ic_filename                = ic_filename
    config%bc_filename                = bc_filename
    config%pft_id                     = plant_functional_type_id
    config%bnf_scheme                 = TRIM(veg_bnf_scheme)
    config%veg_dynamics_scheme        = TRIM(veg_dynamics_scheme)
    config%biomass_alloc_scheme       = TRIM(biomass_alloc_scheme)
    config%leaf_stoichom_scheme       = TRIM(leaf_stoichom_scheme)
    config%read_veg_state             = read_veg_state
    config%flag_dynamic_roots         = flag_dynamic_roots
    config%flag_dynroots_h2o_n_limit  = flag_dynroots_h2o_n_limit
    config%flag_herbivory             = flag_herbivory
    config%flag_log_negative_vegpool  = flag_log_negative_vegpool
    config%flag_log_c_conservation    = flag_log_c_conservation
    config%flag_veg_interactive_n     = flag_veg_interactive_n
    config%flag_veg_interactive_p     = flag_veg_interactive_p
    config%l_use_product_pools        = l_use_product_pools
    config%flag_apply_nml_parameters      = flag_apply_nml_parameters
    config%nml_fmaint_rate_base           = nml_fmaint_rate_base
    config%nml_fstore_sap_wood_max        = nml_fstore_sap_wood_max
    config%nml_k_f_demand                 = nml_k_f_demand
    config%nml_omega_nutrient_demand      = nml_omega_nutrient_demand
    config%nml_background_mort_rate_tree  = nml_background_mort_rate_tree
    config%nml_k_herbivory_grass          = nml_k_herbivory_grass
    config%nml_k_herbivory_pasture        = nml_k_herbivory_pasture
    config%nml_sm2lm_grass                = nml_sm2lm_grass
    config%nml_k1_root_alloc              = nml_k1_root_alloc

  END SUBROUTINE Init_veg_config

  ! ====================================================================================================== !
  !>
  !> Returns the number of veg compartments defined in the enumerator
  !>
  FUNCTION Get_number_of_veg_compartments() RESULT(last_type_id)
    INTEGER :: last_type_id

    last_type_id = LAST_VEG_PART_IDX - 1
  END FUNCTION Get_number_of_veg_compartments

#endif
END MODULE mo_veg_config_class
