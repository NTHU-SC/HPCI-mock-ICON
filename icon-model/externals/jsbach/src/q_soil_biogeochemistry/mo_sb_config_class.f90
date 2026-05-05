!> QUINCY soil-biogeochemistry process config
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
!>#### define soil-biogeochemistry config structure, read soil-biogeochemistry namelist and init configuration parameters
!>
MODULE mo_sb_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters
  USE mo_jsb_config_class,    ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_config, get_number_of_sb_compartments, simple_1d_compartments
  PUBLIC :: SB_PART_DOM_IDX , SB_PART_DOM_ASSOC_IDX, SB_PART_SOLUBLE_LITTER_IDX,                   &
    &       SB_PART_POLYMERIC_LITTER_IDX, SB_PART_WOODY_LITTER_IDX
  PUBLIC :: SB_PART_MYCORRHIZA_IDX, SB_PART_MICROBIAL_IDX, SB_PART_RESIDUE_IDX, &
    &       SB_PART_RESIDUE_ASSOC_IDX
  PUBLIC :: CONST_DEP, REFYEAR_DEP, TRANS_DEP

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the soil_biogeochemistry process, derived from t_jsb_config
  !!
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_sb_config
    !> options
    !!
    CHARACTER(15)         :: sb_model_scheme
    CHARACTER(15)         :: sb_nloss_scheme
    CHARACTER(15)         :: bnf_scheme                 !< select scheme to simulate asymbiotic BNF (biological nitrogen fixation)
                                                        !! identical with the symbiotic N fixation in VEG_
                                                        !! none:      asymbiotic N fixation = 0
                                                        !! fixed:     identical with 'dynamic scheme' code
                                                        !! dynamic:   symbiotic N fixation as a dynamic trade-off of carbon and nitrogen opportunity costs;
                                                        !!              based on Rastaetter et al. 2001, Meyerholt et al. 2016, Kern, 2021
                                                        !! unlimited: asymbiotic N fixation set to satify plant growth demands at any point in time
    CHARACTER(15)         :: sb_adsorp_scheme
    !!
    INTEGER :: deposition_scheme !< one of: CONST_DEP=1, REFYEAR_DEP=2, TRANS_DEP=3
                                 !< - in case of REFYEAR_DEP and TRANS_DEP the data is read from an external file in IQ
    !> flags
    !!
    LOGICAL               :: read_sb_state          !< TRUE: call sb_read_states (mo_sb_init) after init/restart from file
    LOGICAL               :: flag_sb_prescribe_nh4  !< TRUE: solute NH4 is kept at a prescribed level
    LOGICAL               :: flag_sb_prescribe_no3  !< TRUE: solute NO3 is kept at a prescribed level
    LOGICAL               :: flag_sb_prescribe_po4  !< TRUE: solute PO4 is kept at a prescribed level
    LOGICAL               :: flag_sb_interactive_n  !< TRUE: N cycle state affects C,P processes via stoichiometry or allocation changes
    LOGICAL               :: flag_sb_interactive_p  !< TRUE: P cycle state affects C,N processes via stoichiometry or allocation changes
    ! mycorrhiza
    LOGICAL               :: flag_mycorrhiza
    LOGICAL               :: flag_mycorrhiza_org
    LOGICAL               :: flag_mycorrhiza_prim
    ! JSM
    LOGICAL               :: flag_sb_jsm_transport
    LOGICAL               :: flag_sb_jsm_litter_input
    LOGICAL               :: flag_sb_jsm_OM_sorption
    ! P sorption
    LOGICAL               :: flag_sb_double_langmuir
    !> parameter
    !!
    INTEGER               :: usda_taxonomy_class
    INTEGER               :: nwrb_taxonomy_class
    REAL(wp)              :: soil_ph                  !< soil pH - site specific from forcing info file
    REAL(wp)              :: soil_p_depth
    REAL(wp)              :: soil_p_labile
    REAL(wp)              :: soil_p_slow
    REAL(wp)              :: soil_p_occluded
    REAL(wp)              :: soil_p_primary
    REAL(wp)              :: qmax_org_fine_particle    !< maximum sorption capacity of OM to fine soil particle, mol C/ kg fine particle

    CHARACTER(len=filename_max) :: bc_quincy_soil_filename    !< IQ (QUINCY) soil input data
    CHARACTER(len=filename_max) :: deposition_filename_prefix !< IQ. Note: should NOT include an ending "_"

    ! parameter calibration
    LOGICAL          :: flag_apply_nml_parameters  !< replace values of sb_constants parameters by the values of the below nml_* variables
    REAL(wp)         :: nml_temp_freeze_thres_bio
    REAL(wp)         :: nml_temp_freeze_window_bio
    REAL(wp)         :: nml_frac_woody2dec_litter
    REAL(wp)         :: nml_frac_litter2fast_som
    REAL(wp)         :: nml_f_nit_noy
    REAL(wp)         :: nml_f_nit_n2o
    REAL(wp)         :: nml_vmax_nitrification
    REAL(wp)         :: nml_vmax_denitrification
    REAL(wp)         :: nml_t_ref_decomposition
    REAL(wp)         :: nml_ea_decomposition
    REAL(wp)         :: nml_ed_decomposition

   CONTAINS
     PROCEDURE :: Init => Init_sb_config
  END type t_sb_config

  !> INDs of soil biogeochemistry compartments (since sb compartments are static, these are also used as IDs)
  ENUM, BIND(C)
    ENUMERATOR ::                   &
      & SB_PART_DOM_IDX = 1,          &  !< dissolved organic matter
      & SB_PART_DOM_ASSOC_IDX,        &  !< minerally associated dissolved organic matter
      & SB_PART_SOLUBLE_LITTER_IDX,   &  !< metabolic litter (sugars, cellulose etc.)
      & SB_PART_POLYMERIC_LITTER_IDX, &  !< structural litter (lignified litter)
      & SB_PART_WOODY_LITTER_IDX,     &  !< physically protected wood litter
      & SB_PART_MYCORRHIZA_IDX,       &  !< mycorrhizal fungal biomass
      & SB_PART_MICROBIAL_IDX,        &  !< microbial biomass (bacteria!)
      & SB_PART_RESIDUE_IDX,          &  !< necromass of fungi, microbial and mycorrhiza
      & SB_PART_RESIDUE_ASSOC_IDX,    &  !< minerally associated microbial necromass
      & LAST_SB_PART_IDX                 !< last of this ENUM, used to determine the max number of compartments
  END ENUM

  !> deposition schemes (possible values for deposition_scheme)
  ENUM, BIND(C)
    ENUMERATOR ::       &
      & CONST_DEP = 1,      &  !< use global constant nhx, noy and p deposition
      & REFYEAR_DEP,        &  !< read monthly maps of nhx, noy and p deposition but only initially for reference year (IQ)
      & TRANS_DEP              !< read monthly maps of nhx, noy and p deposition for each simulation year (IQ)
  END ENUM

  INTEGER, PARAMETER :: simple_1d_compartments(6) = [ SB_PART_SOLUBLE_LITTER_IDX, SB_PART_POLYMERIC_LITTER_IDX, &
    & SB_PART_WOODY_LITTER_IDX, SB_PART_MYCORRHIZA_IDX, SB_PART_MICROBIAL_IDX, SB_PART_RESIDUE_IDX ]

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_sb_config
  !!
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_util_string,           ONLY: tolower
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_GENERIC

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_sb_config), INTENT(inout) :: config    !< config type for sb
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    LOGICAL                     :: lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename
    CHARACTER(len=filename_max) :: deposition_filename_prefix
    CHARACTER(15)               :: sb_model_scheme
    CHARACTER(15)               :: sb_nloss_scheme
    CHARACTER(15)               :: sb_bnf_scheme
    CHARACTER(15)               :: sb_adsorp_scheme
    CHARACTER(9)                :: deposition_scheme         !< one of constant, refyear, transient
    LOGICAL                     :: read_sb_state
    LOGICAL                     :: flag_sb_prescribe_nh4
    LOGICAL                     :: flag_sb_prescribe_no3
    LOGICAL                     :: flag_sb_prescribe_po4
    LOGICAL                     :: flag_sb_interactive_n
    LOGICAL                     :: flag_sb_interactive_p
    LOGICAL                     :: flag_mycorrhiza
    LOGICAL                     :: flag_mycorrhiza_org
    LOGICAL                     :: flag_mycorrhiza_prim
    LOGICAL                     :: flag_sb_jsm_transport
    LOGICAL                     :: flag_sb_jsm_litter_input
    LOGICAL                     :: flag_sb_jsm_OM_sorption
    LOGICAL                     :: flag_sb_double_langmuir
    INTEGER                     :: usda_taxonomy_class
    INTEGER                     :: nwrb_taxonomy_class
    REAL(wp)                    :: soil_ph
    REAL(wp)                    :: soil_p_depth,        &
                                   soil_p_labile,       &
                                   soil_p_slow,         &
                                   soil_p_occluded,     &
                                   soil_p_primary
    REAL(wp)                    :: qmax_org_fine_particle
    LOGICAL                     :: flag_apply_nml_parameters
    REAL(wp)                    :: nml_temp_freeze_thres_bio
    REAL(wp)                    :: nml_temp_freeze_window_bio
    REAL(wp)                    :: nml_frac_woody2dec_litter
    REAL(wp)                    :: nml_frac_litter2fast_som
    REAL(wp)                    :: nml_f_nit_noy
    REAL(wp)                    :: nml_f_nit_n2o
    REAL(wp)                    :: nml_vmax_nitrification
    REAL(wp)                    :: nml_vmax_denitrification
    REAL(wp)                    :: nml_t_ref_decomposition
    REAL(wp)                    :: nml_ea_decomposition
    REAL(wp)                    :: nml_ed_decomposition

    NAMELIST /lnd_sb_nml/                     &
                      active                    , &
                      lrestart_cont             , &
                      ic_filename               , &
                      bc_filename               , &
                      bc_quincy_soil_filename   , &
                      deposition_filename_prefix, &
                      sb_model_scheme           , &
                      sb_nloss_scheme           , &
                      sb_bnf_scheme             , &
                      sb_adsorp_scheme          , &
                      deposition_scheme         , &
                      read_sb_state             , &
                      flag_sb_prescribe_nh4     , &
                      flag_sb_prescribe_no3     , &
                      flag_sb_prescribe_po4     , &
                      flag_sb_interactive_n     , &
                      flag_sb_interactive_p     , &
                      flag_mycorrhiza           , &
                      flag_mycorrhiza_org       , &
                      flag_mycorrhiza_prim      , &
                      flag_sb_jsm_transport     , &
                      flag_sb_jsm_litter_input  , &
                      flag_sb_jsm_OM_sorption   , &
                      flag_sb_double_langmuir   , &
                      usda_taxonomy_class       , &
                      nwrb_taxonomy_class       , &
                      soil_ph                   , &
                      soil_p_depth              , &
                      soil_p_labile             , &
                      soil_p_slow               , &
                      soil_p_occluded           , &
                      soil_p_primary            , &
                      qmax_org_fine_particle, &
      &   flag_apply_nml_parameters, &
      &   nml_temp_freeze_thres_bio, &
      &   nml_temp_freeze_window_bio, &
      &   nml_frac_woody2dec_litter, &
      &   nml_frac_litter2fast_som, &
      &   nml_f_nit_noy, &
      &   nml_f_nit_n2o, &
      &   nml_vmax_nitrification, &
      &   nml_vmax_denitrification, &
      &   nml_t_ref_decomposition, &
      &   nml_ea_decomposition, &
      &   nml_ed_decomposition

    TYPE(t_jsb_vgrid), POINTER  :: time_months_vgrid
    INTEGER :: i_month
    INTEGER :: nml_handler, nml_unit, istat      ! variables for reading model-options from namelist
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_config'

    CALL message(TRIM(routine), 'Starting soil_biogeochemistry configuration')

    ! Set defaults
    active                      = .TRUE.
    lrestart_cont               = .FALSE.           ! TRUE: Continue even if SB variables are missing in restart file
    ic_filename                 = 'ic_land_sb.nc'
    bc_filename                 = 'bc_land_srf.nc'
    bc_quincy_soil_filename     = 'bc_quincy_soil.nc'
    deposition_filename_prefix  = 'bc_land_npdep_quincy'
    sb_model_scheme             = "simple_1d"       ! simple_1d jsm
    sb_nloss_scheme             = "dynamic"         ! none fixed dynamic
    sb_bnf_scheme               = "dynamic"         ! none dynamic unlimited (fixed == dynamic)
    sb_adsorp_scheme            = "eca_part"        ! eca_full; eca_part (only uptake, not for adsorption); eca_none
    deposition_scheme           = 'constant'
    read_sb_state               = .FALSE.
    flag_sb_prescribe_nh4       = .FALSE.
    flag_sb_prescribe_no3       = .FALSE.
    flag_sb_prescribe_po4       = .TRUE.
    flag_sb_interactive_n       = .TRUE.
    flag_sb_interactive_p       = .FALSE.
    flag_mycorrhiza             = .FALSE.
    flag_mycorrhiza_org         = .FALSE.
    flag_mycorrhiza_prim        = .FALSE.
    flag_sb_jsm_transport       = .FALSE.
    flag_sb_jsm_litter_input    = .FALSE.
    flag_sb_jsm_OM_sorption     = .FALSE.
    flag_sb_double_langmuir     = .FALSE.
    usda_taxonomy_class         = 1
    nwrb_taxonomy_class         = 1
    soil_ph                     = 6.5_wp
    soil_p_depth                = 0.5_wp        ! these soil_p* values are valid for DE-Hainich ?! and are used as default
    soil_p_labile               = 50._wp
    soil_p_slow                 = 70._wp
    soil_p_occluded             = 150._wp
    soil_p_primary              = 100._wp
    qmax_org_fine_particle      = 6.5537_wp
    flag_apply_nml_parameters   = .FALSE.
    nml_temp_freeze_thres_bio   = def_parameters
    nml_temp_freeze_window_bio  = def_parameters
    nml_frac_woody2dec_litter   = def_parameters
    nml_frac_litter2fast_som    = def_parameters
    nml_f_nit_noy               = def_parameters
    nml_f_nit_n2o               = def_parameters
    nml_vmax_nitrification      = def_parameters
    nml_vmax_denitrification    = def_parameters
    nml_t_ref_decomposition     = def_parameters
    nml_ea_decomposition        = def_parameters
    nml_ed_decomposition        = def_parameters

    ! read the namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('lnd_sb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_sb_nml)

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active                          = active
    config%lrestart_cont                   = lrestart_cont
    config%ic_filename                     = ic_filename
    config%bc_filename                     = bc_filename
    config%bc_quincy_soil_filename         = bc_quincy_soil_filename
    config%deposition_filename_prefix      = deposition_filename_prefix
    config%sb_model_scheme                 = TRIM(sb_model_scheme)
    config%sb_nloss_scheme                 = TRIM(sb_nloss_scheme)
    config%bnf_scheme                      = TRIM(sb_bnf_scheme)
    config%sb_adsorp_scheme                = TRIM(sb_adsorp_scheme)
    config%read_sb_state                   = read_sb_state
    config%flag_sb_prescribe_nh4           = flag_sb_prescribe_nh4
    config%flag_sb_prescribe_no3           = flag_sb_prescribe_no3
    config%flag_sb_prescribe_po4           = flag_sb_prescribe_po4
    config%flag_sb_interactive_n           = flag_sb_interactive_n
    config%flag_sb_interactive_p           = flag_sb_interactive_p
    config%flag_mycorrhiza                 = flag_mycorrhiza
    config%flag_mycorrhiza_org             = flag_mycorrhiza_org
    config%flag_mycorrhiza_prim            = flag_mycorrhiza_prim
    config%flag_sb_jsm_transport           = flag_sb_jsm_transport
    config%flag_sb_jsm_litter_input        = flag_sb_jsm_litter_input
    config%flag_sb_jsm_OM_sorption         = flag_sb_jsm_OM_sorption
    config%flag_sb_double_langmuir         = flag_sb_double_langmuir
    config%usda_taxonomy_class             = usda_taxonomy_class
    config%nwrb_taxonomy_class             = nwrb_taxonomy_class
    config%soil_ph                         = soil_ph
    config%soil_p_depth                    = soil_p_depth
    config%soil_p_labile                   = soil_p_labile
    config%soil_p_slow                     = soil_p_slow
    config%soil_p_occluded                 = soil_p_occluded
    config%soil_p_primary                  = soil_p_primary
    config%qmax_org_fine_particle          = qmax_org_fine_particle
    config%flag_apply_nml_parameters       = flag_apply_nml_parameters
    config%nml_temp_freeze_thres_bio       = nml_temp_freeze_thres_bio
    config%nml_temp_freeze_window_bio      = nml_temp_freeze_window_bio
    config%nml_frac_woody2dec_litter       = nml_frac_woody2dec_litter
    config%nml_frac_litter2fast_som        = nml_frac_litter2fast_som
    config%nml_f_nit_noy                   = nml_f_nit_noy
    config%nml_f_nit_n2o                   = nml_f_nit_n2o
    config%nml_vmax_nitrification          = nml_vmax_nitrification
    config%nml_vmax_denitrification        = nml_vmax_denitrification
    config%nml_t_ref_decomposition         = nml_t_ref_decomposition
    config%nml_ea_decomposition            = nml_ea_decomposition
    config%nml_ed_decomposition            = nml_ed_decomposition

    ! test sb_model_scheme
    SELECT CASE(TRIM(sb_model_scheme))
    CASE ("simple_1d")
      CALL message(TRIM(routine), 'Running the soil model: simple_1d (simple soil model)')
    CASE ("jsm")
      CALL message(TRIM(routine), 'Running the soil model: JSM (Jena Soil Model)')
    CASE DEFAULT
      CALL finish(TRIM(routine),  'No such soil model available: '//TRIM(sb_model_scheme))
    END SELECT

    ! Set deposition scheme
    SELECT CASE (tolower(TRIM(deposition_scheme)))
    CASE ("constant")
      CALL message(TRIM(routine), 'Running with constant deposition.')
      config%deposition_scheme = CONST_DEP
    CASE ("refyear")
      CALL message(TRIM(routine), 'Running with deposition from a refyear.')
      config%deposition_scheme = REFYEAR_DEP
    CASE ("transient")
      CALL message(TRIM(routine), 'Running with transient deposition.')
      config%deposition_scheme = TRANS_DEP
    CASE DEFAULT
      CALL finish(TRIM(routine), 'deposition_scheme == '//tolower(TRIM(deposition_scheme))//' not available.')
    END SELECT

    ! Create vertical axis
    ! - this grid might be useful in other contexts, too?
    time_months_vgrid => new_vgrid('time_months', ZAXIS_GENERIC, 12, &
      & levels=(/ (REAL(i_month, kind=wp), i_month = 1,12) /), units='')
    CALL register_vgrid(time_months_vgrid)

  END SUBROUTINE Init_sb_config

  ! ====================================================================================================== !
  !
  !> Returns the number of sb compartments defined in the enumerator
  !
  FUNCTION get_number_of_sb_compartments() RESULT(last_type_id)
    INTEGER :: last_type_id

    last_type_id = LAST_SB_PART_IDX - 1
  END FUNCTION get_number_of_sb_compartments

#endif
END MODULE mo_sb_config_class
