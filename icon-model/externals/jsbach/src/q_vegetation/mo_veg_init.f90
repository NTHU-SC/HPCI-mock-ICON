!> vegetation variables init (QUINCY)
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
!>#### initialization of vegetation memory variables using, e.g., ic & bc input files
!>
MODULE mo_veg_init
#ifndef __NO_QUINCY__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_jsb_math_constants,      ONLY: eps8
  USE mo_jsb_impl_constants,      ONLY: def_parameters

  USE mo_jsb_process_class,       ONLY: VEG_, Q_AGR_, Q_SYL_, ALCC_
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_class,               ONLY: get_model

  USE mo_jsb_grid_class,          ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,                ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_SEED_BED_POOL_ID

  dsl4jsb_Use_config(VEG_)

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: veg_init
#ifndef __QUINCY_STANDALONE__
  PUBLIC :: veg_read_states
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  ! ======================================================================================================= !
  !>Intialize vegetation process (after memory has been set up)
  !>
  !>  differtiates between the various QUINCY models: land, plant, soil, canopy, test_canopy
  !>
  SUBROUTINE veg_init(tile)
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_quincy_model_config,             ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_process_class,               ONLY: HYDRO_
    USE mo_q_assimi_parameters,             ONLY: CiCa_default_C3, CiCa_default_C4
    USE mo_q_assimi_constants,              ONLY: ic3phot, ic4phot
    USE mo_atmland_constants,               ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
      &                                           standard_press_srf, def_co2_deltaC13, def_co2_deltaC14
    USE mo_jsb_physical_constants,          ONLY: molar_mass_N
    USE mo_isotope_util,                    ONLY: calc_fractionation,calc_mixing_ratio_N15N14,calc_mixing_ratio_C14C
    USE mo_q_assimi_process,                ONLY: discrimination_ps_c13, discrimination_ps_c14
    USE mo_veg_constants,                   ONLY: root2leaf_cn, root2leaf_np, ITREE, &
      &   background_mort_rate_grass, k1_root_alloc, rfr_ratio_toc, fresp_growth, &
      &   fmaint_rate_base, fmaint_rate_wood, fstore_sap_wood_max, k_f_demand, omega_nutrient_demand, &
      &   background_mort_rate_tree, k_herbivory_grass, k_herbivory_pasture, sm2lm_grass
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),              POINTER :: model                !< the model
    TYPE(t_lnd_bgcm_store),         POINTER :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_jsb_grid),               POINTER :: hgrid                !< Horizontal grid
    TYPE(t_jsb_vgrid),              POINTER :: vgrid_soil_w         !< Vertical grid
    TYPE(t_lctlib_element),         POINTER :: lctlib               !< land-cover-type library - parameter across pft's
    REAL(wp)                                :: hlp1                 !< helper variable
    REAL(wp)                                :: ref_lai              !< LAI                                 | for canopy model init
    REAL(wp)                                :: ref_nleaf            !< N leaf per unit leaf area in gN/m2  | for canopy model init
    INTEGER                                 :: iblk, ic, is         !< loop dimensions
    INTEGER                                 :: nproma, nblks        !< dimensions
    INTEGER                                 :: nsoil_w              !< number of soil layers (water)
    INTEGER                                 :: bgcm_idx             !< index of the veg/seed bed pool within the bgcm store
    REAL(wp), POINTER           :: veg_pool_mt_domain(:,:,:,:)      !< dim: compartments, elements, nc, nblks -- on domain!
    REAL(wp), POINTER           :: seed_bed_pool_mt_domain(:,:,:)   !< dim: elements, nc, nblks -- on domain!
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !

    ! HYDRO_ 2D
    dsl4jsb_Real2D_onDomain    :: num_sl_above_bedrock
    dsl4jsb_Real2D_onDomain    :: root_depth
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    ! VEG_ 2D
    dsl4jsb_Real2D_onDomain      :: target_cn_leaf
    dsl4jsb_Real2D_onDomain      :: target_np_leaf
    dsl4jsb_Real2D_onDomain      :: target_cn_fine_root
    dsl4jsb_Real2D_onDomain      :: target_np_fine_root
    dsl4jsb_Real2D_onDomain      :: dens_ind
    dsl4jsb_Real2D_onDomain      :: height
    dsl4jsb_Real2D_onDomain      :: lai
    dsl4jsb_Real2D_onDomain      :: dphi
    dsl4jsb_Real2D_onDomain      :: npp_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: unit_uptake_p_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: growth_np_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: growth_p_limit_based_on_n_mavg
    dsl4jsb_Real2D_onDomain      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: uptake_n_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: growth_cn_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: growth_np_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: npp_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: beta_sinklim_ps
    dsl4jsb_Real2D_onDomain      :: leaf_cn_direction
    dsl4jsb_Real2D_onDomain      :: t_air_daytime
    dsl4jsb_Real2D_onDomain      :: press_srf_daytime
    dsl4jsb_Real2D_onDomain      :: co2_mixing_ratio_daytime
    dsl4jsb_Real2D_onDomain      :: ga_daytime
    dsl4jsb_Real2D_onDomain      :: beta_sinklim_ps_daytime
    dsl4jsb_Real2D_onDomain      :: t_jmax_opt_daytime
    dsl4jsb_Real2D_onDomain      :: t_jmax_opt_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: t_air_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: press_srf_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: co2_mixing_ratio_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: ga_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_sinklim_ps_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onDomain      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onDomain      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onDomain      :: t_soil_srf_tphen_mavg
    dsl4jsb_Real2D_onDomain      :: npp_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: ga_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onDomain      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain      :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain      :: lai_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain      :: beta_sinklim_ps_tacclim_mavg
    dsl4jsb_Real2D_onDomain      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onDomain      :: t_soil_root_tacclim_mavg
    dsl4jsb_Real2D_onDomain      :: t_air_week_mavg
    dsl4jsb_Real2D_onDomain      :: t_air_month_mavg
    dsl4jsb_Real2D_onDomain      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onDomain      :: wtr_rootzone_rel_talloc_mavg
    dsl4jsb_Real2D_onDomain      :: exudation_c_tmyc_mavg
    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc
    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc_tvegdyn_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain    :: root_fraction_sl
    dsl4jsb_Real3D_onDomain    :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain    :: fleaf_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onDomain    :: fn_oth_cl
    dsl4jsb_Real3D_onDomain    :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onDomain    :: lai_cl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    hgrid  => Get_grid(model%grid_id)
    nblks  =  hgrid%nblks
    nproma =  hgrid%nproma
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       = vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    bgcm_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_POOL_ID, tile%name, routine)
    veg_pool_mt_domain => bgcm_store%store_2l_2d_bgcms(bgcm_store%idx_in_store(bgcm_idx))%mt_2l_2d_bgcm(:,:,:,:)
    bgcm_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_POOL_ID, tile%name, routine)
    seed_bed_pool_mt_domain => bgcm_store%store_1l_2d_bgcms(bgcm_store%idx_in_store(bgcm_idx))%mt_1l_2d_bgcm(:,:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   num_sl_above_bedrock)        ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   root_depth)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_ubound_sl)
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onDomain(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onDomain(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onDomain(VEG_, target_cn_fine_root)
    dsl4jsb_Get_var2D_onDomain(VEG_, target_np_fine_root)
    dsl4jsb_Get_var2D_onDomain(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onDomain(VEG_, height)
    dsl4jsb_Get_var2D_onDomain(VEG_, lai)
    dsl4jsb_Get_var2D_onDomain(VEG_, dphi)
    dsl4jsb_Get_var2D_onDomain(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_uptake_n_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_uptake_p_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_np_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_p_limit_based_on_n_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, uptake_n_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_cn_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_np_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, npp_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps)
    dsl4jsb_Get_var2D_onDomain(VEG_, leaf_cn_direction)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_soil_srf_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, npp_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, an_boc_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, net_growth_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, lai_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_soil_root_tacclim_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_tacclim_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_week_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_month_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, wtr_rootzone_rel_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, exudation_c_tmyc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, rfr_ratio_boc)
    dsl4jsb_Get_var2D_onDomain(VEG_, rfr_ratio_boc_tvegdyn_mavg)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onDomain(VEG_, root_fraction_sl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fleaf_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fleaf_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, leaf_nitrogen_cl)
    SELECT CASE(model%config%qmodel_id)
    CASE(QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION)
      dsl4jsb_Get_var3D_onDomain(VEG_, lai_cl)
    END SELECT
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.8 apply calibration of vegetation constants if enabled
    !>
    IF (dsl4jsb_Config(VEG_)%flag_apply_nml_parameters) THEN
      ! fmaint_rate_base
      IF (dsl4jsb_Config(VEG_)%nml_fmaint_rate_base > def_parameters) THEN
        fmaint_rate_base = dsl4jsb_Config(VEG_)%nml_fmaint_rate_base
        fmaint_rate_wood = 0.25_wp * fmaint_rate_base
        !$ACC UPDATE DEVICE(fmaint_rate_base, fmaint_rate_wood) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: fmaint_rate_base & fmaint_rate_wood')
      END IF
      ! fstore_sap_wood_max
      IF (dsl4jsb_Config(VEG_)%nml_fstore_sap_wood_max > def_parameters) THEN
        fstore_sap_wood_max = dsl4jsb_Config(VEG_)%nml_fstore_sap_wood_max
        !$ACC UPDATE DEVICE(fstore_sap_wood_max) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: fstore_sap_wood_max')
      END IF
      ! k_f_demand
      IF (dsl4jsb_Config(VEG_)%nml_k_f_demand > def_parameters) THEN
        k_f_demand = dsl4jsb_Config(VEG_)%nml_k_f_demand
        !$ACC UPDATE DEVICE(k_f_demand) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: k_f_demand')
      END IF
      ! omega_nutrient_demand
      IF (dsl4jsb_Config(VEG_)%nml_omega_nutrient_demand > def_parameters) THEN
        omega_nutrient_demand = dsl4jsb_Config(VEG_)%nml_omega_nutrient_demand
        !$ACC UPDATE DEVICE(omega_nutrient_demand) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: omega_nutrient_demand')
      END IF
      ! background_mort_rate_tree
      IF (dsl4jsb_Config(VEG_)%nml_background_mort_rate_tree > def_parameters) THEN
        background_mort_rate_tree = dsl4jsb_Config(VEG_)%nml_background_mort_rate_tree
        !$ACC UPDATE DEVICE(background_mort_rate_tree) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: background_mort_rate_tree')
      END IF
      ! k_herbivory_grass
      IF (dsl4jsb_Config(VEG_)%nml_k_herbivory_grass > def_parameters) THEN
        k_herbivory_grass = dsl4jsb_Config(VEG_)%nml_k_herbivory_grass
        !$ACC UPDATE DEVICE(k_herbivory_grass) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: k_herbivory_grass')
      END IF
      ! k_herbivory_pasture
      IF (dsl4jsb_Config(VEG_)%nml_k_herbivory_pasture > def_parameters) THEN
        k_herbivory_pasture = dsl4jsb_Config(VEG_)%nml_k_herbivory_pasture
        !$ACC UPDATE DEVICE(k_herbivory_pasture) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: k_herbivory_pasture')
      END IF
      ! sm2lm_grass
      IF (dsl4jsb_Config(VEG_)%nml_sm2lm_grass > def_parameters) THEN
        sm2lm_grass = dsl4jsb_Config(VEG_)%nml_sm2lm_grass
        !$ACC UPDATE DEVICE(sm2lm_grass) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: sm2lm_grass')
      END IF
      ! k1_root_alloc
      IF (dsl4jsb_Config(VEG_)%nml_k1_root_alloc > def_parameters) THEN
        k1_root_alloc = dsl4jsb_Config(VEG_)%nml_k1_root_alloc
        !$ACC UPDATE DEVICE(k1_root_alloc) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: k1_root_alloc')
      END IF
    END IF

    !>0.9 init veg_pool_mt_domain with zero
    !>
    veg_pool_mt_domain(:,:,:,:) = 0.0_wp
    seed_bed_pool_mt_domain(:,:,:) = 0.0_wp

    !>1.0 vegetation init specific to the SLM model
    !>
    ! QUINCY models: land, plant,soil (case 1) and canopy (case 2)
    SELECT CASE(model%config%qmodel_id)
    !>  1.1 Land & Plant & Soil
    !>
    CASE(QLAND, QPLANT, QSOIL)
      target_cn_leaf(:,:)      = lctlib%cn_leaf
      target_np_leaf(:,:)      = lctlib%np_leaf
      target_cn_fine_root(:,:) = target_cn_leaf(:,:) / root2leaf_cn
      target_np_fine_root(:,:) = target_np_leaf(:,:) / root2leaf_np
      ! differ between vegetation dynamics schemes
      SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
      ! start bareground with vegetation dynamics switched on: population scheme
      CASE ("population")
        seed_bed_pool_mt_domain(ixC, :, :) = lctlib%seed_size * 10.0_wp
        dens_ind(:,:)                       = 0.0_wp
        seed_bed_pool_mt_domain(ixN, :, :) &
          & = seed_bed_pool_mt_domain(ixC, :, :) / (1._wp + fresp_growth) / lctlib%cn_leaf
        seed_bed_pool_mt_domain(ixP, :, :) = seed_bed_pool_mt_domain(ixN, :, :) / lctlib%np_leaf
        ! pool%seed_bed%carbon13
        IF (lctlib%ps_pathway == ic3phot) THEN
          hlp1 = discrimination_ps_c13(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ELSE
          hlp1 = discrimination_ps_c13(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ENDIF
        DO ic = 1,nproma
          DO iblk = 1,nblks
            seed_bed_pool_mt_domain(ixC13, ic, iblk) = seed_bed_pool_mt_domain( ixC, ic, iblk) &
              &                                        * calc_fractionation(def_co2_mixing_ratio,     &
              &                                                             def_co2_mixing_ratio_C13, hlp1)
          ENDDO
        ENDDO
        ! pool%seed_bed%carbon14
        IF (lctlib%ps_pathway == ic3phot) THEN
          hlp1 = discrimination_ps_c14(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ELSE
          hlp1 = discrimination_ps_c14(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ENDIF
        DO ic = 1,nproma
          DO iblk = 1,nblks
            seed_bed_pool_mt_domain(ixC14, ic, iblk)  = calc_mixing_ratio_C14C(def_co2_deltaC13, &
              &                                                             def_co2_deltaC14 - hlp1 )    &
              &                                         * seed_bed_pool_mt_domain(ixC, ic, iblk)
            seed_bed_pool_mt_domain(ixN15, ic, iblk)  = seed_bed_pool_mt_domain(ixN, ic, iblk) &
              &                                         / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ENDDO
        ENDDO
      ! no vegetation dynamics and not starting with bare ground
      CASE ("none")
        veg_pool_mt_domain(ix_reserve, ixC, :, :) = 5.0_wp
        dens_ind(:,:)                      = 0.4_wp
        veg_pool_mt_domain(ix_reserve, ixN, :, :) &
          & = veg_pool_mt_domain(ix_reserve, ixC, :, :) / (1._wp + fresp_growth) / lctlib%cn_leaf
        veg_pool_mt_domain(ix_reserve, ixP, :, :) = veg_pool_mt_domain(ix_reserve, ixN, :, :) / lctlib%np_leaf
        ! pool%reserve%carbon13
        IF (lctlib%ps_pathway == ic3phot) THEN
          hlp1 = discrimination_ps_c13(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ELSE
          hlp1 = discrimination_ps_c13(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ENDIF
        DO ic = 1,nproma
          DO iblk = 1,nblks
            veg_pool_mt_domain(ix_reserve, ixC13, ic, iblk) = veg_pool_mt_domain(ix_reserve, ixC, ic, iblk) &
              &                                               * calc_fractionation(def_co2_mixing_ratio,    &
              &                                                         def_co2_mixing_ratio_C13, hlp1)
          ENDDO
        ENDDO
        ! pool%reserve%carbon14
        IF (lctlib%ps_pathway == ic3phot) THEN
          hlp1 = discrimination_ps_c14(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ELSE
          hlp1 = discrimination_ps_c14(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, lctlib%ps_pathway)
        ENDIF
        DO ic = 1,nproma
          DO iblk = 1,nblks
            veg_pool_mt_domain(ix_reserve, ixC14, ic, iblk) = calc_mixing_ratio_C14C(def_co2_deltaC13, &
              &                                                           def_co2_deltaC14 - hlp1)     &
              &                                               * veg_pool_mt_domain(ix_reserve, ixC, ic, iblk)
            veg_pool_mt_domain(ix_reserve, ixN15, ic, iblk) = veg_pool_mt_domain(ix_reserve, ixN, ic, iblk) &
              &                                    / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ENDDO
        ENDDO
      CASE DEFAULT
        WRITE(message_text,'(2a)') 'not a valid veg_dynamics_scheme: ',TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme)
        CALL finish(TRIM(routine), message_text)
      ENDSELECT
    !>  1.2 Canopy
    !>
    CASE(QCANOPY)
      height(:,:) = dsl4jsb_Lctlib_param(vegetation_height)
    !>  1.3 Test_Canopy & Test_Radiation
    !>
    CASE(Q_TEST_CANOPY, Q_TEST_RADIATION)
      ref_lai   = 6.0_wp  !! the ref_lai of 6.0 works great for all sites that are not too dry; but use 0.5_wp for dry sites like ZA-Kru
      ref_nleaf = 2.0_wp
      lai(:,:)  = ref_lai
      !leaf_nitrogen
      IF (dsl4jsb_Config(VEG_)%pft_id == 0) THEN
        DO ic = 1,nproma
          DO iblk = 1,nblks
            veg_pool_mt_domain(ix_leaf, ixN, ic, iblk) = ref_nleaf / molar_mass_N * lai(ic, iblk)
          ENDDO
        ENDDO
        height                          = 10.0_wp
      ELSE
        DO ic = 1,nproma
          DO iblk = 1,nblks
            veg_pool_mt_domain(ix_leaf, ixN, ic, iblk) = lai(ic, iblk) / lctlib%sla / lctlib%cn_leaf
          ENDDO
        ENDDO
        veg_pool_mt_domain(ix_leaf, ixP, :, :) = veg_pool_mt_domain(ix_leaf, ixN, :, :) / lctlib%np_leaf
        IF (lctlib%growthform == ITREE) THEN
          height(:,:) = 10.0_wp
        ELSE
          height(:,:) = 1.0_wp
        ENDIF
      ENDIF
      leaf_nitrogen_cl(:,:,:) = 0.0_wp  ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)
      lai_cl(:,:,:)           = 0.0_wp  ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)
    !>  1.4 Default (error catching)
    !>
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid QUINCY model name. Use either land, plant, soil, canopy, test_canopy or test_radiation.')
    ENDSELECT

    !> 2.0 root fraction
    !>

    ! QSOIL, QCANOPY, Q_TEST_CANOPY and Q_TEST_RADIATION need fixed roots.
    ! (QSOIL cannot be properly initialised with dyn roots and
    ! the other modes would be subject to water limitations with dyn roots.)
    ! We thus switch off dynamic roots for these modes.
    SELECT CASE(model%config%qmodel_id)
    CASE(QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION)
      dsl4jsb_Config(VEG_)%flag_dynamic_roots = .FALSE.
      IF (iblk==1) CALL message(TRIM(routine), &
          & ' set "flag_dynamic_roots = .FALSE." (canopy, soil only, test_canopy or test_radiation model)')
    END SELECT

    ! init root fraction
    IF (dsl4jsb_Config(VEG_)%flag_dynamic_roots) THEN
      !>     - dynamic root fractions (see mo_q_veg_growth)
      !>
      ! all roots are in the first soil layer
      root_fraction_sl(:, 1, :)          = 1.0_wp
      root_fraction_sl(:, 2:nsoil_w, :)  = 0.0_wp
    ELSE
      !>     - fixed root fractions (no changes after init) depending on root_depth (SPQ variable)
      !>
      ! root fractions assuming an exponential distribution of roots \n
      !   Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth) \n
      ! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) \n
      !   thus for each layer
      !
      ! NOTE: upper bound of a layer is the larger value compared to lower bound !
      !
      DO iblk = 1,nblks
        DO ic = 1,nproma
          ! root depth larger than thickness of first soil layer
          IF (root_depth(ic, iblk) > soil_lay_depth_ubound_sl(ic, 1, iblk)) THEN
            root_fraction_sl(ic, 1, iblk) = (1._wp - EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(ic, 1, iblk))) &
              &                             / (1._wp - EXP(-lctlib%k_root_dist * root_depth(ic, iblk)))
            ! soil layers 2 to x
            !
            DO is = 2,INT(num_sl_above_bedrock(ic, iblk))
              IF (root_depth(ic, iblk) > soil_lay_depth_ubound_sl(ic, is-1, iblk)) THEN
                IF (root_depth(ic, iblk) > soil_lay_depth_ubound_sl(ic, is, iblk)) THEN
                  root_fraction_sl(ic, is, iblk) = (EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(ic, is-1, iblk)) &
                    &                                - EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(ic, is, iblk))) &
                    &                                / (1._wp - EXP(-lctlib%k_root_dist * root_depth(ic, iblk)))
                ELSE
                  root_fraction_sl(ic, is, iblk) = (EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(ic, is-1, iblk)) &
                    &                                - EXP(-lctlib%k_root_dist * root_depth(ic, iblk))) &
                    &                                / (1._wp - EXP(-lctlib%k_root_dist * root_depth(ic, iblk)))
                END IF
              ELSE
                root_fraction_sl(ic, is, iblk) = 0.0_wp
              END IF
            END DO
          ! root depth lower than / equal thickness of first soil layer
          ! all roots are in the first soil layer
          ELSE
            root_fraction_sl(ic, :, iblk) = 0.0_wp
            root_fraction_sl(ic, 1, iblk) = 1.0_wp
          END IF
        END DO
      END DO
    END IF  ! flag_dynamic_roots

    !>   2.1 set root fraction to zero for tiles such as bare soil and urban area
    !>
    IF (lctlib%BareSoilFlag) THEN
      root_fraction_sl(:, :, :) = 0.0_wp
    END IF

    !> 3.0 set veg_pool bgcm to zero for tiles such as bare soil and urban area
    !>
    IF (lctlib%BareSoilFlag) THEN
      veg_pool_mt_domain(:,:,:,:) = 0.0_wp
    END IF

    !> 4.0 vegetation init for all QUINCY models
    !>
    ! plant variables
    dphi(:,:)                             = -1.0_wp * lctlib%phi_leaf_min
    beta_sinklim_ps(:,:)                  = 1.0_wp
    leaf_cn_direction(:,:)                = 0.0_wp
    fn_oth_cl(:,:,:)                      = 1.0_wp  ! it is necessary to start with fn_oth = 1.0 leaving the rest of the fraction 0.0
                                                    ! fn_chl_cl, fn_rub_cl, fn_et_cl, fn_pepc_cl are set to 0.0_wp by default with the add_var()
    ! initialising moving averages
    ! allocation
    npp_talloc_mavg(:,:)                  = 0.0_wp
    n_fixation_talloc_mavg(:,:)           = 0.0_wp
    transpiration_talloc_mavg(:,:)        = eps8
    unit_transpiration_talloc_mavg(:,:)   = 1.e-6_wp
    unit_uptake_n_talloc_mavg(:,:)        = eps8
    unit_uptake_p_talloc_mavg(:,:)        = 4._wp * 1.e-10_wp
    dphi_talloc_mavg(:,:)                 = -1.0_wp * lctlib%phi_leaf_min
    growth_cn_talloc_mavg(:,:)            = 60.0_wp
    unit_npp_talloc_mavg(:,:)             = unit_uptake_n_talloc_mavg(:,:) * growth_cn_talloc_mavg(:,:)
    growth_cp_talloc_mavg(:,:)            = 2000.0_wp
    growth_np_talloc_mavg(:,:)            = growth_cp_talloc_mavg(:,:) / growth_cn_talloc_mavg(:,:)
    labile_carbon_talloc_mavg(:,:)        = 0.0_wp
    labile_nitrogen_talloc_mavg(:,:)      = 0.0_wp
    labile_phosphorus_talloc_mavg(:,:)    = 0.0_wp
    growth_p_limit_based_on_n_mavg(:,:)   = 1.0_wp
    reserve_carbon_talloc_mavg(:,:)       = 0.0_wp
    reserve_nitrogen_talloc_mavg(:,:)     = 0.0_wp
    reserve_phosphorus_talloc_mavg(:,:)   = 0.0_wp
    wtr_rootzone_rel_talloc_mavg(:,:)     = k1_root_alloc
    exudation_c_tmyc_mavg(:,:)            = 0.0_wp
    ! N fraction optimisation
    fleaf_sunlit_tfrac_mavg_cl(:,:,:)     = 0.5_wp
    t_air_tfrac_mavg(:,:)                 = 283.15_wp
    press_srf_tfrac_mavg(:,:)             = standard_press_srf
    co2_mixing_ratio_tfrac_mavg(:,:)      = def_co2_mixing_ratio
    ga_tfrac_mavg(:,:)                    = 0.2_wp
    beta_sinklim_ps_tfrac_mavg(:,:)       = 1.0_wp
    ! leaf N optimisation
    fleaf_sunlit_tcnl_mavg_cl(:,:,:)      = 0.5_wp
    t_air_tcnl_mavg(:,:)                  = 283.15_wp
    press_srf_tcnl_mavg(:,:)              = standard_press_srf
    co2_mixing_ratio_tcnl_mavg(:,:)       = def_co2_mixing_ratio
    growth_cn_tcnl_mavg(:,:)              = 60.0_wp
    growth_np_tcnl_mavg(:,:)              = 30.0_wp
    uptake_n_tcnl_mavg(:,:)               = eps8
    npp_tcnl_mavg(:,:)                    = uptake_n_tcnl_mavg(:,:) * growth_cn_tcnl_mavg(:,:)
    ! daytime averages plus accumulation var
    t_air_daytime(:,:)                    = 283.15_wp
    press_srf_daytime(:,:)                = standard_press_srf
    co2_mixing_ratio_daytime(:,:)         = def_co2_mixing_ratio
    ga_daytime(:,:)                       = 0.2_wp
    beta_sinklim_ps_daytime(:,:)          = 1.0_wp
    t_jmax_opt_daytime(:,:)               = 17.0_wp
    t_jmax_opt_daytime_dacc(:,:)          = 17.0_wp
    t_air_daytime_dacc(:,:)               = 283.15_wp
    press_srf_daytime_dacc(:,:)           = standard_press_srf
    co2_mixing_ratio_daytime_dacc(:,:)    = def_co2_mixing_ratio
    ga_daytime_dacc(:,:)                  = 0.2_wp
    beta_sinklim_ps_daytime_dacc(:,:)     = 1.0_wp
    ! moving averages
    growth_req_n_tlabile_mavg(:,:)        = 1._wp / 25._wp
    growth_req_p_tlabile_mavg(:,:)        = 1._wp / 14._wp
    t_air_tphen_mavg(:,:)                 = 283.15_wp
    t_soil_srf_tphen_mavg(:,:)            = 283.15_wp
    npp_tuptake_mavg(:,:)                 = 1.0_wp
    demand_uptake_n_tuptake_mavg(:,:)     = 0.1_wp
    demand_uptake_p_tuptake_mavg(:,:)     = 0.1_wp
    growth_req_n_tuptake_mavg(:,:)        = 1._wp / 25._wp
    growth_req_p_tuptake_mavg(:,:)        = 1._wp / 14._wp
    ga_tcnl_mavg(:,:)                     = 0.2_wp
    leaf2root_troot_mavg(:,:)             = 1.0_wp
    an_boc_tvegdyn_mavg(:,:)              = 5.0_wp
    net_growth_tvegdyn_mavg(:,:)          = 5.0_wp
    lai_tvegdyn_mavg(:,:)                 = 0.1_wp
    beta_sinklim_ps_tacclim_mavg(:,:)     = 1.0_wp
    t_air_tacclim_mavg(:,:)               = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_q_veg_respiration:temperature_response_respiration
    t_soil_root_tacclim_mavg(:,:)         = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_q_veg_respiration:temperature_response_respiration
    t_air_week_mavg(:,:)                  = 283.15_wp
    t_air_month_mavg(:,:)                 = 283.15_wp
    t_jmax_opt_mavg(:,:)                  = 17.0_wp
    ! bottom of the canopy radiation
    rfr_ratio_boc(:,:)                    = rfr_ratio_toc
    rfr_ratio_boc_tvegdyn_mavg(:,:)       = rfr_ratio_boc(:,:)

    !$ACC UPDATE DEVICE(root_fraction_sl(:,:,:)) &
    !$ACC   DEVICE(height(:,:), leaf_nitrogen_cl(:,:,:), lai(:,:), lai_cl(:,:,:)) &
    !$ACC   DEVICE(dphi(:,:), beta_sinklim_ps(:,:), leaf_cn_direction(:,:), fn_oth_cl(:,:,:), npp_talloc_mavg(:,:)) &
    !$ACC   DEVICE(n_fixation_talloc_mavg(:,:), transpiration_talloc_mavg(:,:), unit_transpiration_talloc_mavg(:,:)) &
    !$ACC   DEVICE(unit_uptake_n_talloc_mavg(:,:), unit_uptake_p_talloc_mavg(:,:), dphi_talloc_mavg(:,:)) &
    !$ACC   DEVICE(growth_cn_talloc_mavg(:,:), unit_npp_talloc_mavg(:,:), growth_cp_talloc_mavg(:,:)) &
    !$ACC   DEVICE(growth_np_talloc_mavg(:,:), labile_carbon_talloc_mavg(:,:), labile_nitrogen_talloc_mavg(:,:)) &
    !$ACC   DEVICE(labile_phosphorus_talloc_mavg(:,:), growth_p_limit_based_on_n_mavg(:,:), reserve_carbon_talloc_mavg(:,:)) &
    !$ACC   DEVICE(reserve_nitrogen_talloc_mavg(:,:), reserve_phosphorus_talloc_mavg(:,:), wtr_rootzone_rel_talloc_mavg(:,:)) &
    !$ACC   DEVICE(exudation_c_tmyc_mavg(:,:), fleaf_sunlit_tfrac_mavg_cl(:,:,:), t_air_tfrac_mavg(:,:)) &
    !$ACC   DEVICE(press_srf_tfrac_mavg(:,:), co2_mixing_ratio_tfrac_mavg(:,:), ga_tfrac_mavg(:,:)) &
    !$ACC   DEVICE(beta_sinklim_ps_tfrac_mavg(:,:), fleaf_sunlit_tcnl_mavg_cl(:,:,:), t_air_tcnl_mavg(:,:)) &
    !$ACC   DEVICE(press_srf_tcnl_mavg(:,:), co2_mixing_ratio_tcnl_mavg(:,:), growth_cn_tcnl_mavg(:,:)) &
    !$ACC   DEVICE(growth_np_tcnl_mavg(:,:), uptake_n_tcnl_mavg(:,:), npp_tcnl_mavg(:,:), t_air_daytime(:,:)) &
    !$ACC   DEVICE(press_srf_daytime(:,:), co2_mixing_ratio_daytime(:,:), ga_daytime(:,:), beta_sinklim_ps_daytime(:,:)) &
    !$ACC   DEVICE(t_jmax_opt_daytime(:,:), t_jmax_opt_daytime_dacc(:,:), t_air_daytime_dacc(:,:)) &
    !$ACC   DEVICE(press_srf_daytime_dacc(:,:), co2_mixing_ratio_daytime_dacc(:,:), ga_daytime_dacc(:,:)) &
    !$ACC   DEVICE(beta_sinklim_ps_daytime_dacc(:,:), growth_req_n_tlabile_mavg(:,:), growth_req_p_tlabile_mavg(:,:)) &
    !$ACC   DEVICE(t_air_tphen_mavg(:,:), t_soil_srf_tphen_mavg(:,:), npp_tuptake_mavg(:,:)) &
    !$ACC   DEVICE(demand_uptake_n_tuptake_mavg(:,:), demand_uptake_p_tuptake_mavg(:,:), growth_req_n_tuptake_mavg(:,:)) &
    !$ACC   DEVICE(growth_req_p_tuptake_mavg(:,:), ga_tcnl_mavg(:,:), leaf2root_troot_mavg(:,:), an_boc_tvegdyn_mavg(:,:)) &
    !$ACC   DEVICE(net_growth_tvegdyn_mavg(:,:), lai_tvegdyn_mavg(:,:), beta_sinklim_ps_tacclim_mavg(:,:)) &
    !$ACC   DEVICE(t_air_tacclim_mavg(:,:), t_soil_root_tacclim_mavg(:,:), t_air_week_mavg(:,:), t_air_month_mavg(:,:)) &
    !$ACC   DEVICE(t_jmax_opt_mavg(:,:), veg_pool_mt_domain(:,:,:,:)) &
    !$ACC   DEVICE(rfr_ratio_boc(:,:), rfr_ratio_boc_tvegdyn_mavg(:,:)) ASYNC(1)
    !$ACC WAIT(1)

  END SUBROUTINE veg_init

#ifndef __QUINCY_STANDALONE__
  ! ====================================================================================================== !
  !
  !> Read vegetation state from file (after init/restart)
  !> (variables should be in sync with those listed in ./scripts/quincy_sb_and_veg_state_from_restart_file.sh)
  !
  ! ======================================================================================================= !
  SUBROUTINE veg_read_states(tile)
    USE mo_jsb_io_netcdf,           ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_util,                    ONLY: read_2D_var, read_3D_var
    USE mo_lnd_bgcm_store,          ONLY: read_bgcm_from_file
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_PP_CROP_ID, VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, &
      & VEG_BGCM_PP_FIBERBOARD_ID, VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID, VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_SEED_BED_POOL_ID
    dsl4jsb_Use_config(Q_AGR_)
    dsl4jsb_Use_config(Q_SYL_)
    dsl4jsb_Use_config(ALCC_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile             !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model            !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store       !< the bgcm store of this tile
    TYPE(t_input_file)              :: input_file       !< file with to be read veg state variables
    TYPE(t_jsb_vgrid),      POINTER :: vgrid_soil_w     !< vertical grid for sl variables
    INTEGER                         :: isoil            !< soillayer index
    INTEGER                         :: nsoil_w          !< number of soil layers (water)
    REAL(wp),               POINTER :: ptr3D(:,:,:)     !< temporary pointers for 3D variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_read_states'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_SYL_)
    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_config(ALCC_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Real2D_onDomain :: dens_ind
    dsl4jsb_Real2D_onDomain :: mean_leaf_age
    dsl4jsb_Real2D_onDomain :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onDomain :: unit_npp_troot_mavg
    dsl4jsb_Real2D_onDomain :: unit_uptake_n_pot_troot_mavg
    dsl4jsb_Real2D_onDomain :: unit_uptake_p_pot_troot_mavg
    dsl4jsb_Real2D_onDomain :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain :: lai_tvegdyn_mavg
    dsl4jsb_Real3D_onDomain :: root_fraction_sl
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    vgrid_soil_w => Get_vgrid('soil_depth_water')
    nsoil_w = vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_SYL_)
    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_var2D_onDomain(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onDomain(VEG_, mean_leaf_age)
    dsl4jsb_Get_var2D_onDomain(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_npp_troot_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_uptake_n_pot_troot_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_uptake_p_pot_troot_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, net_growth_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, lai_tvegdyn_mavg)
    dsl4jsb_Get_var3D_onDomain(VEG_, root_fraction_sl)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), ' for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(VEG_)%ic_filename), model%grid_id)

    IF (ASSOCIATED(tile%bgcm_store)) THEN
      bgcm_store => tile%bgcm_store
      CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_POOL_ID, "veg", tile%name, routine, input_file)
      CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_SEED_BED_POOL_ID, "veg", tile%name, routine, input_file)

      IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        ! NOTE: here it is assumed that the state file only contains crop products on crop tiles and
        !       wood products on forst tiles. This would not be correct if the states stem from a simulation
        !       with alcc! Because in a simulation with alcc the product pools are passively relocated and
        !       it is thus possible that crop products end up in grass or forest tiles and wood products can
        !       in turn end up on grass or crop tiles, respectively.

        IF (dsl4jsb_Config(Q_AGR_)%active .OR. dsl4jsb_Config(ALCC_)%active) THEN
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_CROP_ID, "veg", tile%name, routine, input_file)
        END IF

        IF (dsl4jsb_Config(Q_SYL_)%active .OR. dsl4jsb_Config(ALCC_)%active) THEN
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_FUEL_ID, "veg", tile%name, routine, input_file)
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_PAPER_ID, "veg", tile%name, routine, input_file)
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_FIBERBOARD_ID, "veg", tile%name, routine, input_file)
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_OIRW_ID, "veg", tile%name, routine, input_file)
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_PV_ID, "veg", tile%name, routine, input_file)
          CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_PP_SAWNWOOD_ID, "veg", tile%name, routine, input_file)
        END IF
      END IF
    END IF

    CALL read_2D_var(input_file, dens_ind, "veg_dens_ind_"//tile%name)
    CALL read_2D_var(input_file, mean_leaf_age, "veg_mean_leaf_age_"//tile%name)
    CALL read_2D_var(input_file, leaf2root_troot_mavg, "veg_leaf2root_troot_mavg_"//tile%name)
    CALL read_2D_var(input_file, unit_npp_troot_mavg, "veg_unit_npp_troot_mavg_"//tile%name)
    CALL read_2D_var(input_file, unit_uptake_n_pot_troot_mavg, "veg_unit_uptake_n_pot_troot_mavg_"//tile%name)
    CALL read_2D_var(input_file, unit_uptake_p_pot_troot_mavg, "veg_unit_uptake_p_pot_troot_mavg_"//tile%name)
    CALL read_2D_var(input_file, net_growth_tvegdyn_mavg, "veg_net_growth_tvegdyn_mavg_"//tile%name)
    CALL read_2D_var(input_file, lai_tvegdyn_mavg, "veg_lai_tvegdyn_mavg_"//tile%name)

    CALL read_3D_var(input_file, root_fraction_sl, 'veg_root_fraction_sl_'//tile%name, 'soillev', nsoil_w)

    CALL input_file%Close()

  END SUBROUTINE veg_read_states

#endif

#endif
END MODULE mo_veg_init
