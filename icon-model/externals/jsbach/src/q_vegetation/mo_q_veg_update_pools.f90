!> QUINCY update vegetation pools
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
!>#### routines for calculating the vegetation pools using the fluxes that have been calculated
!>
MODULE mo_q_veg_update_pools
#ifndef __NO_QUINCY__

  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, finish, message_text

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_SEED_BED_POOL_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_SEED_BED_GROWTH_ID, VEG_BGCM_LITTERFALL_ID, &
    &                                   VEG_BGCM_SEED_BED_LITTERFALL_ID, VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, &
    &                                   VEG_BGCM_RESERVE_USE_ID, VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, &
    &                                   VEG_BGCM_PP_FIBERBOARD_ID, VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID,  &
    &                                   VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_PROD_DECAY_ID, VEG_BGCM_PP_CROP_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_MYCO_EXPORT_ID

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_update_pools'

CONTAINS

  ! ======================================================================================================= !
  !>update vegetation pools
  !>
  SUBROUTINE update_veg_pools(tile, options)
    USE mo_kind,                            ONLY: wp
    USE mo_jsb_impl_constants,              ONLY: test_false_true
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_quincy_model_config,             ONLY: QLAND, QPLANT, QCANOPY
    USE mo_jsb_process_class,               ONLY: HYDRO_, VEG_, Q_ASSIMI_, Q_PHENO_, Q_AGR_, Q_SYL_, ALCC_
    USE mo_jsb_grid_class,                  ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_grid,                        ONLY: Get_vgrid, Get_grid
    USE mo_jsb_physical_constants,          ONLY: lambda_C14
    USE mo_jsb_math_constants,              ONLY: one_day, one_year, eps4, eps8, eps12
    USE mo_isotope_util,                    ONLY: calc_mixing_ratio_N15N14
    USE mo_veg_constants,                   ONLY: eta_nfixation, max_leaf_shedding_rate
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(ALCC_)
    dsl4jsb_Use_config(Q_SYL_)
    dsl4jsb_Use_config(Q_AGR_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_jsb_grid),         POINTER :: grid                   !< horizontal grid
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_soil_w           !< Vertical grid
    INTEGER                           :: nsoil_w                !< number of soil layers (water)
    REAL(wp)                          :: dtime                  !< timestep length
    INTEGER                           :: iblk, ics, ice, nc     !< grid dimensions
    INTEGER                           :: ic, is, icanopy        !< looping indices
    INTEGER                           :: id_elem, ix_elem       !< id and index of bgcm elements
    INTEGER                           :: ix_comp                !< currently id and index of bgcm compartments
    INTEGER                           :: nr_of_veg_bgcm_comp    !< number of compartments in veg bgcms
    REAL(wp) :: lctlib_sla, lctlib_tau_leaf, lctlib_cn_leaf, lctlib_np_leaf
    REAL(wp) :: hlp
    INTEGER  :: lctlib_phenology_type
    LOGICAL  :: config_flag_log_negative_vegpool                !< debug option of veg config
    LOGICAL  :: config_l_use_product_pools, is_active_SYL_or_ALCC, is_active_AGR_or_ALCC
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_pool_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_growth_mt
    dsl4jsb_Def_mt1L2D :: veg_exudation_mt
    dsl4jsb_Def_mt1L2D :: veg_establishment_mt
    dsl4jsb_Def_mt1L2D :: veg_reserve_use_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt

    dsl4jsb_Def_mt1L2D :: veg_pp_crop_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_fuel_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_paper_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_fiberboard_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_oirw_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_pv_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_sawnwood_mt
    dsl4jsb_Def_mt1L2D :: prod_decay_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(ALCC_)
    dsl4jsb_Def_config(Q_SYL_)
    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: lai_max
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: net_growth
    dsl4jsb_Real2D_onChunk      :: npp_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: npp_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: n_fixation_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp_c13
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp_c14
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp_c13
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp_c14
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: net_biosphere_production
    dsl4jsb_Real2D_onChunk      :: net_biosphere_production_c13
    dsl4jsb_Real2D_onChunk      :: net_biosphere_production_c14
    dsl4jsb_Real2D_onChunk      :: biological_n_fixation
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_p
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_p
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_p
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_p
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_n
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_p
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_n
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_p
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_products_total_c
    dsl4jsb_Real2D_onChunk      :: veg_products_total_n
    dsl4jsb_Real2D_onChunk      :: veg_products_total_p
    dsl4jsb_Real2D_onChunk      :: veg_products_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_products_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_products_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_products_decay_c
    dsl4jsb_Real2D_onChunk      :: veg_products_decay_n
    dsl4jsb_Real2D_onChunk      :: veg_products_decay_p
    dsl4jsb_Real2D_onChunk      :: veg_products_decay_c13
    dsl4jsb_Real2D_onChunk      :: veg_products_decay_c14
    dsl4jsb_Real2D_onChunk      :: veg_products_decay_n15
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    dtime     = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    grid                  => Get_grid(model%grid_id)
    vgrid_soil_w          => Get_vgrid('soil_depth_water')
    nsoil_w               =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_config(Q_SYL_)
    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    config_flag_log_negative_vegpool = dsl4jsb_Config(VEG_)%flag_log_negative_vegpool
    config_l_use_product_pools = dsl4jsb_Config(VEG_)%l_use_product_pools

    is_active_AGR_or_ALCC = .FALSE.
    is_active_SYL_or_ALCC = .FALSE.
    IF (model%config%qmodel_id == QLAND .OR. model%config%qmodel_id == QPLANT) THEN
      IF (dsl4jsb_Config(Q_AGR_)%active .OR. dsl4jsb_Config(ALCC_)%active) THEN
        is_active_AGR_or_ALCC = .TRUE.
      END IF
      IF (dsl4jsb_Config(Q_SYL_)%active .OR. dsl4jsb_Config(ALCC_)%active) THEN
        is_active_SYL_or_ALCC = .TRUE.
      END IF
    END IF
    ! ----------------------------------------------------------------------------------------------------- !

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_POOL_ID, seed_bed_pool_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_GROWTH_ID, seed_bed_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_EXUDATION_ID, veg_exudation_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_ESTABLISHMENT_ID, veg_establishment_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_RESERVE_USE_ID, veg_reserve_use_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)

    IF (config_l_use_product_pools) THEN
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PROD_DECAY_ID, prod_decay_mt)

      IF (is_active_AGR_or_ALCC) THEN
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_CROP_ID, veg_pp_crop_mt)
      END IF

      IF (is_active_SYL_or_ALCC) THEN
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt)
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PAPER_ID, veg_pp_paper_mt)
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FIBERBOARD_ID, veg_pp_fiberboard_mt)
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_OIRW_ID, veg_pp_oirw_mt)
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PV_ID, veg_pp_pv_mt)
        dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_SAWNWOOD_ID, veg_pp_sawnwood_mt)
      END IF
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation)        ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C13)    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C14)    ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)             ! in
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, lai_max)                    ! in
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)                ! in
    ! VEG 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)                                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)                           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c13)                              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c14)                              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp_c13)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp_c14)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp_c13)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp_c14)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_c13)         ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_c14)         ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, biological_n_fixation)                ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_p)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c13)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c14)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n15)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_c)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_n)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_p)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_c)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_n)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_p)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_c)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_n)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_p)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_n)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_p)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c13)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c14)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_n15)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_n)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_p)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c13)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c14)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_n15)             ! out
    IF (config_l_use_product_pools) THEN
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_c)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_n)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_p)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_c13)             ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_c14)             ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_n15)             ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_decay_c)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_decay_n)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_decay_p)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_decay_c13)             ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_decay_c14)             ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_decay_n15)             ! out
    END IF
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.8 get number of bgcm compartments
    !>
    nr_of_veg_bgcm_comp = SIZE(veg_pool_mt,1)

    lctlib_sla = dsl4jsb_Lctlib_param(sla)
    lctlib_phenology_type = dsl4jsb_Lctlib_param(phenology_type)
    lctlib_tau_leaf = dsl4jsb_Lctlib_param(tau_leaf)
    lctlib_cn_leaf = dsl4jsb_Lctlib_param(cn_leaf)
    lctlib_np_leaf = dsl4jsb_Lctlib_param(np_leaf)

    !>1.0 for plant and land model, update vegetation pools with fluxes
    !>
    SELECT CASE(model%config%qmodel_id)
    CASE(QPLANT, QLAND)
#ifdef _OPENACC
      CALL finish(routine, 'Code blocks for QPLANT & QLAND modes not ported to GPU, yet. Stop.')
#endif

      DO ic = 1, nc

       !>  1.1 net balance of the labile carbon pool due to photosynthesis and respiration
       !>
       npp(ic)        = gross_assimilation(ic)     - &
                        growth_respiration(ic)     - maint_respiration(ic)     - n_processing_respiration(ic)
       npp_c13(ic)    = gross_assimilation_C13(ic) - &
                        growth_respiration_c13(ic) - maint_respiration_c13(ic) - n_processing_respiration_c13(ic)
       npp_c14(ic)    = gross_assimilation_C14(ic) - &
                        growth_respiration_c14(ic) - maint_respiration_c14(ic) - n_processing_respiration_c14(ic)
       net_growth(ic) = net_growth(ic) + npp(ic)

       veg_growth_mt(ix_labile, ixC, ic)   = veg_growth_mt(ix_labile, ixC, ic)   + npp(ic)     * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixC13, ic) = veg_growth_mt(ix_labile, ixC13, ic) + npp_c13(ic) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixC14, ic) = veg_growth_mt(ix_labile, ixC14, ic) + npp_c14(ic) * dtime / 1000000.0_wp

       !>  1.2 Add nutrients from uptake calculations
       !>
       n_fixation_n15(ic)                  = n_fixation(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))
       veg_growth_mt(ix_labile, ixN, ic)   = veg_growth_mt(ix_labile, ixN, ic) + &
                                             ( uptake_nh4(ic) + uptake_no3(ic) + n_fixation(ic)) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixN15, ic) = veg_growth_mt(ix_labile, ixN15, ic) + &
                                             ( uptake_nh4_n15(ic) + uptake_no3_n15(ic) + n_fixation_n15(ic)) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixP, ic)   = veg_growth_mt(ix_labile, ixP, ic) + &
                                             uptake_po4(ic) * dtime / 1000000.0_wp

       !>  1.3 Add mycorrhiza export
       !>
       DO is = 1, nsoil_w
        veg_growth_mt(ix_labile, ixC, ic)   = veg_growth_mt(ix_labile, ixC, ic)   + &
          &                                   sb_mycorrhiza_export_mt(ixC, ic, is)    * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixC13, ic) = veg_growth_mt(ix_labile, ixC13, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixC13, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixC14, ic) = veg_growth_mt(ix_labile, ixC14, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixC14, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixN, ic)   = veg_growth_mt(ix_labile, ixN, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixN, ic, is)    * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixN15, ic) = veg_growth_mt(ix_labile, ixN15, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixN15, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixP, ic)   = veg_growth_mt(ix_labile, ixP, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixP, ic, is)    * soil_depth_sl(ic,is)
       END DO

      END DO

      !>  1.4 update labile pool and reserve
      !>
      veg_pool_mt(ix_labile, :, :)  = veg_pool_mt(ix_labile, :, :)  + veg_reserve_use_mt(:,:)
      veg_pool_mt(ix_reserve, :, :) = veg_pool_mt(ix_reserve, :, :) - veg_reserve_use_mt(:,:)

      DO ic = 1, nc
        !>  1.5 remove respiration flux associated with herbivory from affected pools
        !>      the remaining matter is already included in litter flux
        !>
        veg_pool_mt(ix_leaf, ixC, ic)    = veg_pool_mt(ix_leaf, ixC, ic)    &
          &                                - herbivory_leaf_resp(ic)      * dtime / 1.e6_wp
        veg_pool_mt(ix_leaf, ixC13, ic)  = veg_pool_mt(ix_leaf, ixC13, ic)  &
          &                                - herbivory_leaf_resp_c13(ic)  * dtime / 1.e6_wp
        veg_pool_mt(ix_leaf, ixC14, ic)  = veg_pool_mt(ix_leaf, ixC14, ic)  &
          &                                - herbivory_leaf_resp_c14(ic)  * dtime / 1.e6_wp
        veg_pool_mt(ix_fruit, ixC, ic)   = veg_pool_mt(ix_fruit, ixC, ic)   &
          &                                - herbivory_fruit_resp(ic)     * dtime / 1.e6_wp
        veg_pool_mt(ix_fruit, ixC13, ic) = veg_pool_mt(ix_fruit, ixC13, ic) &
          &                                - herbivory_fruit_resp_c13(ic) * dtime / 1.e6_wp
        veg_pool_mt(ix_fruit, ixC14, ic) = veg_pool_mt(ix_fruit, ixC14, ic) &
          &                                - herbivory_fruit_resp_c14(ic) * dtime / 1.e6_wp

        !>  1.6 update fluxes from continuous nutrient turnover and resorption
        !>
        ! N
        veg_pool_mt(ix_labile, ixN, ic)     = veg_pool_mt(ix_labile, ixN, ic) &
          &                                   + recycling_leaf_n(ic) + recycling_fine_root_n(ic) + recycling_heart_wood_n(ic)
        veg_pool_mt(ix_leaf, ixN, ic)       = veg_pool_mt(ix_leaf, ixN, ic)       - recycling_leaf_n(ic)
        veg_pool_mt(ix_fine_root, ixN, ic)  = veg_pool_mt(ix_fine_root, ixN, ic)  - recycling_fine_root_n(ic)
        veg_pool_mt(ix_heart_wood, ixN, ic) = veg_pool_mt(ix_heart_wood, ixN, ic) - recycling_heart_wood_n(ic)
        ! P
        veg_pool_mt(ix_labile, ixP, ic)     = veg_pool_mt(ix_labile, ixP, ic) &
          &                                   + recycling_leaf_p(ic) + recycling_fine_root_p(ic) + recycling_heart_wood_p(ic)
        veg_pool_mt(ix_leaf, ixP, ic)       = veg_pool_mt(ix_leaf, ixP, ic)       - recycling_leaf_p(ic)
        veg_pool_mt(ix_fine_root, ixP, ic)  = veg_pool_mt(ix_fine_root, ixP, ic)  - recycling_fine_root_p(ic)
        veg_pool_mt(ix_heart_wood, ixP, ic) = veg_pool_mt(ix_heart_wood, ixP, ic) - recycling_heart_wood_p(ic)
        ! N15
        veg_pool_mt(ix_labile, ixN15, ic)     = veg_pool_mt(ix_labile, ixN15, ic) + recycling_leaf_n15(ic) &
          &                                     + recycling_fine_root_n15(ic) + recycling_heart_wood_n15(ic)
        veg_pool_mt(ix_leaf, ixN15, ic)       = veg_pool_mt(ix_leaf, ixN15, ic)       - recycling_leaf_n15(ic)
        veg_pool_mt(ix_fine_root, ixN15, ic)  = veg_pool_mt(ix_fine_root, ixN15, ic)  - recycling_fine_root_n15(ic)
        veg_pool_mt(ix_heart_wood, ixN15, ic) = veg_pool_mt(ix_heart_wood, ixN15, ic) - recycling_heart_wood_n15(ic)

        !>  1.7 update labile pool with fluxes to mycorrhiza
        !>
        veg_pool_mt(ix_labile, ixC, ic)   = veg_pool_mt(ix_labile, ixC, ic)   - veg_exudation_mt(ixC, ic)
        veg_pool_mt(ix_labile, ixN, ic)   = veg_pool_mt(ix_labile, ixN, ic)   - veg_exudation_mt(ixN, ic)
        veg_pool_mt(ix_labile, ixP, ic)   = veg_pool_mt(ix_labile, ixP, ic)   - veg_exudation_mt(ixP, ic)
        veg_pool_mt(ix_labile, ixC13, ic) = veg_pool_mt(ix_labile, ixC13, ic) - veg_exudation_mt(ixC13, ic)
        veg_pool_mt(ix_labile, ixC14, ic) = veg_pool_mt(ix_labile, ixC14, ic) - veg_exudation_mt(ixC14, ic)
        veg_pool_mt(ix_labile, ixN15, ic) = veg_pool_mt(ix_labile, ixN15, ic) - veg_exudation_mt(ixN15, ic)
      END DO

      !>  1.8 Update plant biogeochemical pools with growth, litterfall due to natural processes and establishment
      !>
      DO ic = 1, nc
        DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (model%config%is_element_used(id_elem)) THEN
            ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
            DO ix_comp = 1, nr_of_veg_bgcm_comp
              veg_pool_mt(ix_comp, ix_elem, ic) = veg_pool_mt(ix_comp, ix_elem, ic) - veg_litterfall_mt(ix_comp, ix_elem, ic)
              veg_pool_mt(ix_comp, ix_elem, ic) = veg_pool_mt(ix_comp, ix_elem, ic) + veg_growth_mt(ix_comp, ix_elem, ic)
            END DO
          END IF
        seed_bed_pool_mt(ix_elem, ic) = seed_bed_pool_mt(ix_elem, ic) - seed_bed_litterfall_mt(ix_elem, ic)
        seed_bed_pool_mt(ix_elem, ic) = seed_bed_pool_mt(ix_elem, ic) + seed_bed_growth_mt(ix_elem, ic)
        END DO
      END DO

      DO ic = 1, nc
        DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (model%config%is_element_used(id_elem)) THEN
            ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
            veg_pool_mt(ix_labile,   ix_elem, ic) = veg_pool_mt(ix_labile,   ix_elem, ic) + veg_establishment_mt(ix_elem, ic)
            seed_bed_pool_mt(ix_elem, ic) = seed_bed_pool_mt(ix_elem, ic) - veg_establishment_mt(ix_elem, ic)
          END IF
        END DO
      END DO

      !>    1.8.1 ensure veg_pool does not become negative - quick-fix and debug
      !>
      DO ic = 1, nc
        DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (model%config%is_element_used(id_elem)) THEN
            ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
            DO ix_comp = 1, nr_of_veg_bgcm_comp
              ! set a negative pool to zero (breaking mass balance !)
              ! and write message to LOG (if enabled)
              ! SZ: temporarily switching off test for negative 14C pools, which occur in rare occasions.
              !     This change does not affect simulation outcomes other than for 14C runs, and for these only few grid cells
              IF (veg_pool_mt(ix_comp, ix_elem, ic) < -eps12 .AND. ix_elem /= ixC14) THEN
                IF (config_flag_log_negative_vegpool) THEN
                  WRITE (message_text,*) 'Negative value of vegetation pool at tile ', tile%name, &
                    &                    ' at (lon:lat) ', grid%lon(ic,iblk), ':', grid%lat(ic,iblk), &
                    &                    ', compartment: ', ix_comp, ' , element: ', ix_elem, &
                    &                    ', value ', veg_pool_mt(ix_comp, ix_elem, ic), ' is set to zero here'
                  CALL message(TRIM(routine), message_text, all_print=.TRUE.)
                END IF
              END IF
              IF (veg_pool_mt(ix_comp, ix_elem, ic) < 0.0_wp) THEN
                veg_pool_mt(ix_comp, ix_elem, ic) = 0.0_wp
              END IF
            END DO
            ! set a negative pool to zero (breaking mass balance !)
            ! and write message to LOG (if enabled)
            ! SZ: temporarily switching off test for negative 14C pools, which occur in rare occasions.
            !     This change does not affect simulation outcomes other than for 14C runs, and for these only few grid cells
            IF (seed_bed_pool_mt(ix_elem, ic) < -eps12 .AND. ix_elem /= ixC14) THEN
              IF (config_flag_log_negative_vegpool) THEN
                WRITE (message_text,*) 'Negative value of seed bed pool at tile ', tile%name, &
                  &                    ' at (lon:lat) ', grid%lon(ic,iblk), ':', grid%lat(ic,iblk), &
                  &                    ', element: ', ix_elem, &
                  &                    ', value ', seed_bed_pool_mt(ix_elem, ic), ' is set to zero here'
                CALL message(TRIM(routine), message_text, all_print=.TRUE.)
              END IF
            END IF
            IF (seed_bed_pool_mt(ix_elem, ic) < 0.0_wp) THEN
              seed_bed_pool_mt(ix_elem, ic) = 0.0_wp
            END IF
          END IF
        END DO
      END DO

      !>    1.8.2 Radioactive decay of C14
      !>
      DO ic = 1, nc
        ! in living plant tissue
        DO ix_comp = 1, nr_of_veg_bgcm_comp
          veg_pool_mt(ix_comp, ixC14, ic) = veg_pool_mt(ix_comp, ixC14, ic) * (1._wp - lambda_C14 * dtime)
        END DO
        seed_bed_pool_mt(ixC14, ic) = seed_bed_pool_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)

        ! and if simulating with product pools also within these
        IF (config_l_use_product_pools) THEN

          IF (is_active_AGR_or_ALCC) THEN
            veg_pp_crop_mt(ixC14, ic) = veg_pp_crop_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          END IF

          IF (is_active_SYL_or_ALCC) THEN
            veg_pp_fuel_mt(ixC14, ic) = veg_pp_fuel_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
            veg_pp_paper_mt(ixC14, ic) = veg_pp_paper_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
            veg_pp_fiberboard_mt(ixC14, ic) = veg_pp_fiberboard_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
            veg_pp_oirw_mt(ixC14, ic) = veg_pp_oirw_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
            veg_pp_pv_mt(ixC14, ic) = veg_pp_pv_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
            veg_pp_sawnwood_mt(ixC14, ic) = veg_pp_sawnwood_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          END IF
        END IF
      END DO

    CASE (QCANOPY)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
      DO ic = 1, nc

       ! update leaf C:N:P pools
       veg_pool_mt(ix_leaf, ixC, ic) = veg_pool_mt(ix_leaf, ixC, ic) + veg_growth_mt(ix_leaf, ixC, ic) &
         &                            - veg_litterfall_mt(ix_leaf, ixC, ic)
       veg_pool_mt(ix_leaf, ixN, ic) = veg_pool_mt(ix_leaf, ixN, ic) + veg_growth_mt(ix_leaf, ixN, ic) &
         &                            - veg_litterfall_mt(ix_leaf, ixN, ic)
       veg_pool_mt(ix_leaf, ixP, ic) = veg_pool_mt(ix_leaf, ixP, ic) + veg_growth_mt(ix_leaf, ixP, ic) &
         &                            - veg_litterfall_mt(ix_leaf, ixP, ic)
  !     veg_pool_mt(ix_leaf, ixN, ic) = veg_pool_mt(ix_leaf, ixC, ic) / lctlib_cn_leaf
  !     veg_pool_mt(ix_leaf, ixP, ic) = veg_pool_mt(ix_leaf, ixN, ic) / lctlib_np_leaf

      END DO
      !$ACC END PARALLEL LOOP

    END SELECT

    !> 1.9 update mean leaf age (advance by one timestep), given growth of new and loss of old leaves
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      IF (veg_pool_mt(ix_leaf, ixC, ic) > eps8) THEN
        mean_leaf_age(ic) = (mean_leaf_age(ic) + dtime / one_day) &
          &                * (1._wp - veg_growth_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_leaf, ixC, ic)) &
          &                - mean_leaf_age(ic) * veg_litterfall_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_leaf, ixC, ic)
      ELSE
        mean_leaf_age(ic) = 0.0_wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !> 2.0 bgc_material diagnostics
    !>     sum up veg_pool, veg_growth and veg_litterfall across components for each particular element
    !>     or simply add value (of one/two/few pool variables) to a memory variable for aggregation and output
    !>

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! init pool totals
      veg_pool_total_c(ic)         = 0.0_wp
      veg_pool_total_n(ic)         = 0.0_wp
      veg_pool_total_p(ic)         = 0.0_wp
      veg_pool_total_c13(ic)       = 0.0_wp
      veg_pool_total_c14(ic)       = 0.0_wp
      veg_pool_total_n15(ic)       = 0.0_wp

      IF (config_l_use_product_pools) THEN
        veg_products_total_c(ic)   = 0.0_wp
        veg_products_total_n(ic)   = 0.0_wp
        veg_products_total_p(ic)   = 0.0_wp
        veg_products_total_c13(ic) = 0.0_wp
        veg_products_total_c14(ic) = 0.0_wp
        veg_products_total_n15(ic) = 0.0_wp
      END IF

      ! init flux totals
      veg_growth_total_c(ic) = 0.0_wp
      veg_growth_total_n(ic) = 0.0_wp
      veg_growth_total_p(ic) = 0.0_wp
      veg_growth_total_c13(ic) = 0.0_wp
      veg_growth_total_c14(ic) = 0.0_wp
      veg_growth_total_n15(ic) = 0.0_wp

      veg_litterfall_total_c(ic) = 0.0_wp
      veg_litterfall_total_n(ic) = 0.0_wp
      veg_litterfall_total_p(ic) = 0.0_wp
      veg_litterfall_total_c13(ic) = 0.0_wp
      veg_litterfall_total_c14(ic) = 0.0_wp
      veg_litterfall_total_n15(ic) = 0.0_wp

      ! store component totals for aggregation
      veg_pool_leaf_c(ic)          = veg_pool_mt(ix_leaf, ixC, ic)
      veg_pool_leaf_n(ic)          = veg_pool_mt(ix_leaf, ixN, ic)
      veg_pool_leaf_p(ic)          = veg_pool_mt(ix_leaf, ixP, ic)
      veg_pool_wood_c(ic)          = veg_pool_mt(ix_sap_wood, ixC, ic) + veg_pool_mt(ix_heart_wood, ixC, ic)
      veg_pool_wood_n(ic)          = veg_pool_mt(ix_sap_wood, ixN, ic) + veg_pool_mt(ix_heart_wood, ixN, ic)
      veg_pool_wood_p(ic)          = veg_pool_mt(ix_sap_wood, ixP, ic) + veg_pool_mt(ix_heart_wood, ixP, ic)
      veg_pool_fine_root_c(ic)     = veg_pool_mt(ix_fine_root, ixC, ic)
      veg_pool_fine_root_n(ic)     = veg_pool_mt(ix_fine_root, ixN, ic)
      veg_pool_fine_root_p(ic)     = veg_pool_mt(ix_fine_root, ixP, ic)
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(hlp)
    DO ic = 1, nc
      hlp = 0._wp
      !$ACC LOOP REDUCTION(+: hlp)
      DO ix_comp = 1, nr_of_veg_bgcm_comp
        hlp = hlp + veg_pool_mt(ix_comp,ixC, ic)
      END DO
      veg_pool_total_c(ic) = hlp
    END DO
    !$ACC END PARALLEL

    !TODO: Consider resolving "Complex loop carried dependence of" as above?
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO ix_comp = 1, nr_of_veg_bgcm_comp
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        ! calculate pool totals
!        veg_pool_total_c(ic)         = veg_pool_total_c(ic) + veg_pool_mt(ix_comp,ixC, ic)
        veg_pool_total_n(ic)         = veg_pool_total_n(ic) + veg_pool_mt(ix_comp,ixN, ic)
        veg_pool_total_p(ic)         = veg_pool_total_p(ic) + veg_pool_mt(ix_comp,ixP, ic)
        veg_pool_total_c13(ic)       = veg_pool_total_c13(ic) + veg_pool_mt(ix_comp,ixC13, ic)
        veg_pool_total_c14(ic)       = veg_pool_total_c14(ic) + veg_pool_mt(ix_comp,ixC14, ic)
        veg_pool_total_n15(ic)       = veg_pool_total_n15(ic) + veg_pool_mt(ix_comp,ixN15, ic)

        IF (config_l_use_product_pools) THEN
          veg_products_total_c(ic)   = 0.0_wp
          veg_products_total_n(ic)   = 0.0_wp
          veg_products_total_p(ic)   = 0.0_wp
          veg_products_total_c13(ic) = 0.0_wp
          veg_products_total_c14(ic) = 0.0_wp
          veg_products_total_n15(ic) = 0.0_wp

          IF (is_active_AGR_or_ALCC) THEN
            veg_products_total_c(ic)   = veg_products_total_c(ic)   + veg_pp_crop_mt(ixC,ic)
            veg_products_total_n(ic)   = veg_products_total_n(ic)   + veg_pp_crop_mt(ixN,ic)
            veg_products_total_p(ic)   = veg_products_total_p(ic)   + veg_pp_crop_mt(ixP,ic)
            veg_products_total_c13(ic) = veg_products_total_c13(ic) + veg_pp_crop_mt(ixC13,ic)
            veg_products_total_c14(ic) = veg_products_total_c14(ic) + veg_pp_crop_mt(ixC14,ic)
            veg_products_total_n15(ic) = veg_products_total_n15(ic) + veg_pp_crop_mt(ixN15,ic)
          END IF

          IF (is_active_SYL_or_ALCC) THEN
            veg_products_total_c(ic) = veg_products_total_c(ic)                                &
              &                       + veg_pp_fuel_mt(ixC,ic) + veg_pp_paper_mt(ixC,ic)       &
              &                       + veg_pp_fiberboard_mt(ixC,ic) + veg_pp_oirw_mt(ixC,ic)  &
              &                       + veg_pp_pv_mt(ixC,ic) + veg_pp_sawnwood_mt(ixC,ic)
            veg_products_total_n(ic) = veg_products_total_n(ic)                                &
              &                       + veg_pp_fuel_mt(ixN,ic) + veg_pp_paper_mt(ixN,ic)       &
              &                       + veg_pp_fiberboard_mt(ixN,ic) + veg_pp_oirw_mt(ixN,ic)  &
              &                       + veg_pp_pv_mt(ixN,ic) + veg_pp_sawnwood_mt(ixN,ic)
            veg_products_total_p(ic) = veg_products_total_p(ic)                                &
              &                       + veg_pp_fuel_mt(ixP,ic) + veg_pp_paper_mt(ixP,ic)       &
              &                       + veg_pp_fiberboard_mt(ixP,ic) + veg_pp_oirw_mt(ixP,ic)  &
              &                       + veg_pp_pv_mt(ixP,ic) + veg_pp_sawnwood_mt(ixP,ic)
            veg_products_total_c13(ic) = veg_products_total_c13(ic)                                &
              &                       + veg_pp_fuel_mt(ixC13,ic) + veg_pp_paper_mt(ixC13,ic)       &
              &                       + veg_pp_fiberboard_mt(ixC13,ic) + veg_pp_oirw_mt(ixC13,ic)  &
              &                       + veg_pp_pv_mt(ixC13,ic) + veg_pp_sawnwood_mt(ixC13,ic)
            veg_products_total_c14(ic) = veg_products_total_c14(ic)                                &
              &                       + veg_pp_fuel_mt(ixC14,ic) + veg_pp_paper_mt(ixC14,ic)       &
              &                       + veg_pp_fiberboard_mt(ixC14,ic) + veg_pp_oirw_mt(ixC14,ic)  &
              &                       + veg_pp_pv_mt(ixC14,ic) + veg_pp_sawnwood_mt(ixC14,ic)
            veg_products_total_n15(ic) = veg_products_total_n15(ic)                                &
              &                       + veg_pp_fuel_mt(ixN15,ic) + veg_pp_paper_mt(ixN15,ic)       &
              &                       + veg_pp_fiberboard_mt(ixN15,ic) + veg_pp_oirw_mt(ixN15,ic)  &
              &                       + veg_pp_pv_mt(ixN15,ic) + veg_pp_sawnwood_mt(ixN15,ic)
          END IF
        END IF

        ! calculate flux totals
        veg_growth_total_c(ic) = veg_growth_total_c(ic) + veg_growth_mt(ix_comp, ixC, ic)
        veg_growth_total_n(ic) = veg_growth_total_n(ic) + veg_growth_mt(ix_comp, ixN, ic)
        veg_growth_total_p(ic) = veg_growth_total_p(ic) + veg_growth_mt(ix_comp, ixP, ic)
        veg_growth_total_c13(ic) = veg_growth_total_c13(ic) + veg_growth_mt(ix_comp, ixC13, ic)
        veg_growth_total_c14(ic) = veg_growth_total_c14(ic) + veg_growth_mt(ix_comp, ixC14, ic)
        veg_growth_total_n15(ic) = veg_growth_total_n15(ic) + veg_growth_mt(ix_comp, ixN15, ic)

        veg_litterfall_total_c(ic) = veg_litterfall_total_c(ic) + veg_litterfall_mt(ix_comp, ixC, ic)
        veg_litterfall_total_n(ic) = veg_litterfall_total_n(ic) + veg_litterfall_mt(ix_comp, ixN, ic)
        veg_litterfall_total_p(ic) = veg_litterfall_total_p(ic) + veg_litterfall_mt(ix_comp, ixP, ic)
        veg_litterfall_total_c13(ic) = veg_litterfall_total_c13(ic) + veg_litterfall_mt(ix_comp, ixC13, ic)
        veg_litterfall_total_c14(ic) = veg_litterfall_total_c14(ic) + veg_litterfall_mt(ix_comp, ixC14, ic)
        veg_litterfall_total_n15(ic) = veg_litterfall_total_n15(ic) + veg_litterfall_mt(ix_comp, ixN15, ic)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! adding seed bed pool for reporting and carbon conservation test to total vegetation pool
      veg_pool_total_c(ic)   = veg_pool_total_c(ic) + seed_bed_pool_mt(ixC, ic)
      veg_pool_total_n(ic)   = veg_pool_total_n(ic) + seed_bed_pool_mt(ixN, ic)
      veg_pool_total_p(ic)   = veg_pool_total_p(ic) + seed_bed_pool_mt(ixP, ic)
      veg_pool_total_c13(ic) = veg_pool_total_c13(ic) + seed_bed_pool_mt(ixC13, ic)
      veg_pool_total_c14(ic) = veg_pool_total_c14(ic) + seed_bed_pool_mt(ixC14, ic)
      veg_pool_total_n15(ic) = veg_pool_total_n15(ic) + seed_bed_pool_mt(ixN15, ic)

      ! adding seed bed growth for reporting and carbon conservation test to total vegetation growth
      veg_growth_total_c(ic)   = veg_growth_total_c(ic) + seed_bed_growth_mt(ixC, ic)
      veg_growth_total_n(ic)   = veg_growth_total_n(ic) + seed_bed_growth_mt(ixN, ic)
      veg_growth_total_p(ic)   = veg_growth_total_p(ic) + seed_bed_growth_mt(ixP, ic)
      veg_growth_total_c13(ic) = veg_growth_total_c13(ic) + seed_bed_growth_mt(ixC13, ic)
      veg_growth_total_c14(ic) = veg_growth_total_c14(ic) + seed_bed_growth_mt(ixC14, ic)
      veg_growth_total_n15(ic) = veg_growth_total_n15(ic) + seed_bed_growth_mt(ixN15, ic)

      ! adding seed bed litterfall for reporting and carbon conservation test to total vegetation litterfall
      veg_litterfall_total_c(ic)   = veg_litterfall_total_c(ic) + seed_bed_litterfall_mt(ixC, ic)
      veg_litterfall_total_n(ic)   = veg_litterfall_total_n(ic) + seed_bed_litterfall_mt(ixN, ic)
      veg_litterfall_total_p(ic)   = veg_litterfall_total_p(ic) + seed_bed_litterfall_mt(ixP, ic)
      veg_litterfall_total_c13(ic) = veg_litterfall_total_c13(ic) + seed_bed_litterfall_mt(ixC13, ic)
      veg_litterfall_total_c14(ic) = veg_litterfall_total_c14(ic) + seed_bed_litterfall_mt(ixC14, ic)
      veg_litterfall_total_n15(ic) = veg_litterfall_total_n15(ic) + seed_bed_litterfall_mt(ixN15, ic)
    END DO
    !$ACC END PARALLEL LOOP

    IF (config_l_use_product_pools) THEN
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
      DO ic = 1, nc
        veg_products_decay_c(ic)   =  prod_decay_mt(ixC, ic)   * 1000000.0_wp / dtime
        veg_products_decay_n(ic)   =  prod_decay_mt(ixN, ic)   * 1000000.0_wp / dtime
        veg_products_decay_p(ic)   =  prod_decay_mt(ixP, ic)   * 1000000.0_wp / dtime
        veg_products_decay_c13(ic) =  prod_decay_mt(ixC13, ic) * 1000000.0_wp / dtime
        veg_products_decay_c14(ic) =  prod_decay_mt(ixC14, ic) * 1000000.0_wp / dtime
        veg_products_decay_n15(ic) =  prod_decay_mt(ixN15, ic) * 1000000.0_wp / dtime
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !> 3.0 biosphere-level diagnostics across multiple processes
    !>
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      ! see also 'update_sb_pools()'
      !
      ! TODO calculation may include 'fFire' and 'fLUC' (land-use change) once available
      !
      ! could also be calculated based on NPP
      ! (n_processing_respiration = n_fixation_respiration + n_transform_respiration)
      net_biosphere_production(ic) = net_biosphere_production(ic) &
        & + gross_assimilation(ic) &
        & - maint_respiration(ic)   - growth_respiration(ic) - n_processing_respiration(ic) &
        & - herbivory_leaf_resp(ic) - herbivory_fruit_resp(ic)
      net_biosphere_production_c13(ic) = net_biosphere_production_c13(ic) &
        & + gross_assimilation_C13(ic) &
        & - maint_respiration_c13(ic)   - growth_respiration_c13(ic) - n_processing_respiration_c13(ic) &
        & - herbivory_leaf_resp_c13(ic) - herbivory_fruit_resp_c13(ic)
      net_biosphere_production_c14(ic) = net_biosphere_production_c14(ic) &
        & + gross_assimilation_C14(ic) &
        & - maint_respiration_c14(ic)   - growth_respiration_c14(ic) - n_processing_respiration_c14(ic) &
        & - herbivory_leaf_resp_c14(ic) - herbivory_fruit_resp_c14(ic)

      ! biological N fixation (VEG_ and SB_ processes)
      biological_n_fixation(ic)    = biological_n_fixation(ic) + n_fixation(ic)

      IF (config_l_use_product_pools) THEN
        ! The product pool decay is a loss to the atmosphere (prod_decay_mt unit is mol m-2 timestep-1)
        net_biosphere_production(ic)     = net_biosphere_production(ic)     - (prod_decay_mt(ixC,ic) * 1000000.0_wp / dtime)
        net_biosphere_production_c13(ic) = net_biosphere_production_c13(ic) - (prod_decay_mt(ixC13,ic) * 1000000.0_wp / dtime)
        net_biosphere_production_c14(ic) = net_biosphere_production_c14(ic) - (prod_decay_mt(ixC14,ic) * 1000000.0_wp / dtime)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE update_veg_pools

#endif
END MODULE mo_q_veg_update_pools
