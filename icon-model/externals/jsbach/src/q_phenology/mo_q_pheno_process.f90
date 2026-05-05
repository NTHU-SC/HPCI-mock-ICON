!> QUINCY calculate vegetation phenology
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
!>#### determine onset and end of the growing season across plant functional types (PFT)
!>
MODULE mo_q_pheno_process
#ifndef __NO_QUINCY__

  USE mo_kind,                 ONLY: wp
  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class, ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_LITTERFALL_ID
  USE mo_jsb_impl_constants,   ONLY: true, false, test_false_true
  USE mo_jsb_control,          ONLY: debug_on
  USE mo_exception,            ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_q_phenology

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_process'

CONTAINS

  ! ======================================================================================================= !
  !> update phenology - QUINCY
  !>   determine onset and end of the growing season for some PFT
  !>
  !>   PFT: evergreen, raingreen, cold-deciduous trees and perennial grasses
  !>
  !>  @TODO
  !>  OBS: currently uncalibrated, in particular rain green phenology
  !>  Dev: Have more cost-benefit in decision to start/stop growing season
  !>      what about seasonality beyond off/on behaviour? -> implies some collaboration between
  !>      turnover rates and phenology -> should the two be linked?
  SUBROUTINE update_q_phenology(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_quincy_model_config,     ONLY: QLAND, QPLANT, QCANOPY
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, Q_PHENO_, Q_AGR_, Q_ASSIMI_, VEG_
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_AGR_)
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                  !< the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER   :: lctlib                 !< land-cover-type library - parameter across pft's
    INTEGER                           :: iblk, ics, ice, nc, ic !< dimensions
    REAL(wp)                          :: dtime                  !< timestep length

    INTEGER  :: lctlib_phenology_type
    REAL(wp) :: lctlib_beta_soil_flush, lctlib_sla, lctlib_cn_leaf, lctlib_np_leaf, lctlib_tau_leaf, &
      &         lctlib_gdd_req_max, lctlib_k_gdd_dormance, lctlib_beta_soil_senescence, lctlib_t_air_senescence, &
      &         lctlib_min_leaf_age

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_q_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk      :: t_air
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: gdd
    dsl4jsb_Real2D_onChunk      :: nd_dormance
    dsl4jsb_Real2D_onChunk      :: root_phenology_type
    dsl4jsb_Real2D_onChunk      :: lai_max
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: t_soil_srf_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_PHENO_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN !< run task only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    dsl4jsb_Get_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%CropFlag .AND. dsl4jsb_Config(Q_AGR_)%active) RETURN !< do not run this routine for crops if
                                                                    !< the crop process is active
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    ! A2L_ 2D
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                            ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)               ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, gdd)                          ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, nd_dormance)                  ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, root_phenology_type)          ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, lai_max)                      ! inout
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tphen_mavg)     ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                              ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_soil_srf_tphen_mavg)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)        ! inout

    lctlib_phenology_type = dsl4jsb_Lctlib_param(phenology_type)
    lctlib_beta_soil_flush = dsl4jsb_Lctlib_param(beta_soil_flush)
    lctlib_sla = dsl4jsb_Lctlib_param(sla)
    lctlib_cn_leaf = dsl4jsb_Lctlib_param(cn_leaf)
    lctlib_np_leaf = dsl4jsb_Lctlib_param(np_leaf)
    lctlib_tau_leaf = dsl4jsb_Lctlib_param(tau_leaf)
    lctlib_gdd_req_max = dsl4jsb_Lctlib_param(gdd_req_max)
    lctlib_k_gdd_dormance = dsl4jsb_Lctlib_param(k_gdd_dormance)
    lctlib_beta_soil_senescence = dsl4jsb_Lctlib_param(beta_soil_senescence)
    lctlib_t_air_senescence = dsl4jsb_Lctlib_param(t_air_senescence)
    lctlib_min_leaf_age = dsl4jsb_Lctlib_param(min_leaf_age)

    !> 1.0 Check for start and end of growing season, and update both gdd and nd_dormance counter
    !>
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      CALL calc_phenology(                      &
        & dtime                               , &
        & lctlib_phenology_type               , &  ! lctlib
        & lctlib_beta_soil_flush              , &
        & lctlib_cn_leaf                      , &
        & lctlib_np_leaf                      , &
        & lctlib_gdd_req_max                  , &
        & lctlib_k_gdd_dormance               , &
        & lctlib_beta_soil_senescence         , &
        & lctlib_t_air_senescence             , &
        & lctlib_min_leaf_age                 , &  ! lctlib
        & t_air(ic)                           , &  ! in
        & beta_soil_gs_tphen_mavg(ic)         , &
        & t_air_tphen_mavg(ic)                , &
        & t_soil_srf_tphen_mavg(ic)           , &
        & t_air_week_mavg(ic)                 , &
        & t_air_month_mavg(ic)                , &
        & mean_leaf_age(ic)                   , &
        & gpp_tlabile_mavg(ic)                , &
        & maint_respiration_tlabile_mavg(ic)  , &  ! in
        & growing_season(ic)                  , &  ! inout
        & growth_req_n_tlabile_mavg(ic)       , &
        & growth_req_p_tlabile_mavg(ic)       , &
        & nd_dormance(ic)                     , &
        & gdd(ic)                             , &
        & root_phenology_type(ic)               )  ! inout
    END DO
    !$ACC END PARALLEL LOOP

    !> 2.0 update growth and litterfall of leaves for CANOPY model only, so that in Q_VEG:mo_q_veg_update_pools
    !>     the leaf pool is updated accordingly, and thus Q_VEG:mo_q_veg_plant_characteristics derives the correct
    !>     LAI
    !>
    SELECT CASE (model%config%qmodel_id)
      !$ACC UPDATE HOST(growing_season(:),lai(:),gdd(:), lai_max(:))&
      !$ACC   HOST(veg_pool_mt(:,:,:),veg_growth_mt(:,:,:),veg_litterfall_mt(:, :, :)) ASYNC(1)
      !$ACC WAIT(1)
      CASE (QCANOPY)
        CALL calc_leaf_mass_change_canopy_mode(             &
          & nc                                , &
          & dtime                             , &
          & lctlib_phenology_type             , & ! lctlib
          & lctlib_sla                        , &
          & lctlib_tau_leaf                   , &
          & lctlib_cn_leaf                    , &
          & lctlib_np_leaf                    , & ! lctlib
          & growing_season(:)                 , & ! in
          & lai(:)                            , &
          & lai_max(:)                        , &
          & veg_pool_mt(ix_leaf, ixC, :)      , & ! in
          & veg_growth_mt(ix_leaf, ixC, :)    , & ! inout
          & veg_growth_mt(ix_leaf, ixN, :)    , &
          & veg_growth_mt(ix_leaf, ixP, :)    , &
          & veg_litterfall_mt(ix_leaf, ixC, :), &
          & veg_litterfall_mt(ix_leaf, ixN, :), &
          & veg_litterfall_mt(ix_leaf, ixP, :)  ) ! inout
      !$ACC UPDATE DEVICE(veg_growth_mt(:,:,:),veg_litterfall_mt(:,:,:)) ASYNC(1)
    END SELECT

  END SUBROUTINE update_q_phenology

  ! ======================================================================================================= !
  !> calculate phenology
  !>
  ELEMENTAL SUBROUTINE calc_phenology(  &
    & dtime                           , &
    & phenology_type                  , &
    & beta_soil_flush                 , &
    & cn_leaf                         , &
    & np_leaf                         , &
    & gdd_req_max                     , &
    & k_gdd_dormance                  , &
    & beta_soil_senescence            , &
    & t_air_senescence                , &
    & min_leaf_age                    , &
    & t_air                           , &
    & beta_soil_gs_tphen_mavg         , &
    & t_air_tphen_mavg                , &
    & t_soil_srf_tphen_mavg           , &
    & t_air_week_mavg                 , &
    & t_air_month_mavg                , &
    & mean_leaf_age                   , &
    & gpp_tlabile_mavg                , &
    & maint_respiration_tlabile_mavg  , &
    & growing_season                  , &
    & growth_req_n_tlabile_mavg       , &
    & growth_req_p_tlabile_mavg       , &
    & nd_dormance                     , &
    & gdd                             , &
    & root_phenology_type)

    !$ACC ROUTINE SEQ

    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: true, false, test_false_true
    USE mo_q_pheno_constants,       ONLY: ievergreen, iraingreen, isummergreen, iperennial, icrop_phenology, &
      &                                   ipheno_type_cold_deciduous, ipheno_type_drought_deciduous, &
      &                                   ipheno_type_cbalance_deciduous
    USE mo_q_pheno_parameters,      ONLY: gdd_t_air_threshold
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),     INTENT(in)    :: dtime                            !< timestep length
    INTEGER,      INTENT(in)    :: phenology_type                   !< lctlib parameter
    REAL(wp),     INTENT(in)    :: beta_soil_flush                  !< lctlib parameter
    REAL(wp),     INTENT(in)    :: cn_leaf                          !< lctlib parameter
    REAL(wp),     INTENT(in)    :: np_leaf                          !< lctlib parameter
    REAL(wp),     INTENT(in)    :: gdd_req_max                      !< lctlib parameter
    REAL(wp),     INTENT(in)    :: k_gdd_dormance                   !< lctlib parameter
    REAL(wp),     INTENT(in)    :: beta_soil_senescence             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: t_air_senescence                 !< lctlib parameter
    REAL(wp),     INTENT(in)    :: min_leaf_age                     !< lctlib parameter
    REAL(wp),     INTENT(in)    :: t_air                            !< air temp
    REAL(wp),     INTENT(in)    :: beta_soil_gs_tphen_mavg          !< soil ..
    REAL(wp),     INTENT(in)    :: t_air_tphen_mavg                 !< moving average air temperature over phenology-relevant period
    REAL(wp),     INTENT(in)    :: t_soil_srf_tphen_mavg            !< moving average soil surface temperature
    REAL(wp),     INTENT(in)    :: t_air_week_mavg                  !< moving average air temperature over a week
    REAL(wp),     INTENT(in)    :: t_air_month_mavg                 !< moving average air temperature over a month
    REAL(wp),     INTENT(in)    :: mean_leaf_age                    !< mean leaf age
    REAL(wp),     INTENT(in)    :: gpp_tlabile_mavg                 !< moving average gpp over tlabile
    REAL(wp),     INTENT(in)    :: maint_respiration_tlabile_mavg   !< moving average maintenance respiration over tlabile
    REAL(wp),     INTENT(inout) :: growing_season                   !< growing season Yes/No
    REAL(wp),     INTENT(inout) :: growth_req_n_tlabile_mavg        !< moles N required for a unit of C growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: growth_req_p_tlabile_mavg        !< moles P required for a unit of N growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: nd_dormance                      !< number of days in dormancy since last growing season
    REAL(wp),     INTENT(inout) :: gdd                              !< growing degree days
    REAL(wp),     INTENT(inout) :: root_phenology_type              !< category for trigger of end of growing season in grasses
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: hlp1                             !< helper variable for tmp values
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phenology'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 check growing season has started or ended
    !>
    ! If not in growing season, check whether growing season has started
    IF(growing_season < test_false_true) THEN
      SELECT CASE (phenology_type)
      CASE (ievergreen)
        ! nothing to do here
      CASE (iraingreen)
        ! if soil moisture is above specified threshold
        IF(beta_soil_gs_tphen_mavg > beta_soil_flush) THEN
          growing_season            = true
          growth_req_n_tlabile_mavg = 1._wp / cn_leaf
          growth_req_p_tlabile_mavg = 1._wp / np_leaf
        ENDIF
      CASE (isummergreen)
        ! if GDD requirement is exceeded. The latter depends on the length of the
        !  dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
        hlp1 = gdd_req_max * EXP( - k_gdd_dormance * nd_dormance )
        IF(gdd > hlp1) THEN
           growing_season            = true
           growth_req_n_tlabile_mavg = 1._wp / cn_leaf
           growth_req_p_tlabile_mavg = 1._wp / np_leaf
        ENDIF
      CASE (iperennial,icrop_phenology)
        ! if temperature and soil moisture are not limiting growth
        hlp1 = gdd_req_max * EXP( - k_gdd_dormance * nd_dormance )
        IF(gdd > hlp1 .AND. beta_soil_gs_tphen_mavg > beta_soil_flush) THEN
          growing_season            = true
          growth_req_n_tlabile_mavg = 1._wp / cn_leaf
          growth_req_p_tlabile_mavg = 1._wp / np_leaf
        ENDIF
      END SELECT
    ! If in the growing season, check whether growing season has ended
    ELSE
      SELECT CASE (phenology_type)
      CASE (ievergreen)
        ! nothing to do here
      CASE (iraingreen)
        ! if soil moisture is below given threshold AND leaves have been out a given number of days
        IF(beta_soil_gs_tphen_mavg < beta_soil_senescence .AND. mean_leaf_age > min_leaf_age ) THEN
          growing_season = false
        ENDIF
      CASE (isummergreen)
        ! if temperature are declining below a certain threshold AND
        !  leaves have been out a given number of days
        IF(t_air_week_mavg < t_air_month_mavg .AND. t_air_tphen_mavg < t_air_senescence .AND. &
          & mean_leaf_age > min_leaf_age) THEN
          growing_season = false
        ENDIF
      CASE (iperennial,icrop_phenology)
        ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND
        !  leaves have been out a given number of days
        IF ((t_air_tphen_mavg < t_air_senescence .OR. beta_soil_gs_tphen_mavg < beta_soil_senescence .OR. &
            & gpp_tlabile_mavg < maint_respiration_tlabile_mavg).AND. mean_leaf_age > min_leaf_age) THEN
          growing_season = false
          ! Determine the specific trigger for the end of growing season. Needed for root turnover.
          IF (t_air_tphen_mavg < t_air_senescence) THEN
            ! end of growing season triggered by low temperature
            root_phenology_type = ipheno_type_cold_deciduous
          ELSE
            ! end of growing season triggered by low soil moisture
            IF (beta_soil_gs_tphen_mavg < beta_soil_senescence) THEN
              root_phenology_type = ipheno_type_drought_deciduous
            ELSE
              ! end of growing season triggered by negative C balance
              root_phenology_type = ipheno_type_cbalance_deciduous
            END IF
          END IF
        END IF
      END SELECT
    ENDIF

    !>2.0 modify phenology counters
    !>

    !>  2.1 gdd counter
    !>
    IF(t_air > gdd_t_air_threshold) THEN
      gdd = gdd + ((t_air - gdd_t_air_threshold) * dtime / one_day)
    ENDIF
    ! vegetation growth doesn't start when there is snow coverage and soil is frozen
    IF(t_air_tphen_mavg < Tzero .AND. t_soil_srf_tphen_mavg <= Tzero) THEN
      gdd = 0.0_wp
    ENDIF

    !>  2.2 nd_dormance counter
    !>
    IF (growing_season > test_false_true) THEN
      gdd         = 0.0_wp
      nd_dormance = 0.0_wp
    ELSE
      nd_dormance = nd_dormance + dtime / one_day
    ENDIF
  END SUBROUTINE calc_phenology

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_phenology
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate leaf growth and litterfall required to simulate LAI change depending on current
  !! phenological state in QCANOPY mode (i.e. if only fast biogeophysical processes are calculated)
  !!
  !!   Input:
  !!     (1) phenological status, current LAI and maximum LAI, leaf carbon pool
  !!
  !!   Output: leaf carbon, nitrogen and phosphorus growth and litter fall
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_leaf_mass_change_canopy_mode(   &
    & nc                            , &
    & dtime                         , &
    & lctlib_phenology_type         , &
    & lctlib_sla                    , &
    & lctlib_tau_leaf               , &
    & lctlib_cn_leaf                , &
    & lctlib_np_leaf                , &
    & growing_season                , &
    & lai                           , &
    & lai_max                       , &
    & veg_pool_leaf_carbon          , &
    & veg_growth_leaf_carbon        , &
    & veg_growth_leaf_nitrogen      , &
    & veg_growth_leaf_phosphorus    , &
    & veg_litterfall_leaf_carbon    , &
    & veg_litterfall_leaf_nitrogen  , &
    & veg_litterfall_leaf_phosphorus)

    USE mo_jsb_math_constants,            ONLY: one_day,one_year,eps4
    USE mo_veg_constants,                 ONLY: max_leaf_shedding_rate
    USE mo_q_pheno_constants,             ONLY: ievergreen
    USE mo_q_pheno_parameters,            ONLY: k_leafon_canopy
    !------------------------------------------------------------------------------------------------------ !
    INTEGER,      INTENT(in)    :: nc                                !< dimensions
    REAL(wp),     INTENT(in)    :: dtime                             !< timestep length
    INTEGER,      INTENT(in)    :: lctlib_phenology_type             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: lctlib_sla                        !< lctlib parameter
    REAL(wp),     INTENT(in)    :: lctlib_tau_leaf                   !< lctlib parameter
    REAL(wp),     INTENT(in)    :: lctlib_cn_leaf                    !< lctlib parameter
    REAL(wp),     INTENT(in)    :: lctlib_np_leaf                    !< lctlib parameter
    REAL(wp),     INTENT(in)    :: growing_season(:)                 !< does the plant intend to grow?
    REAL(wp),     INTENT(in)    :: lai(:)                            !< current LAI (m2/m2)
    REAL(wp),     INTENT(in)    :: lai_max(:)                        !< maximum LAI (m2/m2)
    REAL(wp),     INTENT(in)    :: veg_pool_leaf_carbon(:)           !< the plant's leaf carbon pool (mol/m2)
    REAL(wp),     INTENT(inout) :: veg_growth_leaf_carbon(:)         !< the plant's current leaf growth rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_growth_leaf_nitrogen(:)       !< the plant's current leaf growth rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_growth_leaf_phosphorus(:)     !< the plant's current leaf growth rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_litterfall_leaf_carbon(:)     !< the plant's current litter fall rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_litterfall_leaf_nitrogen(:)   !< the plant's current litter fall rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_litterfall_leaf_phosphorus(:) !< the plant's current litter fall rate (mol/time step)
    !------------------------------------------------------------------------------------------------------ !
    INTEGER                     :: ic                       !< loop over point of the chunk
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_leaf_bgcm_change'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1, nc

      ! diagnose required change in leaf C given phenology to attain prescribed LAI max
       IF (growing_season(ic) > test_false_true) THEN
         veg_growth_leaf_carbon(ic) = lai_max(ic) / lctlib_sla &
           &                          * (1._wp - exp(-k_leafon_canopy &
           &                          * (lai_max(ic) - lai(ic)) / lai_max(ic))) &
           &                          * dtime / one_day
       ELSE
         veg_litterfall_leaf_carbon(ic) = MIN(veg_pool_leaf_carbon(ic), &
           &                                  lai_max(ic)/lctlib_sla * max_leaf_shedding_rate &
           &                                  * dtime / one_day)
       END IF

       ! calc leaf turnover for evergreen trees (needed, for e.g. mean_leaf_age, as evergreens are always within growing_season)
       IF (lctlib_phenology_type == ievergreen) THEN
         veg_litterfall_leaf_carbon(ic) = veg_pool_leaf_carbon(ic) / lctlib_tau_leaf / one_day / one_year * dtime
       END IF

       ! update growth and litterfall of leaf N:P pools
       veg_growth_leaf_nitrogen(ic) = veg_growth_leaf_carbon(ic) / lctlib_cn_leaf
       veg_growth_leaf_phosphorus(ic) = veg_growth_leaf_nitrogen(ic) / lctlib_np_leaf

       veg_litterfall_leaf_nitrogen(ic) = veg_litterfall_leaf_carbon(ic) / lctlib_cn_leaf
       veg_litterfall_leaf_phosphorus(ic) = veg_litterfall_leaf_nitrogen(ic) / lctlib_np_leaf

    END DO

  END SUBROUTINE calc_leaf_mass_change_canopy_mode

#endif
END MODULE mo_q_pheno_process
