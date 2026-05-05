!> QUINCY calculate vegetation dynamics
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
!>#### calculate vegetation dynamics, e.g., establishment and mortality
!>
MODULE mo_q_veg_dynamics
#ifndef __NO_QUINCY__

  USE mo_kind,                 ONLY: wp
  USE mo_jsb_impl_constants,   ONLY: true, false, test_false_true
  USE mo_jsb_control,          ONLY: debug_on
  USE mo_exception,            ONLY: message
  USE mo_jsb_math_constants,   ONLY: one_year, one_day, eps8, eps4

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_LITTERFALL_ID, VEG_BGCM_ESTABLISHMENT_ID, &
    &                                   VEG_BGCM_SEED_BED_POOL_ID, VEG_BGCM_SEED_BED_GROWTH_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_veg_dynamics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_dynamics'

CONTAINS

  ! ======================================================================================================= !
  !>Calculate establishment and mortality of a population
  !>
  !>  (currently static mortality rate only)
  !>  a background_mort_rate exists for each trees and grasses, defined as mortality rate per year
  !>
  SUBROUTINE update_veg_dynamics(tile, options)

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: VEG_, Q_AGR_, Q_PHENO_, HYDRO_, Q_SYL_
    USE mo_veg_constants,         ONLY: itree, rfr_ratio_toc, k_r2fr_chl

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(Q_AGR_)
    dsl4jsb_Use_config(Q_SYL_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_lctlib_element),   POINTER :: lctlib                 !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)   :: fpc                    !< current foliage projective cover
    REAL(wp), DIMENSION(options%nc)   :: hlp1                   !< dummy variable
    REAL(wp), DIMENSION(options%nc)   :: fract_herbivory        !< fraction of leaf and fruit (carbon) lost to herbivory
    REAL(wp)                          :: dtime                  !< timestep length
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    INTEGER                           :: ic                     !< loop over chunk
    INTEGER                           :: ielem, ix_elem         !< id and index of elements
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_dynamics'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_pool_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_growth_mt
    dsl4jsb_Def_mt1L2D :: veg_establishment_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_config(Q_SYL_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onChunk      :: wtr_plant_avail_rel
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: mortality_rate
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: lai_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp_c13
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp_c14
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp_c13
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp_c14
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_config(Q_SYL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_POOL_ID, seed_bed_pool_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_GROWTH_ID, seed_bed_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_ESTABLISHMENT_ID, veg_establishment_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)         ! in
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_plant_avail_rel)      ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth_tvegdyn_mavg)    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, lai_tvegdyn_mavg)           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, rfr_ratio_boc)              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, rfr_ratio_boc_tvegdyn_mavg) ! in
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)           ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                  ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)           ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)      ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)     ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)           ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)      ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)     ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)         ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)    ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)   ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp)        ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp_c13)    ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp_c14)    ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp)       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp_c13)   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp_c14)   ! out
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.9 init local var
    !>
    fract_herbivory(:)          = 0.0_wp

    !>1.0 Calculate foliage projective cover if vegetation dynamics are simulated, i.e., not needed for constant mortality
    !>
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
      CASE ("population")
        fpc(:) = calc_foliage_projective_cover_tile(nc, &
          &                                         lctlib%growthform, &
          &                                         lctlib%wood_density, &
          &                                         lctlib%k_latosa, &
          &                                         veg_pool_mt(ix_sap_wood, ixC, :), &
          &                                         lai(:), &
          &                                         dens_ind(:), &
          &                                         height(:), &
          &                                         diameter(:))
    ENDSELECT

    !>2.0 Calculate mortality rate and veg litterfall
    !>
    ! returns constant mortality_rate if veg_dynamics_scheme == none \n
    ! but calcs mortality_rate dynamically if veg_dynamics_scheme == population
    ! does not apply to croplands
    IF (lctlib%CropFlag .AND. dsl4jsb_Config(Q_AGR_)%active) THEN
      mortality_rate(:) = 0.0_wp
    ELSE
      CALL calc_veg_mortality( &
        & nc                                              , & ! in
        & dtime                                           , &
        & model%config%elements_index_map(:)              , &
        & model%config%is_element_used(:)                 , &
        & lctlib%growthform                               , &
        & lctlib%k1_mort_greff                            , &
        & TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme)  , &
        & an_boc_tvegdyn_mavg(:)                          , &
        & net_growth_tvegdyn_mavg(:)                      , &
        & lai_tvegdyn_mavg(:)                             , &
        & fpc(:)                                          , &
        & veg_pool_mt(:,:,:)                              , & ! in
        & veg_litterfall_mt(:,:,:)                        , & ! inout
        & mortality_rate(:)                               )   ! out
    END IF

    !> 2.1 Calculate herbivory and associated respiration and litter fall
    !>
    IF (dsl4jsb_Config(VEG_)%flag_herbivory) THEN
      CALL calc_veg_herbivory( &
        & nc, &                                 ! in
        & dtime, &                              ! in
        & model%config%elements_index_map(:), & ! in
        & model%config%is_element_used(:), &    ! in
        & lctlib%GrassFlag, &                   ! in
        & lctlib%PastureFlag, &                 ! in
        & lctlib%cn_leaf_min, &                 ! in
        & lai(:), &                             ! in
        & net_growth_tvegdyn_mavg(:), &         ! in
        & mortality_rate(:), &                  ! in
        & veg_pool_mt(:,:,:), &                 ! in
        & herbivory_leaf_resp(:), &             ! inout
        & herbivory_leaf_resp_c13(:), &         ! inout
        & herbivory_leaf_resp_c14(:), &         ! inout
        & herbivory_fruit_resp(:), &            ! inout
        & herbivory_fruit_resp_c13(:), &        ! inout
        & herbivory_fruit_resp_c14(:), &        ! inout
        & fract_herbivory(:), &                 ! inout
        & veg_litterfall_mt(:,:,:))             ! inout
    END IF

    !> 2.2 Adjust tissue turnover that generates new tissue by mortality
    !>
    !> Tissues that turnover to become new tissues (sapwood -> heartwood; fruit -> seed bed)
    !> are recorded as growth (or negative growth) by the turnover routine
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (model%config%is_element_used(ielem)) THEN
        ix_elem = model%config%elements_index_map(ielem)    ! get element index in bgcm
        veg_growth_mt(ix_sap_wood,ix_elem,:)   = veg_growth_mt(ix_sap_wood,ix_elem,:)   * (1._wp - mortality_rate(:))
        veg_growth_mt(ix_heart_wood,ix_elem,:) = veg_growth_mt(ix_heart_wood,ix_elem,:) * (1._wp - mortality_rate(:))
        veg_growth_mt(ix_fruit,ix_elem,:)      = veg_growth_mt(ix_fruit,ix_elem,:)      * (1._wp - mortality_rate(:))
        seed_bed_growth_mt(ix_elem,:)          = seed_bed_growth_mt(ix_elem,:)          * (1._wp - mortality_rate(:))
      END IF
    END DO
    !>
    !> Retranslocation of foliar nutrients is recorded separately, as it applies only to nutrients
    !> These have to be reduced by the fraction of individuals that have died (i.e. no heartwood, seedpool
    !> generation or retranslocation by dying individuals)
    recycling_leaf_n(:)         = recycling_leaf_n(:)         * (1._wp - mortality_rate(:) - fract_herbivory(:))
    recycling_fine_root_n(:)    = recycling_fine_root_n(:)    * (1._wp - mortality_rate(:))
    recycling_heart_wood_n(:)   = recycling_heart_wood_n(:)   * (1._wp - mortality_rate(:))
    recycling_leaf_p(:)         = recycling_leaf_p(:)         * (1._wp - mortality_rate(:) - fract_herbivory(:))
    recycling_fine_root_p(:)    = recycling_fine_root_p(:)    * (1._wp - mortality_rate(:))
    recycling_heart_wood_p(:)   = recycling_heart_wood_p(:)   * (1._wp - mortality_rate(:))
    recycling_leaf_n15(:)       = recycling_leaf_n15(:)       * (1._wp - mortality_rate(:) - fract_herbivory(:))
    recycling_fine_root_n15(:)  = recycling_fine_root_n15(:)  * (1._wp - mortality_rate(:))
    recycling_heart_wood_n15(:) = recycling_heart_wood_n15(:) * (1._wp - mortality_rate(:))

    !>
    !> catch numerical issues arising at leaf or fine-root shedding of the order of < E-12
    DO ic = 1, nc
      IF(veg_pool_mt(ix_leaf,ixN,ic) < recycling_leaf_n(ic) + veg_litterfall_mt(ix_leaf,ixN,ic)) THEN
        recycling_leaf_n(ic) = veg_pool_mt(ix_leaf,ixN,ic) - veg_litterfall_mt(ix_leaf,ixN,ic)
      END IF
      IF(veg_pool_mt(ix_leaf,ixP,ic) < recycling_leaf_p(ic) + veg_litterfall_mt(ix_leaf,ixP,ic)) THEN
        recycling_leaf_p(ic) = veg_pool_mt(ix_leaf,ixP,ic) - veg_litterfall_mt(ix_leaf,ixP,ic)
      END IF
      IF(veg_pool_mt(ix_leaf,ixN15,ic) < recycling_leaf_n15(ic) + veg_litterfall_mt(ix_leaf,ixN15,ic)) THEN
        recycling_leaf_n15(ic) = veg_pool_mt(ix_leaf,ixN15,ic) - veg_litterfall_mt(ix_leaf,ixN15,ic)
      END IF
      IF(veg_pool_mt(ix_fine_root,ixN,ic) < recycling_fine_root_n(ic) + veg_litterfall_mt(ix_fine_root,ixN,ic)) THEN
        recycling_fine_root_n(ic) = veg_pool_mt(ix_fine_root,ixN,ic) - veg_litterfall_mt(ix_fine_root,ixN,ic)
      END IF
      IF(veg_pool_mt(ix_fine_root,ixP,ic) < recycling_fine_root_p(ic) + veg_litterfall_mt(ix_fine_root,ixP,ic)) THEN
        recycling_fine_root_p(ic) = veg_pool_mt(ix_fine_root,ixP,ic) - veg_litterfall_mt(ix_fine_root,ixP,ic)
      END IF
      IF(veg_pool_mt(ix_fine_root,ixN15,ic) < recycling_fine_root_n15(ic) + veg_litterfall_mt(ix_fine_root,ixN15,ic)) THEN
        recycling_fine_root_n15(ic) = veg_pool_mt(ix_fine_root,ixN15,ic) - veg_litterfall_mt(ix_fine_root,ixN15,ic)
      END IF
    END DO

    !>3.0 Calculate establishment_rate and veg establishment
    !>
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
    CASE ("population")
      ! Determine establishment flux:

      ! Standard case: calculate establishment from seed pool
      ! (Exceptions: no establishment on cropland or if we just had a stand replacing harvest running QS)
      IF (.NOT. (lctlib%CropFlag .AND. dsl4jsb_Config(Q_AGR_)%active)) THEN
        IF (.NOT. dsl4jsb_Config(Q_SYL_)%flag_stand_harvest_event) THEN
          CALL calc_veg_establishment( &
            & nc                                  , & ! in
            & dtime                               , &
            & model%config%elements_index_map(:)  , &
            & model%config%is_element_used(:)     , &
            & lctlib%tau_seed_est                 , &
            & wtr_plant_avail_rel(:)              , &
            & growing_season(:)                   , &
            & t_air_week_mavg(:)                  , &
            & fpc(:)                              , &
            & rfr_ratio_boc_tvegdyn_mavg(:)       , &
            & seed_bed_pool_mt(:,:)               , & ! in
            & veg_establishment_mt(:,:)           )   ! inout
        ELSE
          ! Stand replacing harvest in QS is an exceptional case in which the establishment is initialised in a particular way
          hlp1(:) = 0.0_wp
          ! set carbon in establishment such that the inital stand density is 1000 / ha (i.e. 0.1/m2)
          ! which seeds the forest floor for natural regeneration
          ! but limit to actual available seed pool to avoid generating mass
          ! NB: This might need to be revised for particular case studies
          DO ic = 1, nc
            IF (seed_bed_pool_mt(ixC, ic) > eps8) THEN
              hlp1(ic) = MIN(0.1_wp * lctlib%seed_size, seed_bed_pool_mt(ixC, ic)) / seed_bed_pool_mt(ixC, ic)
            END IF
          END DO
          DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
            IF (model%config%is_element_used(ielem)) THEN
              ix_elem = model%config%elements_index_map(ielem)    ! get element index in bgcm
              veg_establishment_mt(ix_elem, :) = hlp1(:) * seed_bed_pool_mt(ix_elem, :)
            END IF
          END DO
        END IF
      END IF
    END SELECT

    !>4.0 Calculate change in individual density from mortality and establishment
    !>
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))

    CASE ("population")
      IF (lctlib%growthform == itree) THEN
        delta_dens_ind(:) = calc_delta_dens_ind(lctlib%seed_size            , &
          &                                     dens_ind(:)                 , &
          &                                     mortality_rate(:)           , &
          &                                     veg_establishment_mt(ixC,:) )
      ELSE
        delta_dens_ind(:) = 0.0_wp
      END IF
    END SELECT

    !> 5.0 calculate forest floor red to far-red ratio
    !>
    rfr_ratio_boc(:) = rfr_ratio_toc * &
      &                EXP(k_r2fr_chl * SUM(leaf_nitrogen_cl(:,:) * fn_chl_cl(:,:) * lai_cl(:,:), DIM=2))

  END SUBROUTINE update_veg_dynamics

  ! ======================================================================================================= !
  !>calculates current foliage projective cover for a tile
  !>
  !>  follows Sitch et al. 2003, eq 8
  !>
  FUNCTION calc_foliage_projective_cover_tile( &
    & nc                      , &
    & lctlib_growthform       , &
    & lctlib_wood_density     , &
    & lctlib_k_latosa         , &
    & veg_pool_sap_wood_carbon, &
    & lai                     , &
    & dens_ind                , &
    & height                  , &
    & diameter)                 &
    & RESULT(fpc)

    USE mo_veg_constants,         ONLY: max_crown_area, min_diameter, k_crown_area, k_rp, k_fpc, k_sai2lai, itree
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in) :: nc                          !< dimensions
    INTEGER,                  INTENT(in) :: lctlib_growthform           !< lctlib paramter
    REAL(wp),                 INTENT(in) :: lctlib_wood_density         !< lctlib paramter
    REAL(wp),                 INTENT(in) :: lctlib_k_latosa             !< lctlib paramter
    REAL(wp), DIMENSION(nc),  INTENT(in) :: veg_pool_sap_wood_carbon, & !< sapwood mass (mol C / m2)
                                            lai                     , & !< lai (m2/m2)
                                            dens_ind                , & !< tree density (#/m2)
                                            height                  , & !< tree height (m)
                                            diameter                    !< tree diameter (m)
    REAL(wp), DIMENSION(nc)              :: fpc                         !< foliage projective cover (fraction)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)                     :: crown_area_tile      ! crown area times individual density
    REAL(wp), DIMENSION(nc)                     :: lai_ind              ! LAI of the individual
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_foliage_projective_cover_tile'
    ! ----------------------------------------------------------------------------------------------------- !

    SELECT CASE (lctlib_growthform)
    CASE (itree)
      WHERE((dens_ind(:) > eps4) .AND. (height(:) > eps4))
        ! crown area derived from diameter following Sitch et al. 2003; eq 4
        crown_area_tile(:) = MIN(max_crown_area, k_crown_area * MAX(min_diameter, diameter(:)) ** k_rp) * &
          &                  dens_ind(:)
        ! maximum annual, individual LAI, derived from the sapwood area of an individual tree, devided by
        ! its crown area
        lai_ind(:)    = veg_pool_sap_wood_carbon(:) / lctlib_wood_density / height(:) * lctlib_k_latosa / &
          &             crown_area_tile(:)
      ELSEWHERE
        crown_area_tile(:) = k_crown_area * min_diameter ** k_rp
        lai_ind(:)         = 0.0_wp
      ENDWHERE
    CASE DEFAULT ! all other growthforms
      crown_area_tile(:) = 1.0_wp ! to satisfy the fpc equation
      lai_ind(:)         = lai(:) ! actually to be divided by dens_ind and crown_area,which are both 1
    END SELECT

    ! foliage projective cover is then the product of density and individual crown area
    fpc(:) = crown_area_tile(:) * (1._wp - EXP(-k_fpc * lai_ind(:) * (1._wp + k_sai2lai)))

  END FUNCTION calc_foliage_projective_cover_tile

  ! ======================================================================================================= !
  !>calculate vegetation mortality
  !>
  SUBROUTINE calc_veg_mortality( &
    & nc                        , &
    & dtime                     , &
    & elements_index_map        , &
    & is_element_used           , &
    & lctlib_growthform         , &
    & lctlib_k1_mort_greff      , &
    & veg_dynamics_scheme       , &
    & an_boc_tvegdyn_mavg       , &
    & net_growth_tvegdyn_mavg   , &
    & lai_tvegdyn_mavg          , &
    & fpc                       , &
    & veg_pool_mt               , &
    & veg_litterfall_mt         , &
    & mortality_rate )

    USE mo_veg_constants,         ONLY: min_greff, k2_mort_greff, k3_mort_greff, fpc_max, &
      &                                 itree, igrass, background_mort_rate_tree, background_mort_rate_grass
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in)    :: nc                         !< dimensions
    REAL(wp),                 INTENT(in)    :: dtime                      !< timestep length
    INTEGER,                  INTENT(in)    :: elements_index_map(:)      !< map bgcm element ID -> IDX
    LOGICAL,                  INTENT(in)    :: is_element_used(:)         !< is element in 'elements_index_map' used
    INTEGER,                  INTENT(in)    :: lctlib_growthform          !< lctlib paramter
    REAL(wp),                 INTENT(in)    :: lctlib_k1_mort_greff       !< lctlib paramter
    CHARACTER(len=*),         INTENT(in)    :: veg_dynamics_scheme        !< vegetation dynamics: none population
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: an_boc_tvegdyn_mavg        !< long-term net C balance at the bottom of the canopy (micro-mol / m2 / s)
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: net_growth_tvegdyn_mavg    !< long-term net growth of the plant (micro-mol / m2 / s)
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: lai_tvegdyn_mavg           !< long-term LAI of the plant
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: fpc                        !< foliage projective cover
    REAL(wp),                 INTENT(in)    :: veg_pool_mt(:,:,:)         !< bgcm veg_pool
    REAL(wp),                 INTENT(inout) :: veg_litterfall_mt(:,:,:)   !< bgcm flux: veg_litterfall
    REAL(wp), DIMENSION(nc),  INTENT(out)   :: mortality_rate             !< current mortality rate (1/timestep)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: n_veg_bgcm_comp                !< number of vegetation bgcm compartments to loop over
    INTEGER                     :: ic                             !< loop over point of chunk
    INTEGER                     :: ielem                          !< loop over bgcm elements
    INTEGER                     :: ix_comp, ix_elem               !< index of compartment or element in bgcm, used for looping
    REAL(wp), DIMENSION(nc)     :: boc_cbalance                   !< C balance at bottom of canopy
    REAL(wp), DIMENSION(nc)     :: greff                          !< long-term growth efficiency
    REAL(wp), DIMENSION(nc)     :: greff_mortality_rate           !< rate of mortality due to low growth efficiency
    REAL(wp), DIMENSION(nc)     :: self_thinning_mortality_rate   !< rate of mortality due to self-thinning
    REAL(wp), DIMENSION(nc)     :: bg_mortality_rate              !< background rate of mortality
    REAL(wp)                    :: hlp1                           !< potential biomass mortality (to be corrected for litterfall)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_mortality'
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.8 get number of bgcm compartments
    !>
    n_veg_bgcm_comp = SIZE(veg_pool_mt, 1)

    !> 0.9 init local var
    !>
    boc_cbalance(:)                 = 0.0_wp
    greff(:)                        = 0.0_wp
    greff_mortality_rate(:)         = 0.0_wp
    self_thinning_mortality_rate(:) = 0.0_wp
    bg_mortality_rate(:)            = 0.0_wp

    ! differ between veg_dynamics_schemes: population & none
    SELECT CASE (veg_dynamics_scheme)
    CASE ("population")
      !>1.0 mortality related to growth efficiency mortality
      !>
      DO ic = 1,nc
        IF (lai_tvegdyn_mavg(ic) > eps8) THEN
          ! growth efficiency per unit leaf area in mol C / m2 LAI / yr
          !   the application of 'MAX()' avoids negative values of greff_mortality_rate(ic) and mortality_rate(ic)
          greff(ic) = MAX(net_growth_tvegdyn_mavg(ic) / lai_tvegdyn_mavg(ic) * one_day * one_year / 1.e6_wp, 0.0_wp)
          !>  1.1 deduce mortality from asymptoptic mortality rate given growth efficiency,
          !>     below a minimum threshold mortality rises to 100%
          !>
          !     @TODO the k1_mort_greff of PFT 3 (TrBR) is modified (increased compared to other PFT)
          !             to reflect disturbance (C loss) due to dry periods
          !             otherwise TrBR trees would die because of too strong maintanaince respiration during dry season,
          !           this "static" value of k1_mort_greff per PFT may be replaced by future implementations
          !             of explicit disturbance regimes
          IF (greff(ic) > min_greff) THEN
            greff_mortality_rate(ic) = lctlib_k1_mort_greff / (k2_mort_greff * greff(ic) + 1.0_wp)
          ELSE
            greff_mortality_rate(ic) = 1.0_wp / (k3_mort_greff * greff(ic) + 1.0_wp)
          END IF
        END IF
      END DO
      !>  1.2 self-thinning related mortality
      !>
      self_thinning_mortality_rate(:) = MAX(fpc(:) - fpc_max, 0.0_wp)
      ! ! depends on the carbon balance of lowest canopy layer, i.e. net photosynthesis minus average construction cost \n
      ! ! in mol / m2 LAI / year \n
      ! ! transformed to a mortality estimate using a Weibull function
      ! boc_cbalance(:) = an_boc_tvegdyn_mavg(:) * one_day * one_year / 1.e6_wp &
      !   &               - ( 1._wp + fresp_growth ) * 1.0_wp/lctlib_tau_leaf/lctlib_sla
      ! WHERE(boc_cbalance(:) < (-eps4))
      !   mortality_rate(:) = mortality_rate(:) &
      !     &                 + ( 1._wp - exp ( - ( - lambda_mort_light * boc_cbalance(:) ) ** k_mort_light))
      ! END WHERE

      !>  1.3 quasi-stochastic background mortality
      !>
      SELECT CASE (lctlib_growthform)
      CASE (itree)
        bg_mortality_rate(:) = background_mort_rate_tree
      CASE (igrass)
        bg_mortality_rate(:) = background_mort_rate_grass
      ENDSELECT
      !>  1.4 total mortality is the sum of all mortality terms, and limited to one
      !>
      mortality_rate(:) = MIN(greff_mortality_rate(:) + self_thinning_mortality_rate(:) + bg_mortality_rate(:), 1.0_wp) &
        &                 / one_year / one_day * dtime
    CASE ("none")
      SELECT CASE (lctlib_growthform)
      CASE (itree)
        mortality_rate(:) = ((background_mort_rate_tree /  one_year) / one_day) * dtime
      CASE (igrass)
        mortality_rate(:) = ((background_mort_rate_grass / one_year) / one_day) * dtime
      ENDSELECT
    ENDSELECT ! veg_dynamics_scheme

    !>2.0 calculate litter fall from mortality
    !>
    DO ic = 1,nc
      DO ix_comp = 1,n_veg_bgcm_comp
        ! loop over bgcm elements and soil layers
        DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (is_element_used(ielem)) THEN
            ix_elem = elements_index_map(ielem)    ! get element index in bgcm
            ! Only the fraction of veg pool not yet lost to litterfall is available for mortality.
            ! If remaining biomass is smaller than litterfall + mortality put remaining biomass to litter,
            ! else reduce biomass by mortality rate
            hlp1 = veg_pool_mt(ix_comp, ix_elem, ic) * mortality_rate(ic)
            IF (veg_pool_mt(ix_comp, ix_elem, ic) - veg_litterfall_mt(ix_comp, ix_elem, ic) - hlp1 <= 0.0_wp) THEN
              veg_litterfall_mt(ix_comp, ix_elem, ic) = veg_litterfall_mt(ix_comp, ix_elem, ic) &
                &                                       + (veg_pool_mt(ix_comp, ix_elem, ic)    &
                &                                         - veg_litterfall_mt(ix_comp, ix_elem, ic))
            ELSE
              veg_litterfall_mt(ix_comp, ix_elem, ic) = veg_litterfall_mt(ix_comp, ix_elem, ic) &
                &                                       + (veg_pool_mt(ix_comp, ix_elem, ic)    &
                &                                         * mortality_rate(ic))
            END IF
            veg_litterfall_mt(ix_comp, ix_elem, ic) = MAX(0.0_wp, veg_litterfall_mt(ix_comp, ix_elem, ic))
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE calc_veg_mortality

  ! ======================================================================================================= !
  !>
  !> Routine to calculate litter production from herbivory
  !>
  !> This routine calculate the respiration loss and litter production caused by naturally occurring
  !>  herbivory. Following McNaughton et al. 1989, Nature, and Cry and Pace, 1993 Nature, the biomass
  !>  removal can be predicted by a power-law on productivity, leading to an increase of the fraction removed
  !>  from <1\% to >50\% with increase in productivity. This is implemented here by fitting to their data
  !>  and imposing a miniumum productivity of 250gC/m2/yr. An additional N limitation factor is applied to
  !>  account for different palatability of leaves.
  !>
  !> According to metabolic scaling, 90% of the energy consumption is assumed to be respired (instantly, as
  !>  the model does not keep track of herbivory biomass). The remaining biomass, including all nutrients
  !>  enter litter, and thereby increase net mineralisation of N and P via the tighter stoichiometry.
  !>
  !> input: long-term net growth, current lai and leaf stoichiometry
  !> output: litter production and herbivory respiration
  !>
  SUBROUTINE calc_veg_herbivory( &
      & nc, &
      & dtime, &
      & elements_index_map, &
      & is_element_used, &
      & lctlib_GrassFlag, &
      & lctlib_PastureFlag, &
      & lctlib_cn_leaf_min, &
      & lai, &
      & net_growth_tvegdyn_mavg, &
      & mortality_rate, &
      & veg_pool_mt, &
      & herbivory_leaf_resp, &
      & herbivory_leaf_resp_c13, &
      & herbivory_leaf_resp_c14, &
      & herbivory_fruit_resp, &
      & herbivory_fruit_resp_c13, &
      & herbivory_fruit_resp_c14, &
      & fract_herbivory, &
      & veg_litterfall_mt)

    USE mo_veg_constants,               ONLY: igrass, itree, min_lai_herbivory, ftroph_loss_herbivory, &
      &                                       min_net_growth_herbivory_grass, k_herbivory_grass, &
      &                                       min_net_growth_herbivory_pasture, k_herbivory_pasture

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                           !< dimensions
    REAL(wp), INTENT(in)    :: dtime                        !< timestep length
    INTEGER,  INTENT(in)    :: elements_index_map(:)        !< map bgcm element ID -> IDX
    LOGICAL,  INTENT(in)    :: is_element_used(:)           !< is element in 'elements_index_map' used
    LOGICAL,  INTENT(in)    :: lctlib_GrassFlag             !< if tile is natural grassland (TRUE) or not
    LOGICAL,  INTENT(in)    :: lctlib_PastureFlag           !< if tile is pasture (TRUE) or not
    REAL(wp), INTENT(in)    :: lctlib_cn_leaf_min           !< lctlib paramter (minimal leaf C:N ratio, molar units)
    REAL(wp), INTENT(in)    :: lai(:)                       !< current LAI (m2/m2)
    REAL(wp), INTENT(in)    :: net_growth_tvegdyn_mavg(:)   !< long-term average net growth (micro-mol / m2 / s)
    REAL(wp), INTENT(in)    :: mortality_rate(:)            !< mortality rate
    REAL(wp), INTENT(in)    :: veg_pool_mt(:,:,:)           !< vegetation pool
    REAL(wp), INTENT(inout) :: herbivory_leaf_resp(:)       !< respiration due to herbivory on leaves (micro-mol / m2 / s)
    REAL(wp), INTENT(inout) :: herbivory_leaf_resp_c13(:)   !< respiration due to herbivory on leaves (micro-mol / m2 / s)
    REAL(wp), INTENT(inout) :: herbivory_leaf_resp_c14(:)   !< respiration due to herbivory on leaves (micro-mol / m2 / s)
    REAL(wp), INTENT(inout) :: herbivory_fruit_resp(:)      !< respiration due to herbivory on fruits (micro-mol / m2 / s)
    REAL(wp), INTENT(inout) :: herbivory_fruit_resp_c13(:)  !< respiration due to herbivory on fruits (micro-mol / m2 / s)
    REAL(wp), INTENT(inout) :: herbivory_fruit_resp_c14(:)  !< respiration due to herbivory on fruits (micro-mol / m2 / s)
    REAL(wp), INTENT(inout) :: fract_herbivory(:)           !< fraction of leaf and fruit (carbon) lost to herbivory
    REAL(wp), INTENT(inout) :: veg_litterfall_mt(:,:,:)     !< vegetation litterfall
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                 :: n_veg_bgcm_comp              !< number of vegetation bgcm compartments
    INTEGER                 :: ic                           !< loop over nc
    INTEGER                 :: ielem                        !< loop over bgcm elements
    INTEGER                 :: ix_comp, ix_elem             !< index of element and compartment in bgcm, used for looping
    REAL(wp)                :: k_herbivory_act              !< asymptotic herbivory fraction of net growth
    REAL(wp)                :: min_net_growth_herbivory_act !< minimal net growth above which herbivory occurs (mol/m2/yr)
    REAL(wp)                :: eff_fract_herbivory          !< updated herbivory fraction
    REAL(wp)                :: fact_lim_carbon              !< carbon limitation of herbivory (unitless)
    REAL(wp)                :: fact_lim_nitrogen            !< nitrogen limitation of herbivory (unitless)
    REAL(wp)                :: hlp1                         !< helper var
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_herbivory'
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.8 get number of bgcm compartments
    !>
    n_veg_bgcm_comp = SIZE(veg_pool_mt, 1)

    !> 1.0 herbivory rate for grass PFT
    !>
    IF (lctlib_GrassFlag .OR. lctlib_PastureFlag) THEN
      DO ic = 1,nc
        ! set minimum net growth threshold for herbivory according to land cover class
        IF (lctlib_GrassFlag) THEN
          min_net_growth_herbivory_act = min_net_growth_herbivory_grass
          k_herbivory_act              = k_herbivory_grass
        END IF
        IF (lctlib_PastureFlag) THEN
          min_net_growth_herbivory_act = min_net_growth_herbivory_pasture
          k_herbivory_act              = k_herbivory_pasture
        END IF
        ! long-term average productivity (mol/m2/year)
        hlp1 = net_growth_tvegdyn_mavg(ic) * one_day * one_year / 1.e6_wp

        ! limit herbivory to only occur beyond minimal LAI and vegetation net growth
        IF (lai(ic) > min_lai_herbivory .AND. hlp1 > min_net_growth_herbivory_act ) THEN
          ! maximum herbivory fraction increases with long-term average productivity according to
          ! power law (Cyr & Pace, 1993). Following McNaughton et al. 1989, a lower limit of herbivory is applied
          ! note that this leads to a quadratic rise in the flux because that is fherbivory * pool
          fact_lim_carbon = (hlp1 - min_net_growth_herbivory_act) ** (4._wp / 3._wp - 1._wp)
          ! herbivory fraction increases with N content, note this is
          ! rearranged from (leafN/leafC)/(maximum leaf N:C)
          fact_lim_nitrogen = veg_pool_mt(ix_leaf, ixN, ic) / veg_pool_mt(ix_leaf, ixC, ic) &
            &                 * lctlib_cn_leaf_min
          ! herbivory fraction is limited to the fraction that has not died yet
          fract_herbivory(ic) = fract_herbivory(ic) &
            &                   + (1.0_wp - mortality_rate(ic)) &
            &                   * fact_lim_carbon * fact_lim_nitrogen &
            &                   * k_herbivory_act / one_day * dtime
        END IF
      END DO
    END IF

    !> 2.0 calculate litter fall from herbivory
    !>
    DO ic = 1,nc
      DO ix_comp = 1,n_veg_bgcm_comp
        ! seed bed pool remains unaffected from mortality
        ! TODO: needs proper implementation of looping over bgcm compartments
        IF (ix_comp == ix_leaf .OR. ix_comp == ix_fruit) THEN
          ! loop over bgcm elements and soil layers
          DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
            IF (is_element_used(ielem)) THEN
              ix_elem = elements_index_map(ielem)    ! get element index in bgcm
              ! Only the fraction of veg pool not yet lost to litterfall is available for herbivory
              ! If remaining biomass is smaller than herbivory rate, reduce herbivory to available biomass
              ! Otherwise put remaining biomass to litter
              IF (veg_pool_mt(ix_comp, ix_elem, ic) > eps8) THEN
                eff_fract_herbivory = MIN(fract_herbivory(ic), &
                  &   (veg_pool_mt(ix_comp, ix_elem, ic) - veg_litterfall_mt(ix_comp, ix_elem, ic)) &
                  &   / veg_pool_mt(ix_comp, ix_elem, ic))
              ELSE
                eff_fract_herbivory = 0.0_wp
              END IF
              ! for carbon and carbon isotopes
              IF (ix_elem == ixC .OR. ix_elem == ixC13 .OR. ix_elem == ixC14) THEN
                veg_litterfall_mt(ix_comp, ix_elem, ic) = veg_litterfall_mt(ix_comp, ix_elem, ic) &
                  &                                       + (veg_pool_mt(ix_comp, ix_elem, ic) &
                  &                                       * (1._wp - ftroph_loss_herbivory) * eff_fract_herbivory)
                hlp1 = ftroph_loss_herbivory * eff_fract_herbivory * veg_pool_mt(ix_comp, ix_elem, ic) * 1.e6_wp / dtime
                ! leaf herbivory
                IF (ix_comp == ix_leaf) THEN
                  IF (ix_elem == ixC) THEN
                    herbivory_leaf_resp(ic)      = herbivory_leaf_resp(ic)      + hlp1
                  END IF
                  IF (ix_elem == ixC13) THEN
                    herbivory_leaf_resp_c13(ic)  = herbivory_leaf_resp_c13(ic)  + hlp1
                  END IF
                  IF (ix_elem == ixC14) THEN
                    herbivory_leaf_resp_c14(ic)  = herbivory_leaf_resp_c14(ic)  + hlp1
                  END IF
                END IF
                ! fruit herbivory
                IF (ix_comp == ix_fruit) THEN
                  IF (ix_elem == ixC) THEN
                    herbivory_fruit_resp(ic)     = herbivory_fruit_resp(ic)     + hlp1
                  END IF
                  IF (ix_elem == ixC13) THEN
                    herbivory_fruit_resp_c13(ic) = herbivory_fruit_resp_c13(ic) + hlp1
                  END IF
                  IF (ix_elem == ixC14) THEN
                    herbivory_fruit_resp_c14(ic) = herbivory_fruit_resp_c14(ic) + hlp1
                  END IF
                END IF
              ! for all elements but carbon and carbon isotopes
              ELSE
                veg_litterfall_mt(ix_comp, ix_elem, ic) = veg_litterfall_mt(ix_comp, ix_elem, ic) &
                  &                                       + (veg_pool_mt(ix_comp, ix_elem, ic) &
                  &                                       * eff_fract_herbivory)
              END IF
              veg_litterfall_mt(ix_comp, ix_elem, ic) = MAX(0.0_wp, veg_litterfall_mt(ix_comp, ix_elem, ic))
            END IF
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE calc_veg_herbivory

  ! ======================================================================================================= !
  !>calculate vegetation establishment
  !>
  SUBROUTINE calc_veg_establishment( &
    & nc                          , &
    & dtime                       , &
    & elements_index_map          , &
    & is_element_used             , &
    & lctlib_tau_seed_est         , &
    & wtr_plant_avail_rel         , &
    & growing_season              , &
    & t_air_week_mavg             , &
    & fpc                         , &
    & rfr_ratio_boc_tvegdyn_mavg  , &
    & seed_bed_pool_mt            , &
    & veg_establishment_mt         )

    USE mo_veg_constants,               ONLY: fpc_max, lambda_est_temp, k_est_temp, lambda_est_moist, k_est_moist, &
      &                                       rfr_ratio_toc
    USE mo_jsb_physical_constants,      ONLY: Tzero
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in)    :: nc                           !< dimensions
    REAL(wp),                 INTENT(in)    :: dtime                        !< timestep length
    INTEGER,                  INTENT(in)    :: elements_index_map(:)        !< map bgcm element ID -> IDX
    LOGICAL,                  INTENT(in)    :: is_element_used(:)           !< is element in 'elements_index_map' used
    REAL(wp),                 INTENT(in)    :: lctlib_tau_seed_est          !< lctlib paramter
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: wtr_plant_avail_rel          !< plant available water [fraction of maximum]
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: growing_season               !< growing season
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: t_air_week_mavg              !< weekly air temperature [K]
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: fpc                          !< folage projective cover
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: rfr_ratio_boc_tvegdyn_mavg   !< red-farred ratio at the bottom of the canopy
    REAL(wp),                 INTENT(in)    :: seed_bed_pool_mt(:,:)        !< bgcm: seed_bed pool
    REAL(wp),                 INTENT(inout) :: veg_establishment_mt(:,:)    !< bgcm flux: vegetation establishment
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem              !< loop over bgcm elements
    INTEGER                     :: ix_elem            !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc)     :: establishment_rate
    REAL(wp), DIMENSION(nc)     :: tc
    REAL(wp), DIMENSION(nc)     :: rfr_ratio
    REAL(wp), DIMENSION(nc)     :: flim_light
    REAL(wp), DIMENSION(nc)     :: flim_moist
    REAL(wp), DIMENSION(nc)     :: flim_temp
    INTEGER                     :: icanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_establishment'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 establishment rate
    !>

    !>  1.1 light limitation
    !>
    flim_light(:) = MAX(fpc_max - fpc(:), 0.0_wp)

    ! depending on the red to far-red ratio of the light at the bottom of the canopy
    ! flim_light(:) = EXP(-(lctlib_lambda_est_light * &
    !   &             (rfr_ratio_toc-rfr_ratio_boc_tvegdyn_mavg(:)))**lctlib_k_est_light)

    !> 1.2 temperature and soil moisture limitation
    !!
    !! @NOTE  OBS: should be top soil moisture!
    WHERE((t_air_week_mavg(:) - Tzero) < 0.0_wp .OR. growing_season(:) < test_false_true)
      flim_temp(:) = 0.0_wp
    ELSEWHERE
      flim_temp(:) = 1.0_wp - EXP(-(lambda_est_temp * (t_air_week_mavg(:) - Tzero)) ** k_est_temp)
    END WHERE
    flim_moist(:)  = 1.0_wp - EXP(-(lambda_est_moist * wtr_plant_avail_rel(:)) ** k_est_moist)

    !>  1.3 actual establishment rate and associated matter flux from seed bed to reserve pool
    !>
    establishment_rate(:) = flim_light(:) * flim_temp(:) * flim_moist(:) &
      &                     * 1.0_wp / lctlib_tau_seed_est / one_day / one_year * dtime
    WHERE(establishment_rate(:) < eps8)
      establishment_rate(:) = 0.0_wp
    ENDWHERE
    ! calc vegetation establishment
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_establishment_mt(ix_elem, :) = veg_establishment_mt(ix_elem, :) &
          &                                + seed_bed_pool_mt(ix_elem, :) * establishment_rate(:)
      END IF
    END DO
  END SUBROUTINE calc_veg_establishment

  ! ======================================================================================================= !
  !>calculate change in tree density
  !>
  PURE ELEMENTAL FUNCTION calc_delta_dens_ind( &
    & lctlib_seed_size         , &
    & dens_ind                 , &
    & mortality_rate           , &
    & veg_establishment_carbon)  &
    & RESULT (delta_dens_ind)

    REAL(wp), INTENT(in) :: lctlib_seed_size            !< seed size to convert C into individuals, lctlib parameter
    REAL(wp), INTENT(in) :: dens_ind                    !< current individuum density (#/m2)
    REAL(wp), INTENT(in) :: mortality_rate              !< mortality per timestep
    REAL(wp), INTENT(in) :: veg_establishment_carbon    !< vegetation establishment measured in C
    REAL(wp)             :: delta_dens_ind              !< change in individual density
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_density_individuals'
    ! ----------------------------------------------------------------------------------------------------- !

    delta_dens_ind = (veg_establishment_carbon / lctlib_seed_size) - (mortality_rate * dens_ind)
  END FUNCTION calc_delta_dens_ind

#endif
END MODULE mo_q_veg_dynamics
