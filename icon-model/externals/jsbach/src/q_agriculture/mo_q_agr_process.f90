!> QUINCY agriculture process calculations
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
!>#### routines for agriculture processes
!>
MODULE mo_q_agr_process
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_math_constants,  ONLY: eps4, eps8

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_crop_phenology, calc_planting_flux, set_crop_n_fixation_status
  PUBLIC :: update_crop_growth_phase, calc_crop_allocation_factors, calc_fertiliser_application
  PUBLIC :: calc_crop_harvest_fraction, calc_crop_leaf_mass_change_canopy_mode

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_process'

  CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to set GDD required to reach maturity given crop type specific minimum and maximum
  !! as well as a scaling factor for some crops
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION get_gdd_mat(ix_ct, gdd_mavg) RESULT (gdd_mat)

    USE mo_q_agr_constants,         ONLY: cttab_fgdd_mat, cttab_gdd_mat_min, cttab_gdd_mat_max

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in)   :: ix_ct                 !< crop type of current tile/gridcell
    REAL(wp), INTENT(in)  :: gdd_mavg              !< growing season average GDD (degree days)
    REAL(wp)              :: gdd_mat
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_gdd_mat'

    gdd_mat = MIN(MAX(cttab_fgdd_mat(ix_ct) * gdd_mavg, &
                      cttab_gdd_mat_min(ix_ct)),cttab_gdd_mat_max(ix_ct))

  END FUNCTION get_gdd_mat

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate current gdd and its long-term average gdd_mavg for crops,
  !!  and determine inital values of growing demand for N and P at beginning of growing season
  !!  Also initialises crop growth phase to 0 at beginning of growing season
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!  temperature thresholds for gdd vary with crop type
  !!  deviations from CLM5:
  !!    - growing season is additionally constrained by soil moisture (same for all crop types)
  !!    - frost is assessed as frozen soil temperature, no other minimum temperatures are used
  !!    - weekly instead of 10day temperatures are used
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_crop_phenology(  &
    & nc                              , &
    & dtime                           , &
    & is_newyear                      , &
    & beta_soil_flush                 , &
    & gdd_req_max                     , &
    & beta_soil_senescence            , &
    & t_air_senescence                , &
    & cn_leaf                         , &
    & np_leaf                         , &
    & crop_type_index                 , &
    & t_air                           , &
    & t_air_week_mavg                 , &
    & t_air_month_mavg                , &
    & beta_soil_gs_tphen_mavg         , &
    & t_air_tphen_mavg                , &
    & t_soil_srf_tphen_mavg           , &
    & gpp_tlabile_mavg                , &
    & maint_respiration_tlabile_mavg  , &
    & daylength_prev_day              , &
    & growing_season                  , &
    & crop_season_per_year            , &
    & crop_season_per_year_mavg       , &
    & growth_req_n_tlabile_mavg       , &
    & growth_req_p_tlabile_mavg       , &
    & gdd                             , &
    & gdd_mavg                        , &
    & crop_growth_phase               , &
    & nd_dormance)

    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: true, false, test_false_true
    USE mo_q_agr_constants,         ONLY: mavg_period_cropseason, gdd_t_thres_dormseason, cttab_gdd_max,   &
      &                                   cttab_t_thres_gdd_gs, cttab_t_thres_planting, cttab_nd_dorm_max, &
      &                                   cttab_gdd_mat_max, cttab_min_daylength, cttab_crop_season_max,   &
      &                                   ix_planting, ix_harvest
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,      INTENT(in)    :: nc                                  !< dimensions
    REAL(wp),     INTENT(in)    :: dtime                               !< timestep length
    LOGICAL,      INTENT(in)    :: is_newyear                          !< logical indicating start of new year
    REAL(wp),     INTENT(in)    :: beta_soil_flush                     !< lctlib parameter
    REAL(wp),     INTENT(in)    :: gdd_req_max                         !< lctlib parameter
    REAL(wp),     INTENT(in)    :: beta_soil_senescence                !< lctlib parameter
    REAL(wp),     INTENT(in)    :: t_air_senescence                    !< lctlib parameter
    REAL(wp),     INTENT(in)    :: cn_leaf                             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: np_leaf                             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: crop_type_index(:)                  !< crop_type_index
    REAL(wp),     INTENT(in)    :: t_air(:)                            !< air temp
    REAL(wp),     INTENT(in)    :: t_air_week_mavg(:)                  !< moving average air temperature over a week
    REAL(wp),     INTENT(in)    :: t_air_month_mavg(:)                 !< moving average air temperature over a month
    REAL(wp),     INTENT(in)    :: beta_soil_gs_tphen_mavg(:)          !< soil water stress indicator
    REAL(wp),     INTENT(in)    :: t_air_tphen_mavg(:)                 !< moving average air temperature over phenology-relevant period
    REAL(wp),     INTENT(in)    :: t_soil_srf_tphen_mavg(:)            !< moving average soil surface temperature
    REAL(wp),     INTENT(in)    :: gpp_tlabile_mavg(:)                 !< moving average gpp over tlabile
    REAL(wp),     INTENT(in)    :: maint_respiration_tlabile_mavg(:)   !< moving average maintenance respiration over tlabile
    REAL(wp),     INTENT(in)    :: daylength_prev_day(:)               !< daylength of the previous day (seconds)
    REAL(wp),     INTENT(inout) :: growing_season(:)                   !< growing season Yes/No
    REAL(wp),     INTENT(inout) :: crop_season_per_year(:)             !< number of crop seasons this year
    REAL(wp),     INTENT(inout) :: crop_season_per_year_mavg(:)        !< long-term average of crop seasons per year
    REAL(wp),     INTENT(inout) :: growth_req_n_tlabile_mavg(:)        !< moles N required for a unit of C growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: growth_req_p_tlabile_mavg(:)        !< moles P required for a unit of N growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: gdd(:)                              !< growing degree days
    REAL(wp),     INTENT(inout) :: gdd_mavg(:)                         !< long-term average growing degree days
    REAL(wp),     INTENT(inout) :: crop_growth_phase(:)                !< category for growth phase of crops
    REAL(wp),     INTENT(inout) :: nd_dormance(:)                      !< number of days of dormancy
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic                                  !< loop over grid cells
    INTEGER                     :: ix_ct                               !< index of crop type
    REAL(wp)                    :: hlp1                                !< helper variable for tmp values
    REAL(wp)                    :: growing_season_old                  !< growing season state of last time step
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_phenology'
    ! ----------------------------------------------------------------------------------------------------- !

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(ix_ct, growing_season_old, hlp1)
    DO ic = 1,nc
      !>  0.9 memorise old growing season state, and update last year's number of growing season if now is
      !>      a new year
      !>
      ix_ct = INT(crop_type_index(ic))
      growing_season_old = growing_season(ic)
      IF (is_newyear) THEN
        crop_season_per_year_mavg(ic) = crop_season_per_year_mavg(ic) * (1._wp - 1._wp/mavg_period_cropseason) + &
                                       crop_season_per_year(ic) * 1._wp/mavg_period_cropseason
        crop_season_per_year(ic) = 0.0_wp
      END IF

      !>  1.0 check growing season has started or ended
      !>
      ! If not in growing season, check whether growing season has started
      IF (growing_season(ic) < test_false_true) THEN
        ! if temperature and soil moisture are not limiting growth
        IF (gdd(ic) > gdd_req_max .AND. beta_soil_gs_tphen_mavg(ic) > beta_soil_flush &
            & .AND. daylength_prev_day(ic) > cttab_min_daylength(ix_ct) &
            & .AND. t_air_week_mavg(ic) > cttab_t_thres_planting(ix_ct) .AND. t_soil_srf_tphen_mavg(ic) > TZero &
            & .AND. nd_dormance(ic) > cttab_nd_dorm_max(ix_ct) &
            & .AND. crop_season_per_year(ic) < cttab_crop_season_max(ix_ct)) THEN
          growing_season(ic) = true
          crop_season_per_year(ic) = crop_season_per_year(ic) + 1._wp
          growth_req_n_tlabile_mavg(ic) = 1._wp / cn_leaf
          growth_req_p_tlabile_mavg(ic) = 1._wp / np_leaf
        END IF
      ! If in the growing season, check whether growing season has ended
      ELSE
        ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND
        ! harvesting has occurred (ie no end of growing season during crop growth)
        IF ((t_air_tphen_mavg(ic) < t_air_senescence .OR. beta_soil_gs_tphen_mavg(ic) < beta_soil_senescence &
            & .OR. gpp_tlabile_mavg(ic) < maint_respiration_tlabile_mavg(ic) &
            & .OR. gdd(ic) > cttab_gdd_mat_max(ix_ct)) .AND. &
            & .NOT.(crop_growth_phase(ic) >= ix_planting .AND. crop_growth_phase(ic) <= ix_harvest) ) THEN
          growing_season(ic) = false
        END IF
      END IF

      !>  2.0 update long-term average growing season length if now is the end of the season
      !>
      IF (growing_season(ic) <= test_false_true .AND. growing_season_old > test_false_true) THEN
        gdd_mavg(ic) = gdd_mavg(ic) * (1._wp - 1._wp / mavg_period_cropseason) &
          &                    + gdd(ic) * 1._wp / mavg_period_cropseason
      END IF

      !>  2.1 calculate current GDD increment
      !>
      !>    follows AgroIBIS in that there is a difference between the temperature
      !>    threshold for growing (crop-type specific) and dormant season
      !>    in the growing season, a maximum increment is applied
      IF (growing_season(ic) >= test_false_true) THEN
        hlp1 = MIN(MAX(t_air(ic) - cttab_t_thres_gdd_gs(ix_ct),0.0_wp),cttab_gdd_max(ix_ct))
      ELSE
        hlp1 = MAX(t_air(ic) - gdd_t_thres_dormseason,0.0_wp)
      END IF
      gdd(ic) = gdd(ic) + hlp1 * dtime / one_day

      !>  2.2 reset GDD counter at start and end of growing season
      !>       set growth phase to zero at start of growing season
      !>
      IF (growing_season(ic) <= test_false_true .AND. growing_season_old > test_false_true) THEN
        gdd(ic) = 0.0_wp
      END IF
      IF (growing_season(ic) > test_false_true .AND. growing_season_old <= test_false_true) THEN
        gdd(ic) = 0.0_wp
        crop_growth_phase(ic) = 0.0_wp
      END IF
      ! reset gdd when there is snow coverage and soil is frozen
      IF (t_air_tphen_mavg(ic) < Tzero .OR. t_soil_srf_tphen_mavg(ic) <= Tzero) THEN
        gdd(ic) = 0.0_wp
      END IF

      !>  2.3 calculate current dormant days increment
      !>
      IF (growing_season(ic) < test_false_true) THEN
        nd_dormance(ic) = nd_dormance(ic) + dtime / one_day
      ELSE
        nd_dormance(ic) = 0.0_wp
      END IF

    END DO
    !$ACC END PARALLEL LOOP
  END SUBROUTINE calc_crop_phenology

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate flux from seed-bed to establishment if this is a planting time step
  !!  (i.e. crop growth phase is 0)
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_planting_flux( &
      & nc, &
      & crop_growth_phase, &
      & seed_bed_pool_mt, &
      & veg_establishment_mt)

    USE mo_lnd_bgcm_idx
    USE mo_q_agr_constants,         ONLY: crop_planting_mass, ix_planting
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                        !< dimensions
    REAL(wp), INTENT(in)    :: crop_growth_phase(:)      !< growth phase (0 = planting)
    REAL(wp), INTENT(in)    :: seed_bed_pool_mt(:,:)     !< bgcm seed-bed pool
    REAL(wp), INTENT(inout) :: veg_establishment_mt(:,:) !< bgcm flux of establishment
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                !< loop over grid cells
    REAL(wp)       :: planting_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_planting_flux'
    ! ----------------------------------------------------------------------------------------------------- !

    ! if growth phase equal zero and there is a seed pool, add transfer initial planting mass to establishment
    ! flux
    DO ic = 1,nc
      IF (INT(crop_growth_phase(ic)) == ix_planting .AND. seed_bed_pool_mt(ixC, ic) > eps8 ) THEN
        planting_rate = MIN(crop_planting_mass, seed_bed_pool_mt(ixC, ic)) / seed_bed_pool_mt(ixC, ic)
        veg_establishment_mt(:,ic) = veg_establishment_mt(:,ic) + planting_rate * seed_bed_pool_mt(:, ic)
      END IF
    END DO

  END SUBROUTINE calc_planting_flux

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to set N fixation status to true (1) or false (0) depending on croptype
  !!  only soybeans are N fixers
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION set_crop_n_fixation_status(crop_type_index) RESULT (active_n_fixation)

    USE mo_q_agr_constants,         ONLY: cttab_active_n_fixation

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: crop_type_index       !< crop type of current tile/gridcell
    REAL(wp)             :: active_n_fixation
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_crop_n_fixation_status'

    active_n_fixation = cttab_active_n_fixation(INT(crop_type_index))

  END FUNCTION set_crop_n_fixation_status

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate current crop growth phase
  !!  This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!   Different to CLM5: crop growth phase 1 starts immediately after planting. This is because
  !!     inital growth has a larger fraction of roots, hence LAI development is slower than in CLM5
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_crop_growth_phase( &
      & nc, &
      & dtime, &
      & crop_type_index, &
      & gdd, &
      & gdd_mavg, &
      & lai, &
      & nd_crop_season, &
      & nd_crop_season_mavg, &
      & crop_growth_phase)

    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_q_agr_constants
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                       !< dimensions
    REAL(wp), INTENT(in)    :: dtime                    !< timestep length
    REAL(wp), INTENT(in)    :: crop_type_index(:)       !< index of crop-type occupying this tile
    REAL(wp), INTENT(in)    :: gdd(:)                   !< current growing degree days (degC day)
    REAL(wp), INTENT(in)    :: gdd_mavg(:)              !< long-term average gdd (degC day)
    REAL(wp), INTENT(in)    :: lai(:)                   !< current leaf area index
    REAL(wp), INTENT(inout) :: nd_crop_season(:)        !< number of days in this growing season
    REAL(wp), INTENT(inout) :: nd_crop_season_mavg(:)   !< long-term average length of growing season
    REAL(wp), INTENT(inout) :: crop_growth_phase(:)     !< growth phase (0 = planting,
                                                        !<                1 = emergence,
                                                        !<                2 = grain fill,
                                                        !<                3 = harvest
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                     !< loop over grid cells
    INTEGER        :: ix_ct                  !< crop type index
    REAL(wp)       :: hlp1
    REAL(wp)       :: gdd_mat                !< GDD at maturity
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_crop_growth_phase'
    ! ----------------------------------------------------------------------------------------------------- !

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(ix_ct, gdd_mat)
    DO ic = 1,nc

      !> 0.9  cropindex and GDD requirement for maturity
      !>
      ix_ct = INT(crop_type_index(ic))
      gdd_mat = get_gdd_mat(ix_ct,gdd_mavg(ic))

      !> 1.0 advance growth phase to next stage if conditions are met
      SELECT CASE (INT(crop_growth_phase(ic)))
        !> 1.1  when the crop is in planting mode
        !>
        CASE (ix_planting)
          ! this only happens at the first time step of the growing season - move to emergence phase
          ! note this differs from AgroIBIS, which requires 5% of GDD_mat to be reached.
          ! this is a result of the labile pool in QUINCY, which ensures that crop growth starts only in
          ! favourable conditions.
          crop_growth_phase(ic) = crop_growth_phase(ic) + 1.0_wp

        !> 1.2  when the crop is in emergence mode
        !>
        CASE (ix_emergence)
          ! if GDD has accumulated to a given fraction of GDD maturity or
          ! LAI has reached crop specific LAI maximum, move to phase 2 (grain filling)
          IF (gdd(ic) > cttab_fract_gdd_mat_gp2to3(ix_ct) * gdd_mat &
              & .OR. lai(ic) > cttab_lai_max(ix_ct)) THEN
            crop_growth_phase(ic) = crop_growth_phase(ic) + 1.0_wp
          END IF

        !> 1.3  when the crop is in grainfilling mode
        !>
        CASE (ix_grainfilling)
          ! if GDD maturity is reached, move to phase 3 (harvest)
          IF (gdd(ic) > gdd_mat .OR. nd_crop_season(ic) >= cttab_nd_gs_max(ix_ct)) THEN
            crop_growth_phase(ic) = crop_growth_phase(ic) + 1.0_wp
          END IF

        !> 1.4  when the crop is in harvesting mode
        !>
        CASE (ix_harvest)
          ! if LAI reaches zero during harvest, enter phase 4 (currently undefined - can be intercrop later)
          ! and update the long-term average crop growing season length
          IF (lai(ic) < eps8) THEN
            crop_growth_phase(ic) = crop_growth_phase(ic) + 1.0_wp
            nd_crop_season_mavg(ic) = nd_crop_season_mavg(ic) * (1._wp - 1._wp / mavg_period_cropseason) &
              &                        + nd_crop_season(ic) * 1._wp / mavg_period_cropseason
          END IF

        CASE DEFAULT ! nothing to do, wait for reintialisation by crop phenology routine

      END SELECT

      !> 2.0  advance nd_crop_season counter, or set to zero outside crop growing season
      !>
      IF (crop_growth_phase(ic) <= ix_harvest) THEN
        nd_crop_season(ic) = nd_crop_season(ic) + dtime / one_day
      ELSE
        nd_crop_season(ic) = 0.0_wp
      END IF

      !> 2.1 if crop growing season has been longer than allowed, set growth phase to harvest
      !>
      IF (nd_crop_season(ic) > cttab_nd_gs_max(ix_ct) .AND. crop_growth_phase(ic) < ix_harvest) THEN
        crop_growth_phase(ic) = ix_harvest
      END IF

    END DO
    !$ACC END PARALLEL LOOP
  END SUBROUTINE update_crop_growth_phase

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate allometry of crops as well as fraction of carbon allocation to fruit
  !!  production. This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation, but has been
  !!  adjusted to match QUINCY's allocation routine, which is based on allometry, not fixed fractions
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_crop_allocation_factors( &
      & nc, &
      & crop_type_index, &
      & crop_growth_phase, &
      & gdd, &
      & gdd_mavg, &
      & lai, &
      & leaf2sapwood_mass_ratio, &
      & leaf2root_mass_ratio, &
      & falloc_fruit_crop)

    USE mo_q_agr_constants
    USE mo_veg_constants,      ONLY: leaf2root_min_ratio
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                          !< dimensions
    REAL(wp), INTENT(in)    :: crop_type_index(:)          !< index of crop-type occupying this tile
    REAL(wp), INTENT(in)    :: crop_growth_phase(:)        !< growth phase (0 = planting, ..., 3 = harvest)
    REAL(wp), INTENT(in)    :: gdd(:)                      !< current growing degree days (degC day)
    REAL(wp), INTENT(in)    :: gdd_mavg(:)                 !< long-term average gdd (degC day)
    REAL(wp), INTENT(in)    :: lai(:)                      !< current leaf area index
    REAL(wp), INTENT(out)   :: leaf2sapwood_mass_ratio(:)  !< leaf to halm mass ratio (unitless)
    REAL(wp), INTENT(out)   :: leaf2root_mass_ratio(:)     !< leaf to fine root mass ratio (unitless)
    REAL(wp), INTENT(out)   :: falloc_fruit_crop(:)        !< fraction of allocation to fruits
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                     !< loop over grid cells
    INTEGER        :: ix_ct                  !< crop type index
    REAL(wp)       :: hlp1, hlp2
    REAL(wp)       :: gdd_mat                !< GDD at maturity
    REAL(wp)       :: heat_index             !< GDD based heat index
    REAL(wp)       :: fheat_index            !< GDD based heat index
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_allocation_factors'

    ! ----------------------------------------------------------------------------------------------------- !

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(ix_ct, gdd_mat, heat_index, fheat_index, hlp1)
    DO ic = 1,nc

      !> 0.9  crop type index and other helper variables
      !>
      ix_ct = INT(crop_type_index(ic))
      ! GDD needed to reach maturity
      gdd_mat = get_gdd_mat(ix_ct,gdd_mavg(ic))
      ! Heat index (defined as the fraction of GDD for maturity that induces grain filling)
      heat_index = cttab_fract_gdd_mat_gp2to3(ix_ct) * gdd_mat
      ! Normalised heat index
      fheat_index = MIN(gdd(ic) / heat_index, 1._wp)

      !> 1.0  allometry of leaf, fine root and stem
      !>
      !> This is adapted from AgroIBIS but reformulated to transform from a fraction based to an
      !> allometry based allocation (consistent with the assumptions in mo_veg_growth:calc_allocation_fraction.
      !> This means that the parameters have been adjusted to match implied outcome
      !> Allocation is responding to cummulated heat sums during emergence (declining allocation to leaves
      !> with rising maturity). Once the heat requirement for grainfilling is met, allocation is more strongly
      !> towards roots.
      !>
      hlp1 = (exp(-k_heatindex) - exp(-k_heatindex*fheat_index))/(exp(-k_heatindex) - 1._wp)
      leaf2root_mass_ratio(ic) = cttab_leaf2root_mass_fin(ix_ct) &
        &                        + (cttab_leaf2root_mass_init(ix_ct) - cttab_leaf2root_mass_fin(ix_ct)) * hlp1
      leaf2sapwood_mass_ratio(ic) = 1._wp / (cttab_sapwood2leaf_mass_init(ix_ct) &
        &                                    + cttab_sapwood2leaf_mass_fin(ix_ct) * (1._wp - hlp1))

      ! if maximum LAI is reached, reduce leaf allocation to minimum, implying that in the allocation
      ! routine the actual allocation to leaves will be zero
      ! this is also applied during harvesting to avoid prolonged harvesting due to a burst in leaf growth
      ! as the LAI drops below the cttab_lai_max(ix_ct)
      IF (lai(ic) > cttab_lai_max(ix_ct) &
          .OR. INT(crop_growth_phase(ic)) == ix_harvest) leaf2root_mass_ratio(ic) = leaf2root_min_ratio

      !> 2.0  calculate grainfill ratio for crop allocation
      !>
      !> this is adapted from AgroIBIS in that instead of decreasing the leaf and shoot allocation fraction
      !> by (1-hlp1)^d, we increase the fruit allocation fraction as accordingly. This is needed to make
      !> the calculation consistent with the assumptions in calc_allocation_fraction is mo_veg_growth. The consequence
      !> of this formulation is nevertheless that leaf, root and shoot allocation decline in the grainfilling period
      !>
      IF (INT(crop_growth_phase(ic)) == ix_grainfilling .OR. INT(crop_growth_phase(ic)) == ix_harvest) THEN
        hlp1 = MIN(1.0_wp,MAX(0.01_wp,(gdd(ic) - heat_index) / (gdd_mat - heat_index)))
        falloc_fruit_crop(ic) = cttab_falloc_fruit_max(ix_ct) * (1._wp - (1._wp - hlp1)**cttab_k_hi_leaf(ix_ct))
      ELSE
        falloc_fruit_crop(ic) = 0.0_wp
      END IF

    END DO
    !$ACC END PARALLEL LOOP
  END SUBROUTINE calc_crop_allocation_factors

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to deal with fertiliser application
  !!  This follows OCN in that at specific days during the crop growing season
  !!  (expressed as a fraction of the number of days the crop growing season lasts) a given rate
  !!  of fertiliser is applied.
  !!  P fertiliser is assumed a time invariant constant fraction of N fertiliser. This is only a first
  !!  order guess
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_fertiliser_application( &
      & nc, &
      & nd_crop_season, &
      & nd_crop_season_mavg, &
      & crop_season_per_year_mavg, &
      & n_fertiliser, &
      & fertiliser_nh4, &
      & fertiliser_nh4_n15, &
      & fertiliser_no3, &
      & fertiliser_no3_n15, &
      & fertiliser_po4)

    USE mo_jsb_math_constants,     ONLY: one_day
    USE mo_jsb_impl_constants,     ONLY: test_false_true
    USE mo_q_agr_constants,        ONLY: fcs_fertappl_day1, fcs_fertappl_day2, fcs_fertappl_day3, &
                                         fcs_fertappl_day4, frac_fertil_nh4, frac_fertil_po4, eta_fertiliser
    USE mo_isotope_util,           ONLY: calc_mixing_ratio_N15N14
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                         !< dimensions
    REAL(wp), INTENT(in)    :: nd_crop_season(:)          !< number of days this growing season has lasted
    REAL(wp), INTENT(in)    :: nd_crop_season_mavg(:)     !< number of days this growing season has lasted
    REAL(wp), INTENT(in)    :: crop_season_per_year_mavg(:) !< average number of growing season per year (#/yr)
    REAL(wp), INTENT(in)    :: n_fertiliser(:)            !< annual N fertiliser rate (mol/m2/yr)
    REAL(wp), INTENT(inout) :: fertiliser_nh4(:),       & !< NH4 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_nh4_n15(:),   & !< 15NH4 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_no3(:),       & !< NO3 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_no3_n15(:),   & !< 15NO3 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_po4(:)          !< PO4 fertiliser application [micro mol / m2 / s]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                !< loop over grid cells
    REAL(wp)       :: hlp1
    LOGICAL        :: is_fertiliser_day
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fertiliser_application'
    ! ----------------------------------------------------------------------------------------------------- !

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(is_fertiliser_day, hlp1)
    DO ic = 1,nc
      !> 1.0 determine if today is a fertiliser day (given that this is a FLOOR, this is true for 24h)
      !>     the first application day is forced to be the first day after planting to avoid fertilisation
      !>     outside the growing season
      !>
      IF (FLOOR(nd_crop_season(ic)) == MAX(1,FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day1)) &
          .OR. FLOOR(nd_crop_season(ic)) == FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day2)   &
          .OR. FLOOR(nd_crop_season(ic)) == FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day3)   &
          .OR. FLOOR(nd_crop_season(ic)) == FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day4)) THEN
        is_fertiliser_day = .TRUE.
      ELSE
        is_fertiliser_day = .FALSE.
      END IF
      !> 2.0 Mineral fertiliser application for croplands
      !>
      IF (is_fertiliser_day) THEN
        IF (crop_season_per_year_mavg(ic) > 1.5_wp) THEN
          hlp1 = n_fertiliser(ic) / one_day * 1.e6_wp / 4._wp / 2._wp ! conversion and split into 4 doses per crop cycle
        ELSE
          hlp1 = n_fertiliser(ic) / one_day * 1.e6_wp / 4._wp ! conversion and split into 4 doses
        END IF
        fertiliser_nh4(ic) = hlp1 * frac_fertil_nh4
        fertiliser_nh4_n15(ic) = fertiliser_nh4(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(eta_fertiliser))
        fertiliser_no3(ic) = hlp1 * (1._wp - frac_fertil_nh4)
        fertiliser_no3_n15(ic) = fertiliser_no3(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(eta_fertiliser))
        fertiliser_po4(ic) = hlp1 * frac_fertil_po4
      END IF
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_fertiliser_application

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics, called by update_veg_pool_on_harvest
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate flux from harvesting if crop has reached maturity (i.e. growth phase is 3)
  !!  This follows the Q_AGR_ calculation of a litter flux, which can be later rerouted to products pools
  !!  The assumption is as the model is applied at grid scale that harvesting takes 1/max_harvesting_rate
  !!  days to complete
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_crop_harvest_fraction(dtime, crop_type_index, &
                                                crop_growth_phase, lai) RESULT (fract_harvest_rel_to_tile)

    USE mo_jsb_math_constants,     ONLY: one_day
    USE mo_q_agr_constants,        ONLY: cttab_lai_max, max_harvesting_rate, ix_harvest
    USE mo_veg_constants,          ONLY: min_lai

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: dtime                 !< time step length (s)
    REAL(wp), INTENT(in) :: crop_type_index       !< crop type of current tile/gridcell
    REAL(wp), INTENT(in) :: crop_growth_phase     !< crop growth phase
    REAL(wp), INTENT(in) :: lai                   !< curent crop lai
    REAL(wp)             :: fract_harvest_rel_to_tile
    ! ---------------------------
    ! 0.2 Local
    INTEGER              :: ix_ct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_harvest_fraction'

    ix_ct = INT(crop_type_index)
    IF (INT(crop_growth_phase) == ix_harvest) THEN
      IF (lai > min_lai) THEN
        fract_harvest_rel_to_tile = MIN(1.0_wp, MAX(eps4, max_harvesting_rate * dtime / one_day &
          &                                             * cttab_lai_max(ix_ct) / lai))
      ELSE
        fract_harvest_rel_to_tile = 1.0_wp
      END IF
    ELSE
      fract_harvest_rel_to_tile = 0.0_wp
    END IF

  END FUNCTION calc_crop_harvest_fraction

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate leaf growth and litterfall required to simulate LAI change depending on current
  !! phenological state in QCANOPY mode (i.e. if only fast biogeophysical processes are calculated)
  !!
  !!   Input:
  !!     (1) crop_growth_phase, GDD, leaf carbon pool
  !!
  !!   Output: leaf carbon, nitrogen and phosphorus growth and litter fall
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_crop_leaf_mass_change_canopy_mode(   &
    & nc                            , &
    & dtime                         , &
    & lctlib_sla                    , &
    & lctlib_cn_leaf                , &
    & lctlib_np_leaf                , &
    & crop_type_index, &
    & crop_growth_phase, &
    & gdd, &
    & gdd_mavg, &
    & veg_pool_leaf_carbon          , &
    & veg_growth_leaf_carbon        , &
    & veg_growth_leaf_nitrogen      , &
    & veg_growth_leaf_phosphorus    , &
    & veg_litterfall_leaf_carbon    , &
    & veg_litterfall_leaf_nitrogen  , &
    & veg_litterfall_leaf_phosphorus)

    USE mo_jsb_math_constants,            ONLY: one_day,one_year
    USE mo_q_agr_constants,               ONLY: cttab_fract_gdd_mat_gp2to3, cttab_lai_max, max_harvesting_rate, &
      &                                         canopy_mode_correction_factor, &
      &                                         ix_planting, ix_emergence, ix_grainfilling, ix_harvest
    !------------------------------------------------------------------------------------------------------ !
    INTEGER,      INTENT(in)    :: nc                                !< dimensions
    REAL(wp),     INTENT(in)    :: dtime                             !< timestep length
    REAL(wp),     INTENT(in)    :: lctlib_sla                        !< lctlib parameter
    REAL(wp),     INTENT(in)    :: lctlib_cn_leaf                    !< lctlib parameter
    REAL(wp),     INTENT(in)    :: lctlib_np_leaf                    !< lctlib parameter
    REAL(wp),     INTENT(in)    :: crop_type_index(:)                !< index of crop-type occupying this tile
    REAL(wp),     INTENT(in)    :: crop_growth_phase(:)              !< current crop growth phase
    REAL(wp),     INTENT(in)    :: gdd(:)                            !< current GDD (degC days)
    REAL(wp),     INTENT(in)    :: gdd_mavg(:)                       !< long-term average of GDD (degC days)
    REAL(wp),     INTENT(in)    :: veg_pool_leaf_carbon(:)           !< the plant's leaf carbon pool (mol/m2)
    REAL(wp),     INTENT(inout) :: veg_growth_leaf_carbon(:)         !< the plant's current leaf growth rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_growth_leaf_nitrogen(:)       !< the plant's current leaf growth rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_growth_leaf_phosphorus(:)     !< the plant's current leaf growth rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_litterfall_leaf_carbon(:)     !< the plant's current litter fall rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_litterfall_leaf_nitrogen(:)   !< the plant's current litter fall rate (mol/time step)
    REAL(wp),     INTENT(inout) :: veg_litterfall_leaf_phosphorus(:) !< the plant's current litter fall rate (mol/time step)
    !------------------------------------------------------------------------------------------------------ !
    INTEGER                     :: ic                                !< loop over point of the chunk
    INTEGER                     :: ix_ct                             !< crop type index
    REAL(wp)                    :: heat_index,fheat_index,lai_max    !< helper variables
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_leaf_mass_change_canopy_mode'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1, nc
      !>
      !> 0.9  crop type index
      !>
      ix_ct = INT(crop_type_index(ic))
      !> adjust LAI max from CTTAB maximum (with dynamic carbon) to calibrated value
      lai_max = canopy_mode_correction_factor * cttab_lai_max(ix_ct)

      !> 1.0 update leaf growth and/or litter fall based on crop growth phase
      SELECT CASE (INT(crop_growth_phase(ic)))
        !> 1.1  when the crop is growing, add leaves to meet target LAI
        !>
        CASE (ix_planting, ix_emergence, ix_grainfilling )
          ! Heat index (defined as the fraction of GDD for maturity that induces grain filling)
          heat_index = cttab_fract_gdd_mat_gp2to3(ix_ct) * get_gdd_mat(ix_ct,gdd_mavg(ic))
          ! Normalised heat index
          fheat_index = MIN(gdd(ic) / heat_index, 1._wp)
          veg_growth_leaf_carbon(ic) = MAX(fheat_index * lai_max / lctlib_sla &
                                           - veg_pool_leaf_carbon(ic),0.0_wp)
          ! update growth and litterfall of leaf N:P pools
          veg_growth_leaf_nitrogen(ic) = veg_growth_leaf_carbon(ic) / lctlib_cn_leaf
          veg_growth_leaf_phosphorus(ic) = veg_growth_leaf_nitrogen(ic) / lctlib_np_leaf

        !> 1.1  when the crop is harvested, remove leaves to meet target LAI
        !>
        CASE (ix_harvest)
          veg_litterfall_leaf_carbon(ic) = MIN(veg_pool_leaf_carbon(ic), &
              &       lai_max/lctlib_sla * max_harvesting_rate * dtime / one_day)
          ! update growth and litterfall of leaf N:P pools
          veg_litterfall_leaf_nitrogen(ic) = veg_litterfall_leaf_carbon(ic) / lctlib_cn_leaf
          veg_litterfall_leaf_phosphorus(ic) = veg_litterfall_leaf_nitrogen(ic) / lctlib_np_leaf

        CASE DEFAULT ! nothing to do, wait for reintialisation by crop phenology routine

      END SELECT

    END DO

  END SUBROUTINE calc_crop_leaf_mass_change_canopy_mode

#endif
END MODULE mo_q_agr_process
