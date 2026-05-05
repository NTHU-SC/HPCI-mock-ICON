!> QUINCY agriculture constants
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
!>#### declare and define agriculture constants
!>
MODULE mo_q_agr_constants
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  !> INDs of crop types  (since crop types are static, these are also used as IDs)
  ENUM, BIND(C)
    ENUMERATOR ::                 &
      & AGR_SPRING_WHEAT_IDX = 1, &   !< spring wheat
      & AGR_TE_CORN_IDX         , &   !< temperate corn
      & AGR_TE_SOYBEAN_IDX      , &   !< temperate soybean
      & AGR_COTTON_IDX          , &   !< cotton
      & AGR_RICE_IDX            , &   !< rice
      & AGR_SUGARCANE_IDX       , &   !< sugarcane
      & AGR_TR_CORN_IDX         , &   !< tropical corn
      & AGR_TR_SOYBEAN_IDX      , &   !< tropical soybean
      & LAST_CROP_IDX ! needs to be the last -- it is used to determine the number of croptypes
  END ENUM

  INTEGER, PARAMETER :: ix_wheat  = AGR_SPRING_WHEAT_IDX
  INTEGER, PARAMETER :: ix_tecorn = AGR_TE_CORN_IDX
  INTEGER, PARAMETER :: ix_tesoy  = AGR_TE_SOYBEAN_IDX
  INTEGER, PARAMETER :: ix_cot    = AGR_COTTON_IDX
  INTEGER, PARAMETER :: ix_rice   = AGR_RICE_IDX
  INTEGER, PARAMETER :: ix_suca   = AGR_SUGARCANE_IDX
  INTEGER, PARAMETER :: ix_trcorn = AGR_TR_CORN_IDX
  INTEGER, PARAMETER :: ix_trsoy  = AGR_TR_SOYBEAN_IDX
  INTEGER, PARAMETER :: ncroptype = LAST_CROP_IDX - 1

  !> INDs of growth stages types  (since growth_stages are static, these are also used as IDs)
  ENUM, BIND(C)
    ENUMERATOR ::               &
      & AGR_PLANTING_IDX = 0  , &   !< planting
      & AGR_EMERGENCE_IDX     , &   !< emergence
      & AGR_GRAINFILLING_IDX  , &   !< grainfilling
      & AGR_HARVEST_IDX             !< harvesting
  END ENUM
  INTEGER, PARAMETER :: ix_planting     = AGR_PLANTING_IDX
  INTEGER, PARAMETER :: ix_emergence    = AGR_EMERGENCE_IDX
  INTEGER, PARAMETER :: ix_grainfilling = AGR_GRAINFILLING_IDX
  INTEGER, PARAMETER :: ix_harvest      = AGR_HARVEST_IDX

  ! generic parameters for agriculture model
  REAL(wp), SAVE :: &
    & qs_def_n_fertiliser             = def_parameters, &  !< QS-only default rate of fertiliser application for crops per year (mol / m2 / year)
    & iq_background_n_fertiliser_rate = def_parameters, &  !< IQ-only minimum entry from natural fertilisation sources (manure etc.)
                                                           !< in extensive agriculture lacking mineral N input (mol / m2 / year)
    & crop_planting_mass           = def_parameters, &  !< initial mass of crops on planting, taken from the seed bed (mol/m2)
    & fstore_seed_min              = def_parameters, &  !< minimal seed pool size as fraction of crop planting mass
    & fstore_seed_max              = def_parameters, &  !< maximal seed pool size as fraction of crop planting mass
    & max_harvesting_rate          = def_parameters, &  !< maximum harvesting fraction (1/days)
    & fract_crop_to_slash          = def_parameters, &  !< initial mass of crops on planting, taken from the seed bed (mol/m2)
    & frac_fertil_nh4              = def_parameters, &  !< fraction of fertiliser applied that is NH4 (unitless)
    & frac_fertil_po4              = def_parameters, &  !< molar ratio of P fertiliser to N fertiliser application
    & fcs_fertappl_day1            = def_parameters, &  !< fraction of crop growing season when 1st fertiliser dose is applied
    & fcs_fertappl_day2            = def_parameters, &  !< fraction of crop growing season when 2nd fertiliser dose is applied
    & fcs_fertappl_day3            = def_parameters, &  !< fraction of crop growing season when 3rd fertiliser dose is applied
    & fcs_fertappl_day4            = def_parameters, &  !< fraction of crop growing season when 4th fertiliser dose is applied
    & canopy_mode_correction_factor= def_parameters, &  !< calibration factor to allow use of cttab_lai_max with all QUINCY modes
    & eta_fertiliser               = def_parameters     !< fractionation of N15 due to NH3 volatilisation [per mill]

  ! phenology parameters for all crop types
  REAL(wp), SAVE :: &
    & mavg_period_cropseason            = def_parameters, &      !< averaging period for long-term cropseason length (years)
    & gdd_t_thres_dormseason            = def_parameters         !< GDD temperature threshold for C3/C4 crops in dormant season (degC days)

  !< phenology parameters for crop types
  REAL(wp), DIMENSION(ncroptype), SAVE :: &
    & cttab_t_thres_planting, &                                  !< weekly air temperature threshold for planting (K)
    & cttab_t_thres_gdd_gs, &                                    !< reference temperature for GDD calculation during growing season (K)
    & cttab_fgdd_mat, &                                          !< fraction of GDD_X needed to reach maturity (unitless)
    & cttab_gdd_mat_min, &                                       !< minimum GDD requirement for maturity (degC days)
    & cttab_gdd_mat_max, &                                       !< maximum GDD requirement for maturity (degC days)
    & cttab_gdd_max, &                                           !< maximum daily GDD increment (degC days)
    & cttab_nd_gs_max, &                                         !< maximum number of days in to maturity (days)
    & cttab_nd_dorm_max, &                                       !< maximum number of days in dormant period (days)
    & cttab_min_daylength, &                                     !< minimum daylength for emergence (seconds)
    & cttab_fract_gdd_mat_gp2to3, &                              !< fraction of GDD maturity required to start grainfilling (unitless)
    & cttab_crop_season_max                                      !< maximum permitted number of crop seasons per year (# yr-1)

  ! allometry parameters for all crop types
  REAL(wp), SAVE :: &
    & k_heatindex, &                                             !< factor in the sensitivity of allometry to GDD (unitless)
    & leaf2sapwood_min_ratio                                     !< minimum leaf to stem mass ration (unitless)

  !< allometry and allocation parameters for crop types
  REAL(wp), DIMENSION(ncroptype), SAVE :: &
    & cttab_k_hi_leaf, &                                         !< shape parameter for leaf allocation response to heat index (unitless)
    & cttab_k_hi_stem, &                                         !< shape parameter for stem allocation response to heat index (unitless)
    & cttab_leaf2root_mass_init, &                               !< initial ratio of leaf to fine root mass at planting (unitless)
    & cttab_leaf2root_mass_fin, &                                !< final ratio of leaf to fine root mass at maturity (unitless)
    & cttab_sapwood2leaf_mass_init, &                            !< initial ratio of leaf to stem mass at planting (unitless)
    & cttab_sapwood2leaf_mass_fin, &                             !< final ratio of leaf to stem mass at maturity (unitless)
    & cttab_falloc_fruit_max, &                                  !< maximum allocation fraction to fruits (unitless)
    & cttab_active_n_fixation, &                                 !< whether crop is N fixing crop or not (0-1)
    & cttab_lai_max                                              !< maximum leaf area of crop (unitless)

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_agr_constants'

  !$ACC DECLARE CREATE(qs_def_n_fertiliser, iq_background_n_fertiliser_rate, crop_planting_mass)
  !$ACC DECLARE CREATE(fstore_seed_min, fstore_seed_max)
  !$ACC DECLARE CREATE(max_harvesting_rate, fract_crop_to_slash, frac_fertil_nh4, frac_fertil_po4, fcs_fertappl_day1)
  !$ACC DECLARE CREATE(fcs_fertappl_day2, fcs_fertappl_day3, fcs_fertappl_day4, canopy_mode_correction_factor, eta_fertiliser)
  !$ACC DECLARE CREATE(mavg_period_cropseason, gdd_t_thres_dormseason, cttab_t_thres_planting)
  !$ACC DECLARE CREATE(cttab_t_thres_gdd_gs, cttab_gdd_mat_min, cttab_gdd_mat_max, cttab_gdd_max)
  !$ACC DECLARE CREATE(cttab_fgdd_mat, cttab_nd_gs_max, cttab_nd_dorm_max, cttab_min_daylength)
  !$ACC DECLARE CREATE(cttab_fract_gdd_mat_gp2to3, cttab_crop_season_max, k_heatindex, cttab_k_hi_leaf)
  !$ACC DECLARE CREATE(cttab_k_hi_stem, cttab_leaf2root_mass_init, cttab_leaf2root_mass_fin)
  !$ACC DECLARE CREATE(cttab_sapwood2leaf_mass_init, cttab_sapwood2leaf_mass_fin, cttab_falloc_fruit_max)
  !$ACC DECLARE CREATE(cttab_active_n_fixation, cttab_lai_max, leaf2sapwood_min_ratio)

CONTAINS

  ! ======================================================================================================= !
  !> initialize parameters for the process: agriculture
  !>
  !>   routine is called in mo_jsb_base
  !>
  SUBROUTINE init_q_agr_constants
    USE mo_jsb_physical_constants,    ONLY: Tzero, molar_mass_C, molar_mass_N, molar_mass_P
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_q_agr_constants'

    ! general crop parameters
    qs_def_n_fertiliser               = 100._wp / 10._wp / molar_mass_N !< 100 kgN/ha/yr -> 1/10 to give g/m2/yr (only used in QS)
    iq_background_n_fertiliser_rate   = 20._wp / 10._wp / molar_mass_N  !< 20 kgN/ha/yr -> 1/10 to give g/m2/yr (only used in IQ)
    crop_planting_mass                = 4.5_wp / molar_mass_C       !< tuned from 3 gC/m2, AgroIBIS, as AgroIBIS does not need roots for N uptake
    fstore_seed_min                   = 2.0_wp                      !< top up seed pool with fruit harvest if seed pool is smaller than 2 time crop planting mass
    fstore_seed_max                   = 20.0_wp                     !< put seed pool to litter if pool is larger than 20 times crop planting mass
    max_harvesting_rate               = 0.1_wp                      !< harvesting season takes 1/par days
    fract_crop_to_slash               = 0.08_wp                     !< from IPCC bookkeeping model
    fcs_fertappl_day1                 = 0.0_wp                      !< at planting
    fcs_fertappl_day2                 = 0.10_wp                     !< intermediate
    fcs_fertappl_day3                 = 0.35_wp                     !< intermediate
    fcs_fertappl_day4                 = 0.55_wp                     !< during grainfilling
    frac_fertil_nh4                   = 7._wp / 8._wp               !< an global average value
    frac_fertil_po4                   = 0.5_wp * 0.4366197_wp * &   !< FAO statistic global N:P2O5 application
                                        molar_mass_N / molar_mass_P
    canopy_mode_correction_factor     = 0.8_wp                      !< tuned to get reasonable maximum LAI globally in CANOPY mode
    eta_fertiliser                    = 0._wp                       !< not clear how to design this now

    ! phenology and growth stage parameters (all AgroIBIS/CLM5)
    mavg_period_cropseason            = 20.0_wp                     !< AgroIBIS/CLM5
    gdd_t_thres_dormseason            = 8.0_wp + Tzero              !< AgroIBIS/CLM5
    cttab_t_thres_planting = &
      & (/ 7._wp, 10._wp, 13._wp, 21._wp, 21._wp, 21._wp, 21._wp, 21._wp /) + TZero
    cttab_t_thres_gdd_gs = &
      & (/ 5._wp, 8._wp, 10._wp, 5._wp, 5._wp, 8._wp, 10._wp, 10._wp /) + TZero
    cttab_fgdd_mat = &
      & (/ 1._wp, 0.85_wp, 1._wp, 0.85_wp, 1._wp, 1._wp, 1._wp, 1._wp /)
    cttab_gdd_mat_min = &
      & (/ 0._wp, 950._wp, 0._wp , 0._wp, 0._wp, 950._wp, 950._wp, 0._wp /)
    cttab_gdd_mat_max = &
      & (/ 1700._wp, 1850._wp, 1900._wp, 1700._wp, 2100._wp, 1850._wp, 1850._wp, 2100._wp/)
    cttab_gdd_max = &
      & (/ 26._wp, 30._wp, 30._wp, 26._wp, 26._wp, 30._wp, 30._wp, 30._wp/)
    cttab_nd_gs_max = &
      & (/ 150._wp, 165._wp, 150._wp, 160._wp, 150._wp, 300._wp, 160._wp, 150._wp /)
    cttab_nd_dorm_max = &
      & (/ 100._wp, 100._wp, 100._wp, 50._wp, 50._wp, 50._wp, 30._wp, 30._wp /)
    cttab_min_daylength = &
      & (/ 12._wp, 12._wp, 12._wp, 8._wp, 6._wp, 8._wp, 8._wp, 8._wp /) * 3600._wp
    cttab_fract_gdd_mat_gp2to3 = &
      & (/ 0.6_wp, 0.65_wp, 0.5_wp, 0.5_wp, 0.4_wp, 0.65_wp, 0.5_wp, 0.5_wp /)
    cttab_crop_season_max = &
      & (/ 1._wp, 1._wp, 1._wp, 1._wp, 2._wp, 2._wp, 2._wp, 2._wp /)

    ! allocation parameters per crop type
    k_heatindex                             = 0.1_wp      ! AgroIBIS/CLM5 parameter b
    leaf2sapwood_min_ratio                  = 3.0_wp      ! calibrated

    cttab_k_hi_leaf = &
      & (/ 3._wp, 5._wp, 2._wp, 2._wp, 3._wp, 5._wp, 5._wp, 2._wp /)
    cttab_k_hi_stem = &
      & (/ 1._wp, 2._wp, 5._wp, 5._wp, 1._wp, 2._wp, 2._wp, 5._wp /)
    cttab_leaf2root_mass_init = &
      & (/ 2.0_wp, 1.5_wp, 1.5_wp, 2.0_wp, 1.7_wp, 1.5_wp, 1.5_wp, 1.5_wp /)
    cttab_leaf2root_mass_fin = &
      & (/ 1.2_wp, 1.0_wp, 1.0_wp, 1.25_wp, 1.2_wp, 1.0_wp, 1.0_wp, 1.2_wp /)
    cttab_sapwood2leaf_mass_init = &
      & (/ 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp /)
    cttab_sapwood2leaf_mass_fin = &
      & (/ 0.1_wp, 0.1_wp, 0.3_wp, 0.3_wp, 0.1_wp, 0.1_wp, 0.1_wp, 0.3_wp /)
    cttab_falloc_fruit_max = &
      & (/ 0.8_wp, 0.7_wp, 0.6_wp, 0.6_wp, 0.7_wp, 0.8_wp, 0.6_wp, 0.6_wp /)
    cttab_active_n_fixation = &
      & (/ 0._wp, 0._wp, 1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 1._wp /)
    cttab_lai_max = &
      & (/ 7._wp, 5._wp, 6._wp, 6._wp, 7._wp, 5._wp, 5._wp, 6._wp /)

    !$ACC UPDATE DEVICE(qs_def_n_fertiliser, iq_background_n_fertiliser_rate) ASYNC(1)
    !$ACC UPDATE DEVICE(fstore_seed_min, fstore_seed_max) ASYNC(1)
    !$ACC UPDATE DEVICE(crop_planting_mass, max_harvesting_rate, canopy_mode_correction_factor) ASYNC(1)
    !$ACC UPDATE DEVICE(fract_crop_to_slash, frac_fertil_nh4, frac_fertil_po4, fcs_fertappl_day1) ASYNC(1)
    !$ACC UPDATE DEVICE(fcs_fertappl_day2, fcs_fertappl_day3, fcs_fertappl_day4, eta_fertiliser) ASYNC(1)
    !$ACC UPDATE DEVICE(mavg_period_cropseason, gdd_t_thres_dormseason, cttab_t_thres_planting) ASYNC(1)
    !$ACC UPDATE DEVICE(cttab_t_thres_gdd_gs, cttab_gdd_mat_min, cttab_gdd_mat_max, cttab_gdd_max) ASYNC(1)
    !$ACC UPDATE DEVICE(cttab_fgdd_mat, cttab_nd_gs_max, cttab_nd_dorm_max, cttab_min_daylength) ASYNC(1)
    !$ACC UPDATE DEVICE(cttab_fract_gdd_mat_gp2to3, cttab_crop_season_max, k_heatindex, cttab_k_hi_leaf) ASYNC(1)
    !$ACC UPDATE DEVICE(cttab_k_hi_stem, cttab_leaf2root_mass_init, cttab_leaf2root_mass_fin) ASYNC(1)
    !$ACC UPDATE DEVICE(cttab_sapwood2leaf_mass_init, cttab_sapwood2leaf_mass_fin, cttab_falloc_fruit_max) ASYNC(1)
    !$ACC UPDATE DEVICE(cttab_active_n_fixation, cttab_lai_max, leaf2sapwood_min_ratio) ASYNC(1)

  END SUBROUTINE init_q_agr_constants

#endif

END MODULE mo_q_agr_constants
