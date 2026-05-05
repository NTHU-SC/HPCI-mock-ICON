!> QUINCY routines for multiple processes of JSM
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
!>#### multiple routines for the processes using with JSM (jena soil model), e.g., turnover, in- / organic matter fluxes
!>
MODULE mo_q_sb_jsm_processes
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_math_constants,  ONLY: one_day, eps8
  USE mo_lnd_bgcm_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_litter_turnover, &
    &       calc_microbial_turnover, calc_asymb_bnf_fraction, calc_microbial_growth, &
    &       calc_bulk_soil_carbon, &
    &       calc_bulk_density_correction, calc_qmax_bulk_density_correction, &
    &       calc_sorption_desorption_of_org_material, calc_depolymerisation_SESAM, &
    &       calc_sinking_flux, calc_sourcing_flux, &
    &       calc_fast_po4, calc_nitri_denitri_partitioning_rate, &
    &       calc_Psorption_parameter

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_jsm_processes'

CONTAINS

  ! ======================================================================================================= !
  !>calculate turnover of litter pools according to first order decay
  !>
  !>  provides a link between the soluable and woody litter pools introduced in QUINCY, and the
  !>    DOM and polymeric pools of JSM
  !>
  !>  Input: temperature and moisture dependence of litter decomposition, litter pool size, and its
  !>         turnover time at standard conditions
  !>  Output: transfer of litter to its sink pool (soluable -> DOM; woody -> polymeric litter)
  !>
  SUBROUTINE calc_litter_turnover( &
    & nc                    , &
    & nsoil_w               , &
    & dtime                 , &
    & tau_litter_pool       , &
    & fc_litter2sink        , &
    & rtm_depolymerisation  , &
    & rmm_depolymerisation  , &
    & pool_mt_litter        , &
    & het_respiration       , &
    & het_respiration_c13   , &
    & het_respiration_c14   , &
    & loss_mt_litter        , &
    & formation_mt_sink)

    INTEGER,                            INTENT(in)    :: nc                       !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_w                  !< number of soil layers
    REAL(wp),                           INTENT(in)    :: dtime                    !< timestep length
    REAL(wp),                           INTENT(in)    :: tau_litter_pool          !< turnover time of the litter pool [1/s]
    REAL(wp),                           INTENT(in)    :: fc_litter2sink           !< fraction of litter C entering sink pool (rest is respired)
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rtm_depolymerisation     !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rmm_depolymerisation     !< moisture rate modifier [-]
    REAL(wp),                           INTENT(in)    :: pool_mt_litter(:,:,:)    !< bgcm: litter pool to decompose [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: het_respiration          !< heterotrophic respiration [micro-mol C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: het_respiration_c13      !< heterotrophic respiration [micro-mol 13C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: het_respiration_c14      !< heterotrophic respiration [pmol 14C/m3/s]
    REAL(wp),                           INTENT(inout) :: loss_mt_litter(:,:,:)    !< bgcm: reduction in litter pool [mol/m3/timestep]
    REAL(wp),                           INTENT(inout) :: formation_mt_sink(:,:,:) !< bgcm: formation of SOM pool [mol/m3/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_w)                  :: hlp1
    REAL(wp), DIMENSION(nc, nsoil_w)                  :: decay_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_turnover'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calulate decay per timestep
    !>
    decay_rate(:,:) = rtm_depolymerisation * rmm_depolymerisation / tau_litter_pool * dtime

    !>2.0 calculate source and sink terms
    !>
    !>  assuming to change in composition or mass of the material to decompose
    !>
    ! C
    hlp1(:,:)                       = decay_rate(:,:) * pool_mt_litter(ixC, :, :)
    loss_mt_litter(ixC, :, :)       = loss_mt_litter(ixC, :, :) + hlp1(:,:)
    formation_mt_sink(ixC, :, :)    = formation_mt_sink(ixC, :, :) + fc_litter2sink * hlp1(:,:)
    het_respiration(:,:)            = het_respiration(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime
    ! C13
    hlp1(:,:)                       = decay_rate(:,:) * pool_mt_litter(ixC13, :, :)
    loss_mt_litter(ixC13, :, :)     = loss_mt_litter(ixC13, :, :) + hlp1(:,:)
    formation_mt_sink(ixC13, :, :)  = formation_mt_sink(ixC13, :, :) + fc_litter2sink * hlp1(:,:)
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime
    ! C14
    hlp1(:,:)                       = decay_rate(:,:) * pool_mt_litter(ixC14, :, :)
    loss_mt_litter(ixC14, :, :)     = loss_mt_litter(ixC14, :, :) + hlp1(:,:)
    formation_mt_sink(ixC14, :, :)  = formation_mt_sink(ixC14, :, :) + fc_litter2sink * hlp1(:,:)
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime
    ! N
    hlp1(:,:)                       = decay_rate(:,:) * pool_mt_litter(ixN, :, :)
    loss_mt_litter(ixN, :, :)       = loss_mt_litter(ixN, :, :) + hlp1(:,:)
    formation_mt_sink(ixN, :, :)    = formation_mt_sink(ixN, :, :) + hlp1(:,:)
    ! N15
    hlp1(:,:)                       = decay_rate(:,:) * pool_mt_litter(ixN15, :, :)
    loss_mt_litter(ixN15, :, :)     = loss_mt_litter(ixN15, :, :) + hlp1(:,:)
    formation_mt_sink(ixN15, :, :)  = formation_mt_sink(ixN15, :, :) + hlp1(:,:)
    ! P
    hlp1(:,:)                       = decay_rate(:,:) * pool_mt_litter(ixP, :, :)
    loss_mt_litter(ixP, :, :)       = loss_mt_litter(ixP, :, :) + hlp1(:,:)
    formation_mt_sink(ixP, :, :)    = formation_mt_sink(ixP, :, :) + hlp1(:,:)
  END SUBROUTINE calc_litter_turnover

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from sb_jsm_main
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate the flux rates of P processes
  !! weathering and occlusion according to Wang et al. 2010;
  !! desorption (slow exchange) rate based on the ANIMO model (Groenendijk et al. 2005)
  !! biocehmical mineralization is calculated using MM-kinetics assuming Mic C is the limiting factor
  !! Input: rate modifiers and source pools
  !!
  !! Output: biochemical P mineralization rate and loss of source pools, gross P mineralization,
  !!         weathering, occlusion and (slow) desorption, microbial uptake
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sourcing_flux( &
    & nc, &
    & nsoil_w, &
    & dtime, &
    & num_sl_above_bedrock, &
    & flag_sb_prescribe_po4, &
    & microbial_carbon, &
    & microbial_nitrogen, &
    & microbial_po4, &
    & dom_carbon, &
    & dom_nitrogen, &
    & dom_phosphorus, &
    & rtm_mic_uptake, &
    & rmm_mic_uptake, &
    & rtm_hsc, &
    & rmm_hsc, &
    & rtm_sorption, &
    & rmm_sorption, &
    & rtm_desorption, &
    & rmm_desorption, &
    & dom_assoc_po4, &
    & residue_assoc_po4, &
    & residue_po4, &
    & dom_assoc_carbon, &
    & residue_assoc_carbon, &
    & residue_carbon, &
    & solute_po4, &
    & solute_nh4, &
    & solute_no3, &
    & assoc_fast, &
    & assoc_slow, &
    & po4_primary, &
    & bulk_dens_sl, &
    & volume_min_sl, &
    & fine_root_carbon_sl, &
    & vmax_weath_mineral_sl, &
    & enzyme_frac_AP, &
    & growth_p_limit_smoothed_sl, &
    & loss_dom_assoc_po4, &
    & loss_residue_assoc_po4, &
    & loss_residue_po4, &
    & loss_dom_po4, &
    & biochem_min_po4, &
    & weathering_po4, &
    & mic_upt_nh4, &
    & mic_upt_no3, &
    & mic_upt_po4, &
    & gross_min_po4, &
    & occlusion_po4, &
    & slow_exchange_po4)

    USE mo_sb_constants

    INTEGER,                                     INTENT(in)        :: nc                           !< dimensions
    INTEGER,                                     INTENT(in)        :: nsoil_w                      !< number of soil layers
    REAL(wp),                                    INTENT(in)        :: dtime                        !< timestep length
    REAL(wp), DIMENSION(nc),                     INTENT(in)        :: num_sl_above_bedrock         !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    LOGICAL,                                     INTENT(in)        :: flag_sb_prescribe_po4        !< whether PO4 is prescribed or not
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: microbial_carbon             !< microbial carbon density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: microbial_nitrogen           !< nitrogen in the microbial biomass pool (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: microbial_po4                !< PO4 in the microbial biomass pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: dom_carbon                   !< dom carbon density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: dom_nitrogen                 !< dom nitrogen density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: dom_phosphorus               !< dom phosphorus density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_mic_uptake               !< temperature rate modifier for microbial (and plant) uptake
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_mic_uptake               !< moisture rate modifier for microbial (and plant) uptake
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_hsc                      !< temperature rate modifier for half_saturation constants
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_hsc                      !< moisture rate modifier for half_saturation constants
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_sorption                 !< temperature rate modifier for sorption
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_sorption                 !< moisture rate modifier for sorption
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_desorption               !< temperature rate modifier for desorption
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_desorption               !< moisture rate modifier for desorption
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: dom_assoc_po4                !< PO4 in the mineral associated dom pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: residue_assoc_po4            !< PO4 in the mineral associated microbial necromass pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: residue_po4                  !< PO4 in the microbial necromass pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: dom_assoc_carbon             !< dom_assoc carbon density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: residue_assoc_carbon         !< residue_assoc carbon density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: residue_carbon               !< residue carbon density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: solute_po4                   !< soluble po4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: solute_nh4                   !< soluble nh4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: solute_no3                   !< soluble no3 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: assoc_slow                   !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: assoc_fast                   !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: po4_primary                  !< primary PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: bulk_dens_sl                 !< soil bulk density [kg/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: volume_min_sl                !< mineral soil volume fraction [m3/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: fine_root_carbon_sl          !< carbon of fine roots per soil layer, [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: vmax_weath_mineral_sl        !< weathering coefficient of mineral soil [mol P /m3 /s]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: enzyme_frac_AP               !< fraction of Acid (alkaline) phosphatase of the total enzyme
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: growth_p_limit_smoothed_sl   !< scalor of plant P uptake demand, used for phosphatase downscaling
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(inout)     :: loss_dom_assoc_po4           !< loss of PO4 in the mineral associated dom pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(inout)     :: loss_residue_assoc_po4       !< loss of PO4 in the mineral associated microbial necromass pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(inout)     :: loss_residue_po4             !< loss of PO4 in the microbial necromass pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(inout)     :: loss_dom_po4                 !< loss of PO4 in the dom pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: biochem_min_po4              !< phosphorus biochemical mineralization rate, [micro-mol P /m3 /s]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: gross_min_po4                !< DOM mineralization rate [micro-mol P /m3 /s]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: weathering_po4               !< phosphorus weathering rate [micro-mol P / m3 / s]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: mic_upt_nh4                  !< microbial nh4 uptake rate [micro-mol N / m3 / s], for eca_none only
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: mic_upt_no3                  !< microbial no3 uptake rate [micro-mol N / m3 / s], for eca_none only
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: mic_upt_po4                  !< microbial po4 uptake rate [micro-mol P / m3 / s], for eca_none only
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: occlusion_po4                !< phosphorus occlusion rate, assoc_slow --> occluded [micro-mol P / m3 / s]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: slow_exchange_po4            !< phosphorus desorption rate, assoc_fast --> assoc_slow [micro-mol P / m3 / s]

    INTEGER                             :: ic, is       !< loop over dimensions
    REAL(wp)                            :: hlp1, hlp2, hlp3, hlp4, hlp5, hlp6, hlp7
    REAL(wp)                            :: uptake_c
    REAL(wp)                            :: f_uptake
    REAL(wp)                            :: uptake_n
    REAL(wp)                            :: uptake_p
    REAL(wp)                            :: mic_upt_n
    REAL(wp)                            :: som_pc
    REAL(wp)                            :: enzyme_AP
    REAL(wp), DIMENSION(nc, nsoil_w)    :: gross_mineralisation_nh4

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sourcing_flux'

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.9 init out var
    !>
    biochem_min_po4(:,:)   = 0.0_wp
    gross_min_po4(:,:)     = 0.0_wp
    weathering_po4(:,:)    = 0.0_wp
    mic_upt_nh4(:,:)       = 0.0_wp
    mic_upt_no3(:,:)       = 0.0_wp
    mic_upt_po4(:,:)       = 0.0_wp
    occlusion_po4(:,:)     = 0.0_wp
    slow_exchange_po4(:,:) = 0.0_wp

    !> 1.0 calculate the gross P mineralization, based on the potential depolymerization of DOM
    !>
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        hlp1     = dom_carbon(ic,is) / (km_mic_uptake * rtm_hsc(ic,is) * rmm_hsc(ic,is) + dom_carbon(ic,is) )
        uptake_c = vmax_mic_uptake * rtm_mic_uptake(ic,is) * rmm_mic_uptake(ic,is) * hlp1 * microbial_carbon(ic,is)
        IF (dom_carbon(ic,is) > eps8) THEN
          f_uptake = MIN(uptake_c / dom_carbon(ic,is), 1._wp * 1.E6_wp / dtime)
        ELSE
          f_uptake = 0.0_wp
        END IF
        uptake_p = f_uptake * dom_phosphorus(ic,is)
        uptake_n = f_uptake * dom_nitrogen(ic,is)
        uptake_c = f_uptake * dom_carbon(ic,is)

        !> 1.1 calculate the potential microbial N, P uptake
        !>
        mic_upt_po4(ic,is) = MAX(0._wp, uptake_c * microbial_cue_max / microbial_cn / microbial_np - uptake_p * microbial_pue)
        mic_upt_n          = MAX(0._wp, uptake_c * microbial_cue_max / microbial_cn - uptake_n * microbial_nue)
        hlp1               = vmax_mic_uptake_nh4 * microbial_carbon(ic,is) * solute_nh4(ic,is) / (solute_nh4(ic,is) &
          &                  + km_mic_uptake_nh4)
        hlp2               = vmax_mic_uptake_no3 * microbial_carbon(ic,is) * solute_no3(ic,is) / (solute_no3(ic,is) &
          &                  + km_mic_uptake_no3)
        IF ((hlp1 + hlp2) > eps8) THEN
          mic_upt_nh4(ic,is) = mic_upt_n * hlp1 / (hlp1 + hlp2)
          mic_upt_no3(ic,is) = mic_upt_n * hlp2 / (hlp1 + hlp2)
        ELSE
          mic_upt_nh4(ic,is) = 0._wp
          mic_upt_no3(ic,is) = 0._wp
        END IF
        gross_mineralisation_nh4(ic,is) = (1._wp - microbial_nue) * uptake_n
      END DO
    END DO

    !> 2.0 Calculates the biochemical P mineralization and P gross mineralization, and the loss of the contributing pools.
    !> P is assumed to clevage out from dom, dom associated, residue and residue associated
    !>
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        IF (.NOT. flag_sb_prescribe_po4) THEN
           hlp1 = (dom_assoc_po4(ic,is) + residue_assoc_po4(ic,is) + residue_po4(ic,is) + dom_phosphorus(ic,is))
           ! Biomineralization only occurs when there is available organic P
           IF (hlp1 > eps8) THEN
              ! biomineralisation is regulated by solute_po4 and enzyme concentration (from root and microbe)
              enzyme_AP = fine_root_carbon_sl(ic,is) * f_plant_carrier_po4 * growth_p_limit_smoothed_sl(ic,is) &
                &         + microbial_carbon(ic,is) * enzyme_frac_total * enzyme_frac_AP(ic,is)
              hlp4      = enzyme_AP / (rtm_hsc(ic,is) * rmm_hsc(ic,is) * km_mic_biochem_po4 + enzyme_AP)
              ! instantaneous rate of biochemical mineralisation 1/s,
              ! needs to multiply by 1.E6_wp for standard flux rate unit
              biochem_min_po4(ic,is) = rtm_mic_uptake(ic,is) * rmm_mic_uptake(ic,is) * vmax_biochem_po4 * hlp4

              ! unit of loss terms: mol/m3/timestep
              ! loss from mineral associated DOM
              IF (dom_assoc_carbon(ic,is) <eps8) THEN
                hlp5 = 0.0_wp
              ELSE
                som_pc = dom_assoc_po4(ic,is) / dom_assoc_carbon(ic,is)
                IF (som_pc < biomin_threshold_pc_ratio_min_assoc_om) THEN
                  hlp5 = 0.0_wp
                ELSE
                  hlp5 = som_pc / ( km_pc_biochem_po4 + som_pc)
                END IF
              END IF
              loss_dom_assoc_po4(ic,is) = loss_dom_assoc_po4(ic,is) + biochem_min_po4(ic,is) * dom_assoc_po4(ic,is) * hlp5 * dtime

              ! loss from mineral associated residue
              IF (residue_assoc_carbon(ic,is) <eps8) THEN
                hlp6 = 0.0_wp
              ELSE
                som_pc = residue_assoc_po4(ic,is) / residue_assoc_carbon(ic,is)
                IF (som_pc < biomin_threshold_pc_ratio_min_assoc_om) THEN
                  hlp6 = 0.0_wp
                ELSE
                  hlp6 = som_pc / ( km_pc_biochem_po4 + som_pc)
                END IF
              END IF
              loss_residue_assoc_po4(ic,is) = loss_residue_assoc_po4(ic,is) + biochem_min_po4(ic,is) * residue_assoc_po4(ic,is) &
                &                             * hlp6 * dtime

              ! loss from microbial residue
              IF (residue_carbon(ic,is) <eps8) THEN
                hlp7 = 0._wp
              ELSE
                som_pc = residue_po4(ic,is) / residue_carbon(ic,is)
                IF (som_pc < biomin_threshold_pc_ratio_mic_res) THEN
                  hlp7 = 0.0_wp
                ELSE
                  hlp7 = som_pc / ( km_pc_biochem_po4 / 5._wp + som_pc)
                END IF
              END IF
              loss_residue_po4(ic,is) = loss_residue_po4(ic,is) + biochem_min_po4(ic,is) * residue_po4(ic,is) * hlp7 * dtime

              ! loss from dom
              IF (dom_phosphorus(ic,is) <eps8) THEN
                gross_min_po4(ic,is) = 0._wp
              ELSE
                gross_min_po4(ic,is) = MAX(0._wp, MIN( (dom_phosphorus(ic,is)-uptake_p) / dtime, biochem_min_po4(ic,is) &
                  &                    * dom_phosphorus(ic,is))) * 1.E6_wp
              END IF
              loss_dom_po4(ic,is)  = loss_dom_po4(ic,is) + gross_min_po4(ic,is) * dtime / 1.E6_wp

              ! biochemical mineralisation in micro-mol P / m3 / s
              biochem_min_po4(ic,is) = biochem_min_po4(ic,is) * (dom_assoc_po4(ic,is) * hlp5 + residue_po4(ic,is) &
                &                      * hlp7 + residue_assoc_po4(ic,is) * hlp6) * 1.E6_wp
           ELSE ! IF (hlp1 > eps8)
              biochem_min_po4(ic,is) = 0.0_wp
           END IF
        ELSE  ! IF (.NOT. flag_sb_prescribe_po4)
          biochem_min_po4(ic,is) = 0.0_wp
          gross_min_po4(ic,is)   = 0._wp
        END IF
      END DO
    END DO

    !> 3.0 Pi fluxes [micro-mol / m3 / s] between different mineral P pools
    !> weathering, occlusion and slow exchange between sorbed and labile P pools
    !>
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        hlp3 = fine_root_carbon_sl(ic,is) * f_plant_carrier_po4 + microbial_carbon(ic,is) * f_mic_carrier_po4
        hlp4 = hlp3 / (rtm_hsc(ic,is) * rmm_hsc(ic,is) * km_root_enzyme_weath_po4 + hlp3)
        ! minimum microbial biomass to exudate bio-weathering enzyme
        hlp5 = 0.1_wp
        IF (fine_root_carbon_sl(ic,is) < eps8 .AND. microbial_carbon(ic,is) < hlp5) THEN
          hlp4 = 0._wp
        END IF
        weathering_po4(ic,is) = rtm_hsc(ic,is) / rmm_hsc(ic,is) * vmax_weath_mineral_sl(ic,is) * po4_primary(ic,is) &
          &                     * (frac_background_weath_mineral + hlp4) * 1.E6_wp
        IF ((weathering_po4(ic,is) / 1.E6_wp) > (po4_primary(ic,is) / dtime)) THEN
          weathering_po4(ic,is) = po4_primary(ic,is) / dtime * 1.E6_wp
        END IF
        occlusion_po4(ic,is) = assoc_slow(ic,is) * k_occlude_po4 * 1.E6_wp

        hlp1 = rtm_desorption(ic,is) * k_desorpt_po4 * assoc_slow(ic,is)
        hlp2 = rtm_sorption(ic,is) * k_adsorpt_po4 * assoc_fast(ic,is)
        slow_exchange_po4(ic,is) = (hlp2 - hlp1) * 1.E6_wp
      END DO
    END DO

  END SUBROUTINE calc_sourcing_flux

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_jsm for: po4, nh4, no3
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate sinks of the bio-available inorganic po4 and nh4 and no3 pools
  !! ECA approach (Equilibrium Chemistry Approximation) modified based on \n
  !! Zhu et al. (2016) and Zhu et al. (2017), fast exchange of po4 and nh4 based on Van der Zee et al. (1989)
  !!
  !! Input: solute concentration, microbial biomass, fine root biomass, rate modifiers, km values, \n
  !!        sorption capacity and available sorption sites (po4, nh4)
  !!
  !! Output: microbial and plant inorganic uptake rates, fast adsorption rates (po4, nh4), \n
  !!         nitrification/denitrification rates (nh4/no3)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sinking_flux(nc, nsoil_w, dtime, num_sl_above_bedrock, &                            ! in
    &                                        lctlib_vmax_uptake, substrate, sb_adsorp_scheme, &       ! in
    &                                        microbial_carbon, fine_root_carbon, solute, &            ! in
    &                                        rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &      ! in
    &                                        rtm_plant_uptake, rmm_plant_uptake, km1_upt, km2_upt, &  ! in
    &                                        km_nitrification, km_denitrification, &                  ! in
    &                                        rtm_nit_denit_fication, rmm_nit_denit_fication, &        ! in
    &                                        anaerobic_volume_fraction, km_adsorpt_sl, &              ! in
    &                                        avail_soil_surface_sorption_capacity_sl, &               ! in
    &                                        solute_available, &                                      ! in
    &                                        vmax_nitri_denitri, &      ! OPTIONAL (in)
    &                                        solute_alternative, &      ! OPTIONAL (in)
    &                                        microbial_uptake, &        ! inout
    &                                        plant_uptake, &            ! out
    &                                        plant_uptake_pot, &        ! out
    &                                        fast_adsorpt, &            ! OPTIONAL (out)
    &                                        flux_nitri_denitri)        ! OPTIONAL (out)

    USE mo_sb_constants
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                     INTENT(in)        :: nc                        !< dimensions
    INTEGER,                                     INTENT(in)        :: nsoil_w                   !< number of soil layers
    REAL(wp),                                    INTENT(in)        :: dtime                     !< timestep length
    REAL(wp), DIMENSION(nc),                     INTENT(in)        :: num_sl_above_bedrock      !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp),                                    INTENT(in)        :: lctlib_vmax_uptake        !< lctlib uptake paramter, for eca_none
    CHARACTER(len=*),                            INTENT(in)        :: substrate                 !< po4, nh4, no3, determine the specific ECA paramaters
    CHARACTER(len=*),                            INTENT(in)        :: sb_adsorp_scheme          !< eca_full; eca_none; eca_part
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: microbial_carbon          !< bgcm sb_pool: microbial carbon density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: fine_root_carbon          !< fine_root density (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_mic_uptake            !< temperature rate modifier for microbial uptake
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_mic_uptake            !< moisture rate modifier for microbial uptake
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_plant_uptake          !< temperature rate modifier for plant uptake
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_plant_uptake          !< moisture rate modifier for plant uptake
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: km1_upt                   !< low-affinity parameter for plant uptake [m3/mol]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: km2_upt                   !< high-affinity parameter for plant uptake [m3/mol]
    REAL(wp),                                    INTENT(in)        :: km_nitrification          !< half-saturation constant for carbon-limited (de-)nitrification [mol C m-3]
    REAL(wp),                                    INTENT(in)        :: km_denitrification        !< half-saturation constant for nitrate-limited denitrifcation [mol N m-3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rtm_nit_denit_fication    !< actual temperature response for de/nitrification
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: rmm_nit_denit_fication    !< actual moisture response for de/nitrification
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: anaerobic_volume_fraction !< actual anaerobic volume fraction
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: km_adsorpt_sl             !< half-saturation of solute sorption to mineral soil [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: solute                    !< concentration of the solvent to be taken up (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: avail_soil_surface_sorption_capacity_sl  !< available sorption sites [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(in)        :: solute_available                         !< available solvent for all the sinking processes (mol / m3)
    REAL(wp), OPTIONAL,                          INTENT(in)        :: vmax_nitri_denitri                       !< maximum nitrification/denitrification rate, for eca_none
    REAL(wp), OPTIONAL, DIMENSION(nc, nsoil_w),  INTENT(in)        :: solute_alternative                       !< concentration of the alternative solvent to be taken up by plant (mol / m3)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(inout)     :: microbial_uptake                         !< microbial uptake flux (micro-mol / m3 / s)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: plant_uptake                             !< plant uptake flux (micro-mol / m3 / s)
    REAL(wp), DIMENSION(nc, nsoil_w),            INTENT(out)       :: plant_uptake_pot                         !< plant potential uptake rate (micro-mol / m3 / s)
    REAL(wp), OPTIONAL, DIMENSION(nc, nsoil_w),  INTENT(out)       :: fast_adsorpt                             !< fast adsorption flux (micro-mol / m3 / s)
    REAL(wp), OPTIONAL, DIMENSION(nc, nsoil_w),  INTENT(out)       :: flux_nitri_denitri                       !< nitrification or denitrification flux (micro-mol / m3 / s)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                             :: ic, is       !< loop over dimensions
    REAL(wp)                            :: mic_hlp, plant_hlp, mineral_hlp, nitri_hlp
    REAL(wp)                            :: hlp1, hlp2, hlp3
    REAL(wp)                            :: f_mic_carrier, km1_mic_uptake, km2_mic_uptake, &
                                           f_plant_carrier, km_plant_uptake, km1_plant_uptake, &
                                           k_mic_uptake, k_plant_uptake, vmax_mic_uptake_sub, vmax_plant_uptake, &
                                           k_desorpt_f, vmax_nit_denit, km_nit_denit, &
                                           k_nit_denit, f_nit_carrier, fast_desorpt
    REAL(wp), DIMENSION(nc, nsoil_w)    :: reaction_volume
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sinking_flux'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.8 init out var
    !>
    !> unit conversion from mol to micro-mol is introduced in the calculation when necessary
    !>
    plant_uptake(:,:)         = 0.0_wp
    plant_uptake_pot(:,:)     = 0.0_wp
    IF (substrate == "po4" .OR. substrate == "nh4") THEN
      fast_adsorpt(:,:)         = 0.0_wp
    END IF
    IF (substrate == "nh4" .OR. substrate == "no3") THEN
      flux_nitri_denitri(:,:)   = 0.0_wp
    END IF

    !>0.9 init local variables
    !>
    !> assign ECA specific parameter values based on substrate
    !>
    SELECT CASE (substrate)
      CASE ("po4")
        f_mic_carrier         = f_mic_carrier_po4
        km1_mic_uptake        = km_mic_uptake_po4
        km_plant_uptake       = km_plant_uptake_po4
        f_plant_carrier       = f_plant_carrier_po4
        k_mic_uptake          = k_mic_uptake_po4
        k_plant_uptake        = k_plant_uptake_po4
        vmax_mic_uptake_sub   = vmax_mic_uptake_po4
        vmax_plant_uptake     = vmax_plant_uptake_po4
        k_desorpt_f           = k_desorpt_f_po4
        km_nit_denit          = 0._wp
        k_nit_denit           = 0._wp
        f_nit_carrier         = 0._wp
        vmax_nit_denit        = 0._wp
        reaction_volume(:,:)  = 0._wp
      CASE ("nh4")
        f_mic_carrier         = f_mic_carrier_n
        km1_mic_uptake        = km_mic_uptake_nh4
        km2_mic_uptake        = km_mic_uptake_no3
        km_plant_uptake       = km_plant_uptake_nh4
        km1_plant_uptake      = km_plant_uptake_no3
        f_plant_carrier       = f_plant_carrier_n
        k_mic_uptake          = k_mic_uptake_nh4
        k_plant_uptake        = k_plant_uptake_nh4
        vmax_mic_uptake_sub   = vmax_mic_uptake_nh4
        vmax_plant_uptake     = vmax_plant_uptake_nh4
        k_desorpt_f           = k_desorpt_nh4
        vmax_nit_denit        = vmax_nitrification_nh4
        km_nit_denit          = km_nitrification_nh4
        f_nit_carrier         = f_nit_carrier_nh4
        k_nit_denit           = k_nit_carrier_n
        reaction_volume(:,:)  = 1._wp - anaerobic_volume_fraction(:,:)
      CASE ("no3")
        f_mic_carrier         = f_mic_carrier_n
        km1_mic_uptake        = km_mic_uptake_no3
        km2_mic_uptake        = km_mic_uptake_nh4
        km_plant_uptake       = km_plant_uptake_no3
        km1_plant_uptake      = km_plant_uptake_nh4
        f_plant_carrier       = f_plant_carrier_n
        k_mic_uptake          = k_mic_uptake_no3
        k_plant_uptake        = k_plant_uptake_no3
        vmax_mic_uptake_sub   = vmax_mic_uptake_no3
        vmax_plant_uptake     = vmax_plant_uptake_no3
        k_desorpt_f           = 0._wp
        vmax_nit_denit        = vmax_denitrification_no3
        km_nit_denit          = km_denit_no3
        f_nit_carrier         = f_denit_carrier_no3
        k_nit_denit           = k_nit_carrier_n
        reaction_volume(:,:)  = anaerobic_volume_fraction(:,:)
    END SELECT


    !> 1.0 uptake
    !>
    !> co-limited rates of substrate uptake for microbes and plants
    !>
    !> define uptake scheme to use, 1. ECA; 2. non-ECA
    !>
    SELECT CASE (sb_adsorp_scheme)
    !> 1.1 ECA approach
    !>
    CASE("eca_full", "eca_part")
      DO ic = 1,nc
        DO is = 1, NINT(num_sl_above_bedrock(ic))
          mic_hlp     = km1_mic_uptake + solute_available(ic,is) + f_mic_carrier * microbial_carbon(ic,is) + &
            &           km1_mic_uptake / km_plant_uptake * f_plant_carrier * fine_root_carbon(ic,is)
          plant_hlp   = km_plant_uptake + solute_available(ic,is) + f_plant_carrier * fine_root_carbon(ic,is) + &
            &           km_plant_uptake / km1_mic_uptake * f_mic_carrier * microbial_carbon(ic,is)
          mineral_hlp = km_adsorpt_sl(ic,is) + solute_available(ic,is) + avail_soil_surface_sorption_capacity_sl(ic,is) + &
            &           km_adsorpt_sl(ic,is) / km_plant_uptake * f_plant_carrier * fine_root_carbon(ic,is) + &
            &           km_adsorpt_sl(ic,is) / km1_mic_uptake * f_mic_carrier * microbial_carbon(ic,is)
          nitri_hlp   = km_nit_denit + solute_available(ic,is) + f_nit_carrier * microbial_carbon(ic,is) + &
            &            km_nit_denit / km1_mic_uptake * f_mic_carrier * microbial_carbon(ic,is) + &
            &            km_nit_denit / km_plant_uptake * f_plant_carrier * fine_root_carbon(ic,is)
          ! nh4 & no3
          ! including the effects of alternative solute and nitrification/denitrification process
          IF (substrate == "nh4" .OR. substrate == "no3") THEN
            plant_hlp = plant_hlp + &
              &         km_plant_uptake / km1_plant_uptake * solute_alternative(ic,is) + &
              &         km_plant_uptake / km_nit_denit * f_nit_carrier * microbial_carbon(ic,is)
            mic_hlp   = mic_hlp + &
              &         km1_mic_uptake / km2_mic_uptake * solute_alternative(ic,is) + &
              &         km1_mic_uptake / km_nit_denit * f_nit_carrier * microbial_carbon(ic,is)
            ! nh4
            ! including the effect of adsorption
            IF (substrate=="nh4") THEN
              mic_hlp      = mic_hlp + &
                &            km1_mic_uptake / km_adsorpt_sl(ic,is) * avail_soil_surface_sorption_capacity_sl(ic,is)
              plant_hlp    = plant_hlp + &
                &            km_plant_uptake / km_adsorpt_sl(ic,is) * avail_soil_surface_sorption_capacity_sl(ic,is)
              mineral_hlp  = mineral_hlp + &
                &            km_adsorpt_sl(ic,is) / km_nit_denit * f_nit_carrier * microbial_carbon(ic,is)
              nitri_hlp    = nitri_hlp + &
                &            km_nit_denit / km_adsorpt_sl(ic,is) * avail_soil_surface_sorption_capacity_sl(ic,is)
              hlp1         = k_desorpt_f / km_adsorpt_sl(ic,is) * solute(ic,is)
              fast_adsorpt(ic,is) = MIN(1._wp / dtime, hlp1 * solute_available(ic,is) / mineral_hlp) &
                &            * avail_soil_surface_sorption_capacity_sl(ic,is) * 1.E6_wp
            END IF
          ! po4
          ! including the effect of adsorption
          ELSE
            mic_hlp      = mic_hlp + &
              &            km1_mic_uptake / km_adsorpt_sl(ic,is) * avail_soil_surface_sorption_capacity_sl(ic,is)
            plant_hlp    = plant_hlp + &
              &            km_plant_uptake / km_adsorpt_sl(ic,is) * avail_soil_surface_sorption_capacity_sl(ic,is)
            hlp1         = k_desorpt_f / km_adsorpt_sl(ic,is) * solute(ic,is)
            fast_adsorpt(ic,is) = MIN(1._wp / dtime, hlp1 * solute_available(ic,is) / mineral_hlp) &
              &            * avail_soil_surface_sorption_capacity_sl(ic,is) * 1.E6_wp
          END IF
          ! eca_full, rates calculation based on enzyme abundance
          ! Calculate plant uptake, microbial uptake and nitrification/denitrification rates
          microbial_uptake(ic,is)   = MIN(rtm_mic_uptake(ic,is) * rmm_mic_uptake(ic,is) * k_mic_uptake &
            &                         * solute_available(ic,is) * f_mic_carrier &
            &                         * microbial_carbon(ic,is) / mic_hlp, microbial_uptake(ic,is))
          plant_uptake(ic,is)       = rtm_plant_uptake(ic,is) * rmm_plant_uptake(ic,is) * k_plant_uptake &
            &                         * solute_available(ic,is) * f_plant_carrier &
            &                         * fine_root_carbon(ic,is) / plant_hlp
          IF (substrate == "nh4" .OR. substrate == "no3") THEN
            flux_nitri_denitri(ic,is) = rtm_nit_denit_fication(ic,is) * rmm_nit_denit_fication(ic,is) &
              &                         * reaction_volume(ic,is) * k_nit_denit &
              &                         * solute_available(ic,is) * f_nit_carrier * microbial_carbon(ic,is) / nitri_hlp
          END IF
          ! eca_part, rates calculation based on vmax (prescribed)
          !! Calculate plant uptake, microbial uptake and nitrification/denitrification rates
          IF (sb_adsorp_scheme == "eca_part") THEN
            microbial_uptake(ic,is)   = rtm_mic_uptake(ic,is) * rmm_mic_uptake(ic,is) * vmax_mic_uptake_sub * &
              &                         solute_available(ic,is) * microbial_carbon(ic,is) / mic_hlp
            plant_uptake(ic,is)       = rtm_plant_uptake(ic,is) * rmm_plant_uptake(ic,is) * vmax_plant_uptake * &
              &                         solute_available(ic,is) * fine_root_carbon(ic,is) / plant_hlp
            IF (substrate == "nh4" .OR. substrate == "no3") THEN
              flux_nitri_denitri(ic,is) = rtm_nit_denit_fication(ic,is) * rmm_nit_denit_fication(ic,is) &
                &                         * reaction_volume(ic,is) * vmax_nit_denit &
                &                         * solute_available(ic,is) * microbial_carbon(ic,is) / nitri_hlp
            END IF
          END IF
        END DO
      END DO

    !> 1.2 non-ECA approach
    !>
    !> same as ssm (simple_1d soil model)
    !>
    CASE("eca_none")
      DO ic = 1,nc
        DO is = 1, NINT(num_sl_above_bedrock(ic))
          plant_uptake(ic,is)   = rtm_plant_uptake(ic,is) * rmm_plant_uptake(ic,is) * fine_root_carbon(ic,is) &
            &                     * lctlib_vmax_uptake * solute(ic,is) &
            &                     * (km1_upt(ic,is) + 1._wp / (km2_upt(ic,is) + solute(ic,is)))
          IF (substrate == "po4" .OR. substrate == "nh4") THEN
            fast_adsorpt(ic,is) = 0._wp
          END IF
          IF (substrate == "nh4" .OR. substrate == "no3") THEN
            flux_nitri_denitri(ic,is) = rtm_nit_denit_fication(ic,is) * rmm_nit_denit_fication(ic,is) &
              &                         * vmax_nitri_denitri * microbial_carbon(ic,is) &
              &                         / (km_nitrification + microbial_carbon(ic,is))  * reaction_volume(ic,is) * solute(ic,is)
          END IF
        END DO
      END DO
    END SELECT

    !> 2.0 check the availability of solute considering all sinking fluxes
    !>
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        fast_desorpt = (solute_available(ic,is) - solute(ic,is)) * k_desorpt_f * 1.E6_wp
        SELECT CASE (substrate)
          CASE("po4")
            hlp1 = microbial_uptake(ic,is) + plant_uptake(ic,is)
            hlp2 = MAX(0._wp, ((solute(ic,is)) / dtime * 1.E6_wp - fast_adsorpt(ic,is) + fast_desorpt))
          CASE("nh4")
            hlp1 = microbial_uptake(ic,is) + plant_uptake(ic,is) + flux_nitri_denitri(ic,is)
            hlp2 = MAX(0._wp, ((solute(ic,is)) / dtime * 1.E6_wp - fast_adsorpt(ic,is) + fast_desorpt))
          CASE("no3")
            hlp1 = microbial_uptake(ic,is) + plant_uptake(ic,is) + flux_nitri_denitri(ic,is)
            hlp2 = MAX(0._wp, (solute(ic,is)) / dtime * 1.E6_wp)
        END SELECT
        IF (hlp1 > hlp2 .AND. hlp1 > 0._wp) THEN
          microbial_uptake(ic,is)   = microbial_uptake(ic,is)   * hlp2 / hlp1
          plant_uptake(ic,is)       = plant_uptake(ic,is)       * hlp2 / hlp1
          SELECT CASE (substrate)
            CASE("po4")
              fast_adsorpt(ic,is)       = fast_adsorpt(ic,is)
            CASE("nh4")
              fast_adsorpt(ic,is)       = fast_adsorpt(ic,is)
              flux_nitri_denitri(ic,is) = flux_nitri_denitri(ic,is) * hlp2 / hlp1
            CASE("no3")
              flux_nitri_denitri(ic,is) = flux_nitri_denitri(ic,is) * hlp2 / hlp1
          END SELECT
        END IF
        plant_uptake_pot(ic,is) = plant_uptake(ic,is)
      END DO
    END DO
  END SUBROUTINE calc_sinking_flux

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_jsm
  !
  !-----------------------------------------------------------------------------------------------------
  !> For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to
  !! net mineralisation \n
  !! For the sb_nloss_scheme = dynamic: calculates the partitioning of nitrification and denitrification
  !! to NO3, NOy, N2O and N2, and diagnoses the amount of leaching below the rooting zone
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_nitri_denitri_partitioning_rate( &
                                    nc, &                         ! in
                                    nsoil_w, &
                                    soil_depth_sl, &
                                    sb_nloss_scheme, &
                                    rtm_denitrification_act, &
                                    nitrification_no3, &          ! inout
                                    denitrification_n2, &
                                    nitrification_noy, &          ! out
                                    nitrification_n2o, &
                                    denitrification_noy, &
                                    denitrification_n2o, &
                                    net_mineralisation_nh4, &
                                    net_mineralisation_nh4_n15, &
                                    nitrification_no3_n15, &
                                    volatilisation_nh4, &
                                    volatilisation_nh4_n15, &
                                    emission_n2, &
                                    emission_n2_n15)

    USE mo_sb_constants,            ONLY: f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o,floss_nmin,fnit_nmin

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                           !< dimensions
                                                                       nsoil_w
    CHARACTER(len=*),                                 INTENT(in)    :: sb_nloss_scheme                 !< sb n loss scheme
    REAL(wp), DIMENSION(nc, nsoil_w),                 INTENT(in)    :: soil_depth_sl, &                !< soil depth per layer [m]
                                                                       net_mineralisation_nh4, &       !< net mineralisation rate of NH4 [micro-mol/m2/s]
                                                                       net_mineralisation_nh4_n15, &   !< net mineralisation rate of 15NH4 [micro-mol/m2/s]
                                                                       rtm_denitrification_act         !< temperature modifier for denitrification
    REAL(wp), DIMENSION(nc, nsoil_w),                 INTENT(inout) :: nitrification_no3, &            !< nitrification to NO3 [micro-mol/m2/s]
                                                                       denitrification_n2, &           !< denitrification to N2 [micro-mol/m2/s]
                                                                       nitrification_noy, &            !< nitrification to NOy [micro-mol/m2/s]
                                                                       nitrification_n2o, &            !< nitrification to N2O [micro-mol/m2/s]
                                                                       denitrification_noy, &          !< denitrification to NOy [micro-mol/m2/s]
                                                                       denitrification_n2o, &          !< denitrification to N2O [micro-mol/m2/s]
                                                                       nitrification_no3_n15, &        !< nitrification to 15NO3 [micro-mol/m2/s]
                                                                       volatilisation_nh4, &           !< volatilisation of NH4 [micro-mol/m2/s]
                                                                       volatilisation_nh4_n15          !< volatilisation of 15NH4 [micro-mol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: emission_n2, &                  !< soil efflux of N2 [micro-mol/m2/s]
                                                                       emission_n2_n15                 !< soil efflux of 15N2 [micro-mol/m2/s]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nitri_denitri_partitioning_rate'


    !> 1.0 case of dynamic N loss calculations, following Zaehle et al. 2011
    !!
    IF (sb_nloss_scheme=="dynamic") THEN
      ! nitrification rate to produce NO3, NOy and N2O from NH4 [mol m-3 timestep-1]
      nitrification_noy(:,:) = f_nit_noy * nitrification_no3(:,:)
      nitrification_n2o(:,:) = f_nit_n2o * nitrification_no3(:,:)
      nitrification_no3(:,:) = nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:)

      ! denitrification rate to produce NOy, N2O and N2 from NO3 [mol m-3 timestep-1]
      denitrification_noy(:,:) = f_denit_noy * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
      denitrification_n2o(:,:) = f_denit_n2o * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
      denitrification_n2(:,:)  = denitrification_n2(:,:) - denitrification_noy(:,:) - denitrification_n2o(:,:)

    !> 2.0 case of fixed N losses: N loss proportional to net mineralisation
    !!
    ELSE
      WHERE(net_mineralisation_nh4(:,:) > 0.0_wp)
        ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
        volatilisation_nh4(:,:)      = net_mineralisation_nh4(:,:)     * floss_nmin
        volatilisation_nh4_n15(:,:)  = net_mineralisation_nh4_n15(:,:) * floss_nmin
        ! NH4 transformed to NO3 [mol m-3 timestep-1]
        nitrification_no3(:,:)     = net_mineralisation_nh4(:,:)     * (1._wp - floss_nmin) * fnit_nmin
        nitrification_no3_n15(:,:) = net_mineralisation_nh4_n15(:,:) * (1._wp - floss_nmin) * fnit_nmin
      ENDWHERE
      ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
      emission_n2(:)     = SUM(volatilisation_nh4(:,:)     * soil_depth_sl(:,:), DIM=2)
      emission_n2_n15(:) = SUM(volatilisation_nh4_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    ENDIF

  END SUBROUTINE calc_nitri_denitri_partitioning_rate

  ! ======================================================================================================= !
  !>calculate turnover of the microbial pool
  !>  corresponds to the term 'pi * Cb' in Ahrens et al. 2015, eq. 7
  !>
  !> Input: microbial biomass
  !>
  !> Output: transfer of microbial biomass to its sink pools (soluable -> DOM; polymeric -> residue)
  !>
  SUBROUTINE calc_microbial_turnover( &
    & nc, &
    & nsoil_w, &
    & dtime, &
    & num_sl_above_bedrock, &
    & fact_n_status_mic_c_growth, &
    & fact_p_status_mic_c_growth, &
    & sb_pool_mt_microbial, &
    & sb_loss_mt_microbial, &
    & sb_formation_mt)

    USE mo_sb_constants,          ONLY: k_mort_microbial, min_mic_biomass, f_soluable_residue, fc_p_res2dom, fc_n_res2dom, &
      &                                 k_n_recyc_corr, k_p_recyc_corr
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                          !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_w                     !< number of soil layers
    REAL(wp),                          INTENT(in)    :: dtime                       !< timestep length
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock        !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: fact_n_status_mic_c_growth  !< scalor of N availability on microbial C growth
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: fact_p_status_mic_c_growth  !< scalor of P availability on microbial C growth
    REAL(wp),                          INTENT(in)    :: sb_pool_mt_microbial(:,:,:) !< bgcm sb_pool: microbial biomass [mol/m3]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt_microbial(:,:,:) !< bgcm sb_loss flux: microbial death [mol/m3/timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)    !< bgcm sb_formation: DOM formation from microbial death [mol/m3/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: hlp1
    REAL(wp)                                         :: hlp2
    REAL(wp)                                         :: hlp3
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: p_res2dom
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: n_res2dom
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: mortality_rate
    INTEGER                                          :: ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_turnover'
    ! ----------------------------------------------------------------------------------------------------- !


    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        !>0.9 init local variables
        !>
        !>  Fractional of N and P in microbial necromass is recycled to DOM,
        !>    regulated by the microbial nutritional scalors following exponential decline curve
        !>
        n_res2dom(ic,is) = fc_n_res2dom * MIN(1._wp, EXP(fact_n_status_mic_c_growth(ic,is) * (-1._wp)) ** k_n_recyc_corr)
        p_res2dom(ic,is) = fc_p_res2dom * MIN(1._wp, EXP(fact_p_status_mic_c_growth(ic,is) * (-1._wp)) ** k_p_recyc_corr)

        !>1.0 calculate microbial mortality rate, limited to minimum microbial biomass density (-> dormant microbes)
        !>
        mortality_rate(ic,is) = k_mort_microbial * dtime
        hlp1(ic,is)           = mortality_rate(ic,is) * sb_pool_mt_microbial(ixC, ic, is)
        IF ((sb_pool_mt_microbial(ixC, ic, is) - hlp1(ic,is)) < min_mic_biomass) THEN
          IF (sb_pool_mt_microbial(ixC, ic, is) > eps8) THEN
            mortality_rate(ic,is) = MAX(0._wp, 1._wp - (min_mic_biomass / sb_pool_mt_microbial(ixC, ic, is)))
          ELSE
            mortality_rate(ic,is) = 0._wp
          END IF
        END IF

        !>2.0 calculate matter fluxes associated with microbial death
        !>  a fraction turns to DOM, the remainder becomes residue
        !>  the assumption is made that the composition of this material is not changed
        !>
        !>  @TODO define C:N:P stoichiometries of DOM and residues from microbial death
        !>
        ! C
        hlp1(ic,is)                              = mortality_rate(ic,is) * sb_pool_mt_microbial(ixC, ic, is)
        sb_loss_mt_microbial(ixC, ic, is)        = sb_loss_mt_microbial(ixC, ic, is) + hlp1(ic,is)
        sb_formation_mt(ix_dom, ixC, ic, is)     = sb_formation_mt(ix_dom, ixC, ic, is) + f_soluable_residue * hlp1(ic,is)
        sb_formation_mt(ix_residue, ixC, ic, is) = sb_formation_mt(ix_residue, ixC, ic, is) + (1._wp - f_soluable_residue) &
                                                   * hlp1(ic,is)
        ! N
        hlp1(ic,is)                              = mortality_rate(ic,is) * sb_pool_mt_microbial(ixN, ic, is)
        sb_loss_mt_microbial(ixN, ic, is)        = sb_loss_mt_microbial(ixN, ic, is) + hlp1(ic,is)
        sb_formation_mt(ix_dom, ixN, ic, is)     = sb_formation_mt(ix_dom, ixN, ic, is) + f_soluable_residue * hlp1(ic,is) &
          &                                      + (1._wp - f_soluable_residue) * hlp1(ic,is) * n_res2dom(ic,is)
        sb_formation_mt(ix_residue, ixN, ic, is) = sb_formation_mt(ix_residue, ixN, ic, is) + (1._wp - f_soluable_residue) &
          &                                      * (1._wp - n_res2dom(ic,is)) * hlp1(ic,is)
        ! P
        hlp1(ic,is)                              = mortality_rate(ic,is) * sb_pool_mt_microbial(ixP, ic, is)
        sb_loss_mt_microbial(ixP, ic, is)        = sb_loss_mt_microbial(ixP, ic, is) + hlp1(ic,is)
        sb_formation_mt(ix_dom, ixP, ic, is)     = sb_formation_mt(ix_dom, ixP, ic, is) + f_soluable_residue * hlp1(ic,is) &
          &                                      + (1._wp - f_soluable_residue) * hlp1(ic,is) * p_res2dom(ic,is)
        sb_formation_mt(ix_residue, ixP, ic, is) = sb_formation_mt(ix_residue, ixP, ic, is) + (1._wp - f_soluable_residue) &
          &                                      * (1._wp - p_res2dom(ic,is)) * hlp1(ic,is)
        ! C13
        hlp1(ic,is)                                = mortality_rate(ic,is) * sb_pool_mt_microbial(ixC13, ic, is)
        sb_loss_mt_microbial(ixC13, ic, is)        = sb_loss_mt_microbial(ixC13, ic, is) + hlp1(ic,is)
        sb_formation_mt(ix_dom, ixC13, ic, is)     = sb_formation_mt(ix_dom, ixC13, ic, is) + f_soluable_residue * hlp1(ic,is)
        sb_formation_mt(ix_residue, ixC13, ic, is) = sb_formation_mt(ix_residue, ixC13, ic, is) + (1._wp - f_soluable_residue) &
                                                     * hlp1(ic,is)
        ! C14
        hlp1(ic,is)                                = mortality_rate(ic,is) * sb_pool_mt_microbial(ixC14, ic, is)
        sb_loss_mt_microbial(ixC14, ic, is)        = sb_loss_mt_microbial(ixC14, ic, is) + hlp1(ic,is)
        sb_formation_mt(ix_dom, ixC14, ic, is)     = sb_formation_mt(ix_dom, ixC14, ic, is) + f_soluable_residue * hlp1(ic,is)
        sb_formation_mt(ix_residue, ixC14, ic, is) = sb_formation_mt(ix_residue, ixC14, ic, is) + (1._wp - f_soluable_residue) &
                                                     * hlp1(ic,is)
        ! N15
        hlp1(ic,is)                                = mortality_rate(ic,is) * sb_pool_mt_microbial(ixN15, ic, is)
        sb_loss_mt_microbial(ixN15, ic, is)        = sb_loss_mt_microbial(ixN15, ic, is) + hlp1(ic,is)
        sb_formation_mt(ix_dom, ixN15, ic, is)     = sb_formation_mt(ix_dom, ixN15, ic, is) + f_soluable_residue * hlp1(ic,is) &
          &                                        + (1._wp - f_soluable_residue) * hlp1(ic,is) * n_res2dom(ic,is)
        sb_formation_mt(ix_residue, ixN15, ic, is) = sb_formation_mt(ix_residue, ixN15, ic, is) + (1._wp - f_soluable_residue) &
          &                                        * (1._wp - n_res2dom(ic,is)) * hlp1(ic,is)
      END DO
    END DO
  END SUBROUTINE calc_microbial_turnover

  ! ======================================================================================================= !
  !>calculates fraction of N deficit of microbial N growth from asymbiotic N fixation
  !>  given limiting functions from soil temperature, moisture and soil nutrient status
  !>
  ELEMENTAL FUNCTION calc_asymb_bnf_fraction( &
    & rtm_asymb_bnf, &
    & rmm_asymb_bnf, &
    & nh4_solute, &
    & no3_solute) &
    & RESULT (asymb_bnf_fraction)

    USE mo_sb_constants,          ONLY: km_asymb_bnf, vmax_asymb_bnf
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(in)             :: rtm_asymb_bnf                 !< temperature rate modifier []
    REAL(wp), INTENT(in)             :: rmm_asymb_bnf                 !< moisture rate modifier []
    REAL(wp), INTENT(in)             :: nh4_solute                    !< NH4 in solution [mol / m3]
    REAL(wp), INTENT(in)             :: no3_solute                    !< NO3 in solution [mol / m3]
    REAL(wp)                         :: asymb_bnf_fraction            !< fraction of N deficit satisfied by N fixation
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                         :: f_nav
    CHARACTER(len=*), PARAMETER      :: routine = TRIM(modname)//':calc_asymb_bnf_fraction'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 rate modifier for mineral soil nitrogen
    !>
    IF (nh4_solute + no3_solute > eps8) THEN
      f_nav = MIN(1.0_wp, 1.0_wp - (nh4_solute + no3_solute) / (nh4_solute + no3_solute + km_asymb_bnf))
    ELSE
      f_nav = 1._wp
    ENDIF
    ! ...
    asymb_bnf_fraction = vmax_asymb_bnf * rtm_asymb_bnf * rmm_asymb_bnf * f_nav
  END FUNCTION calc_asymb_bnf_fraction

  ! ======================================================================================================= !
  !>calculate growth of the microbial pool
  !>
  !>  corresponds to Ahrens et al. 2015, eq. 2
  !>  extended to account for N and P limitation of microbial growth given a prescribed stoichometry
  !>  enzyme partitioning added according to Wutzler et al. 2017, by Lin
  !>
  !>  Input: DOM and microbial pools, microbial nutrient uptake rate, temperature and moisture rate modifiers
  !>
  !>  Output: DOM consumption, microbial growth, hetertrophic respiration and net mineralisation
  !>
  SUBROUTINE calc_microbial_growth( &
    & nc, &
    & nsoil_w, &
    & dtime, &
    & num_sl_above_bedrock, &
    & bnf_scheme, &
    & rtm_mic_uptake, &
    & rmm_mic_uptake, &
    & rtm_hsc, &
    & rmm_hsc, &
    & microbial_cue_eff_tmic_mavg, &
    & asymb_n_fixation_rel_rate, &
    & sb_pool_mt, &
    & microbial_uptake_nh4_sl, &
    & microbial_uptake_nh4_n15_sl, &
    & microbial_uptake_no3_sl, &
    & microbial_uptake_no3_n15_sl, &
    & microbial_uptake_po4_sl, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & net_mineralisation_nh4, &
    & net_mineralisation_nh4_n15, &
    & net_mineralisation_no3, &
    & net_mineralisation_no3_n15, &
    & gross_mineralisation_po4, &
    & net_mineralisation_po4, &
    & microbial_cue_eff, &
    & microbial_nue_eff, &
    & microbial_pue_eff, &
    & asymb_n_fixation, &
    & fact_n_status_mic_c_growth, &
    & fact_p_status_mic_c_growth, &
    & sb_loss_mt, &
    & sb_formation_mt )

    USE mo_isotope_util,          ONLY: calc_fractionation, calc_mixing_ratio_N15N14
    USE mo_sb_constants,          ONLY: vmax_mic_uptake, km_mic_uptake, microbial_cue, microbial_nue, microbial_pue, &
      &                                 min_mic_biomass, microbial_cn, microbial_np, microbial_cue_min, microbial_cue_max, &
      &                                 eta_ammonification
    USE mo_veg_constants,         ONLY: eta_nfixation
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc                           !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_w                      !< number of soil layers
    REAL(wp),                           INTENT(in)    :: dtime                        !< timestep length
    REAL(wp), DIMENSION(nc),            INTENT(in)    :: num_sl_above_bedrock         !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    CHARACTER(len=*),                   INTENT(in)    :: bnf_scheme                   !< SB_ asymbiotic N fixation scheme
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rtm_mic_uptake               !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rmm_mic_uptake               !< moisture rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rtm_hsc                      !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rmm_hsc                      !< moisture rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: microbial_cue_eff_tmic_mavg  !< microbial effective CUE given NP constraints
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: asymb_n_fixation_rel_rate    !< asymb N fixation rate relative to demand []
    REAL(wp),                           INTENT(in)    :: sb_pool_mt(:,:,:,:)          !< bgcm sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_uptake_nh4_sl      !< microbial uptake rate [micro-mol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_uptake_nh4_n15_sl  !< microbial uptake rate [micro-mol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_uptake_no3_sl      !< microbial uptake rate [micro-mol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_uptake_no3_n15_sl  !< microbial uptake rate [micro-mol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_uptake_po4_sl      !< microbial uptake rate [micro-mol P/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: het_respiration              !< respiration rate [micro-mol C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: het_respiration_c13          !< respiration rate [micro-mol 13C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: het_respiration_c14          !< respiration rate [micro-mol 14C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: net_mineralisation_nh4       !< net mineralisation rate [micro-mol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: net_mineralisation_nh4_n15   !< net mineralisation rate [micro-mol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: net_mineralisation_no3       !< net mineralisation rate [micro-mol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: net_mineralisation_no3_n15   !< net mineralisation rate [micro-mol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: gross_mineralisation_po4     !< gross mienralization rate by DOM , [micro-mol P/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: net_mineralisation_po4       !< net mineralisation rate [micro-mol P/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_cue_eff            !< current effective CUE given NP constraints
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_nue_eff            !< current effective NUE given NP constraints
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: microbial_pue_eff            !< current effective PUE given NP constraints
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: asymb_n_fixation             !< current effective CUE given NP constraints
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: fact_n_status_mic_c_growth   !< scalor of N availability on microbial C growth
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(inout) :: fact_p_status_mic_c_growth   !< scalor of P availability on microbial C growth
    REAL(wp),                           INTENT(inout) :: sb_loss_mt(:,:,:,:)          !< bgcm flux: sb_loss [mol/m3/timestep]
    REAL(wp),                           INTENT(inout) :: sb_formation_mt(:,:,:,:)     !< bgcm flux: sb_formation [mol/m3/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: hlp1
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: hlp2
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: hlp3
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: growth
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: growth_c
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: growth_n
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: growth_p
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: f_uptake
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: uptake_dom
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: uptake_c
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: uptake_n
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: uptake_p
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: scal
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: gross_mineralisation_nh4
    REAL(wp), DIMENSION(nc, nsoil_w)                 :: gross_mineralisation_nh4_n15
    INTEGER                                          :: ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_growth'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        !>1.0 potential microbial DOM consumption [micro-mol / m3 / s]
        !>
        hlp1(ic,is)     = sb_pool_mt(ix_dom, ixC, ic, is) / (km_mic_uptake * rtm_hsc(ic,is) &
          &               * rmm_hsc(ic,is) + sb_pool_mt(ix_dom, ixC, ic, is) )
        uptake_c(ic,is) = vmax_mic_uptake * rtm_mic_uptake(ic,is) * rmm_mic_uptake(ic,is) &
          &               * hlp1(ic,is) * sb_pool_mt(ix_microbial, ixC, ic, is)
        IF (sb_pool_mt(ix_dom, ixC, ic, is) > eps8) THEN
          f_uptake(ic,is) = MIN(uptake_c(ic,is) / sb_pool_mt(ix_dom, ixC, ic, is), 1._wp * 1.E6_wp / dtime)
        ELSE
          f_uptake(ic,is) = 0.0_wp
        END IF
        uptake_n(ic,is) = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixN, ic, is)
        uptake_p(ic,is) = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixP, ic, is)

        !>  1.1 free-living N fixation rate supporting microbial growth [micro-mol / m3 / s]
        !>
        SELECT CASE(TRIM(bnf_scheme))
        ! for the SB_ asymbiotic N fixation no difference is made between "dynamic" and "fixed"
        !  ("fixed" is only available in the VEG_ symbiotic N fixation)
        CASE("dynamic", "fixed")
          asymb_n_fixation(ic,is) = asymb_n_fixation_rel_rate(ic,is) * &
            &                     MAX(0.0_wp, microbial_cue_max * uptake_c(ic,is) / microbial_cn - &
            &                     microbial_nue * f_uptake(ic,is) * sb_pool_mt(ix_dom, ixN, ic, is))
        CASE("unlimited")
          asymb_n_fixation(ic,is) = MAX(0.0_wp, microbial_cue_max * uptake_c(ic,is) / microbial_cn - &
            &                     microbial_nue * f_uptake(ic,is) * sb_pool_mt(ix_dom, ixN, ic, is) - &
            &                     microbial_uptake_nh4_sl(ic,is) - microbial_uptake_no3_sl(ic,is))
        CASE("none")
          asymb_n_fixation(ic,is) = 0.0_wp
        CASE DEFAULT
          CALL finish(TRIM(routine),'No such bnf_scheme available here.')
        END SELECT

        !>  1.2 microbial growth rates in [micro-mol C / m3 / s], given by C, N and P availability, and the rate-limiting growth rate
        !>
        growth_c(ic,is) = microbial_cue_max * uptake_c(ic,is)
        growth_n(ic,is) = (microbial_nue * uptake_n(ic,is) + asymb_n_fixation(ic,is) + &
          &               microbial_uptake_nh4_sl(ic,is) + microbial_uptake_no3_sl(ic,is)) * microbial_cn
        growth_p(ic,is) = (microbial_pue * uptake_p(ic,is) + &
          &               microbial_uptake_po4_sl(ic,is)) * microbial_np * microbial_cn

        IF (growth_c(ic,is) > eps8 .AND. sb_pool_mt(ix_microbial, ixC, ic, is) > (min_mic_biomass + eps8)) THEN
          fact_n_status_mic_c_growth(ic,is) = growth_n(ic,is) / growth_c(ic,is)
          fact_p_status_mic_c_growth(ic,is) = growth_p(ic,is) / growth_c(ic,is)
        ELSE
          fact_n_status_mic_c_growth(ic,is) = 1._wp
          fact_p_status_mic_c_growth(ic,is) = 1._wp
        END IF
        IF (growth_c(ic,is) <= growth_n(ic,is) .AND. growth_c(ic,is) <= growth_p(ic,is)) THEN
          ! carbon uptake is limiting
          growth(ic,is) = growth_c(ic,is)
        ELSE IF (growth_n(ic,is) < growth_c(ic,is) .AND. growth_n(ic,is) <= growth_p(ic,is)) THEN
          ! nitrogen uptake is limiting
          growth(ic,is) = growth_n(ic,is)
        ELSE IF (growth_p(ic,is) < growth_c(ic,is) .AND. growth_p(ic,is) < growth_n(ic,is)) THEN
          ! phosphorus uptake is limiting
          growth(ic,is) = growth_p(ic,is)
        END IF

        !>2.0 Adjust carbon-use efficiency such that the growth is stoichiometrically balanced
        !>
        microbial_cue_eff(ic,is) = microbial_cue_max
        IF (growth(ic,is) < growth_c(ic,is) .AND. uptake_c(ic,is) > eps8) THEN
          ! excess carbon is going to be respired
          microbial_cue_eff(ic,is) = MAX(microbial_cue_min,growth(ic,is) / uptake_c(ic,is))
        ELSE IF (growth(ic,is) < growth_c(ic,is)) THEN
          microbial_cue_eff(ic,is) = microbial_cue_min
        END IF

        !>  2.1 actual growth given actual CUE, adjust uptake rate of DOM to the most
        !>     limiting nutrient
        !>
        ! growth_n and growth_p are unchanged, only growth_c needs to be adjusted to current CUE
        growth_c(ic,is) = microbial_cue_eff_tmic_mavg(ic,is) * uptake_c(ic,is)
        ! ensure asymb_n_fixation does not exceed uptake demand
        IF (asymb_n_fixation(ic,is) > growth_c(ic,is) / microbial_cn - microbial_nue * uptake_n(ic,is)) THEN
          asymb_n_fixation(ic,is) = MAX(0.0_wp, growth_c(ic,is) / microbial_cn - microbial_nue * uptake_n(ic,is))
        END IF

        IF (growth_c(ic,is) <= growth_n(ic,is) .AND. growth_c(ic,is) <= growth_p(ic,is)) THEN
          ! carbon uptake is limiting, full uptake rate is realised
          f_uptake(ic,is) = f_uptake(ic,is)
        ELSE IF (growth_n(ic,is) < growth_c(ic,is) .AND. growth_n(ic,is) <= growth_p(ic,is) .AND. &
                 (microbial_cue_eff_tmic_mavg(ic,is) * sb_pool_mt(ix_dom, ixC, ic, is) - microbial_nue &
                  * sb_pool_mt(ix_dom, ixN, ic, is) * microbial_cn) > eps8) THEN
          ! nitrogen uptake is limiting
          ! calculate uptake rate that satisfies N constraint and CUE
          f_uptake(ic,is) = MAX(0.0_wp,MIN(f_uptake(ic,is), &
            &   ((asymb_n_fixation(ic,is) + microbial_uptake_nh4_sl(ic,is) + microbial_uptake_no3_sl(ic,is)) * microbial_cn) / &
            &   (microbial_cue_eff_tmic_mavg(ic,is) * sb_pool_mt(ix_dom, ixC, ic, is) - &
            &    microbial_nue * sb_pool_mt(ix_dom, ixN, ic, is) * microbial_cn)))
        ELSE IF (growth_p(ic,is) < growth_c(ic,is) .AND. growth_p(ic,is) < growth_n(ic,is) .AND. &
                 (microbial_cue_eff_tmic_mavg(ic,is) * sb_pool_mt(ix_dom, ixC, ic, is) - microbial_pue &
                 * sb_pool_mt(ix_dom, ixP, ic, is) * microbial_np * microbial_cn) > eps8) THEN
          ! phosphorus uptake is limiting
          ! calculate uptake rate that satisfies P constraint and CUE
          f_uptake(ic,is) = MAX(0.0_wp,MIN(f_uptake(ic,is), &
            &             (microbial_uptake_po4_sl(ic,is) * microbial_np * microbial_cn) / &
            &             (microbial_cue_eff_tmic_mavg(ic,is) * sb_pool_mt(ix_dom, ixC, ic, is) - &
            &              microbial_pue * sb_pool_mt(ix_dom, ixP, ic, is) * microbial_np * microbial_cn)))
        ELSE
          f_uptake(ic,is) = 0.0_wp
        END IF

        growth_c(ic,is) = microbial_cue_eff_tmic_mavg(ic,is) * f_uptake(ic,is) * sb_pool_mt(ix_dom, ixC, ic, is)
        uptake_n(ic,is) = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixN, ic, is)
        uptake_p(ic,is) = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixP, ic, is)
        growth_n(ic,is) = (microbial_nue * uptake_n(ic,is) + asymb_n_fixation(ic,is) + &
                         microbial_uptake_nh4_sl(ic,is) + microbial_uptake_no3_sl(ic,is)) * microbial_cn
        growth_p(ic,is) = (microbial_pue * uptake_p(ic,is) + &
                         microbial_uptake_po4_sl(ic,is)) * microbial_np * microbial_cn

        !>  2.2 Adjust nitrogen-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
        !>     and calculate gross mineralisation rate (and its isotopic signature)
        !>
        microbial_nue_eff(ic,is) = microbial_nue
        IF (growth_c(ic,is) < growth_n(ic,is)) THEN
          ! inorganic N immobilisation flux in agreement with C growth, taking into account current fixation rate
          ! by design, N fixation is not larger than the deficit
          hlp1(ic,is) = MAX(growth_c(ic,is) / microbial_cn - microbial_nue * uptake_n(ic,is) - asymb_n_fixation(ic,is),0.0_wp)
          IF ((microbial_uptake_nh4_sl(ic,is) + microbial_uptake_no3_sl(ic,is)) > eps8) THEN
            scal(ic,is) = hlp1(ic,is) / (microbial_uptake_nh4_sl(ic,is) + microbial_uptake_no3_sl(ic,is))
          ELSE
            scal(ic,is) = 0.0_wp
          END IF
          microbial_uptake_nh4_sl(ic,is)     = microbial_uptake_nh4_sl(ic,is)     * scal(ic,is)
          microbial_uptake_nh4_n15_sl(ic,is) = microbial_uptake_nh4_n15_sl(ic,is) * scal(ic,is)
          microbial_uptake_no3_sl(ic,is)     = microbial_uptake_no3_sl(ic,is)     * scal(ic,is)
          microbial_uptake_no3_n15_sl(ic,is) = microbial_uptake_no3_n15_sl(ic,is) * scal(ic,is)
          IF (uptake_n(ic,is) > eps8) THEN
            microbial_nue_eff(ic,is) = MIN(microbial_nue, &
              &                      ( growth_c(ic,is) / &
              &                        microbial_cn - microbial_uptake_nh4_sl(ic,is) - microbial_uptake_no3_sl(ic,is) ) / &
              &                      uptake_n(ic,is))
          END IF
        END IF
        gross_mineralisation_nh4(ic,is)     = (1._wp - microbial_nue_eff(ic,is)) * uptake_n(ic,is)
        gross_mineralisation_nh4_n15(ic,is) = gross_mineralisation_nh4(ic,is) * calc_fractionation( &
          &                                                                                sb_pool_mt(ix_dom, ixN, ic, is), &
          &                                                                                sb_pool_mt(ix_dom, ixN15, ic, is), &
          &                                                                                eta_ammonification)

        !>  2.3 Adjust phosphorus-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
        !>     and calculate gross mineralisation rate
        !>
        microbial_pue_eff(ic,is) = microbial_pue
        IF (growth_c(ic,is) < growth_p(ic,is)) THEN
          ! inorganic P immobilisation flux, in agreement with C and N uptake
          microbial_uptake_po4_sl(ic,is) = MAX(growth_c(ic,is) / microbial_cn / microbial_np &
                                           - microbial_pue * uptake_p(ic,is),0.0_wp)
          ! In the case that the source C:N:P is too tight for microbes, cause gross mineralisation
          IF (uptake_p(ic,is) > eps8) THEN
            microbial_pue_eff(ic,is) = ( growth_c(ic,is) / microbial_cn / microbial_np - microbial_uptake_po4_sl(ic,is) ) &
                                       / uptake_p(ic,is)
          END IF
        END IF
        gross_mineralisation_po4(ic,is) = gross_mineralisation_po4(ic,is) + (1._wp - microbial_pue_eff(ic,is)) * uptake_p(ic,is)

        !>3.0 update fluxes
        !>  note that loss and formation of organic matter is in [mol / m3 / timestep]
        !>  whereas fluxes are in [micro-mol / m3 / s]
        !>
        ! carbon
        uptake_dom(ic,is)                          = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixC, ic, is)
        sb_loss_mt(ix_dom, ixC, ic, is)            = sb_loss_mt(ix_dom, ixC, ic, is) + uptake_dom(ic,is) * dtime / 1000000._wp
        sb_formation_mt(ix_microbial, ixC, ic, is) = sb_formation_mt(ix_microbial, ixC, ic, is) &
          &                                        + microbial_cue_eff_tmic_mavg(ic,is) * uptake_dom(ic,is) * dtime / 1000000._wp
        het_respiration(ic,is)                     = het_respiration(ic,is) &
          &                                        + (1._wp - microbial_cue_eff_tmic_mavg(ic,is)) * uptake_dom(ic,is)
        ! nitrogen
        net_mineralisation_nh4(ic,is)              = gross_mineralisation_nh4(ic,is) - microbial_uptake_nh4_sl(ic,is)
        net_mineralisation_no3(ic,is)              = - microbial_uptake_no3_sl(ic,is)
        uptake_dom(ic,is)                          = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixN, ic, is)
        sb_loss_mt(ix_dom, ixN, ic, is)            = sb_loss_mt(ix_dom, ixN, ic, is) + uptake_dom(ic,is) * dtime / 1000000.0_wp
        sb_formation_mt(ix_microbial, ixN, ic, is) = sb_formation_mt(ix_microbial, ixN, ic, is) &
          &                                        + (microbial_nue_eff(ic,is) * uptake_dom(ic,is) &
          &                                        + asymb_n_fixation(ic,is) &
          &                                        + microbial_uptake_nh4_sl(ic,is) + microbial_uptake_no3_sl(ic,is)) &
          &                                        * dtime / 1000000.0_wp
        ! phosphorus
        net_mineralisation_po4(ic,is)              = gross_mineralisation_po4(ic,is) - microbial_uptake_po4_sl(ic,is)
        uptake_dom(ic,is)                          = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixP, ic, is)
        sb_loss_mt(ix_dom, ixP, ic, is)            = sb_loss_mt(ix_dom, ixP, ic, is) + uptake_dom(ic,is) * dtime / 1000000.0_wp
        sb_formation_mt(ix_microbial, ixP, ic, is) = sb_formation_mt(ix_microbial, ixP, ic, is) &
          &                                        + (microbial_pue_eff(ic,is) * uptake_dom(ic,is) &
          &                                        + microbial_uptake_po4_sl(ic,is)) * dtime / 1000000.0_wp
        ! carbon13
        uptake_dom(ic,is)                            = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixC13, ic, is)
        sb_loss_mt(ix_dom, ixC13, ic, is)            = sb_loss_mt(ix_dom, ixC13, ic, is) + uptake_dom(ic,is) * dtime / 1000000._wp
        sb_formation_mt(ix_microbial, ixC13, ic, is) = sb_formation_mt(ix_microbial, ixC13, ic, is) &
          &                                          + microbial_cue_eff_tmic_mavg(ic,is) * uptake_dom(ic,is) &
          &                                          * dtime / 1000000.0_wp
        het_respiration_c13(ic,is)                   = het_respiration_c13(ic,is) + (1._wp - microbial_cue_eff_tmic_mavg(ic,is)) &
          &                                          * uptake_dom(ic,is)
        ! carbon14
        uptake_dom(ic,is)                            = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixC14, ic, is)
        sb_loss_mt(ix_dom, ixC14, ic, is)            = sb_loss_mt(ix_dom, ixC14, ic, is) + uptake_dom(ic,is) * dtime / 1000000.0_wp
        sb_formation_mt(ix_microbial, ixC14, ic, is) = sb_formation_mt(ix_microbial, ixC14, ic, is) &
          &                                            + microbial_cue_eff_tmic_mavg(ic,is) * uptake_dom(ic,is) &
          &                                            * dtime / 1000000.0_wp
        het_respiration_c14(ic,is)                   = het_respiration_c14(ic,is) + (1._wp - microbial_cue_eff_tmic_mavg(ic,is) ) &
          &                                          * uptake_dom(ic,is)
        ! nitrogen15
        net_mineralisation_nh4_n15(ic,is)            = gross_mineralisation_nh4_n15(ic,is) - microbial_uptake_nh4_n15_sl(ic,is)
        net_mineralisation_no3_n15(ic,is)            = - microbial_uptake_no3_n15_sl(ic,is)
        uptake_dom(ic,is)                            = f_uptake(ic,is) * sb_pool_mt(ix_dom, ixN15, ic, is)
        sb_loss_mt(ix_dom, ixN15, ic, is)            = sb_loss_mt(ix_dom, ixN15, ic, is) + uptake_dom(ic,is) * dtime / 1000000.0_wp
        sb_formation_mt(ix_microbial, ixN15, ic, is) = sb_formation_mt(ix_microbial, ixN15, ic, is) &
          &                                          + (microbial_nue * uptake_dom(ic,is) &
          &                                          + asymb_n_fixation(ic,is) &
          &                                          / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(- eta_nfixation)) &
          &                                          + microbial_uptake_nh4_n15_sl(ic,is) + microbial_uptake_no3_n15_sl(ic,is)) &
          &                                          * dtime / 1000000.0_wp
      END DO
    END DO
  END SUBROUTINE calc_microbial_growth

  ! ======================================================================================================= !
  !>calculate depolymerisation rate of polymeric litter and microbial residue
  !>  corresponds Ahrens et al. 2015, eq. 1
  !>
  !>  improves by including microbial residue as alternative polymeric pool
  !>  enzyme partitioning simulated using SESAM approach, by Wutzler et al. 2017
  !>
  !>  Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for
  !>         depolymerisation and half-saturation, polymeric/residue pool and microbial biomass
  !>
  !>  Output: transfer of residue/polymeric pool to DOM
  !>
  SUBROUTINE calc_depolymerisation_SESAM( &
    & nc, nsoil_w, dtime, num_sl_above_bedrock, &
    & rtm_hsc,rmm_hsc,rtm_depolymerisation,rmm_depolymerisation, &
    & enzyme_depoly, microbial_cue_mavg, microbial_nue_mavg, microbial_pue_mavg, &
    & fact_n_status_mic_c_growth, fact_p_status_mic_c_growth, &
    & sb_pool_mt, &
    & enzyme_frac_poly, enzyme_frac_residue, enzyme_frac_AP, &
    & enzyme_frac_poly_c, enzyme_frac_poly_n, enzyme_frac_poly_p, &
    & enzyme_frac_poly_c_mavg, enzyme_frac_poly_n_mavg, enzyme_frac_poly_p_mavg, &
    & sb_loss_mt, &
    & sb_formation_mt)

    USE mo_sb_constants
    USE mo_lnd_time_averages,               ONLY: mavg_period_tenzyme
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                           !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_w                      !< number of soil layers
    REAL(wp),                          INTENT(in)    :: dtime                        !< timestep length
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock         !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: rtm_hsc                      !< temperature rate modifier for half-saturation constants
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: rmm_hsc                      !< moisture rate modifier for half-saturation constants
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: rtm_depolymerisation         !< temperature rate modifier for depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: rmm_depolymerisation         !< moisture rate modifier for depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: enzyme_depoly                !< total enzyme for depolymerisation [mol C/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: microbial_cue_mavg           !< microbial cue moving average
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: microbial_nue_mavg           !< microbial nue moving average
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: microbial_pue_mavg           !< microbial pue moving average
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: fact_n_status_mic_c_growth   !< scalor of N availability on microbial C growth
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: fact_p_status_mic_c_growth   !< scalor of P availability on microbial C growth
    REAL(wp),                          INTENT(in)    :: sb_pool_mt(:,:,:,:)          !< bgcm sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly             !< fraction of enzyme for poly litter depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_residue          !< fraction of enzyme for microbial residue depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_AP               !< fraction of enzyme for phosphatase
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly_c           !< temporal fraction of enzyme for poly litter depolymerisation, C
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly_n           !< temporal fraction of enzyme for poly litter depolymerisation, N
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly_p           !< temporal fraction of enzyme for poly litter depolymerisation, P
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly_c_mavg      !< mavg of enzyme for poly litter depolymerisation, C
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly_n_mavg      !< mavg of enzyme for poly litter depolymerisation, N
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(inout) :: enzyme_frac_poly_p_mavg      !< mavg of enzyme for poly litter depolymerisation, P
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)          !< bgcm sb_loss [mol/m3/timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)     !< bgcm sb_formation [mol/m3/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_w)    :: hlp1
    REAL(wp), DIMENSION(nc, nsoil_w)    :: hlp2
    REAL(wp), DIMENSION(nc, nsoil_w)    :: hlp3
    REAL(wp), DIMENSION(nc, nsoil_w)    :: poly_depolymerisation_rate
    REAL(wp), DIMENSION(nc, nsoil_w)    :: residue_depolymerisation_rate
    REAL(wp), DIMENSION(nc, nsoil_w)    :: depoly_c
    REAL(wp), DIMENSION(nc, nsoil_w)    :: depoly_n
    REAL(wp), DIMENSION(nc, nsoil_w)    :: depoly_p
    REAL(wp), DIMENSION(nc, nsoil_w)    :: hlp_depoly
    REAL(wp), DIMENSION(nc, nsoil_w)    :: enz_AP_direction
    INTEGER                             :: ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation_SESAM'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init local variables
    !>
    depoly_c(:,:) = 0.0_wp
    depoly_n(:,:) = 0.0_wp
    depoly_p(:,:) = 0.0_wp

    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        !>1.0 calculate polymeric litter depolymerisation rate
        !>
        hlp1(ic,is) = 1._wp
        hlp1(ic,is) = enzyme_depoly(ic,is) * enzyme_frac_poly(ic,is) &
          &         / (km_depolymerisation * rtm_hsc(ic,is) * rmm_hsc(ic,is) + enzyme_depoly(ic,is) * enzyme_frac_poly(ic,is))
        hlp2(ic,is) = vmax_depolymerisation_litter * rtm_depolymerisation(ic,is) * rmm_depolymerisation(ic,is) &
          &           * hlp1(ic,is) * dtime &
          &           * sb_pool_mt(ix_polymeric_litter, ixC, ic, is)

        IF (sb_pool_mt(ix_polymeric_litter, ixC, ic, is) > eps8) THEN
          poly_depolymerisation_rate(ic,is) = hlp2(ic,is) / sb_pool_mt(ix_polymeric_litter, ixC, ic, is)
        ELSE
          poly_depolymerisation_rate(ic,is) = 0.0_wp
        END IF

        !>2.0 calculate microbial residue depolymerisation rate
        !>
        hlp1(ic,is) = 1._wp
        hlp1(ic,is) = enzyme_depoly(ic,is) * enzyme_frac_residue(ic,is) &
          &         / (km_depolymerisation * rtm_hsc(ic,is) * rmm_hsc(ic,is) + enzyme_depoly(ic,is) * enzyme_frac_residue(ic,is))
        hlp2(:,:) = vmax_depolymerisation_residue * rtm_depolymerisation(ic,is) &
          &         * rmm_depolymerisation(ic,is) * hlp1(ic,is) * dtime &
          &         * sb_pool_mt(ix_residue, ixC, ic, is)

        IF (sb_pool_mt(ix_residue, ixC, ic, is) > eps8) THEN
          residue_depolymerisation_rate(ic,is) = hlp2(ic,is) / sb_pool_mt(ix_residue, ixC, ic, is)
        ELSE
          residue_depolymerisation_rate(ic,is) = 0.0_wp
        END IF

        !>3.0 calculate the loss and formation of all elements in polymeric litter pool
        !>
        ! C
        hlp2(ic,is)                                  = poly_depolymerisation_rate(ic,is) &
                                                       * sb_pool_mt(ix_polymeric_litter, ixC, ic, is)
        sb_loss_mt(ix_polymeric_litter, ixC, ic, is) = sb_loss_mt(ix_polymeric_litter, ixC, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixC, ic, is)         = sb_formation_mt(ix_dom, ixC, ic, is) + hlp2(ic,is)
        depoly_c(ic,is)                              = depoly_c(ic,is) + hlp2(ic,is)
        ! N
        hlp2(ic,is)                                  = poly_depolymerisation_rate(ic,is) &
                                                       * sb_pool_mt(ix_polymeric_litter, ixN, ic, is)
        sb_loss_mt(ix_polymeric_litter, ixN, ic, is) = sb_loss_mt(ix_polymeric_litter, ixN, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixN, ic, is)         = sb_formation_mt(ix_dom, ixN, ic, is) + hlp2(ic,is)
        depoly_n(ic,is)                              = depoly_n(ic,is) + hlp2(ic,is)
        ! P
        hlp2(ic,is)                                  = poly_depolymerisation_rate(ic,is) &
                                                       * sb_pool_mt(ix_polymeric_litter, ixP, ic, is)
        sb_loss_mt(ix_polymeric_litter, ixP, ic, is) = sb_loss_mt(ix_polymeric_litter, ixP, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixP, ic, is)         = sb_formation_mt(ix_dom, ixP, ic, is) + hlp2(ic,is)
        depoly_p(ic,is)                              = depoly_p(ic,is) + hlp2(ic,is)
        ! C13
        hlp2(ic,is)                                    = poly_depolymerisation_rate(ic,is) &
                                                         * sb_pool_mt(ix_polymeric_litter, ixC13, ic, is)
        sb_loss_mt(ix_polymeric_litter, ixC13, ic, is) = sb_loss_mt(ix_polymeric_litter, ixC13, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixC13, ic, is)         = sb_formation_mt(ix_dom, ixC13, ic, is) + hlp2(ic,is)
        ! C14
        hlp2(ic,is)                                    = poly_depolymerisation_rate(ic,is) &
                                                         * sb_pool_mt(ix_polymeric_litter, ixC14, ic, is)
        sb_loss_mt(ix_polymeric_litter, ixC14, ic, is) = sb_loss_mt(ix_polymeric_litter, ixC14, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixC14, ic, is)         = sb_formation_mt(ix_dom, ixC14, ic, is) + hlp2(ic,is)
        ! N15
        hlp2(ic,is)                                    = poly_depolymerisation_rate(ic,is) &
                                                         * sb_pool_mt(ix_polymeric_litter, ixN15, ic, is)
        sb_loss_mt(ix_polymeric_litter, ixN15, ic, is) = sb_loss_mt(ix_polymeric_litter, ixN15, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixN15, ic, is)         = sb_formation_mt(ix_dom, ixN15, ic, is) + hlp2(ic,is)

        !>4.0 calculate the loss and formation of all elements in microbial residue pool
        !>
        ! C
        hlp2(ic,is)                          = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixC, ic, is)
        sb_loss_mt(ix_residue, ixC, ic, is)  = sb_loss_mt(ix_residue, ixC, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixC, ic, is) = sb_formation_mt(ix_dom, ixC, ic, is) + hlp2(ic,is)
        depoly_c(ic,is)                      = depoly_c(ic,is) + hlp2(ic,is)
        ! N
        hlp2(ic,is)                          = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixN, ic, is)
        sb_loss_mt(ix_residue, ixN, ic, is)  = sb_loss_mt(ix_residue, ixN, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixN, ic, is) = sb_formation_mt(ix_dom, ixN, ic, is) + hlp2(ic,is)
        depoly_n(ic,is)                      = depoly_n(ic,is) + hlp2(ic,is)
        ! P
        hlp2(ic,is)                          = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixP, ic, is)
        sb_loss_mt(ix_residue, ixP, ic, is)  = sb_loss_mt(ix_residue, ixP, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixP, ic, is) = sb_formation_mt(ix_dom, ixP, ic, is) + hlp2(ic,is)
        depoly_p(ic,is)                      = depoly_p(ic,is) + hlp2(ic,is)
        ! C13
        hlp2(ic,is)                            = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixC13, ic, is)
        sb_loss_mt(ix_residue, ixC13, ic, is)  = sb_loss_mt(ix_residue, ixC13, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixC13, ic, is) = sb_formation_mt(ix_dom, ixC13, ic, is) + hlp2(ic,is)
        ! C14
        hlp2(ic,is)                            = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixC14, ic, is)
        sb_loss_mt(ix_residue, ixC14, ic, is)  = sb_loss_mt(ix_residue, ixC14, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixC14, ic, is) = sb_formation_mt(ix_dom, ixC14, ic, is) + hlp2(ic,is)
        ! N15
        hlp2(ic,is)                            = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixN15, ic, is)
        sb_loss_mt(ix_residue, ixN15, ic, is)  = sb_loss_mt(ix_residue, ixN15, ic, is) + hlp2(ic,is)
        sb_formation_mt(ix_dom, ixN15, ic, is) = sb_formation_mt(ix_dom, ixN15, ic, is) + hlp2(ic,is)

        !>5.0 calculate the temporal enzyme partitioning for poly litter based on C, N, and P
        !>
        !> following the algorithm described in Wutzler et al. 2025 (DOI:10.5194/gmd-15-8377-2022)
        !>
        IF (poly_depolymerisation_rate(ic,is) > 0.0_wp .AND. residue_depolymerisation_rate(ic,is)> 0.0_wp) THEN
          !! alpha_C, poly enzyme allocation coefficient based on carbon
          hlp1(ic,is)               = poly_depolymerisation_rate(ic,is) * sb_pool_mt(ix_polymeric_litter, ixC, ic, is)
          hlp2(ic,is)               = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixC, ic, is)
          hlp3(ic,is)               = km_depolymerisation * rtm_hsc(ic,is) * rmm_hsc(ic,is)
          enzyme_frac_poly_c(ic,is) = calc_enzyme_allocation(hlp1(ic,is), hlp2(ic,is), hlp3(ic,is), enzyme_depoly(ic,is))
          !! alpha_N, poly enzyme allocation coefficient based on nitrogen
          hlp1(ic,is)               = poly_depolymerisation_rate(ic,is) * sb_pool_mt(ix_polymeric_litter, ixN, ic, is)
          hlp2(ic,is)               = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixN, ic, is)
          enzyme_frac_poly_n(ic,is) = calc_enzyme_allocation(hlp1(ic,is), hlp2(ic,is), hlp3(ic,is), enzyme_depoly(ic,is))
          !! alpha_P, poly enzyme allocation coefficient based on phosphorus
          hlp1(ic,is)               = poly_depolymerisation_rate(ic,is) * sb_pool_mt(ix_polymeric_litter, ixP, ic, is)
          hlp2(ic,is)               = residue_depolymerisation_rate(ic,is) * sb_pool_mt(ix_residue, ixP, ic, is)
          enzyme_frac_poly_p(ic,is) = calc_enzyme_allocation(hlp1(ic,is), hlp2(ic,is), hlp3(ic,is), enzyme_depoly(ic,is))
        END IF

        !>6.0 calculate the enzyme partitioning based on microbial nutrition status
        !>
        hlp1(ic,is) = depoly_c(ic,is) * microbial_cue_mavg(ic,is)
        hlp2(ic,is) = depoly_n(ic,is) * microbial_cn * microbial_nue_mavg(ic,is)
        hlp3(ic,is) = depoly_p(ic,is) * microbial_cn * microbial_np * microbial_pue_mavg(ic,is)
        enz_AP_direction(ic,is) = 1._wp
        IF (hlp1(ic,is) > 0.0_wp .AND. hlp2(ic,is) > 0.0_wp .AND. hlp3(ic,is) > 0.0_wp) THEN
          IF (fact_n_status_mic_c_growth(ic,is) > 1._wp .AND. fact_p_status_mic_c_growth(ic,is) > 1._wp) THEN
            hlp_depoly(ic,is)       = enzyme_frac_poly_c_mavg(ic,is)
          ELSE IF (fact_n_status_mic_c_growth(ic,is) > 1._wp .AND. fact_p_status_mic_c_growth(ic,is) <= 1._wp) THEN
            hlp_depoly(ic,is)       = enzyme_frac_poly_c_mavg(ic,is)
            enz_AP_direction(ic,is) = -1._wp
          ELSE IF (fact_n_status_mic_c_growth(ic,is) <= 1._wp .AND. fact_p_status_mic_c_growth(ic,is) > 1._wp) THEN
            hlp_depoly(ic,is)       = enzyme_frac_poly_n_mavg(ic,is)
          ELSE
            hlp_depoly(ic,is)       = enzyme_frac_poly_n_mavg(ic,is)
            IF (fact_n_status_mic_c_growth(ic,is) > fact_p_status_mic_c_growth(ic,is)) THEN
              enz_AP_direction(ic,is) = -1._wp
            END IF
          END IF
        ELSE
          hlp_depoly(ic,is) = enzyme_frac_poly_c_mavg(ic,is)
        END IF
        enzyme_frac_poly(ic,is)    = MIN( MAX(enzyme_frac_poly(ic,is) + (hlp_depoly(ic,is) - enzyme_frac_poly(ic,is)) / &
          &                        mavg_period_tenzyme / one_day * dtime, min_enzyme_fraction), (1._wp - min_enzyme_fraction))
        enzyme_frac_residue(ic,is) = 1._wp - enzyme_frac_poly(ic,is)
        enzyme_frac_AP(ic,is)      = MIN( MAX(enzyme_frac_AP(ic,is) * (1._wp - enz_AP_direction(ic,is) &
          &                          / mavg_period_tenzyme / one_day * dtime), &
          &                          min_enzyme_fraction), max_enzyme_fraction)
      END DO
    END DO
  END SUBROUTINE calc_depolymerisation_SESAM

  ! ======================================================================================================= !
  !>calculates the bulk density of the soil corrected by its organic component
  !>  derived from Federer et al. 1993, Can J. For. Res. 23 1026-1032
  !>
  !>  Input: carbon content by layer
  !>
  !>  Output: corrected bulk density of the soil
  !>
  ELEMENTAL FUNCTION calc_bulk_density_correction( &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & bulk_dens_sl) &
    & RESULT (bulk_dens_corr_sl)

    USE mo_jsb_physical_constants, ONLY: molar_mass_C
    USE mo_sb_constants,           ONLY: rho_bulk_org, carbon_per_dryweight_SOM
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(in)    :: bulk_soil_carbon_sl     !< soil non-litter carbon contributing to density [mol C / m3]
    REAL(wp), INTENT(in)    :: soil_litter_carbon_sl   !< soil litter carbon contributing to density [mol C / m3]
    REAL(wp), INTENT(in)    :: bulk_dens_sl            !< bulk density of the mineral soil [kg / m3]
    REAL(wp)                :: bulk_dens_corr_sl       !< combined bulk density of the mineral and organic soil [kg / m3]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: om_dens_sl       !< organic matter density [kg / m3]
    REAL(wp)                    :: lit_dens_sl      !< litter density [kg / m3]
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_density_correction'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate organic matter density and litter density per soil layer [kg / m3]
    !>
    om_dens_sl  = bulk_soil_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp
    lit_dens_sl = soil_litter_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp

    !>2.0 actual bulk density of the soil
    !>
    bulk_dens_corr_sl = MAX(rho_bulk_org, om_dens_sl + lit_dens_sl + bulk_dens_sl &
      &                 - om_dens_sl * bulk_dens_sl / rho_bulk_org &
      &                 - lit_dens_sl * bulk_dens_sl / rho_bulk_org)
  END FUNCTION calc_bulk_density_correction

  ! ======================================================================================================= !
  !>calculates the maximum sorption capacity for organic material of a soil layer
  !>  derived from COMISSION, the implict assumption here is that no association happens on organic material
  !>
  !>  Input: bulk density of the soil (mineral and organic combined)
  !>
  !>  Output: sorption capacity for organic material
  !>
  ELEMENTAL SUBROUTINE calc_qmax_bulk_density_correction( &
    & bulk_soil_carbon_sl, &
    & volume_min_sl, &
    & bulk_dens_sl, &
    & qmax_om, &
    & qmax_min, &
    & qmax, &
    & qmax_fast, &
    & qmax_slow)

    USE mo_jsb_physical_constants, ONLY: molar_mass_C
    USE mo_sb_constants,           ONLY: carbon_per_dryweight_SOM, rho_bulk_org
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),           INTENT(in)  :: bulk_soil_carbon_sl   !< non-litter organic carbon content of the soil layer [mol C / m3]
    REAL(wp),           INTENT(in)  :: volume_min_sl         !< mineral soil volume fraction the soil layer [m3/ m3]
    REAL(wp),           INTENT(in)  :: bulk_dens_sl          !< bulk density of the mineral soil [kg / m3]
    REAL(wp),           INTENT(in)  :: qmax_om               !< maximium sorption capacity on OM [mol / kg]
    REAL(wp),           INTENT(in)  :: qmax_min              !< maximium sorption capacity on  mineral soil [mol / kg]
    REAL(wp),           INTENT(out) :: qmax                  !< maximum sorption capacity of the soil [mol / m3]
    REAL(wp), OPTIONAL, INTENT(out) :: qmax_fast             !< maximum sorption capacity of fast pool [mol / m3]
    REAL(wp), OPTIONAL, INTENT(out) :: qmax_slow             !< maximum sorption capacity of slow pool [mol / m3]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: vol_fraction_om
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_bulk_density_correction'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate volume fraction filled by organic material
    !>  given assumed density and carbon content of OM and litter
    !>
    vol_fraction_om = MIN(bulk_soil_carbon_sl * molar_mass_C / 1000._wp &
      &               / carbon_per_dryweight_SOM / rho_bulk_org, 1.0_wp)

    !>2.0 maximum sorption capacity per soil volume as the mix of OM and mineral sorption (mol/m3)
    !>
    IF (PRESENT(qmax_fast)) THEN
      qmax_fast = vol_fraction_om * qmax_om  * rho_bulk_org
    ENDIF
    IF (PRESENT(qmax_slow)) THEN
      qmax_slow = volume_min_sl   * qmax_min * bulk_dens_sl
    ENDIF
    qmax = vol_fraction_om * qmax_om  * rho_bulk_org &
      &    + volume_min_sl * qmax_min * bulk_dens_sl
  END SUBROUTINE calc_qmax_bulk_density_correction

  ! ======================================================================================================= !
  !>calculates the km for phosphate of a soil layer
  !>
  !>  Input: default km, solute concentration, qmax of different sorption pool (mineral and organic combined)
  !>
  !>  Output: km of the combined sorption pool
  !>
  ELEMENTAL FUNCTION calc_km_qmax_correction( &
    & solute, &
    & qmax1, &
    & qmax2, &
    & km1, &
    & km2, &
    & w_soil) &
    & RESULT (km)

    REAL(wp),           INTENT(in) :: solute          !< po4 concentration of the soil layer [mol P / m3]
    REAL(wp),           INTENT(in) :: qmax1           !< maximum po4 sorption capacity of the 1st pool [mol P/ m3]
    REAL(wp),           INTENT(in) :: qmax2           !< maximum po4 sorption capacity of the 2nd pool [mol P/ m3]
    REAL(wp),           INTENT(in) :: km1             !< km of the 1st pool [mol P/ m3]
    REAL(wp),           INTENT(in) :: km2             !< km of the 2nd pool [mol P/ m3]
    REAL(wp), OPTIONAL, INTENT(in) :: w_soil          !< water content of the soil layer [L / m3]
    REAL(wp)                       :: km              !< km of the combined pool [mol P/ m3]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: hlp
    REAL(wp)                    :: dq1
    REAL(wp)                    :: dq2
    REAL(wp)                    :: sol1
    REAL(wp)                    :: sol2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_km_qmax_correction'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate the change of sorbed P pool given the solute concentration
    !>
    dq1  = qmax1 * km1 / (km1 + solute) / (km1 + solute)
    dq2  = qmax2 * km2 / (km2 + solute) / (km2 + solute)
    hlp  = dq1 + dq2

    !>2.0 calculate km for the combined P pool
    !>
    ! Qmax*km/(km+solute)^2 = dq1 + dq2, hlp = dq1 + dq2
    ! km = (Qmax/hlp - 2*solute +/- sqrt((Qmax/hlp)^2-4 * Qmax/hlp))/2
    IF (hlp < eps8) THEN
      km = MAX(km1, km2)
    ELSE IF ((km1 + solute) < eps8 .OR. (km2 + solute) < eps8) THEN
      km = MAX(km1, km2)
    ELSE
      hlp  = (qmax1 + qmax2) / hlp
      sol1 = (hlp - 2._wp * solute + SQRT(hlp * hlp - 4._wp * hlp * solute)) / 2._wp
      sol2 = (hlp - 2._wp * solute - SQRT(hlp * hlp - 4._wp * hlp * solute)) / 2._wp
      IF (sol1 /= sol1 .OR. sol2 /= sol2) THEN
        km = MAX(km1, km2)
      ELSE IF ((km1 - sol1) * (km2 - sol1) <= 0._wp) THEN
        km = sol1
      ELSE
        km = sol2
      ENDIF
    ENDIF
    ! consider soil water content
    IF (PRESENT(w_soil)) THEN
      km = km * w_soil
    END IF
  END FUNCTION calc_km_qmax_correction

  ! ======================================================================================================= !
  !>calculates bulk soil carbon used for the calculation of soil bulk density and sorption capacity
  !>
  !>
  !>  Input: soil pools
  !>
  !>  Output: bulk_soil_carbon_sl, soil_litter_carbon_sl
  !>
  SUBROUTINE calc_bulk_soil_carbon( &
    & nc, &
    & nsoil_w, &
    & num_sl_above_bedrock, &
    & sb_model_scheme, &
    & sb_pool_mt, &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & volume_min_sl)

    USE mo_jsb_physical_constants, ONLY: molar_mass_C
    USE mo_sb_constants,           ONLY: rho_bulk_org, carbon_per_dryweight_SOM
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,          INTENT(in)    :: nc                         !<
    INTEGER,          INTENT(in)    :: nsoil_w                    !<
    REAL(wp),         INTENT(in)    :: num_sl_above_bedrock(:)    !<
    CHARACTER(len=*), INTENT(in)    :: sb_model_scheme            !< soil model scheme: simple_1d or jsm
    REAL(wp),         INTENT(in)    :: sb_pool_mt(:,:,:,:)        !< pools [mol/m3]
    REAL(wp),         INTENT(inout) :: bulk_soil_carbon_sl(:,:)   !< mol/m3
    REAL(wp),         INTENT(inout) :: soil_litter_carbon_sl(:,:) !< mol/m3
    REAL(wp),         INTENT(inout) :: volume_min_sl(:,:)         !< m3/m3
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                             :: isoil, ic
    REAL(wp), DIMENSION(nc, nsoil_w)                    :: om_dens_sl   !< organic matter density [kg / m3]
    REAL(wp), DIMENSION(nc, nsoil_w)                    :: lit_dens_sl  !< litter density [kg / m3]
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_soil_carbon'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate the content of organic matter [mol/m3] and litter [mol/m3]
    !>
    bulk_soil_carbon_sl(:,:) = sb_pool_mt(ix_microbial, ixC, :, :)       &
      &                      + sb_pool_mt(ix_residue, ixC, :, :)         &
      &                      + sb_pool_mt(ix_residue_assoc, ixC, :, :)   &
      &                      + sb_pool_mt(ix_dom, ixC, :, :)             &
      &                      + sb_pool_mt(ix_dom_assoc, ixC, :, :)

    ! Note: woody litter in the top-layer is assumed to be intact trunks of wood on the ground. They do not
    !   mix with the mineral soil and therefore a correction term for the soil bulk density is not needed. ...
    soil_litter_carbon_sl(:,1) = sb_pool_mt(ix_soluable_litter, ixC, :, 1)   &
      &                        + sb_pool_mt(ix_polymeric_litter, ixC, :, 1)

    ! ... Woody litter in the lower soil layers are from dying course roots. They disintegrate rather quickly
    ! and therefore mix with the mineral soil, requiring a correction of the soil bulk density term (mixing
    ! OM and mineral soils.)
    DO ic = 1,nc
      DO isoil = 2, NINT(num_sl_above_bedrock(ic))
        soil_litter_carbon_sl(ic,isoil) = sb_pool_mt(ix_woody_litter, ixC, ic, isoil)    &
          &                            + sb_pool_mt(ix_soluable_litter, ixC, ic, isoil)  &
          &                            + sb_pool_mt(ix_polymeric_litter, ixC, ic, isoil)
      END DO
    END DO

    !>  1.1 Only jsm distinguish litter C from SOM C
    !>
    IF (sb_model_scheme == 'simple_1d') THEN
      bulk_soil_carbon_sl(:,:)   = bulk_soil_carbon_sl(:,:) + soil_litter_carbon_sl(:,:)
      soil_litter_carbon_sl(:,:) = 0._wp
    ENDIF

    !>2.0 calculate the volume of mineral soil [m3/m3]
    !>  assuming the bulk soil is composed of OM, litter and soil mineral
    !>
    om_dens_sl  = bulk_soil_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp    ! OM density [kg/m3]
    lit_dens_sl = soil_litter_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp  ! litter density [kg/m3]
    DO ic = 1,nc
      DO isoil = 1, NINT(num_sl_above_bedrock(ic))
        volume_min_sl(ic,isoil) = MAX(0.01_wp, &
          &   MIN( 0.99_wp, 1._wp - om_dens_sl(ic,isoil) / rho_bulk_org - lit_dens_sl(ic,isoil) / rho_bulk_org))
      END DO
    END DO
  END SUBROUTINE calc_bulk_soil_carbon

  ! ======================================================================================================= !
  !>calculates the fast_assoc_po4 based on Langmuir equilibrium and Labile P pool
  !>
  !>  Input: Labile P pool, Qmax, Km, all in the unit [mol m-3]
  !>
  !>  Output: fast_assoc_po4
  !>
  FUNCTION calc_fast_po4(Labile_p, Qmax, km) RESULT (fast_po4)

    REAL(wp), INTENT(in) :: Labile_p                    !< labile P pool of the soil layer [mol P / m3]
    REAL(wp), INTENT(in) :: Qmax                        !< maximium sorption capacity of po4 [mol P/ m3]
    REAL(wp), INTENT(in) :: km                          !< sorption affinity of po4 [mol P/ m3]
    REAL(wp)             :: fast_po4                    !< fast associated po4 of the soil layer [mol P / m3]
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fast_po4'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate fast_assoc_po4
    !>
    ! based on the assumption that labile P is the sum of both,
    ! and they follow the Langmuir equilibrium
    ! Labile_p - fast_po4 = fast_po4 * km / (Labile_p * (Qmax - fast_po4))
    fast_po4 = (Labile_p + Qmax + km - SQRT((Labile_p + Qmax + km) ** 2._wp - 4._wp * Labile_p * Qmax)) / 2._wp
  END FUNCTION calc_fast_po4

  ! ======================================================================================================= !
  !>calculates potential enzyme allocation coefficient
  !>  based on SESAM concept, by Thomas Wutzler
  !>
  !> Input: potential depolymerisation rate (C, N, P), microbial biomass, half-saturation coefficient
  !>
  !> Output: potential enzyme coefficient for C, N, P
  !>
  ELEMENTAL FUNCTION calc_enzyme_allocation(dL, dR, km, Cmic) RESULT (alpha)

    REAL(wp), INTENT(in) :: dL          !< depolymerization from litter [mol / m3]
    REAL(wp), INTENT(in) :: dR          !< depolymerization from residue [mol / m3]
    REAL(wp), INTENT(in) :: km          !< half-saturation concentration of depolymerization [mol / m3]
    REAL(wp), INTENT(in) :: Cmic        !< microbial C biomass [mol C/ m3]
    REAL(wp)             :: alpha       !< partitioning coefficient for enzyme allocation [--]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: hlp1    !< b**2-4*a*c in the quadratic equation
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_enzyme_allocation'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculates potential enzyme allocation coefficient
    !>
    ! Solution of the equation:
    !   dL/dR=alpha(km+alpha*Cmic)/((1-alpha)*(km+Cmic*(1-alpha)))
    !   only one of the solution fits in the range [0,1]
    hlp1  = 4._wp * Cmic ** 2._wp * dL * dR &
      &     + 8._wp * Cmic * km * dL * dR &
      &     + km ** 2._wp * dL ** 2._wp &
      &     + 2._wp * km ** 2._wp * dL * dR &
      &     + km ** 2._wp * dR ** 2._wp
    alpha = (dL * km + 2._wp * dL * Cmic + dR * km - SQRT(hlp1)) / (2._wp * Cmic * (dL - dR))
  END FUNCTION calc_enzyme_allocation

  ! ======================================================================================================= !
  !>calculates sorption and desorption of organic material
  !>  corresponds to the sorption and desorption terms in Ahrens et al. 2015, eq. 5 & 6
  !>
  !> Input: sorption temperature and moisture modifiers, sorption pools
  !>
  !> Output: sorption and desorption rate (mol/m2/timestep)
  !>
  SUBROUTINE calc_sorption_desorption_of_org_material( &
    & nc, &
    & nsoil_w, &
    & dtime, &
    & elements_index_map, &
    & is_element_used, &
    & qmax_org, &
    & rtm_sorption, &
    & rmm_sorption, &
    & rtm_desorption, &
    & rmm_desorption, &
    & sb_pool_mt, &
    & sb_formation_mt, &
    & sb_loss_mt )

    USE mo_sb_constants,          ONLY: k_adsorpt_dom, k_desorpt_dom, k_adsorpt_det, k_desorpt_det
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc                       !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_w                  !< number of soil layers
    REAL(wp),                           INTENT(in)    :: dtime                    !< timestep length
    INTEGER,                            INTENT(in)    :: elements_index_map(:)    !< map bgcm element ID -> IDX
    LOGICAL,                            INTENT(in)    :: is_element_used(:)       !< is element in 'elements_index_map' used
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: qmax_org                 !< ...
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rtm_sorption             !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rmm_sorption             !< moisture rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rtm_desorption           !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_w),   INTENT(in)    :: rmm_desorption           !< moisture rate modifier [-]
    REAL(wp),                           INTENT(in)    :: sb_pool_mt(:,:,:,:)      !< bgcm sb_pool: [mol/m3]
    REAL(wp),                           INTENT(inout) :: sb_formation_mt(:,:,:,:) !< bgcm flux: transport [mol/m3/timestep]
    REAL(wp),                           INTENT(inout) :: sb_loss_mt(:,:,:,:)      !< bgcm flux: transport [mol/m3/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                               :: ielem                !< loop over bgcm elements
    INTEGER                               :: ix_elem              !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc, nsoil_w)      :: hlp1
    REAL(wp), DIMENSION(nc, nsoil_w)      :: hlp2
    REAL(wp), DIMENSION(nc, nsoil_w)      :: avail_soil_surface_sorption_capacity_sl
    REAL(wp), DIMENSION(nc, nsoil_w)      :: sorption_rate
    REAL(wp), DIMENSION(nc, nsoil_w)      :: desorption_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sorption_desorption_of_org_material'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 Determine available sorption sites [mol C m-3]
    !>
    avail_soil_surface_sorption_capacity_sl(:,:) = MAX(0._wp, qmax_org(:,:) &
      &                                            - (sb_pool_mt(ix_dom_assoc, ixC, :, :) &
      &                                            + sb_pool_mt(ix_residue_assoc, ixC, :, :)))

    !>2.0 mineral association of DOM
    !>
    sorption_rate(:,:)   = k_adsorpt_dom * rtm_sorption(:,:)   * rmm_sorption(:,:) &
      &                    * avail_soil_surface_sorption_capacity_sl(:,:) * dtime
    desorption_rate(:,:) = k_desorpt_dom * rtm_desorption(:,:) * rmm_desorption(:,:) * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        hlp1(:,:)                                    = sorption_rate(:,:) * sb_pool_mt(ix_dom, ix_elem, :, :)
        sb_formation_mt(ix_dom_assoc, ix_elem, :, :) = sb_formation_mt(ix_dom_assoc, ix_elem, :, :) + hlp1(:,:)
        sb_loss_mt(ix_dom, ix_elem, :, :)            = sb_loss_mt(ix_dom, ix_elem, :, :) + hlp1(:,:)
        hlp2(:,:)                                    = desorption_rate(:,:) * sb_pool_mt(ix_dom_assoc, ix_elem, :, :)
        sb_formation_mt(ix_dom, ix_elem, :, :)       = sb_formation_mt(ix_dom, ix_elem, :, :) + hlp2(:,:)
        sb_loss_mt(ix_dom_assoc, ix_elem, :, :)      = sb_loss_mt(ix_dom_assoc, ix_elem, :, :) + hlp2(:,:)
      END IF
    END DO

    !>3.0 mineral association of microbial residue
    !>
    sorption_rate(:,:)   = k_adsorpt_det * rtm_sorption(:,:)   * rmm_sorption(:,:) &
      &                    * avail_soil_surface_sorption_capacity_sl(:,:) * dtime
    desorption_rate(:,:) = k_desorpt_det * rtm_desorption(:,:) * rmm_desorption(:,:) * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        hlp1(:,:)                                        = sorption_rate(:,:) * sb_pool_mt(ix_residue, ix_elem, :, :)
        sb_formation_mt(ix_residue_assoc, ix_elem, :, :) = sb_formation_mt(ix_residue_assoc, ix_elem, :, :) + hlp1(:,:)
        sb_loss_mt(ix_residue, ix_elem, :, :)            = sb_loss_mt(ix_residue, ix_elem, :, :) + hlp1(:,:)
        hlp2(:,:)                                        = desorption_rate(:,:) * sb_pool_mt(ix_residue_assoc, ix_elem, :, :)
        sb_formation_mt(ix_residue, ix_elem, :, :)       = sb_formation_mt(ix_residue, ix_elem, :, :) + hlp2(:,:)
        sb_loss_mt(ix_residue_assoc, ix_elem, :, :)      = sb_loss_mt(ix_residue_assoc, ix_elem, :, :) + hlp2(:,:)
      END IF
    END DO
  END SUBROUTINE calc_sorption_desorption_of_org_material

  ! ======================================================================================================= !
  !>calculate P sorption parameters, including both single and double Langmuir
  !>
  SUBROUTINE calc_Psorption_parameter( &
    & flag_sb_double_langmuir, &
    & clay_sl, silt_sl, soil_water, ph_sl, &
    & bulk_soil_carbon_sl, soil_litter_carbon_sl, woody_litter_carbon_sl, &
    & bulk_dens_sl, volume_min_sl, po4_solute, &
    & qmax_po4_min_sl, qmax_po4_om_sl, qmax_po4, km_adsorpt_po4_sl, &             ! inout
    & qmax_fast_po4, qmax_slow_po4, km_fast_po4, km_slow_po4, Qmax_AlFe_cor)      ! inout

    USE mo_sb_constants,   ONLY: qmax_po4_silt, qmax_po4_sand, qmax_po4_clay, qmax_po4_OM, qmax_po4_mineral, &
      & km_adsorpt_silt_po4, km_adsorpt_clay_po4, km_adsorpt_sand_po4, km_adsorpt_OM_po4, &
      & km_adsorpt_mineral_po4, km_po4_ph, min_fraction_default_sb
    USE mo_spq_util,       ONLY: calc_qmax_texture
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL,      INTENT(in)    :: flag_sb_double_langmuir      !< conventional or double Langmuir
    REAL(wp),     INTENT(in)    :: clay_sl                      !< clay content of soil layer [-]
    REAL(wp),     INTENT(in)    :: silt_sl                      !< silt content of soil layer [-]
    REAL(wp),     INTENT(in)    :: soil_water                   !< soil water content of soil layer [L/m3]
    REAL(wp),     INTENT(in)    :: ph_sl                        !< ph of soil layer [-]
    REAL(wp),     INTENT(in)    :: bulk_soil_carbon_sl          !< bulk soil carbon content of soil layer [mol C/m3]
    REAL(wp),     INTENT(in)    :: soil_litter_carbon_sl        !< soil litter carbon content of soil layer [mol C/m3]
    REAL(wp),     INTENT(in)    :: woody_litter_carbon_sl       !< woody litter carbon content of soil layer [mol C/m3]
    REAL(wp),     INTENT(in)    :: bulk_dens_sl                 !< soil bulk density [kg/m3]
    REAL(wp),     INTENT(in)    :: volume_min_sl                !< soil mineral volumetric content [m3/m3]
    REAL(wp),     INTENT(in)    :: po4_solute                   !< soluble P concentration of soil layer [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_po4_min_sl              !< maximium po4 sorption capacity for mineral soil [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_po4_om_sl               !< maximium po4 sorption capacity for OM soil [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_po4                     !< maximium po4 sorption capacity for soil layer [mol P/m3]
    REAL(wp),     INTENT(inout) :: km_adsorpt_po4_sl            !< half-saturation po4 concentration [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_fast_po4                !< maximium po4 sorption capacity of fast associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_slow_po4                !< maximium po4 sorption capacity of slow associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: km_fast_po4                  !< half-saturation po4 concentration of fast associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: km_slow_po4                  !< half-saturation po4 concentration of slow associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: Qmax_AlFe_cor                !< qmax correlation factor of Al/Fe content [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)  :: hlp1, hlp2, hlp3, hlp4
    REAL(wp)  :: sand_sl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_Psorption_parameter'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.8 initialize local variable
    !>
    sand_sl = 1._wp - clay_sl - silt_sl

    !>0.9 set Qmax_AlFe_cor zero if lower eps8
    !>
    IF (Qmax_AlFe_cor < eps8) THEN
      Qmax_AlFe_cor = 0._wp
    END IF

    !>1.0 double-surface Langmuir
    !>
    IF (flag_sb_double_langmuir) THEN
      !>  1.1 calculate Qmax based on soil texture
      !>
      ! fine soil + sandy soil [mol P /kg soil]
      qmax_po4_min_sl = (calc_qmax_texture(qmax_po4_silt, clay_sl, silt_sl) + &
        & calc_qmax_texture(qmax_po4_sand, 0._wp, sand_sl)) * Qmax_AlFe_cor
      ! soil OM [mol P /kg OM]
      qmax_po4_om_sl = qmax_po4_OM
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl, &
        & volume_min_sl, &
        & bulk_dens_sl, &
        & qmax_po4_om_sl, &
        & qmax_po4_min_sl, &
        & qmax_po4, &
        & qmax_fast=qmax_fast_po4, &
        & qmax_slow=qmax_slow_po4)
      ! Fe/Al in OM [mol P/m3]
      hlp1 = calc_qmax_texture(qmax_po4_clay, 0._wp, clay_sl) * volume_min_sl * bulk_dens_sl * Qmax_AlFe_cor
      ! update qmax_fast_po4 and qmax_po4
      qmax_fast_po4 = qmax_fast_po4 + hlp1
      qmax_po4      = qmax_po4 + hlp1

      !>  1.2 calculate km based on soil texture
      !>
      ! km for the OM-Fe/Al pool, fast pool
      hlp2 = qmax_fast_po4 - hlp1                        ! qmax for non-adsorbed Fe/Al [mol P/m3]
      hlp3 = km_po4_ph * ph_sl * km_adsorpt_clay_po4     ! km for Fe/Al, mol P/m3
      hlp4 = km_po4_ph * ph_sl * km_adsorpt_OM_po4       ! km for OM, mol P/m3
      km_fast_po4 = calc_km_qmax_correction(po4_solute, hlp1, hlp2, hlp3, hlp4, w_soil=soil_water)
      ! km for the clay-silt-sand pool, slow pool
      IF (qmax_slow_po4 < eps8) THEN
        qmax_slow_po4 = 0._wp
        hlp3        = km_po4_ph * ph_sl * km_adsorpt_silt_po4    ! km for clay/silt [mol P/L]
        hlp4        = km_po4_ph * ph_sl * km_adsorpt_sand_po4    ! km for sand [mol P/L]
        km_slow_po4 = MAX(hlp3, hlp4)
      ELSE
        hlp1        = qmax_slow_po4 * MIN(calc_qmax_texture(qmax_po4_silt, clay_sl, silt_sl) * &
          & Qmax_AlFe_cor / qmax_po4_min_sl, (1._wp - min_fraction_default_sb)) ! qmax for clay/silt
        hlp2        = qmax_slow_po4 - hlp1                                      ! qmax for sand
        hlp3        = km_po4_ph * ph_sl * km_adsorpt_silt_po4                   ! km for clay/silt [mol P/L]
        hlp4        = km_po4_ph * ph_sl * km_adsorpt_sand_po4                   ! km for sand [mol P/L]
        km_slow_po4 = calc_km_qmax_correction(po4_solute, hlp1, hlp2, hlp3, hlp4, w_soil=soil_water)
      END IF
      ! km for the whole soil layer
      km_adsorpt_po4_sl = calc_km_qmax_correction(po4_solute, qmax_fast_po4,qmax_slow_po4, km_fast_po4, km_slow_po4)

    !>2.0 single-surface Langmuir
    !>
    ! Qmax and km are based on texture and calculated both as the sum of OM and mineral soil
    ELSE
      qmax_po4_min_sl   = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
      qmax_po4_om_sl    = qmax_po4_min_sl * min_fraction_default_sb
      CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, volume_min_sl, bulk_dens_sl, &
        &                                    qmax_po4_om_sl, qmax_po4_min_sl, &
        &                                    qmax_po4, &
        &                                    qmax_fast=qmax_fast_po4, &
        &                                    qmax_slow=qmax_slow_po4)
      qmax_po4          = qmax_po4 * soil_water / 1000._wp ! corrected by water content
      CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, volume_min_sl, bulk_dens_sl, &
        &                                    km_adsorpt_OM_po4, km_adsorpt_mineral_po4, &
        &                                    km_adsorpt_po4_sl, &
        &                                    qmax_fast=km_fast_po4, &
        &                                    qmax_slow=km_slow_po4)
    END IF  ! double/single Langmuir
  END SUBROUTINE calc_Psorption_parameter

#endif
END MODULE mo_q_sb_jsm_processes
