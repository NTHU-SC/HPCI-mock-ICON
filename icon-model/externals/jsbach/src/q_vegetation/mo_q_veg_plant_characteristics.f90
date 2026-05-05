!> QUINCY update plant characteristics
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
!>#### routine for updating the characteristics of the plants on a tile after changes in the vegetation pools
!>
MODULE mo_q_veg_plant_characteristics
#ifndef __NO_QUINCY__

  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, finish

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_plant_characteristics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_plant_characteristics'

CONTAINS

  ! ======================================================================================================= !
  !>update vegetation pools
  !>
  SUBROUTINE update_plant_characteristics(tile, options)
    USE mo_kind,                            ONLY: wp
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_quincy_model_config,             ONLY: QLAND, QPLANT, QCANOPY
    USE mo_jsb_process_class,               ONLY: VEG_, RAD_, Q_ASSIMI_, Q_PHENO_, TURB_
    USE mo_jsb_grid_class,                  ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                        ONLY: Get_vgrid
    USE mo_jsb_physical_constants,          ONLY: rhoh2o, grav
    USE mo_jsb_math_constants,              ONLY: one_day, eps8
    USE mo_q_veg_canopy,                    ONLY: calc_canopy_layers
    USE mo_q_veg_growth,                    ONLY: calc_diameter_from_woody_biomass, calc_height_from_diameter
    USE mo_veg_constants,                   ONLY: itree, k_fpc, k_sai2lai
    USE mo_turb_constants,                  ONLY: veg_height_to_rough_momentum, veg_min_roughness
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(Q_ASSIMI_)
    dsl4jsb_Use_config(TURB_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_canopy_q_assimi  !< Vertical grid
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_soil_w           !< Vertical grid
    INTEGER                           :: ncanopy                !< number of canopy layers
    INTEGER                           :: nsoil_w                !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(options%nc)   :: leaf_nitrogen          !< helper array for leaf nitrogen
    REAL(wp), ALLOCATABLE             :: vgrid_canopy_dz(:)     !< canopy grid helper variable dz
    REAL(wp), ALLOCATABLE             :: vgrid_canopy_lbounds(:)!< canopy grid helper variable lbounds
    REAL(wp), ALLOCATABLE             :: vgrid_canopy_ubounds(:)!< canopy grid helper variable ubounds
    REAL(wp)                          :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    REAL(wp)                          :: dtime                  !< timestep length
    INTEGER                           :: iblk, ics, ice, nc     !< grid dimensions
    INTEGER                           :: ic, is, icanopy        !< looping indices
    CHARACTER(len=1024)               :: canopy_cond_scheme     !< canopy_conductance_scheme: medlyn / ballberry - q_assimi config
    LOGICAL                           :: flag_optimal_Nfraction !< on/off optimise leaf internal N allocation - q_assimi config
    REAL(wp) :: config_blending_height
    REAL(wp) :: lctlib_sla, lctlib_np_leaf, lctlib_allom_k1, lctlib_allom_k2
    REAL(wp) :: lctlib_wood_density, lctlib_phi_leaf_min, lctlib_k_latosa, lctlib_k0_fn_struc, lctlib_fn_oth_min
    REAL(wp) :: lctlib_gmin, lctlib_g0, lctlib_t_jmax_omega
    INTEGER  :: lctlib_growthform, lctlib_ps_pathway
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_plant_characteristics'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: wtr_soil_root_pot
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: lai_max
    ! RAD_ 3D
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: sai
    dsl4jsb_Real2D_onChunk      :: fract_fpc
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: rough_veg_star
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    dsl4jsb_Real3D_onChunk      :: delta_root_fraction_sl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
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
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    vgrid_soil_w          => Get_vgrid('soil_depth_water')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    nsoil_w               =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (dsl4jsb_Lctlib_param(BareSoilFlag)) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_config(TURB_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_tfrac_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soa_tphen_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_tfrac_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tfrac_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, wtr_soil_root_pot)         ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, lai_max)                    ! in
    ! RAD_ 3D
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! VEG 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                             ! inout (only calculated for trees)
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                             ! inout (only calculated for trees)
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                               ! inout (only calculated for trees)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, sai)                                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, fract_fpc)                            ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)                                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)                        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, rough_veg_star)                       ! out
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)               ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                      ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                       ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                     ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                      ! inout (in case of optimality, else out)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)                      ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)               ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, delta_root_fraction_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                         ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)                    ! out
    ! ----------------------------------------------------------------------------------------------------- !

    ALLOCATE(vgrid_canopy_dz(ncanopy), vgrid_canopy_lbounds(ncanopy), vgrid_canopy_ubounds(ncanopy))
    !$ACC ENTER DATA CREATE(vgrid_canopy_dz(:), vgrid_canopy_lbounds(:), vgrid_canopy_ubounds(:))

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO icanopy = 1, ncanopy
      vgrid_canopy_dz(icanopy) = vgrid_canopy_q_assimi%dz(icanopy)
      vgrid_canopy_lbounds(icanopy) = vgrid_canopy_q_assimi%lbounds(icanopy)
      vgrid_canopy_ubounds(icanopy) = vgrid_canopy_q_assimi%ubounds(icanopy)
    END DO
    !$ACC END PARALLEL LOOP

    config_blending_height = dsl4jsb_Config(TURB_)%blending_height

    lctlib_sla = dsl4jsb_Lctlib_param(sla)
    lctlib_np_leaf = dsl4jsb_Lctlib_param(np_leaf)
    lctlib_growthform = dsl4jsb_Lctlib_param(growthform)
    lctlib_allom_k1 = dsl4jsb_Lctlib_param(allom_k1)
    lctlib_allom_k2 = dsl4jsb_Lctlib_param(allom_k2)
    lctlib_wood_density = dsl4jsb_Lctlib_param(wood_density)
    lctlib_phi_leaf_min = dsl4jsb_Lctlib_param(phi_leaf_min)
    lctlib_k_latosa = dsl4jsb_Lctlib_param(k_latosa)
    lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)

    lctlib_k0_fn_struc = dsl4jsb_Lctlib_param(k0_fn_struc)
    lctlib_fn_oth_min = dsl4jsb_Lctlib_param(fn_oth_min)
    lctlib_gmin = dsl4jsb_Lctlib_param(gmin)
    lctlib_g0 = dsl4jsb_Lctlib_param(g0)
    lctlib_t_jmax_omega = dsl4jsb_Lctlib_param(t_jmax_omega)

    canopy_cond_scheme = TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme)
    flag_optimal_Nfraction = dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction

    !>0.9 set g1 according to canopy_conductance_scheme
    !>
    SELECT CASE(canopy_cond_scheme)
      CASE ("medlyn")
        lctlib_g1 = dsl4jsb_Lctlib_param(g1_medlyn)
      CASE ("ballberry")
        lctlib_g1 = dsl4jsb_Lctlib_param(g1_bberry)
    END SELECT

    !>1.0 Update plant diagnostics given the above fluxes
    !>

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      lai(ic)      = veg_pool_mt(ix_leaf, ixC, ic) * lctlib_sla

      IF (lctlib_growthform == itree) THEN
        dens_ind(ic) = MAX(0.0_wp, dens_ind(ic) + delta_dens_ind(ic))
        ! reset delta_dens_ind to zero in cases plant characteristics are updated several times in the timestep
        delta_dens_ind(ic) = 0.0_wp

        ! In CANOPY mode dens_ind is zero and height is a pft specific lctlib constant
        IF (dens_ind(ic) > eps8) THEN
          diameter(ic)  = calc_diameter_from_woody_biomass(lctlib_allom_k1                     , &
            &                                              lctlib_allom_k2                     , &
            &                                              lctlib_wood_density                 , &
            &                                              veg_pool_mt(ix_sap_wood, ixC, ic)   , &
            &                                              veg_pool_mt(ix_heart_wood, ixC, ic) , &
            &                                              dens_ind(ic))
          height(ic)    = calc_height_from_diameter(lctlib_allom_k1, lctlib_allom_k2, diameter(ic))
        END IF
      END IF

      dphi(ic) = wtr_soil_root_pot(ic) - lctlib_phi_leaf_min - grav * rhoh2o * height(ic) * 1.e-6_wp
    END DO
    !$ACC END PARALLEL LOOP

    IF (lctlib_growthform == itree) THEN
      SELECT CASE(model%config%qmodel_id)
      CASE(QPLANT, QLAND)
        DO ic = 1, nc
          IF (height(ic) > eps8) THEN
            sai(ic) = k_sai2lai * veg_pool_mt(ix_sap_wood, ixC, ic) &
              &       * lctlib_k_latosa / lctlib_wood_density / height(ic)
          ELSE
            sai(ic) = 0.0_wp
          END IF
        END DO
      CASE(QCANOPY)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
        DO ic = 1, nc
          sai(ic) = k_sai2lai * lai_max(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END SELECT
    END IF

    !>  1.1 calculate foliage projected cover fraction and rough_veg_star for turbulence calculations
    !>
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      fract_fpc(ic) = 1.0_wp - EXP(-k_fpc * (lai(ic) + sai(ic)))
      ! NOTE use as minimum value for all PFT veg_min_roughness, because
      !      vegetation density and height can change and
      !      therefore a tile occupied by a woody PFT can have a height lower than 10m, implying a rough_m less than 1;
      !      jsbach is using the lctlib parameter MinVegRoughness instead
      rough_veg_star(ic) = 1._wp / LOG(config_blending_height / MAX((height(ic) / veg_height_to_rough_momentum), &
        &                  veg_min_roughness)) ** 2
    END DO
    !$ACC END PARALLEL LOOP

    !>  1.2 implied change of root fraction given root growth
    !>
    IF (dsl4jsb_Config(VEG_)%flag_dynamic_roots) THEN
#ifdef _OPENACC
      CALL finish(routine, 'Code block ported to GPU but not tested, yet. Stop.')
#endif

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT)
      DO is = 1, nsoil_w
        DO ic = 1, nc
          root_fraction_sl(ic,is) = root_fraction_sl(ic,is) + delta_root_fraction_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !> 2.0 Given current vegetation pools, update canopy layers
    !>
    !$ACC DATA CREATE(leaf_nitrogen)
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      leaf_nitrogen(ic) = veg_pool_mt(ix_leaf, ixN, ic)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    CALL calc_canopy_layers( nc                                                       , & ! in
                             ncanopy                                                  , &
                             dtime                                                    , &
                             vgrid_canopy_dz(:)                                       , &
                             vgrid_canopy_lbounds(:)                                  , &
                             vgrid_canopy_ubounds(:)                                  , &
                             lctlib_ps_pathway                                        , &
                             lctlib_k0_fn_struc                                       , &
                             lctlib_fn_oth_min                                        , &
                             lctlib_sla                                               , &
                             lctlib_np_leaf                                           , &
                             lctlib_gmin                                              , &
                             lctlib_g0                                                , &
                             lctlib_g1                                                , &
                             lctlib_t_jmax_omega                                      , &
                             flag_optimal_Nfraction                                   , & ! Q_ASSIMI_ config
                             canopy_cond_scheme                                       , & ! Q_ASSIMI_ config (medlyn/ballberry)
                             leaf_nitrogen(:)                                         , & ! in
                             lai(:)                                                   , &
                             ppfd_sunlit_tfrac_mavg_cl(:,:)                           , & ! inout
                             ppfd_shaded_tfrac_mavg_cl(:,:)                           , &
                             fleaf_sunlit_tfrac_mavg_cl(:,:)                          , &
                             fn_rub_cl(:,:)                                           , & ! out
                             fn_et_cl(:,:)                                            , &
                             fn_pepc_cl(:,:)                                          , &
                             fn_chl_cl(:,:)                                           , & ! inout (in case of optimality, else out)
                             fn_oth_cl(:,:)                                           , & ! out
                             lai_cl(:,:)                                              , & ! out
                             cumm_lai_cl(:,:)                                         , & ! out
                             leaf_nitrogen_cl(:,:)                                    , & ! out
                             t_air            = t_air_tfrac_mavg(:), &                    ! optional in
                             t_acclim         = t_air_tacclim_mavg(:), &
                             press_srf        = press_srf_tfrac_mavg(:), &
                             co2_mixing_ratio = co2_mixing_ratio_tfrac_mavg(:), &
                             aerodyn_cond     = ga_tfrac_mavg(:), &
                             beta_air         = beta_air_tfrac_mavg(:), &
                             beta_soa         = beta_soa_tphen_mavg(:), &
                             beta_soil_ps     = beta_soil_ps_tfrac_mavg(:), &
                             beta_sinklim_ps  = beta_sinklim_ps_tfrac_mavg(:), &
                             beta_soil_gs     = beta_soil_gs_tfrac_mavg(:), &
                             t_jmax_opt       = t_jmax_opt_mavg(:) )                      ! optional in
    !$ACC WAIT(1)
    !$ACC END DATA
    !$ACC EXIT DATA DELETE(vgrid_canopy_dz, vgrid_canopy_lbounds, vgrid_canopy_ubounds)
    DEALLOCATE(vgrid_canopy_dz, vgrid_canopy_lbounds, vgrid_canopy_ubounds)

  END SUBROUTINE update_plant_characteristics

#endif
END MODULE mo_q_veg_plant_characteristics
