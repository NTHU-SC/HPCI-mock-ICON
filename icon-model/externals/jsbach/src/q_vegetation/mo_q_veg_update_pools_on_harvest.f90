!> QUINCY update veg and product pools on harvest
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
!>#### update veg and product pools and call litter formation on harvest
!>
MODULE mo_q_veg_update_pools_on_harvest
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message
  USE mo_jsb_time,            ONLY: is_newyear, is_newday

  USE mo_jsb_math_constants,      ONLY: eps8, eps4, one_day
  USE mo_jsb_process_class,       ONLY: VEG_, Q_SYL_, Q_AGR_, SB_, HYDRO_
  USE mo_q_sb_litter_processes,   ONLY: calc_litter_partitioning

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, VEG_BGCM_PP_FIBERBOARD_ID, &
    &       VEG_BGCM_PP_CROP_ID, VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID, VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_HARVEST_TO_PROD_ID, &
    &       VEG_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, VEG_BGCM_HARVEST_LITTER_ID, VEG_BGCM_SEED_BED_POOL_ID


  dsl4jsb_use_memory(HYDRO_)
  dsl4jsb_use_memory(Q_SYL_)
  dsl4jsb_use_memory(VEG_)
  dsl4jsb_use_memory(Q_AGR_)

  dsl4jsb_use_config(VEG_)
  dsl4jsb_use_config(SB_)
  dsl4jsb_use_config(Q_SYL_)

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_pools_on_harvest

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_update_pools_on_harvest'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Update quincy pools upon harvest
  !> Note: within this routine the elements within the compartments of the veg pool of the given tile are
  !>       already reduced proportionally to the harvested fractions -- i.e. OUTSIDE of update veg pools!
  !>       The formation of litter, in contrast, remains part of sb update pools.
  !>
  SUBROUTINE update_pools_on_harvest(tile, options)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_veg_config_class,      ONLY: get_number_of_veg_compartments
    USE mo_q_agr_process,         ONLY: calc_crop_harvest_fraction
    USE mo_veg_constants,         ONLY: min_lai
    USE mo_q_syl_constants,       ONLY: fract_wood_to_pp_fuel, fract_wood_to_pp_paper, fract_wood_to_pp_fiberboard,   &
      &                                 fract_wood_to_pp_oirw, fract_wood_to_pp_pv, fract_wood_to_pp_sawnwood
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: box_tile                   !< pointer to the box tile to get the harvest slash fraction
    TYPE(t_jsb_model),          POINTER :: model                      !< the model
    TYPE(t_lnd_bgcm_store),     POINTER :: bgcm_store                 !< the bgcm store of this tile
    TYPE(t_lctlib_element),     POINTER :: lctlib                     !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_soil_w               !< Vertical grid
    INTEGER                             :: nsoil_w                    !< number of soil layers (water)
    INTEGER                             :: elem_idx_map(LAST_ELEM_ID) !< element mapper ID -> IX
    LOGICAL                             :: is_elem_used(LAST_ELEM_ID) !< indicates which elements are used in this simulation
    REAL(wp)                            :: zero_mt(LAST_ELEM_ID, options%nc) !< helper array to pass zero seedbed litter
    REAL(wp), DIMENSION(options%nc)     :: fract_harvest_rel_to_tile  !< harvested area fraction relative to the tile area
    REAL(wp), DIMENSION(options%nc)     :: cover_fraction             !< current tile area
    REAL(wp)                            :: dtime                      !< timestep length
    INTEGER                             :: ix_ct                      !< index of crop type (if crop)
    INTEGER                             :: ic, iblk, ics, ice, nc     !< dimensions and loop counter
    INTEGER                             :: i_compartment, nr_of_veg_bgcm_comp !< loop counter and dim for veg compartments
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':update_pools_on_harvest'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(Q_SYL_)

    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onChunk :: crop_type_index
    dsl4jsb_Real2D_onChunk :: crop_growth_phase

    dsl4jsb_Def_memory(Q_SYL_)
    dsl4jsb_Real2D_onChunk :: fract_forest_harvest

    dsl4jsb_Def_memory_tile(Q_SYL_, box_tile)
    dsl4jsb_Real2D_onChunk :: fract_wood_to_slash

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: soil_depth_sl

    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: dens_ind
    dsl4jsb_Real2D_onChunk :: delta_dens_ind
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L2D :: seed_bed_pool_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_crop_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_fuel_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_paper_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_fiberboard_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_oirw_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_pv_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_sawnwood_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_flux_harvest_mt
    dsl4jsb_Def_mt2L2D :: veg_litter_flux_harvest_mt
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! currently this only needs to be called for a tile if either Q_SYL_ or Q_AGR_ is running on that tile
    ! later, with alcc, this will be different
    IF (.NOT. (tile%Is_process_calculated(Q_SYL_) .OR. tile%Is_process_calculated(Q_AGR_))) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    nr_of_veg_bgcm_comp = get_number_of_veg_compartments()
    CALL model%Get_top_tile(box_tile)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_SYL_)
    ! ----------------------------------------------------------------------------------------------------- !
    zero_mt(:,:) = 0.0_wp

    IF (tile%Is_process_calculated(Q_SYL_)) THEN
      IF (dsl4jsb_Config(Q_SYL_)%l_daily_harvest) THEN
        ! In case of daily harvest we have to check if this is a new day
        IF (.NOT. is_newday(options%current_datetime, dtime)) THEN
          ! If not then we do not have a harvest event and can return
          RETURN
        END IF
      ELSE
        ! In case of annual harvest we have to check if this is a new year
        IF (.NOT. is_newyear(options%current_datetime, dtime)) THEN
          ! If not then we do not have a harvest event and can return
          RETURN
        END IF
      END IF
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)

    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)

    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)

    IF (tile%Is_process_calculated(Q_AGR_)) THEN
      dsl4jsb_Get_memory(Q_AGR_)
      dsl4jsb_Get_var2D_onChunk(Q_AGR_, crop_type_index)
      dsl4jsb_Get_var2D_onChunk(Q_AGR_, crop_growth_phase)

    ELSE IF (tile%Is_process_calculated(Q_SYL_)) THEN
      dsl4jsb_Get_memory(Q_SYL_)
      dsl4jsb_Get_var2D_onChunk(Q_SYL_, fract_forest_harvest)

      IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        dsl4jsb_Get_memory_tile(Q_SYL_, box_tile)
        dsl4jsb_Get_var2D_onChunk_tile(Q_SYL_, fract_wood_to_slash, box_tile)
      END IF
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_HARVEST_LITTER_ID, veg_litter_flux_harvest_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)

    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      IF (tile%Is_process_calculated(Q_SYL_) .OR. tile%Is_process_calculated(Q_AGR_)) THEN
        dsl4jsb_Get_mt1L2D(VEG_BGCM_HARVEST_TO_PROD_ID, veg_pp_flux_harvest_mt)

        IF (tile%Is_process_calculated(Q_AGR_)) THEN
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_CROP_ID, veg_pp_crop_mt)
          dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_POOL_ID, seed_bed_pool_mt)

        ELSE IF (tile%Is_process_calculated(Q_SYL_)) THEN
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt)
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PAPER_ID, veg_pp_paper_mt)
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FIBERBOARD_ID, veg_pp_fiberboard_mt)
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_OIRW_ID, veg_pp_oirw_mt)
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PV_ID, veg_pp_pv_mt)
          dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_SAWNWOOD_ID, veg_pp_sawnwood_mt)
        END IF
      END IF
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    elem_idx_map(:) = model%config%elements_index_map(:)
    is_elem_used(:) = model%config%is_element_used(:)


    ! Harvest fraction depends on the kind of harvest - currently Q_AGR_ or Q_SYL_
    fract_harvest_rel_to_tile(:) = 0.0_wp
    IF (tile%Is_process_calculated(Q_AGR_)) THEN
      ! calculate flux from harvesting if crop has reached maturity (i.e. growth phase is 3)
      fract_harvest_rel_to_tile(:) = calc_crop_harvest_fraction(dtime, crop_type_index(:),   &
                                                                crop_growth_phase(:), lai(:))

    ELSE IF (tile%Is_process_calculated(Q_SYL_)) THEN
      ! ----------------------------------------------------------------------------------------------------- !
      ! Derive the harvested fraction in proportion to the tile area given the absolute harvested fraction
      CALL tile%Get_fraction(ics, ice, iblk, fract=cover_fraction)
      DO ic = 1,nc
        IF (cover_fraction(ic) > eps8) THEN
          fract_harvest_rel_to_tile(ic) = MIN(1.0_wp, fract_forest_harvest(ic) / cover_fraction(ic))
        ELSE
          fract_harvest_rel_to_tile(ic) = 0.0_wp
        END IF
      END DO
    ELSE
      ! Currently only the processes Q_SYL_ and Q_AGR_ lead to harvest
      ! fract_harvest_rel_to_tile(:) = 0.0_wp
    END IF

    IF (ANY(fract_harvest_rel_to_tile(:) /= 0.0_wp)) THEN
      ! ----------------------------------------------------------------------------------------------------- !
      !>
      !> - Note: later we might care for establishment here to represent re-planting
      !>
      ! veg_establishment_mt
      ! Note: in mo_q_agr_interface this is approached - suggestion: moved this here

      ! ----------------------------------------------------------------------------------------------------- !
      ! - Transfer harvested amount of each element of each compartment to harvest litter flux (but leave the seeds)
      DO ic = 1,nc
        DO i_compartment = 1,nr_of_veg_bgcm_comp
          veg_litter_flux_harvest_mt(i_compartment,:,ic) = veg_pool_mt(i_compartment,:,ic) * fract_harvest_rel_to_tile(ic)
          veg_pool_mt(i_compartment,:,ic) = veg_pool_mt(i_compartment,:,ic) - veg_litter_flux_harvest_mt(i_compartment,:,ic)
        END DO
      END DO

      ! ----------------------------------------------------------------------------------------------------- !
      !>
      !> - Adapt the change of individual density to account for the harvested fraction (only for tiles with forest type)
      !>
      IF (dsl4jsb_Lctlib_param(ForestFlag)) THEN
        DO ic = 1,nc
          delta_dens_ind(ic) = - (fract_harvest_rel_to_tile(ic) * dens_ind(ic))
        END DO
      END IF

      ! ----------------------------------------------------------------------------------------------------- !
      !>
      !> - If product pools are used, the litter flux is decreased by the fluxes to the product pools
      !>
      IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        ! leaf, fine_root, fruit, coarse_root, labile, reserve and seeds (if any): 100% to litter

        IF (tile%Is_process_calculated(Q_AGR_)) THEN
          CALL calc_crop_harvest_to_products_flux(nc, &                                 ! in
            &                                     veg_litter_flux_harvest_mt(:,:,:), & ! inout
            &                                     veg_pool_mt(:,:,:), &
            &                                     seed_bed_pool_mt(:,:), &
            &                                     veg_pp_crop_mt(:,:), &
            &                                     veg_pp_flux_harvest_mt(:,:))

        ELSE IF (tile%Is_process_calculated(Q_SYL_)) THEN
          DO ic = 1,nc
            ! sapwood: lctlib_frac_sapwood_branch goes to litter, rest to product pool -> but mind: slash fraction
            veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
              & + (veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
              &    * (1._wp - lctlib%frac_sapwood_branch) * (1._wp - fract_wood_to_slash(ic)))
            veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) = veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
              & - (veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
              &    * (1._wp - lctlib%frac_sapwood_branch) * (1._wp - fract_wood_to_slash(ic)))

            ! heartwood: goes to product pool -> but mind: slash fraction
            veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
              & + veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) * (1._wp - fract_wood_to_slash(ic))
            veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) = veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) &
              & - veg_litter_flux_harvest_mt(ix_heart_wood,:,ic) * (1._wp - fract_wood_to_slash(ic))
          END DO

          ! distribute to the different product pools
          CALL add_flux_to_product(nc, is_elem_used, elem_idx_map, &
            &                      fract_wood_to_pp_fuel, veg_pp_flux_harvest_mt, veg_pp_fuel_mt)
          CALL add_flux_to_product(nc, is_elem_used, elem_idx_map, &
            &                      fract_wood_to_pp_paper, veg_pp_flux_harvest_mt, veg_pp_paper_mt)
          CALL add_flux_to_product(nc, is_elem_used, elem_idx_map, &
            &                      fract_wood_to_pp_fiberboard, veg_pp_flux_harvest_mt, veg_pp_fiberboard_mt)
          CALL add_flux_to_product(nc, is_elem_used, elem_idx_map, &
            &                      fract_wood_to_pp_oirw, veg_pp_flux_harvest_mt, veg_pp_oirw_mt)
          CALL add_flux_to_product(nc, is_elem_used, elem_idx_map, &
            &                      fract_wood_to_pp_pv, veg_pp_flux_harvest_mt, veg_pp_pv_mt)
          CALL add_flux_to_product(nc, is_elem_used, elem_idx_map, &
            &                      fract_wood_to_pp_sawnwood, veg_pp_flux_harvest_mt, veg_pp_sawnwood_mt)
        END IF ! IF (tile%Is_process_calculated(Q_SYL_))
      END IF ! IF dsl4jsb_Config(VEG_)%l_use_product_pools

      ! ----------------------------------------------------------------------------------------------------- !
      !>
      !> - After potential reduction due to product usage, the litter flux is put into the sb formation flux
      !>
      CALL calc_litter_partitioning( &
        & nc, &                                         ! in
        & nsoil_w, &
        & num_sl_above_bedrock(:), &
        & lctlib%sla, &
        & lctlib%growthform, &
        & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
        & soil_depth_sl(:,:), &
        & root_fraction_sl(:,:), &
        & veg_litter_flux_harvest_mt(:,:,:), &          ! in
        & zero_mt(:,:), &
        & sb_formation_mt(:,:,:,:) )                    ! inout

    END IF

  END SUBROUTINE update_pools_on_harvest


  ! ======================================================================================================= !
  !>
  !> add flux to product pool
  !>
  !> simply takes the given fraction that should go into the given product pool and moves it there...
  !>
  SUBROUTINE add_flux_to_product(nc, is_elem_used, elem_idx_map, fract_wood_to_pp, veg_pp_flux_harvest_mt, veg_pp_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,      INTENT(in) :: nc                          !< block dimension
    LOGICAL,      INTENT(in) :: is_elem_used(:)             !< indicates which elements are in use
    INTEGER,      INTENT(in) :: elem_idx_map(:)             !< map bgcm element ID -> IX
    REAL(wp),     INTENT(in) :: fract_wood_to_pp            !< wood fraction to go into this product pool
    REAL(wp),     INTENT(in) :: veg_pp_flux_harvest_mt(:,:) !< harvest flux to product pools
    REAL(wp),  INTENT(inout) :: veg_pp_mt(:,:)              !< this product pool
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic, id_elem              !< loop counter
    INTEGER                     :: ix_elem                  !< element index
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_flux_to_product'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc
      DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (is_elem_used(id_elem)) THEN
          ix_elem = elem_idx_map(id_elem)
          veg_pp_mt(ix_elem,ic) = veg_pp_mt(ix_elem,ic) + veg_pp_flux_harvest_mt(ix_elem,ic) * fract_wood_to_pp
        END IF
      END DO
    END DO

  END SUBROUTINE add_flux_to_product

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate partitioning of harvesting flux into litter and products
  !!  This follows the Q_AGR calculation of a litter flux
  !!  leaves and stems (minus a slash fraction) as well as fruits are removed, the remainder is
  !!  returned to litter. In case the seed bed is too low to allow for planting next growing season
  !!  a fraction of the fruit harvested is retained and transferred to seed bed
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_crop_harvest_to_products_flux( &
    & nc, &
    & veg_litter_flux_harvest_mt, &
    & veg_pool_mt, &
    & seed_bed_pool_mt, &
    & veg_pp_crop_mt, &
    & veg_pp_flux_harvest_mt)

    USE mo_lnd_bgcm_idx
    USE mo_q_agr_constants,         ONLY: crop_planting_mass, fstore_seed_min, fract_crop_to_slash
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                                !< dimensions
    REAL(wp), INTENT(inout) :: veg_litter_flux_harvest_mt(:,:,:) !< bgcm flux from harvest to litter
    REAL(wp), INTENT(inout) :: veg_pool_mt(:,:,:)                !< bgcm vegetation pool
    REAL(wp), INTENT(inout) :: seed_bed_pool_mt(:,:)             !< bgcm seed bed pool
    REAL(wp), INTENT(inout) :: veg_pp_crop_mt(:,:)               !< bgcm crop products pool
    REAL(wp), INTENT(inout) :: veg_pp_flux_harvest_mt(:,:)       !< bgcm harvest flux to products
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                !< loop over grid cells
    REAL(wp)       :: hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_harvest_to_products_flux'
    ! ----------------------------------------------------------------------------------------------------- !

    ! fine_root, coarse_root, labile, reserve and seeds (if any): 100% to litter
    DO ic = 1,nc

      ! leaves: slash fraction goes to litter, rest to product pool
      veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
        & + veg_litter_flux_harvest_mt(ix_leaf,:,ic) * (1._wp - fract_crop_to_slash)
      veg_litter_flux_harvest_mt(ix_leaf,:,ic) = veg_litter_flux_harvest_mt(ix_leaf,:,ic) &
        & - veg_litter_flux_harvest_mt(ix_leaf,:,ic) * (1._wp - fract_crop_to_slash)

      ! stems: slash fraction goes to litter, rest to product pool
      veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
        & + veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) * (1._wp - fract_crop_to_slash)
      veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) = veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) &
        & - veg_litter_flux_harvest_mt(ix_sap_wood,:,ic) * (1._wp - fract_crop_to_slash)

      ! fruits: first ensure that the seed pool does not decrease below twice planting level to have planting
      !         material for next growing season
      IF (seed_bed_pool_mt(ixC,ic) < fstore_seed_min * crop_planting_mass &
          & .AND. veg_litter_flux_harvest_mt(ix_fruit,ixC,ic) > eps8) THEN
        hlp1 = MIN(MAX( &
          &             (fstore_seed_min * crop_planting_mass - seed_bed_pool_mt(ixC,ic)) &
          &               / veg_litter_flux_harvest_mt(ix_fruit,ixC,ic),0.0_wp),1.0_wp)
        seed_bed_pool_mt(:,ic) = seed_bed_pool_mt(:,ic) &
          & + hlp1 * veg_litter_flux_harvest_mt(ix_fruit,:,ic)
        veg_litter_flux_harvest_mt(ix_fruit,:,ic) = veg_litter_flux_harvest_mt(ix_fruit,:,ic) &
          & - hlp1 * veg_litter_flux_harvest_mt(ix_fruit,:,ic)
      END IF

      veg_pp_flux_harvest_mt(:,ic) = veg_pp_flux_harvest_mt(:,ic) &
        & + veg_litter_flux_harvest_mt(ix_fruit,:,ic)
      veg_litter_flux_harvest_mt(ix_fruit,:,ic) = veg_litter_flux_harvest_mt(ix_fruit,:,ic) &
        & - veg_litter_flux_harvest_mt(ix_fruit,:,ic)

      veg_pp_crop_mt(:,ic) = veg_pp_crop_mt(:,ic) + veg_pp_flux_harvest_mt(:,ic)

    END DO

  END SUBROUTINE calc_crop_harvest_to_products_flux

#endif
END MODULE mo_q_veg_update_pools_on_harvest
