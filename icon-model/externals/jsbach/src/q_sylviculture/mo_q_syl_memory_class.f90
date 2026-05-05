!> QUINCY sylviculture process memory
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
!>#### definition and init of (memory) variables for the quincy sylviculture process
!>
MODULE mo_q_syl_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  dsl4jsb_Use_processes VEG_
  dsl4jsb_Use_config(VEG_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_syl_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  ! ======================================================================================================= !
  !> Type definition for q_syl memory
  !>
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_syl_memory

    TYPE(t_jsb_var_real2d) ::        &
      & fract_forest_harvest_y_read, & !< Grid-cell fraction of forest area to be harvested - only available on box memory (unitless)
      & fract_forest_harvest,        & !< Grid-cell fraction harvested from this forest tile (unitless)
      & fract_wood_to_slash            !< Fraction of harvested wood allocated to slash (unitless)

  CONTAINS
    PROCEDURE :: Init => Init_q_syl_memory
  END TYPE t_q_syl_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_syl_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: q_syl
  !>
  !>
  SUBROUTINE Init_q_syl_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_jsb_class,           ONLY: Get_model
    USE mo_jsb_varlist,         ONLY: BASIC, MEDIUM, FULL, NONE
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_quincy_output_class, ONLY: unitless
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_syl_memory),    INTENT(inout), TARGET :: mem
    CHARACTER(len=*),        INTENT(in)            :: prefix
    CHARACTER(len=*),        INTENT(in)            :: suffix
    INTEGER,                 INTENT(in)            :: lct_ids(:)
    INTEGER,                 INTENT(in)            :: lib_id     !< id of primary lct in lctlib
    INTEGER,                 INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),  POINTER :: model                      !< model
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    INTEGER                     :: table                        ! ...
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_syl_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(model_id)
    table                 = tables(1)
    hgrid                 => Get_grid(mem%grid_id)
    surface               => Get_vgrid('surface')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)

    ! created memory depends on tile type
    IF ( .NOT. ASSOCIATED(mem%parent)) THEN
      ! ... memory to only create on the box tile
      CALL mem%Add_var('fract_forest_harvest_y_read', mem%fract_forest_harvest_y_read, &
        & hgrid, surface, &
        & t_cf('fract_forest_harvest_y_read', unitless, 'Grid-cell forest area to be harvested in this year, read from file.'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = NONE, &
        & lrestart = .TRUE., &
        & lrestart_cont = .TRUE., &
        & initval_r = 0.0_wp)

      IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        CALL mem%Add_var('fract_wood_to_slash', mem%fract_wood_to_slash, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_slash', unitless, 'Fraction of harvested wood allocated to slash'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)
      END IF ! IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
    END IF

    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN
      ! ... memory to create on tiles of or containing LAND_TYPE or VEG_TYPE
      CALL mem%Add_var('fract_forest_harvest', mem%fract_forest_harvest, &
        & hgrid, surface, &
        & t_cf('fract_forest_harvest', unitless, 'Grid-cell area harvested from this forest tile'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & lrestart = .TRUE., &
        & lrestart_cont = .TRUE., &
        & initval_r = 0.0_wp)
    END IF
  END SUBROUTINE Init_q_syl_memory

#endif
END MODULE mo_q_syl_memory_class
