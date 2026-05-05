!> QUINCY agriculture process memory
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
!>#### definition and init of (memory) variables for the quincy agriculture process
!>
MODULE mo_q_agr_memory_class
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

  PUBLIC :: t_q_agr_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  ! ======================================================================================================= !
  !> Type definition for q_agr memory
  !>
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_agr_memory

    ! fertiliser application rate
    TYPE(t_jsb_var_real2d) ::   &
      & crop_growth_phase, &             !< growth phase of crop [0 = Planting, 1 = Emergence, 2 = GrainFill,3 = Harvest]
      & gdd_mavg, &                 !< long-term average growing degree days [degC day]
      & nd_crop_season, &           !< number of days in this crop growing season [# days]
      & nd_crop_season_mavg, &      !< long-term average crop growing season length [# days]
      & n_fertiliser_c3, &          !< N fertiliser application rate for C3 crops [mol m-2 yr-1]
      & n_fertiliser_c4, &          !< N fertiliser application rate for C4 crops [mol m-2 yr-1]
      & crop_season_per_year, &      !< number of crop seasons (agricultural growing seasons) in the current year
      & crop_season_per_year_mavg, & !< long-term average number of crop seasons in a year
      & crop_type_index             !< crop type index (1-ncroptypes)

  CONTAINS
    PROCEDURE :: Init => Init_q_agr_memory
  END TYPE t_q_agr_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: q_agr
  !>
  !>
  SUBROUTINE Init_q_agr_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_jsb_class,           ONLY: Get_model
    USE mo_jsb_varlist,         ONLY: BASIC, MEDIUM, FULL, NONE
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_quincy_output_class, ONLY: unitless
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_agr_memory),    INTENT(inout), TARGET :: mem
    CHARACTER(len=*),        INTENT(in)            :: prefix
    CHARACTER(len=*),        INTENT(in)            :: suffix
    INTEGER,                 INTENT(in)            :: lct_ids(:)
    INTEGER,                 INTENT(in)            :: lib_id          !< id of primary lct in lctlib
    INTEGER,                 INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),  POINTER :: model                        !< model
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    INTEGER                     :: table                        ! ...
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_agr_memory'
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
      ! memory to only create on the box tile
      CALL mem%Add_var('n_fertiliser_c3', mem%n_fertiliser_c3, &
        & hgrid, surface, &
        & t_cf('n_fertiliser_c3', 'mol m-2 yr-1', 'N fertiliser application rate for C3 crops'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & lrestart = .TRUE., &
        & output_level = NONE, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_fertiliser_c4', mem%n_fertiliser_c4, &
        & hgrid, surface, &
        & t_cf('n_fertiliser_c4', 'mol m-2 yr-1', 'N fertiliser application rate for C4 crops'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & lrestart = .TRUE., &
        & output_level = NONE, &
        & initval_r = 0.0_wp)
    END IF

    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN
      ! memory to create on tiles of or containing LAND_TYPE or VEG_TYPE
      CALL mem%Add_var('crop_growth_phase', mem%crop_growth_phase, &
        & hgrid, surface, &
        & t_cf('crop_growth_phase', '[0 - 3]', 'growth phase of crop from planting to harvest'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = NONE, &
        & lrestart = .TRUE., &
        & initval_r = 5.0_wp)

      CALL mem%Add_var('gdd_mavg', mem%gdd_mavg, &
        & hgrid, surface, &
        & t_cf('gdd_mavg', 'degC days', 'long-term average growing degree days'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = NONE, &
        & lrestart = .TRUE., &
        & initval_r = 1700.0_wp)

      CALL mem%Add_var('nd_crop_season', mem%nd_crop_season, &
        & hgrid, surface, &
        & t_cf('nd_crop_season', '# days', 'number of days in this crop growing season'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = NONE, &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nd_crop_season_mavg', mem%nd_crop_season_mavg, &
        & hgrid, surface, &
        & t_cf('nd_crop_season_mavg', '# days', 'long-term average length of crop growing season'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = NONE, &
        & lrestart = .TRUE., &
        & initval_r = 150.0_wp)

      CALL mem%Add_var('crop_season_per_year', mem%crop_season_per_year, &
        & hgrid, surface, &
        & t_cf('crop_season_per_year', '# crop seasons yr-1', 'number of crop seasons in this year'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & lrestart = .TRUE., &
        & output_level = NONE, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('crop_season_per_year_mavg', mem%crop_season_per_year_mavg, &
        & hgrid, surface, &
        & t_cf('crop_season_per_year_mavg', '# crop seasons yr-1', 'long-term average of number of crop seasons per year'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & lrestart = .TRUE., &
        & output_level = NONE, &
        & initval_r = 1.0_wp)

      CALL mem%Add_var('crop_type_index', mem%crop_type_index, &
        & hgrid, surface, &
        & t_cf('crop_type_index', '[1-8]', 'index of croptype'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & lrestart = .TRUE., &
        & output_level = NONE, &
        & initval_r = 0.0_wp)

    END IF

  END SUBROUTINE Init_q_agr_memory

#endif
END MODULE mo_q_agr_memory_class
