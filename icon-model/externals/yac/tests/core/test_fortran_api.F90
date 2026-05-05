! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

!> \file test_fortran_api.F90
!! \test
!! This example tests the usage of the Fortran interface of the
!! YAC core.

#define NOP(x) associate( x => x ); end associate

PROGRAM main

  USE, INTRINSIC :: iso_c_binding
  USE utest
  USE yac_core
  USE mpi

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE c_free ( ptr ) BIND ( c, name='free' )

      USE, INTRINSIC :: iso_c_binding, only : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE c_free

  END INTERFACE

  INTEGER :: ierror

  REAL(kind=c_double), PARAMETER :: &
    YAC_RAD = 0.017453292519943295769_c_double ! M_PI / 180

  TYPE(c_ptr) :: basic_grid_empty
  TYPE(c_ptr) :: basic_grid_reg_2d
  TYPE(c_ptr) :: basic_grid_reg_2d_deg
  TYPE(c_ptr) :: basic_grid_curve_2d
  TYPE(c_ptr) :: basic_grid_curve_2d_deg
  TYPE(c_ptr) :: basic_grid_unstruct
  TYPE(c_ptr) :: basic_grid_unstruct_deg
  TYPE(c_ptr) :: basic_grid_unstruct_ll
  TYPE(c_ptr) :: basic_grid_unstruct_ll_deg
  TYPE(c_ptr) :: basic_grid_unstruct_edge
  TYPE(c_ptr) :: basic_grid_unstruct_edge_deg
  TYPE(c_ptr) :: basic_grid_unstruct_edge_ll
  TYPE(c_ptr) :: basic_grid_unstruct_edge_ll_deg
  TYPE(c_ptr) :: basic_grid_cloud
  TYPE(c_ptr) :: basic_grid_cloud_deg
  TYPE(c_ptr) :: basic_grid_reg_2d_rot
  TYPE(c_ptr) :: basic_grid_reg_2d_rot_deg
  INTEGER(kind=c_size_t) :: num_cells
  INTEGER(kind=c_size_t) :: num_corners
  INTEGER(kind=c_size_t) :: num_edges
  REAL(kind=c_double), ALLOCATABLE :: cell_areas(:)

  REAL(kind=c_double) :: &
    cell_coords(3) = (/1.0_c_double, 0.0_c_double, 0.0_c_double/)
  INTEGER(kind=c_int) :: cell_mask(1) = (/1_c_int/)
  INTEGER(kind=c_size_t) :: cell_coord_idx, cell_coord_ll_idx
  INTEGER(kind=c_size_t) :: cell_mask_idx, cell_mask_ll_idx

  TYPE(c_ptr) :: grid_pair
  TYPE(c_ptr) :: interp_grid
  TYPE(c_ptr) :: interp_stack_config
  TYPE(c_ptr) :: interp_stack_config_cpy

  ! data types required to generate extended spmap configuration
  TYPE(c_ptr) :: scale_config
  TYPE(c_ptr) :: spmap_config
  TYPE(c_ptr) :: src_point_selection
  TYPE(c_ptr) :: source_cell_area_config, target_cell_area_config
  TYPE(c_ptr) :: overwrite_spmap_config
  TYPE(c_ptr) :: overwrite_config
  TYPE(c_ptr), TARGET :: overwrite_configs(2)

  INTEGER(kind=c_int) :: compare_value, test_value
  LOGICAL :: ltest_value

  INTEGER, TARGET :: do_search_call_count
  INTEGER, TARGET :: constructor_call_count
  TYPE(c_ptr) :: interp_method
  TYPE(c_ptr) :: interp_weights
  TYPE(c_ptr) :: interpolation
  TYPE(c_ptr) :: interpolation_frac

  REAL(kind=c_double), TARGET :: src_cell_data(1)
  TYPE(c_ptr), TARGET :: src_field_(1)
  TYPE(c_ptr), TARGET :: src_field_collection(1)
  REAL(kind=c_double), TARGET :: src_frac_mask(1)
  TYPE(c_ptr), TARGET :: src_frac_mask_(1)
  TYPE(c_ptr), TARGET :: src_frac_mask_collection(1)
  REAL(kind=c_double), TARGET :: tgt_cell_data(1)
  TYPE(c_ptr), TARGET :: tgt_field_collection(1)

  INTERFACE
    SUBROUTINE C_UNLINK ( path ) BIND ( c, name='unlink' )
      USE, INTRINSIC :: iso_c_binding, only : c_char
      CHARACTER(KIND=c_char), DIMENSION(*) :: path
    END SUBROUTINE C_UNLINK
  END INTERFACE

  ! ===================================================================

  CALL start_test('fortran_api')

  CALL yac_mpi_init_c()
  CALL yac_yaxt_init_c(MPI_COMM_WORLD)

  CALL test(yac_mpi_is_initialised_c() /= 0_c_int)

  ! empty grids
  basic_grid_empty = &
    yac_basic_grid_empty_new_c("empty_grid" // c_null_char)
  num_cells = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_empty, INT(YAC_LOC_CELL, c_int))
  CALL test(num_cells == 0_c_int)

  ! regular 2d grid (rad coords)
  basic_grid_reg_2d = &
    yac_basic_grid_reg_2d_new_c( &
      "reg_2d_grid" // c_null_char, &
      (/2_c_size_t, 2_c_size_t/), (/0_c_int, 0_c_int/), &
      (/-0.1_c_double, 0.1_c_double/), (/-0.1_c_double, 0.1_c_double/))
  num_corners = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_reg_2d, INT(YAC_LOC_CORNER, c_int))
  CALL test(num_corners == 4_c_int)

  ALLOCATE( &
    cell_areas( &
      INT( &
        yac_basic_grid_get_data_size_c( &
          basic_grid_reg_2d, INT(YAC_LOC_CELL, c_int)))))
  cell_areas(:) = -1.0_c_double
  CALL yac_basic_grid_compute_cell_areas_c(basic_grid_reg_2d, cell_areas)
  CALL test(ALL(cell_areas(:) /= -1.0_c_double))
  DEALLOCATE(cell_areas)

  ! regular 2d grid (deg coords)
  basic_grid_reg_2d_deg = &
    yac_basic_grid_reg_2d_deg_new_c( &
      "reg_2d_grid_deg" // c_null_char, &
      (/2_c_size_t, 2_c_size_t/), (/0_c_int, 0_c_int/), &
      (/-0.1_c_double, 0.1_c_double/), (/-0.1_c_double, 0.1_c_double/))
  num_corners = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_reg_2d_deg, INT(YAC_LOC_CORNER, c_int))
  CALL test(num_corners == 4_c_int)

  ! curvilinear 2d grid (rad coords)
  basic_grid_curve_2d = &
    yac_basic_grid_curve_2d_new_c( &
      "curve_2d_grid" // c_null_char, &
      (/2_c_size_t, 2_c_size_t/), (/0_c_int, 0_c_int/), &
      (/-0.1_c_double, 0.1_c_double, -0.1_c_double, 0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_curve_2d, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! curvilinear 2d grid (deg coords)
  basic_grid_curve_2d_deg = &
    yac_basic_grid_curve_2d_deg_new_c( &
      "curve_2d_grid_deg" // c_null_char, &
      (/2_c_size_t, 2_c_size_t/), (/0_c_int, 0_c_int/), &
      (/-0.1_c_double, 0.1_c_double, -0.1_c_double, 0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_curve_2d_deg, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (rad coords)
  basic_grid_unstruct = &
    yac_basic_grid_unstruct_new_c( &
      "unstruct_grid" // c_null_char, &
      4_c_size_t, 1_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (deg coords)
  basic_grid_unstruct_deg = &
    yac_basic_grid_unstruct_deg_new_c( &
      "unstruct_grid_deg" // c_null_char, &
      4_c_size_t, 1_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_deg, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (rad coords)
  basic_grid_unstruct_ll = &
    yac_basic_grid_unstruct_ll_new_c( &
      "unstruct_grid_ll" // c_null_char, &
      4_c_size_t, 1_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_ll, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (deg coords)
  basic_grid_unstruct_ll_deg = &
    yac_basic_grid_unstruct_ll_deg_new_c( &
      "unstruct_grid_ll_deg" // c_null_char, &
      4_c_size_t, 1_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_ll_deg, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (rad coords)
  basic_grid_unstruct_edge = &
    yac_basic_grid_unstruct_edge_new_c( &
      "unstruct_edge_grid" // c_null_char, &
      4_c_size_t, 1_c_size_t, 4_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/), &
      (/0_c_int, 1_c_int, 1_c_int, 2_c_int, &
        2_c_int, 3_c_int, 3_c_int, 0_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_edge, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (deg coords)
  basic_grid_unstruct_edge_deg = &
    yac_basic_grid_unstruct_edge_deg_new_c( &
      "unstruct_grid_edge_deg" // c_null_char, &
      4_c_size_t, 1_c_size_t, 4_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/), &
      (/0_c_int, 1_c_int, 1_c_int, 2_c_int, &
        2_c_int, 3_c_int, 3_c_int, 0_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_edge_deg, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (rad coords)
  basic_grid_unstruct_edge_ll = &
    yac_basic_grid_unstruct_edge_ll_new_c( &
      "unstruct_grid_edge_ll" // c_null_char, &
      4_c_size_t, 1_c_size_t, 4_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/), &
      (/0_c_int, 1_c_int, 1_c_int, 2_c_int, &
        2_c_int, 3_c_int, 3_c_int, 0_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_edge_ll, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! unstructured grid (deg coords)
  basic_grid_unstruct_edge_ll_deg = &
    yac_basic_grid_unstruct_edge_ll_deg_new_c( &
      "unstruct_grid_edge_ll_deg" // c_null_char, &
      4_c_size_t, 1_c_size_t, 4_c_size_t, (/4_c_int/), &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/), &
      (/0_c_int, 1_c_int, 2_c_int, 3_c_int/), &
      (/0_c_int, 1_c_int, 1_c_int, 2_c_int, &
        2_c_int, 3_c_int, 3_c_int, 0_c_int/))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_unstruct_edge_ll_deg, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_edges == 4_c_int)

  ! cloud grid (rad coords)
  basic_grid_cloud = &
    yac_basic_grid_cloud_new_c( &
      "cloud_grid" // c_null_char, 4_c_size_t, &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/))
  num_cells = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_cloud, INT(YAC_LOC_CELL, c_int))
  num_corners = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_cloud, INT(YAC_LOC_CORNER, c_int))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_cloud, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_cells == 0_c_int)
  CALL test(num_corners == 4_c_int)
  CALL test(num_edges == 0_c_int)

  ! cloud grid (deg coords)
  basic_grid_cloud_deg = &
    yac_basic_grid_cloud_deg_new_c( &
      "cloud_grid" // c_null_char, 4_c_size_t, &
      (/-0.1_c_double, 0.1_c_double, 0.1_c_double, -0.1_c_double/), &
      (/-0.1_c_double, -0.1_c_double, 0.1_c_double, 0.1_c_double/))
  num_cells = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_cloud_deg, INT(YAC_LOC_CELL, c_int))
  num_corners = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_cloud_deg, INT(YAC_LOC_CORNER, c_int))
  num_edges = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_cloud_deg, INT(YAC_LOC_EDGE, c_int))
  CALL test(num_cells == 0_c_int)
  CALL test(num_corners == 4_c_int)
  CALL test(num_edges == 0_c_int)

  ! rotated regular 2d grid (rad coords)
  basic_grid_reg_2d_rot = &
    yac_basic_grid_reg_2d_rot_new_c( &
      "reg_2d_grid_deg" // c_null_char, &
      (/2_c_size_t, 2_c_size_t/), (/0_c_int, 0_c_int/), &
      (/-0.5_c_double, 0.5_c_double/) * YAC_RAD, &
      (/-0.5_c_double, 0.5_c_double/) * YAC_RAD, &
      180.0_c_double * YAC_RAD, 45.0_c_double * YAC_RAD)
  num_corners = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_reg_2d_rot, INT(YAC_LOC_CORNER, c_int))
  CALL test(num_corners == 4_c_int)

  ! rotated regular 2d grid (deg coords)
  basic_grid_reg_2d_rot_deg = &
    yac_basic_grid_reg_2d_rot_deg_new_c( &
      "reg_2d_grid_deg" // c_null_char, &
      (/2_c_size_t, 2_c_size_t/), (/0_c_int, 0_c_int/), &
      (/-0.5_c_double, 0.5_c_double/), &
      (/-0.5_c_double, 0.5_c_double/), &
      180.0_c_double, 45.0_c_double)
  num_corners = &
    yac_basic_grid_get_data_size_c( &
      basic_grid_reg_2d_rot_deg, INT(YAC_LOC_CORNER, c_int))
  CALL test(num_corners == 4_c_int)

  ! adding field coordinates and masks
  cell_coord_idx = &
    yac_basic_grid_add_coordinates_c( &
      basic_grid_unstruct, YAC_LOC_CELL, cell_coords, 1_c_size_t)
  cell_coord_ll_idx = &
    yac_basic_grid_add_coordinates_c( &
      basic_grid_unstruct_ll, YAC_LOC_CELL, cell_coords, 1_c_size_t)
  cell_mask_idx = &
    yac_basic_grid_add_mask_c( &
      basic_grid_unstruct, YAC_LOC_CELL, cell_mask, 1_c_size_t, &
      "cell_mask" // c_null_char)
  cell_mask_ll_idx = &
    yac_basic_grid_add_mask_c( &
      basic_grid_unstruct_ll, YAC_LOC_CELL, cell_mask, 1_c_size_t, &
      "cell_mask_ll" // c_null_char)

  ! distributed grid
  grid_pair = &
    yac_dist_grid_pair_new_c( &
      basic_grid_unstruct, basic_grid_unstruct_ll, MPI_COMM_WORLD)

  ! interpolation grid
  interp_grid = &
    yac_interp_grid_new_c( &
      grid_pair, "unstruct_grid" // c_null_char, &
      "unstruct_grid_ll" // c_null_char, &
      1_c_size_t, (/YAC_LOC_CELL/), (/cell_coord_idx/), (/cell_mask_idx/), &
      YAC_LOC_CELL, cell_coord_ll_idx, cell_mask_ll_idx)

  ! interpolation check
  CALL yac_interp_method_check_add_constructor_callback_c( &
    c_funloc(constructor_callback), c_loc(constructor_call_count), &
    "constructor_callback" // c_null_char)
  CALL yac_interp_method_check_add_do_search_callback_c( &
    c_funloc(do_search_callback), c_loc(do_search_call_count), &
    "do_search_callback" // c_null_char)

  ! setup structures for extended spmap configuration
  source_cell_area_config = &
    yac_spmap_cell_area_config_yac_new_c( &
      YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT_F)
  target_cell_area_config = &
    yac_spmap_cell_area_config_file_new_c( &
      "area.nc", "cell_areas", YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT_F)
  scale_config = &
    yac_spmap_scale_config_new_c( &
      YAC_INTERP_SPMAP_SCALE_TYPE_DEFAULT_F, &
      source_cell_area_config, target_cell_area_config)
  src_point_selection = &
    yac_point_selection_bnd_circle_new_c( &
      0.0_c_double, 0.0_c_double, 0.1_c_double)
  spmap_config = &
    yac_interp_spmap_config_new_c( &
      YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT_F, &
      YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT_F, &
      YAC_INTERP_SPMAP_WEIGHTED_DEFAULT_F, &
      scale_config)
  overwrite_spmap_config = &
    yac_interp_spmap_config_new_c( &
      0.1_c_double, &
      YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT_F, &
      YAC_INTERP_SPMAP_WEIGHTED_DEFAULT_F, &
      scale_config)
  overwrite_config = &
    yac_spmap_overwrite_config_new_c( &
      src_point_selection, overwrite_spmap_config)
  overwrite_configs = (/overwrite_config, c_null_ptr/)

  ! interpolation stack configuration
  constructor_call_count = 0
  do_search_call_count = 0
  interp_stack_config = yac_interp_stack_config_new_c()
  CALL yac_interp_stack_config_add_check_c( &
    interp_stack_config, "constructor_callback" // c_null_char, &
    "do_search_callback" // c_null_char)
  CALL yac_interp_stack_config_add_nnn_c( &
    interp_stack_config, YAC_INTERP_NNN_WEIGHTED_DEFAULT_F, &
    YAC_INTERP_NNN_N_DEFAULT_F, YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT_F, &
    YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT_F)
  interp_stack_config_cpy = &
    yac_interp_stack_config_copy_c(interp_stack_config)
  CALL yac_interp_stack_config_add_average_c( &
    interp_stack_config, YAC_INTERP_AVG_WEIGHT_TYPE_DEFAULT_F, &
    YAC_INTERP_AVG_PARTIAL_COVERAGE_DEFAULT_F)
  CALL yac_interp_stack_config_add_ncc_c( &
    interp_stack_config, YAC_INTERP_NCC_WEIGHT_TYPE_DEFAULT_F, &
    YAC_INTERP_NCC_PARTIAL_COVERAGE_DEFAULT_F)
  CALL yac_interp_stack_config_add_nnn_c( &
    interp_stack_config, YAC_INTERP_NNN_WEIGHTED_DEFAULT_F, &
    YAC_INTERP_NNN_N_DEFAULT_F, YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT_F, &
    YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT_F)
  CALL yac_interp_stack_config_add_rbf_c( &
    interp_stack_config, YAC_INTERP_RBF_N_DEFAULT_F, &
    YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT_F, &
    YAC_INTERP_RBF_SCALE_DEFAULT_F)
  CALL yac_interp_stack_config_add_conservative_c( &
    interp_stack_config, YAC_INTERP_CONSERV_ORDER_DEFAULT_F, &
    YAC_INTERP_CONSERV_ENFORCED_CONSERV_DEFAULT_F, &
    YAC_INTERP_CONSERV_PARTIAL_COVERAGE_DEFAULT_F, &
    YAC_INTERP_CONSERV_NORMALISATION_DEFAULT_F)
  CALL yac_interp_stack_config_add_spmap_c( &
    interp_stack_config, YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT_F, &
    YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT_F, &
    YAC_INTERP_SPMAP_WEIGHTED_DEFAULT_F, YAC_INTERP_SPMAP_SCALE_TYPE_DEFAULT_F, &
    YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT_F, &
    YAC_INTERP_SPMAP_FILENAME_DEFAULT_F, YAC_INTERP_SPMAP_VARNAME_DEFAULT_F, &
    YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT_F, &
    YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT_F, &
    YAC_INTERP_SPMAP_FILENAME_DEFAULT_F, YAC_INTERP_SPMAP_VARNAME_DEFAULT_F, &
    YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT_F)
  CALL yac_interp_stack_config_add_spmap_ext_c( &
    interp_stack_config, YAC_INTERP_SPMAP_CONFIG_DEFAULT_F, &
    YAC_INTERP_SPMAP_OVERWRITE_DEFAULT_F )
  CALL yac_interp_stack_config_add_spmap_ext_c( &
    interp_stack_config, spmap_config, C_LOC(overwrite_configs(1)))
  CALL yac_interp_stack_config_add_hcsbb_c(interp_stack_config)
  CALL yac_interp_stack_config_add_user_file_c( &
    interp_stack_config, "test_fortran_api_weights.nc" // c_null_char, &
    YAC_INTERP_FILE_ON_MISSING_FILE_DEFAULT_F, &
    YAC_INTERP_FILE_ON_SUCCESS_DEFAULT_F)
  CALL yac_interp_stack_config_add_fixed_c( &
    interp_stack_config, -1.0_c_double)
  CALL yac_interp_stack_config_add_creep_c( &
    interp_stack_config, YAC_INTERP_CREEP_DISTANCE_DEFAULT_F)
  CALL yac_interp_stack_config_add_check_c( &
    interp_stack_config, YAC_INTERP_CHECK_CONSTRUCTOR_KEY_DEFAULT_F, &
    YAC_INTERP_CHECK_DO_SEARCH_KEY_DEFAULT_F)

  ! clean up data structures for extended spmap configuration
  CALL yac_spmap_overwrite_config_delete_c(overwrite_config)
  CALL yac_interp_spmap_config_delete_c(overwrite_spmap_config)
  CALL yac_point_selection_delete_c(src_point_selection)
  CALL yac_interp_spmap_config_delete_c(spmap_config)
  CALL yac_spmap_scale_config_delete_c(scale_config);
  CALL yac_spmap_cell_area_config_delete_c(target_cell_area_config)
  CALL yac_spmap_cell_area_config_delete_c(source_cell_area_config)

  compare_value = &
    yac_interp_stack_config_compare_c( &
      interp_stack_config, interp_stack_config_cpy)
  CALL test(compare_value /= 0_c_int)

  ! interpolation method stack
  CALL test(constructor_call_count == 0)
  interp_method = &
    yac_interp_stack_config_generate_c(interp_stack_config_cpy)
  CALL test(constructor_call_count == 1)

  ! interpolation method
  CALL test(do_search_call_count == 0)
  interp_weights = &
    yac_interp_method_do_search_c(interp_method, interp_grid)
  CALL yac_interp_weights_write_to_file_c( &
    interp_weights, "test_fortran_api_weights.nc" // c_null_char, &
    "unstruct_grid" // c_null_char, "unstruct_grid_ll" // c_null_char, &
    0_c_size_t, 0_c_size_t, INT(YAC_WEIGHT_FILE_OVERWRITE, c_int))
  CALL test(do_search_call_count == 1)

  ! interpolation
  interpolation = &
    yac_interp_weights_get_interpolation_c( &
      interp_weights, YAC_MAPPING_ON_SRC, 1_c_size_t, &
      yac_interpolation_get_const_frac_mask_no_value_c(), &
      1.0_c_double, 0.0_c_double, c_null_char, 1_c_int, 1_c_int)
  interpolation_frac = &
    yac_interp_weights_get_interpolation_c( &
      interp_weights, YAC_MAPPING_ON_SRC, 1_c_size_t, &
      0.0_c_double, 1.0_c_double, 0.0_c_double, c_null_char, 1_c_int, 1_c_int)
  ltest_value = &
    yac_interpolation_get_const_frac_mask_no_value_c() /= &
    yac_interpolation_get_const_frac_mask_undef_c()
  CALL test(ltest_value)

  src_field_(1) = c_loc(src_cell_data(1))
  src_field_collection(1) = c_loc(src_field_(1))
  src_frac_mask_(1) = c_loc(src_frac_mask)
  src_frac_mask_collection(1) = c_loc(src_frac_mask_(1))
  tgt_field_collection(1) = c_loc(tgt_cell_data(1))
  src_cell_data(1) = 1.0_c_double
  src_frac_mask(1) = 0.0_c_double

  tgt_cell_data(1) = -1.0_c_double
  CALL yac_interpolation_execute_c( &
    interpolation, c_loc(src_field_collection), c_loc(tgt_field_collection))
  CALL test(tgt_cell_data(1) == src_cell_data(1))

  tgt_cell_data(1) = -1.0_c_double
  CALL yac_interpolation_execute_frac_c( &
    interpolation_frac, c_loc(src_field_collection), &
    c_loc(src_frac_mask_collection), c_loc(tgt_field_collection))
  CALL test(tgt_cell_data(1) == 0.0_c_double)

  tgt_cell_data(1) = -1.0_c_double
  test_value = yac_interpolation_execute_put_test_c(interpolation)
  CALL test(test_value == 1_c_int)
  CALL yac_interpolation_execute_put_c( &
    interpolation, c_loc(src_field_collection))
  CALL yac_interpolation_execute_get_c( &
    interpolation, c_loc(tgt_field_collection))
  CALL test(tgt_cell_data(1) == src_cell_data(1))

  tgt_cell_data(1) = -1.0_c_double
  CALL yac_interpolation_execute_get_async_c( &
    interpolation_frac, c_loc(tgt_field_collection))
  test_value = yac_interpolation_execute_get_test_c(interpolation_frac)
  CALL test(test_value == 0_c_int)
  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
  CALL yac_interpolation_execute_put_frac_c( &
    interpolation_frac, c_loc(src_field_collection), &
    c_loc(src_frac_mask_collection))
  CALL yac_interpolation_execute_wait_c(interpolation_frac)
  CALL test(tgt_cell_data(1) == 0.0_c_double)
  tgt_cell_data(1) = -1.0_c_double

  ! cleanup
  CALL yac_interpolation_delete_c(interpolation_frac)
  CALL yac_interpolation_delete_c(interpolation)
  CALL C_UNLINK("test_fortran_api_weights.nc" // c_null_char)
  CALL yac_interp_weights_delete_c(interp_weights)
  CALL yac_interp_method_delete_c(interp_method)
  CALL c_free(interp_method)
  CALL yac_interp_stack_config_delete_c(interp_stack_config_cpy)
  CALL yac_interp_stack_config_delete_c(interp_stack_config)
  CALL yac_interp_method_cleanup_c()
  CALL yac_interp_grid_delete_c(interp_grid)
  CALL yac_dist_grid_pair_delete_c(grid_pair)
  CALL yac_basic_grid_delete_c(basic_grid_reg_2d_rot_deg)
  CALL yac_basic_grid_delete_c(basic_grid_reg_2d_rot)
  CALL yac_basic_grid_delete_c(basic_grid_cloud_deg)
  CALL yac_basic_grid_delete_c(basic_grid_cloud)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_edge_ll_deg)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_edge_ll)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_edge_deg)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_edge)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_ll_deg)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_ll)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct_deg)
  CALL yac_basic_grid_delete_c(basic_grid_unstruct)
  CALL yac_basic_grid_delete_c(basic_grid_curve_2d_deg)
  CALL yac_basic_grid_delete_c(basic_grid_curve_2d)
  CALL yac_basic_grid_delete_c(basic_grid_reg_2d_deg)
  CALL yac_basic_grid_delete_c(basic_grid_reg_2d)
  CALL yac_basic_grid_delete_c(basic_grid_empty)

  CALL yac_mpi_cleanup_c()
  CALL yac_mpi_finalize_c()

  CALL stop_test
  CALL exit_tests

CONTAINS

  SUBROUTINE constructor_callback(user_data) BIND(C)

    TYPE(c_ptr), value :: user_data

    INTEGER, POINTER :: constructor_call_count

    CALL c_f_pointer(user_data, constructor_call_count)

    constructor_call_count = constructor_call_count + 1

  END SUBROUTINE constructor_callback

  SUBROUTINE do_search_callback( &
    global_ids, coordinates_xyz, count, user_data) BIND(C)

    TYPE(c_ptr), value :: global_ids
    TYPE(c_ptr), value :: coordinates_xyz
    INTEGER(kind=c_size_t), value :: count
    TYPE(c_ptr), value :: user_data

    INTEGER, POINTER :: do_search_call_count

    NOP(global_ids)
    NOP(coordinates_xyz)
    NOP(count)

    CALL c_f_pointer(user_data, do_search_call_count)

    do_search_call_count = do_search_call_count + 1

  END SUBROUTINE

END PROGRAM main
