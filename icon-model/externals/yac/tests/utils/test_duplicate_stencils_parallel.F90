! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_duplicate_stencils_parallel.F90
!! \test
!! Fortran test for stencil duplication.

#include "test_macros.inc"

PROGRAM main

  USE, INTRINSIC :: iso_c_binding
  USE utest
  USE yac_core
  USE yac_utils
  USE mpi

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE write_dummy_scrip_grid_file_c ( &
      grid_name, grid_filename, mask_filename, &
      with_corner, num_lon, num_lat, lon_range, lat_range ) &
      BIND ( c, name='write_dummy_scrip_grid_file' )
      USE, INTRINSIC :: iso_c_binding, only : c_char, c_int, c_size_t, c_double
      CHARACTER(KIND=c_char), DIMENSION(*) :: grid_name
      CHARACTER(KIND=c_char), DIMENSION(*) :: grid_filename
      CHARACTER(KIND=c_char), DIMENSION(*) :: mask_filename
      INTEGER(KIND=c_int), value           :: with_corner
      INTEGER(KIND=c_size_t), value        :: num_lon
      INTEGER(KIND=c_size_t), value        :: num_lat
      REAL(KIND=c_double), DIMENSION(*)    :: lat_range
      REAL(KIND=c_double), DIMENSION(*)    :: lon_range
    END SUBROUTINE write_dummy_scrip_grid_file_c
  END INTERFACE

  INTERFACE
    SUBROUTINE C_UNLINK ( path ) BIND ( c, name='unlink' )
      USE, INTRINSIC :: iso_c_binding, only : c_char
      CHARACTER(KIND=c_char), DIMENSION(*) :: path
    END SUBROUTINE C_UNLINK
  END INTERFACE

  INTEGER :: comm_size, comm_rank
  INTEGER :: ierror

  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: src_grid_name = &
    "source_grid" // c_null_char
  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: src_grid_filename = &
    "test_duplicate_stencils_parallel_src_grids.nc" // c_null_char
  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: src_mask_filename = &
    "test_duplicate_stencils_parallel_src_masks.nc" // c_null_char
  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: tgt_grid_name = &
    "target_grid" // c_null_char
  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: tgt_grid_filename = &
    "test_duplicate_stencils_parallel_tgt_grids.nc" // c_null_char
  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: tgt_mask_filename = &
    "test_duplicate_stencils_parallel_tgt_masks.nc" // c_null_char

  INTEGER(KIND=C_INT), PARAMETER :: with_corner = 1_c_int
  INTEGER(KIND=C_INT), PARAMETER :: valid_mask_value = 0_c_int
  INTEGER(KIND=C_INT), PARAMETER :: use_ll_edges = 0_c_int

  INTEGER(KIND=C_SIZE_T) :: src_cell_coord_idx
  TYPE(C_PTR)            :: src_duplicated_cell_idx
  TYPE(C_PTR)            :: src_orig_cell_global_id
  INTEGER(KIND=C_SIZE_T) :: src_nbr_duplicated_cells

  INTEGER(KIND=C_SIZE_T) :: tgt_cell_coord_idx
  TYPE(C_PTR)            :: tgt_duplicated_cell_idx
  TYPE(C_PTR)            :: tgt_orig_cell_global_id
  INTEGER(KIND=C_SIZE_T) :: tgt_nbr_duplicated_cells

  TYPE(c_ptr) :: src_grid, tgt_grid
  TYPE(c_ptr) :: dist_grid_pair
  TYPE(c_ptr) :: interp_grid
  TYPE(c_ptr) :: interp_stack_config
  TYPE(c_ptr) :: interp_method_stack
  TYPE(c_ptr) :: interp_weights

  ! ===================================================================

  CALL start_test('duplicate_stencils_parallel')

  CALL yac_mpi_init_c()
  CALL yac_yaxt_init_c(MPI_COMM_WORLD)

  CALL MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierror)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierror)

  CALL test(comm_size == 3)

  ! write source and target scrip grid files
  IF (comm_rank == 0) THEN
    CALL write_dummy_scrip_grid_file_c( &
      src_grid_name, src_grid_filename, src_mask_filename, with_corner, &
      360_c_size_t, 180_c_size_t, (/0.0_c_double, 360.0_c_double/), &
      (/-90.0_c_double, 90.0_c_double/))
    CALL write_dummy_scrip_grid_file_c( &
      tgt_grid_name, tgt_grid_filename, tgt_mask_filename, with_corner, &
      380_c_size_t, 180_c_size_t, (/0.0_c_double, 380.0_c_double/), &
      (/-90.0_c_double, 90.0_c_double/))
  END IF
  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! read in source and target grid on all processes
  src_grid = &
    yac_read_scrip_basic_grid_parallel_c( &
      src_grid_filename, src_mask_filename, MPI_COMM_WORLD, src_grid_name, &
      valid_mask_value, src_grid_name, use_ll_edges, src_cell_coord_idx, &
      src_duplicated_cell_idx, src_orig_cell_global_id, &
      src_nbr_duplicated_cells)
  tgt_grid = &
    yac_read_scrip_basic_grid_parallel_c( &
      tgt_grid_filename, tgt_mask_filename, MPI_COMM_WORLD, tgt_grid_name, &
      valid_mask_value, tgt_grid_name, use_ll_edges, tgt_cell_coord_idx, &
      tgt_duplicated_cell_idx, tgt_orig_cell_global_id, &
      tgt_nbr_duplicated_cells)

  ! generate distributed grid pair
  dist_grid_pair = &
    yac_dist_grid_pair_new_c(src_grid, tgt_grid, MPI_COMM_WORLD)

  ! generate interpolation grid
  interp_grid = &
    yac_interp_grid_new_c( &
      dist_grid_pair, TRIM(src_grid_name) // c_null_char, &
      TRIM(tgt_grid_name) // c_null_char,  INT(1, c_size_t), &
      (/INT(YAC_LOC_CELL, c_int)/), (/src_cell_coord_idx/), (/-1_c_size_t/), &
      INT(YAC_LOC_CELL, c_int), tgt_cell_coord_idx, -1_c_size_t)

  ! configure the interpolation stack
  interp_stack_config = yac_interp_stack_config_new_c()
  CALL yac_interp_stack_config_add_conservative_c( &
    interp_stack_config, 1, 0, 1, YAC_INTERP_CONSERV_DESTAREA)

  ! generate the actual interpolation stack
  interp_method_stack = &
    yac_interp_stack_config_generate_c(interp_stack_config)

  ! execute the interpolation stack and generate the weights
  !   YAC starts by extracting all non-masked target points, which
  !   are then passed to the interpolation stack.
  !   The resulting interpolation weights contains the interpolation
  !   stencils, which are distributed across all processes.
  !   (this operation is collective)
  interp_weights = &
    yac_interp_method_do_search_c(interp_method_stack, interp_grid)

  CALL yac_duplicate_stencils_c( &
    interp_weights, tgt_grid, tgt_orig_cell_global_id, &
    tgt_duplicated_cell_idx, tgt_nbr_duplicated_cells, YAC_LOC_CELL)

  ! cleanup
  CALL yac_interp_weights_delete_c(interp_weights)
  CALL yac_interp_method_delete_c(interp_method_stack)
  CALL yac_free_c(interp_method_stack)
  CALL yac_interp_stack_config_delete_c(interp_stack_config)
  CALL yac_interp_grid_delete_c(interp_grid)
  CALL yac_dist_grid_pair_delete_c(dist_grid_pair)
  CALL yac_free_c(tgt_duplicated_cell_idx)
  CALL yac_free_c(tgt_orig_cell_global_id)
  CALL yac_basic_grid_delete_c(tgt_grid)
  CALL yac_free_c(src_duplicated_cell_idx)
  CALL yac_free_c(src_orig_cell_global_id)
  CALL yac_basic_grid_delete_c(src_grid)

  ! delete grid files scrip grid file
  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
  IF (comm_rank == 0) THEN
    CALL C_UNLINK(src_grid_filename)
    CALL C_UNLINK(src_mask_filename)
    CALL C_UNLINK(tgt_grid_filename)
    CALL C_UNLINK(tgt_mask_filename)
  END IF

  CALL yac_mpi_finalize_c()

  CALL stop_test
  CALL exit_tests

END PROGRAM main
