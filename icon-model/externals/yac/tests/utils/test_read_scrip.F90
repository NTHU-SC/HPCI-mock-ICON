! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

!> \file test_read_scrip.F90
!! \test
!! This contains Fortran examples for yac_read_scrip_basic_grid

PROGRAM main

  USE, INTRINSIC :: iso_c_binding
  USE utest
  USE yac_core
  USE yac_utils

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

  ! ===================================================================

  CALL start_test('read_scrip')

  ! ===================================================================
  ! yac_read_scrip_basic_grid
  ! ===================================================================
  CALL test_yac_read_scrip_basic_grid()

  ! ===================================================================
  ! yac_read_scrip_cloud_basic_grid
  ! ===================================================================
  CALL test_yac_read_scrip_cloud_basic_grid()

  CALL stop_test
  CALL exit_tests

CONTAINS

  SUBROUTINE test_yac_read_scrip_basic_grid()

    IMPLICIT NONE

    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_name = &
      "dummy_grid" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_filename = &
      "test_read_scrip_grids.nc" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: mask_filename = &
      "test_read_scrip_masks.nc" // c_null_char

    INTEGER(KIND=C_INT), PARAMETER :: with_corner = 1_c_int
    INTEGER(KIND=C_INT), PARAMETER :: valid_mask_value = 0_c_int
    INTEGER(KIND=C_INT), PARAMETER :: use_ll_edges = 0_c_int

    INTEGER(KIND=C_SIZE_T) :: cell_coord_idx
    TYPE(C_PTR)            :: duplicated_cell_idx
    TYPE(C_PTR)            :: orig_cell_global_id
    INTEGER(KIND=C_SIZE_T) :: nbr_duplicated_cells

    TYPE(C_PTR) :: scrip_grid

    ! write dummy scrip grid file
    CALL write_dummy_scrip_grid_file_c( &
      grid_name, grid_filename, mask_filename, with_corner, &
      360_c_size_t, 10_c_size_t, (/0.0_c_double, 360.0_c_double/), &
      (/0.0_c_double, 10.0_c_double/))

    scrip_grid = &
      yac_read_scrip_basic_grid_c( &
        grid_filename, mask_filename, grid_name, &
        valid_mask_value, grid_name, use_ll_edges, cell_coord_idx, &
        duplicated_cell_idx, orig_cell_global_id, nbr_duplicated_cells)

    CALL yac_free_c(duplicated_cell_idx)
    CALL yac_free_c(orig_cell_global_id)

    CALL yac_basic_grid_delete_c(scrip_grid)

    ! delete dummy scrip grid file
    CALL C_UNLINK(grid_filename)
    CALL C_UNLINK(mask_filename)

  END SUBROUTINE test_yac_read_scrip_basic_grid

  SUBROUTINE test_yac_read_scrip_cloud_basic_grid()

    IMPLICIT NONE

    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_name = &
      "dummy_grid" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_filename = &
      "test_read_scrip_cloud_grids.nc" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: mask_filename = &
      "test_read_scrip_cloud_masks.nc" // c_null_char

    INTEGER(KIND=C_INT), PARAMETER :: with_corner = 0_c_int
    INTEGER(KIND=C_INT), PARAMETER :: valid_mask_value = 0_c_int

    INTEGER(KIND=C_SIZE_T) :: vertex_coord_idx
    TYPE(C_PTR)            :: duplicated_vertex_idx
    TYPE(C_PTR)            :: orig_vertex_global_id
    INTEGER(KIND=C_SIZE_T) :: nbr_duplicated_vertices

    TYPE(C_PTR) :: scrip_cloud_grid

    ! write dummy scrip grid file
    CALL write_dummy_scrip_grid_file_c( &
      grid_name, grid_filename, mask_filename, with_corner, &
      380_c_size_t, 180_c_size_t, (/0.0_c_double, 380.0_c_double/), &
      (/-90.0_c_double, 90.0_c_double/))

    scrip_cloud_grid = &
      yac_read_scrip_cloud_basic_grid_c( &
        grid_filename, mask_filename, grid_name, &
        valid_mask_value, grid_name, vertex_coord_idx, &
        duplicated_vertex_idx, orig_vertex_global_id, nbr_duplicated_vertices)

    CALL test(nbr_duplicated_vertices == 20 * 180)

    CALL yac_free_c(duplicated_vertex_idx)
    CALL yac_free_c(orig_vertex_global_id)

    CALL yac_basic_grid_delete_c(scrip_cloud_grid)

    ! delete dummy scrip grid file
    CALL C_UNLINK(grid_filename)
    CALL C_UNLINK(mask_filename)

  END SUBROUTINE test_yac_read_scrip_cloud_basic_grid

END PROGRAM main
