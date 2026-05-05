! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

!> \file test_read_scrip_parallel.F90
!! \test
!! This contains Fortran examples for yac_read_scrip_basic_grid_parallel

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

  ! ===================================================================

  CALL start_test('read_scrip_parallel')

  CALL MPI_Init(ierror)


  CALL MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierror)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierror)

  CALL test(comm_size == 4)

  ! ===================================================================
  ! yac_read_scrip_basic_grid_parallel
  ! ===================================================================
  CALL test_yac_read_scrip_basic_grid_parallel()

  ! ===================================================================
  ! yac_read_scrip_cloud_basic_grid_parallel
  ! ===================================================================
  CALL test_yac_read_scrip_cloud_basic_grid_parallel()

  ! ===================================================================
  ! yac_read_scrip_generic_basic_grid_parallel
  ! ===================================================================
  CALL test_yac_read_scrip_generic_basic_grid_parallel()

  CALL MPI_Finalize(ierror)

  CALL stop_test
  CALL exit_tests

CONTAINS

  SUBROUTINE test_yac_read_scrip_basic_grid_parallel()

    IMPLICIT NONE

    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_name = &
      "dummy_grid" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_filename = &
      "test_read_scrip_parallel_grids.nc" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: mask_filename = &
      "test_read_scrip_parallel_masks.nc" // c_null_char

    INTEGER(KIND=C_INT), PARAMETER :: with_corner = 1_c_int
    INTEGER(KIND=C_INT), PARAMETER :: valid_mask_value = 0_c_int
    INTEGER(KIND=C_INT), PARAMETER :: use_ll_edges = 0_c_int

    INTEGER(KIND=C_SIZE_T) :: cell_coord_idx
    TYPE(C_PTR)            :: duplicated_cell_idx
    TYPE(C_PTR)            :: orig_cell_global_id
    INTEGER(KIND=C_SIZE_T) :: nbr_duplicated_cells

    TYPE(C_PTR) :: scrip_grid

    ! write dummy scrip grid file
    IF (comm_rank == 0) &
      CALL write_dummy_scrip_grid_file_c( &
        grid_name, grid_filename, mask_filename, with_corner, &
        360_c_size_t, 180_c_size_t, (/0.0_c_double, 360.0_c_double/), &
        (/-90.0_c_double, 90.0_c_double/))
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

    scrip_grid = &
      yac_read_scrip_basic_grid_parallel_c( &
        grid_filename, mask_filename, MPI_COMM_WORLD, grid_name, &
        valid_mask_value, grid_name, use_ll_edges, cell_coord_idx, &
        duplicated_cell_idx, orig_cell_global_id, nbr_duplicated_cells)

    CALL yac_free_c(duplicated_cell_idx)
    CALL yac_free_c(orig_cell_global_id)

    CALL yac_basic_grid_delete_c(scrip_grid)

    ! delete dummy scrip grid file
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
    IF (comm_rank == 0) THEN
      CALL C_UNLINK(grid_filename)
      CALL C_UNLINK(mask_filename)
    END IF

  END SUBROUTINE test_yac_read_scrip_basic_grid_parallel

  SUBROUTINE test_yac_read_scrip_cloud_basic_grid_parallel()

    IMPLICIT NONE

    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_name = &
      "dummy_grid" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_filename = &
      "test_read_scrip_cloud_parallel_grids.nc" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: mask_filename = &
      "test_read_scrip_cloud_parallel_masks.nc" // c_null_char

    INTEGER(KIND=C_INT), PARAMETER :: with_corner = 0_c_int
    INTEGER(KIND=C_INT), PARAMETER :: valid_mask_value = 0_c_int

    INTEGER(KIND=C_SIZE_T) :: vertex_coord_idx
    TYPE(C_PTR)            :: duplicated_vertex_idx
    TYPE(C_PTR)            :: orig_vertex_global_id
    INTEGER(KIND=C_SIZE_T) :: nbr_duplicated_vertices

    TYPE(C_PTR) :: scrip_cloud_grid

    INTEGER :: local_nbr_duplicated_vertices(1)
    INTEGER :: global_nbr_duplicated_vertices(1)

    ! write dummy scrip grid file
    IF (comm_rank == 0) &
      CALL write_dummy_scrip_grid_file_c( &
        grid_name, grid_filename, mask_filename, with_corner, &
        380_c_size_t, 180_c_size_t, (/0.0_c_double, 380.0_c_double/), &
        (/-90.0_c_double, 90.0_c_double/))
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

    scrip_cloud_grid = &
      yac_read_scrip_cloud_basic_grid_parallel_c( &
        grid_filename, mask_filename, MPI_COMM_WORLD, grid_name, &
        valid_mask_value, grid_name, vertex_coord_idx, &
        duplicated_vertex_idx, orig_vertex_global_id, nbr_duplicated_vertices)

    local_nbr_duplicated_vertices(1) = INT(nbr_duplicated_vertices)
    CALL MPI_Allreduce( &
      local_nbr_duplicated_vertices, global_nbr_duplicated_vertices, &
      1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    CALL test(global_nbr_duplicated_vertices(1) == 20 * 180)

    CALL yac_free_c(duplicated_vertex_idx)
    CALL yac_free_c(orig_vertex_global_id)

    CALL yac_basic_grid_delete_c(scrip_cloud_grid)

    ! delete dummy scrip grid file
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
    IF (comm_rank == 0) THEN
      CALL C_UNLINK(grid_filename)
      CALL C_UNLINK(mask_filename)
    END IF

  END SUBROUTINE test_yac_read_scrip_cloud_basic_grid_parallel

  SUBROUTINE test_yac_read_scrip_generic_basic_grid_parallel()

    IMPLICIT NONE

    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_name = &
      "dummy_grid" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_filename = &
      "test_read_scrip_parallel_2_grids.nc" // c_null_char
    CHARACTER(KIND=c_char,LEN=*), PARAMETER :: mask_filename = &
      "test_read_scrip_parallel_2_masks.nc" // c_null_char

    INTEGER(KIND=C_INT) :: with_corner
    INTEGER(KIND=C_INT), PARAMETER :: valid_mask_value = 0_c_int
    INTEGER(KIND=C_INT), PARAMETER :: use_ll_edges = 0_c_int

    INTEGER(KIND=C_SIZE_T) :: cell_coord_idx
    TYPE(C_PTR)            :: duplicated_cell_idx
    TYPE(C_PTR)            :: orig_cell_global_id
    INTEGER(KIND=C_SIZE_T) :: nbr_duplicated_cells

    INTEGER(KIND=C_SIZE_T) :: vertex_coord_idx
    TYPE(C_PTR)            :: duplicated_vertex_idx
    TYPE(C_PTR)            :: orig_vertex_global_id
    INTEGER(KIND=C_SIZE_T) :: nbr_duplicated_vertices

    INTEGER(KIND=C_INT)    :: point_location

    TYPE(C_PTR) :: scrip_grid
    TYPE(C_PTR) :: scrip_cloud_grid

    ! write dummy scrip grid file
    IF (comm_rank == 0) THEN
      with_corner = 1_c_int
      CALL write_dummy_scrip_grid_file_c( &
        grid_name, grid_filename, mask_filename, with_corner, &
        360_c_size_t, 180_c_size_t, (/0.0_c_double, 360.0_c_double/), &
        (/-90.0_c_double, 90.0_c_double/))
    END IF
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

    scrip_grid = &
      yac_read_scrip_generic_basic_grid_parallel_c( &
        grid_filename, mask_filename, MPI_COMM_WORLD, grid_name, &
        valid_mask_value, grid_name, use_ll_edges, cell_coord_idx, &
        duplicated_cell_idx, orig_cell_global_id, nbr_duplicated_cells, &
        point_location)

    CALL test(point_location == INT(YAC_LOC_CELL, C_INT))

    CALL yac_free_c(duplicated_cell_idx)
    CALL yac_free_c(orig_cell_global_id)

    CALL yac_basic_grid_delete_c(scrip_grid)

    ! delete dummy scrip grid file
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
    IF (comm_rank == 0) THEN
      CALL C_UNLINK(grid_filename)
      CALL C_UNLINK(mask_filename)
    END IF

    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

    ! write dummy scrip cloud grid file
    IF (comm_rank == 0) THEN
      with_corner = 0_c_int
      CALL write_dummy_scrip_grid_file_c( &
        grid_name, grid_filename, mask_filename, with_corner, &
        380_c_size_t, 180_c_size_t, (/0.0_c_double, 380.0_c_double/), &
        (/-90.0_c_double, 90.0_c_double/))
    END IF
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

    scrip_cloud_grid = &
      yac_read_scrip_generic_basic_grid_parallel_c( &
        grid_filename, mask_filename, MPI_COMM_WORLD, grid_name, &
        valid_mask_value, grid_name, use_ll_edges, vertex_coord_idx, &
        duplicated_vertex_idx, orig_vertex_global_id, nbr_duplicated_vertices, &
        point_location)

    CALL test(point_location == INT(YAC_LOC_CORNER, C_INT))

    CALL yac_free_c(duplicated_vertex_idx)
    CALL yac_free_c(orig_vertex_global_id)

    CALL yac_basic_grid_delete_c(scrip_cloud_grid)

    ! delete dummy scrip grid file
    CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
    IF (comm_rank == 0) THEN
      CALL C_UNLINK(grid_filename)
      CALL C_UNLINK(mask_filename)
    END IF

  END SUBROUTINE test_yac_read_scrip_generic_basic_grid_parallel

END PROGRAM main
