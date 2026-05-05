! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

!> \file test_read_icon_parallel.F90
!! \test
!! This contains examples for parallel use of yac_read_icon_basic_grid_data.

PROGRAM main

  USE, INTRINSIC :: iso_c_binding
  USE utest
  USE yac_core
  USE yac_utils
  USE mpi

  IMPLICIT NONE

  CHARACTER(KIND=c_char,LEN=*), PARAMETER :: grid_filename = &
    "test_read_icon_parallel_grid.nc" // c_null_char

  INTERFACE
    SUBROUTINE write_test_grid_file_c ( grid_filename, coord_unit ) &
      BIND ( c, name='write_test_grid_file_f2c' )
      USE, INTRINSIC :: iso_c_binding, only : c_char, c_int
      CHARACTER(KIND=c_char), DIMENSION(*) :: grid_filename
      INTEGER(KIND=c_int), value           :: coord_unit
    END SUBROUTINE write_test_grid_file_c
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

  CALL start_test('read_icon_parallel')

  CALL MPI_Init(ierror)


  CALL MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierror)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierror)

  CALL test(comm_size == 4)

  CALL write_test_grid_file_c(grid_filename, 0_c_int)

  ! ===================================================================
  ! yac_read_icon_grid_information_parallel (incl. halos)
  ! ===================================================================
  BLOCK

    INTEGER(KIND=c_int), PARAMETER :: ref_num_cells(4) = (/11,9,14,9/)
    INTEGER(KIND=c_int), PARAMETER :: ref_num_vertices(4) = (/11,10,13,10/)

    INTEGER(KIND=C_INT) :: num_vertices
    INTEGER(KIND=C_INT) :: num_cells
    TYPE(C_PTR)         :: num_vertices_per_cell ! int **
    TYPE(C_PTR)         :: cell_to_vertex ! int **
    TYPE(C_PTR)         :: global_cell_id ! int **
    TYPE(C_PTR)         :: cell_owner ! int **
    TYPE(C_PTR)         :: global_vertex_ids ! int **
    TYPE(C_PTR)         :: vertex_owner ! int **
    TYPE(C_PTR)         :: x_vertices ! double **
    TYPE(C_PTR)         :: y_vertices ! double **
    TYPE(C_PTR)         :: x_cells ! double **
    TYPE(C_PTR)         :: y_cells ! double **
    TYPE(C_PTR)         :: cell_mask ! int **

  CALL yac_read_icon_grid_information_parallel_c( &
    grid_filename, MPI_COMM_WORLD, num_vertices, num_cells, &
    num_vertices_per_cell, cell_to_vertex, global_cell_id, &
    cell_owner, global_vertex_ids, vertex_owner, &
    x_vertices, y_vertices, x_cells, y_cells, cell_mask)

  CALL test(num_vertices == ref_num_vertices(comm_rank + 1))
  CALL test(num_cells == ref_num_cells(comm_rank + 1))

  CALL yac_free_c(num_vertices_per_cell)
  CALL yac_free_c(cell_to_vertex)
  CALL yac_free_c(global_cell_id)
  CALL yac_free_c(cell_owner)
  CALL yac_free_c(global_vertex_ids)
  CALL yac_free_c(vertex_owner)
  CALL yac_free_c(x_vertices)
  CALL yac_free_c(y_vertices)
  CALL yac_free_c(x_cells)
  CALL yac_free_c(y_cells)
  CALL yac_free_c(cell_mask)
  END BLOCK

  ! ===================================================================
  ! yac_read_icon_basic_grid_parallel
  ! ===================================================================
  BLOCK

    TYPE(C_PTR) :: icon_grid
    icon_grid = &
      yac_read_icon_basic_grid_parallel_c( &
        grid_filename, TRIM("icon_grid")//c_null_char, MPI_COMM_WORLD)
    CALL test(yac_basic_grid_get_data_size_c(icon_grid, YAC_LOC_CELL) == 4)
    CALL test(yac_basic_grid_get_data_size_c(icon_grid, YAC_LOC_CORNER) == 6)
    CALL yac_basic_grid_delete_c(icon_grid)

  END BLOCK

  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
  IF (comm_rank == 0) THEN
    CALL C_UNLINK(grid_filename)
  END IF

  CALL MPI_Finalize(ierror)

  CALL stop_test
  CALL exit_tests

END PROGRAM main
