! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_read_icon.F90
!! \test
!! This contains examples for yac_read_icon_basic_grid_data.

#include "test_macros.inc"

PROGRAM main

  USE, INTRINSIC :: iso_c_binding
  USE utest
  USE yac_core
  USE yac_utils

  CHARACTER(LEN=1024) :: grid_dir

  INTERFACE
    FUNCTION yac_file_exists_c ( path ) BIND ( c, name='yac_file_exists' )
      USE, INTRINSIC :: iso_c_binding, only : c_char, c_int
      CHARACTER(KIND=c_char), DIMENSION(*) :: path
      INTEGER(KIND=c_int) :: yac_file_exists_c
    END FUNCTION yac_file_exists_c
  END INTERFACE

  ! ===================================================================

  CALL start_test('read_icon')

  CALL test(command_argument_count() == 1)

  IF (command_argument_count() /= 1) THEN
    PRINT *, "ERROR: missing grid file directory"
    CALL stop_test
    CALL exit_tests
  ELSE
    CALL get_command_argument(1, grid_dir)
  END IF

  IF (0_c_int == &
    yac_file_exists_c( &
      TRIM(grid_dir)//"icon_grid_0030_R02B03_G.nc"//c_null_char)) THEN

    STOP EXIT_SKIP_TEST
  END IF

  ! ===================================================================
  ! yac_read_icon_basic_grid
  ! ===================================================================
  BLOCK
    TYPE(c_ptr) :: icon_grid

    icon_grid = &
      yac_read_icon_basic_grid_c( &
        TRIM(grid_dir)//"icon_grid_0030_R02B03_G.nc"//c_null_char, &
        "icon_grid")

    CALL test(yac_basic_grid_get_data_size_c(icon_grid, YAC_LOC_CELL) == 5120)

    CALL yac_basic_grid_delete_c(icon_grid)
  END BLOCK

  CALL stop_test
  CALL exit_tests

END PROGRAM main
