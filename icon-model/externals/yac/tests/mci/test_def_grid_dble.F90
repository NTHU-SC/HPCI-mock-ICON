! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_def_grid_dble.F90
!! \test
!! This shows how to use yac_fdef_grid

program test_def_grid_dble

#define TEST_PRECISION dp
#define YAC_PTR_TYPE yac_dble_ptr

#include "test_def_grid.inc"

end program test_def_grid_dble
