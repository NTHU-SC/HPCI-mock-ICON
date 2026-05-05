! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_def_grid_real.F90
!! \test
!! This shows how to use yac_fdef_grid

program test_def_grid_real

#define TEST_PRECISION sp
#define YAC_PTR_TYPE yac_real_ptr

#include "test_def_grid.inc"

end program test_def_grid_real
