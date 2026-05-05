! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_def_points_dble.F90
!! \test
!! This shows how to use yac_fdef_points.

program test_def_points_dble

#define TEST_PRECISION dp
#define YAC_PTR_TYPE yac_dble_ptr

#include "test_def_points.inc"

end program test_def_points_dble
