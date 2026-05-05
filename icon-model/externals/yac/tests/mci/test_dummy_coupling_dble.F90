! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#define TEST_PRECISION SELECTED_REAL_KIND(15, 307)
#define TEST_GET_ASYNC

!> \file test_dummy_coupling_dble.F90
!! \test
!! This example simulates a whole model setup with three components (ocean,
!! atmosphere, io). It uses one process for each component.

PROGRAM dummy_coupling_dble

#include "test_dummy_coupling.inc"

END PROGRAM dummy_coupling_dble
