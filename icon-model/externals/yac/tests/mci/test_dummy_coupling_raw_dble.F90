! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_dummy_coupling_raw_dble.F90
!! \test
!! This test checks the raw data exchange feature.

#define TEST_PRECISION SELECTED_REAL_KIND(15, 307)
#define YAC_PTR_TYPE yac_dble_ptr
#define TEST_GET_ASYNC

PROGRAM dummy_coupling_raw_dble

#include "test_dummy_coupling_raw.inc"

END PROGRAM dummy_coupling_raw_dble
