! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_dummy_coupling_raw_real.F90
!! \test
!! This test checks the raw data exchange feature.

#define TEST_PRECISION SELECTED_REAL_KIND(6, 37)
#define YAC_PTR_TYPE yac_real_ptr
!#define TEST_GET_ASYNC

PROGRAM dummy_coupling_raw_real

#include "test_dummy_coupling_raw.inc"

END PROGRAM dummy_coupling_raw_real
