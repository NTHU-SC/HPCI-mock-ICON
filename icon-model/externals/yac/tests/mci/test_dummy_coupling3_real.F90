! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_dummy_coupling3_real.F90
!! \test
!! Fortran version of \ref test_dummy_coupling3_c.c

#define TEST_PRECISION SELECTED_REAL_KIND(6, 37)

PROGRAM test_dummy_coupling3_real

#include "test_dummy_coupling3.inc"

END PROGRAM test_dummy_coupling3_real
