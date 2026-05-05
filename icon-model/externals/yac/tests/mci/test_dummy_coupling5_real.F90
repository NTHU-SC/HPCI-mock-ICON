! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_dummy_coupling5_real.F90
!! \test
!! Fortran version of \ref test_dummy_coupling5_c.c

#define TEST_PRECISION sp
#define YAC_PTR_TYPE yac_real_ptr

PROGRAM test_dummy_coupling5_real

#include "test_dummy_coupling5.inc"

END PROGRAM test_dummy_coupling5_real
