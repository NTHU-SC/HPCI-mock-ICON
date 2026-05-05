! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE TEST_mo_util_backtrace
  USE FORTUTF
  USE fortran_support, ONLY: ftn_util_backtrace

CONTAINS

  SUBROUTINE TEST_ftn_util_backtrace
    CALL TAG_TEST("TEST_ftn_util_backtrace")
    CALL ftn_util_backtrace()
  END SUBROUTINE

END MODULE TEST_mo_util_backtrace
