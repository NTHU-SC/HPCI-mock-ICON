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

!>
!!   Contains real and integer kinds for libfortran-support and libiconmath
!!
MODULE mo_iconlib_kind

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: real64, real32, int32, int64

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wp !< Selected working precision
  PUBLIC :: vp !< Selected variable precision
  PUBLIC :: dp !< 8 byte real (double precision)
  PUBLIC :: sp !< 4 byte real (single precision)
  PUBLIC :: i8 !< 8 byte integer
  PUBLIC :: i4 !< 4 byte integer

  INTEGER, PARAMETER :: dp = real64
  INTEGER, PARAMETER :: sp = real32
  INTEGER, PARAMETER :: i8 = int64
  INTEGER, PARAMETER :: i4 = int32

#ifdef __SINGLE_PRECISION
  INTEGER, PARAMETER :: wp = sp
#else
  INTEGER, PARAMETER :: wp = dp
#endif

#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif

END MODULE mo_iconlib_kind
