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

MODULE helpers
  USE mo_iconlib_kind, ONLY: wp

CONTAINS

  FUNCTION calculate_mean_wp(array) RESULT(mean)
    REAL(wp), INTENT(IN) :: array(:)
    REAL(wp)             :: mean
    INTEGER :: i, n

    mean = 0.0_wp
    variance = 0.0_wp

    n = SIZE(array)

    DO i = 1, SIZE(array)
      mean = mean + array(i)
    END DO
    mean = mean/n

  END FUNCTION

  FUNCTION calculate_variance_wp(array) RESULT(variance)
    REAL(wp), INTENT(IN) :: array(:)
    REAL(wp)             :: mean, variance
    INTEGER :: i, n

    mean = 0.0_wp
    variance = 0.0_wp

    n = SIZE(array)
    mean = calculate_mean_wp(array)

    DO i = 1, SIZE(array)
      variance = variance + (array(i) - mean)**2
    END DO
    variance = variance/n

  END FUNCTION

  SUBROUTINE assert_statistics(test_name, array, mean, variance, tol, max, min)
    USE FORTUTF
    CHARACTER(LEN=*), INTENT(IN) :: test_name
    REAL(wp), INTENT(IN) :: array(:)
    REAL(wp), INTENT(IN) :: mean, variance, tol, max, min

    CALL TAG_TEST(test_name//"__mean")
    CALL ASSERT_ALMOST_EQUAL(calculate_mean_wp(array), mean, tol)
    CALL TAG_TEST(test_name//"__variance")
    CALL ASSERT_ALMOST_EQUAL(calculate_variance_wp(array), variance, tol)
    CALL TAG_TEST(test_name//"__max")
    CALL ASSERT_LESS_THAN_EQUAL(MAXVAL(array), max)
    CALL TAG_TEST(test_name//"__min")
    CALL ASSERT_GREATER_THAN_EQUAL(MINVAL(array), min)

  END SUBROUTINE

END MODULE helpers
