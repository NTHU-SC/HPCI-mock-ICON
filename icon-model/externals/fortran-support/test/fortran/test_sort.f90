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

MODULE TEST_mo_util_sort
  USE FORTUTF
  USE mo_iconlib_kind, ONLY: wp, dp, sp
  USE fortran_support, ONLY: quicksort, insertion_sort

CONTAINS
  SUBROUTINE TEST_quicksort_real_dp
    REAL(dp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    CALL TAG_TEST("TEST_quicksort_real_dp_before")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_dp_after")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real_dp2
    REAL(dp) :: to_sort(6) = (/144.4, 11.0, 4.3, 58.6, 10.0, 7.8/)
    CALL TAG_TEST("TEST_quicksort_real_dp2_before")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_dp2_after")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real_dp_random
    REAL(dp) :: to_sort(6)
    ! Generate random numbers geq 0.0 and < 256.0
    CALL RANDOM_NUMBER(to_sort)
    to_sort = to_sort*256

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_dp_random")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_permutation_real_dp
    REAL(dp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    INTEGER :: idx_permutation(6) = (/1, 2, 3, 4, 5, 6/)
    CALL TAG_TEST("TEST_quicksort_permutation_real_before")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .FALSE.)

    CALL quicksort(to_sort, idx_permutation)

    CALL TAG_TEST("TEST_quicksort_permutation_real_dp_after")
    CALL ASSERT_EQUAL(is_sorted_real_dp(to_sort), .TRUE.)
    CALL TAG_TEST("TEST_quicksort_permutation_real_dp_permutation")
    CALL ASSERT_EQUAL(has_same_values_int(idx_permutation, (/3, 4, 5, 6, 2, 1/)), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real_sp
    REAL(sp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    CALL TAG_TEST("TEST_quicksort_real_sp_before")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_sp_after")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real_sp2
    REAL(sp) :: to_sort(6) = (/144.4, 11.0, 4.3, 58.6, 10.0, 7.8/)
    CALL TAG_TEST("TEST_quicksort_real_sp2_before")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_sp2_after")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real_sp_random
    REAL(sp) :: to_sort(6)
    ! Generate random numbers geq 0.0 and < 256.0
    CALL RANDOM_NUMBER(to_sort)
    to_sort = to_sort*256

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_sp_random")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_permutation_real_sp
    REAL(sp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    INTEGER :: idx_permutation(6) = (/1, 2, 3, 4, 5, 6/)
    CALL TAG_TEST("TEST_quicksort_permutation_real_before")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .FALSE.)

    CALL quicksort(to_sort, idx_permutation)

    CALL TAG_TEST("TEST_quicksort_permutation_real_sp_after")
    CALL ASSERT_EQUAL(is_sorted_real_sp(to_sort), .TRUE.)
    CALL TAG_TEST("TEST_quicksort_permutation_real_sp_permutation")
    CALL ASSERT_EQUAL(has_same_values_int(idx_permutation, (/3, 4, 5, 6, 2, 1/)), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    CALL TAG_TEST("TEST_quicksort_int_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_int_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int2
    INTEGER :: to_sort(6) = (/58, 4, 144, 10, 7, 11/)
    CALL TAG_TEST("TEST_quicksort_int2_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_int2_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int_random
    INTEGER :: to_sort(6)
    REAL(wp) :: random_wp(6)
    ! Generate random numbers between 0 and 255
    CALL RANDOM_NUMBER(random_wp)
    to_sort = FLOOR(random_wp*256)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_int_random")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_permutation_int
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    INTEGER :: idx_permutation(6) = (/1, 2, 3, 4, 5, 6/)
    CALL TAG_TEST("TEST_quicksort_permutation_int_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL quicksort(to_sort, idx_permutation)

    CALL TAG_TEST("TEST_quicksort_permutation_int_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
    CALL TAG_TEST("TEST_quicksort_permutation_int_permutation")
    CALL ASSERT_EQUAL(has_same_values_int(idx_permutation, (/3, 4, 5, 6, 2, 1/)), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string
    CHARACTER :: to_sort(6) = (/'A', 'C', 'Y', 'E', 'S', 'H'/)
    CALL TAG_TEST("TEST_quicksort_string_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string2
    CHARACTER :: to_sort(6) = (/'Y', 'H', 'A', 'S', 'E', 'C'/)
    CALL TAG_TEST("TEST_quicksort_string2_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string2_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string3
    CHARACTER :: to_sort(6) = (/'P', 'M', 'W', 'G', 'K', 'D'/)
    CALL TAG_TEST("TEST_quicksort_string3_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string3_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string4
    CHARACTER :: to_sort(6) = (/'B', 'L', 'Q', 'S', 'Z', 'T'/)
    CALL TAG_TEST("TEST_quicksort_string4_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string4_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_insertion_sort_int
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    CALL TAG_TEST("TEST_insertion_sort_int_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL insertion_sort(to_sort)

    CALL TAG_TEST("TEST_insertion_sort_int_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_insertion_sort_int_random
    INTEGER :: to_sort(6)
    REAL(wp) :: random_wp(6)
    ! Generate random numbers between 0 and 255
    CALL RANDOM_NUMBER(random_wp)
    to_sort = FLOOR(random_wp*256)

    CALL insertion_sort(to_sort)

    CALL TAG_TEST("TEST_insertion_sort_int_random")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  LOGICAL FUNCTION has_same_values_int(array, ref)
    INTEGER, INTENT(IN) :: array(:), ref(:)
    INTEGER :: i

    has_same_values_int = .TRUE.
    DO i = 1, SIZE(array)
      IF (array(i) /= ref(i)) THEN
        has_same_values_int = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION has_same_values_int

  LOGICAL FUNCTION is_sorted_real_dp(array)
    REAL(dp), INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_real_dp = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_real_dp = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_real_dp

  LOGICAL FUNCTION is_sorted_real_sp(array)
    REAL(sp), INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_real_sp = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_real_sp = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_real_sp

  LOGICAL FUNCTION is_sorted_int(array)
    INTEGER, INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_int = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_int = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_int

  LOGICAL FUNCTION is_sorted_string(array)
    CHARACTER, INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_string = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_string = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_string

END MODULE TEST_mo_util_sort
