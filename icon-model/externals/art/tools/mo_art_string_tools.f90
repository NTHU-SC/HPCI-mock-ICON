!
! mo_art_string_tools
! This module provides routines for parsing arrays and production elements from
! XML strings
!
!
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

MODULE mo_art_string_tools
  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_key_value_store,               ONLY: t_key_value_store
  ! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_get_no_elements_in_string
  PUBLIC :: art_split_string_in_array
  PUBLIC :: art_parse_production
  PUBLIC :: key_value_storage_as_string
  PUBLIC :: print_integer
  PUBLIC :: get_number_of_digits

CONTAINS
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_get_no_elements_in_string(num_elem,string,sep,positions_sep)
!<
! SUBROUTINE art_get_no_elements_in_string
! Get number of elements if string is separated via sep
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: &
    &  string,                    &
    &  sep
  INTEGER, INTENT(out) ::             &
    &  num_elem
  INTEGER, INTENT(out), ALLOCATABLE, OPTIONAL :: &
    &  positions_sep(:)

  ! local variables
  INTEGER ::       &
    & number_of_seps_in_str, &
    & i,                     &
    & pos(50)


  IF (LEN(sep) > 1) THEN
    CALL finish('mo_art_string_tools:art_get_number_of_elements',  &
           &    'parameter sep is too long (must have length 1).')
  END IF

  number_of_seps_in_str = 0

  DO i = 1,LEN_TRIM(string)
    IF (string(i:i) == sep) THEN
      number_of_seps_in_str = number_of_seps_in_str + 1
      IF (number_of_seps_in_str <= SIZE(pos)) THEN
        pos(number_of_seps_in_str) = i
      END IF
    END IF
  END DO


  num_elem = number_of_seps_in_str + 1

  IF (PRESENT(positions_sep)) THEN
    IF (number_of_seps_in_str > 50) THEN
      CALL finish('mo_art_string_tools:art_get_number_of_elements', &
             &    'number of '//sep//' in string must not exceed 50.')
    ELSE IF (number_of_seps_in_str >= 1) THEN
      ALLOCATE(positions_sep(number_of_seps_in_str))
      positions_sep(:) = pos(1:number_of_seps_in_str)
    END IF
  END IF


END SUBROUTINE art_get_no_elements_in_string
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_split_string_in_array(str_arr,string,sep)
!<
! SUBROUTINE art_split_string_in_array
! Returns a string array of string separated via sep
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  CHARACTER(LEN=*), ALLOCATABLE, INTENT(out) :: &
    &  str_arr(:)
  CHARACTER(LEN=*), INTENT(in) :: &
    &  string,                    &
    &  sep

  ! local variables
  INTEGER :: &
    &  num_elem, i
  INTEGER, ALLOCATABLE :: &
    &  positions_sep(:)

  CALL art_get_no_elements_in_string(num_elem,TRIM(string),sep,positions_sep)


  ALLOCATE(str_arr(num_elem))

  IF (num_elem == 1) THEN
    str_arr(1) = TRIM(ADJUSTL(string))
  ELSE
    DO i = 1,num_elem
      IF (i == 1) THEN
        str_arr(i) = string(1:positions_sep(i)-1)
      ELSE IF (i < num_elem) THEN
        str_arr(i) = string(positions_sep(i-1)+1:positions_sep(i)-1)
      ELSE
        str_arr(i) = string(positions_sep(i-1)+1:)
      END IF
    END DO
  END IF
END SUBROUTINE art_split_string_in_array
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_parse_production(factors,tracer_names,string)
!<
! SUBROUTINE art_parse_production
! Returns an array of factors and tracer names from a string of format:
! 0.007*TRSO2;0.993*TROCS
!
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  REAL(wp), ALLOCATABLE, INTENT(out) :: &
    &  factors(:)
  CHARACTER(LEN=*), ALLOCATABLE, INTENT(out) :: &
    &  tracer_names(:)
  CHARACTER(LEN=*), INTENT(in) :: &
    &  string

  ! local variables
  CHARACTER(LEN = 30), ALLOCATABLE :: &
    &  str_arr(:)
  CHARACTER(LEN=30) :: &
    &  prod_str
  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &  tracer_name
  INTEGER ::      &
    &  ioerr,     &
    &  num_prods, &
    &  i,         &
    &  index_star


  CALL art_split_string_in_array(str_arr,TRIM(string),';')

  num_prods = SIZE(str_arr)
  ALLOCATE(factors(num_prods))
  ALLOCATE(tracer_names(num_prods))

  DO i = 1,SIZE(str_arr)
    prod_str = str_arr(i)
    index_star = INDEX(prod_str,'*')

    IF (index_star > 0) THEN
      READ(prod_str(1:index_star-1),*,iostat=ioerr) factors(i)

      IF (ioerr /= 0) THEN
        CALL finish('mo_art_string_tools:art_parse_production',  &
               &    'Could not read factor from '//prod_str(1:index_star-1)  &
               &  //' which is part of '//TRIM(string)//'.')
      END IF
    ELSE
      factors(i) = 1
    END IF

    tracer_names(i) = TRIM(ADJUSTL(prod_str(index_star+1:)))
    tracer_name = tracer_names(i)
    IF (tracer_name(1:2) /= 'TR') THEN
      CALL finish('mo_art_string_tools:art_parse_production',  &
             &    'parsed tracer name '//TRIM(tracer_names(i))  &
             &  //' does not begin with TR.')
    END IF
  END DO

END SUBROUTINE art_parse_production

SUBROUTINE key_value_storage_as_string(key_value_store, key, val, ierror)
  IMPLICIT NONE
  TYPE(t_key_value_store), INTENT(in) :: &
    &  key_value_store
  CHARACTER(LEN = *), INTENT(in) :: &
    &  key
  CHARACTER(:), ALLOCATABLE, INTENT(out) :: &
    &  val
  INTEGER, INTENT(out), OPTIONAL :: &
    &  ierror

  IF (PRESENT(ierror)) THEN
    CALL key_value_store%get(key,val,ierror)
  ELSE
    CALL key_value_store%get(key,val)
  END IF

END SUBROUTINE key_value_storage_as_string


FUNCTION print_integer(number, number_of_digits) RESULT(number_string)
    ! Print an integer in a string variable
    !
    ! Take an integer and convert it to a string without any practical limit but the range
    ! of the Fortran `INTEGER` type
    !
    ! Parameters
    ! ----------
    ! number : INTEGER, INTENT(IN)
    !     Integer that will be printed
    ! number_of_digits : INTEGER, INTENT(IN)
    !     Number of digits of `number` as calculated by `get_number_of_digits`
    !
    ! Returns
    ! -------
    ! number_string : CHARACTER(LEN=number_of_digits)
    !     The string equivalent of the input integer

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: number
    INTEGER, INTENT(IN) :: number_of_digits
    CHARACTER(LEN=number_of_digits) :: number_string

    WRITE(number_string, "(I0)") number

END FUNCTION print_integer


FUNCTION get_number_of_digits(number) RESULT(number_of_digits)
    ! Get the number of digits of an integer
    !
    ! This function uses the proposition that
    !
    ! .. math: k = \lfloor \log_{10}(n) \rfloor + 1
    !
    ! if :math:`n \in \mathbb{Z}^{+}` has :math:`k \in \mathbb{Z}^{+}` digits.
    ! :math:`\lfloor x\rfloor: \mathbb{R} \longrightarrow \mathbb{R}` is the floor
    ! function that results in the nearest integer to :math:`x` from below. In the case of
    ! :math:`n \in \mathbb{Z}^{-}` one uses :math:`\lvert n\rvert \in \mathbb{Z}^{+}`. The
    ! special case of zero is defined as 1, given that the logarithm functions are
    ! undefined at zero.
    !
    ! Parameters
    ! ----------
    ! number : INTEGER, INTENT(IN)
    !     Integer to which we will count the digits
    !
    ! Returns
    ! -------
    ! number_of_digits : INTEGER
    !     The number of digits of the input integer
    !
    ! Notes
    ! -----
    ! Here we show a brief proof of the formula. The basis is the following inequality
    !
    ! .. math:: 10^{k - 1} \leq n < 10^{k}
    !
    ! for a :math:`n \in \mathbb{Z}^{+}` that has :math:`k \in \mathbb{Z}^{+}` digits in
    ! base 10. This fact means that the number can be written as a sum of powers of 10.
    !
    ! .. math:: n = \sum_{i=0}^{k - 1} n_{i}\, 10^{i}
    !
    ! where :math:`n_{i} \in \{0, 1,\dots, 9\}` for :math:`i\neq k - 1` and
    ! :math:`n_{k - 1} \in \{1,\dots, 9\}` are the digits. :math:`n_{k - 1} \neq 0`: if it
    ! were, then :math:`n` would have :math:`k - 1` digits. Given that :math:`0\leq n_{i}`
    ! for all :math:`i \neq k - 1`, we can write
    !
    ! .. math::
    !
    !    0                       &\leq \sum_{i=0}^{k-2} n_{i}\, 10^{i}
    !    n_{k - 1}\, 10^{k - 1}  &\leq \sum_{i=0}^{k-1} n_{i}\, 10^{i} = n
    !    10^{k - 1}              &\leq n,\, 1\leq n_{k - 1}
    !
    ! Now, we will prove :math:`n < 10^{k}` by contradiction. Suppose that
    ! :math:`10^{k} \leq n`. We remember that :math:`n_{i} \leq 9 \forall i`. Therefore,
    !
    ! .. math::
    !
    !    n      &\leq 9\sum_{i=0}^{k-1}10^{i}
    !    10^{k} &\leq 9\sum_{i=0}^{k-1}10^{i}
    !    1      &\leq 9\sum_{i=0}^{k-1}10^{i - k}
    !           &\leq 9\sum_{j=1}^{k}\left(\dfrac{1}{10}\right)^{j}
    !           &\leq 1 - \left(\dfrac{1}{10}\right)^{k}
    !    0      &\geq \left(\dfrac{1}{10}\right)^{k}
    !
    ! which is a contradiction because the power of a positive integer is always positive.
    ! We conclude that it should be :math:`10^{k-1}\leq n < 10^{k}`. QED
    !
    ! Now, logarithms are monotonic increasing functions defined for all
    ! :math:`x \in \mathbb{R}^{+}`, therefore, they preserve inequalities in this region.
    ! If we apply the logarithm base 10 (:math:`\log_{10}`) to both sides of the former
    ! inequality, we obtain
    !
    ! .. math:: k - 1 \leq \log_{10}(n) < k
    !
    ! Given that between :math:`k-1` and :math:`k` there are no other integers, the
    ! nearest integer to the logarithm base 10 of :math:`n` from below is :math:`k-1` and
    ! :math:`k = \lfloor\log_{10}(n)\rfloor + 1`. QED
    !
    ! Note on the implementation: Because the intrinsic `LOG10` only accepts floating
    ! point numbers, we convert the absolute value of `number` to a Fortran `REAL` with
    ! working precision. `FLOOR` returns an INTEGER of default kind.

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: number
    INTEGER :: number_of_digits

    number_of_digits = 1
    IF ( number /= 0 ) THEN
        number_of_digits = FLOOR(LOG10(REAL(ABS(number), wp))) + 1
    END IF

END FUNCTION get_number_of_digits
!!
!!-------------------------------------------------------------------
!!
END MODULE mo_art_string_tools
