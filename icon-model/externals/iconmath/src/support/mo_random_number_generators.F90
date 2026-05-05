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

MODULE mo_random_number_generators

  USE mo_iconlib_kind, ONLY: wp, dp, i8, i4

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: random_normal_values, clip, initialize_random_number_generator, generate_uniform_random_number

  INTERFACE generate_uniform_random_number
    MODULE PROCEDURE generate_uniform_random_number_wp
  END INTERFACE generate_uniform_random_number

  INTERFACE random_normal_values
    MODULE PROCEDURE random_normal_values_wp
  END INTERFACE random_normal_values

  INTERFACE clip
    MODULE PROCEDURE clip_wp
  END INTERFACE clip

  ! integer parameters for using revised minstd parameters as in
  ! https://en.wikipedia.org/wiki/Lehmer_random_number_generator#Parameters_in_common_use
  INTEGER(KIND=i8), PARAMETER :: minstd_A0 = 16807
  INTEGER(KIND=i8), PARAMETER :: minstd_A = 48271

  INTEGER(KIND=i8), PARAMETER :: minstd_M = 2147483647
  REAL(KIND=dp), PARAMETER   :: minstd_Scale = 1.0_dp/REAL(minstd_M, dp)

CONTAINS

  ! initializes Lehmer random number generator
  PURE FUNCTION initialize_random_number_generator(iseed, jseed) RESULT(minstd_state)

    INTEGER(KIND=i4), INTENT(IN)      :: iseed
    INTEGER(KIND=i4), INTENT(IN)      :: jseed

    INTEGER(KIND=i4) :: minstd_state

    ! doing kind of an xorshift32 see https://en.wikipedia.org/wiki/Xorshift
    minstd_state = iseed
    minstd_state = IEOR(minstd_state, jseed*2**13)
    minstd_state = IEOR(minstd_state, jseed/2**17)
    minstd_state = IEOR(minstd_state, jseed*2**5)

    minstd_state = MOD(minstd_A0*minstd_state, minstd_M)

    minstd_state = MOD(minstd_A*minstd_state, minstd_M)

  END FUNCTION initialize_random_number_generator

  ! generated Lehmer random number generator https://en.wikipedia.org/wiki/Lehmer_random_number_generator
  FUNCTION generate_uniform_random_number_wp(minstd_state) RESULT(random_number)

    INTEGER(KIND=i4), INTENT(INOUT) :: minstd_state
    REAL(KIND=wp)                :: random_number

    minstd_state = MOD(minstd_A*minstd_state, minstd_M)

    random_number = ABS(minstd_Scale*minstd_state)

  END FUNCTION generate_uniform_random_number_wp

  SUBROUTINE clip_wp(values_range, values)

    ! Subroutine arguments (in/out/inout)
    REAL(wp), INTENT(IN)                           :: values_range ! allowed range of random numbers
    REAL(wp), DIMENSION(:), INTENT(INOUT)          :: values ! array of value to be filled with random numbers

    ! Limit random number range
    WHERE (values(:) > values_range)
      values(:) = values_range
    ELSEWHERE(values(:) < -values_range)
      values(:) = -values_range
    END WHERE

  END SUBROUTINE clip_wp

  SUBROUTINE random_normal_values_wp(seed, values_range, values, mean, stdv)

    ! Subroutine arguments (in/out/inout)
    INTEGER(KIND=i4), INTENT(IN)                  :: seed ! seed to be used
    REAL(wp), INTENT(IN)                          :: values_range ! allowed range of random numbers
    REAL(wp), DIMENSION(:), INTENT(INOUT)         :: values ! array of value to be filled with random numbers
    REAL(wp), INTENT(IN), OPTIONAL                :: mean
    REAL(wp), INTENT(IN), OPTIONAL                :: stdv

    REAL(wp), PARAMETER :: pi = 3.141592653589793238_wp ! replace with math constants from math-support once available

    ! Local variables
    INTEGER(KIND=i4)                         :: i, rng_state, n
    REAL(KIND=wp)                            :: u1, u2, z0, z1, magnitude, theta, mu, sigma

    IF (PRESENT(mean)) THEN
      mu = mean
    ELSE
      mu = 0.0_wp
    END IF

    IF (PRESENT(stdv)) THEN
      sigma = stdv
    ELSE
      sigma = 1.0_wp
    END IF

    n = SIZE(values)

    ! -----------------------------------------------------------------------
    ! Initialize a random number generator, by using the MINSTD linear
    ! congruential generator (LCG). "nmaxstreams" indicates that random
    ! numbers will be requested in blocks of this length. The generator
    ! is seeded with "seed".
    !-------------------------------------------------------------------------

    DO i = 1, n/2 + 1
      rng_state = initialize_random_number_generator(seed, i)

      u1 = generate_uniform_random_number(rng_state)
      DO WHILE (u1 <= EPSILON(0.0_wp))
        u1 = generate_uniform_random_number(rng_state)
      END DO

      u2 = generate_uniform_random_number(rng_state)

      ! using the Box-Muller transform https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
      magnitude = sigma*SQRT(-2.0_wp*LOG(u1)); 
      theta = 2.0_wp*pi*u2

      z0 = magnitude*COS(theta) + mu; 
      z1 = magnitude*SIN(theta) + mu; 
      values(2*i - 1) = z0
      IF (2*i <= n) values(2*i) = z1

    END DO

    CALL clip(values_range, values)

  END SUBROUTINE random_normal_values_wp

END MODULE mo_random_number_generators
