!> Contains constants for the turbulence processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_turb_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! vegetation parameters for turbulence calculations
  REAL(wp), SAVE :: &
    & veg_height_to_rough_momentum, &   !< factor for vegetation height to momentum [unitless]
    & veg_min_roughness                 !< minimum vegetation roughness, similar to 'lctlib_MinVegRoughness' of jsbach


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_turb_constants'

  !$ACC DECLARE CREATE(veg_height_to_rough_momentum, veg_min_roughness)

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used of the turbulence process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_turb_constants
    ! ---------------------------
    ! 0.1 InOut
    !
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_turb_constants'

    ! turbulence calculations
    veg_height_to_rough_momentum  = 10.0_wp           !< expert knowledge
    veg_min_roughness             = 0.1_wp            !< expert knowledge


    !$ACC UPDATE DEVICE(veg_height_to_rough_momentum, veg_min_roughness)
  END SUBROUTINE init_turb_constants

#endif
END MODULE mo_turb_constants
