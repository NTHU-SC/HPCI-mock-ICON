!> QUINCY sylviculture constants
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and define sylviculture constants
!>
MODULE mo_q_syl_constants
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  ! for now generic parameters for sylviculture model - these might get spatially and temporally varrying values later
  REAL(wp), SAVE :: &
    & fract_wood_to_pp_fuel       = def_parameters, &  !< Fraction of harvested wood allocated to the fuel product pool (unitless)
    & fract_wood_to_pp_paper      = def_parameters, &  !< Fraction of harvested wood allocated to the paper product pool (unitless)
    & fract_wood_to_pp_fiberboard = def_parameters, &  !< Fraction of harvested wood allocated to the fiberboard product pool (unitless)
    & fract_wood_to_pp_oirw       = def_parameters, &  !< Fraction of harvested wood allocated to the other industrial roundwood pp (unitless)
    & fract_wood_to_pp_pv         = def_parameters, &  !< Fraction of harvested wood allocated to the plywood and veneer product pool (unitless)
    & fract_wood_to_pp_sawnwood   = def_parameters     !< Fraction of harvested wood allocated to the sawnwood product pool (unitless)

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_syl_constants'

CONTAINS

  ! ======================================================================================================= !
  !> initialize parameters for the process: sylviculture
  !>
  !>   routine is called in mo_jsb_base
  !>
  SUBROUTINE init_q_syl_constants
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_q_syl_constants'

    ! Note: allocation fractions to product pools need to add up to one (asserted in update_pools_on_harvest)
    ! Below some European longterm averages inspired by numbers from Nuetzel (2021, LMU master thesis)
    ! possibly to be read from spatially and temporally varrying maps later
    fract_wood_to_pp_fuel       = 0.4_wp  !< European longterm average inspired by numbers from Nuetzel (2021, LMU master thesis)
    fract_wood_to_pp_paper      = 0.25_wp !< European longterm average inspired by numbers from Nuetzel (2021, LMU master thesis)
    fract_wood_to_pp_fiberboard = 0.05_wp !< European longterm average inspired by numbers from Nuetzel (2021, LMU master thesis)
    fract_wood_to_pp_oirw       = 0.05_wp !< European longterm average inspired by numbers from Nuetzel (2021, LMU master thesis)
    fract_wood_to_pp_pv         = 0.05_wp !< European longterm average inspired by numbers from Nuetzel (2021, LMU master thesis)
    fract_wood_to_pp_sawnwood   = 0.2_wp  !< European longterm average inspired by numbers from Nuetzel (2021, LMU master thesis)

  END SUBROUTINE init_q_syl_constants

#endif

END MODULE mo_q_syl_constants
