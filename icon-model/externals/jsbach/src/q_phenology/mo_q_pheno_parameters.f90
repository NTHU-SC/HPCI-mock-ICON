!> QUINCY phenology parameters
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
!>#### declare and init q_phenology parameters
!>
MODULE mo_q_pheno_parameters
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  REAL(wp), SAVE :: &
    & gdd_t_air_threshold = def_parameters, &   !< temperature threshold for the accumulation of growing degree days [deg K]
    & k_leafon_canopy = def_parameters          !< scaling coefficient for leaf development rate in CANOPY-only model

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_pheno_parameters'

  !$ACC DECLARE CREATE(gdd_t_air_threshold, k_leafon_canopy)

CONTAINS

  ! ======================================================================================================= !
  !> initialize parameters for the process: phenology
  !>
  !>   routine is called in mo_jsb_base
  !>
  SUBROUTINE init_q_pheno_parameters
    USE mo_jsb_physical_constants,    ONLY: Tzero
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_q_pheno_parameters'

    gdd_t_air_threshold       = 5.0_wp + Tzero              !< citation or source ....  [ ? ]
    k_leafon_canopy           = 0.2_wp

    !$ACC UPDATE DEVICE(gdd_t_air_threshold, k_leafon_canopy) ASYNC(1)

  END SUBROUTINE init_q_pheno_parameters

#endif
END MODULE mo_q_pheno_parameters
