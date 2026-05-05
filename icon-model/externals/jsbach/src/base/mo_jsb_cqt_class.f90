!> Contains conserved quantity types, definitions and methods for to-be-conserved quantities
!>
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
!>#### Contains conserved quantity types (CQTs) and the type for cqt collection (used on tiles)
!>
!> NOTE: 2D and 3D variables cannot be of the same CQT
!>
MODULE mo_jsb_cqt_class
#ifndef __NO_JSBACH__

  USE mo_util,                ONLY: int2string
  USE mo_exception,           ONLY: finish
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WATER_CQ_TYPE
  PUBLIC :: FLUX_C_CQ_TYPE, LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE
  PUBLIC :: IQ_2L2D_POOL_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE, IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE
  PUBLIC :: Get_number_of_types, Get_cqt_name, t_jsb_consQuan, t_jsb_consQuan_p, max_cqt_name_length

  ENUM, BIND(C)
    ENUMERATOR ::                 &
      & WATER_CQ_TYPE = 0,        &
      & LIVE_CARBON_CQ_TYPE,      &
      & AG_DEAD_C_CQ_TYPE,        &
      & BG_DEAD_C_CQ_TYPE,        &
      & PRODUCT_CARBON_CQ_TYPE,   &
      & FLUX_C_CQ_TYPE,           &
      & IQ_2L2D_POOL_CQ_TYPE,     & !< conserved quincy quantity type for two layer 2D pool type
      & IQ_1L2D_POOL_CQ_TYPE,     & !< conserved quincy quantity type for one layer 2D pool type
      & IQ_FLUX_CQ_TYPE,          & !< conserved quincy quantity type for fluxes
      & IQ_SL_POOL_CQ_TYPE,       & !< conserved quincy quantity of soil layered pool type
      & IQ_SL_FLUX_CQ_TYPE,       & !< conserved quincy quantity of soil layered flux type
      & LAST_CQ_TYPE ! needs to always be the last -- it is only used to determine the max number of types
  END ENUM

  !> Type used on tiles to collect conserved quantities from the process memories of the tile
  !
  TYPE :: t_jsb_consQuan
    INTEGER :: type_id = -1           !< one of the CQ_TYPEs found in mo_jsb_cqt_class
    INTEGER :: var_type = 0           !< var type (supported: REAL2D or REAL3D) -- all cq vars of one cqt should have the same var type
    LOGICAL :: transfer_all = .FALSE. !< Flag indicating if all matter should be collected (true) or if matter
                                      !< should only be collected proportionally to the lost area (false, default)
                                      !< -- all cq vars of one cqt should have the same transfer_all setting
    INTEGER :: no_of_vars = 0         !< number of vars - required for allocation of cq_vars
    INTEGER :: last_index_used = 0    !< last index used - required for filling of cq_vars
    INTEGER, ALLOCATABLE :: associated_process(:) !< process to which each CQ in cq_vars belongs
    TYPE(t_jsb_var_p), ALLOCATABLE :: cq_vars(:)
        !< collection of all variables of a certain CQT potentially from different process memories
  END TYPE t_jsb_consQuan
  TYPE :: t_jsb_consQuan_p
    TYPE(t_jsb_consQuan), POINTER :: p => NULL()
  END TYPE t_jsb_consQuan_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_cqt_class'

  INTEGER, PARAMETER :: max_cqt_name_length = 25

CONTAINS

  ! ====================================================================================================== !
  !
  !> Returns the number of CQTs defined in the enumerator
  !
  FUNCTION Get_number_of_types() RESULT(last_type_id)
    INTEGER :: last_type_id

    last_type_id = LAST_CQ_TYPE

  END FUNCTION Get_number_of_types

  ! ====================================================================================================== !
  !
  !> Returns the name for this id - restricted length (max_cqt_name_length)
  !
  FUNCTION Get_cqt_name(id) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)  :: id
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_cqt_name'
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id)
      CASE (WATER_CQ_TYPE)
        return_value = 'water'
      CASE (LIVE_CARBON_CQ_TYPE)
        return_value = 'living_carbon'
      CASE (AG_DEAD_C_CQ_TYPE)
        return_value = 'ag_dead_carbon'
      CASE (BG_DEAD_C_CQ_TYPE)
        return_value = 'bg_dead_carbon'
      CASE (PRODUCT_CARBON_CQ_TYPE)
        return_value = 'product_carbon'
      CASE (FLUX_C_CQ_TYPE)
        return_value = 'carbon_flux'
      CASE (IQ_2L2D_POOL_CQ_TYPE)
        return_value = 'quincy_bgcm_2L2D_pool'
      CASE (IQ_1L2D_POOL_CQ_TYPE)
        return_value = 'quincy_bgcm_1L2D_pool'
      CASE (IQ_FLUX_CQ_TYPE)
        return_value = 'quincy_bgcm_flux'
      CASE (IQ_SL_POOL_CQ_TYPE)
        return_value = 'quincy_sl_bgcm_pool'
      CASE (IQ_SL_FLUX_CQ_TYPE)
        return_value = 'quincy_sl_bgcm_flux'
      CASE DEFAULT
        CALL finish(TRIM(routine), 'No name specified for cq type of id '//int2string(id))
    END SELECT

  END FUNCTION Get_cqt_name

#endif
END MODULE mo_jsb_cqt_class
