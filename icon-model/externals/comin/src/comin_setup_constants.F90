!> @file comin_setup_constants.F90
!! @brief ComIn utilities, containing named constants for entry points and
!!        named constants for variable state (READ, WRITE, ...).
!
!  @authors 01/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_setup_constants
  USE ISO_C_BINDING,        ONLY : C_DOUBLE, C_FLOAT

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: wp, dp, sp
  PUBLIC :: EP_SECONDARY_CONSTRUCTOR,     &
    &       EP_ATM_YAC_DEFCOMP_BEFORE,    &
    &       EP_ATM_YAC_DEFCOMP_AFTER,     &
    &       EP_ATM_YAC_SYNCDEF_BEFORE,    &
    &       EP_ATM_YAC_SYNCDEF_AFTER,     &
    &       EP_ATM_YAC_ENDDEF_BEFORE,     &
    &       EP_ATM_YAC_ENDDEF_AFTER,      &
    &       EP_ATM_INIT_FINALIZE,         &
    &       EP_ATM_TIMELOOP_BEFORE,       &
    &       EP_ATM_TIMELOOP_START,        &
    &       EP_ATM_TIMELOOP_END,          &
    &       EP_ATM_TIMELOOP_AFTER,        &
    &       EP_ATM_INTEGRATE_BEFORE,      &
    &       EP_ATM_INTEGRATE_START,       &
    &       EP_ATM_INTEGRATE_END,         &
    &       EP_ATM_INTEGRATE_AFTER,       &
    &       EP_ATM_WRITE_OUTPUT_BEFORE,   &
    &       EP_ATM_WRITE_OUTPUT_AFTER,    &
    &       EP_ATM_CHECKPOINT_BEFORE,     &
    &       EP_ATM_CHECKPOINT_AFTER,      &
    &       EP_ATM_ADVECTION_BEFORE,      &
    &       EP_ATM_ADVECTION_AFTER,       &
    &       EP_ATM_PHYSICS_BEFORE,        &
    &       EP_ATM_PHYSICS_AFTER,         &
    &       EP_ATM_NUDGING_BEFORE,        &
    &       EP_ATM_NUDGING_AFTER,         &
    &       EP_ATM_SURFACE_BEFORE,        &
    &       EP_ATM_SURFACE_AFTER,         &
    &       EP_ATM_TURBULENCE_BEFORE,     &
    &       EP_ATM_TURBULENCE_AFTER,      &
    &       EP_ATM_MICROPHYSICS_BEFORE,   &
    &       EP_ATM_MICROPHYSICS_AFTER,    &
    &       EP_ATM_CONVECTION_BEFORE,     &
    &       EP_ATM_CONVECTION_AFTER,      &
    &       EP_ATM_RADIATION_BEFORE,      &
    &       EP_ATM_RADIATION_AFTER,       &
    &       EP_ATM_RADHEAT_BEFORE,        &
    &       EP_ATM_RADHEAT_AFTER,         &
    &       EP_ATM_GWDRAG_BEFORE,         &
    &       EP_ATM_GWDRAG_AFTER,          &
    &       EP_FINISH,                    &
    &       EP_DESTRUCTOR
  PUBLIC :: COMIN_FLAG_NONE, COMIN_FLAG_READ, COMIN_FLAG_WRITE, COMIN_FLAG_SYNC_HALO, COMIN_FLAG_DEVICE
  PUBLIC :: COMIN_ZAXIS_NONE, COMIN_ZAXIS_2D, COMIN_ZAXIS_3D, COMIN_ZAXIS_3D_HALF, COMIN_ZAXIS_UNDEF
  PUBLIC :: DOMAIN_UNDEFINED, COMIN_DOMAIN_OUTSIDE_LOOP
  PUBLIC :: COMIN_HGRID_UNSTRUCTURED_CELL, COMIN_HGRID_UNSTRUCTURED_EDGE, COMIN_HGRID_UNSTRUCTURED_VERTEX
  PUBLIC :: EP_NAME
  PUBLIC :: COMIN_METADATA_TYPEID_UNDEFINED, &
    &       COMIN_METADATA_TYPEID_INTEGER,   &
    &       COMIN_METADATA_TYPEID_REAL,      &
    &       COMIN_METADATA_TYPEID_CHARACTER, &
    &       COMIN_METADATA_TYPEID_LOGICAL
  PUBLIC :: COMIN_DIM_SEMANTICS_UNDEF, &
    &       COMIN_DIM_SEMANTICS_NPROMA, &
    &       COMIN_DIM_SEMANTICS_BLOCK, &
    &       COMIN_DIM_SEMANTICS_UNBLOCK, &
    &       COMIN_DIM_SEMANTICS_LEVEL, &
    &       COMIN_DIM_SEMANTICS_CONTAINER, &
    &       COMIN_DIM_SEMANTICS_OTHER, &
    &       COMIN_DIM_SEMANTICS_UNUSED
  PUBLIC :: COMIN_VAR_DATATYPE_DOUBLE, &
    &       COMIN_VAR_DATATYPE_FLOAT,  &
    &       COMIN_VAR_DATATYPE_INT

#include "comin_global.inc"

  !> List of entry points, named constants and accessor functions that
  !> are exposed to both, the host interface and the plugin interface.
  !!
  !! @defgroup constants Constants
  !! @defgroup common Common
  !! @{
  !! @}

  !> working precision
  !! @ingroup common
  INTEGER, PARAMETER   :: wp = C_DOUBLE
  INTEGER, PARAMETER   :: dp = C_DOUBLE
  INTEGER, PARAMETER   :: sp = C_FLOAT

  !> id of current domain, two states possible if not in domain loop
  !! @ingroup common
  INTEGER, PARAMETER   :: DOMAIN_UNDEFINED    = -2

  !> id of current domain, two states possible if not in domain loop
  !! @ingroup common
  INTEGER, PARAMETER   :: COMIN_DOMAIN_OUTSIDE_LOOP = -1

  !> List of entry points
  !!
  !! @ingroup common
  !! Note: EP_DESTRUCTOR should always be the last entry
  ENUM, BIND(C)
    ENUMERATOR :: EP_SECONDARY_CONSTRUCTOR = 1, &
      &           EP_ATM_YAC_DEFCOMP_BEFORE,    &
      &           EP_ATM_YAC_DEFCOMP_AFTER,     &
      &           EP_ATM_YAC_SYNCDEF_BEFORE,    &
      &           EP_ATM_YAC_SYNCDEF_AFTER,     &
      &           EP_ATM_YAC_ENDDEF_BEFORE,     &
      &           EP_ATM_YAC_ENDDEF_AFTER,      &
      &           EP_ATM_INIT_FINALIZE,         &
      &           EP_ATM_TIMELOOP_BEFORE,       &
      &           EP_ATM_TIMELOOP_START,        &
      &           EP_ATM_TIMELOOP_END,          &
      &           EP_ATM_TIMELOOP_AFTER,        &
      &           EP_ATM_INTEGRATE_BEFORE,      &
      &           EP_ATM_INTEGRATE_START,       &
      &           EP_ATM_INTEGRATE_END,         &
      &           EP_ATM_INTEGRATE_AFTER,       &
      &           EP_ATM_WRITE_OUTPUT_BEFORE,   &
      &           EP_ATM_WRITE_OUTPUT_AFTER,    &
      &           EP_ATM_CHECKPOINT_BEFORE,     &
      &           EP_ATM_CHECKPOINT_AFTER,      &
      &           EP_ATM_ADVECTION_BEFORE,      &
      &           EP_ATM_ADVECTION_AFTER,       &
      &           EP_ATM_PHYSICS_BEFORE,        &
      &           EP_ATM_PHYSICS_AFTER,         &
      &           EP_ATM_NUDGING_BEFORE,        &
      &           EP_ATM_NUDGING_AFTER,         &
      &           EP_ATM_SURFACE_BEFORE,        &
      &           EP_ATM_SURFACE_AFTER,         &
      &           EP_ATM_TURBULENCE_BEFORE,     &
      &           EP_ATM_TURBULENCE_AFTER,      &
      &           EP_ATM_MICROPHYSICS_BEFORE,   &
      &           EP_ATM_MICROPHYSICS_AFTER,    &
      &           EP_ATM_CONVECTION_BEFORE,     &
      &           EP_ATM_CONVECTION_AFTER,      &
      &           EP_ATM_RADIATION_BEFORE,      &
      &           EP_ATM_RADIATION_AFTER,       &
      &           EP_ATM_RADHEAT_BEFORE,        &
      &           EP_ATM_RADHEAT_AFTER,         &
      &           EP_ATM_GWDRAG_BEFORE,         &
      &           EP_ATM_GWDRAG_AFTER,          &
      &           EP_FINISH,                    &
      &           EP_DESTRUCTOR
  END ENUM

  !> Entry point names (character strings)
   !! @ingroup common
  CHARACTER(LEN=COMIN_MAX_LEN_EP_NAME), PARAMETER :: EP_NAME(EP_DESTRUCTOR) = [ &
    &  "EP_SECONDARY_CONSTRUCTOR    ",  &
    &  "EP_ATM_YAC_DEFCOMP_BEFORE   ",  &
    &  "EP_ATM_YAC_DEFCOMP_AFTER    ",  &
    &  "EP_ATM_YAC_SYNCDEF_BEFORE   ",  &
    &  "EP_ATM_YAC_SYNCDEF_AFTER    ",  &
    &  "EP_ATM_YAC_ENDDEF_BEFORE    ",  &
    &  "EP_ATM_YAC_ENDDEF_AFTER     ",  &
    &  "EP_ATM_INIT_FINALIZE        ",  &
    &  "EP_ATM_TIMELOOP_BEFORE      ",  &
    &  "EP_ATM_TIMELOOP_START       ",  &
    &  "EP_ATM_TIMELOOP_END         ",  &
    &  "EP_ATM_TIMELOOP_AFTER       ",  &
    &  "EP_ATM_INTEGRATE_BEFORE     ",  &
    &  "EP_ATM_INTEGRATE_START      ",  &
    &  "EP_ATM_INTEGRATE_END        ",  &
    &  "EP_ATM_INTEGRATE_AFTER      ",  &
    &  "EP_ATM_WRITE_OUTPUT_BEFORE  ",  &
    &  "EP_ATM_WRITE_OUTPUT_AFTER   ",  &
    &  "EP_ATM_CHECKPOINT_BEFORE    ",  &
    &  "EP_ATM_CHECKPOINT_AFTER     ",  &
    &  "EP_ATM_ADVECTION_BEFORE     ",  &
    &  "EP_ATM_ADVECTION_AFTER      ",  &
    &  "EP_ATM_PHYSICS_BEFORE       ",  &
    &  "EP_ATM_PHYSICS_AFTER        ",  &
    &  "EP_ATM_NUDGING_BEFORE       ",  &
    &  "EP_ATM_NUDGING_AFTER        ",  &
    &  "EP_ATM_SURFACE_BEFORE       ",  &
    &  "EP_ATM_SURFACE_AFTER        ",  &
    &  "EP_ATM_TURBULENCE_BEFORE    ",  &
    &  "EP_ATM_TURBULENCE_AFTER     ",  &
    &  "EP_ATM_MICROPHYSICS_BEFORE  ",  &
    &  "EP_ATM_MICROPHYSICS_AFTER   ",  &
    &  "EP_ATM_CONVECTION_BEFORE    ",  &
    &  "EP_ATM_CONVECTION_AFTER     ",  &
    &  "EP_ATM_RADIATION_BEFORE     ",  &
    &  "EP_ATM_RADIATION_AFTER      ",  &
    &  "EP_ATM_RADHEAT_BEFORE       ",  &
    &  "EP_ATM_RADHEAT_AFTER        ",  &
    &  "EP_ATM_GWDRAG_BEFORE        ",  &
    &  "EP_ATM_GWDRAG_AFTER         ",  &
    &  "EP_FINISH                   ",  &
    &  "EP_DESTRUCTOR               " ]

  !> Variable access flags.
  !! @ingroup plugin_interface
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_FLAG_NONE         = 0,          &
      &           COMIN_FLAG_READ         = IBSET(0,1), &
      &           COMIN_FLAG_WRITE        = IBSET(0,2), &
      &           COMIN_FLAG_SYNC_HALO    = IBSET(0,3), &
      &           COMIN_FLAG_DEVICE       = IBSET(0,4)
  END ENUM

  !> integer constant, which gives an interpretation of the horizontal
  !> grid location (cell, edge, vertex).
  !! @ingroup plugin_interface
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_HGRID_UNSTRUCTURED_CELL   = 1,  &
     &            COMIN_HGRID_UNSTRUCTURED_EDGE   = 2,  &
     &            COMIN_HGRID_UNSTRUCTURED_VERTEX = 3
  END ENUM

  !> Integer constants, giving an interpretation of the vertical axis
  !> (2D, atmospheric levels, ...)
  !! @ingroup plugin_interface
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_ZAXIS_UNDEF         = -1,   &
      &           COMIN_ZAXIS_NONE          =  0,   &
      &           COMIN_ZAXIS_2D            =  1,   &
      &           COMIN_ZAXIS_3D            =  2,   &
      &           COMIN_ZAXIS_3D_HALF       =  3
  END ENUM

  ENUM, BIND(C)
    ENUMERATOR :: COMIN_METADATA_TYPEID_UNDEFINED = -1, &
         &        COMIN_METADATA_TYPEID_INTEGER   =  0, &
         &        COMIN_METADATA_TYPEID_REAL      =  1, &
         &        COMIN_METADATA_TYPEID_CHARACTER =  2, &
         &        COMIN_METADATA_TYPEID_LOGICAL   =  3
  END ENUM

  !> Integer constants, giving an interpretation of variable data types
  !! @ingroup plugin_interface
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_VAR_DATATYPE_DOUBLE = 1,   &
       &          COMIN_VAR_DATATYPE_FLOAT,        &
       &          COMIN_VAR_DATATYPE_INT
  END ENUM

  !> Dimension semantics
  !! @ingroup plugin_interface
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_DIM_SEMANTICS_UNDEF           =  1, & ! not explicitly defined
       &          COMIN_DIM_SEMANTICS_NPROMA    =  2, & ! nproma
       &          COMIN_DIM_SEMANTICS_BLOCK     =  3, & ! blocked data layout
       &          COMIN_DIM_SEMANTICS_UNBLOCK   =  4, & ! unblocked data layout
       &          COMIN_DIM_SEMANTICS_LEVEL     =  5, & ! vertical axis dimension
       &          COMIN_DIM_SEMANTICS_CONTAINER =  6, & ! container
       &          COMIN_DIM_SEMANTICS_OTHER     =  7, & ! special cases (tracers, tiles, etc)
       &          COMIN_DIM_SEMANTICS_UNUSED    =  8    ! unused dimension
  END ENUM

END MODULE comin_setup_constants
