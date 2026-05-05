! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! This module contains the I/O routines for initicon

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_initicono

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish
  USE mo_impl_constants,          ONLY: SUCCESS, MODE_IAU_OCE, max_ntracer
  USE mo_impl_constants,          ONLY: max_ntracer
  USE mo_ocean_types,             ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,           ONLY: t_sea_ice
  USE mo_ocean_nml,               ONLY: init_mode_oce
  USE mo_initicon_types,          ONLY: t_pi_tracer
  USE mtime,                      ONLY: datetime
  USE mo_run_config,              ONLY: ntracer
  USE mo_input_instructions,      ONLY: t_readInstructionListPtr, &
  &                               kInputSourceFg, kInputSourceAna, kInputSourceBoth
  USE mo_input_request_list,      ONLY: t_InputRequestList
  USE mo_initicon_types,          ONLY: t_init_state_const
  USE mo_dynamics_config,         ONLY: nold
  USE mo_fortran_tools,           ONLY: DO_DEALLOCATE, DO_PTR_DEALLOCATE
  USE mo_initicon_io,             ONLY: fetch3d, fetch3d_with_status, fetchSurface, t_fetchParams
  IMPLICIT NONE
  PRIVATE

  !
  !variables

  PUBLIC :: t_inito_state
  PUBLIC :: t_initicono_state  !> state vector for initicon-o
  PUBLIC :: t_pi_oce_in
  PUBLIC :: t_pi_seaice_in
  PUBLIC :: t_pi_oce
  PUBLIC :: t_pi_seaice
  PUBLIC :: t_pi_seaice_inc
  PUBLIC :: t_initicono_read
  PUBLIC :: read_initicono

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_ocean_initicono'


  ! functions
  PUBLIC :: t_inito_state_finalize
  PUBLIC :: fetch_dwdfg_oce, fetch_dwdfg_seaice
  PUBLIC :: fetch_dwdana_oce, fetch_dwdana_seaice


  ! ocean input variables
  TYPE :: t_pi_oce_in

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! vertical dimension of 3D input fields
    INTEGER :: nlev

    REAL(wp), POINTER, DIMENSION(:,:,:) :: to    => NULL(), &
    &                                    so      => NULL(), &
    &                                    u       => NULL(), &
    &                                    v       => NULL(), &
    &                                    vn       => NULL(), &
    &                                    zos     => NULL(), &
    &                                    depth   => NULL()
    REAL(wp), POINTER, DIMENSION(:,:) :: stretch_c => NULL()

    TYPE (t_pi_tracer), DIMENSION(max_ntracer) :: tracer

    CONTAINS
      PROCEDURE :: finalize => t_pi_oce_in_finalize   !< destructor
  END TYPE t_pi_oce_in

  ! surface input variables
  TYPE :: t_pi_seaice_in

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! number of ice levels (should be one)
    INTEGER :: lev

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: hs, hi, conc

    CONTAINS
      PROCEDURE :: finalize => t_pi_seaice_in_finalize !< destructor
  END TYPE t_pi_seaice_in


  TYPE :: t_pi_oce

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! vertical dimension of 3D input fields
    INTEGER :: nlev

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: u, v, vn, to, so, zos, depth
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: stretch_c

    TYPE (t_pi_tracer), DIMENSION(max_ntracer) :: tracer

    CONTAINS
      PROCEDURE :: finalize => t_pi_oce_finalize   !< destructor
  END TYPE t_pi_oce

  !
  TYPE :: t_pi_seaice

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: hi, hs, conc

    CONTAINS
      PROCEDURE :: finalize => t_pi_seaice_finalize   !< destructor
  END TYPE t_pi_seaice


  ! seaice field increments
  TYPE :: t_pi_seaice_inc
    !
    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:)   :: hi
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:)   :: hs
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:)   :: conc

    CONTAINS
      PROCEDURE :: finalize => t_pi_seaice_inc_finalize   !< destructor
  END TYPE t_pi_seaice_inc

  ! state vector type: base class
  !
  TYPE :: t_inito_state

    TYPE (t_pi_oce_in)     :: oce_in
    TYPE (t_pi_seaice_in)  :: seaice_in
    TYPE (t_pi_oce)        :: oce
    TYPE (t_pi_seaice)     :: seaice

    TYPE (t_init_state_const), POINTER :: const => NULL()

    TYPE (datetime)        :: vDateTime  ! validity DateTime
    CONTAINS
      PROCEDURE, PUBLIC :: finalize => t_inito_state_finalize
  END TYPE t_inito_state

  ! Adapt the types of mo_initicon_types for the ocean case. Name it initicono.
  ! complete state vector type
  !
  TYPE, EXTENDS(t_inito_state) :: t_initicono_state

    TYPE (t_pi_oce)      :: oce_inc
    TYPE (t_pi_seaice)   :: seaice_inc

    CONTAINS
      PROCEDURE, PUBLIC :: finalize => t_initicono_state_finalize
  END TYPE t_initicono_state


  !This type contains all variables that can be read
  !with initicon-o. For now this needs to be extended
  !manually.
  TYPE :: t_initicono_read

    LOGICAL :: u = .FALSE., v = .FALSE., vn = .FALSE.
    LOGICAL :: to= .FALSE., so = .FALSE.
    LOGICAL :: zos = .FALSE., stretch_c = .FALSE.
    LOGICAL :: hi = .FALSE., hs = .FALSE., conc = .FALSE.

  END TYPE t_initicono_read

  TYPE(t_initicono_read) :: read_initicono

  CHARACTER(LEN=4) :: to_var


  CONTAINS

  SUBROUTINE t_pi_oce_in_finalize(oce_in)
    CLASS(t_pi_oce_in), INTENT(INOUT) :: oce_in

    INTEGER :: idx

    oce_in%linitialized = .FALSE.
    CALL DO_PTR_DEALLOCATE(oce_in%to)
    CALL DO_PTR_DEALLOCATE(oce_in%so)
    CALL DO_PTR_DEALLOCATE(oce_in%u)
    CALL DO_PTR_DEALLOCATE(oce_in%v)
    CALL DO_PTR_DEALLOCATE(oce_in%vn)
    CALL DO_PTR_DEALLOCATE(oce_in%zos)
    CALL DO_PTR_DEALLOCATE(oce_in%depth)
    CALL DO_PTR_DEALLOCATE(oce_in%stretch_c)

    DO idx=1, ntracer
      CALL oce_in%tracer(idx)%finalize()
    END DO
  END SUBROUTINE t_pi_oce_in_finalize

  SUBROUTINE t_pi_seaice_in_finalize(seaice_in)
    CLASS(t_pi_seaice_in), INTENT(INOUT) :: seaice_in

    seaice_in%linitialized = .FALSE.
    CALL DO_DEALLOCATE(seaice_in%hi)
    CALL DO_DEALLOCATE(seaice_in%hs)
    CALL DO_DEALLOCATE(seaice_in%conc)
  END SUBROUTINE t_pi_seaice_in_finalize


  SUBROUTINE t_pi_oce_finalize(oce)
    CLASS(t_pi_oce), INTENT(INOUT) :: oce

    INTEGER :: idx

    oce%linitialized = .FALSE.
    CALL DO_DEALLOCATE(oce%to)
    CALL DO_DEALLOCATE(oce%so)
    CALL DO_DEALLOCATE(oce%u)
    CALL DO_DEALLOCATE(oce%v)
    CALL DO_DEALLOCATE(oce%vn)
    CALL DO_DEALLOCATE(oce%zos)
    CALL DO_DEALLOCATE(oce%depth)
    CALL DO_DEALLOCATE(oce%stretch_c)

    DO idx=1, ntracer
      CALL oce%tracer(idx)%finalize()
    END DO
  END SUBROUTINE t_pi_oce_finalize


  SUBROUTINE t_pi_seaice_finalize(seaice)
    CLASS(t_pi_seaice), INTENT(INOUT) :: seaice

    seaice%linitialized = .FALSE.
    CALL DO_DEALLOCATE(seaice%hs)
    CALL DO_DEALLOCATE(seaice%hi)
    CALL DO_DEALLOCATE(seaice%conc)

  END SUBROUTINE t_pi_seaice_finalize


  SUBROUTINE t_pi_seaice_inc_finalize(seaice_inc)
    CLASS(t_pi_seaice_inc), INTENT(INOUT) :: seaice_inc

    seaice_inc%linitialized = .FALSE.
    CALL DO_DEALLOCATE(seaice_inc%hi)
    CALL DO_DEALLOCATE(seaice_inc%hs)
    CALL DO_DEALLOCATE(seaice_inc%conc)
  END SUBROUTINE t_pi_seaice_inc_finalize


  SUBROUTINE t_inito_state_finalize(inito_data)
    CLASS(t_inito_state), INTENT(INOUT) :: inito_data

    CALL inito_data%oce_in%finalize()
    CALL inito_data%seaice_in%finalize()
    CALL inito_data%oce%finalize()
    CALL inito_data%seaice%finalize()

  END SUBROUTINE t_inito_state_finalize


  SUBROUTINE t_initicono_state_finalize(inito_data)
    CLASS(t_initicono_state), INTENT(INOUT) :: inito_data

    ! call base class destructor
    CALL t_inito_state_finalize(inito_data)

    CALL inito_data%oce_inc%finalize()
    CALL inito_data%seaice_inc%finalize()
  END SUBROUTINE t_initicono_state_finalize

  !>
  !! Fetch the DWD first guess from the request list (ocean only)
  !! First guess (FG) is read for to, so, u, v, zos, depth
  !! whereas DA output is read for to, so, u, v, depth
  SUBROUTINE fetch_dwdfg_oce(requestList, ocean_state, inputInstructions, read_initicono)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_initicono_read) :: read_initicono

    CHARACTER(*), PARAMETER :: routine = modname//':fetch_dwdfg_oce'
    TYPE(t_fetchParams)      :: params
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    LOGICAL :: lfound_u, lfound_v, lfound_vn, lfound_to, lfound_so

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

      !request the first guess fields (ocean only)
      IF(read_initicono%zos) CALL fetchSurface(params, 'zos', 1, ocean_state(1)%p_prog(nold(1))%h)
      IF(read_initicono%stretch_c) CALL fetchSurface(params, 'stretch_c', 1, ocean_state(1)%p_prog(nold(1))%stretch_c)

      !The following variables are read in even though they are in
      !the diagnostic group.

      IF(read_initicono%u) CALL fetch3d_with_status(routine, 'dwdfg file', params, 'u', 1, ocean_state(1)%p_diag%u, lfound_u)
      IF(read_initicono%v) CALL fetch3d_with_status(routine, 'dwdfg file', params, 'v', 1, ocean_state(1)%p_diag%v, lfound_v)
      IF(read_initicono%vn) CALL fetch3d_with_status(routine, 'dwdfg file', params, 'vn', 1, ocean_state(1)%p_prog(nold(1))%vn, &
                                                     & lfound_vn)

      IF(read_initicono%to) THEN
        my_ptr3d => ocean_state(1)%p_diag%SWPT
        CALL fetch3d_with_status(routine, 'dwdfg file', params, 'to', 1, my_ptr3d, lfound_to)

        IF(lfound_to) THEN
          to_var = 'to'
        ELSE
          CALL fetch3d(params, 'SWPT', 1, my_ptr3d)
          to_var = 'SWPT'
        ENDIF
        !GRIB files contain temperature in Kelvin but ICON-O wants degrees Celsius, conversion is done with a post-op
        IF(MAXVAL(ocean_state(1)%p_diag%SWPT) >= 100._wp) THEN
          write(0,*) MAXVAL(ocean_state(1)%p_diag%SWPT)
          call finish(routine, "Wrong temperature unit. Please check the naming of your temperature field in your input file.")
        ENDIF
        ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,1) = ocean_state(1)%p_diag%SWPT
      ENDIF

      IF(read_initicono%so) THEN
        my_ptr3d => ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,2)
        CALL fetch3d_with_status(routine, 'dwdfg file', params, 'so', 1, my_ptr3d, lfound_so)
      ENDIF

  END SUBROUTINE fetch_dwdfg_oce

  !>
  !! Fetch DA-analysis DATA from the request list (ocean only)
  !!
  !! Depending on the initialization mode, either full fields or increments
  !! are read (atmosphere only). The following full fields are read, if available:
  !!     u, v, to, so, zos, stretch_c
  !!
  SUBROUTINE fetch_dwdana_oce(requestList, ocean_state, initicono, inputInstructions, read_initicono)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_initicono_state), INTENT(INOUT), TARGET :: initicono(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_initicono_read) :: read_initicono

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_oce'
    TYPE(t_pi_oce), POINTER :: my_ptr
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg, lfound_to, lfound_so, lfound_u, lfound_v

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    ! Depending on the initialization mode chosen (incremental vs. non-incremental)
    ! input fields are stored in different locations.
    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      my_ptr => initicono(1)%oce_inc
    ELSE
      my_ptr => initicono(1)%oce
    ENDIF

    ! start reading DA output (ocean only)
    ! The dynamical variables temp, salinity, u and v can be directly taken from the analysis and are thus written
    ! to the prognostic state
    IF ( ( init_mode_oce == MODE_IAU_OCE ) .AND. read_initicono%to) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar(to_var) == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, to_var, 1, my_ptr%to, lfound_to)
      ! check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar(to_var) == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource(to_var, kInputSourceBoth)
      END IF
    ELSEIF(read_initicono%to) THEN
      my_ptr3d => ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,1)
      CALL fetch3d_with_status(routine, 'dwdana file', params, to_var, 1, my_ptr3d, lfound_to)
    ENDIF

    IF ( ( init_mode_oce == MODE_IAU_OCE ) .AND. read_initicono%so) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('so') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'so', 1, my_ptr%so, lfound_so)
      ! check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('so') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('so', kInputSourceBoth)
      END IF
    ELSEIF(read_initicono%so) THEN
      my_ptr3d => ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,2)
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'so', 1, my_ptr3d, lfound_so)
    ENDIF


    IF ( init_mode_oce == MODE_IAU_OCE .AND. read_initicono%u) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('u') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'u', 1, my_ptr%u, lfound_u)
      !check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('u') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('u', kInputSourceBoth)
      END IF
    ELSEIF(read_initicono%u) THEN
      my_ptr3d => ocean_state(1)%p_diag%u
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'u', 1, my_ptr3d, lfound_u)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE .AND. read_initicono%v) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('v') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'v', 1, my_ptr%v, lfound_v)
      !check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('v') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('v', kInputSourceBoth)
      END IF
    ELSEIF (read_initicono%v) THEN
      my_ptr3d => ocean_state(1)%p_diag%v
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'v', 1, my_ptr3d, lfound_v)
    ENDIF

  END SUBROUTINE fetch_dwdana_oce

  !>
  !! Fetch DWD first guess DATA from request list (ocean only)
  SUBROUTINE fetch_dwdfg_seaice(requestList, p_sea_ice, inputInstructions, read_initicono)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_sea_ice
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_initicono_read) :: read_initicono
    TYPE(t_fetchParams) :: params
    LOGICAL :: lfound_hi, lfound_hs, lfound_conc
    INTEGER :: error

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_seaice'

    ! Workaround for the intel compiler botching implicit allocation.
    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    ! sea-ice related fields
    IF(read_initicono%hi) CALL fetch3d_with_status(routine, 'dwdfg file', params, 'hi', 1, p_sea_ice%hi(:,:,:), lfound_hi)
    IF(read_initicono%hs) CALL fetch3d_with_status(routine, 'dwdfg file', params, 'hs', 1, p_sea_ice%hs(:,:,:), lfound_hs)
    IF(read_initicono%conc) CALL fetch3d_with_status(routine, 'dwdfg file', params, 'conc', 1, p_sea_ice%conc(:,:,:), lfound_conc)

  END SUBROUTINE fetch_dwdfg_seaice

  !>
  !! Fetch DWD analysis DATA from the request list (seaice only)
  !!
  !! Analysis is read for:
  !! hi, hs, conc
  !!
  SUBROUTINE fetch_dwdana_seaice(requestList, p_sea_ice, initicono, inputInstructions, read_initicono)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_sea_ice), TARGET, INTENT(INOUT) :: p_sea_ice
    TYPE(t_initicono_state), INTENT(INOUT), TARGET :: initicono(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_pi_seaice), POINTER :: my_ptr
    TYPE(t_initicono_read) :: read_initicono

    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg, lfound_hi, lfound_hs, lfound_conc

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_seaice'

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      my_ptr => initicono(1)%seaice_inc
    ENDIF

    ! hi
    IF ( init_mode_oce == MODE_IAU_OCE .AND. read_initicono%hi) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('hi') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hi', 1, my_ptr%hi, lfound_hi)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg.AND.inputInstructions(1)%ptr%sourceOfVar('hi') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('hi', kInputSourceBoth)
      END IF
    ELSEIF(read_initicono%hi) THEN
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hi', 1, p_sea_ice%hi(:,:,:), lfound_hi)
    ENDIF

    ! hs
    IF ( init_mode_oce == MODE_IAU_OCE .AND. read_initicono%hs) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('hs') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hs', 1, my_ptr%hs, lfound_hs)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg.AND.inputInstructions(1)%ptr%sourceOfVar('hs') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('hs', kInputSourceBoth)
      END IF
    ELSEIF(read_initicono%hs) THEN
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hs', 1, p_sea_ice%hs(:,:,:), lfound_hs)
    ENDIF

    ! conc
    IF ( init_mode_oce == MODE_IAU_OCE .AND. read_initicono%conc) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('conc') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'conc', 1, my_ptr%conc, lfound_conc)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg.AND.inputInstructions(1)%ptr%sourceOfVar('conc') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('conc', kInputSourceBoth)
      END IF
    ELSEIF(read_initicono%conc) THEN
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'conc', 1, p_sea_ice%conc(:,:,:), lfound_conc)
    ENDIF

  END SUBROUTINE fetch_dwdana_seaice

END MODULE mo_ocean_initicono
