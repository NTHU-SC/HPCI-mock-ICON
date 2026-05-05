! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_dummy_coupling_raw.inc
!! \test
!! This test checks the raw data exchange feature.

#ifndef TEST_PRECISION
#error "TEST_PRECISION is not defined"
#endif

#ifndef YAC_PTR_TYPE
#error "YAC_PTR_TYPE is not defined"
#endif

! ------------------------------------------------------------------------

#include "test_macros.inc"

! ------------------------------------------------------------------------

  USE mpi
  USE, INTRINSIC :: iso_c_binding
  USE yac
  USE utest

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = TEST_PRECISION !< selected working precision

  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: src_comp_name = "source_component"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: tgt_comp_name = "target_component"

  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: src_grid_name = "source_grid"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: tgt_grid_name = "target_grid"
  INTEGER, PARAMETER :: cyclic(2) = (/0, 0/)
  INTEGER, PARAMETER :: src_nbr_vertices(2) = (/2,3/)
  INTEGER, PARAMETER :: tgt_nbr_vertices(2) = (/4,3/)
  REAL(KIND=wp), PARAMETER :: RAD = 0.01745329251994329576923690768489_wp
  REAL(KIND=wp), PARAMETER :: FRAC_MASK_TOL = 1e-12_wp
  REAL(KIND=wp), PARAMETER :: src_x_vertices(2,2) = &
    RESHAPE((/-1.0_wp,0.0_wp,0.0_wp,1.0_wp/),(/2,2/)) * RAD
  REAL(KIND=wp), PARAMETER :: src_y_vertices(3) = &
    (/-1.0_wp,0.0_wp,1.0_wp/) * RAD
  REAL(KIND=wp), PARAMETER :: tgt_x_vertices(4) = &
    (/-1.5_wp,-0.5_wp,0.5_wp,1.5_wp/) * RAD
  REAL(KIND=wp), PARAMETER :: tgt_y_vertices(3,2) = &
    RESHAPE((/-1.5_wp,-0.5_wp,0.5_wp,-0.5_wp,0.5_wp,1.5_wp/),(/3,2/)) * RAD
  LOGICAL, PARAMETER :: src_cell_is_valid(2,2) = &
    RESHAPE((/.FALSE., .TRUE., .TRUE., .TRUE./),(/2,2/))

  INTEGER, PARAMETER :: src_num_cells(2) = (/1,2/)
  INTEGER, PARAMETER :: tgt_num_cells(2) = (/3,2/)
  REAL(KIND=wp), PARAMETER :: src_x_cells(1,2) = &
    RESHAPE((/-0.5_wp,0.5_wp/),(/1,2/)) * RAD
  REAL(KIND=wp), PARAMETER :: src_y_cells(2) = (/-0.5_wp,0.5_wp/) * RAD
  REAL(KIND=wp), PARAMETER :: tgt_x_cells(3) = (/-1.0_wp,0.0_wp,1.0_wp/) * RAD
  REAL(KIND=wp), PARAMETER :: tgt_y_cells(2,2) = &
    RESHAPE((/-1.0_wp,0.0_wp,0.0_wp,1.0_wp/),(/2,2/)) * RAD

  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: src_field_name = "source_field"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: tgt_field_name = "target_field"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: src_field_csr_name = "source_field_csr"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: tgt_field_csr_name = "target_field_csr"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: src_field_frac_name = "source_field_frac"
  CHARACTER(LEN=YAC_MAX_CHARLEN), PARAMETER :: tgt_field_frac_name = "target_field_frac"

  TYPE t_interp_weights_data
    DOUBLE PRECISION :: frac_mask_fallback_value
    DOUBLE PRECISION :: scaling_factor
    DOUBLE PRECISION :: scaling_summand
    INTEGER :: num_fixed_values
    DOUBLE PRECISION, ALLOCATABLE :: fixed_values(:)
    INTEGER, ALLOCATABLE :: num_tgt_per_fixed_value(:)
    INTEGER, ALLOCATABLE :: tgt_idx_fixed(:)
    INTEGER :: num_wgt_tgt
    INTEGER, ALLOCATABLE :: wgt_tgt_idx(:)
    INTEGER, ALLOCATABLE :: num_src_per_tgt(:)
    INTEGER, ALLOCATABLE :: src_indptr(:)
    DOUBLE PRECISION, ALLOCATABLE :: weights(:)
    INTEGER, ALLOCATABLE :: src_field_idx(:)
    INTEGER, ALLOCATABLE :: src_idx(:)
    INTEGER :: num_src_fields
    INTEGER, ALLOCATABLE :: src_field_buffer_sizes(:)
  END TYPE t_interp_weights_data

  INTEGER :: src_info, tgt_info, ierror, npes, rank

  LOGICAL :: is_src, is_tgt;
  INTEGER :: src_rank, tgt_rank
  INTEGER :: comp_ids(2), src_comp_id, tgt_comp_id, num_comps
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: comp_names(2)
  INTEGER :: src_grid_id, tgt_grid_id
  INTEGER :: src_cell_mask_id
  INTEGER :: src_cell_point_id, tgt_cell_point_id
  INTEGER :: src_field_id, tgt_field_id
  INTEGER :: src_field_csr_id, tgt_field_csr_id
  INTEGER :: src_field_frac_id, tgt_field_frac_id
  INTEGER :: interp_stack_config

  INTEGER :: collection_size
  INTEGER, PARAMETER :: MAX_COLLECTION_SIZE = 1
  INTEGER, PARAMETER :: NUM_TGT_CELLS = 6
  INTEGER, PARAMETER :: NUM_TESTS_FRAC = 2
  INTEGER, PARAMETER :: NUM_SRC_FIELDS = 1
  REAL(KIND=wp), PARAMETER :: tol = 1e-2_wp
  DOUBLE PRECISION, PARAMETER :: fixed_value = -2.0
  DOUBLE PRECISION, PARAMETER :: frac_mask_fallback_value = 3.0

  TYPE (t_interp_weights_data) :: interp_weights_data, &
                                  interp_weights_data_csr, &
                                  interp_weights_data_frac

  REAL(KIND=wp), ALLOCATABLE, TARGET :: src_field(:,:,:)
  REAL(KIND=wp), ALLOCATABLE, TARGET :: src_frac_mask(:,:,:,:)
  REAL(KIND=wp), ALLOCATABLE, TARGET :: src_field_buffer(:,:)
  REAL(KIND=wp), ALLOCATABLE, TARGET :: src_frac_mask_buffer(:,:)
  TYPE(YAC_PTR_TYPE), ALLOCATABLE :: src_field_ptr(:,:)
  TYPE(YAC_PTR_TYPE), ALLOCATABLE :: src_frac_mask_ptr(:,:,:)
  TYPE(YAC_PTR_TYPE), ALLOCATABLE :: src_field_buffer_ptr(:,:)
  TYPE(YAC_PTR_TYPE), ALLOCATABLE :: src_frac_mask_buffer_ptr(:,:)
  REAL(KIND=wp), ALLOCATABLE :: tgt_field(:,:)
  REAL(KIND=wp), ALLOCATABLE :: ref_tgt_field(:,:)
  REAL(KIND=wp), ALLOCATABLE :: ref_tgt_field_frac(:,:,:)
  INTEGER :: src_field_buffer_size
  INTEGER :: num_cells
  INTEGER :: i, t

  CALL start_test("dummy_coupling_raw")

  ! initialise YAC
  CALL yac_finit()
  CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  CALL yac_fdef_datetime("2000-01-01T00:00:00", "2000-01-02T00:00:00")

  CALL mpi_comm_rank ( MPI_COMM_WORLD, rank, ierror )
  CALL mpi_comm_size ( MPI_COMM_WORLD, npes, ierror )

  IF (npes /= 4) THEN
    WRITE ( * , * ) "Invalid number of processes"
    CALL test ( .FALSE. )
    CALL stop_test
    CALL exit_tests
    CALL mpi_abort (MPI_COMM_WORLD, 999, ierror)
  END IF

  is_src = rank <= 1
  is_tgt = (rank >= 1) .AND. (rank < 3)
  src_rank = rank
  tgt_rank = rank - 1

  ! define local components
  num_comps = 0
  IF (is_src) THEN
    num_comps = num_comps + 1
    comp_names(num_comps) = src_comp_name
  END IF
  IF (is_tgt) THEN
    num_comps = num_comps + 1
    comp_names(num_comps) = tgt_comp_name
  END IF
  IF (num_comps > 0) THEN
    CALL yac_fdef_comps( &
      comp_names(1:num_comps), num_comps, comp_ids(1:num_comps))
  ELSE
    CALL yac_fdef_comp_dummy()
  END IF
  IF (is_tgt) THEN
    tgt_comp_id = comp_ids(num_comps)
    num_comps = num_comps - 1
  ELSE
    tgt_comp_id = -1
  END IF
  IF (is_src) THEN
    src_comp_id = comp_ids(num_comps)
    num_comps = num_comps - 1
  ELSE
    src_comp_id = -1
  END IF

  ! define local grids
  IF (is_src) THEN
    CALL yac_fdef_grid( &
      src_grid_name, src_nbr_vertices, cyclic, &
      src_x_vertices(:,src_rank+1), src_y_vertices, src_grid_id)
    CALL yac_fdef_mask( &
      src_grid_id, (src_nbr_vertices(1) - 1) * (src_nbr_vertices(2) - 1), &
      YAC_LOCATION_CELL, src_cell_is_valid(:, src_rank+1), src_cell_mask_id)
  ELSE
    src_grid_id = -1
    src_cell_mask_id = -1
  END IF
  IF (is_tgt) THEN
    CALL yac_fdef_grid( &
      tgt_grid_name, tgt_nbr_vertices, cyclic, &
      tgt_x_vertices, tgt_y_vertices(:,tgt_rank+1), tgt_grid_id)
  ELSE
    tgt_grid_id = -1
  END IF

  ! define cell points
  IF (is_src) THEN
    CALL yac_fdef_points( &
      src_grid_id, src_num_cells, YAC_LOCATION_CELL, &
      src_x_cells(:,src_rank+1), src_y_cells, src_cell_point_id)
  ELSE
    src_cell_point_id = -1
  END IF
  IF (is_tgt) THEN
    CALL yac_fdef_points( &
      tgt_grid_id, tgt_num_cells, YAC_LOCATION_CELL, &
      tgt_x_cells, tgt_y_cells(:,tgt_rank+1), tgt_cell_point_id)
  ELSE
    tgt_cell_point_id = -1
  END IF

  ! define fields
  IF (is_src) THEN
    CALL yac_fdef_field_mask( &
      src_field_name, src_comp_id, (/src_cell_point_id/), &
      (/src_cell_mask_id/), 1, 1, &
      "1", YAC_TIME_UNIT_SECOND, src_field_id)
    CALL yac_fdef_field_mask( &
      src_field_csr_name, src_comp_id, (/src_cell_point_id/), &
      (/src_cell_mask_id/), 1, 1, &
      "1", YAC_TIME_UNIT_SECOND, src_field_csr_id)
    CALL yac_fdef_field_mask( &
      src_field_frac_name, src_comp_id, (/src_cell_point_id/), &
      (/src_cell_mask_id/), 1, 1, &
      "1", YAC_TIME_UNIT_SECOND, src_field_frac_id)
  ELSE
    src_field_id = -1
    src_field_csr_id = -1
    src_field_frac_id = -1
  END IF
  IF (is_tgt) THEN
    CALL yac_fdef_field( &
      tgt_field_name, tgt_comp_id, (/tgt_cell_point_id/), 1, 1, &
      "1", YAC_TIME_UNIT_SECOND, tgt_field_id)
    CALL yac_fdef_field( &
      tgt_field_csr_name, tgt_comp_id, (/tgt_cell_point_id/), 1, 1, &
      "1", YAC_TIME_UNIT_SECOND, tgt_field_csr_id)
    CALL yac_fdef_field( &
      tgt_field_frac_name, tgt_comp_id, (/tgt_cell_point_id/), 1, 1, &
      "1", YAC_TIME_UNIT_SECOND, tgt_field_frac_id)
  ELSE
    tgt_field_id = -1
    tgt_field_csr_id = -1
    tgt_field_frac_id = -1
  END IF

  ! define couple
  CALL yac_fget_interp_stack_config(interp_stack_config)
  CALL yac_fadd_interp_stack_config_conservative( &
    interp_stack_config, 1, 0, 1, YAC_CONSERV_DESTAREA)
  CALL yac_fadd_interp_stack_config_fixed( &
    interp_stack_config, fixed_value)
  CALL yac_fdef_couple( &
    src_comp_name, src_grid_name, src_field_name, &
    tgt_comp_name, tgt_grid_name, tgt_field_name, &
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
    interp_stack_config, use_raw_exchange=.TRUE.)
  CALL yac_fdef_couple( &
    src_comp_name, src_grid_name, src_field_csr_name, &
    tgt_comp_name, tgt_grid_name, tgt_field_csr_name, &
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
    interp_stack_config, use_raw_exchange=.TRUE.)
  CALL yac_fdef_couple( &
    src_comp_name, src_grid_name, src_field_frac_name, &
    tgt_comp_name, tgt_grid_name, tgt_field_frac_name, &
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
    interp_stack_config, use_raw_exchange=.TRUE.)
  CALL yac_ffree_interp_stack_config(interp_stack_config)

  ! enable fractional masking
  CALL yac_fenable_field_frac_mask( &
    src_comp_name, src_grid_name, src_field_frac_name, &
    frac_mask_fallback_value)

  ! end definition phase
  CALL yac_fenddef()

  ! extract interpolation weight data
  IF (is_tgt) THEN
    CALL yac_fget_raw_interp_weights_data( &
      tgt_field_id, interp_weights_data%frac_mask_fallback_value, &
      interp_weights_data%scaling_factor, interp_weights_data%scaling_summand, &
      interp_weights_data%num_fixed_values, interp_weights_data%fixed_values, &
      interp_weights_data%num_tgt_per_fixed_value, &
      interp_weights_data%tgt_idx_fixed, interp_weights_data%num_wgt_tgt, &
      interp_weights_data%wgt_tgt_idx, interp_weights_data%num_src_per_tgt, &
      interp_weights_data%weights, interp_weights_data%src_field_idx, &
      interp_weights_data%src_idx, interp_weights_data%num_src_fields, &
      interp_weights_data%src_field_buffer_sizes)
    CALL yac_fget_raw_interp_weights_data_csr( &
      tgt_field_csr_id, &
      interp_weights_data_csr%frac_mask_fallback_value, &
      interp_weights_data_csr%scaling_factor, &
      interp_weights_data_csr%scaling_summand, &
      interp_weights_data_csr%num_fixed_values, &
      interp_weights_data_csr%fixed_values, &
      interp_weights_data_csr%num_tgt_per_fixed_value, &
      interp_weights_data_csr%tgt_idx_fixed, &
      interp_weights_data_csr%src_indptr, &
      interp_weights_data_csr%weights, &
      interp_weights_data_csr%src_field_idx, &
      interp_weights_data_csr%src_idx, &
      interp_weights_data_csr%num_src_fields, &
      interp_weights_data_csr%src_field_buffer_sizes)
    CALL yac_fget_raw_interp_weights_data( &
      tgt_field_frac_id, &
      interp_weights_data_frac%frac_mask_fallback_value, &
      interp_weights_data_frac%scaling_factor, &
      interp_weights_data_frac%scaling_summand, &
      interp_weights_data_frac%num_fixed_values, &
      interp_weights_data_frac%fixed_values, &
      interp_weights_data_frac%num_tgt_per_fixed_value, &
      interp_weights_data_frac%tgt_idx_fixed, &
      interp_weights_data_frac%num_wgt_tgt, &
      interp_weights_data_frac%wgt_tgt_idx, &
      interp_weights_data_frac%num_src_per_tgt, &
      interp_weights_data_frac%weights, &
      interp_weights_data_frac%src_field_idx, &
      interp_weights_data_frac%src_idx, &
      interp_weights_data_frac%num_src_fields, &
      interp_weights_data_frac%src_field_buffer_sizes)
  ELSE
    interp_weights_data%num_src_fields = 0
    interp_weights_data_csr%num_src_fields = 0
    interp_weights_data_frac%num_src_fields = 0
  END IF

  collection_size = 1

  SELECT CASE ( MERGE(1,0,is_src) + MERGE(2,0,is_tgt) )

    CASE ( 0 ) ! neither source nor target

      ALLOCATE( &
        src_field(0,0,0), src_frac_mask(0,0,0,0), &
        src_field_buffer(0,0), src_frac_mask_buffer(0,0), &
        src_field_ptr(0,0), src_frac_mask_ptr(0,0,0), &
        src_field_buffer_ptr(0,0), src_frac_mask_buffer_ptr(0,0), &
        tgt_field(0,0), ref_tgt_field(0,0), ref_tgt_field_frac(0,0,0))

    CASE ( 1 ) ! only source

      ALLOCATE( &
        src_field_buffer(0,0), src_field_buffer_ptr(0,0), &
        src_frac_mask_buffer(0,0), src_frac_mask_buffer_ptr(0,0), &
        tgt_field(0,0), ref_tgt_field(0,0), ref_tgt_field_frac(0,0,0))
      num_cells = PRODUCT(src_num_cells)
      ALLOCATE( &
        src_field, &
        SOURCE = &
          RESHAPE( &
            REAL((/0+1,2+1/), wp), &
            (/num_cells,NUM_SRC_FIELDS,collection_size/)))
      ALLOCATE( &
        src_frac_mask, &
        SOURCE = &
          RESHAPE( &
            REAL(&
              (/1.0,1.0,0.5,0.5/), wp), &
              (/num_cells, NUM_SRC_FIELDS, collection_size, NUM_TESTS_FRAC/)))
      ALLOCATE( &
        src_field_ptr(NUM_SRC_FIELDS,collection_size), &
        src_frac_mask_ptr(NUM_SRC_FIELDS,collection_size, NUM_TESTS_FRAC))
      DO i = 1, collection_size
        src_field_ptr(1,i)%p => src_field(:,1,i)
        DO t = 1, NUM_TESTS_FRAC
          src_frac_mask_ptr(1,i,t)%p => src_frac_mask(:,1,i,t)
        END DO
      END DO

      ! send source data
      CALL yac_fput( &
        src_field_id, num_cells, collection_size, src_field(:,1,:), &
        src_info, ierror)
      CALL yac_fput( &
        src_field_id, num_cells, NUM_SRC_FIELDS, collection_size, src_field, &
        src_info, ierror)
      CALL yac_fput( &
        src_field_id, NUM_SRC_FIELDS, collection_size, src_field_ptr, &
        src_info, ierror)
#ifdef TEST_GET_ASYNC
      CALL yac_fput( &
        src_field_id, NUM_SRC_FIELDS, collection_size, src_field_ptr, &
        src_info, ierror)
#endif
      CALL yac_fput( &
        src_field_csr_id, num_cells, collection_size, src_field(:,1,:), &
        src_info, ierror)

      DO t = 1, NUM_TESTS_FRAC

        ! send source data and fractional mask
        CALL yac_fput( &
          src_field_frac_id, num_cells, collection_size, &
          src_field(:,1,:), src_frac_mask(:,1,:,t), src_info, ierror)
        CALL yac_fput( &
          src_field_frac_id, num_cells, NUM_SRC_FIELDS, collection_size, &
          src_field, src_frac_mask(:,:,:,t), src_info, ierror)
        CALL yac_fput( &
          src_field_frac_id, NUM_SRC_FIELDS, collection_size, &
          src_field_ptr, src_frac_mask_ptr(:,:,t), src_info, ierror)
#ifdef TEST_GET_ASYNC
        CALL yac_fput( &
          src_field_frac_id, NUM_SRC_FIELDS, collection_size, &
          src_field_ptr, src_frac_mask_ptr(:,:,t), src_info, ierror)
#endif

      END DO

    CASE ( 2 ) ! only target

      ALLOCATE( &
        src_field(0,0,0), src_frac_mask(0,0,0,0), &
        src_field_ptr(0,0), src_frac_mask_ptr(0,0,0))
      src_field_buffer_size = SUM(interp_weights_data%src_field_buffer_sizes)
      ALLOCATE( &
        src_field_buffer(src_field_buffer_size,collection_size), &
        src_frac_mask_buffer(src_field_buffer_size,collection_size))
      ALLOCATE( &
        src_field_buffer_ptr(NUM_SRC_FIELDS,collection_size), &
        src_frac_mask_buffer_ptr(NUM_SRC_FIELDS,collection_size))
      num_cells = PRODUCT(tgt_num_cells)
      ALLOCATE(tgt_field(num_cells, collection_size), SOURCE=-1.0_wp)
      ALLOCATE( &
        ref_tgt_field, &
        SOURCE = &
          RESHAPE( &
              REAL( &
                (/0.25*(3),0.25*(2+3+4),0.25*(2+4), &
                  0.25*(3),0.25*(3+4),0.25*(4)/), wp), &
              (/num_cells, collection_size/)))
      ALLOCATE( &
        ref_tgt_field_frac, &
        SOURCE = &
          RESHAPE( &
              REAL( &
                (/0.25*(3),0.25*(2+3+4),0.25*(2+4), &
                  0.25*(3),0.25*(3+4),0.25*(4), &
                  (0.25*(1.0*3))/((0.25*1.0)/(0.25)), &
                  (0.25*(0.0*2+0.5*3+1.0*4))/((0.25*(0.0+0.5+1.0))/(0.25+0.25+0.25)), &
                  (0.25*(0.0*2+1.0*4))/((0.25*(0.0+1.0))/(0.25+0.25)), &
                  (0.25*(0.5*3))/((0.25*0.5)/(0.25)), &
                  (0.25*(0.5*3+1.0*4))/((0.25*(0.5+1.0))/(0.25+0.25)), &
                  (0.25*(1.0*4))/((0.25*1.0)/(0.25))/), wp), &
              (/num_cells, collection_size, NUM_TESTS_FRAC/)))
      DO i = 1, collection_size
        src_field_buffer_ptr(1,i)%p => src_field_buffer(:,i)
        src_frac_mask_buffer_ptr(1,i)%p => src_frac_mask_buffer(:,i)
      END DO

      ! receive source field buffer
      CALL yac_fget_raw( &
        tgt_field_id, src_field_buffer_size, collection_size, &
        src_field_buffer, tgt_info, ierror)

      ! compute target field
      num_cells = PRODUCT(tgt_num_cells)
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

      ! receive source field buffer
      CALL yac_fget_raw( &
        tgt_field_id, src_field_buffer_size, collection_size, &
        src_field_buffer, tgt_info, ierror)

      ! compute target field
      num_cells = PRODUCT(tgt_num_cells)
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

      ! receive source field buffer
      CALL yac_fget_raw( &
        tgt_field_id, NUM_SRC_FIELDS, collection_size, src_field_buffer_ptr, &
        tgt_info, ierror)

      ! compute target field
      num_cells = PRODUCT(tgt_num_cells)
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

#ifdef TEST_GET_ASYNC
      ! receive source field buffer
      CALL yac_fget_raw_async( &
        tgt_field_id, NUM_SRC_FIELDS, collection_size, src_field_buffer_ptr, &
        tgt_info, ierror)
      CALL yac_fwait(tgt_field_id)

      ! compute target field
      num_cells = PRODUCT(tgt_num_cells)
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))
#endif

      ! receive source field buffer
      CALL yac_fget_raw( &
        tgt_field_csr_id, src_field_buffer_size, collection_size, &
        src_field_buffer, tgt_info, ierror)

      ! compute target field
      num_cells = PRODUCT(tgt_num_cells)
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data_csr, collection_size, &
        src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

      DO t = 1, NUM_TESTS_FRAC

        ! receive source field buffer and source fractional mask buffer
        CALL yac_fget_raw( &
          tgt_field_frac_id, src_field_buffer_size, collection_size, &
          src_field_buffer, src_frac_mask_buffer, tgt_info, ierror)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))

        ! receive source field buffer and source fractional mask buffer
        CALL yac_fget_raw( &
          tgt_field_frac_id, src_field_buffer_size, collection_size, &
          src_field_buffer, src_frac_mask_buffer, tgt_info, ierror)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))

        ! receive source field buffer and source fractional mask buffer
        CALL yac_fget_raw( &
          tgt_field_frac_id, NUM_SRC_FIELDS, collection_size, &
          src_field_buffer_ptr, src_frac_mask_buffer_ptr, tgt_info, ierror)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))

#ifdef TEST_GET_ASYNC
        ! receive source field buffer and source fractional mask buffer
        CALL yac_fget_raw_async( &
          tgt_field_frac_id, NUM_SRC_FIELDS, collection_size, &
          src_field_buffer_ptr, src_frac_mask_buffer_ptr, tgt_info, ierror)
        CALL yac_fwait(tgt_field_frac_id)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))
#endif

      END DO

    CASE ( 3 ) ! source and target

      num_cells = PRODUCT(src_num_cells)
      ALLOCATE( &
        src_field, &
        SOURCE = &
          RESHAPE( &
            REAL((/1+1,3+1/), wp), &
            (/num_cells,NUM_SRC_FIELDS,collection_size/)))
      ALLOCATE( &
        src_frac_mask, &
        SOURCE = &
          RESHAPE( &
            REAL(&
              (/1.0,1.0,0.0,1.0/), wp), &
              (/num_cells, NUM_SRC_FIELDS, collection_size, NUM_TESTS_FRAC/)))
      ALLOCATE( &
        src_field_ptr(NUM_SRC_FIELDS,collection_size), &
        src_frac_mask_ptr(NUM_SRC_FIELDS,collection_size, NUM_TESTS_FRAC))
      num_cells = PRODUCT(tgt_num_cells)
      src_field_buffer_size = SUM(interp_weights_data%src_field_buffer_sizes)
      ALLOCATE( &
        src_field_buffer(src_field_buffer_size,collection_size), &
        src_frac_mask_buffer(src_field_buffer_size,collection_size), &
        src_field_buffer_ptr(NUM_SRC_FIELDS,collection_size), &
        src_frac_mask_buffer_ptr(NUM_SRC_FIELDS,collection_size))
      ALLOCATE(tgt_field(num_cells, collection_size), SOURCE=-1.0_wp)
      ALLOCATE( &
        ref_tgt_field, &
        SOURCE = &
          RESHAPE( &
              REAL( &
                (/REAL(fixed_value),0.25*(2),0.25*(2), &
                  0.25*(3),0.25*(2+3+4),0.25*(2+4)/), wp), &
              (/num_cells, collection_size/)))
      ALLOCATE( &
        ref_tgt_field_frac, &
        SOURCE = &
          RESHAPE( &
              REAL( &
                (/REAL(fixed_value),0.25*(2),0.25*(2), &
                  0.25*(3),0.25*(2+3+4),0.25*(2+4), &
                  REAL(fixed_value), &
                  REAL(frac_mask_fallback_value), &
                  REAL(frac_mask_fallback_value), &
                  (0.25*(0.5*3))/((0.25*0.5)/(0.25)), &
                  (0.25*(0.0*2+0.5*3+1.0*4))/((0.25*(0.0+0.5+1.0))/(0.25+0.25+0.25)), &
                  (0.25*(0.0*2+1.0*4))/((0.25*(0.0+1.0))/(0.25+0.25))/), wp), &
              (/num_cells, collection_size, NUM_TESTS_FRAC/)))
      DO i = 1, collection_size
        src_field_ptr(1,i)%p => src_field(:,1,i)
        DO t = 1, NUM_TESTS_FRAC
          src_frac_mask_ptr(1,i,t)%p => src_frac_mask(:,1,i,t)
        END DO
      END DO
      DO i = 1, collection_size
        src_field_buffer_ptr(1,i)%p => src_field_buffer(:,i)
        src_frac_mask_buffer_ptr(1,i)%p => src_frac_mask_buffer(:,i)
      END DO

      ! exchange source field buffer
      CALL yac_fexchange_raw( &
        src_field_id, tgt_field_id, PRODUCT(src_num_cells), &
        src_field_buffer_size, collection_size, &
        src_field(:,1,:), src_field_buffer, src_info, tgt_info, ierror)

      ! compute target field
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

      ! exchange source field buffer
      CALL yac_fexchange_raw( &
        src_field_id, tgt_field_id, PRODUCT(src_num_cells), &
        NUM_SRC_FIELDS, src_field_buffer_size, collection_size, &
        src_field, src_field_buffer, src_info, tgt_info, ierror)

      ! compute target field
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

      ! exchange source field buffer
      CALL yac_fexchange_raw( &
        src_field_id, tgt_field_id, NUM_SRC_FIELDS, collection_size, &
        src_field_ptr, src_field_buffer_ptr, &
        src_info, tgt_info, ierror)

      ! compute target field
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

#ifdef TEST_GET_ASYNC
      CALL yac_fexchange_raw( &
        src_field_id, tgt_field_id, NUM_SRC_FIELDS, collection_size, &
        src_field_ptr, src_field_buffer_ptr, &
        src_info, tgt_info, ierror)

      ! compute target field
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data, collection_size, src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))
#endif

      ! exchange source field buffer
      CALL yac_fexchange_raw( &
        src_field_csr_id, tgt_field_csr_id, PRODUCT(src_num_cells), &
        src_field_buffer_size, collection_size, &
        src_field(:,1,:), src_field_buffer, src_info, tgt_info, ierror)

      ! compute target field
      CALL compute_tgt_field( &
        tgt_field, interp_weights_data_csr, collection_size, &
        src_field_buffer)

      ! check results
      CALL test(ALL(ABS(tgt_field - ref_tgt_field) < tol))

      DO t = 1, NUM_TESTS_FRAC

        ! exchange source field buffer and source fractional mask buffer
        CALL yac_fexchange_raw( &
          src_field_frac_id, tgt_field_frac_id, PRODUCT(src_num_cells), &
          src_field_buffer_size, collection_size, &
          src_field(:,1,:), src_frac_mask(:,1,:,t), &
          src_field_buffer, src_frac_mask_buffer, &
          src_info, tgt_info, ierror)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))

        ! exchange source field buffer and source fractional mask buffer
        CALL yac_fexchange_raw( &
          src_field_frac_id, tgt_field_frac_id, &
          PRODUCT(src_num_cells), NUM_SRC_FIELDS, &
          src_field_buffer_size, collection_size, &
          src_field, src_frac_mask(:,:,:,t), &
          src_field_buffer, src_frac_mask_buffer, &
          src_info, tgt_info, ierror)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))

        ! exchange source field buffer and source fractional mask buffer
        CALL yac_fexchange_raw( &
          src_field_frac_id, tgt_field_frac_id, NUM_SRC_FIELDS, &
          collection_size, src_field_ptr, src_frac_mask_ptr(:,:,t), &
          src_field_buffer_ptr, src_frac_mask_buffer_ptr, &
          src_info, tgt_info, ierror)

        ! compute target field
        tgt_field = -1.0_wp
        CALL compute_tgt_field( &
          tgt_field, interp_weights_data_frac, collection_size, &
          src_field_buffer, src_frac_mask_buffer)

        ! check results
        CALL test(ALL(ABS(tgt_field - ref_tgt_field_frac(:,:,t)) < tol))

#ifdef TEST_GET_ASYNC
        ! exchange source field buffer and source fractional mask buffer
        CALL yac_fexchange_raw( &
          src_field_frac_id, tgt_field_frac_id, NUM_SRC_FIELDS, &
          collection_size, src_field_ptr, src_frac_mask_ptr(:,:,t), &
          src_field_buffer_ptr, src_frac_mask_buffer_ptr, &
          src_info, tgt_info, ierror)
#endif

      END DO

    CASE DEFAULT
      ALLOCATE( &
        src_field(0,0,0), src_frac_mask(0,0,0,0), &
        src_field_buffer(0,0), src_frac_mask_buffer(0,0), &
        src_field_ptr(0,0), src_frac_mask_ptr(0,0,0), &
        src_field_buffer_ptr(0,0), src_frac_mask_buffer_ptr(0,0), &
        tgt_field(0,0), ref_tgt_field(0,0), ref_tgt_field_frac(0,0,0))
      WRITE ( * , * ) "Invalid process type", is_src, is_tgt, rank, npes
      CALL test ( .FALSE. )
      CALL stop_test
      CALL exit_tests
      CALL mpi_abort (MPI_COMM_WORLD, 999, ierror)

  END SELECT

  DEALLOCATE( &
    src_field, src_frac_mask, src_field_buffer, src_frac_mask_buffer, &
    src_field_ptr, src_frac_mask_ptr, &
    src_field_buffer_ptr, src_frac_mask_buffer_ptr, &
    tgt_field, ref_tgt_field, ref_tgt_field_frac)
  IF (is_tgt) THEN
    CALL interp_weights_data_free(interp_weights_data)
    CALL interp_weights_data_free(interp_weights_data_csr)
    CALL interp_weights_data_free(interp_weights_data_frac)
  END IF

  CALL yac_ffinalize

  CALL stop_test
  CALL exit_tests

CONTAINS

  SUBROUTINE compute_tgt_field( &
    tgt_field, interp_weights_data, collection_size, &
    src_field_buffer, src_frac_mask_buffer)

    REAL(KIND=wp), INTENT(OUT) :: tgt_field(:,:)
    TYPE(t_interp_weights_data), INTENT(IN) :: interp_weights_data
    INTEGER, INTENT(IN) :: collection_size
    REAL(KIND=wp), INTENT(IN) :: src_field_buffer(:,:)
    REAL(KIND=wp), OPTIONAL, INTENT(IN) :: src_frac_mask_buffer(:,:)

    INTEGER :: collection_idx, i, j, k
    LOGICAL :: with_frac_mask
    DOUBLE PRECISION :: fixed_value, tgt_value, frac_weight_sum, weight_sum
    INTEGER :: num_src
    INTEGER :: accu
    INTEGER :: src_field_buffer_offsets(interp_weights_data%num_src_fields)
    INTEGER :: tgt_idx, num_wgt_tgt
    LOGICAL :: use_csr

    use_csr = ALLOCATED(interp_weights_data%src_indptr)

    IF (use_csr) THEN
      num_wgt_tgt = SIZE(interp_weights_data%src_indptr) - 1
    ELSE
      num_wgt_tgt = interp_weights_data%num_wgt_tgt
    ENDIF

    with_frac_mask = PRESENT(src_frac_mask_buffer)

    tgt_field = -1.0

    accu = 0
    DO i = 1, interp_weights_data%num_src_fields
      src_field_buffer_offsets(i) = accu
      accu = accu + interp_weights_data%src_field_buffer_sizes(i)
    END DO

    DO collection_idx = 1, collection_size

      ! set fixed targets
      k = 1
      DO i = 1, interp_weights_data%num_fixed_values

        fixed_value = interp_weights_data%fixed_values(i)

        DO j = 1, interp_weights_data%num_tgt_per_fixed_value(i)

          tgt_field(interp_weights_data%tgt_idx_fixed(k), collection_idx) = &
            REAL(fixed_value, wp)
          k = k + 1

        END DO
      END DO

      IF (with_frac_mask) THEN

        ! set weighted targets
        k = 1
        DO i = 1, num_wgt_tgt

          IF (use_csr) THEN
            num_src = &
              interp_weights_data%src_indptr(i+1) - &
              interp_weights_data%src_indptr(i)
            tgt_idx = i
          ELSE
            num_src = interp_weights_data%num_src_per_tgt(i)
            tgt_idx = interp_weights_data%wgt_tgt_idx(i)
          END IF

          IF (num_src == 0) CYCLE

          tgt_value = 0.0_wp
          frac_weight_sum = 0.0_wp
          weight_sum = 0.0_wp

          DO j = 1, num_src

            tgt_value = &
              tgt_value + interp_weights_data%weights(k) * &
              src_field_buffer( &
                interp_weights_data%src_idx(k) + &
                src_field_buffer_offsets(interp_weights_data%src_field_idx(k)), &
                collection_idx)
            frac_weight_sum = &
              frac_weight_sum + interp_weights_data%weights(k) * &
              src_frac_mask_buffer( &
                interp_weights_data%src_idx(k) + &
                src_field_buffer_offsets(interp_weights_data%src_field_idx(k)), &
                collection_idx)
            weight_sum = weight_sum + interp_weights_data%weights(k)
            k = k + 1

          END DO

          IF (ABS(frac_weight_sum) > FRAC_MASK_TOL) THEN
            tgt_field(tgt_idx,collection_idx) = &
              REAL( &
                interp_weights_data%scaling_factor * &
                (tgt_value / frac_weight_sum) * weight_sum + &
                interp_weights_data%scaling_summand, wp)
          ELSE
            tgt_field(tgt_idx,collection_idx) = &
              REAL(interp_weights_data%frac_mask_fallback_value, wp)
          END IF

        END DO

      ELSE

        ! set weighted targets
        k = 1
        DO i = 1, num_wgt_tgt

          IF (use_csr) THEN
            num_src = &
              interp_weights_data%src_indptr(i+1) - &
              interp_weights_data%src_indptr(i)
            tgt_idx = i
          ELSE
            num_src = interp_weights_data%num_src_per_tgt(i)
            tgt_idx = interp_weights_data%wgt_tgt_idx(i)
          END IF

          IF (num_src == 0) CYCLE

          tgt_value = 0.0_wp

          DO j = 1, num_src

            tgt_value = &
              tgt_value + interp_weights_data%weights(k) * &
              src_field_buffer( &
                interp_weights_data%src_idx(k) + &
                src_field_buffer_offsets(interp_weights_data%src_field_idx(k)), &
                collection_idx)
            k = k + 1

          END DO

          tgt_field(tgt_idx,collection_idx) = &
            REAL( &
              interp_weights_data%scaling_factor * tgt_value + &
              interp_weights_data%scaling_summand, wp)

        END DO

      END IF

    END DO

  END SUBROUTINE compute_tgt_field

  SUBROUTINE interp_weights_data_free(interp_weights_data)

    TYPE(t_interp_weights_data), INTENT(INOUT) :: interp_weights_data

    DEALLOCATE( &
      interp_weights_data%fixed_values, &
      interp_weights_data%num_tgt_per_fixed_value, &
      interp_weights_data%tgt_idx_fixed, &
      interp_weights_data%weights, &
      interp_weights_data%src_field_idx, &
      interp_weights_data%src_idx, &
      interp_weights_data%src_field_buffer_sizes)
    IF (ALLOCATED(interp_weights_data%wgt_tgt_idx)) &
      DEALLOCATE(interp_weights_data%wgt_tgt_idx)
    IF (ALLOCATED(interp_weights_data%num_src_per_tgt)) &
      DEALLOCATE(interp_weights_data%num_src_per_tgt)
    IF (ALLOCATED(interp_weights_data%src_indptr)) &
      DEALLOCATE(interp_weights_data%src_indptr)
  END SUBROUTINE interp_weights_data_free
