!> @file comin_variable_types.F90
!! @brief Data types for variable definition
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_variable_types

  USE iso_c_binding,           ONLY: c_int, c_int64_t, c_char, c_bool, c_ptr, c_null_ptr, c_f_pointer
  USE comin_setup_constants,   ONLY: wp, COMIN_ZAXIS_3D, &
    &                                COMIN_HGRID_UNSTRUCTURED_CELL, dp, sp, COMIN_DIM_SEMANTICS_CONTAINER, &
    &                                COMIN_DIM_SEMANTICS_UNDEF
  USE comin_metadata_types,    ONLY: t_comin_var_metadata
  USE comin_c_utils,           ONLY: convert_c_string

  IMPLICIT NONE

  PUBLIC

#include "comin_global.inc"

  EXTERNAL :: comin_plugin_finish_external

  ! ------------------------------------
  ! data types for variable definition
  ! ------------------------------------

  !> Variable descriptor.
  !> identifies (uniquely) a variable. Do not confuse with meta-data
  !! @ingroup common
  TYPE :: t_comin_var_descriptor
    CHARACTER(LEN=:), ALLOCATABLE :: name
    ! domain id
    INTEGER                       :: id
  END TYPE t_comin_var_descriptor

  INTERFACE t_comin_var_descriptor
    MODULE PROCEDURE create_comin_var_descriptor
    MODULE PROCEDURE create_comin_var_descriptor_from_c
  END INTERFACE

  TYPE, BIND(C) :: t_comin_var_descriptor_c
    CHARACTER(KIND=c_char) :: name(COMIN_MAX_LEN_VAR_NAME+1)
    ! domain id
    INTEGER(kind=c_int)    :: id
  END TYPE t_comin_var_descriptor_c

  !> Variable pointer. Fortran interface for accessing variables
  !! @ingroup common
  TYPE :: t_comin_var_handle
    !> pointer to the internal data structure, private (not part of the Fortran interface)
    TYPE(t_comin_var_item), PRIVATE, POINTER :: var_item

  CONTAINS
    GENERIC, PUBLIC :: get_ptr       => get_ptr_dp, get_ptr_sp, get_ptr_i
    PROCEDURE, PRIVATE :: get_ptr_dp => comin_var_get_ptr_dp
    PROCEDURE, PRIVATE :: get_ptr_sp => comin_var_get_ptr_sp
    PROCEDURE, PRIVATE :: get_ptr_i  => comin_var_get_ptr_i
    PROCEDURE, PUBLIC  :: array_shape => comin_var_get_array_shape
    PROCEDURE, PUBLIC  :: descriptor => comin_var_get_descriptor
    PROCEDURE, PUBLIC  :: lcontainer => comin_var_get_lcontainer
    PROCEDURE, PUBLIC  :: ncontained => comin_var_get_ncontained
    PROCEDURE, PUBLIC  :: dim_semantics => comin_var_get_dim_semantics
    PROCEDURE, PUBLIC  :: valid      => comin_var_get_valid
    GENERIC, PUBLIC    :: to_3d      => to_3d_dp, to_3d_sp, to_3d_i
    PROCEDURE, PRIVATE :: to_3d_dp   => comin_var_to_3d_dp
    PROCEDURE, PRIVATE :: to_3d_sp   => comin_var_to_3d_sp
    PROCEDURE, PRIVATE :: to_3d_i   => comin_var_to_3d_i
  END TYPE t_comin_var_handle

  !> Variable item
  TYPE :: t_comin_var_item
    !> the var_descriptor
    TYPE(t_comin_var_descriptor) :: descriptor

    !> the (current) pointer to the data
    ! REAL(wp), POINTER :: ptr(:,:,:,:,:) => NULL()
    TYPE(c_ptr) :: cptr = c_null_ptr

    !> type id for the array cptr is pointing to
    INTEGER :: type_id

    !> shape for the array cptr is pointing to
    INTEGER :: array_shape(5)

    !> the (current) device ptr to the data (if any)
    TYPE(c_ptr) :: device_ptr = c_null_ptr

    ! index positions in the 5D array.
    INTEGER, DIMENSION(5) :: dim_semantics

    !> if (tracer==.TRUE.) and (ncontained > 0), then the variable
    !  pointer refers to an array slice pointer
    !  ptr(:,:,:,:,ncontained)
    INTEGER :: ncontained = 0
    !> LOGICAL flag. TRUE, if this is a container (contains variables)
    LOGICAL(kind=c_bool) :: lcontainer = .FALSE.

    !> metadata store
    TYPE(t_comin_var_metadata)           :: metadata
  END TYPE t_comin_var_item

  !> Variable list for context access
  TYPE :: t_comin_var_context_item
    TYPE(t_comin_var_item), POINTER :: var_item
    INTEGER :: access_flag
  END TYPE t_comin_var_context_item

  !> Information on requested variables
  TYPE :: t_comin_request_item
    TYPE(t_comin_var_descriptor)    :: descriptor
    TYPE(t_comin_var_metadata)      :: metadata
    INTEGER, ALLOCATABLE            :: moduleID(:)

    !> LOGICAL flag. TRUE, if this variable is intended to be used
    !> exclusively by a particular 3rd party plugin:
    LOGICAL(kind=c_bool) :: lmodexclusive = .FALSE.
  END TYPE t_comin_request_item

  INTERFACE
    SUBROUTINE comin_var_sync_device_mem_fct(var_ptr, direction)
      IMPORT t_comin_var_handle
      TYPE(t_comin_var_handle), INTENT(IN) :: var_ptr
      LOGICAL, INTENT(IN) :: direction
    END SUBROUTINE comin_var_sync_device_mem_fct
  END INTERFACE
  INTERFACE
    SUBROUTINE comin_var_sync_halo_fct(var_ptr, halo_sync_mode)
      IMPORT t_comin_var_handle
      TYPE(t_comin_var_handle), INTENT(IN) :: var_ptr
      INTEGER, INTENT(IN) :: halo_sync_mode
    END SUBROUTINE comin_var_sync_halo_fct
  END INTERFACE

  ! ------------------------------------
  ! lists of exposed variables
  ! ------------------------------------

  !> Array of variable lists (array of pointer lists) each entry
  !  stores the lists of variables registered for the context
  !  (dimension of array) points to the first element of the variable
  !  list
  !  - contains TYPE(t_comin_var_item)
  TYPE :: t_comin_var_list_context
    TYPE(C_PTR) :: var_list
  END TYPE t_comin_var_list_context

CONTAINS

  !> Create a variable descriptor from a variable name and id.
  FUNCTION create_comin_var_descriptor(name, id) RESULT(desc)
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: id

    TYPE(t_comin_var_descriptor) :: desc

    desc%name = TRIM(ADJUSTL(name))
    desc%id = id

  END FUNCTION create_comin_var_descriptor

  !> Create a variable descriptor from a C structure.
  FUNCTION create_comin_var_descriptor_from_c(desc_c) RESULT(desc)
    TYPE(t_comin_var_descriptor_c), INTENT(IN) :: desc_c

    TYPE(t_comin_var_descriptor) :: desc

    desc = create_comin_var_descriptor(convert_c_string(desc_c%name), desc_c%id)
  END FUNCTION create_comin_var_descriptor_from_c

  !> compare two variable descriptors.
  FUNCTION comin_var_descr_match(var_descriptor1, var_descriptor2)
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descriptor1, var_descriptor2
    LOGICAL :: comin_var_descr_match
    ! local
    LOGICAL :: l_name, l_domain

    l_domain = (var_descriptor1%id == var_descriptor2%id)

    comin_var_descr_match = l_domain
    IF (.NOT. comin_var_descr_match) RETURN

    l_name = (var_descriptor1%name == var_descriptor2%name)

    comin_var_descr_match = l_name
  END FUNCTION comin_var_descr_match

  FUNCTION comin_var_ptr_init(var_item) &
    RESULT(var_ptr)
    TYPE(t_comin_var_item), POINTER, INTENT(IN) :: var_item
    TYPE(t_comin_var_handle) :: var_ptr
    var_ptr = t_comin_var_handle(var_item = var_item)
  END FUNCTION comin_var_ptr_init

  SUBROUTINE comin_var_get_ptr_dp(this, ptr)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    REAL(dp), POINTER, INTENT(INOUT)        :: ptr(:,:,:,:,:)

    CALL C_F_POINTER(this%var_item%cptr, ptr, this%var_item%array_shape)
  END SUBROUTINE comin_var_get_ptr_dp

  SUBROUTINE comin_var_get_ptr_sp(this, ptr)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    REAL(sp), POINTER, INTENT(INOUT)        :: ptr(:,:,:,:,:)

    CALL C_F_POINTER(this%var_item%cptr, ptr, this%var_item%array_shape)
  END SUBROUTINE comin_var_get_ptr_sp

  SUBROUTINE comin_var_get_ptr_i(this, ptr)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    INTEGER(C_INT), POINTER, INTENT(INOUT)        :: ptr(:,:,:,:,:)

    CALL C_F_POINTER(this%var_item%cptr, ptr, this%var_item%array_shape)
  END SUBROUTINE comin_var_get_ptr_i

  FUNCTION comin_var_get_array_shape(this) &
    RESULT(array_shape)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    INTEGER :: array_shape(5)
    array_shape = this%var_item%array_shape
  END FUNCTION comin_var_get_array_shape

  FUNCTION comin_var_get_descriptor(this) &
    RESULT(descriptor)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    TYPE(t_comin_var_descriptor) :: descriptor
    descriptor = this%var_item%descriptor
  END FUNCTION comin_var_get_descriptor

  FUNCTION comin_var_get_lcontainer(this) &
    RESULT(lcontainer)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    LOGICAL :: lcontainer
    lcontainer = this%var_item%lcontainer
  END FUNCTION comin_var_get_lcontainer

  FUNCTION comin_var_get_ncontained(this) &
    RESULT(ncontained)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    INTEGER :: ncontained
    ncontained = this%var_item%ncontained
  END FUNCTION comin_var_get_ncontained

  FUNCTION comin_var_get_dim_semantics(this) &
    RESULT(dim_semantics)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    INTEGER :: dim_semantics(5)
    dim_semantics = this%var_item%dim_semantics
  END FUNCTION comin_var_get_dim_semantics

  FUNCTION comin_var_get_valid(this) &
    RESULT(valid)
    CLASS(t_comin_var_handle), INTENT(IN) :: this
    LOGICAL :: valid
    valid = ASSOCIATED(this%var_item)
  END FUNCTION comin_var_get_valid

  !> Convenience operation for accessing 2D/3D fields.
  !! @ingroup plugin_interface
  !!
  !! Assumes that the last dimension is not used!
  SUBROUTINE comin_var_to_3d_dp(var, slice)
    CLASS(t_comin_var_handle), INTENT(IN)  :: var
    REAL(dp), POINTER :: slice(:,:,:)
    REAL(dp), POINTER :: tmp_ptr(:,:,:,:,:)
    INTEGER :: pos_jn

    ! this operation is invalid if the field is a container
    IF (var%lcontainer()) THEN
      CALL comin_plugin_finish_external("comin_var_to_3d",                      &
        &         " ERROR: Attempt to convert container variable into 3D field.")
    END IF

    CALL var%get_ptr(tmp_ptr)

    pos_jn = FINDLOC(var%dim_semantics(), COMIN_DIM_SEMANTICS_CONTAINER, DIM=1)
    SELECT CASE (pos_jn)
    CASE(1)
      slice => tmp_ptr(1, :, :, :, 1)
    CASE(2)
      slice => tmp_ptr(:, 1, :, :, 1)
    CASE(3)
      slice => tmp_ptr(:, :, 1, :, 1)
    CASE DEFAULT
      slice => tmp_ptr(:, :, :, 1, 1)
    END SELECT
  END SUBROUTINE comin_var_to_3d_dp

  !> Convenience operation for accessing 2D/3D fields.
  !! @ingroup plugin_interface
  !!
  !! Assumes that the last dimension is not used!
  SUBROUTINE comin_var_to_3d_sp(var, slice)
    CLASS(t_comin_var_handle), INTENT(IN)  :: var
    REAL(sp), POINTER :: slice(:,:,:)
    REAL(sp), POINTER :: tmp_ptr(:,:,:,:,:)
    INTEGER :: pos_jn

    ! this operation is invalid if the field is a container
    IF (var%lcontainer()) THEN
      CALL comin_plugin_finish_external("comin_var_to_3d",                      &
        &         " ERROR: Attempt to convert container variable into 3D field.")
    END IF

    CALL var%get_ptr(tmp_ptr)

    pos_jn = FINDLOC(var%dim_semantics(), COMIN_DIM_SEMANTICS_CONTAINER, DIM=1)
    SELECT CASE (pos_jn)
    CASE(1)
      slice => tmp_ptr(1, :, :, :, 1)
    CASE(2)
      slice => tmp_ptr(:, 1, :, :, 1)
    CASE(3)
      slice => tmp_ptr(:, :, 1, :, 1)
    CASE DEFAULT
      slice => tmp_ptr(:, :, :, 1, 1)
    END SELECT
  END SUBROUTINE comin_var_to_3d_sp

  !> Convenience operation for accessing 2D/3D fields.
  !! @ingroup plugin_interface
  !!
  !! Assumes that the last dimension is not used!
  SUBROUTINE comin_var_to_3d_i(var, slice)
    CLASS(t_comin_var_handle), INTENT(IN)  :: var
    INTEGER(C_INT), POINTER :: slice(:,:,:)
    INTEGER(C_INT), POINTER :: tmp_ptr(:,:,:,:,:)
    INTEGER :: pos_jn

    ! this operation is invalid if the field is a container
    IF (var%lcontainer()) THEN
      CALL comin_plugin_finish_external("comin_var_to_3d",                      &
        &         " ERROR: Attempt to convert container variable into 3D field.")
    END IF

    CALL var%get_ptr(tmp_ptr)

    pos_jn = FINDLOC(var%dim_semantics(), COMIN_DIM_SEMANTICS_CONTAINER, DIM=1)
    SELECT CASE (pos_jn)
    CASE(1)
      slice => tmp_ptr(1, :, :, :, 1)
    CASE(2)
      slice => tmp_ptr(:, 1, :, :, 1)
    CASE(3)
      slice => tmp_ptr(:, :, 1, :, 1)
    CASE DEFAULT
      slice => tmp_ptr(:, :, :, 1, 1)
    END SELECT
  END SUBROUTINE comin_var_to_3d_i

  SUBROUTINE comin_var_handle_set_cptr(var, cptr)
    TYPE(t_comin_var_handle), INTENT(INOUT) :: var
    TYPE(c_ptr), INTENT(IN) :: cptr

    var%var_item%cptr = cptr
#if defined(OPENACC)
    var%var_item%device_ptr = acc_deviceptr(cptr)
#endif

  END SUBROUTINE

END MODULE comin_variable_types
