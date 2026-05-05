!> @file comin_variable.F90
!! @brief Functions to modify and retrieve Variable definition
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_variable

  USE iso_c_binding,           ONLY: C_INT, C_PTR, C_LOC, C_NULL_PTR, C_F_POINTER, C_BOOL
  USE comin_ftnlist_ifc,       ONLY: comin_ftnlist_new, comin_ftnlist_push_back, comin_ftnlist_iterator_begin, &
    &                                comin_ftnlist_iterator_next, comin_ftnlist_iterator_value,                &
    &                                comin_ftnlist_iterator_delete, comin_ftnlist_is_end, comin_ftnlist_delete
  USE comin_keyval_ifc,        ONLY: comin_varmap_delete, comin_varmap_get, comin_varmap_put, &
    &                                comin_varmap_iterator_begin, comin_varmap_iterator_value, &
    &                                comin_varmap_iterator_next, comin_varmap_iterator_delete, &
    &                                comin_varmap_iterator_is_end
  USE comin_errhandler_constants, ONLY: COMIN_ERROR_POINTER_NOT_ASSOCIATED,                     &
    &                                   COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR,       &
    &                                   COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE,        &
    &                                   COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE,   &
    &                                   COMIN_ERROR_FIELD_NOT_ALLOCATED,                        &
    &                                   COMIN_ERROR_VAR_SYNC_DEVICE_MEM_NOT_ASSOCIATED,         &
    &                                   COMIN_ERROR_VAR_GET_NO_DEVICE,                          &
    &                                   COMIN_ERROR_VAR_SYNC_HALO_NOT_ASSOCIATED,               &
    &                                   COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR,      &
    &                                   COMIN_ERROR_VAR_GET_VARIABLE_NOT_FOUND,                 &
    &                                   COMIN_ERROR_VAR_GET_CONTAINER_CAN_NOT_HALO_SYNCHRONIZED, &
    &                                   COMIN_ERROR_VAR_GET_IRREGULAR_VAR_CAN_NOT_HALO_SYNCHRONIZED, &
    &                                   COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE
  USE comin_errhandler,        ONLY: comin_plugin_finish, comin_error_set
  USE comin_setup_constants,   ONLY: wp, EP_SECONDARY_CONSTRUCTOR,                              &
    &                                EP_DESTRUCTOR, COMIN_FLAG_DEVICE, COMIN_FLAG_SYNC_HALO, &
    &                                COMIN_VAR_DATATYPE_DOUBLE, COMIN_VAR_DATATYPE_FLOAT, &
    &                                COMIN_VAR_DATATYPE_INT, COMIN_DIM_SEMANTICS_UNDEF
  USE comin_state,             ONLY: state
  USE comin_c_utils,           ONLY: convert_f_string
  USE comin_variable_types,    ONLY: t_comin_var_descriptor, t_comin_var_item,               &
    &                                t_comin_var_descriptor_c, t_comin_var_context_item,     &
    &                                t_comin_var_handle, t_comin_request_item,                  &
    &                                comin_var_descr_match, comin_var_sync_device_mem_fct,   &
    &                                comin_var_sync_halo_fct, comin_var_ptr_init, comin_var_handle_set_cptr
  USE comin_descrdata_types,   ONLY: t_comin_descrdata_global
  USE comin_descrdata,         ONLY: comin_descrdata_get_global

  IMPLICIT NONE

  PRIVATE

  ! Public procedures, intention: called by host
  PUBLIC :: comin_var_list_finalize, comin_var_list_append
  PUBLIC :: comin_request_get_list, comin_var_request_list_finalize
  PUBLIC :: comin_var_descr_list_finalize
  ! Public procedures, intention: called by host and plugin
  PUBLIC :: comin_var_get_descr_list_head
  PUBLIC :: comin_var_descr_list_iterator_delete
  ! Public procedures, intention: called by plugin
  PUBLIC :: comin_var_request_add
  PUBLIC :: comin_var_get
  PUBLIC :: comin_var_get_descr_list_next, comin_var_get_descr_list_var_desc
  ! Public procedures, intention: for internal use
  PUBLIC :: comin_var_complete
  PUBLIC :: comin_var_get_from_exposed
  PUBLIC :: comin_var_set_sync_device_mem
  PUBLIC :: comin_var_set_sync_halo
  ! PUBLIC procedures only exposed to host model
  PUBLIC :: comin_var_set_cptr
  PUBLIC :: comin_var_is_used

#include "comin_global.inc"

CONTAINS

  !> Get first element of variable descriptor list.
  !! Returns a C-pointer that can be evaluated with
  !! the auxiliary function `comin_var_get_descr_list_var_desc`.
  FUNCTION comin_var_get_descr_list_head() RESULT(ptr_c) BIND(C)
    TYPE(c_ptr) :: ptr_c
    CALL comin_ftnlist_iterator_begin(state%comin_var_descr_list, ptr_c)
    IF (comin_ftnlist_is_end(state%comin_var_descr_list, ptr_c)) THEN
      ptr_c = C_NULL_PTR
    END IF
  END FUNCTION comin_var_get_descr_list_head

  !> Get next element of variable descriptor list.
  !! Returns a C-pointer that can be evaluated with
  !! the auxiliary function `comin_var_get_descr_list_var_desc`.
  !! Returns null pointer if end of list has been reached.
  FUNCTION comin_var_get_descr_list_next(current) RESULT(ptr_c) BIND(C)
    TYPE(C_PTR), INTENT(IN), VALUE :: current
    TYPE(C_PTR) :: ptr_c, tmp
    tmp = current
    CALL comin_ftnlist_iterator_next(tmp)
    IF (comin_ftnlist_is_end(state%comin_var_descr_list, tmp)) THEN
      CALL comin_ftnlist_iterator_delete(tmp)
      tmp = C_NULL_PTR
    END IF
    ptr_c = tmp
  END FUNCTION comin_var_get_descr_list_next

  !> @return head of variable list.
  !! @ingroup host_interface
  FUNCTION comin_request_get_list()  RESULT(ptr) BIND(C, name="comin_request_get_list")
    TYPE(c_ptr) :: ptr
    ptr = state%comin_var_request_list
  END FUNCTION comin_request_get_list

  !> Auxiliary function: Evaluates a list iterator of the
  !! variable descriptor list and returns the corresponding
  !! variable descriptor.
  !!
  !! This is the C-variant of the subroutine.
  SUBROUTINE comin_var_get_descr_list_var_desc_c(it, var_desc_out) BIND(C, name="comin_var_get_descr_list_var_desc")
    TYPE(C_PTR), INTENT(IN), VALUE :: it
    TYPE(t_comin_var_descriptor_c), INTENT(INOUT) :: var_desc_out
    !
    TYPE(t_comin_var_descriptor) :: var_desc_ftn

    CALL comin_var_get_descr_list_var_desc(it, var_desc_ftn)
    var_desc_out%id = var_desc_ftn%id
    CALL convert_f_string(var_desc_ftn%name, var_desc_out%name)
  END SUBROUTINE comin_var_get_descr_list_var_desc_c

  !> Auxiliary function: Evaluates a list iterator of the
  !! variable descriptor list and returns the corresponding
  !! variable descriptor.
  SUBROUTINE comin_var_get_descr_list_var_desc(it, var_desc_out)
    TYPE(C_PTR), INTENT(IN), VALUE :: it
    TYPE(t_comin_var_descriptor), INTENT(INOUT) :: var_desc_out
    !
    TYPE(c_ptr) :: cptr
    TYPE(t_comin_var_descriptor), POINTER :: item => NULL()

    CALL comin_ftnlist_iterator_value(it, cptr)
    CALL c_f_POINTER(cptr, item)
    IF (.NOT. ASSOCIATED(item)) THEN
      CALL comin_error_set(COMIN_ERROR_POINTER_NOT_ASSOCIATED); RETURN
    END IF
    var_desc_out = item
  END SUBROUTINE comin_var_get_descr_list_var_desc

  !> Delete list iterator.
  SUBROUTINE comin_var_descr_list_iterator_delete(it) BIND(C)
    TYPE(C_PTR), INTENT(INOUT) :: it
    CALL comin_ftnlist_iterator_delete(it)
  END SUBROUTINE comin_var_descr_list_iterator_delete

  ! destructor.
  SUBROUTINE comin_request_item_finalize(this)
    TYPE(t_comin_request_item), INTENT(INOUT) :: this
    CALL this%metadata%delete()
  END SUBROUTINE comin_request_item_finalize

  !> Append item to variable list.
  !! @ingroup host_interface
  SUBROUTINE comin_var_list_append(var_descr, cptr, device_ptr, &
       & array_shape, type_id, &
       & dim_semantics, &
       & lcontainer, ncontained, &
       & var_handle)
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descr
    TYPE(C_PTR), INTENT(IN)                  :: cptr, device_ptr
    INTEGER, INTENT(IN)                      :: array_shape(5), type_id
    INTEGER, INTENT(IN)                      :: dim_semantics(5)
    LOGICAL, INTENT(IN)                      :: lcontainer
    INTEGER, INTENT(IN)                      :: ncontained
    TYPE(t_comin_var_handle), INTENT(OUT)    :: var_handle
    !
    TYPE(t_comin_var_item),       POINTER :: var_item
    TYPE(t_comin_var_descriptor), POINTER :: var_descr_item

    ! first, add the descriptor to a separate list
    ! (the one that is also exposed to the plugins):
    ALLOCATE(var_descr_item)
    var_descr_item = var_descr
    CALL comin_ftnlist_push_back(state%comin_var_descr_list, c_loc(var_descr_item))

    ! add an entry to the other (internal) list of variables
    ! which contains a pointer to the above descriptor
    ALLOCATE(var_item)
    var_item%descriptor = var_descr
    var_item%cptr = cptr
    var_item%device_ptr = device_ptr
    var_item%array_shape = array_shape
    var_item%type_id = type_id
    var_item%dim_semantics = dim_semantics
    var_item%lcontainer = lcontainer
    var_item%ncontained = ncontained

    CALL var_item%metadata%create()
    CALL comin_varmap_put(state%comin_var_list, var_descr%name, var_descr%id, c_loc(var_item))

    var_handle = comin_var_ptr_init(var_item)

  END SUBROUTINE comin_var_list_append

  !> Destruct variable descriptor list, deallocate memory.
  SUBROUTINE comin_var_descr_list_finalize()
    ! local
    TYPE(c_ptr) :: it, cptr
    TYPE(t_comin_var_descriptor),  POINTER :: item

    CALL comin_ftnlist_iterator_begin(state%comin_var_descr_list, it)
    DO WHILE (.NOT. comin_ftnlist_is_end(state%comin_var_descr_list,it))
      CALL comin_ftnlist_iterator_value(it, cptr)
      CALL c_f_POINTER(cptr, item)
      DEALLOCATE(item)
      CALL comin_ftnlist_iterator_next(it)
    END DO
    CALL comin_ftnlist_iterator_delete(it)
    CALL comin_ftnlist_delete(state%comin_var_descr_list)
  END SUBROUTINE comin_var_descr_list_finalize

  !> Destruct variable list, deallocate memory.
  !! @ingroup host_interface
  SUBROUTINE comin_var_list_finalize()
    ! local
    TYPE(c_ptr) :: it, cptr
    TYPE(t_comin_var_item), POINTER :: item

    CALL comin_varmap_iterator_begin(state%comin_var_list, it)
    DO WHILE (.NOT. comin_varmap_iterator_is_end(state%comin_var_list, it))
      CALL comin_varmap_iterator_value(it, cptr)
      CALL c_f_POINTER(cptr, item)
      CALL item%metadata%delete()
      DEALLOCATE(item)
      CALL comin_varmap_iterator_next(it)
    END DO
    CALL comin_varmap_iterator_delete(it)
    CALL comin_varmap_delete(state%comin_var_list)
  END SUBROUTINE comin_var_list_finalize

  !> Destruct variable request list, deallocate memory.
  !! @ingroup host_interface
  SUBROUTINE comin_var_request_list_finalize()
    ! local
    TYPE(c_ptr) :: it, cptr
    TYPE(t_comin_request_item), POINTER :: item

    CALL comin_ftnlist_iterator_begin(state%comin_var_request_list, it)
    DO WHILE (.NOT. comin_ftnlist_is_end(state%comin_var_request_list,it))
      CALL comin_ftnlist_iterator_value(it, cptr)
      CALL c_f_POINTER(cptr, item)
      CALL comin_request_item_finalize(item)
      CALL comin_ftnlist_iterator_next(it)
    END DO
    CALL comin_ftnlist_iterator_delete(it)
    CALL comin_ftnlist_delete(state%comin_var_request_list)
  END SUBROUTINE comin_var_request_list_finalize

  FUNCTION comin_var_get_c(context_len, context, var_descriptor, flag) &
       & RESULT(var_pointer) &
       & BIND(C, name="comin_var_get")
    INTEGER(c_int),VALUE,           INTENT(IN) :: context_len
    INTEGER(c_int),                 INTENT(IN) :: context(context_len)
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN) :: var_descriptor
    INTEGER(c_int), VALUE,          INTENT(IN) :: flag
    TYPE(c_ptr)                                :: var_pointer
    TYPE(t_comin_var_item), POINTER            :: var_item => NULL()

    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran

    var_pointer = C_NULL_PTR

    var_descriptor_fortran = t_comin_var_descriptor(var_descriptor)

    CALL comin_var_get_internal(context, var_descriptor_fortran, flag, var_item)
    IF(ASSOCIATED(var_item)) var_pointer = C_LOC(var_item)
  END FUNCTION comin_var_get_c

  FUNCTION comin_var_get_ptr(handle) &
       & RESULT(dataptr)                          &
       & BIND(C, NAME="comin_var_get_ptr")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: dataptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
      dataptr = C_NULL_PTR
    ELSE
      dataptr  = p%cptr
    END IF
  END FUNCTION comin_var_get_ptr

  FUNCTION comin_var_get_ptr_double(handle) &
       & RESULT(dataptr)                    &
       & BIND(C, NAME="comin_var_get_ptr_double")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: dataptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    dataptr = C_NULL_PTR
    CALL C_F_POINTER(handle, p)
    IF (ASSOCIATED(p)) THEN
      IF(p%type_id /= COMIN_VAR_DATATYPE_DOUBLE) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE); RETURN
      ENDIF
      dataptr  = p%cptr
    END IF
  END FUNCTION comin_var_get_ptr_double

  FUNCTION comin_var_get_ptr_float(handle) &
       & RESULT(dataptr)                 &
       & BIND(C, NAME="comin_var_get_ptr_float")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: dataptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    dataptr = C_NULL_PTR
    CALL C_F_POINTER(handle, p)
    IF (ASSOCIATED(p)) THEN
      IF(p%type_id /= COMIN_VAR_DATATYPE_FLOAT) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE); RETURN
      ENDIF
      dataptr  = p%cptr
    END IF
  END FUNCTION comin_var_get_ptr_float

  FUNCTION comin_var_get_ptr_int(handle) &
       & RESULT(dataptr)                 &
       & BIND(C, NAME="comin_var_get_ptr_int")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: dataptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    dataptr = C_NULL_PTR
    CALL C_F_POINTER(handle, p)
    IF (ASSOCIATED(p)) THEN
      IF(p%type_id /= COMIN_VAR_DATATYPE_INT) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE); RETURN
      ENDIF
      dataptr  = p%cptr
    END IF
  END FUNCTION comin_var_get_ptr_int

  FUNCTION comin_var_get_device_ptr(handle)          &
       & RESULT(device_ptr)                          &
       & BIND(C, NAME="comin_var_get_device_ptr")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: device_ptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    device_ptr = p%device_ptr
  END FUNCTION comin_var_get_device_ptr

  FUNCTION comin_var_get_device_ptr_double(handle)          &
       & RESULT(device_ptr)                          &
       & BIND(C, NAME="comin_var_get_device_ptr_double")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: device_ptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()

    device_ptr = C_NULL_PTR
    CALL C_F_POINTER(handle, p)
    IF (ASSOCIATED(p)) THEN
      IF(p%type_id /= COMIN_VAR_DATATYPE_DOUBLE) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE); RETURN
      ENDIF
      device_ptr = p%device_ptr
    ENDIF
  END FUNCTION comin_var_get_device_ptr_double

  FUNCTION comin_var_get_device_ptr_float(handle)          &
       & RESULT(device_ptr)                          &
       & BIND(C, NAME="comin_var_get_device_ptr_float")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: device_ptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()

    device_ptr = C_NULL_PTR
    CALL C_F_POINTER(handle, p)
    IF (ASSOCIATED(p)) THEN
      IF(p%type_id /= COMIN_VAR_DATATYPE_FLOAT) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE); RETURN
      ENDIF
      device_ptr = p%device_ptr
    ENDIF
  END FUNCTION comin_var_get_device_ptr_float

  FUNCTION comin_var_get_device_ptr_int(handle)          &
       & RESULT(device_ptr)                          &
       & BIND(C, NAME="comin_var_get_device_ptr_int")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: device_ptr
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()

    device_ptr = C_NULL_PTR
    CALL C_F_POINTER(handle, p)
    IF (ASSOCIATED(p)) THEN
      IF(p%type_id /= COMIN_VAR_DATATYPE_INT) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_WRONG_TYPE); RETURN
      ENDIF
      device_ptr = p%device_ptr
    ENDIF
  END FUNCTION comin_var_get_device_ptr_int

  SUBROUTINE comin_var_get_shape(handle, data_shape) &
       & BIND(C, NAME="comin_var_get_shape")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    INTEGER(C_INT), INTENT(INOUT)       :: data_shape(5)
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
      CALL comin_error_set(COMIN_ERROR_POINTER_NOT_ASSOCIATED); RETURN
    ELSE
      data_shape = p%array_shape
    END IF
  END SUBROUTINE comin_var_get_shape

  SUBROUTINE comin_var_get_dim_semantics(handle, dim_semantics) &
       & BIND(C, NAME="comin_var_get_dim_semantics")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT)       :: dim_semantics(5)
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
      CALL comin_error_set(COMIN_ERROR_POINTER_NOT_ASSOCIATED); RETURN
    ELSE
      dim_semantics = p%dim_semantics
    END IF
  END SUBROUTINE comin_var_get_dim_semantics

  SUBROUTINE comin_var_get_ncontained(handle, ncontained) &
       & BIND(C, NAME="comin_var_get_ncontained")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT)       :: ncontained
    !
    TYPE(t_comin_var_item), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
      CALL comin_error_set(COMIN_ERROR_POINTER_NOT_ASSOCIATED); RETURN
    ELSE
      ! Convert to C dimension index
      ncontained = p%ncontained - 1
    END IF
  END SUBROUTINE comin_var_get_ncontained

  SUBROUTINE comin_var_get_descriptor(handle, descr) &
       & BIND(C, NAME="comin_var_get_descriptor")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(t_comin_var_descriptor_c), INTENT(INOUT) :: descr

    TYPE(t_comin_var_item), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
      CALL comin_error_set(COMIN_ERROR_POINTER_NOT_ASSOCIATED); RETURN
    ELSE
      CALL convert_f_string(p%descriptor%name, descr%name)
      descr%id   = p%descriptor%id
    END IF
  END SUBROUTINE comin_var_get_descriptor

  !> Request a pointer to an ICON variable in context(s).
  !! @ingroup plugin_interface
  SUBROUTINE comin_var_get(context, var_descriptor, flag, var_ptr)
    INTEGER,                      INTENT(IN)  :: context(:)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor
    INTEGER,                      INTENT(IN)  :: flag
    TYPE(t_comin_var_handle),        INTENT(OUT) :: var_ptr
    ! local
    TYPE(t_comin_var_item), POINTER          :: var_item => NULL()

    CALL comin_var_get_internal(context, var_descriptor, flag, var_item)
    var_ptr = comin_var_ptr_init(var_item)
  END SUBROUTINE comin_var_get

  !> get pointer to a variable exposed by ICON
  FUNCTION comin_var_get_from_exposed(var_descriptor)  RESULT(comin_get_var)
    TYPE(t_comin_var_item), POINTER :: comin_get_var
    TYPE (t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    !

    CALL C_F_POINTER(comin_varmap_get(state%comin_var_list, var_descriptor%name, var_descriptor%id), comin_get_var)
  END FUNCTION comin_var_get_from_exposed

  !> get pointer to a variable according to context and descriptor
  FUNCTION comin_var_get_by_context(context, plugin_id, var_descriptor)  RESULT(comin_get_var)
    TYPE(t_comin_var_context_item), POINTER   :: comin_get_var
    INTEGER, INTENT(IN)                       :: context, plugin_id
    TYPE (t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    ! local
    TYPE(c_ptr) :: it, cptr
    TYPE(t_comin_var_context_item), POINTER :: item

    comin_get_var => NULL()
    IF (.NOT. ALLOCATED(state%comin_var_list_context)) RETURN
    ASSOCIATE(var_list => state%comin_var_list_context(context, plugin_id)%var_list)
      CALL comin_ftnlist_iterator_begin(var_list, it)
      DO WHILE (.not. comin_ftnlist_is_end(var_list,it))
        ! test if already registered for context
        CALL comin_ftnlist_iterator_value(it, cptr)
        CALL c_f_pointer(cptr, item)
        IF (comin_var_descr_match(item%var_item%descriptor, var_descriptor)) THEN
          comin_get_var => item
          EXIT
        END IF
        CALL comin_ftnlist_iterator_next(it)
      END DO
      CALL comin_ftnlist_iterator_delete(it)
    END ASSOCIATE
  END FUNCTION comin_var_get_by_context

  SUBROUTINE comin_var_request_add_c(var_descriptor, lmodexclusive) &
    &  BIND(C, name="comin_var_request_add")
    TYPE (t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor
    LOGICAL(C_BOOL), VALUE,          INTENT(IN)  :: lmodexclusive
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran = t_comin_var_descriptor(var_descriptor)
    CALL comin_var_request_add(var_descriptor_fortran, LOGICAL(lmodexclusive))
  END SUBROUTINE comin_var_request_add_c

  !> By calling this subroutine inside the primary constructor, 3rd
  !> party plugins may request the creation of additional variables.
  !! @ingroup plugin_interface
  !!
  !!  Note: The lmodexclusive argument provides the information if this
  !!        variable is exclusive to the calling plugin.
  !!
  !!  Note: If a 3rd party plugin requests the creation of a variable
  !!        through this subroutine, it is still not guaranteed that
  !!        this variable is actually created! It might be skipped
  !!        due to inconsistencies, it could be a duplicate
  !!        etc. Therefore, 3rd party plugins still have to evaluate
  !!        the return code of `comin_var_request_add`.
  !!
  SUBROUTINE comin_var_request_add(var_descriptor, lmodexclusive)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor
    LOGICAL,                       INTENT(IN)  :: lmodexclusive
    ! local
    TYPE(t_comin_descrdata_global), POINTER :: comin_global
    TYPE (t_comin_var_descriptor) :: var_descriptor_domain
    INTEGER :: domain_id

    comin_global => comin_descrdata_get_global()

    IF (state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR); RETURN
    ENDIF

    IF (var_descriptor%id == -1) THEN
      comin_global => comin_descrdata_get_global()
      IF (.NOT. ASSOCIATED(comin_global)) CALL comin_plugin_finish("variable ", "global data missing")

      DO domain_id = 1, comin_global%n_dom
        var_descriptor_domain    = var_descriptor
        var_descriptor_domain%id = domain_id
        CALL comin_var_request_add_element(var_descriptor_domain, lmodexclusive)
      ENDDO
    ELSE
      CALL comin_var_request_add_element(var_descriptor, lmodexclusive)
    ENDIF

  CONTAINS

    SUBROUTINE comin_var_request_add_element(var_descriptor, lmodexclusive)
      TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor
      LOGICAL,                       INTENT(IN)  :: lmodexclusive
      !
      TYPE(c_ptr) :: it, cptr
      TYPE(t_comin_request_item), POINTER :: var_list_request_element, comin_request_item

      ! check if requested variable already requested or if modexlusive conflicts exist
      ! first find the variable in list of all ICON variables and set the pointer
      CALL comin_ftnlist_iterator_begin(state%comin_var_request_list, it)
      DO WHILE (.NOT. comin_ftnlist_is_end(state%comin_var_request_list,it))
        CALL comin_ftnlist_iterator_value(it, cptr)
        CALL c_f_POINTER(cptr, var_list_request_element)
        IF (comin_var_descr_match(var_list_request_element%descriptor, var_descriptor)) THEN
          !> first criterion for abort: variable exists and was requested exclusively
          IF (var_list_request_element%lmodexclusive) THEN
            CALL comin_error_set(COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE); RETURN
            !> second criterion for abort: variable exists and now requested exclusively
          ELSEIF (lmodexclusive) THEN
            CALL comin_error_set(COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE); RETURN
            !> if existing but no conflicts with exclusiveness: expand moduleID information
          ELSE
            IF (.NOT. ALLOCATED(var_list_request_element%moduleID)) THEN
              ! if not allocated something went wrong before (should not happen)
              CALL comin_error_set(COMIN_ERROR_FIELD_NOT_ALLOCATED); RETURN
            ELSE
              var_list_request_element%moduleID = [var_list_request_element%moduleID(:), &
                &                                  state%current_plugin%id]
            END IF
            RETURN
          END IF
        END IF
        CALL comin_ftnlist_iterator_next(it)
      END DO
      CALL comin_ftnlist_iterator_delete(it)

      ! register new variable request
      ASSOCIATE( var_list => state%comin_var_request_list)
        ALLOCATE(comin_request_item)
        comin_request_item%descriptor =  var_descriptor
        comin_request_item%lmodexclusive = lmodexclusive
        comin_request_item%moduleID = [state%current_plugin%id]
        CALL comin_request_item%metadata%create()
        CALL comin_ftnlist_push_back(var_list, c_loc(comin_request_item))
      END ASSOCIATE
    END SUBROUTINE comin_var_request_add_element

  END SUBROUTINE comin_var_request_add

  ! Internal subroutine. Consistency checks and similar operations,
  ! done after primary constructors.
  SUBROUTINE comin_var_complete()
    INTEGER :: i,j

    ALLOCATE(state%comin_var_list_context(EP_DESTRUCTOR, state%num_plugins))
    DO i=1,size(state%comin_var_list_context,1)
      DO j=1,size(state%comin_var_list_context,2)
        CALL comin_ftnlist_new(state%comin_var_list_context(i,j)%var_list)
      END DO
    END DO

  END SUBROUTINE comin_var_complete

  SUBROUTINE comin_var_set_sync_device_mem(sync_device_mem)
    PROCEDURE(comin_var_sync_device_mem_fct) :: sync_device_mem

    state%sync_device_mem => sync_device_mem
    IF (.NOT. ASSOCIATED(state%sync_device_mem)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_SYNC_DEVICE_MEM_NOT_ASSOCIATED); RETURN
    END IF
  END SUBROUTINE comin_var_set_sync_device_mem

  SUBROUTINE comin_var_set_sync_halo(sync_halo)
    PROCEDURE(comin_var_sync_halo_fct)  :: sync_halo
    state%sync_halo => sync_halo
    IF (.NOT. ASSOCIATED(state%sync_halo)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_SYNC_HALO_NOT_ASSOCIATED); RETURN
    END IF
  END SUBROUTINE comin_var_set_sync_halo

  SUBROUTINE comin_var_get_internal(context, var_descriptor, flag, var_item)
    INTEGER,                      INTENT(IN) :: context(:)
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    INTEGER,                      INTENT(IN) :: flag
    TYPE(t_comin_var_item), POINTER :: var_item
    ! local
    TYPE(t_comin_var_context_item), POINTER :: item
    INTEGER :: ic

    ! Routine should only be called during secondary constructor
    IF ((.NOT. state%l_primary_done) .OR. &
     &   state%current_ep > EP_SECONDARY_CONSTRUCTOR) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR); RETURN
    END IF

    ! device pointers can only be accessed if a device is available
    IF ((.NOT. state%comin_descrdata_global%has_device) .AND. &
         & IAND(flag, COMIN_FLAG_DEVICE) /= 0) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_GET_NO_DEVICE); RETURN
    ENDIF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_GET_VARIABLE_NOT_FOUND); RETURN
    ENDIF
    ! a container can not halo synchronized
    IF ((var_item%lcontainer) .AND. &
         & IAND(flag, COMIN_FLAG_SYNC_HALO) /= 0) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_GET_CONTAINER_CAN_NOT_HALO_SYNCHRONIZED); RETURN
    ENDIF
    ! an irregular var can not be halo synchronized
    IF ((ANY(var_item%dim_semantics == COMIN_DIM_SEMANTICS_UNDEF)) .AND. &
         & IAND(flag, COMIN_FLAG_SYNC_HALO) /= 0) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_GET_IRREGULAR_VAR_CAN_NOT_HALO_SYNCHRONIZED); RETURN
    ENDIF

    DO ic = 1, SIZE(context)
      ! ignore EP_SECONDARY_CONSTRUCTOR for var_list
      IF (context(ic) == EP_SECONDARY_CONSTRUCTOR) CYCLE
      item => comin_var_get_by_context(context(ic), state%current_plugin%id, var_item%descriptor)
      IF (.NOT. ASSOCIATED(item)) THEN
        ! not in context list: register variable, set access flag
        ASSOCIATE(var_list => state%comin_var_list_context(context(ic) , state%current_plugin%id)%var_list)
          ALLOCATE(item)
          item = t_comin_var_context_item( var_item = var_item, &
          &                                access_flag = flag)
          CALL comin_ftnlist_push_back(var_list, c_loc(item))
        END ASSOCIATE
      END IF
    END DO
  END SUBROUTINE comin_var_get_internal

  SUBROUTINE comin_var_set_cptr(var, cptr)
    TYPE(t_comin_var_handle), INTENT(INOUT) :: var
    TYPE(c_ptr), INTENT(IN)               :: cptr
    !

    CALL comin_var_handle_set_cptr(var, cptr)

  END SUBROUTINE comin_var_set_cptr

  !> Check if a variable is actually used by any comin plugin.
  LOGICAL FUNCTION comin_var_is_used(var_descriptor)
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descriptor

    TYPE(t_comin_var_context_item), POINTER :: item
    TYPE(c_ptr) :: list, it, p
    INTEGER :: ep, id

    comin_var_is_used = .FALSE.

    DO ep = 1, SIZE(state%comin_var_list_context, 1)
      DO id = 1, SIZE(state%comin_var_list_context, 2)
        list = state%comin_var_list_context(ep, id)%var_list
        CALL comin_ftnlist_iterator_begin(list, it)
        DO WHILE(.NOT. comin_ftnlist_is_end(list, it))
          CALL comin_ftnlist_iterator_value(it, p)
          CALL c_f_pointer(p, item)

          IF (comin_var_descr_match(item%var_item%descriptor, var_descriptor)) THEN
            comin_var_is_used = .TRUE.
            CALL comin_ftnlist_iterator_delete(it)
            RETURN
          END IF

          CALL comin_ftnlist_iterator_next(it)
        END DO
        CALL comin_ftnlist_iterator_delete(it)
      END DO
    END DO
  END FUNCTION

END MODULE comin_variable
