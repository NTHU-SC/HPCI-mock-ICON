!> @file comin_keyval_ifc.F90
!! @brief Interfaces for key-value maps based on std::unordered_map.
!
!  @authors 06/2025 :: ICON Community Interface  <icon@dwd.de>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_keyval_ifc
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_int64_t, c_ptr, c_char, c_bool, c_double

  PRIVATE

  PUBLIC :: comin_keyval_create_c
  PUBLIC :: comin_keyval_delete_c
  PUBLIC :: comin_keyval_query_c
  PUBLIC :: comin_keyval_set_int_c
  PUBLIC :: comin_keyval_get_int_c
  PUBLIC :: comin_keyval_set_double_c
  PUBLIC :: comin_keyval_get_double_c
  PUBLIC :: comin_keyval_set_char_c
  PUBLIC :: comin_keyval_get_char_c
  PUBLIC :: comin_keyval_set_bool_c
  PUBLIC :: comin_keyval_get_bool_c
  PUBLIC :: comin_keyval_iterator_begin_c
  PUBLIC :: comin_keyval_iterator_end_c
  PUBLIC :: comin_keyval_iterator_get_key_c
  PUBLIC :: comin_keyval_iterator_next_c
  PUBLIC :: comin_keyval_iterator_delete_c
  PUBLIC :: comin_keyval_iterator_compare_c

  PUBLIC :: comin_varmap_new
  PUBLIC :: comin_varmap_delete

  PUBLIC :: comin_varmap_get
  PUBLIC :: comin_varmap_put

  PUBLIC :: comin_varmap_iterator_begin
  PUBLIC :: comin_varmap_iterator_delete
  PUBLIC :: comin_varmap_iterator_next
  PUBLIC :: comin_varmap_iterator_is_end
  PUBLIC :: comin_varmap_iterator_value

  INTERFACE
    SUBROUTINE comin_keyval_create_c(comin_keyval_c) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr) :: comin_keyval_c
    END SUBROUTINE comin_keyval_create_c

    SUBROUTINE comin_keyval_delete_c(comin_keyval_c) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: comin_keyval_c
    END SUBROUTINE comin_keyval_delete_c

    SUBROUTINE comin_keyval_query_c(ckey, idx, keyval_c) BIND(C)
      IMPORT c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      INTEGER(KIND=c_int), INTENT(out) :: idx
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_query_c

    SUBROUTINE comin_keyval_set_int_c(ckey,  val, keyval_c) BIND(C)
      IMPORT c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      INTEGER(KIND=c_int), VALUE, INTENT(in) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_int_c

    SUBROUTINE comin_keyval_get_int_c(ckey, val, keyval_c) BIND(C)
      IMPORT c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      INTEGER(KIND=c_int), INTENT(out) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_int_c

    SUBROUTINE comin_keyval_set_double_c(ckey,  val, keyval_c) BIND(C)
      IMPORT c_ptr, c_double
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      REAL(KIND=c_double), VALUE, INTENT(in) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_double_c

    SUBROUTINE comin_keyval_get_double_c(ckey, val, keyval_c) BIND(C)
      IMPORT c_ptr, c_double
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      REAL(KIND=c_double), INTENT(out) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_double_c

    SUBROUTINE comin_keyval_set_char_c(ckey,  cval, keyval_c) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      TYPE(c_ptr), VALUE, INTENT(in) :: cval
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_char_c

    SUBROUTINE comin_keyval_get_char_c(ckey,  cval, keyval_c) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      TYPE(c_ptr) :: cval
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_char_c

    SUBROUTINE comin_keyval_set_bool_c(ckey,  val, keyval_c) BIND(C)
      IMPORT c_ptr, c_bool
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_bool_c

    SUBROUTINE comin_keyval_get_bool_c(ckey, val, keyval_c) BIND(C)
      IMPORT c_ptr, c_bool
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      LOGICAL(KIND=c_bool), INTENT(out) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_bool_c

    SUBROUTINE comin_keyval_iterator_begin_c(keyval_c, iterator) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
      TYPE(c_ptr) :: iterator
    END SUBROUTINE comin_keyval_iterator_begin_c

    SUBROUTINE comin_keyval_iterator_end_c(keyval_c, iterator) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
      TYPE(c_ptr) :: iterator
    END SUBROUTINE comin_keyval_iterator_end_c

    FUNCTION comin_keyval_iterator_get_key_c(iterator) RESULT(ckey) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator
      TYPE(c_ptr) :: ckey
    END FUNCTION comin_keyval_iterator_get_key_c

    SUBROUTINE comin_keyval_iterator_next_c(iterator) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator
    END SUBROUTINE comin_keyval_iterator_next_c

    SUBROUTINE comin_keyval_iterator_delete_c(iterator) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator
    END SUBROUTINE comin_keyval_iterator_delete_c

    FUNCTION comin_keyval_iterator_compare_c(iterator1, iterator2) BIND(C)
      IMPORT c_ptr, c_bool
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator1, iterator2
      LOGICAL(KIND=c_bool) :: comin_keyval_iterator_compare_c
    END FUNCTION comin_keyval_iterator_compare_c
  END INTERFACE

  ! Varmap functions
  INTERFACE
    SUBROUTINE comin_varmap_new(map) BIND(C, name='comin_varmap_new_c')
      IMPORT c_ptr
      TYPE(c_ptr), INTENT(OUT) :: map
    END SUBROUTINE
    SUBROUTINE comin_varmap_delete(map) BIND(C, name='comin_varmap_delete_c')
      IMPORT c_ptr
      TYPE(c_ptr), VALUE :: map
    END SUBROUTINE
    TYPE(c_ptr) FUNCTION comin_varmap_get_c (map, name, len, id) BIND(C)
      IMPORT c_ptr, c_int64_t, c_char, c_int
      TYPE(c_ptr), VALUE, INTENT(IN) :: map
      INTEGER(c_int64_t), VALUE, INTENT(IN) :: len
      CHARACTER(kind=c_char), INTENT(IN) :: name(len)
      INTEGER(c_int), VALUE, INTENT(IN) :: id
    END FUNCTION
    SUBROUTINE comin_varmap_put_c (map, name, len, id, ptr) BIND(C)
      IMPORT c_ptr, c_int64_t, c_char, c_int
      TYPE(c_ptr), VALUE :: map
      INTEGER(c_int64_t), VALUE, INTENT(IN) :: len
      CHARACTER(kind=c_char), INTENT(IN) :: name(len)
      INTEGER(c_int), VALUE, INTENT(IN) :: id
      TYPE(c_ptr), VALUE, INTENT(IN) :: ptr
    END SUBROUTINE
    SUBROUTINE comin_varmap_iterator_begin (map, it) BIND(C, name='comin_varmap_iterator_begin_c')
      IMPORT c_ptr
      TYPE(c_ptr), VALUE :: map
      TYPE(c_ptr), INTENT(OUT) :: it
    END SUBROUTINE
    SUBROUTINE comin_varmap_iterator_delete (it) BIND(C, name='comin_varmap_iterator_delete_c')
      IMPORT c_ptr
      TYPE(c_ptr), VALUE :: it
    END SUBROUTINE
    SUBROUTINE comin_varmap_iterator_next (it) BIND(C, name='comin_varmap_iterator_next_c')
      IMPORT c_ptr
      TYPE(c_ptr), VALUE :: it
    END SUBROUTINE
    LOGICAL(c_bool) FUNCTION comin_varmap_iterator_is_end (map, it) BIND(C, name='comin_varmap_iterator_is_end_c')
      IMPORT c_ptr, c_bool
      TYPE(c_ptr), VALUE, INTENT(IN) :: map
      TYPE(c_ptr), VALUE, INTENT(IN) :: it
    END FUNCTION
    SUBROUTINE comin_varmap_iterator_value (it, ptr) BIND(C, name='comin_varmap_iterator_value_c')
      IMPORT c_ptr
      TYPE(c_ptr), VALUE, INTENT(IN) :: it
      TYPE(c_ptr), INTENT(OUT) :: ptr
    END SUBROUTINE
  END INTERFACE

CONTAINS

  SUBROUTINE comin_varmap_put (map, name, id, ptr)
    TYPE(c_ptr), INTENT(INOUT) :: map
    CHARACTER(kind=c_char, len=*), INTENT(IN) :: name
    INTEGER(c_int), INTENT(IN) :: id
    TYPE(c_ptr), INTENT(IN) :: ptr

    CALL comin_varmap_put_c(map, name, LEN(name, kind=c_int64_t), id, ptr)
  END SUBROUTINE

  TYPE(c_ptr) FUNCTION comin_varmap_get (map, name, id)
    TYPE(c_ptr), INTENT(IN) :: map
    CHARACTER(kind=c_char, len=*), INTENT(IN) :: name
    INTEGER(c_int), INTENT(IN) :: id

    comin_varmap_get = comin_varmap_get_c(map, name, LEN(name, kind=c_int64_t), id)
  END FUNCTION
END MODULE
