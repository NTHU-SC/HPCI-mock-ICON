!> @file comin_ftnlist_ifc.F90
!! @brief Interfaces for list implementation based on std::list.
!
!  @authors 10/2023 :: ICON Community Interface  <icon@dwd.de>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
module comin_ftnlist_ifc
  use iso_c_binding
  implicit none

  PUBLIC

  interface
    subroutine comin_ftnlist_new(ptr) bind(c, name="comin_ftnlist_new")
      use iso_c_binding
      type(c_ptr), intent(inout) :: ptr
    end subroutine comin_ftnlist_new

    subroutine comin_ftnlist_delete(ptr) bind(c, name="comin_ftnlist_delete")
      use iso_c_binding
      type(c_ptr), intent(inout) :: ptr
    end subroutine comin_ftnlist_delete

    subroutine comin_ftnlist_push_back(listptr, ptr) bind(c, name="comin_ftnlist_push_back")
      use iso_c_binding
      type(c_ptr), intent(in), value :: listptr, ptr
    end subroutine comin_ftnlist_push_back

    subroutine comin_ftnlist_iterator_begin(listptr, itptr) bind(c, name="comin_ftnlist_iterator_begin")
      use iso_c_binding
      type(c_ptr), intent(in), value :: listptr
      type(c_ptr), intent(inout) :: itptr
    end subroutine comin_ftnlist_iterator_begin

    subroutine comin_ftnlist_iterator_next(ptr) bind(c, name="comin_ftnlist_iterator_next")
      use iso_c_binding
      type(c_ptr), intent(in), value :: ptr
    end subroutine comin_ftnlist_iterator_next

    subroutine comin_ftnlist_iterator_delete(ptr) bind(c, name="comin_ftnlist_iterator_delete")
      use iso_c_binding
      type(c_ptr), intent(inout) :: ptr
    end subroutine comin_ftnlist_iterator_delete

    subroutine comin_ftnlist_iterator_value(listptr, ptr) bind(c, name="comin_ftnlist_iterator_value")
      use iso_c_binding
      type(c_ptr), intent(in), value :: listptr
      type(c_ptr), intent(inout) :: ptr
    end subroutine comin_ftnlist_iterator_value

    logical(c_bool) function comin_ftnlist_is_end(listptr, ptr) bind(c, name="comin_ftnlist_is_end")
      use iso_c_binding
      type(c_ptr), intent(in), value :: listptr, ptr
    end function comin_ftnlist_is_end

  end interface
end module comin_ftnlist_ifc
