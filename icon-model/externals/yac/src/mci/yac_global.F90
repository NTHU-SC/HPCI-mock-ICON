! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

module mo_yac_iso_c_helpers

  use, intrinsic :: iso_c_binding, only: c_null_char

  private
  public :: yac_internal_cptr2char
  public :: yac_dble2cptr

  contains

#ifndef YAC_CODE_COVERAGE_TEST
#define YAC_FASSERT(exp, msg) IF (.NOT. exp) call yac_abort_message(TRIM(msg) // c_null_char, TRIM(__FILE__) // c_null_char, __LINE__)
#define YAC_CHECK_STRING_LEN(routine, str) YAC_FASSERT(LEN_TRIM(str) < YAC_MAX_CHARLEN, "ERROR(" // TRIM(routine) // "): string '" // TRIM(str) // "' exceeds length of YAC_MAX_CHARLEN")
#else
#define YAC_FASSERT(exp, msg)
#define YAC_CHECK_STRING_LEN(routine, str)
#endif

  !> Convertes a C-pointer to a C-string to a Fortran string
  !! @param[in] cptr C-cointer to C-string
  !! @returns Fortran string
  function yac_internal_cptr2char( cptr ) result (string)

    USE, intrinsic :: iso_c_binding, only: c_ptr, c_char, &
        c_f_pointer,c_size_t

    implicit none

    TYPE(c_ptr), intent(in) :: cptr
    CHARACTER(len=:), allocatable :: string
    CHARACTER(kind=c_char), dimension(:), pointer :: chars
    INTEGER(kind=c_size_t) :: i, strlen

    interface
      function strlen_c(str_ptr) bind ( C, name = "strlen" ) result(len)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value      :: str_ptr
        integer(kind=c_size_t)  :: len
      end function strlen_c
    end interface

    strlen = strlen_c(cptr)
    CALL c_f_pointer(cptr, chars, [ strlen ])
    ALLOCATE(character(len=strlen) :: string)
    DO i=1,strlen
      string(i:i) = chars(i)
    END DO
  end function  yac_internal_cptr2char

  function yac_dble2cptr(routine, ptr_name, dble_ptr)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_ptr

    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: ptr_name
    type(yac_dble_ptr), intent(in) :: dble_ptr
    type(c_ptr) :: yac_dble2cptr

    if (SIZE(dble_ptr%p) > 0) then
        YAC_FASSERT(is_contiguous(dble_ptr%p), "ERROR(" // TRIM(routine) // "): " // TRIM(ptr_name) // " is not contiguous")
        yac_dble2cptr = c_loc(dble_ptr%p(1))
    else
        yac_dble2cptr = c_null_ptr
    endif
  end function yac_dble2cptr

end module mo_yac_iso_c_helpers

module mo_yac_real_to_dble_utils

  public :: send_field_to_dble, &
            send_field_to_dble_single, &
            send_field_to_dble_ptr, &
            recv_field_to_dble, &
            recv_field_to_dble_ptr, &
            recv_field_from_dble, &
            recv_field_from_dble_ptr

contains

  subroutine send_field_to_dble(field_id,        &
                                nbr_hor_points,  &
                                nbr_pointsets,   &
                                collection_size, &
                                send_field,      &
                                send_field_dble, &
                                send_frac_mask,  &
                                send_frac_mask_dble)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_put_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_put_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_put_mask_c2f_c

      end function yac_get_field_put_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: nbr_pointsets
    integer, intent (in)  :: collection_size
    real, intent (in)     :: send_field(nbr_hor_points, &
                                        nbr_pointsets,  &
                                        collection_size)
    double precision, intent (out) :: send_field_dble(nbr_hor_points, &
                                                      nbr_pointsets,  &
                                                      collection_size)
    real, optional, intent (in) :: send_frac_mask(nbr_hor_points, &
                                                  nbr_pointsets,  &
                                                  collection_size)
    double precision, optional, intent (out) :: send_frac_mask_dble(nbr_hor_points, &
                                                                    nbr_pointsets,  &
                                                                    collection_size)

    integer :: i, j, k
    type(c_ptr) :: put_mask_
    type(c_ptr), pointer :: put_mask(:)
    integer(kind=c_int), pointer :: pointset_put_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_SOURCE) then

      put_mask_ = yac_get_field_put_mask_c2f_c(field_id)
      if (c_associated(put_mask_)) then
        call c_f_pointer(put_mask_, put_mask, shape=[nbr_pointsets])
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            call c_f_pointer( &
                  put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
            do k = 1, nbr_hor_points
              if (pointset_put_mask(k) /= 0) then
                send_field_dble(k, j, i) = dble(send_field(k, j, i))
              else
                send_field_dble(k, j, i) = 0d0
              end if
            end do
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            do j = 1, nbr_pointsets
              call c_f_pointer( &
                    put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
              do k = 1, nbr_hor_points
                if (pointset_put_mask(k) /= 0) then
                  send_frac_mask_dble(k, j, i) = dble(send_frac_mask(k, j, i))
                else
                  send_frac_mask_dble(k, j, i) = 0d0
                end if
              end do
            end do
          end do
        end if
      else
        send_field_dble = dble(send_field)
        if (present(send_frac_mask)) then
          send_frac_mask_dble = dble(send_frac_mask)
        end if
      end if
    else
      send_field_dble = 0d0
      if (present(send_frac_mask)) then
        send_frac_mask_dble = 0d0
      end if
    end if
  end subroutine send_field_to_dble

  subroutine send_field_to_dble_single(field_id,        &
                                      nbr_hor_points,  &
                                      collection_size, &
                                      send_field,      &
                                      send_field_dble, &
                                      send_frac_mask,  &
                                      send_frac_mask_dble)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_put_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_put_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_put_mask_c2f_c

      end function yac_get_field_put_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: collection_size
    real, intent (in)     :: send_field(nbr_hor_points, &
                                        collection_size)
    double precision, intent (out) :: send_field_dble(nbr_hor_points, &
                                                      collection_size)
    real, optional, intent (in)     :: send_frac_mask(nbr_hor_points, &
                                                      collection_size)
    double precision, optional, intent (out) :: send_frac_mask_dble(nbr_hor_points, &
                                                                    collection_size)

    integer :: i, j
    type(c_ptr) :: put_mask_
    type(c_ptr), pointer :: put_mask(:)
    integer(kind=c_int), pointer :: pointset_put_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_SOURCE) then

      put_mask_ = yac_get_field_put_mask_c2f_c(field_id)
      if (c_associated(put_mask_)) then
        call c_f_pointer(put_mask_, put_mask, shape=[1])
        do i = 1, collection_size
          call c_f_pointer( &
                  put_mask(1), pointset_put_mask, shape=[nbr_hor_points])
          do j = 1, nbr_hor_points
            if (pointset_put_mask(j) /= 0) then
              send_field_dble(j, i) = dble(send_field(j, i))
            else
              send_field_dble(j, i) = 0d0
            end if
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            call c_f_pointer( &
                    put_mask(1), pointset_put_mask, shape=[nbr_hor_points])
            do j = 1, nbr_hor_points
              if (pointset_put_mask(j) /= 0) then
                send_frac_mask_dble(j, i) = dble(send_frac_mask(j, i))
              else
                send_frac_mask_dble(j, i) = 0d0
              end if
            end do
          end do
        end if
      else
        send_field_dble = dble(send_field)
        if (present(send_frac_mask)) then
          send_frac_mask_dble = dble(send_frac_mask)
        end if
      end if
    else
      send_field_dble = 0d0
      if (present(send_frac_mask)) then
        send_frac_mask_dble = 0d0
      end if
    end if
  end subroutine send_field_to_dble_single

  subroutine send_field_to_dble_ptr(field_id,        &
                                    nbr_pointsets,   &
                                    collection_size, &
                                    send_field,      &
                                    send_field_dble, &
                                    send_frac_mask,  &
                                    send_frac_mask_dble)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_put_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_put_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_put_mask_c2f_c

      end function yac_get_field_put_mask_c2f_c

    end interface

    integer, intent (in)            :: field_id
    integer, intent (in)            :: nbr_pointsets
    integer, intent (in)            :: collection_size
    type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, &
                                                  collection_size)
    type(yac_dble_ptr), intent (out) :: send_field_dble(nbr_pointsets, &
                                                        collection_size)
    type(yac_real_ptr), optional, intent (in) :: send_frac_mask(nbr_pointsets, &
                                                                collection_size)
    type(yac_dble_ptr), optional, intent (out) :: send_frac_mask_dble(nbr_pointsets, &
                                                                      collection_size)

    integer :: i, j, k, nbr_hor_points
    type(c_ptr) :: put_mask_
    type(c_ptr), pointer :: put_mask(:)
    integer(kind=c_int), pointer :: pointset_put_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_SOURCE) then

      put_mask_ = yac_get_field_put_mask_c2f_c(field_id)
      if (c_associated(put_mask_)) then
        call c_f_pointer(put_mask_, put_mask, shape=[nbr_pointsets])
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            nbr_hor_points = size(send_field(j,i)%p)
            allocate(send_field_dble(j,i)%p(nbr_hor_points))
            call c_f_pointer( &
                  put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
            do k = 1, nbr_hor_points
              if (pointset_put_mask(k) /= 0) then
                send_field_dble(j, i)%p(k) = dble(send_field(j, i)%p(k))
              else
                send_field_dble(j, i)%p(k) = 0d0
              end if
            end do
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            do j = 1, nbr_pointsets
              nbr_hor_points = size(send_frac_mask(j,i)%p)
              allocate(send_frac_mask_dble(j,i)%p(nbr_hor_points))
              call c_f_pointer( &
                    put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
              do k = 1, nbr_hor_points
                if (pointset_put_mask(k) /= 0) then
                  send_frac_mask_dble(j, i)%p(k) = dble(send_frac_mask(j, i)%p(k))
                else
                  send_frac_mask_dble(j, i)%p(k) = 0d0
                end if
              end do
            end do
          end do
        end if
      else
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            nbr_hor_points = size(send_field(j,i)%p)
            allocate(send_field_dble(j,i)%p(nbr_hor_points))
            send_field_dble(j,i)%p = dble(send_field(j,i)%p)
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            do j = 1, nbr_pointsets
              nbr_hor_points = size(send_frac_mask(j,i)%p)
              allocate(send_frac_mask_dble(j,i)%p(nbr_hor_points))
              send_frac_mask_dble(j,i)%p = dble(send_frac_mask(j,i)%p)
            end do
          end do
        end if
      end if
    else
      do i = 1, collection_size
        do j = 1, nbr_pointsets
          nbr_hor_points = size(send_field(j,i)%p)
          allocate(send_field_dble(j,i)%p(nbr_hor_points))
          send_field_dble(j,i)%p = 0d0
        end do
      end do
      if (present(send_frac_mask)) then
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            nbr_hor_points = size(send_frac_mask(j,i)%p)
            allocate(send_frac_mask_dble(j,i)%p(nbr_hor_points))
            send_frac_mask_dble(j,i)%p = 0d0
          end do
        end do
      end if
    end if
  end subroutine send_field_to_dble_ptr

  ! -----------------------------------------------------------------------

  subroutine recv_field_to_dble(field_id,        &
                                nbr_hor_points,  &
                                collection_size, &
                                recv_field,      &
                                recv_field_dble)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: collection_size
    real, intent (in)     :: recv_field(nbr_hor_points, &
                                        collection_size)
    double precision, intent (out) :: recv_field_dble(nbr_hor_points, &
                                                      collection_size)

    integer :: i, j
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_)) then
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field_dble(j, i) = dble(recv_field(j, i))
            else
              recv_field_dble(j, i) = 0d0
            end if
          end do
        end do
      else
        recv_field_dble = dble(recv_field)
      end if
    else
      recv_field_dble = 0d0
    end if
  end subroutine recv_field_to_dble

  subroutine recv_field_to_dble_ptr(field_id,        &
                                    collection_size, &
                                    recv_field,      &
                                    recv_field_dble)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)             :: field_id
    integer, intent (in)             :: collection_size
    type(yac_real_ptr), intent (in)  :: recv_field(collection_size)
    type(yac_dble_ptr), intent (out) :: recv_field_dble(collection_size)

    integer :: i, j, nbr_hor_points
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_) .and. (collection_size > 0)) then
        nbr_hor_points = size(recv_field(1)%p)
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          nbr_hor_points = size(recv_field(i)%p)
          allocate(recv_field_dble(i)%p(nbr_hor_points))
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field_dble(i)%p(j) = dble(recv_field(i)%p(j))
            else
              recv_field_dble(i)%p(j) = 0d0
            end if
          end do
        end do
      else
        do i = 1, collection_size
          nbr_hor_points = size(recv_field(i)%p)
          allocate(recv_field_dble(i)%p(nbr_hor_points))
          recv_field_dble(i)%p = dble(recv_field(i)%p)
        end do
      end if
    else
      do i = 1, collection_size
        nbr_hor_points = size(recv_field(i)%p)
        allocate(recv_field_dble(i)%p(nbr_hor_points))
        recv_field_dble(i)%p = 0d0
      end do
    end if
  end subroutine recv_field_to_dble_ptr

  ! -----------------------------------------------------------------------

  subroutine recv_field_from_dble(field_id,        &
                                  nbr_hor_points,  &
                                  collection_size, &
                                  recv_field_dble, &
                                  recv_field)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: collection_size
    double precision, intent (in) :: recv_field_dble(nbr_hor_points, &
                                                    collection_size)
    real, intent (inout)    :: recv_field(nbr_hor_points, &
                                          collection_size)

    integer :: i, j
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_)) then
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field(j, i) = real(recv_field_dble(j, i))
            end if
          end do
        end do
      else
        recv_field = real(recv_field_dble)
      end if
    end if
  end subroutine recv_field_from_dble

  subroutine recv_field_from_dble_ptr(field_id,        &
                                      collection_size, &
                                      recv_field_dble, &
                                      recv_field)

    use yac
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)               :: field_id
    integer, intent (in)               :: collection_size
    type(yac_dble_ptr), intent (inout) :: recv_field_dble(collection_size)
    type(yac_real_ptr), intent (inout) :: recv_field(collection_size)

    integer :: i, j, nbr_hor_points
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_) .and. (collection_size > 0)) then
        nbr_hor_points = size(recv_field(1)%p)
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          nbr_hor_points = size(recv_field(i)%p)
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field(i)%p(j) = real(recv_field_dble(i)%p(j))
            end if
          end do
          deallocate(recv_field_dble(i)%p)
        end do
      else
        do i = 1, collection_size
          recv_field(i)%p = real(recv_field_dble(i)%p)
          deallocate(recv_field_dble(i)%p)
        end do
      end if
    end if
  end subroutine recv_field_from_dble_ptr

end module mo_yac_real_to_dble_utils

! -------------------------- init -------------------------------------

subroutine yac_fmpi_handshake ( comm, group_names, group_comms )
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_char, c_loc
  use yac, dummy => yac_fmpi_handshake

  implicit none

  interface
     subroutine yac_cmpi_handshake_c (comm, n, group_names, group_comms) &
          bind ( c, name='yac_cmpi_handshake_f2c')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use yac, only : YAC_MPI_FINT_KIND
       integer (kind = YAC_MPI_FINT_KIND ), intent(in), value :: comm
       integer(c_int), intent(in), value :: n
       type (c_ptr) , intent(in) :: group_names(n)
       integer (kind = YAC_MPI_FINT_KIND ), intent(out) :: group_comms(n)
     end subroutine yac_cmpi_handshake_c
  end interface

  integer, intent(in) :: comm
  character(len=YAC_MAX_CHARLEN), intent(in) :: group_names(:)
  integer, intent(out) :: group_comms(SIZE(group_names))

  CHARACTER (kind=c_char, len=YAC_MAX_CHARLEN+1), TARGET :: &
    group_names_cpy(SIZE(group_names))
  type( c_ptr ) :: group_name_ptr(SIZE(group_names))
  integer :: i
  DO i=1,SIZE(group_names)
     group_names_cpy(i) = TRIM(group_names(i)) // c_null_char
     group_name_ptr(i) = c_loc(group_names_cpy(i))
  END DO

  call yac_cmpi_handshake_c( &
    comm, SIZE(group_names), group_name_ptr, group_comms)

end subroutine yac_fmpi_handshake

subroutine yac_finit_emitter_flags()

  use yac, only : YAC_YAML_EMITTER_DEFAULT_F, &
                  YAC_YAML_EMITTER_JSON_F

  implicit none

  interface

      function yac_cyaml_get_emitter_flag_default () &
        bind ( c, name='yac_cyaml_get_emitter_flag_default_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int) :: yac_cyaml_get_emitter_flag_default

      end function yac_cyaml_get_emitter_flag_default

      function yac_cyaml_get_emitter_flag_json () &
        bind ( c, name='yac_cyaml_get_emitter_flag_json_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int) :: yac_cyaml_get_emitter_flag_json

      end function yac_cyaml_get_emitter_flag_json

  end interface

  YAC_YAML_EMITTER_DEFAULT_F = yac_cyaml_get_emitter_flag_default()
  YAC_YAML_EMITTER_JSON_F = yac_cyaml_get_emitter_flag_json()

end subroutine

subroutine yac_finit_comm ( comm )

  use yac, dummy => yac_finit_comm

  implicit none

  interface

      subroutine yac_cinit_comm_c ( comm ) &
        bind ( c, name='yac_cinit_comm_f2c' )

        use yac, only : YAC_MPI_FINT_KIND

        integer ( kind=YAC_MPI_FINT_KIND ), value :: comm

      end subroutine yac_cinit_comm_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(in) :: comm

  call yac_finit_emitter_flags ( )
  call yac_cinit_comm_c ( comm )

end subroutine yac_finit_comm

subroutine yac_finit_comm_instance( comm, yac_instance_id)

  use yac, dummy => yac_finit_comm_instance

  implicit none

  interface

      subroutine yac_cinit_comm_instance_c ( comm, yac_instance_id) &
        bind ( c, name='yac_cinit_comm_instance_f2c' )

        use, intrinsic :: iso_c_binding, only : c_int
        use yac, only : YAC_MPI_FINT_KIND

        integer ( kind=YAC_MPI_FINT_KIND ), value :: comm
        integer (kind=c_int)                      :: yac_instance_id

      end subroutine yac_cinit_comm_instance_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(in)  :: comm
  integer, intent(out) :: yac_instance_id

  call yac_finit_emitter_flags ( )
  call yac_cinit_comm_instance_c ( comm, yac_instance_id )

end subroutine yac_finit_comm_instance

subroutine yac_finit (  )

  use yac, dummy => yac_finit

  implicit none

  interface

    subroutine yac_cinit_c (  ) &
      bind ( c, name='yac_cinit' )

    end subroutine yac_cinit_c

    subroutine yac_finit_emitter_flags()
    end subroutine yac_finit_emitter_flags

  end interface

  call yac_finit_emitter_flags ( )
  call yac_cinit_c ( )

end subroutine yac_finit

subroutine yac_finit_instance ( yac_instance_id )

  use yac, dummy => yac_finit_instance

  implicit none

  interface

      subroutine yac_cinit_instance_c ( yac_instance_id ) &
        bind ( c, name='yac_cinit_instance' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer (kind=c_int)                  :: yac_instance_id

      end subroutine yac_cinit_instance_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(out)         :: yac_instance_id !< [OUT] returned handle to the YAC instance

  call yac_finit_emitter_flags ( )
  call yac_cinit_instance_c ( yac_instance_id )

end subroutine yac_finit_instance

subroutine yac_finit_comm_dummy ( world_comm )

  use yac, dummy => yac_finit_comm_dummy

  implicit none

  interface

      subroutine yac_cinit_comm_dummy_c ( world_comm ) &
        bind ( c, name='yac_cinit_comm_dummy_f2c' )

        use yac, only : YAC_MPI_FINT_KIND

        integer ( kind=YAC_MPI_FINT_KIND ), value :: world_comm

      end subroutine yac_cinit_comm_dummy_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(in) :: world_comm !< [IN] MPI world communicator (optional)

  call yac_finit_emitter_flags ( )
  call yac_cinit_comm_dummy_c ( world_comm )

end subroutine yac_finit_comm_dummy

subroutine yac_finit_dummy ( )

  use yac, dummy => yac_finit_dummy

  implicit none

  interface

      subroutine yac_cinit_dummy_c ( ) &
        bind ( c, name='yac_cinit_dummy' )

      end subroutine yac_cinit_dummy_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  call yac_finit_emitter_flags ( )
  call yac_cinit_dummy_c ( )

end subroutine yac_finit_dummy

! ----------------------- getting default instance id ------------------------

function yac_fget_default_instance_id()

  use yac, dummy => yac_fget_default_instance_id

  implicit none

  interface

      function yac_cget_default_instance_id_c ( ) &
          bind ( c, name='yac_cget_default_instance_id' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer(kind=c_int) :: yac_cget_default_instance_id_c

      end function yac_cget_default_instance_id_c

  end interface

  integer :: yac_fget_default_instance_id

  yac_fget_default_instance_id = yac_cget_default_instance_id_c ( )

end function yac_fget_default_instance_id

! -------------------------- reading config file ---------------------------

subroutine yac_fread_config_yaml (yaml_filename)
  use yac, dummy => yac_fread_config_yaml
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_yaml_c(yaml_filename) &
          bind ( c, name='yac_cread_config_yaml' )
       use, intrinsic :: iso_c_binding, only : c_char
       character (kind=c_char), dimension(*) :: yaml_filename
     end subroutine yac_cread_config_yaml_c
  end interface

  character(len=*), intent(in) :: yaml_filename

  call yac_cread_config_yaml_c(TRIM(yaml_filename) // c_null_char)

end subroutine yac_fread_config_yaml

subroutine yac_fread_config_yaml_instance(yac_instance_id, yaml_filename)
  use yac, dummy => yac_fread_config_yaml_instance
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_yaml_instance_c( &
          yac_instance_id, yaml_filename) &
          bind ( c, name='yac_cread_config_yaml_instance' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       integer (kind=c_int), value           :: yac_instance_id
       character (kind=c_char), dimension(*) :: yaml_filename
     end subroutine yac_cread_config_yaml_instance_c
  end interface

  integer, intent(in)          :: yac_instance_id
  character(len=*), intent(in) :: yaml_filename

  call yac_cread_config_yaml_instance_c(yac_instance_id, &
       & TRIM(yaml_filename) // c_null_char)

end subroutine yac_fread_config_yaml_instance

subroutine yac_fread_config_json (json_filename)
  use yac, dummy => yac_fread_config_json
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_json_c(json_filename) &
          bind ( c, name='yac_cread_config_json' )
       use, intrinsic :: iso_c_binding, only : c_char
       character (kind=c_char), dimension(*) :: json_filename
     end subroutine yac_cread_config_json_c
  end interface

  character(len=*), intent(in) :: json_filename

  call yac_cread_config_json_c(TRIM(json_filename) // c_null_char)

end subroutine yac_fread_config_json

subroutine yac_fread_config_json_instance(yac_instance_id, json_filename)
  use yac, dummy => yac_fread_config_json_instance
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_json_instance_c( &
          yac_instance_id, json_filename) &
          bind ( c, name='yac_cread_config_json_instance' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       integer (kind=c_int), value           :: yac_instance_id
       character (kind=c_char), dimension(*) :: json_filename
     end subroutine yac_cread_config_json_instance_c
  end interface

  integer, intent(in)          :: yac_instance_id
  character(len=*), intent(in) :: json_filename

  call yac_cread_config_json_instance_c(yac_instance_id, &
       & TRIM(json_filename) // c_null_char)

end subroutine yac_fread_config_json_instance

! -------------------------- writing config file ---------------------------

subroutine yac_fset_config_output_file ( &
  filename, fileformat, sync_location, include_definitions)
  use yac, dummy => yac_fset_config_output_file
  use, intrinsic :: iso_c_binding, only : c_null_char, c_int
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cset_config_output_file_c(&
          filename, fileformat, sync_location, include_definitions) &
          bind ( c, name='yac_cset_config_output_file' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       character (kind=c_char), dimension(*) :: filename
       integer (kind=c_int), value           :: fileformat
       integer (kind=c_int), value           :: sync_location
       integer (kind=c_int), value           :: include_definitions
     end subroutine yac_cset_config_output_file_c
  end interface

  character(len=*), intent(in) :: filename
  integer, intent(in) :: fileformat
  integer, intent(in) :: sync_location
  logical, intent(in), optional :: include_definitions

  integer(kind=c_int) :: c_include_definitions

  if (present(include_definitions)) then
    c_include_definitions = MERGE(1_c_int, 0_c_int, include_definitions)
  else
    c_include_definitions = 0_c_int
  end if

  call yac_cset_config_output_file_c( &
    TRIM(filename) // c_null_char, fileformat, sync_location, &
    c_include_definitions)

end subroutine yac_fset_config_output_file

subroutine yac_fset_config_output_file_instance ( &
  yac_instance_id, filename, fileformat, sync_location, include_definitions)
  use yac, dummy => yac_fset_config_output_file_instance
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cset_config_output_file_instance_c(&
          yac_instance_id, filename, fileformat, sync_location, &
          include_definitions) &
          bind ( c, name='yac_cset_config_output_file_instance' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       integer (kind=c_int), value           :: yac_instance_id
       character (kind=c_char), dimension(*) :: filename
       integer (kind=c_int), value           :: fileformat
       integer (kind=c_int), value           :: sync_location
       integer (kind=c_int), value           :: include_definitions
     end subroutine yac_cset_config_output_file_instance_c
  end interface

  integer, intent(in) :: yac_instance_id
  character(len=*), intent(in) :: filename
  integer, intent(in) :: fileformat
  integer, intent(in) :: sync_location
  logical, intent(in), optional :: include_definitions

  integer(kind=c_int) :: c_include_definitions

  if (present(include_definitions)) then
    c_include_definitions = MERGE(1_c_int, 0_c_int, include_definitions)
  else
    c_include_definitions = 0_c_int
  end if

  call yac_cset_config_output_file_instance_c( &
    yac_instance_id, TRIM(filename) // c_null_char, fileformat, sync_location, &
    c_include_definitions)

end subroutine yac_fset_config_output_file_instance

! -------------------------- writing grid file ---------------------------

subroutine yac_fset_grid_output_file ( &
  gridname, filename)
  use yac, dummy => yac_fset_grid_output_file
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cset_grid_output_file_c(&
          gridname, filename) &
          bind ( c, name='yac_cset_grid_output_file' )
       use, intrinsic :: iso_c_binding, only : c_char
       character (kind=c_char), dimension(*) :: gridname
       character (kind=c_char), dimension(*) :: filename
     end subroutine yac_cset_grid_output_file_c
  end interface

  character(len=*), intent(in) :: gridname
  character(len=*), intent(in) :: filename

  call yac_cset_grid_output_file_c( &
    TRIM(gridname) // c_null_char, TRIM(filename) // c_null_char)

end subroutine yac_fset_grid_output_file

subroutine yac_fset_grid_output_file_instance ( &
  yac_instance_id, gridname, filename)
  use yac, dummy => yac_fset_grid_output_file_instance
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cset_grid_output_file_instance_c(&
          yac_instance_id, gridname, filename) &
          bind ( c, name='yac_cset_grid_output_file_instance' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       integer (kind=c_int), value           :: yac_instance_id
       character (kind=c_char), dimension(*) :: gridname
       character (kind=c_char), dimension(*) :: filename
     end subroutine yac_cset_grid_output_file_instance_c
  end interface

  integer, intent(in) :: yac_instance_id
  character(len=*), intent(in) :: gridname
  character(len=*), intent(in) :: filename

  call yac_cset_grid_output_file_instance_c( &
    yac_instance_id, TRIM(gridname) // c_null_char, &
    TRIM(filename) // c_null_char)

end subroutine yac_fset_grid_output_file_instance

! -------------------------- cleanup -----------------------------------

subroutine yac_fcleanup ( )

  use yac, dummy => yac_fcleanup

  implicit none

  interface

      subroutine yac_ccleanup_c () bind ( c, name='yac_ccleanup' )
      end subroutine yac_ccleanup_c

  end interface

  call yac_ccleanup_c ( )

end subroutine yac_fcleanup

subroutine yac_fcleanup_instance ( yac_instance_id )

  use yac, dummy => yac_fcleanup_instance

  implicit none

  interface

      subroutine yac_ccleanup_instance_c ( yac_instance_id ) &
        bind ( c, name='yac_ccleanup_instance' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer (kind=c_int), value :: yac_instance_id

      end subroutine yac_ccleanup_instance_c

  end interface

  integer, intent(in) :: yac_instance_id !< [IN] YAC instance identifier

  call yac_ccleanup_instance_c ( yac_instance_id )

end subroutine yac_fcleanup_instance

! -------------------------- final -------------------------------------

subroutine yac_ffinalize ( )

  use yac, dummy => yac_ffinalize

  implicit none

  interface
     subroutine yac_cfinalize_c () bind ( c, name='yac_cfinalize' )
     end subroutine yac_cfinalize_c
  end interface

  call yac_cfinalize_c ( )

end subroutine yac_ffinalize

subroutine yac_ffinalize_instance ( yac_instance_id )

  use yac, dummy => yac_ffinalize_instance

  implicit none

  interface
     subroutine yac_cfinalize_instance_c ( yac_instance_id ) &
        bind ( c, name='yac_cfinalize_instance' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer (kind=c_int), value :: yac_instance_id

     end subroutine yac_cfinalize_instance_c
  end interface

  integer, intent(in) :: yac_instance_id !< [IN] YAC instance identifier

  call yac_cfinalize_instance_c ( yac_instance_id )

end subroutine yac_ffinalize_instance

! -------------------------- version ----------------------------------

function yac_fget_version () result (version_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use yac, dummy => yac_fget_version
  use mo_yac_iso_c_helpers

  implicit none

  interface
     function yac_cget_version_c () bind ( c, name='yac_cget_version' )

     use, intrinsic :: iso_c_binding, only : c_ptr
     type(c_ptr) :: yac_cget_version_c

     end function yac_cget_version_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: version_string

  c_string_ptr = yac_cget_version_c()
  version_string = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_version

! -------------------------- mpi_handshake group name -----------------

function yac_fget_mpi_handshake_group_name () result (mpi_handshake_group_name)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use yac, dummy => yac_fget_mpi_handshake_group_name
  use mo_yac_iso_c_helpers

  implicit none

  interface
    function yac_cget_mpi_handshake_group_name_c () &
      bind ( c, name='yac_cget_mpi_handshake_group_name' )

    use, intrinsic :: iso_c_binding, only : c_ptr
    type(c_ptr) :: yac_cget_mpi_handshake_group_name_c

    end function yac_cget_mpi_handshake_group_name_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: mpi_handshake_group_name

  c_string_ptr = yac_cget_mpi_handshake_group_name_c()
  mpi_handshake_group_name = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_mpi_handshake_group_name

! -------------------------- dates ------------------------------------

subroutine yac_fdef_datetime ( start_datetime, end_datetime )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_datetime
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_datetime_c ( start_datetime, end_datetime ) &
    &      bind ( c, name='yac_cdef_datetime' )

       use, intrinsic :: iso_c_binding, only : c_char

       character ( kind=c_char), dimension(*) :: start_datetime
       character ( kind=c_char), dimension(*) :: end_datetime

     end subroutine yac_cdef_datetime_c

  end interface

  character(len=*), intent(in), optional :: start_datetime !< [IN] start datetime of job
  character(len=*), intent(in), optional :: end_datetime   !< [IN] end datetime of job

  integer :: index

  index = 0


  if (present(start_datetime)) then
    YAC_CHECK_STRING_LEN ( "yac_fdef_datetime", start_datetime )
    index = index + 1
  end if

  if (present(end_datetime)) then
    YAC_CHECK_STRING_LEN ( "yac_fdef_datetime", end_datetime )
    index = index + 2
  end if

  select case ( index )

    case ( 3 )
      call yac_cdef_datetime_c ( TRIM(start_datetime) // c_null_char, &
                                 TRIM(end_datetime)   // c_null_char )
    case ( 2 )
      call yac_cdef_datetime_c ( c_null_char, &
                                 TRIM(end_datetime)   // c_null_char )
    case ( 1 )
      call yac_cdef_datetime_c ( TRIM(start_datetime) // c_null_char, &
                                 c_null_char )
    end select

end subroutine yac_fdef_datetime

subroutine yac_fdef_datetime_instance ( &
  yac_instance_id, start_datetime, end_datetime )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_datetime_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_datetime_instance_c ( yac_instance_id, &
                                               start_datetime,  &
                                               end_datetime )   &
       bind ( c, name='yac_cdef_datetime_instance' )

       use, intrinsic :: iso_c_binding, only : c_char, c_int

       integer (kind=c_int), value            :: yac_instance_id
       character ( kind=c_char), dimension(*) :: start_datetime
       character ( kind=c_char), dimension(*) :: end_datetime

     end subroutine yac_cdef_datetime_instance_c

  end interface

  integer, intent(in)                    :: yac_instance_id !< [IN] YAC instance identifier
  character(len=*), intent(in), optional :: start_datetime  !< [IN] start datetime of job
  character(len=*), intent(in), optional :: end_datetime    !< [IN] end datetime of job

  integer :: index

  index = 0


  if (present(start_datetime)) then
    YAC_CHECK_STRING_LEN ( "yac_fdef_datetime_instance", start_datetime )
    index = index + 1
  end if

  if (present(end_datetime)) then
    YAC_CHECK_STRING_LEN ( "yac_fdef_datetime_instance", end_datetime )
    index = index + 2
  end if

  select case ( index )

    case ( 3 )
      call yac_cdef_datetime_instance_c ( yac_instance_id, &
                                          TRIM(start_datetime) // c_null_char, &
                                          TRIM(end_datetime)   // c_null_char )
    case ( 2 )
      call yac_cdef_datetime_instance_c ( yac_instance_id, &
                                          c_null_char, &
                                          TRIM(end_datetime)   // c_null_char )
    case ( 1 )
      call yac_cdef_datetime_instance_c ( yac_instance_id, &
                                          TRIM(start_datetime) // c_null_char, &
                                           c_null_char )
    end select

end subroutine yac_fdef_datetime_instance

subroutine yac_fdef_calendar ( calendar )

  use yac, dummy => yac_fdef_calendar

  implicit none

  interface

      subroutine yac_cdef_calendar_c ( calendar ) &
        bind ( c, name='yac_cdef_calendar' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: calendar

      end subroutine yac_cdef_calendar_c

   end interface

  integer, intent(in) :: calendar

  call yac_cdef_calendar_c ( calendar )

end subroutine yac_fdef_calendar

subroutine yac_fget_calendar ( calendar )

  use yac, dummy => yac_fget_calendar

  implicit none

  interface

      function yac_cget_calendar_c ( ) &
        bind ( c, name='yac_cget_calendar' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ) :: yac_cget_calendar_c

      end function yac_cget_calendar_c

   end interface

  integer, intent(out) :: calendar

  calendar = yac_cget_calendar_c()

end subroutine yac_fget_calendar

! ------------------------ predef_comp ------------------------------------

SUBROUTINE yac_fpredef_comp ( comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fpredef_comp
  use mo_yac_iso_c_helpers

  implicit none

  INTERFACE

     SUBROUTINE yac_cpredef_comp_c ( comp_name, comp_id ) &
          bind ( c, name='yac_cpredef_comp' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     END SUBROUTINE yac_cpredef_comp_c

  END INTERFACE

  character(len=*), intent(in) :: comp_name !< [IN]  component name
  integer, intent(out)         :: comp_id   !< [OUT] returned handle to the component

  YAC_CHECK_STRING_LEN ( "yac_fpredef_comp", comp_name )

  call yac_cpredef_comp_c ( TRIM(comp_name) // c_null_char, comp_id )

END SUBROUTINE yac_fpredef_comp

SUBROUTINE yac_fpredef_comp_instance ( yac_instance_id, comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fpredef_comp_instance
  use mo_yac_iso_c_helpers

  implicit none

  INTERFACE

     SUBROUTINE yac_cpredef_comp_instance_c ( yac_instance_id, &
          comp_name,       &
          comp_id )        &
          bind ( c, name='yac_cpredef_comp_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer (kind=c_int), value            :: yac_instance_id
       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     END SUBROUTINE yac_cpredef_comp_instance_c

  END INTERFACE

  integer, intent(in)          :: yac_instance_id !< [IN]  YAC instance identifier
  character(len=*), intent(in) :: comp_name       !< [IN]  component name
  integer, intent(out)         :: comp_id         !< [OUT] returned handle to the component

  YAC_CHECK_STRING_LEN ( "yac_fpredef_comp_instance", comp_name )

  call yac_cpredef_comp_instance_c( yac_instance_id,                &
       TRIM(comp_name) // c_null_char, &
       comp_id )

END SUBROUTINE yac_fpredef_comp_instance

! ------------------------ def_comp ------------------------------------

subroutine yac_fdef_comp ( comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_comp
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comp_c ( comp_name, comp_id ) &
       bind ( c, name='yac_cdef_comp' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     end subroutine yac_cdef_comp_c

  end interface

  character(len=*), intent(in) :: comp_name !< [IN]  component name
  integer, intent(out)         :: comp_id   !< [OUT] returned handle to the component

  YAC_CHECK_STRING_LEN ( "yac_fdef_comp", comp_name )

  call yac_cdef_comp_c ( TRIM(comp_name) // c_null_char, comp_id )

end subroutine yac_fdef_comp

subroutine yac_fdef_comp_instance ( yac_instance_id, comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_comp_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comp_instance_c ( yac_instance_id, &
                                           comp_name,       &
                                           comp_id )        &
       bind ( c, name='yac_cdef_comp_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer (kind=c_int), value            :: yac_instance_id
       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     end subroutine yac_cdef_comp_instance_c

  end interface

  integer, intent(in)          :: yac_instance_id !< [IN]  YAC instance identifier
  character(len=*), intent(in) :: comp_name       !< [IN]  component name
  integer, intent(out)         :: comp_id         !< [OUT] returned handle to the component

  YAC_CHECK_STRING_LEN ( "yac_fdef_comp_instance", comp_name )

  call yac_cdef_comp_instance_c( yac_instance_id,                &
                                 TRIM(comp_name) // c_null_char, &
                                 comp_id )

end subroutine yac_fdef_comp_instance

! ------------------------ def_comps ------------------------------------

subroutine yac_fdef_comps ( comp_names, num_comps, comp_ids )

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
  use yac, dummy => yac_fdef_comps
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comps_c ( comp_names, num_comps, comp_ids ) &
           bind ( c, name='yac_cdef_comps' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       type ( c_ptr ), value         :: comp_names
       integer ( kind=c_int ), value :: num_comps
       integer ( kind=c_int )        :: comp_ids(*)

     end subroutine yac_cdef_comps_c

  end interface

  integer, intent(in)   :: num_comps !< [IN]  number of components
  character(kind=c_char, len=*), intent(in) :: &
                           comp_names(num_comps) !< [IN]  component names
  integer, intent(out)  :: comp_ids(num_comps) !< [OUT] returned handle to the components

  integer :: i, j
  character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
  type(c_ptr), target :: comp_name_ptrs(num_comps)

  comp_names_cpy = c_null_char

  do i = 1, num_comps
     YAC_CHECK_STRING_LEN ( "yac_fdef_comps", comp_names(i))
    do j = 1, len_trim(comp_names(i))
      comp_names_cpy(j,i) = comp_names(i)(j:j)
    end do
    comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
  end do

  call yac_cdef_comps_c ( c_loc(comp_name_ptrs), num_comps, comp_ids )

end subroutine yac_fdef_comps

subroutine yac_fdef_comps_instance ( yac_instance_id, &
                                    comp_names,       &
                                    num_comps,        &
                                    comp_ids )

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
  use yac, dummy => yac_fdef_comps_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comps_instance_c ( yac_instance_id, &
                                            comp_names,      &
                                            num_comps,       &
                                            comp_ids )       &
       bind ( c, name='yac_cdef_comps_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer (kind=c_int), value   :: yac_instance_id
       type ( c_ptr ), value         :: comp_names
       integer ( kind=c_int ), value :: num_comps
       integer ( kind=c_int )        :: comp_ids(*)

     end subroutine yac_cdef_comps_instance_c

  end interface

  integer, intent(in)   :: yac_instance_id       !< [IN]  YAC instance identifier
  integer, intent(in)   :: num_comps             !< [IN]  number of components
  character(kind=c_char, len=*), intent(in) :: &
                           comp_names(num_comps) !< [IN]  component names
  integer, intent(out)  :: comp_ids(num_comps)   !< [OUT] returned handle to the components

  integer :: i, j
  character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
  type(c_ptr), target :: comp_name_ptrs(num_comps)

  comp_names_cpy = c_null_char

  do i = 1, num_comps
     YAC_CHECK_STRING_LEN ( "yac_fdef_comps_instance", comp_names(i))
    do j = 1, len_trim(comp_names(i))
      comp_names_cpy(j,i) = comp_names(i)(j:j)
    end do
    comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
  end do

  call yac_cdef_comps_instance_c ( yac_instance_id, &
                                   c_loc(comp_name_ptrs),  &
                                   num_comps,       &
                                   comp_ids )

end subroutine yac_fdef_comps_instance


! ------------------------- def_comp_dummy ------------------------------

subroutine yac_fdef_comp_dummy ( )

  use, intrinsic :: iso_c_binding, only : c_null_ptr, c_int
  use yac, dummy => yac_fdef_comp_dummy

  implicit none

  interface

     subroutine yac_cdef_comps_dummy_c ( comp_names, num_comps, comp_ids ) &
           bind ( c, name='yac_cdef_comps' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       type ( c_ptr ), value         :: comp_names
       integer ( kind=c_int ), value :: num_comps
       integer ( kind=c_int )        :: comp_ids(*)

     end subroutine yac_cdef_comps_dummy_c

  end interface

  call yac_cdef_comps_dummy_c ( c_null_ptr, 0_c_int, [ integer( kind=c_int ) :: ] )

end subroutine yac_fdef_comp_dummy

subroutine yac_fdef_comp_dummy_instance ( yac_instance_id )

  use, intrinsic :: iso_c_binding, only : c_null_ptr, c_int
  use yac, dummy => yac_fdef_comp_dummy_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cdef_comps_dummy_instance_c ( &
           yac_instance_id, comp_names, num_comps, comp_ids ) &
           bind ( c, name='yac_cdef_comps_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer (kind=c_int), value   :: yac_instance_id
       type ( c_ptr ), value         :: comp_names
       integer ( kind=c_int ), value :: num_comps
       integer ( kind=c_int )        :: comp_ids(*)

     end subroutine yac_cdef_comps_dummy_instance_c

  end interface

  integer, intent(in)          :: yac_instance_id !< [IN]  YAC instance identifier

  call yac_cdef_comps_dummy_instance_c( yac_instance_id, &
                                        c_null_ptr,      &
                                        0_c_int,         &
                                        [ integer( kind=c_int ) :: ] )

end subroutine yac_fdef_comp_dummy_instance

! ------------------------- def_points ----------------------------------

subroutine yac_fdef_points_reg2d_real ( grid_id,       &
                                        nbr_points,    &
                                        location,      &
                                        x_points_real, &
                                        y_points_real, &
                                        point_id )

  use yac, dummy => yac_fdef_points_reg2d_real

  implicit none

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                !< [IN]  number of points in x and y
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
  real, intent(in)     :: x_points_real(nbr_points(1)) !< [IN]  longitudes of points
  real, intent(in)     :: y_points_real(nbr_points(2)) !< [IN]  latitudes of points
  integer, intent(out) :: point_id                     !< [OUT] point identifier

  double precision    ::  x_points(nbr_points(1))
  double precision    ::  y_points(nbr_points(2))

  x_points(:) = dble(x_points_real(:))
  y_points(:) = dble(y_points_real(:))

  call yac_fdef_points_reg2d_dble ( grid_id,    &
                                    nbr_points, &
                                    location,   &
                                    x_points,   &
                                    y_points,   &
                                    point_id )

end subroutine yac_fdef_points_reg2d_real

subroutine yac_fdef_points_reg2d_dble ( grid_id,    &
                                        nbr_points, &
                                        location,   &
                                        x_points,   &
                                        y_points,   &
                                        point_id )

  use yac, dummy => yac_fdef_points_reg2d_dble

  implicit none

  interface

     subroutine yac_cdef_points_reg2d_c ( grid_id,      &
                                          nbr_points,   &
                                          location,     &
                                          x_points,     &
                                          y_points,     &
                                          point_id )    &
    &      bind ( c, name='yac_cdef_points_reg2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_int )        :: nbr_points(2)
       integer ( kind=c_int ), value :: location

       real    ( kind=c_double )     :: x_points(nbr_points(1))
       real    ( kind=c_double )     :: y_points(nbr_points(2))

       integer ( kind=c_int )        :: point_id

     end subroutine yac_cdef_points_reg2d_c

  end interface

  integer, intent(in)  :: grid_id                         !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                   !< [IN]  number of points in x and y
  integer, intent(in)  :: location                        !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: x_points(nbr_points(1)) !< [IN]  longitudes of points
  double precision, intent(in) :: y_points(nbr_points(2)) !< [IN]  latitudes of points

  integer, intent(out) :: point_id                        !< [OUT] point identifier

  call yac_cdef_points_reg2d_c ( grid_id,    &
                                 nbr_points, &
                                 location,   &
                                 x_points,   &
                                 y_points,   &
                                 point_id )

end subroutine yac_fdef_points_reg2d_dble

subroutine yac_fdef_points_curve2d_real ( grid_id,       &
                                          nbr_points,    &
                                          location,      &
                                          x_points_real, &
                                          y_points_real, &
                                          point_id )

  use yac, dummy => yac_fdef_points_curve2d_real

  implicit none

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                !< [IN]  number of points in x and y
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
  real, intent(in)     :: &
    x_points_real(nbr_points(1),nbr_points(2))         !< [IN]  longitudes of points
  real, intent(in)     :: &
    y_points_real(nbr_points(1),nbr_points(2))         !< [IN]  latitudes of points
  integer, intent(out) :: point_id                     !< [OUT] point identifier

  double precision    ::  x_points(nbr_points(1),nbr_points(2))
  double precision    ::  y_points(nbr_points(1),nbr_points(2))

  x_points(:,:) = dble(x_points_real(:,:))
  y_points(:,:) = dble(y_points_real(:,:))

  call yac_fdef_points_curve2d_dble ( grid_id,    &
                                      nbr_points, &
                                      location,   &
                                      x_points,   &
                                      y_points,   &
                                      point_id )

end subroutine yac_fdef_points_curve2d_real

subroutine yac_fdef_points_curve2d_dble ( grid_id,    &
                                          nbr_points, &
                                          location,   &
                                          x_points,   &
                                          y_points,   &
                                          point_id )

  use yac, dummy => yac_fdef_points_curve2d_dble

  implicit none

  interface

     subroutine yac_cdef_points_curve2d_c ( grid_id,    &
                                            nbr_points, &
                                            location,   &
                                            x_points,   &
                                            y_points,   &
                                            point_id )  &
           bind ( c, name='yac_cdef_points_curve2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_int )        :: nbr_points(2)
       integer ( kind=c_int ), value :: location

       real    ( kind=c_double )     :: x_points(nbr_points(1),nbr_points(2))
       real    ( kind=c_double )     :: y_points(nbr_points(1),nbr_points(2))

       integer ( kind=c_int )        :: point_id

     end subroutine yac_cdef_points_curve2d_c

  end interface

  integer, intent(in)  :: grid_id                         !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                   !< [IN]  number of points in x and y
  integer, intent(in)  :: location                        !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: &
    x_points(nbr_points(1),nbr_points(2))                 !< [IN]  longitudes of points
  double precision, intent(in) :: &
    y_points(nbr_points(1),nbr_points(2))                 !< [IN]  latitudes of points

  integer, intent(out) :: point_id                        !< [OUT] point identifier

  call yac_cdef_points_curve2d_c ( grid_id,    &
                                   nbr_points, &
                                   location,   &
                                   x_points,   &
                                   y_points,   &
                                   point_id )

end subroutine yac_fdef_points_curve2d_dble

subroutine yac_fdef_points_unstruct_real ( grid_id,       &
                                           nbr_points,    &
                                           location,      &
                                           x_points_real, &
                                           y_points_real, &
                                           point_id )

  use yac, dummy => yac_fdef_points_unstruct_real

  implicit none

  integer, intent(in)  :: grid_id                   !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points                !< [IN]  number of points
  integer, intent(in)  :: location                  !< [IN]  location, one of center/edge/vertex

  real, intent(in)     :: x_points_real(nbr_points) !< [IN]  longitudes of points
  real, intent(in)     :: y_points_real(nbr_points) !< [IN]  latitudes of points

  integer, intent(out) :: point_id                  !< [OUT] point identifier

  double precision    ::  x_points(nbr_points)
  double precision    ::  y_points(nbr_points)

  x_points(:) = dble(x_points_real(:))
  y_points(:) = dble(y_points_real(:))

  call yac_fdef_points_unstruct_dble ( grid_id,    &
                                       nbr_points, &
                                       location,   &
                                       x_points,   &
                                       y_points,   &
                                       point_id )

end subroutine yac_fdef_points_unstruct_real

subroutine yac_fdef_points_unstruct_dble ( grid_id,    &
                                           nbr_points, &
                                           location,   &
                                           x_points,   &
                                           y_points,   &
                                           point_id )

  use yac, dummy => yac_fdef_points_unstruct_dble

  implicit none

  interface

     subroutine yac_cdef_points_unstruct_c ( grid_id,      &
                                             nbr_points,   &
                                             location,     &
                                             x_points,     &
                                             y_points,     &
                                             point_id )    &
           bind ( c, name='yac_cdef_points_unstruct' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer (kind=c_int), value :: grid_id
       integer (kind=c_int), value :: nbr_points
       integer (kind=c_int), value :: location

       real    (kind=c_double)     :: x_points(nbr_points)
       real    (kind=c_double)     :: y_points(nbr_points)

       integer (kind=c_int)        :: point_id

     end subroutine yac_cdef_points_unstruct_c

  end interface

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points                   !< [IN]  number of points
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: x_points(nbr_points) !< [IN]  longitudes of points
  double precision, intent(in) :: y_points(nbr_points) !< [IN]  latitudes of points

  integer, intent(out) :: point_id                     !< [OUT] grid identifier

  call yac_cdef_points_unstruct_c  ( grid_id,    &
                                     nbr_points, &
                                     location,   &
                                     x_points,   &
                                     y_points,   &
                                     point_id )

end subroutine yac_fdef_points_unstruct_dble

subroutine yac_fdef_points_reg2d_rot_real ( grid_id,           &
                                            nbr_points,        &
                                            location,          &
                                            x_points_real,     &
                                            y_points_real,     &
                                            x_north_pole_real, &
                                            y_north_pole_real, &
                                            point_id )

  use yac, dummy => yac_fdef_points_reg2d_rot_real

  implicit none

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                !< [IN]  number of points in x and y
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
  real, intent(in)     :: x_points_real(nbr_points(1)) !< [IN]  longitudes of points in radians
  real, intent(in)     :: y_points_real(nbr_points(2)) !< [IN]  latitudes of points in radians
  real, intent(in)     :: x_north_pole_real            !< [IN]  longitude of north pole in radians
  real, intent(in)     :: y_north_pole_real            !< [IN]  latitude of north pole in radians
  integer, intent(out) :: point_id                     !< [OUT] point identifier

  double precision    ::  x_points(nbr_points(1))
  double precision    ::  y_points(nbr_points(2))
  double precision    ::  x_north_pole
  double precision    ::  y_north_pole

  x_points(:) = dble(x_points_real(:))
  y_points(:) = dble(y_points_real(:))

  x_north_pole = dble(x_north_pole_real)
  y_north_pole = dble(y_north_pole_real)

  call yac_fdef_points_reg2d_rot_dble ( grid_id,      &
                                        nbr_points,   &
                                        location,     &
                                        x_points,     &
                                        y_points,     &
                                        x_north_pole, &
                                        y_north_pole, &
                                        point_id )

end subroutine yac_fdef_points_reg2d_rot_real

subroutine yac_fdef_points_reg2d_rot_dble ( grid_id,      &
                                            nbr_points,   &
                                            location,     &
                                            x_points,     &
                                            y_points,     &
                                            x_north_pole, &
                                            y_north_pole, &
                                            point_id )

  use yac, dummy => yac_fdef_points_reg2d_rot_dble

  implicit none

  interface

     subroutine yac_cdef_points_reg2d_rot_c ( grid_id,      &
                                              nbr_points,   &
                                              location,     &
                                              x_points,     &
                                              y_points,     &
                                              x_north_pole, &
                                              y_north_pole, &
                                              point_id )    &
    &      bind ( c, name='yac_cdef_points_reg2d_rot' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value    :: grid_id
       integer ( kind=c_int )           :: nbr_points(2)
       integer ( kind=c_int ), value    :: location

       real    ( kind=c_double )        :: x_points(nbr_points(1))
       real    ( kind=c_double )        :: y_points(nbr_points(2))

       real    ( kind=c_double ), value :: x_north_pole
       real    ( kind=c_double ), value :: y_north_pole

       integer ( kind=c_int )           :: point_id

     end subroutine yac_cdef_points_reg2d_rot_c

  end interface

  integer, intent(in)  :: grid_id                         !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                   !< [IN]  number of points in x and y
  integer, intent(in)  :: location                        !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: x_points(nbr_points(1)) !< [IN]  longitudes of points in radians
  double precision, intent(in) :: y_points(nbr_points(2)) !< [IN]  latitudes of points in radians
  double precision, intent(in) :: x_north_pole            !< [IN]  longitude of north pole in radians
  double precision, intent(in) :: y_north_pole            !< [IN]  latitude of north pole in radians

  integer, intent(out) :: point_id                        !< [OUT] point identifier

  call yac_cdef_points_reg2d_rot_c ( grid_id,      &
                                     nbr_points,   &
                                     location,     &
                                     x_points,     &
                                     y_points,     &
                                     x_north_pole, &
                                     y_north_pole, &
                                     point_id )

end subroutine yac_fdef_points_reg2d_rot_dble

! ------------------------- def_grid -------------------------------

!> Definition of a non-uniform unstructured grid (cells have
!! varying numbers of vertices)
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_connections       total size of cell_to_vertex
!! @param[in]  nbr_vertices_per_cell number of vertices for each cell
!! @param[in]  x_vertices_real       longitudes of vertices
!! @param[in]  y_vertices_real       latitudes of vertices
!! @param[in]  cell_to_vertex_in     connectivity between vertices and cells\n
!!                                   (the vertex indices per cell have to be
!!                                   in clockwise or counterclockwise ordering)
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_nonuniform_real ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_connections,       &
                                           nbr_vertices_per_cell, &
                                           x_vertices_real,       &
                                           y_vertices_real,       &
                                           cell_to_vertex_in,     &
                                           grid_id,               &
                                           use_ll_edges)

  use yac, dummy => yac_fdef_grid_nonuniform_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices
  integer, intent(in)  :: nbr_cells
  integer, intent(in)  :: nbr_connections
  integer, intent(in)  :: nbr_vertices_per_cell(nbr_cells)

  real, intent(in)     :: x_vertices_real(nbr_vertices)
  real, intent(in)     :: y_vertices_real(nbr_vertices)

  integer, intent(in)  :: cell_to_vertex_in(nbr_connections)

  integer, intent(out) :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  double precision     :: x_vertices(nbr_vertices)
  double precision     :: y_vertices(nbr_vertices)

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_nonuniform_dble ( grid_name,             &
                                       nbr_vertices,          &
                                       nbr_cells,             &
                                       nbr_connections,       &
                                       nbr_vertices_per_cell, &
                                       x_vertices,            &
                                       y_vertices,            &
                                       cell_to_vertex_in,     &
                                       grid_id,               &
                                       use_ll_edges )

end subroutine yac_fdef_grid_nonuniform_real

!> Definition of a non-uniform unstructured grid (cells have
!! varying numbers of vertices)
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_connections       total size of cell_to_vertex
!! @param[in]  nbr_vertices_per_cell number of vertices for each cell
!! @param[in]  x_vertices            longitudes of vertices
!! @param[in]  y_vertices            latitudes of vertices
!! @param[in]  cell_to_vertex_in     connectivity between vertices and cells\n
!!                                   (the vertex indices per cell have to be in
!!                                   clockwise or counterclockwise ordering)
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_nonuniform_dble ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_connections,       &
                                           nbr_vertices_per_cell, &
                                           x_vertices,            &
                                           y_vertices,            &
                                           cell_to_vertex_in,     &
                                           grid_id,               &
                                           use_ll_edges )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_nonuniform_dble

  implicit none

  interface

     subroutine yac_cdef_grid_unstruct_c ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_vertices_per_cell, &
                                           x_vertices,            &
                                           y_vertices,            &
                                           cell_to_vertex,        &
                                           grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(*)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_c

     subroutine yac_cdef_grid_unstruct_ll_c ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_vertices_per_cell, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_vertex,        &
                                              grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct_ll' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(*)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_ll_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices
  integer, intent(in)          :: nbr_cells
  integer, intent(in)          :: nbr_connections
  integer, intent(in)          :: nbr_vertices_per_cell(nbr_cells)

  double precision, intent(in) :: x_vertices(nbr_vertices)
  double precision, intent(in) :: y_vertices(nbr_vertices)

  integer, intent(in)          :: cell_to_vertex_in(nbr_connections)

  integer, intent(out)         :: grid_id

  logical, optional, intent(in) :: use_ll_edges


  integer                      :: cell_to_vertex(nbr_connections)

  logical :: use_ll_edges_

  YAC_FASSERT(ALL(cell_to_vertex_in > 0), "ERROR(yac_fdef_grid_nonuniform_dble): all entries of cell_to_vertex have to be > 0")

  YAC_FASSERT(ALL(cell_to_vertex_in <= nbr_vertices), "ERROR(yac_fdef_grid_nonuniform_dble): all entries of cell_to_vertex have to be <= nbr_vertices")

  cell_to_vertex(:) = cell_to_vertex_in(:) - 1

  if (present(use_ll_edges)) then
    use_ll_edges_ = use_ll_edges
  else
    use_ll_edges_ = .false.
  end if

  if (use_ll_edges_) then
    call yac_cdef_grid_unstruct_ll_c  ( TRIM(grid_name) // c_null_char, &
                                        nbr_vertices,                   &
                                        nbr_cells,                      &
                                        nbr_vertices_per_cell,          &
                                        x_vertices,                     &
                                        y_vertices,                     &
                                        cell_to_vertex,                 &
                                        grid_id )
  else
    call yac_cdef_grid_unstruct_c  ( TRIM(grid_name) // c_null_char, &
                                     nbr_vertices,                   &
                                     nbr_cells,                      &
                                     nbr_vertices_per_cell,          &
                                     x_vertices,                     &
                                     y_vertices,                     &
                                     cell_to_vertex,                 &
                                     grid_id )
  end if

end subroutine yac_fdef_grid_nonuniform_dble

!> Definition of a uniform unstructured grid (all cells have the
!! number of vertices)
!! @param[in]  grid_name                grid name
!! @param[in]  nbr_vertices             number of vertices
!! @param[in]  nbr_cells                number of cells
!! @param[in]  nbr_vertices_per_cell_in number of vertices for each cell
!! @param[in]  x_vertices_real          longitudes of vertices
!! @param[in]  y_vertices_real          latitudes of vertices
!! @param[in]  cell_to_vertex_in        connectivity between vertices and cells\n
!!                                      (the vertex indices per cell have to be
!!                                      in clockwise or counterclockwise ordering)
!! @param[out] grid_id                  grid identifier
!! @param[in]  use_ll_edges             use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_unstruct_real ( grid_name,                &
                                         nbr_vertices,             &
                                         nbr_cells,                &
                                         nbr_vertices_per_cell_in, &
                                         x_vertices_real,          &
                                         y_vertices_real,          &
                                         cell_to_vertex_in,        &
                                         grid_id,                  &
                                         use_ll_edges )

  use yac, dummy => yac_fdef_grid_unstruct_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices
  integer, intent(in)  :: nbr_cells
  integer, intent(in)  :: nbr_vertices_per_cell_in

  real, intent(in)     :: x_vertices_real(nbr_vertices)
  real, intent(in)     :: y_vertices_real(nbr_vertices)

  integer, intent(in)  :: cell_to_vertex_in(nbr_vertices_per_cell_in,nbr_cells)

  integer, intent(out)   :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  double precision       :: x_vertices(nbr_vertices)
  double precision       :: y_vertices(nbr_vertices)

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_unstruct_dble ( grid_name,                &
                                     nbr_vertices,             &
                                     nbr_cells,                &
                                     nbr_vertices_per_cell_in, &
                                     x_vertices,               &
                                     y_vertices,               &
                                     cell_to_vertex_in,        &
                                     grid_id,                  &
                                     use_ll_edges )

end subroutine yac_fdef_grid_unstruct_real

!> Definition of a uniform unstructured grid (all cells have the
!! number of vertices)
!! @param[in]  grid_name                grid name
!! @param[in]  nbr_vertices             number of vertices
!! @param[in]  nbr_cells                number of cells
!! @param[in]  nbr_vertices_per_cell_in number of vertices for each cell
!! @param[in]  x_vertices               longitudes of vertices
!! @param[in]  y_vertices               latitudes of vertices
!! @param[in]  cell_to_vertex_in        connectivity between vertices and cells\n
!!                                      (the vertex indices per cell have to be
!!                                      in clockwise or counterclockwise ordering)
!! @param[out] grid_id                  grid identifier
!! @param[in]  use_ll_edges             use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_unstruct_dble ( grid_name,                &
                                         nbr_vertices,             &
                                         nbr_cells,                &
                                         nbr_vertices_per_cell_in, &
                                         x_vertices,               &
                                         y_vertices,               &
                                         cell_to_vertex_in,        &
                                         grid_id,                  &
                                         use_ll_edges )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_unstruct_dble

  implicit none

  interface

     subroutine yac_cdef_grid_unstruct_c ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_vertices_per_cell, &
                                           x_vertices,            &
                                           y_vertices,            &
                                           cell_to_vertex,        &
                                           grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(*)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_c

     subroutine yac_cdef_grid_unstruct_ll_c ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_vertices_per_cell, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_vertex,        &
                                              grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct_ll' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(*)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_ll_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices
  integer, intent(in)          :: nbr_cells
  integer, intent(in)          :: nbr_vertices_per_cell_in

  double precision, intent(in) :: x_vertices(nbr_vertices)
  double precision, intent(in) :: y_vertices(nbr_vertices)

  integer, intent(in)  :: cell_to_vertex_in( &
                            nbr_vertices_per_cell_in, &
                            nbr_cells)

  integer, intent(out)         :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  integer                      :: nbr_vertices_per_cell(nbr_cells)
  integer                      :: cell_to_vertex(nbr_vertices_per_cell_in,nbr_cells)
  logical                      :: use_ll_edges_

  nbr_vertices_per_cell(:) = nbr_vertices_per_cell_in

  YAC_FASSERT(ALL(cell_to_vertex_in > 0), "ERROR(yac_fdef_grid_unstruct_dble): all entries of cell_to_vertex have to be > 0")

  YAC_FASSERT(ALL(cell_to_vertex_in <= nbr_vertices), "ERROR(yac_fdef_grid_unstruct_dble): all entries of cell_to_vertex have to be <= nbr_vertices")

  cell_to_vertex(:,:) = cell_to_vertex_in(:,:) - 1

  if (present(use_ll_edges)) then
    use_ll_edges_ = use_ll_edges
  else
    use_ll_edges_ = .false.
  end if

  if (use_ll_edges_) then
    call yac_cdef_grid_unstruct_ll_c  ( TRIM(grid_name) // c_null_char, &
                                        nbr_vertices,                   &
                                        nbr_cells,                      &
                                        nbr_vertices_per_cell,          &
                                        x_vertices,                     &
                                        y_vertices,                     &
                                        cell_to_vertex,                 &
                                        grid_id )
  else
    call yac_cdef_grid_unstruct_c  ( TRIM(grid_name) // c_null_char, &
                                     nbr_vertices,                   &
                                     nbr_cells,                      &
                                     nbr_vertices_per_cell,          &
                                     x_vertices,                     &
                                     y_vertices,                     &
                                     cell_to_vertex,                 &
                                     grid_id )
  end if

end subroutine yac_fdef_grid_unstruct_dble

!> Definition of a non-uniform unstructured grid (cells have
!! varying numbers of vertices) with explicit edge definition
!! @param[in]  grid_name          grid name
!! @param[in]  nbr_vertices       number of vertices
!! @param[in]  nbr_cells          number of cells
!! @param[in]  nbr_edges          number of edges
!! @param[in]  nbr_connections    total size of cell_to_edge
!! @param[in]  nbr_edges_per_cell number of edges for each cell
!! @param[in]  x_vertices_real    longitudes of vertices
!! @param[in]  y_vertices_real    latitudes of vertices
!! @param[in]  cell_to_edge       connectivity between edges and cells\n
!!                                (the edge indices per cell have to be
!!                                in clockwise or counterclockwise ordering)
!! @param[in]  edge_to_vertex     connectivity between edges and vertices
!! @param[out] grid_id            grid identifier
!! @param[in]  use_ll_edges       use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_nonuniform_edge_real ( grid_name,          &
                                                nbr_vertices,       &
                                                nbr_cells,          &
                                                nbr_edges,          &
                                                nbr_connections,    &
                                                nbr_edges_per_cell, &
                                                x_vertices_real,    &
                                                y_vertices_real,    &
                                                cell_to_edge,       &
                                                edge_to_vertex,     &
                                                grid_id,            &
                                                use_ll_edges)

  use yac, dummy => yac_fdef_grid_nonuniform_edge_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices
  integer, intent(in)  :: nbr_cells
  integer, intent(in)  :: nbr_edges
  integer, intent(in)  :: nbr_connections
  integer, intent(in)  :: nbr_edges_per_cell(nbr_cells)

  real, intent(in)     :: x_vertices_real(nbr_vertices)
  real, intent(in)     :: y_vertices_real(nbr_vertices)

  integer, intent(in)  :: cell_to_edge(nbr_connections)
  integer, intent(in)  :: edge_to_vertex(2,nbr_edges)

  integer, intent(out) :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  double precision     :: x_vertices(nbr_vertices)
  double precision     :: y_vertices(nbr_vertices)

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_nonuniform_edge_dble ( grid_name,          &
                                            nbr_vertices,       &
                                            nbr_cells,          &
                                            nbr_edges,          &
                                            nbr_connections,    &
                                            nbr_edges_per_cell, &
                                            x_vertices,         &
                                            y_vertices,         &
                                            cell_to_edge,       &
                                            edge_to_vertex,     &
                                            grid_id,            &
                                            use_ll_edges )

end subroutine yac_fdef_grid_nonuniform_edge_real

!> Definition of a non-uniform unstructured grid (cells have
!! varying numbers of vertices) with explicit edge definition
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_edges             number of edges
!! @param[in]  nbr_connections       total size of cell_to_edge
!! @param[in]  nbr_edges_per_cell number of edges for each cell
!! @param[in]  x_vertices            longitudes of vertices
!! @param[in]  y_vertices            latitudes of vertices
!! @param[in]  cell_to_edge_in       connectivity between edges and cells\n
!!                                   (the edge indices per cell have to be
!!                                   in clockwise or counterclockwise ordering)
!! @param[in]  edge_to_vertex_in     connectivity between edges and vertices
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_nonuniform_edge_dble ( grid_name,          &
                                                nbr_vertices,       &
                                                nbr_cells,          &
                                                nbr_edges,          &
                                                nbr_connections,    &
                                                nbr_edges_per_cell, &
                                                x_vertices,         &
                                                y_vertices,         &
                                                cell_to_edge_in,    &
                                                edge_to_vertex_in,  &
                                                grid_id,            &
                                                use_ll_edges)

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_nonuniform_edge_dble

  implicit none

  interface

     subroutine yac_cdef_grid_unstruct_edge_c ( grid_name,          &
                                                nbr_vertices,       &
                                                nbr_cells,          &
                                                nbr_edges,          &
                                                nbr_edges_per_cell, &
                                                x_vertices,         &
                                                y_vertices,         &
                                                cell_to_edge,       &
                                                edge_to_vertex,     &
                                                grid_id )           &
           bind ( c, name='yac_cdef_grid_unstruct_edge' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    ), value        :: nbr_edges
       integer ( kind=c_int    )               :: nbr_edges_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_edge(*)
       integer ( kind=c_int    )               :: edge_to_vertex(2,nbr_edges)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_edge_c

     subroutine yac_cdef_grid_unstruct_edge_ll_c ( grid_name,             &
                                                   nbr_vertices,          &
                                                   nbr_cells,             &
                                                   nbr_edges,             &
                                                   nbr_edges_per_cell,    &
                                                   x_vertices,            &
                                                   y_vertices,            &
                                                   cell_to_edge,          &
                                                   edge_to_vertex,        &
                                                   grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct_edge_ll' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    ), value        :: nbr_edges
       integer ( kind=c_int    )               :: nbr_edges_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_edge(*)
       integer ( kind=c_int    )               :: edge_to_vertex(2,nbr_edges)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_edge_ll_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices
  integer, intent(in)          :: nbr_cells
  integer, intent(in)          :: nbr_edges
  integer, intent(in)          :: nbr_connections
  integer, intent(in)          :: nbr_edges_per_cell(nbr_cells)

  double precision, intent(in) :: x_vertices(nbr_vertices)
  double precision, intent(in) :: y_vertices(nbr_vertices)

  integer, intent(in)          :: cell_to_edge_in(nbr_connections)
  integer, intent(in)          :: edge_to_vertex_in(2,nbr_edges)

  integer, intent(out)         :: grid_id

  logical, optional, intent(in) :: use_ll_edges


  integer                      :: cell_to_edge(nbr_connections)
  integer                      :: edge_to_vertex(2,nbr_edges)

  logical :: use_ll_edges_

  YAC_FASSERT(ALL(cell_to_edge_in > 0), "ERROR(yac_fdef_grid_nonuniform_edge_dble): all entries of cell_to_edge have to be > 0")

  YAC_FASSERT(ALL(cell_to_edge_in <= nbr_edges), "ERROR(yac_fdef_grid_nonuniform_edge_dble): all entries of cell_to_edge have to be <= nbr_edges")

  YAC_FASSERT(ALL(edge_to_vertex_in > 0), "ERROR(yac_fdef_grid_nonuniform_edge_dble): all entries of edge_to_vertex have to be > 0")

  YAC_FASSERT(ALL(edge_to_vertex_in <= nbr_vertices), "ERROR(yac_fdef_grid_nonuniform_edge_dble): all entries of edge_to_vertex have to be <= nbr_vertices")

  cell_to_edge(:) = cell_to_edge_in(:) - 1
  edge_to_vertex(:,:) = edge_to_vertex_in(:,:) - 1

  if (present(use_ll_edges)) then
    use_ll_edges_ = use_ll_edges
  else
    use_ll_edges_ = .false.
  end if

  if (use_ll_edges_) then
    call yac_cdef_grid_unstruct_edge_ll_c  ( TRIM(grid_name) // c_null_char, &
                                             nbr_vertices,                   &
                                             nbr_cells,                      &
                                             nbr_edges,                      &
                                             nbr_edges_per_cell,             &
                                             x_vertices,                     &
                                             y_vertices,                     &
                                             cell_to_edge,                   &
                                             edge_to_vertex,                 &
                                             grid_id )
  else
    call yac_cdef_grid_unstruct_edge_c  ( TRIM(grid_name) // c_null_char, &
                                          nbr_vertices,                   &
                                          nbr_cells,                      &
                                          nbr_edges,                      &
                                          nbr_edges_per_cell,             &
                                          x_vertices,                     &
                                          y_vertices,                     &
                                          cell_to_edge,                   &
                                          edge_to_vertex,                 &
                                          grid_id )
  end if

end subroutine yac_fdef_grid_nonuniform_edge_dble

!> Definition of a uniform unstructured grid (all cells have the
!! number of vertices) with explicit edge definition
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_edges             number of edges
!! @param[in]  nbr_edges_per_cell_in number of edges for each cell
!! @param[in]  x_vertices_real       longitudes of vertices
!! @param[in]  y_vertices_real       latitudes of vertices
!! @param[in]  cell_to_edge_in       connectivity between edges and cells\n
!!                                   (the edge indices per cell have to be
!!                                   in clockwise or counterclockwise ordering)
!! @param[in]  edge_to_vertex_in     connectivity between edges and vertices
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_unstruct_edge_real ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_edges,             &
                                              nbr_edges_per_cell_in, &
                                              x_vertices_real,       &
                                              y_vertices_real,       &
                                              cell_to_edge_in,       &
                                              edge_to_vertex_in,     &
                                              grid_id,               &
                                              use_ll_edges)

  use yac, dummy => yac_fdef_grid_unstruct_edge_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices
  integer, intent(in)  :: nbr_cells
  integer, intent(in)  :: nbr_edges
  integer, intent(in)  :: nbr_edges_per_cell_in

  real, intent(in)     :: x_vertices_real(nbr_vertices)
  real, intent(in)     :: y_vertices_real(nbr_vertices)

  integer, intent(in)  :: cell_to_edge_in(nbr_edges_per_cell_in,nbr_cells)
  integer, intent(in)  :: edge_to_vertex_in(2,nbr_edges)

  integer, intent(out)   :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  double precision       :: x_vertices(nbr_vertices)
  double precision       :: y_vertices(nbr_vertices)

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_unstruct_edge_dble ( grid_name,             &
                                          nbr_vertices,          &
                                          nbr_cells,             &
                                          nbr_edges,             &
                                          nbr_edges_per_cell_in, &
                                          x_vertices,            &
                                          y_vertices,            &
                                          cell_to_edge_in,       &
                                          edge_to_vertex_in,     &
                                          grid_id,               &
                                          use_ll_edges )

end subroutine yac_fdef_grid_unstruct_edge_real

!> Definition of a uniform unstructured grid (all cells have the
!! number of vertices) with explicit edge definition
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_edges             number of edges
!! @param[in]  nbr_edges_per_cell_in number of edges for each cell
!! @param[in]  x_vertices            longitudes of vertices
!! @param[in]  y_vertices            latitudes of vertices
!! @param[in]  cell_to_edge_in       connectivity between edges and cells\n
!!                                   (the edge indices per cell have to be
!!                                   in clockwise or counterclockwise ordering)
!! @param[in]  edge_to_vertex_in     connectivity between edges and vertices
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_unstruct_edge_dble ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_edges,             &
                                              nbr_edges_per_cell_in, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_edge_in,       &
                                              edge_to_vertex_in,     &
                                              grid_id,               &
                                              use_ll_edges)

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_unstruct_edge_dble

  implicit none

  interface

     subroutine yac_cdef_grid_unstruct_edge_c ( grid_name,          &
                                                nbr_vertices,       &
                                                nbr_cells,          &
                                                nbr_edges,          &
                                                nbr_edges_per_cell, &
                                                x_vertices,         &
                                                y_vertices,         &
                                                cell_to_edge,       &
                                                edge_to_vertex,     &
                                                grid_id )           &
           bind ( c, name='yac_cdef_grid_unstruct_edge' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    ), value        :: nbr_edges
       integer ( kind=c_int    )               :: nbr_edges_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_edge(*)
       integer ( kind=c_int    )               :: edge_to_vertex(2,nbr_edges)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_edge_c

     subroutine yac_cdef_grid_unstruct_edge_ll_c ( grid_name,          &
                                                   nbr_vertices,       &
                                                   nbr_cells,          &
                                                   nbr_edges,          &
                                                   nbr_edges_per_cell, &
                                                   x_vertices,         &
                                                   y_vertices,         &
                                                   cell_to_edge,       &
                                                   edge_to_vertex,     &
                                                   grid_id )           &
           bind ( c, name='yac_cdef_grid_unstruct_edge_ll' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    ), value        :: nbr_edges
       integer ( kind=c_int    )               :: nbr_edges_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_edge(*)
       integer ( kind=c_int    )               :: edge_to_vertex(2,nbr_edges)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_edge_ll_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices
  integer, intent(in)          :: nbr_cells
  integer, intent(in)          :: nbr_edges
  integer, intent(in)          :: nbr_edges_per_cell_in

  double precision, intent(in) :: x_vertices(nbr_vertices)
  double precision, intent(in) :: y_vertices(nbr_vertices)

  integer, intent(in)  :: cell_to_edge_in(nbr_edges_per_cell_in, nbr_cells)
  integer, intent(in)  :: edge_to_vertex_in(2, nbr_edges)

  integer, intent(out)         :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  integer                      :: nbr_edges_per_cell(nbr_cells)
  integer                      :: cell_to_edge(nbr_edges_per_cell_in,nbr_cells)
  integer                      :: edge_to_vertex(2,nbr_edges)
  logical                      :: use_ll_edges_

  nbr_edges_per_cell(:) = nbr_edges_per_cell_in

  YAC_FASSERT(ALL(cell_to_edge_in > 0), "ERROR(yac_fdef_grid_unstruct_edge_dble): all entries of cell_to_edge have to be > 0")

  YAC_FASSERT(ALL(cell_to_edge_in <= nbr_edges), "ERROR(yac_fdef_grid_unstruct_edge_dble): all entries of cell_to_edge have to be <= nbr_edges")

  YAC_FASSERT(ALL(edge_to_vertex_in > 0), "ERROR(yac_fdef_grid_unstruct_edge_dble): all entries of edge_to_vertex have to be > 0")

  YAC_FASSERT(ALL(edge_to_vertex_in <= nbr_vertices), "ERROR(yac_fdef_grid_unstruct_edge_dble): all entries of edge_to_vertex have to be <= nbr_vertices")

  cell_to_edge(:,:) = cell_to_edge_in(:,:) - 1
  edge_to_vertex(:,:) = edge_to_vertex_in(:,:) - 1

  if (present(use_ll_edges)) then
    use_ll_edges_ = use_ll_edges
  else
    use_ll_edges_ = .false.
  end if

  if (use_ll_edges_) then
    call yac_cdef_grid_unstruct_edge_ll_c  ( TRIM(grid_name) // c_null_char, &
                                             nbr_vertices,                   &
                                             nbr_cells,                      &
                                             nbr_edges,                      &
                                             nbr_edges_per_cell,             &
                                             x_vertices,                     &
                                             y_vertices,                     &
                                             cell_to_edge,                   &
                                             edge_to_vertex,                 &
                                             grid_id )
  else
    call yac_cdef_grid_unstruct_edge_c  ( TRIM(grid_name) // c_null_char, &
                                          nbr_vertices,                   &
                                          nbr_cells,                      &
                                          nbr_edges,                      &
                                          nbr_edges_per_cell,             &
                                          x_vertices,                     &
                                          y_vertices,                     &
                                          cell_to_edge,                   &
                                          edge_to_vertex,                 &
                                          grid_id )
  end if

end subroutine yac_fdef_grid_unstruct_edge_dble

!> Definition of a 2d curvilinear grid
!! @param[in]  grid_name       grid name
!! @param[in]  nbr_vertices    number of cells in each dimension
!! @param[in]  cyclic          cyclic behavior of cells in each dimension
!! @param[in]  x_vertices_real longitudes of vertices
!! @param[in]  y_vertices_real latitudes of vertices
!! @param[out] grid_id         grid identifier
subroutine yac_fdef_grid_curve2d_real ( grid_name,       &
                                        nbr_vertices,    &
                                        cyclic,          &
                                        x_vertices_real, &
                                        y_vertices_real, &
                                        grid_id )

  use yac, dummy => yac_fdef_grid_curve2d_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices(2)
  integer, intent(in)  :: cyclic(2)
  real, intent(in)     :: &
    x_vertices_real(nbr_vertices(1),nbr_vertices(2))
  real, intent(in)     :: &
    y_vertices_real(nbr_vertices(1),nbr_vertices(2))
  integer, intent(out) :: grid_id

  double precision     ::  x_vertices(nbr_vertices(1),nbr_vertices(2))
  double precision     ::  y_vertices(nbr_vertices(1),nbr_vertices(2))

  x_vertices(:,:) = dble(x_vertices_real(:,:))
  y_vertices(:,:) = dble(y_vertices_real(:,:))

  call yac_fdef_grid_curve2d_dble ( grid_name,    &
                                    nbr_vertices, &
                                    cyclic,       &
                                    x_vertices,   &
                                    y_vertices,   &
                                    grid_id )

end subroutine yac_fdef_grid_curve2d_real

!> Definition of a 2d curvilinear grid
!! @param[in]  grid_name    grid name
!! @param[in]  nbr_vertices number of cells in each dimension
!! @param[in]  cyclic       cyclic behavior of cells in each dimension
!! @param[in]  x_vertices   longitudes of vertices
!! @param[in]  y_vertices   latitudes of vertices
!! @param[out] grid_id      grid identifier
subroutine yac_fdef_grid_curve2d_dble ( grid_name,    &
                                        nbr_vertices, &
                                        cyclic,       &
                                        x_vertices,   &
                                        y_vertices,   &
                                        grid_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_curve2d_dble

  implicit none

  interface

     subroutine yac_cdef_grid_curve2d_c ( grid_name,    &
                                          nbr_vertices, &
                                          cyclic,       &
                                          x_vertices,   &
                                          y_vertices,   &
                                          grid_id )     &
           bind ( c, name='yac_cdef_grid_curve2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    )               :: nbr_vertices(2)
       integer ( kind=c_int    )               :: cyclic(2)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices(1),nbr_vertices(2))
       real    ( kind=c_double )               :: y_vertices(nbr_vertices(1),nbr_vertices(2))
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_curve2d_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices(2)
  integer, intent(in)          :: cyclic(2)
  double precision, intent(in) :: &
    x_vertices(nbr_vertices(1),nbr_vertices(2))
  double precision, intent(in) :: &
    y_vertices(nbr_vertices(1),nbr_vertices(2))
  integer, intent(out)         :: grid_id

  call yac_cdef_grid_curve2d_c ( TRIM(grid_name) // c_null_char, &
                                 nbr_vertices,                   &
                                 cyclic,                         &
                                 x_vertices,                     &
                                 y_vertices,                     &
                                 grid_id )

end subroutine yac_fdef_grid_curve2d_dble

!> Definition of a 2d regular grid
!! @param[in]  grid_name       grid name
!! @param[in]  nbr_vertices    number of cells in each dimension
!! @param[in]  cyclic          cyclic behavior of cells in each dimension
!! @param[in]  x_vertices_real longitudes of vertices
!! @param[in]  y_vertices_real latitudes of vertices
!! @param[out] grid_id         grid identifier
subroutine yac_fdef_grid_reg2d_real ( grid_name,       &
                                      nbr_vertices,    &
                                      cyclic,          &
                                      x_vertices_real, &
                                      y_vertices_real, &
                                      grid_id )

  use yac, dummy => yac_fdef_grid_reg2d_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices(2)
  integer, intent(in)  :: cyclic(2)
  real, intent(in)     :: x_vertices_real(nbr_vertices(1))
  real, intent(in)     :: y_vertices_real(nbr_vertices(2))
  integer, intent(out) :: grid_id

  double precision     ::  x_vertices(nbr_vertices(1))
  double precision     ::  y_vertices(nbr_vertices(2))

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_reg2d_dble ( grid_name,    &
                                  nbr_vertices, &
                                  cyclic,       &
                                  x_vertices,   &
                                  y_vertices,   &
                                  grid_id )

end subroutine yac_fdef_grid_reg2d_real

!> Definition of a 2d regular grid
!! @param[in]  grid_name    grid name
!! @param[in]  nbr_vertices number of cells in each dimension
!! @param[in]  cyclic       cyclic behavior of cells in each dimension
!! @param[in]  x_vertices   longitudes of vertices
!! @param[in]  y_vertices   latitudes of vertices
!! @param[out] grid_id      grid identifier
subroutine yac_fdef_grid_reg2d_dble ( grid_name,    &
                                      nbr_vertices, &
                                      cyclic,       &
                                      x_vertices,   &
                                      y_vertices,   &
                                      grid_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_reg2d_dble

  implicit none

  interface

     subroutine yac_cdef_grid_reg2d_c ( grid_name,    &
                                        nbr_vertices, &
                                        cyclic,       &
                                        x_vertices,   &
                                        y_vertices,   &
                                        grid_id )     &
           bind ( c, name='yac_cdef_grid_reg2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    )               :: nbr_vertices(2)
       integer ( kind=c_int    )               :: cyclic(2)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices(1))
       real    ( kind=c_double )               :: y_vertices(nbr_vertices(2))
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_reg2d_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices(2)
  integer, intent(in)          :: cyclic(2)
  double precision, intent(in) :: x_vertices(nbr_vertices(1))
  double precision, intent(in) :: y_vertices(nbr_vertices(2))
  integer, intent(out)         :: grid_id

  call yac_cdef_grid_reg2d_c ( TRIM(grid_name) // c_null_char, &
                               nbr_vertices,                   &
                               cyclic,                         &
                               x_vertices,                     &
                               y_vertices,                     &
                               grid_id )

end subroutine yac_fdef_grid_reg2d_dble

!> Definition of a grid consisting of a cloud of points
!! @param[in]  grid_name     grid name
!! @param[in]  nbr_points    number of points in each dimension
!! @param[in]  x_points_real longitudes of points
!! @param[in]  y_points_real latitudes of points
!! @param[out] grid_id       grid identifier
subroutine yac_fdef_grid_cloud_real ( grid_name,     &
                                      nbr_points,    &
                                      x_points_real, &
                                      y_points_real, &
                                      grid_id)

  use yac, dummy => yac_fdef_grid_cloud_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_points
  real, intent(in)     :: x_points_real(nbr_points)
  real, intent(in)     :: y_points_real(nbr_points)
  integer, intent(out) :: grid_id

  double precision     ::  x_points(nbr_points)
  double precision     ::  y_points(nbr_points)

  x_points(:) = dble(x_points_real(:))
  y_points(:) = dble(y_points_real(:))

  call yac_fdef_grid_cloud_dble ( grid_name,  &
                                  nbr_points, &
                                  x_points,   &
                                  y_points,   &
                                  grid_id )

end subroutine yac_fdef_grid_cloud_real

!> Definition of a grid consisting of a cloud of points
!! @param[in]  grid_name  grid name
!! @param[in]  nbr_points number of points in each dimension
!! @param[in]  x_points   longitudes of points
!! @param[in]  y_points   latitudes of points
!! @param[out] grid_id    grid identifier
subroutine yac_fdef_grid_cloud_dble ( grid_name,  &
                                      nbr_points, &
                                      x_points,   &
                                      y_points,   &
                                      grid_id)

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_cloud_dble

  implicit none

  interface

     subroutine yac_fdef_grid_cloud_c ( grid_name,  &
                                        nbr_points, &
                                        x_points,   &
                                        y_points,   &
                                        grid_id )   &
           bind ( c, name='yac_cdef_grid_cloud' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_points
       real    ( kind=c_double )               :: x_points(nbr_points)
       real    ( kind=c_double )               :: y_points(nbr_points)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_fdef_grid_cloud_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_points
  double precision, intent(in) :: x_points(nbr_points)
  double precision, intent(in) :: y_points(nbr_points)
  integer, intent(out)         :: grid_id

  call yac_fdef_grid_cloud_c ( TRIM(grid_name) // c_null_char, &
                               nbr_points,                     &
                               x_points,                       &
                               y_points,                       &
                               grid_id )

end subroutine yac_fdef_grid_cloud_dble

!> Definition of a 2d regular rotated grid
!! @param[in]  grid_name         grid name
!! @param[in]  nbr_vertices      number of cells in each dimension
!! @param[in]  cyclic            cyclic behavior of cells in each dimension
!! @param[in]  x_vertices_real   longitudes of vertices
!! @param[in]  y_vertices_real   latitudes of vertices
!! @param[in]  x_north_pole_real longitude of north pole in radians
!! @param[in]  y_north_pole_real latitude of north pole in radians
!! @param[out] grid_id           grid identifier
subroutine yac_fdef_grid_reg2d_rot_real ( grid_name,         &
                                          nbr_vertices,      &
                                          cyclic,            &
                                          x_vertices_real,   &
                                          y_vertices_real,   &
                                          x_north_pole_real, &
                                          y_north_pole_real, &
                                          grid_id )

  use yac, dummy => yac_fdef_grid_reg2d_rot_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices(2)
  integer, intent(in)  :: cyclic(2)
  real, intent(in)     :: x_vertices_real(nbr_vertices(1))
  real, intent(in)     :: y_vertices_real(nbr_vertices(2))
  real, intent(in)     :: x_north_pole_real
  real, intent(in)     :: y_north_pole_real
  integer, intent(out) :: grid_id

  double precision     ::  x_vertices(nbr_vertices(1))
  double precision     ::  y_vertices(nbr_vertices(2))
  double precision     ::  x_north_pole
  double precision     ::  y_north_pole

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  x_north_pole = dble(x_north_pole_real)
  y_north_pole = dble(y_north_pole_real)

  call yac_fdef_grid_reg2d_rot_dble ( grid_name,    &
                                      nbr_vertices, &
                                      cyclic,       &
                                      x_vertices,   &
                                      y_vertices,   &
                                      x_north_pole, &
                                      y_north_pole, &
                                      grid_id )

end subroutine yac_fdef_grid_reg2d_rot_real

!> Definition of a 2d regular rotated grid
!! @param[in]  grid_name    grid name
!! @param[in]  nbr_vertices number of cells in each dimension
!! @param[in]  cyclic       cyclic behavior of cells in each dimension
!! @param[in]  x_vertices   longitudes of vertices
!! @param[in]  y_vertices   latitudes of vertices
!! @param[in]  x_north_pole longitude of north pole in radians
!! @param[in]  y_north_pole latitude of north pole in radians
!! @param[out] grid_id      grid identifier
subroutine yac_fdef_grid_reg2d_rot_dble ( grid_name,    &
                                          nbr_vertices, &
                                          cyclic,       &
                                          x_vertices,   &
                                          y_vertices,   &
                                          x_north_pole, &
                                          y_north_pole, &
                                          grid_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_grid_reg2d_rot_dble

  implicit none

  interface

     subroutine yac_cdef_grid_reg2d_rot_c ( grid_name,    &
                                            nbr_vertices, &
                                            cyclic,       &
                                            x_vertices,   &
                                            y_vertices,   &
                                            x_north_pole, &
                                            y_north_pole, &
                                            grid_id )     &
           bind ( c, name='yac_cdef_grid_reg2d_rot' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    )               :: nbr_vertices(2)
       integer ( kind=c_int    )               :: cyclic(2)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices(1))
       real    ( kind=c_double )               :: y_vertices(nbr_vertices(2))
       real    ( kind=c_double ), value        :: x_north_pole
       real    ( kind=c_double ), value        :: y_north_pole
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_reg2d_rot_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices(2)
  integer, intent(in)          :: cyclic(2)
  double precision, intent(in) :: x_vertices(nbr_vertices(1))
  double precision, intent(in) :: y_vertices(nbr_vertices(2))
  double precision, intent(in) :: x_north_pole
  double precision, intent(in) :: y_north_pole
  integer, intent(out)         :: grid_id

  call yac_cdef_grid_reg2d_rot_c ( TRIM(grid_name) // c_null_char, &
                                   nbr_vertices,                   &
                                   cyclic,                         &
                                   x_vertices,                     &
                                   y_vertices,                     &
                                   x_north_pole,                   &
                                   y_north_pole,                   &
                                   grid_id )

end subroutine yac_fdef_grid_reg2d_rot_dble

! ------------------------- set global_index ------------------------------

subroutine yac_fset_global_index ( global_index, &
                                   location,     &
                                   grid_id )

  use yac, dummy => yac_fset_global_index

  implicit none

  interface

     subroutine yac_cset_global_index_c ( global_index, &
                                          location,     &
                                          grid_id )     &
           bind ( c, name='yac_cset_global_index' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int )        :: global_index(*)
       integer ( kind=c_int ), value :: location
       integer ( kind=c_int ), value :: grid_id

     end subroutine yac_cset_global_index_c

  end interface

  integer, intent(in)  :: global_index(*) !< [IN] global index
  integer, intent(in)  :: location        !< [IN] location
  integer, intent(in)  :: grid_id         !< [IN] point set indentifier

  call yac_cset_global_index_c ( global_index, &
                                 location,     &
                                 grid_id )

end subroutine yac_fset_global_index

! ------------------------- set core_mask ------------------------------

subroutine yac_fset_core_lmask ( is_core,  &
                                 location, &
                                 grid_id )

  use yac, dummy => yac_fset_core_lmask
  use, intrinsic :: iso_c_binding, only : c_size_t

  implicit none

  logical, intent(in)  :: is_core(*) !< [IN] core flag
                                     !< false, cell/vertex/edge is halo
                                     !< true, cell/vertex/edge is core
  integer, intent(in)  :: location   !< [IN] location
  integer, intent(in)  :: grid_id    !< [IN] point set indentifier

  integer (kind=c_size_t) :: i, count
  integer, allocatable :: int_is_core(:)

  interface

     function yac_cget_grid_size_c ( location, &
                                    grid_id ) &
           result ( grid_size )               &
           bind ( c, name='yac_cget_grid_size' )

       use, intrinsic :: iso_c_binding, only : c_int, c_size_t

       integer ( kind=c_int ), value :: location
       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_size_t) :: grid_size

     end function yac_cget_grid_size_c

  end interface

  count = yac_cget_grid_size_c(location, grid_id)
  allocate(int_is_core(count))

  do i = 1, count
     if ( is_core(i) ) then
        int_is_core(i) = 1
     else
        int_is_core(i) = 0
     endif
  enddo

  call yac_fset_core_imask ( int_is_core, &
                             location,    &
                             grid_id )

end subroutine yac_fset_core_lmask

subroutine yac_fset_core_imask ( is_core,  &
                                 location, &
                                 grid_id )

  use yac, dummy => yac_fset_core_imask

  implicit none

  interface

     subroutine yac_cset_core_mask_c ( mask,        &
                                       location,    &
                                       grid_id )    &
           bind ( c, name='yac_cset_core_mask' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int )        :: mask(*)
       integer ( kind=c_int ), value :: location
       integer ( kind=c_int ), value :: grid_id

     end subroutine yac_cset_core_mask_c

  end interface


  integer, intent(in)  :: is_core(*) !< [IN] core flag
                                     !< 0, cell/vertex/edge is halo
                                     !< 1, cell/vertex/edge is core
  integer, intent(in)  :: location   !< [IN] location
  integer, intent(in)  :: grid_id    !< [IN] point set indentifier

  call yac_cset_core_mask_c ( is_core,  &
                              location, &
                              grid_id )

end subroutine yac_fset_core_imask

! ------------------------- set_mask ------------------------------

subroutine yac_fset_lmask ( is_valid,     &
                            points_id )

  use yac, dummy => yac_fset_lmask
  use, intrinsic :: iso_c_binding, only : c_size_t

  implicit none

  logical, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                       !< false, point is masked out
                                       !< true, point is valid
  integer, intent(in)  :: points_id    !< [IN] point set indentifier

  integer ( kind=c_size_t) :: i, count
  integer, allocatable :: int_is_valid(:)

  interface

     function yac_cget_points_size_c ( points_id ) &
           result ( points_size )                 &
           bind ( c, name='yac_cget_points_size' )

       use, intrinsic :: iso_c_binding, only : c_int, c_size_t

       integer ( kind=c_int ), value :: points_id
       integer ( kind=c_size_t) :: points_size

     end function yac_cget_points_size_c

  end interface

  count = yac_cget_points_size_c(points_id)
  allocate(int_is_valid(count))

  do i = 1, count
    int_is_valid(i) = MERGE(1,0,is_valid(i))
  enddo

  call yac_fset_imask ( int_is_valid, points_id )

end subroutine yac_fset_lmask

subroutine yac_fset_imask ( is_valid, &
                            points_id )

  use yac, dummy => yac_fset_imask

  implicit none

  interface

     subroutine yac_cset_mask_c ( is_valid,   &
                                  points_id ) &
           bind ( c, name='yac_cset_mask' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int )        :: is_valid(*)
       integer ( kind=c_int ), value :: points_id

     end subroutine yac_cset_mask_c

  end interface

  integer, intent(in)  :: is_valid(*)  !< [IN] integer mask
                                       !< 0, point is masked out
                                       !< 1, point is valid
  integer, intent(in)  :: points_id    !< [IN] point set indentifier

  call yac_cset_mask_c ( is_valid, points_id )

end subroutine yac_fset_imask

! ------------------------- def_mask ------------------------------

subroutine yac_fdef_lmask ( grid_id,    &
                            nbr_points, &
                            location,   &
                            is_valid,   &
                            mask_id )

  use yac, dummy => yac_fdef_lmask

  implicit none

  integer, intent(in)  :: grid_id     !< [IN] grid identifier
  integer, intent(in)  :: nbr_points  !< [IN] number of points
  integer, intent(in)  :: location    !< [IN] location, one of center/edge/vertex
  logical, intent(in)  :: is_valid(*) !< [IN] logical mask
                                      !< false, point is masked out
                                      !< true, point is valid
  integer, intent(out) :: mask_id     !< [OUT] mask identifier

  integer :: i
  integer, allocatable :: int_is_valid(:)

  allocate(int_is_valid(nbr_points))

  do i = 1, nbr_points
     int_is_valid(i) = MERGE(1,0,is_valid(i))
  enddo

  call yac_fdef_imask ( grid_id,      &
                        nbr_points,   &
                        location,     &
                        int_is_valid, &
                        mask_id )

end subroutine yac_fdef_lmask

subroutine yac_fdef_imask ( grid_id,    &
                            nbr_points, &
                            location,   &
                            is_valid,   &
                            mask_id )

  use yac, dummy => yac_fdef_imask

  implicit none

  interface

     subroutine yac_cdef_mask_c ( grid_id,    &
                                  nbr_points, &
                                  location,   &
                                  is_valid,   &
                                  mask_id )   &
           bind ( c, name='yac_cdef_mask' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_int ), value :: nbr_points
       integer ( kind=c_int ), value :: location
       integer ( kind=c_int )        :: is_valid(*)
       integer ( kind=c_int )        :: mask_id

     end subroutine yac_cdef_mask_c

  end interface

  integer, intent(in)  :: grid_id     !< [IN] grid identifier
  integer, intent(in)  :: nbr_points  !< [IN] number of points
  integer, intent(in)  :: location    !< [IN] location, one of center/edge/vertex
  integer, intent(in)  :: is_valid(*) !< [IN] logical mask
                                      !< false, point is masked out
                                      !< true, point is valid
  integer, intent(out) :: mask_id     !< [OUT] mask identifier

  call yac_cdef_mask_c ( grid_id,    &
                         nbr_points, &
                         location,   &
                         is_valid,   &
                         mask_id )

end subroutine yac_fdef_imask

! ------------------------- def_mask_named ------------------------

subroutine yac_fdef_lmask_named ( grid_id,    &
                                  nbr_points, &
                                  location,   &
                                  is_valid,   &
                                  name,       &
                                  mask_id )

  use yac, dummy => yac_fdef_lmask_named

  implicit none

  integer, intent(in)  :: grid_id      !< [IN] grid identifier
  integer, intent(in)  :: nbr_points   !< [IN] number of points
  integer, intent(in)  :: location     !< [IN] location, one of center/edge/vertex
  logical, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                       !< false, point is masked out
                                       !< true, point is valid
  character(len=*), intent(in) :: name !< [IN] name of the mask
  integer, intent(out) :: mask_id      !< [OUT] mask identifier

  integer :: i
  integer, allocatable :: int_is_valid(:)

  allocate(int_is_valid(nbr_points))

  do i = 1, nbr_points
     int_is_valid(i) = MERGE(1,0,is_valid(i))
  enddo

  call yac_fdef_imask_named ( grid_id,      &
                              nbr_points,   &
                              location,     &
                              int_is_valid, &
                              name,         &
                              mask_id )

end subroutine yac_fdef_lmask_named

subroutine yac_fdef_imask_named ( grid_id,    &
                                  nbr_points, &
                                  location,   &
                                  is_valid,   &
                                  name,       &
                                  mask_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_imask_named
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_mask_named_c ( grid_id,    &
                                        nbr_points, &
                                        location,   &
                                        is_valid,   &
                                        name,       &
                                        mask_id )   &
           bind ( c, name='yac_cdef_mask_named' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer ( kind=c_int ), value           :: grid_id
       integer ( kind=c_int ), value           :: nbr_points
       integer ( kind=c_int ), value           :: location
       integer ( kind=c_int )                  :: is_valid(*)
       character ( kind=c_char ), dimension(*) :: name
       integer ( kind=c_int )                  :: mask_id

     end subroutine yac_cdef_mask_named_c

  end interface

  integer, intent(in)  :: grid_id      !< [IN] grid identifier
  integer, intent(in)  :: nbr_points   !< [IN] number of points
  integer, intent(in)  :: location     !< [IN] location, one of center/edge/vertex
  integer, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                       !< false, point is masked out
                                       !< true, point is valid
  character(len=*), intent(in) :: name !< [IN] name of the mask
  integer, intent(out) :: mask_id      !< [OUT] mask identifier

  YAC_CHECK_STRING_LEN ( "yac_fdef_imask_named", name )

  call yac_cdef_mask_named_c ( grid_id,                   &
                               nbr_points,                &
                               location,                  &
                               is_valid,                  &
                               TRIM(name) // c_null_char, &
                               mask_id )

end subroutine yac_fdef_imask_named

! ----------------------------- def_field -------------------------------

subroutine yac_fdef_field ( field_name,      &
                            component_id,    &
                            point_ids,       &
                            num_pointsets,   &
                            collection_size, &
                            timestep,        &
                            time_unit,       &
                            field_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_field
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_field_c ( field_name,      &
                                   component_id,    &
                                   point_ids,       &
                                   num_pointsets,   &
                                   collection_size, &
                                   timestep,        &
                                   time_unit,       &
                                   field_id )       &
           bind ( c, name='yac_cdef_field' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: field_name !< [IN] short name of the field
       integer ( kind=c_int ), value :: component_id         !< [IN] component identifier
       integer ( kind=c_int )        :: point_ids(*)         !< [IN] point identifier
       integer ( kind=c_int ), value :: num_pointsets        !< [IN] number of pointsets per grid
       integer ( kind=c_int ), value :: collection_size      !< [IN] collection_size
       character ( kind=c_char ), dimension(*) :: timestep   !< [IN] timestep
       integer ( kind=c_int ), value :: time_unit            !< [IN] unit of timestep
       integer ( kind=c_int )        :: field_id             !< [OUT] returned field handle

     end subroutine yac_cdef_field_c

  end interface

  !
  ! Definition of coupling fields
  !
  character(len=*), intent (in) :: field_name      !< [IN] short name of the field
  integer, intent (in)          :: component_id    !< [IN] component identifier
  integer, intent (in)          :: point_ids(*)    !< [IN] point identifier
  integer, intent (in)          :: num_pointsets   !< [IN] number of pointsets per grid
  integer, intent (in)          :: collection_size !< [IN] collection size
  character(len=*), intent (in) :: timestep        !< [IN] timestep
  integer, intent (in)          :: time_unit       !< [IN] unit of timestep
  integer, intent (out)         :: field_id        !< [OUT] returned field handle

  YAC_CHECK_STRING_LEN ( "yac_fdef_field", field_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_field", timestep )


  call yac_cdef_field_c ( TRIM(field_name) // c_null_char, &
                          component_id,                    &
                          point_ids,                       &
                          num_pointsets,                   &
                          collection_size,                 &
                          TRIM(timestep) // c_null_char,   &
                          time_unit,                       &
                          field_id )

end subroutine yac_fdef_field

! ----------------------------- def_field_mask---------------------------

subroutine yac_fdef_field_mask ( field_name,      &
                                 component_id,    &
                                 point_ids,       &
                                 mask_ids,        &
                                 num_pointsets,   &
                                 collection_size, &
                                 timestep,        &
                                 time_unit,       &
                                 field_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fdef_field_mask
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_field_mask_c ( field_name,      &
                                        component_id,    &
                                        point_ids,       &
                                        mask_ids,        &
                                        num_pointsets,   &
                                        collection_size, &
                                        timestep,        &
                                        time_unit,       &
                                        field_id )       &
           bind ( c, name='yac_cdef_field_mask' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int ), value :: component_id
       integer ( kind=c_int )        :: point_ids(*)
       integer ( kind=c_int )        :: mask_ids(*)
       integer ( kind=c_int ), value :: num_pointsets
       integer ( kind=c_int ), value :: collection_size
       character ( kind=c_char ), dimension(*) :: timestep
       integer ( kind=c_int ), value :: time_unit
       integer ( kind=c_int )        :: field_id

     end subroutine yac_cdef_field_mask_c

  end interface


  !
  ! Definition of coupling fields
  !
  character(len=*), intent (in) :: field_name      !< [IN] short name of the field
  integer, intent (in)          :: component_id    !< [IN] component identifier
  integer, intent (in)          :: point_ids(*)    !< [IN] point identifier
  integer, intent (in)          :: mask_ids(*)     !< [IN] mask identifier
  integer, intent (in)          :: num_pointsets   !< [IN] number of pointsets per grid
  integer, intent (in)          :: collection_size !< [IN] collection size
  character(len=*), intent (in) :: timestep        !< [IN] timestep
  integer, intent (in)          :: time_unit       !< [IN] unit of timestep
  integer, intent (out)         :: field_id        !< [OUT] returned field handle

  YAC_CHECK_STRING_LEN ( "yac_fdef_field_mask", field_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_field_mask", timestep )

  call yac_cdef_field_mask_c ( TRIM(field_name) // c_null_char, &
                               component_id,                    &
                               point_ids,                       &
                               mask_ids,                        &
                               num_pointsets,                   &
                               collection_size,                 &
                               TRIM(timestep) // c_null_char,   &
                               time_unit,                       &
                               field_id )

end subroutine yac_fdef_field_mask

! -----------------------------------------------------------------------

subroutine yac_fcheck_field_dimensions( field_id,          &
                                        collection_size,   &
                                        num_interp_fields, &
                                        interp_field_sizes )

  use yac, dummy => yac_fcheck_field_dimensions

  implicit none

  interface

     subroutine yac_ccheck_field_dimensions_c ( field_id,            &
                                                collection_size,     &
                                                num_interp_fields,   &
                                                interp_field_sizes ) &
       bind ( c, name='yac_ccheck_field_dimensions' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       integer ( kind=c_int ), value :: num_interp_fields
       integer ( kind=c_int ), dimension(*) :: interp_field_sizes

     end subroutine yac_ccheck_field_dimensions_c

  end interface

  integer, intent (in) :: field_id                              !<[IN] field handle
  integer, intent (in) :: collection_size                       !<[IN] collection size
  integer, intent (in) :: num_interp_fields                     !<[IN] number of interpolation fields
                                                                !!     (number of pointsets)
  integer, intent (in) :: interp_field_sizes(num_interp_fields) !<[IN] data size of each
                                                                !!     interpolation field

  call yac_ccheck_field_dimensions_c(field_id,          &
                                     collection_size,   &
                                     num_interp_fields, &
                                     interp_field_sizes)

end subroutine yac_fcheck_field_dimensions

subroutine yac_fcheck_src_field_buffer_size( field_id,          &
                                             collection_size,   &
                                             src_field_buffer_size )

  use yac, dummy => yac_fcheck_src_field_buffer_size

  implicit none

  interface

     subroutine yac_ccheck_src_field_buffer_size_c ( field_id,               &
                                                     collection_size,        &
                                                     src_field_buffer_size ) &
       bind ( c, name='yac_ccheck_src_field_buffer_size' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       integer ( kind=c_int ), value :: src_field_buffer_size

     end subroutine yac_ccheck_src_field_buffer_size_c

  end interface

  integer, intent (in) :: field_id              !<[IN] field handle
  integer, intent (in) :: collection_size       !<[IN] collection size
  integer, intent (in) :: src_field_buffer_size !<[IN] source field buffer size
                                                !!     (SUM(src_field_buffer_sizes(:))

  call yac_ccheck_src_field_buffer_size_c(field_id,          &
                                          collection_size,   &
                                          src_field_buffer_size)

end subroutine yac_fcheck_src_field_buffer_size

subroutine yac_fcheck_src_field_buffer_sizes( field_id,        &
                                              num_src_fields,  &
                                              collection_size, &
                                              src_field_buffer_sizes )

  use yac, dummy => yac_fcheck_src_field_buffer_sizes

  implicit none

  interface

     subroutine yac_ccheck_src_field_buffer_sizes_c ( field_id,                &
                                                      num_src_fields,          &
                                                      collection_size,         &
                                                      src_field_buffer_sizes ) &
       bind ( c, name='yac_ccheck_src_field_buffer_sizes' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: num_src_fields
       integer ( kind=c_int ), value :: collection_size
       integer ( kind=c_int )        :: src_field_buffer_sizes(*)

     end subroutine yac_ccheck_src_field_buffer_sizes_c

  end interface

  integer, intent (in) :: field_id        !<[IN] field handle
  integer, intent (in) :: num_src_fields  !<[IN] number of source fields
  integer, intent (in) :: collection_size !<[IN] collection size
  integer, intent (in) :: src_field_buffer_sizes(num_src_fields)
                                          !<[IN] source field buffer size
                                          !!     (SUM(src_field_buffer_sizes(:))

  call yac_ccheck_src_field_buffer_sizes_c(field_id,        &
                                           num_src_fields,  &
                                           collection_size, &
                                           src_field_buffer_sizes)

end subroutine yac_fcheck_src_field_buffer_sizes

! -----------------------------------------------------------------------

subroutine yac_fget_raw_interp_weights_data (  &
  field_id,                                    &
  frac_mask_fallback_value,                    &
  scaling_factor,                              &
  scaling_summand,                             &
  num_fixed_values,                            &
  fixed_values,                                &
  num_tgt_per_fixed_value,                     &
  tgt_idx_fixed,                               &
  num_wgt_tgt,                                 &
  wgt_tgt_idx,                                 &
  num_src_per_tgt,                             &
  weights,                                     &
  src_field_idx,                               &
  src_idx,                                     &
  num_src_fields,                              &
  src_field_buffer_size )

  use yac, dummy => yac_fget_raw_interp_weights_data

  use, intrinsic :: iso_c_binding, only : &
    c_int, c_size_t, c_double, c_ptr, c_f_pointer

  implicit none

  integer, intent (in)                       :: field_id
  double precision, intent(out)              :: frac_mask_fallback_value
  double precision, intent(out)              :: scaling_factor
  double precision, intent(out)              :: scaling_summand
  integer, intent(out)                       :: num_fixed_values
  double precision, allocatable, intent(out) :: fixed_values(:)
  integer, allocatable, intent(out)          :: num_tgt_per_fixed_value(:)
  integer, allocatable, intent(out)          :: tgt_idx_fixed(:)
  integer, intent(out)                       :: num_wgt_tgt
  integer, allocatable, intent(out)          :: wgt_tgt_idx(:)
  integer, allocatable, intent(out)          :: num_src_per_tgt(:)
  double precision, allocatable, intent(out) :: weights(:)
  integer, allocatable, intent(out)          :: src_field_idx(:)
  integer, allocatable, intent(out)          :: src_idx(:)
  integer, intent(out)                       :: num_src_fields
  integer, allocatable, intent(out)          :: src_field_buffer_size(:)

  integer :: total_num_fixed_tgt, total_num_weights

  real (kind=c_double)               :: frac_mask_fallback_value_c
  real (kind=c_double)               :: scaling_factor_c
  real (kind=c_double)               :: scaling_summand_c
  integer ( kind=c_size_t )          :: num_fixed_values_c
  real (kind=c_double), pointer      :: fixed_values_c(:)
  integer ( kind=c_size_t ), pointer :: num_tgt_per_fixed_value_c(:)
  integer ( kind=c_size_t ), pointer :: tgt_idx_fixed_c(:)
  integer(kind=c_size_t)             :: num_wgt_tgt_c
  integer ( kind=c_size_t ), pointer :: wgt_tgt_idx_c(:)
  integer ( kind=c_size_t ), pointer :: num_src_per_tgt_c(:)
  real (kind=c_double), pointer      :: weights_c(:)
  integer ( kind=c_size_t ), pointer :: src_field_idx_c(:)
  integer ( kind=c_size_t ), pointer :: src_idx_c(:)
  integer ( kind=c_size_t )          :: num_src_fields_c
  integer ( kind=c_size_t ), pointer :: src_field_buffer_size_c(:)

  type(c_ptr) :: fixed_values_c_ptr
  type(c_ptr) :: num_tgt_per_fixed_value_c_ptr
  type(c_ptr) :: tgt_idx_fixed_c_ptr
  type(c_ptr) :: wgt_tgt_idx_c_ptr
  type(c_ptr) :: num_src_per_tgt_c_ptr
  type(c_ptr) :: weights_c_ptr
  type(c_ptr) :: src_field_idx_c_ptr
  type(c_ptr) :: src_idx_c_ptr
  type(c_ptr) :: src_field_buffer_size_c_ptr

  interface

    subroutine yac_cget_raw_interp_weights_data_c ( &
      field_id,                                     &
      frac_mask_fallback_value,                     &
      scaling_factor,                               &
      scaling_summand,                              &
      num_fixed_values,                             &
      fixed_values,                                 &
      num_tgt_per_fixed_value,                      &
      tgt_idx_fixed,                                &
      num_wgt_tgt,                                  &
      wgt_tgt_idx,                                  &
      num_src_per_tgt,                              &
      weights,                                      &
      src_field_idx,                                &
      src_idx,                                      &
      num_src_fields,                               &
      src_field_buffer_size )                       &
      bind ( c, name='yac_cget_raw_interp_weights_data' )

      use, intrinsic :: iso_c_binding, only : c_int, c_size_t, c_double, c_ptr

      integer ( kind=c_int ), value       :: field_id
      real (kind=c_double), intent(out)   :: frac_mask_fallback_value
      real (kind=c_double), intent(out)   :: scaling_factor
      real (kind=c_double), intent(out)   :: scaling_summand
      integer(kind=c_size_t), intent(out) :: num_fixed_values
      type(c_ptr), intent(out)            :: fixed_values
      type(c_ptr), intent(out)            :: num_tgt_per_fixed_value
      type(c_ptr), intent(out)            :: tgt_idx_fixed
      integer(kind=c_size_t), intent(out) :: num_wgt_tgt
      type(c_ptr), intent(out)            :: wgt_tgt_idx
      type(c_ptr), intent(out)            :: num_src_per_tgt
      type(c_ptr), intent(out)            :: weights
      type(c_ptr), intent(out)            :: src_field_idx
      type(c_ptr), intent(out)            :: src_idx
      integer(kind=c_size_t), intent(out) :: num_src_fields
      type(c_ptr), intent(out)            :: src_field_buffer_size

    end subroutine yac_cget_raw_interp_weights_data_c

    subroutine free_c ( ptr ) BIND ( c, NAME='free' )

      use, intrinsic :: iso_c_binding, only : c_ptr

      type ( c_ptr ), intent(in), value :: ptr

    end subroutine free_c

  end interface

  ! get interpolation weight data through C interface
  CALL yac_cget_raw_interp_weights_data_c ( &
      INT(field_id, c_int),                 &
      frac_mask_fallback_value_c,           &
      scaling_factor_c,                     &
      scaling_summand_c,                    &
      num_fixed_values_c,                   &
      fixed_values_c_ptr,                   &
      num_tgt_per_fixed_value_c_ptr,        &
      tgt_idx_fixed_c_ptr,                  &
      num_wgt_tgt_c,                        &
      wgt_tgt_idx_c_ptr,                    &
      num_src_per_tgt_c_ptr,                &
      weights_c_ptr,                        &
      src_field_idx_c_ptr,                  &
      src_idx_c_ptr,                        &
      num_src_fields_c,                     &
      src_field_buffer_size_c_ptr )

  call c_f_pointer( &
    fixed_values_c_ptr, fixed_values_c, [int(num_fixed_values_c)])
  call c_f_pointer( &
    num_tgt_per_fixed_value_c_ptr, num_tgt_per_fixed_value_c, &
    [int(num_fixed_values_c)])
  total_num_fixed_tgt = int(sum(num_tgt_per_fixed_value_c))
  call c_f_pointer( &
    tgt_idx_fixed_c_ptr, tgt_idx_fixed_c, [int(total_num_fixed_tgt)])

  call c_f_pointer( &
    wgt_tgt_idx_c_ptr, wgt_tgt_idx_c, [int(num_wgt_tgt_c)])
  call c_f_pointer( &
    num_src_per_tgt_c_ptr, num_src_per_tgt_c, [int(num_wgt_tgt_c)])
  total_num_weights = int(sum(num_src_per_tgt_c))
  call c_f_pointer( &
    weights_c_ptr, weights_c, [int(total_num_weights)])
  call c_f_pointer( &
    src_field_idx_c_ptr, src_field_idx_c, [int(total_num_weights)])
  call c_f_pointer( &
    src_idx_c_ptr, src_idx_c, [int(total_num_weights)])

  call c_f_pointer( &
    src_field_buffer_size_c_ptr, src_field_buffer_size_c, &
    [int(num_src_fields_c)])

  ! convert C to Fortran data
  frac_mask_fallback_value = dble(frac_mask_fallback_value_c)
  scaling_factor = dble(scaling_factor_c)
  scaling_summand = dble(scaling_summand_c)
  num_fixed_values = int(num_fixed_values_c)
  allocate(fixed_values, source=dble(fixed_values_c))
  allocate(num_tgt_per_fixed_value, source=int(num_tgt_per_fixed_value_c))
  allocate(tgt_idx_fixed, source=int(tgt_idx_fixed_c)+1)
  num_wgt_tgt = int(num_wgt_tgt_c)
  allocate(wgt_tgt_idx, source=int(wgt_tgt_idx_c)+1)
  allocate(num_src_per_tgt, source=int(num_src_per_tgt_c))
  allocate(weights, source=dble(weights_c))
  allocate(src_field_idx, source=int(src_field_idx_c)+1)
  allocate(src_idx, source=int(src_idx_c)+1)
  num_src_fields = int(num_src_fields_c)
  allocate(src_field_buffer_size, source=int(src_field_buffer_size_c))

  ! free C arrays
  CALL free_c(fixed_values_c_ptr)
  CALL free_c(num_tgt_per_fixed_value_c_ptr)
  CALL free_c(tgt_idx_fixed_c_ptr)
  CALL free_c(wgt_tgt_idx_c_ptr)
  CALL free_c(num_src_per_tgt_c_ptr)
  CALL free_c(weights_c_ptr)
  CALL free_c(src_field_idx_c_ptr)
  CALL free_c(src_idx_c_ptr)
  CALL free_c(src_field_buffer_size_c_ptr)

end subroutine yac_fget_raw_interp_weights_data

subroutine yac_fget_raw_interp_weights_data_csr (  &
  field_id,                                        &
  frac_mask_fallback_value,                        &
  scaling_factor,                                  &
  scaling_summand,                                 &
  num_fixed_values,                                &
  fixed_values,                                    &
  num_tgt_per_fixed_value,                         &
  tgt_idx_fixed,                                   &
  src_indptr,                                      &
  weights,                                         &
  src_field_idx,                                   &
  src_idx,                                         &
  num_src_fields,                                  &
  src_field_buffer_size )

  use yac, dummy => yac_fget_raw_interp_weights_data_csr

  use, intrinsic :: iso_c_binding, only : &
    c_int, c_size_t, c_double, c_ptr, c_f_pointer

  implicit none

  integer, intent (in)                       :: field_id
  double precision, intent(out)              :: frac_mask_fallback_value
  double precision, intent(out)              :: scaling_factor
  double precision, intent(out)              :: scaling_summand
  integer, intent(out)                       :: num_fixed_values
  double precision, allocatable, intent(out) :: fixed_values(:)
  integer, allocatable, intent(out)          :: num_tgt_per_fixed_value(:)
  integer, allocatable, intent(out)          :: tgt_idx_fixed(:)
  integer, allocatable, intent(out)          :: src_indptr(:)
  double precision, allocatable, intent(out) :: weights(:)
  integer, allocatable, intent(out)          :: src_field_idx(:)
  integer, allocatable, intent(out)          :: src_idx(:)
  integer, intent(out)                       :: num_src_fields
  integer, allocatable, intent(out)          :: src_field_buffer_size(:)

  integer :: total_num_fixed_tgt, total_num_weights

  real (kind=c_double)               :: frac_mask_fallback_value_c
  real (kind=c_double)               :: scaling_factor_c
  real (kind=c_double)               :: scaling_summand_c
  integer ( kind=c_size_t )          :: num_fixed_values_c
  real (kind=c_double), pointer      :: fixed_values_c(:)
  integer ( kind=c_size_t ), pointer :: num_tgt_per_fixed_value_c(:)
  integer ( kind=c_size_t ), pointer :: tgt_idx_fixed_c(:)
  integer ( kind=c_size_t ), pointer :: src_indptr_c(:)
  real (kind=c_double), pointer      :: weights_c(:)
  integer ( kind=c_size_t ), pointer :: src_field_idx_c(:)
  integer ( kind=c_size_t ), pointer :: src_idx_c(:)
  integer ( kind=c_size_t )          :: num_src_fields_c
  integer ( kind=c_size_t ), pointer :: src_field_buffer_size_c(:)
  integer ( kind=c_size_t )          :: tgt_field_data_size_c

  type(c_ptr) :: fixed_values_c_ptr
  type(c_ptr) :: num_tgt_per_fixed_value_c_ptr
  type(c_ptr) :: tgt_idx_fixed_c_ptr
  type(c_ptr) :: src_indptr_c_ptr
  type(c_ptr) :: weights_c_ptr
  type(c_ptr) :: src_field_idx_c_ptr
  type(c_ptr) :: src_idx_c_ptr
  type(c_ptr) :: src_field_buffer_size_c_ptr

  interface

    subroutine yac_cget_raw_interp_weights_data_csr_c2f_c ( &
      field_id,                                             &
      frac_mask_fallback_value,                             &
      scaling_factor,                                       &
      scaling_summand,                                      &
      num_fixed_values,                                     &
      fixed_values,                                         &
      num_tgt_per_fixed_value,                              &
      tgt_idx_fixed,                                        &
      src_indptr,                                           &
      weights,                                              &
      src_field_idx,                                        &
      src_idx,                                              &
      num_src_fields,                                       &
      src_field_buffer_size,                                &
      tgt_field_data_size)                                  &
      bind ( c, name='yac_cget_raw_interp_weights_data_csr_c2f' )

      use, intrinsic :: iso_c_binding, only : c_int, c_size_t, c_double, c_ptr

      integer ( kind=c_int ), value       :: field_id
      real (kind=c_double), intent(out)   :: frac_mask_fallback_value
      real (kind=c_double), intent(out)   :: scaling_factor
      real (kind=c_double), intent(out)   :: scaling_summand
      integer(kind=c_size_t), intent(out) :: num_fixed_values
      type(c_ptr), intent(out)            :: fixed_values
      type(c_ptr), intent(out)            :: num_tgt_per_fixed_value
      type(c_ptr), intent(out)            :: tgt_idx_fixed
      type(c_ptr), intent(out)            :: src_indptr
      type(c_ptr), intent(out)            :: weights
      type(c_ptr), intent(out)            :: src_field_idx
      type(c_ptr), intent(out)            :: src_idx
      integer(kind=c_size_t), intent(out) :: num_src_fields
      type(c_ptr), intent(out)            :: src_field_buffer_size
      integer(kind=c_size_t), intent(out) :: tgt_field_data_size

    end subroutine yac_cget_raw_interp_weights_data_csr_c2f_c

    subroutine free_c ( ptr ) BIND ( c, NAME='free' )

      use, intrinsic :: iso_c_binding, only : c_ptr

      type ( c_ptr ), intent(in), value :: ptr

    end subroutine free_c

  end interface

  ! get interpolation weight data through C interface
  CALL yac_cget_raw_interp_weights_data_csr_c2f_c ( &
      INT(field_id, c_int),                         &
      frac_mask_fallback_value_c,                   &
      scaling_factor_c,                             &
      scaling_summand_c,                            &
      num_fixed_values_c,                           &
      fixed_values_c_ptr,                           &
      num_tgt_per_fixed_value_c_ptr,                &
      tgt_idx_fixed_c_ptr,                          &
      src_indptr_c_ptr,                             &
      weights_c_ptr,                                &
      src_field_idx_c_ptr,                          &
      src_idx_c_ptr,                                &
      num_src_fields_c,                             &
      src_field_buffer_size_c_ptr,                  &
      tgt_field_data_size_c)

  call c_f_pointer( &
    fixed_values_c_ptr, fixed_values_c, [int(num_fixed_values_c)])
  call c_f_pointer( &
    num_tgt_per_fixed_value_c_ptr, num_tgt_per_fixed_value_c, &
    [int(num_fixed_values_c)])
  total_num_fixed_tgt = int(sum(num_tgt_per_fixed_value_c))
  call c_f_pointer( &
    tgt_idx_fixed_c_ptr, tgt_idx_fixed_c, [int(total_num_fixed_tgt)])

  call c_f_pointer( &
    src_indptr_c_ptr, src_indptr_c, [int(tgt_field_data_size_c+1)])
  total_num_weights = int(src_indptr_c(tgt_field_data_size_c+1))
  call c_f_pointer( &
    weights_c_ptr, weights_c, [int(total_num_weights)])
  call c_f_pointer( &
    src_field_idx_c_ptr, src_field_idx_c, [int(total_num_weights)])
  call c_f_pointer( &
    src_idx_c_ptr, src_idx_c, [int(total_num_weights)])

  call c_f_pointer( &
    src_field_buffer_size_c_ptr, src_field_buffer_size_c, &
    [int(num_src_fields_c)])

  ! convert C to Fortran data
  frac_mask_fallback_value = dble(frac_mask_fallback_value_c)
  scaling_factor = dble(scaling_factor_c)
  scaling_summand = dble(scaling_summand_c)
  num_fixed_values = int(num_fixed_values_c)
  allocate(fixed_values, source=dble(fixed_values_c))
  allocate(num_tgt_per_fixed_value, source=int(num_tgt_per_fixed_value_c))
  allocate(tgt_idx_fixed, source=int(tgt_idx_fixed_c)+1)
  allocate(src_indptr, source=int(src_indptr_c)+1)
  allocate(weights, source=dble(weights_c))
  allocate(src_field_idx, source=int(src_field_idx_c)+1)
  allocate(src_idx, source=int(src_idx_c)+1)
  num_src_fields = int(num_src_fields_c)
  allocate(src_field_buffer_size, source=int(src_field_buffer_size_c))

  ! free C arrays
  CALL free_c(fixed_values_c_ptr)
  CALL free_c(num_tgt_per_fixed_value_c_ptr)
  CALL free_c(tgt_idx_fixed_c_ptr)
  CALL free_c(src_indptr_c_ptr)
  CALL free_c(weights_c_ptr)
  CALL free_c(src_field_idx_c_ptr)
  CALL free_c(src_idx_c_ptr)
  CALL free_c(src_field_buffer_size_c_ptr)

end subroutine yac_fget_raw_interp_weights_data_csr

! ---------------------------------- put --------------------------------
!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      returned field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_real ( field_id,         &
                           nbr_hor_points,   &
                           nbr_pointsets,    &
                           collection_size,  &
                           send_field,       &
                           info,             &
                           ierror )

  use yac, dummy => yac_fput_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id
  integer, intent (in)  :: nbr_hor_points
  integer, intent (in)  :: nbr_pointsets
  integer, intent (in)  :: collection_size
  real, intent (in)     :: send_field(nbr_hor_points, &
                                      nbr_pointsets,  &
                                      collection_size)
  integer, intent (out) :: info
  integer, intent (out) :: ierror

  double precision      :: send_field_dble(nbr_hor_points, &
                                           nbr_pointsets,  &
                                           collection_size)
  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )

  call send_field_to_dble(field_id,        &
                          nbr_hor_points,  &
                          nbr_pointsets,   &
                          collection_size, &
                          send_field,      &
                          send_field_dble)

  call yac_fput_dble ( field_id,        &
                       nbr_hor_points,  &
                       nbr_pointsets,   &
                       collection_size, &
                       send_field_dble, &
                       info,            &
                       ierror )

end subroutine yac_fput_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_real ( field_id,         &
                                nbr_hor_points,   &
                                nbr_pointsets,    &
                                collection_size,  &
                                send_field,       &
                                send_frac_mask,   &
                                info,             &
                                ierror )

  use yac, dummy => yac_fput_frac_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id
  integer, intent (in)  :: nbr_hor_points
  integer, intent (in)  :: nbr_pointsets
  integer, intent (in)  :: collection_size
  real, intent (in)     :: send_field(nbr_hor_points, &
                                      nbr_pointsets,  &
                                      collection_size)
  real, intent (in)     :: send_frac_mask(nbr_hor_points, &
                                          nbr_pointsets,  &
                                          collection_size)
  integer, intent (out) :: info
  integer, intent (out) :: ierror

  double precision      :: send_field_dble(nbr_hor_points, &
                                           nbr_pointsets,  &
                                           collection_size)
  double precision      :: send_frac_mask_dble(nbr_hor_points, &
                                               nbr_pointsets,  &
                                               collection_size)
  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )

  call send_field_to_dble(field_id,        &
                          nbr_hor_points,  &
                          nbr_pointsets,   &
                          collection_size, &
                          send_field,      &
                          send_field_dble, &
                          send_frac_mask,  &
                          send_frac_mask_dble)

  call yac_fput_frac_dble ( field_id,            &
                            nbr_hor_points,      &
                            nbr_pointsets,       &
                            collection_size,     &
                            send_field_dble,     &
                            send_frac_mask_dble, &
                            info,                &
                            ierror )

end subroutine yac_fput_frac_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_real_ptr ( field_id,        &
                               nbr_pointsets,   &
                               collection_size, &
                               send_field,      &
                               info,            &
                               ierror )

  use yac, dummy => yac_fput_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)            :: field_id
  integer, intent (in)            :: nbr_pointsets
  integer, intent (in)            :: collection_size
  type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, &
                                                collection_size)
  integer, intent (out)           :: info
  integer, intent (out)           :: ierror

  integer :: i, j
  type(yac_dble_ptr)              :: send_field_dble(nbr_pointsets, &
                                                     collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  call send_field_to_dble_ptr(field_id,        &
                              nbr_pointsets,   &
                              collection_size, &
                              send_field,      &
                              send_field_dble)

  call yac_fput_dble_ptr ( field_id,        &
                           nbr_pointsets,   &
                           collection_size, &
                           send_field_dble, &
                           info,            &
                           ierror )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
    end do
  end do

end subroutine yac_fput_real_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_real_ptr ( field_id,        &
                                    nbr_pointsets,   &
                                    collection_size, &
                                    send_field,      &
                                    send_frac_mask,  &
                                    info,            &
                                    ierror )

  use yac, dummy => yac_fput_frac_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)            :: field_id
  integer, intent (in)            :: nbr_pointsets
  integer, intent (in)            :: collection_size
  type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, &
                                                collection_size)
  type(yac_real_ptr), intent (in) :: send_frac_mask(nbr_pointsets, &
                                                    collection_size)
  integer, intent (out)           :: info
  integer, intent (out)           :: ierror

  integer :: i, j
  type(yac_dble_ptr)              :: send_field_dble(nbr_pointsets, &
                                                     collection_size)
  type(yac_dble_ptr)              :: send_frac_mask_dble(nbr_pointsets, &
                                                         collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  call send_field_to_dble_ptr(field_id,        &
                              nbr_pointsets,   &
                              collection_size, &
                              send_field,      &
                              send_field_dble, &
                              send_frac_mask,  &
                              send_frac_mask_dble)

  call yac_fput_frac_dble_ptr ( field_id,            &
                                nbr_pointsets,       &
                                collection_size,     &
                                send_field_dble,     &
                                send_frac_mask_dble, &
                                info,                &
                                ierror )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
      deallocate(send_frac_mask_dble(j, i)%p)
    end do
  end do

end subroutine yac_fput_frac_real_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_single_pointset_real ( field_id,        &
                                           nbr_hor_points,  &
                                           collection_size, &
                                           send_field,      &
                                           info,            &
                                           ierror )

  use yac, dummy => yac_fput_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  real, intent (in)             :: send_field(nbr_hor_points, &
                                              collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  double precision              :: send_field_dble(nbr_hor_points, &
                                                   collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call send_field_to_dble_single(field_id,        &
                                 nbr_hor_points,  &
                                 collection_size, &
                                 send_field,      &
                                 send_field_dble)

  call yac_fput_single_pointset_dble ( field_id,        &
                                       nbr_hor_points,  &
                                       collection_size, &
                                       send_field_dble, &
                                       info,            &
                                       ierror )

end subroutine yac_fput_single_pointset_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_single_pointset_real ( field_id,        &
                                                nbr_hor_points,  &
                                                collection_size, &
                                                send_field,      &
                                                send_frac_mask,  &
                                                info,            &
                                                ierror )

  use yac, dummy => yac_fput_frac_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  real, intent (in)             :: send_field(nbr_hor_points, &
                                              collection_size)
  real, intent (in)             :: send_frac_mask(nbr_hor_points, &
                                                  collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  double precision              :: send_field_dble(nbr_hor_points, &
                                                   collection_size)
  double precision              :: send_frac_mask_dble(nbr_hor_points, &
                                                       collection_size)
  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call send_field_to_dble_single(field_id,        &
                                 nbr_hor_points,  &
                                 collection_size, &
                                 send_field,      &
                                 send_field_dble, &
                                 send_frac_mask,  &
                                 send_frac_mask_dble)

  call yac_fput_frac_single_pointset_dble ( field_id,            &
                                            nbr_hor_points,      &
                                            collection_size,     &
                                            send_field_dble,     &
                                            send_frac_mask_dble, &
                                            info,                &
                                            ierror )

end subroutine yac_fput_frac_single_pointset_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_dble ( field_id,        &
                           nbr_hor_points,  &
                           nbr_pointsets,   &
                           collection_size, &
                           send_field,      &
                           info,            &
                           ierror )

  use yac, dummy => yac_fput_dble

  implicit none

  interface

     subroutine yac_cput__c ( field_id,         &
                              collection_size,  &
                              send_field,       &
                              info,             &
                              ierror ) bind ( c, name='yac_cput_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: nbr_pointsets
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              nbr_pointsets,  &
                                              collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )

  call yac_cput__c ( field_id,         &
                     collection_size,  &
                     send_field,       &
                     info,             &
                     ierror )

end subroutine yac_fput_dble

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_dble ( field_id,        &
                                nbr_hor_points,  &
                                nbr_pointsets,   &
                                collection_size, &
                                send_field,      &
                                send_frac_mask,  &
                                info,            &
                                ierror )

  use yac, dummy => yac_fput_frac_dble

  implicit none

  interface

     subroutine yac_cput_frac__c ( field_id,         &
                                   collection_size,  &
                                   send_field,       &
                                   send_frac_mask,   &
                                   info,             &
                                   ierror ) &
       bind ( c, name='yac_cput_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_frac__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: nbr_pointsets
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              nbr_pointsets,  &
                                              collection_size)
  double precision, intent (in) :: send_frac_mask(nbr_hor_points, &
                                                  nbr_pointsets,  &
                                                  collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )

  call yac_cput_frac__c ( field_id,         &
                          collection_size,  &
                          send_field,       &
                          send_frac_mask,   &
                          info,             &
                          ierror )

end subroutine yac_fput_frac_dble

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_dble_ptr ( field_id,        &
                               nbr_pointsets,   &
                               collection_size, &
                               send_field,      &
                               info,            &
                               ierror )

  use yac, dummy => yac_fput_dble_ptr
  use mo_yac_iso_c_helpers
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none

  interface

     subroutine yac_cput_ptr__c ( field_id,        &
                                  collection_size, &
                                  send_field,      &
                                  info,            &
                                  ierror )         &
         bind ( c, name='yac_cput_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_ptr__c

  end interface

  integer, intent (in)     :: field_id
  integer, intent (in)     :: nbr_pointsets
  integer, intent (in)     :: collection_size
  type(yac_dble_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
  integer, intent (out)    :: info
  integer, intent (out)    :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(nbr_pointsets, collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      send_field_(j, i) = &
        yac_dble2cptr("yac_fput_dble_ptr", "send_field", send_field(j, i))
    end do
  end do

  call yac_cput_ptr__c ( field_id,        &
                         collection_size, &
                         send_field_,     &
                         info,            &
                         ierror )

end subroutine yac_fput_dble_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_dble_ptr ( field_id,        &
                                    nbr_pointsets,   &
                                    collection_size, &
                                    send_field,      &
                                    send_frac_mask,  &
                                    info,            &
                                    ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fput_frac_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none

  interface

     subroutine yac_cput_frac_ptr__c ( field_id,        &
                                       collection_size, &
                                       send_field,      &
                                       send_frac_mask,  &
                                       info,            &
                                       ierror )         &
         bind ( c, name='yac_cput_frac_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: send_frac_mask(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_frac_ptr__c

  end interface

  integer, intent (in)     :: field_id
  integer, intent (in)     :: nbr_pointsets
  integer, intent (in)     :: collection_size
  type(yac_dble_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
  type(yac_dble_ptr), intent (in) :: send_frac_mask(nbr_pointsets, collection_size)
  integer, intent (out)    :: info
  integer, intent (out)    :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(nbr_pointsets, collection_size)
  type(c_ptr) :: send_frac_mask_(nbr_pointsets, collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  do i = 1, collection_size
     do j = 1, nbr_pointsets
        send_field_(j, i) = &
          yac_dble2cptr(    &
            "yac_fput_frac_dble_ptr", "send_field", send_field(j, i))
        send_frac_mask_(j, i) = &
          yac_dble2cptr(        &
            "yac_fput_frac_dble_ptr", "send_frac_mask", send_frac_mask(j, i))
    end do
  end do

  call yac_cput_frac_ptr__c ( field_id,        &
                              collection_size, &
                              send_field_,     &
                              send_frac_mask_, &
                              info,            &
                              ierror )

end subroutine yac_fput_frac_dble_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_single_pointset_dble ( field_id,        &
                                           nbr_hor_points,  &
                                           collection_size, &
                                           send_field,      &
                                           info,            &
                                           ierror )

  use yac, dummy => yac_fput_single_pointset_dble

  implicit none

  interface

     subroutine yac_cput__c ( field_id,        &
                              collection_size, &
                              send_field,      &
                              info,            &
                              ierror ) bind ( c, name='yac_cput_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call yac_cput__c ( field_id,        &
                     collection_size, &
                     send_field,      &
                     info,            &
                     ierror )

end subroutine yac_fput_single_pointset_dble

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_single_pointset_dble ( field_id,        &
                                                nbr_hor_points,  &
                                                collection_size, &
                                                send_field,      &
                                                send_frac_mask,  &
                                                info,            &
                                                ierror )

  use yac, dummy => yac_fput_frac_single_pointset_dble

  implicit none

  interface

     subroutine yac_cput_frac__c ( field_id,        &
                                   collection_size, &
                                   send_field,      &
                                   send_frac_mask,  &
                                   info,            &
                                   ierror )         &
       bind ( c, name='yac_cput_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_frac__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              collection_size)
  double precision, intent (in) :: send_frac_mask(nbr_hor_points, &
                                                  collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call yac_cput_frac__c ( field_id,        &
                          collection_size, &
                          send_field,      &
                          send_frac_mask,  &
                          info,            &
                          ierror )

end subroutine yac_fput_frac_single_pointset_dble

! ---------------------------------- get -------------------------------

subroutine yac_fget_real ( field_id,        &
                           nbr_hor_points,  &
                           collection_size, &
                           recv_field,      &
                           info,            &
                           ierror )

  use yac, dummy => yac_fget_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: nbr_hor_points  !< [IN] number of horizontal points
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  real, intent (inout)  :: recv_field(nbr_hor_points, collection_size)
                                           !< [INOUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error

  double precision :: recv_field_dble(nbr_hor_points, collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call recv_field_to_dble(field_id,        &
                          nbr_hor_points,  &
                          collection_size, &
                          recv_field,      &
                          recv_field_dble)

  call yac_fget_dble ( field_id,        &
                       nbr_hor_points,  &
                       collection_size, &
                       recv_field_dble, &
                       info,            &
                       ierror )

  call recv_field_from_dble(field_id,        &
                            nbr_hor_points,  &
                            collection_size, &
                            recv_field_dble, &
                            recv_field)

end subroutine yac_fget_real

subroutine yac_fget_real_ptr ( field_id,        &
                               collection_size, &
                               recv_field,      &
                               info,            &
                               ierror )

  use yac, dummy => yac_fget_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_real_ptr)    :: recv_field(collection_size) !< [INOUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  type(yac_dble_ptr)    :: recv_field_dble(collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size,  1, (/SIZE(recv_field(1)%p, 1)/) )

  call recv_field_to_dble_ptr(field_id,        &
                              collection_size, &
                              recv_field,      &
                              recv_field_dble)

  call yac_fget_dble_ptr ( field_id,        &
                           collection_size, &
                           recv_field_dble, &
                           info,            &
                           ierror )

  call recv_field_from_dble_ptr(field_id,        &
                                collection_size, &
                                recv_field_dble, &
                                recv_field)

end subroutine yac_fget_real_ptr

subroutine yac_fget_dble ( field_id,        &
                           nbr_hor_points,  &
                           collection_size, &
                           recv_field,      &
                           info,            &
                           ierror )

  use yac, dummy => yac_fget_dble

  implicit none

  interface

     subroutine yac_cget__c ( field_id,        &
                              collection_size, &
                              recv_field,      &
                              info,            &
                              ierror ) bind ( c, name='yac_cget_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       real    ( kind=c_double )     :: recv_field(*)       !< [OUT] returned field
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget__c

  end interface

  integer, intent (in)          :: field_id        !< [IN] field identifier
  integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
  integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
  double precision, intent (inout):: recv_field(nbr_hor_points, collection_size)
                                                  !< [INOUT] returned field
  integer, intent (out)         :: info           !< [OUT] returned info
  integer, intent (out)         :: ierror         !< [OUT] returned error

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call yac_cget__c ( field_id,        &
                     collection_size, &
                     recv_field,      &
                     info,            &
                     ierror )

end subroutine yac_fget_dble

subroutine yac_fget_dble_ptr ( field_id,        &
                               collection_size, &
                               recv_field,      &
                               info,            &
                               ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fget_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none

  interface

     subroutine yac_cget_c ( field_id,        &
                             collection_size, &
                             recv_field,      &
                             info,            &
                             ierror ) bind ( c, name='yac_cget' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       type(c_ptr)                   :: recv_field(*)       !< [OUT] returned field
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget_c

  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr)    :: recv_field(collection_size) !< [OUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  integer :: i
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    recv_field_(i) = yac_dble2cptr("yac_fget_dble_ptr", "recv_field", recv_field(i))
  end do

  call yac_cget_c ( field_id,        &
                    collection_size, &
                    recv_field_,     &
                    info,            &
                    ierror )

end subroutine yac_fget_dble_ptr

! ---------------------------------- get_raw -------------------------------

subroutine yac_fget_raw_real ( field_id,              &
                               src_field_buffer_size, &
                               collection_size,       &
                               src_field_buffer,      &
                               info,                  &
                               ierror )

  use yac, dummy => yac_fget_raw_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: src_field_buffer_size
                                           !< [IN] source field buffer size
                                           !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  real, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                           !< [OUT] returned source field buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  double precision :: &
    src_field_buffer_dble(src_field_buffer_size, collection_size)

  call yac_fget_raw_dble ( field_id,              &
                           src_field_buffer_size, &
                           collection_size,       &
                           src_field_buffer_dble, &
                           info,                  &
                           ierror )

  src_field_buffer = real(src_field_buffer_dble)

end subroutine yac_fget_raw_real

subroutine yac_fget_raw_real_ptr ( field_id,         &
                                   num_src_fields,   &
                                   collection_size,  &
                                   src_field_buffer, &
                                   info,             &
                                   ierror )

  use yac, dummy => yac_fget_raw_real_ptr
  use, intrinsic :: iso_c_binding, only: c_null_char

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  integer :: i, j
  integer :: src_field_buffer_sizes(num_src_fields)
  type(yac_dble_ptr) :: src_field_buffer_dble(num_src_fields, collection_size)

  do i = 1, num_src_fields
    src_field_buffer_sizes(i) = size(src_field_buffer(i,1)%p(:))
    do j = 1, collection_size
      YAC_FASSERT(src_field_buffer_sizes(i) == size(src_field_buffer(i,j)%p), "ERROR(yac_fget_raw_real_ptr): inconsistent source buffer sizes")
      allocate(src_field_buffer_dble(i,j)%p(src_field_buffer_sizes(i)))
    end do
  end do

  call yac_fcheck_src_field_buffer_sizes( &
    field_id, num_src_fields, collection_size, src_field_buffer_sizes)

  call yac_fget_raw_dble_ptr ( field_id,              &
                               num_src_fields,        &
                               collection_size,       &
                               src_field_buffer_dble, &
                               info,                  &
                               ierror )

  do i = 1, num_src_fields
    do j = 1, collection_size
      src_field_buffer(i,j)%p = &
        real(src_field_buffer_dble(i,j)%p)
      deallocate(src_field_buffer_dble(i,j)%p)
    end do
  end do

end subroutine yac_fget_raw_real_ptr

subroutine yac_fget_raw_dble ( field_id,              &
                               src_field_buffer_size, &
                               collection_size,       &
                               src_field_buffer,      &
                               info,                  &
                               ierror )

  use yac, dummy => yac_fget_raw_dble

  implicit none

  interface

     subroutine yac_cget_raw__c ( field_id,         &
                                  collection_size,  &
                                  src_field_buffer, &
                                  info,             &
                                  ierror )          &
       bind ( c, name='yac_cget_raw_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       real    ( kind=c_double )     :: src_field_buffer(*) !< [OUT] returned source field buffer
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget_raw__c

  end interface

  integer, intent (in)  :: field_id              !< [IN] field identifier
  integer, intent (in)  :: src_field_buffer_size !< [IN] source buffer size
                                                 !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  double precision, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] returned source field buffer
  integer, intent (out) :: info                  !< [OUT] returned info
  integer, intent (out) :: ierror                !< [OUT] returned error handler

  call yac_fcheck_src_field_buffer_size( &
    field_id, collection_size, src_field_buffer_size)

  call yac_cget_raw__c ( field_id,         &
                         collection_size,  &
                         src_field_buffer, &
                         info,             &
                         ierror )

end subroutine yac_fget_raw_dble

subroutine yac_get_raw_dble_ptr ( field_id,         &
                                  num_src_fields,   &
                                  collection_size,  &
                                  src_field_buffer, &
                                  info,             &
                                  ierror,           &
                                  use_async )

  use mo_yac_iso_c_helpers
  use yac
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

  implicit none

  interface

     subroutine yac_cget_raw_ptr__c ( field_id,         &
                                      collection_size,  &
                                      src_field_buffer, &
                                      info,             &
                                      ierror )          &
        bind ( c, name='yac_cget_raw_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       type(c_ptr)                   :: src_field_buffer(*) !< [OUT] source field buffer
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget_raw_ptr__c

     subroutine yac_cget_raw_async_ptr__c ( field_id,         &
                                            collection_size,  &
                                            src_field_buffer, &
                                            info,             &
                                            ierror )          &
        bind ( c, name='yac_cget_raw_async_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       type(c_ptr)                   :: src_field_buffer(*) !< [OUT] source field buffer
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget_raw_async_ptr__c

  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number of source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler
  logical, intent(in) :: use_async         !< [IN] use asynchronous get

  integer :: i, j
  integer :: src_field_buffer_sizes(num_src_fields)
  type(c_ptr) :: src_field_buffer_(num_src_fields, collection_size)

  src_field_buffer_sizes = &
    (/(SIZE(src_field_buffer(i, 1)%p),i=1,num_src_fields)/)

  call yac_fcheck_src_field_buffer_sizes( &
    field_id, num_src_fields, collection_size, src_field_buffer_sizes)

  do i = 1, collection_size
    do j = 1, num_src_fields
      YAC_FASSERT(src_field_buffer_sizes(j) == size(src_field_buffer(j,i)%p), "ERROR(yac_get_raw_dble_ptr): inconsistent source buffer sizes")
      src_field_buffer_(j, i) = &
        yac_dble2cptr( &
          "yac_get_raw_dble_ptr", "src_field_buffer", src_field_buffer(j, i))
    end do
  end do

  if (use_async) then
    call yac_cget_raw_async_ptr__c ( field_id,          &
                                     collection_size,   &
                                     src_field_buffer_, &
                                     info,              &
                                     ierror )
  else
    call yac_cget_raw_ptr__c ( field_id,          &
                               collection_size,   &
                               src_field_buffer_, &
                               info,              &
                               ierror )
  end if

end subroutine yac_get_raw_dble_ptr

subroutine yac_fget_raw_dble_ptr ( field_id,         &
                                   num_src_fields,   &
                                   collection_size,  &
                                   src_field_buffer, &
                                   info,             &
                                   ierror )

  use yac, dummy => yac_fget_raw_dble_ptr

  implicit none

  interface
    subroutine yac_get_raw_dble_ptr ( field_id,             &
                                      num_src_fields,       &
                                      collection_size,      &
                                      src_field_buffer,     &
                                      info,                 &
                                      ierror,               &
                                      use_async )
      import :: yac_dble_ptr
      integer, intent (in)  :: field_id
      integer, intent (in)  :: num_src_fields
      integer, intent (in)  :: collection_size
      type(yac_dble_ptr), intent(inout) :: &
        src_field_buffer(num_src_fields, collection_size)
      integer, intent (out) :: info
      integer, intent (out) :: ierror
      logical, intent(in) :: use_async
    end subroutine yac_get_raw_dble_ptr
  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number of source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  call yac_get_raw_dble_ptr( &
    field_id, num_src_fields, collection_size, &
    src_field_buffer, info, ierror, .false. )

end subroutine yac_fget_raw_dble_ptr

subroutine yac_fget_raw_frac_real ( field_id,              &
                                    src_field_buffer_size, &
                                    collection_size,       &
                                    src_field_buffer,      &
                                    src_frac_mask_buffer,  &
                                    info,                  &
                                    ierror )

  use yac, dummy => yac_fget_raw_frac_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id              !< [IN] field identifier
  integer, intent (in)  :: src_field_buffer_size !< [IN] source buffer size
                                                 !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  real, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] returned source field buffer
  real, intent (out)    :: &
    src_frac_mask_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] returned source fractional mask buffer
  integer, intent (out) :: info                  !< [OUT] returned info
  integer, intent (out) :: ierror                !< [OUT] returned error handler

  double precision :: &
    src_field_buffer_dble(src_field_buffer_size, collection_size)
  double precision :: &
    src_frac_mask_buffer_dble(src_field_buffer_size, collection_size)

  call yac_fget_raw_frac_dble ( field_id,                  &
                                src_field_buffer_size,     &
                                collection_size,           &
                                src_field_buffer_dble,     &
                                src_frac_mask_buffer_dble, &
                                info,                      &
                                ierror )

  src_field_buffer = real(src_field_buffer_dble)
  src_frac_mask_buffer = real(src_frac_mask_buffer_dble)

end subroutine yac_fget_raw_frac_real

subroutine yac_fget_raw_frac_real_ptr ( field_id,             &
                                        num_src_fields,       &
                                        collection_size,      &
                                        src_field_buffer,     &
                                        src_frac_mask_buffer, &
                                        info,                 &
                                        ierror )

  use yac, dummy => yac_fget_raw_frac_real_ptr
  use, intrinsic :: iso_c_binding, only: c_null_char

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  type(yac_real_ptr), intent(inout) :: &
    src_frac_mask_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source fractional mask buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  integer :: i, j
  integer :: src_field_buffer_sizes(num_src_fields)
  type(yac_dble_ptr) :: src_field_buffer_dble(num_src_fields, collection_size)
  type(yac_dble_ptr) :: src_frac_mask_buffer_dble(num_src_fields, collection_size)

  do i = 1, num_src_fields
    src_field_buffer_sizes(i) = size(src_field_buffer(i,1)%p(:))
    do j = 1, collection_size
      YAC_FASSERT(src_field_buffer_sizes(i) == size(src_field_buffer(i,j)%p), "ERROR(yac_fget_raw_frac_real_ptr): inconsistent source buffer sizes")
      allocate(src_field_buffer_dble(i,j)%p(src_field_buffer_sizes(i)))
      allocate(src_frac_mask_buffer_dble(i,j)%p(src_field_buffer_sizes(i)))
    end do
  end do

  call yac_fcheck_src_field_buffer_sizes( &
    field_id, num_src_fields, collection_size, src_field_buffer_sizes)

  call yac_fget_raw_frac_dble_ptr ( field_id,                  &
                                    num_src_fields,            &
                                    collection_size,           &
                                    src_field_buffer_dble,     &
                                    src_frac_mask_buffer_dble, &
                                    info,                      &
                                    ierror )

  do i = 1, num_src_fields
    do j = 1, collection_size
      src_field_buffer(i,j)%p = &
        real(src_field_buffer_dble(i,j)%p)
      deallocate(src_field_buffer_dble(i,j)%p)
      src_frac_mask_buffer(i,j)%p = &
        real(src_frac_mask_buffer_dble(i,j)%p)
      deallocate(src_frac_mask_buffer_dble(i,j)%p)
    end do
  end do

end subroutine yac_fget_raw_frac_real_ptr

subroutine yac_fget_raw_frac_dble ( field_id,              &
                                    src_field_buffer_size, &
                                    collection_size,       &
                                    src_field_buffer,      &
                                    src_frac_mask_buffer,  &
                                    info,                  &
                                    ierror )

  use yac, dummy => yac_fget_raw_frac_dble

  implicit none

  interface

     subroutine yac_cget_raw_frac__c ( field_id,             &
                                       collection_size,      &
                                       src_field_buffer,     &
                                       src_frac_mask_buffer, &
                                       info,                 &
                                       ierror )              &
       bind ( c, name='yac_cget_raw_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id                !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size         !< [IN] collection size
       real    ( kind=c_double )     :: src_field_buffer(*)     !< [OUT] returned source field buffer
       real    ( kind=c_double )     :: src_frac_mask_buffer(*) !< [OUT] returned source fractional mask buffer
       integer ( kind=c_int )        :: info                    !< [OUT] returned info
       integer ( kind=c_int )        :: ierror                  !< [OUT] returned error handler

     end subroutine yac_cget_raw_frac__c

  end interface

  integer, intent (in)  :: field_id              !< [IN] field identifier
  integer, intent (in)  :: src_field_buffer_size !< [IN] source buffer size
                                                !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  double precision, intent (out) :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                !< [OUT] returned source field buffer
  double precision, intent (out) :: &
    src_frac_mask_buffer(src_field_buffer_size, collection_size)
                                                !< [OUT] returned source fractional mask buffer
  integer, intent (out) :: info                  !< [OUT] returned info
  integer, intent (out) :: ierror                !< [OUT] returned error handler

  call yac_fcheck_src_field_buffer_size( &
    field_id, collection_size, src_field_buffer_size)

  call yac_cget_raw_frac__c ( field_id,             &
                              collection_size,      &
                              src_field_buffer,     &
                              src_frac_mask_buffer, &
                              info,                 &
                              ierror )

end subroutine yac_fget_raw_frac_dble

subroutine yac_get_raw_frac_dble_ptr ( field_id,             &
                                       num_src_fields,       &
                                       collection_size,      &
                                       src_field_buffer,     &
                                       src_frac_mask_buffer, &
                                       info,                 &
                                       ierror,               &
                                       use_async )

  use mo_yac_iso_c_helpers
  use yac
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

  implicit none

  interface

     subroutine yac_cget_raw_frac_ptr__c ( field_id,             &
                                           collection_size,      &
                                           src_field_buffer,     &
                                           src_frac_mask_buffer, &
                                           info,                 &
                                           ierror )              &
        bind ( c, name='yac_cget_raw_frac_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id                !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size         !< [IN] collection size
       type(c_ptr)                   :: src_field_buffer(*)     !< [OUT] source field buffer
       type(c_ptr)                   :: src_frac_mask_buffer(*) !< [OUT] source fractional mask buffer
       integer ( kind=c_int )        :: info                    !< [OUT] returned info
       integer ( kind=c_int )        :: ierror                  !< [OUT] returned error handler

     end subroutine yac_cget_raw_frac_ptr__c

     subroutine yac_cget_raw_frac_async_ptr__c ( field_id,             &
                                                 collection_size,      &
                                                 src_field_buffer,     &
                                                 src_frac_mask_buffer, &
                                                 info,                 &
                                                 ierror )              &
        bind ( c, name='yac_cget_raw_frac_async_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id                !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size         !< [IN] collection size
       type(c_ptr)                   :: src_field_buffer(*)     !< [OUT] source field buffer
       type(c_ptr)                   :: src_frac_mask_buffer(*) !< [OUT] source fractional mask buffer
       integer ( kind=c_int )        :: info                    !< [OUT] returned info
       integer ( kind=c_int )        :: ierror                  !< [OUT] returned error handler

     end subroutine yac_cget_raw_frac_async_ptr__c

  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number of source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  type(yac_dble_ptr), intent(inout) :: &
    src_frac_mask_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source fractional mask buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler
  logical, intent(in) :: use_async         !< [IN] use asynchronous get

  integer :: i, j
  integer :: src_field_buffer_sizes(num_src_fields)
  type(c_ptr) :: src_field_buffer_(num_src_fields, collection_size)
  type(c_ptr) :: src_frac_mask_buffer_(num_src_fields, collection_size)

  src_field_buffer_sizes = &
    (/(SIZE(src_field_buffer(i, 1)%p),i=1,num_src_fields)/)

  call yac_fcheck_src_field_buffer_sizes( &
    field_id, num_src_fields, collection_size, src_field_buffer_sizes)

  do i = 1, collection_size
    do j = 1, num_src_fields
      YAC_FASSERT(src_field_buffer_sizes(j) == size(src_field_buffer(j,i)%p), "ERROR(yac_get_raw_frac_dble_ptr): inconsistent source buffer sizes")
      YAC_FASSERT(src_field_buffer_sizes(j) == size(src_frac_mask_buffer(j,i)%p), "ERROR(yac_get_raw_frac_dble_ptr): inconsistent source buffer sizes")
      src_field_buffer_(j, i) = &
        yac_dble2cptr( &
          "yac_get_raw_frac_dble_ptr", "src_field_buffer", src_field_buffer(j, i))
      src_frac_mask_buffer_(j, i) = &
        yac_dble2cptr( &
          "yac_get_raw_frac_dble_ptr", "src_frac_mask_buffer", src_frac_mask_buffer(j, i))
    end do
  end do

  if (use_async) then
    call yac_cget_raw_frac_async_ptr__c ( field_id,              &
                                          collection_size,       &
                                          src_field_buffer_,     &
                                          src_frac_mask_buffer_, &
                                          info,                  &
                                          ierror )
  else
    call yac_cget_raw_frac_ptr__c ( field_id,              &
                                    collection_size,       &
                                    src_field_buffer_,     &
                                    src_frac_mask_buffer_, &
                                    info,                  &
                                    ierror )
  end if

end subroutine yac_get_raw_frac_dble_ptr

subroutine yac_fget_raw_frac_dble_ptr ( field_id,             &
                                        num_src_fields,       &
                                        collection_size,      &
                                        src_field_buffer,     &
                                        src_frac_mask_buffer, &
                                        info,                 &
                                        ierror )

  use yac, dummy => yac_fget_raw_frac_dble_ptr

  implicit none

  interface
    subroutine yac_get_raw_frac_dble_ptr ( field_id,             &
                                           num_src_fields,       &
                                           collection_size,      &
                                           src_field_buffer,     &
                                           src_frac_mask_buffer, &
                                           info,                 &
                                           ierror,               &
                                           use_async )
      import :: yac_dble_ptr
      integer, intent (in)  :: field_id
      integer, intent (in)  :: num_src_fields
      integer, intent (in)  :: collection_size
      type(yac_dble_ptr), intent(inout) :: &
        src_field_buffer(num_src_fields, collection_size)
      type(yac_dble_ptr), intent(inout) :: &
        src_frac_mask_buffer(num_src_fields, collection_size)
      integer, intent (out) :: info
      integer, intent (out) :: ierror
      logical, intent(in) :: use_async
    end subroutine yac_get_raw_frac_dble_ptr
  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number of source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  type(yac_dble_ptr), intent(inout) :: &
    src_frac_mask_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source fractional mask buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  call yac_get_raw_frac_dble_ptr( &
    field_id, num_src_fields, collection_size, &
    src_field_buffer, src_frac_mask_buffer, &
    info, ierror, .false.)

end subroutine yac_fget_raw_frac_dble_ptr

! ---------------------------------- get_async -------------------------------

subroutine yac_fget_async_dble_ptr ( field_id,        &
                                     collection_size, &
                                     recv_field,      &
                                     info,            &
                                     ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fget_async_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none

  interface

     subroutine yac_cget_async_c ( field_id,        &
                                   collection_size, &
                                   recv_field,      &
                                   info,            &
                                   ierror ) bind ( c, name='yac_cget_async' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       type(c_ptr)                   :: recv_field(*)       !< [OUT] returned field
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget_async_c

  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr)    :: recv_field(collection_size) !< [OUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  integer :: i
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    recv_field_(i) = yac_dble2cptr("yac_fget_async_dble_ptr", "recv_field", recv_field(i))
  end do

  call yac_cget_async_c ( field_id,        &
                          collection_size, &
                          recv_field_,     &
                          info,            &
                          ierror )

end subroutine yac_fget_async_dble_ptr

! ---------------------------------- get_raw_async -------------------------------

subroutine yac_fget_raw_async_dble_ptr ( field_id,             &
                                         num_src_fields,       &
                                         collection_size,      &
                                         src_field_buffer,     &
                                         info,                 &
                                         ierror )

  use yac, dummy => yac_fget_raw_async_dble_ptr

  implicit none

  interface
    subroutine yac_get_raw_dble_ptr ( field_id,             &
                                      num_src_fields,       &
                                      collection_size,      &
                                      src_field_buffer,     &
                                      info,                 &
                                      ierror,               &
                                      use_async )
      import :: yac_dble_ptr
      integer, intent (in)  :: field_id
      integer, intent (in)  :: num_src_fields
      integer, intent (in)  :: collection_size
      type(yac_dble_ptr), intent(inout) :: &
        src_field_buffer(num_src_fields, collection_size)
      integer, intent (out) :: info
      integer, intent (out) :: ierror
      logical, intent(in) :: use_async
    end subroutine yac_get_raw_dble_ptr
  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number of source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  call yac_get_raw_dble_ptr( &
    field_id, num_src_fields, collection_size, src_field_buffer, &
    info, ierror, .true.)

end subroutine yac_fget_raw_async_dble_ptr

subroutine yac_fget_raw_frac_async_dble_ptr ( field_id,             &
                                              num_src_fields,       &
                                              collection_size,      &
                                              src_field_buffer,     &
                                              src_frac_mask_buffer, &
                                              info,                 &
                                              ierror )

  use yac, dummy => yac_fget_raw_frac_async_dble_ptr

  implicit none

  interface
    subroutine yac_get_raw_frac_dble_ptr ( field_id,             &
                                           num_src_fields,       &
                                           collection_size,      &
                                           src_field_buffer,     &
                                           src_frac_mask_buffer, &
                                           info,                 &
                                           ierror,               &
                                           use_async )
      import :: yac_dble_ptr
      integer, intent (in)  :: field_id
      integer, intent (in)  :: num_src_fields
      integer, intent (in)  :: collection_size
      type(yac_dble_ptr), intent(inout) :: &
        src_field_buffer(num_src_fields, collection_size)
      type(yac_dble_ptr), intent(inout) :: &
        src_frac_mask_buffer(num_src_fields, collection_size)
      integer, intent (out) :: info
      integer, intent (out) :: ierror
      logical, intent(in) :: use_async
    end subroutine yac_get_raw_frac_dble_ptr
  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: num_src_fields  !< [IN] number of source fields
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source field buffer
  type(yac_dble_ptr), intent(inout) :: &
    src_frac_mask_buffer(num_src_fields, collection_size)
                                           !< [INOUT] returned source fractional mask buffer
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  call yac_get_raw_frac_dble_ptr( &
    field_id, num_src_fields, collection_size, &
    src_field_buffer, src_frac_mask_buffer, &
    info, ierror, .true.)

end subroutine yac_fget_raw_frac_async_dble_ptr

! ---------------------------------- exchange --------------------------------

subroutine yac_fexchange_real ( send_field_id,       &
                                recv_field_id,       &
                                send_nbr_hor_points, &
                                send_nbr_pointsets,  &
                                recv_nbr_hor_points, &
                                collection_size,     &
                                send_field,          &
                                recv_field,          &
                                send_info,           &
                                recv_info,           &
                                ierror )

  use yac, dummy => yac_fexchange_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] field identifier
  integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      send_nbr_pointsets,  &
                                      collection_size)
                                               !< [IN] send field
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, &
                                  send_nbr_pointsets,  &
                                  collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, &
                                  collection_size)

  integer :: i

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_field,          &
                          send_buffer)
  call recv_field_to_dble(recv_field_id,        &
                          recv_nbr_hor_points,  &
                          collection_size,      &
                          recv_field,           &
                          recv_buffer)

  call yac_fexchange_dble ( send_field_id,       &
                            recv_field_id,       &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            send_buffer,         &
                            recv_buffer,         &
                            send_info,           &
                            recv_info,           &
                            ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_real

subroutine yac_fexchange_frac_real ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_hor_points, &
                                     send_nbr_pointsets,  &
                                     recv_nbr_hor_points, &
                                     collection_size,     &
                                     send_field,          &
                                     send_frac_mask,      &
                                     recv_field,          &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

  use yac, dummy => yac_fexchange_frac_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] field identifier
  integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      send_nbr_pointsets,  &
                                      collection_size)
                                               !< [IN] send field
  real, intent (in)     :: send_frac_mask(send_nbr_hor_points, &
                                          send_nbr_pointsets,  &
                                          collection_size)
                                               !< [IN] fractional mask
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, &
                                  send_nbr_pointsets,  &
                                  collection_size)
  double precision :: send_frac_mask_buffer(send_nbr_hor_points, &
                                            send_nbr_pointsets,  &
                                            collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, &
                                  collection_size)

  integer :: i

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_field,          &
                          send_buffer,         &
                          send_frac_mask,      &
                          send_frac_mask_buffer)
  call recv_field_to_dble(recv_field_id,        &
                          recv_nbr_hor_points,  &
                          collection_size,      &
                          recv_field,           &
                          recv_buffer)

  call yac_fexchange_frac_dble ( send_field_id,         &
                                 recv_field_id,         &
                                 send_nbr_hor_points,   &
                                 send_nbr_pointsets,    &
                                 recv_nbr_hor_points,   &
                                 collection_size,       &
                                 send_buffer,           &
                                 send_frac_mask_buffer, &
                                 recv_buffer,           &
                                 send_info,             &
                                 recv_info,             &
                                 ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_frac_real

subroutine yac_fexchange_raw_real ( send_field_id,         &
                                    recv_field_id,         &
                                    send_nbr_hor_points,   &
                                    send_nbr_pointsets,    &
                                    src_field_buffer_size, &
                                    collection_size,       &
                                    send_field,            &
                                    src_field_buffer,      &
                                    send_info,             &
                                    recv_info,             &
                                    ierror )

  use yac, dummy => yac_fexchange_raw_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id         !< [IN] field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets    !< [IN] number of send point sets
  integer, intent (in)  :: src_field_buffer_size !< [IN] source field buffer size
                                                 !!      (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      send_nbr_pointsets,  &
                                      collection_size)
                                                 !< [IN] send field
  real, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] source field buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, &
                                  send_nbr_pointsets,  &
                                  collection_size)
  double precision :: &
    src_field_buffer_dble(src_field_buffer_size, collection_size)

  integer :: i

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_field,          &
                          send_buffer)

  call yac_fexchange_raw_dble ( send_field_id,         &
                                recv_field_id,         &
                                send_nbr_hor_points,   &
                                send_nbr_pointsets,    &
                                src_field_buffer_size, &
                                collection_size,       &
                                send_buffer,           &
                                src_field_buffer_dble, &
                                send_info,             &
                                recv_info,             &
                                ierror )

  src_field_buffer = real(src_field_buffer_dble)

end subroutine yac_fexchange_raw_real

subroutine yac_fexchange_raw_frac_real ( send_field_id,         &
                                          recv_field_id,         &
                                          send_nbr_hor_points,   &
                                          send_nbr_pointsets,    &
                                          src_field_buffer_size, &
                                          collection_size,       &
                                          send_field,            &
                                          send_frac_mask,        &
                                          src_field_buffer,      &
                                          src_frac_mask_buffer,  &
                                          send_info,             &
                                          recv_info,             &
                                          ierror )

  use yac, dummy => yac_fexchange_raw_frac_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id         !< [IN] field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets    !< [IN] number of send point sets
  integer, intent (in)  :: src_field_buffer_size !< [IN] source buffer size
                                                  !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  real, intent (in)     :: &
    send_field(send_nbr_hor_points, send_nbr_pointsets, collection_size)
                                                  !< [IN] send field
  real, intent (in)     :: &
    send_frac_mask(send_nbr_hor_points, send_nbr_pointsets, collection_size)
                                                  !< [IN] fractional mask
  real, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                  !< [OUT] returned source field buffer
  real, intent (out)    :: &
    src_frac_mask_buffer(src_field_buffer_size, collection_size)
                                                  !< [OUT] returned source fractional mask buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  double precision :: send_field_dble(send_nbr_hor_points, &
                                      send_nbr_pointsets,  &
                                      collection_size)
  double precision :: send_frac_mask_dble(send_nbr_hor_points, &
                                          send_nbr_pointsets,  &
                                          collection_size)
  double precision :: &
    src_field_buffer_dble(src_field_buffer_size, collection_size)
  double precision :: &
    src_frac_mask_buffer_dble(src_field_buffer_size, collection_size)

  integer :: i

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_field,          &
                          send_field_dble)
  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_frac_mask,      &
                          send_frac_mask_dble)

  call yac_fexchange_raw_frac_dble ( send_field_id,             &
                                     recv_field_id,             &
                                     send_nbr_hor_points,       &
                                     send_nbr_pointsets,        &
                                     src_field_buffer_size,     &
                                     collection_size,           &
                                     send_field_dble,           &
                                     send_frac_mask_dble,       &
                                     src_field_buffer_dble,     &
                                     src_frac_mask_buffer_dble, &
                                     send_info,                 &
                                     recv_info,                 &
                                     ierror )

  src_field_buffer = real(src_field_buffer_dble)
  src_frac_mask_buffer = real(src_frac_mask_buffer_dble)

end subroutine yac_fexchange_raw_frac_real

subroutine yac_fexchange_real_ptr ( send_field_id,      &
                                    recv_field_id,      &
                                    send_nbr_pointsets, &
                                    collection_size,    &
                                    send_field,         &
                                    recv_field,         &
                                    send_info,          &
                                    recv_info,          &
                                    ierror )

  use yac, dummy => yac_fexchange_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
                                               !< [IN] send field
  type(yac_real_ptr)    :: recv_field(collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  integer :: i, j
  type(yac_dble_ptr) :: send_field_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: recv_field_dble(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  call send_field_to_dble_ptr(send_field_id,      &
                              send_nbr_pointsets, &
                              collection_size,    &
                              send_field,         &
                              send_field_dble)
  call recv_field_to_dble_ptr(recv_field_id,   &
                              collection_size, &
                              recv_field,      &
                              recv_field_dble)

  call yac_fexchange_dble_ptr ( send_field_id,      &
                                recv_field_id,      &
                                send_nbr_pointsets, &
                                collection_size,    &
                                send_field_dble,    &
                                recv_field_dble,    &
                                send_info,          &
                                recv_info,          &
                                ierror )

  call recv_field_from_dble_ptr(recv_field_id,   &
                                collection_size, &
                                recv_field_dble, &
                                recv_field)
  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
    end do
  end do

end subroutine yac_fexchange_real_ptr

subroutine yac_fexchange_frac_real_ptr ( send_field_id,      &
                                         recv_field_id,      &
                                         send_nbr_pointsets, &
                                         collection_size,    &
                                         send_field,         &
                                         send_frac_mask,     &
                                         recv_field,         &
                                         send_info,          &
                                         recv_info,          &
                                         ierror )

  use yac, dummy => yac_fexchange_frac_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
                                               !< [IN] send field
  type(yac_real_ptr), intent (in) ::                      &
                           send_frac_mask(send_nbr_pointsets, &
                                          collection_size)
                                               !< [IN] fractional mask
  type(yac_real_ptr)    :: recv_field(collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  integer :: i, j
  type(yac_dble_ptr) :: send_field_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: send_frac_mask_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: recv_field_dble(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  call send_field_to_dble_ptr(send_field_id,      &
                              send_nbr_pointsets, &
                              collection_size,    &
                              send_field,         &
                              send_field_dble,    &
                              send_frac_mask,     &
                              send_frac_mask_dble)
  call recv_field_to_dble_ptr(recv_field_id,   &
                              collection_size, &
                              recv_field,      &
                              recv_field_dble)

  call yac_fexchange_frac_dble_ptr ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_pointsets,  &
                                     collection_size,     &
                                     send_field_dble,     &
                                     send_frac_mask_dble, &
                                     recv_field_dble,     &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

  call recv_field_from_dble_ptr(recv_field_id,   &
                                collection_size, &
                                recv_field_dble, &
                                recv_field)
  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
      deallocate(send_frac_mask_dble(j, i)%p)
    end do
  end do

end subroutine yac_fexchange_frac_real_ptr

subroutine yac_fexchange_raw_real_ptr ( send_field_id,      &
                                        recv_field_id,      &
                                        send_nbr_pointsets, &
                                        collection_size,    &
                                        send_field,         &
                                        src_field_buffer,   &
                                        send_info,          &
                                        recv_info,          &
                                        ierror )

  use yac, dummy => yac_fexchange_raw_real_ptr
  use mo_yac_real_to_dble_utils
  use, intrinsic :: iso_c_binding, only: c_null_char

  implicit none

  integer, intent (in)               :: send_field_id      !< [IN] send field identifier
  integer, intent (in)               :: recv_field_id      !< [IN] recv field identifier
  integer, intent (in)               :: send_nbr_pointsets !< [IN] number of send point sets
  integer, intent (in)               :: collection_size    !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent (in)    :: &
    send_field(send_nbr_pointsets, collection_size)        !< [IN] send field
  type(yac_real_ptr), intent (inout) :: &
    src_field_buffer(send_nbr_pointsets, collection_size)  !< [INOUT] returned source field buffer
  integer, intent (out)              :: send_info          !< [OUT] returned send info
  integer, intent (out)              :: recv_info          !< [OUT] returned recv info
  integer, intent (out)              :: ierror             !< [OUT] returned error

  integer :: i, j
  integer :: src_field_buffer_sizes(send_nbr_pointsets)
  type(yac_dble_ptr) :: send_field_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: src_field_buffer_dble(send_nbr_pointsets, collection_size)

  src_field_buffer_sizes = &
    (/(SIZE(src_field_buffer(i, 1)%p),i=1,send_nbr_pointsets)/)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_sizes( &
    recv_field_id, send_nbr_pointsets, collection_size, src_field_buffer_sizes)

  call send_field_to_dble_ptr(send_field_id,      &
                              send_nbr_pointsets, &
                              collection_size,    &
                              send_field,         &
                              send_field_dble)

  do i = 1, send_nbr_pointsets
    do j = 1, collection_size
      YAC_FASSERT(src_field_buffer_sizes(i) == size(src_field_buffer(i,j)%p), "ERROR(yac_fexchange_raw_real_ptr): inconsistent source buffer sizes")
      allocate(src_field_buffer_dble(i,j)%p(src_field_buffer_sizes(i)))
    end do
  end do

  call yac_fcheck_src_field_buffer_sizes( &
    recv_field_id, send_nbr_pointsets, collection_size, src_field_buffer_sizes)

  call yac_fexchange_raw_dble_ptr ( send_field_id,         &
                                    recv_field_id,         &
                                    send_nbr_pointsets,    &
                                    collection_size,       &
                                    send_field_dble,       &
                                    src_field_buffer_dble, &
                                    send_info,             &
                                    recv_info,             &
                                    ierror )

  do i = 1, send_nbr_pointsets
    do j = 1, collection_size
      src_field_buffer(i,j)%p = &
        real(src_field_buffer_dble(i,j)%p)
      deallocate(src_field_buffer_dble(i,j)%p)
    end do
  end do

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
    end do
  end do

end subroutine yac_fexchange_raw_real_ptr

subroutine yac_fexchange_raw_frac_real_ptr ( send_field_id,        &
                                             recv_field_id,        &
                                             send_nbr_pointsets,   &
                                             collection_size,      &
                                             send_field,           &
                                             send_frac_mask,       &
                                             src_field_buffer,     &
                                             src_frac_mask_buffer, &
                                             send_info,            &
                                             recv_info,            &
                                             ierror )

  use yac, dummy => yac_fexchange_raw_frac_real_ptr
  use mo_yac_real_to_dble_utils
  use, intrinsic :: iso_c_binding, only: c_null_char

  implicit none

  integer, intent (in)               :: send_field_id       !< [IN] send field identifier
  integer, intent (in)               :: recv_field_id       !< [IN] recv field identifier
  integer, intent (in)               :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)               :: collection_size     !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent (in)    :: &
    send_field(send_nbr_pointsets, collection_size)         !< [IN] send field
  type(yac_real_ptr), intent (in)    :: &
    send_frac_mask(send_nbr_pointsets, collection_size)     !< [IN] fractional mask
  type(yac_real_ptr), intent(inout) :: &
    src_field_buffer(send_nbr_pointsets, collection_size)   !< [INOUT] returned source field buffer
  type(yac_real_ptr), intent(inout) :: &
    src_frac_mask_buffer(send_nbr_pointsets, collection_size)
                                                            !< [INOUT] returned source fractional mask buffer
  integer, intent (out)              :: send_info           !< [OUT] returned send info
  integer, intent (out)              :: recv_info           !< [OUT] returned recv info
  integer, intent (out)              :: ierror              !< [OUT] returned error

  integer :: i, j
  integer :: src_field_buffer_sizes(send_nbr_pointsets)
  type(yac_dble_ptr) :: send_field_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: send_frac_mask_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: src_field_buffer_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: src_frac_mask_buffer_dble(send_nbr_pointsets, collection_size)

  src_field_buffer_sizes = &
    (/(SIZE(src_field_buffer(i, 1)%p),i=1,send_nbr_pointsets)/)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_sizes( &
    recv_field_id, send_nbr_pointsets, collection_size, src_field_buffer_sizes)

  call send_field_to_dble_ptr(send_field_id,      &
                              send_nbr_pointsets, &
                              collection_size,    &
                              send_field,         &
                              send_field_dble,    &
                              send_frac_mask,     &
                              send_frac_mask_dble)

  do i = 1, send_nbr_pointsets
    do j = 1, collection_size
      YAC_FASSERT(src_field_buffer_sizes(i) == size(src_field_buffer(i,j)%p), "ERROR(yac_fexchange_raw_frac_real_ptr): inconsistent source buffer sizes")
      allocate(src_field_buffer_dble(i,j)%p(src_field_buffer_sizes(i)))
      allocate(src_frac_mask_buffer_dble(i,j)%p(src_field_buffer_sizes(i)))
    end do
  end do

  call yac_fexchange_raw_frac_dble_ptr ( send_field_id,             &
                                         recv_field_id,             &
                                         send_nbr_pointsets,        &
                                         collection_size,           &
                                         send_field_dble,           &
                                         send_frac_mask_dble,       &
                                         src_field_buffer_dble,     &
                                         src_frac_mask_buffer_dble, &
                                         send_info,                 &
                                         recv_info,                 &
                                         ierror )

  do i = 1, send_nbr_pointsets
    do j = 1, collection_size
      src_field_buffer(i,j)%p = &
        real(src_field_buffer_dble(i,j)%p)
      deallocate(src_field_buffer_dble(i,j)%p)
      src_frac_mask_buffer(i,j)%p = &
        real(src_frac_mask_buffer_dble(i,j)%p)
      deallocate(src_frac_mask_buffer_dble(i,j)%p)
    end do
  end do

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
      deallocate(send_frac_mask_dble(j, i)%p)
    end do
  end do

end subroutine yac_fexchange_raw_frac_real_ptr

subroutine yac_fexchange_single_pointset_real ( send_field_id,       &
                                                recv_field_id,       &
                                                send_nbr_hor_points, &
                                                recv_nbr_hor_points, &
                                                collection_size,     &
                                                send_field,          &
                                                recv_field,          &
                                                send_info,           &
                                                recv_info,           &
                                                ierror )

  use yac, dummy => yac_fexchange_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id        !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id        !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points  !< [IN] number of horizontal send points
  integer, intent (in)  :: recv_nbr_hor_points  !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size      !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      collection_size)
                                                !< [IN] send field
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                                !< [INOUT] returned recv field
  integer, intent (out) :: send_info            !< [OUT] returned send info
  integer, intent (out) :: recv_info            !< [OUT] returned recv info
  integer, intent (out) :: ierror               !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, collection_size)

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble_single(send_field_id,        &
                                 send_nbr_hor_points,  &
                                 collection_size,      &
                                 send_field,           &
                                 send_buffer)

  call recv_field_to_dble(recv_field_id,       &
                          recv_nbr_hor_points, &
                          collection_size,     &
                          recv_field,          &
                          recv_buffer)

  call yac_fexchange_single_pointset_dble ( send_field_id,       &
                                            recv_field_id,       &
                                            send_nbr_hor_points, &
                                            recv_nbr_hor_points, &
                                            collection_size,     &
                                            send_buffer,         &
                                            recv_buffer,         &
                                            send_info,           &
                                            recv_info,           &
                                            ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_single_pointset_real

subroutine yac_fexchange_frac_single_pointset_real ( send_field_id,       &
                                                     recv_field_id,       &
                                                     send_nbr_hor_points, &
                                                     recv_nbr_hor_points, &
                                                     collection_size,     &
                                                     send_field,          &
                                                     send_frac_mask,      &
                                                     recv_field,          &
                                                     send_info,           &
                                                     recv_info,           &
                                                     ierror )

  use yac, dummy => yac_fexchange_frac_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id        !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id        !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points  !< [IN] number of horizontal send points
  integer, intent (in)  :: recv_nbr_hor_points  !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size      !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      collection_size)
                                                !< [IN] send field
  real, intent (in)     :: send_frac_mask(send_nbr_hor_points, &
                                          collection_size)
                                                !< [IN] fractional mask
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                                !< [INOUT] returned recv field
  integer, intent (out) :: send_info            !< [OUT] returned send info
  integer, intent (out) :: recv_info            !< [OUT] returned recv info
  integer, intent (out) :: ierror               !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, collection_size)
  double precision :: send_frac_mask_buffer(send_nbr_hor_points, collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, collection_size)

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble_single(send_field_id,        &
                                 send_nbr_hor_points,  &
                                 collection_size,      &
                                 send_field,           &
                                 send_buffer,          &
                                 send_frac_mask,       &
                                 send_frac_mask_buffer)

  call recv_field_to_dble(recv_field_id,       &
                          recv_nbr_hor_points, &
                          collection_size,     &
                          recv_field,          &
                          recv_buffer)

  call yac_fexchange_frac_single_pointset_dble ( send_field_id,         &
                                                 recv_field_id,         &
                                                 send_nbr_hor_points,   &
                                                 recv_nbr_hor_points,   &
                                                 collection_size,       &
                                                 send_buffer,           &
                                                 send_frac_mask_buffer, &
                                                 recv_buffer,           &
                                                 send_info,             &
                                                 recv_info,             &
                                                 ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_frac_single_pointset_real

subroutine yac_fexchange_raw_single_pointset_real ( send_field_id,         &
                                                    recv_field_id,         &
                                                    send_nbr_hor_points,   &
                                                    src_field_buffer_size, &
                                                    collection_size,       &
                                                    send_field,            &
                                                    src_field_buffer,      &
                                                    send_info,             &
                                                    recv_info,             &
                                                    ierror )

  use yac, dummy => yac_fexchange_raw_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id         !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: src_field_buffer_size !< [IN] source field buffer size
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  real, intent (in)     :: &
    send_field(send_nbr_hor_points, collection_size)
                                                 !< [IN] send field
  real, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] source field buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, collection_size)
  double precision :: &
    src_field_buffer_dble(src_field_buffer_size, collection_size)

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call send_field_to_dble_single(send_field_id,       &
                                 send_nbr_hor_points, &
                                 collection_size,     &
                                 send_field,          &
                                 send_buffer)

  call yac_fexchange_raw_single_pointset_dble ( send_field_id,         &
                                                recv_field_id,         &
                                                send_nbr_hor_points,   &
                                                src_field_buffer_size, &
                                                collection_size,       &
                                                send_buffer,           &
                                                src_field_buffer_dble, &
                                                send_info,             &
                                                recv_info,             &
                                                ierror )

  src_field_buffer = real(src_field_buffer_dble)

end subroutine yac_fexchange_raw_single_pointset_real

subroutine yac_fexchange_raw_frac_single_pointset_real ( send_field_id,         &
                                                          recv_field_id,         &
                                                          send_nbr_hor_points,   &
                                                          src_field_buffer_size, &
                                                          collection_size,       &
                                                          send_field,            &
                                                          send_frac_mask,        &
                                                          src_field_buffer,      &
                                                          src_frac_mask_buffer,  &
                                                          send_info,             &
                                                          recv_info,             &
                                                          ierror )

  use yac, dummy => yac_fexchange_raw_frac_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id         !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: src_field_buffer_size !< [IN] source buffer size
                                                  !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  real, intent (in)     :: &
    send_field(send_nbr_hor_points, collection_size)
                                                  !< [IN] send field
  real, intent (in)     :: &
    send_frac_mask(send_nbr_hor_points, collection_size)
                                                  !< [IN] fractional mask
  real, intent (out)    :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                  !< [OUT] returned source field buffer
  real, intent (out)    :: &
    src_frac_mask_buffer(src_field_buffer_size, collection_size)
                                                  !< [OUT] returned source fractional mask buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  double precision :: send_field_dble(send_nbr_hor_points, collection_size)
  double precision :: send_frac_mask_dble(send_nbr_hor_points, collection_size)
  double precision :: &
    src_field_buffer_dble(src_field_buffer_size, collection_size)
  double precision :: &
    src_frac_mask_buffer_dble(src_field_buffer_size, collection_size)

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call send_field_to_dble_single(send_field_id,       &
                                 send_nbr_hor_points, &
                                 collection_size,     &
                                 send_field,          &
                                 send_field_dble)
  call send_field_to_dble_single(send_field_id,       &
                                 send_nbr_hor_points, &
                                 collection_size,     &
                                 send_frac_mask,      &
                                 send_frac_mask_dble)

  call yac_fexchange_raw_frac_single_pointset_dble ( send_field_id,             &
                                                     recv_field_id,             &
                                                     send_nbr_hor_points,       &
                                                     src_field_buffer_size,     &
                                                     collection_size,           &
                                                     send_field_dble,           &
                                                     send_frac_mask_dble,       &
                                                     src_field_buffer_dble,     &
                                                     src_frac_mask_buffer_dble, &
                                                     send_info,                 &
                                                     recv_info,                 &
                                                     ierror )

  src_field_buffer = real(src_field_buffer_dble)
  src_frac_mask_buffer = real(src_frac_mask_buffer_dble)

end subroutine yac_fexchange_raw_frac_single_pointset_real

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_dble ( send_field_id,       &
                                recv_field_id,       &
                                send_nbr_hor_points, &
                                send_nbr_pointsets,  &
                                recv_nbr_hor_points, &
                                collection_size,     &
                                send_field,          &
                                recv_field,          &
                                send_info,           &
                                recv_info,           &
                                ierror )

  use yac, dummy => yac_fexchange_dble

  implicit none

  interface

     subroutine yac_cexchange__c ( send_field_id,   &
                                   recv_field_id,   &
                                   collection_size, &
                                   send_field,      &
                                   recv_field,      &
                                   send_info,       &
                                   recv_info,       &
                                   ierror )         &
         bind ( c, name='yac_cexchange_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::               &
                           send_field(           &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            collection_size)
  double precision, intent (inout)::             &
                           recv_field(           &
                            recv_nbr_hor_points, &
                            collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call yac_cexchange__c ( send_field_id,   &
                          recv_field_id,   &
                          collection_size, &
                          send_field,      &
                          recv_field,      &
                          send_info,       &
                          recv_info,       &
                          ierror )

end subroutine yac_fexchange_dble

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[in]    send_frac_mask      fractional mask
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_frac_dble ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_hor_points, &
                                     send_nbr_pointsets,  &
                                     recv_nbr_hor_points, &
                                     collection_size,     &
                                     send_field,          &
                                     send_frac_mask,      &
                                     recv_field,          &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

  use yac, dummy => yac_fexchange_frac_dble

  implicit none

  interface

     subroutine yac_cexchange_frac__c ( send_field_id,   &
                                        recv_field_id,   &
                                        collection_size, &
                                        send_field,      &
                                        send_frac_mask,  &
                                        recv_field,      &
                                        send_info,       &
                                        recv_info,       &
                                        ierror )         &
         bind ( c, name='yac_cexchange_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_frac__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::               &
                           send_field(           &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            collection_size)
  double precision, intent (in) ::               &
                           send_frac_mask(           &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            collection_size)
  double precision, intent (inout)::             &
                           recv_field(           &
                            recv_nbr_hor_points, &
                            collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call yac_cexchange_frac__c ( send_field_id,   &
                               recv_field_id,   &
                               collection_size, &
                               send_field,      &
                               send_frac_mask,  &
                               recv_field,      &
                               send_info,       &
                               recv_info,       &
                               ierror )

end subroutine yac_fexchange_frac_dble

subroutine yac_fexchange_raw_dble ( send_field_id,         &
                                    recv_field_id,         &
                                    send_nbr_hor_points,   &
                                    send_nbr_pointsets,    &
                                    src_field_buffer_size, &
                                    collection_size,       &
                                    send_field,            &
                                    src_field_buffer,      &
                                    send_info,             &
                                    recv_info,             &
                                    ierror )

  use yac, dummy => yac_fexchange_raw_dble

  implicit none

  interface

     subroutine yac_cexchange_raw__c ( send_field_id,    &
                                       recv_field_id,    &
                                       collection_size,  &
                                       send_field,       &
                                       src_field_buffer, &
                                       send_info,        &
                                       recv_info,        &
                                       ierror )          &
         bind ( c, name='yac_cexchange_raw_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: src_field_buffer(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_raw__c

  end interface

  integer, intent (in)  :: send_field_id         !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets    !< [IN] number of send point sets
  integer, intent (in)  :: src_field_buffer_size !< [IN] source field buffer size
                                                 !!      (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  double precision, intent (in)  :: &
    send_field(send_nbr_hor_points, send_nbr_pointsets, collection_size)
                                                 !< [IN] send field
  double precision, intent (out) :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] source field buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  integer :: i

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call yac_cexchange_raw__c ( send_field_id,    &
                              recv_field_id,    &
                              collection_size,  &
                              send_field,       &
                              src_field_buffer, &
                              send_info,        &
                              recv_info,        &
                              ierror )

end subroutine yac_fexchange_raw_dble

subroutine yac_fexchange_raw_frac_dble ( send_field_id,         &
                                         recv_field_id,         &
                                         send_nbr_hor_points,   &
                                         send_nbr_pointsets,    &
                                         src_field_buffer_size, &
                                         collection_size,       &
                                         send_field,            &
                                         send_frac_mask,        &
                                         src_field_buffer,      &
                                         src_frac_mask_buffer,  &
                                         send_info,             &
                                         recv_info,             &
                                         ierror )

  use yac, dummy => yac_fexchange_raw_frac_dble

  implicit none

  interface

     subroutine yac_cexchange_raw_frac__c ( send_field_id,        &
                                            recv_field_id,        &
                                            collection_size,      &
                                            send_field,           &
                                            send_frac_mask,       &
                                            src_field_buffer,     &
                                            src_frac_mask_buffer, &
                                            send_info,            &
                                            recv_info,            &
                                            ierror )              &
         bind ( c, name='yac_cexchange_raw_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       real    ( kind=c_double )     :: src_field_buffer(*)
       real    ( kind=c_double )     :: src_frac_mask_buffer(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_raw_frac__c

  end interface

  integer, intent (in)  :: send_field_id         !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets    !< [IN] number of send point sets
  integer, intent (in)  :: src_field_buffer_size !< [IN] source buffer size
                                                 !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  double precision, intent (in) :: &
    send_field(send_nbr_hor_points, send_nbr_pointsets, collection_size)
                                                 !< [IN] send field
  double precision, intent (in) :: &
    send_frac_mask(send_nbr_hor_points, send_nbr_pointsets, collection_size)
                                                 !< [IN] fractional mask
  double precision, intent (out):: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] returned source field buffer
  double precision, intent (out):: &
    src_frac_mask_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] returned source field buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  integer :: i

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call yac_cexchange_raw_frac__c ( send_field_id,        &
                                   recv_field_id,        &
                                   collection_size,      &
                                   send_field,           &
                                   send_frac_mask,       &
                                   src_field_buffer,     &
                                   src_frac_mask_buffer, &
                                   send_info,            &
                                   recv_info,            &
                                   ierror )

end subroutine yac_fexchange_raw_frac_dble

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field handle
!! @param[inout] recv_field          returned recv field handle
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_dble_ptr ( send_field_id,      &
                                    recv_field_id,      &
                                    send_nbr_pointsets, &
                                    collection_size,    &
                                    send_field,         &
                                    recv_field,         &
                                    send_info,          &
                                    recv_info,          &
                                    ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fexchange_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none

  interface

     subroutine yac_cexchange_ptr__c ( send_field_id,   &
                                       recv_field_id,   &
                                       collection_size, &
                                       send_field,      &
                                       recv_field,      &
                                       send_info,       &
                                       recv_info,       &
                                       ierror )         &
         bind ( c, name='yac_cexchange_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_ptr__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: collection_size
  type(yac_dble_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
  type(yac_dble_ptr)      :: recv_field(collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      send_field_(j, i) = yac_dble2cptr("yac_fexchange_dble_ptr", "send_field", send_field(j, i))
    end do
  end do
  do i = 1, collection_size
    recv_field_(i) = yac_dble2cptr("yac_fexchange_dble_ptr", "recv_field", recv_field(i))
  end do

  call yac_cexchange_ptr__c ( send_field_id,   &
                              recv_field_id,   &
                              collection_size, &
                              send_field_,     &
                              recv_field_,     &
                              send_info,       &
                              recv_info,       &
                              ierror )

end subroutine yac_fexchange_dble_ptr

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field handle
!! @param[in]    send_frac_mask      fractional mask handle
!! @param[inout] recv_field          returned recv field handle
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_frac_dble_ptr ( send_field_id,      &
                                         recv_field_id,      &
                                         send_nbr_pointsets, &
                                         collection_size,    &
                                         send_field,         &
                                         send_frac_mask,     &
                                         recv_field,         &
                                         send_info,          &
                                         recv_info,          &
                                         ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fexchange_frac_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none

  interface

     subroutine yac_cexchange_frac_ptr__c ( send_field_id,   &
                                            recv_field_id,   &
                                            collection_size, &
                                            send_field,      &
                                            send_frac_mask,  &
                                            recv_field,      &
                                            send_info,       &
                                            recv_info,       &
                                            ierror )         &
         bind ( c, name='yac_cexchange_frac_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: send_frac_mask(*)
       type(c_ptr)                   :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_frac_ptr__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: collection_size
  type(yac_dble_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
  type(yac_dble_ptr), intent (in) ::                      &
                           send_frac_mask(send_nbr_pointsets, &
                                          collection_size)
  type(yac_dble_ptr)      :: recv_field(collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: send_frac_mask_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_frac_mask(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      send_field_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_frac_dble_ptr", "send_field", send_field(j, i))
      send_frac_mask_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_frac_dble_ptr", "send_frac_mask", send_frac_mask(j, i))
    end do
  end do
  do i = 1, collection_size
    recv_field_(i) = yac_dble2cptr("yac_fexchange_frac_dble_ptr", "recv_field", recv_field(i))
  end do

  call yac_cexchange_frac_ptr__c ( send_field_id,   &
                                   recv_field_id,   &
                                   collection_size, &
                                   send_field_,     &
                                   send_frac_mask_, &
                                   recv_field_,     &
                                   send_info,       &
                                   recv_info,       &
                                   ierror )

end subroutine yac_fexchange_frac_dble_ptr

subroutine yac_fexchange_raw_dble_ptr ( send_field_id,      &
                                        recv_field_id,      &
                                        send_nbr_pointsets, &
                                        collection_size,    &
                                        send_field,         &
                                        src_field_buffer,   &
                                        send_info,          &
                                        recv_info,          &
                                        ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fexchange_raw_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

  implicit none

  interface

     subroutine yac_cexchange_raw_ptr__c ( send_field_id,    &
                                           recv_field_id,    &
                                           collection_size,  &
                                           send_field,       &
                                           src_field_buffer, &
                                           send_info,        &
                                           recv_info,        &
                                           ierror )          &
         bind ( c, name='yac_cexchange_raw_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: src_field_buffer(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_raw_ptr__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: collection_size
  type(yac_dble_ptr), intent (in) :: &
    send_field(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(send_nbr_pointsets, collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i, j
  integer :: src_field_buffer_sizes(send_nbr_pointsets)
  type(c_ptr) :: send_field_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: src_field_buffer_(send_nbr_pointsets, collection_size)

  src_field_buffer_sizes = &
    (/(SIZE(src_field_buffer(i, 1)%p),i=1,send_nbr_pointsets)/)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_sizes( &
    recv_field_id, send_nbr_pointsets, collection_size, src_field_buffer_sizes)

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      send_field_(j, i) = yac_dble2cptr("yac_fexchange_raw_dble_ptr", "send_field", send_field(j, i))
    end do
  end do
  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      YAC_FASSERT(src_field_buffer_sizes(j) == size(src_field_buffer(j,i)%p), "ERROR(yac_fexchange_raw_dble_ptr): inconsistent source buffer sizes")
      src_field_buffer_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_raw_dble_ptr", "src_field_buffer", src_field_buffer(j, i))
    end do
  end do

  call yac_cexchange_raw_ptr__c ( send_field_id,     &
                                  recv_field_id,     &
                                  collection_size,   &
                                  send_field_,       &
                                  src_field_buffer_, &
                                  send_info,         &
                                  recv_info,         &
                                  ierror )

end subroutine yac_fexchange_raw_dble_ptr


subroutine yac_fexchange_raw_frac_dble_ptr ( send_field_id,        &
                                             recv_field_id,        &
                                             send_nbr_pointsets,   &
                                             collection_size,      &
                                             send_field,           &
                                             send_frac_mask,       &
                                             src_field_buffer,     &
                                             src_frac_mask_buffer, &
                                             send_info,            &
                                             recv_info,            &
                                             ierror )

  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fexchange_raw_frac_dble_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

  implicit none

  interface

     subroutine yac_cexchange_raw_frac_ptr__c ( send_field_id,        &
                                                recv_field_id,        &
                                                collection_size,      &
                                                send_field,           &
                                                send_frac_mask,       &
                                                src_field_buffer,     &
                                                src_frac_mask_buffer, &
                                                send_info,            &
                                                recv_info,            &
                                                ierror )              &
         bind ( c, name='yac_cexchange_raw_frac_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: send_frac_mask(*)
       type(c_ptr)                   :: src_field_buffer(*)
       type(c_ptr)                   :: src_frac_mask_buffer(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_raw_frac_ptr__c

  end interface

  integer, intent (in)              :: send_field_id      !< [IN] send field identifier
  integer, intent (in)              :: recv_field_id      !< [IN] recv field identifier
  integer, intent (in)              :: send_nbr_pointsets !< [IN] number of send point sets
  integer, intent (in)              :: collection_size    !< [IN] number of vertical level or bundles
  type(yac_dble_ptr), intent (in)   :: &
    send_field(send_nbr_pointsets, collection_size)       !< [IN] send field
  type(yac_dble_ptr), intent (in)   :: &
    send_frac_mask(send_nbr_pointsets, collection_size)   !< [IN] fractional mask
  type(yac_dble_ptr), intent(inout) :: &
    src_field_buffer(send_nbr_pointsets, collection_size) !< [INOUT] returned source field buffer
  type(yac_dble_ptr), intent(inout) :: &
    src_frac_mask_buffer(send_nbr_pointsets, collection_size)
                                                          !< [INOUT] returned source fractional mask buffer
  integer, intent (out)             :: send_info          !< [OUT] returned send info
  integer, intent (out)             :: recv_info          !< [OUT] returned recv info
  integer, intent (out)             :: ierror             !< [OUT] returned error

  integer :: i, j
  integer :: src_field_buffer_sizes(send_nbr_pointsets)
  type(c_ptr) :: send_field_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: send_frac_mask_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: src_field_buffer_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: src_frac_mask_buffer_(send_nbr_pointsets, collection_size)

  src_field_buffer_sizes = &
    (/(SIZE(src_field_buffer(i, 1)%p),i=1,send_nbr_pointsets)/)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_frac_mask(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_src_field_buffer_sizes( &
    recv_field_id, send_nbr_pointsets, collection_size, src_field_buffer_sizes)

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      send_field_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_raw_frac_dble_ptr", "send_field", send_field(j, i))
      send_frac_mask_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_raw_frac_dble_ptr", "send_frac_mask", send_frac_mask(j, i))
    end do
  end do

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      YAC_FASSERT(src_field_buffer_sizes(j) == size(src_field_buffer(j,i)%p), "ERROR(yac_fexchange_raw_frac_dble_ptr): inconsistent source buffer sizes")
      YAC_FASSERT(src_field_buffer_sizes(j) == size(src_frac_mask_buffer(j,i)%p), "ERROR(yac_fexchange_raw_frac_dble_ptr): inconsistent source buffer sizes")
      src_field_buffer_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_raw_frac_dble_ptr", "src_field_buffer", src_field_buffer(j, i))
      src_frac_mask_buffer_(j, i) = &
        yac_dble2cptr( &
          "yac_fexchange_raw_frac_dble_ptr", "src_frac_mask_buffer", src_frac_mask_buffer(j, i))
    end do
  end do

  call yac_cexchange_raw_frac_ptr__c ( send_field_id,         &
                                       recv_field_id,         &
                                       collection_size,       &
                                       send_field_,           &
                                       send_frac_mask_,       &
                                       src_field_buffer_,     &
                                       src_frac_mask_buffer_, &
                                       send_info,             &
                                       recv_info,             &
                                       ierror )

end subroutine yac_fexchange_raw_frac_dble_ptr

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_single_pointset_dble ( send_field_id,       &
                                                recv_field_id,       &
                                                send_nbr_hor_points, &
                                                recv_nbr_hor_points, &
                                                collection_size,     &
                                                send_field,          &
                                                recv_field,          &
                                                send_info,           &
                                                recv_info,           &
                                                ierror )

  use yac, dummy => yac_fexchange_single_pointset_dble

  implicit none

  interface

     subroutine yac_cexchange__c ( send_field_id,   &
                                   recv_field_id,   &
                                   collection_size, &
                                   send_field,      &
                                   recv_field,      &
                                   send_info,       &
                                   recv_info,       &
                                   ierror )         &
         bind ( c, name='yac_cexchange_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::                         &
                           send_field(send_nbr_hor_points, &
                                      collection_size)
  double precision, intent (inout)::                       &
                           recv_field(recv_nbr_hor_points, &
                                      collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call yac_cexchange__c ( send_field_id,   &
                          recv_field_id,   &
                          collection_size, &
                          send_field,      &
                          recv_field,      &
                          send_info,       &
                          recv_info,       &
                          ierror )

end subroutine yac_fexchange_single_pointset_dble

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[in]    send_frac_mask      fractional mask
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_frac_single_pointset_dble ( send_field_id,       &
                                                     recv_field_id,       &
                                                     send_nbr_hor_points, &
                                                     recv_nbr_hor_points, &
                                                     collection_size,     &
                                                     send_field,          &
                                                     send_frac_mask,      &
                                                     recv_field,          &
                                                     send_info,           &
                                                     recv_info,           &
                                                     ierror )

  use yac, dummy => yac_fexchange_frac_single_pointset_dble

  implicit none

  interface

     subroutine yac_cexchange_frac__c ( send_field_id,   &
                                        recv_field_id,   &
                                        collection_size, &
                                        send_field,      &
                                        send_frac_mask,  &
                                        recv_field,      &
                                        send_info,       &
                                        recv_info,       &
                                        ierror )         &
         bind ( c, name='yac_cexchange_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_frac__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::                         &
                           send_field(send_nbr_hor_points, &
                                      collection_size)
  double precision, intent (in) ::                         &
                           send_frac_mask(send_nbr_hor_points, &
                                          collection_size)
  double precision, intent (inout)::                       &
                           recv_field(recv_nbr_hor_points, &
                                      collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call yac_cexchange_frac__c ( send_field_id,   &
                               recv_field_id,   &
                               collection_size, &
                               send_field,      &
                               send_frac_mask,  &
                               recv_field,      &
                               send_info,       &
                               recv_info,       &
                               ierror )

end subroutine yac_fexchange_frac_single_pointset_dble

subroutine yac_fexchange_raw_single_pointset_dble ( send_field_id,         &
                                                    recv_field_id,         &
                                                    send_nbr_hor_points,   &
                                                    src_field_buffer_size, &
                                                    collection_size,       &
                                                    send_field,            &
                                                    src_field_buffer,      &
                                                    send_info,             &
                                                    recv_info,             &
                                                    ierror )

  use yac, dummy => yac_fexchange_raw_single_pointset_dble

  implicit none

  interface

     subroutine yac_cexchange_raw__c ( send_field_id,    &
                                       recv_field_id,    &
                                       collection_size,  &
                                       send_field,       &
                                       src_field_buffer, &
                                       send_info,        &
                                       recv_info,        &
                                       ierror )          &
         bind ( c, name='yac_cexchange_raw_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: src_field_buffer(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_raw__c

  end interface

  integer, intent (in)  :: send_field_id         !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id         !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)  :: src_field_buffer_size !< [IN] source field buffer size
  integer, intent (in)  :: collection_size       !< [IN] number of vertical level or bundles
  double precision, intent (in) :: &
    send_field(send_nbr_hor_points, collection_size)
                                                 !< [IN] send field
  double precision, intent(out) :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                 !< [OUT] source field buffer
  integer, intent (out) :: send_info             !< [OUT] returned send info
  integer, intent (out) :: recv_info             !< [OUT] returned recv info
  integer, intent (out) :: ierror                !< [OUT] returned error

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call yac_cexchange_raw__c ( send_field_id,    &
                              recv_field_id,    &
                              collection_size,  &
                              send_field,       &
                              src_field_buffer, &
                              send_info,        &
                              recv_info,        &
                              ierror )

end subroutine yac_fexchange_raw_single_pointset_dble

subroutine yac_fexchange_raw_frac_single_pointset_dble ( send_field_id,         &
                                                          recv_field_id,         &
                                                          send_nbr_hor_points,   &
                                                          src_field_buffer_size, &
                                                          collection_size,       &
                                                          send_field,            &
                                                          send_frac_mask,        &
                                                          src_field_buffer,      &
                                                          src_frac_mask_buffer,  &
                                                          send_info,             &
                                                          recv_info,             &
                                                          ierror )

  use yac, dummy => yac_fexchange_raw_frac_single_pointset_dble

  implicit none

  interface

     subroutine yac_cexchange_raw_frac__c ( send_field_id,        &
                                            recv_field_id,        &
                                            collection_size,      &
                                            send_field,           &
                                            send_frac_mask,       &
                                            src_field_buffer,     &
                                            src_frac_mask_buffer, &
                                            send_info,            &
                                            recv_info,            &
                                            ierror )              &
         bind ( c, name='yac_cexchange_raw_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       real    ( kind=c_double )     :: src_field_buffer(*)
       real    ( kind=c_double )     :: src_frac_mask_buffer(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_raw_frac__c

  end interface

  integer, intent (in)           :: send_field_id         !< [IN] send field identifier
  integer, intent (in)           :: recv_field_id         !< [IN] recv field identifier
  integer, intent (in)           :: send_nbr_hor_points   !< [IN] number of horizontal send points
  integer, intent (in)           :: src_field_buffer_size !< [IN] source buffer size
                                                          !! (SUM(src_field_buffer_sizes(:)))
  integer, intent (in)           :: collection_size       !< [IN] number of vertical level or bundles
  double precision, intent (in)  :: &
    send_field(send_nbr_hor_points, collection_size)      !< [IN] send field
  double precision, intent (in)  :: &
    send_frac_mask(send_nbr_hor_points, collection_size)  !< [IN] fractional mask
  double precision, intent (out) :: &
    src_field_buffer(src_field_buffer_size, collection_size)
                                                          !< [OUT] returned source field buffer
  double precision, intent (out) :: &
    src_frac_mask_buffer(src_field_buffer_size, collection_size)
                                                          !< [OUT] returned source fractional mask buffer
  integer, intent (out)          :: send_info             !< [OUT] returned send info
  integer, intent (out)          :: recv_info             !< [OUT] returned recv info
  integer, intent (out)          :: ierror                !< [OUT] returned error

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_src_field_buffer_size( &
    recv_field_id, collection_size, src_field_buffer_size)

  call yac_cexchange_raw_frac__c ( send_field_id,        &
                                   recv_field_id,        &
                                   collection_size,      &
                                   send_field,           &
                                   send_frac_mask,       &
                                   src_field_buffer,     &
                                   src_frac_mask_buffer, &
                                   send_info,            &
                                   recv_info,            &
                                   ierror )

end subroutine yac_fexchange_raw_frac_single_pointset_dble

! ----------------------------------------------------------------------

subroutine yac_ftest_i ( field_id, flag )

  use yac, dummy => yac_ftest_i

  implicit none

  interface

     subroutine yac_ctest_c ( field_id, flag ) &
         bind ( c, name='yac_ctest' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int )        :: flag

     end subroutine yac_ctest_c

  end interface

  integer, intent (in)  :: field_id
  integer, intent (out) :: flag

  call yac_ctest_c ( field_id, flag )

end subroutine yac_ftest_i

subroutine yac_ftest_l ( field_id, flag )

  use yac, dummy => yac_ftest_l

  implicit none

  interface

     subroutine yac_ctest_c ( field_id, flag ) &
         bind ( c, name='yac_ctest' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int )        :: flag

     end subroutine yac_ctest_c

  end interface

  integer, intent (in)  :: field_id
  logical, intent (out) :: flag

  integer :: iflag

  call yac_ctest_c ( field_id, iflag )

  flag = iflag /= 0

end subroutine yac_ftest_l

! ----------------------------------------------------------------------

subroutine yac_fwait ( field_id )

  use yac, dummy => yac_fwait

  implicit none

  interface

     subroutine yac_cwait_c ( field_id ) &
         bind ( c, name='yac_cwait' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id

     end subroutine yac_cwait_c

  end interface

  integer, intent (in)  :: field_id

  call yac_cwait_c ( field_id )

end subroutine yac_fwait

! ----------------------------------------------------------------------

subroutine yac_fget_comp_comm ( comp_id, comp_comm )

   use yac, dummy => yac_fget_comp_comm

   implicit none

   interface

     subroutine yac_get_comp_comm_c ( comp_id, comp_comm ) &
       bind ( c, name='yac_get_comp_comm_f2c' )

       use, intrinsic :: iso_c_binding, only : c_int
       use yac, only : YAC_MPI_FINT_KIND

       integer ( kind=c_int ), value      :: comp_id
       integer ( kind=YAC_MPI_FINT_KIND ) :: comp_comm

     end subroutine yac_get_comp_comm_c

   end interface

   integer, intent (in)  :: comp_id   !< [IN] component identifier
   integer, intent (out) :: comp_comm !< [OUT] component MPI communicator

   call yac_get_comp_comm_c ( comp_id, comp_comm )

end subroutine yac_fget_comp_comm

! ----------------------------------------------------------------------

subroutine yac_fget_comps_comm ( comp_names, num_comps, comps_comm )

   use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
   use yac, dummy => yac_fget_comps_comm
   use mo_yac_iso_c_helpers

   implicit none

   interface

     subroutine yac_cget_comps_comm_c ( comp_names, &
                                        num_comps,  &
                                        comps_comm) &
           bind ( c, name='yac_cget_comps_comm_f2c' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr
       use yac, only : YAC_MPI_FINT_KIND

       type ( c_ptr )                     :: comp_names(*)
       integer ( kind=c_int ), value      :: num_comps
       integer ( kind=YAC_MPI_FINT_KIND ) :: comps_comm

     end subroutine yac_cget_comps_comm_c

   end interface

   integer, intent(in)   :: num_comps             !< [IN]  number of components
   character(kind=c_char, len=*), intent(in) :: &
                            comp_names(num_comps) !< [IN]  component names
   integer, intent (out) :: comps_comm            !< [OUT] components MPI communicator

   integer :: i, j
   character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
   type(c_ptr) :: comp_name_ptrs(num_comps)

   comp_names_cpy = c_null_char

   do i = 1, num_comps
     YAC_CHECK_STRING_LEN ( "yac_fget_comps_comm", comp_names(i))
     do j = 1, len_trim(comp_names(i))
       comp_names_cpy(j,i) = comp_names(i)(j:j)
     end do
     comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
   end do

   call yac_cget_comps_comm_c ( comp_name_ptrs, &
                                num_comps,      &
                                comps_comm )

end subroutine yac_fget_comps_comm

subroutine yac_fget_comps_comm_instance ( yac_instance_id, &
                                          comp_names,      &
                                          num_comps,       &
                                          comps_comm )

   use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
   use yac, dummy => yac_fget_comps_comm_instance
   use mo_yac_iso_c_helpers

   implicit none

   interface

     subroutine yac_cget_comps_comm_instance_c ( yac_instance_id, &
                                                 comp_names,      &
                                                 num_comps,       &
                                                 comps_comm)      &
           bind ( c, name='yac_cget_comps_comm_instance_f2c' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr
       use yac, only : YAC_MPI_FINT_KIND

       integer   ( kind=c_int ), value    :: yac_instance_id
       type ( c_ptr )                     :: comp_names(*)
       integer ( kind=c_int ), value      :: num_comps
       integer ( kind=YAC_MPI_FINT_KIND ) :: comps_comm

     end subroutine yac_cget_comps_comm_instance_c

   end interface

   integer, intent(in)   :: yac_instance_id       !< [IN]  YAC instance identifier
   integer, intent(in)   :: num_comps             !< [IN]  number of components
   character(kind=c_char, len=*), intent(in) :: &
                            comp_names(num_comps) !< [IN]  component names
   integer, intent (out) :: comps_comm            !< [OUT] components MPI communicator

   integer :: i, j
   character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
   type(c_ptr) :: comp_name_ptrs(num_comps)

   comp_names_cpy = c_null_char

   do i = 1, num_comps
     YAC_CHECK_STRING_LEN ( "yac_fget_comps_comm_instance", comp_names(i))
     do j = 1, len_trim(comp_names(i))
       comp_names_cpy(j,i) = comp_names(i)(j:j)
     end do
     comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
   end do

   call yac_cget_comps_comm_instance_c ( yac_instance_id, &
                                         comp_name_ptrs,  &
                                         num_comps,       &
                                         comps_comm )

end subroutine yac_fget_comps_comm_instance

! ------------------- search/end of definition -------------------------

subroutine yac_fsync_def ( )

   use yac, dummy => yac_fsync_def

   implicit none

   interface

     subroutine yac_csync_def_c ( ) bind ( c, name='yac_csync_def' )

     end subroutine yac_csync_def_c

   end interface

   call yac_csync_def_c ( )

end subroutine yac_fsync_def

subroutine yac_fsync_def_instance ( yac_instance_id )

   use yac, dummy => yac_fsync_def_instance

   implicit none

   interface

     subroutine yac_csync_def_instance_c ( yac_instance_id ) &
       bind ( c, name='yac_csync_def_instance' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: yac_instance_id

     end subroutine yac_csync_def_instance_c

   end interface

   integer, intent(in)  :: yac_instance_id !< [IN]  YAC instance identifier

   call yac_csync_def_instance_c ( yac_instance_id )

end subroutine yac_fsync_def_instance


subroutine yac_fenddef ( )

   use yac, dummy => yac_fenddef

   implicit none

   interface

     subroutine yac_cenddef_c ( ) bind ( c, name='yac_cenddef' )

     end subroutine yac_cenddef_c

   end interface

   call yac_cenddef_c ( )

end subroutine yac_fenddef

subroutine yac_fenddef_instance ( yac_instance_id )

   use yac, dummy => yac_fenddef_instance

   implicit none

   interface

     subroutine yac_cenddef_instance_c ( yac_instance_id ) &
       bind ( c, name='yac_cenddef_instance' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: yac_instance_id

     end subroutine yac_cenddef_instance_c

   end interface

   integer, intent(in)  :: yac_instance_id !< [IN]  YAC instance identifier

   call yac_cenddef_instance_c ( yac_instance_id )

end subroutine yac_fenddef_instance

subroutine yac_fenddef_and_emit_config(emit_flags, config)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use yac, dummy => yac_fenddef_and_emit_config
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_fenddef_and_emit_config_c ( &
       emit_flags, config)                    &
       bind ( c, name='yac_cenddef_and_emit_config' )

     use, intrinsic :: iso_c_binding, only : c_int, c_ptr
     integer ( kind=c_int ), value :: emit_flags
     type(c_ptr) :: config

     end subroutine yac_fenddef_and_emit_config_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  integer, intent (in)           :: emit_flags !< [IN] flags for emitting the config
  character (len=:), ALLOCATABLE :: config     !< [IN,OUT] configuration string

  type (c_ptr)                   :: c_string_ptr

  call yac_fenddef_and_emit_config_c(emit_flags, c_string_ptr)
  config = yac_internal_cptr2char(c_string_ptr)
  call free_c(c_string_ptr)

end subroutine yac_fenddef_and_emit_config

subroutine yac_fenddef_and_emit_config_instance( &
  yac_instance_id, emit_flags, config)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use yac, dummy => yac_fenddef_and_emit_config_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_fenddef_and_emit_config_instance_c ( &
       yac_instance_id, emit_flags, config)            &
       bind ( c, name='yac_cenddef_and_emit_config_instance' )

     use, intrinsic :: iso_c_binding, only : c_int, c_ptr
     integer ( kind=c_int ), value :: yac_instance_id
     integer ( kind=c_int ), value :: emit_flags
     type(c_ptr) :: config

     end subroutine yac_fenddef_and_emit_config_instance_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  integer, intent (in)           :: yac_instance_id !< [IN] YAC instance identifier
  integer, intent (in)           :: emit_flags      !< [IN] flags for emitting the config
  character (len=:), ALLOCATABLE :: config          !< [IN,OUT] configuration string

  type (c_ptr)                   :: c_string_ptr

  call yac_fenddef_and_emit_config_instance_c( &
    yac_instance_id, emit_flags, c_string_ptr)
  config = yac_internal_cptr2char(c_string_ptr)
  call free_c(c_string_ptr)

end subroutine yac_fenddef_and_emit_config_instance

! ------------------------ query routines -----------------------------

function yac_fget_grid_size( location, grid_id ) result (grid_size)

  use yac, dummy => yac_fget_grid_size

  use, intrinsic :: iso_c_binding, ONLY: c_size_t, c_int, c_null_char
  implicit none

  interface
    function yac_cget_grid_size_c(location, grid_id) result (grid_size) &
        bind(c, name='yac_cget_grid_size')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int
      integer(kind=c_int), value :: location
      integer(kind=c_int), value :: grid_id
      integer(kind=c_size_t) :: grid_size
    end function yac_cget_grid_size_c
  end interface

  integer, intent(in) :: location
  integer, intent(in) :: grid_id
  integer :: grid_size

  integer(kind=c_size_t) :: c_grid_size

  c_grid_size = &
    yac_cget_grid_size_c(int(location, c_int), int(grid_id, c_int))

  YAC_FASSERT(INT(HUGE(grid_size), c_size_t) >= c_grid_size, "ERROR(yac_fget_grid_size): grid size exceeds HUGE(grid_size)")

  grid_size = INT(c_grid_size)

end function yac_fget_grid_size

! ---------------------------------------------------------------------

subroutine yac_fcompute_grid_cell_areas_real( &
  grid_id, nbr_cells, cell_areas )

  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_null_char
  use yac, dummy => yac_fcompute_grid_cell_areas_real

  implicit none

  interface

    subroutine yac_ccompute_grid_cell_areas_c ( grid_id, cell_areas )&
        bind ( c, name='yac_ccompute_grid_cell_areas' )

      use, intrinsic :: iso_c_binding, only : c_int, c_double

      integer ( kind=c_int), value :: grid_id
      real    ( kind=c_double)     :: cell_areas(*)

    end subroutine yac_ccompute_grid_cell_areas_c

  end interface

  integer, intent(in) :: grid_id
  integer, intent(in) :: nbr_cells
  real, intent(out)   :: cell_areas(nbr_cells)

  real(kind=c_double), allocatable :: c_cell_areas(:)
  integer :: ref_nbr_cells

  ref_nbr_cells = yac_fget_grid_size(YAC_LOCATION_CELL, grid_id)

  YAC_FASSERT(nbr_cells == ref_nbr_cells, "ERROR(yac_fcompute_grid_cell_areas_real): wrong number of cells for provided grid")

  allocate(c_cell_areas(nbr_cells))

  CALL yac_ccompute_grid_cell_areas_c(int(grid_id, c_int), c_cell_areas)

  cell_areas(:) = REAL(c_cell_areas(:))

  deallocate(c_cell_areas)

end subroutine yac_fcompute_grid_cell_areas_real


subroutine yac_fcompute_grid_cell_areas_dble( &
  grid_id, nbr_cells, cell_areas )

  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_null_char
  use yac, dummy => yac_fcompute_grid_cell_areas_dble

  implicit none

  interface

    subroutine yac_ccompute_grid_cell_areas_c ( grid_id, cell_areas )&
        bind ( c, name='yac_ccompute_grid_cell_areas' )

      use, intrinsic :: iso_c_binding, only : c_int, c_double

      integer ( kind=c_int), value :: grid_id
      real    ( kind=c_double)     :: cell_areas(*)

    end subroutine yac_ccompute_grid_cell_areas_c

  end interface

  integer, intent(in) :: grid_id
  integer, intent(in) :: nbr_cells
  double precision, intent(out)   :: cell_areas(nbr_cells)

  real(kind=c_double), allocatable :: c_cell_areas(:)
  integer :: ref_nbr_cells

  ref_nbr_cells = yac_fget_grid_size(YAC_LOCATION_CELL, grid_id)

  YAC_FASSERT(nbr_cells == ref_nbr_cells, "ERROR(yac_fcompute_grid_cell_areas_dble): wrong number of cells for provided grid")

  allocate(c_cell_areas(nbr_cells))

  CALL yac_ccompute_grid_cell_areas_c(int(grid_id, c_int), c_cell_areas)

  cell_areas(:) = DBLE(c_cell_areas(:))

  deallocate(c_cell_areas)

end subroutine yac_fcompute_grid_cell_areas_dble

! ---------------------------------------------------------------------

function yac_fget_points_size( point_id ) result (points_size)

  use yac, dummy => yac_fget_points_size

  use, intrinsic :: iso_c_binding, ONLY: c_size_t, c_int, c_null_char
  implicit none

  interface
    function yac_cget_points_size_c(point_id) result (points_size) &
        bind(c, name='yac_cget_points_size')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int
      integer(kind=c_int), value :: point_id
      integer(kind=c_size_t) :: points_size
    end function yac_cget_points_size_c
  end interface

  integer, intent(in) :: point_id
  integer :: points_size

  integer(kind=c_size_t) :: c_points_size

  c_points_size = yac_cget_points_size_c(int(point_id, c_int))

  YAC_FASSERT(INT(HUGE(points_size), c_size_t) >= c_points_size, "ERROR(yac_fget_point_size): point size exceeds HUGE(point_size)")

  points_size = INT(c_points_size)

end function yac_fget_points_size

! ---------------------------------------------------------------------

 function yac_fget_comp_names ( ) result( comp_names )

   use yac, dummy => yac_fget_comp_names
   use mo_yac_iso_c_helpers

   use, intrinsic :: iso_c_binding, ONLY: c_ptr
   implicit none

   interface
      function yac_cget_nbr_comps_c() result( nbr_comps ) &
           bind(c, name='yac_cget_nbr_comps')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int) :: nbr_comps

      end function yac_cget_nbr_comps_c

      subroutine yac_cget_comp_names_c( nbr_comps, comp_names ) &
           bind(c, name='yac_cget_comp_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: nbr_comps
        TYPE(c_ptr), intent(out) :: comp_names(nbr_comps)
      end subroutine yac_cget_comp_names_c
   end interface

   type(yac_string), allocatable :: comp_names(:)
   integer :: nbr_comps
   INTEGER :: i
   TYPE(c_ptr), allocatable :: comp_ptr(:)

   nbr_comps = yac_cget_nbr_comps_c()
   allocate(comp_ptr(nbr_comps))
   allocate(comp_names(nbr_comps))
   CALL yac_cget_comp_names_c(nbr_comps, comp_ptr)
   DO i=1,nbr_comps
      comp_names(i)%string = yac_internal_cptr2char(comp_ptr(i))
   END DO
 end function yac_fget_comp_names

 function yac_fget_comp_names_instance ( yac_instance_id ) result ( comp_names )

   use yac, dummy => yac_fget_comp_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr
   implicit none

   interface
      function yac_cget_nbr_comps_instance_c( yac_instance_id ) &
          result( nbr_comps )                                   &
           bind(c, name='yac_cget_nbr_comps_instance')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        integer(kind=c_int) :: nbr_comps
      end function yac_cget_nbr_comps_instance_c

      subroutine yac_cget_comp_names_instance_c( yac_instance_id, &
                                                 nbr_comps,       &
                                                 comp_names )     &
           bind(c, name='yac_cget_comp_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: yac_instance_id
        integer(kind=c_int), intent(in), value :: nbr_comps
        TYPE(c_ptr), intent(out) :: comp_names(nbr_comps)
      end subroutine yac_cget_comp_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   type(yac_string), allocatable :: comp_names(:)
   integer :: nbr_comps
   INTEGER :: i
   TYPE(c_ptr), allocatable :: comp_ptr(:)

   nbr_comps = yac_cget_nbr_comps_instance_c(yac_instance_id)
   allocate(comp_names(nbr_comps))
   allocate(comp_ptr(nbr_comps))
   CALL yac_cget_comp_names_instance_c(yac_instance_id, &
                                       nbr_comps,       &
                                       comp_ptr)
   DO i=1,nbr_comps
      comp_names(i)%string = yac_internal_cptr2char(comp_ptr(i))
   END DO
 end function  yac_fget_comp_names_instance

! ---------------------------------------------------------------------

 function yac_fget_grid_names ( ) result ( grid_names )

   use yac, dummy => yac_fget_grid_names
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr

   implicit none

   interface
      function yac_cget_nbr_grids_c() result( nbr_grids ) &
           bind(c, name='yac_cget_nbr_grids')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int) :: nbr_grids
      end function yac_cget_nbr_grids_c

      subroutine yac_cget_grid_names_c( nbr_grids, grid_names ) &
           bind(c, name='yac_cget_grid_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_grid_names_c
   end interface

   type(yac_string), allocatable :: grid_names(:)
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids = yac_cget_nbr_grids_c()
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_grid_names_c(nbr_grids, grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_grid_names

 function yac_fget_grid_names_instance ( yac_instance_id ) result ( grid_names )

   use yac, dummy => yac_fget_grid_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr

   implicit none

   interface
      function yac_cget_nbr_grids_instance_c( yac_instance_id ) &
           result( nbr_grids )                                  &
           bind(c, name='yac_cget_nbr_grids_instance')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        integer(kind=c_int) :: nbr_grids

      end function yac_cget_nbr_grids_instance_c

      subroutine yac_cget_grid_names_instance_c( yac_instance_id, &
                                                 nbr_grids,       &
                                                 grid_names )     &
           bind(c, name='yac_cget_grid_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: yac_instance_id
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_grid_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   type(yac_string), allocatable :: grid_names(:)
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids = yac_cget_nbr_grids_instance_c(yac_instance_id)
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_grid_names_instance_c(yac_instance_id, &
                                       nbr_grids,       &
                                       grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_grid_names_instance

 ! ---------------------------------------------------------------------

 function yac_fget_comp_grid_names ( comp_name ) result ( grid_names )

   use yac, dummy => yac_fget_comp_grid_names
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_comp_nbr_grids_c( comp_name ) result( nbr_grids ) &
           bind(c, name='yac_cget_comp_nbr_grids')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int) :: nbr_grids

      end function yac_cget_comp_nbr_grids_c

      subroutine yac_cget_comp_grid_names_c( comp_name,   &
                                             nbr_grids,   &
                                             grid_names ) &
           bind(c, name='yac_cget_comp_grid_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_comp_grid_names_c
   end interface

   type(yac_string), allocatable :: grid_names(:)
   CHARACTER(len=*), intent(in)  :: comp_name
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids = yac_cget_comp_nbr_grids_c(TRIM(comp_name) // c_null_char)
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_comp_grid_names_c(TRIM(comp_name) // c_null_char, &
                                   nbr_grids,                      &
                                   grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_comp_grid_names

 function yac_fget_comp_grid_names_instance ( yac_instance_id, comp_name) result ( grid_names )

   use yac, dummy => yac_fget_comp_grid_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_comp_nbr_grids_instance_c( yac_instance_id, &
                                                   comp_name )      &
           result( nbr_grids )                                      &
           bind(c, name='yac_cget_comp_nbr_grids_instance')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        integer(kind=c_int) :: nbr_grids
      end function yac_cget_comp_nbr_grids_instance_c

      subroutine yac_cget_comp_grid_names_instance_c( yac_instance_id, &
                                                      comp_name,       &
                                                      nbr_grids,       &
                                                      grid_names )     &
           bind(c, name='yac_cget_comp_grid_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int), intent(in), value :: yac_instance_id
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_comp_grid_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   type(yac_string), allocatable :: grid_names(:)
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids =                                          &
    yac_cget_comp_nbr_grids_instance_c(yac_instance_id, &
                                       TRIM(comp_name) // c_null_char)
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_comp_grid_names_instance_c(yac_instance_id,                &
                                            TRIM(comp_name) // c_null_char, &
                                            nbr_grids,                      &
                                            grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_comp_grid_names_instance

 ! ---------------------------------------------------------------------

 function yac_fget_field_names ( comp_name, grid_name ) result( field_names )

   use yac, dummy => yac_fget_field_names
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_nbr_fields_c(comp_name, grid_name) &
           result( nbr_fields )                            &
           bind(c, name='yac_cget_nbr_fields')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int) :: nbr_fields
      end function yac_cget_nbr_fields_c

      subroutine yac_cget_field_names_c( comp_name,   &
                                         grid_name,   &
                                        nbr_fields,   &
                                        field_names ) &
           bind(c, name='yac_cget_field_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int), intent(in), value :: nbr_fields
        TYPE(c_ptr), intent(out) :: field_names(nbr_fields)
      end subroutine yac_cget_field_names_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   type(yac_string), allocatable :: field_names(:)
   integer :: nbr_fields
   INTEGER :: i
   TYPE(c_ptr), allocatable :: field_ptr(:)

   nbr_fields = yac_cget_nbr_fields_c(TRIM(comp_name)//c_null_char, &
                                      TRIM(grid_name)//c_null_char)
   allocate(field_ptr(nbr_fields))
   CALL yac_cget_field_names_c(TRIM(comp_name)//c_null_char, &
                               TRIM(grid_name)//c_null_char, &
                               nbr_fields,                   &
                               field_ptr)
   allocate(field_names(nbr_fields))
   DO i=1,nbr_fields
      field_names(i)%string = yac_internal_cptr2char(field_ptr(i))
   END DO
 end function yac_fget_field_names

 function yac_fget_field_names_instance ( yac_instance_id, &
                                          comp_name,       &
                                          grid_name )      &
      result( field_names )

   use yac, dummy => yac_fget_field_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_nbr_fields_instance_c( yac_instance_id, &
                                               comp_name,       &
                                               grid_name )      &
           result( nbr_fields ) bind(c, name='yac_cget_nbr_fields_instance')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int) :: nbr_fields
      end function yac_cget_nbr_fields_instance_c

      subroutine yac_cget_field_names_instance_c( yac_instance_id, &
                                                  comp_name,       &
                                                  grid_name,       &
                                                  nbr_fields,      &
                                                  field_names )    &
           bind(c, name='yac_cget_field_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int), intent(in), value :: nbr_fields
        TYPE(c_ptr), intent(out) :: field_names(nbr_fields)
      end subroutine yac_cget_field_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   type(yac_string), allocatable :: field_names(:)
   integer :: nbr_fields
   INTEGER :: i
   TYPE(c_ptr), allocatable :: field_ptr(:)

   nbr_fields =                                                   &
     yac_cget_nbr_fields_instance_c(yac_instance_id,              &
                                    TRIM(comp_name)//c_null_char, &
                                    TRIM(grid_name)//c_null_char)
   allocate(field_ptr(nbr_fields))
   CALL yac_cget_field_names_instance_c(yac_instance_id,              &
                                        TRIM(comp_name)//c_null_char, &
                                        TRIM(grid_name)//c_null_char, &
                                        nbr_fields,                   &
                                        field_ptr)
   allocate(field_names(nbr_fields))
   DO i=1,nbr_fields
      field_names(i)%string = yac_internal_cptr2char(field_ptr(i))
   END DO
 end function yac_fget_field_names_instance

! ---------------------------------------------------------------------

function yac_fget_field_is_defined ( comp_name, grid_name, field_name ) &
     result(field_is_defined)

  use yac, dummy => yac_fget_field_is_defined
  use, intrinsic :: iso_c_binding, only : c_null_char, c_int

  implicit none

  interface

     function yac_cget_field_is_defined_c ( comp_name, grid_name, field_name ) &
          bind ( c, name='yac_cget_field_is_defined' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: comp_name
       character ( kind=c_char ), dimension(*) :: grid_name
       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int ) :: yac_cget_field_is_defined_c

     end function yac_cget_field_is_defined_c

  end interface

  character(len=*), intent (in)  :: comp_name
  character(len=*), intent (in)  :: grid_name
  character(len=*), intent (in)  :: field_name
  logical :: field_is_defined

  field_is_defined =                                                       &
    0_c_int /= yac_cget_field_is_defined_c ( TRIM(comp_name)//c_null_char, &
                                             TRIM(grid_name)//c_null_char, &
                                             TRIM(field_name)//c_null_char )

end function yac_fget_field_is_defined

function yac_fget_field_is_defined_instance ( yac_id,      &
                                              comp_name,   &
                                              grid_name,   &
                                              field_name ) &
     result(field_is_defined)

  use yac, dummy => yac_fget_field_is_defined_instance
  use, intrinsic :: iso_c_binding, only : c_null_char, c_int

  implicit none

  interface

     function yac_cget_field_is_defined_instance_c ( yac_id,      &
                                                     comp_name,   &
                                                     grid_name,   &
                                                     field_name ) &
          bind ( c, name='yac_cget_field_is_defined_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), value, intent(in) :: yac_id
       character ( kind=c_char ), dimension(*) :: comp_name
       character ( kind=c_char ), dimension(*) :: grid_name
       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int )  :: yac_cget_field_is_defined_instance_c

     end function yac_cget_field_is_defined_instance_c

  end interface

  integer, intent(in) :: yac_id
  character(len=*), intent (in)  :: comp_name
  character(len=*), intent (in)  :: grid_name
  character(len=*), intent (in)  :: field_name
  logical :: field_is_defined

  field_is_defined =                                                         &
    0_c_int /= yac_cget_field_is_defined_instance_c ( yac_id,                &
                                               TRIM(comp_name)//c_null_char, &
                                               TRIM(grid_name)//c_null_char, &
                                               TRIM(field_name)//c_null_char )

end function yac_fget_field_is_defined_instance

! ---------------------------------------------------------------------

function yac_fget_field_id ( comp_name, grid_name, field_name ) &
     result(field_id)

  use yac, dummy => yac_fget_field_id
  use, intrinsic :: iso_c_binding, only : c_null_char

  implicit none

  interface

     function yac_cget_field_id_c ( comp_name, grid_name, field_name ) &
          result(field_id) &
          bind ( c, name='yac_cget_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: comp_name
       character ( kind=c_char ), dimension(*) :: grid_name
       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int ) :: field_id

     end function yac_cget_field_id_c

  end interface

  character(len=*), intent (in)  :: comp_name
  character(len=*), intent (in)  :: grid_name
  character(len=*), intent (in)  :: field_name
  integer :: field_id

  field_id = yac_cget_field_id_c ( TRIM(comp_name)//c_null_char, &
                                   TRIM(grid_name)//c_null_char, &
                                   TRIM(field_name)//c_null_char )

end function yac_fget_field_id

function yac_fget_field_id_instance ( yac_id,      &
                                      comp_name,   &
                                      grid_name,   &
                                      field_name ) &
     result(field_id)

  use yac, dummy => yac_fget_field_id_instance
  use, intrinsic :: iso_c_binding, only : c_null_char

  implicit none

  interface

     function yac_cget_field_id_instance_c ( yac_id,      &
                                             comp_name,   &
                                             grid_name,   &
                                             field_name ) &
          result(field_id)                                &
          bind ( c, name='yac_cget_field_id_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), value, intent(in) :: yac_id
       character ( kind=c_char ), dimension(*) :: comp_name
       character ( kind=c_char ), dimension(*) :: grid_name
       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int )  :: field_id

     end function yac_cget_field_id_instance_c

  end interface

  integer, intent(in) :: yac_id
  character(len=*), intent (in)  :: comp_name
  character(len=*), intent (in)  :: grid_name
  character(len=*), intent (in)  :: field_name
  integer :: field_id

  field_id =                                                     &
    yac_cget_field_id_instance_c ( yac_id,                       &
                                   TRIM(comp_name)//c_null_char, &
                                   TRIM(grid_name)//c_null_char, &
                                   TRIM(field_name)//c_null_char )

end function yac_fget_field_id_instance

! ---------------------------------------------------------------------

subroutine yac_fget_action ( field_id, action )

   use yac, dummy => yac_fget_action

   implicit none

   interface

    subroutine yac_cget_action_c ( field_id, action ) &
      bind ( c, name='yac_cget_action' )

      use, intrinsic :: iso_c_binding, only : c_int

      integer ( kind=c_int ), value          :: field_id
      integer ( kind=c_int)                  :: action

    end subroutine yac_cget_action_c

   end interface

  integer, intent (in)  :: field_id !< [IN]  field identifier
  integer, intent (out) :: action   !< [OUT] action for the current timestep\n
                                    !!       (\ref YAC_ACTION_NONE,
                                    !!        \ref YAC_ACTION_COUPLING,
                                    !!        \ref YAC_ACTION_GET_FOR_RESTART,
                                    !!        \ref YAC_ACTION_PUT_FOR_RESTART, or
                                    !!        \ref YAC_ACTION_OUT_OF_BOUND)

  call yac_cget_action_c(field_id, action)

end subroutine yac_fget_action

! ---------------------------------------------------------------------

subroutine yac_fupdate ( field_id )

   use yac, dummy => yac_fupdate

   implicit none

   interface

    subroutine yac_cupdate_c ( field_id ) &
      bind ( c, name='yac_cupdate' )

      use, intrinsic :: iso_c_binding, only : c_int

      integer ( kind=c_int ), value :: field_id

    end subroutine yac_cupdate_c

   end interface

  integer, intent (in)  :: field_id !< [IN]  field identifier

  call yac_cupdate_c(field_id)

end subroutine yac_fupdate

subroutine yac_fdef_couple( src_comp_name, src_grid_name, src_field_name, &
  tgt_comp_name, tgt_grid_name, tgt_field_name, &
  coupling_timestep, time_unit, time_reduction, interp_stack_config_id, &
  src_lag, tgt_lag, weight_file, weight_file_on_existing, mapping_side, &
  scale_factor, scale_summand, src_mask_names, tgt_mask_name, &
  yaxt_exchanger_name, use_raw_exchange)

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_null_ptr, c_loc
  use yac, dummy => yac_fdef_couple
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_couple__c( src_comp_name,           &
                                    src_grid_name,           &
                                    src_field_name,          &
                                    tgt_comp_name,           &
                                    tgt_grid_name,           &
                                    tgt_field_name,          &
                                    coupling_timestep,       &
                                    time_unit,               &
                                    time_reduction,          &
                                    interp_stack_config_id,  &
                                    src_lag,                 &
                                    tgt_lag,                 &
                                    weight_file,             &
                                    weight_file_on_existing, &
                                    mapping_side,            &
                                    scale_factor,            &
                                    scale_summand,           &
                                    num_src_mask_names,      &
                                    src_mask_names,          &
                                    tgt_mask_name,           &
                                    yaxt_exchanger_name,     &
                                    use_raw_exchange)        &
       bind ( c, name='yac_cdef_couple_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr, c_double

       character ( kind=c_char ), dimension(*) :: src_comp_name
       character ( kind=c_char ), dimension(*) :: src_grid_name
       character ( kind=c_char ), dimension(*) :: src_field_name
       character ( kind=c_char ), dimension(*) :: tgt_comp_name
       character ( kind=c_char ), dimension(*) :: tgt_grid_name
       character ( kind=c_char ), dimension(*) :: tgt_field_name
       character ( kind=c_char ), dimension(*) :: coupling_timestep
       integer ( kind=c_int ), value           :: time_unit
       integer ( kind=c_int ), value           :: time_reduction
       integer ( kind=c_int ), value           :: interp_stack_config_id
       integer ( kind=c_int ), value           :: src_lag
       integer ( kind=c_int ), value           :: tgt_lag
       type ( c_ptr ), value                   :: weight_file
       integer ( kind=c_int ), value           :: weight_file_on_existing
       integer ( kind=c_int ), value           :: mapping_side
       real ( kind=c_double ), value           :: scale_factor
       real ( kind=c_double ), value           :: scale_summand
       integer ( kind=c_int ), value           :: num_src_mask_names
       type ( c_ptr )                          :: src_mask_names(*)
       type ( c_ptr ), value                   :: tgt_mask_name
       type ( c_ptr ), value                   :: yaxt_exchanger_name
       integer ( kind=c_int ), value           :: use_raw_exchange
     end subroutine yac_cdef_couple__c

  end interface

  character ( len=* ), intent(in)           :: src_comp_name
  character ( len=* ), intent(in)           :: src_grid_name
  character ( len=* ), intent(in)           :: src_field_name
  character ( len=* ), intent(in)           :: tgt_comp_name
  character ( len=* ), intent(in)           :: tgt_grid_name
  character ( len=* ), intent(in)           :: tgt_field_name
  character ( len=* ), intent(in)           :: coupling_timestep
  integer, intent(in)                       :: time_unit
  integer, intent(in)                       :: time_reduction
  integer, intent(in)                       :: interp_stack_config_id
  integer, intent(in), optional             :: src_lag
  integer, intent(in), optional             :: tgt_lag
  character ( len=* ), intent(in), optional :: weight_file
  integer, intent(in), optional             :: weight_file_on_existing
  integer, intent(in), optional             :: mapping_side
  double precision, intent(in), optional    :: scale_factor
  double precision, intent(in), optional    :: scale_summand
  type(yac_string), intent(in), optional    :: src_mask_names(:)
  character ( len=* ), intent(in), optional :: tgt_mask_name
  character ( len=* ), intent(in), optional :: yaxt_exchanger_name
  logical, intent(in), optional             :: use_raw_exchange

  integer :: i, j
  integer :: src_lag_cpy, tgt_lag_cpy, weight_file_on_existing_cpy, &
             mapping_side_cpy
  character(kind=c_char), target :: weight_file_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: weight_file_ptr
  double precision :: scale_factor_cpy, scale_summand_cpy
  integer :: num_src_mask_names
  character(kind=c_char), allocatable, target :: src_mask_names_cpy(:,:)
  type(c_ptr), allocatable :: src_mask_names_ptr(:)
  character(kind=c_char), target :: tgt_mask_name_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: tgt_mask_name_ptr
  character(kind=c_char), target :: yaxt_exchanger_name_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: yaxt_exchanger_name_ptr
  integer(kind=c_int) :: use_raw_exchange_cpy

  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", src_comp_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", src_grid_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", src_field_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", tgt_comp_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", tgt_grid_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", tgt_field_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple", coupling_timestep )
  if ( present(src_lag) ) then
    src_lag_cpy = src_lag
  else
    src_lag_cpy = 0
  end if
  if ( present(tgt_lag) ) then
    tgt_lag_cpy = tgt_lag
  else
    tgt_lag_cpy = 0
  end if
  if ( present(weight_file) ) then
     YAC_CHECK_STRING_LEN ( "yac_fdef_couple", weight_file )
     weight_file_cpy = c_null_char
     do i = 1, len_trim(weight_file)
       weight_file_cpy(i) = weight_file(i:i)
     end do
     weight_file_ptr = c_loc(weight_file_cpy(1))
  else
     weight_file_ptr = c_null_ptr
  end if
  if ( present(weight_file_on_existing) ) then
    weight_file_on_existing_cpy = weight_file_on_existing
  else
    weight_file_on_existing_cpy = YAC_WGT_ON_EXISTING_OVERWRITE
  end if
  if ( present(mapping_side) ) then
    mapping_side_cpy = mapping_side
  else
    mapping_side_cpy = 1
  end if
  if ( present(scale_factor) ) then
    scale_factor_cpy = scale_factor
  else
    scale_factor_cpy = 1.0
  end if
  if ( present(scale_summand) ) then
    scale_summand_cpy = scale_summand
  else
    scale_summand_cpy = 0.0
  end if
  if ( present(src_mask_names) ) then
    num_src_mask_names = size(src_mask_names)
    allocate(src_mask_names_ptr(num_src_mask_names))
    allocate(src_mask_names_cpy(YAC_MAX_CHARLEN+1,num_src_mask_names))
    src_mask_names_cpy = c_null_char
    do i = 1, num_src_mask_names
      YAC_FASSERT(allocated(src_mask_names(i)%string), "ERROR(yac_fdef_couple): source mask name not allocated")
      YAC_CHECK_STRING_LEN ( "yac_fdef_couple", src_mask_names(i)%string )
      do j = 1, len_trim(src_mask_names(i)%string)
        src_mask_names_cpy(j, i) = src_mask_names(i)%string(j:j)
      end do
      src_mask_names_ptr(i) = c_loc(src_mask_names_cpy(1,i))
    end do
  else
    num_src_mask_names = 0
    allocate(src_mask_names_ptr(0))
  end if
  if ( present(tgt_mask_name) ) then
     YAC_CHECK_STRING_LEN ( "yac_fdef_couple", tgt_mask_name )
     tgt_mask_name_cpy = c_null_char
     do i = 1, len_trim(tgt_mask_name)
       tgt_mask_name_cpy(i) = tgt_mask_name(i:i)
     end do
     tgt_mask_name_ptr = c_loc(tgt_mask_name_cpy(1))
  else
     tgt_mask_name_ptr = c_null_ptr
  end if
  if ( present(yaxt_exchanger_name) ) then
     YAC_CHECK_STRING_LEN ( "yac_fdef_couple", yaxt_exchanger_name )
     yaxt_exchanger_name_cpy = c_null_char
     do i = 1, len_trim(yaxt_exchanger_name)
       yaxt_exchanger_name_cpy(i) = yaxt_exchanger_name(i:i)
     end do
     yaxt_exchanger_name_ptr = c_loc(yaxt_exchanger_name_cpy(1))
  else
     yaxt_exchanger_name_ptr = c_null_ptr
  end if
  if ( present(use_raw_exchange) ) then
    use_raw_exchange_cpy = MERGE(1_c_int, 0_c_int, use_raw_exchange)
  else
    use_raw_exchange_cpy = 0
  end if

  call yac_cdef_couple__c( TRIM(src_comp_name) // c_null_char,     &
                           TRIM(src_grid_name) // c_null_char,     &
                           TRIM(src_field_name) // c_null_char,    &
                           TRIM(tgt_comp_name) // c_null_char,     &
                           TRIM(tgt_grid_name) // c_null_char,     &
                           TRIM(tgt_field_name) // c_null_char,    &
                           TRIM(coupling_timestep) // c_null_char, &
                           time_unit,                              &
                           time_reduction,                         &
                           interp_stack_config_id,                 &
                           src_lag_cpy,                            &
                           tgt_lag_cpy,                            &
                           weight_file_ptr,                        &
                           weight_file_on_existing_cpy,            &
                           mapping_side_cpy,                       &
                           scale_factor_cpy,                       &
                           scale_summand_cpy,                      &
                           num_src_mask_names,                     &
                           src_mask_names_ptr,                     &
                           tgt_mask_name_ptr,                      &
                           yaxt_exchanger_name_ptr,                &
                           use_raw_exchange_cpy)

end subroutine yac_fdef_couple

subroutine yac_fdef_couple_instance( instance_id,  &
     src_comp_name, src_grid_name, src_field_name, &
     tgt_comp_name, tgt_grid_name, tgt_field_name, &
     coupling_timestep, time_unit, time_reduction, &
     interp_stack_config_id, src_lag, tgt_lag,     &
     weight_file, weight_file_on_existing,         &
     mapping_side, scale_factor, scale_summand,    &
     src_mask_names, tgt_mask_name,                &
     yaxt_exchanger_name, use_raw_exchange )

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_null_ptr, c_loc
  use yac, dummy => yac_fdef_couple_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_couple_instance__c( instance_id,             &
                                             src_comp_name,           &
                                             src_grid_name,           &
                                             src_field_name,          &
                                             tgt_comp_name,           &
                                             tgt_grid_name,           &
                                             tgt_field_name,          &
                                             coupling_timestep,       &
                                             time_unit,               &
                                             time_reduction,          &
                                             interp_stack_config_id,  &
                                             src_lag,                 &
                                             tgt_lag,                 &
                                             weight_file,             &
                                             weight_file_on_existing, &
                                             mapping_side,            &
                                             scale_factor,            &
                                             scale_summand,           &
                                             num_src_mask_names,      &
                                             src_mask_names,          &
                                             tgt_mask_name,           &
                                             yaxt_exchanger_name,     &
                                             use_raw_exchange)        &
          bind ( c, name='yac_cdef_couple_instance_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr, c_double

       integer ( kind=c_int ), value           :: instance_id
       character ( kind=c_char ), dimension(*) :: src_comp_name
       character ( kind=c_char ), dimension(*) :: src_grid_name
       character ( kind=c_char ), dimension(*) :: src_field_name
       character ( kind=c_char ), dimension(*) :: tgt_comp_name
       character ( kind=c_char ), dimension(*) :: tgt_grid_name
       character ( kind=c_char ), dimension(*) :: tgt_field_name
       character ( kind=c_char ), dimension(*) :: coupling_timestep
       integer ( kind=c_int ), value           :: time_unit
       integer ( kind=c_int ), value           :: time_reduction
       integer ( kind=c_int ), value           :: interp_stack_config_id
       integer ( kind=c_int ), value           :: src_lag
       integer ( kind=c_int ), value           :: tgt_lag
       type ( c_ptr ), value                   :: weight_file
       integer ( kind=c_int ), value           :: weight_file_on_existing
       integer ( kind=c_int ), value           :: mapping_side
       real ( kind=c_double ), value           :: scale_factor
       real ( kind=c_double ), value           :: scale_summand
       integer ( kind=c_int ), value           :: num_src_mask_names
       type ( c_ptr )                          :: src_mask_names(*)
       type ( c_ptr ), value                   :: tgt_mask_name
       type ( c_ptr ), value                   :: yaxt_exchanger_name
       integer ( kind=c_int ), value           :: use_raw_exchange
     end subroutine yac_cdef_couple_instance__c

  end interface

  integer, intent(in)                       :: instance_id
  character ( len=* ), intent(in)           :: src_comp_name
  character ( len=* ), intent(in)           :: src_grid_name
  character ( len=* ), intent(in)           :: src_field_name
  character ( len=* ), intent(in)           :: tgt_comp_name
  character ( len=* ), intent(in)           :: tgt_grid_name
  character ( len=* ), intent(in)           :: tgt_field_name
  character ( len=* ), intent(in)           :: coupling_timestep
  integer, intent(in)                       :: time_unit
  integer, intent(in)                       :: time_reduction
  integer, intent(in)                       :: interp_stack_config_id
  integer, intent(in), optional             :: src_lag
  integer, intent(in), optional             :: tgt_lag
  character ( len=* ), intent(in), optional :: weight_file
  integer, intent(in), optional             :: weight_file_on_existing
  integer, intent(in), optional             :: mapping_side
  double precision, intent(in), optional    :: scale_factor
  double precision, intent(in), optional    :: scale_summand
  type(yac_string), intent(in), optional    :: src_mask_names(:)
  character ( len=* ), intent(in), optional :: tgt_mask_name
  character ( len=* ), intent(in), optional :: yaxt_exchanger_name
  logical, intent(in), optional             :: use_raw_exchange

  integer :: i, j
  integer :: src_lag_cpy, tgt_lag_cpy, weight_file_on_existing_cpy, &
             mapping_side_cpy
  character(kind=c_char), target :: weight_file_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: weight_file_ptr
  double precision :: scale_factor_cpy, scale_summand_cpy
  integer :: num_src_mask_names
  character(kind=c_char), allocatable, target :: src_mask_names_cpy(:,:)
  type(c_ptr), allocatable :: src_mask_names_ptr(:)
  character(kind=c_char), target :: tgt_mask_name_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: tgt_mask_name_ptr
  character(kind=c_char), target :: yaxt_exchanger_name_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: yaxt_exchanger_name_ptr
  integer(kind=c_int) :: use_raw_exchange_cpy

  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", src_comp_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", src_grid_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", src_field_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", tgt_comp_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", tgt_grid_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", tgt_field_name )
  YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", coupling_timestep )
  if ( present(src_lag) ) then
    src_lag_cpy = src_lag
  else
    src_lag_cpy = 0
  end if
  if ( present(tgt_lag) ) then
    tgt_lag_cpy = tgt_lag
  else
    tgt_lag_cpy = 0
  end if
  if ( present(weight_file) ) then
     YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", weight_file )
     weight_file_cpy = c_null_char
     do i = 1, len_trim(weight_file)
       weight_file_cpy(i) = weight_file(i:i)
     end do
     weight_file_ptr = c_loc(weight_file_cpy(1))
  else
     weight_file_ptr = c_null_ptr
  end if
  if ( present(weight_file_on_existing) ) then
    weight_file_on_existing_cpy = weight_file_on_existing
  else
    weight_file_on_existing_cpy = YAC_WGT_ON_EXISTING_OVERWRITE
  end if
  if ( present(mapping_side) ) then
    mapping_side_cpy = mapping_side
  else
    mapping_side_cpy = 1
  end if
  if ( present(scale_factor) ) then
    scale_factor_cpy = scale_factor
  else
    scale_factor_cpy = 1.0
  end if
  if ( present(scale_summand) ) then
    scale_summand_cpy = scale_summand
  else
    scale_summand_cpy = 0.0
  end if
  if ( present(src_mask_names) ) then
    num_src_mask_names = size(src_mask_names)
    allocate(src_mask_names_ptr(num_src_mask_names))
    allocate(src_mask_names_cpy(YAC_MAX_CHARLEN+1,num_src_mask_names))
    src_mask_names_cpy = c_null_char
    do i = 1, num_src_mask_names
      YAC_FASSERT(allocated(src_mask_names(i)%string), "ERROR(yac_fdef_couple): source mask name not allocated")
      YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", src_mask_names(i)%string )
      do j = 1, len_trim(src_mask_names(i)%string)
        src_mask_names_cpy(j, i) = src_mask_names(i)%string(j:j)
      end do
      src_mask_names_ptr(i) = c_loc(src_mask_names_cpy(1,i))
    end do
  else
    num_src_mask_names = 0
    allocate(src_mask_names_ptr(0))
  end if
  if ( present(tgt_mask_name) ) then
     YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", tgt_mask_name )
     tgt_mask_name_cpy = c_null_char
     do i = 1, len_trim(tgt_mask_name)
       tgt_mask_name_cpy(i) = tgt_mask_name(i:i)
     end do
     tgt_mask_name_ptr = c_loc(tgt_mask_name_cpy(1))
  else
     tgt_mask_name_ptr = c_null_ptr
  end if
  if ( present(yaxt_exchanger_name) ) then
     YAC_CHECK_STRING_LEN ( "yac_fdef_couple_instance", yaxt_exchanger_name )
     yaxt_exchanger_name_cpy = c_null_char
     do i = 1, len_trim(yaxt_exchanger_name)
       yaxt_exchanger_name_cpy(i) = yaxt_exchanger_name(i:i)
     end do
     yaxt_exchanger_name_ptr = c_loc(yaxt_exchanger_name_cpy(1))
  else
     yaxt_exchanger_name_ptr = c_null_ptr
  end if
  if ( present(use_raw_exchange) ) then
    use_raw_exchange_cpy = MERGE(1_c_int, 0_c_int, use_raw_exchange)
  else
    use_raw_exchange_cpy = 0
  end if

  call yac_cdef_couple_instance__c( instance_id, &
                                    TRIM(src_comp_name) // c_null_char,     &
                                    TRIM(src_grid_name) // c_null_char,     &
                                    TRIM(src_field_name) // c_null_char,    &
                                    TRIM(tgt_comp_name) // c_null_char,     &
                                    TRIM(tgt_grid_name) // c_null_char,     &
                                    TRIM(tgt_field_name) // c_null_char,    &
                                    TRIM(coupling_timestep) // c_null_char, &
                                    time_unit,                              &
                                    time_reduction,                         &
                                    interp_stack_config_id,                 &
                                    src_lag_cpy,                            &
                                    tgt_lag_cpy,                            &
                                    weight_file_ptr,                        &
                                    weight_file_on_existing_cpy,            &
                                    mapping_side_cpy,                       &
                                    scale_factor_cpy,                       &
                                    scale_summand_cpy,                      &
                                    num_src_mask_names,                     &
                                    src_mask_names_ptr,                     &
                                    tgt_mask_name_ptr,                      &
                                    yaxt_exchanger_name_ptr,                &
                                    use_raw_exchange_cpy )

end subroutine yac_fdef_couple_instance

! ---------------------------------------------------------------------

function yac_fget_component_name_from_field_id ( field_id ) &
     result(comp_name)

   use yac, dummy => yac_fget_component_name_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_component_name_from_field_id_c ( field_id ) &
           result(comp_name)                                        &
       bind ( c, name='yac_cget_component_name_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id !< [IN]  field ID
       TYPE(c_ptr) :: comp_name

     end function yac_cget_component_name_from_field_id_c

   end interface

   integer, intent (in)           :: field_id  !< [IN]  field identifier
   character (len=:), allocatable :: comp_name !< [OUT] field name

   comp_name = yac_internal_cptr2char( &
        yac_cget_component_name_from_field_id_c ( field_id ))

 end function yac_fget_component_name_from_field_id

 ! ---------------------------------------------------------------------

function yac_fget_grid_name_from_field_id ( field_id ) &
     result(grid_name)

   use yac, dummy => yac_fget_grid_name_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_grid_name_from_field_id_c ( field_id ) &
           result(grid_name)                                   &
       bind ( c, name='yac_cget_grid_name_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value  :: field_id !< [IN]  field ID
       TYPE(c_ptr) :: grid_name

     end function yac_cget_grid_name_from_field_id_c

   end interface

   integer, intent (in)           :: field_id  !< [IN]  field identifier
   character (len=:), ALLOCATABLE :: grid_name !< [OUT] field name

   grid_name = yac_internal_cptr2char( &
        yac_cget_grid_name_from_field_id_c ( field_id ))

 end function yac_fget_grid_name_from_field_id

! ---------------------------------------------------------------------

function yac_fget_field_name_from_field_id ( field_id ) &
     result(field_name)

   use yac, dummy => yac_fget_field_name_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_field_name_from_field_id_c ( field_id ) &
           result(field_name)                                   &
       bind ( c, name='yac_cget_field_name_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id !< [IN]  field ID
       TYPE(c_ptr)                   :: field_name

     end function yac_cget_field_name_from_field_id_c

   end interface

   integer, intent (in)           :: field_id   !< [IN]  field identifier
   character (len=:), ALLOCATABLE :: field_name !< [OUT] field name

   field_name = yac_internal_cptr2char( &
        yac_cget_field_name_from_field_id_c ( field_id ))

 end function yac_fget_field_name_from_field_id

! ---------------------------------------------------------------------

 function yac_fget_role_from_field_id ( field_id )

   use yac, dummy => yac_fget_role_from_field_id

   implicit none

   interface

     function yac_cget_role_from_field_id_c ( field_id ) &
       bind ( c, name='yac_cget_role_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id   !< [IN]  field ID
       integer ( kind=c_int )        :: yac_cget_role_from_field_id_c

     end function  yac_cget_role_from_field_id_c

   end interface

   integer, intent (in) :: field_id   !< [IN]  field identifier
   integer              :: yac_fget_role_from_field_id

   yac_fget_role_from_field_id = &
     yac_cget_role_from_field_id_c ( field_id )

 end function yac_fget_role_from_field_id

 ! Note that in contrast to most of the other functions in this file, we have to
 ! introduce a separate result variable below. Otherwise, NVHPC fails to compile
 ! the file due to the fact that the name of the function is also a name of the
 ! interface in the YAC module.
 function yac_fget_field_role ( comp_name, grid_name, field_name ) result( res )

   use yac, dummy => yac_fget_field_role
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_role_c ( comp_name, grid_name, field_name ) &
       bind ( c, name='yac_cget_field_role' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_role_c

     end function  yac_cget_field_role_c

   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer                      :: res

   res = yac_cget_field_role_c (     &
     TRIM(comp_name) // c_null_char, &
     TRIM(grid_name) // c_null_char, &
     TRIM(field_name) // c_null_char )

 end function yac_fget_field_role

 function yac_fget_field_role_instance ( yac_instance_id, comp_name, grid_name, field_name )

   use yac, dummy => yac_fget_field_role_instance
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_role_instance_c ( yac_instance_id, &
                                               comp_name,       &
                                               grid_name,       &
                                               field_name )     &
       bind ( c, name='yac_cget_field_role_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_role_instance_c

     end function  yac_cget_field_role_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer :: yac_fget_field_role_instance

   yac_fget_field_role_instance =       &
     yac_cget_field_role_instance_c (   &
        yac_instance_id,                &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char )

 end function yac_fget_field_role_instance

! ---------------------------------------------------------------------

subroutine yac_fget_field_source ( tgt_comp_name,   &
                                  tgt_grid_name,   &
                                  tgt_field_name,  &
                                  src_comp_name,   &
                                  src_grid_name,   &
                                  src_field_name )

  use yac, dummy => yac_fget_field_source
  use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cget_field_source_c ( tgt_comp_name,  &
                                        tgt_grid_name,  &
                                        tgt_field_name, &
                                        src_comp_name,  &
                                        src_grid_name,  &
                                        src_field_name) &
       bind ( c, name='yac_cget_field_source' )

       use, intrinsic :: iso_c_binding, only : c_ptr, c_char

       character( kind=c_char), dimension(*), intent(in)  :: tgt_comp_name
       character( kind=c_char), dimension(*), intent(in)  :: tgt_grid_name
       character( kind=c_char), dimension(*), intent(in)  :: tgt_field_name
       type(c_ptr)  :: src_comp_name
       type(c_ptr)  :: src_grid_name
       type(c_ptr)  :: src_field_name

     end subroutine  yac_cget_field_source_c

  end interface

  character(len=*), intent(in) :: tgt_comp_name
  character(len=*), intent(in) :: tgt_grid_name
  character(len=*), intent(in) :: tgt_field_name
  character(len=:), ALLOCATABLE :: src_comp_name
  character(len=:), ALLOCATABLE :: src_grid_name
  character(len=:), ALLOCATABLE :: src_field_name

  type(c_ptr)  :: src_comp_name_c
  type(c_ptr)  :: src_grid_name_c
  type(c_ptr)  :: src_field_name_c

  CALL yac_cget_field_source_c ( TRIM(tgt_comp_name) // c_null_char,  &
                                 TRIM(tgt_grid_name) // c_null_char,  &
                                 TRIM(tgt_field_name) // c_null_char, &
                                 src_comp_name_c,                     &
                                 src_grid_name_c,                     &
                                 src_field_name_c)

  src_comp_name = yac_internal_cptr2char(src_comp_name_c)
  src_grid_name = yac_internal_cptr2char(src_grid_name_c)
  src_field_name = yac_internal_cptr2char(src_field_name_c)

end subroutine yac_fget_field_source

subroutine yac_fget_field_source_instance ( yac_instance_id, &
                                            tgt_comp_name,   &
                                            tgt_grid_name,   &
                                            tgt_field_name,  &
                                            src_comp_name,   &
                                            src_grid_name,   &
                                            src_field_name )

  use yac, dummy => yac_fget_field_source_instance
  use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cget_field_source_instance_c ( yac_instance_id, &
                                                   tgt_comp_name,  &
                                                   tgt_grid_name,  &
                                                   tgt_field_name, &
                                                   src_comp_name,  &
                                                   src_grid_name,  &
                                                   src_field_name) &
       bind ( c, name='yac_cget_field_source_instance' )

       use, intrinsic :: iso_c_binding, only : c_ptr, c_char, c_int

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: tgt_comp_name
       character( kind=c_char), dimension(*), intent(in)  :: tgt_grid_name
       character( kind=c_char), dimension(*), intent(in)  :: tgt_field_name
       type(c_ptr)  :: src_comp_name
       type(c_ptr)  :: src_grid_name
       type(c_ptr)  :: src_field_name

     end subroutine  yac_cget_field_source_instance_c

  end interface

  integer, intent(in) :: yac_instance_id
  character(len=*), intent(in) :: tgt_comp_name
  character(len=*), intent(in) :: tgt_grid_name
  character(len=*), intent(in) :: tgt_field_name
  character(len=:), ALLOCATABLE :: src_comp_name
  character(len=:), ALLOCATABLE :: src_grid_name
  character(len=:), ALLOCATABLE :: src_field_name

  type(c_ptr)  :: src_comp_name_c
  type(c_ptr)  :: src_grid_name_c
  type(c_ptr)  :: src_field_name_c

  CALL yac_cget_field_source_instance_c ( yac_instance_id,                     &
                                          TRIM(tgt_comp_name) // c_null_char,  &
                                          TRIM(tgt_grid_name) // c_null_char,  &
                                          TRIM(tgt_field_name) // c_null_char, &
                                          src_comp_name_c,                     &
                                          src_grid_name_c,                     &
                                          src_field_name_c)

  src_comp_name = yac_internal_cptr2char(src_comp_name_c)
  src_grid_name = yac_internal_cptr2char(src_grid_name_c)
  src_field_name = yac_internal_cptr2char(src_field_name_c)

end subroutine yac_fget_field_source_instance

! ---------------------------------------------------------------------

function yac_fget_timestep_from_field_id ( field_id ) result(string)

   use yac, dummy => yac_fget_timestep_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_timestep_from_field_id_c ( field_id ) &
           result(string)                                     &
           bind ( c, name='yac_cget_timestep_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id
       type(c_ptr)                   :: string

     end function  yac_cget_timestep_from_field_id_c

    end interface

    integer, intent (in)           :: field_id !< [IN]  field identifier
    character (len=:), ALLOCATABLE :: string   !< [OUT] timestep in iso format

    string =                  &
      yac_internal_cptr2char( &
        yac_cget_timestep_from_field_id_c(field_id))

  end function yac_fget_timestep_from_field_id

  function yac_fget_field_timestep ( comp_name, grid_name, field_name ) &
       result( timestep )

    use yac, dummy => yac_fget_field_timestep
    use mo_yac_iso_c_helpers
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

   implicit none

   interface

      function yac_cget_field_timestep_c ( comp_name, grid_name, field_name) &
           result( timestep ) &
           bind ( c, name='yac_cget_field_timestep' )

        use, intrinsic :: iso_c_binding, only : c_char, c_ptr

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       type(c_ptr)  :: timestep
     end function  yac_cget_field_timestep_c

   end interface

   character(len=*), intent(in)  :: comp_name
   character(len=*), intent(in)  :: grid_name
   character(len=*), intent(in)  :: field_name
   character(len=:), ALLOCATABLE :: timestep
   TYPE(c_ptr)                   :: c_char_ptr

   c_char_ptr =                         &
     yac_cget_field_timestep_c (        &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char)

   timestep = yac_internal_cptr2char(c_char_ptr)
 end function  yac_fget_field_timestep

 function yac_fget_field_timestep_instance ( yac_instance_id, comp_name, grid_name, field_name ) &
      result( timestep )

   use yac, dummy => yac_fget_field_timestep_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

   implicit none

   interface

      function yac_cget_field_timestep_instance_c ( yac_instance_id, &
                                                    comp_name,       &
                                                    grid_name,       &
                                                    field_name )     &
           result( timestep )                                        &
           bind ( c, name='yac_cget_field_timestep_instance' )

        use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr

        integer( kind=c_int ), intent(in), value :: yac_instance_id
        character( kind=c_char), dimension(*), intent(in) :: comp_name
        character( kind=c_char), dimension(*), intent(in) :: grid_name
        character( kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr)  :: timestep

     end function  yac_cget_field_timestep_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=:), ALLOCATABLE :: timestep
   type(c_ptr)                  :: c_char_ptr

   c_char_ptr =                           &
     yac_cget_field_timestep_instance_c ( &
        yac_instance_id,                  &
        TRIM(comp_name) // c_null_char,   &
        TRIM(grid_name) // c_null_char,   &
        TRIM(field_name) // c_null_char)

   timestep = yac_internal_cptr2char(c_char_ptr)

 end function   yac_fget_field_timestep_instance

! ---------------------------------------------------------------------

function yac_fget_field_datetime(field_id) &
     result(datetime)

  use yac, dummy => yac_fget_field_datetime
  use mo_yac_iso_c_helpers

  implicit none

  interface
     function yac_cget_field_datetime_c( field_id ) &
          result(datetime) &
          bind ( c, name="yac_cget_field_datetime" )
       use, intrinsic :: iso_c_binding, only : c_int, c_ptr
       integer (kind=c_int ), intent(in), value :: field_id
       type(c_ptr) :: datetime
     end function yac_cget_field_datetime_c
  end interface

  integer, intent(in) :: field_id
  character(len=:), allocatable :: datetime
  datetime = yac_internal_cptr2char(yac_cget_field_datetime_c(field_id))
end function yac_fget_field_datetime

! ---------------------------------------------------------------------

 subroutine yac_fenable_field_frac_mask( &
   comp_name, grid_name, field_name, frac_mask_fallback_value)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fenable_field_frac_mask

   implicit none

   interface
      subroutine yac_cenable_field_frac_mask_c (comp_name,                &
                                                grid_name,                &
                                                field_name,               &
                                                frac_mask_fallback_value) &
           bind( c, name="yac_cenable_field_frac_mask" )
        use, intrinsic :: iso_c_binding, only: c_char, c_double
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        real(kind=c_double), value :: frac_mask_fallback_value
      end subroutine yac_cenable_field_frac_mask_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision, intent(in) :: frac_mask_fallback_value

   CALL yac_cenable_field_frac_mask_c (  &
        TRIM(comp_name) // c_null_char,  &
        TRIM(grid_name) // c_null_char,  &
        TRIM(field_name) // c_null_char, &
        frac_mask_fallback_value)
 end subroutine yac_fenable_field_frac_mask

 subroutine yac_fenable_field_frac_mask_instance(yac_instance_id, comp_name, &
      grid_name, field_name, frac_mask_fallback_value)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fenable_field_frac_mask_instance

   implicit none

   interface
      subroutine yac_cenable_field_frac_mask_instance_c ( &
        yac_instance_id, comp_name, grid_name, field_name, &
        frac_mask_fallback_value) &
          bind( c, name="yac_cenable_field_frac_mask_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        real(kind=c_double), value :: frac_mask_fallback_value
      end subroutine yac_cenable_field_frac_mask_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision, intent(in) :: frac_mask_fallback_value

   CALL yac_cenable_field_frac_mask_instance_c (  &
        yac_instance_id,                          &
        TRIM(comp_name) // c_null_char,           &
        TRIM(grid_name) // c_null_char,           &
        TRIM(field_name) // c_null_char,          &
        frac_mask_fallback_value)
 end subroutine yac_fenable_field_frac_mask_instance

! ---------------------------------------------------------------------

 ! Note that in contrast to most of the other functions in this file, we have to
 ! introduce a separate result variable below. Otherwise, NVHPC fails to compile
 ! the file due to the fact that the name of the function is also a name of the
 ! interface in the YAC module.
 function yac_fget_field_frac_mask_fallback_value ( &
   comp_name, grid_name, field_name ) result( res )

   use yac, dummy => yac_fget_field_frac_mask_fallback_value
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_frac_mask_fallback_value_c ( &
       comp_name, grid_name, field_name ) &
       bind ( c, name='yac_cget_field_frac_mask_fallback_value' )

       use, intrinsic :: iso_c_binding, only : c_double, c_char

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       real ( kind=c_double ) :: yac_cget_field_frac_mask_fallback_value_c

     end function  yac_cget_field_frac_mask_fallback_value_c

   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision             :: res

   res = yac_cget_field_frac_mask_fallback_value_c ( &
     TRIM(comp_name) // c_null_char,                 &
     TRIM(grid_name) // c_null_char,                 &
     TRIM(field_name) // c_null_char )

 end function yac_fget_field_frac_mask_fallback_value

 function yac_fget_field_frac_mask_fallback_value_instance ( &
   yac_instance_id, comp_name, grid_name, field_name )

   use yac, dummy => yac_fget_field_frac_mask_fallback_value_instance
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_frac_mask_fallback_value_instance_c ( &
       yac_instance_id, comp_name, grid_name, field_name ) &
       bind ( c, name='yac_cget_field_frac_mask_fallback_value_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       real ( kind=c_double ) :: yac_cget_field_frac_mask_fallback_value_instance_c

     end function  yac_cget_field_frac_mask_fallback_value_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision :: yac_fget_field_frac_mask_fallback_value_instance

   yac_fget_field_frac_mask_fallback_value_instance =     &
     yac_cget_field_frac_mask_fallback_value_instance_c ( &
        yac_instance_id,                                  &
        TRIM(comp_name) // c_null_char,                   &
        TRIM(grid_name) // c_null_char,                   &
        TRIM(field_name) // c_null_char )

 end function yac_fget_field_frac_mask_fallback_value_instance

 ! ---------------------------------------------------------------------

function yac_fget_collection_size_from_field_id ( field_id ) &
     RESULT(collection_size)

   use yac, dummy => yac_fget_collection_size_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_collection_size_from_field_id_c ( field_id ) &
           RESULT(collection_size) &
           bind ( c, name='yac_cget_collection_size_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int )        :: collection_size

     end function  yac_cget_collection_size_from_field_id_c

   end interface

   integer, intent (in)  :: field_id        !< [IN]  field identifier
   integer               :: collection_size !< [OUT] collection size

   collection_size = &
     yac_cget_collection_size_from_field_id_c ( field_id )

 end function  yac_fget_collection_size_from_field_id

 ! Note that in contrast to most of the other functions in this file, we have to
 ! introduce a separate result variable below. Otherwise, NVHPC fails to compile
 ! the file due to the fact that the name of the function is also a name of the
 ! interface in the YAC module.
 function yac_fget_field_collection_size ( comp_name, grid_name, field_name ) result( res )

   use yac, dummy => yac_fget_field_collection_size
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_collection_size_c ( comp_name,  &
                                                grid_name,   &
                                                field_name ) &
       bind ( c, name='yac_cget_field_collection_size' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_collection_size_c

     end function  yac_cget_field_collection_size_c

   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer                      :: res

   res = yac_cget_field_collection_size_c ( &
     TRIM(comp_name) // c_null_char,        &
     TRIM(grid_name) // c_null_char,        &
     TRIM(field_name) // c_null_char )

 end function yac_fget_field_collection_size

 function yac_fget_field_collection_size_instance ( yac_instance_id, comp_name, grid_name, field_name )

   use yac, dummy => yac_fget_field_collection_size_instance
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_collection_size_instance_c ( yac_instance_id, &
                                                          comp_name,       &
                                                          grid_name,       &
                                                          field_name )     &
       bind ( c, name='yac_cget_field_collection_size_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_collection_size_instance_c

     end function  yac_cget_field_collection_size_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer :: yac_fget_field_collection_size_instance

   yac_fget_field_collection_size_instance =     &
     yac_cget_field_collection_size_instance_c ( &
        yac_instance_id,                         &
        TRIM(comp_name) // c_null_char,          &
        TRIM(grid_name) // c_null_char,          &
        TRIM(field_name) // c_null_char )

 end function yac_fget_field_collection_size_instance

 ! ---------------------------------------------------------------------

 subroutine yac_fdef_component_metadata(comp_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fdef_component_metadata

   implicit none

   interface
      subroutine yac_cdef_component_metadata_c (comp_name, metadata) &
           bind( c, name="yac_cdef_component_metadata" )
        use, intrinsic :: iso_c_binding, only: c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_component_metadata_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_component_metadata_c (TRIM(comp_name) // c_null_char, &
                                       TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_component_metadata

 subroutine yac_fdef_component_metadata_instance(yac_instance_id, comp_name, &
      metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fdef_component_metadata_instance

   implicit none

   interface
      subroutine yac_cdef_component_metadata_instance_c (yac_instance_id, &
                                                         comp_name,       &
                                                         metadata)        &
           bind( c, name="yac_cdef_component_metadata_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_component_metadata_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_component_metadata_instance_c ( &
        yac_instance_id,                         &
        TRIM(comp_name) // c_null_char,          &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_component_metadata_instance

 subroutine yac_fdef_grid_metadata(grid_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fdef_grid_metadata

   implicit none

   interface
      subroutine yac_cdef_grid_metadata_c (grid_name, metadata) &
           bind( c, name="yac_cdef_grid_metadata" )
        use, intrinsic :: iso_c_binding, only: c_char
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_grid_metadata_c
   end interface

   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_grid_metadata_c ( TRIM(grid_name) // c_null_char, &
                                   TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_grid_metadata

 subroutine yac_fdef_grid_metadata_instance(yac_instance_id, grid_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fdef_grid_metadata_instance

   implicit none

   interface
      subroutine yac_cdef_grid_metadata_instance_c (yac_instance_id, &
                                                    grid_name,       &
                                                    metadata)        &
           bind( c, name="yac_cdef_grid_metadata_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_grid_metadata_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_grid_metadata_instance_c ( &
        yac_instance_id,                    &
        TRIM(grid_name) // c_null_char,     &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_grid_metadata_instance

 subroutine yac_fdef_field_metadata(comp_name, grid_name, field_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fdef_field_metadata

   implicit none

   interface
      subroutine yac_cdef_field_metadata_c (comp_name,  &
                                            grid_name,  &
                                            field_name, &
                                            metadata)   &
           bind( c, name="yac_cdef_field_metadata" )
        use, intrinsic :: iso_c_binding, only: c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_field_metadata_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_field_metadata_c (      &
        TRIM(comp_name) // c_null_char,  &
        TRIM(grid_name) // c_null_char,  &
        TRIM(field_name) // c_null_char, &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_field_metadata

 subroutine yac_fdef_field_metadata_instance(yac_instance_id, comp_name, &
      grid_name, field_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use yac, dummy => yac_fdef_field_metadata_instance

   implicit none

   interface
      subroutine yac_cdef_field_metadata_instance_c (yac_instance_id, &
                                                     comp_name,       &
                                                     grid_name,       &
                                                     field_name,      &
                                                     metadata)        &
           bind( c, name="yac_cdef_field_metadata_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_field_metadata_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_field_metadata_instance_c ( &
        yac_instance_id,                     &
        TRIM(comp_name) // c_null_char,      &
        TRIM(grid_name) // c_null_char,      &
        TRIM(field_name) // c_null_char,     &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_field_metadata_instance

 function yac_fcomponent_has_metadata(comp_name) result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use yac, dummy => yac_fcomponent_has_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_c (comp_name) &
           result(metadata)                              &
           bind( c, name="yac_cget_component_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   logical :: has_metadata
   has_metadata = &
    c_associated( &
      yac_cget_component_metadata_c(TRIM(comp_name) // c_null_char))
 end function yac_fcomponent_has_metadata

 function yac_fcomponent_has_metadata_instance(yac_instance_id, comp_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use yac, dummy => yac_fcomponent_has_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_instance_c(yac_instance_id, &
                                                      comp_name)       &
           result(metadata)                                            &
           bind( c, name="yac_cget_component_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   logical :: has_metadata
   has_metadata =                              &
     c_associated(                             &
       yac_cget_component_metadata_instance_c( &
         yac_instance_id, TRIM(comp_name) // c_null_char))
 end function yac_fcomponent_has_metadata_instance

 function yac_fgrid_has_metadata(grid_name) result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use yac, dummy => yac_fgrid_has_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_c(grid_name) result(metadata) &
           bind( c, name="yac_cget_grid_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_c
   end interface
   character(len=*), intent(in) :: grid_name
   logical :: has_metadata
   has_metadata = &
    c_associated(yac_cget_grid_metadata_c(TRIM(grid_name) // c_null_char))
 end function yac_fgrid_has_metadata

 function yac_fgrid_has_metadata_instance(yac_instance_id, grid_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use yac, dummy => yac_fgrid_has_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_instance_c(yac_instance_id, grid_name) &
           result(metadata) &
           bind( c, name="yac_cget_grid_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: grid_name
   logical :: has_metadata
   has_metadata =                         &
     c_associated(                        &
       yac_cget_grid_metadata_instance_c( &
         yac_instance_id, TRIM(grid_name) // c_null_char))
 end function yac_fgrid_has_metadata_instance

 function yac_ffield_has_metadata(comp_name, grid_name, field_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use yac, dummy => yac_ffield_has_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_c(comp_name, grid_name, field_name) &
           result(metadata) &
           bind( c, name="yac_cget_field_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   logical :: has_metadata
   has_metadata =                        &
     c_associated(                       &
       yac_cget_field_metadata_c(        &
         TRIM(comp_name) // c_null_char, &
         TRIM(grid_name) // c_null_char, &
         TRIM(field_name) // c_null_char))
 end function yac_ffield_has_metadata

 function yac_ffield_has_metadata_instance(              &
      yac_instance_id, comp_name, grid_name, field_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use yac, dummy => yac_ffield_has_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_instance_c(yac_instance_id, &
                                                  comp_name,       &
                                                  grid_name,       &
                                                  field_name)      &
           result(metadata)                                        &
           bind( c, name="yac_cget_field_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   logical :: has_metadata
   has_metadata =                          &
     c_associated(                         &
       yac_cget_field_metadata_instance_c( &
         yac_instance_id,                  &
         TRIM(comp_name) // c_null_char,   &
         TRIM(grid_name) // c_null_char,   &
         TRIM(field_name) // c_null_char))
 end function yac_ffield_has_metadata_instance

 function yac_fget_component_metadata(comp_name) result( metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use yac, dummy => yac_fget_component_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_c(comp_name) result(metadata) &
           bind( c, name="yac_cget_component_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata = yac_cget_component_metadata_c(TRIM(comp_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_component_metadata): no metadata defined for component " // TRIM(comp_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_component_metadata

 function yac_fget_component_metadata_instance(yac_instance_id, comp_name) &
      result( metadata )

   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use yac, dummy => yac_fget_component_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_instance_c(yac_instance_id, &
                                                      comp_name)       &
          result(metadata)                                             &
           bind( c, name="yac_cget_component_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                              &
     yac_cget_component_metadata_instance_c( &
       yac_instance_id, TRIM(comp_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_component_metadata_instance): no metadata defined for component " // TRIM(comp_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_component_metadata_instance

 function yac_fget_grid_metadata(grid_name) result( metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use yac, dummy => yac_fget_grid_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_c(grid_name) result(metadata) &
           bind( c, name="yac_cget_grid_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_c
   end interface
   character(len=*), intent(in) :: grid_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata = yac_cget_grid_metadata_c(TRIM(grid_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_grid_metadata): no metadata defined for grid " // TRIM(grid_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_grid_metadata

 function yac_fget_grid_metadata_instance(yac_instance_id, grid_name) &
      result( metadata )

   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use yac, dummy => yac_fget_grid_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_instance_c(yac_instance_id, &
                                                 grid_name)       &
           result(metadata)                                       &
           bind( c, name="yac_cget_grid_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: grid_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                         &
     yac_cget_grid_metadata_instance_c( &
       yac_instance_id, TRIM(grid_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_grid_metadata_instance): no metadata defined for grid " // TRIM(grid_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_grid_metadata_instance

 function yac_fget_field_metadata(comp_name, grid_name, field_name) &
      result( metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use yac, dummy => yac_fget_field_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_c(comp_name, grid_name, field_name) &
           result(metadata) &
           bind( c, name="yac_cget_field_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                         &
    yac_cget_field_metadata_c(          &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_field_metadata): no metadata defined for field " // TRIM(comp_name) // "::" // TRIM(grid_name) // "::" // TRIM(field_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_field_metadata

 function yac_fget_field_metadata_instance(yac_instance_id, comp_name, &
      grid_name, field_name) result( metadata )

   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use yac, dummy => yac_fget_field_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_instance_c(yac_instance_id, &
                                                  comp_name,       &
                                                  grid_name,       &
                                                  field_name)      &
           result(metadata)                                        &
           bind( c, name="yac_cget_field_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                         &
    yac_cget_field_metadata_instance_c( &
        yac_instance_id,                &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_field_metadata_instance): no metadata defined for field " // TRIM(comp_name) // "::" // TRIM(grid_name) // "::" // TRIM(field_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_field_metadata_instance

 ! ---------------------------------------------------------------------

function yac_fget_start_datetime () result (start_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fget_start_datetime

  implicit none

  interface
     function yac_cget_start_datetime_c() &
        bind ( c, name='yac_cget_start_datetime' )

     use, intrinsic :: iso_c_binding, only : c_ptr
     type(c_ptr) :: yac_cget_start_datetime_c

     end function yac_cget_start_datetime_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: start_datetime_string

  c_string_ptr = yac_cget_start_datetime_c()
  start_datetime_string = yac_internal_cptr2char(c_string_ptr)
  CALL free_c(c_string_ptr)

end function yac_fget_start_datetime

function yac_fget_start_datetime_instance (yac_instance_id) &
  result (start_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fget_start_datetime_instance

  implicit none

  interface
     function yac_cget_start_datetime_instance_c(yac_instance_id) &
        bind ( c, name='yac_cget_start_datetime_instance' )

     use, intrinsic :: iso_c_binding, only : c_ptr, c_int
     integer ( kind=c_int ), value :: yac_instance_id
     type(c_ptr)                   :: yac_cget_start_datetime_instance_c

     end function yac_cget_start_datetime_instance_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  integer, intent(in) :: yac_instance_id !< [IN]  YAC instance identifier

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: start_datetime_string

  c_string_ptr = yac_cget_start_datetime_instance_c(yac_instance_id)
  start_datetime_string = yac_internal_cptr2char(c_string_ptr)
  CALL free_c(c_string_ptr)

end function yac_fget_start_datetime_instance

! ---------------------------------------------------------------------

function yac_fget_end_datetime () result (end_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fget_end_datetime

  implicit none

  interface
     function yac_cget_end_datetime_c() &
        bind ( c, name='yac_cget_end_datetime' )

     use, intrinsic :: iso_c_binding, only : c_ptr
     type(c_ptr) :: yac_cget_end_datetime_c

     end function yac_cget_end_datetime_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: end_datetime_string

  c_string_ptr = yac_cget_end_datetime_c()
  end_datetime_string = yac_internal_cptr2char(c_string_ptr)
  CALL free_c(c_string_ptr)

end function yac_fget_end_datetime

function yac_fget_end_datetime_instance (yac_instance_id) &
  result (end_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use yac, dummy => yac_fget_end_datetime_instance

  implicit none

  interface
     function yac_cget_end_datetime_instance_c(yac_instance_id) &
        bind ( c, name='yac_cget_end_datetime_instance' )

     use, intrinsic :: iso_c_binding, only : c_ptr, c_int
     integer ( kind=c_int ), value :: yac_instance_id
     type(c_ptr)                   :: yac_cget_end_datetime_instance_c

     end function yac_cget_end_datetime_instance_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  integer, intent(in) :: yac_instance_id !< [IN]  YAC instance identifier

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: end_datetime_string

  c_string_ptr = yac_cget_end_datetime_instance_c(yac_instance_id)
  end_datetime_string = yac_internal_cptr2char(c_string_ptr)
  CALL free_c(c_string_ptr)

end function yac_fget_end_datetime_instance

! ---------------------------------------------------------------------

subroutine yac_fget_interp_stack_config ( interp_stack_config_id )

  use yac, dummy => yac_fget_interp_stack_config

  implicit none

  interface

      subroutine yac_cget_interp_stack_config_c ( interp_stack_config_id ) &
        bind ( c, name='yac_cget_interp_stack_config' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ) :: interp_stack_config_id

      end subroutine yac_cget_interp_stack_config_c

   end interface

  integer, intent(out) :: interp_stack_config_id

  call yac_cget_interp_stack_config_c ( interp_stack_config_id )

end subroutine yac_fget_interp_stack_config

subroutine yac_ffree_interp_stack_config ( interp_stack_config_id )

  use yac, dummy => yac_ffree_interp_stack_config

  implicit none

  interface

      subroutine yac_cfree_interp_stack_config_c ( interp_stack_config_id ) &
        bind ( c, name='yac_cfree_interp_stack_config' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: interp_stack_config_id

      end subroutine yac_cfree_interp_stack_config_c

   end interface

  integer, intent(in) :: interp_stack_config_id

  call yac_cfree_interp_stack_config_c ( interp_stack_config_id )

end subroutine yac_ffree_interp_stack_config

subroutine yac_fadd_interp_stack_config_average ( interp_stack_config_id, &
     reduction_type, partial_coverage)

  use yac, dummy => yac_fadd_interp_stack_config_average

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_average_c (             &
          interp_stack_config_id, reduction_type, partial_coverage ) &
        bind ( c, name='yac_cadd_interp_stack_config_average' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: interp_stack_config_id
        integer ( kind=c_int ), value :: reduction_type
        integer ( kind=c_int ), value :: partial_coverage

      end subroutine yac_cadd_interp_stack_config_average_c

   end interface

   integer, intent(in)          :: interp_stack_config_id
   integer, intent(in)          :: reduction_type
   integer, intent(in)          :: partial_coverage

  call yac_cadd_interp_stack_config_average_c ( &
         interp_stack_config_id, reduction_type, partial_coverage )

end subroutine yac_fadd_interp_stack_config_average

subroutine yac_fadd_interp_stack_config_ncc ( interp_stack_config_id, &
     weight_type, partial_coverage)

  use yac, dummy => yac_fadd_interp_stack_config_ncc

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_ncc_c (              &
          interp_stack_config_id, weight_type, partial_coverage ) &
        bind ( c, name='yac_cadd_interp_stack_config_ncc' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: interp_stack_config_id
        integer ( kind=c_int ), value :: weight_type
        integer ( kind=c_int ), value :: partial_coverage

      end subroutine yac_cadd_interp_stack_config_ncc_c

   end interface

   integer, intent(in)          :: interp_stack_config_id
   integer, intent(in)          :: weight_type
   integer, intent(in)          :: partial_coverage

  call yac_cadd_interp_stack_config_ncc_c ( &
         interp_stack_config_id, weight_type, partial_coverage )

end subroutine yac_fadd_interp_stack_config_ncc

subroutine yac_fadd_interp_stack_config_nnn(interp_stack_config_id, &
     type, n, max_search_distance, scale)

  use yac, dummy => yac_fadd_interp_stack_config_nnn
  use, intrinsic :: iso_c_binding, only : c_size_t

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_nnn_c( &
          interp_stack_config_id, type, n, max_search_distance, scale)   &
        bind ( c, name='yac_cadd_interp_stack_config_nnn' )
        use, intrinsic :: iso_c_binding, only : c_int, c_double, c_size_t

        integer ( kind=c_int ), value    :: interp_stack_config_id
        integer ( kind=c_int ), value    :: type
        integer ( kind=c_size_t ), value :: n
        real ( kind=c_double ), value    :: max_search_distance
        real ( kind=c_double ), value    :: scale

      end subroutine yac_cadd_interp_stack_config_nnn_c

   end interface

   integer, intent(in)          :: interp_stack_config_id
   integer, intent(in)          :: type
   integer, intent(in)          :: n
   double precision, intent(in) :: max_search_distance
   double precision, intent(in) :: scale

  call yac_cadd_interp_stack_config_nnn_c( &
         interp_stack_config_id, type, int(n, c_size_t), &
         max_search_distance, scale)

end subroutine yac_fadd_interp_stack_config_nnn

subroutine yac_fadd_interp_stack_config_rbf(interp_stack_config_id, &
     n, max_search_distance, scale)

  use yac, dummy => yac_fadd_interp_stack_config_rbf
  use, intrinsic :: iso_c_binding, only : c_size_t
  use yac_core, only: YAC_INTERP_RBF_N_DEFAULT_F, &
                      YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT_F, &
                      YAC_INTERP_RBF_SCALE_DEFAULT_F

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_rbf_c( &
          interp_stack_config_id, n, max_search_distance, scale)   &
        bind ( c, name='yac_cadd_interp_stack_config_rbf' )
        use, intrinsic :: iso_c_binding, only : c_int, c_double, c_size_t

        integer ( kind=c_int ), value    :: interp_stack_config_id
        integer ( kind=c_size_t ), value :: n
        real ( kind=c_double ), value    :: max_search_distance
        real ( kind=c_double ), value    :: scale

      end subroutine yac_cadd_interp_stack_config_rbf_c

   end interface

   integer, intent(in)                    :: interp_stack_config_id
   integer, optional, intent(in)          :: n
   double precision, optional, intent(in) :: max_search_distance
   double precision, optional, intent(in) :: scale

   integer          :: n_
   double precision :: max_search_distance_
   double precision :: scale_

   if (present(n)) then
     n_ = n
   else
     n_ = YAC_INTERP_RBF_N_DEFAULT_F
   endif
   if (present(max_search_distance)) then
     max_search_distance_ = max_search_distance
   else
     max_search_distance_ = YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT_F
   endif
   if (present(scale)) then
     scale_ = scale
   else
     scale_ = YAC_INTERP_RBF_SCALE_DEFAULT_F
   endif

  call yac_cadd_interp_stack_config_rbf_c( &
         interp_stack_config_id, int(n_, c_size_t), &
         max_search_distance_, scale_)

end subroutine yac_fadd_interp_stack_config_rbf

subroutine yac_fadd_interp_stack_config_conservative( &
     interp_stack_config_id, order, enforced_conserv, &
     partial_coverage, normalization)

  use yac, dummy => yac_fadd_interp_stack_config_conservative

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_conservative_c( &
         interp_stack_config_id, order, enforced_conserv,    &
         partial_coverage, normalization)                    &
        bind ( c, name='yac_cadd_interp_stack_config_conservative' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: interp_stack_config_id
        integer ( kind=c_int ), value :: order
        integer ( kind=c_int ), value :: enforced_conserv
        integer ( kind=c_int ), value :: partial_coverage
        integer ( kind=c_int ), value :: normalization

      end subroutine yac_cadd_interp_stack_config_conservative_c

   end interface

   integer, intent(in)     :: interp_stack_config_id
   integer, intent(in)     :: order
   integer, intent(in)     :: enforced_conserv
   integer, intent(in)     :: partial_coverage
   integer, intent(in)     :: normalization

  call yac_cadd_interp_stack_config_conservative_c(       &
         interp_stack_config_id, order, enforced_conserv, &
         partial_coverage, normalization)

end subroutine yac_fadd_interp_stack_config_conservative

subroutine yac_fget_spmap_overwrite_config_bnd_circle ( &
    center_lon, center_lat, inc_angle, overwrite_config_id, &
    spread_distance, max_search_distance, weight_type)

  use yac, dummy => yac_fget_spmap_overwrite_config_bnd_circle

  implicit none

  interface

     subroutine yac_cget_spmap_overwrite_config_c(overwrite_config_id) &
          bind ( c, name='yac_cget_spmap_overwrite_config' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ) :: overwrite_config_id

     end subroutine yac_cget_spmap_overwrite_config_c

     subroutine yac_cset_spmap_overwrite_config_src_point_selection_c( &
          overwrite_config_id, center_lon, center_lat, inc_angle) &
          bind ( c, name='yac_cset_spmap_overwrite_config_src_point_selection_bnd_circle' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: overwrite_config_id
       real ( kind=c_double ), value :: center_lon
       real ( kind=c_double ), value :: center_lat
       real ( kind=c_double ), value :: inc_angle

     end subroutine yac_cset_spmap_overwrite_config_src_point_selection_c

     subroutine yac_cset_spmap_overwrite_config_spread_distance_c( &
          overwrite_config_id, spread_distance) &
          bind ( c, name='yac_cset_spmap_overwrite_config_spread_distance' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: overwrite_config_id
       real ( kind=c_double ), value :: spread_distance

     end subroutine yac_cset_spmap_overwrite_config_spread_distance_c

     subroutine yac_cset_spmap_overwrite_config_max_search_distance_c( &
          overwrite_config_id, max_search_distance) &
          bind ( c, name='yac_cset_spmap_overwrite_config_max_search_distance' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: overwrite_config_id
       real ( kind=c_double ), value :: max_search_distance

     end subroutine yac_cset_spmap_overwrite_config_max_search_distance_c

     subroutine yac_cset_spmap_overwrite_config_weight_type_c( &
          overwrite_config_id, weight_type) &
          bind ( c, name='yac_cset_spmap_overwrite_config_weight_type' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: overwrite_config_id
       integer ( kind=c_int ), value :: weight_type

     end subroutine yac_cset_spmap_overwrite_config_weight_type_c

  end interface

  double precision, intent(in)           :: center_lon
  double precision, intent(in)           :: center_lat
  double precision, intent(in)           :: inc_angle
  integer, intent(out)                   :: overwrite_config_id
  double precision, intent(in), optional :: spread_distance
  double precision, intent(in), optional :: max_search_distance
  integer, intent(in), optional          :: weight_type

  integer ( kind = c_int ) :: c_overwrite_config_id

  call yac_cget_spmap_overwrite_config_c(c_overwrite_config_id)

  call  yac_cset_spmap_overwrite_config_src_point_selection_c( &
    c_overwrite_config_id, center_lon, center_lat, inc_angle)
  if (present(spread_distance)) then
    call yac_cset_spmap_overwrite_config_spread_distance_c( &
      c_overwrite_config_id, spread_distance)
  end if
  if (present(max_search_distance)) then
    call yac_cset_spmap_overwrite_config_max_search_distance_c( &
      c_overwrite_config_id, max_search_distance)
  end if
  if (present(weight_type)) then
    call yac_cset_spmap_overwrite_config_weight_type_c( &
      c_overwrite_config_id, weight_type)
  end if

  overwrite_config_id = int(c_overwrite_config_id)

end subroutine yac_fget_spmap_overwrite_config_bnd_circle

subroutine yac_ffree_spmap_overwrite_config ( &
    overwrite_config_id )

  use yac, dummy => yac_ffree_spmap_overwrite_config

  implicit none

  interface

     subroutine yac_cfree_spmap_overwrite_config_c(overwrite_config_id) &
          bind ( c, name='yac_cfree_spmap_overwrite_config' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: overwrite_config_id

     end subroutine yac_cfree_spmap_overwrite_config_c

  end interface

  integer, intent(in) :: overwrite_config_id

  CALL yac_cfree_spmap_overwrite_config_c(overwrite_config_id)

end subroutine yac_ffree_spmap_overwrite_config

subroutine yac_fadd_interp_stack_config_spmap(interp_stack_config_id, &
     spread_distance, max_search_distance, weight_type, scale_type,   &
     src_sphere_radius, src_filename, src_varname, src_min_global_id, &
     tgt_sphere_radius, tgt_filename, tgt_varname, tgt_min_global_id, &
     overwrite_config_ids)

  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_loc
  use yac, dummy => yac_fadd_interp_stack_config_spmap

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_spmap_c(                      &
          interp_stack_config_id, spread_distance,                         &
          max_search_distance, weight_type, scale_type,                    &
          src_sphere_radius, src_filename, src_varname, src_min_global_id, &
          tgt_sphere_radius, tgt_filename, tgt_varname, tgt_min_global_id, &
          overwrite_config_ids, overwrite_config_count) &
          bind ( c, name='yac_cadd_interp_stack_config_spmap_f2c' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double, c_char, c_ptr

       integer ( kind=c_int ), value         :: interp_stack_config_id
       real ( kind=c_double ), value         :: spread_distance
       real ( kind=c_double ), value         :: max_search_distance
       integer ( kind=c_int ), value         :: weight_type
       integer ( kind=c_int ), value         :: scale_type
       real ( kind=c_double ), value         :: src_sphere_radius
       character (kind=c_char), dimension(*) :: src_filename
       character (kind=c_char), dimension(*) :: src_varname
       integer ( kind=c_int ), value         :: src_min_global_id
       real ( kind=c_double ), value         :: tgt_sphere_radius
       character (kind=c_char), dimension(*) :: tgt_filename
       character (kind=c_char), dimension(*) :: tgt_varname
       integer ( kind=c_int ), value         :: tgt_min_global_id
       type ( c_ptr ), value                 :: overwrite_config_ids
       integer ( kind=c_int ), value         :: overwrite_config_count

     end subroutine yac_cadd_interp_stack_config_spmap_c

  end interface

  integer, intent(in)          :: interp_stack_config_id
  double precision, intent(in) :: spread_distance
  double precision, intent(in) :: max_search_distance
  integer, intent(in)          :: weight_type
  integer, intent(in)          :: scale_type
  double precision, intent(in)  :: src_sphere_radius
  character (len=*), intent(in) :: src_filename
  character (len=*), intent(in) :: src_varname
  integer, intent(in)           :: src_min_global_id
  double precision, intent(in)  :: tgt_sphere_radius
  character (len=*), intent(in) :: tgt_filename
  character (len=*), intent(in) :: tgt_varname
  integer, intent(in)           :: tgt_min_global_id
  integer, optional, intent(in) :: overwrite_config_ids(:)

  integer (kind=c_int), allocatable, target :: c_overwrite_config_ids(:)
  type ( c_ptr ) :: ptr_overwrite_config_ids
  integer :: overwrite_config_count

  if (present(overwrite_config_ids)) then
    overwrite_config_count = SIZE(overwrite_config_ids)
    allocate(c_overwrite_config_ids(overwrite_config_count))
    c_overwrite_config_ids = int(overwrite_config_ids, c_int)
  else
    overwrite_config_count = 0;
  end if

  if (overwrite_config_count > 0) then
    ptr_overwrite_config_ids = c_loc(c_overwrite_config_ids(1))
  else
    ptr_overwrite_config_ids = c_null_ptr
  end if

  call yac_cadd_interp_stack_config_spmap_c(                              &
         interp_stack_config_id, spread_distance,                         &
         max_search_distance, weight_type, scale_type,                    &
         src_sphere_radius, src_filename, src_varname, src_min_global_id, &
         tgt_sphere_radius, tgt_filename, tgt_varname, tgt_min_global_id, &
         ptr_overwrite_config_ids, overwrite_config_count)

end subroutine yac_fadd_interp_stack_config_spmap

subroutine yac_fadd_interp_stack_config_hcsbb(interp_stack_config_id)

  use yac, dummy => yac_fadd_interp_stack_config_hcsbb

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_hcsbb_c( &
          interp_stack_config_id)                     &
          bind ( c, name='yac_cadd_interp_stack_config_hcsbb' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: interp_stack_config_id

     end subroutine yac_cadd_interp_stack_config_hcsbb_c

  end interface

  integer, intent(in) :: interp_stack_config_id

  call yac_cadd_interp_stack_config_hcsbb_c(interp_stack_config_id)

end subroutine yac_fadd_interp_stack_config_hcsbb

subroutine yac_fadd_interp_stack_config_user_file( &
     interp_stack_config_id, filename, on_missing_file, on_success)

  use yac, dummy => yac_fadd_interp_stack_config_user_file
  use, intrinsic :: iso_c_binding, only : c_null_char

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_user_file_2_c(              &
          interp_stack_config_id, filename, on_missing_file, on_success) &
          bind ( c, name='yac_cadd_interp_stack_config_user_file_2' )
       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer ( kind=c_int ), value         :: interp_stack_config_id
       character (kind=c_char), dimension(*) :: filename
       integer ( kind=c_int ), value         :: on_missing_file
       integer ( kind=c_int ), value         :: on_success

     end subroutine yac_cadd_interp_stack_config_user_file_2_c

  end interface

  integer, intent(in)           :: interp_stack_config_id
  character (len=*), intent(in) :: filename
  integer, optional, intent(in) :: on_missing_file
  integer, optional, intent(in) :: on_success

  integer :: on_missing_file_cpy, on_success_cpy

  if (present(on_missing_file)) then
    on_missing_file_cpy = on_missing_file
  else
    on_missing_file_cpy = YAC_FILE_MISSING_ERROR
  end if

  if (present(on_success)) then
    on_success_cpy = on_success
  else
    on_success_cpy = YAC_FILE_SUCCESS_CONT
  end if

  call yac_cadd_interp_stack_config_user_file_2_c( &
         interp_stack_config_id,                   &
         TRIM(filename) // c_null_char,            &
         on_missing_file_cpy, on_success_cpy)

end subroutine yac_fadd_interp_stack_config_user_file

subroutine yac_fadd_interp_stack_config_fixed(interp_stack_config_id, &
     val)

  use yac, dummy => yac_fadd_interp_stack_config_fixed

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_fixed_c( &
          interp_stack_config_id, val)                &
          bind ( c, name='yac_cadd_interp_stack_config_fixed' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: interp_stack_config_id
       real ( kind=c_double ), value :: val

     end subroutine yac_cadd_interp_stack_config_fixed_c

  end interface

  integer, intent(in)          :: interp_stack_config_id
  double precision, intent(in) :: val

  call yac_cadd_interp_stack_config_fixed_c( &
         interp_stack_config_id, val)

end subroutine yac_fadd_interp_stack_config_fixed

subroutine yac_fadd_interp_stack_config_creep(interp_stack_config_id, &
     creep_distance)

  use yac, dummy => yac_fadd_interp_stack_config_creep

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_creep_c( &
          interp_stack_config_id, creep_distance)     &
          bind ( c, name='yac_cadd_interp_stack_config_creep' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: interp_stack_config_id
       integer ( kind=c_int ), value :: creep_distance

     end subroutine yac_cadd_interp_stack_config_creep_c

  end interface

  integer, intent(in) :: interp_stack_config_id
  integer, intent(in) :: creep_distance

  call yac_cadd_interp_stack_config_creep_c( &
         interp_stack_config_id, creep_distance)

end subroutine yac_fadd_interp_stack_config_creep

subroutine yac_fget_interp_stack_config_from_string_yaml ( &
    interp_stack_config, interp_stack_config_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fget_interp_stack_config_from_string_yaml

  implicit none

  interface

      subroutine yac_cget_interp_stack_config_from_string_yaml_c ( &
        interp_stack_config, interp_stack_config_id ) &
        bind ( c, name='yac_cget_interp_stack_config_from_string_yaml' )
        use, intrinsic :: iso_c_binding, only : c_char, c_int

        character (kind=c_char), dimension(*) :: interp_stack_config
        integer ( kind=c_int ) :: interp_stack_config_id

      end subroutine yac_cget_interp_stack_config_from_string_yaml_c

   end interface

  character ( len=* ), intent(in) :: interp_stack_config
  integer, intent(out) :: interp_stack_config_id

  call yac_cget_interp_stack_config_from_string_yaml_c ( &
         TRIM(interp_stack_config) // c_null_char, &
         interp_stack_config_id )

end subroutine yac_fget_interp_stack_config_from_string_yaml

subroutine yac_fget_interp_stack_config_from_string_json ( &
    interp_stack_config, interp_stack_config_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use yac, dummy => yac_fget_interp_stack_config_from_string_json

  implicit none

  interface

      subroutine yac_cget_interp_stack_config_from_string_json_c ( &
        interp_stack_config, interp_stack_config_id ) &
        bind ( c, name='yac_cget_interp_stack_config_from_string_json' )
        use, intrinsic :: iso_c_binding, only : c_char, c_int

        character (kind=c_char), dimension(*) :: interp_stack_config
        integer ( kind=c_int ) :: interp_stack_config_id

      end subroutine yac_cget_interp_stack_config_from_string_json_c

   end interface

  character ( len=* ), intent(in) :: interp_stack_config
  integer, intent(out) :: interp_stack_config_id

  call yac_cget_interp_stack_config_from_string_json_c ( &
         TRIM(interp_stack_config) // c_null_char, &
         interp_stack_config_id )

end subroutine yac_fget_interp_stack_config_from_string_json
