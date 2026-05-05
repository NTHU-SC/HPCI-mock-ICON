! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_abort.F90
!! \test
!! This contains an example of how to use yac_abort_message.

program test_abort

  use yac, only : yac_abort_message
  implicit none

  CALL yac_abort_message("some error message", __FILE__, __LINE__)

end program test_abort

