!> @file comin_print.F90
!! @brief Utility functions for printing log messages.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_print

  USE comin_state,                ONLY: state
  USE comin_c_utils,              ONLY: convert_c_string
  USE iso_c_binding,              ONLY: C_PTR

  PUBLIC :: comin_print_debug
  PUBLIC :: comin_print_info
  PUBLIC :: comin_print_warning

CONTAINS

  !> Prints a debug message if the plugin has set log_debug (disabled per default)
  !> The message will only be printed on process 0.
  !! @ingroup plugin_interface
  SUBROUTINE comin_print_debug(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg

    IF (state%lstdout .AND. state%current_plugin%log_debug) THEN
      WRITE(state%output_unit, *) "DEBUG(" // state%current_plugin%name // "): " // TRIM(msg)
    END IF
  END SUBROUTINE comin_print_debug

  ! C-wrapper
  SUBROUTINE comin_print_debug_c(msg) &
    BIND(C, name="comin_print_debug")
    TYPE(c_ptr), VALUE, INTENT(IN) :: msg

    CALL comin_print_debug(convert_c_string(msg))
  END SUBROUTINE comin_print_debug_c

  !> Prints a info message if the plugin has set log_info (enabled per default)
  !> The message will only be printed on process 0.
  !! @ingroup plugin_interface
  SUBROUTINE comin_print_info(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg

    IF (state%lstdout .AND. state%current_plugin%log_info) THEN
      WRITE(state%output_unit, *) "INFO(" // state%current_plugin%name // "): " // TRIM(msg)
    END IF
  END SUBROUTINE comin_print_info

  ! C-wrapper
  SUBROUTINE comin_print_info_c(msg) &
    BIND(C, name="comin_print_info")
    TYPE(c_ptr), VALUE, INTENT(IN) :: msg

    CALL comin_print_info(convert_c_string(msg))
  END SUBROUTINE comin_print_info_c

  !> Prints a warning if the plugin has set log_warning (enabled per default)
  !> The message will only be printed on process 0.
  !! @ingroup plugin_interface
  SUBROUTINE comin_print_warning(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg

    IF (state%lstdout .AND. state%current_plugin%log_warning) THEN
      WRITE(state%output_unit, *) "WARNING(" // state%current_plugin%name // "): " // TRIM(msg)
    END IF
  END SUBROUTINE comin_print_warning

  ! C-wrapper
  SUBROUTINE comin_print_warning_c(msg) &
    BIND(C, name="comin_print_warning")
    TYPE(c_ptr), VALUE, INTENT(IN) :: msg

    CALL comin_print_warning(convert_c_string(msg))
  END SUBROUTINE comin_print_warning_c

END MODULE comin_print
