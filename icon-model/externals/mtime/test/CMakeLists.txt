# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

find_package(Check REQUIRED)

add_library(
  mtime_check STATIC
  mtime_calendar_test.c
  mtime_calendar_test.h
  mtime_date_test.c
  mtime_date_test.h
  mtime_datetime_test.c
  mtime_datetime_test.h
  mtime_julianDay_test.c
  mtime_julianDay_test.h
  mtime_time_test.c
  mtime_time_test.h
  mtime_timedelta_test.c
  mtime_timedelta_test.h
)
target_link_libraries(mtime_check PUBLIC mtime::mtime Check::Interface)

add_executable(test_runner test_runner.c)
target_compile_definitions(test_runner PRIVATE ENABLE_CHECK)
target_link_libraries(
  test_runner PRIVATE mtime_check mtime::mtime Check::Interface
)
add_test(NAME test_runner COMMAND test_runner)

if(MTIME_ENABLE_PYTHON)
  find_package(UnixCommands)
  if(NOT BASH)
    message(
      FATAL_ERROR
        "Could not find the Bash interpreter required for the Python tests"
    )
  endif()
  include(AutoconfConfigure)
  set(test_python_tmp
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/mtime/test_python.sh"
  )
  autoconf_configure(
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/test_python.sh.in
    OUTPUTS ${test_python_tmp}
    VARIABLES SHELL ENABLE_PYTHON_FALSE abs_top_builddir PYTHON
    VALUES "${BASH}" "#" "${PROJECT_BINARY_DIR}" "${Python_EXECUTABLE}"
  )
  file(
    COPY ${test_python_tmp}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FILE_PERMISSIONS
      OWNER_READ
      OWNER_WRITE
      OWNER_EXECUTE
      GROUP_READ
      GROUP_EXECUTE
      WORLD_READ
      WORLD_EXECUTE
  )
  add_test(NAME test_python.sh COMMAND test_python.sh)
endif()
