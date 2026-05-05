# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

cmake_minimum_required(VERSION 3.18)

project(
  mtime
  VERSION 1.3.0
  LANGUAGES C Fortran
)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build tests" ON)

option(MTIME_ENABLE_PYTHON "Enable the Python interface" OFF)
option(MTIME_BUILD_DOCS "Build documentation" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Check whether the Fortran compiler in use is NEC:
  include(CheckCompilerIsNEC)
  check_compiler_is_nec(MTIME_Fortran_is_NEC Fortran)
  if(MTIME_Fortran_is_NEC)
    if(NOT CMAKE_Fortran_COMPILE_OPTIONS_PREPROCESS_ON)
      set(CMAKE_Fortran_COMPILE_OPTIONS_PREPROCESS_ON -fpp)
    endif()
  endif()
endif()

include(GNUInstallDirs)

add_subdirectory(src)

if(MTIME_ENABLE_PYTHON)
  # Technically, we can produce the Python interface without the interpreter but
  # we REQUIRE it for now for the consistency with the Autotools-based building:
  find_package(Python REQUIRED COMPONENTS Interpreter)
  include(GNUInstallPythonDirs)
  add_subdirectory(python)
endif()

if(MTIME_BUILD_DOCS)
  add_subdirectory(doc)
endif()

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(examples)
  add_subdirectory(test)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/config.cmake.in
  ${PROJECT_BINARY_DIR}/mtime-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/mtime/cmake
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/mtime-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES ${PROJECT_BINARY_DIR}/mtime-config.cmake
        ${PROJECT_BINARY_DIR}/mtime-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/mtime/cmake
)

export(EXPORT mtime-targets FILE ${PROJECT_BINARY_DIR}/mtime-targets.cmake)

install(EXPORT mtime-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/mtime/cmake)
