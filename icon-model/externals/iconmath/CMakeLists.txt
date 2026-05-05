# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

cmake_minimum_required(VERSION 3.18)

project(
  iconmath
  VERSION 1.0.0
  LANGUAGES Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build tests" ON)
option(BUILD_ICONMATH_INTERPOLATION "Build interpolation library" ON)
option(BUILD_ICONMATH_HORIZONTAL "Build horizontal library" ON)

option(IM_ENABLE_MIXED_PRECISION "Enable mixed precision" OFF)
option(IM_ENABLE_LOOP_EXCHANGE "Enable loop exchange" OFF)
option(IM_ENABLE_DIM_SWAP "Enable dimension swap" OFF)
option(IM_ENABLE_OPENACC "Enable OpenACC support" OFF)
option(IM_ENABLE_OPENMP "Enable OpenMP support" OFF)

# GNUInstallDirs issues a warning if CMAKE_SIZEOF_VOID_P is not defined, which
# is the case with NAG. One way to circumvent that is to enable C language for
# the project but that introduces a redundant and potentially misleading
# dependency on the C compiler. For the time being, we set CMAKE_SIZEOF_VOID_P
# to 8, which results in the lib64 being an installation subdirectory for
# libraries on some systems. To avoid any other side effects, we unset the
# variable right after running the GNUInstallDirs script.
if(NOT DEFINED CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P 8)
  set(unset_CMAKE_SIZEOF_VOID_P 1)
else()
  set(unset_CMAKE_SIZEOF_VOID_P 0)
endif()
include(GNUInstallDirs)
if(unset_CMAKE_SIZEOF_VOID_P)
  unset(CMAKE_SIZEOF_VOID_P)
endif()
unset(unset_CMAKE_SIZEOF_VOID_P)

# The contents of the INTERFACE_INCLUDE_DIRECTORIES of IMPORTED targets are
# treated as SYSTEM includes by default (see
# https://cmake.org/cmake/help/v3.0/prop_tgt/NO_SYSTEM_FROM_IMPORTED.html). This
# happens if CMake knows how to provide system includes to the compiler in use.
# This knowledge for the NVHPC compilers was add in CMake 3.22.0 (see
# https://gitlab.kitware.com/cmake/cmake/-/merge_requests/6691). Unfortunately,
# the knowledge was incomplete because NVHPC Fortran 21.7 and older do not
# support the '-isystem' flag. This was fixed only in CMake 3.27.0 (see
# https://gitlab.kitware.com/cmake/cmake/-/merge_requests/8588). The following
# circumvents the issue for older NVHPC and older CMake:
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" AND CMAKE_Fortran_COMPILER_VERSION
                                                  VERSION_LESS "21.7")
  unset(CMAKE_INCLUDE_SYSTEM_FLAG_Fortran)
endif()

# Dependencies

include("${PROJECT_SOURCE_DIR}/cmake/get_library.cmake")
# build fortran-support
get_library("fortran-support" 0f57ee7a57094dfaf6d7d9276683ca6ec6e8cf34)

include("${PROJECT_SOURCE_DIR}/cmake/check_macro.cmake")
check_macro_defined(HAS_OPENMP_MACRO _OPENMP QUIET)
if(IM_ENABLE_OPENMP)
  if (NOT HAS_OPENMP_MACRO)
    find_package(OpenMP QUIET)
    # all math-libraries need OpenMP
    if(OpenMP_Fortran_FOUND)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    else()
      message(FATAL_ERROR "Could NOT find OpenMP_Fortran.")
    endif()
  else()
    if(HAS_OPENMP_MACRO)
      message(
        FATAL_ERROR
        "The OpenMP support is disabled but the _OPENMP macro is defined")
    endif()
  endif()
endif()

check_macro_defined(HAS_OPENACC_MACRO _OPENACC QUIET)

if(IM_ENABLE_OPENACC)
  # If _OPENACC is defined, assume that the required compiler flags are already
  # provided, e.g. in CMAKE_Fortran_FLAGS:
  if(NOT HAS_OPENACC_MACRO)
    # On LUMI, we only have OpenACC_Fortran, but no OpenACC_C
    find_package(OpenACC QUIET)
    if(NOT OpenACC_Fortran_FOUND)
      message(FATAL_ERROR "Could NOT find OpenACC_Fortran.")
    endif()
  endif()
else()
  if(HAS_OPENACC_MACRO)
    message(
      FATAL_ERROR
        "The OpenACC support is disabled but the _OPENACC macro is defined")
  endif()
endif()

add_subdirectory(src)

# Allow for 'make test' even if the tests are disabled:
enable_testing()
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# Install

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO)

export(EXPORT "${PROJECT_NAME}-targets"
       FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake")

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)

install(EXPORT "${PROJECT_NAME}-targets"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")

# The following targets are only added if this is the parent project
if("${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")

  include("${PROJECT_SOURCE_DIR}/cmake/list_sources.cmake")
  # Collect source files of the project:
  list_sources(format_sources)
  # Exclude source files of the external projects:
  list(FILTER format_sources EXCLUDE REGEX "/_deps/")
  # Collect Fortran source files:
  set(fortran_format_sources "${format_sources}")
  list(FILTER fortran_format_sources INCLUDE REGEX "\.(f|F|f90|F90)$")
  # fprettify fails when formatting long mantissa numbers
  list(REMOVE_ITEM fortran_format_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/support/mo_math_constants.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/mo_lib_interpolation_scalar.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/mo_lib_interpolation_vector.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/mo_lib_intp_rbf.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/horizontal/mo_lib_gradients.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/horizontal/mo_lib_divrot.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/horizontal/mo_lib_laplace.F90"
  )

  message(STATUS "Fortran source files for formatting: ${fortran_format_sources}")
  add_custom_target(
    format
    COMMAND # We do not use a config file but specify the formatting arguments
            # on the command line because the feature requires an extra Python
            # package and does not work for certain arguments even in the latest
            # release (e.g., see
            # https://github.com/pseewald/fprettify/issues/94):
            fprettify --case 2 2 2 1 --indent 2 --whitespace 2 --strip-comments
            ${fortran_format_sources}
    COMMENT "Format Fortran codes")

endif()
