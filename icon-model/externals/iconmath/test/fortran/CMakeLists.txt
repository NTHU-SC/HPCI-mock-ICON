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

message(CHECK_START "Fetching external FortUTF")
if("${CMAKE_VERSION}" VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW)
endif()
include(FetchContent)
FetchContent_Declare(
  fortutf
  URL https://github.com/artemis-beta/FortUTF/archive/ab4ceadf0ab07648ef60b71ae721cd1fe2955521.tar.gz
  URL_HASH MD5=6526aa1498e601131b727f1cfd6613ac)
FetchContent_MakeAvailable(fortutf)
message(CHECK_PASS "done")

set(FORTUTF_PROJECT_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(FORTUTF_PROJECT_SRC_LIBRARY iconmath-support)
set(FORTUTF_PROJECT_SRC_FILES helpers.f90)

include("${fortutf_SOURCE_DIR}/cmake/fortutf.cmake")

set(EXTRA_LIBRARIES fortran-support::fortran-support)

# We have to prevent FortUTF from generating invalid Fortran code:
set(name_holder ${PROJECT_NAME})
string(REPLACE "-" "_" PROJECT_NAME ${PROJECT_NAME})

FortUTF_Find_Tests(
  LIBRARIES ${EXTRA_LIBRARIES}
  MODULE_DIRS ${fortutf_BINARY_DIR}/modules)

set(PROJECT_NAME ${name_holder})
add_test(NAME FortUTF_UnitTest COMMAND iconmath_Tests)
