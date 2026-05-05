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

# cmake-format: off
# fs_add_c_test(<test_name>
#               [SOURCES <sources>]
#               [ARGS <args>])
# cmake-format: on
# -----------------------------------------------------------------------------
# Compiles a test executable with the name <test_name> using the source code
# <source>. Specify ctest arguments in <args> if necessary. The googletest and
# libfortran-support libraries will be linked automatically.
#
# The C++ standard is set to C++17.
#
function(fs_add_c_test test_name)

  cmake_parse_arguments(PARSE_ARGV 1 ARG "" "" "SOURCES;ARGS")
  if(NOT ARG_SOURCES)
    set(ARG_SOURCES "${ARG_UNPARSED_ARGUMENTS}")
  endif()

  add_executable("CTest_${test_name}" ${ARG_SOURCES})
  target_link_libraries(
    "CTest_${test_name}" PRIVATE fortran-support::fortran-support
                                 GTest::gtest_main)
  add_test(NAME "CTest_${test_name}" COMMAND "CTest_${test_name}" ${ARG_ARGS})
  set_property(TEST "CTest_${test_name}" PROPERTY LABELS C)
  set_target_properties("CTest_${test_name}"
                        PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

endfunction()
