# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

if(NOT TARGET Check::Interface)
  include(FindPackageHandleStandardArgs)
  find_package(Check CONFIG QUIET)
  if(Check_FOUND AND TARGET Check::checkShared)
    add_library(Check::Interface ALIAS Check::checkShared)
    find_package_handle_standard_args(Check CONFIG_MODE)
  else()
    find_library(Check_LIBRARY NAMES check)
    mark_as_advanced(Check_LIBRARY)
    find_path(Check_INCLUDE_DIR NAMES check.h)
    mark_as_advanced(Check_INCLUDE_DIR)
    find_package_handle_standard_args(
      Check REQUIRED_VARS Check_LIBRARY Check_INCLUDE_DIR
    )
    if(Check_FOUND)
      add_library(Check::Interface UNKNOWN IMPORTED)
      set_target_properties(
        Check::Interface
        PROPERTIES IMPORTED_LOCATION ${Check_LIBRARY}
                   INTERFACE_INCLUDE_DIRECTORIES ${Check_INCLUDE_DIR}
      )
    endif()
  endif()
endif()
