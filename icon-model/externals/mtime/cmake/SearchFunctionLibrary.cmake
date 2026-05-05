# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

# ~~~
# search_function_library(<function>,
#                         <variable-found>,
#                         <variable-library>,
#                         [OPTIONS <options>])
# ~~~
# Searches for a library defining the <function> if it is not already available.
# Calls the standard macro check_function_exists with the <function> and
# <variable-found> arguments. First with no modifications to the value of the
# CMAKE_REQUIRED_LIBRARIES variable and then with each of the libraries listed
# in the <options> prepended to CMAKE_REQUIRED_LIBRARIES. If the result of the
# macro is positive, the internal cache <variable-library> is set to the first
# library found to contain the <function> (can be an empty string if no extra
# library is required). Otherwise, the value of the <variable-library> is not
# modified.
#
function(search_function_library function var_found var_library)
  cmake_parse_arguments(PARSE_ARGV 3 ARG "" "" "OPTIONS")

  if(DEFINED ${var_found})
    return()
  endif()

  if(NOT CMAKE_REQUIRED_QUIET)
    message(CHECK_START "Looking for ${function}")
  endif()
  set(save_CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET}")
  set(CMAKE_REQUIRED_QUIET True)
  set(save_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  if(DEFINED ${var_library})
    set(options "${${var_library}}")
  else()
    set(options ";${ARG_OPTIONS}")
  endif()

  include(CheckFunctionExists)
  foreach(option IN LISTS options)
    unset(${var_found} CACHE)
    set(CMAKE_REQUIRED_LIBRARIES "${option};${save_CMAKE_REQUIRED_LIBRARIES}")
    check_function_exists(${function} ${var_found})
    if(${var_found})
      set(${var_library}
          "${option}"
          CACHE INTERNAL "Library providing function ${function}"
      )
      break()
    endif()
  endforeach()

  if(NOT save_CMAKE_REQUIRED_QUIET)
    if(${var_found})
      message(CHECK_PASS "found")
    else()
      message(CHECK_FAIL "not found")
    endif()
  endif()
endfunction()
