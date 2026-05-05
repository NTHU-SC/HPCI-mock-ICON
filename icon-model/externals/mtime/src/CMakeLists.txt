# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

set(include_directory ${PROJECT_SOURCE_DIR}/include)
set(public_headers
    # cmake-format: sort
    ${include_directory}/mtime_calendar.h
    ${include_directory}/mtime_calendar360day.h
    ${include_directory}/mtime_calendar365day.h
    ${include_directory}/mtime_calendarGregorian.h
    ${include_directory}/mtime_date.h
    ${include_directory}/mtime_datetime.h
    ${include_directory}/mtime_eventHandling.h
    ${include_directory}/mtime_eventList.h
    ${include_directory}/mtime_iso8601.h
    ${include_directory}/mtime_julianDay.h
    ${include_directory}/mtime_time.h
    ${include_directory}/mtime_timedelta.h
    ${include_directory}/mtime_utilities.h
    ${include_directory}/orbit.h
)

add_library(
  mtime_c OBJECT
  # cmake-format: sort
  ${public_headers}
  kepler.c
  mtime_calendar.c
  mtime_calendar360day.c
  mtime_calendar365day.c
  mtime_calendarGregorian.c
  mtime_date.c
  mtime_datetime.c
  mtime_eventHandling.c
  mtime_eventList.c
  mtime_iso8601.c
  mtime_julianDay.c
  mtime_time.c
  mtime_timedelta.c
  mtime_utilities.c
  orbit.c
  vsop87.c
)

add_library(mtime::mtime_c ALIAS mtime_c)

set_target_properties(
  mtime_c
  PROPERTIES C_STANDARD 99
             C_STANDARD_REQUIRED ON
             PUBLIC_HEADER "${public_headers}"
             EXPORT_NAME mtime::mtime_c
)

if(BUILD_SHARED_LIBS OR MTIME_ENABLE_PYTHON)
  set_target_properties(mtime_c PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(
  mtime_c
  PUBLIC
    $<BUILD_INTERFACE:$<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:${include_directory}>>
    $<INSTALL_INTERFACE:$<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>>
)

include(SearchFunctionLibrary)
search_function_library(roundf MTIME_HAVE_roundf MTIME_LIBM OPTIONS m)

if(NOT MTIME_HAVE_roundf)
  message(
    FATAL_ERROR "Could not find libraries that enable the roundf function"
  )
endif()

target_link_libraries(mtime_c PRIVATE ${MTIME_LIBM})

add_library(
  mtime
  libmtime.f90
  mtime_c_bindings.f90
  mtime_constants.f90
  mtime_error_handling.f90
)

add_library(mtime::mtime ALIAS mtime)

set(Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
set_target_properties(
  mtime
  PROPERTIES Fortran_PREPROCESS ON
             Fortran_MODULE_DIRECTORY ${Fortran_MODULE_DIRECTORY}
             EXPORT_NAME mtime::mtime
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(
  mtime
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
    $<INSTALL_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>>
)

target_link_libraries(mtime PUBLIC mtime_c)

install(TARGETS mtime_c mtime EXPORT mtime-targets)

install(DIRECTORY ${Fortran_MODULE_DIRECTORY}/ TYPE INCLUDE)
