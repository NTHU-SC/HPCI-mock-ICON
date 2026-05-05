# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

set(examples_Fortran
    # cmake-format: sort
    callback_test
    comp_weights
    duration
    example
    iconatm
    iconoce
    output_control
    repetitor
    simulate_iau
    tas
    test_cf_timeaxis
    test_dace
    test_jd_logic
    uniq
)

foreach(example IN LISTS examples_Fortran)
  add_executable(${example} ${example}.f90)
  add_test(NAME ${example} COMMAND ${example})
endforeach()

set(examples_C
    # cmake-format: sort
    model_integration
    modulo
    recurrence
)

foreach(example IN LISTS examples_C)
  add_executable(${example} ${example}.c)
  add_test(NAME ${example} COMMAND ${example})
endforeach()

# Additional sources of the example executables:
target_sources(duration PRIVATE mo_kind.f90)
# Avoid multiple rules for the Fortran module file:
add_library(event_manager OBJECT mo_event_manager.f90)
target_sources(iconatm PRIVATE $<TARGET_OBJECTS:event_manager>)
target_sources(iconoce PRIVATE $<TARGET_OBJECTS:event_manager>)
target_sources(test_dace PRIVATE $<TARGET_OBJECTS:event_manager>)

foreach(
  target IN
  LISTS examples_Fortran
  ITEMS event_manager
)
  set_target_properties(${target} PROPERTIES Fortran_PREPROCESS ON)
endforeach()

foreach(target IN LISTS examples_C)
  set_target_properties(
    ${target} PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON
  )
endforeach()

foreach(
  target IN
  LISTS examples_Fortran examples_C
  ITEMS event_manager
)
  target_link_libraries(${target} PRIVATE mtime::mtime)
endforeach()

foreach(namelist IN ITEMS iconatm.nml iconoce.nml)
  file(
    CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/${namelist}
    ${CMAKE_CURRENT_BINARY_DIR}/${namelist}
    RESULT
      dummy # just don't fail if something goes wrong (e.g. building in-source)
    COPY_ON_ERROR SYMBOLIC
  )
endforeach()
