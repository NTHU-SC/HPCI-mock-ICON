# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

set(module_dir_name mtime)
set(module_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${module_dir_name})
set(module_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/${module_dir_name})

add_library(mtime_py MODULE)

set_target_properties(
  mtime_py
  PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${module_binary_dir}
             OUTPUT_NAME __mtime
             PREFIX ""
             SUFFIX "${CMAKE_SHARED_MODULE_SUFFIX}"
)

target_link_libraries(mtime_py PRIVATE mtime_c)

include(AutoconfConfigure)
# Also create ${module_binary_dir}:
autoconf_configure(
  INPUTS ${module_source_dir}/_mtime.py.in
  OUTPUTS ${module_binary_dir}/_mtime.py
  VARIABLES acx_modext
  VALUES "${CMAKE_SHARED_MODULE_SUFFIX}"
)

file(
  CREATE_LINK ${module_source_dir}/__init__.py ${module_binary_dir}/__init__.py
  RESULT
    dummy # just don't fail if something goes wrong (e.g. building in-source)
  COPY_ON_ERROR SYMBOLIC
)

install(
  TARGETS mtime_py
  LIBRARY DESTINATION ${CMAKE_INSTALL_PYTHON_PLATLIBDIR}/${module_dir_name}
)

install(
  FILES ${module_source_dir}/__init__.py ${module_binary_dir}/_mtime.py
  DESTINATION ${CMAKE_INSTALL_PYTHON_PLATLIBDIR}/${module_dir_name}
)
