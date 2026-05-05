# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

if(NOT DEFINED Python_FOUND)
  message(
    AUTHOR_WARNING
      "The macro requires Python: you either need to find_package(Python) or "
      "make sure that the following variables are set (e.g. you can get valid "
      "values under different variable names from find_package(Python3)): "
      "Python_FOUND, Python_EXECUTABLE, "
      "Python_VERSION_MAJOR, Python_VERSION_MINOR"
  )
endif()

# ~~~
# _GNUInstallPythonDirs_get_sitedir(<plat_specific>,
#                                   <fallback>,
#                                   <variable>,
#                                   <description>)
# ~~~
# Requests ${Python_EXECUTABLE} for the subdirectory for either the general
# (when <plat_specific> is set to 0) or platform-dependent (when <plat_specific>
# is set to 1) library installation. Set the cash <variable> with the provided
# <description> to the result of the request. If the request fails, the
# <variable> is set to the <fallback> value.
#
macro(
  _GNUInstallPythonDirs_get_sitedir plat_specific fallback variable description
)
  unset(_sitedir)
  if(Python_FOUND)
    execute_process(
      COMMAND
        ${Python_EXECUTABLE} -c
        "from distutils import sysconfig as sc;print(sc.get_python_lib(${plat_specific},0,''))"
      RESULT_VARIABLE _success
      OUTPUT_VARIABLE _sitedir
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT _success EQUAL 0)
      unset(_sitedir)
    endif()
    unset(_success)
  endif()
  if(NOT _sitedir)
    set(_sitedir "${fallback}")
  endif()
  set(${variable}
      "${_sitedir}"
      CACHE PATH "${description}"
  )
  unset(_sitedir)
endmacro()

if(NOT CMAKE_INSTALL_PYTHON_PURELIBDIR)
  _gnuinstallpythondirs_get_sitedir(
    0
    "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages"
    CMAKE_INSTALL_PYTHON_PURELIBDIR
    "Python platform-independent libraries"
  )
endif()

if(NOT CMAKE_INSTALL_PYTHON_PLATLIB)
  _gnuinstallpythondirs_get_sitedir(
    1
    "${CMAKE_INSTALL_PYTHON_PURELIBDIR}"
    CMAKE_INSTALL_PYTHON_PLATLIBDIR
    "Python platform-specific libraries"
  )
endif()

foreach(dir PYTHON_PURELIBDIR PYTHON_PLATLIBDIR)
  if(COMMAND GNUInstallDirs_get_absolute_install_dir)
    GNUInstallDirs_get_absolute_install_dir(
      CMAKE_INSTALL_FULL_${dir} CMAKE_INSTALL_${dir} ${dir}
    )
  else()
    set(CMAKE_INSTALL_FULL_${dir}
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_${dir}}"
    )
  endif()
endforeach()
