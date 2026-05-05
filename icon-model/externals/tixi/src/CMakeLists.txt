# SPDX-FileCopyrightText: 2025 DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# SPDX-License-Identifier: BSD-3-Clause

# Function to avoid namespace pollution:
function(configure_version_file)
  set(TIXI_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
  set(TIXI_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
  set(TIXI_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
  set(TIXI_VERSION "${PROJECT_VERSION}")
  configure_file(tixi_version.h.in tixi_version.h @ONLY)
endfunction()

configure_version_file()

add_library(
  icon_tixi
  ${CMAKE_CURRENT_BINARY_DIR}/tixi_version.h
  tixi.f90
  tixi.h
  tixiData.h
  tixiImpl.c
  tixiInternal.c
  tixiInternal.h
  tixiUtils.c
  tixiUtils.h
  uidHelper.c
  uidHelper.h
  xpathFunctions.c
  xpathFunctions.h)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_library(IconTixi::icon_tixi ALIAS icon_tixi)

set_target_properties(
  icon_tixi
  PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
             EXPORT_NAME IconTixi::icon_tixi
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

target_compile_features(icon_tixi PRIVATE c_std_99)

target_include_directories(
  icon_tixi
  INTERFACE
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
    $<INSTALL_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_INSTALL_INCLUDEDIR}>>
  PRIVATE
    $<BUILD_INTERFACE:$<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:${CMAKE_CURRENT_BINARY_DIR}>>
)

target_link_libraries(icon_tixi PRIVATE LibXml2::LibXml2)

install(TARGETS icon_tixi EXPORT IconTixiTargets)

install(DIRECTORY "${Fortran_MODULE_DIRECTORY}/" TYPE INCLUDE)
