# SPDX-FileCopyrightText: 2025 DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.18)

project(
  IconTixi
  VERSION 2.2.2
  LANGUAGES C Fortran
)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

find_package(LibXml2 REQUIRED)

include(GNUInstallDirs)

add_subdirectory(src)

export(
  EXPORT "IconTixiTargets" FILE "${PROJECT_BINARY_DIR}/IconTixiTargets.cmake"
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/Config.cmake.in"
  "${PROJECT_BINARY_DIR}/IconTixiConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/IconTixiConfigVersion.cmake"
  VERSION
    "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  COMPATIBILITY AnyNewerVersion
)

install(
  EXPORT "IconTixiTargets" DESTINATION "${CMAKE_INSTALL_LIBDIR}/IconTixi/cmake"
)

install(
  FILES "${PROJECT_BINARY_DIR}/IconTixiConfig.cmake"
        "${PROJECT_BINARY_DIR}/IconTixiConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/IconTixi/cmake"
)
