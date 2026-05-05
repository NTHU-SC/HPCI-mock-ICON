cmake_minimum_required(VERSION 3.17)

project(ComIn
  VERSION 0.4.0
  LANGUAGES Fortran C CXX
  DESCRIPTION "Community Interface")

include(GNUInstallDirs)
list( APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR} )
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Scripts/")
add_subdirectory(cmake)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build Testing" OFF) # turn testing off by default as it introduces a lot of dependencies
option(COMIN_ENABLE_EXAMPLES "Enable the example plugins" ON)
option(COMIN_ENABLE_PYTHON_ADAPTER "Enable the python adapter" OFF) # introduces python as a dependency
option(COMIN_ENABLE_STATIC_LINKING_TEST "Enable the static linking test" OFF) # may yield problems in non-default environments
option(COMIN_ENABLE_REPLAY_TOOL "Enable the replay tool" OFF)
option(COMIN_ENABLE_YAC "Enable YAC for building examples and the replay tool" OFF)
option(COMIN_ENABLE_YAXT "Enable YAXT for building exmples and the replay tool" OFF)

if( BUILD_TESTING AND NOT COMIN_ENABLE_REPLAY_TOOL)
  message(STATUS "Activating COMIN_ENABLE_REPLAY_TOOL for BUILD_TESTING")
  set(COMIN_ENABLE_REPLAY_TOOL ON CACHE BOOL "Enable the replay tool" FORCE)
endif()

include(cmake/Scripts/SetDLLibs.cmake)

set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI REQUIRED COMPONENTS C Fortran)

if( COMIN_ENABLE_YAXT )
  find_package(YAXT REQUIRED)
endif()

if( COMIN_ENABLE_YAC )
  find_package(YAC REQUIRED)
endif()

# If it is a out-of-source build throw an error if comin_version.inc exists in source tree
if( NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  if( EXISTS "${CMAKE_SOURCE_DIR}/include/comin_version.inc" )
  message( FATAL_ERROR
    "Source tree is already configured. Out-of-source builds require a clean source tree. Please clean the source tree before proceeding." )
  endif()
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/comin_version.inc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/comin_version.inc" @ONLY)

set(SOURCE_FILES
  src/comin_keyval.cpp
  src/comin_ftnlist.cpp
  src/mpi_handshake.c
  src/mo_mpi_handshake.F90
  src/comin_c_utils.F90
  src/comin_keyval_ifc.F90
  src/comin_ftnlist_ifc.F90
  src/comin_metadata_types.F90
  src/comin_setup_utils.F90
  src/comin_setup.F90
  src/comin_plugin_interface.F90
  src/comin_parallel_types.F90
  src/comin_parallel.F90
  src/comin_host_interface.F90
  src/comin_errhandler_constants.F90
  src/comin_errhandler_types.F90
  src/comin_errhandler.F90
  src/comin_callback_types.F90
  src/comin_callback.F90
  src/comin_variable_types.F90
  src/comin_metadata.F90
  src/comin_variable.F90
  src/comin_descrdata_types.F90
  src/comin_descrdata.F90
  src/comin_plugin_types.F90
  src/comin_state.F90
  src/comin_setup_constants.F90
  src/comin_print.F90
  src/comin_c_utils.c
)

# Create the library
add_library(ComIn ${SOURCE_FILES})
add_library(${PROJECT_NAME}::ComIn ALIAS ComIn)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(ComIn PROPERTIES
  OUTPUT_NAME "comin"
  EXPORT_NAME "${PROJECT_NAME}::ComIn"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
  POSITION_INDEPENDENT_CODE 1
  CXX_STANDARD 17
)

target_include_directories(ComIn
PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  "$<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
INTERFACE
  $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILER_INCLUDE_DIRS}>
PRIVATE
  "$<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_CURRENT_SOURCE_DIR}/src>>"
)

target_link_libraries(ComIn PRIVATE MPI::MPI_Fortran ${COMIN_DL_LIBS})

install(TARGETS ComIn EXPORT "${PROJECT_NAME}Targets")
install(DIRECTORY
  "${Fortran_MODULE_DIRECTORY}/"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  TYPE INCLUDE
  PATTERN "*.in" EXCLUDE
)

include(CTest)
if(BUILD_TESTING)
  include(ComInAddTest)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()

if(COMIN_ENABLE_REPLAY_TOOL)
  add_subdirectory(replay_tool)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

# Make the Config available to subprojects:
list(PREPEND CMAKE_PREFIX_PATH "${PROJECT_BINARY_DIR}")

add_subdirectory(plugins)

add_subdirectory(doc)

export(
  EXPORT "${PROJECT_NAME}Targets"
  FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
)

install(
  EXPORT "${PROJECT_NAME}Targets"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake"
)

install(
  FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake"
)
