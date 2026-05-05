find_package(NetCDF REQUIRED COMPONENTS C Fortran)

add_subdirectory(src)

add_library(comin_run_recorder_plugin MODULE comin_run_recorder_plugin.F90 src/comin_descrdata_save.F90)
target_link_libraries(comin_run_recorder_plugin MPI::MPI_Fortran ComIn NetCDF::NetCDF_Fortran utils)

add_executable(comin_replay
  comin_replay.F90
  src/comin_descrdata_load.F90
  src/global_index_lookup.F90
)
add_executable(${PROJECT_NAME}::comin_replay ALIAS comin_replay)

set_target_properties(comin_replay PROPERTIES
  EXPORT_NAME "${PROJECT_NAME}::comin_replay"
)
if(YAXT_FOUND)
  target_link_libraries(comin_replay YAXT::YAXT_Fortran)
  target_compile_definitions(comin_replay PRIVATE -DENABLE_YAXT)
endif()

if(YAC_FOUND)
  # We need to link YAC before NetCDF because CMake does not resolve
  # the linking order correctly as both YAC and NetCDF_Fortran depend
  # on NetCDF_C
  target_link_libraries(comin_replay YAC::YAC)
  target_compile_definitions(comin_replay PRIVATE -DENABLE_YAC)
endif()
target_link_libraries(comin_replay MPI::MPI_Fortran ComIn NetCDF::NetCDF_Fortran utils)

install(TARGETS comin_replay EXPORT "${PROJECT_NAME}Targets")

add_library(comin_var_recorder_plugin MODULE comin_var_recorder_plugin.F90)
target_link_libraries(comin_var_recorder_plugin ComIn NetCDF::NetCDF_Fortran utils)

add_library(comin_var_replay_plugin MODULE comin_var_replay_plugin.F90)
add_library(${PROJECT_NAME}::comin_var_replay_plugin ALIAS comin_var_replay_plugin)

set_target_properties(comin_var_replay_plugin PROPERTIES
  EXPORT_NAME "${PROJECT_NAME}::comin_var_replay_plugin"
)

target_link_libraries(comin_var_replay_plugin
  PUBLIC ComIn NetCDF::NetCDF_Fortran
  PRIVATE utils
)

install(
  TARGETS comin_var_replay_plugin
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

### testing ###
# check if the output of the recorder plugins can be used as replay_data

if(BUILD_TESTING)

  comin_add_replay_test(NAME record
    REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
    NUM_PROCS 3)
  comin_test_add_plugin(TEST record
    NAME "comin_run_recorder_plugin"
    PLUGIN_LIBRARY "$<TARGET_FILE:comin_run_recorder_plugin>"
  )
  comin_test_add_plugin(TEST record
    NAME comin_var_recorder_plugin
    PLUGIN_LIBRARY "$<TARGET_FILE:comin_var_recorder_plugin>"
    OPTIONS "pres,u,v,qr,qi,qv,qs,qc")

  set_tests_properties(record PROPERTIES
    FIXTURES_SETUP record
  )
  get_test_property(record WORKING_DIRECTORY RECORD_WD)

  comin_add_replay_test(NAME replay_record
    REPLAY_DATA_PATH ${RECORD_WD}
    NUM_PROCS 3)

  # TODO: Compare output of record test with the input
  # needs some netcdf comparison tool

endif(BUILD_TESTING)
