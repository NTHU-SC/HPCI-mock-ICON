if(YAC_FOUND)
  add_library(yac_input_plugin MODULE yac_input_plugin.c)
  target_link_libraries(yac_input_plugin ComIn YAC::YAC MPI::MPI_C)

  add_executable(yac_input_external yac_input_external.c)
  target_link_libraries(yac_input_external YAC::YAC)

  if(BUILD_TESTING)

    comin_add_replay_test(NAME yac_input
      REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
      NUM_PROCS 3
      REFERENCE_OUTPUT yac_input_reference
    )
    comin_test_add_plugin(TEST yac_input
      NAME "yac_input_plugin"
      PLUGIN_LIBRARY "$<TARGET_FILE:yac_input_plugin>"
      OPTIONS "yac_input_var"
      COMM "yac_input"
    )
    comin_test_add_external_process(TEST yac_input
      NUM_PROCS 1
      COMMAND "$<TARGET_FILE:yac_input_external> yac_input_var"
    )

  endif(BUILD_TESTING)
else()
  message(STATUS "YAC not found yac_plugin disabled")
endif()
