if(BUILD_TESTING)
  find_package(MPI)
  if(MPI_FOUND)
    comin_add_replay_test(NAME dummy_plugin_test
      REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
      NUM_PROCS 3
    )

    comin_test_add_plugin(TEST dummy_plugin_test
      NAME "dummy"
      PLUGIN_LIBRARY "$<TARGET_FILE:dummy_plugin>"
    )
  endif()
endif()
