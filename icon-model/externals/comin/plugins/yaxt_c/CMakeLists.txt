if(YAXT_FOUND)
  add_library(yaxt_c_plugin MODULE yaxt_c_plugin.c)
  target_link_libraries(yaxt_c_plugin ComIn YAXT::YAXT_C MPI::MPI_C)

  add_executable(yaxt_c_external yaxt_c_external.c ${PROJECT_SOURCE_DIR}/src/mpi_handshake.c)
  target_link_libraries(yaxt_c_external YAXT::YAXT_C MPI::MPI_C)
  target_include_directories(yaxt_c_external PRIVATE  ${PROJECT_SOURCE_DIR}/include)

  if(BUILD_TESTING)

    comin_add_replay_test(NAME yaxt_c
      REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
      NUM_PROCS 3
    )
    comin_test_add_plugin(TEST yaxt_c
      NAME yaxt_c_plugin
      PLUGIN_LIBRARY "$<TARGET_FILE:yaxt_c_plugin>"
      COMM "yaxt_c"
    )
    comin_test_add_external_process(TEST yaxt_c
      COMMAND "$<TARGET_FILE:yaxt_c_external> 5126 6"
    )

  endif(BUILD_TESTING)
endif()
