if(YAXT_FOUND)
  add_library(yaxt_fortran_plugin MODULE yaxt_fortran_plugin.F90)
  target_link_libraries(yaxt_fortran_plugin ComIn YAXT::YAXT_Fortran MPI::MPI_Fortran)

  if(BUILD_TESTING)

    comin_add_replay_test(NAME yaxt_fortran_plugin
      REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
      NUM_PROCS 3
      REFERENCE_OUTPUT yaxt_fortran_plugin_reference
    )
    comin_test_add_plugin(TEST yaxt_fortran_plugin
      PLUGIN_LIBRARY "$<TARGET_FILE:yaxt_fortran_plugin>"
      COMM "yaxt_fortran_plugin_comm"
    )

  endif(BUILD_TESTING)
else()
  message(STATUS "disabling yaxt_fortran_plugin")
endif()
