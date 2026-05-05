# link example directory to make relative references to the python files
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/../examples ${CMAKE_CURRENT_BINARY_DIR}/examples SYMBOLIC)

comin_add_replay_test(NAME simple_python_plugin
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
  REFERENCE_OUTPUT simple_python_plugin_reference)
comin_test_add_plugin(TEST simple_python_plugin
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../examples/simple_python_plugin.py --hello=world 'argument with spaces'")

if( YAC_FOUND )
  comin_add_replay_test(NAME yac_python_example
    REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
    NUM_PROCS 3
    REFERENCE_OUTPUT yac_python_example_reference
  )
  comin_test_add_plugin(TEST yac_python_example
    NAME "python_adapter"
    PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
    OPTIONS "../examples/yac_example.py")
endif()

if(YAXT_FOUND)
  comin_add_replay_test(NAME halo_sync
    REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
    NUM_PROCS 3
  )
  comin_test_add_plugin(TEST halo_sync
    NAME "python_adapter"
    PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
    OPTIONS "${CMAKE_CURRENT_SOURCE_DIR}/halo_sync_test.py")
endif()

if( NOT MPI4PY_FOUND )
  # check if a working mpi4py module exists
  execute_process(COMMAND
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${Python3_EXECUTABLE} -c
    "from mpi4py import MPI; print(MPI.COMM_WORLD.size if MPI.COMM_WORLD.rank == 0 else '')"
    TIMEOUT 5
    OUTPUT_VARIABLE mpi4py_config_test_output
    RESULT_VARIABLE mpi4py_config_test_result
    ERROR_VARIABLE mpi4py_config_test_error
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(mpi4py_config_test_result)
    message(WARNING "mpi4py is not working properly. Tests with mpi4py are disabled")
    set(MPI4PY_FOUND FALSE)
  else()
    if(mpi4py_config_test_output MATCHES "2")
      set(MPI4PY_FOUND TRUE)
    else()
      message(WARNING "mpi4py is not working properly. Tests with mpi4py are disabled")
      set(MPI4PY_FOUND FALSE)
    endif()
  endif()
endif()

if(MPI4PY_FOUND)
  file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/test_mpi.py ${CMAKE_CURRENT_BINARY_DIR}/test_mpi.py SYMBOLIC)
  comin_add_replay_test(NAME python_test_mpi
    REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
    NUM_PROCS 3
    REFERENCE_OUTPUT test_mpi_reference)
  comin_test_add_plugin(TEST python_test_mpi
    NAME "python_adapter"
    PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
    OPTIONS "../test_mpi.py"
    COMM "python_mpi_test_comm")
  comin_test_add_external_process(TEST python_test_mpi
    COMMAND "$<TARGET_FILE:mpi_communicator_test_external> python_mpi_test_comm")
endif()

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/variables.py" "${CMAKE_CURRENT_BINARY_DIR}/variables.py" SYMBOLIC)
comin_add_replay_test(NAME python_variables
  REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
  NUM_PROCS 3
  REFERENCE_OUTPUT python_variables_reference)
comin_test_add_plugin(TEST python_variables
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../variables.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/descrdata.py" "${CMAKE_CURRENT_BINARY_DIR}/descrdata.py" SYMBOLIC)
comin_add_replay_test(NAME python_descrdata
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
  REFERENCE_OUTPUT python_descrdata_reference
)
comin_test_add_plugin(TEST python_descrdata
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../descrdata.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/generate_python_api_doc.py" "${CMAKE_CURRENT_BINARY_DIR}/generate_python_api_doc.py" SYMBOLIC)
comin_add_replay_test(NAME python_generate_python_api_doc
  REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
  NUM_PROCS 3)
comin_test_add_plugin(TEST python_generate_python_api_doc
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../generate_python_api_doc.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/finish.py" "${CMAKE_CURRENT_BINARY_DIR}/finish.py" SYMBOLIC)
comin_add_replay_test(NAME python_finish
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
)
comin_test_add_plugin(TEST python_finish
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../finish.py")
set_tests_properties(python_finish PROPERTIES WILL_FAIL ON)

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/test_print.py" "${CMAKE_CURRENT_BINARY_DIR}/test_print.py" SYMBOLIC)
comin_add_replay_test(NAME python_test_print
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
  REFERENCE_OUTPUT python_test_print_reference)
comin_test_add_plugin(TEST python_test_print
  NAME "test_print(debug)"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_print.py"
  LOG_DEBUG
)
comin_test_add_plugin(TEST python_test_print
  NAME "test_print(no info)"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_print.py"
  NO_LOG_INFO
)
comin_test_add_plugin(TEST python_test_print
  NAME "test_print(no info/warning)"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_print.py"
  NO_LOG_INFO
  NO_LOG_WARNING
)

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/exceptions.py" "${CMAKE_CURRENT_BINARY_DIR}/exceptions.py" SYMBOLIC)
comin_add_replay_test(NAME python_exceptions
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
  REFERENCE_OUTPUT python_exceptions_reference
)
comin_test_add_plugin(TEST python_exceptions
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../exceptions.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/uncaught_exception.py" "${CMAKE_CURRENT_BINARY_DIR}/uncaught_exception.py" SYMBOLIC)
comin_add_replay_test(NAME python_uncaught_exception
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
)
comin_test_add_plugin(TEST python_uncaught_exception
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../uncaught_exception.py")
set_tests_properties(python_uncaught_exception
  PROPERTIES
  PASS_REGULAR_EXPRESSION "python_adapter: Uncaught exception in callback")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/test_var_shapes.py" "${CMAKE_CURRENT_BINARY_DIR}/test_var_shapes.py" SYMBOLIC)
comin_add_replay_test(NAME python_test_var_shapes
  REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
  NUM_PROCS 3
)
comin_test_add_plugin(TEST python_test_var_shapes
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_var_shapes.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/test_datatypes.py" "${CMAKE_CURRENT_BINARY_DIR}/test_datatypes.py" SYMBOLIC)
comin_add_replay_test(NAME python_test_datatypes
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1
  REFERENCE_OUTPUT python_test_datatypes_reference)
comin_test_add_plugin(TEST python_test_datatypes
  NAME "test_datatypes"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_datatypes.py"
  LOG_DEBUG
)

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/gpu_test.py" "${CMAKE_CURRENT_BINARY_DIR}/gpu_test.py" SYMBOLIC)
comin_add_replay_test(NAME python_gpu_test
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1)
comin_test_add_plugin(TEST python_gpu_test
  NAME "gpu_test"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../gpu_test.py"
  LOG_DEBUG
)

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/test_sys_path.py" "${CMAKE_CURRENT_BINARY_DIR}/test_sys_path.py" SYMBOLIC)
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/dummy_module.py" "${CMAKE_CURRENT_BINARY_DIR}/dummy_module.py" SYMBOLIC)
comin_add_replay_test(NAME python_sys_path
  REPLAY_DATA "test_nwp_R02B04"
  NUM_PROCS 1)
comin_test_add_plugin(TEST python_sys_path
  NAME "sys_path"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_sys_path.py"
  LOG_DEBUG
)
comin_test_add_plugin(TEST python_sys_path
  NAME "sys_path2"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../test_sys_path.py"
  LOG_DEBUG
)
