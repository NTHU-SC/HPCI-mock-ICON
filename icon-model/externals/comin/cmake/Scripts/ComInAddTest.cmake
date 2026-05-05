if(COMMAND comin_add_test)
    return()
endif()

if(NOT TARGET update_test_references)
  add_custom_target(update_test_references
    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red --bold "CAUTION: Please check if the updated references contain only changes that are intended!")
endif()

# this is a target to collect test_data targets that are used in tests
if(NOT TARGET download_test_data)
  add_custom_target(download_test_data)
  add_test(NAME "download_test_data"
    COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target "download_test_data"
  )
  set_tests_properties("download_test_data" PROPERTIES
    FIXTURES_SETUP download_test_data)
endif()

if (TARGET ComIn::comin_replay OR COMIN_ENABLE_REPLAY_TOOL)
  # DOCUMENTATION CAN BE FOUND IN doc/cmake.md
  # PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
  function(comin_add_replay_data)
    set(options)
    set(oneValueArgs NAME URL MD5HASH)
    set(multiValueArgs)
    cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(ARGS_MD5HASH)
      set(EH " EXPECTED_MD5 ${ARGS_MD5HASH}")
    else()
      set(EH "")
    endif()
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}.cmake
      "file(DOWNLOAD ${ARGS_URL} ${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}.tar.gz SHOW_PROGRESS ${EH})\n"
    )
    add_custom_target("download_replay_data_${ARGS_NAME}"
      BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}.tar.gz"
      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}.cmake
      COMMAND ${CMAKE_COMMAND} -E tar xf "${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}.tar.gz"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}
    )
    set_target_properties("download_replay_data_${ARGS_NAME}" PROPERTIES
      DATA_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/replay_data/${ARGS_NAME}")
  endfunction()

  # DOCUMENTATION CAN BE FOUND IN doc/cmake.md
  # PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
  function(comin_add_replay_test)
    set(options)
    set(oneValueArgs NAME REPLAY_DATA REPLAY_DATA_PATH NUM_PROCS REFERENCE_OUTPUT)
    set(multiValueArgs)
    cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(NOT ARGS_NUM_PROCS)
      set(ARGS_NUM_PROCS 1)
    endif()

    set(test_dir ${CMAKE_CURRENT_BINARY_DIR}/${ARGS_NAME})
    file(MAKE_DIRECTORY ${test_dir})

    add_test(NAME ${ARGS_NAME}
      COMMAND sh ${test_dir}/run.sh
      WORKING_DIRECTORY ${test_dir})

    if(NOT ARGS_REPLAY_DATA_PATH)
      if(ARGS_REPLAY_DATA)
        IF(NOT TARGET "download_replay_data_${ARGS_REPLAY_DATA}")
          message(FATAL_ERROR "Undefined replay_data: ${ARGS_REPLAY_DATA}. You can define it with comin_add_replay_data.")
        endif()
        get_target_property(REPLAY_DATA_PATH "download_replay_data_${ARGS_REPLAY_DATA}" DATA_DIRECTORY)
        IF(NOT REPLAY_DATA_PATH)
          message(FATAL_ERROR "Undefined replay_data: ${ARGS_REPLAY_DATA}")
        endif()
        add_dependencies("download_test_data" "download_replay_data_${ARGS_REPLAY_DATA}")
      else()
        message(FATAL_ERROR "No replay data provided")
      endif()
    else()
      set(REPLAY_DATA_PATH ${ARGS_REPLAY_DATA_PATH})
    endif()
    file(RELATIVE_PATH REPLAY_DATA_PATH "${test_dir}" "${REPLAY_DATA_PATH}")

    # add a dummy target to collect properties
    add_custom_target(_comin_test_${ARGS_NAME})
    set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY PLUGIN_COUNT 0)
    set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY NAMELIST "&comin_nml\n")
    set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY NUM_PROCS ${ARGS_NUM_PROCS})
    set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY EXTERNAL_PROCESSES "")

    if(NOT BUILD_TESTING)
      return()
    endif()

    # generate run.sh
    set(run_sh "#!/bin/sh\n")
    string(APPEND run_sh
      "cd ${test_dir}\n"
      "rm -rf output\n"
      "mkdir -p output\n"
      "set -e\n")
    if (MPI_FOUND)
      string(APPEND run_sh "${MPIEXEC}")

      string(APPEND run_sh " ${MPIEXEC_NUMPROC_FLAG} $<TARGET_PROPERTY:_comin_test_${ARGS_NAME},NUM_PROCS> ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ComIn::comin_replay> ${MPIEXEC_POSTFLAGS} ${ARGS_NAME}.nml ")

      # add external processes
      string(APPEND run_sh "$<TARGET_GENEX_EVAL:_comin_test_${ARGS_NAME},$<TARGET_PROPERTY:_comin_test_${ARGS_NAME},EXTERNAL_PROCESSES>>\n")
    else()
      if(${ARGS_NUM_PROCS} GREATER 1)
        message(FATAL_ERROR "MPI is required to run tests with more than 1 process.")
      endif()
      string(APPEND run_sh "$<TARGET_FILE:ComIn::comin_replay> ${ARGS_NAME}.nml\n")
    endif()

    # check the references
    if(ARGS_REFERENCE_OUTPUT)
      get_filename_component(ref_path ${ARGS_REFERENCE_OUTPUT} REALPATH BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
      string(APPEND run_sh "diff -w -r ${test_dir}/output ${ref_path} \n")

      add_custom_target(update_test_references_${ARGS_NAME})
      add_custom_command(TARGET update_test_references_${ARGS_NAME}
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${ref_path}
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${test_dir}/output ${ref_path}
      )
      add_dependencies(update_test_references update_test_references_${ARGS_NAME})
    endif()

    file(GENERATE OUTPUT ${test_dir}/run.sh
      CONTENT ${run_sh}
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ)

    # generate master.nml
    file(GENERATE OUTPUT ${test_dir}/${ARGS_NAME}.nml
      CONTENT
      "!
&replay_tool_nml
  replay_data_path=\"${REPLAY_DATA_PATH}/\"
$<$<BOOL:${ARGS_REFERENCE_OUTPUT}>:  output_dir=\"output/\"\n>/
$<TARGET_GENEX_EVAL:_comin_test_${ARGS_NAME},$<TARGET_PROPERTY:_comin_test_${ARGS_NAME},NAMELIST>>/\n"
    )

    comin_test_add_plugin(TEST ${ARGS_NAME}
      NAME "comin_var_replay_plugin"
      PLUGIN_LIBRARY "$<TARGET_FILE:ComIn::comin_var_replay_plugin>"
      OPTIONS "${REPLAY_DATA_PATH}/vars_"
    )

  endfunction()

  # define test datasets
  include(${CMAKE_CURRENT_LIST_DIR}/predef_replay_data.cmake)
endif()

# DOCUMENTATION CAN BE FOUND IN doc/cmake.md
# PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
function(comin_test_add_plugin)
  if(NOT BUILD_TESTING)
    return()
  endif()
  set(options LOG_DEBUG NO_LOG_INFO NO_LOG_WARNING)
  set(oneValueArgs TEST NAME PLUGIN_LIBRARY PRIMARY_CONSTRUCTOR OPTIONS COMM)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(ARGS_UNPARSED_ARGUMENTS)
    message(WARNING "Unknown arguments: ${ARGS_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT ARGS_TEST)
    message(FATAL_ERROR "No test provided")
  endif()

  if( NOT ARGS_PLUGIN_LIBRARY AND NOT ARGS_PRIMARY_CONSTRUCTOR)
    message(FATAL_ERROR "Either PLUGIN_LIBRARY or PRIMARY_CONSTRUCTOR must be provided")
  endif()

  get_property(PLUGIN_COUNT TARGET _comin_test_${ARGS_TEST} PROPERTY PLUGIN_COUNT)
  MATH(EXPR PLUGIN_COUNT "${PLUGIN_COUNT}+1")
  set_property(TARGET _comin_test_${ARGS_TEST} PROPERTY PLUGIN_COUNT ${PLUGIN_COUNT})
  set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
    "  plugin_list(${PLUGIN_COUNT})%name = \"${ARGS_NAME}\"\n")
  if(ARGS_PLUGIN_LIBRARY) # optional
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%plugin_library = \"${ARGS_PLUGIN_LIBRARY}\"\n")
  endif()
  if(ARGS_PRIMARY_CONSTRUCTOR)
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%primary_constructor = \"${ARGS_PRIMARY_CONSTRUCTOR}\"\n")
  endif()
  if(ARGS_OPTIONS) # optional
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%options = \"${ARGS_OPTIONS}\"\n")
  endif()
  if(ARGS_COMM) # optional
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%comm = \"${ARGS_COMM}\"\n")
  endif()
  if(ARGS_LOG_DEBUG)
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%log_debug = .TRUE.\n")
  endif()
  if(ARGS_NO_LOG_INFO)
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%log_info = .FALSE.\n")
  endif()
  if(ARGS_NO_LOG_WARNING)
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%log_warning = .FALSE.\n")
  endif()
endfunction()

# DOCUMENTATION CAN BE FOUND IN doc/cmake.md
# PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
function(comin_test_add_external_process)
  if(NOT BUILD_TESTING)
    return()
  endif()
  if(NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is required to add tests with external processes.")
  endif()
  set(options)
  set(oneValueArgs TEST NUM_PROCS COMMAND)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(ARGS_UNPARSED_ARGUMENTS)
    message(WARNING "Unknown arguments: ${ARGS_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT ARGS_TEST)
    message(FATAL_ERROR "No test provided")
  endif()
  if(NOT ARGS_COMMAND)
    message(FATAL_ERROR "No command provided")
  endif()
  if(NOT ARGS_NUM_PROCS)
    set(ARGS_NUM_PROCS 1)
  endif()
  set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY EXTERNAL_PROCESSES
    "  : ${MPIEXEC_NUMPROC_FLAG} ${ARGS_NUM_PROCS} ${ARGS_COMMAND}")
endfunction()
