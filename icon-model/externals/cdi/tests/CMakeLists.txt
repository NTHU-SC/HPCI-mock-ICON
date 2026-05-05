
add_executable(calendar_test1.run calendar_test1.c)
add_executable(cksum_read.run cksum_read.c )
add_executable(cksum_verify.run cksum_verify.c)
add_executable(cksum_write.run cksum_write.c)
add_executable(cksum_write_chunk.run cksum_write_chunk.c)
add_executable(deco2d_model.run deco2d_model.c)
add_executable(ensure_array_size.run ensure_array_size.c)
add_executable(pio_write.run pio_write.c)
add_executable(pio_write_setup_grid.run pio_write_setup_grid.c)
add_executable(simple_model.run simple_model.c)
add_executable(simple_model_helper.run simple_model_helper.c)
add_executable(stream_cksum.run stream_cksum.c)
add_executable(test_byteswap.run test_byteswap.c)
add_executable(test_cdf_read.run test_cdf_read.c)
add_executable(test_cdf_write.run test_cdf_write.c)
add_executable(test_grib.run test_grib.c)
add_executable(test_month_adjust.run test_month_adjust.c)
add_executable(test_resource_copy.run test_resource_copy.c)
add_executable(test_table.run test_table.c)

get_property(current_targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
list(APPEND test_lib_src
  var_cksum.c
  ensure_array_size.c
  stream_cksum.c
  simple_model.c
  simple_model_helper.c
  pio_write_setup_grid.c
)

add_library(test_lib ${test_lib_src})
target_include_directories(test_lib PRIVATE ../src/)
target_link_libraries(test_lib PRIVATE cdilib ${cdi_linked_libs})

foreach(target ${current_targets})
  target_include_directories(${target} PRIVATE ../src/)
  target_link_libraries(${target} PRIVATE cdilib test_lib ${cdi_linked_libs})
  add_test(NAME ${target} COMMAND ${target})
  if(labels_${target})
    set_property(TEST ${target} PROPERTY LABELS labels_${target})
  endif()
endforeach()
