add_executable(cdi
  cdi.c printinfo.c
)

target_include_directories(cdi PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include
)

target_link_libraries(cdi PRIVATE cdilib)
target_include_directories(cdi PRIVATE "${PROJECT_BINARY_DIR}/src/config.h")
target_compile_definitions(cdi PRIVATE  ${cdi_netcdf_definitions})
