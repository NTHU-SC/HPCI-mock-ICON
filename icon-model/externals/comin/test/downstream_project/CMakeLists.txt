cmake_minimum_required(VERSION 3.20)
project(comin_downstream
  VERSION 0.1.0
  LANGUAGES C
  DESCRIPTION "Dummy ComIn downstream Plugin")


include(CTest)

find_package(ComIn REQUIRED)

add_subdirectory(src/)
add_subdirectory(test/)
