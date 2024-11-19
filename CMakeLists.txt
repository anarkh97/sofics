cmake_minimum_required(VERSION 3.10)

# set the project name
project(SOFICS VERSION 1.0)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# compile packages
option(BUILD_AEROS "Build local copy of Aero-S" OFF)
option(BUILD_M2C "Build local copy of M2C" ON)

# add aeros
if(BUILD_AEROS)
  add_subdirectory(packages/aeros)
endif()

# add m2c
if(BUILD_M2C)
  add_subdirectory(packages/m2c)
endif()

# add sofics source
add_subdirectory(src)
