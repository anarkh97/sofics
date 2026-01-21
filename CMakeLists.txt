cmake_minimum_required(VERSION 3.31)

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
option(BUILD_AEROS "Build local copy of Aero-S" ON)
option(BUILD_M2C "Build local copy of M2C" ON)

# add aeros
if(BUILD_AEROS)

  # fortran flags for aeros
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")

  #-----------------------------------------
  # Aero-S variables which can be used to specify your custom libraries.
  #
  # E.g.,
  #
  # Zoltan library
  # set(ZOLTAN_INCLUDE_PATH "/path/to/ZOLTAN/include/")
  # set(ZOLTAN_SRC_PATH "/path/to/ZOLTAN/src")
  # set(ZOLTAN_LIBRARY "/path/to/libzoltan.a") 
  #
  # Spooles library
  # set(SPOOLES_INCLUDE_PATH "/path/to/SPOOLES/include/")
  # set(SPOOLES_spooles_LIBRARY "/path/to/spooles.a")
  # set(SPOOLES_spoolesMT_LIBRARY "/path/to/spoolesMT.a")
  #
  # BLACS and SCALAPACK library
  # set(BLACS_blacs_LIBRARY "/path/to/libblacs.a")
  # set(BLACS_Cinit_LIBRARY "/path/to/libcinit.a")
  # set(SCALAPACK_LIBRARY "/path/to/libscalapack.a")
  #-----------------------------------------

  add_subdirectory(packages/aeros)
  
endif()

# add m2c
if(BUILD_M2C)
  add_subdirectory(packages/m2c)
endif()

# add sofics source
add_subdirectory(src)
