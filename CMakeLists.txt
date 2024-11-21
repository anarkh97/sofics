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
option(BUILD_AEROS "Build local copy of Aero-S" ON)
option(BUILD_M2C "Build local copy of M2C" ON)

# add aeros
if(BUILD_AEROS)
  # fortran flags -- for aeros
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")

  # Zoltan library -- might be required by aero-s
  set(ZOLTAN_INCLUDE_PATH "/projects/wang_aoe_lab/KevinWang/Software/Trilinos/MyBuild/include/")
  set(ZOLTAN_SRC_PATH "/projects/wang_aoe_lab/KevinWang/Software/Trilinos/packages/zoltan/src/")
  set(ZOLTAN_LIBRARY "/projects/wang_aoe_lab/KevinWang/Software/Trilinos/MyBuild/lib/libzoltan.a") 

  # Spooles library -- might be required by aero-s
  set(SPOOLES_INCLUDE_PATH "/projects/wang_aoe_lab/KevinWang/Software/SPOOLES/")
  set(SPOOLES_spooles_LIBRARY "/projects/wang_aoe_lab/KevinWang/Software/SPOOLES/spooles.a")
  set(SPOOLES_spoolesMT_LIBRARY "/projects/wang_aoe_lab/KevinWang/Software/SPOOLES/MT/src/spoolesMT.a")

  # BLACS and SCALAPACK library
  #set(BLACS_blacs_LIBRARY "")
  #set(BLACS_Cinit_LIBRARY "")
  #set(SCALAPACK_LIBRARY "")

  add_subdirectory(packages/aeros)
endif()

# add m2c
if(BUILD_M2C)
  add_subdirectory(packages/m2c)
endif()

# add sofics source
add_subdirectory(src)
