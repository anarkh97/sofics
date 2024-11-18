cmake_minimum_required(VERSION 3.10)

# set the project name
project(SOFICS VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# compiler flags (turn on all the warnings)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# add the executables
add_executable(gmsh2aeros 
gmsh2aeros.cpp
)
add_executable(postprocessor
postprocessor
)

# copy bash scripts
# Note: It might be better to use configure_file() here instead of 
# file(COPY ...).The advantage of configure_file() is that it can 
# automatically update the files being copied, ensuring that the 
# destination files are always up-to-date. This is especially useful
# if there are any placeholders or variables within the files that 
# need to be replaced.
# TODO: will update in future.
# Collect all .sh files in the source directory into a list
file(GLOB BASH_FILES "${CMAKE_SOURCE_DIR}/*.sh")

# Copy the collected files to the desired destination
if(BASH_FILES)
  foreach(file ${BASH_FILES})
    file(COPY ${file} DESTINATION ${CMAKE_BINARY_DIR})
  endforeach()
else()
  message(ERROR "The SOFICS bash files not found.")
endif()
