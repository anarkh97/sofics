# set the project name
project(UTILITIES VERSION 1.0)

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
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# add the executables
add_executable(gmsh2aeros 
gmsh2aeros.cpp
)
add_executable(postprocessor
postprocessor
)

# copy bash scripts
add_custom_target(copy_scripts ALL)

# Collect all .sh files in the source directory into a list
file(GLOB BASH_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.sh")

# Copy the collected files to the desired destination
if(BASH_FILES)
  foreach(shfile ${BASH_FILES})
    get_filename_component(fname ${shfile} NAME)
    add_custom_command(
        TARGET copy_scripts
        COMMAND ${CMAKE_COMMAND} -E copy ${shfile} ${CMAKE_CURRENT_BINARY_DIR}/${fname}
        DEPENDS ${shfile}
    )
  endforeach()
else()
  message(ERROR "The SOFICS bash files not found.")
endif()
