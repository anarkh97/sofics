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
postprocessor.cpp
)

# Collect all .sh files in the source directory into a list
file(GLOB BASH_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.sh")

# Copy the collected files to the desired destination
if(BASH_FILES)

  set(SOFICS_SCRIPT_OUTPUTS "")

  foreach(shfile ${BASH_FILES})
    get_filename_component(fname ${shfile} NAME)
    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${fname}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${shfile} "${CMAKE_CURRENT_BINARY_DIR}/${fname}"
      DEPENDS ${shfile}
    )

    list(APPEND SOFICS_SCRIPT_OUTPUTS "${CMAKE_CURRENT_BINARY_DIR}/${fname}")
  endforeach()

  add_custom_target(sofics_scripts ALL
    DEPENDS ${SOFICS_SCRIPT_OUTPUTS}
    COMMENT ""
  )

else()
  message(ERROR "The SOFICS bash files not found.")
endif()
