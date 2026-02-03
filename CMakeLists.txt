cmake_minimum_required(VERSION 3.31)

# set the project name
project(SOFICS VERSION 1.0)

include(GNUInstallDirs)

# -----------------------------------------------------------------------------
# Build configuration 
# -----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_AEROS "Build Aero-S" ON)
option(BUILD_M2C "Build M2C" ON)
option(BUILD_TOOLS "Build Gmsh2Aeros and Postprocessor" ON)
option(UPDATE_SUBMODULES "Check submodules during build" ON)

# -----------------------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------------------
# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# -----------------------------------------------------------------------------
# CMake Modules
# -----------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Submodule management (Aero-S/M2C)
include(Submodules)

# Find/fetch dependencies (related to Aero-S/M2C)
include(Dependencies)

# -----------------------------------------------------------------------------
# Build submodules (Aero-S and M2C)
# -----------------------------------------------------------------------------
# add aeros
if(BUILD_AEROS)
  message(STATUS "Building Aero-S ...")
  include(aeros)
endif()

# add m2c
if(BUILD_M2C)
  message(STATUS "Building M2C ...")
  include(m2c)
endif()

# -----------------------------------------------------------------------------
# Build internal tools
# -----------------------------------------------------------------------------
if(BUILD_TOOLS)
  message(STATUS "Building Gmsh2Aeros and Postprocessor ...")
  add_subdirectory(tools)
endif()

# -----------------------------------------------------------------------------
# SOFICS environment setup
# -----------------------------------------------------------------------------
add_subdirectory(scripts)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
message("")
message("========================================")
message("SOFICS ${PROJECT_VERSION} Configuration Summary")
message("========================================")
message("Build type:          ${CMAKE_BUILD_TYPE}")
message("Install prefix:      ${CMAKE_INSTALL_PREFIX}")
message("")
message("Components:")
message("  AERO-S:            ${BUILD_AEROS}")
message("  M2C:               ${BUILD_M2C}")
message("  Internal tools:    ${BUILD_TOOLS}")
message("")
message("========================================")
message("")
