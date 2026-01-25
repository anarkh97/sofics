# Aero-S build file

# -----------------------------------------------------------------------------
# Fortran flags
# -----------------------------------------------------------------------------
set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran-10")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")

# -----------------------------------------------------------------------------
# Aero-S variables which can be used to specify your custom libraries.
# -----------------------------------------------------------------------------
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

set(USE_CHOLMOD OFF)
set(USE_PARDISO OFF)
set(TRY_ACME OFF)

# -----------------------------------------------------------------------------
# Call Aero-S CMakeLists
# -----------------------------------------------------------------------------
add_subdirectory(
  "${PROJECT_SOURCE_DIR}/packages/aeros" 
  "${CMAKE_BINARY_DIR}/aeros-build"
)

# -----------------------------------------------------------------------------
# Change output directory and install target 
# -----------------------------------------------------------------------------
set_target_properties(femexecutable
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

#add_executable(aeros ALIAS femexecutable)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(AEROS_OUTPUT_NAME "aeros.debug")
else()
  set(AEROS_OUTPUT_NAME "aeros")
endif()

install(PROGRAMS
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${AEROS_OUTPUT_NAME}
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT aeros
  RENAME aeros
)
