# Aero-S build file

# -----------------------------------------------------------------------------
# Fortran flags
# -----------------------------------------------------------------------------
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
set_target_properties(aeros
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

install(TARGETS aeros
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT aeros
)
