# M2C build file

# -----------------------------------------------------------------------------
# Call M2C CMakeLists
# -----------------------------------------------------------------------------
add_subdirectory(
  "${PROJECT_SOURCE_DIR}/packages/m2c" 
  "${CMAKE_BINARY_DIR}/m2c-build"
)

# -----------------------------------------------------------------------------
# Change output directory and install target 
# -----------------------------------------------------------------------------
set_target_properties(m2c
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

install(TARGETS m2c 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT m2c 
)
