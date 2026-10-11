# M2C build file

# -----------------------------------------------------------------------------
# Temporary patch: M2C frees unset/aliased communicators on exit.
# Drop once fixed upstream.
# -----------------------------------------------------------------------------
set(M2C_PATCH "${PROJECT_SOURCE_DIR}/cmake/patches/M2CCommFree.patch")

# A clean reverse-apply means the patch is already in the working tree.
execute_process(COMMAND git apply --reverse --check "${M2C_PATCH}"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/src/solvers/m2c"
  RESULT_VARIABLE applied OUTPUT_QUIET ERROR_QUIET)

if(NOT applied EQUAL "0")
  execute_process(COMMAND git apply "${M2C_PATCH}"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/src/solvers/m2c"
    COMMAND_ERROR_IS_FATAL ANY)
  message(STATUS "Applied the communicator error patch.")
else()
  message(STATUS "Communicator error patch already applied.")
endif()

# -----------------------------------------------------------------------------
# Call M2C CMakeLists
# -----------------------------------------------------------------------------
add_subdirectory(
  "${PROJECT_SOURCE_DIR}/src/solvers/m2c" 
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
