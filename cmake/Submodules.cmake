# Submodule management script

if(UPDATE_SUBMODULES)

  find_package(GIT QUIET REQUIRED)
  
  # PROJECT_SOURCE_DIR is the set by the last project() call
  # and is unaffected by add_subdirectory() changes
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(STATUS "Updating Submodules ...")

    # perform recursive submodule update in the main SOFICS directory
    # i.e., where SOFICS CMakeLists is placed.
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
      OUTPUT_VARIABLE GIT_SUBMOD_OUTPUT
      ERROR_VARIABLE GIT_SUBMOD_ERROR
    )
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(WARNING "git submodule update failed with ${GIT_SUBMOD_RESULT}")
      message(WARNING "stdout: ${GIT_SUBMOD_OUTPUT}")
      message(WARNING "stderr: ${GIT_SUBMOD_ERROR}")
    endif()
  endif()

endif()

# -----------------------------------------------------------------------------
# Verify Required Submodules
# -----------------------------------------------------------------------------
if(BUILD_AEROS)
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/packages/aeros/.git")
    message(FATAL_ERROR 
    "Aero-S submodule not found!\n"
    "Please run: git submodule update --init --recursive")
  endif()
  message(STATUS "Aero-S submodule found")
endif()

if(BUILD_M2C)
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/packages/m2c/.git")
    message(FATAL_ERROR 
    "M2C submodule not found!\n"
    "Please run: git submodule update --init --recursive")
  endif()
  message(STATUS "M2C submodule found")
endif()
