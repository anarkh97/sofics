# Make a environment script with current build details
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/driver.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/driver.sh
  @ONLY
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
      COMMENT "Copying ${fname} to build directory"
    )

    list(APPEND SOFICS_SCRIPT_OUTPUTS "${CMAKE_CURRENT_BINARY_DIR}/${fname}")
  endforeach()

  # target for copying scripts during the build
  add_custom_target(sofics_scripts ALL
    DEPENDS ${SOFICS_SCRIPT_OUTPUTS}
    COMMENT "Building SOFICS bash scripts"
  )

  # install from build
  foreach(shfile ${BASH_FILES})
    get_filename_component(fname ${shfile} NAME)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${fname}
      DESTINATION ${CMAKE_INSTALL_BINDIR}
      COMPONENT scripts
    )
  endforeach()
  
else()
  message(ERROR "The SOFICS bash scripts not found.")
endif()

# Install configured environment setup
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/driver.sh
  DESTINATION ${CMAKE_INSTALL_BINDIR} 
  COMPONENT scripts
)
