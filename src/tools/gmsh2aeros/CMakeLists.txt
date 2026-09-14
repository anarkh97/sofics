add_executable(gmsh2aeros
  gmsh2aeros.cpp
)

install(TARGETS gmsh2aeros
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} 
  COMPONENT tools
)
