add_executable(aeros2dakota
  aeros2dakota.cpp
)

install(TARGETS aeros2dakota
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} 
  COMPONENT tools
)
