add_subdirectory(gmsh2aeros)
add_subdirectory(aeros2dakota)

add_custom_target(internal-tools)
add_dependencies(internal-tools
  gmsh2aeros
  aeros2dakota
)
