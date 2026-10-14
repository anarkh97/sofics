set(MESH_DATA "${CMAKE_CURRENT_SOURCE_DIR}/data/one_tet.msh")
set(MASS_DATA "${CMAKE_CURRENT_SOURCE_DIR}/data/mass.out")
set(DISP_DATA "${CMAKE_CURRENT_SOURCE_DIR}/data/disp")

# One tetrahedron in, NODES and TOPOLOGY out.
add_test(NAME gmsh2aeros_convert COMMAND sh -c
  "$<TARGET_FILE:gmsh2aeros> ${MESH_DATA} one_tet.include &&
   grep -q '^NODES' one_tet.include &&
   grep -q '^TOPOLOGY' one_tet.include")

# A scalar Aero-S result
add_test(NAME aeros2dakota_scalar COMMAND sh -c
  "$<TARGET_FILE:aeros2dakota> -w -s ${MASS_DATA} -d scalar_results.out &&
   grep -q '1.312199e+03' scalar_results.out")

# A vector Aero-S result
add_test(NAME aeros2dakota_disp COMMAND sh -c
  "$<TARGET_FILE:aeros2dakota> -w -a ${DISP_DATA} -d disp_results.out &&
   grep -q '5.000000e-02' disp_results.out")
