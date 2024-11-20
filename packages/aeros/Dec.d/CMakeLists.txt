ADD_LIBRARY(Dec Decomp.d/Greedy.C Decomp.d/DecOpti.C 
Geom.d/IPoints.C Geom.d/Mesh.C Geom.d/PolyObj.C
dec.C ElemMFCheck.C AddE.C)

SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
TARGET_LINK_LIBRARIES(Dec Eigen3::Eigen)
