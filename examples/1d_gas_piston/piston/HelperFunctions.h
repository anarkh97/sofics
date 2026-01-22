#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include<vector>
#include<IoData.h>
#include<Vector3D.h>
#include<GeomState.h>
#include<TriangulatedSurface.h>
#include<Utils.h>

namespace HelperFunctions {

  //! Problem initialization methods
  void SetupEmbeddedSurface(GeomState &geom_state,
                                   TriangulatedSurface &surf);

  void InitializeAcceleration(GeomState &geom_state,
                                     std::vector<Vec3D> &force);

  //! Time stepping methods
  void AdvancePartialTimeStep(GeomState &geom_state, 
                                     TriangulatedSurface &surf,
                                     double dt);

  void AdvancePartialTimeStep(GeomState &geom_state, 
                                     TriangulatedSurface &surf,
                                     std::vector<Vec3D> &force,
                                     double dt);

  //! Output functions
  void OutputTriangulatedMesh(const OutputData &iod_output, 
                              const std::vector<Vec3D> &X0,
                              const std::vector<Int3> &elems);
  
  void OutputSolution(const OutputData &iod_output,
                      const TriangulatedSurface &surface,
                      const std::vector<Vec3D> &force,
                      double dt, double time, int time_step,
                      bool write_header=false);

};

#endif
