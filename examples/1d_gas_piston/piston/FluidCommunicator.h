#ifndef _FLUID_COMMUNICATOR_H_
#define _FLUID_COMMUNICATOR_H_

#include<mpi.h>
#include<memory>
#include<utility>
#include<TriangulatedSurface.h>
#include<Vector3D.h>
#include<GeomState.h>

/************************************************
 * Simple communicator class for exchanging 
 * pressure and velocities between the fluid
 * solver (M2C/FEST) and the structural solver
 * (PISTON). The class combines implementation
 * of M2C's ConcurrentProgramsHandler,
 * AerosMessenger, and M2CTwinMessenger.
 ***********************************************/

class FluidCommunicator {

  struct SendInfoData {
    int fluid_rank;
    int start_index;
    std::vector<int> node_ids;
  };

  int fluid_color; //!< the id ("color") of M2C/FEST in the split
  int piston_color; //!< the id ("color") of PISTON in split
  int max_color; //!< the total number of "colors" (must be same in all concurrent programs)

  //! Communicators
  MPI_Comm global_comm; //!< the global communicator
  int global_size, global_rank;

  MPI_Comm piston_comm; //!< the communicator for PISTON
  int piston_size, piston_rank;

  MPI_Comm joint_comm; //!< joint communicator between PISTON and M2C/FEST
  std::vector<SendInfoData> send_table;

  //! vector with all communicators
  std::vector<MPI_Comm> all_comm;

  //! Embedded surface
  std::shared_ptr<TriangulatedSurface> surf;
  //! Forces on the embedded boundary
  std::shared_ptr<std::vector<Vec3D>> force;

  //! Time stepping
  double dt, tmax;

public:

  FluidCommunicator(MPI_Comm global_comm_, MPI_Comm &comm_);
  ~FluidCommunicator();

  bool IsCoupled(); 

  //! functions that set info to be sent
  void ComputeAndSetTimeInfo(double dt_, double tmax_);
  void InitializeMessengers(TriangulatedSurface *surf_,
                            std::vector<Vec3D> *F_);

  //! main functions that handle communications
  void CommunicateBeforeTimeStepping();
  void Exchange();
  void FinalExchange();

  void Destroy();

private:

  void SetupCommunicators(); //!< called by the constructure similar to M2C 

  void Negotiate();
  void SendInfo();
  void SendEmbeddedWetSurface();
  void SendEmbeddedWetSurfaceInfo();
  void SendDisplacementAndVelocity();
  void GetForce();

};

#endif
