#include<IoData.h>
#include<FluidCommunicator.h>
#include<HelperFunctions.h>

int verbose;
double start_time;
MPI_Comm piston_comm;

/*************************************
 * Main Function
 ************************************/
int main(int argc, char* argv[])
{

  //! Initialize MPI
  MPI_Init(NULL, NULL);
  start_time = walltime();

  //! Read user's input file (read the parameters)
  IoData iod(argc, argv);
  verbose = iod.output_data.verbose;

  //! Partition MPI
  piston_comm = MPI_COMM_WORLD;  
  MPI_Comm comm;
  FluidCommunicator communicator(MPI_COMM_WORLD, comm);
  piston_comm = comm; // correct the communicator

  //! Finalize IoData (read additional files and check for errors)
  iod.finalize();

  if(!communicator.IsCoupled()) {
    print_error("*** Error: PISTON needs to be coupled with a fluid solver.\n");
    exit_mpi();
  }





  //! Initialize the spring state 
  GeomState geom_state(iod.spring_data);

  //! Embedded surface and force
  TriangulatedSurface surface;
  std::vector<Vec3D> force;

  HelperFunctions::SetupEmbeddedSurface(geom_state, surface);
  force.resize(surface.active_nodes, Vec3D(0.0));

  /*************************************
   * Main Loop 
   ************************************/
  double t_n = 0.0; //!< simulation (i.e. physical) time
  int time_step = 0;
  double dt = iod.time_data.dt;
  double tmax = iod.time_data.tmax;




  communicator.ComputeAndSetTimeInfo(dt, tmax);
  communicator.InitializeMessengers(&surface, &force);
  communicator.CommunicateBeforeTimeStepping();

  HelperFunctions::InitializeAcceleration(geom_state, force);

  //HelperFunctions::UpdateSurfaceStaticLocation(geom_state, surface, force);

  HelperFunctions::OutputTriangulatedMesh(iod.output_data, surface.X0, surface.elems);
  HelperFunctions::OutputSolution(iod.output_data, surface, force, dt, t_n, time_step, true);

  while(t_n < tmax+0.01*dt) {

    HelperFunctions::AdvancePartialTimeStep(geom_state, surface, dt);

    // C0: Send predicted displacements and get updated forces
    communicator.Exchange();

    // Update solution
    HelperFunctions::AdvancePartialTimeStep(geom_state, surface, force, dt);

    time_step++;
    t_n += dt;

    HelperFunctions::OutputSolution(iod.output_data, surface, force, dt, t_n, time_step);

  }
  
  //! finalize
  communicator.Destroy();
  MPI_Finalize();

  return 0;

}
