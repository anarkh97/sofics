#include<HelperFunctions.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<cstdio>

extern MPI_Comm piston_comm;

//------------------------------------------------------------------------------
// Helper methods local this cpp file.
//! Used for writing output
static double last_snapshot_time;
static std::vector<double> last_probe_time;

static void OutputSurfaceSolution(const char *filename, 
                                  const std::vector<Vec3D> &X0, 
                                  const std::vector<Vec3D> &X,
                                  double t, 
                                  bool write_header=false);

static void OutputProbeSolution(const char *filename,
                                const std::vector<Vec3D> &X0,
                                const std::vector<Vec3D> &X,
                                int id,
                                double t,
                                bool write_header=false);

static void OutputProbeSolution(const char *filename,
                                const std::vector<Vec3D> &V,
                                int id,
                                double t,
                                bool write_header=false);
                                  
//-----------------------------------------------------------------------------

void HelperFunctions::SetupEmbeddedSurface(GeomState &state, TriangulatedSurface &surf)
{

  std::vector<Vec3D> &X0 = surf.X0;
  std::vector<Vec3D> &X = surf.X;
  std::vector<Vec3D> &Udot = surf.Udot;
  std::vector<Int3> &Es = surf.elems;

  //! Surface x-coordinate.
  X0.resize(4, Vec3D(0.0)); 
  X0[0] = Vec3D(1000.0,  0.5,  0.5);
  X0[1] = Vec3D(1000.0, -0.5,  0.5);
  X0[2] = Vec3D(1000.0, -0.5, -0.5);
  X0[3] = Vec3D(1000.0,  0.5, -0.5);
  //X0[4] = Vec3D( 1000.0,    0.0,    0.0);

  X.resize(4, Vec3D(0.0));
  for(int i=0; i<(int)X.size(); ++i) {
    X[i][0] = X0[i][0] + state.disp;
    X[i][1] = X0[i][1];
    X[i][2] = X0[i][2];
  }

  Udot.resize(4, Vec3D(0.0));
  for(int i=0; i<(int)Udot.size(); ++i)
    Udot[i][0] = state.velo;

  Es.resize(2, Int3(-1));
  Es[0] = Int3(0, 2, 1);
  Es[1] = Int3(0, 3, 2);
  //Es[0] = Int3(1, 0, 4);
  //Es[1] = Int3(2, 1, 4);
  //Es[2] = Int3(3, 2, 4);
  //Es[3] = Int3(0, 3, 4);

  surf.active_nodes = X0.size();
  surf.active_elems = Es.size();
  surf.BuildConnectivities();
  surf.CalculateNormalsAndAreas();

}

//-----------------------------------------------------------------------------

void HelperFunctions::InitializeAcceleration(GeomState &state, std::vector<Vec3D> &force)
{

  double &u0 = state.disp;
  double &ms = state.mass;
  double &ks = state.stiffness;

  double &a0 = state.accl; 

  //! Calculate total force on the spring
  double f = 0.0;
  for(int i=0; i<(int)force.size(); ++i)
    f += force[i].norm();

  a0 = (f - ks*u0)/ms;

}

//-----------------------------------------------------------------------------

void HelperFunctions::AdvancePartialTimeStep(GeomState &state, TriangulatedSurface &surf,
                                             double dt)
{

  //! Calculate partial update velocities
  double v_n_h = state.velo + 0.5*dt*state.accl;
  
  //! Update displacement
  double d_n_p = state.disp + dt*v_n_h;

  //! Update surface and state
  for(int i=0; i<(int)surf.X.size(); ++i) {
    surf.X[i] = surf.X0[i] + Vec3D(d_n_p,0.0,0.0); 
    surf.Udot[i] = Vec3D(v_n_h,0.0,0.0);
  }

  state.disp = d_n_p;
  state.velo = v_n_h;

}

//-----------------------------------------------------------------------------

void HelperFunctions::AdvancePartialTimeStep(GeomState &state, TriangulatedSurface &surf,
                                             std::vector<Vec3D> &force, double dt)
{

  assert(surf.active_nodes == (int)force.size());

  double &d_n_p = state.disp;
  double &v_n_h = state.velo;

  double &ms = state.mass;
  double &ks = state.stiffness;

  //! Sum forces on the surface
  double f = 0.0;
  for(int i=0; i<(int)force.size(); ++i)
    f += force[i][0];

  //! Update surface and state
  state.accl = (f - ks*d_n_p)/ms;
  state.velo = v_n_h + 0.5*dt*state.accl;

  for(int i=0; i<(int)surf.X.size(); ++i) {
    surf.Udot[i] = Vec3D(state.velo,0.0,0.0);
  }

}

//-----------------------------------------------------------------------------

void HelperFunctions::OutputTriangulatedMesh(const OutputData &iod_output,
                                             const std::vector<Vec3D> &X0, 
                                             const std::vector<Int3> &elems)
{

  if(!strcmp(iod_output.surf_file,""))
    return; 

  //! Repurposed from M2C LagrangianOutput::OutputTriangulatedMesh
  int mpi_rank = -1;
  MPI_Comm_rank(piston_comm, &mpi_rank);

  //Only Proc #0 writes
  if(mpi_rank == 0) {

    char filename[512];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.surf_file);

    std::fstream out;
    out.open(filename, std::fstream::out);

    if(!out.is_open()) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", filename);
      exit(-1);
    }

    out << "Nodes MyNodes" << std::endl;

    for(int i=0; i<(int)X0.size(); i++) 
      out << std::setw(10) << i+1 
          << std::setw(14) << std::scientific << X0[i][0]
          << std::setw(14) << std::scientific << X0[i][1] 
          << std::setw(14) << std::scientific << X0[i][2] << "\n";

    out << "Elements MyElems using MyNodes" << std::endl;

    for(int i=0; i<(int)elems.size(); i++)
      out << std::setw(10) << i+1 << "  4  "  //"4" for triangle elements
          << std::setw(10) << elems[i][0]+1
          << std::setw(10) << elems[i][1]+1
          << std::setw(10) << elems[i][2]+1 << "\n";

    out.flush();
    out.close();
  }

  MPI_Barrier(piston_comm);

}

//------------------------------------------------------------------------------

void HelperFunctions::OutputSolution(const OutputData &iod_output,
                                     const TriangulatedSurface &surface,
                                     const std::vector<Vec3D> &force,
                                     double dt, double time, int time_step,
                                     bool write_header)
{

  //! reset counters when header is requested
  if(write_header) {
    last_snapshot_time = -1;
    last_probe_time.resize(iod_output.probes.dataMap.size(), -1);
  }

  // first output probes 
  for(auto it=iod_output.probes.dataMap.begin();
      it != iod_output.probes.dataMap.end(); ++it) {

    const ProbeData &probe = *(it->second);

    // checks
    if(probe.node_id < 0 || probe.node_id >= surface.active_nodes) continue;

    if(probe.frequency<=0.0 && probe.frequency_dt<=0.0) continue;

    if(!isTimeToWrite(time, dt, time_step, probe.frequency_dt, 
                      probe.frequency, last_probe_time[it->first], false)) {
      continue;
    };


    // output
    char outname[512];
    if(strcmp(probe.disp_file,"")) {
      sprintf(outname, "%s%s", iod_output.prefix, probe.disp_file);
      OutputProbeSolution(outname, surface.X0, surface.X, probe.node_id, time, write_header); 
    }
    if(strcmp(probe.velo_file, "")) {
      sprintf(outname, "%s%s", iod_output.prefix, probe.velo_file);
      OutputProbeSolution(outname, surface.Udot, probe.node_id, time, write_header); 
    }
    if(strcmp(probe.forc_file, "")) {
      sprintf(outname, "%s%s", iod_output.prefix, probe.forc_file);
      OutputProbeSolution(outname, force, probe.node_id, time, write_header); 
    }

    last_probe_time[it->first] = time;

  }


  //! Output surface solutions ...
  if(iod_output.frequency <= 0.0 && iod_output.frequency_dt <= 0.0)
    return;

  if(!strcmp(iod_output.disp_file, ""))
    return;

  if(!isTimeToWrite(time, dt, time_step, iod_output.frequency_dt, iod_output.frequency, last_snapshot_time, false))
    return;

  char outname[512];
  sprintf(outname, "%s%s", iod_output.prefix, iod_output.disp_file);
  OutputSurfaceSolution(outname, surface.X0, surface.X, time, write_header);
  last_snapshot_time = time;

}

//-----------------------------------------------------------------------------

static void OutputSurfaceSolution(const char *filename, const std::vector<Vec3D>& X0, 
                                  const std::vector<Vec3D>& X, double t, 
                                  bool write_header)
{

  int mpi_rank = -1;
  MPI_Comm_rank(piston_comm, &mpi_rank);

  //Only Proc #0 writes
  if(mpi_rank == 0) {
    FILE *disp_file = (write_header) ? fopen(filename, "w") : fopen(filename, "a");
    if(disp_file == NULL) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", filename);
      exit(-1);
    }

    if(write_header) {
      fprintf(disp_file, "Vector DISP under NLDynamic for MyNodes\n");
      fprintf(disp_file, "%d\n", (int)X0.size());
    }

    fprintf(disp_file,"%e\n", t);
    for(int i=0; i<(int)X0.size(); i++)
      fprintf(disp_file,"%12.8e    %12.8e    %12.8e\n", X[i][0]-X0[i][0], X[i][1]-X0[i][1], X[i][2]-X0[i][2]);

    fclose(disp_file);
    print("- Wrote solution on a Lagrangian mesh at %e.\n", t);
  }

  MPI_Barrier(piston_comm);

}

//-----------------------------------------------------------------------------

static void OutputProbeSolution(const char *filename, const std::vector<Vec3D> &X0, 
                                const std::vector<Vec3D> &X, int id, double t,
                                bool write_header)
{


  int mpi_rank = -1;
  MPI_Comm_rank(piston_comm, &mpi_rank);

  //Only Proc #0 writes
  if(mpi_rank == 0) {
    FILE *probe_file = (write_header) ? fopen(filename, "w") : fopen(filename, "a");
    if(probe_file == NULL) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", filename);
      exit(-1);
    }

    if(write_header)
      fprintf(probe_file, "# node %d\n", id);

    fprintf(probe_file,"%e    ", t);
    fprintf(probe_file,"%12.8e    %12.8e    %12.8e\n", X[id][0]-X0[id][0], X[id][1]-X0[id][1], X[id][2]-X0[id][2]);

    fclose(probe_file);
  }

  MPI_Barrier(piston_comm);
  
}

//-----------------------------------------------------------------------------

static void OutputProbeSolution(const char *filename, const std::vector<Vec3D> &V, 
                                int id, double t, bool write_header)
{

  int mpi_rank = -1;
  MPI_Comm_rank(piston_comm, &mpi_rank);

  //Only Proc #0 writes
  if(mpi_rank == 0) {
    FILE *probe_file = (write_header) ? fopen(filename, "w") : fopen(filename, "a");
    if(probe_file == NULL) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", filename);
      exit(-1);
    }

    if(write_header)
      fprintf(probe_file, "# node %d\n", id);

    fprintf(probe_file,"%e    ", t);
    fprintf(probe_file,"%12.8e    %12.8e    %12.8e\n", V[id][0], V[id][1], V[id][2]);

    fclose(probe_file);
  }

  MPI_Barrier(piston_comm);
  
}
//-----------------------------------------------------------------------------

