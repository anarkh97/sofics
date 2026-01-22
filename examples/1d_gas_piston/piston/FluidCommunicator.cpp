#include<FluidCommunicator.h>
#include<Utils.h>
#include<cassert>

#define SUGGEST_DT_TAG 444
#define WET_SURF_TAG1 555
#define WET_SURF_TAG2 666
#define WET_SURF_TAG3 888
#define WET_SURF_TAG4 999
#define SUBCYCLING_TAG 777

#define NEGO_NUM_TAG 10000
#define NEGO_BUF_TAG 10001

#define FORCE_TAG 1000
#define DISP_TAG 2000
#define INFO_TAG 3000

static int bufsize = 6;

using std::vector;
using std::pair;
using std::make_pair;

//---------------------------------------------------------

FluidCommunicator::FluidCommunicator(MPI_Comm global_comm_, MPI_Comm &comm_)
                 : global_comm(global_comm_), piston_comm(global_comm_), 
                   joint_comm(global_comm_), surf(nullptr), force(nullptr),
                   dt(0.0), tmax(0.0)
{
  //The following parameters are the same as "FLUID_ID" and "MAX_CODES" in AERO-S and AERO-F
  fluid_color = 0; //"FLUID_ID" in M2C  
  piston_color = 1; //"STRUCT_ID" in AERO-S
  max_color = 4; 

  // simultaneous operations w/ other programs 
  SetupCommunicators();

  // create inter-communicators 
  joint_comm = all_comm[fluid_color];

  comm_ = piston_comm;
}

//---------------------------------------------------------

FluidCommunicator::~FluidCommunicator()
{

}

//---------------------------------------------------------

void
FluidCommunicator::Destroy()
{

}

//---------------------------------------------------------

bool
FluidCommunicator::IsCoupled() {
  int joint_size = -1;
  MPI_Comm_remote_size(joint_comm, &joint_size);
  return joint_size > 0;
}

//---------------------------------------------------------

void
FluidCommunicator::ComputeAndSetTimeInfo(double dt_, double tmax_)
{
  dt = dt_;
  tmax = tmax_;

  int max_step = (int)( (tmax + 0.49*dt)/dt );
  double remainder = tmax - max_step*dt;

  if(std::abs(remainder) > 0.01*dt) {
    tmax = max_step*dt;
    if(piston_rank == 0)
      fprintf(stdout, " Warning: Total time is being changed to : %e\n", tmax);
  }
}

//---------------------------------------------------------

//for M2C messengers
void
FluidCommunicator::InitializeMessengers(TriangulatedSurface *surf_, 
                                        vector<Vec3D> *F_)
{

  assert(surf_); //cannot be NULL
  assert(F_); //cannot be NULL

  surf.reset(surf_);
  force.reset(F_);

  // send info first
  SendEmbeddedWetSurfaceInfo();

  SendEmbeddedWetSurface();

  Negotiate();  

  SendInfo();

  // send sub-cycling info, required for C0.
  // M2C currently does not support this feature.
  double info = 0;
  if(piston_rank == 0) {
    MPI_Send(&info, 1, MPI_DOUBLE, 0, SUBCYCLING_TAG, joint_comm);
  }

}

//---------------------------------------------------------

void
FluidCommunicator::SetupCommunicators()
{

  MPI_Comm_rank(global_comm, &global_rank);
  MPI_Comm_size(global_comm, &global_size);

  //! split the global communicator into groups. This split call
  //! gives us the processes assigned to piston. A similar call is
  //! undertaken in m2c as well.
  MPI_Comm_split(global_comm, piston_color + 1, global_rank, &piston_comm);
  MPI_Comm_rank(piston_comm, &piston_rank);
  MPI_Comm_size(piston_comm, &piston_size);

  if(piston_size > 1 and piston_rank == 0) {
    fprintf(stderr, "\033[;31m*** Error: PISTON is a single processor "
                    "code.\033[0m\n");
    exit(-1);
  }

  all_comm.resize(max_color);

  all_comm[piston_color] = piston_comm;

  vector<int> leaders(max_color, -1);
  vector<int> newleaders(max_color, -1);

  if(piston_rank == 0) { // check not required but included for clarity.
    leaders[piston_color] = global_rank;
  }

  //! collects leaders from all programs.
  MPI_Allreduce(leaders.data(), newleaders.data(), max_color, MPI_INTEGER,
    MPI_MAX, global_comm);

  for(int i=0; i<max_color; i++) {
    if(i != piston_color && newleaders[i] >= 0) {
      // create a communicator between piston and program i
      int tag;
      if(piston_color < i)
        tag = max_color * (piston_color + 1) + i + 1;
      else
        tag = max_color * (i + 1) + piston_color + 1;

      MPI_Intercomm_create(piston_comm, 0, global_comm, newleaders[i], tag, &all_comm[i]);
    }
  }

/*
  int fluid_size = 0;
  MPI_Comm_remote_size(all_comm[fluid_color], &fluid_size);

  if(fluid_size < 1 and piston_rank == 0) {
    fprintf(stderr, "\033[;31m*** Error: PISTON is expected to run in a "
                    "coupled sense. Try running with M2C/FEST.\033[0m\n");
    exit(-1);
  }
*/

}

//---------------------------------------------------------

void
FluidCommunicator::Negotiate()
{

  int num_fluid_procs_comm = -1;
  MPI_Comm_remote_size(joint_comm, &num_fluid_procs_comm);
  for(int i_fluid=0; i_fluid<num_fluid_procs_comm; ++i_fluid) {
  
    // receive the number of matched nodes that the fluid 
    // rank claims it has.
    int fluid_matched_nodes = 0;
    MPI_Status status;
    MPI_Recv(&fluid_matched_nodes, 1, MPI_INT, MPI_ANY_SOURCE,
             NEGO_NUM_TAG, joint_comm, &status);

    if(fluid_matched_nodes == 0) continue;

    // this processor is involved in communication

    // rank of the fluid processor that sent the previous
    // message.
    int fluid_rank = status.MPI_SOURCE;

    // receive the list of nodes ids that the this fluid
    // rank claims to have.
    vector<int> node_index(fluid_matched_nodes, -1);
    MPI_Recv(node_index.data(), fluid_matched_nodes, MPI_INT,
             fluid_rank, NEGO_BUF_TAG, joint_comm, MPI_STATUS_IGNORE); 

    // update send table
    SendInfoData send_data;
    send_data.fluid_rank = fluid_rank;
    send_data.node_ids = node_index;

    send_table.push_back(send_data);

  }

  // calculate index displacements
  send_table[0].start_index = 0;
  for(int i=1; i<send_table.size(); ++i) {
    send_table[i].start_index =
      send_table[i-1].start_index + send_table[i-1].node_ids.size();
  }

  // In our case, the wetted surface is "not" split across multiple
  // processors. Hence, we send back whatever was received.
  vector<MPI_Request> send_requests;
  for(int i_fluid=0; i_fluid<send_table.size(); ++i_fluid) {

    send_requests.push_back(MPI_Request());

    int send_rank = send_table[i_fluid].fluid_rank;
    int send_size = send_table[i_fluid].node_ids.size();
    vector<int> &send_ids = send_table[i_fluid].node_ids;

    // send the actual matched nodes back to the fluid solver.
    MPI_Isend(&send_size, 1, MPI_INT, send_rank, NEGO_NUM_TAG,
      joint_comm, &(send_requests.back()));

    MPI_Isend(send_ids.data(), send_ids.size(), MPI_INT, send_rank, 
      NEGO_BUF_TAG, joint_comm, &(send_requests.back()));

  }
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

  //print("... [E] Negotiation successful ...\n");

}

//---------------------------------------------------------

void
FluidCommunicator::SendInfo()
{

  double info[5];

  info[0] = 22; // C0 --- (AN) Explicit-Explicit 
  info[1] = dt;
  info[2] = tmax;
  info[3] = 0; // Note: not used by M2C/FEST
  info[4] = 0; // Note: not used by M2C/FEST

  if(piston_rank == 0) {
    MPI_Send(info, 5, MPI_DOUBLE, 0, INFO_TAG, joint_comm);
  }

  //print(stderr, " ... [S] Sending info ...\n"); 

}

//---------------------------------------------------------

void
FluidCommunicator::SendEmbeddedWetSurfaceInfo()
{

  assert(surf);

  int info[4];

  info[0] = 3; // Element type fixed to triangles.
  info[1] = 0; // cracking is disabled here.
  info[2] = surf->active_nodes;
  info[3] = surf->active_elems;

  if(piston_rank == 0)
    MPI_Send(info, 4, MPI_INT, 0, WET_SURF_TAG1, joint_comm); 

  //print(stderr, " ... [S] Sending Embedded surface info ...\n"); 

}

//---------------------------------------------------------

void
FluidCommunicator::SendEmbeddedWetSurface()
{

  // send nodes
  if(piston_rank == 0)
    MPI_Send((double*)surf->X0.data(), surf->active_nodes*3, MPI_DOUBLE,
             0, WET_SURF_TAG2, joint_comm);

  // send elements --- element type fixed to triangles (See Info)
  if(piston_rank == 0)
    MPI_Send((int*)surf->elems.data(), surf->active_elems*3, MPI_INT,
             0, WET_SURF_TAG3, joint_comm);

  //print(stderr, " ... [S] Sending embedded surface ...\n"); 

}

//---------------------------------------------------------

void
FluidCommunicator::SendDisplacementAndVelocity()
{

  if(piston_rank == 0) {
    
    // Calculate the surface displacement and velocity
    int active_nodes = surf->active_nodes;
    vector<double> disp_velo_buffer(active_nodes*bufsize, 0.0);

    for(int i=0; i<active_nodes; ++i) {
      for(int j=0; j<3; ++j) {
        disp_velo_buffer[bufsize*i+j] = surf->X[i][j]-surf->X0[i][j];
        disp_velo_buffer[bufsize*i+3+j] = surf->Udot[i][j];
      }
    }

    vector<MPI_Request> send_requests;
    for(int i=0; i<send_table.size(); ++i) {
      int size = bufsize*send_table[i].node_ids.size();
      int first_index = bufsize*send_table[i].start_index;

      double *local_buffer = disp_velo_buffer.data() + first_index;

/*
      vector<double> local_buffer(bufsize*local_size, 0.0);
      for(int n=0; n<local_size; ++n) {
        int c_index = send_table[i].node_ids[n]; 
        for(int j=0; j<3; ++j) {
          local_buffer[bufsize*n+j] = disp_velo_buffer[bufsize*c_index+j];
          local_buffer[bufsize*n+3+j] = disp_velo_buffer[bufsize*c_index+3+j];
        }
      }
*/

      send_requests.push_back(MPI_Request());
      MPI_Isend(local_buffer, size, MPI_DOUBLE, send_table[i].fluid_rank,
        DISP_TAG, joint_comm, &(send_requests.back()));
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

  }

  //print(stderr, " ... [S] Sending displacements ...\n"); 

}

//---------------------------------------------------------

void
FluidCommunicator::GetForce()
{
  if(piston_rank == 0) {

    for(int i=0; i<send_table.size(); ++i) {
      vector<double> local_buffer(3*send_table[i].node_ids.size(), 0.0);
      MPI_Recv(local_buffer.data(), local_buffer.size(), MPI_DOUBLE, 
        send_table[i].fluid_rank, FORCE_TAG, joint_comm, MPI_STATUS_IGNORE);

      //assemble into our vector
      for(int n=0; n<send_table[i].node_ids.size(); ++n) {
        int c_id = send_table[i].node_ids[n];
        for(int j=0; j<3; ++j) {
          (*force)[c_id][j] = local_buffer[3*n+j];
        }
      }
    }

  }
  
  //print(stderr, " ... [S] Recieved fluid load ...\n"); 

}

//---------------------------------------------------------

void
FluidCommunicator::CommunicateBeforeTimeStepping()
{

  SendDisplacementAndVelocity();

  GetForce();

  SendInfo(); // send dt, tmax

}

//---------------------------------------------------------

void
FluidCommunicator::Exchange()
{

  SendInfo(); // send dt, tmax

  SendDisplacementAndVelocity();

  GetForce();

}

//---------------------------------------------------------

void
FluidCommunicator::FinalExchange()
{

  GetForce();

}

//---------------------------------------------------------



