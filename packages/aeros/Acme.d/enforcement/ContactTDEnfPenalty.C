// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#include "ContactTDEnfModel.h"
#include "ContactTDTied.h"
#include "ContactErrors.h"
#include "ContactTDEnfPenalty.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "ContactShellNode.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactScratchManager.h"
#include "cstring"
#include <cmath>
#include <cstdio>
#include <unistd.h>

// Internal constants
#undef  NDIM
#define NDIM 3
#undef  ZERO_TOL
#define ZERO_TOL 1.0E-14
#undef  SMALL
#define SMALL 1.0e-6
#undef  CONVERGENCE_TOL
#define CONVERGENCE_TOL 1.0e-10
#undef  FAIL_TOL
#define FAIL_TOL 0.01

// Internal flags 
// LOCAL_PRINT_FLAG = 0 turns off all internal debugging prints
#undef  LOCAL_PRINT_FLAG 
#define LOCAL_PRINT_FLAG 0

extern int contactPrintFlag;

#if LOCAL_PRINT_FLAG
#include "ContactParOStream.h"
#else
#ifdef CONTACT_DEBUG_NODE
#include "ContactParOStream.h"
#endif
#endif

/* SYMBOLS ______________________________________________________
 * x = position
 * u = displacement
 * v = velocity
 * a = acceleration
 * f = force
 * m = mass
 * nm = outward unit normal vector ( f. nm < 0 for compression )
 *      (nm is typically outward for the master face and the
 *       force is typically stored with the slave node)
 * uj = jump in displacement i.e. relative displacement
 * vj = jump in velocity i.e. relative velocity
 * gap = uj . nm  , where  gap > 0 for separation
 * slip = vj * dt
 * kp ? = kinematic partition
 *     ... subscripts ....
 * f_n = the normal component of force
 * f_t = the tangential force
 * x_sn = slave node position
 * x_mn = master node position
 * x_cp = contact point position
 * x_cur = "CURRENT" position ?? from last time step after contact correction
 * x_new = "NEW" position ?? from current iteration
 * x_pre = "PREDICTED" position ?? from central difference predictor w/o contact *         force.
 * _____________________________________________________________*/

extern "C" {

    void FORTRAN(kappa_3d)(
	  int& npair,     // pair number
          int& ndm,       // number of dimensions
          int& ndf,       // number of dofs per node
          Real* x,        // coordinates
          Real* u,        // dofs (displacements)
          int* ix1,       // original vertice on slave face
          int* rx1,       // reordered vertice on slave face
          Real* kappa,    // scale factor for mortar matrix
          int& ns,        // # of node on slave surface
          int& nope1,     // number of vertice on slave element
          int& selem);    // number od slave elements
  
    void FORTRAN(mortar_matrix)(
          int& npair,       // pair number
          int& ndm,         // number of dimensions
          int& ndf,         // number of dofs per node
          Real* x,	       // coordinates
          Real* u,	       // dofs (displacements)
          int* ix1, int* ix2,    // nodes on slave and master
          int* rx1, int* rx2,    // reordered node on slave and master 
          Real* mortar,      // mortar matrice
          int* seginfo,     // segments information
          int& nseg,        // # of segments
          Real* xsl,         // natural coord for projected node 
          Real* xml,         // natural coord for projected node 
          Real* xsp,         // projection of slave nodes onto the help surface
          Real* xmp,         // projection of master nodes onto the help surface
          int& ns,	       // # of node on slave surface
          int& nm,	       // # of node on master surface
          int& ng,	       // # of gauss quadrature points
          int& nope1, int& nope2 );
    
    void FORTRAN(mortar_gap) (
     	 int& npair,      //pair number
     	 int& ndm,	     //number of dimensions
     	 int& ndf,	     //number of dofs per node
     	 Real* x,	     //coordinates
     	 Real* u,	     //dofs (displacements)
     	 int* ix1, int* ix2,   // nodes on slave and master
     	 int* rx1,
     	 Real* mortar,     // mortar matrice
     	 Real* normals,    // normal vectors at slave nodes
     	 int& ns,	     // # of node on slave surface
     	 int& nm,	     // # of node on master surface
     	 int* inn1, int* inn2, // node ID mapping to reordered one (slave, master)
     	 Real* kappa,      // make gap function to unit value
     	 int* bcontact,   // whether the slave node is incontact with master surface.
     	 Real* gapvec,     // the gap vector 
     	 Real* gap,	     // gap at the slave nodes
     	 Real* traction,   // traction at the slave nodes
     	 Real& penn,
     	 int& contact_sign,
     	 int& nope1, int& nope2,
     	 int& selem);
}


ContactTDEnfPenalty::ContactTDEnfPenalty( const Real* enf_data, 
					              ContactSearch* Search,
			                              ContactSearch::ContactErrorCode& error,
                                                      bool get_cvars,
                                                      bool get_plot_force)
  : ContactTDEnforcement( enf_data, Search, error, get_cvars, get_plot_force),
    cnfi(NULL),cnei(NULL),penalty_scale(1.0)
{  
  if (error != ContactSearch::NO_ERROR) return;
  num_iterations = 1;
  convergence_tolerance = CONVERGENCE_TOL;
}

ContactTDEnfPenalty::ContactTDEnfPenalty( ContactSearch* Search,
			  		              const Real* restart_data,
						      ContactSearch::ContactErrorCode& error,
                                                      bool get_cvars,
                                                      bool get_plot_force)
  : ContactTDEnforcement( Search, restart_data, error, get_cvars, get_plot_force ),
    cnfi(NULL),cnei(NULL),penalty_scale(1.0)
{  
  if (error != ContactSearch::NO_ERROR) return;
  num_iterations = 1;
}

ContactTDEnfPenalty::~ContactTDEnfPenalty()
{
}

ContactSearch::ContactErrorCode
ContactTDEnfPenalty::Set_Penalty_Scale( double scale )
{
  if( scale > 0 ){
    penalty_scale = scale;
    return ContactSearch::NO_ERROR;
  }
  errors->Add_Error_Message( "Penalty scale must be positive" );
  return ContactSearch::INVALID_DATA;
}

void compute_location_on_face(ContactFace<Real> *face, 
                              Real * face_node_coords, 
                              const Real local_coord_x, 
                              const Real local_coord_y,
                              Real * shape_functions, 
                              Real *global_coords) {
  Real local_coords[2] = {local_coord_x, local_coord_y};
  face->Evaluate_Shape_Functions(local_coords, shape_functions);
  PRECONDITION(face->Nodes_Per_Face() >=  1);
  global_coords[0] = face_node_coords[0] * shape_functions[0];
  global_coords[1] = face_node_coords[1] * shape_functions[0];
  global_coords[2] = face_node_coords[2] * shape_functions[0];
  for(int inode = 1; inode < face->Nodes_Per_Face(); ++inode) {
    global_coords[0] += face_node_coords[inode * 3 + 0] * shape_functions[inode];
    global_coords[1] += face_node_coords[inode * 3 + 1] * shape_functions[inode];
    global_coords[2] += face_node_coords[inode * 3 + 2] * shape_functions[inode];
  }
}

ContactSearch::ContactErrorCode 
ContactTDEnfPenalty::Compute_Contact_Force( Real DT_old, Real DT,
						  Real* Mass,
						  Real* Density,
						  Real* Wavespeed,
						  Real* Force  )
{
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Start_Timer( enforcement_time );
#endif
#if LOCAL_PRINT_FLAG>0
  ContactParOStream& postream = ParOStream();
#endif
#if LOCAL_PRINT_FLAG > 9
  postream << "TDPenalty::Compute_Contact_Force START \n";
  postream.flush();
#endif

  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;

  dt = DT;
  dt_old = DT_old;
  dt2 = 1.0/(0.5*(dt+dt_old)*dt);
  // global scaling of inertia based penalty
  Real penalty = 0.5*penalty_scale*dt2; 

  // Call base class Set_Up function to prepare for enforcement
  error_code = Set_Up();
  error_code = (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, communicator );

  if( error_code ){
    Clean_Up();
    return error_code;
  }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Start_Timer( initialization_time );
#endif

  if(SAVE_CVARS) Initialize_CVARS(); // PJSA even though CVARS are not saved for penalty, at least they should be initialized to zero

  // search data
  //
  // Get and initialize scratch memory 
  //
  NODAL_MASS.       Allocate_Scratch(number_of_total_nodes, 1, ScratchVariable::ZERO_SCRATCH);
  INC_FORCE.        Allocate_Scratch(number_of_total_nodes, 3, ScratchVariable::ZERO_SCRATCH);
  TOTAL_FORCE.      Allocate_Scratch(number_of_total_nodes, 3, ScratchVariable::ZERO_SCRATCH);
  NEW_POSITION.     Allocate_Scratch(number_of_total_nodes, 3, ScratchVariable::ZERO_SCRATCH);

  KIN_CONSTR_VECTOR.Allocate_Scratch(number_of_total_nodes, 3);
  NUM_KIN_CONSTR.   Allocate_Scratch(number_of_total_nodes, 1);
  PREDICTED_POS.    Allocate_Scratch(number_of_total_nodes, 3);
  CURRENT_POS.      Allocate_Scratch(number_of_total_nodes, 3);
  //
  //  Copy permanante nodal variables to scratch variables
  //
  {
    VariableHandle var_handles[4];
    var_handles[0] = topology->Variable_Handle( ContactTopology::Current_Position  );
    var_handles[1] = topology->Variable_Handle( ContactTopology::Predicted_Position  );
    var_handles[2] = topology->Variable_Handle( ContactTopology::Kin_Constr_Vector );
    var_handles[3] = topology->Variable_Handle( ContactTopology::Num_Kin_Constr );
    ScratchVariable* scratch_arrays[4];  
    scratch_arrays[0] = &CURRENT_POS;
    scratch_arrays[1] = &PREDICTED_POS;
    scratch_arrays[2] = &KIN_CONSTR_VECTOR;
    scratch_arrays[3] = &NUM_KIN_CONSTR;
    Copy_Variable_to_Scratch(4, var_handles, scratch_arrays);
  }



  // Initialize the scratch memory 
  std::memset( Force,0, number_host_code_nodes*sizeof(Real)*NDIM );
  NEW_POSITION.Duplicate_Scratch(PREDICTED_POS);
  
  // copy host code mass array to local scratch
  Real* host_arrays[1];
  host_arrays[0] = Mass;
  ScratchVariable* scratch_arrays[2];
  scratch_arrays[0] = &NODAL_MASS;
  Copy_Host_Scalar_Arrays_to_Node_Scratch( 1, host_arrays, scratch_arrays );

  // mass & predicted positions of the phantom nodes
  // REJ : communicating the NEW_POSITION=PREDICTED POSITION isn't necessary?
  scratch_arrays[0] = &NODAL_MASS;
  scratch_arrays[1] = &NEW_POSITION;
#ifndef CONTACT_NO_MPI
  if(contact_number_of_processors( communicator ) > 1 ){
    contact_import_scratch_vars(communicator,
                                *Node_AsymComm,
                                *(search->Get_Comm_Buffer()),
                                2,
                                scratch_arrays);
  }
#endif
  //  Compact out duplicate or invalid interactions
  //Compact_Node_Entity_List();
  //  Set the fast lookup methods for the entity list
  node_entity_list.Finalize();

  // Build lists grouping the nodes with single and mulitple interactions
  error_code = (ContactSearch::ContactErrorCode) contact_global_error_check( error_code, communicator );

#if !defined (CONTACT_NO_MPI) & defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( initialization_time );
#endif

  if(  error_code ){
    Clean_Up();
    return error_code;
  }

  for(int i=0 ; i<number_of_total_nodes ; ++i){
    Real* f_tot = TOTAL_FORCE.Get_Scratch(i);
    for(int j=0 ; j<NDIM ; ++j){ f_tot[j] = 0.0;}
  }
  
#ifdef CONTACT_TD_FACE_FACE_ENF
  int num_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  int num_faces = topology->Number_of_Faces();
  ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
#endif

  // set to a value that won't trigger an exit on the 1st iteration
  int iteration;
  for( iteration=0 ; iteration<num_iterations ; ++iteration ){

    if (contactPrintFlag > 1 && contact_processor_number(communicator)==0) {
      std::cout << "  Iteration "<<iteration<<"\n";
    }

    // Initialize for iteration
    INC_FORCE.Zero_Scratch();
#if LOCAL_PRINT_FLAG > 9
    postream << "TDPenalty : initialization done\n";
    postream.flush();
#endif
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( penalty_time );
#endif

#ifdef CONTACT_TD_FACE_FACE_ENF
    { // begin FACE_FACE include
    // *NOTE* - there is some hardcoding in here for quads
    for(int i = 0; i < num_faces; ++i) {
      Faces[i]->temp_tag=0;
    }
    int num_slave_faces  = 0;
    int num_master_faces = 0;
    int global_num_face_face_interactions = 0;
    int global_num_vertices = 0;
    for(int i = 0; i < num_faces; ++i) {
      ContactFace<Real> *face = Faces[i];
      int num_face_face_interactions = face->Number_FaceFace_Interactions();
      global_num_face_face_interactions += face->Number_FaceFace_Interactions();
      for(int j = 0; j < num_face_face_interactions; ++j) {
        ContactFaceFaceInteraction<Real> *ffi = face->Get_FaceFace_Interaction(j);
        ContactFace<Real> *master_face = ffi->MasterFace();
        ContactFace<Real> *slave_face  = ffi->SlaveFace();
        POSTCONDITION(face==slave_face);
	if (master_face->temp_tag == 0) {
	  master_face->temp_tag = -1;
	  ++num_master_faces;
          for (int k=0; k<master_face->Nodes_Per_Face(); ++k) {
            master_face->Node(k)->temp_tag  = -1;
            master_face->Node(k)->temp_tag1 = -1;
          }
	}
	if (slave_face->temp_tag == 0) {
	  slave_face->temp_tag = -2;
	  ++num_slave_faces;
          for (int k=0; k<slave_face->Nodes_Per_Face(); ++k) {
            slave_face->Node(k)->temp_tag  = -2;
            slave_face->Node(k)->temp_tag1 = -2;
          }
	}
      }
    }
    int* ix1 = new int[4*num_slave_faces];
    int* ix2 = new int[4*num_master_faces];     //nodes on slave and master
    int* rx1 = new int[4*num_slave_faces];
    int* rx2 = new int[4*num_master_faces];     //nodes on slave and master
    num_slave_faces  = 0;
    num_master_faces = 0;
    int num_slave_face_nodes  = 0;
    int num_master_face_nodes = 0;
    for(int i = 0; i < num_faces; ++i) {
      ContactFace<Real> *face = Faces[i];
      if (face->temp_tag == -1) {
        for (int k=0; k<face->Nodes_Per_Face(); ++k) {
          if (face->Node(k)->temp_tag==-1) {
            face->Node(k)->temp_tag = ++num_master_face_nodes;
          }
          ix2[4*num_master_faces+k] = face->Node(k)->EnfArrayIndex()+1; //FORTRAN indexing
          rx2[4*num_master_faces+k] = face->Node(k)->temp_tag;
        }
        num_master_faces++;
        face->temp_tag = num_master_faces;
      }
      if (face->temp_tag == -2) {
        for (int k=0; k<face->Nodes_Per_Face(); ++k) {
          if (face->Node(k)->temp_tag==-2) {
            face->Node(k)->temp_tag = ++num_slave_face_nodes;
          }
          ix1[4*num_slave_faces+k] = face->Node(k)->EnfArrayIndex()+1; //FORTRAN indexing
          rx1[4*num_slave_faces+k] = face->Node(k)->temp_tag;
        }
        num_slave_faces++;
        face->temp_tag = num_slave_faces;
      }
    }
    
    int   ns      = num_slave_face_nodes;
    int   nm      = num_master_face_nodes;
    int*  inn1    = new int [num_slave_face_nodes];
    int*  inn2    = new int [num_master_face_nodes];
    Real* normals = new Real [num_slave_face_nodes*3];
    num_slave_face_nodes = 0;
    num_master_face_nodes = 0;
    for (int i=0; i<number_of_total_nodes; ++i) {
      if (Nodes[i]->temp_tag1 == -1) {
	inn2[Nodes[i]->temp_tag-1] = i+1;
	num_master_face_nodes++;
      }
      if (Nodes[i]->temp_tag1 == -2)  {
	inn1[Nodes[i]->temp_tag-1] = i+1;
	Real* node_normal = Nodes[i]->Variable(ContactTopology::Node_Normal);
	normals[3*num_slave_face_nodes  ] = node_normal[0];
	normals[3*num_slave_face_nodes+1] = node_normal[1];
	normals[3*num_slave_face_nodes+2] = node_normal[2];
	num_slave_face_nodes++;
      }
    }
    
    int  n = 0;
    int* seginfo = new int[27*global_num_face_face_interactions]; //segments information
    for(int i = 0; i < num_faces; ++i) {
      ContactFace<Real> *face = Faces[i];
      int num_face_face_interactions = face->Number_FaceFace_Interactions();
      for(int j = 0; j < num_face_face_interactions; ++j) {
        ContactFaceFaceInteraction<Real> *ffi = face->Get_FaceFace_Interaction(j);
        ContactFace<Real> *master_face = ffi->MasterFace();
        ContactFace<Real> *slave_face  = ffi->SlaveFace();
	seginfo[27*n  ] = slave_face->temp_tag;
	seginfo[27*n+1] = master_face->temp_tag;
	seginfo[27*n+2] = ffi->NumEdges();
	++n;
      }
    }
    int npair;         // not really used but needed in the arg list
    int contact_sign;  //output
    int selem      = num_slave_faces;
    int ng         = 3;  // note: this needs to be input from the host-code
    Real penn      = 0.0;  // note: this needs to be input from the host-code
    int ndm        = 3;
    int ndf        = 3;
    int nope1      = 4;
    int nope2      = 4;
    int nseg       = global_num_face_face_interactions;
    int num_pair   = global_num_face_face_interactions;
    Real* x        = CURRENT_POS.Get_Data();      // coordinates
    Real* u        = new Real[3*number_of_total_nodes];      // dofs (displacements)
    Real* xsl      = new Real[16*global_num_face_face_interactions];    // 
    Real* xml      = new Real[16*global_num_face_face_interactions];    //
    Real* xsp      = new Real[12*global_num_face_face_interactions];    //projection of slave nodes onto the help surface
    Real* xmp      = new Real[12*global_num_face_face_interactions];    //projection of master nodes onto the help surface
    Real* mortar   = new Real[ns*(ns+nm)];
    Real* kappa    = new Real[ns];
    int* bcontact  = new int[ns];
    Real* gapvec   = new Real [3*ns];
    Real* gap      = new Real [ns];
    Real* traction = new Real [ns];
    memset(u, 0, number_of_total_nodes*3*sizeof(Real));
    
    int num_seg = 0;
    for(int i = 0; i < num_faces; ++i) {
      ContactFace<Real> *face = Faces[i];
      int num_face_face_interactions = face->Number_FaceFace_Interactions();
      if(num_face_face_interactions == 0) continue;
      Real* normal = face->Variable(ContactTopology::Face_Normal);
      //
      // Each face may overlap any number of other faces, loop over all the 
      // defined face-face overlaps and process each independantly
      //
      for(int j = 0; j < num_face_face_interactions; ++j) {
        ContactFaceFaceInteraction<Real> *ffi = face->Get_FaceFace_Interaction(j);
	//
	//  Extract pointers to the two faces involved in the interaction.  Save basic information such as the
	//  number of nodes in each face and the number of vertexes in the overlapping polygon.
	//
        Real master_shape_functions[MAX_NODES_PER_FACE];
        Real slave_shape_functions[MAX_NODES_PER_FACE];
        ContactFace<Real> *master_face = ffi->MasterFace();
        ContactFace<Real> *slave_face  = ffi->SlaveFace();
        int num_master_nodes     = master_face->Nodes_Per_Face();
        int num_slave_nodes      = slave_face->Nodes_Per_Face();
        int num_vertex           = ffi->NumEdges();
	//
	// Extract face nodal coordinates, this will be looped over several 
        // times so are saved to temporary arrays.  Also, calculate the 
        // projection of the slave and master face nodes to the flat plane.
	//
	Real xso[3] = {0.0, 0.0, 0.0};
        Real master_coords[MAX_NODES_PER_FACE * 3];
        Real slave_coords[MAX_NODES_PER_FACE * 3];
        for(int inode = 0; inode < master_face->Nodes_Per_Face(); ++inode) {
          Real *x_cur = CURRENT_POS.Get_Scratch(master_face->Node(inode));
          master_coords[inode*3+0] = x_cur[0];
          master_coords[inode*3+1] = x_cur[1];
          master_coords[inode*3+2] = x_cur[2];
	}

        for(int inode = 0; inode < slave_face->Nodes_Per_Face(); ++inode) {
          Real *x_cur = CURRENT_POS.Get_Scratch(slave_face->Node(inode));
          slave_coords[inode*3+0] = x_cur[0];
          slave_coords[inode*3+1] = x_cur[1];
          slave_coords[inode*3+2] = x_cur[2];
	  xso[0] += x_cur[0];
	  xso[1] += x_cur[1];
	  xso[2] += x_cur[2];
	}
	xso[0] *= 0.25;
	xso[1] *= 0.25;
	xso[2] *= 0.25;

        Real xdotn1 = (slave_coords[0*3+0]-xso[0])*normal[0]
                     +(slave_coords[0*3+1]-xso[1])*normal[1]
                     +(slave_coords[0*3+2]-xso[2])*normal[2];

        Real xdotn2 = (slave_coords[1*3+0]-xso[0])*normal[0]
                     +(slave_coords[1*3+1]-xso[1])*normal[1]
                     +(slave_coords[1*3+2]-xso[2])*normal[2];

        Real xdotn3 = (slave_coords[2*3+0]-xso[0])*normal[0]
                     +(slave_coords[2*3+1]-xso[1])*normal[1]
                     +(slave_coords[2*3+2]-xso[2])*normal[2];

        Real xdotn4 = (slave_coords[3*3+0]-xso[0])*normal[0]
                     +(slave_coords[3*3+1]-xso[1])*normal[1]
                     +(slave_coords[3*3+2]-xso[2])*normal[2];

        Real xdotn5 = (master_coords[0*3+0]-xso[0])*normal[0]
                     +(master_coords[0*3+1]-xso[1])*normal[1]
                     +(master_coords[0*3+2]-xso[2])*normal[2];

        Real xdotn6 = (master_coords[1*3+0]-xso[0])*normal[0]
                     +(master_coords[1*3+1]-xso[1])*normal[1]
                     +(master_coords[1*3+2]-xso[2])*normal[2];

        Real xdotn7 = (master_coords[2*3+0]-xso[0])*normal[0]
                     +(master_coords[2*3+1]-xso[1])*normal[1]
                     +(master_coords[2*3+2]-xso[2])*normal[2];

        Real xdotn8 = (master_coords[3*3+0]-xso[0])*normal[0]
                     +(master_coords[3*3+1]-xso[1])*normal[1]
                     +(master_coords[3*3+2]-xso[2])*normal[2];

        for (int ii=0; ii<3; ++ii) {
          xsp[num_seg*12+0*3+ii]=slave_coords[0*3+ii]-xdotn1*normal[ii];
          xsp[num_seg*12+1*3+ii]=slave_coords[1*3+ii]-xdotn2*normal[ii];
          xsp[num_seg*12+2*3+ii]=slave_coords[2*3+ii]-xdotn3*normal[ii];
          xsp[num_seg*12+3*3+ii]=slave_coords[3*3+ii]-xdotn4*normal[ii];

          xmp[num_seg*12+0*3+ii]=master_coords[0*3+ii]-xdotn5*normal[ii];
          xmp[num_seg*12+1*3+ii]=master_coords[1*3+ii]-xdotn6*normal[ii];
          xmp[num_seg*12+2*3+ii]=master_coords[2*3+ii]-xdotn7*normal[ii];
          xmp[num_seg*12+3*3+ii]=master_coords[3*3+ii]-xdotn8*normal[ii];
        }
	
	//
	// Compute the area of contact on the master face, note, this assumes the 
        // contact patch is convex, I belive this will always be the case for a 
        // well posed problem, but that is something that should probably be double 
        // checked.....  
	//
	//  The area is computed by summing the areas of sub triangles using cross products
	//
        Real overlap_area = 0.0;

        Real sub_tri_origin[3];
        Real sub_tri_side_a[3];
        Real sub_tri_side_b[3];
        for(int ivert = 0; ivert < num_vertex; ++ivert) {
          ContactFaceFaceVertex<Real> *vertex = ffi->Get_Vertex(ivert);
	  //
	  // load local coordinates for ecah interaction vertex for both
          // the slave and master faces
	  //
	  xsl[num_seg*16+2*ivert]   = vertex->slave_x;
	  xsl[num_seg*16+2*ivert+1] = vertex->slave_y;
	  xml[num_seg*16+2*ivert]   = vertex->master_x;
	  xml[num_seg*16+2*ivert+1] = vertex->master_y;

	  //
	  //  Compute the true location of the current vertex on the master face
	  //
          Real master_global_coords[3];
          compute_location_on_face(master_face, master_coords, 
                                   vertex->master_x, vertex->master_y, 
                                   master_shape_functions,
                                   master_global_coords);
				   
	  if(ivert == 0) {
            sub_tri_origin[0] = master_global_coords[0];
            sub_tri_origin[1] = master_global_coords[1];
            sub_tri_origin[2] = master_global_coords[2];
	  } else if(ivert == 1) {
            sub_tri_side_a[0] = master_global_coords[0] - sub_tri_origin[0];
            sub_tri_side_a[1] = master_global_coords[1] - sub_tri_origin[1];
            sub_tri_side_a[2] = master_global_coords[2] - sub_tri_origin[2];
	  } else {
            sub_tri_side_b[0] = master_global_coords[0] - sub_tri_origin[0];
            sub_tri_side_b[1] = master_global_coords[1] - sub_tri_origin[1];
            sub_tri_side_b[2] = master_global_coords[2] - sub_tri_origin[2];
	    //
	    // sub_tri_side_vec a and b are two sides of the current sub triangle.  
            // Take half the magnitudeof the cross product to compute the current 
            // sub triange area.
	    //
            Real cross_prod[3] = {sub_tri_side_a[1]*sub_tri_side_b[2] - 
                                  sub_tri_side_a[2] * sub_tri_side_b[1],
                                  sub_tri_side_a[0]*sub_tri_side_b[2] - 
                                  sub_tri_side_a[2] * sub_tri_side_b[0],
				  sub_tri_side_a[0]*sub_tri_side_b[1] - 
                                  sub_tri_side_a[1] * sub_tri_side_b[0]};
            overlap_area += 0.5*std::sqrt(cross_prod[0]*cross_prod[0]+
                                          cross_prod[1]*cross_prod[1]+
					  cross_prod[2]*cross_prod[2]);
            sub_tri_side_a[0] = sub_tri_side_b[0];
            sub_tri_side_a[1] = sub_tri_side_b[1];
            sub_tri_side_a[2] = sub_tri_side_b[2];
	  }
	}
	num_seg++;
	//
	//  Compute kinematic penalty factors
	//
	Real  kp = enforcement_data.Get_Data( KINEMATIC_PARTITION, 
                                              master_face->Entity_Key(), 
                                              slave_face->Entity_Key());
        Real pen = kp * dt2 * overlap_area / num_vertex;
	penn = std::max(penn, kp);
      }
    }
    if(iteration == 0) {
      FORTRAN(kappa_3d)(
	  npair,      // pair number
          ndm,        // number of dimensions
          ndf,        // number of dofs per node
          x,          // coordinates
          u,          // dofs (displacements)
          ix1,        // original vertice on slave face
          rx1,        // reordered vertice on slave face
          kappa,      // scale factor for mortar matrix
          ns,         // # of node on slave surface
          nope1,      // number of vertice on slave element
          selem       // number od slave elements
	  );
      
    }
    FORTRAN(mortar_matrix)(
          npair,       // pair number, doesn't seem to be used
          ndm,         // number of dimensions
          ndf,         // number of dofs per node
          x,	       // coordinates
          u,	       // dofs (displacements)
          ix1, ix2,    // nodes on slave and master
          rx1, rx2,    // reordered node on slave and master 
          mortar,      // mortar matrice
          seginfo,     // segments information
          nseg,        // # of segments
          xsl,         // natural coord for projected node 
          xml,         // natural coord for projected node 
          xsp,         // projection of slave nodes onto the help surface
          xmp,         // projection of master nodes onto the help surface
          ns,	       // # of node on slave surface
          nm,	       // # of node on master surface
          ng,	       // # of gauss quadrature points
          nope1, nope2 );
    FORTRAN(mortar_gap) (
     	 npair,      //pair number
     	 ndm,	     //number of dimensions
     	 ndf,	     //number of dofs per node
     	 x,	     //coordinates
     	 u,	     //dofs (displacements)
     	 ix1, ix2,   // nodes on slave and master
     	 rx1,
     	 mortar,     // mortar matrice
     	 normals,    // normal vectors at slave nodes
     	 ns,	     // # of node on slave surface
     	 nm,	     // # of node on master surface
     	 inn1, inn2, // node ID mapping to reordered one (slave, master)
     	 kappa,      // make gap function to unit value
     	 bcontact,   // whether the slave node is incontact with master surface.
     	 gapvec,     // the gap vector 
     	 gap,	     // gap at the slave nodes
     	 traction,   // traction at the slave nodes
     	 penn,
     	 contact_sign,
     	 nope1, nope2,
     	 selem);
    
    // now do a vector-matrix multiply to compute the forces, f = traction*mortar
    
    // slave nodes
    for (int i=0; i<ns; ++i) {
      PRECONDITION(Nodes[inn1[i]-1]);
      Real* f_K = INC_FORCE.Get_Scratch(Nodes[inn1[i]-1]);
      for (int j=0; j<ns; ++j) {
        int index = i*ns+j;
        f_K[0] += mortar[index]*traction[j]*normals[3*j+0];
        f_K[1] += mortar[index]*traction[j]*normals[3*j+1];
        f_K[2] += mortar[index]*traction[j]*normals[3*j+2];
      }
    }
    
    // master nodes
    for (int i=0; i<nm; ++i) {
      PRECONDITION(Nodes[inn2[i]-1]);
      Real* f_K = INC_FORCE.Get_Scratch(Nodes[inn2[i]-1]);
      for (int j=0; j<ns; ++j) {
        int index = (i+ns)*ns+j;
        f_K[0] -= mortar[index]*traction[j]*normals[3*j+0];
        f_K[1] -= mortar[index]*traction[j]*normals[3*j+1];
        f_K[2] -= mortar[index]*traction[j]*normals[3*j+2];
      }
    }
	 
    } // end FACE_FACE include
#endif

    //
    // Compute penalty forces for node-face interactions
    //
    for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {

      const int node_index = cnei_group.Get_Constrained_Node_Index();
      for( int i=0 ; i< cnei_group.Num_Interactions() ; ++i ){
	cnei = cnei_group.Get_Interaction(i);
        if(cnei != NULL){
	  int e_key = cnei->Get_Entity_Key();
	  int n_key = cnei->Get_Node_Key();
	  Real  kp = enforcement_data.Get_Data( KINEMATIC_PARTITION, e_key, n_key );
          cnfi = dynamic_cast<ContactNodeFaceInteraction*> (cnei);
	  if (cnfi && kp > 0.0) {
	    Real  m_sn = *(NODAL_MASS.Get_Scratch(node_index));
	    Real pen = kp*penalty*m_sn;
	    Real* nm = cnei->Get_Normal();

            Real shape_functions[MAX_NODES_PER_FACE];
	    Real uj_cur[3] = {0.0,0.0,0.0};
	    Real uj_new[3] = {0.0,0.0,0.0};
	    Real* x_sn_cur = CURRENT_POS.Get_Scratch(node_index);
	    Real* x_sn_new = NEW_POSITION.Get_Scratch(node_index);
	    for(int j=0 ; j<NDIM ; ++j){
	      uj_cur[j] = x_sn_cur[j];
	      uj_new[j] = x_sn_new[j];
	    }

	    // compute contact point
	    ContactFace<Real>* face = cnfi->Face();
            const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(i);
	    PRECONDITION( num_face_nodes <= MAX_NODES_PER_FACE );
	    // face distribution factors
	    face->Evaluate_Shape_Functions(cnfi->Get_Coordinates(), shape_functions );
	    for(int k=0 ; k<num_face_nodes  ; ++k){
	      Real* x_mn_cur = CURRENT_POS.Get_Scratch(face->Node(k));
	      Real* x_mn_new = NEW_POSITION.Get_Scratch(face->Node(k));
	      for(int j=0 ; j<NDIM ; ++j){
	        uj_cur[j] -=  shape_functions[k]*x_mn_cur[j] ;
	        uj_new[j] -=  shape_functions[k]*x_mn_new[j] ;
	      }
	    }

            Real f_sn[3] = {0.0,0.0,0.0};
            ContactEnfModel *model      = TDEnfModel(cnfi);
            ContactTDTied   *tied_model = dynamic_cast<ContactTDTied*>(model);
            if (tied_model != NULL) { // tied
              Real gap0 = cnei->Get_Gap_Initial();
              Real n_dir[3] = {0,0,0};
              Real* x_face[MAX_NODES_PER_FACE];
              PRECONDITION( MAX_NODES_PER_FACE>= num_face_nodes);
              for(int ii=0 ; ii<num_face_nodes ; ++ii ){
                const int node_face_index = cnei_group.Get_Face_Node_Index(i, ii);
                 x_face[ii] = NEW_POSITION.Get_Scratch(node_face_index);
              }
              cnfi->Face()->Compute_Normal(x_face,cnfi->Get_Coordinates(),n_dir);
              for(int j=0 ; j<NDIM ; ++j ) uj_new[j] -=  gap0*n_dir[j]; 
	      for(int j=0 ; j<NDIM ; ++j ) f_sn[j] = -pen*uj_new[j];
#if LOCAL_PRINT_FLAG > 9
	      postream << "node: " << node_index
		<< " uj: " << uj_new[0] << " " << uj_new[1] << " " << uj_new[2] 
		<< "\n";
#endif
            }
            else { // all others -> frictionless 

	      // compute penalty acceleration
	      Real gap = 0.0;
	      for (int j=0 ; j<NDIM ; ++j) gap += uj_new[j]*nm[j];
	      Real gapdot = gap;
	      for (int j=0 ; j<NDIM ; ++j) gapdot -= uj_cur[j]*nm[j];
#if LOCAL_PRINT_FLAG > 9
	      postream << "node: " << node_index
		     << " gap: "<< gap  << " gapdot:" << gapdot << "\n";
#endif
              if (gap < 0.0 && gapdot < 0.0)  {
	        // Penalty force
	        for(int j=0 ; j<NDIM ; ++j) f_sn[j] = -pen*gap*nm[j];
	      } // g < 0, g_dot < 0
            } // tied or else
#if LOCAL_PRINT_FLAG > 9
	    postream << "node: " << node_index
		<< " f_pen: " << f_sn[0] << " " << f_sn[1] << " " << f_sn[2] 
		<< "\n";
#endif
	    // Assemble to nodes
	    Real* f_K = INC_FORCE.Get_Scratch(node_index);
	    for(int j=0 ; j<NDIM ; ++j) f_K[j] += f_sn[j];
	    // Assemble to face nodes
	    for(int k=0 ; k<face->Nodes_Per_Face() ; ++k){
	      f_K = INC_FORCE.Get_Scratch(face->Node(k));
	      // Force on the face is equal & opposite the slave node force
	      for(int j=0 ; j<NDIM ; ++j){
                f_K[j] -= shape_functions[k]*f_sn[j];
              }
	    } 
	  } // cnfi &&  kp > 0.0
        } // cnei exists
      } // max interactions
    }

    // SWAPADD the assembled_force
    ScratchVariable* vars[1];
    vars[0] = &INC_FORCE;
    swapadd_node_scratch_vars( 1, vars);

    if( topology->Have_Shells() ) Assemble_Shell_Forces(); // PJSA

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( penalty_time );
    Timer().Start_Timer( update_config_time );
#endif

    // make force compatible with prescibed boundary conditions
    Make_Scratch_Vector_Consistent_with_BCs( INC_FORCE );

#if LOCAL_PRINT_FLAG > 9
    postream << "TDPenalty : contact force done\n";
    postream.flush();
#endif

    // Compute total force & new coordinates 
    for( int inode=0 ; inode<number_of_nodes ; inode++ ){
      Real* x_new =   NEW_POSITION.Get_Scratch(inode);
      Real* f_inc =   INC_FORCE.   Get_Scratch(inode);
      Real m_sn   = *(NODAL_MASS.  Get_Scratch(inode));
      Real* f_tot =   TOTAL_FORCE. Get_Scratch(inode);
      Real scale = 1.0/(m_sn*dt2);
      for(int j=0 ; j<NDIM ; ++j){ 
         f_tot[j] += f_inc[j];
         x_new[j] += f_inc[j]*scale; 
      }
    }

    // We need to get the predicted positions of the phantom nodes
    scratch_arrays[0] = &NEW_POSITION;
#ifndef CONTACT_NO_MPI
    if(contact_number_of_processors( communicator ) > 1 ){
      contact_import_scratch_vars(communicator,
                                  *Node_AsymComm,
                                  *(search->Get_Comm_Buffer()),
                                  1,
                                  scratch_arrays);
    }
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( update_config_time );
#endif

    // calculate norm of the incremental force (multiplier)
    if ( convergence_tolerance > 0.0) {
      Real inc_force_norm = 0.0;
      for(int i=0 ; i<number_of_nodes ; ++i ){
	Real* f_inc = INC_FORCE.Get_Scratch(i);
	for(int j=0 ; j<NDIM ; ++j ){ inc_force_norm += f_inc[j]*f_inc[j]; }
      }
      global_inc_force_norm
	= std::sqrt(contact_global_sum(inc_force_norm, communicator));
      if (iteration == 0) {
	initial_global_inc_force_norm = global_inc_force_norm;

        if (contactPrintFlag > 1 && contact_processor_number(communicator)==0) {
          std::cout << "    alt convergence tol = "<<convergence_tolerance*initial_global_inc_force_norm<<"\n";
        }
      }
#if LOCAL_PRINT_FLAG > 0
      postream << iteration << " >> convg. norm.: "
	       << global_inc_force_norm << "\n";
      postream.flush();
#endif

      if (contactPrintFlag > 1 && contact_processor_number(communicator)==0) {
        std::cout << "    convergence norm    = "<<global_inc_force_norm<<"\n";
      }

      if( global_inc_force_norm
          < convergence_tolerance*initial_global_inc_force_norm
	  || global_inc_force_norm < convergence_tolerance ) break;
    }
#if LOCAL_PRINT_FLAG > 3 
  postream << "\n FORCES____iteration = " << iteration << "\n";
  for(int i=0 ; i<number_of_nodes ; ++i ){
    ContactNode<Real>* node = enforcement_node_list[i];
    int ExoID = node->Exodus_ID();
    Real* f_tot =   TOTAL_FORCE. Get_Scratch(i);
    postream << ExoID << " f_tot: "
	     << f_tot[0] <<" "<< f_tot[1] <<" "<<f_tot[2] << "\n";
  }
  postream << "_______________________________________\n";
  postream.flush();
#endif
  } // end of iteration loop

  if(iteration == num_iterations && num_iterations > 1 && global_inc_force_norm > std::max(initial_global_inc_force_norm,FAIL_TOL)) {
    if(contact_processor_number(communicator)==0)
      std::cerr << "\rerror: contact enforcement failed (inc. force = " << global_inc_force_norm
                << ", target = " << std::max(convergence_tolerance*initial_global_inc_force_norm, convergence_tolerance) << ")\n";
    if(global_inc_force_norm > initial_global_inc_force_norm) {
      // Set final updated position to the current predicted positiion
      NEW_POSITION.Duplicate_Scratch(PREDICTED_POS);
    }
  }

  Compute_Extra_Ghosting_Length();

  // for Dynamic_Search
  for(int i=0 ; i<number_of_nodes ; ++i){
    Real* x_i =   PREDICTED_POS.Get_Scratch(i);
    Real* x_f =   NEW_POSITION.Get_Scratch(i);
    for(int j=0;j<NDIM; ++j){
      disp_old_aug[j*number_of_nodes+i] = x_f[j]-x_i[j];
    }
  }

  scratch_arrays[0] = &TOTAL_FORCE;
  Copy_Node_Vector_Scratch_to_Host_Arrays( 1, &Force, scratch_arrays );

  if( topology->Have_Shells() ) Store_Shell_Final_Lofted_Positions(); // PJSA

#if LOCAL_PRINT_FLAG > 3 
  postream << "\n FORCES___total_iterations = " << iteration << "\n";
  double f_sum[3] = {0.0,0.0,0.0};
  double f_abs[3] = {0.0,0.0,0.0};
  double m_sum = 0.0;
  for(int i=0 ; i<number_of_nodes ; ++i ){
    ContactNode<Real>* node = enforcement_node_list[i];
    int ExoID = node->Exodus_ID();
    Real* f_tot =   TOTAL_FORCE. Get_Scratch(i);
    postream << ExoID << " f_tot "
	     << f_tot[0] <<" "<< f_tot[1] <<" "<<f_tot[2] << "\n";
    for(int j=0;j<NDIM; ++j){ f_sum[j] += f_tot[j]; }
    for(int j=0;j<NDIM; ++j){ f_abs[j] += std::fabs(f_tot[j]); }
    Real m   = *(NODAL_MASS.  Get_Scratch(i));
    m_sum += m;
  }
  postream << "=======================================\n";
  postream << iteration << " sum f   "
	     << f_sum[0] <<" "<< f_sum[1] <<" "<<f_sum[2] << "\n";
  postream << iteration << " sum |f| "
	     << f_abs[0] <<" "<< f_abs[1] <<" "<<f_abs[2] << "\n";
  postream << iteration << " sum m " << m_sum << "\n";
  postream << "_______________________________________\n";
  postream.flush();
#endif

  // 2DO ifdef this?
  // for regression testing
  if(CALC_PLOT_FORCE) {
    for(int i=0 ; i<number_of_nodes; ++i){
      Real* f_t = TOTAL_FORCE.Get_Scratch(i);
      for(int j=0 ; j<NDIM ; ++j) plot_force[j*number_of_nodes+i] = f_t[j];
    }
  }

  // Compute Old timestep multiplier for augmented search
  time_mult_old_aug = dt*(0.5*(dt+dt_old));

  // Release scratch memory & Call the base class clean up function
  TDPenalty_Release_Scratch();
  Clean_Up();

#if !defined(CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( enforcement_time );
#endif
#if LOCAL_PRINT_FLAG > 9
  postream << "TDPenalty::Compute_Contact_Force FINISH \n";
  postream.flush();
#endif

  return (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, communicator );
}

void ContactTDEnfPenalty::TDPenalty_Release_Scratch(void)
{
  NODAL_MASS.       Clear_Scratch();
  INC_FORCE.        Clear_Scratch();
  TOTAL_FORCE.      Clear_Scratch();
  NEW_POSITION.     Clear_Scratch();
  KIN_CONSTR_VECTOR.Clear_Scratch();
  NUM_KIN_CONSTR.   Clear_Scratch();
  PREDICTED_POS.    Clear_Scratch();
  CURRENT_POS.      Clear_Scratch();
  
  node_entity_list.Clear();
}
