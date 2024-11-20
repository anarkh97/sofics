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


#include "ContactErrors.h"
#include "ContactTopology.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactMPCs.h"
#include "ContactAsymComm.h"
#include "ContactSymComm.h"
#include <cstring>
#include <cstdio>
#include <cmath>

ContactMPCs::ContactMPCs( const Real* enf_data,
			  ContactSearch* Search,
			  ContactSearch::ContactErrorCode& error )
  : ContactEnforcement( error, Search, ContactEnforcement::TiedKinematics,
			NSIZED, enf_data, false )
{
  if (error != ContactSearch::NO_ERROR) return;
  // plot vars
  number_global_plot_vars = 0;
  number_element_plot_vars = 0;
  number_nodal_plot_vars = 0;

#ifndef CONTACT_NO_MPI
  node_asymcomm = new ContactAsymComm( *topology->Node_Sym_Comm() );
#endif

  if( number_of_nodes ){
    slave_mpcs = new Contact_MPC_Equation*[3*number_of_nodes];
    for( int i=0 ; i<3*number_of_nodes ; ++i) slave_mpcs[i] = NULL;
  } else
    slave_mpcs = NULL;

  num_face_mpcs = 0;
  face_mpcs = NULL;

  error = ContactSearch::NO_ERROR;
}

ContactMPCs::ContactMPCs( ContactSearch* Search,
			  const Real* restart_data,
                          ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::TiedKinematics,
                        restart_data )
{
  if (error != ContactSearch::NO_ERROR) return;
  number_global_plot_vars = 0;
  number_nodal_plot_vars = 0;
  number_element_plot_vars = 0;

#ifndef CONTACT_NO_MPI
  node_asymcomm = new ContactAsymComm( *topology->Node_Sym_Comm() );
#endif

  if( number_of_nodes ){
    slave_mpcs = new Contact_MPC_Equation*[3*number_of_nodes];
    for( int i=0 ; i<3*number_of_nodes ; ++i) slave_mpcs[i] = NULL;
  } else
    slave_mpcs = NULL;

  num_face_mpcs = 0;
  face_mpcs = NULL;

  error = ContactSearch::NO_ERROR;
}


ContactMPCs::~ContactMPCs()
{
  Remove_MPCs();
  delete [] slave_mpcs;
  delete [] face_mpcs;
#ifndef CONTACT_NO_MPI
  delete node_asymcomm;
#endif
}

ContactSearch::ContactErrorCode
ContactMPCs::Extract_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Extract_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Implant_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Implant_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}


ContactSearch::ContactErrorCode
ContactMPCs::Implant_Nodal_Restart_Variable( int n, Real* data )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}
ContactSearch::ContactErrorCode
ContactMPCs::Extract_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Implant_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Extract_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Implant_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Extract_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Implant_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Compute_MPCs( )
{
  int j;
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;
  
  // Clear out old MPCs
  Remove_MPCs();

  // Call base class Set_Up
  error_code = Set_Up();
  if( error_code ) return error_code;
  
  // Define the MPCs on the slave nodes

    for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
         cnfi != NULL;
         cnfi = node_entity_list.Node_Face_Iterator_Next()) {
    ContactNode<Real>* sn_node = cnfi->Node();
    ContactFace<Real>* ms_face = cnfi->Face();
    int sn_index = sn_node->ProcArrayIndex();
    Real shape_functions[MAX_NODES_PER_FACE];
    Real* coords = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
    PRECONDITION(cnfi->Face()->Nodes_Per_Face() <= 8 );
    cnfi->Face()->Evaluate_Shape_Functions( coords, shape_functions );
    for( j=0 ; j<3 ; ++j){
      if( !slave_mpcs[3*sn_index+j] ){
	Contact_MPC_Equation* mpc = new Contact_MPC_Equation;
	slave_mpcs[3*sn_index+j] = mpc;
        mpc->slave_node_proc_id = sn_node->Global_ID().HiInt();
        mpc->slave_node_local_id = sn_node->Global_ID().LoInt();
	mpc->master_face_proc_id = ms_face->Global_ID().HiInt();
	mpc->master_face_local_id = ms_face->Global_ID().LoInt();
	int num_face_nodes = cnfi->Face()->Nodes_Per_Face();
	mpc->num_face_nodes = num_face_nodes;
	for( int k=0 ; k<num_face_nodes ; ++k){
	  ContactNode<Real>* face_node = cnfi->Face()->Node(k);
	  mpc->face_node_proc_id[k] = face_node->Global_ID().HiInt();
	  mpc->face_node_local_id[k] = face_node->Global_ID().LoInt();
	  mpc->face_node_coeff[k] = -shape_functions[k];
	}
	break;
      }
    }
  }
    
#ifndef CONTACT_NO_MPI
  Communicate_MPCs( );
#endif

#ifndef CONTACT_NO_EXODUS_OUTPUT
  Plot_MPC_Vars( );
#endif

  // Call the base class clean up functions
  Clean_Up();

  return (ContactSearch::ContactErrorCode)
    contact_global_error_check( error_code, communicator );
}



void ContactMPCs::Remove_MPCs()
{
  int i;
  for( i=0 ; i<3*number_of_nodes ; ++i) delete slave_mpcs[i];
  for( i=0 ; i<num_face_mpcs ; ++i) delete face_mpcs[i];
}


ContactSearch::ContactErrorCode
ContactMPCs::Update_For_Topology_Change( int new_number_of_nodes,
						   int* old_to_new_map )
{
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactMPCs::Number_of_MPC_Equations(int& num_mpcs)
{
  ContactSearch::ContactErrorCode error = ContactSearch::NO_ERROR;

  num_mpcs = num_face_mpcs;
  for( int i=0 ; i<number_of_nodes ; ++i){
    for( int j=0 ; j<3 ; ++j){
      if( slave_mpcs[3*i+j] ) num_mpcs++;
    }
  }
  return error;
}

ContactSearch::ContactErrorCode
ContactMPCs::Get_MPC_Equations(int num_mpcs,
			       int* snode_pid,
			       int* snode_lid,
			       int* mface_pid,
			       int* mface_lid,
			       int* nface_nodes,
			       int* fnode_pid,
			       int* fnode_lid,
			       Real* fnode_coefs)
{
  int i,j;
  const int max_num_fnodes = 8;
  std::memset( snode_pid, 0, 
	  num_mpcs*sizeof(int));
  std::memset( snode_lid, 0,
	  num_mpcs*sizeof(int));
  std::memset( mface_pid, 0,
	  num_mpcs*sizeof(int));
  std::memset( mface_lid, 0,
	  num_mpcs*sizeof(int));
  std::memset( nface_nodes, 0,
	  num_mpcs*sizeof(int));
  std::memset( fnode_pid, 0,
	  num_mpcs*max_num_fnodes*sizeof(int));
  std::memset( fnode_lid, 0,
	  num_mpcs*max_num_fnodes*sizeof(int));
  std::memset( fnode_coefs, 0,
	  num_mpcs*max_num_fnodes*sizeof(int));

  int mpc_count = 0;
  for( i=0 ; i<num_face_mpcs ; ++i) { //load up face mpcs
    if( face_mpcs[i] ) {
      Contact_MPC_Equation* mpc = face_mpcs[i];
      snode_pid[mpc_count] = mpc->slave_node_proc_id;
      snode_lid[mpc_count] = mpc->slave_node_local_id;
      mface_pid[mpc_count] = mpc->master_face_proc_id;
      mface_lid[mpc_count] = mpc->master_face_local_id;
      int number_of_face_nodes =  mpc->num_face_nodes;
      nface_nodes[mpc_count] = number_of_face_nodes;
      for(int k=0; k<number_of_face_nodes; ++k) {
	fnode_pid[8*mpc_count+k] = mpc->face_node_proc_id[k];
	fnode_lid[8*mpc_count+k] = mpc->face_node_local_id[k];
	fnode_coefs[8*mpc_count+k] = mpc->face_node_coeff[k];
      }
      mpc_count++;
    }
  }
  for( i=0 ; i<number_of_nodes ; ++i) { //load up sn mpcs 
    for( j=0 ; j<3 ; ++j){
      if( slave_mpcs[3*i+j]) {
	Contact_MPC_Equation* mpc = slave_mpcs[3*i+j];
	snode_pid[mpc_count] = mpc->slave_node_proc_id;
	snode_lid[mpc_count] = mpc->slave_node_local_id;
	mface_pid[mpc_count] = mpc->master_face_proc_id;
	mface_lid[mpc_count] = mpc->master_face_local_id;
	int number_of_face_nodes =  mpc->num_face_nodes;
	nface_nodes[mpc_count] = number_of_face_nodes;
	for(int k=0; k<number_of_face_nodes; ++k) {
	  fnode_pid[8*mpc_count+k] = mpc->face_node_proc_id[k];
	  fnode_lid[8*mpc_count+k] = mpc->face_node_local_id[k];
	  fnode_coefs[8*mpc_count+k] = mpc->face_node_coeff[k];
	}
	mpc_count++;
      }
    }
  }

  POSTCONDITION(mpc_count == num_mpcs);

  return ContactSearch::NO_ERROR;
}


void ContactMPCs::Communicate_MPCs( )
{
#ifndef CONTACT_NO_MPI
  int i,j,k;

  // ***********************************************************************
  //
  // Communicate the MPCs to ghost copies of the slave nodes
  //
  // ***********************************************************************
  //
  // At this point, we have the MPC's on the owning processor for each
  // slave node.  What we need to do is to communicate the MPCs to all
  // the "ghost" copies of the slave node.  We have a communication plan
  // in node_symcomm.  We are back to the same old issue that I know what
  // I need to send but I don't know what I need to receive.  I do know
  // who might have data for my but I don't know how much (i.e., how many
  // MPCS per shared node).
  //
  // I see two solutions (both with pros and cons)
  //   1) Do an initial communication to give the number of MPCs per node.
  //      Followed by communicating the data.
  //          Pros:  Communicate the minimum amount of data
  //          Cons:  Two communication steps
  //   2) Send a buffer that is long enough to hold 3 MPCs per node even
  //      though some of the data is all zeros and not used.
  //          Pros:  One communication step.
  //          Cons:  Send a lot of unneccesary data.
  //
  // I'm going to use option 1.  My reasoning is that I will face the same
  // issue with communicating to the faces but I can't use option 2 because
  // there is no bound on how many nodes a face may constrain (each of which
  // generates an MPC equation).
  //
  // The steps are then
  //   1) Communicate how many mpc equations I will send for each node
  //   2) Send the mpcs for each node

  if( contact_number_of_processors( communicator ) == 1 ) return;

  // Output the comm plan
  ContactSymComm* node_symcomm = topology->Node_Sym_Comm();

  // Create the asymmetric comm plan I need to push the data
  ContactAsymComm Node_AsymComm_Plan( *node_symcomm );

  // Get the buffers for communicating
  int size_send_buf = Node_AsymComm_Plan.Size_Export()*sizeof(int);
  int size_recv_buf = Node_AsymComm_Plan.Size_Import()*sizeof(int);
  int num_request_handles = Node_AsymComm_Plan.Num_Import_Comm_Partners();
  char* send_buf;
  char* recv_buf;
  RequestHandle* request_handles;
  CommBuffer()->Buffers( size_send_buf, size_recv_buf, num_request_handles,
			 &send_buf, &recv_buf, &request_handles );

  // Put the number of mpcs per node into the buffers & communicate
  int buf_offset = 0;
  int* isb = (int*) send_buf;
  for( i=0 ; i<Node_AsymComm_Plan.Num_Export_Comm_Partners() ; ++i){
    int num_export_to_proc = Node_AsymComm_Plan.Num_Export_to_Proc( i );
    ContactTopologyEntity<Real>** export_entity_list = 
      Node_AsymComm_Plan.Export_Entity_List(i);
    for( j=0 ; j<num_export_to_proc ; ++j){
      int node_index = export_entity_list[j]->ProcArrayIndex();
      int num_mpcs = 0;
      for( k=0 ; k<3 ; ++k){ 
	if( slave_mpcs[3*node_index+k] ) 
	  num_mpcs++;
      }
      isb[buf_offset++] = num_mpcs;
    }
  }

  contact_communicate_packed_buffers( communicator, Node_AsymComm_Plan, 
				      sizeof(int), send_buf, 
				      recv_buf, request_handles );

  // Use the information just received to build a comm plan 
  // a) Build the export part
  int num_export_partners = 0;
  int num_export_entities = 0;
  for( i=0 ; i<Node_AsymComm_Plan.Num_Export_Comm_Partners() ; ++i){
    int num_export_to_proc = Node_AsymComm_Plan.Num_Export_to_Proc(i);
    ContactTopologyEntity<Real>** export_entity_list = 
      Node_AsymComm_Plan.Export_Entity_List(i);
    int num_mpcs = 0;
    for( j=0 ; j<num_export_to_proc ; ++j){
      int node_index = export_entity_list[j]->ProcArrayIndex();
      for( k=0 ; k<3 ; ++k) if( slave_mpcs[3*node_index+k] ) num_mpcs++;
    }
    if( num_mpcs ){
      num_export_partners++;
      num_export_entities += num_mpcs;
    }
  }
  int* export_comm_procs = new int[num_export_partners];
  int* num_export_to_proc = new int[num_export_partners];
  int eoffset = 0;
  int ioffset = 0;
  ContactTopologyEntity<Real>** export_entities = new ContactTopologyEntity<Real>*[num_export_entities];
  for( i=0 ; i<Node_AsymComm_Plan.Num_Export_Comm_Partners() ; ++i){
    int num_export = Node_AsymComm_Plan.Num_Export_to_Proc(i);
    int num_export_to_partner = 0;
    ContactTopologyEntity<Real>** export_entity_list = 
      Node_AsymComm_Plan.Export_Entity_List(i);
    for( j=0 ; j<num_export ; ++j){
      int node_index = export_entity_list[j]->ProcArrayIndex();
      for( k=0 ; k<3 ; ++k){
	if( slave_mpcs[3*node_index+k] ){
	  export_entities[eoffset++] = export_entity_list[j];
	  num_export_to_partner++;
	}
      }
    }
    if( num_export_to_partner ){
      export_comm_procs[ioffset] = Node_AsymComm_Plan.Export_Comm_Proc_ID( i );
      num_export_to_proc[ioffset] = num_export_to_partner;
      ioffset++;
    }    
  }
  // b) Build the import part from the data just received
  int num_import_partners = 0;
  int num_import_entities = 0;
  int* irb = (int*) recv_buf;
  int offset = 0;
  for( i=0 ; i<Node_AsymComm_Plan.Num_Import_Comm_Partners() ; ++i){
    int num_to_partner = 0;
    for( j=0 ; j<Node_AsymComm_Plan.Num_Import_from_Proc( i ) ; ++j)
      num_to_partner += irb[offset++];
    if( num_to_partner ){
      num_import_partners++;
      num_import_entities += num_to_partner;
    }
  }
  int* import_proc_ids = NULL;
  int* num_import_from_proc = NULL;
  ContactTopologyEntity<Real>** import_entity_list = NULL;
  if( num_import_partners ){
    import_proc_ids = new int[num_import_partners];
    num_import_from_proc = new int[num_import_partners];
    import_entity_list = new ContactTopologyEntity<Real>*[num_import_entities];
  }
  offset = 0;
  eoffset = 0;
  ioffset = 0;
  for( i=0 ; i<Node_AsymComm_Plan.Num_Import_Comm_Partners() ; ++i){
    int num_to_partner = 0;
    ContactTopologyEntity<Real>** elist = Node_AsymComm_Plan.Import_Entity_List( i );
    for( j=0 ; j<Node_AsymComm_Plan.Num_Import_from_Proc( i ) ; ++j){
      int num_mpcs = irb[offset++];
      for( k=0 ; k<num_mpcs ; ++k)
	import_entity_list[eoffset++] = elist[j];
      num_to_partner += num_mpcs;
    }
    if( num_to_partner ){
      import_proc_ids[ioffset] = Node_AsymComm_Plan.Import_Comm_Proc_ID( i );
      num_import_from_proc[ioffset] = num_to_partner;
      ioffset++;
    }
  }

  ContactAsymComm Export_NodeMPCs( num_export_partners,
				   num_import_partners,
				   num_export_to_proc,
				   num_import_from_proc,
				   export_comm_procs,
				   import_proc_ids,
				   export_entities,
				   import_entity_list );

  delete [] export_entities;
  delete [] export_comm_procs;
  delete [] num_export_to_proc;
  delete [] num_import_from_proc;
  delete [] import_proc_ids;
  delete [] import_entity_list;
  
  // Get Buffers to send slave node mpcs
  int size_mpc = sizeof(Contact_MPC_Equation);
  size_send_buf = Export_NodeMPCs.Size_Export()*size_mpc;
  size_recv_buf = Export_NodeMPCs.Size_Import()*size_mpc;
  num_request_handles = Export_NodeMPCs.Num_Import_Comm_Partners();
  CommBuffer()->Buffers( size_send_buf, size_recv_buf, num_request_handles,
			 &send_buf, &recv_buf, &request_handles );
  char* sb = send_buf;
  // Load up the send buffers
  for( i=0 ; i<Export_NodeMPCs.Num_Export_Comm_Partners() ; ++i){
    ContactTopologyEntity<Real>** export_entity_list = Export_NodeMPCs.Export_Entity_List(i);
    // Zero temp tag so we know which MPC to load
    for( j=0 ; j<Export_NodeMPCs.Num_Export_to_Proc( i ) ; ++j)
      export_entity_list[j]->temp_tag = 0;
    for( j=0 ; j<Export_NodeMPCs.Num_Export_to_Proc( i ) ; ++j){
      int sn_index = export_entity_list[j]->ProcArrayIndex();
      Contact_MPC_Equation* mpc = 
	slave_mpcs[3*sn_index+export_entity_list[j]->temp_tag];
      export_entity_list[j]->temp_tag++;
      std::memcpy( sb, mpc, size_mpc );
      sb += size_mpc;
    }
  }
      
  contact_communicate_packed_buffers( communicator, Export_NodeMPCs, size_mpc, 
				      send_buf, recv_buf, request_handles );

  // Create MPC equations for my ghost nodes and store them
  char* rb = recv_buf;
  for( i=0 ; i<Export_NodeMPCs.Num_Import_Comm_Partners() ; ++i){
    import_entity_list = Export_NodeMPCs.Import_Entity_List(i);
    for( j=0 ; j<Export_NodeMPCs.Num_Import_from_Proc( i ) ; ++j){
      int sn_index = import_entity_list[j]->ProcArrayIndex();
      for( k=0 ; k<3 ; ++k){
	if( !slave_mpcs[3*sn_index+k] ){
	  Contact_MPC_Equation* mpc = new Contact_MPC_Equation;
	  slave_mpcs[3*sn_index+k] = mpc;
	  std::memcpy( mpc, rb, size_mpc );
	  rb += size_mpc;
	  break;
	}
      }
    }
  }

  // -----------------------------------------------------------------
  // Communicate the MPCs to the processor that has the face
  // -----------------------------------------------------------------
  //
  // Build a comm plan to send data to the faces.
  //   Note: This is the inverse map of Face_AsymComm since we are 
  //         sending back data that we had been importing
  int* import_length_per_proc = 
    new int[Face_AsymComm->Num_Export_Comm_Partners()];
  for( i=0 ; i<Face_AsymComm->Num_Export_Comm_Partners() ; ++i)
    import_length_per_proc[i] = 1;
  int* export_length_per_proc = 
    new int[Face_AsymComm->Num_Import_Comm_Partners()];
  for( i=0 ; i<Face_AsymComm->Num_Import_Comm_Partners() ; ++i)
    export_length_per_proc[i] = 1;

  ContactAsymComm face_asymcomm1( Face_AsymComm->Num_Import_Comm_Partners(),
				  Face_AsymComm->Num_Export_Comm_Partners(),
				  export_length_per_proc,
				  import_length_per_proc,
				  Face_AsymComm->Import_Comm_Proc_IDs(),
				  Face_AsymComm->Export_Comm_Proc_IDs() );
  
  // Zero temp tag for all phatnom faces and then loop over the MPCs
  // and mark each face with the number of MPCs it has.
  int my_proc = contact_processor_number( communicator );
  ContactFace<Real>** phantoms = Phantom_Faces();
  for( i=0 ; i<Number_Imported_Phantom_Faces() ; ++i)
    phantoms[i]->temp_tag = 0;


  for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
       cnfi != NULL;
       cnfi = node_entity_list.Node_Face_Iterator_Next()) {
    if( cnfi->Face()->Owner() != my_proc ){
      cnfi->Face()->temp_tag += 1;
    }
  }

  //
  // Get comm buffers to send the number of MPCs to each processor
  size_send_buf = face_asymcomm1.Num_Export_Comm_Partners()*sizeof(int);
  size_recv_buf = face_asymcomm1.Num_Import_Comm_Partners()*sizeof(int);
  num_request_handles = face_asymcomm1.Num_Import_Comm_Partners();
  CommBuffer()->Buffers( size_send_buf, size_recv_buf, num_request_handles,
			 &send_buf, &recv_buf, &request_handles );
  //
  // Load the buffers
  int* sbi = (int*) send_buf;
  for( i=0 ; i<face_asymcomm1.Num_Export_Comm_Partners() ; ++i){
    int proc = face_asymcomm1.Export_Comm_Proc_ID( i );
    sbi[i] = 0;
    for( j=0 ; j<Number_Imported_Phantom_Faces() ; ++j){
      if( phantoms[j]->Owner() == proc ){
	sbi[i] += phantoms[j]->temp_tag;
      }
    }
  }

  contact_communicate_packed_buffers( communicator, face_asymcomm1, sizeof(int),
				      send_buf,  recv_buf, request_handles );

  // Using the data just received (and that from the send buffer),
  // build a comm plan.
  num_import_partners = 0;
  num_export_partners = 0;
  num_face_mpcs = 0;
  int* ex_comm_proc_id = new int[face_asymcomm1.Num_Export_Comm_Partners()];
  int* im_comm_proc_id = new int[face_asymcomm1.Num_Import_Comm_Partners()];
  int* rbi = (int*) recv_buf;
  for( i=0 ; i<face_asymcomm1.Num_Import_Comm_Partners() ; ++i){
    if( rbi[i] ){
      im_comm_proc_id[num_import_partners] = 
	face_asymcomm1.Import_Comm_Proc_ID( i );
      import_length_per_proc[num_import_partners] = rbi[i];
      num_import_partners++;
      num_face_mpcs += rbi[i];
    }
  }
  for( i=0 ; i<face_asymcomm1.Num_Export_Comm_Partners() ; ++i){
    if( sbi[i] ){
      ex_comm_proc_id[num_export_partners] = 
	face_asymcomm1.Export_Comm_Proc_ID( i );
      export_length_per_proc[num_export_partners] = sbi[i];
      num_export_partners++;
    }
  }

  ContactAsymComm face_asymcomm2( num_export_partners,
				  num_import_partners,
				  export_length_per_proc,
				  import_length_per_proc,
				  ex_comm_proc_id,
				  im_comm_proc_id );
  
  // Allocate memory to hold the face mpcs
  if( num_face_mpcs )
    face_mpcs = new Contact_MPC_Equation*[num_face_mpcs];
  else
    face_mpcs = NULL;

  // Get comm buffers to send the MPCs to each processor
  //   Note: We can NOT use the Size_Export or Size_Import functions
  //         because we didn't construct the entity lists.  Instead
  //         loop to build this values.
  size_send_buf = 0;
  for( i=0 ; i<face_asymcomm2.Num_Export_Comm_Partners() ; ++i)
    size_send_buf += size_mpc*face_asymcomm2.Num_Export_to_Proc(i);
  size_recv_buf = 0;
  for( i=0 ; i<face_asymcomm2.Num_Import_Comm_Partners() ; ++i)
    size_recv_buf += size_mpc*face_asymcomm2.Num_Import_from_Proc(i);
  num_request_handles = face_asymcomm2.Num_Import_Comm_Partners();
  CommBuffer()->Buffers( size_send_buf, size_recv_buf, num_request_handles,
			 &send_buf, &recv_buf, &request_handles );

  // Load the data into the buffers
  buf_offset = 0;
  for( i=0 ; i<face_asymcomm2.Num_Export_Comm_Partners() ; ++i){
    int export_proc = face_asymcomm2.Export_Comm_Proc_ID( i );




    for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
         cnfi != NULL;
         cnfi = node_entity_list.Node_Face_Iterator_Next()) {
      if( export_proc == cnfi->Face()->Owner() ){
	// We now know the node which has the MPC, we must now find the
	// correct MPC for this face (since the node may be constrained
	// by more than one face).
	for( k=0 ; k<3 ; ++k){
	  PRECONDITION( slave_mpcs[3*cnfi->Node()->ProcArrayIndex()+k] );
            if( cnfi->Face()->Global_ID().LoInt() ==
              slave_mpcs[3*cnfi->Node()->ProcArrayIndex()+k]->
              master_face_local_id ) break;
	}
      
	POSTCONDITION( k<3 );
	PRECONDITION( buf_offset <= size_send_buf-size_mpc );
	Contact_MPC_Equation* mpc = 
	  slave_mpcs[3*cnfi->Node()->ProcArrayIndex()+k];
	std::memcpy( send_buf+buf_offset, mpc, size_mpc );
	buf_offset += size_mpc;
      }
    }
  }
  POSTCONDITION( buf_offset == size_send_buf );
  
  contact_communicate_packed_buffers( communicator, face_asymcomm2, size_mpc, 
				      send_buf, recv_buf, request_handles );

  // Create the face_mpcs
  rb = recv_buf;
  for( i=0 ; i<num_face_mpcs ; ++i){
    Contact_MPC_Equation* mpc = new Contact_MPC_Equation;
    face_mpcs[i] = mpc;
    std::memcpy( mpc, rb, size_mpc );
    rb += size_mpc;
  }

  delete [] import_length_per_proc;
  delete [] export_length_per_proc;
  delete [] im_comm_proc_id;
  delete [] ex_comm_proc_id;

#endif
}


void ContactMPCs::Display_MPCs()
{
  int i,j,k;

  ContactParOStream& postream = ParOStream();
  postream << "Entering Display_MPCs()\n";
  for( i=0 ; i<number_of_nodes ; ++i){
    for( k=0 ; k<3 ; ++k){
      Contact_MPC_Equation* mpc = slave_mpcs[3*i+k];
      if( mpc){
	postream << "MPC Equation (Slave Node) for node " 
		 << enforcement_node_list[i]->Exodus_ID()
		 << " (" << mpc->slave_node_proc_id << ","
		 << mpc->slave_node_local_id << ") "
		 << " and face (" << mpc->master_face_proc_id << ","
		 << mpc->master_face_local_id << ")\n";
	for( j=0 ; j<mpc->num_face_nodes ; ++j){
	  postream << "  Coefficient for node (" << mpc->face_node_proc_id[j]
		   << "," << mpc->face_node_local_id[j] << ") = " 
		   << mpc->face_node_coeff[j] << "\n";
	}
      }
    }
  }
  for( i=0 ; i<num_face_mpcs ; ++i){
    Contact_MPC_Equation* mpc = face_mpcs[i];
    postream << "MPC Equation (Master Face) for node "
	     << " (" << mpc->slave_node_proc_id << ","
	     << mpc->slave_node_local_id << ") "
	     << " and face (" << mpc->master_face_proc_id << ","
	     << mpc->master_face_local_id << ")\n";
    for( j=0 ; j<mpc->num_face_nodes ; ++j){
      postream << "  Coefficient for node (" << mpc->face_node_proc_id[j]
	       << "," << mpc->face_node_local_id[j] << ") = " 
	       << mpc->face_node_coeff[j] << "\n";
    }
  }
  postream.flush();
}


void ContactMPCs::Plot_MPC_Vars( )
{
  // Dump slave node mpc coefficients to exodus file
  // first node variable will have all mpc_coeff[0]'s
  // written at the slave node to which they apply.
  // Similarly, the second node var will have all 
  // mpc_coeff[1]'s etc...
  Real* plot_mpc = new Real[8*number_of_nodes];
  for( int i=0 ; i<3 ; ++i){
    std::memset( plot_mpc, 0, 8*number_of_nodes*sizeof(Real));
    for( int j=0 ; j<number_of_nodes ; ++j){
      if( slave_mpcs[3*j+i] ){
	Contact_MPC_Equation* mpc = slave_mpcs[3*j+i];
	for( int k=0; k<mpc->num_face_nodes; ++k) {
	  plot_mpc[j*8+k] = mpc->face_node_coeff[k];
	}
      }
    }
    Set_NVARS( 8, plot_mpc );
  }
  delete [] plot_mpc;
}


void ContactMPCs::Set_NVARS( int num_node_vars, Real* nodal_data )
{
  PRECONDITION( num_node_vars );
  int i, j;

  int old_data_length = 0;
  // if already exists resize and add most recent data to end
  if( nodal_plot_vars ) {
    old_data_length = number_nodal_plot_vars*number_of_nodes;
    Real *temp_vars = NULL;
    temp_vars = new Real[old_data_length];
    std::memcpy(temp_vars, nodal_plot_vars, old_data_length*sizeof(Real));
    delete [] nodal_plot_vars;
    nodal_plot_vars = NULL;
    number_nodal_plot_vars += num_node_vars;
    int new_data_length = number_nodal_plot_vars*number_of_nodes;
    nodal_plot_vars = new Real[new_data_length];
    std::memcpy(nodal_plot_vars, temp_vars, old_data_length*sizeof(Real));
    delete [] temp_vars;
  } else {
    old_data_length = 0;
    number_nodal_plot_vars = num_node_vars;
    nodal_plot_vars = new Real[num_node_vars*number_of_nodes];
  }

  // add new data  
  int var_count = old_data_length;
  for( i=0; i<num_node_vars; ++i) {
    for( j=0; j<number_of_nodes; ++j) {
      nodal_plot_vars[var_count++] = 
	nodal_data[j*num_node_vars+i];
    }
  }

  return;
}

void ContactMPCs::Get_Nodal_Plot_Variable( int var_num, Real* data )
{
  // variable numbering starts at 0
#ifndef CONTACT_NO_EXODUS_OUTPUT
  PRECONDITION( var_num < number_nodal_plot_vars );
  std::memcpy( data, nodal_plot_vars+var_num*number_of_nodes, 
	  number_of_nodes*sizeof(Real) );
#endif
}
