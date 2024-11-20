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

//#define SHELL_DEBUG

#include "ContactUtilities.h"
#include "ContactShellHandler.h"
#include "ContactNode.h"
#include "ContactSearch.h"
#include "ContactShellNode.h"
#include "Contact_Communication.h"
#include "ContactErrors.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include "ContactTDEnforcement.h"
#include "ContactParOStream.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactHostGlobalID.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#include "ContactAsymComm.h"
#include "ContactSymComm.h"
#include "ContactCommBuffer.h"
#include "Contact_Communication.h"
#endif

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

using acme::Dot;
using acme::Cross;
using acme::Magnitude;
using acme::Normalize;
using acme::Invert_3x3Matrix;

namespace {
  std::vector<int> faces_to_process;   // shell_face
}

ContactShellHandler::ContactShellHandler( ContactErrors* Errors,
					  int Dimensionality, 
					  int number_node_blocks,
      const ContactSearch::ContactNode_Type* node_block_types, 
					  const int* host_number_nodes_per_block,
					  int*& number_nodes_per_block,
					  const int* host_exodus_node_ids,
					  int*& exodus_node_ids,
					  const int* host_global_node_ids,
					  int*& global_node_ids,
	     ContactType** node_entity_types,
					  const double* coords,
					  int number_face_blocks, 
      const ContactSearch::ContactFace_Type* face_block_types,
					  const int* number_faces_per_block, 
					  const int* global_face_ids,
					  const int* host_face_connectivity,
					  int*& face_connectivity,
					  const double* face_lofting_factors,
					  int number_comm_partners, 
					  const int* comm_proc_ids,
					  const int* host_number_nodes_to_partner, 
					  int*& number_nodes_to_partner, 
					  const int* host_comm_nodes,
					  int*& comm_nodes,
					  MPI_Comm& communicator,
					  ContactSearch* search,
					  ContactTopology* Topology,
					  ContactCommBuffer* CommBuffer,
					  ShellID_Resolution res_method,
					  ContactShellHandler * old_handler,
                           ContactSearch::ContactErrorCode& error_code )
  : postream(communicator), dimensionality(Dimensionality), lofting_computed(false), 
    topology(Topology), comm_buffer(CommBuffer), SearchComm(communicator)
{

  //-----------------------------------------
  // Here is a list of the fundamental variables that we will use
  // for the mapping from user shell faces to the lofted surface.
  //
  //   number_host_code_nodes (Int) -- number of nodes supplied from 
  //        host code, both shell and not
  //   inv_conn_count[number_host_code_nodes] (Int) -- stores for
  //        ALL nodes the number of shell faces that are attached
  //        to it
  //   conn_offsets (Int) [number_face_blocks] -- stores the offsets into
  //        the face connectivity array for ALL face blocks 
  //   num_shell_faces (Int) -- number of user-supplied shell faces. Also
  //        includes faces that are not on shells but which are attached
  //        to a shell node.
  //   num_shell_nodes (Int) -- number of user-supplied nodes which are
  //        connected to original shell faces (<= number_host_code_nodes)
  //   tot_inv_conn_count (Int) -- total number of entries in the inv_conn
  //        array
  //   offset_inv_conn[number_host_code_nodes] (Int) -- stores the offsets
  //        into the inv_conn array for each node
  //   offset_counters[number_host_code_nodes] (Int) -- temporary array
  //        used for filling inv_conn arrays; initialized to same values
  //        as offset_inv_conn
  //   inv_conn[tot_inv_conn_count] (Int) -- this is the inverse connectivity
  //        going from nodes to the shell faces they are connected to. The
  //        number stored is the shell face that the node is connected to.
  //        To get the face id and face block passed in, must dereference
  //        using shell_face_id and shell_face_blk.
  //   shell_face_id[num_shell_faces] -- array which maps from the 
  //        shell face to the local face id of the face relative to
  //        its user-supplied face block. This list also includes
  //        faces that are not on shells but which are attached
  //        to a shell node.
  //   shell_face_blk[num_shell_faces] -- array which maps from the 
  //        shell face to the face block id of the user-supplied face.
  //        This list also includes faces that are not on shells but 
  //        which are attached to a shell node.
  //   opposite_face[num_shell_faces] -- array which defines what faces 
  //        are opposite sides of the shell. The numbering is in
  //        relation to only the faces related to shells. If the face
  //        is in the shell_face lists but is not a shell face (i.e. it
  //        is attached to a shell node), then the value in this
  //        array is -1.
  //   face_process_flag[4*num_shell_faces] (Int) -- flag array for 
  //        indicating if the face has been processed yet.
  //   face_node_process_flag[4*num_shell_faces] (Int) -- flag array for 
  //        indicating if the node on the face has been processed yet.
  //   num_created_shell_nodes[num_host_code_nodes] (Int) -- number of acme
  //        shell nodes created for each host code shell node
  //   host_to_acme_node_map[num_acme_nodes] (Int) -- map from 
  //        host code ids to acme ids. The first entry is the same as 
  //        the original number, and the second through N are the 
  //        acme numbers of the other nodes created for the host code
  //        shell node.  This array is jagged and requires going through
  //        the offset array to get the start for a given node.
  //   acme_to_host_node_map[num_shell_nodes] (Int) -- map from 
  //        acme ids to host code ids. 
  //   face_lofting_factors[num_faces] (double) -- has the lofting factors
  //        for all faces. Needed in topology creation.
  //-----------------------------------------

#ifndef CONTACT_NO_MPI
  my_proc_id = contact_processor_number(SearchComm);
#else
  my_proc_id = 0;
#endif
  
  int i,j,k,blk;
  int have_an_error = 0;
  
  //only written to handle shells in 3-D, not for 2-D analyses
  PRECONDITION(dimensionality == 3);

  // Do some initialization so if things fail, the code won't seg fault
  // in the destructor
  lofted_nodes = NULL;
  host_to_acme_node_map = NULL;
  acme_to_host_node_map = NULL;
  num_created_shell_nodes = NULL;
  num_acme_nodes_for_host_code_node = NULL;
  offset_for_node = NULL;
  block_offset_for_nodes = NULL;
  host_nodes_per_block = NULL;
  opposite_face = NULL;
  shell_face_id = NULL;
  shell_face_blk = NULL;
  shell_loft_factors = NULL;
  orig_node_block_types = NULL;
  orig_number_nodes_per_block = NULL;
  orig_exodus_node_ids = NULL;
  orig_global_node_ids = NULL;
  orig_face_connectivity = NULL;
  orig_comm_proc_ids = NULL;
  orig_number_nodes_to_partner = NULL;
  orig_comm_nodes = NULL;
  new_number_nodes_per_block = NULL;
  new_exodus_node_ids = NULL;
  new_global_node_ids = NULL;
  new_face_connectivity = NULL;
  new_number_nodes_to_partner = NULL;
  new_comm_nodes = NULL;
  max_num_acme_nodes_for_host_code_node = 0;
#ifndef CONTACT_NO_MPI
  comm_plan = NULL;
  comm_plan_to_owner = NULL;
  comm_plan_to_ghost = NULL;
  import_loft_face_data_comm_plan = NULL;
  number_lofting_ghosts = NULL;
  offset_to_buffer_offsets = NULL;
  buffer_offsets = NULL;
#endif
  is_a_tab_node = NULL;
  
#ifdef SHELL_DEBUG
  //-- use following 3 lines to convert postream to an std::ofstream
  //-- also look in header, and comment out initialization in 
  //-- constructor above.
  //   char ook[100];
  //   std::sprintf(ook,"dump_%d",my_proc_id);
  //   postream.open(ook);
  postream << my_proc_id << ":Hey there, I am in ContactShellHandler::constructor." << "\n";
#endif


  //---------------------------------------------------------------------
  // store away the nodal data that was passed in and allocate some 
  // basic data structures we need.  Also allocate arrays for the
  // host code data that we need to manipulate in this routine.

  int number_faces = 0;
  for ( i = 0; i < number_face_blocks; ++i) 
    number_faces += number_faces_per_block[i];
  m_number_faces = number_faces;

  number_host_code_nodes = 0;
  host_nodes_per_block = new int[number_node_blocks];
  for ( i = 0; i < number_node_blocks; ++i){ 
    number_host_code_nodes += host_number_nodes_per_block[i];
    host_nodes_per_block[i] = host_number_nodes_per_block[i];
  }

  // make sure coords and lofting factors were passed in if nodes and/or
  // faces were passed in on this processor
  PRECONDITION(number_host_code_nodes == 0 || coords);
  PRECONDITION(number_faces == 0 || face_lofting_factors);

  orig_exodus_node_ids = new int[number_host_code_nodes];
  orig_global_node_ids = new int[2*number_host_code_nodes];
  for ( i = 0; i < number_host_code_nodes; ++i){
    orig_exodus_node_ids[i] = host_exodus_node_ids[i];
    orig_global_node_ids[2*i] = host_global_node_ids[2*i];
    orig_global_node_ids[2*i + 1] = host_global_node_ids[2*i + 1];
  }

  Initialize_Node_Maps();

  orig_node_block_types = new int[number_node_blocks];
  orig_number_nodes_per_block = new int[number_node_blocks];
  new_number_nodes_per_block = new int[number_node_blocks];
  number_nodes_per_block = new_number_nodes_per_block;
  for (i = 0; i < number_node_blocks; ++i){
    orig_node_block_types[i] = node_block_types[i];
    orig_number_nodes_per_block[i] = host_number_nodes_per_block[i];
    new_number_nodes_per_block[i] = host_number_nodes_per_block[i];
  }

  // create connectivity offsets
  int * conn_offsets = new int[number_face_blocks];
  int num_nodes_in_orig_conn = 0;
  for ( blk = 0; blk < number_face_blocks; blk ++) {
    conn_offsets[blk]=num_nodes_in_orig_conn;
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk]);
    POSTCONDITION(num_face_nodes>0);
    num_nodes_in_orig_conn += num_face_nodes*number_faces_per_block[blk];
  }
  
  orig_face_connectivity = new int[num_nodes_in_orig_conn];
  new_face_connectivity = new int[num_nodes_in_orig_conn];
  face_connectivity = new_face_connectivity;
  for ( i = 0; i < num_nodes_in_orig_conn; ++i) {
    orig_face_connectivity[i] = host_face_connectivity[i];
    new_face_connectivity[i] = host_face_connectivity[i];
  }

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> number_host_code_nodes: " << number_host_code_nodes << "\n";
  for( i=0 ; i<number_host_code_nodes ; ++i){
    postream << "  Node " << i << " Coords: " << coords[3*i] << " " 
	     << coords[3*i+1] << " " << coords[3*i+2] << "\n";
  }
#endif
  
#ifdef SHELL_DEBUG
  {
    postream << my_proc_id << ":HEY!!!! Here is original connectivity." << "\n";
    int print_count = 0;
    for ( blk = 0; blk < number_face_blocks; blk ++) {
      postream << my_proc_id << ":   >> Block " << blk << "\n";
      int num_face_nodes = 
	search->Number_Nodes_Per_Face(face_block_types[blk]);
      POSTCONDITION(num_face_nodes>0);
      for ( int face = 0; face < number_faces_per_block[blk]; face ++) {
	char ook[80];
	std::sprintf(ook,"%d:      ->Face %d :",my_proc_id,face);
	//foo << my_proc_id << ":     -> Face " << face << " :";
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int node = face_connectivity[print_count]-1; 
	  char fishy[10];
	  std::sprintf(fishy," %d",node);
	  std::strcat(ook,fishy);
	  //foo << " " << node;
	  print_count++;
	}
	postream << ook << "\n";
      }
    }
    postream << my_proc_id << ":HEY!!!! Done with original connectivity." << "\n";
  }
  {
    postream << my_proc_id << ":HEY!!!! Here is original face ids:" << "\n";
    int print_count = 0;
    for ( blk = 0; blk < number_face_blocks; blk ++) {
      postream << my_proc_id << ":   >> Block " << blk << "\n";
      for ( int face = 0; face < number_faces_per_block[blk]; face ++) {
	postream << my_proc_id << ":      >> face " << face
	     << " id: " << global_face_ids[print_count++]
	     << " " << global_face_ids[print_count++] << "\n";
      }
    }
    postream << my_proc_id << ":HEY!!!! Done with original face ids." << "\n";
  }
  {
    postream << my_proc_id << ":HEY!!!! Here is original node ids." << "\n";
    for ( i = 0; i < number_host_code_nodes; ++i){
      postream << my_proc_id << ":TEST>      node: " << i << " HID: "
	       << host_global_node_ids[2*i] << "," << host_global_node_ids[2*i+1]
	       << "  EXID: " << host_exodus_node_ids[i] << "\n";
    }
    postream << my_proc_id << ":HEY!!!! Done with original node ids." << "\n";
    postream.flush();
  }
#endif

  // initialize the number of acme nodes created per host code shell node
  num_created_shell_nodes = new int[number_host_code_nodes];
  for( i=0 ; i<number_host_code_nodes ; ++i) num_created_shell_nodes[i] = 0;
    
  // Create the first node from each host node (we will only create 
  // additional lofted nodes from here.
  for( i=0 ; i<number_host_code_nodes ; ++i){
    Create_Node( i, i);
  }
  next_new_id = number_host_code_nodes;
  

  //---------------------------------------------------------------------
  // get inverse connectivities 
  int num_shell_faces(0);
  int num_shell_nodes(0);
  int* inv_conn = NULL;
  int *offset_inv_conn = NULL;
  int *inv_conn_count = NULL;
  compute_inv_conn(number_host_code_nodes,
		   number_face_blocks,
		   face_connectivity,
		   conn_offsets,
		   face_block_types,
		   number_faces_per_block,
		   face_lofting_factors,
		   number_comm_partners,
		   comm_proc_ids,
		   host_number_nodes_to_partner,
		   host_comm_nodes,
		   inv_conn,
		   offset_inv_conn,
		   inv_conn_count,
		   shell_face_blk,
		   shell_face_id,
		   shell_loft_factors,
		   num_shell_faces,
		   num_shell_nodes,
		   true,
		   search);
  POSTCONDITION(inv_conn != NULL);
  POSTCONDITION(offset_inv_conn != NULL);
  POSTCONDITION(inv_conn_count != NULL);
  POSTCONDITION(shell_face_id != NULL);
  POSTCONDITION(shell_face_blk != NULL);

#ifdef SHELL_DEBUG
  {
    postream << my_proc_id << ":HEY!!!! Here is local inv_conn." << "\n";
    for ( i = 0; i < number_host_code_nodes ; i ++) {
      char ook[80];
      std::sprintf(ook,"%d:   Node %d :",my_proc_id,i);
      int off = offset_inv_conn[i];
      for ( j = 0; j < inv_conn_count[i]; ++j){
	char fishy[10];
	std::sprintf(fishy," %d",inv_conn[off+j]);
	std::strcat(ook, fishy);
      }
      postream << ook << "\n";
    }
    postream << my_proc_id << ":HEY!!!! Done with inv_conn." << "\n";
    postream.flush();
  }
#endif
  
  //---------------------------------------------------------------------
  // for parallel, we want to ghost all of the faces attached to shared
  // nodes to the sharing processors. That way all processors can
  // compute their topology independently.
  //
  // The data structures that have both local and ghosted face connectivites
  // start with "full_", for example full_face_connectivities. For serial
  // cases, we just make the "full_" pointers point to the original
  // data structures.

  int * full_face_connectivity = NULL;
  int * full_face_conn_update = NULL;
  int * full_conn_offsets = NULL;
  ContactSearch::ContactFace_Type * full_face_block_types = NULL;
  int * full_number_faces_per_block = NULL;
  int full_number_face_blocks = 0;
  int * full_inv_conn = NULL;
  int * full_inv_conn_count = NULL;
  int * full_offset_inv_conn = NULL;
  int * full_shell_face_id = NULL;
  int * full_shell_face_blk = NULL;
  int * full_global_face_ids = NULL;
  double * full_coords = NULL;
  int full_num_shell_faces = 0;
  double * full_shell_loft_factors = NULL;
#ifndef CONTACT_NO_MPI
  int full_num_shell_nodes = 0;
  int full_num_host_code_nodes = 0;
#endif

  if ( contact_number_of_processors(SearchComm) == 1 ) {
    // serial case -- just make "full_" data structures point to 
    // original data structures
    full_face_connectivity = orig_face_connectivity;
    full_face_conn_update = face_connectivity;
    full_conn_offsets = conn_offsets;
    full_face_block_types = 
      new ContactSearch::ContactFace_Type[number_face_blocks];
    full_number_faces_per_block = new int[number_face_blocks];
    for ( i = 0; i < number_face_blocks; ++i){
      full_face_block_types[i] = face_block_types[i];
      full_number_faces_per_block[i] = number_faces_per_block[i];
    }
    full_number_face_blocks = number_face_blocks;
    full_inv_conn = inv_conn;
    full_inv_conn_count = inv_conn_count;
    full_offset_inv_conn = offset_inv_conn;
    full_shell_face_id = shell_face_id;
    full_shell_face_blk = shell_face_blk;
    full_global_face_ids = new int[number_faces*2];
    for ( i = 0; i < number_faces*2; ++i)
      full_global_face_ids[i] = global_face_ids[i];
    full_num_shell_faces = num_shell_faces;
    full_coords = new double[number_host_code_nodes*3];
    for ( i = 0; i < number_host_code_nodes*3; ++i)
      full_coords[i] = coords[i]; 
    full_shell_loft_factors = shell_loft_factors;
#ifndef CONTACT_NO_MPI
    full_num_shell_nodes = num_shell_nodes;
    full_num_host_code_nodes = number_host_code_nodes;
#endif
  }
#ifndef CONTACT_NO_MPI
  else {
    // parallel case.
    
    // compute and communicate the face connectivities for faces attached
    // to shell nodes.
    error_code = ghost_shell_faces(number_comm_partners,
				   comm_proc_ids,
				   host_number_nodes_to_partner,
				   host_comm_nodes,
				   num_shell_faces,
				   number_face_blocks,
				   inv_conn,
				   inv_conn_count,
				   offset_inv_conn,
				   face_connectivity,
				   conn_offsets,
				   face_block_types,
				   number_faces_per_block,
				   global_face_ids,
				   coords,
				   shell_face_blk,
				   shell_face_id,
				   shell_loft_factors,
				   full_number_face_blocks,
				   full_face_connectivity,
				   full_conn_offsets,
				   full_face_block_types,
				   full_number_faces_per_block,
				   full_global_face_ids,
				   full_num_host_code_nodes,
				   full_coords,
				   full_shell_loft_factors,
				   search );
    error_code = (ContactSearch::ContactErrorCode) 
      contact_global_error_check( error_code, SearchComm );
    if( error_code != ContactSearch::NO_ERROR ) return;

    POSTCONDITION(full_face_connectivity != NULL);
    POSTCONDITION(full_conn_offsets != NULL);
    POSTCONDITION(full_face_block_types != NULL);
    POSTCONDITION(full_number_faces_per_block != NULL);
  
    // make copy of host code connectivity into scratch array. This version 
    // we will permit to be changed.
    int last_blk = full_number_face_blocks-1;
    int full_conn_size = full_conn_offsets[last_blk]
      + search->Number_Nodes_Per_Face(full_face_block_types[last_blk])
      * full_number_faces_per_block[last_blk];
    full_face_conn_update = new int[full_conn_size];
    std::memset (full_face_conn_update, 0, full_conn_size*sizeof(int));
    // copy full connectivity data into copy, but don't set values
    // for nodes that are not on this processor
    int idx = 0;
    for ( i = 0; i < full_number_face_blocks; ++i){
      int num_face_nodes = 
	search->Number_Nodes_Per_Face(full_face_block_types[i]);
      POSTCONDITION(num_face_nodes>0);
      for ( j = 0; j < full_number_faces_per_block[i]*num_face_nodes; ++j){
	int node = full_face_connectivity[idx];
	if (node <= number_host_code_nodes) // fortran #ering
	  full_face_conn_update[idx] = node;
	idx++;
      }
    }
    
#ifdef SHELL_DEBUG
    {
      postream << my_proc_id << ":HEY!!!! Here is full connectivity." << "\n";
      int print_count = 0;
      for ( blk = 0; blk < full_number_face_blocks; blk ++) {
	postream << my_proc_id << ":   >> Block " << blk << "\n";
	int num_face_nodes = 
	  search->Number_Nodes_Per_Face(full_face_block_types[blk]);
	POSTCONDITION(num_face_nodes>0);
	for ( int face = 0; face < full_number_faces_per_block[blk]; face ++) {
	  char ook[80];
	  std::sprintf(ook,"%d:      ->Face %d :",my_proc_id,face);
	  //foo << my_proc_id << ":     -> Face " << face << " :";
	  for ( int inode = 0; inode < num_face_nodes; inode ++) {
	    int node = full_face_connectivity[print_count]-1; 
	    char fishy[10];
	    std::sprintf(fishy," %d",node);
	    std::strcat(ook,fishy);
	    //foo << " " << node;
	    print_count++;
	  }
	  postream << ook << "\n";
	}
      }
      postream << my_proc_id << ":HEY!!!! Done with full connectivity." << "\n";
      
    }
#endif
    
    // use new connectivity to compute a new inverse connectivity
    double * dummy = NULL;
    compute_inv_conn(full_num_host_code_nodes,
		     full_number_face_blocks,
		     full_face_connectivity,
		     full_conn_offsets,
		     full_face_block_types,
		     full_number_faces_per_block,
		     NULL, 
		     0,
		     NULL,
		     NULL,
		     NULL,
		     full_inv_conn,
		     full_offset_inv_conn,
		     full_inv_conn_count,
		     full_shell_face_blk,
		     full_shell_face_id,
		     dummy,
		     full_num_shell_faces,
		     full_num_shell_nodes,
		     false,
		     search);

    POSTCONDITION( full_face_block_types != NULL);
    POSTCONDITION ( full_num_shell_faces >= num_shell_faces);
#ifdef SHELL_DEBUG
    {
      postream << my_proc_id << ":HEY!!!! Here is full inv_conn." << "\n";
      for ( i = 0; i < full_num_host_code_nodes ; i ++) {
	char ook[80];
	std::sprintf(ook,"%d:   Node %d :",my_proc_id,i);
	int off = full_offset_inv_conn[i];
	for ( j = 0; j < full_inv_conn_count[i]; ++j){
	  char fishy[10];
	  std::sprintf(fishy," %d",full_inv_conn[off+j]);
	  std::strcat(ook, fishy);
	}
	postream << ook << "\n";
      }
      postream << my_proc_id << ":HEY!!!! Done with full inv_conn." << "\n";
      postream.flush();
    }
    {
      postream << my_proc_id << ":HEY!!!! Here is full face ids -- only "
	       << "for shell (or shell-attached) faces:" << "\n";
      for ( i = 0; i < full_num_shell_faces; i ++) {
	int blk = full_shell_face_blk[i];
	int id = full_shell_face_id[i];
	postream << my_proc_id << ":   >> face blk " << blk 
		 << " id " << id << " gid " <<  full_global_face_ids[2*i]
		 << "|" <<  full_global_face_ids[2*i+1] 
		 << " loft fact:" << full_shell_loft_factors[i] << "\n";
      }
      postream << my_proc_id << ":HEY!!!! Done with full face ids." << "\n";
      postream.flush();
    }
#endif
    
  } 
#endif // CONTACT_NO_MPI

  // compute full number of faces & store it on class. This is for 
  // asserts in the code.
  m_full_number_faces = 0;
  for ( i = 0; i < full_number_face_blocks; ++i)
    m_full_number_faces += full_number_faces_per_block[i];

  //---------------------------------------------------------------------
  // Now we have completed the inverse connectivities. We can now use these to 
  // construct the shell topology. 
  //

  // allocate data needed
  //   NOTE: Here we assume largest number of edges per face = 4
  //   NOTE: Here we assume largest number of nodes per face = 4
  opposite_face = new int[full_num_shell_faces];
  bool * face_process_flag = new bool[full_num_shell_faces];
  bool * face_node_process_flag = new bool[full_num_shell_faces*4];
  for (i = 0; i < full_num_shell_faces; ++i) {
    opposite_face[i]=-1;
    face_process_flag[i] = false;
    for (j =0; j < 4; ++j)
      face_node_process_flag[i*4 + j] = false;
  }

  // set the node-on-face process flag to true if the node is
  // part of a ghosted face and is not shared by this processor.
  // This will permit us to properly process the connectivities
  // of ghosted faces (important for n-way sections) while 
  // avoiding the creation of acme nodes for host code nodes
  // that don't exist on this processor.
  // Also mark as processed all non-shell nodes, which arise from
  // non-shell faces that are attached to shell nodes. 
  for (i = 0; i < full_num_shell_faces; ++i) {
    // get block and id of shell face
    int block = full_shell_face_blk[i];
    int id = full_shell_face_id[i];
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(full_face_block_types[block]);
    POSTCONDITION(num_face_nodes>0);
    bool face_is_shell = ContactSearch::Is_a_Shell_Face(full_face_block_types[block]);
    int conn_base_idx = full_conn_offsets[block] + id*num_face_nodes;
    for (j =0; j < num_face_nodes; ++j) {
      int node_id = full_face_connectivity[conn_base_idx + j] -1;
      if ( face_is_shell ) {
	if (node_id < number_host_code_nodes) continue;
	//   NOTE: Here we assume largest number of nodes per face = 4
	face_node_process_flag[i*4 + j] = true;
      } else {
	//   NOTE: Here we assume largest number of nodes per face = 4
	if (node_id >= number_host_code_nodes) 
	  face_node_process_flag[i*4 + j] = true;
	else if ( inv_conn_count[node_id] == 0 )
	  face_node_process_flag[i*4 + j] = true;
      }
    }
  }

  // loop over the shell faces -- process them
  for ( int shell_face = 0; shell_face < full_num_shell_faces; shell_face++){

#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>>>Process face " << shell_face << "\n";
#endif

    // we now have a face to process. Call face process routine.
    faces_to_process.clear();
    faces_to_process.push_back(shell_face);
    for (int faces_it = 0; faces_it < faces_to_process.size(); ++faces_it) {

      if (face_process_flag[faces_to_process[faces_it]]) {
        continue;
      }

      bool face_processed = process_face(faces_to_process[faces_it],
  				       face_process_flag,
  				       face_node_process_flag,
  				       full_face_block_types,
  				       full_face_conn_update,  // changable
  				       full_face_connectivity, // reference
  				       full_conn_offsets,
  				       full_inv_conn,
  				       full_offset_inv_conn,
  				       full_inv_conn_count,
  				       full_shell_face_blk,
  				       full_shell_face_id,
  				       full_shell_loft_factors,
  				       number_nodes_per_block,
  				       full_coords,
  				       search ); 
      if (! face_processed) {
        topology->errors->Add_Error_Message("Internal Error Processing Faces (1)");
        have_an_error = 1;
        break;
      }
    }
  } // end loop on faces

#ifdef SHELL_DEBUG
    postream.flush();
#endif

  if( contact_global_error_check( have_an_error, SearchComm ) >= 1 ){
    error_code = ContactSearch::INTERNAL_ERROR;
    return;
  }
  
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> Done with all faces." << "\n";
  postream.flush();
#endif

  //---------------------------------------------------------------------
  // we have completed the construction of the shell topology. now make data
  // structures consistent and create correct ids for the new nodes.

  // renumber std::internal node numbering to be consistent with multiple node 
  // blocks
#ifdef SHELL_DEBUG
  {
  postream<<"face connectivity prior to renumber\n";
  int cnt = 0;
  int idx = 0;
  for ( i = 0; i < number_face_blocks; ++i){
    int nnodes = search->Number_Nodes_Per_Face(face_block_types[i]);
    for ( j = 0; j < number_faces_per_block[i]; ++j){
      postream<<"  face "<<cnt<<" ("<<global_face_ids[2*cnt+0]
                              <<", "<<global_face_ids[2*cnt+1]<<"):";
      for ( k = 0; k < nnodes; ++k){
        postream<<"  "<<full_face_conn_update[idx]; idx++;
      }
      postream<<"\n"; cnt++;
    }
  }
  postream.flush();
  }
#endif
  Renumber_Nodes( number_node_blocks,
		  host_number_nodes_per_block,
		  number_nodes_per_block,
		  full_number_face_blocks,
		  full_face_block_types,
		  full_number_faces_per_block,
		  full_face_conn_update );
#ifdef SHELL_DEBUG
  {
  postream<<"face connectivity after renumber\n";
  int cnt = 0;
  int idx = 0;
  for ( i = 0; i < number_face_blocks; ++i){
    int nnodes = search->Number_Nodes_Per_Face(face_block_types[i]);
    for ( j = 0; j < number_faces_per_block[i]; ++j){
      postream<<"  face "<<cnt<<" ("<<global_face_ids[2*cnt+0]
                              <<", "<<global_face_ids[2*cnt+1]<<"):";
      for ( k = 0; k < nnodes; ++k){
        postream<<"  "<<full_face_conn_update[idx]; idx++;
      }
      postream<<"\n"; cnt++;
    }
  }
  postream.flush();
  }
#endif

#ifndef CONTACT_NO_MPI
  // if in parallel, we changed a copy of the full connectivity as we 
  // processed the faces. Now copy this changed connectivity data into
  // the array passed in by the user
  int idx = 0;
  for ( i = 0; i < number_face_blocks; ++i){
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[i]);
    POSTCONDITION(num_face_nodes>0);
    for ( j = 0; j < number_faces_per_block[i]*num_face_nodes; ++j){
      face_connectivity[idx] = full_face_conn_update[idx];
      idx++;
    }
  }
#endif

#ifdef SHELL_DEBUG
  {
    postream << my_proc_id << ":HEY!!!! Here is full updated conn." << "\n";
    int print_count = 0;
    for ( blk = 0; blk < full_number_face_blocks; blk ++) {
      postream << my_proc_id << ":   >> Block " << blk << "\n";
      int num_face_nodes = 
	search->Number_Nodes_Per_Face(full_face_block_types[blk]);
      POSTCONDITION(num_face_nodes>0);
      for ( int face = 0; face < full_number_faces_per_block[blk]; face ++) {
	char ook[80];
	std::sprintf(ook,"%d:      ->Face %d :",my_proc_id,face);
	//foo << my_proc_id << ":     -> Face " << face << " :";
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int node = full_face_conn_update[print_count]-1; 
	  char fishy[10];
	  std::sprintf(fishy," %d",node);
	  std::strcat(ook,fishy);
	  //foo << " " << node;
	  print_count++;
	}
	postream << ook << "\n";
      }
    }
    postream << my_proc_id << ":HEY!!!! Done with full updated conn." << "\n";
    postream.flush();
  }
#endif

  Complete_Node_Maps();

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> *** The global ids as passed in.." << "\n";
  for ( i = 0; i < number_host_code_nodes; ++i){
    postream << my_proc_id << ":TEST>      node: " << i << " HID: "
	 << host_global_node_ids[2*i] << "," << host_global_node_ids[2*i+1]
	 << "  EXID: " << host_exodus_node_ids[i] << "\n";
  }
  postream.flush();
#endif

  // allocate arrays for new node ids, and copy old ids to new ids for one
  // of the nodes at each shell node.
  new_exodus_node_ids = new int[number_acme_nodes];
  exodus_node_ids = new_exodus_node_ids;
  new_global_node_ids = new int[2*number_acme_nodes];
  global_node_ids = new_global_node_ids;
  for ( i = 0; i < number_host_code_nodes; ++i){
    new_exodus_node_ids[i] = host_exodus_node_ids[i];
    new_global_node_ids[2*i] = host_global_node_ids[2*i];
    new_global_node_ids[2*i + 1] = host_global_node_ids[2*i + 1];
  }
		   
  // allocate new nodal comm maps to include newly created nodes from
  // shell nodes.
  orig_comm_proc_ids = new int[number_comm_partners];
  orig_number_nodes_to_partner = new int[number_comm_partners];
  new_number_nodes_to_partner = new int[number_comm_partners];
  number_nodes_to_partner = new_number_nodes_to_partner;
  orig_number_comm_nodes = 0;
  int new_number_comm_nodes = 0;
  int comm_offset = 0;
  for ( i = 0; i < number_comm_partners; ++i){
    orig_comm_proc_ids[i] = comm_proc_ids[i];
    orig_number_nodes_to_partner[i] = host_number_nodes_to_partner[i];
    new_number_nodes_to_partner[i] = host_number_nodes_to_partner[i];
    orig_number_comm_nodes += number_nodes_to_partner[i];
    for( j=0 ; j<host_number_nodes_to_partner[i] ; ++j){
      int host_comm_node = host_comm_nodes[comm_offset++];
      int num_acme_nodes = num_acme_nodes_for_host_code_node[host_comm_node-1];
      new_number_comm_nodes += num_acme_nodes;
    }
  }
  orig_comm_nodes = new int[orig_number_comm_nodes];
  new_comm_nodes = new int[new_number_comm_nodes];
  comm_nodes = new_comm_nodes;
  for ( i = 0; i < orig_number_comm_nodes; ++i) {
    orig_comm_nodes[i] = host_comm_nodes[i];
    new_comm_nodes[i] = host_comm_nodes[i];
  }

  // compute new inverse connectivity based on the new updated
  // full topology, i.e. all owned and ghosted faces. This is needed
  // for both id resolution and for looking for nodes at tab corners.  
  int * full_updated_inv_conn(NULL);
  int * full_updated_offset_inv_conn(NULL);
  int * full_updated_inv_conn_count(NULL);
  int * full_updated_shell_face_blk(NULL);
  int * full_updated_shell_face_id(NULL);
  int full_updated_num_shell_faces;
  int full_updated_num_shell_nodes;
  double * dummy = NULL;
  compute_inv_conn(number_acme_nodes,
		   full_number_face_blocks,
		   full_face_conn_update,
		   full_conn_offsets,
		   full_face_block_types,
		   full_number_faces_per_block,
		   NULL,
		   0,
		   NULL,
		   NULL,
		   NULL,
		   full_updated_inv_conn,
		   full_updated_offset_inv_conn,
		   full_updated_inv_conn_count,
		   full_updated_shell_face_blk,
		   full_updated_shell_face_id,
		   dummy,
		   full_updated_num_shell_faces,
		   full_updated_num_shell_nodes,
		   false,
		   search);
  
  POSTCONDITION(full_updated_inv_conn);
  POSTCONDITION(full_updated_offset_inv_conn);
  POSTCONDITION(full_updated_inv_conn_count );
  POSTCONDITION(full_updated_shell_face_blk );
  POSTCONDITION(full_updated_shell_face_id  );
  
#ifdef SHELL_DEBUG
  {
    postream << my_proc_id << ":HEY!!!! Here is full_updated_inv_conn" << "\n";
    for ( i = 0; i < number_acme_nodes ; i ++) {
      char ook[80];
      std::sprintf(ook,"%d:   Node %d :",my_proc_id,i);
      int off = full_updated_offset_inv_conn[i];
      for ( j = 0; j < full_updated_inv_conn_count[i]; ++j){
	char fishy[10];
	std::sprintf(fishy," %d",full_updated_inv_conn[off+j]);
	std::strcat(ook, fishy);
      }
      postream << ook << "\n";
    }
    postream << my_proc_id << ":HEY!!!! Done with inv_conn in par_res" << "\n";
  }
#endif

  // find any tab nodes -- i.e. acme nodes that are at the corner of a tab
  // construction
#ifdef SHELL_DEBUG
  postream << my_proc_id << ": Find tab nodes\n";
#endif
  is_a_tab_node = new bool[number_acme_nodes];
  for ( i = 0; i < number_acme_nodes; ++i) {
    bool tab_node = false;
    for ( j = 0; j < full_updated_inv_conn_count[i]-1; ++j) {
      for ( k = j+1; k < full_updated_inv_conn_count[i]; ++k) {
	int face_j = 
	  full_updated_inv_conn[full_updated_offset_inv_conn[i]+j];
	int face_k = 
	  full_updated_inv_conn[full_updated_offset_inv_conn[i]+k];
	// see if faces are opposites -- if so, this is a tab node
#ifdef SHELL_DEBUG
	postream << my_proc_id << ":         node " << i 
		 << " j,k: " << j << " " << k
		 << " faces " << face_j << " " << face_k 
		 << " opposite of face_j: " << opposite_face[face_j]
		 << "\n";
#endif
	if ( opposite_face[face_j] == face_k ) {
	  tab_node = true;
	  break;
	}
      }
      if ( tab_node ) break;
    }
    is_a_tab_node[i] = tab_node;
#ifdef SHELL_DEBUG
    if (tab_node) postream << my_proc_id << ":    ****** node " << i 
			   << " is a tab node\n";
#endif
  }
#ifdef SHELL_DEBUG
  postream << my_proc_id << ": Done finding tab nodes\n";
  postream.flush();
#endif

  // if single processor, then finish node ids
  if( contact_number_of_processors(SearchComm) == 1 ){
    if ( res_method ==  SEQUENTIAL_IDS ) {
      int max_id = 0;
      for ( i = 0; i < number_host_code_nodes; ++i)
	max_id = std::max(max_id,exodus_node_ids[i]);
      for( i=number_host_code_nodes ; i<number_acme_nodes ; ++i){
        exodus_node_ids[i]     = max_id+i+1;
        global_node_ids[2*i]   = 0;
        global_node_ids[2*i+1] = max_id+i+1;
      }
    } else {
      number_nodes_using_handler( 0, 0, NULL,
  			  	  full_updated_inv_conn,
  				  full_updated_inv_conn_count,
				  full_updated_offset_inv_conn,
				  full_updated_shell_face_blk,
				  full_updated_shell_face_id,
				  full_number_faces_per_block,
				  full_global_face_ids,
				  old_handler,
				  exodus_node_ids, 
				  global_node_ids);

    } 
  }
#ifndef CONTACT_NO_MPI
  else {
  //---------------------------------------------------------------------
  // now do a final commuication step to finalize all of the shell 
  // topology construction that we have been doing in serial.
    error_code = 
      parallel_node_resolution(full_number_faces_per_block,
			       full_updated_inv_conn,
			       full_updated_offset_inv_conn,
			       full_updated_inv_conn_count,
			       full_updated_shell_face_blk,
			       full_updated_shell_face_id,
			       full_global_face_ids,
			       number_comm_partners,
			       orig_comm_proc_ids,
			       orig_number_nodes_to_partner,
			       orig_comm_nodes,
			       number_nodes_to_partner, // changeable version
			       comm_nodes,              // changeable version
			       exodus_node_ids,         // changeable version
			       global_node_ids,         // changeable version
			       res_method,
			       old_handler,
			       search);
    
    error_code = (ContactSearch::ContactErrorCode) 
      contact_global_error_check( error_code, SearchComm );
    if( error_code != ContactSearch::NO_ERROR ) return;
  }
#endif
  
  delete [] full_updated_inv_conn;
  delete [] full_updated_offset_inv_conn;
  delete [] full_updated_inv_conn_count;
  delete [] full_updated_shell_face_blk;
  delete [] full_updated_shell_face_id;
  
  // create node map between acme and host code
  number_lofted_nodes = 0;
  acme_to_host_node_map = new int[number_acme_nodes];
  *node_entity_types = new ContactType[number_acme_nodes];
  for( i=0 ; i<number_host_code_nodes ; ++i){
    if ( num_created_shell_nodes[i] > 1) 
      number_lofted_nodes += num_created_shell_nodes[i];
    for( j=0 ; j<num_acme_nodes_for_host_code_node[i] ; ++j)
      acme_to_host_node_map[Acme_Node_for_Host_Node( i, j )] = i;
    if( Num_Acme_Nodes_for_Host_Node( i ) == 1 ){
      (*node_entity_types)[Acme_Node_for_Host_Node( i, 0)] = CT_NODE;
    } else {
      for( j=0 ; j<Num_Acme_Nodes_for_Host_Node(i) ; ++j)
	(*node_entity_types)[Acme_Node_for_Host_Node( i, j )] = 
	  CT_SHELL_NODE;
    }
  }
  
  // allocate space to hold pointers to lofted nodes  
  lofted_nodes = new ContactShellNode*[number_lofted_nodes];
  
  // delete temporary data
  delete [] face_process_flag;
  delete [] face_node_process_flag;
  delete [] conn_offsets;
  delete [] inv_conn;
  delete [] offset_inv_conn;
  delete [] inv_conn_count;

  delete [] full_face_block_types;
  delete [] full_number_faces_per_block;
  delete [] full_global_face_ids;
  delete [] full_coords;
  if ( contact_number_of_processors(SearchComm) > 1 ) {
    delete [] full_face_connectivity;
    delete [] full_face_conn_update;
    delete [] full_conn_offsets;
    delete [] full_inv_conn;
    delete [] full_inv_conn_count;
    delete [] full_offset_inv_conn;
    delete [] full_shell_face_id;
    delete [] full_shell_face_blk;
    delete [] full_shell_loft_factors;
  }
  
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":HEY!!!! Here is returned connectivity." << "\n";
  int print_count = 0;
  for ( blk = 0; blk < number_face_blocks; blk ++) {
    postream << my_proc_id << ":   >> Block " << blk << "\n";
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk]);
    POSTCONDITION(num_face_nodes>0);
    for ( int face = 0; face < number_faces_per_block[blk]; face ++) {
      char foo[80];
      std::sprintf(foo,"%d:      ->Face %d :",my_proc_id,face);
      for ( int inode = 0; inode < num_face_nodes; inode ++) {
	int node = face_connectivity[print_count]-1; 
	char fishy[10];
	std::sprintf(fishy," %d",node);
	std::strcat(foo,fishy);
	print_count++;
      }
      postream << foo << "\n";
    }
  }

  postream << my_proc_id << ":HEY!!!! Here is acme_to_host_node_map." << "\n";
  for ( i = 0; i < number_acme_nodes ; ++i) { 
    postream << my_proc_id << ":   acme id:" << i << " host code node (c#): " 
	 << acme_to_host_node_map[i] << "\n";
  }

  postream << my_proc_id << ":TEST> *** The new global node ids..." << "\n";
  for ( i = 0; i < number_acme_nodes; ++i){
    postream << my_proc_id << ":TEST>      node: " << i << " HID: "
	 << new_global_node_ids[2*i] << "," << new_global_node_ids[2*i+1]
	 << "  EXID: " << new_exodus_node_ids[i] << "\n";
  }
  postream.flush();

  // Here we want to check the node map and make sure that it is consistsent.
  // Do this by using it to convert new connectivity to old connectivity 
  // and seeing if it all works out
#ifdef KHB_NEEDS_FIXING
  // This won't actually work if we have nodes in multiple node blocks
  //   The check (conn_node >= number_host_code_nodes) is not std::right in this
  //   case.  Arne and I agreed to shut this off for now.
  postream << my_proc_id << ":HEY!!!! Checking node map." << "\n";
  int conn_count = 0;
  for ( blk = 0; blk < number_face_blocks; blk ++) {
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk]);
    POSTCONDITION(num_face_nodes>0);
    if (ContactSearch::Is_a_Shell_Face(face_block_types[blk]) ) {
      for ( int face = 0; face < number_faces_per_block[blk]; face ++) {
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int conn_node = face_connectivity[conn_count]-1; // -1 for fortran
	  if ( conn_node >= number_host_code_nodes) 
	    conn_node = acme_to_host_node_map[conn_node];
	  int orig_conn_node = orig_face_connectivity[conn_count]-1;
	  if ( conn_node != orig_conn_node) {
	    have_an_error++;
	    contact_global_error_check( have_an_error, SearchComm );
	    topology->errors->Add_Error_Message("Shell Connectivity Problems");
	    error_code = ContactSearch::INTERNAL_ERROR;
	    return;
	  }
	  conn_count++;
	}
      }
    }
    else {
      conn_count += num_face_nodes*number_faces_per_block[blk];
    }
  }
  postream.flush();
  if( contact_global_error_check( have_an_error, SearchComm)  >= 1 ){
    error_code = ContactSearch::INTERNAL_ERROR;
    return;
  }
#endif
#endif

#ifndef CONTACT_NO_MPI
  comm_plan = NULL;
  comm_plan_to_owner = NULL;
  comm_plan_to_ghost = NULL;
  import_loft_face_data_comm_plan = NULL;
  number_lofting_ghosts = NULL;
  offset_to_buffer_offsets = NULL;
  buffer_offsets = NULL;
#endif

}


//------------------------------------------------------------------
//
//    compute_inv_conn -- this method computes the inverse
//        connectivity, i.e. the elements connected to each
//        node, given a face_connectivity and its offsets.
//        Also set up arrays to hold the block and face ids for
//        all shell faces.
//
//------------------------------------------------------------------
void ContactShellHandler::compute_inv_conn(
				     const int & num_host_code_nodes,
				     const int & number_face_blocks,
				     const int * face_connectivity,
				     const int * conn_offsets,
    const ContactSearch::ContactFace_Type* face_block_types,
				     const int * number_faces_per_block,
				     const double * face_lofting_factors,
				     const int& number_comm_partners,
				     const int* comm_proc_ids,
				     const int* number_nodes_to_partner,
				     const int* comm_nodes,
					   int *& inv_conn,
					   int *& offset_inv_conn,
					   int *& inv_conn_count,
					   int *& shell_face_blk_list,
					   int *& shell_face_id_list,
					   double *& shell_loft_factors,
					   int & num_shell_faces,
					   int & num_shell_nodes,
				           bool parallel_consistent,
				           ContactSearch * search ) 
{
  //---------------------------------------------------------------------
  // we need to create the inverse connectivity from the
  // nodes to the shell faces. First count the number of faces attached to
  // each of the shell nodes, allocate space, then store the inverse 
  // connectivity. 
 
  int i, blk;
   
  // count number of shell faces connected to each node 
  inv_conn_count = new int[num_host_code_nodes];
  for ( i = 0; i < num_host_code_nodes; ++i) inv_conn_count[i] = 0;
  
  int conn_count = 0;
  num_shell_faces = 0;
  for ( blk = 0; blk < number_face_blocks; blk ++) {
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk]);
    POSTCONDITION(num_face_nodes>0);
    if (ContactSearch::Is_a_Shell_Face(face_block_types[blk]) ) {
      for ( int face = 0; face < number_faces_per_block[blk]; face ++) {
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int node = face_connectivity[conn_count]-1; // -1 for fortran
	  conn_count ++;
	  if (node < 0 ) continue; // skip nodes < 0 
	  inv_conn_count[node]++;
	}
      }
      num_shell_faces += number_faces_per_block[blk];
    }
    else {
      conn_count += num_face_nodes*number_faces_per_block[blk];
    }
  }

  // if parallel_consistent is set, then we want to communicate to all
  // processors which shared nodes are shell nodes. Once we know this, 
  // we can properly include all non-shell faces that are attached to shell
  // nodes
  int * shared_shell_node_flag = new int[num_host_code_nodes];
  for (i = 0; i < num_host_code_nodes; ++i) shared_shell_node_flag[i] = 0;
#ifndef CONTACT_NO_MPI
  int j;
  if ( parallel_consistent ) {
    int num_comm_nodes = 0;
    for (i = 0; i < number_comm_partners; ++i) {
      for (j = 0; j < number_nodes_to_partner[i]; ++j) {
	num_comm_nodes ++;
      }
    }
    
    // set flag for shared shell nodes & setup comm structures
    RequestHandle * recv_handles = new RequestHandle[number_comm_partners];
    char * recv_shell_nodes_buff = new char[num_comm_nodes*sizeof(int)];
    char * send_shell_nodes_buff = new char[num_comm_nodes*sizeof(int)];
    int *recv_shell_nodes, *send_shell_nodes;
    recv_shell_nodes = reinterpret_cast<int *>(recv_shell_nodes_buff);
    send_shell_nodes = reinterpret_cast<int *>(send_shell_nodes_buff);
    int idx = 0;
    for (i = 0; i < num_comm_nodes; ++i) send_shell_nodes[i] = 0;
    for (i = 0; i < number_comm_partners; ++i) {
      for (j = 0; j < number_nodes_to_partner[i]; ++j) {
	int node = comm_nodes[idx]-1; // -1 for fortran
	if ( inv_conn_count[node] > 0 ) {
	  shared_shell_node_flag[node] = 1;
	  send_shell_nodes[idx] = 1;
	}
	idx++;
      }
    }
    
    // do communication
    ContactAsymComm node_comm(number_comm_partners,
			      number_comm_partners,
			      number_nodes_to_partner,
			      number_nodes_to_partner,
			      comm_proc_ids,
			      comm_proc_ids);
    
    contact_communicate_packed_buffers(SearchComm,
				       node_comm,
				       sizeof(int),
				       send_shell_nodes_buff,
				       recv_shell_nodes_buff,
				       recv_handles);
    // resolve all shared nodes
    idx = 0;
    for (i = 0; i < number_comm_partners; ++i) {
      for (j = 0; j < number_nodes_to_partner[i]; ++j) {
	int node = comm_nodes[idx]-1; // -1 for fortran
	if ( recv_shell_nodes[idx] > 0 ) shared_shell_node_flag[node] = 1;
	idx++;
      }
    }

    delete [] recv_handles;
    delete [] recv_shell_nodes_buff;
    delete [] send_shell_nodes_buff;
  }
#endif
  
  // loop over non-shell faces. If any have a node which is a shell node,
  // then include the face in the inv_conn list.   Note that this brings up
  // a tricky situation in parallel -- may not know the non-shell face is
  // attached to a shell face. Will have to resolve this later.
  conn_count = 0;
  for ( blk = 0; blk < number_face_blocks; blk ++) {
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk]);
    POSTCONDITION(num_face_nodes>0);
    if (! ContactSearch::Is_a_Shell_Face(face_block_types[blk]) ) {
      for ( int face = 0; face < number_faces_per_block[blk]; face ++) {
	bool face_attached_to_shell = false;
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int node = face_connectivity[conn_count]-1; // -1 for fortran
	  conn_count ++;
	  if (node < 0 ) continue; // skip nodes < 0 
	  if (inv_conn_count[node] > 0 || shared_shell_node_flag[node] > 0) {
	    inv_conn_count[node]++;
	    face_attached_to_shell = true;
	  }
	}
	if (face_attached_to_shell) num_shell_faces++;
      }
    }
    else {
      conn_count += num_face_nodes*number_faces_per_block[blk];
    }
  }

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> number of local shell faces: " 
       << num_shell_faces << "\n";
  postream << my_proc_id << ":TEST> inv_conn_count:" << "\n";
  for ( i = 0; i < num_host_code_nodes; ++i) {
    postream << my_proc_id << ":TEST>   node " << i << " count " 
	 << inv_conn_count[i] << "\n";
  }
#endif

  // sum up all the node->face connectivities,  This is the size
  // required for the inverse connectivity arrays. Also compute an 
  // offset array so we know where to put data for a specific node.
  num_shell_nodes = 0;
  int tot_inv_conn_count = 0;
  offset_inv_conn = new int[num_host_code_nodes];
  for ( i = 0; i < num_host_code_nodes; ++i)
    offset_inv_conn[i] = -1;
  int * offset_counters = new int[num_host_code_nodes];
  for ( int inode = 0; inode < num_host_code_nodes; inode++) {
    if ( inv_conn_count[inode] > 0 ) num_shell_nodes++;
    offset_inv_conn[inode] = tot_inv_conn_count;
    offset_counters[inode] = tot_inv_conn_count;
    tot_inv_conn_count += inv_conn_count[inode];
  }

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> num_shell_nodes: " << num_shell_nodes << "\n";
#endif

  // allocate arrays for the inverse connectivity and for the mapping
  // from user-input face to shell face (neglecting non-shell faces).
  // if we passed in face_lofting_factors, then we want the mapping
  // for the face loft factors as well (only for the original call
  // to compute_inv_conn -- the lofting factors are dealt with 
  // differently for the ghosted versions). 
  inv_conn = new int[tot_inv_conn_count];
  shell_face_id_list  = new int[num_shell_faces];
  shell_face_blk_list = new int[num_shell_faces];
  if (NULL != face_lofting_factors)
    shell_loft_factors  = new double[num_shell_faces];
  
  // fill inverse connectivity.  The face number is the number of the
  // face in concatinated shell numbering (ignoring non-shell faces).
  //  postream << "**** do inv conn" << "\n";
  int shell_face_count = 0;
  int host_face_count = 0;
  conn_count = 0;
  for ( blk = 0; blk < number_face_blocks; blk ++) {
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk]);
    POSTCONDITION(num_face_nodes>0);
    //    postream << " = blk " << blk << " num face nodes:" << num_face_nodes << "\n";
    if (ContactSearch::Is_a_Shell_Face(face_block_types[blk]) ) {
      //      postream << "   * its a shell blk " << "\n"; 
      // these faces are all shell faces -- add to inv_conn
      for ( int face = 0; 
	    face < number_faces_per_block[blk]; 
	    face++, host_face_count++) {
	//	postream << "     - face " << face << "\n"; 
	shell_face_id_list[shell_face_count] = face;
	shell_face_blk_list[shell_face_count] = blk;
	// only process lofting factors if they were passed in
	if (NULL != face_lofting_factors)
	  shell_loft_factors[shell_face_count] = face_lofting_factors[host_face_count];
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int node = face_connectivity[conn_count] - 1; // -1 for fortran
	  //	  postream << "       ~ node " << node << "\n"; 
	  conn_count ++;
	  if (node < 0 ) continue; // skip nodes < 0 
	  PRECONDITION( offset_counters[node]-offset_inv_conn[node] 
			< inv_conn_count[node]);
	  inv_conn[offset_counters[node]]= shell_face_count;
	  //	  postream << "         inv_conn(" << offset_counters[node]
	  //		   << ")=" << shell_face_count << "\n"; 
	  offset_counters[node]++;
	}
	shell_face_count ++;
      }
    }
    else {
      //      postream << "   * its a non-shell blk " << "\n"; 
      // these are non-shell faces. Only add to inv_conn if one of their
      // nodes is a shell node.
      for ( int face = 0; 
	    face < number_faces_per_block[blk]; 
	    face ++, host_face_count++) {
	bool face_attached_to_shell_node = false;
	for ( int inode = 0; inode < num_face_nodes; inode ++) {
	  int node = face_connectivity[conn_count+inode] - 1; // -1 for fortran
	  if (node < 0 ) continue; // skip nodes < 0 
	  PRECONDITION (node < num_host_code_nodes && node >= 0);
	  if ( inv_conn_count[node] > 0 || shared_shell_node_flag[node] > 0) {
	    face_attached_to_shell_node = true;
	    break;
	  }
	}
	if (face_attached_to_shell_node) {
	  //	  postream << "     - face " << face << "\n"; 
	  shell_face_id_list[shell_face_count] = face;
	  shell_face_blk_list[shell_face_count] = blk;
	  // only process lofting factors if they were passed in
	  if (NULL != face_lofting_factors)
	    shell_loft_factors[shell_face_count] = face_lofting_factors[host_face_count];
	  for ( int inode = 0; inode < num_face_nodes; inode ++) {
	    int node = face_connectivity[conn_count] - 1; // -1 for fortran
	    //	    postream << "       ~ node " << node << "\n"; 
	    conn_count ++;
	    if (node < 0 ) continue; // skip nodes < 0 
	    if (inv_conn_count[node] <= 0 && shared_shell_node_flag[node] <= 0)
	      continue; // skip non-shell nodes
	    PRECONDITION( offset_counters[node]-offset_inv_conn[node] 
			  < inv_conn_count[node]);
	    inv_conn[offset_counters[node]]= shell_face_count;
	    //	    postream << "         inv_conn(" << offset_counters[node]
	    //		     << ")=" << shell_face_count << "\n"; 
	    offset_counters[node]++;
	  }
	  shell_face_count ++;
	} else {
	  //	  postream << "     - SKIPPED face " << face << "\n"; 
	  conn_count += num_face_nodes;
	}
      }
    }
  }

  delete [] offset_counters;
  delete [] shared_shell_node_flag;
}

//------------------------------------------------------------------
//
//    process_face -- this method processes a shell face. It 
//         identifies all of its nodes and creates any that are
//         needed. It also completes all of the shell node data
//         structures needed by contact for shells.
//         This routine also loops over all of the faces connected
//         to the nodes of this face and makes sure that they also
//         connect to the proper nodes.
//         Finally, it loops over all adjacent faces and calls this
//         routine on them, thus recursively covering all faces
//         on this contact surface.
//
//------------------------------------------------------------------
bool ContactShellHandler::process_face(const int shell_face,
				       bool * face_process_flag,
				       bool * face_node_process_flag,
      const ContactSearch::ContactFace_Type* face_block_types,
				       int * face_conn_changeable,
				       int * face_conn_reference,
				 const int * conn_offsets,
				 const int * inv_conn,
				 const int * offset_inv_conn,
				 const int * inv_conn_count,
				 const int * shell_face_blk_list,
				 const int * shell_face_id_list,
				 const double * shell_loft_factors,
				       int * num_nds_per_blk,
				 const double * coords,
				       ContactSearch * search ) 
{
  // this routine processes a shell face

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>  >> Processing face " << shell_face << "\n";
#endif
  
  PRECONDITION (shell_face < m_full_number_faces); //bad face id in process_face
  
  // NOTE: ASG HACK
  //   assumes maximum of 4 nodes and edges per face.
  
  int i, j, k,num_adj_faces(0);
  int adj_faces[4];
  
  // if this face was already processed, return
  if (face_process_flag[shell_face]) {
#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>   -> skipping face cuz its been processed." << "\n";
#endif
    return true;
  }

  // get id and blk of shell face in host-code provided data
  int face_id = shell_face_id_list[shell_face];
  int face_blk = shell_face_blk_list[shell_face];
  int num_face_nodes = 
    search->Number_Nodes_Per_Face(face_block_types[face_blk]);
  POSTCONDITION(num_face_nodes>0);
  int conn_base_idx = conn_offsets[face_blk] + face_id*num_face_nodes;
  
  //-----------------------------------------------------
  // loop over the nodes on the face and resolve its nodes
  for ( i = 0; i < num_face_nodes; ++i) {
    
    int new_node_id;
    int node = face_conn_reference[conn_base_idx+i]-1; // -1 for fortran
    
#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>    >> loop over face node" << i << ":" << node << "\n";
#endif

    // skip node if already processed on this face
    // ASG HACK -- max 4 nodes/face assumption
    if (face_node_process_flag[shell_face*4 + i]) continue;
    
    //------------------------------------
    // setup node data structures and create node if needed
    if (num_created_shell_nodes[node] == 0){
      // node not touched by other face.
      
#ifdef SHELL_DEBUG
      postream << my_proc_id << ":TEST>     -> use original node" << "\n";
#endif

      num_created_shell_nodes[node]++;
      new_node_id = node;
    } 
    else {
      // node has been touched by other face. Create a new one
      
      num_created_shell_nodes[node]++;
      new_node_id = next_new_id++;
#ifdef SHELL_DEBUG
      postream << my_proc_id << ":TEST>     -> ** create node " 
               << new_node_id << " from " << node << "\n";
#endif
      Create_Node( node, new_node_id );
      face_conn_changeable[conn_base_idx + i] 
	= new_node_id + 1;// +1 for fortran
    }
    // ASG HACK -- max 4 edges/face assumption
    face_node_process_flag[shell_face*4 + i] = true;
    
    //------------------------------------
    // for each edge on this face that contains the current node we are
    // resolving (hereby denoted as the anchor node), 

    // loop over edges on face
    // ASG HACK -- same number edges/face assumption
    for (j = 0; j < num_face_nodes; ++j){
      int node1, node2;
      if ( j == num_face_nodes-1 ) {
        node1 = face_conn_reference[conn_base_idx + num_face_nodes-1] - 1; 
	                                 // -1 for fortran
        node2 = face_conn_reference[conn_base_idx] - 1; // -1 for fortran
      } else {		 
        node1 = face_conn_reference[conn_base_idx+j] - 1; // -1 for fortran
        node2 = face_conn_reference[conn_base_idx+j+1] -1; // -1 for fortran
      }

      // skip edge if no anchor
      if ( node1 != node && node2 != node) {
	//	postream << my_proc_id << ":                  ->skip | 1:" << node1
	//	     << "  2:" << node2 << "  anchor:" << node << "\n";
        continue; 
      }
      
      // call routine to process edge -- returns -1 if no faces are connected
      int adj_face = process_edge(shell_face,
				  node,
				  new_node_id,
				  node1, 
				  node2,
				  inv_conn,
				  offset_inv_conn,
				  inv_conn_count,
				  conn_offsets,
				  shell_face_blk_list,
				  shell_face_id_list,
				  shell_loft_factors,
				  face_conn_changeable,
				  face_conn_reference,
				  face_block_types,
				  face_node_process_flag,
				  coords,
				  search ); 

#ifdef SHELL_DEBUG
      // std::flush to prevent postream from getting too big and crashing
      postream.flush();
#endif
      
      if ( adj_face >= 0) {
        // make sure adj_face is not already in list of adjacent
        // faces
        bool found = false;
        for (k = 0; k < num_adj_faces; ++k) {
          if ( adj_faces[k] == adj_face) {
            found = true;
            break;
          }
        }
        if (!found) {
          PRECONDITION(num_adj_faces < 4);
          adj_faces[num_adj_faces] = adj_face;
          num_adj_faces++;
        }
      } else if ( adj_face == -2 ){
        return(false);
      }
    } // end loop on edges of face anchored at node

  } // end loop on nodes of face

  // mark face as processed.
  face_process_flag[shell_face] = true;

  // print out resulting connectivity
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>   *** connectivity for face " << shell_face << "\n";
  char foo[80];
  std::sprintf(foo,"%d:TEST>     -> ",my_proc_id);
  for ( i = 0; i < num_face_nodes; ++i) {
    int node = face_conn_changeable[conn_base_idx+i]-1;
    char fishy[10];
    std::sprintf(fishy," %d",node);
    std::strcat(foo,fishy);
  }
  postream << foo << "\n";
  postream.flush();
#endif    

  //-----------------------------------------------------
  // loop over faces we have touched while resolving the edges. Process
  // each one.
  for ( i =0; i < num_adj_faces; ++i) {
    faces_to_process.push_back(adj_faces[i]);
  }

  return true;
}

//------------------------------------------------------------------
//
//    process_edge -- this routine takes a face and an edge and finds
//         its adjacent face on the contact surface. It marks the 
//         anchor node of the edge as resolved on the new face. It then 
//         finds the edge of the adjacent face that is connected to one 
//         of the nodes on the original edge and calls process_edge to 
//         find its adjacent face. This routine thus drives the process
//         of walking the faces attached to the anchor node.
//
//         NOTE: a return value of -2 indicates a problem
//------------------------------------------------------------------
int ContactShellHandler::process_edge(
				const int shell_face,
				const int anchor_node,
				const int shell_node,
				const int node1, 
				const int node2,
				const int * inv_conn,
				const int * offset_inv_conn,
				const int * inv_conn_count,
				const int * conn_offsets,
				const int * shell_face_blk_list,
				const int * shell_face_id_list,
	         	        const double * shell_loft_factors,
				      int * face_conn_changeable,
				      int * face_conn_reference,
     const ContactSearch::ContactFace_Type* face_block_types,
				      bool * face_node_process_flag,
				const double * coords,
				ContactSearch * search ) 
{
  PRECONDITION( shell_face >= 0);
  PRECONDITION( shell_node >= 0);
  PRECONDITION( node1 >= 0);
  PRECONDITION( node2 >= 0);
  PRECONDITION( anchor_node == node1 || anchor_node == node2 );
  int i, j;

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>     >> process edge " << node1 << " & " << node2 
       << " on face " << shell_face << "\n"
       << "              anchor: " << anchor_node << " currentid:" 
       << shell_node << "\n";
#endif

  int blk = shell_face_blk_list[shell_face];
  // if face passed in does not yet have its opposite determined, then 
  // find it. If it is a non-shell face, just leave its opposite as -1
  if (! ContactSearch::Is_a_Shell_Face(face_block_types[blk]) ) {
    opposite_face[shell_face] = -1;
  } else {
    if (opposite_face[shell_face] < 0) {
      bool found_opposite = false;
      found_opposite = find_opposite(shell_face,
		      face_block_types,
		      conn_offsets,
		      face_conn_reference,
		      offset_inv_conn,
		      inv_conn_count,
		      inv_conn,
		      shell_face_blk_list,
		      shell_face_id_list,
		      search);
      if (! found_opposite ) {
	// if shell doesn't have an opposite, just mark it as -1.
	opposite_face[shell_face] = -1;
      }
    }
  }
  int opposite = opposite_face[shell_face];
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>      -> opposite face: " << opposite << "\n";
#endif

  int num_faces_on_anchor = inv_conn_count[anchor_node];
  int inv_conn_idx = offset_inv_conn[anchor_node];
  int num_possible_adj_faces(0);
  int * adj_faces = new int [num_faces_on_anchor];
  int * adj_face_idxes = new int [num_faces_on_anchor];
  int * adj_edges = new int [num_faces_on_anchor];
  int * adj_num_face_nodes_for_face = new int [num_faces_on_anchor];
  int * anchor_idxes = new int [num_faces_on_anchor];

  // loop over the faces connected to the anchor node
  for ( i = 0; i < num_faces_on_anchor; ++i) {
    int attached_face = inv_conn[inv_conn_idx+i];
    if ( attached_face == shell_face) continue; // same as face we passed in
    if ( attached_face == opposite) continue; // same as face we passed in
    int attached_face_id = shell_face_id_list[attached_face];
    int attached_face_blk = shell_face_blk_list[attached_face];
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[attached_face_blk]);
    POSTCONDITION(num_face_nodes>0);
    int attached_face_idx = 
      conn_offsets[attached_face_blk] + attached_face_id*num_face_nodes;
    
    // loop over edges of attached face
    // ASG HACK -- same number edges/face assumption
    for ( j = 0; j < num_face_nodes; ++j) {
      int edge_node1, edge_node2;
      if ( j == num_face_nodes-1 ) {
	edge_node1 
	  = face_conn_reference[attached_face_idx+num_face_nodes-1] - 1;
	                                      // -1 for fortran
	edge_node2 
	  = face_conn_reference[attached_face_idx]-1; // -1 for fortran
      } else {	     	      
	edge_node1 
	  = face_conn_reference[attached_face_idx+j]-1; // -1 for fortran
	edge_node2 
	  = face_conn_reference[attached_face_idx+j+1]-1; //-1 for fortran
      }

      // skip edges which do not have the opposite direction
      if ( edge_node1 != node2 || edge_node2 != node1) continue;
 
#ifdef SHELL_DEBUG
      postream << my_proc_id << ":TEST>      -> found connected face: " << attached_face << "\n";
#endif

      // found attached edge. set flags and keep looking.

      if ( edge_node1 == anchor_node ) {
	anchor_idxes[num_possible_adj_faces] = j;
      } else if ( edge_node2 == anchor_node ) {
	if ( j == num_face_nodes-1 ) anchor_idxes[num_possible_adj_faces]=0;
	else anchor_idxes[num_possible_adj_faces] = j+1;
      }
      POSTCONDITION(anchor_idxes[num_possible_adj_faces] < num_face_nodes);

      adj_faces[num_possible_adj_faces] = attached_face;
      adj_face_idxes[num_possible_adj_faces] = attached_face_idx;
      adj_edges[num_possible_adj_faces] = j;
      adj_num_face_nodes_for_face[num_possible_adj_faces] = num_face_nodes;
      num_possible_adj_faces ++;
      break;
    } // end loop over edges
    
    
  } // end loop over faces attached to anchor_node

  // if no adjacent faces were found, return out of routine
  if ( num_possible_adj_faces == 0 ){
    delete [] anchor_idxes;
    delete [] adj_num_face_nodes_for_face;
    delete [] adj_edges;
    delete [] adj_face_idxes;
    delete [] adj_faces;
    return -1;
  }

  int joining_edge;
  if ( num_possible_adj_faces == 1 ) {
    // only one edge matches -- choose it
    joining_edge = 0;
  } else {
    // More than one face found that shared the edge -- we have an N-way 
    // section. Need to call a routine to determine which face is 
    // closest
    joining_edge = find_best_face(node1, 
				  node2,
				  shell_face,
				  face_block_types,
				  conn_offsets,
				  face_conn_reference,
				  shell_face_blk_list,
				  shell_face_id_list,
				  shell_loft_factors,
				  num_possible_adj_faces,
				  adj_faces,
				  coords,
				  search );
    if ( joining_edge == -2 ){
      delete [] anchor_idxes;
      delete [] adj_num_face_nodes_for_face;
      delete [] adj_edges;
      delete [] adj_face_idxes;
      delete [] adj_faces;
      return -2;
    }
  }

  // copy identified edge and 
  int adj_face = adj_faces[joining_edge];
  int adj_face_idx = adj_face_idxes[joining_edge];
  int adj_edge = adj_edges[joining_edge];
  int adj_num_face_nodes = adj_num_face_nodes_for_face[joining_edge];
  int anchor_idx = anchor_idxes[joining_edge];

  delete [] anchor_idxes;
  delete [] adj_num_face_nodes_for_face;
  delete [] adj_edges;
  delete [] adj_face_idxes;
  delete [] adj_faces;
  
  // some asserts for valid return values if we found an adjacent face      
  POSTCONDITION(adj_edge >=0);
  POSTCONDITION(adj_num_face_nodes >=0);
  POSTCONDITION(adj_face >=0);
  POSTCONDITION(adj_face_idx >=0);
  POSTCONDITION(anchor_idx >=0);

  // if anchor node on found face was already processed, then 
  // we have circled back to face that started processing this
  // node. return.
  // ASG HACK -- max 4 nodes/face assumption
  if ( face_node_process_flag[adj_face*4+anchor_idx] )
    return adj_face;

  // set anchor node on found face 
  // ASG HACK -- max 4 nodes/face assumption
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>      -> fixing local node idx: " << anchor_idx << "\n";
#endif
  face_node_process_flag[adj_face*4+anchor_idx] = true;
  if ( shell_node != anchor_node ) 
    face_conn_changeable[adj_face_idx + anchor_idx] 
      = shell_node + 1; // +1 for fortran

  // we found the face which shares an edge with original face passed in 
  // (shell_face). Now find the other edge on this face that has the same
  // anchor node.
  for ( i = 0; i < adj_num_face_nodes; ++i) {
    if ( i == adj_edge) continue; // skip edge adjacent to passed-in face
    int edge_node1, edge_node2;
    if ( i == adj_num_face_nodes-1 ) {
      edge_node1 = face_conn_reference[adj_face_idx+adj_num_face_nodes-1] - 1;
                                // -1 for fortran
      edge_node2 = face_conn_reference[adj_face_idx]-1; // -1 for fortran
    } else {	   	    
      edge_node1 = face_conn_reference[adj_face_idx+i]-1; // -1 for fortran
      edge_node2 = face_conn_reference[adj_face_idx+i+1]-1; //-1 for fortran
    }
    
    // skip edges which do not have anchor_node
    if ( edge_node1 != anchor_node && edge_node2 != anchor_node) continue;
    
    // found the other edges with this anchor node. process it.

#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>      -> opposite edge: " << edge_node1 
	 << " & " << edge_node2 << "\n";
#endif

    int return_value;
    return_value = process_edge(adj_face,
				anchor_node,
				shell_node,
				edge_node1, 
				edge_node2,
				inv_conn,
				offset_inv_conn,
				inv_conn_count,
				conn_offsets,
				shell_face_blk_list,
				shell_face_id_list,
				shell_loft_factors,
				face_conn_changeable,
				face_conn_reference,
				face_block_types,
				face_node_process_flag,
				coords,
				search ); 
    // if error was found in processing, return
    if (return_value == -2) {
      return -2;
    }
  } // end loop on adj_face edges
    
  return adj_face;
}


//------------------------------------------------------------------
//
//   find_best_face -- 
//
//------------------------------------------------------------------
int ContactShellHandler::find_best_face(int node1,
					int node2,
					int face,
       const ContactSearch::ContactFace_Type* face_block_types,
					const int * conn_offsets,
					const int * face_connectivity,
					const int * shell_face_blk_list,
					const int * shell_face_id_list,
		                        const double * shell_loft_factors,
					int num_adj_faces,
					const int * adj_faces,
					const double * coords,
					ContactSearch * search)
{  
  int i,j;

  double pi = 2.0 * std::asin(1.0);

#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>          -> in best face with : " 
	     << node1 << " & " << node2 <<  " on face " << face << "\n";
#endif
    
  // compute midpoint of edge
  double edge_midpoint[3];
  edge_midpoint[0] = (coords[node1*3]   + coords[node2*3])/2.0;
  edge_midpoint[1] = (coords[node1*3+1] + coords[node2*3+1])/2.0;
  edge_midpoint[2] = (coords[node1*3+2] + coords[node2*3+2])/2.0;

  // compute vector direction of edge + its magnitude
  double edge_dir[3];
  edge_dir[0] = coords[node2*3]   - coords[node1*3];
  edge_dir[1] = coords[node2*3+1] - coords[node1*3+1];
  edge_dir[2] = coords[node2*3+2] - coords[node1*3+2];
  double edge_dir_mag = Magnitude(edge_dir);
  
  // compute midpoint of current face
  int face_id = shell_face_id_list[face];
  int face_blk = shell_face_blk_list[face];
  int num_face_nodes = 
    search->Number_Nodes_Per_Face(face_block_types[face_blk]);
  POSTCONDITION(num_face_nodes>0);
  int face_idx = 
    conn_offsets[face_blk] + face_id*num_face_nodes;
  double face_midpoint[3];
  face_midpoint[0] = 0.0;
  face_midpoint[1] = 0.0;
  face_midpoint[2] = 0.0;
  for ( j = 0; j < num_face_nodes; ++j) {
    int node = face_connectivity[face_idx+j] -1; // -1 for fortran
    face_midpoint[0] += coords[3*node]/num_face_nodes;
    face_midpoint[1] += coords[3*node + 1]/num_face_nodes;
    face_midpoint[2] += coords[3*node + 2]/num_face_nodes;
  }
  
  // compute base vector -- midpoint of face - midpoint of edge
  double base_vector[3];
  base_vector[0] = face_midpoint[0] - edge_midpoint[0];
  base_vector[1] = face_midpoint[1] - edge_midpoint[1];
  base_vector[2] = face_midpoint[2] - edge_midpoint[2];

  // make base vector perpendicular to edge vector
  double base_edge_dot = Dot(base_vector, edge_dir);

  double scale = base_edge_dot / (edge_dir_mag * edge_dir_mag);
  base_vector[0] = base_vector[0] - scale * edge_dir[0];
  base_vector[1] = base_vector[1] - scale * edge_dir[1];
  base_vector[2] = base_vector[2] - scale * edge_dir[2];

  // compute mag of base vector
  double base_vec_mag = Magnitude(base_vector);
  
  POSTCONDITION( Dot(base_vector,edge_dir) / (base_vec_mag * edge_dir_mag) < 1.0e-4);

  // setup best values
  int best_face = -1;
  double best_angle = 3.0 * pi;
  
  // loop over adjacent faces
  for ( i = 0; i < num_adj_faces; i ++ ){

    // get data for adj face
    int adj_face = adj_faces[i];
    int adj_face_id = shell_face_id_list[adj_face];
    int adj_face_blk = shell_face_blk_list[adj_face];
    num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[adj_face_blk]);
    POSTCONDITION(num_face_nodes>0);
    int adj_face_idx = 
      conn_offsets[adj_face_blk] + adj_face_id*num_face_nodes;
    
#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>            -> adj face: " 
	     << i << " is " << adj_face << "\n";
#endif
    
    // compute midpoint of adjacent face
    double adj_face_midpoint[3];
    adj_face_midpoint[0] = 0.0;
    adj_face_midpoint[1] = 0.0;
    adj_face_midpoint[2] = 0.0;
    for ( j = 0; j < num_face_nodes; ++j) {
      int node = face_connectivity[adj_face_idx+j] -1; // -1 for fortran
      adj_face_midpoint[0] += coords[3*node]/num_face_nodes;
      adj_face_midpoint[1] += coords[3*node + 1]/num_face_nodes;
      adj_face_midpoint[2] += coords[3*node + 2]/num_face_nodes;
    }
    
    // compute vector -- midpoint of adjacent face - midpoint of edge
    double vector[3];
    vector[0] = adj_face_midpoint[0] - edge_midpoint[0];
    vector[1] = adj_face_midpoint[1] - edge_midpoint[1];
    vector[2] = adj_face_midpoint[2] - edge_midpoint[2];

    // make vector perpendicular to edge vector
    double vec_edge_dot = Dot(vector, edge_dir);

    scale = vec_edge_dot / (edge_dir_mag * edge_dir_mag);
    vector[0] = vector[0] - scale * edge_dir[0];
    vector[1] = vector[1] - scale * edge_dir[1];
    vector[2] = vector[2] - scale * edge_dir[2];
    
    // compute mag of vector
    double vec_mag = Magnitude(vector);

    POSTCONDITION( Dot(vector,edge_dir) / (vec_mag * edge_dir_mag) < 1.0e-4);
    
    // get dot product between base vector and vector
    double ang_dot = Dot(base_vector,vector);
    double ang_dot_normed = ang_dot / (base_vec_mag * vec_mag);

    double angle;
    if ( 1.0 - std::fabs(ang_dot_normed) > 1.0e-8 ) {
      // there is enough of an angle between base vector and vector to 
      // compute.
      
      // get cross product between base vector and vector
      double cross_vec[3];
      Cross(base_vector, vector, cross_vec);
 
      double cross_vec_mag = Magnitude(cross_vec);

      // get dot product between cross and edge dir. Note that 
      // the normalization done here just makes it nondimensional
      // which is good when comparing to constant tolerances.
      double dir_dot = Dot(cross_vec, edge_dir) / (cross_vec_mag * edge_dir_mag);

      POSTCONDITION(std::fabs( 1.0 - std::fabs(dir_dot) ) < 1.0e-6 );

      // compute angle
      if ( dir_dot > 0.0) {
	// cross is positive -- angle between 0 and 180 degrees
	angle = std::acos( ang_dot_normed );
      } else {
	// cross is negative -- angle between 180 and 360 degrees
	angle = 2.0 * pi - std::acos( ang_dot_normed );
      }    
    } else {
      // to avoid roundoff issues, catch numbers close to 1.0
      if ( ang_dot_normed > 0.0) angle = 0.0;
      else angle = pi;
    }

    // if angle is near zero and lofting factors are not correctly
    // opposed, then the zero angle is not optimal -- set to 2 pi
    if ( angle < 1.0e-10) {
      double this_face_loft_fact = shell_loft_factors[face];
      double comp_face_loft_fact = shell_loft_factors[adj_face];

// #ifdef SHELL_DEBUG
//     postream << my_proc_id << ":TEST>            -> angle near zero:" 
// 	     << "this loft fact:" << this_face_loft_fact
// 	     << " new:" << comp_face_loft_fact << "\n";
// #endif

      if ( this_face_loft_fact < 1.0e-10 ){
	// this face loft factor is near zero. if the other face in the
	// same space is also near zero, then the zero angle is o.k. If the
	// other loft factor is not near zero, then we error -- problem in
	// mesh construction.
	if (comp_face_loft_fact > 1.0e-10 ) {
	  topology->errors->
	    Add_Error_Message("Two faces with identical position and ambiguous lofting factors");
	  char msg[81];
	  std::sprintf(msg,"  Faces connected to nodes %d and %d",
		       orig_exodus_node_ids[node1], 
		       orig_exodus_node_ids[node2]);
	  topology->errors->Add_Error_Message(msg);
	  topology->Search()->append_invalid_edge(orig_exodus_node_ids[node1],
						  orig_exodus_node_ids[node2]);
	  return -2;
	}
      } else if ( this_face_loft_fact > 1.0 - 1.0e-10 ){
	// this face loft factor is near one. if the other face in the
	// same space is also near one, then the zero angle should
	// actually be 2 pi -- they are opposite faces. If the
	// other loft factor is not near one, then we error -- problem in
	// mesh construction.
	if (comp_face_loft_fact < 1.0 - 1.0e-10 ) {
	  topology->errors->
	    Add_Error_Message("Two faces have identical position and ambiguous lofting factors");
	  char msg[81];
	  std::sprintf(msg,"  Faces connected to nodes %d and %d",
		       orig_exodus_node_ids[node1], 
		       orig_exodus_node_ids[node2]);
	  topology->errors->Add_Error_Message(msg);
	  topology->Search()->append_invalid_edge(orig_exodus_node_ids[node1],
						  orig_exodus_node_ids[node2]);
	  return -2;
	}
	angle = 2.0 * pi;
#ifdef SHELL_DEBUG
	postream << my_proc_id << ":TEST>              -> new angle 2 pi.\n";
#endif
      }
    }
    
    if ( best_angle - angle > 1.0e-10 ) {

#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>            -> a clear winner:" 
	     << "prev best ang:" << best_angle << " new:" << angle << "\n";
#endif

      // new angle is better than best angle. Replace it.
      best_angle = angle;
      best_face = i;

    } else if ( angle - best_angle < 1.0e-10) {

      // angles are basically identical -- faces are lying on 
      // top of each other. Use lofting factors to determine best match
      if (best_face < 0) {
	topology->errors->
	  Add_Error_Message("Internal error in find_best_face");
	topology->Search()->append_invalid_edge(orig_exodus_node_ids[node1],
						orig_exodus_node_ids[node2]);
	return -2;
      }
      double best_face_loft_fact = shell_loft_factors[adj_faces[best_face]];
      double new_face_loft_fact = shell_loft_factors[adj_face];
      if ( new_face_loft_fact - best_face_loft_fact < 1.0e-5 && 
	   new_face_loft_fact - best_face_loft_fact > -1.0e-5) {
	// two faces are lying on top of each other with no differentiation.
	// std::abort.
	topology->errors->
	  Add_Error_Message("Two faces have identical position and ambiguous lofting factors");
	  char msg[81];
	  std::sprintf(msg,"  Faces connected to nodes %d and %d",
		       orig_exodus_node_ids[node1], 
		       orig_exodus_node_ids[node2]);
	  topology->errors->Add_Error_Message(msg);
	  topology->Search()->append_invalid_edge(orig_exodus_node_ids[node1],
						  orig_exodus_node_ids[node2]);
	return -2;
      } else if (new_face_loft_fact > best_face_loft_fact) {
#ifdef SHELL_DEBUG
    postream << my_proc_id << ":TEST>            -> not clear winner:" 
	     << "prev best ang:" << best_angle << " new:" << angle << "\n";
    postream << my_proc_id << ":TEST>                                " 
	     << "prev best fact:" << best_face_loft_fact 
	     << " new:" << new_face_loft_fact << "\n";
#endif

	// new face has larger lofting factor -- use it as best.
	best_angle = angle;
	best_face = i;
	best_face_loft_fact = new_face_loft_fact;
      }
    }
  }
  
  // return face with minimum angle

  return best_face;
}


//------------------------------------------------------------------
//
//   find_opposite -- given a shell face, finds its opposite.
//       The opposite is the shell face that has the same
//       nodes in its connectivity, but which are in the opposite
//       order. Puts all the opposite data into the member data
//       array "opposite_face".
//
//------------------------------------------------------------------
bool ContactShellHandler::find_opposite(const int shell_face,
       const ContactSearch::ContactFace_Type* face_block_types,
					const int * conn_offsets,
					const int * face_connectivity,
					const int * offset_inv_conn,
					const int * inv_conn_count,
					const int * inv_conn,
					const int * shell_face_blk_list,
					const int * shell_face_id_list,
					ContactSearch * search)
{  
  int j,k;

  // skip if we already found the opposite
  if (opposite_face[shell_face]!=-1) return true;

  //  postream << my_proc_id << ": ***** finding opposite..." << "\n";

  // get host-code face and block from shell face
  int face_id = shell_face_id_list[shell_face];
  int face_blk = shell_face_blk_list[shell_face];
  
  // loop over the faces attached to the first node, check if all
  // its nodes are same as this one
  int num_face_nodes = 
    search->Number_Nodes_Per_Face(face_block_types[face_blk]);
  POSTCONDITION(num_face_nodes>0);
  int conn_idx    = conn_offsets[face_blk] + num_face_nodes*face_id;
  int first_Node  = face_connectivity[conn_idx]-1; // -1 for fortran
  bool face_match = false;
  for (int inode_face = offset_inv_conn[first_Node];
       inode_face < offset_inv_conn[first_Node] + inv_conn_count[first_Node];
       inode_face++ ) {
    face_match = true;
    int shell_face_from_node = inv_conn[inode_face];
    //    postream << my_proc_id << ":  >> check face " << shell_face_from_node << "\n";
    if ( shell_face_from_node == shell_face) continue; // skip if same face
    int face_blk_from_node = shell_face_blk_list[shell_face_from_node]; 
    int face_id_from_node = shell_face_id_list[shell_face_from_node]; 
    //    postream << my_proc_id << ":  >>    with blk,id " 
    //	 << face_blk_from_node << "," << face_id_from_node << "\n";
    if (num_face_nodes != 
	search->Number_Nodes_Per_Face(face_block_types[face_blk_from_node]))
      continue;
    
    //    postream << my_proc_id << ":     -> not same or diff type" << "\n";

    int conn_idx_from_node = conn_offsets[face_blk_from_node] + 
      num_face_nodes*face_id_from_node;

//      postream << my_proc_id << ":      -> check conn:"
//  	 << " (" << conn_idx << ") "
//  	 << face_connectivity[conn_idx+0] << " "
//  	 << face_connectivity[conn_idx+1] << " "
//  	 << face_connectivity[conn_idx+2] << " "
//  	 << face_connectivity[conn_idx+3] << "\n";
//      postream << my_proc_id << ":      -> vs. conn:"
//  	 << " (" << conn_idx_from_node << ") "
//  	 << face_connectivity[conn_idx_from_node+0] << " "
//  	 << face_connectivity[conn_idx_from_node+1] << " "
//  	 << face_connectivity[conn_idx_from_node+2] << " "
//  	 << face_connectivity[conn_idx_from_node+3] << "\n";



    for (j = 0; j < num_face_nodes; ++j) {
      bool node_match = false;
      for (k = 0; k < num_face_nodes; ++k) {
	if (face_connectivity[conn_idx+j] ==
	    face_connectivity[conn_idx_from_node+k]) {
	  node_match = true;
	  break;
	}
      }
      if (!node_match) {
	face_match = false;
	break;
      } 
    }
    if (face_match) {
      opposite_face[shell_face] = shell_face_from_node;
      opposite_face[shell_face_from_node] = shell_face;
      break;
    }
  } // end loop on inv_conn for first node of shell face

  //  postream << my_proc_id << ": ***** done searching for opposite..." << "\n";
  return face_match;
} 

ContactShellHandler::~ContactShellHandler()
{
  delete [] lofted_nodes;
  delete [] host_to_acme_node_map;
  delete [] acme_to_host_node_map;
  delete [] num_created_shell_nodes;
  delete [] num_acme_nodes_for_host_code_node;
  delete [] offset_for_node;
  delete [] block_offset_for_nodes;
  delete [] host_nodes_per_block;

  delete [] opposite_face;
  delete [] shell_face_id;
  delete [] shell_face_blk;
  delete [] shell_loft_factors;
  delete [] is_a_tab_node;

  delete [] orig_node_block_types;
  delete [] orig_number_nodes_per_block;
  delete [] orig_exodus_node_ids;
  delete [] orig_global_node_ids;
  delete [] orig_face_connectivity;
  delete [] orig_comm_proc_ids;
  delete [] orig_number_nodes_to_partner;
  delete [] orig_comm_nodes;

  delete [] new_number_nodes_per_block;
  delete [] new_exodus_node_ids;
  delete [] new_global_node_ids;
  delete [] new_face_connectivity;
  delete [] new_number_nodes_to_partner;
  delete [] new_comm_nodes;

#ifndef CONTACT_NO_MPI
  delete comm_plan;
  delete comm_plan_to_owner;
  delete comm_plan_to_ghost;
  delete import_loft_face_data_comm_plan;
  delete [] number_lofting_ghosts;
  delete [] offset_to_buffer_offsets;
  delete [] buffer_offsets;
#endif

} 


#ifndef CONTACT_NO_MPI
void ContactShellHandler::Build_Shell_Comm_Plans()
{
  // This function creates comm plans that only involve the shell nodes.
  // A total of 3 nodal comm plans are created.
  //   1) A symmetric comm plan that can be used for swapadds.
  //   2) An asymmetric comm plan to send data from the ghosts to the owner
  //   3) An asymmetric comm plan to send data from the owner to the ghosts

  // Mark temp_tag with 1 for shell nodes and 0 for non-shell nodes
  int i;
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    if( Nodes[i]->Is_a_Shell_Node() )
      Nodes[i]->temp_tag = 1;
    else
      Nodes[i]->temp_tag = 0;
  }

#ifdef SHELL_DEBUG
  int j;
  ContactSymComm* node_sym_comm = topology->Node_Sym_Comm();
  if( contact_processor_number(SearchComm) == 0 )
    postream << "Total Comm Plan\n==========================================\n";
  postream << "Number Comm Partners = " << node_sym_comm->Num_Comm_Partners()
	   << "\n";
  for( i=0 ; i<node_sym_comm->Num_Comm_Partners() ; ++i){
    postream << "  Communicating to Processor " 
	     << node_sym_comm->Comm_Proc_ID(i)
	     << "\n";
    for( j=0 ; j<node_sym_comm->Num_to_Proc(i) ; ++j){
      postream << "    Communicating node " 
	       << ((ContactNode<Real>*) node_sym_comm->Entity_List(i)[j])->Exodus_ID()
	       << " " 
	       << ((ContactNode<Real>*) node_sym_comm->Entity_List(i)[j])->Global_ID()
	       << "\n";
    }
  }
  postream.flush();
#endif

  // Build the symmetric comm plan
  comm_plan = new ContactSymComm();
  comm_plan->Build_Subset_Comm_Plan_Using_Temp_Tag(*topology->Node_Sym_Comm());

#ifdef SHELL_DEBUG
  if( contact_processor_number(SearchComm) == 0 )
    postream << "\n\nShell Comm Plan\n======================================\n";
  postream << "Number Comm Partners = " << comm_plan->Num_Comm_Partners()
	   << "\n";
  for( i=0 ; i<comm_plan->Num_Comm_Partners() ; ++i){
    postream << "  Communicating to Processor " << comm_plan->Comm_Proc_ID(i)
	     << "\n";
    for( j=0 ; j<comm_plan->Num_to_Proc(i) ; ++j){
      postream << "    Communicating node " 
	       << ((ContactNode<Real>*) comm_plan->Entity_List(i)[j])->Exodus_ID()
	       << " " 
	       << ((ContactNode<Real>*) comm_plan->Entity_List(i)[j])->Global_ID()
	       << "\n";
    }
  }
  postream.flush();
#endif

  // Build the asymmetric comm plans
  // Note: comm_plan_to_owner is the inverse of comm_plan_to_ghosts 
  //       (i.e., import is export, etc.)
  comm_plan_to_ghost = new ContactAsymComm( *comm_plan );
  comm_plan_to_owner = new 
    ContactAsymComm( comm_plan_to_ghost->Num_Import_Comm_Partners(),
		     comm_plan_to_ghost->Num_Export_Comm_Partners(),
		     comm_plan_to_ghost->Num_Import_from_Procs(),
		     comm_plan_to_ghost->Num_Export_to_Procs(),
		     comm_plan_to_ghost->Import_Comm_Proc_IDs(),
		     comm_plan_to_ghost->Export_Comm_Proc_IDs(),
		     comm_plan_to_ghost->Import_Entity_Lists(),
		     comm_plan_to_ghost->Export_Entity_Lists() );

#ifdef SHELL_DEBUG
  if( contact_processor_number(SearchComm) == 0 ){
    postream << "\n\nComm Plan to Ghost\n";
    postream << "==========================================\n";
  }
  postream << "Export Information\n";
  postream << "   Num Comm Partners = " 
	   << comm_plan_to_ghost->Num_Export_Comm_Partners() << "\n";
  for( i=0 ; i<comm_plan_to_ghost->Num_Export_Comm_Partners() ; ++i){
    postream << "    Exporting to Proc " 
	     << comm_plan_to_ghost->Export_Comm_Proc_ID(i) << "\n";
    for( j=0 ; j<comm_plan_to_ghost->Num_Export_to_Proc(i) ; ++j){
      postream << "       Communicating Node " 
	       << ((ContactNode<Real>*) comm_plan_to_ghost->Export_Entity_List(i)[j])->Exodus_ID()
	       << " " 
	       << ((ContactNode<Real>*) comm_plan_to_ghost->Export_Entity_List(i)[j])->Global_ID()
	       << "\n";
    }
  }
  postream << "Import Information\n";
  postream << "   Num Comm Partners = " 
	   << comm_plan_to_ghost->Num_Import_Comm_Partners() << "\n";
  for( i=0 ; i<comm_plan_to_ghost->Num_Import_Comm_Partners() ; ++i){
    postream << "    Importing from Proc " 
	     << comm_plan_to_ghost->Import_Comm_Proc_ID(i) << "\n";
    for( j=0 ; j<comm_plan_to_ghost->Num_Import_from_Proc(i) ; ++j){
      postream << "       Communicating Node " 
	       << ((ContactNode<Real>*) comm_plan_to_ghost->Import_Entity_List(i)[j])->Exodus_ID()
	       << " " 
	       << ((ContactNode<Real>*) comm_plan_to_ghost->Import_Entity_List(i)[j])->Global_ID()
	       << "\n";
    }
  }
  postream.flush();


  if( contact_processor_number(SearchComm) == 0 ){
    postream << "\n\nComm Plan to Owner\n";
    postream << "==========================================\n";
  }
  postream << "Export Information\n";
  postream << "   Num Comm Partners = " 
	   << comm_plan_to_owner->Num_Export_Comm_Partners() << "\n";
  for( i=0 ; i<comm_plan_to_owner->Num_Export_Comm_Partners() ; ++i){
    postream << "    Exporting to Proc " 
	     << comm_plan_to_owner->Export_Comm_Proc_ID(i) << "\n";
    for( j=0 ; j<comm_plan_to_owner->Num_Export_to_Proc(i) ; ++j){
      postream << "       Communicating Node " 
	       << ((ContactNode<Real>*) comm_plan_to_owner->Export_Entity_List(i)[j])->Exodus_ID()
	       << " " 
	       << ((ContactNode<Real>*) comm_plan_to_owner->Export_Entity_List(i)[j])->Global_ID()
	       << "\n";
    }
  }
  postream << "Import Information\n";
  postream << "   Num Comm Partners = " 
	   << comm_plan_to_owner->Num_Import_Comm_Partners() << "\n";
  for( i=0 ; i<comm_plan_to_owner->Num_Import_Comm_Partners() ; ++i){
    postream << "    Importing from Proc " 
	     << comm_plan_to_owner->Import_Comm_Proc_ID(i) << "\n";
    for( j=0 ; j<comm_plan_to_owner->Num_Import_from_Proc(i) ; ++j){
      postream << "       Communicating Node " 
	       << ((ContactNode<Real>*) comm_plan_to_owner->Import_Entity_List(i)[j])->Exodus_ID()
	       << " " 
	       << ((ContactNode<Real>*) comm_plan_to_owner->Import_Entity_List(i)[j])->Global_ID()
	       << "\n";
    }
  }
  postream.flush();
#endif // SHELL_DEBUG
}
#endif

#ifndef CONTACT_NO_MPI

ContactSearch::ContactErrorCode
ContactShellHandler::ghost_shell_faces(const int& number_comm_partners,
				       const int* comm_proc_ids,
				       const int* number_nodes_to_partner,
				       const int* comm_nodes,
				       const int& num_orig_shell_faces,
				       const int& num_orig_face_blks,
				       const int* inv_conn,
				       const int* inv_conn_count,
				       const int* offset_inv_conn,
				       const int* face_connectivity,
				       const int* conn_offsets,
     const ContactSearch::ContactFace_Type* face_block_types,
				       const int* num_orig_faces_per_blk,
				       const int* global_face_ids,
				       const double* coords,
				       const int* shell_face_blk_list,
				       const int* shell_face_id_list,
				       const double* shell_loft_factors,
				       int & full_number_face_blocks,
				       int*& full_face_connectivity,
				       int*& full_conn_offsets,
           ContactSearch::ContactFace_Type*& full_face_block_types,
				       int*& full_num_faces_per_blk,
				       int*& full_global_face_ids,
				       int&  full_num_host_code_nodes,
				       double *& full_coords,
				       double *& full_shell_loft_factors,
				       ContactSearch* search
				       )
{
//------------------------------------------------------------------
//
//    ghost_shell_faces -- this method does parallel communication to
//        find all faces that attach to shared nodes and commuicate
//        their connectivities to all processors that share that node.
//        This essentially ghosts the connectivities for all faces 
//        connected to shared nodes to all sharing processors. The
//        individual processors can then construct their shell toplogies
//        independently and still have the std::right behavior at shared nodes.
//
//------------------------------------------------------------------
  int i,j,k;

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":----------> ghost_shell_faces <-----------" << "\n";
  postream.flush();
#endif

  // we should only be here if we have more than 1 processor  
  PRECONDITION(contact_number_of_processors(SearchComm) > 1);

  // Parallel case -- the "full" versions of the conn and inv_conn
  //    data structures will include both the local faces and
  //    any faces that need to be ghosted in from other processors.
  
  bool* send_face_flag = new bool[num_orig_shell_faces];
  for (i = 0; i < num_orig_shell_faces; ++i) send_face_flag[i]=false; 
  int* num_send_faces = new int[number_comm_partners];
  int** send_face_list = new int*[number_comm_partners];
  for (i = 0; i < number_comm_partners; ++i) 
    num_send_faces[i] = 0;
  
  //--------------------------------------------------
  // compute number faces connected to shared nodes and store them in a 
  // temporary list
  int node_idx = 0;
  int tot_num_faces_to_send = 0;
  for (i = 0; i < number_comm_partners; ++i) {
#ifdef SHELL_DEBUG
	postream << my_proc_id << ":    -> for proc " 
	     << comm_proc_ids[i] << " num shared nodes " 
	     << number_nodes_to_partner[i]
	     << "\n";
#endif
    for (j = 0; j < number_nodes_to_partner[i]; ++j) {
      int node = comm_nodes[node_idx++]-1; // -1 for fortran
      //postream << my_proc_id << ":      -> comm node " << node << "\n"; 
      for ( k = 0; k < inv_conn_count[node]; ++k) {
	int face = inv_conn[offset_inv_conn[node]+k];
	send_face_flag[face] = true;
#ifdef SHELL_DEBUG
	postream << my_proc_id << ":        -> send " << face << " to proc " 
	     << comm_proc_ids[i] << "\n";
#endif
      }
    }
    for (j = 0; j < num_orig_shell_faces; ++j) 
      if ( send_face_flag[j] ) num_send_faces[i]++;
    tot_num_faces_to_send += num_send_faces[i];
    send_face_list[i] = new int[num_send_faces[i]];
    int face_idx = 0;
    for (j = 0; j < num_orig_shell_faces; ++j) {
      if ( send_face_flag[j] ){
	send_face_flag[j]=false; 
	send_face_list[i][face_idx++] = j;
      }
    }
  }
  
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":   tot num faces to send:" 
       << tot_num_faces_to_send << "\n";
#endif


  //--------------------------------------------------
  // communicate number of faces connected to shared nodes to sharing
  // processors

  RequestHandle * recv_handles = new RequestHandle[number_comm_partners];
  char * recv_num_face_buff = new char[number_comm_partners*sizeof(int)];
  char * send_num_face_buff = new char[number_comm_partners*sizeof(int)];
  int *recv_num_faces, *send_num_faces;
  recv_num_faces = reinterpret_cast<int *>(recv_num_face_buff);
  send_num_faces = reinterpret_cast<int *>(send_num_face_buff);
  for (i=0; i<number_comm_partners; ++i) {
    send_num_faces[i] = num_send_faces[i];
#ifdef SHELL_DEBUG
    postream << my_proc_id << ":       + faces to send to" 
	 << comm_proc_ids[i] << " : " << send_num_faces[i] << "\n"; 
#endif
  }

  int * one_for_all = new int[number_comm_partners];
  for ( i = 0; i < number_comm_partners; ++i) one_for_all[i] = 1;

  ContactAsymComm num_ghost_comm(number_comm_partners,
				 number_comm_partners,
				 one_for_all,
				 one_for_all,
				 comm_proc_ids,
				 comm_proc_ids);

  delete [] one_for_all;
  
  contact_communicate_packed_buffers(SearchComm,
				     num_ghost_comm,
				     sizeof(int),
				     send_num_face_buff,
				     recv_num_face_buff,
				     recv_handles);

  // count total number of faces to recv
  int tot_num_faces_to_recv = 0;
  for (i=0; i<number_comm_partners; ++i) {
    tot_num_faces_to_recv += recv_num_faces[i];
#ifdef SHELL_DEBUG
    postream << my_proc_id << ":       + faces to recv from" 
	 << comm_proc_ids[i] << " : " << recv_num_faces[i] << "\n"; 
#endif
  }
  
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":   tot num faces to recv:" 
       << tot_num_faces_to_recv << "\n";
#endif

  //--------------------------------------------------
  // communicate faces connected to shared nodes to sharing
  // processors
  // ASG HACK -- assuming max of 4 nodes/face
  // real comm to push data
  // arrangement of comm data for each face:
  //   Int data:
  //      0: face blk #
  //      1: face blk id
  //      2-9: host code ids of face connectivity
  //      10-11: global id of face
  //      12: 1 if a shell face, 0 otherwise
  //   Real data (following int data in memory)
  //      0-11: coords of nodes in connectivity
  //      12:   shell face lofting factor
  int size_int_comm_per_face = 13;
  int size_real_comm_per_face = 13;
  int face_buff_size = size_int_comm_per_face*sizeof(int) 
    + size_real_comm_per_face*sizeof(double);

  int * send_face_ints = 
    new int[tot_num_faces_to_send*size_int_comm_per_face];
  double * send_face_reals = 
    new double[tot_num_faces_to_send*size_real_comm_per_face];
  int * recv_face_ints = 
    new int[tot_num_faces_to_recv*size_int_comm_per_face];
  double * recv_face_reals = 
    new double[tot_num_faces_to_recv*size_real_comm_per_face];

  RequestHandle * recv_handles2;
  char * recv_face_buff, * send_face_buff;
  topology->comm_buffer->Buffers( tot_num_faces_to_send*face_buff_size,
				  tot_num_faces_to_recv*face_buff_size,
				  number_comm_partners,
				  &send_face_buff, &recv_face_buff, 
				  &recv_handles2 );

  // store face data to send
  int send_faces_int_idx = 0;
  int send_faces_real_idx = 0;
  int buff_mem_addr =0;
  int face_int_count = 0;
  int face_real_count = 0;
  for( i=0 ; i<number_comm_partners ; ++i){
    for( j=0 ; j<send_num_faces[i] ; ++j){
      int face = send_face_list[i][j];
      //get blk and id of face
      int blk = shell_face_blk_list[face];
      int id = shell_face_id_list[face];
      int num_face_nodes = 
	search->Number_Nodes_Per_Face(face_block_types[blk]);
      POSTCONDITION(num_face_nodes>0);
      send_face_ints[send_faces_int_idx++] = blk;
      send_face_ints[send_faces_int_idx++] = id;
      // get connectivity of face as node global ids
      for ( k = 0; k < num_face_nodes; ++k) {
	int local_node = 
	  face_connectivity[conn_offsets[blk]+id*num_face_nodes+k]-1;
	send_face_ints[send_faces_int_idx + 2*k] = 
	  orig_global_node_ids[2*local_node];
	send_face_ints[send_faces_int_idx + 2*k+1] = 
	  orig_global_node_ids[2*local_node+1];
	send_face_reals[send_faces_real_idx + 3*k] = coords[3*local_node];
	send_face_reals[send_faces_real_idx + 3*k+1] = coords[3*local_node+1];
	send_face_reals[send_faces_real_idx + 3*k+2] = coords[3*local_node+2];
      }
      send_faces_int_idx += 8;
      send_faces_real_idx += 12;
      // get lofting factor
      send_face_reals[send_faces_real_idx++] = shell_loft_factors[face];
      // get face global id
      int face_num = 0;
      for ( k = 0; k < blk; ++k) 
	face_num += num_orig_faces_per_blk[k];
      face_num += id;
      send_face_ints[send_faces_int_idx++] = global_face_ids[2*face_num];
      send_face_ints[send_faces_int_idx++] = global_face_ids[2*face_num+1];
      if (ContactSearch::Is_a_Shell_Face(face_block_types[blk]))
	send_face_ints[send_faces_int_idx++] = 1;
      else
	send_face_ints[send_faces_int_idx++] = 0;
      
#ifdef SHELL_DEBUG
      postream << my_proc_id << ":        -> sending ints to proc " 
	   << comm_proc_ids[i] << ": B"
	   << send_face_ints[send_faces_int_idx - 13] << " ID"
	   << send_face_ints[send_faces_int_idx - 12] << " | "
	   << send_face_ints[send_faces_int_idx - 11] << " "
	   << send_face_ints[send_faces_int_idx - 10] << "/" 
	   << send_face_ints[send_faces_int_idx - 9] << " "
	   << send_face_ints[send_faces_int_idx - 8] << "/"
	   << send_face_ints[send_faces_int_idx - 7] << " "
	   << send_face_ints[send_faces_int_idx - 6] << "/"
	   << send_face_ints[send_faces_int_idx - 5] << " "
	   << send_face_ints[send_faces_int_idx - 4] << " fid: "
	   << send_face_ints[send_faces_int_idx - 3] << " "
	   << send_face_ints[send_faces_int_idx - 2] << " shell?: "
	   << send_face_ints[send_faces_int_idx - 1] << " "
	   << "\n";
      postream << my_proc_id << ":        -> sending reals to proc " 
	   << comm_proc_ids[i] << ": 1) "
	   << send_face_reals[send_faces_real_idx - 13] << ","
	   << send_face_reals[send_faces_real_idx - 12] << ","
	   << send_face_reals[send_faces_real_idx - 11] << " 2) "
	   << send_face_reals[send_faces_real_idx - 10] << "," 
	   << send_face_reals[send_faces_real_idx - 9] << ","
	   << send_face_reals[send_faces_real_idx - 8] << " 3) "
	   << send_face_reals[send_faces_real_idx - 7] << ","
	   << send_face_reals[send_faces_real_idx - 6] << ","
	   << send_face_reals[send_faces_real_idx - 5] << " 4) "
	   << send_face_reals[send_faces_real_idx - 4] << ","
	   << send_face_reals[send_faces_real_idx - 3] << ","
           << send_face_reals[send_faces_real_idx - 2] << "; LF "
	   << send_face_reals[send_faces_real_idx - 1]
	   << "\n";
#endif
    }
    std::memcpy(send_face_buff + buff_mem_addr, &send_face_ints[face_int_count],
	   send_num_faces[i]*sizeof(int)*size_int_comm_per_face);
    buff_mem_addr += send_num_faces[i]*sizeof(int)*size_int_comm_per_face;
    face_int_count += send_num_faces[i]*size_int_comm_per_face;
    std::memcpy(send_face_buff + buff_mem_addr, &send_face_reals[face_real_count],
	   send_num_faces[i]*sizeof(double)*size_real_comm_per_face);
    buff_mem_addr += send_num_faces[i]*sizeof(double)*size_real_comm_per_face;
    face_real_count += send_num_faces[i]*size_real_comm_per_face;
  }
  
#ifdef SHELL_DEBUG
  {
    postream << my_proc_id << ": ***** here is send_face_buff" << "\n";
    PRECONDITION(size_int_comm_per_face == 13);
    PRECONDITION(size_real_comm_per_face == 13);
    int buff_off = 0;
    for( i=0 ; i<number_comm_partners ; ++i){
      int * send_face_int = reinterpret_cast<int *>(send_face_buff + buff_off);
      buff_off += send_num_faces[i]*size_int_comm_per_face*sizeof(int);
      double * send_face_real = 
	reinterpret_cast<double *>(send_face_buff + buff_off);
      buff_off += send_num_faces[i]*size_real_comm_per_face*sizeof(double);
      int send_int_idx = 0;
      int send_real_idx = 0;
      for( j=0 ; j<send_num_faces[i] ; ++j){
	send_int_idx += size_int_comm_per_face;
	postream << my_proc_id << ":        -> sending int to proc " 
	     << comm_proc_ids[i] << ": B"
	     << send_face_int[send_int_idx - 13] << " ID"
	     << send_face_int[send_int_idx - 12] << " | "
	     << send_face_int[send_int_idx - 11] << " "
	     << send_face_int[send_int_idx - 10] << "/" 
	     << send_face_int[send_int_idx - 9] << " "
	     << send_face_int[send_int_idx - 8] << "/"
	     << send_face_int[send_int_idx - 7] << " "
	     << send_face_int[send_int_idx - 6] << "/"
	     << send_face_int[send_int_idx - 5] << " "
	     << send_face_int[send_int_idx - 4] << " fid: "
	     << send_face_int[send_int_idx - 3] << " "
	     << send_face_int[send_int_idx - 2] << " shell?: "
	     << send_face_int[send_int_idx - 1] << " "
	     << "\n";
	send_real_idx += size_real_comm_per_face;
	postream << my_proc_id << ":        -> sending real to proc " 
	     << comm_proc_ids[i] << ": 1) "
	     << send_face_real[send_real_idx - 13] << ","
	     << send_face_real[send_real_idx - 12] << ","
	     << send_face_real[send_real_idx - 11] << " 2) "
	     << send_face_real[send_real_idx - 10]  << "," 
	     << send_face_real[send_real_idx - 9]  << ","
	     << send_face_real[send_real_idx - 8]  << " 3) "
	     << send_face_real[send_real_idx - 7]  << ","
	     << send_face_real[send_real_idx - 6]  << ","
	     << send_face_real[send_real_idx - 5]  << " 4) "
	     << send_face_real[send_real_idx - 4]  << ","
	     << send_face_real[send_real_idx - 3]  << ","
       	     << send_face_real[send_real_idx - 2]  << ":LF " 
	     << send_face_real[send_real_idx - 1] 
	     << "\n";
      }
    }
  }
#endif

  // do communication
  ContactAsymComm ghost_comm(number_comm_partners,
			     number_comm_partners,
			     send_num_faces,
			     recv_num_faces,
			     comm_proc_ids,
			     comm_proc_ids);

  contact_communicate_packed_buffers(SearchComm,
				     ghost_comm,
				     face_buff_size,
				     send_face_buff,
				     recv_face_buff,
				     recv_handles2);
#ifdef SHELL_DEBUG
  {
    postream << my_proc_id << ": ***** here is recv_face_buff" << "\n";
    PRECONDITION(size_int_comm_per_face == 13);
    PRECONDITION(size_real_comm_per_face == 13);
    int buff_off = 0;
    for( i=0 ; i<number_comm_partners ; ++i){
      int * recv_face_int = reinterpret_cast<int *>(recv_face_buff + buff_off);
      buff_off += recv_num_faces[i]*size_int_comm_per_face*sizeof(int);
      double * recv_face_real = 
	reinterpret_cast<double *>(recv_face_buff + buff_off);
      buff_off += recv_num_faces[i]*size_real_comm_per_face*sizeof(double);
      int recv_int_idx = 0;
      int recv_real_idx = 0;
      for( j=0 ; j<recv_num_faces[i] ; ++j){
	recv_int_idx += size_int_comm_per_face;
	postream << my_proc_id << ":        -> recv int from proc " 
	     << comm_proc_ids[i] << ": B"
	     << recv_face_int[recv_int_idx - 13] << " ID"
	     << recv_face_int[recv_int_idx - 12] << " | "
	     << recv_face_int[recv_int_idx - 11] << " "
	     << recv_face_int[recv_int_idx - 10] << "/" 
	     << recv_face_int[recv_int_idx - 9] << " "
	     << recv_face_int[recv_int_idx - 8] << "/"
	     << recv_face_int[recv_int_idx - 7] << " "
	     << recv_face_int[recv_int_idx - 6] << "/"
	     << recv_face_int[recv_int_idx - 5] << " "
	     << recv_face_int[recv_int_idx - 4] << " fid: "
	     << recv_face_int[recv_int_idx - 3] << " "
	     << recv_face_int[recv_int_idx - 2] << " shell?: "
	     << recv_face_int[recv_int_idx - 1] << " "
	     << "\n";
	recv_real_idx += size_real_comm_per_face;
	postream << my_proc_id << ":        -> recv real from proc " 
	     << comm_proc_ids[i] << ": 1) "
	     << recv_face_real[recv_real_idx - 13] << ","
	     << recv_face_real[recv_real_idx - 12] << ","
	     << recv_face_real[recv_real_idx - 11] << " 2) "
	     << recv_face_real[recv_real_idx - 10]  << "," 
	     << recv_face_real[recv_real_idx - 9]  << ","
	     << recv_face_real[recv_real_idx - 8]  << " 3) "
	     << recv_face_real[recv_real_idx - 7]  << ","
	     << recv_face_real[recv_real_idx - 6]  << ","
	     << recv_face_real[recv_real_idx - 5]  << " 4) "
	     << recv_face_real[recv_real_idx - 4]  << ","
	     << recv_face_real[recv_real_idx - 3]  << ","
	     << recv_face_real[recv_real_idx - 2]  << ":LF " 
	     << recv_face_real[recv_real_idx - 1] 
	     << "\n";
      }
    }
  }
#endif


  buff_mem_addr = 0;
  face_int_count = 0;
  face_real_count = 0;
  for( i=0 ; i<number_comm_partners ; ++i){
    std::memcpy(&recv_face_ints[face_int_count] , recv_face_buff + buff_mem_addr,
	   recv_num_faces[i]*sizeof(int)*size_int_comm_per_face);
    buff_mem_addr += recv_num_faces[i]*sizeof(int)*size_int_comm_per_face;
    face_int_count += recv_num_faces[i]*size_int_comm_per_face;
    std::memcpy(&recv_face_reals[face_real_count] , recv_face_buff + buff_mem_addr,
	   recv_num_faces[i]*sizeof(double)*size_real_comm_per_face);
    buff_mem_addr += recv_num_faces[i]*sizeof(double)*size_real_comm_per_face;
    face_real_count += recv_num_faces[i]*size_real_comm_per_face;
  }

  //--------------------------------------------------
  // all the communications have completed -- we now have a complete set
  // of face connectivities for all shell faces that are adjoined to
  // the shared shell nodes. Now create the "full" conn and inv_conn
  // data structures; i.e. including both local and ghosted faces
  //--------------------------------------------------

  //--------------------------------------------------
  // find out number of tri and quad faces that were sent
  // and allocate data structures for constructing full connectivity 
  // ASG HACK -- this assumes only 3-node and 4-node shells
  int num_tri(0), num_tri_shell(0), num_quad(0), num_quad_shell(0),
    num_ghost_face_nodes(0);
  for ( i = 0; i < tot_num_faces_to_recv; ++i) {
    int blk_id = recv_face_ints[i*size_int_comm_per_face];
    POSTCONDITION(search->Number_Nodes_Per_Face(face_block_types[blk_id])>0);
    if ( face_block_types[blk_id] == ContactSearch::SHELLQUADFACEL4 ) {
      num_quad_shell++;
      num_ghost_face_nodes +=4;
    }
    else if ( face_block_types[blk_id] == ContactSearch::SHELLTRIFACEL3 ) {
      num_tri_shell++;
      num_ghost_face_nodes +=3;
    }
    else if ( face_block_types[blk_id] == ContactSearch::QUADFACEL4 ) {    
      num_quad++;
      num_ghost_face_nodes +=4;
    }
    else if ( face_block_types[blk_id] == ContactSearch::TRIFACEL3 ) {     
      num_tri++;
      num_ghost_face_nodes +=3;
    }
    else {
      topology->errors->Add_Error_Message( "Unrecognized face" );
      return ContactSearch::INVALID_DATA;
    }
  }    
  
  // count number of original faces
  int num_orig_faces = 0;
  int num_orig_face_nodes = 0;
  for ( i = 0; i < num_orig_face_blks; ++i) {
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[i]);
    POSTCONDITION(num_face_nodes>0);
    num_orig_faces +=  num_orig_faces_per_blk[i];
    num_orig_face_nodes +=  num_orig_faces_per_blk[i]*num_face_nodes;
  }

  // allocate new conn and id arrays
  full_face_connectivity = 
    new int[num_orig_face_nodes + num_ghost_face_nodes];
  full_conn_offsets = new int[num_orig_face_blks+4];
  full_face_block_types = 
    new ContactSearch::ContactFace_Type[num_orig_face_blks+4];
  full_num_faces_per_blk = new int[num_orig_face_blks+4];
  int * full_global_node_ids = 
    new int[(number_host_code_nodes + num_ghost_face_nodes)*2];
  full_global_face_ids = 
    new int[(num_orig_faces + tot_num_faces_to_recv)*2];
  full_coords = 
    new double[(number_host_code_nodes + num_ghost_face_nodes)*3];
  full_shell_loft_factors = 
    new double [num_orig_faces + tot_num_faces_to_recv];

  //--------------------------------------------------
  // fill full connectivity data structures
  //  -- copy old node ids and coords into full node ids array
  full_num_host_code_nodes = number_host_code_nodes;
  for ( i = 0; i < number_host_code_nodes; ++i) {
    full_global_node_ids[2*i] = orig_global_node_ids[2*i];
    full_global_node_ids[2*i+1] = orig_global_node_ids[2*i+1];
    full_coords[3*i] = coords[3*i];
    full_coords[3*i+1] = coords[3*i+1];
    full_coords[3*i+2] = coords[3*i+2];
  }

  //  -- copy local connectivity into full structure
  int orig_fidx(0);
  int cidx(0);
  for ( i = 0; i < num_orig_face_blks; ++i) {
    
      full_conn_offsets[i] = conn_offsets[i];
    full_face_block_types[i] = face_block_types[i];
    full_num_faces_per_blk[i] = num_orig_faces_per_blk[i];
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[i]);
    POSTCONDITION(num_face_nodes>0);
    for ( j = 0; j < num_orig_faces_per_blk[i]; ++j) {
      orig_fidx++;
      for ( k = 0; k < num_face_nodes; ++k) {
	full_face_connectivity[cidx] = face_connectivity[cidx];
	cidx++;
      }
    }
  }

  //  -- copy global face ids for local faces considered as shell faces
  for ( i = 0; i < num_orig_shell_faces; ++i) {
    int blk = shell_face_blk_list[i];
    int id = shell_face_id_list[i];
    int idx = 0;
    for (j = 0; j < blk; ++j) idx += num_orig_faces_per_blk[j];
    idx += id;
    full_global_face_ids[2*i] = global_face_ids[2*idx];
    full_global_face_ids[2*i+1] = global_face_ids[2*idx+1];
  }

  //  -- copy local loft factors into full array
  for ( i =0; i < num_orig_shell_faces; ++i)
    full_shell_loft_factors[i] = shell_loft_factors[i];

  // -- add ghosted faces
  int blk = num_orig_face_blks;
  full_face_block_types[blk] = ContactSearch::SHELLQUADFACEL4;
  full_conn_offsets[blk] = num_orig_face_nodes;
  full_num_faces_per_blk[blk] = num_quad_shell;
  blk++;
  full_face_block_types[blk] = ContactSearch::SHELLTRIFACEL3;
  full_conn_offsets[blk] = full_conn_offsets[blk-1] + num_quad_shell*4;
  full_num_faces_per_blk[blk] = num_tri_shell;
  blk++;
  full_face_block_types[blk] = ContactSearch::QUADFACEL4;
  full_conn_offsets[blk] = full_conn_offsets[blk-1] + num_tri_shell*3;
  full_num_faces_per_blk[blk] = num_quad;
  blk++;
  full_face_block_types[blk] = ContactSearch::TRIFACEL3;
  full_conn_offsets[blk] = full_conn_offsets[blk-1] + num_quad*4;
  full_num_faces_per_blk[blk] = num_tri;

  int connoff[4];
  connoff[0] = cidx;
  connoff[1] = connoff[0] + num_quad_shell*4;
  connoff[2] = connoff[1] + num_tri_shell*3;
  connoff[3] = connoff[2] + num_quad*4;
  int fidx[4];
  fidx[0] = num_orig_shell_faces;
  fidx[1] = fidx[0] + num_quad_shell;
  fidx[2] = fidx[1] + num_tri_shell;
  fidx[3] = fidx[2] + num_quad;

  for ( i = 0; i < tot_num_faces_to_recv; ++i) {
    int bidx = -1;
    int blk_id = recv_face_ints[i*size_int_comm_per_face];
    int num_face_nodes = 
      search->Number_Nodes_Per_Face(face_block_types[blk_id]);
    POSTCONDITION(num_face_nodes>0);
    if ( face_block_types[blk_id] == ContactSearch::SHELLQUADFACEL4 ) bidx = 0;
    if ( face_block_types[blk_id] == ContactSearch::SHELLTRIFACEL3 ) bidx = 1;
    if ( face_block_types[blk_id] == ContactSearch::QUADFACEL4 ) bidx = 2;
    if ( face_block_types[blk_id] == ContactSearch::TRIFACEL3 ) bidx = 3;
    POSTCONDITION(bidx>=0);
    for ( j = 0; j < num_face_nodes; ++j) {
      int node_id[2];
      node_id[0] = recv_face_ints[i*size_int_comm_per_face+2+2*j];
      node_id[1] = recv_face_ints[i*size_int_comm_per_face+2+2*j+1];
      // find node id in processor local numbering
      bool found_node = false;
      int found_node_idx(-1);
      for ( k = 0; k < full_num_host_code_nodes; ++k) {
	if (node_id[0] == full_global_node_ids[k*2] &&
	    node_id[1] == full_global_node_ids[k*2+1]) {
	  found_node_idx = k;
	  found_node = true;
	  break;
	}
      }
      if ( !found_node ) {
	// node is not in local numbering. Add it.
	full_global_node_ids[full_num_host_code_nodes*2] = node_id[0];
	full_global_node_ids[full_num_host_code_nodes*2+1] = node_id[1];
	full_coords[full_num_host_code_nodes*3] = 
	  recv_face_reals[i*size_real_comm_per_face+3*j];
	full_coords[full_num_host_code_nodes*3+1] = 
	  recv_face_reals[i*size_real_comm_per_face+3*j+1];
	full_coords[full_num_host_code_nodes*3+2] = 
	  recv_face_reals[i*size_real_comm_per_face+3*j+2];
	found_node_idx = full_num_host_code_nodes;
	full_num_host_code_nodes ++;
      }
      POSTCONDITION( found_node_idx >= 0);
      full_face_connectivity[connoff[bidx]++] = found_node_idx+1;// +1 for fortran
    }
    // set face id
    full_global_face_ids[2*fidx[bidx]] = 
      recv_face_ints[i*size_int_comm_per_face+10];
    full_global_face_ids[2*fidx[bidx]+1] = 
      recv_face_ints[i*size_int_comm_per_face+11];
    // set loft factor
    full_shell_loft_factors[fidx[bidx]] = 
      recv_face_reals[i*size_real_comm_per_face + 12];
    fidx[bidx]++;
  }
  
  postream.flush();
  POSTCONDITION ( fidx[3] == num_orig_shell_faces + tot_num_faces_to_recv);
  POSTCONDITION ( connoff[3] == num_orig_face_nodes + num_ghost_face_nodes);

  full_number_face_blocks = num_orig_face_blks+4;
  
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> *** The full global ids:" << "\n";
  for ( i = 0; i < full_num_host_code_nodes; ++i){
    postream << my_proc_id << ":TEST>      node: " << i << " ID: "
	 << full_global_node_ids[2*i] << "," 
	 << full_global_node_ids[2*i+1] << "\n";
  }
#endif

  //----------------------
  // clean up memory
  delete [] num_send_faces;
  delete [] send_face_flag;
  delete [] recv_num_face_buff;
  delete [] send_num_face_buff;
  delete [] recv_handles;
  for ( i = 0; i < number_comm_partners; ++i)
    delete [] send_face_list[i];
  delete [] send_face_list;
  delete [] recv_face_ints;
  delete [] send_face_ints;
  delete [] recv_face_reals;
  delete [] send_face_reals;
  delete [] full_global_node_ids;
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":<------ std::exit ghost_shell_faces ------->" << "\n";
  postream.flush();
#endif

  return ContactSearch::NO_ERROR;
};

//------------------------------------------------------------------
//
//    parallel_node_resolution -- this method takes the shell topology
//        that has been completed in serial by all processors and
//        does the final communcation steps required to finalize all the
//        created nodes and their IDs.
//
//------------------------------------------------------------------
ContactSearch::ContactErrorCode
ContactShellHandler::
parallel_node_resolution(const int* full_number_faces_per_block,
			 const int* full_updated_inv_conn,
			 const int* full_updated_offset_inv_conn,
			 const int* full_updated_inv_conn_count,
			 const int* full_updated_shell_face_blk,
			 const int* full_updated_shell_face_id,
			 const int* full_face_global_ids,
			 const int& number_comm_partners,
			 const int* comm_proc_ids,
			 const int* host_number_nodes_to_partner,
			 const int* host_comm_nodes,
			 int* number_nodes_to_partner,
			 int* comm_nodes,
			 int* exodus_node_ids,
			 int* global_node_ids,
			 ShellID_Resolution res_method,
			 ContactShellHandler * old_handler,
			 ContactSearch* search) 
{

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>=================== Finallizing parallel node data" << "\n";
  if (res_method == USE_OLD_HANDLER) 
    postream << my_proc_id << ":TEST>     *** use old handler for id seed \n";
  postream.flush();
#endif

  int i,j, k, l, m;

  // to do face-id based comparison of shared nodes, we need an inverse
  // connectivity based on the modified full face connectivity. 
  
//   int * full_inv_conn(NULL);
//   int * full_offset_inv_conn(NULL);
//   int * full_inv_conn_count(NULL);
//   int * full_shell_face_blk(NULL);
//   int * full_shell_face_id(NULL);
//   int full_num_shell_faces;
//   int full_num_shell_nodes;
  
//   compute_inv_conn(number_acme_nodes,
// 		   full_number_face_blocks,
// 		   full_face_connectivity,
// 		   full_conn_offsets,
// 		   full_face_block_types,
// 		   full_number_faces_per_block,
// 		   0,
// 		   NULL,
// 		   NULL,
// 		   NULL,
// 		   full_inv_conn,
// 		   full_offset_inv_conn,
// 		   full_inv_conn_count,
// 		   full_shell_face_blk,
// 		   full_shell_face_id,
// 		   full_num_shell_faces,
// 		   full_num_shell_nodes,
// 		   false,
// 		   search);

//   POSTCONDITION(full_inv_conn);
//   POSTCONDITION(full_offset_inv_conn);
//   POSTCONDITION(full_inv_conn_count );
//   POSTCONDITION(full_shell_face_blk );
//   POSTCONDITION(full_shell_face_id  );
// #ifdef SHELL_DEBUG
//   {
//     postream << my_proc_id << ":HEY!!!! Here is inv_conn in par_res" << "\n";
//     for ( i = 0; i < number_acme_nodes ; i ++) {
//       char ook[80];
//       std::sprintf(ook,"%d:   Node %d :",my_proc_id,i);
//       int off = full_offset_inv_conn[i];
//       for ( j = 0; j < full_inv_conn_count[i]; j++){
// 	char fishy[10];
// 	std::sprintf(fishy," %d",full_inv_conn[off+j]);
// 	std::strcat(ook, fishy);
//       }
//       postream << ook << "\n";
//     }
//     postream << my_proc_id << ":HEY!!!! Done with inv_conn in par_res" << "\n";
//   }
// #endif
  
  //--------------------------------------------  
  // To compute the ids of the nodes, we must first count up how many nodes
  // that each processor owns.  This includes all nodes which are local to 
  // a processor and are not shared, as well as nodes that are shared but
  // for which this processor has the lowest proc id of those which share it.
  
  //   we want total numbering base, and we want the local numbering base
  //   the total numbering base is the total number of host code nodes
  //   the local numbering is the lower processor sum on the number of 
  //   nodes that that processor created and which it owns.

  NodeStatus * node_status = new NodeStatus[number_host_code_nodes];
  int * node_owner = new int[number_host_code_nodes];
  for ( i = 0; i < number_host_code_nodes; ++i) {
    node_status[i] = OWNED_NOT_SHARED;
    node_owner[i] = my_proc_id;
  }
  
  int num_owned_newly_created_nodes =
    number_acme_nodes - number_host_code_nodes;
  int num_owned_host_code_nodes = number_host_code_nodes;
  int idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    int proc_id = comm_proc_ids[i];
    for (j = 0; j < host_number_nodes_to_partner[i]; ++j){
      int node = host_comm_nodes[idx++]-1; // -1 for fortran
      if ( node_status[node] == OWNED_NOT_SHARED || 
	   node_status[node] == OWNED_SHARED) {
	if (proc_id < my_proc_id) {
	  node_status[node] = NOT_OWNED_SHARED;
	  num_owned_host_code_nodes --;
	  if (num_created_shell_nodes[node] > 0 )
	    num_owned_newly_created_nodes -= 
	      Num_Acme_Nodes_for_Host_Node(node) - 1; // only counting
	                                              // created nodes
	                                              // more than original
	  if (node_owner[node] > proc_id) node_owner[node] = proc_id;
	}
	else {
	  node_status[node] = OWNED_SHARED;
	}
      }
    }
  }

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> number_host_code_nodes: " 
       << number_host_code_nodes
       << "  num_owned_host_code_nodes: " 
       << num_owned_host_code_nodes << "\n";
  postream << my_proc_id << ":TEST> number_acme_nodes: " << number_acme_nodes
       << "  num_owned_newly_created_nodes: " 
       << num_owned_newly_created_nodes << "\n";
  postream << my_proc_id << ":TEST> node status:" << "\n";
  for ( i =0 ; i < number_host_code_nodes; ++i) {
    postream << my_proc_id << ":TEST>    node:" << i << " status:" 
	 << node_status[i] << "\n";
  }
#endif

  if (res_method == SEQUENTIAL_IDS) 
    number_nodes_sequentially( num_owned_host_code_nodes,
			       num_owned_newly_created_nodes,
			       node_status,
			       exodus_node_ids, 
			       global_node_ids);
  else
    number_nodes_using_handler( num_owned_host_code_nodes,
				num_owned_newly_created_nodes,
				node_status,
				full_updated_inv_conn,
				full_updated_inv_conn_count,
				full_updated_offset_inv_conn,
				full_updated_shell_face_blk,
				full_updated_shell_face_id,
				full_number_faces_per_block,
				full_face_global_ids,
				old_handler,
				exodus_node_ids, 
				global_node_ids);
#ifdef SHELL_DEBUG
  postream << ":Done with new total local numbering:\n";
#endif

  //--------------------------------------------  
  // At this point, all nodes have been properly numbered on the owning 
  // processors. Now the owning processors need to communicate their
  // determined node numberings to the shared nodes so they can all be
  // syncronized.
  //
  // if we are creating the original shell mesh, then we need only
  // communicate data for the shell nodes. If we are updating the 
  // shell mesh, then we may have had to renumber non-shell nodes
  // as well, so we have to communicate all shared nodes.
  
  // each processor already knows who is the owner for the shared nodes and
  // now much data will be set. 

  // compute total number to send and recv

  int tot_num_nodes_to_recv = 0;
  int tot_num_host_nodes_to_recv = 0;
  int tot_num_nodes_to_send = 0;
  int * num_nodes_to_recv = new int[number_comm_partners];
  int * num_nodes_to_send = new int[number_comm_partners];
  for (i = 0; i < number_comm_partners; ++i) num_nodes_to_recv[i] = 0;
  for (i = 0; i < number_comm_partners; ++i) num_nodes_to_send[i] = 0;
    
  idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    int proc_id = comm_proc_ids[i];
    for (j = 0; j < host_number_nodes_to_partner[i]; ++j){
      int node = host_comm_nodes[idx++]-1; // -1 for fortran
      if ( num_created_shell_nodes[node] <= 0 &&
	   res_method == SEQUENTIAL_IDS) 
	continue; // skip non-shell nodes when numbering sequentially
      if (node_status[node] == OWNED_SHARED) {
	tot_num_nodes_to_send += Num_Acme_Nodes_for_Host_Node(node);
	num_nodes_to_send[i] += Num_Acme_Nodes_for_Host_Node(node);
      }
      else {
	if (node_owner[node] == proc_id) {
	  tot_num_host_nodes_to_recv ++;
	  tot_num_nodes_to_recv += Num_Acme_Nodes_for_Host_Node(node);
	  num_nodes_to_recv[i] += Num_Acme_Nodes_for_Host_Node(node);
	}
      }
    }
  }
  
  
  //--------------------------------------------------
  // communicate ids of shared nodes to sharing processors
  // node info to send:
  //   0-1: host code id
  //   2: exodus id
  //   3-4: global id of a face attached to node
  int node_info_size = 5;
  RequestHandle * recv_handles;
  char * recv_node_buff, * send_node_buff;
  topology->comm_buffer->
    Buffers( tot_num_nodes_to_send*node_info_size*sizeof(int),
	     tot_num_nodes_to_recv*node_info_size*sizeof(int),
	     number_comm_partners,
	     &send_node_buff, &recv_node_buff, 
	     &recv_handles );
  int *recv_nodes, *send_nodes;
  recv_nodes = reinterpret_cast<int *>(recv_node_buff);
  send_nodes = reinterpret_cast<int *>(send_node_buff);
  for (i=0;i<tot_num_nodes_to_send*node_info_size; ++i) send_nodes[i] = 0;

  // compute nodes to send
#ifdef SHELL_DEBUG
  postream << "Computing Nodes to Send\n";
#endif
  idx = 0;
  int node_idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    for (j = 0; j < host_number_nodes_to_partner[i]; ++j){
      int node = host_comm_nodes[idx++]-1 ; // -1 for fortran
      if ( num_created_shell_nodes[node] <= 0 &&
	   res_method == SEQUENTIAL_IDS) 
	continue; // skip non-shell nodes when numbering sequentially
      if (node_status[node] == OWNED_SHARED) {
	for (k = 0; k < Num_Acme_Nodes_for_Host_Node(node); ++k) {
	  int created_node = Acme_Node_for_Host_Node( node, k );
#ifdef SHELL_DEBUG
	  postream << "Created node = " << created_node << "\n";
#endif
	  send_nodes[node_idx*node_info_size] = 
	    global_node_ids[2*created_node];
	  send_nodes[node_idx*node_info_size + 1] = 
	    global_node_ids[2*created_node + 1];
	  send_nodes[node_idx*node_info_size + 2] = 
	    exodus_node_ids[created_node];
	  // get global id of an attached face if the node is a shell node.
	  // if it is a non-shell node, then we don't need to send a
	  // face id
	  if ( num_created_shell_nodes[node] > 0 ) { // shell node
	    PRECONDITION(full_updated_inv_conn_count[created_node] > 0);
	    int face = full_updated_inv_conn[full_updated_offset_inv_conn[created_node]];
	    send_nodes[node_idx*node_info_size + 3] = 
	      full_face_global_ids[2*face];
	    send_nodes[node_idx*node_info_size + 4] = 
	      full_face_global_ids[2*face+1];
	  } else { // non-shell node
	    send_nodes[node_idx*node_info_size + 3] = -1;
	    send_nodes[node_idx*node_info_size + 4] = -1;
	  }
	  node_idx++;
#ifdef SHELL_DEBUG
	int proc_id = comm_proc_ids[i];
	postream << my_proc_id << ":       for node " << node
	     << " created node #:" << k 
	     << " send data to " << proc_id << ": gid " 
	     << send_nodes[node_idx*node_info_size -5] << " " 
	     << send_nodes[node_idx*node_info_size -4] << " xid " 
	     << send_nodes[node_idx*node_info_size -3] << " face " 
	     << send_nodes[node_idx*node_info_size -2] << " " 
	     << send_nodes[node_idx*node_info_size -1] << "\n"; 
#endif
	}
      }
    }
  }

  // do communication
  ContactAsymComm node_comm(number_comm_partners,
			    number_comm_partners,
			    num_nodes_to_send,
			    num_nodes_to_recv,
			    comm_proc_ids,
			    comm_proc_ids);

  contact_communicate_packed_buffers(SearchComm,
				     node_comm,
				     node_info_size*sizeof(int),
				     send_node_buff,
				     recv_node_buff,
				     recv_handles);

#ifdef SHELL_DEBUG
  postream.flush();
#endif
  //-------------------------------------------------------------------
  // communications for extra node ids is complete; now fill in the exodus
  // and host code id maps, and data needed to finish the communication maps 

  int found_an_error = 0;
  node_idx = 0;
  idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    int proc_id = comm_proc_ids[i];
    for (j = 0; j < host_number_nodes_to_partner[i]; ++j){
      int node = host_comm_nodes[idx++] -1 ; // -1 for fortran
      if ( num_created_shell_nodes[node] <= 0 &&
	   res_method == SEQUENTIAL_IDS) 
	continue; // skip non-shell nodes when numbering sequentially
      if (node_status[node] == OWNED_SHARED) continue; // skip owned nodes
      if (node_owner[node] == proc_id) {
	if ( num_created_shell_nodes[node] > 0 ) { // shell node
	     
	  // loop over acme nodes sent for this host node
	  for ( k = 0; k < Num_Acme_Nodes_for_Host_Node(node); ++k) {
#ifdef SHELL_DEBUG
	    postream << my_proc_id << ":   Proc " << i << " is " << proc_id << "\n";
	    postream << my_proc_id << ":     comm node " << j << " is " << node << "\n";
	    postream << my_proc_id << ":       recv data: gid " 
		     << recv_nodes[node_idx*node_info_size +0] << " " 
		     << recv_nodes[node_idx*node_info_size +1] << " xid " 
		     << recv_nodes[node_idx*node_info_size +2] << " face " 
		     << recv_nodes[node_idx*node_info_size +3] << " " 
		     << recv_nodes[node_idx*node_info_size +4] << "\n";
#endif
	    int face_id_1 = recv_nodes[node_idx*node_info_size + 3];
	    int face_id_2 = recv_nodes[node_idx*node_info_size + 4];
	    // loop over created nodes from host code node
	    bool found_node = false;
	    for ( l = 0; l < Num_Acme_Nodes_for_Host_Node(node); ++l) {
	      int created_node = Acme_Node_for_Host_Node( node, l );
	      for ( m = 0; m < full_updated_inv_conn_count[created_node]; ++m) {
		int face_2 = 
		  full_updated_inv_conn[full_updated_offset_inv_conn[created_node]+m];
		if (face_id_1 == full_face_global_ids[2*face_2] && 
		    face_id_2 == full_face_global_ids[2*face_2 + 1] ){
		  found_node = true;
		  global_node_ids[2*created_node] = 
		    recv_nodes[node_idx*node_info_size]; 
		  global_node_ids[2*created_node + 1] =
		    recv_nodes[node_idx*node_info_size + 1]; 
		  exodus_node_ids[created_node] = 
		    recv_nodes[node_idx*node_info_size + 2];
		  break;
		}
	      } // end loop on faces connected to created node
	      if ( found_node ) break;
	    } // end loop on nodes created from host code node
	    if ( ! found_node) {
	      topology->errors->Add_Error_Message("Can't resolve node");
	      found_an_error = 1;
#ifdef SHELL_DEBUG
	      postream << "**** can't resolve node!\n"
		       << "   node is: " << node << "\n";
	      postream.flush();
#endif
	      found_an_error = contact_global_error_check( found_an_error, SearchComm );
	      return ContactSearch::INTERNAL_ERROR;
	    }
	    node_idx++;
	  } // end loop on nodal ids sent from owning proc
	} else { // non-shell node
	  PRECONDITION ( Num_Acme_Nodes_for_Host_Node(node) == 1 );
	  int created_node = Acme_Node_for_Host_Node( node, 0 );
	  global_node_ids[2*created_node] = 
	    recv_nodes[node_idx*node_info_size]; 
	  global_node_ids[2*created_node + 1] =
	    recv_nodes[node_idx*node_info_size + 1]; 
	  node_idx++;
	}
      }
    } // end loop on partner nodes
  } // end loop on communication partners

#ifdef SHELL_DEBUG
  postream.flush();
#endif
  if( contact_global_error_check( found_an_error, SearchComm ) ) 
    return ContactSearch::INTERNAL_ERROR;

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> WAITING...." << "\n";
  contact_global_sync(SearchComm);
  postream << my_proc_id << ":TEST>  >> Whew! here are the resolved parallel"
       << " ids:" << "\n";
  for (i =0 ; i < number_acme_nodes; ++i){
    postream << my_proc_id << ":TEST>      node: " << i << " HID: "
	 << global_node_ids[2*i] << "," << global_node_ids[2*i+1]
	 << "  EXID: " << exodus_node_ids[i] << "\n";
  }
#endif

  //-------------------------------------------------------
  // now we have all global ids for the acme nodes from
  // the host code nodes. Now re-order the
  // host->acme node map so that the acme nodes are
  // in order of globally-resolved exodus_id.
  for ( i = 0 ; i < number_host_code_nodes; i ++) {
    for ( k = 0; k < Num_Acme_Nodes_for_Host_Node(i) - 1; ++k) {
      int created_node1 = Acme_Node_for_Host_Node( i, k );
      for ( l = k+1; l < Num_Acme_Nodes_for_Host_Node(i); ++l) {
	int created_node2 = Acme_Node_for_Host_Node( i, l );
	if (exodus_node_ids[created_node2] < 
	    exodus_node_ids[created_node1] ){
	  int tmp = host_to_acme_node_map[offset_for_node[i] + k];
	  host_to_acme_node_map[offset_for_node[i] + k] = 
	    host_to_acme_node_map[offset_for_node[i] + l];
	  host_to_acme_node_map[offset_for_node[i] + l] = tmp;
	  created_node1 = created_node2;
	}
      }
    }
  }


  //-------------------------------------------------------
  // All of the ids are now syncronized across processors. Now we
  // need to create new communication maps for the created shell nodes
  idx = 0;
  int new_idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    number_nodes_to_partner[i] = 0;
    for (j = 0; j < host_number_nodes_to_partner[i]; ++j){
      int node = host_comm_nodes[idx++] -1 ; // -1 for fortran
      if (num_created_shell_nodes[node] > 0) {
	for (l = 0; l < Num_Acme_Nodes_for_Host_Node(node); ++l)
	  comm_nodes[new_idx++] =
	    Acme_Node_for_Host_Node( node, l ) + 1;// +1 for fortran
	number_nodes_to_partner[i] += Num_Acme_Nodes_for_Host_Node(node);
      } else {
	// comm node is not a shell node.
	comm_nodes[new_idx++] = node + 1; // +1 for fortran
	number_nodes_to_partner[i]++;
      }
    } // end loop on nodes to comm partner
  } // end loop on comm partners

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> Here is original comm map." << "\n";
  idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    int proc_id = comm_proc_ids[i];
    postream << my_proc_id << ":TEST>   partner " << i << " is proc " 
	 << proc_id << "\n";
    for (j = 0; j < host_number_nodes_to_partner[i]; ++j){
      int node = host_comm_nodes[idx++]-1;
      postream << my_proc_id << ":TEST>        node " << node << " (orig id: " 
	   << orig_global_node_ids[2*node] << " " 
	   << orig_global_node_ids[2*node+1]
	   << ")"<< "\n";
    }
  }
  postream << my_proc_id << ":TEST> Here is updated comm map." << "\n";
  idx = 0;
  for (i = 0; i < number_comm_partners; ++i) {
    int proc_id = comm_proc_ids[i];
    postream << my_proc_id << ":TEST>   partner " << i << " is proc " 
	 << proc_id << "\n";
    for (j = 0; j < number_nodes_to_partner[i]; ++j){
      int node = comm_nodes[idx++]-1;
      postream << my_proc_id << ":TEST>        node " << node << " (id: " 
	   << global_node_ids[2*node] << " " <<  global_node_ids[2*node+1]
	   << ")"<< "\n";
    }
  }
#endif

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST>=============Done finalizing parallel node data" << "\n";
  postream.flush();
#endif

  delete [] node_status;
  delete [] node_owner;
  delete [] num_nodes_to_recv;
  delete [] num_nodes_to_send;

  return ContactSearch::NO_ERROR;
}

#endif  // CONTACT_NO_MPI

void ContactShellHandler::number_nodes_sequentially
               ( 
		const int num_owned_host_code_nodes,
		const int num_owned_newly_created_nodes,
		const NodeStatus * node_status,
		int* exodus_node_ids,
		int* global_node_ids
		)
{
  int i,j;

  // Find the highest GID number from existing nodes and the new host 
  // code nodes and communicate this in parallel.
  int starting_gid[2], top_gid[2];
  top_gid[0] = -1; 
  top_gid[1] = -1; 
  for (i = 0; i < number_host_code_nodes; ++i) {
    top_gid[0] = std::max(top_gid[0], orig_global_node_ids[2*i]);
    top_gid[1] = std::max(top_gid[1], orig_global_node_ids[2*i+1]);
  }
  contact_global_maximum(top_gid, starting_gid, 2, SearchComm);
  starting_gid[1]++;
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> starting gid: " 
	   << starting_gid[0] << " | " << starting_gid[1] << "\n";
#endif

  // now do an all-to-all broadcast of the number of local nodes on 
  // each processor so that each processor can compute where to start
  // numbering the nodes
  
  int num_procs = contact_number_of_processors(SearchComm);
  int * num_new_nodes_local_to_proc = new int[num_procs];
  num_new_nodes_local_to_proc[my_proc_id] = num_owned_newly_created_nodes;

  // crummy communication, but couldn't std::right off think of a better way..
  for (i = 0; i < num_procs; ++i)
    contact_broadcast(&num_new_nodes_local_to_proc[i], 1, i, SearchComm);

  // now use this info to compute the local numbers for this processor
  int numbering_base = 0;
  for (i = 0; i < my_proc_id; ++i)
    numbering_base += num_new_nodes_local_to_proc[i];

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> numbering_base: " << numbering_base << "\n";
#endif
  
  starting_gid[1] += numbering_base;

  delete [] num_new_nodes_local_to_proc;

  // recompute the global ids and exodus ids for all nodes. In general, we 
  // want all non-shell nodes and one of the created nodes for each shell
  // node to still have the same exodus id and global id as the host code
  // node. This just makes it a little easier to interpret the results. 
  // Other nodes created from the shell nodes then must also have a 
  // unique id assigned to them, following the numbering_base as computed
  // above. Also, since multiple node blocks may have caused all of 
  // the local ids for the nodes to be renumbered, we map each nodes back to 
  // the host code node to make sure that our numbering scheme above is
  // adhered to.
#ifdef SHELL_DEBUG
  postream << ":New total local numbering:\n";
#endif
  int created_node;
  for (i = 0; i < number_host_code_nodes; ++i) {
    if ( node_status[i] == NOT_OWNED_SHARED ) continue; // skip non owned nodes
#ifdef SHELL_DEBUG
    postream << "  hcnode " << i << "\n";
#endif
    created_node = Acme_Node_for_Host_Node( i, 0 );
    exodus_node_ids[created_node] = orig_exodus_node_ids[i];
    global_node_ids[2*created_node] = orig_global_node_ids[2*i];
    global_node_ids[2*created_node+1] = orig_global_node_ids[2*i+1];

#ifdef SHELL_DEBUG
    postream << "    created #0: " << created_node 
	     << " exid" <<  exodus_node_ids[created_node]
	     << " gid: " <<  global_node_ids[2*created_node]
	     << " " << global_node_ids[2*created_node+1] << "\n";
#endif
    
    if ( num_created_shell_nodes[i] <= 0 ) continue; // skip non shell nodes

    for ( j = 1; j < Num_Acme_Nodes_for_Host_Node(i); ++j) {
      created_node = Acme_Node_for_Host_Node( i, j );
      exodus_node_ids[created_node] = starting_gid[1];
      global_node_ids[2*created_node] = starting_gid[0];
      global_node_ids[2*created_node+1] = starting_gid[1];
      starting_gid[1]++;
#ifdef SHELL_DEBUG
      postream << "    created #" << i << ": " << created_node 
	       << " exid" <<  exodus_node_ids[created_node]
	       << " gid: " <<  global_node_ids[2*created_node]
	       << " " << global_node_ids[2*created_node+1] << "\n";
#endif
    }
  }

}

void ContactShellHandler::number_nodes_using_handler
               ( 
		const int num_owned_host_code_nodes,
		const int num_owned_newly_created_nodes,
		const NodeStatus * node_status,
		const int* inv_conn,
		const int* inv_conn_count,
		const int* offset_inv_conn,
		const int* full_shell_face_blk,
		const int* full_shell_face_id,
		const int* number_faces_per_block,
		const int* global_face_ids,
		ContactShellHandler * old_handler,
		int* exodus_node_ids,
		int* global_node_ids
		)
{
  int i,j,k,l;

#ifdef SHELL_DEBUG
  postream << "In number_nodes_using_handler\n";
  postream.flush();
#endif

  // first, we need to find correspondence between new nodes and old nodes
  // if we can match them up, then have new nodes acquire the old nodes'
  // GID and EXID

  int num_new_created_nodes = 0;

  // loop over all new host code nodes
  for ( i = 0; i < number_host_code_nodes; ++i) {

#ifdef SHELL_DEBUG
    postream << "  * process host code node " << i << "\n";
#endif
    // look for node in old handler that has same gid as new node
    int new_gid[2], old_gid[2];
    GID_of_Host_Code_Node(i,new_gid[0],new_gid[1]);
   
    bool found_node = false;
    int old_host_node = -1;
    for ( j = 0; j < old_handler->Number_Host_Code_Nodes(); ++j) {
      old_handler->GID_of_Host_Code_Node(j,old_gid[0],old_gid[1]);
      if (old_gid[0] == new_gid[0] && old_gid[1] == new_gid[1]) {
	found_node = true;
	old_host_node = j;
	break;
      }
    }

    if (found_node) {
      // found the GID from old handler -- carry over IDS if possible
      
#ifdef SHELL_DEBUG
      postream << "    * matched node with gid " << new_gid[0] << " | "
            << new_gid[1] << "\n";
#endif

      // if only one acme node, then its a std::hex node. Assume we found it.
      if ( Num_Acme_Nodes_for_Host_Node(i) == 1 ) {
	int created_node = Acme_Node_for_Host_Node( i, 0 );
	global_node_ids[2*created_node] = new_gid[0];
	global_node_ids[2*created_node+1] = new_gid[1];
	exodus_node_ids[created_node] = new_gid[1];
	continue;
      }
      
      // flag GIDS as needing to be set
      for ( j = 0; j < Num_Acme_Nodes_for_Host_Node(i); ++j) {
	int created_node = Acme_Node_for_Host_Node( i, j );
	global_node_ids[2*created_node] = -1;
	global_node_ids[2*created_node+1] = -1;
	exodus_node_ids[created_node] = -1;
      }

      // loop over the new nodes for this host code node      
      for ( j = 0; j < Num_Acme_Nodes_for_Host_Node(i); ++j) {
	int created_node = Acme_Node_for_Host_Node( i, j );
	
#ifdef SHELL_DEBUG
	postream << "      * looping over node " << j+1 << " of "
		 << Num_Acme_Nodes_for_Host_Node(i) << ": "
		 << created_node << "\n";
#endif

	bool found_face = false;
	int resolved_node_gid[2];
	// loop over faces connected to this node
	for ( k = 0; k< inv_conn_count[created_node]; ++k){
	  int face = inv_conn[offset_inv_conn[created_node]+k];

#ifdef SHELL_DEBUG
	  int face_blk = full_shell_face_blk[face];
	  int face_id = full_shell_face_id[face];
	  postream << "         * looping over face " << k+1 << " of "
	    << inv_conn_count[created_node] << ": blk "
	    << face_blk << " id " << face_id << "\n";
#endif

 	  int face_gid[2], old_face_gid[2];
 	  face_gid[0] = global_face_ids[2*face];
 	  face_gid[1] = global_face_ids[2*face+1];
	  
#ifdef SHELL_DEBUG
	  postream << "           * face gid: " << face_gid[0] 
	    << " | " << face_gid[1] << "\n";
#endif

	  for ( l = 0; 
		l < old_handler->Num_Acme_Nodes_for_Host_Node(old_host_node); 
		++l) {
	    int old_acme_node = old_handler->Acme_Node_for_Host_Node( old_host_node, l );
	    ContactNode<Real>* old_node_obj = static_cast<ContactNode<Real>*>(topology->NodeList()->Find(old_acme_node));
            if(old_node_obj->ConnectedToFace(ContactHostGlobalID(old_face_gid[0], old_face_gid[1]))) {
              found_face = true;
              resolved_node_gid[0] = old_node_obj->Global_ID().HiInt();
	      resolved_node_gid[1] = old_node_obj->Global_ID().LoInt();
	    }
	    if (found_face) break;
	  }

	  if (found_face) break;
	} // end loop on connected faces
	    
	if ( found_face) {
	  // loop over previous created nodes for this new host code node to
	  // make sure we have not already used this resolved node id
	  bool already_used_node_gid = false;
	  for (k = 0; k < j; ++k){
	    int prev_node = Acme_Node_for_Host_Node( i, j );
	    if (global_node_ids[2*prev_node]   == resolved_node_gid[0] &&
		global_node_ids[2*prev_node+1] == resolved_node_gid[1] ) {
	      already_used_node_gid = true;
	      break;
	    }
	  }	    

	  if ( already_used_node_gid ) {
	    num_new_created_nodes ++;
	  } else {
	    global_node_ids[2*created_node] = resolved_node_gid[0];
	    global_node_ids[2*created_node+1] = resolved_node_gid[1];
	    exodus_node_ids[created_node] = resolved_node_gid[1];
	  }
	} else {
	  num_new_created_nodes ++;
	}
	
      } // end loop on created nodes from host code node
    
    } else {

#ifdef SHELL_DEBUG
      postream << "    * did not match node.\n";
#endif
      // did not find the GID from old handler -- this must be a new node.
      // mark ids as needing a new id number
      num_new_created_nodes += Num_Acme_Nodes_for_Host_Node(i);
      for ( j = 0; j < Num_Acme_Nodes_for_Host_Node(i); ++j) {
	int created_node = Acme_Node_for_Host_Node( i, j );
	global_node_ids[2*created_node] = -1;
	global_node_ids[2*created_node+1] = -1;
	exodus_node_ids[created_node] = -1;
      }
    }
  } // end loop on host code nodes 

#ifdef SHELL_DEBUG
  postream << " Here are the new gids:\n";
  for (i = 0; i < number_host_code_nodes; ++i) {
    int id[2];
    GID_of_Host_Code_Node(i,id[0],id[1]);
    postream << "  -> host code node " << i << " gid:"
	     << id[0] << " | " << id[1] << "\n";
    for ( j = 0; j < Num_Acme_Nodes_for_Host_Node(i); ++j) {
      int created_node = Acme_Node_for_Host_Node( i, j );
      id[0] = global_node_ids[2*created_node];
      id[1] = global_node_ids[2*created_node+1];
      postream << "    -> created node #" << created_node << " gid:"
	       << id[0] << " | " << id[1] << "\n";
    }
  }
#endif

  //-------------------------
  // Resolve ids of any new nodes created from update
  //
  // Find the highest GID number from existing nodes and the new host 
  // code nodes and communicate this in parallel.
  int starting_gid[2], top_gid[2];
  top_gid[0] = -1; 
  top_gid[1] = -1; 
  for (i = 0; i < number_acme_nodes; ++i) {
    top_gid[0] = std::max(top_gid[0], new_global_node_ids[2*i]);
    top_gid[1] = std::max(top_gid[1], new_global_node_ids[2*i+1]);
  }
  contact_global_maximum(top_gid, starting_gid, 2, SearchComm);
  starting_gid[1]++;
#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> starting gid: " 
	   << starting_gid[0] << " | " << starting_gid[1] << "\n";
#endif

  // now do an all-to-all broadcast of the number of local nodes on 
  // each processor so that each processor can compute where to start
  // numbering the nodes
  
  int num_procs = contact_number_of_processors(SearchComm);
  int * num_new_nodes_local_to_proc = new int[num_procs];
  num_new_nodes_local_to_proc[my_proc_id] = num_new_created_nodes;

  // crummy communication, but couldn't std::right off think of a better way..
  for (i = 0; i < num_procs; ++i)
    contact_broadcast(&num_new_nodes_local_to_proc[i], 1, i, SearchComm);

  // now use this info to compute the local numbers for this processor
  int numbering_base = 0;
  for (i = 0; i < my_proc_id; ++i)
    numbering_base += num_new_nodes_local_to_proc[i];

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":TEST> numbering_base: " << numbering_base << "\n";
#endif
  
  starting_gid[1] += numbering_base;

  delete [] num_new_nodes_local_to_proc;
  
  // compute the global ids and exodus ids for new nodes. Because we don't
  // have any idea about the ids that were used before, we can't preserve
  // the exodus or gids for non-shell nodes or one of the nodes created from 
  // a shell node, as we did in the initial construction.
#ifdef SHELL_DEBUG
  postream << ":New total local numbering:\n";
#endif
  for (i = 0; i < number_host_code_nodes; ++i) {
    if ( num_procs > 1) { 
      if ( node_status[i] == NOT_OWNED_SHARED ) continue; // skip non owned nodes
    }
#ifdef SHELL_DEBUG
    postream << "  hcnode " << i << "\n";
#endif
    for ( j = 0; j < Num_Acme_Nodes_for_Host_Node(i); ++j) {
      int created_node = Acme_Node_for_Host_Node( i, j );
      if (global_node_ids[2*created_node+1]>=0) continue;
      exodus_node_ids[created_node] = starting_gid[1];
      global_node_ids[2*created_node] = starting_gid[0];
      global_node_ids[2*created_node+1] = starting_gid[1];
      starting_gid[1]++;
#ifdef SHELL_DEBUG
      postream << "    created #" << i << ": " << created_node 
	       << " exid" <<  exodus_node_ids[created_node]
	       << " gid: " <<  global_node_ids[2*created_node]
	       << " " << global_node_ids[2*created_node+1] << "\n";
#endif
    }
  }

}


void ContactShellHandler::Set_NodeBlk_Position( int NodeBlk_ID,
                                                ContactBlockEntityList* nodes,
						VariableHandle POSITION,
						const Real* positions )
{
  int i,j,k;

  // This function takes the nodal coordinates from the host code numbering
  // and maps it through Host->Acme node maps to put the positions in the
  // actual topology.  This does NOT do the lofting (we want that done later
  // after we can build physical faces, etc.)
  int num_nodes_in_block = host_nodes_per_block[NodeBlk_ID-1];
  int node_offset = block_offset_for_nodes[NodeBlk_ID-1];
  for( i=0 ; i<num_nodes_in_block ; ++i){
    for ( k =0; k < Num_Acme_Nodes_for_Host_Node(node_offset+i); ++k) {
      int hi, lo;
      Acme_NodeGID_for_Host_Node(node_offset+i,k,hi,lo);
      ContactHostGlobalID gid(hi, lo);
      ContactNode<Real>* node = 
        static_cast<ContactNode<Real>*>(nodes->Find(gid));
      POSTCONDITION(node);
      Real* position = node->Variable(POSITION);
      for( j=0 ; j<dimensionality ; ++j) {
	position[j] = positions[dimensionality*i + j];
      }
    }
  }
}

void 
ContactShellHandler::Set_NodeBlk_KinConstr( int NodeBlk_ID,
                                            ContactBlockEntityList* nodes,
					    VariableHandle NUM_KIN_CONSTR,
					    VariableHandle KIN_CONSTR_VECTOR,
					    const int* num_kcs,
					    const Real* kc_vectors)
{
  int i,j,k;

  // This function takes the nodal constraints from the host code numbering
  // and maps it through Host->Acme node maps to put the constraints in the
  // actual topology. 
  int num_nodes_in_block = host_nodes_per_block[NodeBlk_ID-1];
  int node_offset = block_offset_for_nodes[NodeBlk_ID-1];
  for( i=0 ; i<num_nodes_in_block ; ++i){
    for ( k =0; k < Num_Acme_Nodes_for_Host_Node(node_offset+i); ++k) {
      int hi, lo;
      Acme_NodeGID_for_Host_Node(node_offset+i,k,hi,lo);
      ContactHostGlobalID gid(hi, lo);
      ContactNode<Real>* node = 
        static_cast<ContactNode<Real>*>(nodes->Find(gid));
      POSTCONDITION(node);
      *node->Variable(NUM_KIN_CONSTR) = num_kcs[i];
      Real* vector = node->Variable(KIN_CONSTR_VECTOR);
      for( j=0 ; j<dimensionality ; ++j)
	vector[j] = kc_vectors[dimensionality*i + j];
    }
  } 
}

void
ContactShellHandler::Set_NodeBlk_RemainingGap( int NodeBlk_ID,
                                               ContactBlockEntityList* nodes,
					       VariableHandle REMAINING_GAP,
					       const Real* gap)
{
  int i,j,k;

  // This function takes the remaining gap from the host code numbering
  // and maps it through Host->Acme node maps to put the constraints in the
  // actual topology. 
  int num_nodes_in_block = host_nodes_per_block[NodeBlk_ID-1];
  int node_offset = block_offset_for_nodes[NodeBlk_ID-1];
  for( i=0 ; i<num_nodes_in_block ; ++i){
    for ( k =0; k < Num_Acme_Nodes_for_Host_Node(node_offset+i); ++k) {
      int hi, lo;
      Acme_NodeGID_for_Host_Node(node_offset+i,k,hi,lo);
      ContactHostGlobalID gid(hi, lo);
      ContactNode<Real>* node = 
        static_cast<ContactNode<Real>*>(nodes->Find(gid));
      POSTCONDITION(node);
      Real* remaining_gap = node->Variable(REMAINING_GAP);
      for( j=0 ; j<dimensionality ; ++j)
	remaining_gap[j] = gap[dimensionality*i + j];
    }
  } 
}

void
ContactShellHandler::Set_NodeBlk_GhostingGap( int NodeBlk_ID,
                                              ContactBlockEntityList* nodes,
					      VariableHandle GHOSTING_GAP,
					      const Real* gap)
{
  int i,j,k;

  // This function takes the remaining gap from the host code numbering
  // and maps it through Host->Acme node maps to put the constraints in the
  // actual topology. 
  int num_nodes_in_block = host_nodes_per_block[NodeBlk_ID-1];
  int node_offset = block_offset_for_nodes[NodeBlk_ID-1];
  for( i=0 ; i<num_nodes_in_block ; ++i){
    for ( k =0; k < Num_Acme_Nodes_for_Host_Node(node_offset+i); ++k) {
      int hi, lo;
      Acme_NodeGID_for_Host_Node(node_offset+i,k,hi,lo);
      ContactHostGlobalID gid(hi, lo);
      ContactNode<Real>* node = 
        static_cast<ContactNode<Real>*>(nodes->Find(gid));
      POSTCONDITION(node);
      Real* ghosting_gap = node->Variable(GHOSTING_GAP);
      for( j=0 ; j<dimensionality ; ++j)
	ghosting_gap[j] = gap[dimensionality*i + j];
    }
  } 
}

void ContactShellHandler::Complete_Initialization()
{
  int i,k;
  int num_ln_set = 0;
  for( i=0 ; i<number_host_code_nodes ; ++i){
    if( num_created_shell_nodes[i] > 1 ){
      for( k =0; k < Num_Acme_Nodes_for_Host_Node(i); ++k) {
        int hi, lo;
        Acme_NodeGID_for_Host_Node(i,k,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactShellNode* node = static_cast<ContactShellNode*> 
	  (topology->NodeList()->Find(gid));
        POSTCONDITION(node);
	PRECONDITION(node->Is_a_Shell_Node());
	lofted_nodes[num_ln_set++] = node;
	POSTCONDITION( num_ln_set <= number_lofted_nodes );
      }
    }
  }
  POSTCONDITION( num_ln_set == number_lofted_nodes );
#ifndef CONTACT_NO_MPI
  Build_Shell_Comm_Plans();
#endif
  Mark_Shell_Nodes();
#ifndef CONTACT_NO_MPI
  Build_Loft_Nodes_Data_Structures();
#endif

}

#ifndef CONTACT_NO_MPI
void ContactShellHandler::Build_Loft_Nodes_Data_Structures()
{
  int i,j,k;
  int number_of_nodes = topology->Number_of_Nodes();

  // First determine how many ghost faces each node will get
  char* send_buff;
  char* recv_buff;
  RequestHandle* recv_handles;
  int size_per_object = 1*sizeof(int);
  topology->comm_buffer->Buffers( *comm_plan_to_owner, size_per_object,
				  &send_buff, &recv_buff, &recv_handles );
  int* sb_i = (int*) send_buff;
  int offset = 0;
  for( i=0 ; i<comm_plan_to_owner->Num_Export_Comm_Partners() ; ++i){
    ContactNode<Real>** node_list_to_proc = (ContactNode<Real>**) 
      comm_plan_to_owner->Export_Entity_List( i );
    for( j=0 ; j<comm_plan_to_owner->Num_Export_to_Proc( i ) ; ++j){
      if ( node_list_to_proc[j]->Physical_Type() == ContactNode<Real>::MIXED_NODE)
	sb_i[offset++] = 0;
      else
	sb_i[offset++] = node_list_to_proc[j]->Number_Face_Connections();
    }
  }
  contact_communicate_packed_buffers( SearchComm, *comm_plan_to_owner, 
				      sizeof(int), send_buff, recv_buff, 
				      recv_handles );

  // Now using this information build another comm plan that accounts for
  // how many ghost faces each node will get.  (Note: I will include the
  // same node N times to a processor; once for each ghost face it will be
  // receiving).

  // Size the import part of the plan
  int* rb_i = (int*) recv_buff;
  offset = 0;
  int num_import_procs = 0;
  int num_imported_faces = 0;
  for( i=0 ; i<comm_plan_to_owner->Num_Import_Comm_Partners() ; ++i){
    bool added_this_proc = false;
    for( j=0 ; j<comm_plan_to_owner->Num_Import_from_Proc(i) ; ++j){
      if( rb_i[offset] > 0 ){
	num_imported_faces += rb_i[offset];
	if( !added_this_proc ){
	  added_this_proc = true;
	  num_import_procs++;
	}
      }
      offset++;
    }
  }

  // Fill the import part of the plan
  rb_i = (int*) recv_buff;
  int* Num_Import_from_Proc = new int[num_import_procs];
  int* Import_Comm_Proc_IDs = new int[num_import_procs];
  ContactTopologyEntity<Real>** Import_Entity_List = new ContactTopologyEntity<Real>*[num_imported_faces];
  std::memset( Num_Import_from_Proc, 0, num_import_procs*sizeof(int) );
  num_import_procs = 0;
  num_imported_faces = 0;
  offset = 0;
  for( i=0 ; i<comm_plan_to_owner->Num_Import_Comm_Partners() ; ++i){
    bool added_this_proc = false;
    ContactTopologyEntity<Real>** node_list = comm_plan_to_owner->Import_Entity_List( i );
    for( j=0 ; j<comm_plan_to_owner->Num_Import_from_Proc(i) ; ++j){
      if( rb_i[offset] > 0 ){
	for( k=0 ; k<rb_i[offset] ; ++k){
	  Import_Entity_List[num_imported_faces++] = node_list[j];
	}
	if( !added_this_proc ){
	  added_this_proc = true;
	  Import_Comm_Proc_IDs[num_import_procs++] = 
	    comm_plan_to_owner->Import_Comm_Proc_ID(i);
	}
	Num_Import_from_Proc[ num_import_procs-1 ] += rb_i[offset];
      }
      offset++;
    }
  }
  
  // Size the export part of the plan
  offset = 0;
  sb_i = (int*) send_buff;
  int num_export_procs = 0;
  int num_exported_faces = 0;
  for( i=0 ; i<comm_plan_to_owner->Num_Export_Comm_Partners() ; ++i){
    bool added_this_proc = false;
    for( j=0 ; j<comm_plan_to_owner->Num_Export_to_Proc(i) ; ++j){
      if( sb_i[offset] > 0 ){
	num_exported_faces += sb_i[offset];
	if( !added_this_proc ){
	  added_this_proc = true;
	  num_export_procs++;
	}
      }
      offset++;
    }
  }

  // Fill the export part of the plan
  sb_i = (int*) send_buff;
  int* Num_Export_to_Proc = new int[num_export_procs];
  int* Export_Comm_Proc_IDs = new int[num_export_procs];
  std::memset( Num_Export_to_Proc, 0, num_export_procs*sizeof(int) );
  ContactTopologyEntity<Real>** Export_Entity_List = new ContactTopologyEntity<Real>*[num_exported_faces];
  num_export_procs = 0;
  num_exported_faces = 0;
  offset = 0;
  for( i=0 ; i<comm_plan_to_owner->Num_Export_Comm_Partners() ; ++i){
    bool added_this_proc = false;
    ContactTopologyEntity<Real>** node_list = comm_plan_to_owner->Export_Entity_List( i );
    for( j=0 ; j<comm_plan_to_owner->Num_Export_to_Proc(i) ; ++j){
      if( sb_i[offset] > 0 ){
	for( k=0 ; k<sb_i[offset] ; ++k){
	  Export_Entity_List[num_exported_faces++] = node_list[j];
	}
	if( !added_this_proc ){
	  added_this_proc = true;
	  Export_Comm_Proc_IDs[num_export_procs++] =
	    comm_plan_to_owner->Export_Comm_Proc_ID(i);
	}
	Num_Export_to_Proc[ num_export_procs-1 ] += sb_i[offset];
      }
      offset++;
    }
  }
  
  import_loft_face_data_comm_plan = new ContactAsymComm( num_export_procs,
							 num_import_procs,
							 Num_Export_to_Proc,
							 Num_Import_from_Proc,
							 Export_Comm_Proc_IDs,
							 Import_Comm_Proc_IDs,
							 Export_Entity_List,
							 Import_Entity_List );

#ifdef SHELL_DEBUG  
  //  ContactParOStream postream(SearchComm);
  if( contact_processor_number(SearchComm) == 0 ){
    postream << "\n\nImport Loft Face Comm Plan\n";
    postream << "==========================================\n";
  }
  postream << "Export Information\n";
  postream << "   Num Comm Partners = " 
	   << import_loft_face_data_comm_plan->Num_Export_Comm_Partners() << "\n";
  for( i=0 ; i<import_loft_face_data_comm_plan->Num_Export_Comm_Partners() ; ++i){
    postream << "    Exporting to Proc " 
	     << import_loft_face_data_comm_plan->Export_Comm_Proc_ID(i) << "\n";
    for( j=0 ; j<import_loft_face_data_comm_plan->Num_Export_to_Proc(i) ; ++j){
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                          (import_loft_face_data_comm_plan->Export_Entity_List(i)[j]);
      postream << "       Communicating Node " 
	       << node->Exodus_ID()
	       << " " 
	       << node->Global_ID()
	       << "\n";
      int nfaces = node->Number_Face_Connections();
      postream << "          Faces:";
      for (int mwg=0; mwg<nfaces; mwg++) {
        postream << node->Face(mwg)->Global_ID();
      }
      postream << "\n";
    }
  }
  postream << "Import Information\n";
  postream << "   Num Comm Partners = " 
	   << import_loft_face_data_comm_plan->Num_Import_Comm_Partners() << "\n";
  for( i=0 ; i<import_loft_face_data_comm_plan->Num_Import_Comm_Partners() ; ++i){
    postream << "    Importing to Proc " 
	     << import_loft_face_data_comm_plan->Import_Comm_Proc_ID(i) << "\n";
    for( j=0 ; j<import_loft_face_data_comm_plan->Num_Import_from_Proc(i) ; ++j){
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                          (import_loft_face_data_comm_plan->Import_Entity_List(i)[j]);
      postream << "       Communicating Node " 
	       << node->Exodus_ID()
	       << " " 
	       << node->Global_ID()
	       << "\n";
      //int nfaces = node->Number_Face_Connections();
      //postream << "          Faces:";
      //for (int mwg=0; mwg<nfaces; mwg++) {
      //  postream << node->Face(mwg)->Global_ID();
      //}
      //postream << "\n";
    }
  }
  postream.flush();
#endif // SHELL_DEBUG

  // Now build the data structures to use in Loft_Nodes
  number_lofting_ghosts = new int[number_of_nodes];
  offset_to_buffer_offsets = new int[number_of_nodes];
  for( i=0 ; i<number_of_nodes ; ++i){
    number_lofting_ghosts[i] = 0;
    offset_to_buffer_offsets[i] = -1;
  }
  int number_of_ghosts = 0;
  for( i=0;i<import_loft_face_data_comm_plan->Num_Import_Comm_Partners(); ++i){
    ContactNode<Real>** node_list = (ContactNode<Real>**) 
      import_loft_face_data_comm_plan->Import_Entity_List(i);
    for(j=0;j<import_loft_face_data_comm_plan->Num_Import_from_Proc(i); ++j){
      number_lofting_ghosts[node_list[j]->ProcArrayIndex()] += 1;
      number_of_ghosts++;
    }
  }
  int next_offset = 0;
  for( i=0 ; i<number_of_nodes ; ++i){
    if( number_lofting_ghosts[i] > 0 ){
      offset_to_buffer_offsets[i] = next_offset;
      next_offset += number_lofting_ghosts[i];
    }
  }

  buffer_offsets = new int[number_of_ghosts];
  for( i=0 ; i<number_of_ghosts ; ++i) 
    buffer_offsets[i] = -1;

  offset = 0;
  for( i=0;i<import_loft_face_data_comm_plan->Num_Import_Comm_Partners(); ++i){
    ContactNode<Real>** node_list = (ContactNode<Real>**) 
      import_loft_face_data_comm_plan->Import_Entity_List(i);
    for( j=0;j<import_loft_face_data_comm_plan->Num_Import_from_Proc(i); ++j){
      int node_index = node_list[j]->ProcArrayIndex();
      int* node_buffer_offsets =
	buffer_offsets + offset_to_buffer_offsets[node_index];
      for( k=0 ; k<number_lofting_ghosts[node_index] ; ++k){
	if( node_buffer_offsets[k] == -1 ){
	  node_buffer_offsets[k] = offset++;
	  break;
	}
      }
      POSTCONDITION( k != number_lofting_ghosts[node_index] );
    }
  }  
#ifdef SHELL_DEBUG
  if( contact_processor_number(SearchComm) == 0 ){
    postream << "\n\nLofted Communication Information\n";
    postream << "=====================================================\n";
  }
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    postream << "Information for Node " << Nodes[i]->Exodus_ID() << "\n";
    postream << "   Number of Lofting Ghosts = " << number_lofting_ghosts[i] 
	     << "\n";
    postream << "    Buffer Offset   = " 
             << offset_to_buffer_offsets[i] << "\n";
    for( j=0 ; j<number_lofting_ghosts[i] ; ++j)
      postream << "    Buffer Location = " 
	       << buffer_offsets[offset_to_buffer_offsets[i]+j] << "\n";
  }
  postream.flush();
#endif

  delete [] Num_Import_from_Proc;
  delete [] Import_Comm_Proc_IDs;
  delete [] Num_Export_to_Proc;
  delete [] Export_Comm_Proc_IDs;
  delete [] Import_Entity_List;
  delete [] Export_Entity_List;
}
#endif

void ContactShellHandler::NumberNodes() {
  for( int i=0 ; i<number_host_code_nodes ; ++i){

    int host_gid_lo;
    int host_gid_hi;

    GID_of_Host_Code_Node(i, host_gid_lo, host_gid_hi);

    if( num_created_shell_nodes[i] > 1 ){
      for( int k =0; k < Num_Acme_Nodes_for_Host_Node(i); ++k) {
        int hi, lo;
        Acme_NodeGID_for_Host_Node(i,k,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactShellNode* node = static_cast<ContactShellNode*> 
	  (topology->NodeList()->Find(gid));
        node->Shell_Node_Base_ID( host_gid_hi );
        POSTCONDITION(node);
	POSTCONDITION(node->Is_a_Shell_Node());
      }
    }
  }
}

void ContactShellHandler::Mark_Shell_Nodes( )
{
  //std::cerr << "here in ContactShellHandler::Mark_Shell_Nodes: number_lofted_nodes = " << number_lofted_nodes << std::endl;
  int i,j;

  // loop over shell nodes. Check connected faces. If any face is
  // not a shell face, then this is a mixed node. set tmp tag to 1
  // and do a parallel swapadd to make sure all processors know that
  // this node is a mixed node
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->temp_tag = 0;
  }

  for( i=0 ; i<number_lofted_nodes ; ++i) {
   if(!lofted_nodes[i]->ConnectedToAllShellFaces()) lofted_nodes[i]->temp_tag = 1;
  }
  
#ifndef CONTACT_NO_MPI
  int * comm_data = new int[topology->Number_of_Nodes()];
  for ( i = 0; i < topology->Node_Sym_Comm()->Num_Comm_Partners(); ++i) {
    ContactTopologyEntity<Real> ** entity_list = topology->Node_Sym_Comm()->Entity_List(i);
    for ( j = 0; j < topology->Node_Sym_Comm()->Num_to_Proc(i); ++j) 
      comm_data[entity_list[j]->ProcArrayIndex()] = entity_list[j]->temp_tag;
  }
  
  contact_swapadd_data_array(SearchComm, *topology->Node_Sym_Comm(), *comm_buffer,
			     comm_data, 1);
  
  for ( i = 0; i < topology->Node_Sym_Comm()->Num_Comm_Partners(); ++i) {
    ContactTopologyEntity<Real> ** entity_list = topology->Node_Sym_Comm()->Entity_List(i);
    for ( j = 0; j < topology->Node_Sym_Comm()->Num_to_Proc(i); ++j) 
      entity_list[j]->temp_tag += comm_data[entity_list[j]->ProcArrayIndex()];
  }
  if (comm_data) delete [] comm_data;
#endif
  
  int lofted_node_idx = 0;
  for( i=0 ; i<number_host_code_nodes ; ++i) {
    // want all shell nodes created from single host code shell node
    // all to share a unique id so that if we have two created shell
    // nodes we can tell if they are related. Use the host code 
    // exodus ID for this number.
    if ( num_created_shell_nodes[i] <= 1 ) continue; 
    for ( int k =0; k < Num_Acme_Nodes_for_Host_Node(i); ++k) {
      int contact_id = Acme_Node_for_Host_Node(i,k);
      int host_id = acme_to_host_node_map[contact_id];
      POSTCONDITION(host_id<number_host_code_nodes);
      int host_ex_id = orig_exodus_node_ids[host_id];
      lofted_nodes[lofted_node_idx]->Shell_Node_Base_ID( host_ex_id );
      if (lofted_nodes[lofted_node_idx]->temp_tag > 0) {
	lofted_nodes[lofted_node_idx]->Physical_Type(ContactNode<Real>::MIXED_NODE);
#ifdef SHELL_DEBUG
	postream << "=> node " << lofted_nodes[lofted_node_idx]->Exodus_ID()
		 << " (CID:" << contact_id << ")" 
		 << " is a MIXED node" << "\n";
#endif
      } else {
	lofted_nodes[lofted_node_idx]->Physical_Type(ContactNode<Real>::SHELL_NODE);
#ifdef SHELL_DEBUG
	postream << "=> node " << lofted_nodes[lofted_node_idx]->Exodus_ID()
		 << " (CID:" << contact_id << ")" 
		 << " is a pure SHELL node" << "\n";
#endif
      }
      lofted_node_idx++;
    }
  }
  POSTCONDITION ( lofted_node_idx == number_lofted_nodes);

  // Now go back and flag all nodes that were detected to be tab nodes.
  // Note that it is currently possible for some nodes to be missed as
  // mixed nodes -- they are just flagged as contiuum nodes. This 
  // traps those nodes if they are tab nodes. Still need to work on 
  // correcting the erroneous flagging.
  ContactNode<Real>** node_list = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for( i=0 ; i<number_acme_nodes ; ++i){
    ContactNode<Real>* node = node_list[i];
    if (is_a_tab_node[node->HostGlobalArrayIndex()]){
      if ( node->Physical_Type() == ContactNode<Real>::SHELL_NODE){
	node->Physical_Type(ContactNode<Real>::SHELL_TAB_NODE);
#ifdef SHELL_DEBUG
	  postream << "*> node " << node->Exodus_ID()
		   << " :: Make that a SHELL TAB node.\n";
#endif
      } else if ( node->Physical_Type() == ContactNode<Real>::MIXED_NODE){
	node->Physical_Type(ContactNode<Real>::MIXED_TAB_NODE);
#ifdef SHELL_DEBUG
	  postream << "*> node " << node->Exodus_ID()
		   << " :: Make that a MIXED TAB node.\n";
#endif
      } else if ( node->Physical_Type() == ContactNode<Real>::CONTINUUM_NODE){
	node->Physical_Type(ContactNode<Real>::MIXED_TAB_NODE);
#ifdef SHELL_DEBUG
	  postream << "*> node " << node->Exodus_ID()
		   << " :: Fix CONTINUUM node to be a MIXED TAB node.\n";
#endif
      }
    }
  }
}

void ContactShellHandler::Loft_Nodes( int num_configs,
				      VariableHandle CURRENT_POSITION,
				      VariableHandle PREDICTED_POSITION,
				      VariableHandle AUGMENTED_POSITION,
				      VariableHandle FACE_NORMAL,
				      VariableHandle LOFTING_VECTOR )
{
//#define SHELL_DEBUG
#ifdef SHELL_DEBUG
  //  ContactParOStream postream(SearchComm);
#endif
	
  int number_of_nodes = topology->Number_of_Nodes();

#ifndef CONTACT_NO_MPI
  Real* sb = NULL;
  Real* rb = NULL;
  Real* rb0 = NULL;
  if( contact_number_of_processors( SearchComm ) > 1 ) {
    char* send_buf=NULL;
    char* recv_buf=NULL;
    // Import the ghost faces and normals (4 reals are stored, t then n)
    RequestHandle* recv_handles;
    topology->comm_buffer->Buffers( *import_loft_face_data_comm_plan,
				    4*sizeof(Real), 
				    &send_buf, &recv_buf, &recv_handles );
    sb = (Real*) send_buf;
    rb = (Real*) recv_buf;
    rb0 = (Real*) recv_buf;
    int offset = 0;
    for( int i=0 ; i<import_loft_face_data_comm_plan->Num_Export_Comm_Partners() ; ++i ){
      ContactNode<Real>** node_list = (ContactNode<Real>**) 
	import_loft_face_data_comm_plan->Export_Entity_List(i);
      ContactTopologyEntity<Real>* old_node = node_list[0];
      int next_face_index = 0;
      for( int j=0 ; j<import_loft_face_data_comm_plan->Num_Export_to_Proc(i) ; ++j ){
        PRECONDITION(!(node_list[j]->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)));
	if( old_node != node_list[j] ) 
	  next_face_index = 0;
	old_node = node_list[j];
        
	ContactFace<Real>* face = node_list[j]->GetFace(next_face_index++);
        while (face->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) {
	  face = node_list[j]->GetFace(next_face_index++);
        }
        
	PRECONDITION( face->FaceType() == ContactSearch::SHELLQUADFACEL4 ||
		      face->FaceType() == ContactSearch::SHELLTRIFACEL3 );
                      
	switch( face->FaceType() ){
	case(ContactSearch::SHELLQUADFACEL4):{
	  ContactShellQuadFaceL4<Real>* sq4 = (ContactShellQuadFaceL4<Real>*) face;
	  sb[offset++] = sq4->Lofting_Factor()*sq4->Thickness();
	  break;
        }
	case(ContactSearch::SHELLTRIFACEL3):{
	  ContactShellTriFaceL3<Real>* st3 = (ContactShellTriFaceL3<Real>*) face; 
	  sb[offset++] = st3->Lofting_Factor()*st3->Thickness();
	  break;
        }
        default:
          POSTCONDITION(0);
          break;
	}
	Real* fn = face->Variable(FACE_NORMAL);
	sb[offset++] = fn[0];
	sb[offset++] = fn[1];
	sb[offset++] = fn[2];
#ifdef SHELL_DEBUG
	Real dbg_t;
	switch( face->FaceType() ){
	case(ContactSearch::SHELLQUADFACEL4):{
	  ContactShellQuadFaceL4<Real>* sq4 = (ContactShellQuadFaceL4<Real>*) face;
	  dbg_t = sq4->Lofting_Factor()*sq4->Thickness();
	  break;
	}
	case(ContactSearch::SHELLTRIFACEL3):{
	  ContactShellTriFaceL3<Real>* st3 = (ContactShellTriFaceL3<Real>*) face; 
	  dbg_t = st3->Lofting_Factor()*st3->Thickness();
	  break;
	}
	}
	postream << "Packing for node " << node_list[j]->Exodus_ID()
                 << " face "<<face->Global_ID()
		 << " t=" << dbg_t << " & n=("
		 << fn[0] << "," << fn[1] << "," << fn[2] << ")\n";
#endif
      }
    }
#ifdef SHELL_DEBUG
    postream.flush();
#endif
    int size_per_object = 4*sizeof(Real);
    contact_communicate_packed_buffers( SearchComm, 
					*import_loft_face_data_comm_plan, 
					size_per_object,
					send_buf, recv_buf, recv_handles );
  }
#endif

  Real Physical_Face_Normal[3][3];
  Real T_Weighted_Normal[3][3];
  int num_faces_in_pf[3];

  ContactNode<Real>** Nodes = 
      reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  int i = -1;
  for( int ii=0 ; ii<number_of_nodes ; ++ii){
    ContactNode<Real>* node = Nodes[ii];
    if(  node->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) continue;
    ++i;
    PRECONDITION(i<topology->Number_of_Primary_Nodes());
    //std::cerr << "here in ShellHandler::Loft_Nodes, ii = " << ii << ", node->Is_a_Shell_Node() = " << node->Is_a_Shell_Node() << std::endl;
    if( !node->Is_a_Shell_Node() ) continue;
    // We have decided for the time being, to only loft shell nodes.
    // We don't loft mixed nodes because that would change the continuum
    //   faces which are already correct.
    // We don't loft tab nodes (either shell or mixed) because the lofting
    //   is ambiguous for this case.
#ifdef SHELL_DEBUG
    postream << "Node " << node->Exodus_ID() << " type = "<<node->Physical_Type()<<"\n";
#endif

    ContactSearch::Search_Option_Status simple_lofting; 
    Real dummy;
    topology->Search()->Get_Search_Option(ContactSearch::SHELL_SIMPLE_LOFTING,simple_lofting,&dummy);
    if (ContactSearch::INACTIVE == simple_lofting) {
     // here is the complex shell lofting strategy. This strategy involves
     // finding all the physical faces and using these to solve for the 
     // exact location of the lofted node. This is very accurate for
     // many corners, but has problems when the angles are very small and
     // when the number of physical faces changes. 
    if( node->Physical_Type() == ContactNode<Real>::SHELL_NODE ){
      Real loft[3];
      loft[0] = 0.0;
      loft[1] = 0.0;
      loft[2] = 0.0;
      Physical_Face_Normal[0][0] = 0.0;
      Physical_Face_Normal[0][1] = 0.0;
      Physical_Face_Normal[0][2] = 0.0;
      Physical_Face_Normal[1][0] = 0.0;
      Physical_Face_Normal[1][1] = 0.0;
      Physical_Face_Normal[1][2] = 0.0;
      Physical_Face_Normal[2][0] = 0.0;
      Physical_Face_Normal[2][1] = 0.0;
      Physical_Face_Normal[2][2] = 0.0;
      T_Weighted_Normal[0][0] = 0.0;
      T_Weighted_Normal[0][1] = 0.0;
      T_Weighted_Normal[0][2] = 0.0;
      T_Weighted_Normal[1][0] = 0.0;
      T_Weighted_Normal[1][1] = 0.0;
      T_Weighted_Normal[1][2] = 0.0;
      T_Weighted_Normal[2][0] = 0.0;
      T_Weighted_Normal[2][1] = 0.0;
      T_Weighted_Normal[2][2] = 0.0;
      num_faces_in_pf[0] = 0;
      num_faces_in_pf[1] = 0;
      num_faces_in_pf[2] = 0;
      int num_face_connections  = node->Number_Face_Connections();
      int num_physical_faces = 0;
      Real dot;
      Real* face_normal;
      Real thickness;
      Real max_thickness = 0.0;
      int loop_size = num_face_connections;
#ifndef CONTACT_NO_MPI
      int offset = offset_to_buffer_offsets[i];
      int* node_buffer_offsets = buffer_offsets + offset;
      loop_size += number_lofting_ghosts[i];
#endif
#ifdef SHELL_DEBUG
      postream << "Processing node " << node->Exodus_ID() << "\n";
#endif
      for( int j=0 ; j<loop_size ; ++j){
	bool process = false;
	if( j<num_face_connections ){
          if (topology->HaveGhosting()) {
            if (node->GetFace(j)->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) continue;
          }
	  /********************************************************************/
	  /*                                                                  */
	  /*    P R O C E S S    O N    P R O C E S S O R    F A C E S        */
	  /*                                                                  */
	  /********************************************************************/
	  ContactFace<Real>* face = node->GetFace(j);
          POSTCONDITION(face);
	  if( ContactSearch::Is_a_Shell_Face( face->FaceType() ) ){
	    face_normal = face->Variable(FACE_NORMAL);
	    PRECONDITION( face->FaceType() == ContactSearch::SHELLQUADFACEL4 ||
			  face->FaceType() == ContactSearch::SHELLTRIFACEL3 );
	    
	    switch(face->FaceType()){
	    case(ContactSearch::SHELLQUADFACEL4):{
	      ContactShellQuadFaceL4<Real>* sq4 = (ContactShellQuadFaceL4<Real>*) face;
	      thickness = sq4->Lofting_Factor()*sq4->Thickness();
	      break;
	    }
	    case(ContactSearch::SHELLTRIFACEL3):{
	      ContactShellTriFaceL3<Real>* st3 = (ContactShellTriFaceL3<Real>*) face;      
	      thickness = st3->Lofting_Factor()*st3->Thickness();
	      break;
	    }
            default:
              POSTCONDITION(0);
              break;
	    }
#ifdef SHELL_DEBUG
	    postream << "   On Proc:  t=" << thickness << " & n=(" 
		     << face_normal[0] << "," << face_normal[1] << ","
		     << face_normal[2] << ")\n";
#endif
	    if ( thickness > max_thickness) max_thickness = thickness;
	    if ( thickness > 0.0 ) process = true;
	  }
#ifndef CONTACT_NO_MPI
	} else {
	  /********************************************************************/
	  /*                                                                  */
	  /*    P R O C E S S    O F F   P R O C E S S O R    F A C E S       */
	  /*                                                                  */
	  /********************************************************************/
	  int jj = j - num_face_connections;
          PRECONDITION(number_lofting_ghosts[i]>0);
	  Real* node_data = rb + 4*node_buffer_offsets[jj];
	  thickness   = *node_data;
	  face_normal = node_data + 1;
#ifdef SHELL_DEBUG
	  postream << "   Off Proc: t=" << thickness << " & n=(" 
		   << face_normal[0] << "," << face_normal[1] << ","
		   << face_normal[2] << ")\n";
#endif
	  if ( thickness > max_thickness) max_thickness = thickness;
	  if( thickness > 0.0 ) process = true;
#endif
	}
#ifdef SHELL_DEBUG
	postream << "    *** max thickness: " << max_thickness << "\n";
#endif
	if( process ){
	  int  closest_physical_face = -1;
	  Real largest_dot_product   = -2.0;
	  // See if a physical face has been defined that fits this face
          int k;
	  for( k=0 ; k<num_physical_faces ; ++k){
	    dot = (Physical_Face_Normal[k][0]*face_normal[0] +
		   Physical_Face_Normal[k][1]*face_normal[1] +
		   Physical_Face_Normal[k][2]*face_normal[2] );
	    // Store the closest physical face
	    if( dot > largest_dot_product ){
	      largest_dot_product   = dot;
	      closest_physical_face = k;
	    }
	    // Assume sharp is 60 degrees
	    if( dot > 0.5 ) break;
	  }
	  if( k < num_physical_faces ){
	    // this normal is part of Physical Face k, add it to it
	    T_Weighted_Normal[k][0] += thickness*face_normal[0];
	    T_Weighted_Normal[k][1] += thickness*face_normal[1];
	    T_Weighted_Normal[k][2] += thickness*face_normal[2];
	    Physical_Face_Normal[k][0] *= num_faces_in_pf[k];
	    Physical_Face_Normal[k][1] *= num_faces_in_pf[k];
	    Physical_Face_Normal[k][2] *= num_faces_in_pf[k];
	    Physical_Face_Normal[k][0] += face_normal[0];
	    Physical_Face_Normal[k][1] += face_normal[1];
	    Physical_Face_Normal[k][2] += face_normal[2];
	    Normalize(Physical_Face_Normal[k]);
	    num_faces_in_pf[k] += 1;
	  } else {
	    if( num_physical_faces < 3 ){
	      // Add a new physical face
	      T_Weighted_Normal[num_physical_faces][0]=thickness*face_normal[0];
	      T_Weighted_Normal[num_physical_faces][1]=thickness*face_normal[1];
	      T_Weighted_Normal[num_physical_faces][2]=thickness*face_normal[2];
	      Physical_Face_Normal[num_physical_faces][0] = face_normal[0];
	      Physical_Face_Normal[num_physical_faces][1] = face_normal[1];
	      Physical_Face_Normal[num_physical_faces][2] = face_normal[2];
	      num_faces_in_pf[num_physical_faces] = 1;
	      num_physical_faces++;
	    } else {
	      // We already have three physical faces, add this to the 
	      // closest one
              PRECONDITION(closest_physical_face>=0);
	      T_Weighted_Normal[closest_physical_face][0] += 
		thickness*face_normal[0];
	      T_Weighted_Normal[closest_physical_face][1] += 
		thickness*face_normal[1];
	      T_Weighted_Normal[closest_physical_face][2] += 
		thickness*face_normal[2];
	      Physical_Face_Normal[closest_physical_face][0] *= 
		num_faces_in_pf[closest_physical_face];
	      Physical_Face_Normal[closest_physical_face][1] *= 
		num_faces_in_pf[closest_physical_face];
	      Physical_Face_Normal[closest_physical_face][2] *= 
		num_faces_in_pf[closest_physical_face];
	      Physical_Face_Normal[closest_physical_face][0] += face_normal[0];
	      Physical_Face_Normal[closest_physical_face][1] += face_normal[1];
	      Physical_Face_Normal[closest_physical_face][2] += face_normal[2];

	      Normalize(Physical_Face_Normal[closest_physical_face]);

	      num_faces_in_pf[closest_physical_face] += 1;
	    }
	  }
	}
      }
      
      /********************************************************************/
      /*                                                                  */
      /*   S O L V E    F O R   T H E    L O F T E D    P O S I T I O N   */
      /*                                                                  */
      /********************************************************************/
      
      Real inverse[3][3];
      Real mag_twn_0, mag_twn_1, mag_twn_2;
      Real rhs[3];
      Real dot_twns;
      switch( num_physical_faces ){
      case(1):
	loft[0] = T_Weighted_Normal[0][0] / num_faces_in_pf[0];
	loft[1] = T_Weighted_Normal[0][1] / num_faces_in_pf[0];
	loft[2] = T_Weighted_Normal[0][2] / num_faces_in_pf[0];
#ifdef SHELL_DEBUG
	postream << "   Lofted Vector (1 Physical Face) = " << loft[0] << " "
		 << loft[1] << " " << loft[2] << "\n";
#endif
	break;
      case(2):
	mag_twn_0 = ( T_Weighted_Normal[0][0]*T_Weighted_Normal[0][0] +
		      T_Weighted_Normal[0][1]*T_Weighted_Normal[0][1] +
		      T_Weighted_Normal[0][2]*T_Weighted_Normal[0][2] ) /
	  num_faces_in_pf[0];
	rhs[0] = mag_twn_0;
	mag_twn_1 = ( T_Weighted_Normal[1][0]*T_Weighted_Normal[1][0] +
		      T_Weighted_Normal[1][1]*T_Weighted_Normal[1][1] +
		      T_Weighted_Normal[1][2]*T_Weighted_Normal[1][2] ) /
	  num_faces_in_pf[1];
	rhs[1] = mag_twn_1;
	rhs[2] = 0.0;
	
	// if twn_0 and twn_1 are coincident, then don't loft at all.
	// Could do something with detecting edge directions, but not sure 
	// how to do this at this point.
	dot_twns = 
	  ( T_Weighted_Normal[0][0]*T_Weighted_Normal[1][0] +
	    T_Weighted_Normal[0][1]*T_Weighted_Normal[1][1] +
	    T_Weighted_Normal[0][2]*T_Weighted_Normal[1][2] )
	  / (std::sqrt(mag_twn_0 * num_faces_in_pf[0]) * 
	     std::sqrt(mag_twn_1 * num_faces_in_pf[1]) );
	if ( std::fabs ( 1.0 - std::fabs( dot_twns ) ) < 1.0e-6) {
	  loft[0] = 0.0;
	  loft[1] = 0.0;
	  loft[2] = 0.0;
#ifdef SHELL_DEBUG
	  postream << "  NOTE: thickness directions are coincident! \n";
	  postream << "  T_Weighted_Normal : " 
		   << T_Weighted_Normal[0][0] << " "
		   << T_Weighted_Normal[0][1] << " "
		   << T_Weighted_Normal[0][2] << "\n "
		   << "                    : " 
		   << T_Weighted_Normal[1][0] << " "
		   << T_Weighted_Normal[1][1] << " "
		   << T_Weighted_Normal[1][2] << "\n ";
	  postream << "  rhs               : " 
		   << rhs[0] << " "
		   << rhs[1] << " "
		   << rhs[2] << "\n ";
#endif
	  break;
	}

	// For this case, the third equation is T * (n1 x n2) = 0;
	T_Weighted_Normal[2][0] = 
	  T_Weighted_Normal[0][1]*T_Weighted_Normal[1][2] -
	  T_Weighted_Normal[1][1]*T_Weighted_Normal[0][2] ;
	T_Weighted_Normal[2][1] = 
	  T_Weighted_Normal[1][0]*T_Weighted_Normal[0][2] -
	  T_Weighted_Normal[0][0]*T_Weighted_Normal[1][2] ;
	T_Weighted_Normal[2][2] = 
	  T_Weighted_Normal[0][0]*T_Weighted_Normal[1][1] -
	  T_Weighted_Normal[1][0]*T_Weighted_Normal[0][1] ;
#ifdef SHELL_DEBUG
	postream << "  T_Weighted_Normal : " 
		 << T_Weighted_Normal[0][0] << " "
		 << T_Weighted_Normal[0][1] << " "
		 << T_Weighted_Normal[0][2] << "\n "
		 << "                    : " 
		 << T_Weighted_Normal[1][0] << " "
		 << T_Weighted_Normal[1][1] << " "
		 << T_Weighted_Normal[1][2] << "\n "
		 << "                    : " 
		 << T_Weighted_Normal[2][0] << " "
		 << T_Weighted_Normal[2][1] << " "
		 << T_Weighted_Normal[2][2] << "\n ";
	postream << "  rhs               : " 
		 << rhs[0] << " "
		 << rhs[1] << " "
		 << rhs[2] << "\n ";
#endif
	Invert_3x3Matrix( T_Weighted_Normal, inverse );
#ifdef SHELL_DEBUG
	postream << "  inverse           : " 
		 << inverse[0][0] << " "
		 << inverse[0][1] << " "
		 << inverse[0][2] << "\n "
		 << "                    : " 
		 << inverse[1][0] << " "
		 << inverse[1][1] << " "
		 << inverse[1][2] << "\n "
		 << "                    : " 
		 << inverse[2][0] << " "
		 << inverse[2][1] << " "
		 << inverse[2][2] << "\n ";
#endif
	loft[0] = inverse[0][0]*rhs[0] + inverse[0][1]*rhs[1] + inverse[0][2]*rhs[2];
	loft[1] = inverse[1][0]*rhs[0] + inverse[1][1]*rhs[1] + inverse[1][2]*rhs[2];
	loft[2] = inverse[2][0]*rhs[0] + inverse[2][1]*rhs[1] + inverse[2][2]*rhs[2];
#ifdef SHELL_DEBUG
	postream << "   Lofted Vector (2 Physical Face) = " << loft[0] << " "
		 << loft[1] << " " << loft[2] << "\n";
#endif
	break;
      case(3):
	Invert_3x3Matrix( T_Weighted_Normal, inverse );
	mag_twn_0 = ( T_Weighted_Normal[0][0]*T_Weighted_Normal[0][0] +
		      T_Weighted_Normal[0][1]*T_Weighted_Normal[0][1] +
		      T_Weighted_Normal[0][2]*T_Weighted_Normal[0][2] ) /
	  num_faces_in_pf[0];
	rhs[0] = mag_twn_0;
	mag_twn_1 = ( T_Weighted_Normal[1][0]*T_Weighted_Normal[1][0] +
		      T_Weighted_Normal[1][1]*T_Weighted_Normal[1][1] +
		      T_Weighted_Normal[1][2]*T_Weighted_Normal[1][2] ) /
	  num_faces_in_pf[1];
	rhs[1] = mag_twn_1;
	mag_twn_2 = ( T_Weighted_Normal[2][0]*T_Weighted_Normal[2][0] +
		      T_Weighted_Normal[2][1]*T_Weighted_Normal[2][1] +
		      T_Weighted_Normal[2][2]*T_Weighted_Normal[2][2] ) 
	  / num_faces_in_pf[2];
	rhs[2] = mag_twn_2;
	loft[0] = inverse[0][0]*rhs[0] + inverse[0][1]*rhs[1] + inverse[0][2]*rhs[2];
	loft[1] = inverse[1][0]*rhs[0] + inverse[1][1]*rhs[1] + inverse[1][2]*rhs[2];
	loft[2] = inverse[2][0]*rhs[0] + inverse[2][1]*rhs[1] + inverse[2][2]*rhs[2];
#ifdef SHELL_DEBUG
	postream << "   Lofted Vector (3 Physical Face) = " << loft[0] << " "
		 << loft[1] << " " << loft[2] << "\n";
#endif
	break;
      }
      // if loft vector is too large, cut down its magnitude.
      Real mag_loft = Magnitude(loft);

      if ( max_thickness > 0.0) {
	if ( mag_loft/max_thickness > 5.0) {
#ifdef SHELL_DEBUG
	  postream << "   Lofted Vector too large, limited from = " << loft[0]
		   << " " << loft[1] << " " << loft[2] << "\n";
#endif
	  loft[0] *= 5.0*max_thickness/mag_loft;
	  loft[1] *= 5.0*max_thickness/mag_loft;
	  loft[2] *= 5.0*max_thickness/mag_loft;
#ifdef SHELL_DEBUG
	  postream << "                                    to = " << loft[0]
		   << " " << loft[1] << " " << loft[2] << "\n";
#endif
	}
      }

      Real* lofting_vector = node->Variable(LOFTING_VECTOR);
      lofting_vector[0] = loft[0];
      lofting_vector[1] = loft[1];
      lofting_vector[2] = loft[2];
    }
    } else {
      // this is a much simplified lofting strategy. The idea is to
      // take the average normal at the node, then multiply that times
      // the average thickness. This should produce something reasonable
      // in corners, and not change when the number of physical faces 
      // changes.
      // I think this may fall apart if you have many many thin angled 
      // shells on one plane and  just one shell on another plane. The
      // many many shells will have their normals weighted more heavily
      // than the big face. But it my turn out o.k. after all.
    if( node->Physical_Type() == ContactNode<Real>::SHELL_NODE ){
      Real avg_norm[3];
      avg_norm[0] = 0.0;
      avg_norm[1] = 0.0;
      avg_norm[2] = 0.0;
      Real avg_thick = 0.0;
      int num_face_connections = node->Number_Face_Connections();
      int num_faces = num_face_connections;
#ifndef CONTACT_NO_MPI
      int offset = offset_to_buffer_offsets[i];
      int* node_buffer_offsets = buffer_offsets + offset;
      num_faces += number_lofting_ghosts[i];
#endif
#ifdef SHELL_DEBUG
      postream << "Processing node " << node->Exodus_ID() << "\n";
#endif
      for( int j=0 ; j<num_faces ; ++j){
        Real* face_normal = NULL;
        Real thickness;
	if( j<num_face_connections ){
          if (topology->HaveGhosting()) {
            if (node->GetFace(j)->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) continue;
          }
	  /********************************************************************/
	  /*                                                                  */
	  /*    P R O C E S S    O N    P R O C E S S O R    F A C E S        */
	  /*                                                                  */
	  /********************************************************************/
	  ContactFace<Real>* face = node->GetFace(j);
	  if( ContactSearch::Is_a_Shell_Face( face->FaceType() ) ){
	    face_normal = face->Variable(FACE_NORMAL);
	    PRECONDITION( face->FaceType() == ContactSearch::SHELLQUADFACEL4 ||
			  face->FaceType() == ContactSearch::SHELLTRIFACEL3 );
	    
	    switch(face->FaceType()){
	    case(ContactSearch::SHELLQUADFACEL4):{
	      ContactShellQuadFaceL4<Real>* sq4 = (ContactShellQuadFaceL4<Real>*) face;
	      thickness = sq4->Lofting_Factor()*sq4->Thickness();
	      break;
	    }
	    case(ContactSearch::SHELLTRIFACEL3):{
	      ContactShellTriFaceL3<Real>* st3 = (ContactShellTriFaceL3<Real>*) face;      
	      thickness = st3->Lofting_Factor()*st3->Thickness();
	      break;
	    }
            default:
              POSTCONDITION(0);
              break;
	    }
#ifdef SHELL_DEBUG
	    postream << "   On Proc:  t=" << thickness << " & n=(" 
		     << face_normal[0] << "," << face_normal[1] << ","
		     << face_normal[2] << ")\n";
#endif
	  }
#ifndef CONTACT_NO_MPI
	} else {
	  /********************************************************************/
	  /*                                                                  */
	  /*    P R O C E S S    O F F   P R O C E S S O R    F A C E S       */
	  /*                                                                  */
	  /********************************************************************/
	  int jj = j - num_face_connections;
	  Real* node_data = rb + 4*node_buffer_offsets[jj];
	  thickness   = *node_data;
	  face_normal = node_data + 1;
#ifdef SHELL_DEBUG
	  postream << "   Off Proc: t=" << thickness << " & n=(" 
		   << face_normal[0] << "," << face_normal[1] << ","
		   << face_normal[2] << ")\n";
#endif
#endif
	}
        avg_norm[0] += face_normal[0];
        avg_norm[1] += face_normal[1];
        avg_norm[2] += face_normal[2];
        avg_thick += thickness;
      }
      // normalize normals and multiply by avg thickness
      avg_thick /= num_faces;
      Normalize(avg_norm);
#ifdef SHELL_DEBUG
      postream << "  average  : t=" << avg_thick << " & n=(" 
               << avg_norm[0] << "," << avg_norm[1] << ","
	       << avg_norm[2] << ")\n";
#endif

      Real* lofting_vector = node->Variable(LOFTING_VECTOR);
      lofting_vector[0] = avg_norm[0] * avg_thick;
      lofting_vector[1] = avg_norm[1] * avg_thick;
      lofting_vector[2] = avg_norm[2] * avg_thick;
    } // end loop on faces
  } // end of simplified lofting algorithm
  } // end loop on nodes
    
#ifdef SHELL_DEBUG
    postream.flush();
#endif

#ifndef CONTACT_NO_MPI
  // Now send the "lofting vector" to the ghosts
  contact_import_reg_data( SearchComm, *comm_plan_to_ghost,
			   *topology->comm_buffer, LOFTING_VECTOR, 3);
#endif

  // Now loft the nodes
  i = -1;
  for( int ii=0 ; ii<number_of_nodes ; ++ii){
    ContactNode<Real>* node = Nodes[ii];
    if(  node->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) continue;
    ++i;
    PRECONDITION(i<topology->Number_of_Primary_Nodes());
    if( !node->Is_a_Shell_Node() ) continue;
    Real* pre_loft = 
      (static_cast<ContactShellNode*>(node))->Previous_Lofting();
    Real* pos = node->Variable(CURRENT_POSITION);
    Real* loft = node->Variable(LOFTING_VECTOR);
#ifdef SHELL_DEBUG
    postream << "Lofting Node " << node->Exodus_ID() << " ("
	     << loft[0] << "," << loft[1] << "," << loft[2] << ")\n";
#endif
    if( num_configs > 1 ){
      if( lofting_computed ){
	pos[0] += pre_loft[0];
	pos[1] += pre_loft[1];
	pos[2] += pre_loft[2];
      } else {
	pos[0] += loft[0];
	pos[1] += loft[1];
	pos[2] += loft[2];
      }
      pos = node->Variable(PREDICTED_POSITION);
      pos[0] += loft[0];
      pos[1] += loft[1];
      pos[2] += loft[2];
      if( num_configs > 2 ){
	pos = node->Variable(AUGMENTED_POSITION);
	pos[0] += loft[0];
	pos[1] += loft[1];
	pos[2] += loft[2];
      }
    } else {
      pos[0] += loft[0];
      pos[1] += loft[1];
      pos[2] += loft[2];
    }
    pre_loft[0] = loft[0];
    pre_loft[1] = loft[1];
    pre_loft[2] = loft[2];
  }
#ifdef SHELL_DEBUG
  i = -1;
  for( int ii=0 ; ii<number_of_nodes ; ++ii){
    ContactNode<Real>* node = Nodes[ii];
    if(  node->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) continue;
    ++i;
    PRECONDITION(i<topology->Number_of_Primary_Nodes());
    if( !node->Is_a_Shell_Node() ) continue
    Real* pos = node->Variable(CURRENT_POSITION);
    postream << "Position for Node " << node->Exodus_ID() << " ("
	     << pos[0] << "," << pos[1] << "," << pos[2] << ")\n";
  }
  postream.flush();
#endif
  lofting_computed = true;
//#undef SHELL_DEBUG
}

void
ContactShellHandler::Acme_NodeGID_for_Host_Node( int node, int num,
                                                 int& hi, int& lo ) {
    PRECONDITION (node < number_host_code_nodes && node >= 0);
    PRECONDITION (num < Num_Acme_Nodes_for_Host_Node(node) && num >= 0);
    int offset = offset_for_node[node];
    PRECONDITION(offset + num < number_acme_nodes);
    int index = host_to_acme_node_map[offset + num];
    hi = new_global_node_ids[2*index+0];
    lo = new_global_node_ids[2*index+1];
  };


void ContactShellHandler::Initialize_Node_Maps()
{
  number_acme_nodes = 0;
  num_acme_nodes_for_host_code_node = new int[number_host_code_nodes];
  std::memset(num_acme_nodes_for_host_code_node, 0, number_host_code_nodes*sizeof(int));
  allocator = new ContactFixedSizeAllocator( sizeof( created_node ), 100 );
  first_node = NULL;
  last_node = NULL;
}

void ContactShellHandler::Create_Node( int host_id, int acme_id )
{
  number_acme_nodes++;
  num_acme_nodes_for_host_code_node[host_id]++;
  created_node* new_node = (created_node*) allocator->New_Frag();
  new_node->host_id = host_id;
  new_node->acme_id = acme_id;
  new_node->next_node = NULL;
  if( first_node ){
    last_node->next_node = new_node;
    last_node = new_node;
  } else {
    first_node = new_node;
    last_node = first_node;
  }
}

void ContactShellHandler::Complete_Node_Maps()
{
  int i,j;

  // First build the offsets and find max number of acme nodes for a
  // host code node
  offset_for_node = new int[number_host_code_nodes];
  offset_for_node[0] = 0;
  max_num_acme_nodes_for_host_code_node = 0;
  for( i=1 ; i<number_host_code_nodes ; ++i){
    offset_for_node[i] = offset_for_node[i-1] + 
                          num_acme_nodes_for_host_code_node[i-1];
    max_num_acme_nodes_for_host_code_node =
      std::max(num_acme_nodes_for_host_code_node[i],
               max_num_acme_nodes_for_host_code_node);
  }
  int parallel_max =
    contact_global_maximum(max_num_acme_nodes_for_host_code_node,SearchComm );
  max_num_acme_nodes_for_host_code_node = parallel_max;
 
  host_to_acme_node_map = new int[number_acme_nodes];
  for( i=0 ; i<number_acme_nodes ; ++i) 
    host_to_acme_node_map[i] = -1;

  // Now loop over all the created nodes and build the host_to_acme node map
  created_node* node = first_node;  
  for( i=0 ; i<number_acme_nodes ; ++i){
    int host_id = node->host_id;
    int acme_id = node->acme_id;
    int offset = offset_for_node[host_id];
    int* acme_ids = host_to_acme_node_map + offset;
    for( j=0 ; j<num_acme_nodes_for_host_code_node[host_id] ; ++j){
      if( acme_ids[j] == -1 ){
	acme_ids[j] = acme_id;
	break;
      }
    }
    POSTCONDITION( j < num_acme_nodes_for_host_code_node[host_id] );
    node = node->next_node;
  }

  // Delete all the temporary memory we used to build these data structures
  node = first_node;
  created_node* next_node(NULL);
  for( i=0 ; i<number_acme_nodes ; ++i){
    next_node = node->next_node;
    allocator->Delete_Frag(node);
    node = next_node;
  }
  POSTCONDITION( next_node == NULL );
  delete allocator;

#ifdef SHELL_DEBUG
  postream << my_proc_id << ":HEY!!!! Here is host_to_acme_node_map." << "\n";
  for ( i = 0; i < number_host_code_nodes ; ++i) { 
    postream << my_proc_id << ":   host code id (c#): " << i << " " 
	     << " shell nodes: ";
    for( j=0 ; j<num_acme_nodes_for_host_code_node[i] ; ++j)
      postream << Acme_Node_for_Host_Node( i, j ) << " ";
    postream << "\n";
  }
#endif

}    

void 
ContactShellHandler::Renumber_Nodes(int number_node_blocks,
				    const int* host_number_nodes_per_block,
				    int* number_nodes_per_block,
				    const int& number_face_blocks,
               const ContactSearch::ContactFace_Type* face_block_types,
				    const int* number_faces_per_block,
				    int* face_connectivity )
{
  // When we created all of the "lofted" nodes, we did it assuming a
  // single node block.  This isn't necessarily the case.  Here we
  // renumber the node lists and modify the connectivities to account
  // for multiple node blocks.  We could not have done this in one pass
  // because we can't set the ids until we know how many nodes were 
  // created in node block 1, node block 2, etc. which of course we don't
  // know until everything is done.
  //
  // The steps for sorting this out are
  //   1) Count the total number of nodes for each node block
  //   2) Build a 'next index' array where we can quickly get a new index
  //      for the reordering.
  //   3) Create a map from the old numbering (the one we just created) to
  //      the final numbering where everything is consistent with the node
  //      blocks.
  //   4) Reset the connectivity array based on the map created in step 3.

  int i,j;
  created_node* node;

#ifdef SHELL_DEBUG
  node = first_node;
  for( i=0 ; i<number_acme_nodes ; ++i){
    postream << "Host Node " << node->host_id << " created acme node "
	     << node->acme_id << "\n";
    node = node->next_node;
    
  }
  postream.flush();
#endif

  block_offset_for_nodes = new int[number_node_blocks];

  if( number_node_blocks == 1 ){
    number_nodes_per_block[0] = number_acme_nodes;
    block_offset_for_nodes[0] = 0;
    return;
  } else {
    block_offset_for_nodes[0] = 0;
    for( i=1 ; i<number_node_blocks ; ++i)
      block_offset_for_nodes[i] = block_offset_for_nodes[i-1] + 
	host_number_nodes_per_block[i-1];
  }
  
  //
  //=============
  // Step 1
  //=============
  //

  for( i=0 ; i<number_node_blocks ; ++i)
    number_nodes_per_block[i] = 0;
  
  // block_index will be the ending index for a node block.
  int* block_index = new int[number_node_blocks];
  block_index[0] = host_number_nodes_per_block[0]-1;
  for( i=1 ; i<number_node_blocks ; ++i)
    block_index[i] = block_index[i-1] + host_number_nodes_per_block[i];

  node = first_node;
  for( i=0 ; i<number_acme_nodes ; ++i){
    int host_id = node->host_id;
    for( j=0 ; j<number_node_blocks ; ++j)
      if( host_id <= block_index[j] ) break;
    number_nodes_per_block[j] += 1;
    node = node->next_node;
  }

  // 
  //=============
  // Step 2
  //=============
  //
  int* next_index = new int[number_node_blocks];
  next_index[0] = 0;
  for( j=1 ; j<number_node_blocks ; ++j)
    next_index[j] = next_index[j-1] + number_nodes_per_block[j-1];
  
  // 
  //=============
  // Step 3
  //=============
  //
  int* map_created_to_final = new int[number_acme_nodes];
  node = first_node;
  for( i=0 ; i<number_acme_nodes ; ++i){
    int host_id = node->host_id;
    for( j=0 ; j<number_node_blocks ; ++j)
      if( host_id <= block_index[j] ) break;
    map_created_to_final[node->acme_id] = next_index[j];
    node->acme_id = next_index[j];
    next_index[j] += 1;
    node = node->next_node;
  }

#ifdef SHELL_DEBUG
  for( i=0 ; i<number_acme_nodes ; ++i){
    postream << "Created Node " << i << " will be final node " 
	 << map_created_to_final[i] << "\n";
  }
  postream.flush();
#endif

  // 
  //=============
  // Step 4
  //=============
  //
  int conn_offset = 0;
  for( i=0 ; i<number_face_blocks ; ++i){
    int nodes_per_face=0;
    switch( face_block_types[i] ){
    case( ContactSearch::QUADFACEL4 ):
    case( ContactSearch::SHELLQUADFACEL4 ):
      nodes_per_face = 4;
      break;
    case( ContactSearch::QUADFACEQ8 ):
      nodes_per_face = 8;
      break;
    case( ContactSearch::QUADFACEQ9 ):
      nodes_per_face = 9;
      break;
    case( ContactSearch::TRIFACEL3 ):
    case( ContactSearch::SHELLTRIFACEL3 ):
      nodes_per_face = 3;
      break;
    case (ContactSearch::TRIFACEQ6 ):
      nodes_per_face = 6;
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    POSTCONDITION(nodes_per_face>0);
    int num_node_conn_for_block = nodes_per_face*number_faces_per_block[i];
    for( j=0 ; j<num_node_conn_for_block ; ++j){
      int old_conn = face_connectivity[conn_offset] - 1;
      int new_conn;
      if ( old_conn < 0 ) new_conn = 0;
      else new_conn = map_created_to_final[old_conn] + 1;
      face_connectivity[conn_offset] = new_conn;
      conn_offset++;
    }
  }

  delete [] block_index;
  delete [] next_index;
  delete [] map_created_to_final;
}
