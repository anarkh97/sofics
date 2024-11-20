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


#include <algorithm>

#include "contact_assert.h"
#include "contact_sorting.h"
#include "ContactTopology.h"
#include "ContactAsymComm.h"
#include "ContactSymComm.h"
#include "ContactErrors.h"
#include "ContactTopologyEntityHash.h"
#include "ContactEntityDataHash.h"
#include "ContactNode.h"
#include "ContactNodeBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactElementBlock.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactShellNode.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactHexElementL8.h"
#include "ContactAnalyticSurface.h"
#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "CString.h"
#include "search_methods.h"
#include "Contact_Communication.h"
#include "ContactCommBuffer.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSequentialAllocator.h"
#include "ContactShellHandler.h"
#include "ContactSearchData.h"
#include "ContactBoundingBox.h"
#include "ContactBoundingBoxHierarchy.h"
#include "ContactBoundingBoxHierarchy_Int.h"
#include "ContactRangeSearch.h"
#include "contact_tolerances.h"

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifndef CONTACT_NO_MPI
//temporary include
#include "mpi.h"
#include "zoltan.h"
#include "ContactZoltan.h"
#include "ContactZoltanComm.h"
#include "ContactZoltanCommUtils.h"
#endif

using namespace std;

ContactTopology::ContactTopology( ContactErrors* Errors,
				  int Dimensionality, 
				  int number_analytic_surfaces,
				  int number_node_blocks,
	     const ContactSearch::ContactNode_Type* node_block_types, 
				  const int* host_number_nodes_per_block,
				  const int* host_exodus_node_ids,
				  const int* host_global_node_ids,
				  const double* coords,
				  int number_face_blocks, 
	     const ContactSearch::ContactFace_Type* face_block_types,
				  const int* number_faces_per_block, 
				  const int* global_face_ids,
				  const int* host_face_connectivity,
				  const Real* face_lofting_factors,
				  int number_element_blocks, 
	     const ContactSearch::ContactElement_Type* element_block_types,
				  const int* number_elements_per_block, 
				  const int* global_element_ids,
				  const int* element_connectivity,
				  int number_comm_partners, 
				  const int* comm_proc_ids,
				  const int* host_number_nodes_to_partner, 
				  const int* host_comm_nodes,
				  ContactCommBuffer* Comm_buffer,
				  MPI_Comm communicator,
                                  ContactSearch* srch,
				  ContactSearch::ContactErrorCode& error_code )
  :  computed_characteristic_length( false ),
     min_characteristic_length(BIGNUM),
     ghosted_node_blocks(NULL),
     ghosted_edge_blocks(NULL),
     ghosted_face_blocks(NULL),
     ghosted_element_blocks(NULL),
     have_global_ghosting(false),
     have_local_ghosting(false),
     have_tied_ghosting(false)
{
  PRECONDITION( number_node_blocks >= 0 );
  PRECONDITION( number_face_blocks+number_element_blocks >= 0 );
  int i,k;

  topology_type = PRIMARY;
  error_code = ContactSearch::NO_ERROR;
  no_parallel_consistency = ContactSearch::INACTIVE;
  SearchComm = communicator;
  comm_buffer = Comm_buffer;
  errors = Errors;
  search = srch;
  dimensionality = Dimensionality;
  node_block_ids = NULL;
  number_of_nodes  = 0;
  number_of_edges  = 0;
  number_of_faces  = 0;
  number_of_elements = 0;
  number_of_node_blocks = number_node_blocks;
  number_of_edge_blocks = 0;
  number_of_face_blocks = number_face_blocks;
  number_of_element_blocks = number_element_blocks;
  number_of_analytic_surfaces = number_analytic_surfaces;
  number_of_added_analytic_surfaces = 0;
  number_of_comm_partners = number_comm_partners;
  two_configurations = false;
  node_blocks = NULL;
  edge_blocks = NULL;
  face_blocks = NULL;
  if (number_analytic_surfaces>0) {
    AnalyticSurfaces = new ContactAnalyticSurface*[number_of_analytic_surfaces];
    for (i=0; i<number_of_analytic_surfaces; ++i) {
      AnalyticSurfaces[i] = NULL;
    }
  } else {
    AnalyticSurfaces = NULL;
  }
#ifndef CONTACT_NO_MPI
  Node_AsymComm = NULL;
  Node_SymComm  = NULL;
  Edge_SymComm  = NULL;
  GhostingCommSpec = NULL;
  GhostFaces_ZoltanComm = NULL;
  GhostOthers_ZoltanComm = NULL;
  TiedCommSpec = NULL;
  TiedFaces_ZoltanComm = NULL;
  TiedNodes_ZoltanComm = NULL;
#endif
  node_list = new ContactTopologyEntityList();
  edge_list = new ContactTopologyEntityList();
  face_list = new ContactTopologyEntityList();
  elem_list = new ContactTopologyEntityList();
  primary_node_list = new ContactTopologyEntityList();
  primary_edge_list = edge_list;
  primary_face_list = new ContactTopologyEntityList();
  primary_elem_list = new ContactTopologyEntityList();
#ifdef CONTACT_DEBUG_NODE
  number_debug_nodes = 0;
  debug_node_global_ids = NULL;
  debug_node_exodus_ids = NULL;
  debug_nodes = NULL;
#endif
  scratch_set = false;
  have_shells = false;
  old_to_new_node_map = NULL;
  
#ifndef CONTACT_NO_MPI
  ghosted_import_gids = NULL;
  ghosted_import_lids = NULL;
  ghosted_import_pids = NULL;
  num_ghosted_import  = 0;
  ghosted_export_gids = NULL;
  ghosted_export_lids = NULL;
  ghosted_export_pids = NULL;
  num_ghosted_export  = 0;
  tied_import_gids = NULL;
  tied_import_lids = NULL;
  tied_import_pids = NULL;
  num_tied_import  = 0;
  tied_export_gids = NULL;
  tied_export_lids = NULL;
  tied_export_pids = NULL;
  num_tied_export  = 0;
#endif

  Set_Up_Variable_Handles();

  // check if there are any shell faces passed in  
  for ( int blk = 0; blk < number_face_blocks; ++blk){
    if ( ContactSearch::Is_a_Shell_Face(face_block_types[blk]) ) {
      have_shells = true;
      break;
    }
  }
  
  int * number_nodes_per_block(NULL);
  int * exodus_node_ids(NULL);
  int * global_node_ids(NULL);
  int * face_connectivity(NULL);
  int * number_nodes_to_partner(NULL);
  int * comm_nodes(NULL);
  ContactType* node_entity_types(NULL);

  if (have_shells && Dimensionality == 3){
    shell_handler = new ContactShellHandler(errors, 
					    Dimensionality,
					    number_node_blocks,
					    node_block_types,
					    host_number_nodes_per_block,
					    number_nodes_per_block,
					    host_exodus_node_ids,
					    exodus_node_ids,
					    host_global_node_ids,
					    global_node_ids,
					    &node_entity_types,
					    coords,
					    number_face_blocks, 
					    face_block_types,
					    number_faces_per_block,
					    global_face_ids, 
					    host_face_connectivity,
					    face_connectivity,
					    face_lofting_factors,
					    number_comm_partners, 
					    comm_proc_ids, 
					    host_number_nodes_to_partner,
					    number_nodes_to_partner,
					    host_comm_nodes,
					    comm_nodes,
					    communicator,
					    srch,
					    this,
					    Comm_buffer,
					    ContactShellHandler::SEQUENTIAL_IDS,
					    NULL,
					    error_code );

    if( error_code != ContactSearch::NO_ERROR ) return;
  } else {
    shell_handler = NULL;

    // need to copy the host code data in order to pass it around. This
    // is all related to the pain of "const".
    int number_faces = 0;
    for ( i = 0; i < number_face_blocks; ++i) 
      number_faces += number_faces_per_block[i];
    
    int number_host_code_nodes = 0;
    for ( i = 0; i < number_node_blocks; ++i) 
      number_host_code_nodes += host_number_nodes_per_block[i];
    
    number_nodes_per_block = new int[number_node_blocks];
    for (i = 0; i < number_node_blocks; ++i){
      number_nodes_per_block[i] = host_number_nodes_per_block[i];
    }
    
    exodus_node_ids = new int[number_host_code_nodes];
    global_node_ids = new int[2*number_host_code_nodes];
    for ( i = 0; i < number_host_code_nodes; ++i){
      exodus_node_ids[i]       = host_exodus_node_ids[i];
      global_node_ids[2*i]     = host_global_node_ids[2*i];
      global_node_ids[2*i + 1] = host_global_node_ids[2*i + 1];
    }
	
    int num_nodes_in_conn = 0;
    for ( i = 0; i < number_face_blocks; ++i) {
      int num_face_nodes = 
	search->Number_Nodes_Per_Face(face_block_types[i]);
      POSTCONDITION(num_face_nodes>0);
      num_nodes_in_conn += num_face_nodes*number_faces_per_block[i];
    }
    face_connectivity = new int[num_nodes_in_conn];
    for ( i = 0; i < num_nodes_in_conn; ++i ) {
      face_connectivity[i] = host_face_connectivity[i];
    }
    
    number_nodes_to_partner = new int[number_comm_partners];
    int number_comm_nodes = 0;
    for ( i = 0; i < number_comm_partners; ++i ){
      number_nodes_to_partner[i] = host_number_nodes_to_partner[i];
      number_comm_nodes += number_nodes_to_partner[i];
    }
    comm_nodes = new int[number_comm_nodes];
    for ( i = 0; i < number_comm_nodes; ++i) {
      comm_nodes[i] = host_comm_nodes[i];
    }
    node_entity_types =  new ContactType[number_host_code_nodes];
    for( i=0 ; i<number_host_code_nodes ; ++i )
      node_entity_types[i] = CT_NODE;
  }
  
  ContactParOStream& postream = search->ParOStream();
  postream.flush();
  int number_host_code_nodes=0;
  for ( i = 0; i < number_of_node_blocks; ++i){
    number_host_code_nodes += number_nodes_per_block[i];
  }
  int number_host_code_faces=0;
  for ( i = 0; i < number_of_face_blocks; ++i){
    number_host_code_faces += number_faces_per_block[i];
  }
  int number_host_code_elements=0;
  for ( i = 0; i < number_of_element_blocks; ++i){
    number_host_code_elements += number_elements_per_block[i];
  }
  if (contact_number_of_processors(SearchComm)>1) {
    if( number_node_blocks ){
      ghosted_node_blocks = new ContactNodeBlock*[number_node_blocks];
      for( i=0 ; i<number_of_node_blocks ; ++i )
        ghosted_node_blocks[i] = new ContactNodeBlock( node_block_types[i], i, this);
    } else {
      ghosted_node_blocks = NULL;
    }
    if( number_face_blocks ){
      ghosted_face_blocks = new ContactFaceBlock*[number_face_blocks];
      for( i=0 ; i<number_of_face_blocks ; ++i )
        ghosted_face_blocks[i] = new ContactFaceBlock( face_block_types[i], i, this);
    } else {
      ghosted_face_blocks = NULL;
    }
    if( number_element_blocks ){
      ghosted_element_blocks = new ContactElementBlock*[number_element_blocks];
      for( i=0 ; i<number_of_element_blocks ; ++i )
        ghosted_element_blocks[i] = new ContactElementBlock( element_block_types[i], i, this);
    } else {
      ghosted_element_blocks = NULL;
    }
  } else {
    ghosted_node_blocks    = NULL;
    ghosted_face_blocks    = NULL;
    ghosted_element_blocks = NULL;
  }
  
  // Count the total number of nodes and construct the node blocks
  // which in turn construct the nodes.
  int  offset = 0;
  int* exo_ids  = (int*)exodus_node_ids;
  int* host_ids = (int*)global_node_ids;
  node_blocks = new ContactNodeBlock*[number_node_blocks];
  for( i=0 ; i<number_of_node_blocks ; ++i ){
#if CONTACT_DEBUG_PRINT_LEVEL>=5
    postream<<"Adding node block "<<i<<" with "
            <<number_nodes_per_block[i]<<" nodes\n";
    for( int j=0; j<number_nodes_per_block[i]; ++j ) {
      postream<<"  "<<j<<":  exo_id = "<<exo_ids[j]
              <<"  host_id = ("<<host_ids[2*j]<<", "<<host_ids[2*j+1]<<")\n";
    }
    postream.flush();
#endif
    node_blocks[i] = new ContactNodeBlock( node_block_types[i], 
					   i,
					   number_nodes_per_block[i],
					   number_of_nodes,
                                           exo_ids,
                                           host_ids,
					   &node_entity_types[number_of_nodes],
					   this );
    ContactTopologyEntity<Real>* entity;
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
      entity->HostGlobalArrayIndex(entity->HostArrayIndex()+offset);
    }
    offset   += number_nodes_per_block[i];
    exo_ids  += number_nodes_per_block[i];
    host_ids += 2*number_nodes_per_block[i];
  }
  node_list->BuildList(node_blocks, number_of_node_blocks,
                       no_parallel_consistency==ContactSearch::INACTIVE);
  primary_node_list->BuildList(node_blocks, number_of_node_blocks,
                               no_parallel_consistency==ContactSearch::INACTIVE);
  number_of_nodes = node_list->NumEntities();
  number_of_primary_nodes = node_list->NumEntities();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->OwnerProcArrayIndex(Nodes[i]->ProcArrayIndex());
    Nodes[i]->PrimaryProcArrayIndex(Nodes[i]->ProcArrayIndex());
    Nodes[i]->fcs_index = i;
  }

  // Set the entity key for nodes in block 1 if needed
  if (number_of_element_blocks + number_of_face_blocks == 0) {
    int nnodes = node_list->BlockNumEntities(0);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(0));
    for (int j=0; j<nnodes; ++j) {
      nodes[j]->Entity_Key( 0 );
    }
  }
  // Set the entity key for nodes in blocks 2-number_of_node_blocks
  int base_key = std::max(1,number_of_element_blocks + number_of_face_blocks);
  for( i=1 ; i<number_of_node_blocks ; ++i ){
    int entity_key = base_key + (i-1);
    int nnodes = node_list->BlockNumEntities(i);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
    for (int j=0; j<nnodes; ++j) {
      nodes[j]->Entity_Key( entity_key );
    }
  }

  int sum_input[2], sum_output[2];
  sum_input[0] = number_face_blocks;
  sum_input[1] = number_element_blocks;
  contact_global_sum(sum_input, sum_output, 2, SearchComm );
  int number_of_global_face_blocks   = sum_output[0];
  int number_of_global_element_blocks = sum_output[1];


  if (number_of_global_face_blocks>0) {
    // Count the total number of faces and construct the face blocks
    // which in turn construct the faces.
#if CONTACT_DEBUG_PRINT_LEVEL>=5
    int* fconn = face_connectivity;
#endif
    ContactTopologyEntity<Real>* entity;
    offset      = 0;
    host_ids    = (int*)global_face_ids;
    face_blocks = new ContactFaceBlock*[number_face_blocks];
    for( i=0 ; i<number_of_face_blocks ; ++i ){
#if CONTACT_DEBUG_PRINT_LEVEL>=5
      postream<<"Adding face block "<<i<<" with "
              <<number_faces_per_block[i]<<" faces\n";
      for( int j=0; j<number_faces_per_block[i]; ++j ) {
        postream<<"  "<<j<<":  host_id = ("<<host_ids[2*j]<<", "<<host_ids[2*j+1]<<")\n";
        int nnodes = search->Number_Nodes_Per_Face(face_block_types[i]);
        postream<<"      conn: ";
        for (int c=0; c<nnodes; ++c) {
          postream<<"  "<<exodus_node_ids[fconn[c]-1];
        }
        fconn += nnodes;
        postream<<"\n";
      }
      postream.flush();
#endif
      face_blocks[i] = new ContactFaceBlock( face_block_types[i], 
                                             i,
                                             number_faces_per_block[i],
                                             number_of_faces,
                                             host_ids,
                                             this );
      face_blocks[i]->FaceList()->IteratorStart();
      while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
        entity->HostGlobalArrayIndex(entity->HostArrayIndex()+offset);
      }
      offset   += number_faces_per_block[i];
      host_ids += 2*number_faces_per_block[i];
    }
    face_list->BuildList(face_blocks, number_of_face_blocks,
                         no_parallel_consistency==ContactSearch::INACTIVE);
    primary_face_list->BuildList(face_blocks, number_of_face_blocks,
                                 no_parallel_consistency==ContactSearch::INACTIVE);
    number_of_faces = face_list->NumEntities();
    number_of_primary_faces = face_list->NumEntities();
    ContactFace<Real>** Faces = 
      reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
    for (i=0; i<number_of_faces; ++i) {
      Faces[i]->OwnerProcArrayIndex(Faces[i]->ProcArrayIndex());
      Faces[i]->PrimaryProcArrayIndex(Faces[i]->ProcArrayIndex());
      Faces[i]->fcs_index = i;
    }

    // Connect the faces to their nodes (from the connectivity array)
    int index = 0;
    for( i=0 ; i<number_of_face_blocks ; ++i ){
      face_blocks[i]->FaceList()->IteratorStart();
      while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        for( k=0 ; k<face->Nodes_Per_Face() ; ++k ){
          int n = face_connectivity[index++]-1; // -1 Fortran->C
          ContactHostGlobalID global_id( global_node_ids[2*n+0], 
                                         global_node_ids[2*n+1] );
          ContactNode<Real>* node = static_cast<ContactNode<Real>*>(node_list->Find(global_id));
          PRECONDITION( node );
          face->ConnectNode(k, node);
        }
      }
    }
    Connect_Faces_to_Nodes();
  } else {
    face_blocks = NULL;
    number_of_faces = 0;
    number_of_primary_faces = 0;
  }


  if (number_of_global_element_blocks>0) {
    // Count the total number of elements and construct the element blocks
    // which in turn construct the elements.
#if CONTACT_DEBUG_PRINT_LEVEL>=5
    const int* econn = element_connectivity;
#endif
    ContactTopologyEntity<Real>* entity;
    offset         = 0;
    host_ids       = (int*)global_element_ids;
    element_blocks = new ContactElementBlock*[number_element_blocks];
    for( i=0 ; i<number_of_element_blocks ; ++i ){
#if CONTACT_DEBUG_PRINT_LEVEL>=5
      postream<<"Adding element block "<<i<<" with "
              <<number_elements_per_block[i]<<" elements\n";
      for( int j=0; j<number_elements_per_block[i]; ++j ) {
        postream<<"  "<<j<<":  host_id = ("<<host_ids[2*j]<<", "<<host_ids[2*j+1]<<")\n";
        int nnodes = 8; //search->Number_Nodes_Per_Element(element_block_types[i]);
        postream<<"      conn: ";
        for (int c=0; c<nnodes; ++c) {
          postream<<"  "<<exodus_node_ids[econn[c]-1];
        }
        econn += nnodes;
        postream<<"\n";
      }
      postream.flush();
#endif
      element_blocks[i] = new ContactElementBlock( element_block_types[i], 
                                                   i+number_of_face_blocks,
                                                   number_elements_per_block[i],
                                                   number_of_elements,
                                                   host_ids,
                                                   this );
      element_blocks[i]->ElemList()->IteratorStart();
      while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
        entity->HostGlobalArrayIndex(entity->HostArrayIndex()+offset);
      }
      offset   += number_elements_per_block[i];
      host_ids += 2*number_elements_per_block[i];
    }
    elem_list->BuildList(element_blocks, number_of_element_blocks,
                         no_parallel_consistency==ContactSearch::INACTIVE);
    primary_elem_list->BuildList(element_blocks, number_of_element_blocks,
                                 no_parallel_consistency==ContactSearch::INACTIVE);
    number_of_elements = elem_list->NumEntities();
    number_of_primary_elements = elem_list->NumEntities();
    ContactElement** Elements = 
      reinterpret_cast<ContactElement**>(elem_list->EntityList());
    for (i=0; i<number_of_elements; ++i) {
      Elements[i]->OwnerProcArrayIndex(Elements[i]->ProcArrayIndex());
      Elements[i]->PrimaryProcArrayIndex(Elements[i]->ProcArrayIndex());
      Elements[i]->fcs_index = i;
    }
    
    // Connect the elements to their nodes (from the connectivity array)
    int index = 0;
    for( i=0 ; i<number_of_element_blocks ; ++i ){
      element_blocks[i]->ElemList()->IteratorStart();
      while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
        ContactElement* element = static_cast<ContactElement*>(entity);
        for( k=0 ; k<element->Nodes_Per_Element() ; ++k ){
          int n = element_connectivity[index++]-1; // -1 Fortran->C
          ContactHostGlobalID global_id( global_node_ids[2*n+0], 
                                         global_node_ids[2*n+1] );
          ContactNode<Real>* node = static_cast<ContactNode<Real>*>(node_list->Find(global_id));
          PRECONDITION( node );
          element->ConnectNode(k, node);
        }
      }
    }
  } else {
    element_blocks = NULL;
    number_of_elements = 0;
    number_of_primary_elements = 0;
  }

  // Construct Edges and create connections

  if(have_shells) shell_handler->NumberNodes();

  if( dimensionality == 3 ) Construct_and_Connect_Edges(error_code);
  number_of_primary_edges = number_of_edges;
  
  // Check to see if the mesh is valid.  If not return
  error_code = (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, SearchComm );
  if( error_code ) return;

  // Create the communication lists (they must exist even 
  // if we don't communicate to any other processors)
#if CONTACT_DEBUG_PRINT_LEVEL>=5
  postream << "Node communication list" << "\n";
  postream << "  number of communication partners: " 
           << number_comm_partners << "\n";
  int* cnodes = comm_nodes;
  for ( i = 0; i < number_comm_partners; ++i ) {
    postream << "    communication partner id: " 
             << comm_proc_ids[i] << "\n";
    postream << "    number of nodes: " 
         << number_nodes_to_partner[i] << "\n";
    for ( int j = 0; j < number_nodes_to_partner[i]; ++j ) {
      postream << "       node " << j+1 
               << ":  " << exodus_node_ids[cnodes[j]-1] << "\n";
    }
    cnodes += number_nodes_to_partner[i];
  }
  postream.flush();
#endif
  int num_node_comm_ent = 0;
  for( i=0 ; i<number_comm_partners ; ++i )
    num_node_comm_ent += number_nodes_to_partner[i];
  ContactTopologyEntity<Real>** node_ent_comm_list = NULL;
  if( num_node_comm_ent ){
    node_ent_comm_list = new ContactTopologyEntity<Real>*[num_node_comm_ent];
    for( i=0 ; i<num_node_comm_ent ; ++i ) {
      int index = comm_nodes[i]-1;
      int base  = 2*index;
      ContactHostGlobalID global_id( global_node_ids[base+0], 
                                     global_node_ids[base+1] );
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(node_list->Find(global_id));
      POSTCONDITION(node);
      node_ent_comm_list[i] = node;
      node->Shared(true);
    }
  }

#ifndef CONTACT_NO_MPI
  Node_SymComm = new ContactSymComm( number_comm_partners,
				     comm_proc_ids,
				     number_nodes_to_partner,
				     node_ent_comm_list );
#endif

  delete [] node_ent_comm_list;
  delete [] node_entity_types;
  if ( ! shell_handler) {
    delete [] number_nodes_per_block;
    delete [] exodus_node_ids;
    delete [] global_node_ids;
    delete [] face_connectivity;
    delete [] number_nodes_to_partner;
    delete [] comm_nodes;
  }

  // now compute the owners for all nodes
  Compute_Owners(error_code);

#ifndef CONTACT_NO_MPI
  Node_AsymComm = new ContactAsymComm(*Node_SymComm );
#endif
  
  for( i=0 ; i<number_of_edge_blocks ; ++i ){
    edge_blocks[i]->EdgeList()->Rehash();
  }
  if (contact_number_of_processors(SearchComm)>1) {
    if( number_of_edge_blocks ){
      ghosted_edge_blocks = new ContactEdgeBlock*[number_of_edge_blocks];
      for( i=0 ; i<number_of_edge_blocks ; ++i )
        ghosted_edge_blocks[i] = new ContactEdgeBlock( edge_blocks[i]->Type(), i, this);
    } else {
      ghosted_edge_blocks = NULL;
    }
  } else {
    ghosted_edge_blocks = NULL;
  }
  
#ifndef CONTACT_NO_MPI
  contact_swap_edge_faces(SearchComm, *Edge_SymComm, *comm_buffer);
#endif

  ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->SetNeighborFacesInfo();
  }
  
  // Check to see if the mesh is valid across processors.  If not return
  error_code = (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, SearchComm );
  if( error_code ) return;

  // Initialize the scratch memory manager
  scratch_set = false;

  if( have_shells ) shell_handler->Complete_Initialization();
  
#if CONTACT_DEBUG_PRINT_LEVEL>=6
  Display_Entities(postream, 1);
  postream.flush();
#endif
  //Display_Entities(postream, 0);
  //postream.flush();
}

ContactTopology::ContactTopology( ContactErrors* Errors,
				  int Dimensionality,
				  int number_analytic_surfaces, 
				  int number_node_blocks,
	     const ContactSearch::ContactNode_Type* node_block_types,
				  int number_edge_blocks,
	     const ContactSearch::ContactEdge_Type* edge_block_types,
				  int number_face_blocks, 
	     const ContactSearch::ContactFace_Type* face_block_types,
		                  int number_element_blocks,
	     const ContactSearch::ContactElement_Type* element_block_types,
				  MPI_Comm communicator,
                                  ContactSearch* srch )
  :  computed_characteristic_length( true ),
     min_characteristic_length(0),
     ghosted_node_blocks(NULL),
     ghosted_edge_blocks(NULL),
     ghosted_face_blocks(NULL),
     ghosted_element_blocks(NULL),
     have_global_ghosting(false),
     have_local_ghosting(false),
     have_tied_ghosting(false)
{
  int i;
  // This constructor is only intended to construct a secondary decomposition
  // where things are done in pieces.
  
  topology_type = SECONDARY;
  SearchComm = communicator;
  comm_buffer = NULL;
  errors = Errors;
  search = srch;
  dimensionality = Dimensionality;
  node_block_ids = NULL;
  number_of_nodes = 0;
  number_of_edges = 0;
  number_of_faces = 0;
  number_of_elements = 0;
  number_of_node_blocks = number_node_blocks;
  number_of_edge_blocks = number_edge_blocks;
  number_of_face_blocks = number_face_blocks;
  number_of_element_blocks = number_element_blocks;
  number_of_analytic_surfaces = number_analytic_surfaces;
  number_of_added_analytic_surfaces = 0;
  two_configurations = false;
  node_blocks = NULL;
  edge_blocks = NULL;
  face_blocks = NULL;
  element_blocks = NULL;
  ghosted_node_blocks    = NULL;
  ghosted_face_blocks    = NULL;
  ghosted_element_blocks = NULL;
  if (number_analytic_surfaces>0) {
    AnalyticSurfaces = new ContactAnalyticSurface*[number_of_analytic_surfaces];
    for (i=0; i<number_of_analytic_surfaces; ++i) {
      AnalyticSurfaces[i] = NULL;
    }
  } else {
    AnalyticSurfaces = NULL;
  }
  number_of_comm_partners = 0;
#ifndef CONTACT_NO_MPI
  Node_AsymComm = NULL;
  Node_SymComm  = NULL;
  Edge_SymComm  = NULL;
  GhostingCommSpec = NULL;
  GhostFaces_ZoltanComm = NULL;
  GhostOthers_ZoltanComm = NULL;
  TiedCommSpec = NULL;
  TiedFaces_ZoltanComm = NULL;
  TiedNodes_ZoltanComm = NULL;
#endif  
  node_list = new ContactTopologyEntityList();
  edge_list = new ContactTopologyEntityList();
  face_list = new ContactTopologyEntityList();
  elem_list = new ContactTopologyEntityList();
  primary_node_list = NULL;
  primary_edge_list = NULL;
  primary_face_list = NULL;
  primary_elem_list = NULL;
#ifdef CONTACT_DEBUG_NODE
  number_debug_nodes = 0;
  debug_node_global_ids = NULL;
  debug_node_exodus_ids = NULL;
  debug_nodes = NULL;
#endif
  scratch_set = false;
  have_shells = false;
  shell_handler = NULL;
  old_to_new_node_map = NULL;
  
#ifndef CONTACT_NO_MPI
  ghosted_import_gids = NULL;
  ghosted_import_lids = NULL;
  ghosted_import_pids = NULL;
  num_ghosted_import  = 0;
  ghosted_export_gids = NULL;
  ghosted_export_lids = NULL;
  ghosted_export_pids = NULL;
  num_ghosted_export  = 0;
  tied_import_gids = NULL;
  tied_import_lids = NULL;
  tied_import_pids = NULL;
  num_tied_import  = 0;
  tied_export_gids = NULL;
  tied_export_lids = NULL;
  tied_export_pids = NULL;
  num_tied_export  = 0;
#endif

  Set_Up_Variable_Handles();

  // Construct Empty Node Blocks
  if( number_node_blocks ){
    node_blocks = new ContactNodeBlock*[number_node_blocks];
    for( i=0 ; i<number_of_node_blocks ; ++i )
      node_blocks[i] = new ContactNodeBlock( node_block_types[i], i, this);
  } else
    node_blocks = NULL;

  // Construct Empty Edge Blocks
  if( number_edge_blocks ){
    edge_blocks = new ContactEdgeBlock*[number_edge_blocks];
    for( i=0 ; i<number_of_edge_blocks ; ++i )
      edge_blocks[i] = new ContactEdgeBlock( edge_block_types[i], i, this);
  } else
    edge_blocks = NULL;

  // Construct Empty Face Blocks
  if( number_face_blocks ){
    face_blocks = new ContactFaceBlock*[number_face_blocks];
    for( i=0 ; i<number_of_face_blocks ; ++i )
      face_blocks[i] = new ContactFaceBlock( face_block_types[i], i, this);
  } else
    face_blocks = NULL;

  // Construct Empty Element Blocks
  if( number_of_element_blocks ){
    element_blocks = new ContactElementBlock*[number_element_blocks];
    for( i=0 ; i<number_of_element_blocks ; ++i )
      element_blocks[i]=new ContactElementBlock(element_block_types[i],i,this);
  } else
    element_blocks = NULL;

  // Initialize the scratch memory manager
  scratch_set = false;
}

ContactTopology::~ContactTopology()
{
#ifndef CONTACT_NO_MPI
  DeleteGhosting();
  DeleteGhostedNFIfaces();
  DeleteTiedFaces();
#endif

  if( node_list ) delete node_list;
  node_list = NULL;
  if( primary_node_list ) delete primary_node_list;
  primary_node_list = NULL;
  if( number_of_node_blocks && node_blocks ){
    for( int i=0 ; i<number_of_node_blocks ; ++i ) {
      if( node_blocks[i] ) delete node_blocks[i];
    }
    delete [] node_blocks;
    node_blocks = NULL;
  }
  if( number_of_node_blocks && ghosted_node_blocks ){
    for( int i=0 ; i<number_of_node_blocks ; ++i ) {
      if( ghosted_node_blocks[i] ) delete ghosted_node_blocks[i];
    }
    delete [] ghosted_node_blocks;
    ghosted_node_blocks = NULL;
  }
  
  if( edge_list ) delete edge_list;
  edge_list = NULL;
  if( number_of_edge_blocks && edge_blocks ){
    for( int i=0 ; i<number_of_edge_blocks ; ++i ) {
      if( edge_blocks[i] ) delete edge_blocks[i];
    }
    delete [] edge_blocks;
    edge_blocks = NULL;
  }
  if( number_of_edge_blocks && ghosted_edge_blocks ){
    for( int i=0 ; i<number_of_edge_blocks ; ++i ) {
      if( ghosted_edge_blocks[i] ) delete ghosted_edge_blocks[i];
    }
    delete [] ghosted_edge_blocks;
    ghosted_edge_blocks = NULL;
  }
  
  if( face_list ) delete face_list;
  face_list = NULL;
  if( primary_face_list ) delete primary_face_list;
  primary_face_list = NULL;
  if( number_of_face_blocks && face_blocks ){
    for( int i=0 ; i<number_of_face_blocks ; ++i ) {
      if( face_blocks[i] ) delete face_blocks[i];
    }
    delete [] face_blocks;
    face_blocks = NULL;
  }
  if( number_of_face_blocks && ghosted_face_blocks ){
    for( int i=0 ; i<number_of_face_blocks ; ++i ) {
      if( ghosted_face_blocks[i] ) delete ghosted_face_blocks[i];
    }
    delete [] ghosted_face_blocks;
    ghosted_face_blocks = NULL;
  }
  
  if( primary_elem_list ) delete primary_elem_list;
  primary_elem_list = NULL;
  if( elem_list ) delete elem_list;
  elem_list = NULL;
  if( number_of_element_blocks && element_blocks ){
    for( int i=0 ; i<number_of_element_blocks ; ++i ) {
      if( element_blocks[i] ) delete element_blocks[i];
    }
    delete [] element_blocks;
    element_blocks = NULL;
  }
  if( number_of_element_blocks && ghosted_element_blocks ){
    for( int i=0 ; i<number_of_element_blocks ; ++i ) {
      if( ghosted_element_blocks[i] ) delete ghosted_element_blocks[i];
    }
    delete [] ghosted_element_blocks;
    ghosted_element_blocks = NULL;
  }
  
  if( number_of_analytic_surfaces && AnalyticSurfaces ){
    for( int i=0 ; i<number_of_analytic_surfaces ; ++i ) {
      delete AnalyticSurfaces[i];
    }
    delete [] AnalyticSurfaces;
    AnalyticSurfaces = NULL;
  }

  delete [] Var_Handles;
  Var_Handles = NULL;

  if (shell_handler) delete shell_handler;
  shell_handler = NULL;
  
  if (old_to_new_node_map) delete [] old_to_new_node_map;

#ifndef CONTACT_NO_MPI  
  if( Node_AsymComm ) delete Node_AsymComm;
  if( Node_SymComm  ) delete Node_SymComm;
  if( Edge_SymComm  ) delete Edge_SymComm;
  if (GhostingCommSpec!=NULL)       delete GhostingCommSpec;
  if (GhostFaces_ZoltanComm!=NULL)  delete GhostFaces_ZoltanComm;
  if (GhostOthers_ZoltanComm!=NULL) delete GhostOthers_ZoltanComm;
#endif
#ifdef CONTACT_DEBUG_NODE
  if( number_debug_nodes ){
    for( int i=0 ; i<number_debug_nodes ; ++i )
      delete debug_node_global_ids[i];
    delete [] debug_nodes;
    delete [] debug_node_exodus_ids;
    delete [] debug_node_global_ids;
  }
#endif
}

#ifdef CONTACT_DEBUG_NODE
ContactSearch::ContactErrorCode ContactTopology::Add_Debug_Node( int exodus_id )
{
  // See if this node is already in our list
  bool already_have_this_node = false;
  for( int i=0 ; i<number_debug_nodes ; ++i ){
    if( debug_node_exodus_ids[i] == exodus_id ){
      already_have_this_node = true;
      break;
    }
  }
  if( already_have_this_node ){
    if( contact_processor_number( SearchComm ) == 0 ){
      std::sprintf(message,"Node %d has already been added",exodus_id);
      errors->Add_Error_Message(message);
    }
    return ContactSearch::INVALID_ID;
  }
    

  // See if I have this node
  // NOTE -- assuming max of 100 shells nodes created from this one
  //         host code node
  ContactNode<Real>* new_debug_nodes[100];
  for ( int i= 0; i < 100; ++i) new_debug_nodes[i] = NULL;
  int owning_proc = -1;
  int on_processor_id = -1;
  int num_debug_acme_nodes = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    bool found_debug_node = false;
    if ( node->Is_a_Shell_Node() ) {
      if( (static_cast<ContactShellNode*>(node))->Shell_Node_Base_ID() 
	  == exodus_id ) {
	found_debug_node = true;
      }
    } else {
      if( node->Exodus_ID() == exodus_id ) {
	found_debug_node = true;
      }
    }
    if ( found_debug_node) {
      new_debug_nodes[num_debug_acme_nodes] = node;
      if ( num_debug_acme_nodes == 0 ) {
	// NOTE -- assuming all created nodes owned by same processor
	owning_proc = new_debug_nodes[0]->Owner();
	on_processor_id = new_debug_nodes[0]->ProcArrayIndex();
      }
      ++num_debug_acme_nodes;
      POSTCONDITION( num_debug_acme_nodes <= 100);
    }
  }
  
  // Communicate to tell every processor about the first node of this set.
  //    We assume that if multiple acme nodes were created from the
  //    host code node (i.e. have same exodus id), then they are all
  //    owned by the same processor.
  //    The id is initialized to -1 so a global max will communicate the data

  int sum_input[3], sum_output[3];
  sum_input[0] = owning_proc;
  sum_input[1] = on_processor_id;
  sum_input[2] = num_debug_acme_nodes;
  contact_global_maximum( sum_input, sum_output, 3, SearchComm );
  owning_proc = sum_output[0];
  on_processor_id = sum_output[1];
  num_debug_acme_nodes = sum_output[2];

  POSTCONDITION( contact_processor_number(SearchComm) != owning_proc ||
		 new_debug_nodes[0] );
  
  if( owning_proc == -1 || on_processor_id == -1 ){
    if( contact_processor_number( SearchComm ) == 0 ){
      std::sprintf(message,"Node %d is not part of the contact surface",
              exodus_id);
      errors->Add_Error_Message(message);
    }
    return ContactSearch::INVALID_ID;
  }
  if( number_debug_nodes ){
    ContactHostGlobalID** old_gids = debug_node_global_ids;
    debug_node_global_ids = 
      new ContactHostGlobalID*[number_debug_nodes+num_debug_acme_nodes];
    std::memcpy( debug_node_global_ids, old_gids, 
	    number_debug_nodes*sizeof(ContactGlobalID*) );
    delete [] old_gids;
    ContactNode<Real>** old_nodes = debug_nodes;
    debug_nodes = new ContactNode<Real>*[number_debug_nodes+num_debug_acme_nodes];
    std::memcpy( debug_nodes, old_nodes, number_debug_nodes*sizeof(ContactNode<Real>*) );
    delete [] old_nodes;
    int* old_ids = debug_node_exodus_ids;
    debug_node_exodus_ids = new int[number_debug_nodes+num_debug_acme_nodes];
    std::memcpy( debug_node_exodus_ids, old_ids, number_debug_nodes*sizeof(int));
    delete [] old_ids;
  } else {
    debug_node_global_ids = new ContactHostGlobalID*[num_debug_acme_nodes];
    debug_nodes = new ContactNode<Real>*[num_debug_acme_nodes];
    debug_node_exodus_ids = new int[num_debug_acme_nodes];
  }
  for ( int i = 0; i < num_debug_acme_nodes; ++i ) {
    debug_node_global_ids[number_debug_nodes] = new 
      ContactHostGlobalID( owning_proc, on_processor_id );
    debug_nodes[number_debug_nodes] = new_debug_nodes[i];
    debug_node_exodus_ids[number_debug_nodes] = exodus_id;
    ++number_debug_nodes;
  }
  return ContactSearch::NO_ERROR;
}

bool ContactTopology::Is_a_Debug_Node( ContactNode<Real>* node )
{
#ifdef CONTACT_DEBUG_NODE
  for( int i=0 ; i<number_debug_nodes ; ++i ){
    if (node->Is_a_Shell_Node() ) {
      // if this is a shell node map back to original host code exodus ids
      if( (static_cast<ContactShellNode*> (node))->Shell_Node_Base_ID() == 
	  debug_node_exodus_ids[i] )
	return true;
    }
    else {
      if( node->Exodus_ID() == debug_node_exodus_ids[i] )
	return true;
    }
  }
#endif
  return false;
}
  
void ContactTopology::Display_Debug_Node_IDs( ContactParOStream& postream )
{
  if( number_debug_nodes == 0 ) return;
  postream << "\n";
  for( int i=0 ; i<number_debug_nodes ; ++i ){
    if ( debug_node_global_ids[i]->HiInt() != 
	 contact_processor_number(SearchComm) ) continue;
    PRECONDITION (debug_nodes[i]);
    postream << "Debug information will be output for Node ( host id:" 
	     << debug_node_exodus_ids[i] << " acme id:"
	     << debug_nodes[i]->Exodus_ID() << ")\n"
	     << "   Current load step is: " << search->StepNumber() << "\n";
  }
  postream.flush();
}

void ContactTopology::Display_Debug_Nodes( ContactParOStream& postream )
{
  int i,j;
  if( number_debug_nodes ){
    for( j=0 ; j<number_debug_nodes ; ++j ){
      if( debug_node_global_ids[j]->HiInt() == 
	  contact_processor_number( SearchComm ) ) {
	ContactNode<Real>* debug_node = debug_nodes[j];
	PRECONDITION( debug_node );
	Real* cur0 = debug_node->Variable(CURRENT_POSITION);
	postream << "\nDebug Information for Node ( host id:" 
		 << debug_node_exodus_ids[j] << " acme id:" 
		 << debug_node->Exodus_ID()
		 << ")\n";
	postream << "Current Position (x,y,z): " << cur0[0] << " " << cur0[1]
		 << " " << cur0[2] << "\n";
	if( two_configurations ){
	  Real* curp = debug_node->Variable(PREDICTED_POSITION);
	  postream << "Predicted Position (x,y,z): " << curp[0] << " " 
		   << curp[1] << " " << curp[2] << "\n";
	}
	Real* normal = debug_node->Variable(NODE_NORMAL);
	postream << "Node Normal (x,y,z): " << normal[0] << " " 
		 << normal[1] << " " << normal[2] << "\n";
	postream << "Number Kinematic Constraints = " 
		 << *debug_node->Variable(NUM_KIN_CONSTR) << "\n";
	if( *debug_node->Variable(NUM_KIN_CONSTR) ){
	  Real* con_vec = debug_node->Variable(KIN_CONSTR_VECTOR);
	  postream << "Constraint Vector (x,y,z): " << con_vec[0] << " " 
		   << con_vec[1] << " " << con_vec[2] << "\n";
	}
        if (debug_node->Number_NodeFace_Interactions()==0) {
	  postream << "Node has no Node-Face Interactions.\n";
	} else {
	  postream << "Node has " 
                   << debug_node->Number_NodeFace_Interactions() 
                   << " interaction(s).\n";
          ContactNodeEntityInteraction** interactions = 
	    debug_node->Get_NodeEntity_Interactions();
          for (i=0; i<debug_node->Number_NodeEntity_Interactions(); ++i) {
            if (interactions[i]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
            ContactNodeFaceInteraction* cnfi = 
              static_cast<ContactNodeFaceInteraction*>(interactions[i]);
	    postream << "   Interaction " << i << ".\n";
	    postream << "      Face: (" 
		     << cnfi->FaceEntityData()->host_gid[0] << "," 
		     << cnfi->FaceEntityData()->host_gid[1] << ")\n";
	    postream << "      Source           = " 
		     << cnfi->Source_Name().data() << "\n";
	    postream << "      Node Entity Key  = "
		     << cnfi->Scalar_Var(ContactNodeFaceInteraction::NODE_ENTITY_KEY )
		     << "\n";
	    postream << "      GapRate*dt       = " 
		     << cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR)
		     << "\n";
	    postream << "      Old Gap          = " 
		     << cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD)
		     << "\n";
	    Real* coord = 
	      cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
	    postream << "      Coordinates      = " << coord[0] << " "
		     << coord[1] << "\n";
	    Real* pb_dir =
	      cnfi->Get_Pushback_Dir();
	    postream << "      Pushback Dir     = " << pb_dir[0] << " " 
		     << pb_dir[1] << " " << pb_dir[2] << "\n";
	    Real* normal_dir = 
	      cnfi->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
	    postream << "      Normal Dir       = " << normal_dir[0] << " "
		     << normal_dir[1] << " " << normal_dir[2] << "\n";
	    Real* pf_normal =
	      cnfi->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
	    postream << "      PF Normal        = " << pf_normal[0] << " "
		     << pf_normal[1] << " " << pf_normal[2] << "\n";
          }
	}
        if (debug_node->Number_NodeSurface_Interactions()==0) {
	  postream << "Node has no Node-Surface Interactions.\n";
	} else {
          ContactNodeEntityInteraction** interactions = 
            debug_node->Get_NodeEntity_Interactions();
          for (i=0; i<debug_node->Number_NodeEntity_Interactions(); ++i) {
            if (interactions[i]->Get_Type()!=ContactNodeEntityInteraction::NODE_SURFACE_INTERACTION) continue;
            ContactNodeSurfaceInteraction* cnsi = 
              static_cast<ContactNodeSurfaceInteraction*>(interactions[i]);
	    postream << "Node Interacts with Analytic Surface "
	             << cnsi->Surface()->Global_ID() << "\n";
	    postream << "  Gap = "
	             << cnsi->Scalar_Var(ContactNodeSurfaceInteraction::GAP_CUR) 
	             << "\n";
	    Real* cpoint = 
	      cnsi->Vector_Var(ContactNodeSurfaceInteraction::CONTACT_POINT);
	    postream << "  Contact Point: " << cpoint[0] << " " << cpoint[1] 
	             << " " << cpoint[2] << "\n";
	    Real* snorm = 
	      cnsi->Get_Normal();
	    postream << "  Surface Normal: " << snorm[0] << " " << snorm[1]
	             << " " << snorm[2] << "\n";
	  }
	}
      }
    }
    postream.flush();
  }
}
#endif

void ContactTopology::Add_Analytic_Surface( ContactAnalyticSurface* surface )
{
  // Add the new surface
  int my_proc   = contact_processor_number(SearchComm);
  surface->Ownership(ContactTopologyEntity<Real>::OWNED);
  surface->Owner(my_proc);
  surface->Secondary_Owner(my_proc);

  AnalyticSurfaces[number_of_added_analytic_surfaces++] = surface;
}

ContactSearch::ContactErrorCode
ContactTopology::Set_Analytic_Surface_Configuration( int entity_key, 
						     const Real* data )
{
  int array_offset = entity_key - number_of_face_blocks;
  if( array_offset < 0 || array_offset >= number_of_analytic_surfaces ){
    std::sprintf(message,"Set_Analytic_Surface_Configuration: Invalid ID = %d",
	    entity_key );
    return ContactSearch::INVALID_ID;
  }
  AnalyticSurfaces[array_offset]->Set_Configuration( data );
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTopology::Set_NodeBlk_RemainingGap( int id, 
					   const Real* gap )
{
  if( id < 1 || id > number_of_node_blocks ){
    std::sprintf(message,"Set_NodeBlk_RemainingGap: Unknown Node Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }
  if( shell_handler ){
    shell_handler->Set_NodeBlk_RemainingGap( id, node_blocks[id-1]->NodeList(),
    					     REMAINING_GAP, gap );
  } else {
    // Note that the ID is in fortran numbering so decrement by 1
    int nnodes = primary_node_list->BlockNumEntities(id-1);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
    for (int i=0; i<nnodes; ++i) {
      ContactNode<Real>* node = nodes[i];
      int index = node->HostArrayIndex();
      Real* remaining_gap = node->Variable(REMAINING_GAP);
      for( int j=0 ; j<dimensionality ; ++j ) {
	remaining_gap[j] = gap[dimensionality*index + j];
      }
    }
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTopology::Set_NodeBlk_GhostingGap( int id, 
					  const Real* gap )
{
  if( id < 1 || id > number_of_node_blocks ){
    std::sprintf(message,"Set_NodeBlk_GhostingGap: Unknown Node Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }
  if( shell_handler ){
    shell_handler->Set_NodeBlk_GhostingGap( id, node_blocks[id-1]->NodeList(),
    					    NODE_GHOST_GAP, gap );
  } else {
    // Note that the ID is in fortran numbering so decrement by 1
    int nnodes = primary_node_list->BlockNumEntities(id-1);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
    for (int i=0; i<nnodes; ++i) {
      ContactNode<Real>* node = nodes[i];
      int index = node->HostArrayIndex();
      Real* ghosting_gap = node->Variable(NODE_GHOST_GAP);
      for( int j=0 ; j<dimensionality ; ++j ) {
	ghosting_gap[j] = gap[dimensionality*index + j];
      }
    }
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTopology::Set_NodeBlk_KinConstr( int id, const int* num_kcs,
					const Real* kc_vectors )
{
  if( id < 1 || id > number_of_node_blocks ){
    std::sprintf(message,"Set_NodeBlk_KinConstr: Unknown Node Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }
  if( shell_handler ){
    shell_handler->Set_NodeBlk_KinConstr( id, node_blocks[id-1]->NodeList(),
					  NUM_KIN_CONSTR, KIN_CONSTR_VECTOR,
					  num_kcs, kc_vectors );
  } else {
    // Note that the ID is in fortran numbering so decrement by 1
    int nnodes = primary_node_list->BlockNumEntities(id-1);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
    for (int i=0; i<nnodes; ++i) {
      ContactNode<Real>* node = nodes[i];
      int index = node->HostArrayIndex();
      Real* num_kc = node->Variable(NUM_KIN_CONSTR);
      num_kc[0] = num_kcs[index];
      Real* kc_vector = node->Variable(KIN_CONSTR_VECTOR);
      for( int j=0 ; j<dimensionality ; ++j ) {
	kc_vector[j] = kc_vectors[dimensionality*index + j];
      }
    }
  }
  return ContactSearch::NO_ERROR;
}


ContactSearch::ContactErrorCode 
ContactTopology::Set_NodeBlk_Positions( int id, const Real* positions )
{
  if( id < 1 || id > number_of_node_blocks ){
    std::sprintf(message,"Set_NodeBlk_Positions: Unknown Node Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }
  if( shell_handler ){
    shell_handler->Set_NodeBlk_Position( id, node_blocks[id-1]->NodeList(),
					 CURRENT_POSITION, positions );
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=10
    ContactParOStream& postream = search->ParOStream();
    postream<<"Set_NodeBlk_Positions() for block "<<id<<"\n";
#endif
    // Note that the ID is in fortran numbering so decrement by 1
    int nnodes = primary_node_list->BlockNumEntities(id-1);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
    for (int i=0; i<nnodes; ++i) {
      ContactNode<Real>* node = nodes[i];
      int index = node->HostArrayIndex();
      Real* position = node->Variable(CURRENT_POSITION);
      for( int j=0 ; j<dimensionality ; ++j ) {
        position[j] = positions[dimensionality*index + j];
      }
      for( int j=dimensionality ; j<3 ; ++j ) {
        position[j] = 0.0;
      }
#if CONTACT_DEBUG_PRINT_LEVEL>=10
      postream<<"  node "<<node->Global_ID()
              <<":  ("<<position[0]<<", "<<position[1]<<", "<<position[2]<<")\n";
#endif
    }
#if CONTACT_DEBUG_PRINT_LEVEL>=10
    postream.flush();
#endif
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactTopology::Set_NodeBlk_Positions_2( int id, const Real* positions )
{
  two_configurations = true;
  // Note that the ID is in fortran numbering so decrement by 1
  if( id < 1 || id > number_of_node_blocks ){
    std::sprintf(message,"Set_NodeBlk_Positions_2: Unknown Node Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }
  if( shell_handler ){
    shell_handler->Set_NodeBlk_Position( id, node_blocks[id-1]->NodeList(),
					 PREDICTED_POSITION, positions );
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=10
    ContactParOStream& postream = search->ParOStream();
    postream<<"Set_NodeBlk_Positions_2() for block "<<id<<"\n";
#endif
    // Note that the ID is in fortran numbering so decrement by 1
    int nnodes = primary_node_list->BlockNumEntities(id-1);
    ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
    for (int i=0; i<nnodes; ++i) {
      ContactNode<Real>* node = nodes[i];
      int index = node->HostArrayIndex();
      Real* position = node->Variable(PREDICTED_POSITION);
      for( int j=0 ; j<dimensionality ; ++j ) {
        position[j] = positions[dimensionality*index + j];
      }
      for( int j=dimensionality ; j<3 ; ++j ) {
        position[j] = 0.0;
      }
#if CONTACT_DEBUG_PRINT_LEVEL>=10
      postream<<"  node "<<node->Global_ID()
              <<":  ("<<position[0]<<", "<<position[1]<<", "<<position[2]<<")\n";
#endif
    }

#if CONTACT_DEBUG_PRINT_LEVEL>=10
    postream.flush();
#endif
  }
  return ContactSearch::NO_ERROR;
}


ContactSearch::ContactErrorCode 
ContactTopology::Set_NodeBlk_Attributes(
			  ContactSearch::Node_Block_Attribute attr,
			  int id,  const Real* attributes )
{
  if( id < 1 || id > number_of_node_blocks ){
    std::sprintf(message,"Set_NodeBlk_Attributes: Unknown Node Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }
  if( shell_handler ){
    errors->Add_Error_Message( 
     "Meshes with Shells and Node Attributes are not supported" );
    return ContactSearch::UNIMPLEMENTED_FUNCTION;
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=10
    ContactParOStream& postream = search->ParOStream();
    postream<<"Set_NodeBlk_Attributes() for block "<<id<<"\n";
#endif
    switch( attr ){
    case ContactSearch::RADIUS:{
#if CONTACT_DEBUG_PRINT_LEVEL>=10
      postream<<"  attribute: RADIUS\n";
#endif
      Real rmax = 0.0;
      // Note that the ID is in fortran numbering so decrement by 1
      int nnodes = primary_node_list->BlockNumEntities(id-1);
      ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
      for (int i=0; i<nnodes; ++i) {
        ContactNode<Real>* node = nodes[i];
        int index = node->HostArrayIndex();
        Real* attribute = node->Variable(NODE_RADIUS);
        *attribute = attributes[index];
        rmax = std::max(rmax,attributes[index]);
#if CONTACT_DEBUG_PRINT_LEVEL>=10
        postream<<"  node "<<node->Global_ID()
                <<":  "<<*attribute<<"\n";
#endif
      }
      rmax = contact_global_maximum( rmax, SearchComm );
      node_blocks[id-1]->Rmax(rmax);
      node_blocks[id-1]->Has_Attributes(true);
      node_blocks[id-1]->Has_Radius_Attributes(true);
      break;
    }
    case ContactSearch::NORMAL:{
#if CONTACT_DEBUG_PRINT_LEVEL>=10
      postream<<"  attribute: NORMAL\n";
#endif
      // Note that the ID is in fortran numbering so decrement by 1
      int nnodes = primary_node_list->BlockNumEntities(id-1);
      ContactNode<Real>** nodes = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(id-1));
      for (int i=0; i<nnodes; ++i) {
        ContactNode<Real>* node = nodes[i];
        int index = node->HostArrayIndex();
        Real* attribute = node->Variable(NODE_NORMAL);
        for( int j=0 ; j<dimensionality ; ++j ) {
          attribute[j] = attributes[dimensionality*index + j];
        }
#if CONTACT_DEBUG_PRINT_LEVEL>=10
        postream<<"  node "<<node->Global_ID()
                <<":  ("<<attribute[0]<<", "<<attribute[1]<<", "<<attribute[2]<<")\n";
#endif
      }
      node_blocks[id-1]->Has_Attributes(true);
      node_blocks[id-1]->Has_Normal_Attributes(true);
      break;
    }
    default: {
      CString msg("Set_NodeBlk_Attributes: Unknown Node Attribute");
      errors->Add_Error_Message( msg.data() );
      break;
    }
    }
#if CONTACT_DEBUG_PRINT_LEVEL>=10
    postream.flush();
#endif
  }
  return ContactSearch::NO_ERROR;
}


ContactSearch::ContactErrorCode 
ContactTopology::Set_FaceBlk_Attributes(
			  ContactSearch::Face_Block_Attribute attr,
			  int id,  const Real* attributes )
{
  if( id < 1 || id > number_of_face_blocks ){
    std::sprintf(message,"Set_FaceBlk_Attributes: Unknown Face Block ID = %d",id);
    errors->Add_Error_Message( message );
    return ContactSearch::ID_NOT_FOUND;
  }

  if( !search->Is_a_Shell_Face(face_blocks[id-1]->Type()) ){
    std::sprintf(message,
	    "Set_FaceBlk_Attributes: Face Block ID %d is not a shell Block",id);
    errors->Add_Error_Message( message );
    return ContactSearch::INVALID_DATA;
  }
  // Note that the ID is in fortran numbering so decrement by 1
  int i;
  ContactFaceBlock* block = face_blocks[id-1];
  int nfaces = primary_face_list->BlockNumEntities(id-1);
  ContactFace<Real>** faces  = reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(id-1));
  switch( attr ){
  case ContactSearch::SHELL_THICKNESS:{
    if( nfaces ){
      PRECONDITION( block->Type() == ContactSearch::SHELLQUADFACEL4 ||
		    block->Type() == ContactSearch::SHELLTRIFACEL3 );
      switch( block->Type() ){
      case( ContactSearch::SHELLQUADFACEL4 ):{
        for (i=0; i<nfaces; ++i) {
	  ContactShellQuadFaceL4<Real>* face = 
            static_cast<ContactShellQuadFaceL4<Real>*>(faces[i]);
	  face->Thickness( attributes[face->HostArrayIndex()] );
	}
	break;
      }
      case( ContactSearch::SHELLTRIFACEL3 ):{
        for (i=0; i<nfaces; ++i) {
	  ContactShellTriFaceL3<Real>* face = 
            static_cast<ContactShellTriFaceL3<Real>*>(faces[i]);
	  face->Thickness( attributes[face->HostArrayIndex()] );
	}
	break;
      }
      default:
        POSTCONDITION(0);
        break;
      }
    }
    break;
  }
  case ContactSearch::LOFTING_FACTOR:{
    if( nfaces ){
      PRECONDITION( block->Type() == ContactSearch::SHELLQUADFACEL4 ||
		    block->Type() == ContactSearch::SHELLTRIFACEL3 );
      switch( block->Type() ){
      case( ContactSearch::SHELLQUADFACEL4 ):{
        for (i=0; i<nfaces; ++i) {
	  ContactShellQuadFaceL4<Real>* face = 
            static_cast<ContactShellQuadFaceL4<Real>*>(faces[i]);
	  face->Lofting_Factor( attributes[face->HostArrayIndex()] );
	}
	break;
      }
      case( ContactSearch::SHELLTRIFACEL3 ):{
        for (i=0; i<nfaces; ++i) {
	  ContactShellTriFaceL3<Real>* face = 
            static_cast<ContactShellTriFaceL3<Real>*>(faces[i]);
	  face->Lofting_Factor( attributes[face->HostArrayIndex()] );
	}
	break;
      }
      default:
        POSTCONDITION(0);
        break;
      }
    }
    break;
  }
  default: {
    CString msg("Set_FaceBlk_Attributes: Unknown Face Attribute");
    errors->Add_Error_Message( msg.data() );
    break;
  }
  }
  return ContactSearch::NO_ERROR;
}

void ContactTopology::Connect_Faces_to_Nodes()
{
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());

  for (int i=0; i<number_of_faces; ++i) {
    for (int k=0; k<Faces[i]->Nodes_Per_Face(); ++k) {
      ContactNode<Real>* node = Faces[i]->Node(k);
      node->Connect_Face( Faces[i] );
    }
  }
  
  for (int i=0; i<number_of_nodes; ++i) {
    Nodes[i]->SortConnectedFaces();
  }
}

void ContactTopology::Construct_and_Connect_Edges(
				   ContactSearch::ContactErrorCode& error_code)
{
  PRECONDITION( number_of_edges == 0 );
  int i,j,k,n;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());

  // set temp_tag to the face number in the list
  int count = 0;
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->temp_tag = count++;
  }

  int max_edges=4;
  int size = number_of_faces*max_edges;
  int * table = new int[size];
  for( i=0 ; i<size ; ++i ) table[i] = 0;
  
  int  nedge_types = ContactSearch::NEDGE_TYPES;
  int* edge_types  = new int[nedge_types];
  int* edge_count  = new int[nedge_types];
  for( i=0 ; i<nedge_types ; ++i ){
    edge_types[i] = 0;
    edge_count[i] = 0;
  }

  for( i=0 ; i<number_of_face_blocks ; ++i ){
    edge_types[face_blocks[i]->EdgeType()]++;
    int nfaces = face_list->BlockNumEntities(i);
    ContactFace<Real>** faces  = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
    for (k=0; k<nfaces; ++k) {
      ContactFace<Real>* face = faces[k];
      int index = face->temp_tag;
      for( j=0 ; j<face->Edges_Per_Face() ; ++j ){
	if( table[index*max_edges+j] == 0 ){
	  ++number_of_edges;
          edge_count[face_blocks[i]->EdgeType()]++;
	  table[index*max_edges+j] = edge_count[face_blocks[i]->EdgeType()];
	  // Get the nodes that make up the edge
	  ContactNode<Real> *node[3];
	  face->Get_Edge_Nodes( j,node );
	  // Get the faces that connect to these two nodes (1 or 2)
	  ContactFace<Real> *face1,*face2;
	  int Number_Of_Faces = Get_Faces_Connected_to_Nodes( node[0],node[1],
							      &face1,&face2,
							      error_code );
	  if( Number_Of_Faces == 2 ){ 
	    ContactFace<Real>* neighbor=NULL;
	    if( face1 != face )
	      neighbor = face1;
	    else if( face2 != face )
	      neighbor = face2;
	    POSTCONDITION(neighbor!=NULL);
            int edge_num = neighbor->Get_Edge_Number( node );
	    if( edge_num != -1 )
	      table[ neighbor->temp_tag*max_edges+edge_num ] = 
		edge_count[face_blocks[i]->EdgeType()];
	  }
	}
      }
    }
  }

  // create the edge blocks
  number_of_edge_blocks = 0;
  for( i=1 ; i<nedge_types ; ++i ){
    if( edge_types[i] > 0 ) ++number_of_edge_blocks;
  }
  int* offset = NULL;
  if (number_of_edge_blocks>0) {
    offset = new int[number_of_edge_blocks];
    edge_blocks = new ContactEdgeBlock*[number_of_edge_blocks];
    int block=0;
    int entity_key = 1;
    int Next_ID = 1;
    for (n=1; n<nedge_types ; ++n ) {
      if( edge_types[n] > 0 ){
        offset[block] = Next_ID-1;
        edge_blocks[block] =
          new ContactEdgeBlock( (ContactSearch::ContactEdge_Type)n,
                                block,entity_key, edge_count[n],
                                Next_ID, this );
        ++block;
      }
    }
  }
  
  if( number_of_edges ){ 
    int myproc = contact_processor_number(SearchComm);
    int block = 0;
    for (n=1; n<nedge_types ; ++n ) {
      if( edge_types[n] > 0 ){
	for( i=0 ; i<number_of_face_blocks ; ++i ){
	  if(face_blocks[i]->EdgeType()==n) {
            int nfaces = face_list->BlockNumEntities(i);
            ContactFace<Real>** faces  = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
            for (k=0; k<nfaces; ++k) {
              ContactFace<Real>* face = faces[k];
	      int ind = face->temp_tag;
	      for( j=0 ; j<face->Edges_Per_Face() ; ++j ){
		ContactNode<Real> *node[3];
		int edge_id = table[ind*max_edges+j]-1;
		PRECONDITION( edge_id>=0 && edge_id<number_of_edges );
                ContactHostGlobalID global_id( myproc, offset[block]+edge_id+1 );
                ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(edge_blocks[block]->EdgeList()->Find(global_id));
                POSTCONDITION(edge);
		face->ConnectEdge( j, edge );
		face->Get_Edge_Nodes( j, node );
		for (int ii=0; ii<edge->Nodes_Per_Edge(); ++ii) {
		  edge->ConnectNode(ii,node[ii]);
		}
	      }
            }
          }
        }
	++block;
      }
    }
    
    edge_list->BuildList(edge_blocks, number_of_edge_blocks,
                         no_parallel_consistency==ContactSearch::INACTIVE);
    edge_list->SortByNodeGID();
    ContactEdge<Real>** Edges = 
      reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
    POSTCONDITION(number_of_edges == edge_list->NumEntities());
    for (i=0; i<number_of_edges; ++i) {
      Edges[i]->OwnerProcArrayIndex(Edges[i]->ProcArrayIndex());
      Edges[i]->PrimaryProcArrayIndex(Edges[i]->ProcArrayIndex());
      Edges[i]->fcs_index = i;
    }
    Connect_Faces_to_Edges();
  }

  delete [] edge_types;
  delete [] edge_count;
  delete [] table;
  if (offset) delete [] offset;

}


void ContactTopology::Connect_Faces_to_Edges()
{
  int i,k;

  // Connect the faces to the edges
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  int nfaces = face_list->NumEntities();
  for( i=0 ; i<nfaces ; ++i ){
    for( k=0 ; k<Faces[i]->Edges_Per_Face() ; ++k ){
      ContactEdge<Real>* edge = Faces[i]->Edge(k);
      edge->ConnectFace( Faces[i] );
    }
  }
}



int ContactTopology::Get_Faces_Connected_to_Nodes( ContactNode<Real>* node0, 
						   ContactNode<Real>* node1,
						   ContactFace<Real>** face1,
						   ContactFace<Real>** face2,
				    ContactSearch::ContactErrorCode& error_code)
{
  PRECONDITION( node0 && node1 );
  int nfaces = 0;
  int num_face_con = 0;
  for( int i=0 ; i<node0->Number_Face_Connections() ; ++i ){
    ContactFace<Real>* face = node0->GetFace(i);
    for( int j=0 ; j<node1->Number_Face_Connections() ; ++j ){
      if( face == node1->GetFace(j) ){
        int mwg = num_face_con;
	switch(mwg){
	case 0:
	  *face1 = face;
	  ++num_face_con;
	  break;
	case 1:
	  *face2 = face;
	  ++num_face_con;
	  break;
	case 2:
	  std::sprintf(message, 
		  "More than two faces connected to an edge at nodes (%d, %d) & (%d, %d)",
		  node0->Global_ID().HiInt(),node0->Global_ID().LoInt(), 
                  node1->Global_ID().HiInt(),node1->Global_ID().LoInt() );
	  errors->Add_Error_Message( message );
	  
	  std::sprintf(message, 
		  "  Face 0 GID = (%d, %d) in block %d",
      	          (*face1)->Global_ID().HiInt(),(*face1)->Global_ID().LoInt(), 
                  (*face1)->BlockID() );
	  errors->Add_Error_Message( message );
	  std::sprintf(message, 
		  "  Face 1 GID = (%d, %d) in block %d",
      	          (*face2)->Global_ID().HiInt(),(*face2)->Global_ID().LoInt(), 
                  (*face2)->BlockID() );
	  errors->Add_Error_Message( message );
	  std::sprintf(message, 
		  "  Face 2 GID = (%d, %d) in block %d",
      	          face->Global_ID().HiInt(),face->Global_ID().LoInt(), 
                  face->BlockID() );
	  errors->Add_Error_Message( message );

          //
          //  Store the invalid edge for later retreival by the host code
          //
          {
            ContactShellNode *ShellNode0 = dynamic_cast<ContactShellNode*>(node0);
            ContactShellNode *ShellNode1 = dynamic_cast<ContactShellNode*>(node1);
            int global_id0;
            int global_id1;
            if(ShellNode0) {
              global_id0 = ShellNode0->Shell_Node_Base_ID();
            } else {
             global_id0 = node0->Global_ID().LoInt();
            }
            if(ShellNode1) {
              global_id1 = ShellNode1->Shell_Node_Base_ID();
            } else {
              global_id1 = node1->Global_ID().LoInt();
            }
            search->append_invalid_edge(global_id0, global_id1);
          }
          nfaces = 3;
	  num_face_con = -1;
	  error_code = ContactSearch::INVALID_DATA;
	  break;
	case -1:
	  std::sprintf(message, 
		  "  Face %d GID = (%d, %d) in block %d",nfaces,
      	          face->Global_ID().HiInt(),face->Global_ID().LoInt(), 
                  face->BlockID() );
	  errors->Add_Error_Message( message );
          ++nfaces;
          num_face_con = -1;
	  error_code = ContactSearch::INVALID_DATA;
	  break;
	}
      }
    }
  }
  return(num_face_con);
}

void ContactTopology::Display(ContactParOStream& postream)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=6
  postream << "\n\n======================================================\n";
  postream << "\n\n  Contact Topology:\n";
  postream << "     Number of Nodes    = " << number_of_nodes << "\n";
  if( dimensionality == 3) 
  postream << "     Number of Edges    = " << number_of_edges << "\n";
  postream << "     Number of Faces    = " << number_of_faces << "\n";
  postream << "     Number of Elements = " << number_of_elements << "\n";
#if CONTACT_DEBUG_PRINT_LEVEL>=7
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=8
  postream << "\n\n  Contact Node     Global ID    Exodus ID         POSITION\n";
  postream <<     " --------------   -----------  -----------        -----------\n" ;

  int oldprecision = std::cout.precision();
  std::cout.precision(16);
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    postream << "	 "  << i << "		   "  
             << node->Global_ID() 
             << "	   "
             << node->Exodus_ID()
             << "	   ";
    for( int k=0 ; k<dimensionality ; ++k )
      postream << node->Variable(CURRENT_POSITION)[k] << "  ";
    postream << "\n";
    if( two_configurations ){
      postream << "				  ";
      for( int k=0 ; k<dimensionality ; ++k )
        postream << node->Variable(PREDICTED_POSITION)[k] << "  ";
      postream << "\n";
    }
  }
  std::cout.precision(oldprecision);
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=7
  postream << "\n\n  Node Connectivity for Faces " << "\n";
  for( int i=0 ; i<number_of_faces ; ++i ){
    ContactFace<Real>* face = Faces[i];
    postream << "    Face " << face->Global_ID() 
           << " has node connectivity ";
    for( int k=0 ; k<face->Nodes_Per_Face() ; ++k )
      postream << face->Node(k)->Global_ID() << " ";
    postream << "\n";
  }
  if( dimensionality == 3){
    postream << "\n\n  Edge Connectivity for Faces " << "\n";
    for( int i=0 ; i<number_of_faces ; ++i ){
        ContactFace<Real>* face = Faces[i];
	postream << "    Face " << face->Global_ID() 
	     << " has edge connectivity ";
	for( int k=0 ; k<face->Edges_Per_Face() ; ++k )
	  postream << face->Edge(k)->Global_ID() << " ";
	postream << "\n";
    }
    postream << "\n\n  Node Connectivity for Edges and Curvature" << "\n";
    for( int i=0 ; i<number_of_edges ; ++i ){
      ContactEdge<Real>* edge = Edges[i];
      postream << "    Edge " << edge->Global_ID() 
             << " has node connectivity ";
      for( int k=0 ; k<edge->Nodes_Per_Edge() ; ++k ){
        postream << edge->Node(k)->Global_ID() << " ";
      }
      postream << "   Curvature = " << *edge->Variable(CURVATURE) <<"\n";
    }
  }
  postream << "\n\n  Face Connectivity for Nodes\n";
  for( int i=0 ; i<number_of_nodes ; ++i ){
    ContactNode<Real>* node = Nodes[i];
    postream << "    Node " << node->Global_ID() << " has connectivity ";
    for( int k=0 ; k<node->Number_Face_Connections() ; ++k ) {
      postream << node->GetFace(k)->Global_ID() << " ";
    }
    postream << "\n";
  }
  if( dimensionality == 3){
    postream << "\n\n  Face Connectivity for Edges\n";
    for( int i=0 ; i<number_of_edges ; ++i ){
      ContactEdge<Real>* edge = Edges[i];
      postream << "    Edge " << edge->Global_ID() << " has connectivity ";
      for( int k=0 ; k<edge->Number_Face_Connections() ; ++k )
        postream << edge->Face(k)->Global_ID() << " ";
      postream << "\n";
    }
  }
  for( int i=0 ; i<number_of_analytic_surfaces ; ++i )
    AnalyticSurfaces[i]->Display(postream);
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=8
  postream << "\n\n Node Normals\n";
  for( int i=0 ; i<number_of_nodes ; ++i ){
    ContactNode<Real>* node = Nodes[i];
    postream << "   Node " << node->Global_ID() << " has normal ";
    for( int k=0 ; k<dimensionality ; ++k ) {
      postream << node->Variable(NODE_NORMAL)[k] << "  ";
    }
    postream << "\n";
  }
  
  postream << "\n\n Face Normals\n";
  for( int i=0 ; i<number_of_faces ; ++i ){
    ContactFace<Real>* face = Faces[i];
    postream << "   Face " << face->Global_ID() << " has normal ";
    for( int k=0 ; k<dimensionality ; ++k )
      postream << face->Variable(FACE_NORMAL)[k] << "  ";
    postream << "\n";
  }
  postream << "\n\n Face Centroids\n";
  for( int i=0 ; i<number_of_faces ; ++i ){
    ContactFace<Real>* face = Faces[i];
    postream << "   Face " << face->Global_ID() << " has centroid ";
    for( int k=0 ; k<dimensionality ; ++k )
      postream << face->Variable(CENTROID)[k] << "  ";
    postream << "\n";
  }
#endif
  postream << "======================================================\n\n";
  postream.flush();
#endif
}

void 
ContactTopology::Display_Entities( ContactParOStream& postream, int comm_flag )
{
  postream.flush();std::cout<<std::flush;
  int myproc = contact_processor_number(SearchComm);

  if (myproc==0) {
    std::cout<<"\n\nNumber of node blocks = "<<number_of_node_blocks<<std::endl;
  }
  postream<<"\n========== Node List (Size = "
          <<node_list->NumEntities()<<") ==========\n";
  node_list->Display(postream);
  postream.flush();
  
  if (myproc==0) {
    std::cout<<"\n\nNumber of edge blocks = "<<number_of_edge_blocks<<std::endl;
  }
  postream<<"\n========== Edge List (Size = "
          <<edge_list->NumEntities()<<") ==========\n";
  edge_list->Display(postream);
  postream.flush();
  
  if (myproc==0) {
    std::cout<<"\n\nNumber of face blocks = "<<number_of_face_blocks<<std::endl;
  }
  postream<<"\n========== Face List (Size = "
          <<face_list->NumEntities()<<") ==========\n";
  face_list->Display(postream);
  postream.flush();
  
  if (myproc==0) {
    std::cout<<"\n\nNumber of element blocks = "<<number_of_element_blocks<<std::endl;
  }
  postream<<"\n========== Element List (Size = "
          <<elem_list->NumEntities()<<") ==========\n";
  elem_list->Display(postream);
  postream.flush();

#ifndef CONTACT_NO_MPI
  if (comm_flag) {
    int number_of_global_nodes = contact_global_sum(number_of_nodes, SearchComm);
    int number_of_global_edges = contact_global_sum(number_of_edges, SearchComm);

    if (number_of_global_nodes && contact_number_of_processors(SearchComm)>1) {
      int i, j;
      if (Node_SymComm) {
        postream << "Node symmetric communication list" << "\n";
        postream << "  number of communication partners: " 
                 << Node_SymComm->Num_Comm_Partners() << "\n";
        for ( i = 0; i < Node_SymComm->Num_Comm_Partners(); ++i ) {
          postream << "    communication partner id: " 
                   << Node_SymComm->Comm_Proc_ID(i) << "\n";
          postream << "    number of nodes: " 
                   << Node_SymComm->Num_to_Proc(i) << "\n";
          ContactTopologyEntity<Real>** nodes_to_print = Node_SymComm->Entity_List(i);
          for ( j = 0; j < Node_SymComm->Num_to_Proc(i); ++j ) {
            ContactNode<Real> * node = 
              static_cast<ContactNode<Real> *>(nodes_to_print[j]);
            postream << "       node " << j+1 
                     << "  " << node->Global_ID() << "\n";
          }
        }
        postream.flush();
      }
      if (Node_AsymComm) {
        postream << "Node asymmetric communication list" << "\n";
        postream << "  number of export communication partners: " 
                 << Node_AsymComm->Num_Export_Comm_Partners() << "\n";
        for ( i = 0; i < Node_AsymComm->Num_Export_Comm_Partners(); ++i ) {
          postream << "    export communication partner id: " 
                   << Node_AsymComm->Export_Comm_Proc_ID(i) << "\n";
          postream << "    number of nodes: " 
                   << Node_AsymComm->Num_Export_to_Proc(i) << "\n";
          ContactTopologyEntity<Real>** nodes_to_print = Node_AsymComm->Export_Entity_List(i);
          for ( j = 0; j < Node_AsymComm->Num_Export_to_Proc(i); ++j ) {
            ContactNode<Real> * node = 
              static_cast<ContactNode<Real> *>(nodes_to_print[j]);
            postream << "       node " << j+1 
                     << "  " << node->Global_ID() << "\n";
          }
        }
        postream << "  number of import communication partners: " 
                 << Node_AsymComm->Num_Import_Comm_Partners() << "\n";
        for ( i = 0; i < Node_AsymComm->Num_Import_Comm_Partners(); ++i ) {
          postream << "    export communication partner id: " 
                   << Node_AsymComm->Import_Comm_Proc_ID(i) << "\n";
          postream << "    number of nodes: " 
                   << Node_AsymComm->Num_Import_from_Proc(i) << "\n";
          ContactTopologyEntity<Real>** nodes_to_print = Node_AsymComm->Import_Entity_List(i);
          for ( j = 0; j < Node_AsymComm->Num_Import_from_Proc(i); ++j ) {
            ContactNode<Real> * node = 
              static_cast<ContactNode<Real> *>(nodes_to_print[j]);
            postream << "       node " << j+1 
                     << "  " << node->Global_ID() << "\n";
          }
        }
        postream.flush();
      }
    }
      
    if (number_of_global_edges && contact_number_of_processors(SearchComm)>1) {
      int i, j, k;
      if (Edge_SymComm) {
        postream << "Edge symmetric communication list" << "\n";
        postream << "  number of communication partners: " 
                 << Edge_SymComm->Num_Comm_Partners() << "\n";
        for ( i = 0; i < Edge_SymComm->Num_Comm_Partners(); ++i ) {
          postream << "    communication partner id: " 
                   << Edge_SymComm->Comm_Proc_ID(i) << "\n";
          postream << "    number of edges: " 
                   << Edge_SymComm->Num_to_Proc(i) << "\n";
          ContactTopologyEntity<Real>** edges_to_print = Edge_SymComm->Entity_List(i);
          for ( j = 0; j < Edge_SymComm->Num_to_Proc(i); ++j ) {
            postream << "       edge " << j+1 << "\n";
            ContactEdge<Real> * edge = 
              static_cast<ContactEdge<Real> *>(edges_to_print[j]);
            int num_nodes = edge->Nodes_Per_Edge();
            for ( k = 0; k < num_nodes; ++k ) {
              postream << "         node " << k+1 << ":";
              postream << edge->Node(k)->Global_ID() << "\n";
            }
          }
        }
        postream.flush();
      }
    }
  }
#endif
}

void 
ContactTopology::Display_Ghosted_Entities( ContactParOStream& postream, int comm_flag )
{
  postream.flush();std::cout<<std::flush;
  int myproc = contact_processor_number(SearchComm);

  if (myproc==0) {
    std::cout<<"\n\nNumber of node blocks = "<<number_of_node_blocks<<std::endl;
  }
  for (int i=0; i<number_of_node_blocks; ++i) {
    postream<<"\n========== Ghosted Node Block "<<i<<", size = "
          <<ghosted_node_blocks[i]->NodeList()->NumEntities()<<") ==========\n";
    ghosted_node_blocks[i]->NodeList()->Display(postream);
  }
  postream.flush();
  
  if (myproc==0) {
    std::cout<<"\n\nNumber of face blocks = "<<number_of_face_blocks<<std::endl;
  }
  for (int i=0; i<number_of_face_blocks; ++i) {
    postream<<"\n========== Ghosted Face Block "<<i<<", size = "
          <<ghosted_face_blocks[i]->FaceList()->NumEntities()<<") ==========\n";
    ghosted_face_blocks[i]->FaceList()->Display(postream);
  }
  postream.flush();
  
  if (myproc==0) {
    std::cout<<"\n\nNumber of element blocks = "<<number_of_element_blocks<<std::endl;
  }
  for (int i=0; i<number_of_element_blocks; ++i) {
    postream<<"\n========== Ghosted Element Block "<<i<<", size = "
          <<ghosted_element_blocks[i]->ElemList()->NumEntities()<<") ==========\n";
    ghosted_element_blocks[i]->ElemList()->Display(postream);
  }
  postream.flush();

#ifndef CONTACT_NO_MPI
  if (comm_flag && GhostingCommSpec) {
    ContactZoltanLID zoltanLID;
    ContactZoltanGID zoltanGID;
    postream << "Ghosting asymmetric communication list" << "\n";
    int       num_import   = GhostingCommSpec->Num_Import();
    LB_ID_PTR import_lids  = GhostingCommSpec->Import_LIDS();
    LB_ID_PTR import_gids  = GhostingCommSpec->Import_GIDS();
    int*      import_procs = GhostingCommSpec->Import_Procs();
    int       num_export   = GhostingCommSpec->Num_Export();
    LB_ID_PTR export_lids  = GhostingCommSpec->Export_LIDS();
    LB_ID_PTR export_gids  = GhostingCommSpec->Export_GIDS();
    int*      export_procs = GhostingCommSpec->Export_Procs();
    if (export_gids!=NULL && export_lids!=NULL && export_procs!=NULL) {
      postream << "GHoSTING EXPORTS:\n";
      for (int i=0; i< num_export; ++i){
	postream << "  " << i << ":  GID = "
		 << "T" << zoltanGID.Type(&export_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "H" << zoltanGID.Hi(&export_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "L" << zoltanGID.Lo(&export_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "    LID = "
		 << "T" << zoltanLID.Type(&export_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "I" << zoltanLID.Index(&export_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "    PROC = "
		 << export_procs[i]<< "\n";
      }
    }
    if (import_gids!=NULL && import_lids!=NULL && import_procs!=NULL) {
      postream << "GHOSTING IMPORTS:\n";
      for (int i=0; i< num_import; ++i){
	postream << "  " << i << ":  GID = "
		 << "T" << zoltanGID.Type(&import_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "H" << zoltanGID.Hi(&import_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "L" << zoltanGID.Lo(&import_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "    LID = "
		 << "T" << zoltanLID.Type(&import_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "H" << zoltanLID.Index(&import_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "    PROC = "
		 << import_procs[i]<< "\n";
      }
    }
    postream.flush();
  }
#endif
}

void 
ContactTopology::Display_NodeNode_Interactions( ContactParOStream& postream, 
                                                int state )
{
  int cnt = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    Nodes[i]->Display_NodeNode_Interactions(postream, state);
    cnt += Nodes[i]->Number_NodeNode_Interactions(state);
  }
  if (!cnt) {
    postream << "No interactions present\n";
  }
}

void 
ContactTopology::Display_NodeNode_Interactions_Summary( ContactParOStream& postream, 
                                                        char* margin, int state )
{
  int cnt = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    cnt += Nodes[i]->Number_NodeNode_Interactions(state);
  }
  postream << margin << "Found "<<cnt<<" node/node interactions\n";
}

void 
ContactTopology::Display_NodeEntity_Interactions( ContactParOStream& postream, 
                                                  int state )
{
  int cnt = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    Nodes[i]->Display_NodeEntity_Interactions(postream, state);
    cnt += Nodes[i]->Number_NodeEntity_Interactions(state);
  }
  if (!cnt) {
    postream << "No interactions present\n";
  }
}

void 
ContactTopology::Display_NodeEntity_Interactions_Summary( ContactParOStream& postream, 
                                                          unsigned int status, char* margin, int state )
{
  int cnt0 = 0;
  int cnt1 = 0;
  int cnt2 = 0;
  int cnt3 = 0;
  int cnt4 = 0;
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    if (node->CheckContext(status) && node->Ownership() == ContactTopologyEntity<Real>::OWNED) {
      //PRECONDITION(node->Ownership() == ContactTopologyEntity<Real>::OWNED);
      cnt0 += node->Number_NodeEntity_Interactions(state);
      cnt1 += node->Number_NodeFace_Interactions(state);
      cnt2 += node->Number_NodeSurface_Interactions(state);
      cnt3 += node->Num_Tracked_Interactions(state);
      cnt4 += node->Num_Tied_Interactions(state);
    }
  }
  if (number_of_analytic_surfaces>0) {
    postream << margin << "Found "<<cnt1<<" node/face and "<<cnt2<<" node/surface interactions\n";
  } else {
    int cnt1_global = contact_global_sum(cnt1,SearchComm);
    postream << margin << "Found "<<cnt1<<" of "<<cnt1_global<<" total node/face interactions\n";
  }
  postream << margin << "      "<<cnt3<<" tracked interactions\n";
  postream << margin << "      "<<cnt4<<" tied interactions\n";
}

void 
ContactTopology::Display0_NodeEntity_Interactions_Summary( unsigned int status, char* margin, int state )
{
  int cnt1 = 0;
  int cnt2 = 0;
  int cnt3 = 0;
  int cnt4 = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
    if (node->Ownership() == ContactTopologyEntity<Real>::OWNED) {
//      PRECONDITION(node->Ownership() == ContactTopologyEntity<Real>::OWNED);
      cnt1 += node->Number_NodeFace_Interactions(state);
      cnt2 += node->Number_NodeSurface_Interactions(state);
      cnt3 += node->Num_Tracked_Interactions(state);
      cnt4 += node->Num_Tied_Interactions(state);
    }
  }
  int cnt_nfi = contact_global_sum(cnt1,SearchComm);
  int cnt_nsi = contact_global_sum(cnt2,SearchComm);
  int cnt_3   = contact_global_sum(cnt3,SearchComm);
  int cnt_4   = contact_global_sum(cnt4,SearchComm);
  if( contact_processor_number(SearchComm) == 0 ) {
    if (number_of_analytic_surfaces>0) {
      std::cout << margin << "Found "<<cnt_nfi<<" node/face and "<<cnt_nsi<<" node/surface interactions\n";
    } else {
      std::cout << margin << "Found "<<cnt_nfi<<" node/face interactions\n";
    }
    std::cout << margin << "      "<<cnt_3<<" tracked interactions\n";
    std::cout << margin << "      "<<cnt_4<<" tied interactions\n";

  }
}

void 
ContactTopology::Display_FaceFace_Interactions( ContactParOStream& postream, 
                                                int state )
{
  int cnt = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    Faces[i]->Display_FaceFace_Interactions(postream, state);
    cnt += Faces[i]->Number_FaceFace_Interactions(state);
  }
  if (!cnt) {
    postream << "No interactions present\n";
  }
}

void 
ContactTopology::Display_FaceFace_Interactions_Summary( ContactParOStream& postream, 
                                                        char* margin, int state )
{
  int cnt = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    cnt += Faces[i]->Number_FaceFace_Interactions();
  }
  postream << margin << "Found "<<cnt<<" face/face interactions\n";
}

void 
ContactTopology::Display_FaceCoverage_Interactions( ContactParOStream& postream, 
                                                    int state )
{
  int cnt = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    Faces[i]->Display_FaceCoverage_Interactions(postream, state);
    cnt += Faces[i]->Number_FaceCoverage_Interactions();
  }
  if (!cnt) {
    postream << "No interactions present\n";
  }
}

void 
ContactTopology::Display_FaceCoverage_Interactions_Summary( ContactParOStream& postream, 
                                                            char* margin, int state )
{
  int cnt = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    cnt += Faces[i]->Number_FaceCoverage_Interactions(state);
  }
  postream << margin << "Found "<<cnt<<" face/coverage interactions\n";
}

void 
ContactTopology::Display_ElementElement_Interactions( ContactParOStream& postream, 
                                                      int state )
{
  int cnt = 0;
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (int i=0; i<number_of_elements; ++i) {
    Elements[i]->Display_ElementElement_Interactions(postream, state);
    cnt += Elements[i]->Number_ElementElement_Interactions(state);
  }
  if (!cnt) {
    postream << "No interactions present\n";
  }
}

void 
ContactTopology::Display_ElementElement_Interactions_Summary( ContactParOStream& postream, 
                                                              char* margin, int state )
{
  int cnt = 0;
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (int i=0; i<number_of_elements; ++i) {
    cnt += Elements[i]->Number_ElementElement_Interactions(state);
  }
  postream << margin << "Found "<<cnt<<" element/element interactions\n";
}

void
ContactTopology::Set_Up_Variable_Handles()
{
  Var_Handles = new VariableHandle[MAX_VARIABLE_ENUM];

  int i;
  int index = 0;
  int var_handle_offset = 0;
  
  // The first locations in Var_Handles are the Node_Scalar_Vars
  for( i=0 ; i<ContactNode<Real>::NUMBER_SCALAR_VARS ; ++i )
    Var_Handles[var_handle_offset+i] = index++;
  var_handle_offset += ContactNode<Real>::NUMBER_SCALAR_VARS;
  // The next locations in Var_Handles are the Node_Vector_Vars
  for( i=0 ; i<ContactNode<Real>::NUMBER_VECTOR_VARS ; ++i ){
    Var_Handles[var_handle_offset+i] = index;
    index += 3;
  }
  var_handle_offset += ContactNode<Real>::NUMBER_VECTOR_VARS;

  // The next locations are the Edge_Scalar_Vars
  index = 0;
  for( i=0 ; i<ContactEdge<Real>::NUMBER_SCALAR_VARS ; ++i )
    Var_Handles[var_handle_offset+i] = index++;
  var_handle_offset += ContactEdge<Real>::NUMBER_SCALAR_VARS;
  // The next locations in Var_Handles are the Edge_Vector_Vars
  for( i=0 ; i<ContactEdge<Real>::NUMBER_VECTOR_VARS ; ++i ){
    Var_Handles[var_handle_offset+i] = index;
    index += 3;
  }
  var_handle_offset += ContactEdge<Real>::NUMBER_VECTOR_VARS;

  // The next locations are the Face_Scalar_Vars
  index = 0;
  for( i=0 ; i<ContactFace<Real>::NUMBER_SCALAR_VARS ; ++i )
    Var_Handles[var_handle_offset+i] = index++;
  var_handle_offset += ContactFace<Real>::NUMBER_SCALAR_VARS;
  // The next locations in Var_Handles are the Face_Vector_Vars
  for( i=0 ; i<ContactFace<Real>::NUMBER_VECTOR_VARS ; ++i ){
    Var_Handles[var_handle_offset+i] = index;
    index += 3;
  }
  var_handle_offset += ContactFace<Real>::NUMBER_VECTOR_VARS;

  // The next location are the Elem_Scalar_Vars
  index = 0;
  for( i=0 ; i<ContactElem<Real>::NUMBER_SCALAR_VARS ; ++i )
    Var_Handles[var_handle_offset+i] = index++;
  var_handle_offset += ContactElem<Real>::NUMBER_SCALAR_VARS;
  // The next locations in Var_Handles are the Elem_Vector_Vars
  for( i=0 ; i<ContactElem<Real>::NUMBER_VECTOR_VARS ; ++i ){
    Var_Handles[var_handle_offset+i] = index;
    index += 3;
  }
  var_handle_offset += ContactElem<Real>::NUMBER_VECTOR_VARS;

  // The next location are the Element_Scalar_Vars
  index = 0;
  for( i=0 ; i<ContactElement::NUMBER_SCALAR_VARS ; ++i )
    Var_Handles[var_handle_offset+i] = index++;
  var_handle_offset += ContactElement::NUMBER_SCALAR_VARS;
  // The next locations in Var_Handles are the Element_Vector_Vars
  for( i=0 ; i<ContactElement::NUMBER_VECTOR_VARS ; ++i ){
    Var_Handles[var_handle_offset+i] = index;
    index += 3;
  }
  var_handle_offset += ContactElement::NUMBER_VECTOR_VARS;
  
  POSTCONDITION( var_handle_offset == MAX_VARIABLE_ENUM );

#define    NODE_SCALAR_VAR( b,a ) a = *(Var_Handles+b);
#define    NODE_VECTOR_VAR( b,a ) a = *(Var_Handles+b);
#define    EDGE_SCALAR_VAR( b,a ) a = *(Var_Handles+b);
#define    EDGE_VECTOR_VAR( b,a ) a = *(Var_Handles+b);
#define    FACE_SCALAR_VAR( b,a ) a = *(Var_Handles+b);
#define    FACE_VECTOR_VAR( b,a ) a = *(Var_Handles+b);
#define    ELEM_SCALAR_VAR( b,a ) a = *(Var_Handles+b);
#define    ELEM_VECTOR_VAR( b,a ) a = *(Var_Handles+b);
#define ELEMENT_SCALAR_VAR( b,a ) a = *(Var_Handles+b);
#define ELEMENT_VECTOR_VAR( b,a ) a = *(Var_Handles+b);
#define     NEI_SCALAR_VAR( b,a )
#define     NEI_VECTOR_VAR( b,a )
#define     NNI_SCALAR_VAR( b,a )
#define     EEI_SCALAR_VAR( b,a )
#include "contact_variables.def"
#undef    NODE_SCALAR_VAR
#undef    NODE_VECTOR_VAR
#undef    EDGE_SCALAR_VAR
#undef    EDGE_VECTOR_VAR
#undef    FACE_SCALAR_VAR
#undef    FACE_VECTOR_VAR
#undef    ELEM_SCALAR_VAR
#undef    ELEM_VECTOR_VAR
#undef ELEMENT_SCALAR_VAR
#undef ELEMENT_VECTOR_VAR
#undef     NEI_SCALAR_VAR
#undef     NEI_VECTOR_VAR
#undef     NNI_SCALAR_VAR
#undef     EEI_SCALAR_VAR
}



bool ContactTopology::Faces_Connected( ContactFace<Real>* face1, ContactFace<Real>* face2 )
{
  PRECONDITION( face1 && face2 );

  for( int i=0 ; i<face1->Nodes_Per_Face() ; ++i ){
    ContactNode<Real>* node = face1->Node(i);
    for( int j=0 ; j<face2->Nodes_Per_Face() ; ++j )
      if( node == face2->Node(j) ) return true;
  }
  return false;
}


void ContactTopology::Update_State()
{
  int i;

  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->Update_Interactions();
  }

  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->Update_Interactions();
  }

  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (i=0; i<number_of_elements; ++i) {
    Elements[i]->Update_Interactions();
  }
}

void ContactTopology::Delete_All_Interactions()
{
  int i, j;
  int number_of_states = 2;
  ContactInteractionEntity<Real>* link;
  ContactInteractionDLL<Real>* interactions;
  
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    for (j=0; j<number_of_states; ++j) {
      node->Delete_NodeEntity_Interactions(j);
      node->Delete_NodeNode_Interactions(j);
    }
  }
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    for (j=0; j<number_of_states; ++j) {
      interactions = face->Get_FaceFace_Interactions(j);
      if(interactions != NULL) {
        interactions->IteratorStart();
        while ((link = interactions->IteratorForward())) {
  	  ContactFaceFaceInteraction<Real>* cffi = 
  	         static_cast<ContactFaceFaceInteraction<Real>*>(link);
  	  cffi->~ContactFaceFaceInteraction();
          search->Get_Allocators()[ContactSearch::ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi); // PJSA
        }
        interactions->Clear();
      }
      interactions = face->Get_FaceCoverage_Interactions(j);
      if(interactions != NULL) {
        interactions->IteratorStart();
        while ((link = interactions->IteratorForward())) {
  	  ContactFaceCoverageInteraction* cfci = 
  	         static_cast<ContactFaceCoverageInteraction*>(link);
  	  cfci->~ContactFaceCoverageInteraction();
        }
        interactions->Clear();
      }
    }
    face->Clear_FaceFace_Interactions(); // PJSA
  }
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for( i=0; i<number_of_elements; ++i) {
    ContactElement* element = Elements[i];
    for (j=0; j<number_of_states; ++j) {
      interactions = element->Get_ElementElement_Interactions(j);
      interactions->IteratorStart();
      while ((link = interactions->IteratorForward())) {
   	ContactElementElementInteraction* ceei = 
   	       static_cast<ContactElementElementInteraction*>(link);
   	ceei->~ContactElementElementInteraction();
      }
      interactions->Clear();
    }
  }
  search->Get_Allocators()[ContactSearch::ALLOC_ContactFaceFaceInteraction].Purge(); // PJSA
}

int ContactTopology::Number_NodeFace_Interactions()
{
  int n = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    n += Nodes[i]->Number_NodeFace_Interactions();
  }
  return n;
}

int ContactTopology::Number_NodeEntity_Interactions()
{
  int n = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    n += Nodes[i]->Number_NodeEntity_Interactions();
  }
  return n;
}

int ContactTopology::Number_NodeNode_Interactions()
{
  int n = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    n += Nodes[i]->Number_NodeNode_Interactions();
  }
  return n;
}

void ContactTopology::Size_NodeFace_Interactions( int& num_interactions,
						  int& data_size )
{
  data_size = SIZE_NODEFACE_INTERACTION_DATA;
  num_interactions = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    num_interactions += Nodes[i]->Number_NodeFace_Interactions();
  }
}

void ContactTopology::Get_NodeFace_Interactions(int* Node_block_ids,
						int* node_indexes_in_block,
						int* node_entity_keys,
						int* face_block_ids,  
						int* face_indexes_in_block,
						int* face_proc,
						Real* interaction_data )
{
  int kk=0;
  const int data_size = SIZE_NODEFACE_INTERACTION_DATA;

  // pack a dense array of interaction data
  // assumes host code gives enough space for data array
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (int j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
      if (interactions[j]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
      ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[j]);
      // all indexes are incremented to get back to FORTRAN indexing
      Node_block_ids[kk] = cnfi->Node()->BlockID()+1;
      node_indexes_in_block[kk] = cnfi->Node()->HostArrayIndex()+1;
      node_entity_keys[kk] = (int)
   	cnfi->Scalar_Var(ContactNodeFaceInteraction::NODE_ENTITY_KEY)+1;
      face_block_ids[kk] = cnfi->FaceEntityData()->block_id+1;
      face_indexes_in_block[kk] = cnfi->FaceEntityData()->index_in_host_array+1;
      face_proc[kk] = cnfi->FaceEntityData()->owner;
      Real* xi = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
      interaction_data[NF_COORD1+kk*data_size] = xi[0];
      interaction_data[NF_COORD2+kk*data_size] = xi[1];
      interaction_data[NF_GAP_CUR+kk*data_size] = 
   	cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) +
   	cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
      // interaction_data[NF_GAP_OLD+kk*data_size] =
      //   cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
      Real* pb_dir = 
   	cnfi->Get_Pushback_Dir();
      interaction_data[NF_PUSHBACK_X+kk*data_size] = pb_dir[0];
      interaction_data[NF_PUSHBACK_Y+kk*data_size] = pb_dir[1];
      interaction_data[NF_PUSHBACK_Z+kk*data_size] = pb_dir[2];
      Real* surf_norm = 
   	cnfi->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
      interaction_data[NF_NORMAL_X+kk*data_size] = surf_norm[0];
      interaction_data[NF_NORMAL_Y+kk*data_size] = surf_norm[1];
      interaction_data[NF_NORMAL_Z+kk*data_size] = surf_norm[2];
      interaction_data[NF_SOURCE+kk*data_size] = 
   	cnfi->Scalar_Var(ContactNodeFaceInteraction::SOURCE);
      interaction_data[NF_NODE_AREA+kk*data_size] =
   	cnfi->Scalar_Var(ContactNodeFaceInteraction::NODE_AREA);
      interaction_data[NF_GAP_INT+kk*data_size] =
   	cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_INT);
      ++kk;
    }
  }
}

void ContactTopology::Size_NodeNode_Interactions( int& num_interactions,
						  int& data_size )
{
  data_size = SIZE_NODENODE_INTERACTION_DATA;
  num_interactions = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    num_interactions += Nodes[i]->Number_NodeNode_Interactions();
  }
}

void ContactTopology::Get_NodeNode_Interactions(int* slave_node_block_ids,
						int* slave_node_indexes_in_block,
						int* master_node_block_ids,  
						int* master_node_indexes_in_block,
						int* master_node_proc,
						Real* interaction_data )
{
  int index = 0;
  const int data_size = SIZE_NODENODE_INTERACTION_DATA;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    ContactInteractionDLL<Real>* interactions = Nodes[i]->Get_NodeNode_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
        ContactNodeNodeInteraction* cnni =
    	   static_cast<ContactNodeNodeInteraction*> (interaction);
        slave_node_block_ids[index] = cnni->SlaveNode()->BlockID()+1;
        slave_node_indexes_in_block[index] = cnni->SlaveNode()->HostArrayIndex()+1;
        master_node_block_ids[index] = cnni->MasterNodeEntityData()->block_id+1;
        master_node_indexes_in_block[index] = cnni->MasterNodeEntityData()->index_in_host_array+1;
        master_node_proc[index] = cnni->MasterNodeEntityData()->owner;
        interaction_data[NN_DIST+index*data_size] = cnni->Scalar_Var(ContactNodeNodeInteraction::DISTANCE);
        ++index;
      }
    }
  }
}

int ContactTopology::Number_NodeSurface_Interactions()
{
  int n = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    n += Nodes[i]->Number_NodeSurface_Interactions();
  }
  return n;
}

void ContactTopology::Size_NodeSurface_Interactions( int& num_interactions,
						     int& data_size )
{
  data_size = SIZE_NODESURFACE_INTERACTION_DATA;
  num_interactions = 0;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    num_interactions += Nodes[i]->Number_NodeSurface_Interactions();
  }
}

void ContactTopology::Get_NodeSurface_Interactions( int* Node_block_ids,
						    int* node_indexes_in_block,
						    int* surface_ID,
						    Real* interaction_data )
{
  const int data_size = SIZE_NODESURFACE_INTERACTION_DATA;
  int index=0;

  // pack a dense array of interaction data
  // assumes host code gives enough space for data array
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (int j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
      if (interactions[j]->Get_Type()!=ContactNodeEntityInteraction::NODE_SURFACE_INTERACTION) continue;
      ContactNodeSurfaceInteraction* cnsi = static_cast<ContactNodeSurfaceInteraction*>(interactions[j]);
      Node_block_ids[index] = cnsi->Node()->BlockID()+1;
      node_indexes_in_block[index] = cnsi->Node()->HostArrayIndex()+1;
      surface_ID[index] = cnsi->Surface()->ProcArrayIndex()+1;
      Real* contact_point = 
  	cnsi->Vector_Var(ContactNodeSurfaceInteraction::CONTACT_POINT);
      interaction_data[NS_COORD_X+index*data_size] = contact_point[0];
      interaction_data[NS_COORD_Y+index*data_size] = contact_point[1];
      interaction_data[NS_COORD_Z+index*data_size] = contact_point[2];
      interaction_data[NS_GAP_CUR+index*data_size] = 
  	cnsi->Scalar_Var(ContactNodeSurfaceInteraction::GAP_CUR);
      // interaction_data[NS_GAP_OLD+index*data_size] = 
      //   cnsi->Scalar_Var(ContactNodeSurfaceInteraction::GAP_OLD);
      Real* surf_normal = 
  	cnsi->Get_Normal();
      interaction_data[NS_NORMAL_X+index*data_size] = surf_normal[0];
      interaction_data[NS_NORMAL_Y+index*data_size] = surf_normal[1];
      interaction_data[NS_NORMAL_Z+index*data_size] = surf_normal[2];
      Real* pf_normal = 
  	cnsi->Vector_Var(ContactNodeSurfaceInteraction::PHYSICAL_FACE_NORMAL);
      interaction_data[NS_PFNORMAL_X+index*data_size] = pf_normal[0];
      interaction_data[NS_PFNORMAL_Y+index*data_size] = pf_normal[1];
      interaction_data[NS_PFNORMAL_Z+index*data_size] = pf_normal[2];
      interaction_data[NS_ENTITY_KEY+index*data_size] = 
  	cnsi->Scalar_Var(ContactNodeSurfaceInteraction::NODE_ENTITY_KEY);
      interaction_data[NS_NDAREA+index*data_size] = 
  	cnsi->Scalar_Var(ContactNodeSurfaceInteraction::NODE_AREA);
      ++index;
    }
  }
}

int ContactTopology::Number_FaceFace_Interactions()
{
  int n = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( int i=0; i<number_of_faces; ++i) {
    n += Faces[i]->Number_FaceFace_Interactions();
  }
  return n;
}

void ContactTopology::Size_FaceFace_Interactions( int& num_interactions,
						  int& data_size )
{
  data_size = 0;
  num_interactions = 0;
  
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions == NULL) continue;
    interactions->IteratorStart();
    while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
      ContactFaceFaceInteraction<Real>* cffi =
  	static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
      data_size += cffi->Data_Size();
      ++num_interactions;
    }
  }
}

void ContactTopology::Get_FaceFace_Interactions( int* slave_face_block_ids,
						 int* slave_face_indexes_in_block,
						 int* slave_face_proc,
						 int* master_face_block_ids,
						 int* master_face_indexes_in_block,
						 int* master_face_proc,
                                                 int* interaction_index,
						 Real* interaction_data )
{
  int index0 = 0;
  int index1 = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  ContactSearch::Search_Option_Status compute_partials;
  Real order;
  search->Get_Search_Option(ContactSearch::COMPUTE_PARTIALS, compute_partials, &order);
  for( int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions == NULL) continue;
    interactions->IteratorStart();
    while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
      ContactFaceFaceInteraction<Real>* cffi =
  	 static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
      slave_face_block_ids[index0] = cffi->SlaveFace()->BlockID()+1;
      slave_face_indexes_in_block[index0] = cffi->SlaveFace()->HostArrayIndex()+1;
      slave_face_proc[index0] = cffi->SlaveFaceEntityData()->owner;
      master_face_block_ids[index0] = cffi->MasterFaceEntityData()->block_id+1;
      master_face_indexes_in_block[index0] = cffi->MasterFaceEntityData()->index_in_host_array+1;
      master_face_proc[index0] = cffi->MasterFaceEntityData()->owner;
      interaction_index[index0] = index1++;
      *interaction_data++ = cffi->NumEdges();
      ContactFaceFaceVertex<Real>* vertices = cffi->Get_Vertices();
      int j;
      for (j=0; j<cffi->NumEdges(); ++j) {
  	if (vertices[j].master_edge_flag) {
  	  *interaction_data++ = j+1;
  	} else {
  	  *interaction_data++ = 0;
  	}
      }
      index1 += cffi->NumEdges();
      for (j=0; j<cffi->NumEdges(); ++j) {
  	if (vertices[j].master_edge_flag) {
  	  *interaction_data++ = 0;
  	} else {
  	  *interaction_data++ = vertices[j].slave_edge_id;
  	}
      }
      index1 += cffi->NumEdges();
      for (j=0; j<cffi->NumEdges(); ++j) {
  	*interaction_data++ = vertices[j].slave_x;
  	*interaction_data++ = vertices[j].slave_y;
  	*interaction_data++ = vertices[j].master_x;
  	*interaction_data++ = vertices[j].master_y;
  	index1 += 4;
      }
      if(compute_partials == ContactSearch::ACTIVE && order > 0) {
        int k;
        for (j=0; j<cffi->NumEdges(); ++j) {
          for (k=0; k<cffi->NumDerivatives(); ++k)
            *interaction_data++ = vertices[j].slave_x_derivatives[k];
          for (k=0; k<cffi->NumDerivatives(); ++k) 
            *interaction_data++ = vertices[j].slave_y_derivatives[k];
          for (k=0; k<cffi->NumDerivatives(); ++k) 
            *interaction_data++ = vertices[j].master_x_derivatives[k];
          for (k=0; k<cffi->NumDerivatives(); ++k) 
            *interaction_data++ = vertices[j].master_y_derivatives[k];
          index1 += 4*cffi->NumDerivatives();
        }
        if(order == 2) {
          for (j=0; j<cffi->NumEdges(); ++j) {
            for (k=0; k<cffi->NumSecondDerivatives(); ++k)
              *interaction_data++ = vertices[j].slave_x_second_derivatives[k];
            for (k=0; k<cffi->NumSecondDerivatives(); ++k)
              *interaction_data++ = vertices[j].slave_y_second_derivatives[k];
            for (k=0; k<cffi->NumSecondDerivatives(); ++k)
              *interaction_data++ = vertices[j].master_x_second_derivatives[k];
            for (k=0; k<cffi->NumSecondDerivatives(); ++k)
              *interaction_data++ = vertices[j].master_y_second_derivatives[k];
            index1 += 4*cffi->NumSecondDerivatives();
          }
        }
      }
      ++index0;
    }
  }
}

int ContactTopology::Number_FaceCoverage_Interactions()
{
  int n = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( int i=0; i<number_of_faces; ++i) {
    n += Faces[i]->Number_FaceCoverage_Interactions();
  }
  return n;
}

void ContactTopology::Size_FaceCoverage_Interactions( int& num_interactions,
						      int& data_size )
{
  data_size = 0;
  num_interactions = 0;
  
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
        ContactFaceCoverageInteraction* cfci =
  	  static_cast<ContactFaceCoverageInteraction*> (interaction);
        data_size += cfci->Data_Size();
        ++num_interactions;
      }
    }
  }
}

void ContactTopology::Get_FaceCoverage_Interactions( int* face_block_ids,
						     int* face_indexes_in_block,
                                                     int* interaction_index,
						     Real* interaction_data )
{
  int index0 = 0;
  int index1 = 0;
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
        ContactFaceCoverageInteraction* cfci =
  	     static_cast<ContactFaceCoverageInteraction*> (interaction);
        face_block_ids[index0] = cfci->SlaveFace()->BlockID()+1;
        face_indexes_in_block[index0] = cfci->SlaveFace()->HostArrayIndex()+1;
        interaction_index[index0] = index1++;
        *interaction_data++ = cfci->NumVertices();
        for( ContactFaceCoverageVertex* ll_node=cfci->Head(); 
  	     ll_node; ll_node=ll_node->next ){
  	  *interaction_data++ = ll_node->slave_x;
  	  *interaction_data++ = ll_node->slave_y;
  	  index1 += 2;
        }
        ++index0;
      }
    }
  }
}

int ContactTopology::Number_ElementElement_Interactions()
{
  int n = 0;
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for( int i=0; i<number_of_elements; ++i) {
    n += Elements[i]->Number_ElementElement_Interactions();
  }
  return n;
}

void ContactTopology::Size_ElementElement_Interactions( int& num_interactions,
						        int& data_size )
{
  data_size = SIZE_ELEMENTELEMENT_INTERACTION_DATA;
  num_interactions = 0;
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for( int i=0; i<number_of_elements; ++i) {
    num_interactions += Elements[i]->Number_Interactions();
  }
}

void ContactTopology::Get_ElementElement_Interactions( int* slave_element_block_ids,
						       int* slave_element_indexes_in_block,
						       int* master_element_block_ids,
						       int* master_element_indexes_in_block,
						       int* master_element_proc,
						       Real* interaction_data )
{
  int index = 0;
  const int data_size = SIZE_ELEMENTELEMENT_INTERACTION_DATA;
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for( int i=0; i<number_of_elements; ++i) {
    ContactElement* element = Elements[i];
    ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
    interactions->IteratorStart();
    while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
      ContactElementElementInteraction* ceei =
    	 static_cast<ContactElementElementInteraction*> (interaction);
      slave_element_block_ids[index] = ceei->SlaveElement()->BlockID()+1;
      slave_element_indexes_in_block[index] = ceei->SlaveElement()->HostArrayIndex()+1;
      master_element_block_ids[index] = ceei->MasterElementEntityData()->block_id+1;
      master_element_indexes_in_block[index] = ceei->MasterElementEntityData()->index_in_host_array+1;
      master_element_proc[index] = ceei->MasterElementEntityData()->owner;
      interaction_data[EE_VOLUME+index*data_size] = ceei->Scalar_Var(ContactElementElementInteraction::VOLUME);
      ++index;
    }
  }
}

void ContactTopology::Compute_Owners(ContactSearch::ContactErrorCode& error_code){
  //
  // this routine computes the owners for all local nodes
  //

  // at first, set all objects to owned and also their global ids 
  int i;
  int my_proc = contact_processor_number(SearchComm);

  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    Nodes[i]->Owner(my_proc);
    Nodes[i]->Secondary_Owner(my_proc);
  }

  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
  for (i=0; i<number_of_edges; ++i) {
    Edges[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    Edges[i]->Owner(my_proc);
    Edges[i]->Secondary_Owner(my_proc);
  }

  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    Faces[i]->Owner(my_proc);
    Faces[i]->Secondary_Owner(my_proc);
  }

  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (i=0; i<number_of_elements; ++i) {
    Elements[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    Elements[i]->Owner(my_proc);
    Elements[i]->Secondary_Owner(my_proc);
  }

#ifndef CONTACT_NO_MPI
  // return if in serial  
  if ( 1 == contact_number_of_processors( SearchComm ) ) return; 

  // compute owners for nodes
  Compute_Owners_For_Entity(Node_SymComm, node_list, 0);

  // compute edge communication lists
  Compute_Edge_Comm_List(error_code);
  
  // compute owners for edges
  Compute_Owners_For_Entity(Edge_SymComm, edge_list, 1);
  edge_list->Rehash();
#endif  
}

#ifndef CONTACT_NO_MPI
void ContactTopology::Compute_Owners_For_Entity(
                      ContactSymComm* comm_list,
                      ContactTopologyEntityList* entity_list,
                      int gid)
{
  int i,j; // counters
  int number_of_entities = entity_list->NumEntities();
  ContactTopologyEntity<Real>** Entity = entity_list->EntityList();
  if ( 0 == comm_list->Num_Comm_Partners() ) {
    // need this to balance the global sync in contact_swapadd_data_array()
    contact_global_sync(SearchComm);
    return;
  }
  
  //-------------------------------------------
  // compute the owner of all entities in entity_list

  // get local processor id

  const int local_proc_id = contact_processor_number(SearchComm);
  
  // first loop over all entities in the entity list and set temp tag to
  // the local processor id
  for (i=0; i<number_of_entities; ++i) {
    Entity[i]->temp_tag = local_proc_id;
  }
  
  // loop over the communication list
  int count = 0;
  for ( i = 0; i < comm_list->Num_Comm_Partners(); ++i ) {
    int proc_id = comm_list->Comm_Proc_ID(i);
    int num_shared_entities = comm_list->Num_to_Proc(i);
    if ( proc_id < local_proc_id ) {
      ContactTopologyEntity<Real> ** comm_entities = comm_list->Entity_List(i);
      for ( j = 0; j < num_shared_entities; ++j ) {
	ContactTopologyEntity<Real> * entity = comm_entities[j];
	if (entity->temp_tag > proc_id) entity->temp_tag = proc_id;
      }
    }
    count = count + num_shared_entities;
  }
  
  int stride = 2;
  if (gid) stride = 4;
  // loop over all entities in the entity list and set the owner flag
  i = 0;
  int * id_data = new int[stride*entity_list->NumEntities()];
  for (j=0; j<number_of_entities; ++j) {
    if ( local_proc_id == Entity[j]->temp_tag ){
      Entity[j]->Ownership(ContactTopologyEntity<Real>::OWNED);
      id_data[i  ] = Entity[j]->Owner();
      id_data[i+1] = Entity[j]->OwnerProcArrayIndex();
      if (gid) {
        id_data[i+2] = Entity[j]->Global_ID().HiInt();
        id_data[i+3] = Entity[j]->Global_ID().LoInt();
      }
    } else {
      Entity[j]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
      id_data[i  ] = 0;
      id_data[i+1] = 0;
      if (gid) {
        id_data[i+2] = 0;
        id_data[i+3] = 0;
      }
    }
    i += stride;
  }
  
  contact_swapadd_data_array(SearchComm,*comm_list,*comm_buffer,id_data,stride);

  i = 0;
  for (j=0; j<number_of_entities; ++j) {
    Entity[j]->Owner(id_data[i]);
    Entity[j]->OwnerProcArrayIndex(id_data[i+1]);
    if (gid) {
      Entity[j]->Global_ID().HiInt(id_data[i+2]);
      Entity[j]->Global_ID().LoInt(id_data[i+3]);
    }
    i += stride;
  }
  
  delete [] id_data;
  
}
#endif

#ifndef CONTACT_NO_MPI
void ContactTopology::Compute_Edge_Comm_List(
				  ContactSearch::ContactErrorCode& error_code)
{
  //allocate linked list to hold shared edges

  std::vector< std::pair<ContactTopologyEntity<Real>*, int> > shared_edge_list;

  //----------------------
  // find all the edges which are shared between processors
  Find_Shared_Edge_Candidates( &shared_edge_list);

  //----------------------
  // now we need to broadcast the edges (just their nodes)
  //   first transmit number of edges that we will communicate
  int mesg = 1003;
  int num_partners = Node_SymComm->Num_Comm_Partners();
  RequestHandle * recv_handles = new RequestHandle[num_partners];
  char * recv_num_edge_buff = new char[num_partners*sizeof(int)*2];
  char * send_num_edge_buff = new char[num_partners*sizeof(int)*2];
  int * rb, *recv_num_edges, *sb, *send_num_edges;
  rb = recv_num_edges = reinterpret_cast<int*>(recv_num_edge_buff);
  sb = send_num_edges = reinterpret_cast<int*>(send_num_edge_buff);
  for (int i=0;i<num_partners*2;++i) send_num_edges[i] = 0;
    
  //     register recvs for # of edges to recv
  for( int i=0 ; i<num_partners ; ++i ){
    recv_handles[i] = 
      contact_nonblocking_receive( mesg, rb, 2,
				   Node_SymComm->Comm_Proc_ID( i ), 
				   SearchComm );
    rb+=2;
  }
  //     sync to make sure all recvs have been posted
  contact_global_sync(SearchComm);
  //
  //     compute what to send -- loop over edges and count how many to 
  //       send to whom
  //

  for(int i = 0; i < shared_edge_list.size(); ++i) {
    std::pair<ContactTopologyEntity<Real>*, int> &entry = shared_edge_list[i];
    int dest = entry.second;
    for (int j = 0; j<num_partners; ++j) {
      int partner = Node_SymComm->Comm_Proc_ID(j);
      if ( partner == dest){
	send_num_edges[j*2]++;
	send_num_edges[j*2+1] += 
	  ((ContactEdge<Real>*)entry.first)->Nodes_Per_Edge();
	break;
      }
    }
  }

  //     post sends
  for (int i=0; i<num_partners; ++i) {
    contact_blocking_send( mesg, sb, 2,
			   Node_SymComm->Comm_Proc_ID( i ), 
			   SearchComm );
    sb+=2;
  }
  //     Wait till all of my messages have arrived
  for(int i=0 ; i<num_partners ; ++i ){
    contact_wait_msg_done( recv_handles[i] );
  }

  //   now transmit the nodes associated with the edges 
  //     NOTE: This assumes that global ids are two integers
  mesg = 1004;
  int len;
  int total_recv_len = 0;
  int * recv_offsets = new int[num_partners];
  for (int i = 0; i < num_partners; ++i){
    recv_offsets[i] = total_recv_len;
    total_recv_len += recv_num_edges[i*2+1]*2+recv_num_edges[i*2]; 
  }
  int total_send_len = 0;
  for (int i = 0; i < num_partners; ++i)
    total_send_len += send_num_edges[i*2+1]*2+send_num_edges[i*2];
  char * recv_edge_buff = new char[sizeof(int)*total_recv_len];
  char * send_edge_buff = new char[sizeof(int)*total_send_len];
  int * recv_edge, * send_edge;
  rb = recv_edge = reinterpret_cast<int*>(recv_edge_buff);
  sb = send_edge = reinterpret_cast<int*>(send_edge_buff);
  for (int i=0;i<total_send_len;++i) send_edge[i] = 0;
    
  // register recvs for edges to recv
  for(int i=0 ; i<num_partners ; ++i ){
    len = recv_num_edges[i*2] + recv_num_edges[i*2+1]*2;
    if ( len > 0) recv_handles[i] = 
		    contact_nonblocking_receive( mesg, rb, len, 
						 Node_SymComm->Comm_Proc_ID(i),
						 SearchComm );
    rb+=len;
  }
  // sync to make sure all recvs have been posted
  contact_global_sync(SearchComm);

  // compute what to send -- loop over edges store the edge nodal ids
  int index = 0;
  for (int i = 0; i<num_partners; ++i) {
    int partner = Node_SymComm->Comm_Proc_ID(i);

    for(int j = 0; j < shared_edge_list.size(); ++j) {
      std::pair<ContactTopologyEntity<Real>*,int> &entry = shared_edge_list[j];

      int dest = entry.second;
      if ( dest == partner) {
	ContactEdge<Real> * edge = (ContactEdge<Real>*) entry.first;
	int num_nodes = edge->Nodes_Per_Edge();
	send_edge[index] = num_nodes;
	++index;
	for (int k = 0; k < num_nodes; ++k) {
	  send_edge[index]   = edge->Node(k)->Global_ID().HiInt();
	  send_edge[index+1] = edge->Node(k)->Global_ID().LoInt();
	  index += 2;
	}
      }
    }
  }
  PRECONDITION(index == total_send_len);

  // post sends
  for (int i=0; i<num_partners; ++i) {
    len = send_num_edges[i*2] + send_num_edges[i*2+1]*2;
    if ( len > 0 ) contact_blocking_send( mesg, sb, len, 
					  Node_SymComm->Comm_Proc_ID( i ), 
					  SearchComm );
    sb +=len;
  }
  // Wait till all of my messages have arrived
  for(int i=0 ; i<num_partners ; ++i ){
    if (recv_num_edges[i*2] >0) contact_wait_msg_done( recv_handles[i] );
  }
  
  //----------------------
  // now everyone has the corresponding edge data. Now we need to compare
  // to build the communication list

  // use temp_tag to record the number of times the edge is linked with
  // faces in parallel. Here we initialize it to the number of faces 
  // locally connected to the edge.

  for(int i = 0; i < shared_edge_list.size(); ++i) {
    std::pair<ContactTopologyEntity<Real> *, int> &entry = shared_edge_list[i];
    ContactEdge<Real>* edge = (ContactEdge<Real>*) entry.first;
    edge->temp_tag = edge->Number_Face_Connections();    
  }

  
  for(int i = 0; i < shared_edge_list.size(); ++i) {
    std::pair<ContactTopologyEntity<Real> *, int> &entry = shared_edge_list[i];
    int dest = entry.second;
    int found = 0;
    ContactEdge<Real> * local_edge = (ContactEdge<Real>*) entry.first;
    for (int j = 0; j < num_partners; ++j){
      int partner = Node_SymComm->Comm_Proc_ID(j);
      if ( dest == partner) {
	int offset = recv_offsets[j];
	int num_edges = recv_num_edges[j*2];
	
	// now compare local edge against recv'ed edges	
	for (int k = 0; k < num_edges; ++k){
	  if ( local_edge->Nodes_Per_Edge() == recv_edge[offset]){
	    int found_node = 0;
	    for (int l = 0; l < local_edge->Nodes_Per_Edge(); ++l) {
	      int hi = local_edge->Node(l)->Global_ID().HiInt();
	      int lo = local_edge->Node(l)->Global_ID().LoInt();
	      for (int m = 0; m < recv_edge[offset]; ++m) {
		if (recv_edge[offset+m*2+1] == hi &&
		    recv_edge[offset+m*2+2] == lo) {
		  ++found_node;
		  break;
		}
	      } // end loop on nodes on recv'd edge
	    } // end loop on nodes on shared candidate
	    if ( found_node == local_edge->Nodes_Per_Edge() ) {
	      found = 1;
	      // increase temp_tag on edge to increase the number of
	      // faces this edge is connected to. Note that if the 
	      // edge is connected to two faces on the other processor
	      // we only count that as 1 sharing. This may not generate
	      // an error on this processor, but the other processor
	      // will still catch it, thus we don't worry about this
	      // small error.
	      local_edge->temp_tag++;
	      break;
	    }
	  }
	  offset += recv_edge[offset]*2 + 1;
	} // end loop on recv edges from chosen processor
	if ( found == 0) {
          entry.second = -1;
	  send_num_edges[j*2]--;
	}
      } 
    } // end loop on communication partners

    // check to make sure that we have only paired this edge with a total of
    // two edges on this and other processors. If we have done more, then
    // we have an invalid mesh.
    if (local_edge->temp_tag > 2 ) {
      // the count of shared edges is more than 1, thus in parallel
      // this edge represents more than 2 faces sharing an edge.
      // flag this as an invalid mesh.
      ContactNode<Real> *Node0 = local_edge->Node(0);
      ContactNode<Real> *Node1 = local_edge->Node(1);

      std::sprintf(message, 
	      "More than two faces connected to an edge at nodes (%d, %d) & (%d, %d)",
	      Node0->Global_ID().HiInt(),Node0->Global_ID().LoInt(), 
	      Node1->Global_ID().HiInt(),Node1->Global_ID().LoInt());
      //
      //  Store the invalid edge for later retreival by the host code
      //
      ContactShellNode *ShellNode0 = dynamic_cast<ContactShellNode*>(Node0);
      ContactShellNode *ShellNode1 = dynamic_cast<ContactShellNode*>(Node1);

      int global_id0;
      int global_id1;

      if(ShellNode0) {
        global_id0 = ShellNode0->Shell_Node_Base_ID();
      } else {
        global_id0 = Node0->Global_ID().LoInt();
      }
      if(ShellNode1) {
        global_id1 = ShellNode1->Shell_Node_Base_ID();
      } else {
        global_id1 = Node1->Global_ID().LoInt();
      }

      search->append_invalid_edge(global_id0, global_id1);

      errors->Add_Error_Message( message );
      error_code = ContactSearch::INVALID_DATA;
    }
  } // end loop on candidate edges for sharing

  //-----------------------------
  // now complete the specification of the Edge communication lists
  Complete_Edge_Comm_List(num_partners, &shared_edge_list,
			  send_num_edges);

  delete[] recv_handles;
  delete[] send_num_edge_buff;
  delete[] recv_num_edge_buff;
  delete[] send_edge_buff;
  delete[] recv_edge_buff;
  delete[] recv_offsets;

}
#endif

#ifndef CONTACT_NO_MPI
void ContactTopology::Find_Shared_Edge_Candidates( 
      std::vector< std::pair<ContactTopologyEntity<Real>*,int> > *shared_edge_list)
{
  // skip all this if no nodes communicate to other processors.
  if (Node_SymComm->Num_Comm_Partners()!=0) {

    //allocate space for list of processors that might share node
    int * proc_list = new int[Node_SymComm->Num_Comm_Partners()];
    
    // find all edges which are shared with another processor
    ContactEdge<Real>** Edges = 
      reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
    for (int i=0; i<number_of_edges; ++i) {
      ContactEdge<Real>* edge = Edges[i];
      
      PRECONDITION (edge->Number_Face_Connections() <= 2);

      // find all processors that share node 1
      int num_procs = 0;
      for (int j = 0; j < Node_SymComm->Num_Comm_Partners(); ++j ) {
	ContactTopologyEntity<Real> ** nodes = Node_SymComm->Entity_List(j);
	for (int k = 0; k < Node_SymComm->Num_to_Proc(j); ++k ) {
	  if ( nodes[k]->Global_ID() == edge->Node(0)->Global_ID()){
	    proc_list[num_procs] = j;
	    ++num_procs;
	  }
	}
      }
      
      // if no procs share node 1, this is not a shared edge
      if (num_procs == 0 ) continue;
      
      // now loop over sharing processors for node 1 and see if other edge
      // nodes are shared by same processors
      for (int j = 0; j < num_procs; ++j) {
	int sharing_proc = proc_list[j];
	ContactTopologyEntity<Real> ** nodes = Node_SymComm->Entity_List(sharing_proc);
	int num_nodes_sharing_proc = Node_SymComm->Num_to_Proc(sharing_proc);
	int found_edge = 1;
	// loop over remaining nodes on edge
	for (int k = 1; k < edge->Nodes_Per_Edge(); ++k) {
	  int found_node = 0;
	  // loop over nodes in comm list for sharing proc
	  for (int l = 0; l < num_nodes_sharing_proc; ++l ) {
	    // if node is same as one on edge, mark that we found a node
	    if ( nodes[l]->Global_ID() == edge->Node(k)->Global_ID()){
	      found_node = 1;
	      break;
	    }
	  }
	  // if we didn't find last node, then edge is not shared
	  if (!found_node) {
	    found_edge = 0;
	    break;
	  }
	}
	// if edge is shared, add to linked list
	if (found_edge) {
	  shared_edge_list->push_back( std::pair<ContactTopologyEntity<Real>*,int> (edge, 
				       Node_SymComm->Comm_Proc_ID(sharing_proc)));
	}
      } // end of loop on procs which share node 1
    } // end of loop on local edges
    
    delete[] proc_list;
  }
}
#endif

#ifndef CONTACT_NO_MPI
void ContactTopology::Complete_Edge_Comm_List( int num_partners,
      std::vector< std::pair<ContactTopologyEntity<Real>*,int> > *shared_edge_list, int* send_num_edges)
{

  // we have now confirmed all edges that are actually shared, and found those
  // that are not actually shared. Now build the comm list for the edges,
  // based on comparisons on the nodes.

  //   first build everything but the entity list
  int * comm_proc_ids = new int[num_partners];
  int * number_edges_to_partner = new int[num_partners];
  int * send_offsets = new int[num_partners];
  int number_comm_partners = 0;
  int total_num_edges = 0;
  for (int i = 0; i < num_partners; ++i){
    int num_edges_on_proc = send_num_edges[i*2];
    if ( num_edges_on_proc > 0 ) {
      comm_proc_ids[number_comm_partners] = Node_SymComm->Comm_Proc_ID(i);
      number_edges_to_partner[number_comm_partners] = num_edges_on_proc;
      send_offsets[number_comm_partners] = total_num_edges;
      total_num_edges += num_edges_on_proc;
      ++number_comm_partners;
    }
  }

  //  now build the entity list
  ContactTopologyEntity<Real> ** comm_edge_list = new ContactTopologyEntity<Real>*[total_num_edges];
  int * count_edges_to_partner = new int[number_comm_partners];
  for (int i = 0; i < number_comm_partners; ++i ) count_edges_to_partner[i]=0;

  for(int i = 0; i < shared_edge_list->size(); ++i) {
    std::pair<ContactTopologyEntity<Real>*, int> &entry = (*shared_edge_list)[i];

    int dest = entry.second;
    if (dest < 0) continue;
    for (int j = 0; j < number_comm_partners; ++j) {
      if (comm_proc_ids[j] == dest){

	// add edge to list, and make list sorted
	PRECONDITION(count_edges_to_partner[j] <= number_edges_to_partner[j]);
	ContactEdge<Real> * new_edge = (ContactEdge<Real>*) entry.first;
        new_edge->Shared(true);
	int assigned = 0;
	int offset = send_offsets[j];
	for (int k = 0; k < count_edges_to_partner[j]; ++k) {
	  ContactEdge<Real> * stored_edge = 
	    static_cast<ContactEdge<Real> *>(comm_edge_list[offset+k]);
	  // Compare_Edges is > 0 if second edge is less than first
	  if ( Compare_Edges(stored_edge,new_edge) > 0) {
	    for (int l = count_edges_to_partner[j]; l > k; --l)
	      comm_edge_list[offset+l] = comm_edge_list[offset+l-1];
	    comm_edge_list[offset+k] = new_edge;
	    assigned = 1;
	    break;
	  }
	}
	if ( assigned == 0 ) {
	  comm_edge_list[offset+count_edges_to_partner[j]] = new_edge;
	}
	count_edges_to_partner[j]++;
      }
    }
  }
  
  // Now create the edge communication list
  Edge_SymComm = new ContactSymComm( number_comm_partners,
				     comm_proc_ids,
				     number_edges_to_partner,
				     comm_edge_list );
  
  delete[] comm_proc_ids;
  delete[] number_edges_to_partner;
  delete[] count_edges_to_partner;
  delete[] send_offsets;
  delete[] comm_edge_list;

}
#endif

#ifndef CONTACT_NO_MPI
int ContactTopology::Compare_Edges(ContactEdge<Real> * edge1, ContactEdge<Real> * edge2){
  // Compare_Edges compares two edges and returns zero if edge1 is less
  //  than edge2, or 1 if edge1 is greater than edge2.
  PRECONDITION(edge1->Nodes_Per_Edge() <= 3);
  PRECONDITION(edge2->Nodes_Per_Edge() <= 3);
  /*
  PRECONDITION(edge1->Nodes_Per_Edge() == edge2->Nodes_Per_Edge());
  */
  if (edge1->Nodes_Per_Edge() < edge2->Nodes_Per_Edge()) return 0;
  if (edge1->Nodes_Per_Edge() > edge2->Nodes_Per_Edge()) return 1;

  int num_nodes = edge1->Nodes_Per_Edge();
  int edge1_node_owner[3]; 
  int edge1_node_id[3]; 
  int edge2_node_owner[3]; 
  int edge2_node_id[3]; 
  int i,j,k;

  // copy nodes into sorted lists
  for( i = 0; i < num_nodes; ++i){
    int new_owner = edge1->Node(i)->Owner();
    int new_id = edge1->Node(i)->ProcArrayIndex();
    int assigned = 0;
    for (j = 0; j < i; ++j) {
      if ( new_owner < edge1_node_owner[j] ||
	   (new_owner == edge1_node_owner[j] && new_id < edge1_node_id[j])) {
	for (k = i; k > j; --k) {
	  edge1_node_owner[k] = edge1_node_owner[k-1];
	  edge1_node_id[k] = edge1_node_id[k-1];
	}
	edge1_node_owner[j] = new_owner;
	edge1_node_id[j] = new_id;
	assigned = 1;
	break;
      }
    }
    if ( assigned == 0 ) {
      edge1_node_owner[i] = new_owner; 
      edge1_node_id[i] = new_id;
    }
  }

  for( i = 0; i < num_nodes; ++i){
    int new_owner = edge2->Node(i)->Owner();
    int new_id = edge2->Node(i)->ProcArrayIndex();
    int assigned = 0;
    for (j = 0; j < i; ++j) {
      if ( (new_owner < edge2_node_owner[j]) ||
	   (new_owner == edge2_node_owner[j] && new_id < edge2_node_id[j])) {
	for (k = i; k > j; --k) {
	  edge2_node_owner[k] = edge2_node_owner[k-1];
	  edge2_node_id[k] = edge2_node_id[k-1];
	}
	edge2_node_owner[j] = new_owner;
	edge2_node_id[j] = new_id;
	assigned = 1;
	break;
      }
    }
    if ( assigned == 0 ) {
      edge2_node_owner[i] = new_owner; 
      edge2_node_id[i] = new_id;
    }
  }

  // now that both lists are sorted, compare the nodes.
  int result = 0;
  for (i = 0; i< num_nodes; ++i) {
    if ( (edge2_node_owner[i] == edge1_node_owner[i]) &&
	 (edge2_node_id[i] == edge1_node_id[i]) ) continue;
    if ( edge2_node_owner[i] < edge1_node_owner[i] ||
	 (edge2_node_owner[i] == edge1_node_owner[i] && 
	  edge2_node_id[i] < edge1_node_id[i]) ) {
      result = 1;
      break;
    }
    else {
      result = 0;
      break;
    }
  }
  return result;
}
#endif

#ifndef CONTACT_NO_MPI
void ContactTopology::Assign_Secondary_Ownership(ContactZoltan* zoltan,
                                                 VariableHandle POSITION)
{
  // Since it doesn't matter which processor owns the edges and faces
  // in the secondary decomposition (as long as its a processor that would
  // have the object anyway), we will simply assign the owning processor as
  // the owning processor for node 0 in its connectivity (since we are 
  // guaranteed this processor MUST at least ghost this object).

  // There are special cases where we can own a face but not own any of
  // its nodes.  For these cases, we will need an interface here to
  // Zoltans point drop to do the assignment of ownership.
  // Talk to KHB for details.
  
  // Elements will be assigned according to which processor would own
  // the elements centroid.

  int i,j,proc_num;

  Real *position;
  
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    position = node->Variable(POSITION);
    zoltan->Point_Assign(position,&proc_num);
    node->Secondary_Owner(proc_num);
  }

  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
  for (i=0; i<number_of_edges; ++i) {
    ContactEdge<Real>* edge = Edges[i];
    if( edge->Ownership() == ContactTopologyEntity<Real>::OWNED ){
      for( j=0 ; j<edge->Nodes_Per_Edge() ; ++j ){
	if( edge->Node(j)->Ownership() == ContactTopologyEntity<Real>::OWNED ){
	  edge->Secondary_Owner( edge->Node(j)->Secondary_Owner() );
	  break;
	}
      }
      if( j==edge->Nodes_Per_Edge() ){
        position = edge->Node(0)->Variable(POSITION);
        zoltan->Point_Assign( position, &proc_num );
        edge->Secondary_Owner( proc_num );
      }
    }
  }

  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    if( face->Ownership() == ContactTopologyEntity<Real>::OWNED ){
      for( j=0 ; j<face->Nodes_Per_Face() ; ++j ){
	if( face->Node(j)->Ownership() == ContactTopologyEntity<Real>::OWNED ){
	  face->Secondary_Owner( face->Node(j)->Secondary_Owner() );
	  break;
	}
      }
      if( j==face->Nodes_Per_Face() ){
        position = face->Node(0)->Variable(POSITION);
        zoltan->Point_Assign( position, &proc_num );
        face->Secondary_Owner( proc_num );
      }
    }
  }

  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (i=0; i<number_of_elements; ++i) {
    ContactElement* element = Elements[i];
    if( element->Ownership() == ContactTopologyEntity<Real>::OWNED ){
      Real* centroid = element->Variable(ELEMENT_CENTROID);
      zoltan->Point_Assign( centroid, &proc_num );
      element->Secondary_Owner( proc_num );
    }
  }			       
}
#endif

void 
ContactTopology::SortEntityList(int cnt, ContactTopologyEntity<Real>** list)
{
  if (cnt>1) {
    ContactTopologyEntity<Real>* entity;
    int k = cnt>>1;
    int n = cnt;
    for (;;) {
      if (k>0) {
        entity = list[--k];
      } else {
        if (--n <= 0) return;
        entity  = list[n];
        list[n] = list[0];
      }
      int i = k;
      int j = (k<<1)+1;
      while (j<n) {
        if ((j+1<n) && (list[j]->Global_ID() < list[j+1]->Global_ID())) ++j;
        if (entity->Global_ID() < list[j]->Global_ID()) {
          list[i] = list[j];
          i = j;
          j = (i<<1)+1;
        } else break;
      }
      list[i] = entity;
    }
  }
}

void 
ContactTopology::SortEntityList1(int cnt, ContactEdge<Real>** list)
{
  if (cnt>1) {
    ContactEdge<Real>* entity;
    int k = cnt>>1;
    int n = cnt;
    for (;;) {
      if (k>0) {
        entity = list[--k];
      } else {
        if (--n <= 0) return;
        entity  = list[n];
        list[n] = list[0];
      }
      int i = k;
      int j = (k<<1)+1;
      while (j<n) {
        if ((j+1<n) && (EdgeLessThan(list[j],list[j+1]))) ++j;
        if (EdgeLessThan(entity,list[j])) {
          list[i] = list[j];
          i = j;
          j = (i<<1)+1;
        } else break;
      }
      list[i] = entity;
    }
  }
}

void ContactTopology::CleanUp()
{
#ifndef CONTACT_NO_MPI
  if (topology_type==SECONDARY) {
    int i;
    
#ifdef CONTACT_ANALYZE_HASH
    ContactParOStream& postream = search->ParOStream();
    postream.flush();
  
    postream<<"Node block hash tables ===================\n";
    for( i=0 ; i<number_of_node_blocks ; ++i ) {
      postream<<"  Node block "<<i<<"\n";
      postream<<*(node_blocks[i]->NodeList());
    }
  
    postream<<"Face block hash tables ===================\n";
    for( i=0 ; i<number_of_face_blocks ; ++i ) {
      postream<<"  Face block "<<i<<"\n";
      postream<<*(face_blocks[i]->FaceList());
    }
    
    postream<<"Node list hash table ===================\n";
    postream<<*node_list;
 
    postream<<"Face list hash table ===================\n";
    postream<<*face_list;
    
    postream.flush();
#endif

    node_list->CleanUp();
    if( number_of_node_blocks ){
      for( i=0 ; i<number_of_node_blocks ; ++i ) {
        node_blocks[i]->Delete_Node_List();
      }
    }
    
    edge_list->CleanUp();
    if( number_of_edge_blocks ){
      for( i=0 ; i<number_of_edge_blocks ; ++i ) {
        edge_blocks[i]->Delete_Edge_List();
      }
    }
    
    face_list->CleanUp();
    if( number_of_face_blocks ){
      for( i=0 ; i<number_of_face_blocks ; ++i ) {
        face_blocks[i]->Delete_Face_List();
      }
    }
    
    elem_list->CleanUp();
    if( number_of_element_blocks ){
      for( i=0 ; i<number_of_element_blocks ; ++i ) {
        element_blocks[i]->Delete_Element_List();
      }
    }
    
    if( number_of_analytic_surfaces && AnalyticSurfaces){
      for (i=0; i<number_of_analytic_surfaces; ++i) {
        if (AnalyticSurfaces[i] != NULL) delete AnalyticSurfaces[i];
        AnalyticSurfaces[i] = NULL;
      }
      number_of_added_analytic_surfaces = 0;
    }

    if( Node_SymComm ) {
      delete Node_SymComm;
      Node_SymComm = NULL;
    }
    if( Edge_SymComm ) {
      delete Edge_SymComm;
      Edge_SymComm = NULL;
    }

  } else {
  
    if( contact_number_of_processors(SearchComm)>1) DeleteGhosting();
    
  }
#endif  
}    

void ContactTopology::Get_Scratch()
{
  PRECONDITION( scratch_set == false );
  scratch_set = true;
  search->Get_Scratch(number_of_nodes);
}

void ContactTopology::Release_Scratch()
{
  PRECONDITION( scratch_set == true );
  search->Clear_Scratch();
  scratch_set = false;
}

void ContactTopology::Compute_Max_Relative_Node_Motion( Real* max_relative_motion )
{
  int i,j;

  Real local_max_motion[MAX_DIMENSIONALITY];
  Real local_min_motion[MAX_DIMENSIONALITY];
  Real global_max_motion[MAX_DIMENSIONALITY];
  Real global_min_motion[MAX_DIMENSIONALITY];

  for( i=0 ; i<dimensionality ; ++i ) {
    local_max_motion[i] = -BIGNUM;
    local_min_motion[i] =  BIGNUM;
  }

  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  //
  //  Determine the maximum and minimum displacement for each node in the three
  //  coordinate axes.
  //
  for (i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    if (node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
    Real* pos_c = node->Variable(CURRENT_POSITION);
    Real* pos_p = node->Variable(PREDICTED_POSITION);
    for( j=0 ; j<dimensionality ; ++j ) {
      Real delta = pos_p[j] - pos_c[j];
      local_max_motion[j] = std::max(local_max_motion[j],delta);
      local_min_motion[j] = std::min(local_min_motion[j],delta);
    }
  }
  //
  // Determine the global maximum and minimum across all processors
  //
  contact_global_maximum(local_max_motion, global_max_motion,
			 dimensionality, SearchComm);
  contact_global_minimum(local_min_motion, global_min_motion,
			 dimensionality, SearchComm);
  //
  //  Determine the maximum relative motion in each cardinal direction
  //
  for(int idim = 0; idim < dimensionality; ++idim) {
    max_relative_motion[idim] = (global_max_motion[idim] - global_min_motion[idim]);
  }
}

#ifndef CONTACT_NO_MPI
void
ContactTopology::GhostTiedFaces()
{
  if( contact_number_of_processors( SearchComm )>1 ){
    PRECONDITION(have_tied_ghosting==false);
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
      ContactParOStream& postream = search->ParOStream();
      postream<<"    Ghosting off-proccessor master faces for tied interactions\n";
    #endif
    
    int my_proc           = contact_processor_number( SearchComm );
    have_tied_ghosting    = true;
    num_tied_import       = 0;
    int num_off_processor = 0;
    ContactNode<Real>** Nodes   = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      node->temp_tag1   = 0;
      if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
        ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(ContactSearch::STATE0);
        for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE0); ++j) {
          ContactNodeEntityInteraction* cnei = interactions[j];
          if( cnei->Is_Tied() || cnei->Is_InfSlip() ){
            node->SetContextBit(ContactTopologyEntity<Real>::TIED);
            if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
              ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
              if (cnfi->FaceEntityData()->owner != my_proc) {
                ++num_off_processor;
                node->temp_tag1 = 1;
              } else {
                cnfi->Face()->SetContextBit(ContactTopologyEntity<Real>::TIED);
              }
            }
          }//end is tied if
        }//end loop over all interactions on node
      }//end status flag if
    }//end loop over all nodes
    
    ContactZoltanLID zoltanLID;
    ContactZoltanGID zoltanGID;
    LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
    LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
    int        zoltan_pid;
    TiedFaces_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT);
    //=============================================================
    // loop through all the node/face interactions on each node 
    // and add the face to the import list if it is off-processor.
    //=============================================================
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      if (node->temp_tag1==0) continue;
      ContactNodeEntityInteraction** interactions = 
        node->Get_NodeEntity_Interactions(ContactSearch::STATE0);
      for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE0); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if( cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued() ){
          if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
            ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
            if (cnfi->FaceEntityData()->owner != my_proc) {
              // Master face is off-processor so add it to the import list
              ContactHostGlobalID Face_GID( cnfi->FaceEntityData()->host_gid[0], 
                                            cnfi->FaceEntityData()->host_gid[1] );
              zoltan_pid = cnfi->FaceEntityData()->owner;
              PRECONDITION(zoltan_pid>=0);
              cnfi->ZoltanFaceLID(zoltan_lid, 1);
              cnfi->ZoltanFaceGID(zoltan_gid);
              TiedFaces_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
            }
          }
        }
      }
    }
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
      postream<<"      ghosting "<<TiedFaces_ZoltanComm->Num_Import()<<" faces\n";
    #endif
      
    if (contact_global_sum(num_off_processor, SearchComm )>0) {

      //================================================================================
      // we're going to keep the present ghosting so construct the data import commspec
      //================================================================================
      num_tied_import  = TiedFaces_ZoltanComm->Num_Import();
      tied_import_gids = new LB_ID_TYPE [num_tied_import*ZOLTAN_GID_SIZE];
      tied_import_lids = new LB_ID_TYPE [num_tied_import*ZOLTAN_LID_SIZE];
      tied_import_pids = new int [num_tied_import];
       
      num_tied_export  = -1;
      tied_export_gids = NULL;
      tied_export_lids = NULL;
      tied_export_pids = NULL;
       
      std::memcpy(tied_import_gids, 
                  TiedFaces_ZoltanComm->Import_GIDS(), 
                  num_tied_import*ZOLTAN_GID_SIZE*sizeof(int));
             
      std::memcpy(tied_import_lids, 
                  TiedFaces_ZoltanComm->Import_LIDS(), 
                  num_tied_import*ZOLTAN_LID_SIZE*sizeof(int));
             
      std::memcpy(tied_import_pids, 
                  TiedFaces_ZoltanComm->Import_Procs(), 
                  num_tied_import*sizeof(int));
  
      search->Get_Zoltan()->Set_UpdateTiedImportCallBacks(search);
      TiedCommSpec = new ContactZoltanCommUtils(ContactZoltanCommUtils::IMPORT,
                                                search->Get_Zoltan()->Get_ZoltanPtr(),
                                                num_tied_import,tied_import_gids,
                                                tied_import_lids, tied_import_pids,
                                                search->ParOStream());
    }
  } else {
    ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(ContactSearch::STATE0);
      for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE0); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if( cnei->Is_Tied() || cnei->Is_InfSlip() ){
          node->SetContextBit(ContactTopologyEntity<Real>::TIED);
          if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
            ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
            cnfi->Face()->SetContextBit(ContactTopologyEntity<Real>::TIED);
          }
        }//end is tied if
      }//end loop over all interactions on node
    }//end loop over all nodes
  }
}

void
ContactTopology::UpdateTiedFaces()
{
  if( contact_number_of_processors( SearchComm )>1 && TiedCommSpec!=NULL){
    ContactParOStream& postream = search->ParOStream();
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"    Updating Ghosted off-proccessor faces\n";
    #endif

    search->Get_Zoltan()->Set_UpdateTiedImportCallBacks(search);
    TiedCommSpec->Migrate(postream);
  }
}

void
ContactTopology::DeleteTiedFaces()
{
  if( contact_number_of_processors( SearchComm )>1 && have_tied_ghosting){
    have_tied_ghosting = false;
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
    ContactParOStream& postream = search->ParOStream();
    postream<<"    Releasing ghosted tied master faces\n";
    #endif
   
    if (TiedCommSpec!=NULL) {
      delete TiedCommSpec;
      TiedCommSpec = NULL;
    }
    if (TiedFaces_ZoltanComm!=NULL) {
      delete TiedFaces_ZoltanComm;
      TiedFaces_ZoltanComm = NULL;
    }
    if (TiedNodes_ZoltanComm!=NULL) {
      delete TiedNodes_ZoltanComm;
      TiedNodes_ZoltanComm = NULL;
    }
    if (tied_import_gids) delete [] tied_import_gids;
    if (tied_import_lids) delete [] tied_import_lids;
    if (tied_import_pids) delete [] tied_import_pids;
    tied_import_gids = NULL;
    tied_import_lids = NULL;
    tied_import_pids = NULL;
    tied_export_gids = NULL;
    tied_export_lids = NULL;
    tied_export_pids = NULL;
  }
}


ContactSearch::ContactErrorCode
ContactTopology::GhostNFIfaces(int flag)
{
  ContactSearch::ContactErrorCode err_code(ContactSearch::NO_ERROR);
  if( contact_number_of_processors( SearchComm )>1 ){
    PRECONDITION(have_local_ghosting==false);
    have_local_ghosting = true;
    num_ghosted_import  = 0;
    #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_ANALYZE_DATA_XFER)
      ContactParOStream& postream = search->ParOStream();
      #ifdef CONTACT_ANALYZE_DATA_XFER
        search->Bytes_For_Nodes(0.0);
        search->Bytes_For_Faces(0.0);
      #endif
      #if CONTACT_DEBUG_PRINT_LEVEL>=2
        postream<<"    Ghosting off-proccessor faces\n";
      #endif
    #endif
    
    bool added_neighbors = false;
    int my_proc = contact_processor_number( SearchComm );
    ContactZoltanLID zoltanLID;
    ContactZoltanGID zoltanGID;
    LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
    LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
    int        zoltan_pid;
    if (GhostFaces_ZoltanComm!=NULL) {
      GhostFaces_ZoltanComm->Initialize( ContactZoltanComm::ZOLTAN_IMPORT);
    } else {
      GhostFaces_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT);
    }
    //=============================================================
    // loop through all the node/face interactions on each node 
    // and add the face to the import list if it is off-processor.
    //
    // node->temp_tag1 = 1 indicates that a node has an interaction
    // where the master face is off-processor that will have to be
    // connected to it's interaction at the end of this function.
    //=============================================================
    ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      node->temp_tag1 = 0;
      //if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
      ContactNodeEntityInteraction** interactions = 
        node->Get_NodeEntity_Interactions(ContactSearch::STATE1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE1); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if (cnei->Is_Tied() || cnei->Is_InfSlip() || 
            cnei->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
        ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
        if (cnfi->FaceEntityData()->owner != my_proc) {
          // Master face is off-processor so add it to the import list
          ContactHostGlobalID Face_GID( cnfi->FaceEntityData()->host_gid[0], 
                                        cnfi->FaceEntityData()->host_gid[1] );
          node->temp_tag1 = 1;
          zoltan_pid = cnfi->FaceEntityData()->owner;
          PRECONDITION(zoltan_pid>=0);
          cnfi->ZoltanFaceLID(zoltan_lid, 1);
          cnfi->ZoltanFaceGID(zoltan_gid);
          GhostFaces_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
        }
        if (cnfi->Is_Tracked()) {
          for (int k=0; k<cnfi->NumSharedFaces(); ++k) {
            ContactTopologyEntity<Real>::connection_data *face_info = cnfi->SharedFaceData(k);
            if (face_info->owner != my_proc) {
              // Neighbor face is off-processor so add it to the import list
              PRECONDITION(face_info->owner>=0);
              zoltan_pid = face_info->owner;
              ContactHostGlobalID GID( face_info->host_gid[0], face_info->host_gid[1] );
              zoltanLID.ZoltanLID(CT_FACE, face_info->owner_proc_array_index, zoltan_lid);
              zoltanGID.ZoltanGID(CT_FACE, &GID, zoltan_gid);
              GhostFaces_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
              added_neighbors = true;
            }
          }
        }
      }
    }
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
      postream<<"      ghosting "<<GhostFaces_ZoltanComm->Num_Import()<<" faces\n";
    #endif

    num_ghosted_import           = GhostFaces_ZoltanComm->Num_Import();
    int global_num_ghosted_faces = contact_global_sum(num_ghosted_import, SearchComm );
      
    if (global_num_ghosted_faces>0) {
  
      #ifdef CONTACT_ANALYZE_DATA_XFER
        search->Bytes_For_Nodes(0.0);
        search->Bytes_For_Faces(0.0);
      #endif

      //=============================================================
      // migrate all the node/face interaction's off-processor faces
      //=============================================================
      for (int i=0; i<number_of_face_blocks; ++i) {
        ghosted_face_blocks[i]->FaceList()->SetupHash(GhostFaces_ZoltanComm->Num_Import());
      }
      int	num_xtra_import   = GhostFaces_ZoltanComm->Num_Import();
      LB_ID_PTR import_xtra_gids  = GhostFaces_ZoltanComm->Import_GIDS();
      LB_ID_PTR import_xtra_lids  = GhostFaces_ZoltanComm->Import_LIDS();
      int*	import_xtra_procs = GhostFaces_ZoltanComm->Import_Procs();
      int	num_xtra_export   = -1;
      LB_ID_PTR export_xtra_gids  = NULL;
      LB_ID_PTR export_xtra_lids  = NULL;
      int*	export_xtra_procs = NULL;
      if (added_neighbors) {
        search->Get_Zoltan()->Set_GhostingImportCallBacks(search);
      } else {
        search->Get_Zoltan()->Set_GhostingExportCallBacks(search);
      }
      search->Get_Zoltan()->Help_Migrate(num_xtra_import,  import_xtra_gids, 
  		                         import_xtra_lids, import_xtra_procs, 
  		                         num_xtra_export,  export_xtra_gids,
  		                         export_xtra_lids, export_xtra_procs);

      if (GhostOthers_ZoltanComm!=NULL) {
        GhostOthers_ZoltanComm->Initialize( ContactZoltanComm::ZOLTAN_IMPORT);
      } else {
        GhostOthers_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT);
      }
                                                                 
      //========================================================
      // loop thru all the ghosted faces and ghost needed nodes 
      //========================================================
      int num_node_import = 0;
      for (int i=0; i<number_of_face_blocks; ++i) {
        ContactBlockEntityList* block_face_list = ghosted_face_blocks[i]->FaceList();
        block_face_list->IteratorStart();
        while (ContactTopologyEntity<Real>* entity=block_face_list->IteratorForward()) {
          ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
          ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
          POSTCONDITION( node_info );
          for(int j=0 ; j<face->Nodes_Per_Face() ; ++j){
            ContactNode<Real>* node = static_cast<ContactNode<Real> *>(node_list->Find(&node_info[j]));
            if (!node) {
              zoltan_pid = node_info[j].owner;
	      ContactHostGlobalID GID( node_info[j].host_gid[0], 
                                       node_info[j].host_gid[1] );
              zoltanLID.ZoltanLID(CT_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
              zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
              GhostOthers_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
              ++num_node_import;
            }
          }
        }
      }
      for (int i=0; i<number_of_node_blocks; ++i) {
        ghosted_node_blocks[i]->NodeList()->SetupHash(num_node_import);
      }
      num_xtra_import   = GhostOthers_ZoltanComm->Num_Import();
      import_xtra_gids  = GhostOthers_ZoltanComm->Import_GIDS();
      import_xtra_lids  = GhostOthers_ZoltanComm->Import_LIDS();
      import_xtra_procs = GhostOthers_ZoltanComm->Import_Procs();
      num_xtra_export   = -1;
      export_xtra_gids  = NULL;
      export_xtra_lids  = NULL;
      export_xtra_procs = NULL;
      search->Get_Zoltan()->Set_GhostingExportCallBacks(search);
      search->Get_Zoltan()->Help_Migrate(num_xtra_import,  import_xtra_gids, 
  		                         import_xtra_lids, import_xtra_procs, 
  		                         num_xtra_export,  export_xtra_gids,
  		                         export_xtra_lids, export_xtra_procs);

      #ifdef CONTACT_ANALYZE_DATA_XFER
        postream<<"    Data Sent by Zoltan when ghosting tracked faces...\n";
        postream<<"      nodes:  "<<search->Bytes_For_Nodes()/1024.0<<" Kb\n";
        postream<<"      faces:  "<<search->Bytes_For_Faces()/1024.0<<" Kb\n";
        postream.flush();
      #endif

      if (flag>0) {   
        //============================================================
        // If we're going to keep the present ghosting for subsequent 
        // tracking steps, then construct the data import commspec
        //============================================================
        int num_faces       = GhostFaces_ZoltanComm->Num_Import();
        int num_nodes       = GhostOthers_ZoltanComm->Num_Import();
        num_ghosted_import  = num_faces+num_nodes;
        ghosted_import_gids = new LB_ID_TYPE [num_ghosted_import*ZOLTAN_GID_SIZE];
        ghosted_import_lids = new LB_ID_TYPE [num_ghosted_import*ZOLTAN_LID_SIZE];
        ghosted_import_pids = new int [num_ghosted_import];
         
  	num_ghosted_export  = -1;
        ghosted_export_gids = NULL;
  	ghosted_export_lids = NULL;
        ghosted_export_pids = NULL;
         
        std::memcpy(ghosted_import_gids, 
                    GhostFaces_ZoltanComm->Import_GIDS(), 
                    num_faces*ZOLTAN_GID_SIZE*sizeof(int));
        std::memcpy(&ghosted_import_gids[num_faces*ZOLTAN_GID_SIZE], 
                    GhostOthers_ZoltanComm->Import_GIDS(), 
                    num_nodes*ZOLTAN_GID_SIZE*sizeof(int));
               
        std::memcpy(ghosted_import_lids, 
                    GhostFaces_ZoltanComm->Import_LIDS(), 
                    num_faces*ZOLTAN_LID_SIZE*sizeof(int));
        std::memcpy(&ghosted_import_lids[num_faces*ZOLTAN_LID_SIZE], 
                    GhostOthers_ZoltanComm->Import_LIDS(), 
                    num_nodes*ZOLTAN_LID_SIZE*sizeof(int));
               
        std::memcpy(ghosted_import_pids, 
                    GhostFaces_ZoltanComm->Import_Procs(), 
                    num_faces*sizeof(int));
        std::memcpy(&ghosted_import_pids[num_faces], 
                    GhostOthers_ZoltanComm->Import_Procs(), 
                    num_nodes*sizeof(int));
  
        search->Get_Zoltan()->Set_UpdateGhostingImportCallBacks(search);
        GhostingCommSpec = new ContactZoltanCommUtils(ContactZoltanCommUtils::IMPORT,
                                                      search->Get_Zoltan()->Get_ZoltanPtr(),
                                                      num_ghosted_import,ghosted_import_gids,
                                                      ghosted_import_lids, ghosted_import_pids,
                                                      search->ParOStream());
      }

    }

    if (num_ghosted_import>0) {
      //====================================================================
      // loop thru all the ghosted faces and connect the nodes to the faces
      //====================================================================
      for (int i=0; i<number_of_face_blocks; ++i) {
        ContactBlockEntityList* entity_list = ghosted_face_blocks[i]->FaceList();
        entity_list->IteratorStart();
        while (ContactTopologyEntity<Real>* entity=entity_list->IteratorForward()) {
          ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
          //connect the nodes to the face
          ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
          POSTCONDITION( node_info );
          for(int j=0 ; j<face->Nodes_Per_Face() ; ++j){
            ContactNode<Real>* node = NULL;
            for (int k=0; k<number_of_node_blocks; ++k) {
              node = static_cast<ContactNode<Real> *>(ghosted_node_blocks[k]->NodeList()->Find(&node_info[j]));
              if (node) break;
            }
            if (!node) {
              node = static_cast<ContactNode<Real> *>(node_list->Find(&node_info[j]));
            }
            POSTCONDITION(node);
            face->ConnectNode( j, node );
          }
        }
      }
      
      //==============================================================
      // connect all the ghosted faces to the appropriate interaction
      //
      // node->temp_tag1 = 1 indicates that a node has an interaction
      // where the master face is off-processor that will have to be
      // connected to it's interaction.
      //==============================================================
      for (int i=0; i<number_of_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        //if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
        if (node->temp_tag1) {
          bool found_invalid = false;
          ContactNodeEntityInteraction** interactions = 
            node->Get_NodeEntity_Interactions(ContactSearch::STATE1);
          for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE1); ++j) {
            ContactNodeEntityInteraction* cnei = interactions[j];
            if (cnei->Is_Tied() || cnei->Is_InfSlip() || 
                cnei->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
            ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
            if (cnfi->FaceEntityData()->owner != my_proc) {
              ContactHostGlobalID GID( cnfi->FaceEntityData()->host_gid[0], 
                                       cnfi->FaceEntityData()->host_gid[1] );
              int block = cnfi->FaceEntityData()->block_id;
              ContactFace<Real>* face = static_cast<ContactFace<Real> *>
                (ghosted_face_blocks[block]->FaceList()->Find( cnfi->FaceEntityData() ));
              if (face!=NULL) {
                cnfi->Connect_Face( face );
              } else {
                // MWG: We really need to return a fatal error code here and abort!
                //      Since we are explicitly ghosting the master faces for each 
                //      interaction, this condition should never occurr.
                Real* position = node->Variable(CURRENT_POSITION);
                Real* rem_gap  = node->Variable(REMAINING_GAP);
		static int printlimit = 0;
		if (printlimit < 5) {
		  std::cout<<"P"<<my_proc<<": In ContactTopology::GhostNFIfaces(), step "<<search->StepNumber()<<"\n"
			   <<"  error with interaction "<<j
			   <<", cannot connect master face "<<"("<<cnfi->FaceEntityData()->host_gid[1]<<", "<<cnfi->FaceEntityData()->host_gid[0]<<")\n"
			   <<"    node "<<node->Global_ID()<<" exodus_id:"<< node->Exodus_ID()<<"\n"
			   <<"    position      = ("<<position[0]<<", "<<position[1]<<", "<<position[2]<<")\n"
			   <<"    remaining_gap = ("<<rem_gap[0]<<", "<<rem_gap[1]<<", "<<rem_gap[2]<<")\n"<<std::flush;
		  printlimit ++;
		}
		node->Delete_NodeEntity_Interaction(cnfi,1);
                found_invalid = true;
		
		// What we should do is die here. Unfortunately, we 
		// risk branch stability with this, so for now 
		// we just report the error and move on.
		// err_code = ContactSearch::INTERNAL_ERROR;
              }
            }
          }
        }
      }
    }
  }
  return err_code;
}

void
ContactTopology::UpdateGhostedNFIfaces()
{
  if( contact_number_of_processors( SearchComm )>1 && GhostingCommSpec!=NULL){
    ContactParOStream& postream = search->ParOStream();
#if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(ANALYZE_DATA_XFER)
    //ContactParOStream& postream = search->ParOStream();
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"    Updating Ghosted off-proccessor faces\n";
#endif
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(ANALYZE_DATA_XFER)
    search->Bytes_For_Nodes(0.0);
    search->Bytes_For_Faces(0.0);
#endif

    search->Get_Zoltan()->Set_UpdateGhostingImportCallBacks(search);
    GhostingCommSpec->Migrate(postream);

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(ANALYZE_DATA_XFER)
    postream<<"    Data Sent by Zoltan when updating ghosted faces for tracking...\n";
    postream<<"      nodes:  "<<search->Bytes_For_Nodes()/1024.0<<" Kb\n";
    postream<<"      faces:  "<<search->Bytes_For_Faces()/1024.0<<" Kb\n";
    postream.flush();
#endif
    
    //================================================================
    // connect all the ghosted faces to the appropriate interaction
    // NOTE: need to do this because the enforcement does it's own
    // ghosting and connects those faces to the interactions, 
    // overwriting any that are left over from the initial ghosting
    // that's done for tracking.  We'll need to work on the enforcment
    // to make it more aware of the persistant ghosting left over 
    // from track or 'no secondary'
    //================================================================
    int my_proc = contact_processor_number( SearchComm );
    ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(ContactSearch::STATE1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE1); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if (cnei->Is_Tied() || cnei->Is_InfSlip() || 
            cnei->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
        ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
        if (cnfi->FaceEntityData()->owner != my_proc) {
          ContactHostGlobalID GID( cnfi->FaceEntityData()->host_gid[0], 
                                   cnfi->FaceEntityData()->host_gid[1] );
          int block = cnfi->FaceEntityData()->block_id;
          ContactFace<Real>* face = static_cast<ContactFace<Real> *>
            (ghosted_face_blocks[block]->FaceList()->Find( cnfi->FaceEntityData() ));
          if (face!=NULL) {
            cnfi->Connect_Face( face );
          } else {
            // MWG: Since we are explicitly ghosting the master faces for each interaction,
            //      this condition should never occur!
            POSTCONDITION(false);
          }
        }
      }
    }
  }
}

void
ContactTopology::DeleteGhostedNFIfaces()
{
  if( contact_number_of_processors( SearchComm )>1 && have_local_ghosting){
    have_local_ghosting = false;
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    ContactParOStream& postream = search->ParOStream();
    postream<<"    Releasing ghosted faces\n";
#endif
    if( num_ghosted_import>0 ){
      if( number_of_node_blocks ){
        for (int i=0 ; i<number_of_node_blocks ; ++i) {
          ghosted_node_blocks[i]->Delete_Node_List();
        }
      }
      if( number_of_face_blocks ){
        for (int i=0 ; i<number_of_face_blocks ; ++i) {
          ghosted_face_blocks[i]->Delete_Face_List();
        }
      }
      // unconnect all the ghosted faces to the appropriate interaction
      int my_proc = contact_processor_number( SearchComm );
      ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
      for (int i=0; i<number_of_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        //if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
        if (node->temp_tag1) {
          ContactNodeEntityInteraction** interactions = 
            node->Get_NodeEntity_Interactions(ContactSearch::STATE0);
          for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE0); ++j) {
            ContactNodeEntityInteraction* cnei = interactions[j];
            if (cnei->Is_Tied() || cnei->Is_InfSlip() || 
                cnei->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
            ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
            if (cnfi->FaceEntityData()->owner != my_proc) {
              cnfi->Connect_Face( (ContactFace<Real>*)NULL );
            }
          }
          interactions = node->Get_NodeEntity_Interactions(ContactSearch::STATE1);
          for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE1); ++j) {
            ContactNodeEntityInteraction* cnei = interactions[j];
            if (cnei->Is_Tied() || cnei->Is_InfSlip() || 
                cnei->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
            ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
            if (cnfi->FaceEntityData()->owner != my_proc) {
              cnfi->Connect_Face( (ContactFace<Real>*)NULL );
            }
          }
        }
      }
    }
    if (GhostingCommSpec!=NULL) {
      delete GhostingCommSpec;
      GhostingCommSpec = NULL;
    }
    GhostFaces_ZoltanComm->CleanUp();
    GhostOthers_ZoltanComm->CleanUp();
    if (ghosted_import_gids) delete [] ghosted_import_gids;
    if (ghosted_import_lids) delete [] ghosted_import_lids;
    if (ghosted_import_pids) delete [] ghosted_import_pids;
    ghosted_import_gids = NULL;
    ghosted_import_lids = NULL;
    ghosted_import_pids = NULL;
    ghosted_export_gids = NULL;
    ghosted_export_lids = NULL;
    ghosted_export_pids = NULL;
  }
}

void
ContactTopology::DoGhosting(VariableHandle POSITION, const Real& reasonable_gap)
{
  if( contact_number_of_processors( SearchComm ) == 1 ) return;
  if (have_global_ghosting==true) DeleteGhosting(); 
  have_global_ghosting = true;
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    ContactParOStream& postream = search->ParOStream();
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
      postream<<"    Performing ghosting for the primary\n";
    #endif
  #endif
  #ifdef CONTACT_TIMINGS
    search->Timer()->Start_Timer( search->secondary_owner_migration_time );
  #endif
  int my_proc_id = contact_processor_number(SearchComm);
  LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
  int        zoltan_pid;
    
  ContactNode<Real>**    nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  ContactFace<Real>**    faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  ContactElement** elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
        
  if (GhostFaces_ZoltanComm!=NULL) {
    GhostFaces_ZoltanComm->Initialize( ContactZoltanComm::ZOLTAN_EXPORT);
  } else {   
    GhostFaces_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_EXPORT);
  }
                                           
  #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
    search->Bytes_For_Nodes(0.0);
    search->Bytes_For_Faces(0.0);
    search->Bytes_For_Elems(0.0);
  #endif

  for (int i=0; i<number_of_node_blocks; ++i) {
    ghosted_node_blocks[i]->NodeList()->UseHash(1000);
  }
  for (int i=0; i<number_of_face_blocks; ++i) {
    ghosted_face_blocks[i]->FaceList()->UseHash(1000);
  }
  for (int i=0; i<number_of_element_blocks; ++i) {
    ghosted_element_blocks[i]->ElemList()->UseHash(1000);
  }
      
  //=========================================================================
  //  T O P _ L E V E L   F A C E S
  //=========================================================================
    
  Real user_search_tol = search->Search_Data()->Max_Search_Tolerance();
  Real box_inflation   = search->BoxInflation();
  int num_config       = (POSITION==CURRENT_POSITION)?1:2;

  //
  //  Determine the actual bounding box of each processors set of elements.  
  //  This may be significantly smaller than box_tol + RCB size.  Bounding 
  //  boxes will be assemebed into a global list which will be reduced and 
  //  copied to all processors.  This will allow culling out ghosts that are 
  //  not strictly nessecary.
  //  
  //  NKC note:  This calculation is not huge, but is non trivial.  Definetly 
  //  worth it for SPH problems or problems where normal smoothing is used and 
  //  the object sizes vary significantly.  There might also may be some work 
  //  reduction for other problem types.  Need to determine if it is significant
  //  enough to warrant peforming the calculation.  Also may need to do similar 
  //  calcs for other portions of tolerance.  For example, huge velocity on one 
  //  node may yeild universally huge ghosting rads, big tolerance on one
  //  interaction could be a problem.  Large adhesion or other enforcement 
  //  tolerance on a small subset of nodes yields very large ghosting radius.
  //
  for (int i=0; i<number_of_nodes; ++i) {
    nodes[i]->Secondary_Owner(nodes[i]->Owner());
  }
  for (int i=0; i<number_of_faces; ++i) {
    faces[i]->Secondary_Owner(faces[i]->Owner());
  }
  for (int i=0; i<number_of_elements; ++i) {
    elems[i]->Secondary_Owner(elems[i]->Owner());
  }
  int total_num_procs = contact_number_of_processors( SearchComm );  
  ObjectBoundingBox *proc_box_array = NULL;
  proc_box_array = new ObjectBoundingBox[total_num_procs];
  search->Create_Processor_Bounding_Boxes(proc_box_array, POSITION, total_num_procs);
  int hierarchy_size = 2 * total_num_procs - 1;
  ObjectBoundingBoxHierarchy *proc_box_hierarchy = NULL;
  if(hierarchy_size > 0) {
    proc_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
    ObjectBoundingBoxHierarchy::create_hierarchy(proc_box_hierarchy, proc_box_array, total_num_procs);
  }
  int* proc_list = new int [total_num_procs];

  //=======================================================================
  // Ghost all faces to determine where they should be sent.  This section 
  // builds the export part of the asymmetric face communication object.
  //=======================================================================

  for (int ii=0; ii<number_of_faces; ++ii) {
    ContactFace<Real>* face = faces[ii];
    // get bounding box
    ContactBoundingBox face_current_box;
    ContactBoundingBox face_predicted_box;
    ContactBoundingBox object_box;
    face->ComputeBoundingBoxForSearch(num_config,
				      CURRENT_POSITION,
				      POSITION, 
                                      auto_tol,
                                      box_inflation,
                                      user_search_tol,
	                              face_current_box,
	                              face_predicted_box,
                                      object_box);
    //
    //  NKC note, really want to search each processor bounding box
    //  vs. the current object box.  The processor bounding boxes can overlap.
    //  However, the RCB boxes defined by zoltan will not overlap.
    //  By expanding each object box by the global object size tolerance it ensures
    //  the object box will overlap any potential zoltan RCB bounding box.
    //  May need to not use zoltan and do these calcs by hand later to 
    //  optimize the search on actual processor box size.  Using the global object tolerance
    //  here has the potential to find many processor overlaps that need to by thrown out later
    //  by the proc_box_array overlap calculation
    // 
    //
    int numprocs = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,object_box, proc_list, numprocs);
    // loop over the proc_list from the box drop -- if the local 
    // processor is in the list, mark the current face as needing to
    // be copied locally from primary to secondary decomposition (set
    // temp_tag to 1)
    for (int i=0; i<numprocs; ++i) {
      //
      //  Double check that the current bounding box actually overlaps the 
      //  true processor bounding box.  If not, do not ghost this object
      //
      int proc_num = proc_list[i];
      if (proc_num != my_proc_id) {
	// off processor - add to communication object
	zoltan_pid = proc_list[i];
	face->ZoltanLID(CT_FACE, zoltan_lid);
	face->ZoltanGID(CT_FACE, zoltan_gid);
	GhostFaces_ZoltanComm->Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
      }
    }
  }
  
  //=========================================================================
  //  T O P _ L E V E L   E L E M E N T S
  //=========================================================================

  //=======================================================================
  // Ghost all elements to determine where they should be sent.  This section 
  // builds the export part of the asymmetric element communication object.
  //=======================================================================
  for (int ii=0; ii<number_of_elements; ++ii) {
    ContactElement* element = elems[ii];
    // get bounding box
    ContactBoundingBox object_box;
    element->ComputeBoundingBoxForSearch(num_config, 
                                         CURRENT_POSITION,
                                         POSITION, 
                                         auto_tol,
                                         box_inflation,
                                         user_search_tol,
                                         object_box);
    int numprocs = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,object_box, proc_list, numprocs);
    // loop over the proc_list from the box drop -- if the local 
    // processor is in the list, mark the current element as needing to
    // be copied locally from primary to secondary decomposition (set
    // temp_tag to 1)
    for (int i=0; i<numprocs; ++i) {
      //
      //  Double check that the current bounding box actually overlaps the true processor bounding box.  
      //  If not, do not ghost this object
      //
      int proc_num = proc_list[i];
      if (proc_num != my_proc_id) {
	// off processor - add to communication object
	zoltan_pid = proc_list[i];
	element->ZoltanLID(CT_ELEMENT, zoltan_lid);
	element->ZoltanGID(CT_ELEMENT, zoltan_gid);
	GhostFaces_ZoltanComm->Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
      }
    }
  }

  //=========================================================================
  //  T O P _ L E V E L   N O D E S   ( B L O C K S   1 - N )
  //=========================================================================
  // In this block, we are only moving "POINTS" not "NODES", 
  // we know they are all ContactNode<Real>'s not ContactShellNodes 
  // so I won't even check
  int base = 0;
  if (number_of_face_blocks+number_of_element_blocks > 0) base=1;
  for(int i=base; i<number_of_node_blocks; ++i) {
    ContactNodeBlock* node_block = node_blocks[i];
    if( node_block->Type() == ContactSearch::POINT ){
      int nnodes = node_list->BlockNumEntities(i);
      ContactNode<Real>** block_nodes  = 
	reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
          
      if (node_block->Has_Radius_Attributes()) {
	for (int j=0; j<nnodes; ++j) {
	  ContactNode<Real>* node = block_nodes[j];
	  if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
	  Real* position =  node->Variable(POSITION);
	  Real  radius   = *node->Variable(NODE_RADIUS);

	  ContactBoundingBox object_box;
	  object_box.add_sphere(position, radius);
	  int numprocs = 0;
	  ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,object_box, proc_list, numprocs);
	  // loop over the proc_list from the box drop -- if the local 
	  // processor is in the list, mark the current node as needing to
	  // be copied locally from primary to secondary decomposition (set
	  // temp_tag to 1)
	  for (int k=0; k<numprocs; ++k) {
	    //
	    //  Double check that the current bounding box actually overlaps the true processor bounding box.  If not, do not
	    //  ghost this object
	    //
	    if (proc_list[k] != my_proc_id) {
	      // off processor - add to communication object
	      zoltan_pid = proc_list[k];
	      node->ZoltanLID(CT_NODE, zoltan_lid);
	      node->ZoltanGID(CT_NODE, zoltan_gid);
	      GhostFaces_ZoltanComm->Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
	    }   
	  }
	}
      } else {
	for (int j=0; j<nnodes; ++j) {
	  ContactNode<Real>* node = block_nodes[j];
	  if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
	  if (node->Secondary_Owner()!=my_proc_id) {
	    // off processor - add to communication object
	    node->ZoltanLID(CT_NODE, zoltan_lid);
	    node->ZoltanGID(CT_NODE, zoltan_gid);
	    zoltan_pid = node->Secondary_Owner();
	    GhostFaces_ZoltanComm->Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
	  }
	}
      }
    }
  }
  if(proc_box_hierarchy) delete [] proc_box_hierarchy;
  if(proc_box_array) delete [] proc_box_array;
  delete [] proc_list;
  //
  //  Have Zoltan transfer the specified export objects
  //
  MigrateExportedData();
}

void
ContactTopology::DeleteGhosting()
{
  if( contact_number_of_processors( SearchComm )>1 && have_global_ghosting){
    have_global_ghosting = false;
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    ContactParOStream postream( SearchComm );
    postream<<"    Deleting ghosting for the secondary\n";
#endif
    int my_proc = contact_processor_number( SearchComm );
    if( number_of_node_blocks ){
      for(int i=0 ; i<number_of_node_blocks ; ++i ) {
        ghosted_node_blocks[i]->Delete_Node_List();
      }
    }
    if( number_of_face_blocks ){
      for(int i=0 ; i<number_of_face_blocks ; ++i ) {
        ghosted_face_blocks[i]->Delete_Face_List();
      }
    }
    if( number_of_element_blocks ){
      for(int i=0 ; i<number_of_element_blocks ; ++i ) {
        ghosted_element_blocks[i]->Delete_Element_List();
      }
    }

    node_list->BuildList(Node_Blocks(), Number_of_Node_Blocks(),
                         no_parallel_consistency==ContactSearch::INACTIVE);
    face_list->BuildList(Face_Blocks(), Number_of_Face_Blocks(),
                         no_parallel_consistency==ContactSearch::INACTIVE);
    elem_list->BuildList(Element_Blocks(), Number_of_Element_Blocks(),
                         no_parallel_consistency==ContactSearch::INACTIVE);
    number_of_nodes    = node_list->NumEntities();
    number_of_faces    = face_list->NumEntities();
    number_of_elements = elem_list->NumEntities();
    ContactNode<Real>**    nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
    ContactFace<Real>**    faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
    ContactElement** elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
    
#if CONTACT_DEBUG_PRINT_LEVEL>=3
    postream<<"After search:\n";
    postream<<"  number_of_nodes    = "<<number_of_nodes<<"\n";
    postream<<"  number_of_faces    = "<<number_of_faces<<"\n";
    postream<<"  number_of_elements = "<<number_of_elements<<"\n";
#endif

    //===============================
    // Reset final enitity ownership 
    //===============================
    int my_proc_id  = contact_processor_number(SearchComm);
    
    for (int i=0; i<number_of_nodes; ++i) {
      if (nodes[i]->Owner()==my_proc_id)  {
        nodes[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
      } else {
        nodes[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
        //PRECONDITION(nodes[i]->Number_Interactions(1)==0);
      }
    }
    
    for (int i=0; i<number_of_faces; ++i) {
      if (faces[i]->Owner()==my_proc_id)  {
        faces[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
      } else {
        faces[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
      }
    }
    
    for (int i=0; i<number_of_elements; ++i) {
      if (elems[i]->Owner()==my_proc_id)  {
        elems[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
      } else {
        elems[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
      }
    }
    
    if (number_of_faces>0) {      
      //====================================================================
      // loop thru all the faces and reconnect the faces to the nodes
      //====================================================================
      for (int i=0; i<number_of_nodes; ++i) {
        nodes[i]->Delete_Face_Connections();
      }
      for (int i=0; i<number_of_faces; ++i) {
        ContactFace<Real>* face = faces[i];
        for(int k=0 ; k<face->Nodes_Per_Face() ; ++k ){
          ContactNode<Real>* node = face->Node(k);
          node->Connect_Face( face );
        }
      }
      for (int i=0; i<number_of_nodes; ++i) {
        nodes[i]->SortConnectedFaces();
      }
    }
    
    if (search->Do_NodeNode_Search()) {
      // unconnect all the ghosted nodes from the appropriate interaction
      ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
      for (int i=0; i<number_of_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        ContactInteractionDLL<Real>* interactions = node->Get_NodeNode_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()) {
            ContactNodeNodeInteraction* cnni = 
              static_cast<ContactNodeNodeInteraction*>(interaction);
            if (cnni->MasterNodeEntityData()->owner != my_proc) {
              cnni->Connect_MasterNode( (ContactNode<Real>*)NULL );
            }
          }
        }
      }
    }
    
    if (search->Do_NodeFace_Search()) {
      // unconnect all the ghosted faces from the appropriate interaction
      ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
      for (int i=0; i<number_of_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
        ContactNodeEntityInteraction** interactions = 
          node->Get_NodeEntity_Interactions(0);
        for (int j=0; j<node->Number_NodeEntity_Interactions(0); ++j) {
          if (interactions[j]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
          ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[j]);
          if (cnfi->FaceEntityData()->owner != my_proc) {
            cnfi->Connect_Face( (ContactFace<Real>*)NULL );
          }
        }
        interactions = node->Get_NodeEntity_Interactions(1);
        for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
          if (interactions[j]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
          ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[j]);
          if (cnfi->FaceEntityData()->owner != my_proc) {
            cnfi->Connect_Face( (ContactFace<Real>*)NULL );
          }
        }
      }
    }
    
    if (search->Do_FaceFace_Search()) {
      // unconnect all the ghosted faces from the appropriate interaction
      ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
      for (int i=0; i<number_of_faces; ++i) {
        ContactFace<Real>* face = Faces[i];
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()) {
          ContactFaceFaceInteraction<Real>* cffi = 
            static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
          if (cffi->MasterFaceEntityData()->owner != my_proc) {
            cffi->Connect_MasterFace( (ContactFace<Real>*)NULL );
          }
        }
      }
    }
    
    if (search->Do_ElemElem_Search()) {
      // unconnect all the ghosted elements from the appropriate interaction
      ContactElement** Elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
      for (int i=0; i<number_of_elements; ++i) {
        ContactElement* element = Elems[i];
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()) {
          ContactElementElementInteraction* ceei = 
            static_cast<ContactElementElementInteraction*>(interaction);
          if (ceei->MasterElementEntityData()->owner != my_proc) {
            ceei->Connect_MasterElement( (ContactElement*)NULL );
          }
        }
      }
    }
    if (GhostingCommSpec!=NULL) {
      delete GhostingCommSpec;
      GhostingCommSpec = NULL;
    }
    GhostFaces_ZoltanComm->CleanUp();
    GhostOthers_ZoltanComm->CleanUp();
    if (ghosted_import_gids) delete [] ghosted_import_gids;
    if (ghosted_import_lids) delete [] ghosted_import_lids;
    if (ghosted_import_pids) delete [] ghosted_import_pids;
    ghosted_import_gids = NULL;
    ghosted_import_lids = NULL;
    ghosted_import_pids = NULL;
    ghosted_export_gids = NULL;
    ghosted_export_lids = NULL;
    ghosted_export_pids = NULL;
  }
}

void
ContactTopology::DoGhosting_New_NodeFace(VariableHandle POSITION, const Real &reasonable_gap) {
  //
  //  If running in serial, no cross processor ghosting is needed, exit this routine early.
  //
  if( contact_number_of_processors( SearchComm ) == 1 ) return;
  if(have_global_ghosting==true) DeleteGhosting();
  have_global_ghosting = true;
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    ContactParOStream& postream = search->ParOStream();
    postream<<"    Performing ghosting for the primary\n";
  #endif
  #ifdef CONTACT_TIMINGS
    search->Timer()->Start_Timer( search->create_secondary_time );
    search->Timer()->Start_Timer( search->secondary_owner_migration_time );
  #endif
  int my_proc_id      = contact_processor_number(SearchComm);
  int total_num_procs = contact_number_of_processors(SearchComm);
  LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
    
  ContactNode<Real>**    nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  ContactFace<Real>**    faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  ContactElement** elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
        
  if (GhostFaces_ZoltanComm!=NULL) {
    GhostFaces_ZoltanComm->Initialize(ContactZoltanComm::ZOLTAN_EXPORT);
  } else {
    GhostFaces_ZoltanComm = new ContactZoltanComm(ContactZoltanComm::ZOLTAN_EXPORT);
  }
                                           
  #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
    search->Bytes_For_Nodes(0.0);
    search->Bytes_For_Faces(0.0);
    search->Bytes_For_Elems(0.0);
  #endif
  for (int i=0; i<number_of_node_blocks; ++i) {
    ghosted_node_blocks[i]->NodeList()->UseHash(1000);
  }
  for (int i=0; i<number_of_face_blocks; ++i) {
    ghosted_face_blocks[i]->FaceList()->UseHash(1000);
  }
  for (int i=0; i<number_of_element_blocks; ++i) {
    ghosted_element_blocks[i]->ElemList()->UseHash(1000);
  }
  //
  //  Create owning topology.  As no secondary decomposition has been created, the secondary owners of 
  //  all entities are the same as the primary owners
  //
  for (int i=0; i<number_of_nodes; ++i) {
    nodes[i]->Secondary_Owner(nodes[i]->Owner());
  }
  for (int i=0; i<number_of_faces; ++i) {
    faces[i]->Secondary_Owner(faces[i]->Owner());
  }
  for (int i=0; i<number_of_elements; ++i) {
    elems[i]->Secondary_Owner(elems[i]->Owner());
  }

  //=========================================================================
  //  T O P _ L E V E L   F A C E S
  //=========================================================================
    
  //===================================================
  // Compute tolerance for the drop box to do ghosting 
  //===================================================
  Real user_search_tol = search->Search_Data()->Max_Search_Tolerance();
  Real box_inflation   = search->BoxInflation();
  int  num_configs     = (POSITION==CURRENT_POSITION)?1:2;
  //    
  //
  //  Do ghosting, follow many steps
  //
  //  Step 1:  Create the data structures to hold the per processor bounding boxes 
  //           of all active nodes and active faces
  //
  ObjectBoundingBox *proc_total_boxes = new ObjectBoundingBox[total_num_procs*2];
  ObjectBoundingBox *proc_total_node_boxes = proc_total_boxes;
  ObjectBoundingBox *proc_total_face_boxes = proc_total_boxes + total_num_procs;
  //
  //  Step 2:  Assemble lists of objects to ghost on each processor.  Each object
  //           must store the object bounding box
  //
  vector<ContactBoundingBox> local_nodes;
  vector<ContactBoundingBox> local_faces;  
  ObjectBoundingBox my_proc_total_node_box;
  ObjectBoundingBox my_proc_total_face_box;
  //
  //  Create bounding boxes for all nodes on the current processor
  //
  local_nodes.reserve(number_of_nodes);
  for(int inode = 0; inode < number_of_nodes; ++inode) {
    ContactNode<Real> *node = nodes[inode];
    if (node->Secondary_Owner() != my_proc_id ||                    
        !node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE) ||
        node->Physical_Type()   == ContactNode<Real>::SHELL_TAB_NODE ||
        node->Physical_Type()   == ContactNode<Real>::MIXED_TAB_NODE) continue;
    ContactBoundingBox node_BB;
    node->ComputeBoundingBoxForSearch(num_configs,
                                      NODE_GHOST_GAP, 
                                      CURRENT_POSITION,
                                      POSITION, 
                                      auto_tol,
                                      box_inflation,
                                      node_BB);
    local_nodes.push_back(node_BB);
    my_proc_total_node_box.add_box(node_BB);
  }
  //
  //  Loop over faces by entity block to ensure that the face list is sorted by entity_key
  //
  int global_face_index = 0;
  local_faces.reserve(number_of_faces);
  const int num_face_blocks = Number_of_Face_Blocks();
  for(int iface_block = 0; iface_block < num_face_blocks; ++iface_block) {
    ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(iface_block));    
    int num_faces = face_list->BlockNumEntities(iface_block);
    //
    //  Add each face to the list
    //
    
    for(int iface = 0; iface < num_faces; ++iface, ++global_face_index) {    
      ContactFace<Real> *face = Faces[iface];
      PRECONDITION(face->Nodes_Per_Face() > 0);
      ContactBoundingBox face_current_box;
      ContactBoundingBox face_predicted_box;
      ContactBoundingBox face_bounding_box;
      face->ComputeBoundingBoxForSearch(num_configs,
				        CURRENT_POSITION,
				        POSITION, 
                                        auto_tol,
                                        box_inflation,
                                        user_search_tol,
	                                face_current_box,
	                                face_predicted_box,
                                        face_bounding_box);
      local_faces.push_back(face_bounding_box);
      my_proc_total_face_box.add_box(face_bounding_box); 
    }
  }
  proc_total_node_boxes[my_proc_id] = my_proc_total_node_box;
  proc_total_face_boxes[my_proc_id] = my_proc_total_face_box;
  //
  //  Step 3:  Perform global reductions so that all processors know the total
  //           bounding boxes of all other processors.
  //
  ObjectBoundingBox::global_box_combine(proc_total_boxes, total_num_procs*2, SearchComm);
  //
  //  'Remove' the bounding boxes for the current processor from the processor lists.  Don't
  //  want to ghost objects back to the processor that currently owns them
  //
  proc_total_node_boxes[my_proc_id].Reset();
  proc_total_face_boxes[my_proc_id].Reset();

  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    proc_total_node_boxes[iproc].set_object_number(iproc);
    proc_total_face_boxes[iproc].set_object_number(iproc);
  }
  //
  //  Step 4:  Create search trees from the processor bounding boxes.  The search tree
  //           will be used to determine if a given box potentially interacts 
  //           with some object on another processor.
  //
  int hierarchy_size = 2 * total_num_procs - 1;
  ObjectBoundingBoxHierarchy * proc_node_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
  ObjectBoundingBoxHierarchy * proc_face_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
  ObjectBoundingBoxHierarchy::create_hierarchy(proc_node_box_hierarchy, proc_total_node_boxes, total_num_procs);
  ObjectBoundingBoxHierarchy::create_hierarchy(proc_face_box_hierarchy, proc_total_face_boxes, total_num_procs);
  //
  //  Step 5:  Determine the potential communication parnters.  If the node box for one processor 
  //           overlaps the face box of another, than that processor is a potential communication partner.  
  //           Note, since the proc_total_boxes arrays are identical on all processors, the comm_partner 
  //           array will be valid and complete on all processors without requiring any additional 
  //           communication.
  //
  vector<int> comm_partners_send(total_num_procs, 0);
  vector<int> comm_partners_recv(total_num_procs, 0);
  int *search_list = new int[total_num_procs];
  int list_size = 0;
  ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_node_box_hierarchy, 
                                                      my_proc_total_face_box, 
                                                      search_list,
                                                      list_size);
  for(int ilist = 0; ilist < list_size; ++ilist) {
    comm_partners_recv[search_list[ilist]] = 1;
  }
  list_size = 0;
  ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_face_box_hierarchy, 
                                                      my_proc_total_node_box, 
                                                      search_list,
                                                      list_size);
  for(int ilist = 0; ilist < list_size; ++ilist) {
    comm_partners_send[search_list[ilist]] = 1;
  }
  //
  //  Step 6:  Loop over all objects to compute the potential objects to ghost.  If a face 
  //           overlaps another processor node box, add it to the face ghost list.  If a node 
  //           overlaps another processors face box, add it to the node ghost list.
  //
  //           Maintain a seperate list that lists the root face index for each sent face, this will
  //           be used later when recording the actual objects to ghost
  //
  vector< vector<ObjectBoundingBox> > node_boxes_to_send(total_num_procs);
  int num_local_nodes = local_nodes.size();
  for(int inode = 0; inode < num_local_nodes; ++inode) {
    list_size = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_face_box_hierarchy,
                                                        local_nodes[inode],
                                                        search_list,
                                                        list_size);
    for(int ilist = 0; ilist < list_size; ++ilist) {
      int send_to_proc_num = search_list[ilist];
      if(send_to_proc_num == my_proc_id) continue;
      node_boxes_to_send[send_to_proc_num].push_back(ObjectBoundingBox(local_nodes[inode], inode));
    } 
  }
  vector< vector<ObjectBoundingBox> > face_boxes_to_send(total_num_procs);
  int num_local_faces = local_faces.size();
  for(int iface = 0; iface < num_local_faces; ++iface) {
    list_size = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_node_box_hierarchy,
                                                        local_faces[iface],
                                                        search_list,
                                                        list_size);
    for(int ilist = 0; ilist < list_size; ++ilist) {
      int send_to_proc_num = search_list[ilist];
      if(send_to_proc_num == my_proc_id) continue;
      face_boxes_to_send[send_to_proc_num].push_back(ObjectBoundingBox(local_faces[iface], iface));
    } 
  }
  delete [] search_list;
  //
  //  Step 7: Transfer the node boxes to the receiving processor
  //
  vector< vector<ObjectBoundingBox> > node_boxes_to_receive(total_num_procs);
  ACME::Parallel_Data_Exchange(node_boxes_to_send, node_boxes_to_receive, comm_partners_send, comm_partners_recv, SearchComm);
  //
  //  Step 8: Create the ghosting lists, these lists describe which nodes and 
  //          which faces must be ghosted from the various processor groupings.  For each 
  //          valid processor grouping, search the receive node lists versus the sending face 
  //          lists.
  //
  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    if(comm_partners_recv[iproc] == 0) continue;
    vector<ObjectBoundingBox> &received_node_boxes  = node_boxes_to_receive[iproc];
    vector<ObjectBoundingBox> &send_face_boxes     = face_boxes_to_send[iproc];
    int num_received_nodes = received_node_boxes.size();
    int num_sent_faces    = send_face_boxes.size();
    if(num_received_nodes > 0 && num_sent_faces > 0) {
      //
      //  Create a bounding box hierarchy out of the received nodes
      //
      ObjectBoundingBoxHierarchy *received_node_hierarchy = new ObjectBoundingBoxHierarchy[2 * num_received_nodes - 1];
      ObjectBoundingBoxHierarchy::create_hierarchy(received_node_hierarchy,
						   &(received_node_boxes[0]),
						   num_received_nodes);
      //
      //  Search faces against nodes.  If a match is found, must ghost the sent face, and the received node
      //
      for(int iface = 0; iface < num_sent_faces; ++iface) {
	if(ObjectBoundingBoxHierarchy::find_any_overlap_loop(received_node_hierarchy, 
							     send_face_boxes[iface])) {
	  const int face_index = send_face_boxes[iface].get_object_number();
	  ContactFace<Real> *face = faces[face_index];
	  face->ZoltanLID(CT_FACE, zoltan_lid);
	  face->ZoltanGID(CT_FACE, zoltan_gid);
	  GhostFaces_ZoltanComm->Add_Export(zoltan_lid,zoltan_gid,iproc);
	}
      }
      delete [] received_node_hierarchy;
    }
  }
  //
  //  Step 11:  Delete all alloated data.
  //
  delete [] proc_node_box_hierarchy;
  delete [] proc_face_box_hierarchy;
  delete [] proc_total_boxes;
  //
  //  Step 12:  Have Zoltan perform the actual data migration
  //    
  MigrateExportedData();
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->create_secondary_time );
#endif
}

void
ContactTopology::DoCaptureGhosting_New_NodeFace(VariableHandle POSITION, const Real &reasonable_gap) {
  //
  //  If running in serial, no cross processor ghosting is needed, exit this routine early.
  //
  if( contact_number_of_processors( SearchComm ) == 1 ) return;
  if (have_global_ghosting==true) DeleteGhosting();
  have_global_ghosting = true;
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    ContactParOStream& postream = search->ParOStream();
    postream<<"    Performing ghosting for the primary\n";
  #endif
  #ifdef CONTACT_TIMINGS
    search->Timer()->Start_Timer( search->create_secondary_time );
    search->Timer()->Start_Timer( search->secondary_owner_migration_time );
  #endif
  int my_proc_id      = contact_processor_number(SearchComm);
  int total_num_procs = contact_number_of_processors(SearchComm);
  LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
    
  ContactNode<Real>**    nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  ContactFace<Real>**    faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  ContactElement** elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
        
  if (GhostFaces_ZoltanComm!=NULL) {
    GhostFaces_ZoltanComm->Initialize(ContactZoltanComm::ZOLTAN_EXPORT);
  } else {
    GhostFaces_ZoltanComm = new ContactZoltanComm(ContactZoltanComm::ZOLTAN_EXPORT);
  }
                                           
  #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
    search->Bytes_For_Nodes(0.0);
    search->Bytes_For_Faces(0.0);
    search->Bytes_For_Elems(0.0);
  #endif
  for (int i=0; i<number_of_node_blocks; ++i) {
    ghosted_node_blocks[i]->NodeList()->UseHash(1000);
  }
  for (int i=0; i<number_of_face_blocks; ++i) {
    ghosted_face_blocks[i]->FaceList()->UseHash(1000);
  }
  for (int i=0; i<number_of_element_blocks; ++i) {
    ghosted_element_blocks[i]->ElemList()->UseHash(1000);
  }
  //
  //  Create owning topology.  As no secondary decomposition has been created, the secondary owners of 
  //  all entities are the same as the primary owners
  //
  for (int i=0; i<number_of_nodes; ++i) {
    nodes[i]->Secondary_Owner(nodes[i]->Owner());
  }
  for (int i=0; i<number_of_faces; ++i) {
    faces[i]->Secondary_Owner(faces[i]->Owner());
  }
  for (int i=0; i<number_of_elements; ++i) {
    elems[i]->Secondary_Owner(elems[i]->Owner());
  }

  //=========================================================================
  //  T O P _ L E V E L   F A C E S
  //=========================================================================
    
  int  num_configs = (POSITION!=CURRENT_POSITION)?2:1;
  
  //===================================================
  // Compute tolerance for the drop box to do ghosting 
  //===================================================
  Real capture_tol     = search->CaptureMotion();
  Real user_search_tol = search->Search_Data()->Max_Search_Tolerance();
  Real box_inflation   = search->BoxInflation();
  Real box_expand      = user_search_tol+capture_tol;
  //    
  //
  //  Do ghosting, follow many steps
  //
  //  Step 1:  Create the data structures to hold the per processor bounding boxes 
  //           of all active nodes and active faces
  //
  ObjectBoundingBox *proc_total_boxes = new ObjectBoundingBox[total_num_procs*2];
  ObjectBoundingBox *proc_total_node_boxes = proc_total_boxes;
  ObjectBoundingBox *proc_total_face_boxes = proc_total_boxes + total_num_procs;
  //
  //  Step 2:  Assemble lists of objects to ghost on each processor.  Each object
  //           must store the object bounding box
  //
  vector<ContactBoundingBox> local_nodes;
  vector<ContactBoundingBox> local_faces;  
  ObjectBoundingBox my_proc_total_node_box;
  ObjectBoundingBox my_proc_total_face_box;
  //
  //  Create bounding boxes for all nodes on the current processor
  //
  local_nodes.reserve(number_of_nodes);
  for(int inode = 0; inode < number_of_nodes; ++inode) {
    ContactNode<Real> *node = nodes[inode];
    if (node->Secondary_Owner() != my_proc_id ||                    
        !node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE) ||
        node->Physical_Type()   == ContactNode<Real>::SHELL_TAB_NODE ||
        node->Physical_Type()   == ContactNode<Real>::MIXED_TAB_NODE) continue;
    ContactBoundingBox node_BB;
    node->ComputeBoundingBoxForSearch(num_configs,
                                      NODE_GHOST_GAP, 
                                      CURRENT_POSITION,
                                      POSITION, 
                                      auto_tol,
                                      box_inflation,
                                      capture_tol,
                                      node_BB);
    local_nodes.push_back(node_BB);
    my_proc_total_node_box.add_box(node_BB);
  }
  //
  //  Loop over faces by entity block to ensure that the face list is sorted by entity_key
  //
  int global_face_index = 0;
  local_faces.reserve(number_of_faces);
  const int num_face_blocks = Number_of_Face_Blocks();
  for(int iface_block = 0; iface_block < num_face_blocks; ++iface_block) {
    ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(iface_block));    
    int num_faces = face_list->BlockNumEntities(iface_block);
    //
    //  Add each face to the list
    //
    
    for(int iface = 0; iface < num_faces; ++iface, ++global_face_index) {    
      ContactFace<Real> *face = Faces[iface];
      PRECONDITION(face->Nodes_Per_Face() > 0);
      ContactBoundingBox face_current_box;
      ContactBoundingBox face_predicted_box;
      ContactBoundingBox face_bounding_box;
      face->ComputeBoundingBoxForSearch(num_configs,
				        CURRENT_POSITION,
				        POSITION, 
                                        auto_tol,
                                        box_inflation,
   					box_expand,
	                                face_current_box,
	                                face_predicted_box,
                                        face_bounding_box);
      local_faces.push_back(face_bounding_box);
      my_proc_total_face_box.add_box(face_bounding_box); 
    }
  }
  proc_total_node_boxes[my_proc_id] = my_proc_total_node_box;
  proc_total_face_boxes[my_proc_id] = my_proc_total_face_box;
  //
  //  Step 3:  Perform global reductions so that all processors know the total
  //           bounding boxes of all other processors.
  //
  ObjectBoundingBox::global_box_combine(proc_total_boxes, total_num_procs*2, SearchComm);
  //
  //  'Remove' the bounding boxes for the current processor from the processor lists.  Don't
  //  want to ghost objects back to the processor that currently owns them
  //
  proc_total_node_boxes[my_proc_id].Reset();
  proc_total_face_boxes[my_proc_id].Reset();

  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    proc_total_node_boxes[iproc].set_object_number(iproc);
    proc_total_face_boxes[iproc].set_object_number(iproc);
  }
  //
  //  Step 4:  Create search trees from the processor bounding boxes.  The search tree
  //           will be used to determine if a given box potentially interacts 
  //           with some object on another processor.
  //
  int hierarchy_size = 2 * total_num_procs - 1;
  ObjectBoundingBoxHierarchy * proc_node_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
  ObjectBoundingBoxHierarchy * proc_face_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
  ObjectBoundingBoxHierarchy::create_hierarchy(proc_node_box_hierarchy, proc_total_node_boxes, total_num_procs);
  ObjectBoundingBoxHierarchy::create_hierarchy(proc_face_box_hierarchy, proc_total_face_boxes, total_num_procs);
  //
  //  Step 5:  Determine the potential communication parnters.  If the node box for one processor 
  //           overlaps the face box of another, than that processor is a potential communication partner.  
  //           Note, since the proc_total_boxes arrays are identical on all processors, the comm_partner 
  //           array will be valid and complete on all processors without requiring any additional 
  //           communication.
  //
  vector<int> comm_partners_send(total_num_procs, 0);
  vector<int> comm_partners_recv(total_num_procs, 0);
  int *search_list = new int[total_num_procs];
  int list_size = 0;
  ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_node_box_hierarchy, 
                                                      my_proc_total_face_box, 
                                                      search_list,
                                                      list_size);
  for(int ilist = 0; ilist < list_size; ++ilist) {
    comm_partners_recv[search_list[ilist]] = 1;
  }
  list_size = 0;
  ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_face_box_hierarchy, 
                                                      my_proc_total_node_box, 
                                                      search_list,
                                                      list_size);
  for(int ilist = 0; ilist < list_size; ++ilist) {
    comm_partners_send[search_list[ilist]] = 1;
  }
  //
  //  Step 6:  Loop over all objects to compute the potential objects to ghost.  If a face 
  //           overlaps another processor node box, add it to the face ghost list.  If a node 
  //           overlaps another processors face box, add it to the node ghost list.
  //
  //           Maintain a seperate list that lists the root face index for each sent face, this will
  //           be used later when recording the actual objects to ghost
  //
  vector< vector<ObjectBoundingBox> > node_boxes_to_send(total_num_procs);
  int num_local_nodes = local_nodes.size();
  for(int inode = 0; inode < num_local_nodes; ++inode) {
    list_size = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_face_box_hierarchy,
                                                        local_nodes[inode],
                                                        search_list,
                                                        list_size);
    for(int ilist = 0; ilist < list_size; ++ilist) {
      int send_to_proc_num = search_list[ilist];
      if(send_to_proc_num == my_proc_id) continue;
      node_boxes_to_send[send_to_proc_num].push_back(ObjectBoundingBox(local_nodes[inode], inode));
    } 
  }
  vector< vector<ObjectBoundingBox> > face_boxes_to_send(total_num_procs);
  int num_local_faces = local_faces.size();
  for(int iface = 0; iface < num_local_faces; ++iface) {
    list_size = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_node_box_hierarchy,
                                                        local_faces[iface],
                                                        search_list,
                                                        list_size);
    for(int ilist = 0; ilist < list_size; ++ilist) {
      int send_to_proc_num = search_list[ilist];
      if(send_to_proc_num == my_proc_id) continue;
      face_boxes_to_send[send_to_proc_num].push_back(ObjectBoundingBox(local_faces[iface], iface));
    } 
  }
  delete [] search_list;
  //
  //  Step 7: Transfer the node boxes to the receiving processor
  //
  vector< vector<ObjectBoundingBox> > node_boxes_to_receive(total_num_procs);
  ACME::Parallel_Data_Exchange(node_boxes_to_send, node_boxes_to_receive, comm_partners_send, comm_partners_recv, SearchComm);
  //
  //  Step 8: Create the ghosting lists, these lists describe which nodes and 
  //          which faces must be ghosted from the various processor groupings.  For each 
  //          valid processor grouping, search the receive node lists versus the sending face 
  //          lists.
  //
  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    if(comm_partners_recv[iproc] == 0) continue;
    vector<ObjectBoundingBox> &received_node_boxes  = node_boxes_to_receive[iproc];
    vector<ObjectBoundingBox> &send_face_boxes     = face_boxes_to_send[iproc];
    int num_received_nodes = received_node_boxes.size();
    int num_sent_faces    = send_face_boxes.size();
    if(num_received_nodes > 0 && num_sent_faces > 0) {
      //
      //  Create a bounding box hierarchy out of the received nodes
      //
      ObjectBoundingBoxHierarchy *received_node_hierarchy = new ObjectBoundingBoxHierarchy[2 * num_received_nodes - 1];
      ObjectBoundingBoxHierarchy::create_hierarchy(received_node_hierarchy,
						   &(received_node_boxes[0]),
						   num_received_nodes);
      //
      //  Search faces against nodes.  If a match is found, must ghost the sent face, and the received node
      //
      for(int iface = 0; iface < num_sent_faces; ++iface) {
	if(ObjectBoundingBoxHierarchy::find_any_overlap_loop(received_node_hierarchy, 
							     send_face_boxes[iface])) {
	  const int face_index = send_face_boxes[iface].get_object_number();
	  ContactFace<Real> *face = faces[face_index];
	  face->ZoltanLID(CT_FACE, zoltan_lid);
	  face->ZoltanGID(CT_FACE, zoltan_gid);
	  GhostFaces_ZoltanComm->Add_Export(zoltan_lid,zoltan_gid,iproc);
	}
      }
      delete [] received_node_hierarchy;
    }
  }
  //
  //  Step 11:  Delete all alloated data.
  //
  delete [] proc_node_box_hierarchy;
  delete [] proc_face_box_hierarchy;
  delete [] proc_total_boxes;
  //
  //  Step 12:  Have Zoltan perform the actual data migration
  //    
  MigrateExportedData();
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->create_secondary_time );
#endif
}

void ContactTopology::UpdateGhostingSetupForNoSecondary()
{
  ContactParOStream& postream = search->ParOStream();
  
  // initialize the flag indicating entities are to be updated from off-processor
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    Nodes[i]->temp_tag1 = 0;
  }
  ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    Faces[i]->temp_tag1 = 0;
  }
  
  // tag all faces that are connected to nodes in proximity
  // (so the physical face calculation has all the info)
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    if (node->in_proximity) {
      node->temp_tag1 = 1;
      int num_faces = node->Number_Face_Connections();
      for (int j=0; j<num_faces; ++j) {
        node->GetFace(j)->temp_tag1 = 1;
      }
    }
  }
  
  // tag all faces that are in proximity
  for (int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    if (face->in_proximity) face->temp_tag1 = 1;
  }
  
  // tag all nodes connected to tagged faces
  for (int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    if (face->temp_tag1>0) {
      int num_nodes = face->Nodes_Per_Face();
      for (int j=0; j<num_nodes; ++j) {
        ContactNode<Real>* node = face->Node(j);
        node->temp_tag1 = 1;
      }
    }
  }
  
  // Loop thru all the ghosted nodes and count  
  // the number that need to be updated
  int num_nodes = 0;
  for (int i=0; i<number_of_node_blocks; ++i) {
    ContactBlockEntityList* ghost_node_list = ghosted_node_blocks[i]->NodeList();
    ghost_node_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=ghost_node_list->IteratorForward()) {
      if (entity->temp_tag1>0) ++num_nodes;
    }
  }
  
  // Loop thru all the ghosted faces and count  
  // the number that need to be updated
  int num_faces = 0;
  for (int i=0; i<number_of_face_blocks; ++i) {
    ContactBlockEntityList* ghost_face_list = ghosted_face_blocks[i]->FaceList();
    ghost_face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=ghost_face_list->IteratorForward()) {
      if (entity->temp_tag1>0) ++num_faces;
    }
  }
      
  num_ghosted_import  = num_faces+num_nodes;
  if (num_ghosted_import>0) {
    ghosted_import_gids = new LB_ID_TYPE [num_ghosted_import*ZOLTAN_GID_SIZE];
    ghosted_import_lids = new LB_ID_TYPE [num_ghosted_import*ZOLTAN_LID_SIZE];
    ghosted_import_pids = new int [num_ghosted_import];
  } else {
    ghosted_import_gids = NULL;
    ghosted_import_lids = NULL;
    ghosted_import_pids = NULL;
  }
   
  num_ghosted_export  = -1;
  ghosted_export_gids = NULL;
  ghosted_export_lids = NULL;
  ghosted_export_pids = NULL;
  
  int num_entity = 0;
  
  // Add the nodes that need to be updated to the import list
  for (int i=0; i<number_of_node_blocks; ++i) {
    ContactBlockEntityList* ghost_node_list = ghosted_node_blocks[i]->NodeList();
    ghost_node_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=ghost_node_list->IteratorForward()) {
      if (entity->temp_tag1>0) {
        entity->ZoltanGID(CT_NODE,&ghosted_import_gids[num_entity*ZOLTAN_GID_SIZE]);
        ghosted_import_pids[num_entity++] = entity->Owner();
      }
    }
  }
  
  // Add the faces that need to be updated to the import list
  for (int i=0; i<number_of_face_blocks; ++i) {
    ContactBlockEntityList* ghost_face_list = ghosted_face_blocks[i]->FaceList();
    ghost_face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=ghost_face_list->IteratorForward()) {
      if (entity->temp_tag1>0) {
        entity->ZoltanGID(CT_FACE,&ghosted_import_gids[num_entity*ZOLTAN_GID_SIZE]);
        ghosted_import_pids[num_entity++] = entity->Owner();
      }
    }
  }
  
  POSTCONDITION(num_entity==num_ghosted_import);
  
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream<<"    Setup ghosting update for "<<num_ghosted_import<<" entities ("<<num_nodes<<" nodes & "<<num_faces<<" faces)\n";
#endif

  // construct the data import commspec for subquent tracking steps
  search->Get_Zoltan()->Set_UpdateGhostingImportCallBacks(search);
  //ContactParOStream& postream = search->ParOStream();
  GhostingCommSpec = new ContactZoltanCommUtils(ContactZoltanCommUtils::IMPORT,
                                                search->Get_Zoltan()->Get_ZoltanPtr(),
                                                num_entity,
                                                ghosted_import_gids,
                                                ghosted_import_lids,
                                                ghosted_import_pids,
                                                postream);
}

void ContactTopology::UpdateGhosting()
{
  if( contact_number_of_processors( SearchComm )>1 && GhostingCommSpec!=NULL){
    ContactParOStream& postream = search->ParOStream();
#if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(ANALYZE_DATA_XFER)
    //ContactParOStream& postream = search->ParOStream();
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"    Updating Ghosted off-proccessor faces\n";
#endif
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(ANALYZE_DATA_XFER)
    search->Bytes_For_Nodes(0.0);
    search->Bytes_For_Faces(0.0);
#endif

    search->Get_Zoltan()->Set_UpdateGhostingImportCallBacks(search);
    GhostingCommSpec->Migrate(postream);

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(ANALYZE_DATA_XFER)
    postream<<"    Data Sent by Zoltan when updating ghosted faces for tracking...\n";
    postream<<"      nodes:  "<<search->Bytes_For_Nodes()/1024.0<<" Kb\n";
    postream<<"      faces:  "<<search->Bytes_For_Faces()/1024.0<<" Kb\n";
    postream.flush();
#endif
    
    //================================================================
    // connect all the ghosted faces to the appropriate interaction
    // NOTE: need to do this because the enforcement does it's own
    // ghosting and connects those faces to the interactions, 
    // overwriting any that are left over from the initial ghosting
    // that's done for tracking.  We'll need to work on the enforcment
    // to make it more aware of the persistant ghosting left over 
    // from track or 'no secondary'
    //================================================================
    int my_proc = contact_processor_number( SearchComm );
    ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
        if (interactions[j]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
        ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[j]);
        if (cnfi->FaceEntityData()->owner != my_proc) {
          ContactHostGlobalID GID( cnfi->FaceEntityData()->host_gid[0], 
                                   cnfi->FaceEntityData()->host_gid[1] );
          int block = cnfi->FaceEntityData()->block_id;
          ContactFace<Real>* face = static_cast<ContactFace<Real> *>
            (ghosted_face_blocks[block]->FaceList()->Find( cnfi->FaceEntityData() ));
          if (face!=NULL) {
            cnfi->Connect_Face( face );
          } else {
            POSTCONDITION(false);
          }
        }
      }
    }
  }
}

void ContactTopology::MigrateExportedData() 
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  ContactParOStream& postream = search->ParOStream();
#endif
  ContactNode<Real>**    nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  ContactFace<Real>**    faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  ContactElement** elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
  int my_proc_id      = contact_processor_number(SearchComm);
  LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
  int        zoltan_pid;

  int       num_toplevel_import   = -1;
  LB_ID_PTR import_toplevel_gids  = NULL;
  LB_ID_PTR import_toplevel_lids  = NULL;
  int*      import_toplevel_procs = NULL;
  int       num_toplevel_export   = GhostFaces_ZoltanComm->Num_Export();
  LB_ID_PTR export_toplevel_gids  = GhostFaces_ZoltanComm->Export_GIDS();
  LB_ID_PTR export_toplevel_lids  = GhostFaces_ZoltanComm->Export_LIDS();
  int*      export_toplevel_procs = GhostFaces_ZoltanComm->Export_Procs();

  //================================================
  // Migrate all the off processor toplevel objects
  //================================================
#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->owner_help_migrate_time );
#endif
  search->Get_Zoltan()->Set_GhostingExportCallBacks(search);
  search->Get_Zoltan()->Help_Migrate(num_toplevel_import,  import_toplevel_gids, 
				     import_toplevel_lids, import_toplevel_procs, 
				     num_toplevel_export,  export_toplevel_gids,
				     export_toplevel_lids, export_toplevel_procs);
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->owner_help_migrate_time );
#endif

  ContactZoltanLID zoltanLID;
  ContactZoltanGID zoltanGID;


  /*
  //
  //  Alternative to Zoltan:  Manually create the export object lists easier for ACME to do this than Zoltan and
  //  creates a more accurate list of objects to ghost.
  //
  vector< vector<int> > import_local_ids(total_num_procs);
  for(int i=0 ; i<number_of_face_blocks ; ++i ){
    ContactBlockEntityList* block_face_list = ghosted_face_blocks[i]->FaceList();
    block_face_list->IteratorStart();
    while(ContactTopologyEntity<Real>* entity = block_face_list->IteratorForward() ){
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
      POSTCONDITION( node_info );
      for(int j=0 ; j<face->Nodes_Per_Face() ; ++j){
	ContactNode<Real>* node = static_cast<ContactNode<Real> *>(NodeList()->Find(&node_info[j]));
	if (!node) {
	  //
	  //  Add the node to the import list
	  //
	  PRECONDITION(node_info[j].owner !=my_proc_id);
          import_local_ids[node_info[j].owner].push_back(node_info[j].owner_proc_array_index);
	}
      }
    }
  }

  for(int i=0 ; i<number_of_element_blocks ; ++i ){
    ContactBlockEntityList* block_element_list = element_blocks[i]->ElemList();
    block_element_list->IteratorStart();
    while(ContactTopologyEntity<Real>* entity=block_element_list->IteratorForward() ){
      ContactElement* element = static_cast<ContactElement*>(entity);
      ContactTopologyEntity<Real>::connection_data *node_info = element->NodeInfo();
      POSTCONDITION( node_info );
      for(int j=0 ; j<element->Nodes_Per_Element() ; ++j){
	ContactNode<Real>* node = static_cast<ContactNode<Real> *>(NodeList()->Find(&node_info[j]));
	if (!node) {
	  //
	  //  Add the node to the import list
	  //
	  PRECONDITION(node_info[j].owner !=my_proc_id);
          import_local_ids[node_info[j].owner].push_back(node_info[j].owner_proc_array_index);
	}
      }
    }
  }
  vector< vector<int> > export_local_ids(total_num_procs);
  //
  //  Exchange the import and export lists to create a comm map
  //
  ACME::Parallel_Data_Exchange(import_local_ids, export_local_ids, SearchComm);

  ContactZoltanComm GhostOthers_ZoltanComm(ContactZoltanComm::ZOLTAN_EXPORT);
  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    vector<int> &current_export_ids = export_local_ids[iproc];
    for(int inode = 0; inode < current_export_ids.size(); ++inode) {
      zoltanLID.ZoltanLID(CT_NODE, current_export_ids[inode], zoltan_lid);
      ContactHostGlobalID GID = nodes[current_export_ids[inode]]->Global_ID();
      zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
      GhostOthers_ZoltanComm->Add_Export( zoltan_lid, zoltan_gid, iproc );  
    }
  }


  num_toplevel_import   = -1;
  import_toplevel_gids  = NULL;
  import_toplevel_lids  = NULL;
  import_toplevel_procs = NULL;
  num_toplevel_export   = GhostOthers_ZoltanComm->Num_Export();
  export_toplevel_gids  = GhostOthers_ZoltanComm->Export_GIDS();
  export_toplevel_lids  = GhostOthers_ZoltanComm->Export_LIDS();
  export_toplevel_procs = GhostOthers_ZoltanComm->Export_Procs();

  //================================================
  // Migrate all the off processor derived objects
  //================================================
#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->owner_help_migrate_time );
#endif
  search->Get_Zoltan()->Set_GhostingExportCallBacks(search);
  search->Get_Zoltan()->Help_Migrate(num_toplevel_import,  import_toplevel_gids, 
				     import_toplevel_lids, import_toplevel_procs, 
				     num_toplevel_export,  export_toplevel_gids,
				     export_toplevel_lids, export_toplevel_procs);
  //
  //  Some point in the future may want to cut out zoltan entiriely, possible to do 
  //  migration operations much faster ourselves.
  //
  //  Perform the nodal sizing and packing operations
  //
  vector< vector<char> > send_buffer;
  const int state1 = 1;  
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    vector<int> &current_export_ids = export_local_ids[iproc];
    int num_to_export = current_export_ids.size();
    //
    //  Determine the size in bytes of the export buffer
    //
    int export_size = 0;
    for(int inode = 0; inode < num_to_export; ++inode) {
      int node_index = current_export_ids[inode];
      ContactNode<Real> *node = nodes[node_index];
      export_size += node->Size(state1);
    }
    //
    //  Allocate the export buffer
    //
    vector<char> &current_send_buffer = send_buffer[iproc];
    current_send_buffer[iproc].resize(export_size);
    char* send_buffer_ptr = &(current_send_buffer[0]);
    //
    //  Pack the export array
    //
    for(int inode = 0; inode < num_to_export; ++inode) {
      int node_index = current_export_ids[inode];
      ContactNode<Real> *node = nodes[node_index];
      node->Pack(send_buffer_ptr, -2);
      buf += node->Size(state1);
    }
  }
  */

  //
  //  Regular Zoltan based method, construct import lists and let Zoltan do all ghosting calcs
  //
  int num_face_import = 0;
  int num_elem_import = 0;
  for (int i=0; i<number_of_face_blocks; ++i) {
    num_face_import += ghosted_face_blocks[i]->Number_of_Faces();
  }
  for (int i=0; i<number_of_element_blocks; ++i) {
    num_elem_import += ghosted_element_blocks[i]->Number_of_Elements();
  }
  if (GhostOthers_ZoltanComm!=NULL) {
    GhostOthers_ZoltanComm->Initialize( ContactZoltanComm::ZOLTAN_IMPORT);
  } else {
    GhostOthers_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT);
  }
                                              
  for(int i=0 ; i<number_of_face_blocks ; ++i ){
    ContactTopologyEntity<Real>* entity;
    ContactBlockEntityList* block_face_list = ghosted_face_blocks[i]->FaceList();
    block_face_list->IteratorStart();
    while( (entity=block_face_list->IteratorForward()) ){
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
      POSTCONDITION( node_info );
      for(int j=0 ; j<face->Nodes_Per_Face() ; ++j){
	ContactNode<Real>* node = static_cast<ContactNode<Real> *>(NodeList()->Find(&node_info[j]));
	if (!node) {
	  if (node_info[j].owner!=my_proc_id) {
	    zoltan_pid = node_info[j].owner;
	    ContactHostGlobalID GID( node_info[j].host_gid[0], 
				     node_info[j].host_gid[1] );
	    zoltanLID.ZoltanLID(CT_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
	    zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
	    GhostOthers_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
	  }
	}
      }
    }
  }

  for(int i=0 ; i<number_of_element_blocks ; ++i ){
    ContactTopologyEntity<Real>* entity;
    ContactBlockEntityList* block_element_list = element_blocks[i]->ElemList();
    block_element_list->IteratorStart();
    while( (entity=block_element_list->IteratorForward()) ){
      ContactElement* element = static_cast<ContactElement*>(entity);
        
      ContactTopologyEntity<Real>::connection_data *node_info = element->NodeInfo();
      POSTCONDITION( node_info );
      for(int j=0 ; j<element->Nodes_Per_Element() ; ++j){
	ContactNode<Real>* node = static_cast<ContactNode<Real> *>(NodeList()->Find(&node_info[j]));
	if (!node) {
	  if (node_info[j].owner!=my_proc_id) {
	    zoltan_pid = node_info[j].owner;
	    ContactHostGlobalID GID( node_info[j].host_gid[0], 
				     node_info[j].host_gid[1] );
	    zoltanLID.ZoltanLID(CT_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
	    zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
	    GhostOthers_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
	  }
	}
      }
    }
  }

  num_toplevel_import   = GhostOthers_ZoltanComm->Num_Import();
  import_toplevel_gids  = GhostOthers_ZoltanComm->Import_GIDS();
  import_toplevel_lids  = GhostOthers_ZoltanComm->Import_LIDS();
  import_toplevel_procs = GhostOthers_ZoltanComm->Import_Procs();
  num_toplevel_export   = -1;
  export_toplevel_gids  = NULL;
  export_toplevel_lids  = NULL;
  export_toplevel_procs = NULL;

  //================================================
  // Migrate all the off processor derived objects
  //================================================

#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->owner_help_migrate_time );
#endif
  search->Get_Zoltan()->Set_GhostingImportCallBacks(search);
  search->Get_Zoltan()->Help_Migrate(num_toplevel_import,  import_toplevel_gids, 
				     import_toplevel_lids, import_toplevel_procs, 
				     num_toplevel_export,  export_toplevel_gids,
				     export_toplevel_lids, export_toplevel_procs);
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->owner_help_migrate_time );
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  postream<<"    Data Sent by Zoltan when ghosting topology...\n";
  postream<<"      nodes:  "<<search->Bytes_For_Nodes()/1024.0<<" Kb\n";
  postream<<"      faces:  "<<search->Bytes_For_Faces()/1024.0<<" Kb\n";
  postream<<"      elems:  "<<search->Bytes_For_Elems()/1024.0<<" Kb\n";
#endif

#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->secondary_owner_migration_time );
#endif

#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->secondary_connect_mesh_time );
  search->Timer()->Start_Timer( search->secondary_connect_mesh_phase2_time );
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream<<"Before ghosting:\n";
  postream<<"  number_of_nodes    = "<<number_of_nodes<<"\n";
  postream<<"  number_of_faces    = "<<number_of_faces<<"\n";
  postream<<"  number_of_elements = "<<number_of_elements<<"\n";
#endif
  
  node_list->BuildList(node_blocks, ghosted_node_blocks, number_of_node_blocks,
		       no_parallel_consistency==ContactSearch::INACTIVE);
  face_list->BuildList(face_blocks, ghosted_face_blocks, number_of_face_blocks,
		       no_parallel_consistency==ContactSearch::INACTIVE);
  elem_list->BuildList(element_blocks, ghosted_element_blocks, number_of_element_blocks,
		       no_parallel_consistency==ContactSearch::INACTIVE);
  number_of_nodes    = node_list->NumEntities();
  number_of_faces    = face_list->NumEntities();
  number_of_elements = elem_list->NumEntities();
  nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  elems = reinterpret_cast<ContactElement**>(elem_list->EntityList());
    
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->secondary_connect_mesh_phase2_time );
  search->Timer()->Start_Timer( search->secondary_connect_mesh_phase3_time );
#endif
    
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream<<"After ghosting:\n";
  postream<<"  number_of_nodes    = "<<number_of_nodes<<"\n";
  postream<<"  number_of_faces    = "<<number_of_faces<<"\n";
  postream<<"  number_of_elements = "<<number_of_elements<<"\n";
#endif

  //===============================
  // Reset final enitity ownership 
  //===============================
    
  for (int i=0; i<number_of_nodes; ++i) {
    if (nodes[i]->Secondary_Owner()==my_proc_id)  {
      nodes[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    } else {
      nodes[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
    }
  }
    
  for (int i=0; i<number_of_faces; ++i) {
    if (faces[i]->Secondary_Owner()==my_proc_id)  {
      faces[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    } else {
      faces[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
    }
  }
    
  for (int i=0; i<number_of_elements; ++i) {
    if (elems[i]->Secondary_Owner()==my_proc_id)  {
      elems[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
    } else {
      elems[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
    }
  }
  
#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->secondary_connect_mesh_phase3_time );
  search->Timer()->Start_Timer( search->secondary_connect_mesh_phase4_time );
#endif

  //====================================================================
  // loop thru all the ghosted faces and connect the nodes to the faces
  //====================================================================
  for (int i=0; i<number_of_face_blocks; ++i) {
    ContactBlockEntityList* block_list = ghosted_face_blocks[i]->FaceList();
    block_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
      POSTCONDITION( node_info );
      for(int j=0 ; j<face->Nodes_Per_Face() ; ++j){
	ContactNode<Real>* node = static_cast<ContactNode<Real> *>(node_list->Find(&node_info[j]));
	POSTCONDITION(node);
	face->ConnectNode( j, node );
      }
    }
  }
  //==========================================================================
  // loop thru all the ghosted elements and connect the nodes to the elements
  //==========================================================================
  for (int i=0; i<number_of_element_blocks; ++i) {
    ContactBlockEntityList* block_list = ghosted_element_blocks[i]->ElemList();
    block_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      ContactTopologyEntity<Real>::connection_data *node_info = element->NodeInfo();
      POSTCONDITION( node_info );
      for(int j=0 ; j<element->Nodes_Per_Element() ; ++j){
	ContactNode<Real>* node = static_cast<ContactNode<Real> *>(node_list->Find(&node_info[j]));
	POSTCONDITION(node);
	element->ConnectNode( j, node );
      }
    }
  }
    
  if (number_of_faces>0) {
    //====================================================================
    // loop thru all the faces and reconnect the faces to the nodes
    //====================================================================
    for (int i=0; i<number_of_nodes; ++i) {
      nodes[i]->Delete_Face_Connections();
    }
    for (int i=0; i<number_of_faces; ++i) {
      ContactFace<Real>* face = faces[i];
      for(int k=0 ; k<face->Nodes_Per_Face() ; ++k ){
        ContactNode<Real>* node = face->Node(k);
        node->Connect_Face(face );
      }
    }
    for (int i=0; i<number_of_nodes; ++i) {
      nodes[i]->SortConnectedFaces();
    }
  }
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->secondary_connect_mesh_phase4_time );
  search->Timer()->Start_Timer( search->secondary_connect_mesh_phase5_time );
#endif
                                        
  if (search->Do_NodeFace_Search()) {
    // connect all the ghosted faces to the appropriate interaction
    for (int i=0; i<number_of_nodes; ++i) {
      bool found_invalid = false;
      ContactNode<Real>* node = nodes[i];
      ContactNodeEntityInteraction** interactions = 
	node->Get_NodeEntity_Interactions(1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
	if (cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued()) {
	  cnei->Connect_Entity(this);
          if (cnei->Entity()==NULL) {
            // MWG: We really need to return a fatal error code here and abort!
            //      Since we are explicitly ghosting the master faces for each interaction,
            //      this condition should never occurr.
            Real* position = node->Variable(CURRENT_POSITION);
            Real* rem_gap  = node->Variable(REMAINING_GAP);
            if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
	      static int printlimit = 0;
	      if (printlimit < 5) {					       
		ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
		std::cout<<"P"<<my_proc_id<<": In ContactTopology::MigrateExportedData(),step "<<search->StepNumber()<<"\n"
			 <<"  error with interaction "<<j
			 <<", cannot connect master  face "<<"("<<cnfi->FaceEntityData()->host_gid[1]<<", "<<cnfi->FaceEntityData()->host_gid[0]<<")\n"
			 <<"    node "<<node->Global_ID()<<" exodus_id:"<< node->Exodus_ID()<<"\n"
			 <<"    position      = ("<<position[0]<<", "<<position[1]<<", "<<position[2]<<")\n"
			 <<"    remaining_gap = ("<<rem_gap[0]<<", "<<rem_gap[1]<<", "<<rem_gap[2]<<")\n"<<std::flush;
		printlimit ++;
	      }
	      node->Delete_NodeEntity_Interaction(cnei,1);
              found_invalid = true;

	      // What we should do is die here. Unfortunately, we 
	      // risk branch stability with this, so for now 
	      // we just report the error and move on.
	      // err_code = ContactSearch::INTERNAL_ERROR;
            }
	  }
        }
      }
    }
  }

  for(int i=0 ; i<number_of_node_blocks ; ++i ){
    ContactTopologyEntity<Real>* entity = NULL;
    ContactBlockEntityList* block_list = ghosted_node_blocks[i]->NodeList();
    block_list->IteratorStart();
    while( (entity=block_list->IteratorForward()) ){
      entity->SetContextBit(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH);
    }
  }                                     
  for(int i=0 ; i<number_of_face_blocks ; ++i ){
    ContactTopologyEntity<Real>* entity = NULL;
    ContactBlockEntityList* block_list = ghosted_face_blocks[i]->FaceList();
    block_list->IteratorStart();
    while( (entity=block_list->IteratorForward()) ){
      entity->SetContextBit(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH);
    }
  }
  for(int i=0 ; i<number_of_element_blocks ; ++i ){
    ContactTopologyEntity<Real>* entity = NULL;
    ContactBlockEntityList* block_list = element_blocks[i]->ElemList();
    block_list->IteratorStart();
    while( (entity=block_list->IteratorForward()) ){
      entity->SetContextBit(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH);
    }
  }
  
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->secondary_connect_mesh_phase5_time );
  search->Timer()->Stop_Timer( search->secondary_connect_mesh_time );
#endif

}

#endif

void 
ContactTopology::RebuildTopologyEntityList( ContactSearch::Search_Option_Status NPC_Status )
{
  no_parallel_consistency = NPC_Status;
  node_list->BuildList(Node_Blocks(), Number_of_Node_Blocks(),
                       no_parallel_consistency==ContactSearch::INACTIVE);
  edge_list->BuildList(Edge_Blocks(), Number_of_Edge_Blocks(),
                       no_parallel_consistency==ContactSearch::INACTIVE);
  face_list->BuildList(Face_Blocks(), Number_of_Face_Blocks(),
                       no_parallel_consistency==ContactSearch::INACTIVE);
  elem_list->BuildList(Element_Blocks(), Number_of_Element_Blocks(),
                       no_parallel_consistency==ContactSearch::INACTIVE);
  edge_list->SortByNodeGID();
}

