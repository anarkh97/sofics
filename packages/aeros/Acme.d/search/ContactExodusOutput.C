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


#ifndef CONTACT_NO_EXODUS_OUTPUT

#include "ContactTopology.h"
#include "ContactSearch.h"
#include "CString.h"
#include "exodusII.h"
#include "ContactElement.h"
#include "ContactErrors.h"
#include "ContactFace.h"
#include "ContactNode.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactFaceBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactSearchData.h"
#include "ContactSymComm.h"
#include "Contact_Communication.h"
#include "ContactEnforcement.h"

#ifndef CONTACT_NO_MPI
#include "ne_nemesisI.h"
#endif

#include <iostream>
#include <cstdio>
#include <cstring>

ContactSearch::ContactErrorCode
ContactTopology::Exodus_Output( int Exodus_ID, Real Time, 
                                ContactErrors* error_handler,
                                ContactSearchData* search_data,
                                const ContactSearch::Search_Option_Status&
                                normal_smoothing,
                                const ContactSearch::Search_Option_Status&
                                multiple_interaction,
                                const Real& sharp_smooth_ang, 
                                const Real& normal_smoothing_distance, 
                                const ContactSearch::Smoothing_Resolution &
                                smoothing_resolution,
				const ContactSearch::Search_Option_Status&
				compute_node_areas,
                                ContactEnforcement* enforcement )
{
  int i,j,k,ierr,index;
  ContactFaceFaceInteraction<Real>* cffi;
  ContactFaceCoverageInteraction* cfci;
  ContactInteractionEntity<Real>* interaction;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_node_list->EntityList());
  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(primary_edge_list->EntityList());
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_face_list->EntityList());
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(primary_elem_list->EntityList());

#ifndef CONTACT_NO_MPI
  bool PARALLEL  = false;
  int  num_procs = contact_number_of_processors( SearchComm );
  int  my_proc   = contact_processor_number( SearchComm );
  if( num_procs > 1 ) PARALLEL = true;
#endif
  
  int num_nni_l = 0;
  int num_nfi_l = 0;
  int num_nsi_l = 0;
  int max_nni_l = 0;
  int max_nfi_l = 0;
  int max_nsi_l = 0;
  for (i=0; i<number_of_primary_nodes; ++i) {
    num_nni_l += Nodes[i]->Number_NodeNode_Interactions();
    num_nfi_l += Nodes[i]->Number_NodeFace_Interactions();
    num_nsi_l += Nodes[i]->Number_NodeSurface_Interactions();
    max_nni_l = std::max(max_nni_l,Nodes[i]->Number_NodeNode_Interactions());
    max_nfi_l = std::max(max_nfi_l,Nodes[i]->Number_NodeFace_Interactions());
    max_nsi_l = std::max(max_nsi_l,Nodes[i]->Number_NodeSurface_Interactions());
  }
  int num_nni = contact_global_sum(num_nni_l, SearchComm );
  int num_nfi = contact_global_sum(num_nfi_l, SearchComm );
  int num_nsi = contact_global_sum(num_nsi_l, SearchComm );
  int max_nni = contact_global_maximum(max_nni_l, SearchComm );
  int max_nfi = contact_global_maximum(max_nfi_l, SearchComm );
  int max_nsi = contact_global_maximum(max_nsi_l, SearchComm );
  
  int num_eei_l = 0;
  int max_eei_l = 0;
  for (i=0; i<number_of_primary_elements; ++i) {
    num_eei_l += Elements[i]->Number_ElementElement_Interactions();
    max_eei_l  = std::max(max_nfi_l,Elements[i]->Number_ElementElement_Interactions());
  }
  int num_eei = contact_global_sum(num_eei_l, SearchComm );
  int max_eei = contact_global_maximum(max_eei_l, SearchComm );
  
  int num_ffi_blocks  = 0;
  int num_ffi_elems   = 0;
  int num_ffi_nodes   = 0;
  int num_ffi_l       = 0;
  int max_ffi_l       = 0;
  int max_ffi_verts_l = 0;
  for (i=0; i<number_of_primary_faces; ++i) {
    num_ffi_l += Faces[i]->Number_FaceFace_Interactions();
    max_ffi_l  = std::max(max_ffi_l,Faces[i]->Number_FaceFace_Interactions());
    ContactInteractionDLL<Real>* interactions = Faces[i]->Get_FaceFace_Interactions();
    if(interactions == NULL) continue;
    interactions->IteratorStart();
    while (interaction=interactions->IteratorForward()){
      cffi = static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
      max_ffi_verts_l = std::max(max_ffi_verts_l,cffi->NumEdges()+1);
      num_ffi_elems += cffi->NumEdges();
      num_ffi_nodes += cffi->NumEdges()+1;
    }
  }
  int num_ffi       = contact_global_sum(num_ffi_l, SearchComm );
  int max_ffi       = contact_global_maximum(max_ffi_l, SearchComm );
  int max_ffi_verts = contact_global_maximum(max_ffi_verts_l, SearchComm );
#ifndef CONTACT_NO_MPI
  int num_ffi_elems_global = contact_global_sum( num_ffi_elems, SearchComm );
#endif
  if (num_ffi>0) num_ffi_blocks = 1;
  
  int num_fci_blocks  = 0;
  int num_fci_elems   = 0;
  int num_fci_nodes   = 0;
  int num_fci_l       = 0;
  int max_fci_l       = 0;
  int max_fci_verts_l = 0;
  for (i=0; i<number_of_primary_faces; ++i) {
    num_fci_l += Faces[i]->Number_FaceCoverage_Interactions();
    max_fci_l  = std::max(max_fci_l,Faces[i]->Number_FaceCoverage_Interactions());
    ContactInteractionDLL<Real>* interactions = Faces[i]->Get_FaceCoverage_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while (interaction=interactions->IteratorForward()){
        cfci = static_cast<ContactFaceCoverageInteraction*>(interaction);
        max_fci_verts_l = std::max(max_fci_verts_l,cfci->NumVertices());
        num_fci_elems += cfci->NumVertices();
        num_fci_nodes += cfci->NumVertices();
      }
    }
  }
  int num_fci       = contact_global_sum(num_fci_l, SearchComm );
  int max_fci       = contact_global_maximum(max_fci_l, SearchComm );
  int max_fci_verts = contact_global_maximum(max_fci_verts_l, SearchComm );
#ifndef CONTACT_NO_MPI
  int num_fci_elems_global = contact_global_sum( num_fci_elems, SearchComm );
#endif
  if (num_fci>0) num_fci_blocks = 1;

  // Exodus won't support a file with no nodes or elements.
  // We have three possible processor node/element configurations
  //   1) No Nodes or Elements
  //   2) Nodes but no Elements (talk to KHB for a description)
  //   3) Nodes & Elements
  //
  // Dont_Have_Nodes_and_Element = 0  Don't add a sphere element on this proc
  //                             = 1  Add a sphere element on this proc
  //
  // Need_Sphere_Element_Block   = 0 Don't need to add the sphere element block
  //                             > 0 # of global nodes/elements to be added
  //
  int Dont_Have_Nodes_and_Elements = 0;
  if( (!number_of_primary_faces && !number_of_primary_elements) || 
      !number_of_primary_nodes ) Dont_Have_Nodes_and_Elements = 1;
  int Need_Sphere_Element_Block = 
    contact_global_sum( Dont_Have_Nodes_and_Elements,SearchComm );

  //============================================================================
  // Calculate the global and local number of nodes, 
  // elements, and element blocks
  //============================================================================
    
  // Count the number of global nodes 
  // (the global sum of owned nodes)
  int number_nodes_owned = 0;
  for (i=0; i<number_of_primary_nodes; ++i) {
    if( Nodes[i]->Ownership() == ContactTopologyEntity<Real>::OWNED )
      number_nodes_owned++;
  }
  number_nodes_owned += num_ffi_nodes;
  number_nodes_owned += num_fci_nodes;
  int number_nodes_global = contact_global_sum(number_nodes_owned, SearchComm);
  number_nodes_global += Need_Sphere_Element_Block;

  // Count the number of global edges 
  // (the global sum of owned edges)
  int number_edges_owned = 0;
  for (i=0; i<number_of_primary_edges; ++i) {
    if( Edges[i]->Ownership() == ContactTopologyEntity<Real>::OWNED )
      number_edges_owned++;
  }
#ifndef CONTACT_NO_MPI
  int number_edges_global = contact_global_sum(number_edges_owned, SearchComm);
  // Count the number of global faces & elements
  // (currently this is just the global sum)
  int number_faces_global    = contact_global_sum(number_of_primary_faces, SearchComm);
  int number_elements_global = contact_global_sum(number_of_primary_elements, SearchComm );
#endif

  // Count the number of global SPH elements 
  // (the global sum of owned sph elements)
  int num_node_sph_elems = 0;
  int num_sph_elem_blks = 0;
  for( i=1 ; i<number_of_node_blocks ; ++i ){
    if( node_blocks[i]->Type() == ContactSearch::POINT ){
      num_sph_elem_blks++;
      ContactNode<Real>** BlockNodes = 
        reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
      for (j=0; j<primary_node_list->BlockNumEntities(i); ++j) {
	if( BlockNodes[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
	  num_node_sph_elems++;
      }
    }
  }
#ifndef CONTACT_NO_MPI
  int num_node_sph_elems_global = contact_global_sum( num_node_sph_elems, SearchComm );
#endif  
  int number_of_elem_blocks = number_of_element_blocks +
                              number_of_face_blocks + 
                              number_of_edge_blocks +
                              num_sph_elem_blks + 
                              num_ffi_blocks +
                              num_fci_blocks;
  int number_of_elems       = number_of_primary_elements +
                              number_of_primary_faces + 
                              number_edges_owned +
                              num_node_sph_elems + 
                              num_ffi_elems + 
                              num_fci_elems;
  
  int num_proc_nodes = number_of_primary_nodes;
  if( Need_Sphere_Element_Block ){
    number_of_elem_blocks++;
    if( Dont_Have_Nodes_and_Elements ){
      num_proc_nodes++;
      number_of_elems++;
    }
  }
  int num_base_nodes = num_proc_nodes;
  num_proc_nodes    += num_ffi_nodes;
  num_proc_nodes    += num_fci_nodes;
  int number_of_node_sets = 0;
  int number_of_side_sets = 0;
  
  //============================================================================
  // I N I T I A L I Z E   T H E   D A T A B A S E
  //============================================================================
  char Title[] = "ACME Library Output";
  ierr = ex_put_init( Exodus_ID, Title, dimensionality, num_proc_nodes,
                      number_of_elems, number_of_elem_blocks,
                      number_of_node_sets, number_of_side_sets );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_init",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  
#ifndef CONTACT_NO_MPI
  if( PARALLEL ){

    contact_global_sync(SearchComm);
  
    char* file_type = (char*)"p";
    ne_put_init_info( Exodus_ID, num_procs, 1, file_type );

    // Output the global counts
    int num_elem_global = number_edges_global + 
                          number_faces_global +
                          number_elements_global +
                          num_ffi_elems_global + 
                          num_fci_elems_global +
                          num_node_sph_elems_global + 
                          Need_Sphere_Element_Block;
    ne_put_init_global( Exodus_ID, number_nodes_global, num_elem_global,
                        number_of_elem_blocks, 0, 0 );
    
    // Compute the global counts for each element block
    int* global_counts = new int[number_of_elem_blocks];
    int* local_counts  = new int[number_of_elem_blocks];
    std::memset( local_counts, 0, number_of_elem_blocks*sizeof(int) );
    for( i=0 ; i<number_of_element_blocks ; ++i )
      local_counts[i] = element_blocks[i]->Number_of_Elements();
    for( i=0 ; i<number_of_face_blocks ; ++i )
      local_counts[i] = face_blocks[i]->Number_of_Faces();
    // for edges we must decide if we are responsible for them so we can't
    // just use the block count
    for( i=0 ; i<number_of_edge_blocks ; ++i ) {
      k = number_of_face_blocks+i;
      local_counts[k] = 0;
      ContactEdge<Real>** BlockEdges = 
        reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
      for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
        if( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
          local_counts[k]++;
      }
    }
    for( i=1 ; i<number_of_node_blocks ; ++i ){
      if( node_blocks[i]->Type() == ContactSearch::POINT ){
	k = number_of_element_blocks+
            number_of_face_blocks+
            number_of_edge_blocks+(i-1);
        ContactNode<Real>** BlockNodes = 
          reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
        for (j=0; j<primary_node_list->BlockNumEntities(i); ++j) {
	  if( BlockNodes[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
	    local_counts[k]++;
	}
      }
    }
    if( Need_Sphere_Element_Block ){
      if( Dont_Have_Nodes_and_Elements )
        local_counts[number_of_element_blocks+
                     number_of_face_blocks+
                     number_of_edge_blocks+
                     num_sph_elem_blks] = 1;
      else 
        local_counts[number_of_element_blocks+
                     number_of_face_blocks+
                     number_of_edge_blocks+
                     num_sph_elem_blks] = 0;     
    }
    for( i=0 ; i<num_ffi_blocks ; ++i ){
      k = number_of_element_blocks+
          number_of_face_blocks+
          number_of_edge_blocks+
          num_sph_elem_blks;
      for( j=0 ; j<num_ffi_elems ; ++j ){
        local_counts[k]++;
      }
    }
    for( i=0 ; i<num_fci_blocks ; ++i ){
      k = number_of_element_blocks+
          number_of_face_blocks+
          number_of_edge_blocks+
          num_sph_elem_blks+
          num_ffi_blocks;
      for( j=0 ; j<num_fci_elems ; ++j ){
        local_counts[k]++;
      }
    }

    contact_global_sum( local_counts, global_counts, number_of_elem_blocks,
                        SearchComm );

    int* block_ids = new int[number_of_elem_blocks];
    for( i=0 ; i<number_of_elem_blocks ; ++i ) block_ids[i] = i+1;

    // Set the sphere block_id to 50000 (if it exists)
    if( Need_Sphere_Element_Block ){
      PRECONDITION( number_of_element_blocks < 50000 );
      block_ids[number_of_elem_blocks-num_ffi_blocks-num_fci_blocks-1] = 50000;
    }
    // Set the ffi block_id to 60000 (if it exists)
    if( num_ffi>0 ){
      PRECONDITION( number_of_elem_blocks < 60000 );
      block_ids[number_of_elem_blocks-num_fci_blocks-1] = 60000;
    }
    // Set the ffi block_id to 70000 (if it exists)
    if( num_fci>0 ){
      PRECONDITION( number_of_elem_blocks < 70000 );
      block_ids[number_of_elem_blocks-1] = 70000;
    }

    ne_put_eb_info_global( Exodus_ID, block_ids, global_counts );

    // Count the number of std::internal, border & external nodes
    ContactTopologyEntity<Real>** node_comm_list;
    for (i=0; i<number_of_primary_nodes; ++i) {
      Nodes[i]->temp_tag = 0;
    }
    for( i=0 ; i<Node_SymComm->Num_Comm_Partners() ; ++i ){
      node_comm_list = Node_SymComm->Entity_List(i);
      for( j=0 ; j<Node_SymComm->Num_to_Proc( i ) ; ++j )
        node_comm_list[j]->temp_tag++;
    }
    int num_external_nodes = 0;
    int num_internal_nodes = 0;
    int num_border_nodes   = 0;
    for (i=0; i<number_of_primary_nodes; ++i) {
      if( Nodes[i]->temp_tag == 0 )
        num_internal_nodes++;
      else
        num_border_nodes++;
    }
    if( Dont_Have_Nodes_and_Elements ) num_internal_nodes++;
    num_internal_nodes += num_ffi_nodes;
    num_internal_nodes += num_fci_nodes;

    // Count the number of std::internal & border elements.
    // "Elements" (ie, faces) are border if one of their edges is in
    // an edge communication list.
    int num_internal_elems = 0;
    int num_border_elems = 0;
    for (i=0; i<number_of_primary_faces; ++i) {
      Faces[i]->temp_tag = 0;
    }
    for (i=0; i<number_of_primary_edges; ++i) {
      Edges[i]->temp_tag = 0;
    }
    ContactTopologyEntity<Real>** edge_comm_list;
    for( i=0 ; i<Edge_SymComm->Num_Comm_Partners() ; ++i ){
      edge_comm_list = Edge_SymComm->Entity_List(i);
      for( j=0 ; j<Edge_SymComm->Num_to_Proc( i ) ; ++j )
        edge_comm_list[j]->temp_tag++;
    }
    for (i=0; i<number_of_primary_edges; ++i) {
      if( Edges[i]->temp_tag != 0){ 
        for( j=0 ; j<Edges[i]->Number_Face_Connections() ; ++j )
          Edges[i]->Face(j)->temp_tag = 1;
      }
    }
    for (i=0; i<number_of_primary_faces; ++i) {
      if( Faces[i]->temp_tag == 0 )
        num_internal_elems++;
      else
        num_border_elems++;
    }
    num_internal_elems += number_edges_owned;
    if( number_of_node_blocks > 1 ){
      for( i=i ; i<number_of_node_blocks ; ++i ){
	if( node_blocks[i]->Type() == ContactSearch::POINT ){
          ContactNode<Real>** BlockNodes = 
            reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
          for (j=0; j<primary_node_list->BlockNumEntities(i); ++j) {
	    if( BlockNodes[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
	      num_internal_elems++;
	  }
	}
      }
    }
    if( Dont_Have_Nodes_and_Elements ) num_internal_elems++;
    num_internal_elems += num_ffi_elems;
    num_internal_elems += num_fci_elems;

    // Output the global parameters
    int num_node_cmaps = Node_SymComm->Num_Comm_Partners();
    int num_elem_cmaps = Edge_SymComm->Num_Comm_Partners();

    ne_put_loadbal_param( Exodus_ID, num_internal_nodes, num_border_nodes,
                          num_external_nodes, num_internal_elems, 
                          num_border_elems, num_node_cmaps, num_elem_cmaps,
                          my_proc );

    // Create the node maps
    int* node_mapi = new int[num_internal_nodes];
    int* node_mapb = new int[num_border_nodes];
    int* node_mape = NULL;
    int icount_i = 0;
    int icount_b = 0;
    for (i=0; i<number_of_primary_nodes; ++i) {
      if( Nodes[i]->temp_tag == 0 )
        //node_mapi[icount_i++] = Nodes[i]->ProcArrayIndex()+1;
        node_mapi[icount_i++] = Nodes[i]->fcs_index+1;
      else
        //node_mapb[icount_b++] = Nodes[i]->ProcArrayIndex()+1;
        node_mapb[icount_b++] = Nodes[i]->fcs_index+1;
    }
    int offset = number_of_primary_nodes+1;
    if( Dont_Have_Nodes_and_Elements ) {
      node_mapi[icount_i++] = offset++;
    }
    for (i=0; i<num_ffi_nodes; ++i) {
      node_mapi[icount_i++] = offset++;
    }
    for (i=0; i<num_fci_nodes; ++i) {
      node_mapi[icount_i++] = offset++;
    }
    ne_put_node_map( Exodus_ID, node_mapi, node_mapb, node_mape, my_proc );
    delete [] node_mapi;
    delete [] node_mapb;

    // Create the element maps
    int* elem_mapi = new int[num_internal_elems];
    int* elem_mapb = new int[num_border_elems];
    icount_i = 0;
    icount_b = 0;
    for (i=0; i<number_of_primary_faces; ++i) {
      if( Faces[i]->temp_tag == 0 )
        //elem_mapi[icount_i++] = Faces[i]->ProcArrayIndex()+1;
        elem_mapi[icount_i++] = Faces[i]->fcs_index+1;
      else
        //elem_mapb[icount_b++] = Faces[i]->ProcArrayIndex()+1;
        elem_mapb[icount_b++] = Faces[i]->fcs_index+1;
    }
    
    // NOTE from asg: This may be a problem. Here we are determining
    //  the mapping id from the contact index. For this processor,
    //  all edges are numbered without any holes in the numbering. 
    //  However, holes in the numbering appear when we exclude those
    //  we don't own.  If there is a need for local numberings in
    //  the exodus elem map to have no numbering holes, then this code
    //  is not sufficient.
    for (i=0; i<number_of_primary_edges; ++i) {
      if ( Edges[i]->Ownership() == ContactTopologyEntity<Real>::OWNED) 
        elem_mapi[icount_i++] = 
          //number_of_primary_faces + Edges[i]->ProcArrayIndex()+1;
          number_of_primary_faces + Edges[i]->fcs_index+1;
    }
    offset = number_of_primary_faces + 
             number_edges_owned +
             num_node_sph_elems + 1;
    if( Dont_Have_Nodes_and_Elements ) {
      elem_mapi[icount_i++] = offset++;
    }
    for (i=0; i<num_ffi_elems; ++i) {
      elem_mapi[icount_i++] = offset++;
    }
    for (i=0; i<num_fci_elems; ++i) {
      elem_mapi[icount_i++] = offset++;
    }
    ne_put_elem_map( Exodus_ID, elem_mapi, elem_mapb, my_proc );
    delete [] elem_mapi;
    delete [] elem_mapb;

    // Create the node & element communication maps
    int* node_cmap_node_cnts = new int[num_node_cmaps];
    int* node_cmap_ids = new int[num_node_cmaps];
    int* elem_cmap_elem_cnts = new int[num_elem_cmaps];
    int* elem_cmap_ids = new int[num_elem_cmaps];
    for( i=0 ; i<num_node_cmaps ; ++i ){
      node_cmap_ids[i] = Node_SymComm->Comm_Proc_ID( i );
      node_cmap_node_cnts[i] = Node_SymComm->Num_to_Proc( i );
    }
    for( i=0 ; i<num_elem_cmaps ; ++i ){
      elem_cmap_ids[i] = Edge_SymComm->Comm_Proc_ID( i );
      elem_cmap_elem_cnts[i] = Edge_SymComm->Num_to_Proc( i );
    }

    ne_put_cmap_params( Exodus_ID, node_cmap_ids, node_cmap_node_cnts, 
                        elem_cmap_ids, elem_cmap_elem_cnts, my_proc );

    int* nodal_ids = new int[Node_SymComm->Size()];
    int* proc_ids = new int[Node_SymComm->Size()];
    for( i=0 ; i<Node_SymComm->Num_Comm_Partners() ; ++i ){
      int map_id = Node_SymComm->Comm_Proc_ID(i);
      node_comm_list = Node_SymComm->Entity_List(i);
      for( j=0 ; j<Node_SymComm->Num_to_Proc( i ) ; ++j ){
        //nodal_ids[j] = node_comm_list[j]->ProcArrayIndex()+1;
        nodal_ids[j] = node_comm_list[j]->fcs_index+1;
        proc_ids[j] = Node_SymComm->Comm_Proc_ID(i);
      }
      ne_put_node_cmap( Exodus_ID, map_id, nodal_ids, proc_ids, my_proc );
    }

    int* elem_ids = new int[Edge_SymComm->Size()];
    if( Edge_SymComm->Size() > Node_SymComm->Size() ){
      delete [] proc_ids;
      proc_ids = new int[Edge_SymComm->Size()];
    }
    int* side_ids = new int[Edge_SymComm->Size()];
    for( i=0 ; i<Edge_SymComm->Num_Comm_Partners() ; ++i ){
      int mapid = Edge_SymComm->Comm_Proc_ID(i);
      edge_comm_list = Edge_SymComm->Entity_List(i);
      for( j=0 ; j<Edge_SymComm->Num_to_Proc( i ) ; ++j ){
        ContactEdge<Real>* e = reinterpret_cast<ContactEdge<Real>*>(edge_comm_list[j]);
        PRECONDITION( e->Number_Face_Connections() == 1 );
        ContactFace<Real>* f = e->Face(0);
        //elem_ids[j] = f->ProcArrayIndex()+1;
        elem_ids[j] = f->fcs_index+1;
        proc_ids[j] = Edge_SymComm->Comm_Proc_ID(i);
        // Now find the "face" of the shell element.  This corresponds
        // to the edge index + 1 for FORTRAN and + 2 to account for the
        // faces not the edges.
        for( k=0 ; k<f->Edges_Per_Face() ; ++k )
          if( f->Edge(k) == e ) break;
        POSTCONDITION( k < f->Edges_Per_Face() );
        side_ids[j] = k+3;
      }
      ne_put_elem_cmap( Exodus_ID,mapid,elem_ids,side_ids,proc_ids,my_proc );
    }
    
    delete [] local_counts;
    delete [] global_counts;
    delete [] block_ids;
    delete [] nodal_ids;
    delete [] proc_ids;
    delete [] side_ids;
    delete [] elem_ids;
    delete [] node_cmap_node_cnts;
    delete [] node_cmap_ids;
    delete [] elem_cmap_elem_cnts;
    delete [] elem_cmap_ids;

    // output global numbers for nodes and elements 
    int* proc_starting_count = new int[num_procs];
    int* proc_starting_cnt   = new int[num_procs];
  
  
    // Construct a unique numbering for all node_block nodes.
    // I will simply number the owned nodes on proc 0 from 1-N and
    // then assign the node numbers to nodes on proc 1 from N+1 to NN, etc.
    std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
    proc_starting_count[my_proc] = number_nodes_owned;
    if( Dont_Have_Nodes_and_Elements ) proc_starting_count[my_proc]++;
    contact_global_sum( proc_starting_count, proc_starting_cnt, num_procs,
                        SearchComm );
    // Put the "global Exodus ID" into node_num_map for owned nodes.
    // Set the value to zero for "ghost" nodes.  We will then swapadd this
    // array to get the values consistently known.
    int* node_num_map = new int[num_proc_nodes];
    int my_start = 1;
    for( i=0 ; i<my_proc ; ++i )
      my_start += proc_starting_cnt[i];
    for (i=0; i<number_of_primary_nodes; ++i) {
      if( Nodes[i]->Ownership() == ContactTopologyEntity<Real>::OWNED ) 
        node_num_map[i] = my_start++;
      else
        node_num_map[i] = 0;
    }
    for( i=number_of_primary_nodes ; i<num_proc_nodes ; ++i ){
        node_num_map[i] = my_start++;
    }
    contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                    node_num_map, 1 );
    //This statement doesn't seem to be necessary due to the loop 
    //above and is causing problems with the element death output
    //if( Dont_Have_Nodes_and_Elements ) node_num_map[0] = my_start;
    ex_put_node_num_map( Exodus_ID, node_num_map );
    delete [] node_num_map;
    
    
    // Construct a unique number of all face elements
    //   If we don't have any faces (and therefore edges) we will be
    //   adding on SPHERE element so use the max to get the std::right size.
    int  elem_num_map_size  = std::max(1,number_of_elems);
    int* elem_num_map       = new int[elem_num_map_size];
    int  local_elem_offset  = 0;
    int  global_elem_offset = 0;
    
    for( i=0 ; i<number_of_element_blocks ; ++i ){
      int nelements_in_block = element_blocks[i]->Number_of_Elements();
      int nglobal_elements_in_block = contact_global_sum( nelements_in_block, 
                                                          SearchComm );
      std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
      proc_starting_count[my_proc] = nelements_in_block;
      contact_global_sum( proc_starting_count, proc_starting_cnt, num_procs,
                          SearchComm );
      for( my_start=global_elem_offset+1, j=0 ; j<my_proc ; ++j ) {
        my_start += proc_starting_cnt[j];
      }
      for( j=0 ; j<nelements_in_block ; ++j ) {
        elem_num_map[local_elem_offset++] = global_elem_offset + my_start++;
      }
      global_elem_offset += nglobal_elements_in_block;
    }
    
    for( i=0 ; i<number_of_face_blocks ; ++i ){
      int nelements_in_block = face_blocks[i]->Number_of_Faces();
      int nglobal_elements_in_block = contact_global_sum( nelements_in_block, 
                                                          SearchComm );
      std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
      proc_starting_count[my_proc] = nelements_in_block;
      contact_global_sum( proc_starting_count, proc_starting_cnt, num_procs,
                          SearchComm );
      for( my_start=global_elem_offset+1, j=0 ; j<my_proc ; ++j ) {
        my_start += proc_starting_cnt[j];
      }
      for( j=0 ; j<nelements_in_block ; ++j ) {
        elem_num_map[local_elem_offset++] = global_elem_offset + my_start++;
      }
      global_elem_offset += nglobal_elements_in_block;
    }
    
    // Construct a unique numbering for all edges and add to face elements
    // I will simply number the owned edges on proc 0 from 1-N and
    // then assign the edge numbers to edges on proc 1 from N+1 to NN, etc.
    for( i=0 ; i<number_of_edge_blocks ; ++i ){
      int nowned_edges_in_block = 0;
      ContactEdge<Real>** BlockEdges = 
        reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
      for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
        if( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
          nowned_edges_in_block++;
      }
      int nglobal_edges_in_block = contact_global_sum( nowned_edges_in_block, 
                                                       SearchComm );
      std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
      proc_starting_count[my_proc] = nowned_edges_in_block;
      contact_global_sum( proc_starting_count, proc_starting_cnt, 
                          num_procs, SearchComm );
      for( my_start=global_elem_offset+1, j=0 ; j<my_proc ; ++j ) {
        my_start += proc_starting_cnt[j];
      }
      for (j=0; j<edge_blocks[i]->Number_of_Edges(); ++j) {
        if( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
          PRECONDITION(local_elem_offset<elem_num_map_size);
          elem_num_map[local_elem_offset++] = global_elem_offset + my_start++;
        }
      }
      global_elem_offset += nglobal_edges_in_block;
    }

    // Construct a unique numbering for the SPH elements (node blocks)
    for( i=1 ; i<number_of_node_blocks ; ++i ){
      if( node_blocks[i]->Type() == ContactSearch::POINT ){
	int nowned_nodes_in_block = 0;
        ContactNode<Real>** BlockNodes = 
          reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
        for (j=0; j<primary_node_list->BlockNumEntities(i); ++j) {
	  if( BlockNodes[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
	    nowned_nodes_in_block++;
	}
	int nglobal_nodes_in_block = contact_global_sum( nowned_nodes_in_block, 
							 SearchComm );
	std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
	proc_starting_count[my_proc] = nowned_nodes_in_block;
	contact_global_sum( proc_starting_count, proc_starting_cnt, 
			    num_procs, SearchComm );
	for( my_start=global_elem_offset+1, j=0 ; j<my_proc ; ++j ) {
	  my_start += proc_starting_cnt[j];
	}
        for (j=0; j<node_blocks[i]->Number_of_Nodes(); ++j) {
	  if( BlockNodes[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
	    elem_num_map[local_elem_offset++] = global_elem_offset + my_start++;
	  }
	}
	global_elem_offset += nglobal_nodes_in_block;
      }
    }

    // Construct a unique numbering for the SPHERE elements (if needed)
    if( Need_Sphere_Element_Block ){
      std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
      if( Dont_Have_Nodes_and_Elements ) proc_starting_count[my_proc] = 1;
      contact_global_sum( proc_starting_count, proc_starting_cnt,
                          num_procs, SearchComm );
      for( my_start=global_elem_offset+1, i=0 ; i<my_proc ; ++i ) {
        my_start += proc_starting_cnt[i];
      }
      if( Dont_Have_Nodes_and_Elements ) {
        elem_num_map[local_elem_offset++] = number_elements_global +
                                            number_faces_global +
                                            number_edges_global + 
                                            my_start;
      }
      global_elem_offset++;
    }
    
    // Construct a unique numbering for the FFI elements (if needed)
    if (num_ffi>0) {
      std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
      proc_starting_count[my_proc] = num_ffi_elems;
      contact_global_sum( proc_starting_count, proc_starting_cnt,
                          num_procs, SearchComm );
      for( my_start=global_elem_offset+1, j=0 ; j<my_proc ; ++j ) {
        my_start += proc_starting_cnt[j];
      }
      for( j=0 ; j<num_ffi_elems ; ++j ) {
        elem_num_map[local_elem_offset++] = global_elem_offset +  my_start++;
      }
      global_elem_offset += num_ffi_elems_global;
    }
    
    // Construct a unique numbering for the FCI elements (if needed)
    if (num_fci>0) { 
      std::memset( proc_starting_count, 0, num_procs*sizeof(int) );
      proc_starting_count[my_proc] = num_fci_elems;
      contact_global_sum( proc_starting_count, proc_starting_cnt,
                          num_procs, SearchComm );
      for( my_start=global_elem_offset+1, j=0 ; j<my_proc ; ++j ) {
        my_start += proc_starting_cnt[j];
      }
      for( j=0 ; j<num_fci_elems ; ++j ) {
        elem_num_map[local_elem_offset++] = global_elem_offset +  my_start++;
      }
    }
   
    ex_put_elem_num_map( Exodus_ID, elem_num_map );
    delete [] elem_num_map;

    delete [] proc_starting_count;
    delete [] proc_starting_cnt;  
  }
#endif

  //============================================================================
  //  C R E A T E   I N F O   R E C O R D S
  //============================================================================
  // first entry is number of face blocks (number of element blocks
  // without the edge blocks), remainder is the SearchData
  int num_search_entities = search_data->Num_Search_Entities();
  int size_search_data = ContactSearch::NSIZSD*num_search_entities*
        num_search_entities;
  int info_size = size_search_data + 1;
  char** info_strings = new char*[info_size];
  int info_string_len = 80;
  for( i=0 ; i<info_size ; ++i )
    info_strings[i] = new char[info_string_len+1];
  std::sprintf( info_strings[0],"%d",number_of_face_blocks);
  index = 1;
  int num_entity_keys = search_data->Num_Search_Entities();
  for( k=0 ; k<num_entity_keys ; ++k ){
    for( j=0 ; j<num_entity_keys ; ++j ){
      for( i=0 ; i<ContactSearch::NSIZSD ; ++i ){
        Real sd = search_data->Get_Search_Data(
                  (ContactSearch::Search_Data_Index)i,j,k );
        std::sprintf( info_strings[index],"%21.16e",sd);
        index ++;
      }
    }
  }
  ierr = ex_put_info( Exodus_ID,info_size,info_strings );
  for( i=0 ; i<info_size ; ++i )
    delete [] info_strings[i];
  delete [] info_strings;
  
  //============================================================================
  // Collect up the coordinates (from CURRENT_POSITION)
  //============================================================================
  int n=0;
  Real* coordinates = new Real[3*num_proc_nodes];
  std::memset( coordinates, 0, 3*num_proc_nodes*sizeof(Real) );
  for (n=0; n<number_of_primary_nodes; ++n) {
    Real* position = Nodes[n]->Variable(CURRENT_POSITION);
    for( int m=0 ; m<dimensionality ; ++m ) {
      coordinates[n+m*num_proc_nodes] = position[m];
    }
  }
  if (num_ffi>0) {
    for( i=0 ; i<number_of_face_blocks ; ++i ){
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
        ContactInteractionDLL<Real>* interactions = BlockFaces[j]->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        interactions->IteratorStart();
        while (interaction=interactions->IteratorForward()){
          cffi = static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
          ContactFaceFaceVertex<Real>* vertices = cffi->Get_Vertices();
          ContactFace<Real>* slave_face = cffi->SlaveFace();
          for (k=0; k<cffi->NumEdges()+1; ++k) {
            Real local_coords[3];
            Real position[3];
            local_coords[0] = vertices[k].slave_x;
            local_coords[1] = vertices[k].slave_y;
            local_coords[2] = 1.0-local_coords[0]-local_coords[1];
            slave_face->Compute_Global_Coordinates(CURRENT_POSITION,
					           local_coords, position );
            for( int m=0 ; m<dimensionality ; ++m ) {
              coordinates[n+m*num_proc_nodes] = position[m];
            }
            n++;
          }
        }
      }
    }
  }
  if (num_fci>0) {
    for( i=0 ; i<number_of_face_blocks ; ++i ){
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
        ContactInteractionDLL<Real>* interactions = BlockFaces[j]->Get_FaceCoverage_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while (interaction=interactions->IteratorForward()){
            cfci = static_cast<ContactFaceCoverageInteraction*>(interaction);
            ContactFaceCoverageVertex* ll_node;
            for( ll_node=cfci->Head(); ll_node; ll_node=ll_node->next ){
              Real local_coords[3];
              Real position[3];
              local_coords[0] = ll_node->slave_x;
              local_coords[1] = ll_node->slave_y;
              local_coords[2] = 1.0-local_coords[0]-local_coords[1];
              BlockFaces[j]->Compute_Global_Coordinates(CURRENT_POSITION,
		  	          		        local_coords, 
                                                        position );
              for( int m=0 ; m<dimensionality ; ++m ) {
                coordinates[n+m*num_proc_nodes] = position[m];
              }
              n++;
	    }
          }
        }
      }
    }
  }
  
  //============================================================================
  //  O U T P U T   T H E   C O O R D I N A T E S
  //============================================================================
  char* coord_names[3];
  coord_names[0] = (char*) "X";
  coord_names[1] = (char*) "Y";
  coord_names[2] = (char*) "Z";
  ierr = ex_put_coord_names( Exodus_ID, coord_names );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_coord_names",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  
  // Output the coordinates
  ierr = ex_put_coord( Exodus_ID, coordinates, coordinates+num_proc_nodes,
                       coordinates+2*num_proc_nodes );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_coord",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }

  delete [] coordinates;
  
  //============================================================================
  //  O U T P U T   T H E   E L E M E N T   B L O C K S
  //============================================================================
  int max_conn_size=1;
  for( i=0 ; i<number_of_face_blocks ; ++i ){
    max_conn_size = std::max(max_conn_size,face_blocks[i]->Number_of_Faces());
  }
  for( i=0 ; i<number_of_element_blocks ; ++i ){
    max_conn_size = std::max(max_conn_size,element_blocks[i]->Number_of_Elements());
  }
  for( i=0 ; i<number_of_edge_blocks ; ++i ){
    int nedges_in_block = 0;
    ContactEdge<Real>** BlockEdges = 
      reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
    for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
      if( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED )
        nedges_in_block++;
    }
    max_conn_size = std::max(max_conn_size,nedges_in_block);
  }
  max_conn_size = std::max(max_conn_size,num_ffi_elems);
  max_conn_size = std::max(max_conn_size,num_fci_elems);
  char *elem_name=NULL;
  int* connectivity = new int[8*max_conn_size];
  int material_id = 1;
  int number_of_nodes_per_element = -1;
  int attributes_per_element = 0;
  int number_of_elements_in_block;
  for( i=0 ; i<number_of_element_blocks ; ++i ){
    number_of_elements_in_block = element_blocks[i]->Number_of_Elements();
    switch( element_blocks[i]->Type() ){
    case( ContactSearch::CARTESIANHEXELEMENTL8 ):
      elem_name = (char*) "HEX8CAR";
      number_of_nodes_per_element = 8;
      break;
    case( ContactSearch::HEXELEMENTL8 ):
      elem_name = (char*) "HEX8";
      number_of_nodes_per_element = 8;
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    POSTCONDITION(number_of_nodes_per_element>=0);
    POSTCONDITION(elem_name!=NULL);
    ierr = ex_put_elem_block( Exodus_ID, material_id, elem_name,
                  number_of_elements_in_block, 
                  number_of_nodes_per_element,
                  attributes_per_element );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    } 
    if (number_of_elements_in_block>0) {
      index = 0;
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
      for (j=0; j<primary_elem_list->BlockNumEntities(i); ++j) {
        for( k=0 ; k<number_of_nodes_per_element ; ++k )
          //connectivity[index++] = BlockElements[j]->Node(k)->ProcArrayIndex()+1;
          connectivity[index++] = BlockElements[j]->Node(k)->fcs_index+1;
      }
      ierr = ex_put_elem_conn( Exodus_ID, material_id, connectivity );
      if( ierr<0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
    }
    material_id++;
  }
    
  for( i=0 ; i<number_of_face_blocks ; ++i ){
    number_of_elements_in_block = face_blocks[i]->Number_of_Faces();
    number_of_nodes_per_element = ContactSearch::Number_Nodes_Per_Face(face_blocks[i]->Type());
    switch( face_blocks[i]->Type() ){
    case( ContactSearch::LINEFACEL2 ):
      elem_name = (char*) "SHELL";   // "SHELL2";
      break;
    case( ContactSearch::QUADFACEL4 ):
      elem_name = (char*) "SHELL";   // "SHELL4";
      break;
    case( ContactSearch::QUADFACEQ8 ):
      elem_name = (char*) "SHELL";   // "SHELL8";
      break;
    case( ContactSearch::QUADFACEQ9 ):
      elem_name = (char*) "SHELL";   // "SHELL9";
      break;
    case( ContactSearch::TRIFACEL3 ):
      elem_name = (char*) "TRI3";    // "SHELL3";
      break;
    case( ContactSearch::TRIFACEQ6 ):
      elem_name = (char*) "TRI6";    // "SHELL6";
      break;
    case( ContactSearch::SHELLQUADFACEL4 ):
      elem_name = (char*) "SHELL";
      break;
    case( ContactSearch::SHELLTRIFACEL3 ):
      elem_name = (char*) "TRI3";
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    ierr = ex_put_elem_block( Exodus_ID, material_id, elem_name,
                  number_of_elements_in_block, 
                  number_of_nodes_per_element,
                  attributes_per_element );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    } 

    if (number_of_elements_in_block>0) {
      index = 0;
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
        for( k=0 ; k<number_of_nodes_per_element ; ++k )
          //connectivity[index++] = BlockFaces[j]->Node(k)->ProcArrayIndex()+1;
          connectivity[index++] = BlockFaces[j]->Node(k)->fcs_index+1;
      }
      ierr = ex_put_elem_conn( Exodus_ID, material_id, connectivity );
      if( ierr<0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
    }
    material_id++;
  }
  if( dimensionality == 3 ){
    for( i=0 ; i<number_of_edge_blocks ; ++i ){
      // Create the barsets for the edges
      int num_edges_in_block = 0;
#ifndef CONTACT_NO_MPI
      if ( PARALLEL ) {
        ContactEdge<Real>** BlockEdges = 
          reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
        for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
          if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) 
            num_edges_in_block++;
        }
      }
      else
        num_edges_in_block = edge_blocks[i]->Number_of_Edges();
#else
      num_edges_in_block = edge_blocks[i]->Number_of_Edges();
#endif
      int number_of_nodes_per_edge = -1;
      switch (edge_blocks[i]->Type()) {
      case ContactSearch::LINEEDGEL2:
        elem_name = (char*) "BAR";
        number_of_nodes_per_edge = 2;
        break;
      case ContactSearch::LINEEDGEQ3:
        elem_name = (char*) "BAR3";
        number_of_nodes_per_edge = 3;
        break;
      default:
        POSTCONDITION(0);
        break;
      }
      POSTCONDITION(number_of_nodes_per_edge>=0);
      ierr = ex_put_elem_block( Exodus_ID, material_id, elem_name,
                                num_edges_in_block, number_of_nodes_per_edge,
                                0 );
      if( ierr!=0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
      if (num_edges_in_block>0) {
        index = 0;
#ifndef CONTACT_NO_MPI
        if ( PARALLEL ) {
          ContactEdge<Real>** BlockEdges = 
            reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
          for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
            if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
              for( k=0 ; k<number_of_nodes_per_edge ; ++k ){
                //connectivity[index++] = BlockEdges[j]->Node(k)->ProcArrayIndex()+1;
                connectivity[index++] = BlockEdges[j]->Node(k)->fcs_index+1;
              }
            }
          }
        } else {
#endif
          ContactEdge<Real>** BlockEdges = 
            reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
          for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
            for( k=0 ; k<number_of_nodes_per_edge ; ++k ){
              //connectivity[index++] = BlockEdges[j]->Node(k)->ProcArrayIndex()+1;
              connectivity[index++] = BlockEdges[j]->Node(k)->fcs_index+1;
            }
          }
#ifndef CONTACT_NO_MPI
        }
#endif
        ierr = ex_put_elem_conn( Exodus_ID, material_id, connectivity );
        if( ierr!=0 ){
          std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
          error_handler->Add_Error_Message(message);
          return ContactSearch::EXODUS_ERROR;
        }
      }
      material_id++;
    }
  }
  if( number_of_node_blocks > 1 ){
    elem_name = (char*) "SPH";
    for( i=1 ; i<number_of_node_blocks ; ++i){
      if( node_blocks[i]->Type() == ContactSearch::POINT ){
	int num_sph = node_blocks[i]->Number_of_Nodes();
	int eb_id = number_of_face_blocks + number_of_edge_blocks + i;
	ierr = ex_put_elem_block( Exodus_ID, eb_id, elem_name,
				  num_sph, 1, 0 );
	if( ierr!=0 ) {
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
	if( num_sph ){
          index = 0;
          ContactNode<Real>** BlockNodes = 
            reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
          for (j=0; j<primary_node_list->BlockNumEntities(i); ++j) {
	    //connectivity[index++] = BlockNodes[j]->ProcArrayIndex()+1;
	    connectivity[index++] = BlockNodes[j]->fcs_index+1;
          }
	  ierr = ex_put_elem_conn( Exodus_ID, eb_id, connectivity );
	}
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
  }
  if( Need_Sphere_Element_Block ){
    elem_name = (char*) "SPHERE";
    int num_spheres = 0;
    if( Dont_Have_Nodes_and_Elements ) num_spheres = 1;
    ierr = ex_put_elem_block( Exodus_ID, 50000, elem_name,
                              num_spheres, 1, 0 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if( num_spheres ){
      connectivity[0] = num_proc_nodes;
      ierr = ex_put_elem_conn( Exodus_ID, 50000, connectivity );
      if( ierr!=0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
    }
  }
  if( num_ffi ){
    n = 0;
    elem_name = (char*) "TRI3";    // "SHELL3";
    ierr = ex_put_elem_block( Exodus_ID, 60000, elem_name,
                              num_ffi_elems, 3, 0 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if (num_ffi_elems>0) {
      int nstart = num_base_nodes+1;
      for( i=0 ; i<number_of_face_blocks ; ++i){
        ContactFace<Real>** BlockFaces = 
          reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
        for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
          ContactInteractionDLL<Real>* interactions = 
            BlockFaces[j]->Get_FaceFace_Interactions();
          if(interactions == NULL) continue;
          interactions->IteratorStart();
          while (interaction=interactions->IteratorForward()){
            cffi = static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
            int k1, k2;
            int k3  = cffi->NumEdges();
            int num = cffi->NumEdges();
            for (k1=0, k2=1; k1<num; k1++, k2=(k1+1)%num) {
              connectivity[n++] = nstart+k1;
              connectivity[n++] = nstart+k2;
              connectivity[n++] = nstart+k3;
            }
            nstart += cffi->NumEdges()+1;
          }
        }
      }
      ierr = ex_put_elem_conn( Exodus_ID, 60000, connectivity );
      if( ierr!=0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
      num_base_nodes += num_ffi_nodes;
    }
  }
  if( num_fci ){
    n = 0;
    elem_name = (char*) "BAR";
    ierr = ex_put_elem_block( Exodus_ID, 70000, elem_name,
                              num_fci_elems, 2, 0 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_elem_block",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if (num_fci_elems>0) {
      int nstart = num_base_nodes+1;
      for( i=0 ; i<number_of_face_blocks ; ++i){
        ContactFace<Real>** BlockFaces = 
          reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
        for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
          ContactInteractionDLL<Real>* interactions = BlockFaces[j]->Get_FaceCoverage_Interactions();
	  if(interactions != NULL) {
            interactions->IteratorStart();
            while (interaction=interactions->IteratorForward()){
              cfci = static_cast<ContactFaceCoverageInteraction*>(interaction);
              int k1, k2;
              int num = cfci->NumVertices();
              for (k1=0, k2=1; k1<num; k1++, k2=(k1+1)%num) {
                connectivity[n++] = nstart+k1;
                connectivity[n++] = nstart+k2;
              }
              nstart += num;
            }
	  }
        }
      }
      ierr = ex_put_elem_conn( Exodus_ID, 70000, connectivity );
      if( ierr!=0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_elem_conn",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
    }
  }
  delete [] connectivity;

  return Exodus_Output_Results( Exodus_ID, Time, errors,
				enforcement,
				num_proc_nodes,
				number_of_elem_blocks, 
				number_of_elems,
				number_edges_owned,
				Need_Sphere_Element_Block,
				num_nni, max_nni, 
                                num_nfi, max_nfi,
				num_nsi, max_nsi,
				num_ffi, max_ffi, max_ffi_verts,
				num_fci, max_fci, max_fci_verts,
				num_eei, max_eei,
				normal_smoothing,
				multiple_interaction,
				sharp_smooth_ang,
				normal_smoothing_distance,
				smoothing_resolution,
				compute_node_areas );
}

      
#else

#include "ContactTopology.h"
#include "ContactErrors.h"
#include <cstdio>

ContactSearch::ContactErrorCode
ContactTopology::Exodus_Output( int , Real , 
                                ContactErrors* error_handler,
                                ContactSearchData* ,
                                const ContactSearch::Search_Option_Status&,
                                const ContactSearch::Search_Option_Status&,
                                const Real& , 
                                const Real& , 
                                const ContactSearch::Smoothing_Resolution &,
                              const ContactSearch::Search_Option_Status&,
                                ContactEnforcement*  )
{
  std::sprintf(message,"Library wasn't compiled with Exodus support");
  error_handler->Add_Error_Message(message);
  return ContactSearch::EXODUS_ERROR;
}

ContactSearch::ContactErrorCode
ContactTopology::ExodusMesh_Output( int Exodus_ID, 
                                ContactErrors* error_handler,
                                ContactSearchData* search_data)
{
  std::sprintf(message,"Library wasn't compiled with Exodus support");
  error_handler->Add_Error_Message(message);
  return ContactSearch::EXODUS_ERROR;
}

#endif
        
