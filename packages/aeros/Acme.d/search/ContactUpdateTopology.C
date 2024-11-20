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


#include "contact_assert.h"
#include "contact_sorting.h"
#include "ContactAnalyticSurface.h"
#include "ContactSymComm.h"
#include "ContactEdgeBlock.h"
#include "ContactTopologyEntityHash.h"
#include "ContactEntityDataHash.h"
#include "ContactErrors.h"
#include "ContactFaceBlock.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineFaceL2.h"
#include "ContactNode.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactQuadFaceL4.h"
#include "ContactTopology.h"
#include "ContactTriFaceL3.h"
#include "CString.h"
#include "search_methods.h"
#include "Contact_Communication.h"
#include "ContactCommBuffer.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSequentialAllocator.h"
#include "ContactShellHandler.h"

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef CONTACT_NO_MPI
#include "zoltan.h"
#include "ContactZoltan.h"
#include "ContactZoltanComm.h"
#endif

void   
ContactTopology::UpdateTopology( ContactErrors* Errors,

				 const int* number_node_deaths_per_block,
				 const int* global_node_death_ids,
				 const int* number_face_deaths_per_block,
				 const int* global_face_death_ids,
				 const int* number_elem_deaths_per_block,
				 const int* global_elem_death_ids,
                                 
				 const int* number_node_births_per_block,
				 const int* global_node_birth_ids,
				 const int* exodus_node_birth_ids,
				 const int* number_face_births_per_block,
				 const int* global_face_birth_ids,
				 const int* face_connectivity,
				 const int* number_elem_births_per_block,
				 const int* global_elem_birth_ids,
				 const int* elem_connectivity,
       
		       	 	 const int  num_node_exports,
                       	 	 const int* node_export_id_list,
                                 const int* node_export_pids,
		       	 	 const int  num_face_exports,
                       	 	 const int* face_export_id_list,
                                 const int* face_export_pids,
		       	 	 const int  num_element_exports,
                       	 	 const int* element_export_id_list,
                                 const int* element_export_pids,
                                 
                                 const int num_nodes,
                                 const int* node_host_ids,
                                 const int num_faces,
                                 const int* face_host_ids,
                                 const int num_elements,
                                 const int* element_host_ids,
                                 
				 const int number_comm_partners, 
				 const int* comm_proc_ids,
				 const int* number_nodes_to_partner, 
				 const int* comm_nodes,
                                 
                                 ContactParOStream& postream,
				 ContactSearch::ContactErrorCode& error_code )
{
  int i,j,k,n;
  ContactTopologyEntity<Real>* entity=NULL;
#ifndef CONTACT_NO_MPI  
  if( Edge_SymComm ) delete Edge_SymComm;
  Edge_SymComm = NULL;
#endif
  if( edge_list ) delete edge_list;
  edge_list = new ContactTopologyEntityList();
  number_of_edges = 0;
  if( number_of_edge_blocks && edge_blocks ){
    for( i=0 ; i<number_of_edge_blocks ; ++i) {
      if( edge_blocks[i] ) delete edge_blocks[i];
    }
    delete [] edge_blocks;
    edge_blocks = NULL;
    number_of_edge_blocks = 0;
  }
  for( i=0 ; i<number_of_node_blocks ; ++i){
    ContactBlockEntityList* nodes = node_blocks[i]->NodeList();
    nodes->IteratorStart();
    while ((entity=nodes->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      node->Delete_Face_Connections();
    }
  }
  
  if( scratch_set ) search->Get_Scratch_Allocator().Purge_Mem();

#ifdef CONTACT_DEBUG_NODE
  if( number_debug_nodes ){
    for( i=0 ; i<number_debug_nodes ; ++i)
      delete debug_node_global_ids[i];
    delete [] debug_nodes;
    delete [] debug_node_exodus_ids;
    delete [] debug_node_global_ids;
  }
  number_debug_nodes    = 0;
  debug_node_global_ids = NULL;
  debug_node_exodus_ids = NULL;
  debug_nodes           = NULL;
#endif

  error_code = ContactSearch::NO_ERROR;
  errors = Errors;
  number_of_comm_partners = number_comm_partners;
  two_configurations = false;
  scratch_set = false;
  
  ContactType* node_birth_entity_types = NULL;
  if (shell_handler) {
  
  } else {
    int nnodes = 0;
    for( i=0 ; i<number_of_node_blocks ; ++i){
      nnodes += number_node_births_per_block[i];
    }
    node_birth_entity_types = new ContactType[nnodes];
    for( i=0 ; i<nnodes ; ++i) {
      node_birth_entity_types[i] = CT_NODE;
    }
  }
                
  TopologyBirth( Errors,
		 number_node_births_per_block,
		 global_node_birth_ids,
		 exodus_node_birth_ids,
                 node_birth_entity_types,
		 number_face_births_per_block,
		 global_face_birth_ids,
		 face_connectivity,
		 number_elem_births_per_block,
		 global_elem_birth_ids,
		 elem_connectivity,
                 postream,
		 error_code ); 

  TopologyDeath( Errors,
		 number_node_deaths_per_block,
		 global_node_death_ids,
		 number_face_deaths_per_block,
		 global_face_death_ids,
		 number_elem_deaths_per_block,
		 global_elem_death_ids,
                 postream,
		 error_code );
                 
  TopologyDLB( Errors,
               num_node_exports,
               node_export_id_list,
               node_export_pids,
               num_face_exports,
               face_export_id_list,
               face_export_pids,
               num_element_exports,
               element_export_id_list,
               element_export_pids,
               postream,
               error_code ); 
               
  if (node_birth_entity_types) delete [] node_birth_entity_types; 
  
  for( n=0, i=0 ; i<number_of_node_blocks ; ++i){
    n += node_blocks[i]->NodeList()->NumEntities();
  }
  POSTCONDITION(num_nodes==n);
  for( n=0, i=0 ; i<number_of_face_blocks ; ++i){
    n += face_blocks[i]->FaceList()->NumEntities();
  }  
  POSTCONDITION(num_faces==n);
  for( n=0, i=0 ; i<number_of_element_blocks ; ++i){
    n += element_blocks[i]->ElemList()->NumEntities();
  }  
  POSTCONDITION(num_elements==n);
     
  for( n=0, i=0 ; i<number_of_node_blocks ; ++i){
    for (k=0, j=0; j<node_blocks[i]->NodeList()->NumEntities(); ++j) {
      ContactHostGlobalID global_id( node_host_ids[2*n+0], 
                                     node_host_ids[2*n+1] );
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                          (node_blocks[i]->NodeList()->Find(global_id));
      PRECONDITION( node );
      node->HostArrayIndex(k);
      node->HostGlobalArrayIndex(n);
      n++;
      k++;
    }
    node_blocks[i]->NodeList()->Rehash();
  } 
  node_list->CleanUp();
  node_list->BuildList(node_blocks, number_of_node_blocks,
                       no_parallel_consistency==ContactSearch::INACTIVE);
  number_of_nodes = node_list->NumEntities();
  POSTCONDITION(num_nodes==number_of_nodes);


  for( n=0, i=0 ; i<number_of_face_blocks ; ++i){
    for (k=0, j=0; j<face_blocks[i]->FaceList()->NumEntities(); ++j) {
      ContactHostGlobalID global_id( face_host_ids[2*n+0], 
                                     face_host_ids[2*n+1] );
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                          (face_blocks[i]->FaceList()->Find(global_id));
      PRECONDITION( face );
      face->HostArrayIndex(k);
      face->HostGlobalArrayIndex(n);
      n++;
      k++;
    }
    face_blocks[i]->FaceList()->Rehash();
  }
  face_list->CleanUp();
  face_list->BuildList(face_blocks, number_of_face_blocks,
                       no_parallel_consistency==ContactSearch::INACTIVE);
  number_of_faces = face_list->NumEntities();
  POSTCONDITION(num_faces==number_of_faces);


  for( n=0, i=0 ; i<number_of_element_blocks ; ++i){
    for (k=0, j=0; j<element_blocks[i]->ElemList()->NumEntities(); ++j) {
      ContactHostGlobalID global_id( element_host_ids[2*n+0], 
                                     element_host_ids[2*n+1] );
      ContactElement* element = static_cast<ContactElement*>
                                (element_blocks[i]->ElemList()->Find(global_id));
      PRECONDITION( element );
      element->HostArrayIndex(k);
      element->HostGlobalArrayIndex(n);
      n++;
      k++;
    }
    element_blocks[i]->ElemList()->Rehash();
  }
  elem_list->CleanUp();
  elem_list->BuildList(element_blocks, number_of_element_blocks,
                       no_parallel_consistency==ContactSearch::INACTIVE);
  number_of_elements = elem_list->NumEntities();
  POSTCONDITION(num_elements==number_of_elements);

  Connect_Faces_to_Nodes();
  // Construct Edges and create connections
  if( dimensionality == 3 ) Construct_and_Connect_Edges(error_code);
  for( i=0 ; i<number_of_edge_blocks ; ++i){
    edge_blocks[i]->EdgeList()->Rehash();
  }
  number_of_edges = edge_list->NumEntities();
         
  // Check to see if the mesh is valid.  If not return
  error_code = (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, SearchComm );
  if( error_code ) return;
  
  UpdateInteractions( Errors, postream, error_code );  
  
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->OwnerProcArrayIndex(Nodes[i]->ProcArrayIndex());
    Nodes[i]->PrimaryProcArrayIndex(Nodes[i]->ProcArrayIndex());
  }
  ContactEdge<Real>** Edges = reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
  for (i=0; i<number_of_edges; ++i) {
    Edges[i]->OwnerProcArrayIndex(Edges[i]->ProcArrayIndex());
    Edges[i]->PrimaryProcArrayIndex(Edges[i]->ProcArrayIndex());
  }
  ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->OwnerProcArrayIndex(Faces[i]->ProcArrayIndex());
    Faces[i]->PrimaryProcArrayIndex(Faces[i]->ProcArrayIndex());
  }
  ContactElement** Elements = reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (i=0; i<number_of_elements; ++i) {
    Elements[i]->OwnerProcArrayIndex(Elements[i]->ProcArrayIndex());
    Elements[i]->PrimaryProcArrayIndex(Elements[i]->ProcArrayIndex());
  }
  
#ifndef CONTACT_NO_MPI  
  if( Node_SymComm ) delete Node_SymComm;
  Node_SymComm = NULL;
#endif

#ifndef CONTACT_NO_MPI
  // Create the communication lists (they must exist even if we don't
  // communicate to any other processors)
  int num_node_comm_ent = 0;
  for( i=0 ; i<number_comm_partners ; ++i) {
    num_node_comm_ent += number_nodes_to_partner[i];
  }
  ContactTopologyEntity<Real>** node_ent_comm_list = NULL;
  if( num_node_comm_ent ){
    node_ent_comm_list = new ContactTopologyEntity<Real>*[num_node_comm_ent];
    for( i=0 ; i<num_node_comm_ent ; ++i) {
      int index = comm_nodes[i]-1;
      ContactHostGlobalID global_id( node_host_ids[2*index+0], 
                                     node_host_ids[2*index+1] );
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(node_list->Find(global_id));
      POSTCONDITION(node);
      node_ent_comm_list[i] = node;
    }
  }
  
  Node_SymComm = new ContactSymComm( number_comm_partners,
				     comm_proc_ids,
				     number_nodes_to_partner,
				     node_ent_comm_list );
  delete [] node_ent_comm_list;

#endif


  // now compute the owners for all nodes
  Compute_Owners(error_code);
  
  for( i=0 ; i<number_of_edge_blocks ; ++i){
    edge_blocks[i]->EdgeList()->Rehash();
  }
  
  //Compute_Edge_Comm_List(error_code);
  
#ifndef CONTACT_NO_MPI
  contact_swap_edge_faces(SearchComm, *Edge_SymComm, *comm_buffer);
#endif

  Faces = reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->SetNeighborFacesInfo();
  }

  // Check to see if the mesh is valid across processors.  If not return
  error_code = (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, SearchComm );
  if( error_code ) return;

  // Initialize the scratch memory manager
  scratch_set = false;
#if CONTACT_DEBUG_PRINT_LEVEL>=5
  Display_Entities(postream, 1);
  postream.flush();
#endif

}

void   
ContactTopology::TopologyDeath( ContactErrors* Errors,

				const int* number_node_deaths_per_block,
				const int* global_node_death_ids,
				const int* number_face_deaths_per_block,
				const int* global_face_death_ids,
				const int* number_elem_deaths_per_block,
				const int* global_elem_death_ids,
                                 
                                ContactParOStream& postream,
				ContactSearch::ContactErrorCode& error_code )
{ 
  
  int i,j,k,n;
  
  // 1) tag all the entites for deletion
  // 2) look thru all the connections and delete the tagged entries
  // 3) delete all tagged entities
  
  //===========================================================================
  // initalize all the entities' temp_tag:  
  //    =0, entity is inactive  
  //    =1, entity is active
  //===========================================================================
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(node_list->EntityList());
  for( i=0 ; i<number_of_nodes ; ++i){
    Nodes[i]->temp_tag = 1;
  }
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for( i=0 ; i<number_of_faces ; ++i){
    Faces[i]->temp_tag = 1;
  }
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for( i=0 ; i<number_of_elements ; ++i){
    Elements[i]->temp_tag = 1;
  }

  //===========================================================================
  // tag all the entites specified for deletion and any connected entities
  //===========================================================================
  for( n=0, i=0 ; i<number_of_node_blocks ; ++i){
    for( j=0 ; j<number_node_deaths_per_block[i] ; ++j){
      ContactHostGlobalID gid(global_node_death_ids[2*n+0],
                              global_node_death_ids[2*n+1]);
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(node_blocks[i]->NodeList()->Find(gid));
      node->temp_tag = 0;
      n++;
    }
  }
  for( n=0, i=0 ; i<number_of_face_blocks ; ++i){
    for( j=0 ; j<number_face_deaths_per_block[i] ; ++j){
      ContactHostGlobalID gid(global_face_death_ids[2*n+0],
                              global_face_death_ids[2*n+1]);
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(face_blocks[i]->FaceList()->Find(gid));
      face->temp_tag = 0;
      // tag all the connected nodes for deletion
      for (k=0; k<face->Nodes_Per_Face(); ++k) {
        face->Node(k)->temp_tag = 0;
      }
      n++;
    }
  }
  for( n=0, i=0 ; i<number_of_element_blocks ; ++i){
    for( j=0 ; j<number_elem_deaths_per_block[i] ; ++j){
      ContactHostGlobalID gid(global_elem_death_ids[2*n+0],
                              global_elem_death_ids[2*n+1]);
      ContactElement* element = static_cast<ContactElement*>(element_blocks[i]->ElemList()->Find(gid));
      element->temp_tag = 0;
      // tag all the connected nodes and edges for deletion
      for (k=0; k<element->Nodes_Per_Element(); ++k) {
        element->Node(k)->temp_tag = 0;
      }
      n++;
    }
  }
  
  //===========================================================================
  // Check all connections to make sure shared entities are marked correctly
  //===========================================================================
  for( i=0 ; i<number_of_faces ; ++i){
    if (Faces[i]->temp_tag==1) {
      // for active faces, tag all the connected nodes as being active
      for (k=0; k<Faces[i]->Nodes_Per_Face(); ++k) {
        Faces[i]->Node(k)->temp_tag = 1;
      }
    }
  }
  for( i=0 ; i<number_of_elements ; ++i){
    if (!Elements[i]->temp_tag) {
      // for active elements, tag all the connected nodes as being active
      for (k=0; k<Elements[i]->Nodes_Per_Element(); ++k) {
        Elements[i]->Node(k)->temp_tag = 1;
      }
    }
  }

  //===========================================================================
  // Now actually delete the tagged entities
  //===========================================================================
  for( i=0 ; i<number_of_nodes ; ++i){
    if (!Nodes[i]->temp_tag) {
      node_blocks[Nodes[i]->BlockID()]->Delete_Node(Nodes[i]);
    }
  }
  for( i=0 ; i<number_of_faces ; ++i){
    if (!Faces[i]->temp_tag) {
      face_blocks[Faces[i]->BlockID()]->Delete_Face(Faces[i]);
    }
  }
  
  for( i=0 ; i<number_of_elements ; ++i){
    if (!Elements[i]->temp_tag) {
      element_blocks[Elements[i]->BlockID()]->Delete_Element(Elements[i]);
    }
  }
}

void   
ContactTopology::TopologyBirth( ContactErrors* Errors,
                                 
				const int* number_node_births_per_block,
				const int* global_node_birth_ids,
				const int* exodus_node_birth_ids,
                                ContactType* node_birth_entity_types,
				const int* number_face_births_per_block,
				const int* global_face_birth_ids,
				const int* face_connectivity,
				const int* number_elem_births_per_block,
				const int* global_elem_birth_ids,
				const int* elem_connectivity,
                                
                                ContactParOStream& postream,
				ContactSearch::ContactErrorCode& error_code )
{
  int i, j, k;
  int num_node_births = 0;
  int num_face_births = 0;
  int num_elem_births = 0;
  ContactTopologyEntity<Real>* entity;
  
  for( k=0, i=0 ; i<number_of_node_blocks ; ++i){
    num_node_births += number_node_births_per_block[i];
  }
  for( i=0 ; i<number_of_face_blocks ; ++i){
    num_face_births += number_face_births_per_block[i];
  }
  for( i=0 ; i<number_of_element_blocks ; ++i){
    num_elem_births += number_elem_births_per_block[i];
  }
  if (num_node_births+num_face_births+num_elem_births==0) return;
  
  //===========================================================================
  // initalize all the existing entities' temp_tag so:  
  //    =0, entity is new  
  //    =1, entity is pre-existing
  //===========================================================================
  for( i=0 ; i<number_of_node_blocks ; ++i){
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
      entity->temp_tag = 1;
    }
  }
  for( i=0 ; i<number_of_face_blocks ; ++i){
    face_blocks[i]->FaceList()->IteratorStart();
    while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
      entity->temp_tag = 1;
    }
  }
  for( i=0 ; i<number_of_element_blocks ; ++i){
    element_blocks[i]->ElemList()->IteratorStart();
    while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
      entity->temp_tag = 1;
    }
  }

  // update the topology's nodes due to birth
  if (num_node_births) {
    int* gids   = (int*)global_node_birth_ids;
    int* exoids = (int*)exodus_node_birth_ids;
    ContactType* entity_types = node_birth_entity_types;
    for( i=0 ; i<number_of_node_blocks ; ++i){
      int nnodes = number_node_births_per_block[i];
      if (nnodes) {
        node_blocks[i]->Add_Nodes(nnodes, gids, exoids, entity_types);
        gids         += 2*nnodes;
        exoids       +=   nnodes;
        entity_types +=   nnodes;
      }
    }
  }
  // update the topology's faces due to birth
  if (num_face_births) {
    int* gids = (int*)global_face_birth_ids;
    for( i=0 ; i<number_of_face_blocks ; ++i){
      int nfaces = number_face_births_per_block[i];
      if (nfaces) {
        face_blocks[i]->Add_Faces(nfaces, gids);
        gids += 2*nfaces;
      }
    }
  }
  // update the topology's elements due to birth
  if (num_elem_births) {
    int* gids = (int*)global_elem_birth_ids;
    for( i=0 ; i<number_of_element_blocks ; ++i){
      int nelems = number_elem_births_per_block[i];
      if (nelems) {
        element_blocks[i]->Add_Elements(nelems, gids);
        gids += 2*nelems;
      }
    }
  }
 
  // Set the entity key for new nodes in blocks 2-number_of_node_blocks
  if (num_node_births) {
    for( i=1 ; i<number_of_node_blocks ; ++i){
      int entity_key = number_of_element_blocks + number_of_face_blocks + (i-1);
      node_blocks[i]->NodeList()->IteratorStart();
      while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
        if (entity->temp_tag == 0) {
          ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
          node->Entity_Key( entity_key );
        }
      }
    }
  }

  // Connect the faces to their nodes (from the connectivity array)
  if (num_face_births) {
    int index = 0;
    for( i=0 ; i<number_of_face_blocks ; ++i){
      face_blocks[i]->FaceList()->IteratorStart();
      while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
        if (entity->temp_tag == 0) {
          ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
          for( k=0 ; k<face->Nodes_Per_Face() ; ++k){
            ContactHostGlobalID global_id( face_connectivity[2*index+0], 
                                           face_connectivity[2*index+1] );
            ContactNode<Real>* node = NULL;
            for( j=0 ; j<number_of_node_blocks ; ++j){
              node = static_cast<ContactNode<Real>*>
                     (node_blocks[i]->NodeList()->Find(global_id));
              if (node) break;
            }
            PRECONDITION( node );
            face->ConnectNode(k, node);
            index++;
          }
        }
      }
    }
  }

  // Connect the elements to their nodes (from the connectivity array)
  if (num_elem_births) {
    int index = 0;
    for( i=0 ; i<number_of_element_blocks ; ++i){
      element_blocks[i]->ElemList()->IteratorStart();
      while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
        if (entity->temp_tag==0) {
          ContactElement* element = static_cast<ContactElement*>(entity);
          for( k=0 ; k<element->Nodes_Per_Element() ; ++k){
            ContactHostGlobalID global_id( elem_connectivity[2*index+0], 
                                           elem_connectivity[2*index+1] );
            ContactNode<Real>* node = NULL;
            for( j=0 ; j<number_of_node_blocks ; ++j){
              node = static_cast<ContactNode<Real>*>
                     (node_blocks[i]->NodeList()->Find(global_id));
              if (node) break;
            }
            PRECONDITION( node );
            element->ConnectNode(k, node);
            index++;
          }
        }
      }
    }
  }                   
} 
  
void   
ContactTopology::TopologyDLB( ContactErrors* Errors,
                                 
		       	      const int  num_node_exports,
                       	      const int* node_export_ids,
                              const int* node_export_pids,
		       	      const int  num_face_exports,
                       	      const int* face_export_ids,
                              const int* face_export_pids,
		       	      const int  num_element_exports,
                       	      const int* element_export_ids,
                              const int* element_export_pids,
                              
                              ContactParOStream& postream,
			      ContactSearch::ContactErrorCode& error_code )
{
#ifndef CONTACT_NO_MPI
  int size = num_node_exports+num_face_exports+num_element_exports;
  int sum  = contact_global_sum(size, SearchComm );
  if (sum==0) return;
  
  ContactTopologyEntity<Real>* entity=NULL;
  
  for( int i=0 ; i<number_of_node_blocks ; ++i){
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
      entity->temp_tag = 0;
    }
  }
  for( int i=0 ; i<number_of_face_blocks ; ++i){
    face_blocks[i]->FaceList()->IteratorStart();
    while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
      entity->temp_tag = 0;
    }
  }
  for( int i=0 ; i<number_of_element_blocks ; ++i){
    element_blocks[i]->ElemList()->IteratorStart();
    while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
      entity->temp_tag = 0;
    }
  }
  
  //====================================================
  // dynamic load balancing step
  //====================================================
  LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
  int        zoltan_pid;
  ContactZoltanComm Zoltan_Comm( ContactZoltanComm::ZOLTAN_EXPORT);
  for(int i=0; i<num_node_exports; ++i) {
    ContactNode<Real>* node = NULL;
    for( int j=0 ; j<number_of_node_blocks ; ++j){
      ContactHostGlobalID global_id( node_export_ids[2*i+0], 
                                     node_export_ids[2*i+1] );
      node = static_cast<ContactNode<Real>*>
             (node_blocks[j]->NodeList()->Find(global_id));
      if ( node ) break;
    }
    POSTCONDITION( node );
    node->ZoltanLID(CT_NODE, zoltan_lid);
    node->ZoltanGID(CT_NODE, zoltan_gid);
    //node->temp_tag = 1;
    zoltan_pid = node_export_pids[i];
    Zoltan_Comm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
  }
  for(int i=0; i<num_face_exports; ++i) {
    ContactFace<Real>* face = NULL;
    for( int j=0 ; j<number_of_face_blocks ; ++j){
      ContactHostGlobalID global_id( face_export_ids[2*i+0], 
                                     face_export_ids[2*i+1] );
      face = static_cast<ContactFace<Real>*>
             (face_blocks[j]->FaceList()->Find(global_id));
      if ( face ) break;
    }
    POSTCONDITION( face ); 
    face->ZoltanLID(CT_FACE, zoltan_lid);
    face->ZoltanGID(CT_FACE, zoltan_gid);
    face->temp_tag = 1;
    zoltan_pid = face_export_pids[i];
    Zoltan_Comm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
  }
  for(int i=0; i<num_element_exports; ++i) {
    ContactElement* element = NULL;
    for( int j=0 ; j<number_of_element_blocks ; ++j){
      ContactHostGlobalID global_id( element_export_ids[2*i+0], 
                                     element_export_ids[2*i+1] );
      element = static_cast<ContactElement*>
                (element_blocks[j]->ElemList()->Find(global_id));
      if ( element ) break;
    }
    POSTCONDITION( element );
    element->ZoltanLID(CT_ELEMENT, zoltan_lid);
    element->ZoltanGID(CT_ELEMENT, zoltan_gid);
    element->temp_tag = 1;
    zoltan_pid = face_export_pids[i];
    Zoltan_Comm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
  }
  
  int       num_import   = -1;
  LB_ID_PTR import_gids  = NULL;
  LB_ID_PTR import_lids  = NULL;
  int*      import_procs = NULL;
  int       num_export   = Zoltan_Comm.Num_Export();
  LB_ID_PTR export_gids  = Zoltan_Comm.Export_GIDS();
  LB_ID_PTR export_lids  = Zoltan_Comm.Export_LIDS();
  int*      export_procs = Zoltan_Comm.Export_Procs();
  search->Get_Zoltan()->Set_DynamicLoadBalanceCallBacks(search);
  search->Get_Zoltan()->Help_Migrate(num_import,  import_gids, 
                                     import_lids, import_procs, 
                                     num_export,  export_gids,
		                     export_lids, export_procs);
                       
  //===========================================================================
  // Now actually delete the exported faces and elements
  //===========================================================================
  node_list->CleanUp();
  node_list->BuildList(node_blocks, number_of_node_blocks,
                       no_parallel_consistency==ContactSearch::INACTIVE);
  for( int i=0 ; i<number_of_face_blocks ; ++i){
    face_blocks[i]->FaceList()->IteratorStart();
    while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = reinterpret_cast<ContactFace<Real> *>(entity);
      int num_nodes = face->Nodes_Per_Face();
      ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
      POSTCONDITION( node_info );
      for( int j=0 ; j<num_nodes ; ++j){
        ContactNode<Real>* node = static_cast<ContactNode<Real> *>(node_list->Find( &node_info[j] ));
        POSTCONDITION( node );
        face->ConnectNode( j,node );
      }
      if (face->temp_tag) {
        face_blocks[i]->Delete_Face(face);
      } else {
        for( int j=0 ; j<num_nodes ; ++j){
          face->Node(j)->temp_tag = 1;
        }
      }
    }
  }
  for( int i=0 ; i<number_of_element_blocks ; ++i){
    element_blocks[i]->ElemList()->IteratorStart();
    while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
      ContactElement* element = reinterpret_cast<ContactElement *>(entity);
      int num_nodes = element->Nodes_Per_Element();
      ContactTopologyEntity<Real>::connection_data *node_info = element->NodeInfo();
      POSTCONDITION( node_info );
      for( int j=0 ; j<num_nodes ; ++j){
        ContactNode<Real>* node = static_cast<ContactNode<Real> *>(node_list->Find( &node_info[j] ));
        POSTCONDITION( node );
        element->ConnectNode( j,node );
      }
      if (element->temp_tag) {
        element_blocks[i]->Delete_Element(element);
      } else {
        for( int j=0 ; j<num_nodes ; ++j){
          element->Node(j)->temp_tag = 1;
        }
      }
    }
  }     
  for( int i=1 ; i<number_of_node_blocks ; ++i){
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) { 
      entity->temp_tag = 1;
    }
  } 
  for(int i=0; i<num_node_exports; ++i) {
    ContactNode<Real>* node = NULL;
    int j;
    for( j=0 ; j<number_of_node_blocks ; ++j){
      ContactHostGlobalID global_id( node_export_ids[2*i+0], 
                                     node_export_ids[2*i+1] );
      node = static_cast<ContactNode<Real>*>
             (node_blocks[j]->NodeList()->Find(global_id));
      if ( node ) break;
    }
    POSTCONDITION( node );
    if (j>0) {
      node->temp_tag = 0;
    }
  }       
  //===========================================================================
  // Now actually delete the hanging nodes
  //===========================================================================
  for( int i=0 ; i<number_of_node_blocks ; ++i){
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
      if (entity->temp_tag==0) {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
        node_blocks[i]->Delete_Node(node);
      }
    }
  }
#endif
}

void   
ContactTopology::UpdateInteractions( ContactErrors* Errors,
                                     ContactParOStream& postream,
				     ContactSearch::ContactErrorCode& error_code )
{ 
  int i, j;
  int my_proc = contact_processor_number(SearchComm);
  int num_query = 0;
  ContactTopologyEntity<Real>* entity=NULL;

  //===========================================================================
  // now check if any master entities for the
  // remaining interactions have been deleted.
  //===========================================================================
  for( i=0 ; i<number_of_node_blocks ; ++i){
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      node->temp_tag = 0;
      ContactNodeEntityInteraction** interactions = 
	node->Get_NodeEntity_Interactions();
      for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
        ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[j]);
        if (cnfi==NULL) continue;
        if (cnfi->FaceEntityData()->owner == my_proc) {
          int blk = cnfi->FaceEntityData()->block_id;
          ContactFace<Real>* face = static_cast<ContactFace<Real>*>
            (face_blocks[blk]->FaceList()->Find(cnfi->FaceEntityData()));
          if (!face) {
            // master face not present, delete interaction
            node->Delete_NodeEntity_Interaction(cnfi);
            node->temp_tag = 1;
          }
        } else {
          num_query++;
        }
      }
    }
  }
  for( i=0 ; i<number_of_face_blocks ; ++i){
    face_blocks[i]->FaceList()->IteratorStart();
    while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      face->temp_tag = 0;
      ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
      if(interactions == NULL) continue;
      ContactInteractionEntity<Real>* interaction;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactFaceFaceInteraction<Real>* cffi =
          static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
        if (cffi->MasterFaceEntityData()->owner == my_proc) {
          int blk = cffi->MasterFaceEntityData()->block_id;
          face = static_cast<ContactFace<Real>*>
            (face_blocks[blk]->FaceList()->Find(cffi->MasterFaceEntityData()));
          if (!face) {
            // master face not present, delete interaction
            face->Delete_FaceFace_Interaction(cffi);
            face->temp_tag = 1;
          }
        } else {
          num_query++;
        }
      }
    }
  }
  for( i=0 ; i<number_of_element_blocks ; ++i){
    element_blocks[i]->ElemList()->IteratorStart();
    while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      element->temp_tag = 0;
      ContactInteractionDLL<Real>* interactions =
        element->Get_ElementElement_Interactions();
      ContactInteractionEntity<Real>* interaction;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactElementElementInteraction* ceei =
          static_cast<ContactElementElementInteraction*>(interaction);
        if (ceei->MasterElementEntityData()->owner == my_proc) {
          int blk = ceei->MasterElementEntityData()->block_id;
          element = static_cast<ContactElement*>
            (element_blocks[blk]->ElemList()->Find(ceei->MasterElementEntityData()));
          if (!element) {
            // master element not present, delete interaction
            element->Delete_ElementElement_Interaction(ceei);
            element->temp_tag = 1;
          }
        } else {
          num_query++;
        }
      }
    }
  }

#ifndef CONTACT_NO_MPI
  LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
  int        zoltan_pid;
  ContactZoltanComm ZoltanHostidQueryComm( ContactZoltanComm::ZOLTAN_IMPORT);
                                           
  //==========================================================================
  // Process Node/Face Interactions...
  //==========================================================================
  for( i=0 ; i<number_of_node_blocks ; ++i){
    node_blocks[i]->NodeList()->IteratorStart();
    while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      ContactNodeEntityInteraction** interactions = 
	node->Get_NodeEntity_Interactions();
      for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
        ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[j]);
        if (cnfi==NULL) continue;
        if (!cnfi->Face()) {
          zoltan_lid[0] = CT_FACE;
          zoltan_lid[1] = cnfi->FaceEntityData()->index_in_owner_proc_array;
          zoltan_gid[0] = CT_FACE;
          zoltan_gid[1] = cnfi->FaceEntityData()->host_gid[0];
          zoltan_gid[2] = cnfi->FaceEntityData()->host_gid[1];
          zoltan_pid    = cnfi->FaceEntityData()->owner;
          ZoltanHostidQueryComm.Add_Import(zoltan_lid,zoltan_gid,zoltan_pid);
        }
      }
    }
  }
  
  //==========================================================================
  // Process Face/Face Interactions...
  //==========================================================================
  for( i=0 ; i<number_of_face_blocks ; ++i){
    face_blocks[i]->FaceList()->IteratorStart();
    while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
      if(interactions == NULL) continue;
      ContactInteractionEntity<Real>* interaction;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactFaceFaceInteraction<Real>* cffi = 
          static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
        if (!cffi->MasterFace()) {
          zoltan_lid[0] = CT_FACE;
          zoltan_lid[1] = cffi->MasterFaceEntityData()->index_in_owner_proc_array;
          zoltan_gid[0] = CT_FACE;
          zoltan_gid[1] = cffi->MasterFaceEntityData()->host_gid[0];
          zoltan_gid[2] = cffi->MasterFaceEntityData()->host_gid[1];
          zoltan_pid    = cffi->MasterFaceEntityData()->owner;
          ZoltanHostidQueryComm.Add_Import(zoltan_lid,zoltan_gid,zoltan_pid);
        }
      }
    }
  }
  
  //==========================================================================
  // Process ElementElement Interactions...
  //==========================================================================
  for( i=0 ; i<number_of_element_blocks ; ++i){
    element_blocks[i]->ElemList()->IteratorStart();
    while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
      ContactInteractionEntity<Real>* interaction;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactElementElementInteraction* ceei = 
          static_cast<ContactElementElementInteraction*>(interaction);
        if (!ceei->MasterElement()) {
          zoltan_lid[0] = CT_ELEMENT;
          zoltan_lid[1] = ceei->MasterElementEntityData()->index_in_owner_proc_array;
          zoltan_gid[0] = CT_ELEMENT;
          zoltan_gid[1] = ceei->MasterElementEntityData()->host_gid[0];
          zoltan_gid[2] = ceei->MasterElementEntityData()->host_gid[1];
          zoltan_pid    = ceei->MasterElementEntityData()->owner;
          ZoltanHostidQueryComm.Add_Import(zoltan_lid,zoltan_gid,zoltan_pid);
        }
      }
    }
  }

  if (contact_number_of_processors( SearchComm )>1) {
    query_linklist = new std::vector<ContactInteractionEntity<Real>::entity_data*>;
    int       num_import   = ZoltanHostidQueryComm.Num_Import();
    LB_ID_PTR import_gids  = ZoltanHostidQueryComm.Import_GIDS();
    LB_ID_PTR import_lids  = ZoltanHostidQueryComm.Import_LIDS();
    int*      import_procs = ZoltanHostidQueryComm.Import_Procs();
    int       num_export=0;
    LB_ID_PTR export_gids=NULL;
    LB_ID_PTR export_lids=NULL;
    int*      export_procs=NULL;
    search->Get_Zoltan()->Compute_Destinations(num_import,   import_gids,
                                               import_lids,  import_procs,
                                               &num_export,  &export_gids,
                                               &export_lids, &export_procs);
   ZoltanHostidQueryComm.Set_Export(num_export, export_gids, export_lids, export_procs);
    
    search->Get_Zoltan()->Set_HostidQueryCallBacks(search);
    search->Get_Zoltan()->Help_Migrate(num_import,  import_gids, 
                                       import_lids, import_procs,
                                       num_export,  export_gids, 
                                       export_lids, export_procs);
    
    ContactEntityDataHash* hostid_query_hash =
        new ContactEntityDataHash(query_linklist);
  
    //==========================================================================
    // Process Node/Face Interactions...
    //==========================================================================
    for( i=0 ; i<number_of_node_blocks ; ++i){
      node_blocks[i]->NodeList()->IteratorStart();
      while ((entity=node_blocks[i]->NodeList()->IteratorForward())) {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
        ContactNodeEntityInteraction** interactions = 
	  node->Get_NodeEntity_Interactions();
        for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
          ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[j]);
          if (cnfi==NULL) continue;
          if (!cnfi->Face()) {
            ContactInteractionEntity<Real>::entity_data*
                face_struct = hostid_query_hash->find( cnfi->FaceEntityData(), 0, NULL, 0 );
            if (!face_struct) {
              // master face not present, delete interaction
              node->Delete_NodeEntity_Interaction(cnfi);
              node->temp_tag = 1;
            }
          }
        }
      }
    }
    //==========================================================================
    // Process Face/Face Interactions...
    //==========================================================================
    for( i=0 ; i<number_of_face_blocks ; ++i){
      face_blocks[i]->FaceList()->IteratorStart();
      while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        ContactInteractionEntity<Real>* interaction;
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())) {
          ContactFaceFaceInteraction<Real>* cffi = 
            static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
          if (!cffi->MasterFace()) {
            ContactInteractionEntity<Real>::entity_data*
                face_struct = hostid_query_hash->find( cffi->MasterFaceEntityData(), 0, NULL, 0 );
            if (!face_struct) {
              // master face not present, delete interaction
              face->Delete_FaceFace_Interaction(cffi);
              face->temp_tag = 1;
            }
          }
        }
      }
    }
 
    //==========================================================================
    // Process Element/Element Interactions...
    //==========================================================================
    for( i=0 ; i<number_of_element_blocks ; ++i){
      element_blocks[i]->ElemList()->IteratorStart();
      while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
        ContactElement* element = static_cast<ContactElement*>(entity);
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        ContactInteractionEntity<Real>* interaction;
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())) {
          ContactElementElementInteraction* ceei = 
            static_cast<ContactElementElementInteraction*>(interaction);
          if (!ceei->MasterElement()) {
            ContactInteractionEntity<Real>::entity_data*
                element_struct = hostid_query_hash->find( ceei->MasterElementEntityData(), 0, NULL, 0 );
            if (!element_struct) {
              // master element not present, delete interaction
              element->Delete_ElementElement_Interaction(ceei);
              element->temp_tag = 1;
            }
          }
        }
      }
    }
    
    for(i = 0; i < query_linklist->size(); ++i) {
      ContactInteractionEntity<Real>::entity_data* entity_data = (*query_linklist)[i];      
      delete [] entity_data;
    }

    if(query_linklist) {
      delete query_linklist;
      query_linklist = NULL;
    }

    if (hostid_query_hash) delete hostid_query_hash;
  }
#endif
  
  //==========================================================================
  // Process Face/Face Interactions for compaction
  //==========================================================================
  for( i=0 ; i<number_of_face_blocks ; ++i){
    face_blocks[i]->FaceList()->IteratorStart();
    while ((entity=face_blocks[i]->FaceList()->IteratorForward())) {
      if (entity->temp_tag) {
        j = 0;
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        ContactInteractionEntity<Real>* interaction;
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())) {
          interaction->Index(j++);
        }
      }
    }
  }
 
  //==========================================================================
  // Process Element/Element Interactions for compaction
  //==========================================================================
  for( i=0 ; i<number_of_element_blocks ; ++i){
    element_blocks[i]->ElemList()->IteratorStart();
    while ((entity=element_blocks[i]->ElemList()->IteratorForward())) {
      if (entity->temp_tag) {
        j = 0;
        ContactElement* element = static_cast<ContactElement*>(entity);
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        ContactInteractionEntity<Real>* interaction;
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())) {
          interaction->Index(j++);
        }
      }
    }
  }
}
