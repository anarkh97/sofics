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


#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactCommBuffer.h"
#include "ContactErrors.h"
#include "ContactTopologyEntity.h"
#include "ContactInteractionEntity.h"
#include "ContactTopologyEntityHash.h"
#include "ContactEntityDataHash.h"
#include "ContactNodeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactSearch.h"
#include "ContactSearchData.h"
#include "ContactTable.h"
#include "ContactTopology.h"
#include "ContactSymComm.h"
#include "ContactFixedSizeAllocator.h"
#include <cstring>
#include "contact_tolerances.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#include "ContactZoltanComm.h"
#endif


int ContactSearch::Restart_Size(){
  int size = 0;
  int i,k;
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;
  
  // Add the initial size parameters
  size += 33; // NOTE this has to be updated for new Search_Options

  // Add the tables if any
  for( i=0 ; i<num_tables ; ++i){
    size += tables[i]->Restart_Size();
  }

#ifndef CONTACT_NO_MPI
  // Add the node communication map
  ContactSymComm* Node_SymComm = primary_topology->Node_Sym_Comm();
  size += 2*Node_SymComm->Num_Comm_Partners()+1;
  for (i=0; i<Node_SymComm->Num_Comm_Partners(); ++i) {
    size += Node_SymComm->Num_to_Proc(i);
  }
#endif
  
  // Add the Type and Size Data for the node blocks
  size += 2*primary_topology->Number_of_Node_Blocks();
  
  // Add the Type and Size Data for the face blocks
  size += 2*primary_topology->Number_of_Face_Blocks();

  // Add the size of the face connectivity array
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      for( k=0 ; k<face->Nodes_Per_Face() ; ++k){
	size += 1;
      }
    }
  }
  
  // Add the Type and Size Data for the element blocks
  size += 2*primary_topology->Number_of_Element_Blocks();

  // Add the size of the element connectivity array
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    primary_topology->Element_Block(i)->ElemList()->IteratorStart();
    while ((entity=primary_topology->Element_Block(i)->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      for( k=0 ; k<element->Nodes_Per_Element() ; ++k){
	size += 1;
      }
    }
  }

  // Add the size for the nodal exodus ids
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i)
    size += primary_topology->Node_Block(i)->Number_of_Nodes();
    
  // Add the size for the node host ordering
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i)
    size += primary_topology->Node_Block(i)->Number_of_Nodes();
    
  // Add the size for the face host ordering
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i)
    size += primary_topology->Face_Block(i)->Number_of_Faces();
    
  // Add the size for the element host ordering
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i)
    size += primary_topology->Element_Block(i)->Number_of_Elements();

  // Add the size of the node data arrays
  primary_topology->NodeList()->IteratorStart();
  while ((entity=primary_topology->NodeList()->IteratorForward())) {
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    size += node->Restart_Size();
  }
  
  // Add the size of the edge data arrays
  primary_topology->EdgeList()->IteratorStart();
  while ((entity=primary_topology->EdgeList()->IteratorForward())) {
    ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity);
    size += edge->Restart_Size();
  }

  // Add the size of the face data arrays
  primary_topology->FaceList()->IteratorStart();
  while ((entity=primary_topology->FaceList()->IteratorForward())) {
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    size += face->Restart_Size();
  }

  // Add the size of the element data arrays
  primary_topology->ElemList()->IteratorStart();
  while ((entity=primary_topology->ElemList()->IteratorForward())) {
    ContactElement* element = static_cast<ContactElement*>(entity);
    size += element->Restart_Size();
  }

  // Add the size of the search data
  int num_entity_keys = search_data->Num_Search_Entities();
  size += NSIZSD*num_entity_keys*num_entity_keys;

  // Add the size for each of the analytic surfaces
  for( i=0 ; i<primary_topology->Number_of_Analytic_Surfaces() ; ++i){
    size += 1;  // Type of analytic surface
    size += primary_topology->Analytic_Surface(i)->Restart_Size();
  }

  // Add the number of node-face interactions (Current State only)
  size += 1;

  // Add the size of node-face interactions  (Current State only)
  primary_topology->NodeList()->IteratorStart();
  while ((entity=primary_topology->NodeList()->IteratorForward())) {
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[i]);
      if (cnfi!=NULL) {
        size += cnfi->Restart_Size();
        cnfi->Set_NodeEntityData();
      }
    }
  }

  // Add the number of node-surface interactions (Current State only)
  size += 1;

  // Add the size of node-surface interactions  (Current State only)
  primary_topology->NodeList()->IteratorStart();
  while ((entity=primary_topology->NodeList()->IteratorForward())) {
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(interactions[i]);
      if (cnsi!=NULL)  {
        size += cnsi->Restart_Size();
        cnsi->Set_NodeEntityData();
      }
    }
  }

  // Add the number of face-face interactions (Current State only)
  size += 1;

  // Add the size of face-face interactions  (Current State only)
  primary_topology->FaceList()->IteratorStart();
  while ((entity=primary_topology->FaceList()->IteratorForward())) {
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      if(interactions == NULL) continue;
      while ((interaction=interactions->IteratorForward())){
        ContactFaceFaceInteraction<Real>* cffi = 
          static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
        size += cffi->Restart_Size();
      }
    }
  }

  // Add the number of face-coverage interactions (Current State only)
  size += 1;

  // Add the size of face-coverage interactions  (Current State only)
  primary_topology->FaceList()->IteratorStart();
  while ((entity=primary_topology->FaceList()->IteratorForward())) {
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())){
        ContactFaceCoverageInteraction* cfci = 
          static_cast<ContactFaceCoverageInteraction*>(interaction);
        size += cfci->Restart_Size();
      }
    }
  }
  
  // Add the number of element-element interactions (Current State only)
  size += 1;

  // Add the size of element-element interactions  (Current State only)
  primary_topology->ElemList()->IteratorStart();
  while ((entity=primary_topology->ElemList()->IteratorForward())) {
    ContactElement* element = static_cast<ContactElement*>(entity);
    ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
    interactions->IteratorStart();
    while ((interaction=interactions->IteratorForward())){
      ContactElementElementInteraction* ceei = 
        static_cast<ContactElementElementInteraction*>(interaction);
      size += ceei->Restart_Size();
    }
  }

  return size;
}

ContactSearch::ContactErrorCode 
ContactSearch::Extract_Restart_Data( Real* buffer )
{
  int i,j,k;
  Real* buf_loc = buffer;
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;

  *buf_loc++ = step_number;
  *buf_loc++ = num_tracked_nodes;
  *buf_loc++ = dimensionality;
  *buf_loc++ = num_states;
  *buf_loc++ = primary_topology->Number_of_Node_Blocks();
  *buf_loc++ = primary_topology->Number_of_Face_Blocks();
  *buf_loc++ = primary_topology->Number_of_Element_Blocks();
  *buf_loc++ = primary_topology->Number_of_Analytic_Surfaces();
  *buf_loc++ = multiple_interaction_status;
  *buf_loc++ = normal_smoothing_status;
  *buf_loc++ = compute_node_areas;
  *buf_loc++ = partition_gap_status;
  *buf_loc++ = old_dynamic_search;
  *buf_loc++ = tracking_type;
  *buf_loc++ = enable_off_face_tracking;
  *buf_loc++ = no_secondary;
  *buf_loc++ = search_cull;
  *buf_loc++ = no_warped_volume;
  *buf_loc++ = no_parallel_consistency;
  *buf_loc++ = orig_sharp_smooth_angle;
  *buf_loc++ = normal_smoothing_distance;
  *buf_loc++ = smoothing_resolution;
  *buf_loc++ = global_search_cull;
  *buf_loc++ = auto_tol;
  *buf_loc++ = aggressive_tolerances;
  *buf_loc++ = skip_physical_faces;
  *buf_loc++ = physical_face_algorithm;
  *buf_loc++ = box_inflation;
  *buf_loc++ = gap_inflation;
  *buf_loc++ = num_tables;
  *buf_loc++ = tracking_step;
  *buf_loc++ = global_tracking_interval;
  *buf_loc++ = initialized_tied;

  // Write out the tables if any
  if( num_tables ){
    for( i=0 ; i<num_tables ; ++i){
      buf_loc += tables[i]->Extract_Restart_Data( buf_loc );
    }
  } else
    tables = NULL;

#ifndef CONTACT_NO_MPI
  // Write the node communication map
  ContactSymComm* Node_SymComm = primary_topology->Node_Sym_Comm();
  *buf_loc++ = Node_SymComm->Num_Comm_Partners();
  for (i=0; i<Node_SymComm->Num_Comm_Partners(); ++i) {
    *buf_loc++ = Node_SymComm->Comm_Proc_ID(i);
  }
  for (i=0; i<Node_SymComm->Num_Comm_Partners(); ++i) {
    *buf_loc++ = Node_SymComm->Num_to_Proc(i);
  }
  for (i=0; i<Node_SymComm->Num_Comm_Partners(); ++i) {
    ContactNode<Real>** nodes = (ContactNode<Real>**)Node_SymComm->Entity_List(i);
    for (j=0; j<Node_SymComm->Num_to_Proc(i); ++j) {
      //*buf_loc++ = nodes[j]->ProcArrayIndex()+1;
      *buf_loc++ = nodes[j]->HostGlobalArrayIndex()+1;
    }
  }
#endif
  
  // Write the node block types and number of nodes in each block
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    *buf_loc++ = primary_topology->Node_Block(i)->Type();
    *buf_loc++ = primary_topology->Node_Block(i)->Number_of_Nodes();
  }

  // Write the face block_types and number of faces in each block
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    *buf_loc++ = primary_topology->Face_Block(i)->Type();
    *buf_loc++ = primary_topology->Face_Block(i)->Number_of_Faces();
  }
  
  // Write the face connectivity
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      for( k=0 ; k<face->Nodes_Per_Face() ; ++k){
	//*buf_loc++ = face->Node(k)->ProcArrayIndex() + 1;
	*buf_loc++ = face->Node(k)->HostGlobalArrayIndex() + 1;
      }
    }
  }

  // Write the element block_types and number of elements in each block
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    *buf_loc++ = primary_topology->Element_Block(i)->Type();
    *buf_loc++ = primary_topology->Element_Block(i)->Number_of_Elements();
  }

  // Write the element connectivity
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    primary_topology->Element_Block(i)->ElemList()->IteratorStart();
    while ((entity=primary_topology->Element_Block(i)->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      for( k=0 ; k<element->Nodes_Per_Element() ; ++k){
	//*buf_loc++ = element->Node(k)->ProcArrayIndex() + 1;
	*buf_loc++ = element->Node(k)->HostGlobalArrayIndex() + 1;
      }
    }
  }

  // Write the Node Exodus IDs
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      *buf_loc++ = node->Exodus_ID();
    }
  }
  
  // Write the Node host ordering
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      *buf_loc++ = entity->HostGlobalArrayIndex();
    }
  }
  
  // Write the Face host ordering
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      *buf_loc++ = entity->HostGlobalArrayIndex();
    }
  }
  
  // Write the Element host ordering
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    primary_topology->Element_Block(i)->ElemList()->IteratorStart();
    while ((entity=primary_topology->Element_Block(i)->ElemList()->IteratorForward())) {
      *buf_loc++ = entity->HostGlobalArrayIndex();
    }
  }

  // Write the node data arrays
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      node->Restart_Pack( buf_loc );
      buf_loc += node->Restart_Size();
    }
  }

  // Write the edge data arrays
  for( i=0 ; i<primary_topology->Number_of_Edge_Blocks() ; ++i){
    primary_topology->Edge_Block(i)->EdgeList()->IteratorStart();
    while ((entity=primary_topology->Edge_Block(i)->EdgeList()->IteratorForward())) {
      ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity);
      edge->Restart_Pack( buf_loc );
      buf_loc += edge->Restart_Size();
    }
  }

  // Write the face data arrays
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      face->Restart_Pack( buf_loc );
      buf_loc += face->Restart_Size();
    }
  }

  // Write the element data arrays
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    primary_topology->Element_Block(i)->ElemList()->IteratorStart();
    while ((entity=primary_topology->Element_Block(i)->ElemList()->IteratorForward())) {
    ContactElement* element = static_cast<ContactElement*>(entity);
      element->Restart_Pack( buf_loc );
      buf_loc += element->Restart_Size();
    }
  }

  // Write the search data
  search_data->SetAll();  //reset original search data
  int num_entity_keys = search_data->Num_Search_Entities();
  std::memcpy( buf_loc, search_data->Search_Data(), 
	  sizeof(Real)*NSIZSD*num_entity_keys*num_entity_keys );
  buf_loc += NSIZSD*num_entity_keys*num_entity_keys;

  // Write the analytic surface data
  for( i=0 ; i<primary_topology->Number_of_Analytic_Surfaces() ; ++i){
    *buf_loc++ = primary_topology->Analytic_Surface(i)->Surface_Type();
    buf_loc += primary_topology->Analytic_Surface(i)->
      Extract_Restart_Data( buf_loc );
  }

  // Write the number of node-face interactions (Current State only)
  int num_nodeface_interactions = 0;
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      num_nodeface_interactions += node->Number_NodeFace_Interactions();
    }
  }
  *buf_loc++ = num_nodeface_interactions;

  // Write the node-face interactions (Current State only)
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
      for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
        ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[j]);
        if (cnfi!=NULL) {
          cnfi->Restart_Pack( buf_loc );
          buf_loc += cnfi->Restart_Size();
        }
      }
    }
  }
  
  // Write the number of node-surface interactions (Current State only)
  int num_nodesurface_interactions = 0;
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      num_nodesurface_interactions += node->Number_NodeSurface_Interactions();
    }
  }
  *buf_loc++ = num_nodesurface_interactions;

  // Write the node-surface interactions (Current State only)
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    primary_topology->Node_Block(i)->NodeList()->IteratorStart();
    while ((entity=primary_topology->Node_Block(i)->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
      for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
        ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(interactions[j]);
        if (cnsi!=NULL) {
          cnsi->Restart_Pack( buf_loc );
          buf_loc += cnsi->Restart_Size();
        }
      }
    }
  }
  
  // Write the number of face-face interactions (Current State only)
  int num_faceface_interactions = 0;
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      num_faceface_interactions += face->Number_FaceFace_Interactions();
    }
  }
  *buf_loc++ = num_faceface_interactions;
  
  // Write the face-face interactions (Current State only)
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
      if(interactions == NULL) continue;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())){
        ContactFaceFaceInteraction<Real>* cffi = 
          static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
        cffi->Restart_Pack( buf_loc );
        buf_loc += cffi->Restart_Size();
      }
    }
  }
  
  // Write the number of face-coverage interactions (Current State only)
  int num_facecoverage_interactions = 0;
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      num_faceface_interactions += face->Number_FaceCoverage_Interactions();
    }
  }
  *buf_loc++ = num_facecoverage_interactions;
  
  // Write the face-coverage interactions (Current State only)
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    primary_topology->Face_Block(i)->FaceList()->IteratorStart();
    while ((entity=primary_topology->Face_Block(i)->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
      if(interactions != NULL) {
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          ContactFaceCoverageInteraction* cfci = 
            static_cast<ContactFaceCoverageInteraction*>(interaction);
          cfci->Restart_Pack( buf_loc );
          buf_loc += cfci->Restart_Size();
        }
      }
    }
  }
  
  // Write the number of element-element interactions (Current State only)
  int num_elementelement_interactions = 0;
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    primary_topology->Element_Block(i)->ElemList()->IteratorStart();
    while ((entity=primary_topology->Element_Block(i)->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      num_elementelement_interactions += element->Number_ElementElement_Interactions();
    }
  }
  *buf_loc++ = num_elementelement_interactions;
  
  // Write the element-element interactions (Current State only)
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    primary_topology->Element_Block(i)->ElemList()->IteratorStart();
    while ((entity=primary_topology->Element_Block(i)->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())){
        ContactElementElementInteraction* ceei = 
          static_cast<ContactElementElementInteraction*>(interaction);
        ceei->Restart_Pack( buf_loc );
        buf_loc += ceei->Restart_Size();
      }
    }
  }

  POSTCONDITION( Restart_Size() == buf_loc-buffer );
  return ContactSearch::NO_ERROR;
}
  
ContactSearch::ContactSearch( Real* buffer, 
                              const int* Node_Host_IDs, 
                              const int* Face_Host_IDs, 
                              const int* Element_Host_IDs, 
                              MPI_Comm& mpi_communicator,
			      ContactErrorCode& error )
  : initialized(false), 
    initialized_tied(false), 
    initialized_context(false),
    initializing_tied(false),
    step_number(0),
    postream( mpi_communicator ),
    timer( mpi_communicator ),
#ifndef CONTACT_NO_MPI
    zoltan(NULL),
#endif
    num_tracked_nodes(0)
{
  
  restart = true;
  // compute initial internal tolerances
  box_inflation = BOX_INFLATION_FACTOR;
  gap_inflation = GAP_INFLATION_FACTOR+1.0;

  //
  //  Initialize static face data arrays for fast lookup of face type
  //  information from the base class
  //
  ContactFace<Real>::Initialize_Lookup_Arrays();
  ContactEdge<Real>::Initialize_Lookup_Arrays();

  error = NO_ERROR;

  errors = new ContactErrors();
  comm_buffer = new ContactCommBuffer();
 
  Set_Up_Allocators();

  SearchComm = mpi_communicator;

  int i,j,k;
  Real* buf_loc = buffer;

  step_number = (int) *buf_loc++;
  num_tracked_nodes = (int) *buf_loc++;
  dimensionality = (int) *buf_loc++;
  num_states = (int) *buf_loc++;
  int num_node_blocks = (int) *buf_loc++;
  int num_face_blocks = (int) *buf_loc++;
  int num_element_blocks = (int) *buf_loc++;
  int num_analytic_surfaces = (int) *buf_loc++;
  multiple_interaction_status = (Search_Option_Status) *buf_loc++;
  normal_smoothing_status = (Search_Option_Status) *buf_loc++;
  compute_node_areas = (Search_Option_Status) *buf_loc++;
  partition_gap_status = (Search_Option_Status) *buf_loc++;
  old_dynamic_search = (Search_Option_Status) *buf_loc++;
  tracking_type = (Track_Type) *buf_loc++;
  enable_off_face_tracking = (int) *buf_loc++;
  no_secondary = (Search_Option_Status) *buf_loc++;
  search_cull = (Search_Option_Status) *buf_loc++;
  no_warped_volume = (Search_Option_Status) *buf_loc++;
  no_parallel_consistency = (Search_Option_Status) *buf_loc++;
  orig_sharp_smooth_angle = *buf_loc++;
  sharp_smooth_curvature = ComputeCurvatureFromAngle(orig_sharp_smooth_angle);
  normal_smoothing_distance = *buf_loc++;
  smoothing_resolution = (Smoothing_Resolution) *buf_loc++;
  global_search_cull = (Search_Cull) *buf_loc++;
  auto_tol = (Search_Option_Status) *buf_loc++;
  aggressive_tolerances = (Search_Option_Status) *buf_loc++;
  skip_physical_faces = (Search_Option_Status) *buf_loc++;
  physical_face_algorithm = (PF_Algorithm) *buf_loc++;
  box_inflation = *buf_loc++;
  gap_inflation = *buf_loc++;
  num_tables = (int) *buf_loc++;
  tracking_step = (int) *buf_loc++;
  global_tracking_interval = (int) *buf_loc++;
  initialized_tied =(bool) *buf_loc++;
  if (tracking_type==NO_TRACKING) {
    enable_tracking = INACTIVE;
  } else {
    enable_tracking = ACTIVE;
    if (tracking_type==GLOBAL_TRACKING) tracking_step=0;
  }

  // this computes the actual tolerances
  Set_Search_Option(AUTO_TOL,auto_tol,NULL);
  Set_Search_Option(AGGRESSIVE_TOLERANCES,aggressive_tolerances,NULL);

  if( num_tables ){
    tables = new ContactTable*[num_tables];
    for( i=0 ; i<num_tables ; ++i){
      tables[i] = new ContactTable();
      buf_loc += tables[i]->Implant_Restart_Data( buf_loc );
    }
  } else
    tables = NULL;

#ifndef CONTACT_NO_MPI
  // Read the node communication maps
  int  nsize=0;
  int* Node_Comm_Proc_IDs = NULL;
  int* Number_Nodes_To_Partner = NULL;
  int* Communication_Nodes = NULL;
  int  Number_of_Nodal_Comm_Partners = (int) *buf_loc++;
  if (Number_of_Nodal_Comm_Partners>0) {
    Node_Comm_Proc_IDs = new int [Number_of_Nodal_Comm_Partners];
    Number_Nodes_To_Partner = new int [Number_of_Nodal_Comm_Partners];
    for (i=0; i<Number_of_Nodal_Comm_Partners; ++i) {
      Node_Comm_Proc_IDs[i] = (int) *buf_loc++;
    }
    for (i=0; i<Number_of_Nodal_Comm_Partners; ++i) {
      Number_Nodes_To_Partner[i] = (int) *buf_loc++;
      nsize += Number_Nodes_To_Partner[i];
    }
    Communication_Nodes = new int [nsize];
    for (k=0, i=0; i<Number_of_Nodal_Comm_Partners; ++i) {
      for (j=0; j<Number_Nodes_To_Partner[i]; ++j) {
        Communication_Nodes[k++] = (int) *buf_loc++;
      }
    }
  }
#else
  int  Number_of_Nodal_Comm_Partners = 0;
  int* Node_Comm_Proc_IDs = NULL;
  int* Number_Nodes_To_Partner = NULL;
  int* Communication_Nodes = NULL;
#endif
  
  // Read the Node Block Types and Number of Nodes in each Block
  int num_nodes = 0;
  ContactNode_Type* Node_Block_Types = new ContactNode_Type[num_node_blocks];
  int* Number_Nodes_in_Blocks = new int[num_node_blocks];
  for( i=0 ; i<num_node_blocks ; ++i){
    Node_Block_Types[i] = (ContactNode_Type) *buf_loc++;
    Number_Nodes_in_Blocks[i] = (int) *buf_loc++;
    num_nodes += Number_Nodes_in_Blocks[i];
  }

  // Read the Face Block Types and Number of Faces in each Block
  int num_faces = 0;
  ContactFace_Type* Face_Block_Types = new ContactFace_Type[num_face_blocks];
  int* Number_Faces_in_Blocks = new int[num_face_blocks];
  int size_face_connectivity = 0;
  for( i=0 ; i<num_face_blocks ; ++i){
    Face_Block_Types[i] = (ContactFace_Type) *buf_loc++; 
    Number_Faces_in_Blocks[i] = (int) *buf_loc++;
    num_faces += Number_Faces_in_Blocks[i];
    int num_nodes_per_face = ContactSearch::Number_Nodes_Per_Face(Face_Block_Types[i]);
    size_face_connectivity += num_nodes_per_face*Number_Faces_in_Blocks[i];
  }   
  
  int* face_connectivity = new int[size_face_connectivity];
  for( i=0 ; i<size_face_connectivity ; ++i)
    face_connectivity[i] = (int) *buf_loc++;  

  // Read the Element Block Types and Number of Elements in each Block
  int num_elements = 0;
  ContactElement_Type* Element_Block_Types = new ContactElement_Type[num_element_blocks];
  int* Number_Elements_in_Blocks = new int[num_element_blocks];
  int size_element_connectivity = 0;
  for( i=0 ; i<num_element_blocks ; ++i){
    Element_Block_Types[i] = (ContactElement_Type) *buf_loc++; 
    Number_Elements_in_Blocks[i] = (int) *buf_loc++;
    num_elements += Number_Elements_in_Blocks[i];
    switch( Element_Block_Types[i] ){
    case( CARTESIANHEXELEMENTL8 ):
    case( HEXELEMENTL8 ):
      size_element_connectivity += 8*Number_Elements_in_Blocks[i];
      break;
    default:
      POSTCONDITION( 0 );
    }
  }   
  
  int* element_connectivity = new int[size_element_connectivity];
  for( i=0 ; i<size_element_connectivity ; ++i)
    element_connectivity[i] = (int) *buf_loc++;
    
  // Read the Node Exodus IDs
  int* Node_Exodus_IDs = new int[num_nodes];
  for( i=0 ; i<num_nodes ; ++i){
    Node_Exodus_IDs[i] = (int) *buf_loc++;
  }
  
  k = 0;
  int* NodeMap = new int[num_nodes];
  int* Mapped_Node_Host_IDs = new int[2*num_nodes];
  for( i=0 ; i<num_node_blocks ; ++i){
    for (j=0; j<Number_Nodes_in_Blocks[i]; ++j) {
      int map = (int) *buf_loc++;
      NodeMap[k] = map;
      Mapped_Node_Host_IDs[2*k  ] = Node_Host_IDs[2*map];
      Mapped_Node_Host_IDs[2*k+1] = Node_Host_IDs[2*map+1];
      k++;
    }
  }
  
  k = 0;
  int* FaceMap = new int[num_faces];
  int* Mapped_Face_Host_IDs = new int[2*num_faces];
  for( i=0 ; i<num_face_blocks ; ++i){
    for (j=0; j<Number_Faces_in_Blocks[i]; ++j) {
      int map = (int) *buf_loc++;
      FaceMap[k] = map;
      Mapped_Face_Host_IDs[2*k  ] = Face_Host_IDs[2*map];
      Mapped_Face_Host_IDs[2*k+1] = Face_Host_IDs[2*map+1];
      k++;
    }
  }
  
  k = 0;
  int* ElementMap = new int[num_elements];
  int* Mapped_Element_Host_IDs = new int[2*num_elements];
  for( i=0 ; i<num_element_blocks ; ++i){
    for (j=0; j<Number_Elements_in_Blocks[i]; ++j) {
      int map = (int) *buf_loc++;
      ElementMap[k] = map;
      Mapped_Element_Host_IDs[2*k  ] = Element_Host_IDs[2*map];
      Mapped_Element_Host_IDs[2*k+1] = Element_Host_IDs[2*map+1];
      k++;
    }
  }

  // Construct the primary_topology
  primary_topology = new ContactTopology( errors, dimensionality,
					  num_analytic_surfaces, 
					  num_node_blocks,
					  Node_Block_Types,
					  Number_Nodes_in_Blocks, 
					  Node_Exodus_IDs,
					  Mapped_Node_Host_IDs,
					  //Node_Host_IDs,
					  // NOTE -- ASG 2/3/03
					  // the next entry is supposed to
					  // be the nodal coordinates, as
					  // needed by shell contact. Right
					  // now we set this to null, since
					  // we can't do shell restart
					  // with binary restart anyway.
					  // This will be caught in a 
					  // PRECONDITION in ShellHandler
					  // if called with this restart.
					  NULL,
					  num_face_blocks,
					  Face_Block_Types,
					  Number_Faces_in_Blocks, 
					  Mapped_Face_Host_IDs,
					  //Face_Host_IDs,
					  face_connectivity,
					  // NOTE -- ASG 2/3/03
					  // the next entry is supposed to
					  // be the shell lofting factors, as
					  // needed by shell contact. Right
					  // now we set this to null, since
					  // we can't do shell restart
					  // with binary restart anyway.
					  // This will be caught in a 
					  // PRECONDITION in ShellHandler
					  // if called with this restart.
					  NULL,
                                          num_element_blocks,
                                          Element_Block_Types,
                                          Number_Elements_in_Blocks,
                                          Mapped_Element_Host_IDs,
                                          //Element_Host_IDs,
					  element_connectivity,
					  Number_of_Nodal_Comm_Partners, 
					  Node_Comm_Proc_IDs, 
					  Number_Nodes_To_Partner,
					  Communication_Nodes,
					  comm_buffer,
					  SearchComm,
                                          this,
					  error );

  primary_topology->topology_type = ContactTopology::PRIMARY;

  // If the topology is bad, simply return
  // This shouldn't be possible since its a restart but check anyway
  error = (ContactErrorCode) contact_global_error_check( error, SearchComm );
  if( error ) return;
    
  ContactTopologyEntity<Real>* entity;
  
  // Read the node data arrays
  j = 0;
  for( i=0 ; i<primary_topology->Number_of_Node_Blocks() ; ++i){
    ContactNodeBlock* block = primary_topology->Node_Block(i);
    block->NodeList()->IteratorStart();
    while ((entity = block->NodeList()->IteratorForward())) {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      node->HostGlobalArrayIndex(NodeMap[j]);
      node->Restart_Unpack( buf_loc );
      buf_loc += node->Restart_Size();
      j++;
    }
  }

  // Read the edge data arrays
  for( i=0 ; i<primary_topology->Number_of_Edge_Blocks() ; ++i){
    primary_topology->Edge_Block(i)->EdgeList()->IteratorStart();
    while ((entity=primary_topology->Edge_Block(i)->EdgeList()->IteratorForward())) {
      ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity);
      edge->Restart_Unpack( buf_loc );
      buf_loc += edge->Restart_Size();
    }
  }

  // Read the face data arrays
  j = 0;
  for( i=0 ; i<primary_topology->Number_of_Face_Blocks() ; ++i){
    ContactFaceBlock* block = primary_topology->Face_Block(i);
    block->FaceList()->IteratorStart();
    while ((entity = block->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      face->HostGlobalArrayIndex(FaceMap[j]);
      face->Restart_Unpack( buf_loc );
      buf_loc += face->Restart_Size();
      j++;
    }
  }

  // Read the element data arrays
  j = 0;
  for( i=0 ; i<primary_topology->Number_of_Element_Blocks() ; ++i){
    ContactElementBlock* block = primary_topology->Element_Block(i);
    block->ElemList()->IteratorStart();
    while ((entity = block->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      element->HostGlobalArrayIndex(ElementMap[j]);
      element->Restart_Unpack( buf_loc );
      buf_loc += element->Restart_Size();
      j++;
    }
  }
  
  delete[] NodeMap;
  delete[] FaceMap;
  delete[] ElementMap;
  delete[] Mapped_Node_Host_IDs;
  delete[] Mapped_Face_Host_IDs;
  delete[] Mapped_Element_Host_IDs;
  
  // Restore the search data
  search_data = new ContactSearchData(primary_topology);
  search_data->Set_Search_Data( buf_loc );
  int num_entity_keys = search_data->Num_Search_Entities();
  buf_loc += NSIZSD*num_entity_keys*num_entity_keys;
  
  // Restore the analytic surfaces
  if( num_analytic_surfaces ){
    for( i=0 ; i<num_analytic_surfaces ; ++i){
      ContactSearch::AnalyticSurface_Type as_type = 
	(ContactSearch::AnalyticSurface_Type) *buf_loc++;
      ContactAnalyticSurface* asurf=NULL;
      switch( as_type ){
      case ContactSearch::PLANE :
	asurf = new ContactAnalyticPlane( i );
	break;
      case ContactSearch::SPHERE :
	asurf = new ContactAnalyticSphere( i );
	break;
      case ContactSearch::CYLINDER_INSIDE :
	asurf = new ContactAnalyticCylinderInside( i );
	break;
      case ContactSearch::CYLINDER_OUTSIDE :
	asurf = new ContactAnalyticCylinderOutside( i );
	break;
      default:
        POSTCONDITION(0);
        break;
      }
      POSTCONDITION(asurf!=NULL);
      buf_loc += asurf->Implant_Restart_Data( buf_loc );
      int surface_key = primary_topology->Number_of_Added_Analytic_Surfaces() +
                        primary_topology->Number_of_Node_Blocks() - 1 +
                        primary_topology->Number_of_Face_Blocks();
      asurf->Entity_Key(surface_key);
      primary_topology->Add_Analytic_Surface( asurf );
    }
  }
  
  int my_proc = contact_processor_number(SearchComm);

  // Read the number of node-face interactions (Current State only)
  int num_nodeface_interactions = (int) *buf_loc++;

  // Read the node-face interactions (Current State only)
  for( i=0 ; i<num_nodeface_interactions ; ++i){
    ContactNodeFaceInteraction* cnfi = 
      ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
         allocators[ALLOC_ContactNodeFaceInteraction]);
    cnfi->Restart_Unpack( buf_loc );
    buf_loc += cnfi->Restart_Size();
    // need to do this because host_ids are not
    // guarrenteed to be the same across a restart.
    int node_index = cnfi->NodeEntityData()->index_in_owner_proc_array;
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                        (primary_topology->NodeList()->Find(node_index));
    cnfi->Connect_Node(node);
    if (cnfi->FaceEntityData()->owner == my_proc) {
      int face_index = cnfi->FaceEntityData()->index_in_owner_proc_array;
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                          (primary_topology->FaceList()->Find(face_index));
      cnfi->Connect_Face(face);
    }
    node->Add_NodeEntity_Interaction( cnfi );
  }

  // Read the number of node-face interactions (Current State only)
  int num_nodesurface_interactions = (int) *buf_loc++;
  
  // Read the node-surface interactions (Current State only)
  for( i=0 ; i<num_nodesurface_interactions ; ++i){
    ContactNodeSurfaceInteraction* cnsi = 
      ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction(
         allocators[ALLOC_ContactNodeSurfaceInteraction]);
    cnsi->Restart_Unpack( buf_loc );
    buf_loc += cnsi->Restart_Size();
    // need to do this because host_ids are not
    // guarrenteed to be the same across a restart.
    int node_index = cnsi->NodeEntityData()->index_in_owner_proc_array;
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                        (primary_topology->NodeList()->Find(node_index));
    cnsi->Connect_Node(node);
    cnsi->Connect_Surface( primary_topology );
    node->Add_NodeEntity_Interaction( cnsi );
  }
  
  // Read the number of face-face interactions (Current State only)
  int num_faceface_interactions = (int) *buf_loc++;

  // Read the face-face interactions (Current State only)
  for( i=0 ; i<num_faceface_interactions ; ++i){
    ContactFaceFaceInteraction<Real>* cffi = 
      ContactFaceFaceInteraction<Real>::new_ContactFaceFaceInteraction(
         allocators[ALLOC_ContactFaceFaceInteraction]);
    cffi->Restart_Unpack( buf_loc );
    buf_loc += cffi->Restart_Size();
    // need to do this because host_ids are not
    // guarenteed to be the same across a restart.
    int sface_index = cffi->SlaveFaceEntityData()->index_in_owner_proc_array;
    ContactFace<Real>* sface = static_cast<ContactFace<Real>*>
                         (primary_topology->FaceList()->Find(sface_index));
    cffi->Connect_SlaveFace(sface);
    if (cffi->MasterFaceEntityData()->owner == my_proc) {
      int mface_index = cffi->MasterFaceEntityData()->index_in_owner_proc_array;
      ContactFace<Real>* mface = static_cast<ContactFace<Real>*>
                           (primary_topology->FaceList()->Find(mface_index));
      cffi->Connect_MasterFace(mface);
    }
    sface->Store_FaceFace_Interaction( cffi );
  }
  
  // Read the number of face-coverage interactions (Current State only)
  int num_facecoverage_interactions = (int) *buf_loc++;
  
  // Read the face-coverage interactions (Current State only)
  for( i=0 ; i<num_facecoverage_interactions ; ++i){
    ContactFaceCoverageInteraction* cfci = 
      ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
         allocators[ALLOC_ContactFaceCoverageInteraction]);
    cfci->Restart_Unpack( buf_loc );
    buf_loc += cfci->Restart_Size();
    // need to do this because host_ids are not
    // guarenteed to be the same across a restart.
    int sface_index = cfci->SlaveFaceEntityData()->index_in_owner_proc_array;
    ContactFace<Real>* sface = static_cast<ContactFace<Real>*>
                         (primary_topology->FaceList()->Find(sface_index));
    cfci->Connect_SlaveFace(sface);
    sface->Store_FaceCoverage_Interaction( cfci );
  }
  
  // Read the number of element-element interactions (Current State only)
  int num_elementelement_interactions = (int) *buf_loc++;
  
  // Read the element-element interactions (Current State only)
  for( i=0 ; i<num_elementelement_interactions ; ++i){
    ContactElementElementInteraction* ceei = 
      ContactElementElementInteraction::new_ContactElementElementInteraction(
         allocators[ALLOC_ContactElementElementInteraction]);
    ceei->Restart_Unpack( buf_loc );
    buf_loc += ceei->Restart_Size();
    // need to do this because host_ids are not
    // guarenteed to be the same across a restart.
    int selement_index = ceei->SlaveElementEntityData()->index_in_owner_proc_array;
    ContactElement* selement = static_cast<ContactElement*>
                               (primary_topology->ElemList()->Find(selement_index));
    ceei->Connect_SlaveElement(selement);
    if (ceei->MasterElementEntityData()->owner == my_proc) {
      int melement_index = ceei->MasterElementEntityData()->index_in_owner_proc_array;
      ContactElement* melement = static_cast<ContactElement*>
                                 (primary_topology->ElemList()->Find(melement_index));
      ceei->Connect_MasterElement(melement);
    }
    selement->Store_ElementElement_Interaction( ceei );
  }

  POSTCONDITION( Restart_Size() == buf_loc-buffer );
  delete [] Number_Nodes_in_Blocks;
  delete [] Number_Faces_in_Blocks;
  delete [] Number_Elements_in_Blocks;
  delete [] face_connectivity;
  delete [] element_connectivity;

  //===================================================================
  // Create a topology object with the current entity block structures
  //===================================================================
  int num_edge_blocks = primary_topology->Number_of_Edge_Blocks();
  ContactEdge_Type* edge_block_types = new ContactEdge_Type[num_edge_blocks];
  for( i=0 ; i<num_edge_blocks ; ++i)
    edge_block_types[i] = primary_topology->Edge_Block(i)->Type();

  secondary_topology = new ContactTopology( errors,
                                            dimensionality,
                                            num_analytic_surfaces,
                                            num_node_blocks, 
                                            Node_Block_Types,
                                            num_edge_blocks,
                                            edge_block_types,
                                            num_face_blocks, 
                                            Face_Block_Types,
                                            num_element_blocks, 
                                            Element_Block_Types,
                                            SearchComm,
                                            this );

  secondary_topology->topology_type = ContactTopology::SECONDARY;

  delete [] edge_block_types;
  
#ifndef CONTACT_NO_MPI
  if (contact_number_of_processors(SearchComm)>1) {
    create_zoltan_object(error);
  }
#endif

  delete [] Node_Block_Types;
  delete [] Face_Block_Types;
  delete [] Element_Block_Types;
  delete [] Node_Exodus_IDs;
  if( Communication_Nodes ) delete [] Communication_Nodes;
  if( Node_Comm_Proc_IDs ) delete [] Node_Comm_Proc_IDs;
  if( Number_Nodes_To_Partner ) delete [] Number_Nodes_To_Partner;
  
  enforcement = NULL;
  number_registered_enforcements = 0;
  ResetInteractionHostIDs( );

#if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream << "Search Object is Constructed\n";
  postream << "  Number of Nodes    = " 
	   << primary_topology->Number_of_Nodes() << "\n";
  postream << "  Number of Edges    = " 
	   << primary_topology->Number_of_Edges() << "\n";
  postream << "  Number of Faces    = " 
	   << primary_topology->Number_of_Faces() << "\n";
  postream << "  Number of Elements = " 
	   << primary_topology->Number_of_Elements() << "\n";
  postream << "  Number of Surfaces = " 
	   << primary_topology->Number_of_Analytic_Surfaces() << "\n";
  postream.flush();
#endif

#ifdef CONTACT_TIMINGS
  Register_Timers();
#endif

  ws_size                = 1;
  max_facets             = 32;
  data_size              = ContactFace<Real>::LENGTH;
  ctrl                   = new Real[24*ws_size];
  ctrcl                  = new Real[data_size*ws_size];
  ctrcl_facets           = new Real[data_size*max_facets*ws_size];
  pushback_dir_flag      = new int[max_facets];
  ms_coordinates_c       = new Real[9*max_facets*ws_size];
  ms_coordinates_a       = new Real[9*max_facets*ws_size];
  ms_coordinates_p       = new Real[9*max_facets*ws_size];
  ms_normals_c           = new Real[3*max_facets*ws_size];
  ms_normals_a           = new Real[3*max_facets*ws_size];
  ms_normals_p           = new Real[3*max_facets*ws_size];
  /*
  dynamic_process_method = NULL;
  list                   = NULL;
  node_entity_keys       = NULL;
  physical_faces         = NULL;
  node_search_status     = NULL;
  face_search_status     = NULL;
  index                  = NULL;
  rank                   = NULL;
  rank2                  = NULL;
  scratch                = NULL;
  position_search        = NULL;
  position               = NULL;
  allocated_node_size    = 0;
  allocated_face_size    = 0;
  allocated_list_size    = 0;
  */
}

void
ContactSearch::ResetInteractionHostIDs( )
{
  // Need to do this step because HostIDs are not guarenteed to be
  // the same upon restart.  The ordering is the same but the IDs
  // themselves are not.  For example, if the host-code is doing 
  // birth/death, it may compact the IDs upon restart.  

  int i;
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;
  int num_query = 0;
  //==========================================================================
  // Process Node/Face Interactions...
  //==========================================================================
  primary_topology->NodeList()->IteratorStart();
  while( (entity=primary_topology->NodeList()->IteratorForward()) ){
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[i]);
      if (cnfi!=NULL) {
        cnfi->Set_NodeEntityData();
        if (!cnfi->Set_FaceEntityData()) num_query++;
      }
    }
  }

  //==========================================================================
  // Process Node/Surface Interactions...
  //==========================================================================
  primary_topology->NodeList()->IteratorStart();
  while( (entity=primary_topology->NodeList()->IteratorForward()) ){
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(interactions[i]);
      if (cnsi!=NULL) cnsi->Set_NodeEntityData();
    }
  }
  
  //==========================================================================
  // Process Face/Face Interactions...
  //==========================================================================
  primary_topology->FaceList()->IteratorStart();
  while( (entity=primary_topology->FaceList()->IteratorForward()) ){
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions == NULL) continue;
    interactions->IteratorStart();
    while ((interaction=interactions->IteratorForward())) {
      ContactFaceFaceInteraction<Real>* cffi = 
        static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
      cffi->Set_SlaveFaceEntityData();
      if (cffi->Set_MasterFaceEntityData()) num_query++;
    }
  }
  
  //==========================================================================
  // Process Face/Coverage Interactions...
  //==========================================================================
  primary_topology->FaceList()->IteratorStart();
  while( (entity=primary_topology->FaceList()->IteratorForward()) ){
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactFaceCoverageInteraction* cfci = 
          static_cast<ContactFaceCoverageInteraction*>(interaction);
        cfci->Set_SlaveFaceEntityData();
      }
    }
  }
  
  //==========================================================================
  // Process ElementElement Interactions...
  //==========================================================================
  primary_topology->ElemList()->IteratorStart();
  while( (entity=primary_topology->ElemList()->IteratorForward()) ){
    ContactElement* element = static_cast<ContactElement*>(entity);
    ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
    interactions->IteratorStart();
    while ((interaction=interactions->IteratorForward())) {
      ContactElementElementInteraction* ceei = 
        static_cast<ContactElementElementInteraction*>(interaction);
      ceei->Set_SlaveElementEntityData();
      if (ceei->Set_MasterElementEntityData()) num_query++;
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
  primary_topology->NodeList()->IteratorStart();
  while( (entity=primary_topology->NodeList()->IteratorForward()) ){
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[i]);
      if (cnfi!=NULL) {
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
  primary_topology->FaceList()->IteratorStart();
  while( (entity=primary_topology->FaceList()->IteratorForward()) ){
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions == NULL) continue;
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
  
  //==========================================================================
  // Process ElementElement Interactions...
  //==========================================================================
  primary_topology->ElemList()->IteratorStart();
  while( (entity=primary_topology->ElemList()->IteratorForward()) ){
    ContactElement* element = static_cast<ContactElement*>(entity);
    ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
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
#endif
 
#ifndef CONTACT_NO_MPI
  if (contact_number_of_processors( SearchComm )>1) {
    std::vector<ContactInteractionEntity<Real>::entity_data*> query_linklist;
    primary_topology->QueryLinkList(&query_linklist);
    int       num_import   = ZoltanHostidQueryComm.Num_Import();
    LB_ID_PTR import_gids  = ZoltanHostidQueryComm.Import_GIDS();
    LB_ID_PTR import_lids  = ZoltanHostidQueryComm.Import_LIDS();
    int*      import_procs = ZoltanHostidQueryComm.Import_Procs();
    int       num_export=0;
    LB_ID_PTR export_gids=NULL;
    LB_ID_PTR export_lids=NULL;
    int*      export_procs=NULL;
    zoltan->Compute_Destinations(num_import,   import_gids,
                                 import_lids,  import_procs,
                                 &num_export,  &export_gids,
                                 &export_lids, &export_procs);

    ZoltanHostidQueryComm.Set_Export(num_export, export_gids, export_lids, export_procs);
    
    zoltan->Set_HostidQueryCallBacks(this);
    zoltan->Help_Migrate(num_import,  import_gids, 
                         import_lids, import_procs,
                         num_export,  export_gids, 
                         export_lids, export_procs);
    
    ContactEntityDataHash* hostid_query_hash =
        new ContactEntityDataHash(&query_linklist);
  
    //==========================================================================
    // Process Node/Face Interactions...
    //==========================================================================
    primary_topology->NodeList()->IteratorStart();
    while( (entity=primary_topology->NodeList()->IteratorForward()) ){
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
      for (i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
        ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(interactions[i]);
        if (cnfi!=NULL) {
          if (!cnfi->Face()) {
            ContactInteractionEntity<Real>::entity_data*
                face_struct = hostid_query_hash->find( cnfi->FaceEntityData(), 0, NULL, 0 );
            PRECONDITION(face_struct);
            cnfi->FaceEntityData()->host_gid[0] = face_struct->host_gid[0];
            cnfi->FaceEntityData()->host_gid[1] = face_struct->host_gid[1];
          }
        }
      }
    }
 
    //==========================================================================
    // Process Face/Face Interactions...
    //==========================================================================
    primary_topology->FaceList()->IteratorStart();
    while( (entity=primary_topology->FaceList()->IteratorForward()) ){
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
      if(interactions == NULL) continue;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactFaceFaceInteraction<Real>* cffi = 
          static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
        if (!cffi->MasterFace()) {
          ContactInteractionEntity<Real>::entity_data*
              face_struct = hostid_query_hash->find( cffi->MasterFaceEntityData(), 0, NULL, 0 );
          PRECONDITION(face_struct);
          cffi->MasterFaceEntityData()->host_gid[0] = face_struct->host_gid[0];
          cffi->MasterFaceEntityData()->host_gid[1] = face_struct->host_gid[1];
        }
      }
    }
 
    //==========================================================================
    // Process Element/Element Interactions...
    //==========================================================================
    primary_topology->ElemList()->IteratorStart();
    while( (entity=primary_topology->ElemList()->IteratorForward()) ){
      ContactElement* element = static_cast<ContactElement*>(entity);
      ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())) {
        ContactElementElementInteraction* ceei = 
          static_cast<ContactElementElementInteraction*>(interaction);
        if (!ceei->MasterElement()) {
          ContactInteractionEntity<Real>::entity_data*
              element_struct = hostid_query_hash->find( ceei->MasterElementEntityData(), 0, NULL, 0 );
          PRECONDITION(element_struct);
          ceei->MasterElementEntityData()->host_gid[0] = element_struct->host_gid[0];
          ceei->MasterElementEntityData()->host_gid[1] = element_struct->host_gid[1];
        }
      }
    }
    
    for(i = 0; i < query_linklist.size(); ++i) {
      ContactInteractionEntity<Real>::entity_data *entity_data = query_linklist[i];
      delete [] entity_data;
    }
    if (hostid_query_hash) delete hostid_query_hash;
  }
#endif
}
