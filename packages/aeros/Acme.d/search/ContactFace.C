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


#ifndef ContactFace_C_
#define ContactFace_C_

#include "ContactFace.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "contact_assert.h"
#include <cstddef>
#include <iostream>
#include <cstring>
#include <cmath>

template<typename DataType>
  bool ContactFace<DataType>::array_init = false;
template<typename DataType>
  int ContactFace<DataType>::NODES_PER_FACE[ContactSearch::NFACE_TYPES] = {0};
template<typename DataType>
  int ContactFace<DataType>::EDGES_PER_FACE[ContactSearch::NFACE_TYPES] = {0};

template<typename DataType>
ContactFace<DataType>::ContactFace(ContactFixedSizeAllocator* alloc,
                         ContactSearch::ContactFace_Type Type, 
			 int Block_ID, int Host_Index_in_Block, int key,
                         ContactNode<DataType> **node_list_,
                         ContactEdge<DataType> **edge_list_,
                         typename ContactTopologyEntity<DataType>::connection_data *node_info_list_,
                         typename ContactTopologyEntity<DataType>::connection_data *edge_info_list_)
  : ContactTopologyEntity<DataType>( Block_ID, Host_Index_in_Block, DataArray, CT_FACE), 
    allocators(alloc),
    node_list(node_list_),
    edge_list(edge_list_),
    node_info_list(node_info_list_),
    edge_info_list(edge_info_list_) {

  face_type       = Type;
  entity_key = key;
  number_of_neighbors = 0;
  neighbor_face_info  = NULL;
}

template<typename DataType>
ContactFace<DataType>::~ContactFace()
{
  for(int i = 0; i < FaceFaceInteractions.size(); ++i) {
    ContactInteractionEntity<DataType>* link=NULL;
    FaceFaceInteractions[i].IteratorStart();
    while ((link=FaceFaceInteractions[i].IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* cffi = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(link);
      cffi->~ContactFaceFaceInteraction<DataType>();
      allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi);      
    }
  }
  for( int i=0 ; i<FaceCoverageInteractions.size() ; ++i ){
    ContactInteractionEntity<Real>* link=NULL;
    FaceCoverageInteractions[i].IteratorStart();
    while ((link=FaceCoverageInteractions[i].IteratorForward())) {
      ContactFaceCoverageInteraction* cfci = 
        static_cast<ContactFaceCoverageInteraction*>(link);
      cfci->~ContactFaceCoverageInteraction();
      allocators[ContactSearch::ALLOC_ContactFaceCoverageInteraction].Delete_Frag(cfci);
    }
  }
  FaceFaceInteractions.clear();
  FaceCoverageInteractions.clear();
  if (neighbor_face_info) delete [] neighbor_face_info;
}

template<typename DataType>
void ContactFace<DataType>::Initialize_Lookup_Arrays() {
  NODES_PER_FACE[ContactSearch::QUADFACEL4     ] = 4;
  NODES_PER_FACE[ContactSearch::QUADFACEQ8     ] = 8;
  NODES_PER_FACE[ContactSearch::QUADFACEQ9     ] = 9;
  NODES_PER_FACE[ContactSearch::TRIFACEL3      ] = 3;
  NODES_PER_FACE[ContactSearch::TRIFACEQ6      ] = 6;
  NODES_PER_FACE[ContactSearch::SHELLQUADFACEL4] = 4;
  NODES_PER_FACE[ContactSearch::SHELLTRIFACEL3 ] = 3;
  NODES_PER_FACE[ContactSearch::LINEFACEL2     ] = 2;
  NODES_PER_FACE[ContactSearch::LINEFACEQ3     ] = 3;

  EDGES_PER_FACE[ContactSearch::QUADFACEL4     ] = 4;
  EDGES_PER_FACE[ContactSearch::QUADFACEQ8     ] = 4;
  EDGES_PER_FACE[ContactSearch::QUADFACEQ9     ] = 4;
  EDGES_PER_FACE[ContactSearch::TRIFACEL3      ] = 3;
  EDGES_PER_FACE[ContactSearch::TRIFACEQ6      ] = 3;
  EDGES_PER_FACE[ContactSearch::SHELLQUADFACEL4] = 4;
  EDGES_PER_FACE[ContactSearch::SHELLTRIFACEL3 ] = 3;
  EDGES_PER_FACE[ContactSearch::LINEFACEL2     ] = 0;
  EDGES_PER_FACE[ContactSearch::LINEFACEQ3     ] = 0;



  array_init = true;
}


template<typename DataType>
void ContactFace<DataType>::SetNeighborFacesInfo()
{
  number_of_neighbors = 0;
  for (int i=0; i<Edges_Per_Face(); ++i) {
    ContactEdge<DataType>* edge = Edges()[i];
    if (edge->Shared()) {
      ++number_of_neighbors;
    } else {
      if (edge->Neighbor_Face(this)) {
        ++number_of_neighbors;
      }
    }
  }
  if (neighbor_face_info) delete [] neighbor_face_info;
  neighbor_face_info = NULL;
  neighbor_face_info = new typename ContactTopologyEntity<DataType>::connection_data[Edges_Per_Face()];
  
  for (int i=0; i<Edges_Per_Face(); ++i) {
    neighbor_face_info[i].owner = -1; 
    neighbor_face_info[i].owner_proc_array_index = -1;
    neighbor_face_info[i].host_gid[0] = -1; 
    neighbor_face_info[i].host_gid[1] = -1;
    ContactEdge<DataType>* edge = Edges()[i];
    if (edge->Shared()) {
      neighbor_face_info[i] = *(edge->FaceInfo());
    } else {
      ContactFace<DataType>* NeighborFace = edge->Neighbor_Face(this);
      if (NeighborFace) {
        neighbor_face_info[i].owner = NeighborFace->Owner(); 
        neighbor_face_info[i].owner_proc_array_index = NeighborFace->BlockID();
        neighbor_face_info[i].host_gid[0] = NeighborFace->Global_ID().HiInt(); 
        neighbor_face_info[i].host_gid[1] = NeighborFace->Global_ID().LoInt();
      }
    }
  }
}

template<typename DataType>
void ContactFace<DataType>::ConnectNode(const int num, ContactNode<DataType>* node ) {
  PRECONDITION( num < Nodes_Per_Face() );
  PRECONDITION( Nodes() );
  Nodes()[num] = node;
  this->PackConnection(node, &NodeInfo()[num]);
}


template<typename DataType>
void ContactFace<DataType>::ConnectEdge(const int num, ContactEdge<DataType>* edge ) {
  PRECONDITION( num < Edges_Per_Face() );
  PRECONDITION( Edges() );
  Edges()[num] = edge;
}

#ifndef CONTACT_NO_MPI
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
template<typename DataType>
int ContactFace<DataType>::Size_Interactions(int state)
{
  int size = sizeof(int);
  
  ContactInteractionDLL<DataType>* interactions;
  ContactInteractionEntity<DataType>* entity;
  
  interactions = Get_FaceFace_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* i = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
      size += i->Size()+sizeof(int);
    }
  }
  interactions = Get_FaceCoverage_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceCoverageInteraction* i = 
        static_cast<ContactFaceCoverageInteraction*>(entity);
      size += i->Size()+sizeof(int);
    }
  }
  
  return size;
}

template<typename DataType>
void ContactFace<DataType>::Pack_Interactions( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  int num_interactions = Number_Interactions(state);
  
  *i_buf = num_interactions;
  
  ContactInteractionDLL<DataType>* interactions;
  ContactInteractionEntity<DataType>* entity;
  buffer += sizeof(int);
  
  interactions = Get_FaceFace_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* i = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
      int* ibuf = reinterpret_cast<int*>(buffer);
      *ibuf = i->Size();
      i->Pack(buffer+sizeof(int));
      buffer += i->Size()+sizeof(int);
    }
  }
  interactions = Get_FaceCoverage_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceCoverageInteraction* i = 
        static_cast<ContactFaceCoverageInteraction*>(entity);
      int* ibuf = reinterpret_cast<int*>(buffer);
      *ibuf = i->Size();
      i->Pack(buffer+sizeof(int));
      buffer += i->Size()+sizeof(int);
    }
  }
}

template<typename DataType>
void ContactFace<DataType>::Unpack_Interactions( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  int num_interactions = *i_buf;
  char* buf = reinterpret_cast<char*>(buffer+sizeof(int));
  for (int n=0; n<num_interactions; ++n) {
    int* ibuf = reinterpret_cast<int*>(buf);
    int  size = ibuf[0];
    int  interaction_type = ibuf[1];
    buf += sizeof(int);
    switch (interaction_type) {
    case CT_FFI:
      {
	ContactFaceFaceInteraction<DataType>* i = 
          ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction( 
              allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction]);
	i->Unpack( buf );
	i->Connect_SlaveFace( this );
	Store_FaceFace_Interaction( i, state );
      }
      break;
    case CT_FCI:
      {
	ContactFaceCoverageInteraction* i =
          ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction( 
              allocators[ContactSearch::ALLOC_ContactFaceCoverageInteraction]);
	i->Unpack( buf );
	i->Connect_SlaveFace( this );
	Store_FaceCoverage_Interaction( i, state );
      }
      break;
    }
    buf += size;
  }
}

template<typename DataType>
void ContactFace<DataType>::Copy_Interactions( ContactFace<DataType>* src, int state )
{
  ContactInteractionEntity<DataType>* entity;
  if(src->FaceFaceInteractions.size() > state && src->FaceFaceInteractions[state].NumEntities() > 0) {
    ContactInteractionDLL<DataType> &interactions = src->FaceFaceInteractions[state];
    interactions.IteratorStart();
    while ((entity=interactions.IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* new_ffi = 
        ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction( 
            allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction]);
      ContactFaceFaceInteraction<DataType>* old_ffi = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
      new_ffi->Copy( old_ffi );
      new_ffi->Connect_SlaveFace( this );
      Store_FaceFace_Interaction( new_ffi, state );
    }
  }

  if (src->FaceCoverageInteractions.size() > state && src->FaceCoverageInteractions[state].NumEntities()>0) {
    ContactInteractionDLL<Real> &interactions = src->FaceCoverageInteractions[state];
    interactions.IteratorStart();
    while ((entity=interactions.IteratorForward())) {
      ContactFaceCoverageInteraction* new_fci = 
        ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction( 
            allocators[ContactSearch::ALLOC_ContactFaceCoverageInteraction]);
      ContactFaceCoverageInteraction* old_fci = 
        static_cast<ContactFaceCoverageInteraction*>(entity);
      new_fci->Copy( old_fci );
      new_fci->Connect_SlaveFace( this );
      Store_FaceCoverage_Interaction( new_fci, state );
    }
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
template<typename DataType>
int ContactFace<DataType>::Size_Interactions_ForSecondary(int state)
{
  int size = sizeof(int);
  
  ContactInteractionDLL<DataType>* interactions;
  ContactInteractionEntity<DataType>* entity;
  
  interactions = Get_FaceFace_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* i = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
      size += i->Size()+sizeof(int);
    }
  }
  interactions = Get_FaceCoverage_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceCoverageInteraction* i = 
        static_cast<ContactFaceCoverageInteraction*>(entity);
      size += i->Size()+sizeof(int);
    }
  }
  
  return size;
}

template<typename DataType>
void ContactFace<DataType>::Pack_Interactions_ForSecondary( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  int num_interactions = Number_Interactions(state);
  
  *i_buf = num_interactions;
  
  ContactInteractionDLL<DataType>* interactions;
  ContactInteractionEntity<DataType>* entity;
  buffer += sizeof(int);
  
  interactions = Get_FaceFace_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* i = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
      int* ibuf = reinterpret_cast<int*> (buffer);
      *ibuf = i->Size();
      i->Pack(buffer+sizeof(int));
      buffer += i->Size()+sizeof(int);
    }
  }
  interactions = Get_FaceCoverage_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactFaceCoverageInteraction* i = 
        static_cast<ContactFaceCoverageInteraction*>(entity);
      int* ibuf = reinterpret_cast<int*> (buffer);
      *ibuf = i->Size();
      i->Pack(buffer+sizeof(int));
      buffer += i->Size()+sizeof(int);
    }
  }
}

template<typename DataType>
void ContactFace<DataType>::Unpack_Interactions_ForSecondary( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  int num_interactions = *i_buf;
  char* buf = reinterpret_cast<char*> (buffer+sizeof(int));
  for (int n=0; n<num_interactions; ++n) {
    int* ibuf = reinterpret_cast<int*> (buf);
    int  size = ibuf[0];
    int  interaction_type = ibuf[1];
    buf += sizeof(int);
    switch (interaction_type) {
    case CT_FFI:
      {
	ContactFaceFaceInteraction<DataType>* i = 
          ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction( 
              allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction]);
	i->Unpack( buf );
	i->Connect_SlaveFace( this );
	Store_FaceFace_Interaction( i, state );
      }
      break;
    case CT_FCI:
      {
	ContactFaceCoverageInteraction* i =
          ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction( 
              allocators[ContactSearch::ALLOC_ContactFaceCoverageInteraction]);
	i->Unpack( buf );
	i->Connect_SlaveFace( this );
	Store_FaceCoverage_Interaction( i, state );
      }
      break;
    }
    buf += size;
  }
}

template<typename DataType>
void ContactFace<DataType>::Copy_Interactions_ForSecondary( ContactFace<DataType>* src, int state )
{
  ContactInteractionEntity<DataType>* entity;
  if(src->FaceFaceInteractions.size() > state && src->FaceFaceInteractions[state].NumEntities() > 0) {
    ContactInteractionDLL<DataType> &interactions = src->FaceFaceInteractions[state];
    interactions.IteratorStart();
    while ((entity=interactions.IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* new_ffi = 
        ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction( 
            allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction]);
      ContactFaceFaceInteraction<DataType>* old_ffi = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
      new_ffi->Copy( old_ffi );
      new_ffi->Connect_SlaveFace( this );
      Store_FaceFace_Interaction( new_ffi, state );
    }
  }
  if (src->FaceCoverageInteractions.size() > state && src->FaceCoverageInteractions[state].NumEntities()>0) {
    ContactInteractionDLL<Real> &interactions = src->FaceCoverageInteractions[state];
    interactions.IteratorStart();
    while ((entity=interactions.IteratorForward())) {
      ContactFaceCoverageInteraction* new_fci = 
        ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction( 
            allocators[ContactSearch::ALLOC_ContactFaceCoverageInteraction]);
      ContactFaceCoverageInteraction* old_fci = 
        static_cast<ContactFaceCoverageInteraction*>(entity);
      new_fci->Copy( old_fci );
      new_fci->Connect_SlaveFace( this );
      Store_FaceCoverage_Interaction( new_fci, state );
    }
  }
}

template<typename DataType>
DataType ContactFace<DataType>::MaxSize( VariableHandle POSITION )
{
  int  i, j;
  DataType size=0.0, node_positions[8][3];
  
  for (i=0; i<Nodes_Per_Face(); ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  for (i=0,j=1; i<Edges_Per_Face(); ++i,j=(i+1)%Edges_Per_Face()) {
    DataType dx = node_positions[i][0]-node_positions[j][0];
    DataType dy = node_positions[i][1]-node_positions[j][1];
    DataType dz = node_positions[i][2]-node_positions[j][2];
    size = std::max(size,std::sqrt(dx*dx+dy*dy+dz*dz));
  }
  return (size);
}
#endif

template<typename DataType>
ContactFaceFaceInteraction<DataType>* 
ContactFace<DataType>::Get_FaceFace_Interaction(int interaction_number, int state )
{
  ContactInteractionEntity<DataType>* entity=NULL;

  if(FaceFaceInteractions.size() > state) {
    FaceFaceInteractions[state].IteratorStart();
    while ((entity=FaceFaceInteractions[state].IteratorForward())) {
      if (entity->Index()==interaction_number) break;
    }
    ContactFaceFaceInteraction<DataType>* ffi=NULL;
    if (entity) ffi = static_cast<ContactFaceFaceInteraction<DataType>*>(entity);
    return ffi;
  }
  return NULL;
}

template<typename DataType>
void 
ContactFace<DataType>::Store_FaceFace_Interaction( 
             ContactFaceFaceInteraction<DataType>* ffi, int state )
{
  PRECONDITION( ffi );

  if(FaceFaceInteractions.size() <= state) FaceFaceInteractions.resize(state+1);
  int interaction_number = FaceFaceInteractions[state].NumEntities();
  ffi->Index(interaction_number);
  FaceFaceInteractions[state].Append(ffi);
}

template<typename DataType>
void 
ContactFace<DataType>::Delete_FaceFace_Interaction( 
             ContactFaceFaceInteraction<DataType>* ffi, int state )
{
  ContactInteractionEntity<DataType>* link=NULL;

  if(FaceFaceInteractions.size() <= state) return;
  FaceFaceInteractions[state].IteratorStart();
  while ((link=FaceFaceInteractions[state].IteratorForward())) {
    if (link->Index()==ffi->Index()) {
      ContactFaceFaceInteraction<DataType>* cffi = static_cast<ContactFaceFaceInteraction<DataType>*>(link);
      cffi->~ContactFaceFaceInteraction<DataType>();
      allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi);
      FaceFaceInteractions[state].DeletePrev();
      break;
    }
  }
}

template<typename DataType>
void 
ContactFace<DataType>::Store_FaceCoverage_Interaction( 
             ContactFaceCoverageInteraction* fci, int state )
{
  PRECONDITION( fci );
  if(FaceCoverageInteractions.size() <= state) FaceCoverageInteractions.resize(state+1);
  int interaction_number = FaceCoverageInteractions[state].NumEntities();
  fci->Index(interaction_number);
  FaceCoverageInteractions[state].Append(fci);
}

template<typename DataType>
void 
ContactFace<DataType>::Display_FaceFace_Interactions( ContactParOStream& postream, int state )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=4

  if(FaceFaceInteractions.size() <= state) return;
  FaceFaceInteractions[state].IteratorStart();
  while (ContactInteractionEntity<DataType>* entity=FaceFaceInteractions[state].IteratorForward()) {
    int n;
    ContactFaceFaceInteraction<DataType>* cffi =
       static_cast<ContactFaceFaceInteraction<DataType>*> (entity);
    ContactHostGlobalID Master_GID( cffi->MasterFaceEntityData()->host_gid[0], 
                                    cffi->MasterFaceEntityData()->host_gid[1] );
    ContactFaceFaceVertex<DataType>* vertices = cffi->Get_Vertices();
    postream << "Slave Face " << cffi->SlaveFace()->Global_ID() << "\n";
    postream << "Master Face " << Master_GID << "\n";
    postream << "\t\t Num Areas = " << cffi->NumEdges() << "\n";
    postream << "\t\t Slave Areas ---\n";
    for (n=0; n<cffi->NumEdges()+1; ++n) {
      postream << "\t\t\t " << n << ":  "
               << vertices[n].slave_x << " "
               << vertices[n].slave_y << "\n";
    }
    postream << "\t\t Master Areas ---\n";
    for (n=0; n<cffi->NumEdges()+1; ++n) {
      postream << "\t\t\t " << n << ":  "
               << vertices[n].master_x << " "
               << vertices[n].master_y << "\n";
    }
    postream << "\t\t Edge IDs---\n";
    for (n=0; n<cffi->NumEdges(); ++n) {
      postream << "\t\t\t " << n << ":  "
               << vertices[n].slave_edge_id << "\n";
    }
    postream << "\t\t Num Edges = " << cffi->NumEdges() << "\n";
    postream << "\t\t Master Edge Indices---\n";
    for (n=0; n<cffi->NumEdges(); ++n) {
      if (vertices[n].master_edge_flag)
      postream << "\t\t\t " << n << ":  "
               << n+1 << "\n";
    }
  }	
#endif
}

template<typename DataType>
void 
ContactFace<DataType>::Display_FaceCoverage_Interactions( ContactParOStream& postream, int state )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=4
  int j=0;

  if(FaceCoverageInteractions.size() <= state) return;
  FaceCoverageInteractions[state].IteratorStart();
  while (ContactInteractionEntity<Real>* entity=FaceCoverageInteractions[state].IteratorForward()) {
    ContactFaceCoverageInteraction* cfci =
         static_cast<ContactFaceCoverageInteraction*> (entity);
    postream << "Slave Face " << cfci->SlaveFace()->Global_ID() << "\n";
    if (cfci->NumVertices()==0) {
      postream << "\tFully covered\n";
    } else {
      postream << "\tPartially covered, number of areas = "
               << FaceCoverageInteractions[state].NumEntities() << "\n";
      postream << "\tArea " << ++j << " -----\n";
      int k = 0;
      ContactFaceCoverageVertex* ll_node;
      for( ll_node=cfci->Head(); ll_node; ll_node=ll_node->next ){
        postream << "\t  " << ++k << ":  ("
                 << ll_node->slave_x << ", " 
                 << ll_node->slave_y << ")\n";
      }
    }
  }	
#endif
}

template<typename DataType>
void
ContactFace<DataType>::Update_Interactions( )
{
  int i;
  ContactInteractionEntity<DataType>* link;
  
  // State 0 is the "new" state, State 1 is the "n-1" state, etc.
  
  if( FaceFaceInteractions.size() >= NumberOfStates()) {
    ContactInteractionDLL<DataType> temp = FaceFaceInteractions[NumberOfStates()-1];
    for( i=NumberOfStates()-1 ; i>0 ; --i ) {
      FaceFaceInteractions[i] = FaceFaceInteractions[i-1];
    }
    FaceFaceInteractions[0] = temp;
    FaceFaceInteractions[0].IteratorStart();
    while ((link=FaceFaceInteractions[0].IteratorForward())) {
      ContactFaceFaceInteraction<DataType>* cffi = 
        static_cast<ContactFaceFaceInteraction<DataType>*>(link);
      cffi->~ContactFaceFaceInteraction<DataType>();
      allocators[ContactSearch::ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi);
    }
    FaceFaceInteractions[0].Clear();
  }

  if( FaceCoverageInteractions.size() >= NumberOfStates()) {
    ContactInteractionDLL<Real> temp = FaceCoverageInteractions[NumberOfStates()-1];
    for( i=NumberOfStates()-1 ; i>0 ; --i ){
      FaceCoverageInteractions[i] = FaceCoverageInteractions[i-1];
    }
    FaceCoverageInteractions[0] = temp;
    FaceCoverageInteractions[0].IteratorStart();
    while ((link=FaceCoverageInteractions[0].IteratorForward())) {
      ContactFaceCoverageInteraction* cfci =
         static_cast<ContactFaceCoverageInteraction*> (link);
      cfci->~ContactFaceCoverageInteraction();
      allocators[ContactSearch::ALLOC_ContactFaceCoverageInteraction].Delete_Frag(cfci);
    }
    FaceCoverageInteractions[0].Clear();
  }
}

template<typename DataType>
void
ContactFace<DataType>::SetEdgeCurvature(VariableHandle CURVATURE)
{
  DataType* curvature = &DataArray[Edge0_Curvature];
  for (int i=0; i<Edges_Per_Face(); ++i) {
    curvature[i] = *Edges()[i]->Variable(CURVATURE);
  }
}


template<typename DataType>
void
ContactFace<DataType>::SetEdgeCurvature(VariableHandle &CURVATURE,
                              ContactEdge<DataType> *edge)
{
  DataType* curvature = &DataArray[Edge0_Curvature];
  const int num_edge = Edges_Per_Face();
  ContactEdge<DataType> **edges = Edges();
  for (int i=0; i<num_edge; ++i) {
    if(edges[i] == edge) {
      curvature[i] = *(edge->Variable(CURVATURE));
      break;
    }
  }
}



template<typename DataType>
DataType
ContactFace<DataType>::GetEdgeCurvature(int i)
{
  DataType* curvature = &DataArray[Edge0_Curvature];
  PRECONDITION(i>=0 && i<Edges_Per_Face());
  return curvature[i];
}

template<typename DataType>
void 
ContactFace<DataType>::GetEdgeInfo(ContactNode<DataType>* node, 
                         ContactNode<DataType>** edge_nodes, 
                         int* edge_nums)
{
  for (int n=0; n<Nodes_Per_Face(); ++n) {
    if (node==node_list[n]) {
      int np1       = n==Nodes_Per_Face()-1?0:n+1;
      int nm1       = n==0?Nodes_Per_Face()-1:n-1;
      int edge0     = n;
      int edge1     = n==0?Nodes_Per_Face()-1:n-1;
      edge_nums[0]  = edge0;
      edge_nums[1]  = edge1;
      edge_nodes[0] = node_list[np1];
      edge_nodes[1] = node_list[nm1];
    }
  }
}

template<typename DataType>
void
ContactFace<DataType>::SetEdgeSmoothedNormal(VariableHandle SMOOTHED_NORMAL)
{
  DataType* smoothed_normal = &DataArray[NUMBER_SCALAR_VARS+3*Edge0_Smooth_Normal];
  for (int i=0; i<Edges_Per_Face(); ++i) {
    DataType* sn = Edges()[i]->Variable(SMOOTHED_NORMAL);
    smoothed_normal[i*3  ] = sn[0];
    smoothed_normal[i*3+1] = sn[1];
    smoothed_normal[i*3+2] = sn[2];
  }
}

template<typename DataType>
void
ContactFace<DataType>::GetEdgeSmoothedNormal(int i, DataType* smoothed_normal)
{
  PRECONDITION(i>=0 && i<Edges_Per_Face());
  DataType* sn = &DataArray[NUMBER_SCALAR_VARS+3*Edge0_Smooth_Normal+3*i];
  smoothed_normal[0] = sn[0];
  smoothed_normal[1] = sn[1];
  smoothed_normal[2] = sn[2];
}

template<typename DataType>
void 
ContactFace<DataType>::ComputeBoundingBoxForSearch(const int num_configs,
                                         const VariableHandle &NODE_COORD_START,
                                         const VariableHandle &NODE_COORD_END,
                                         const int  auto_tol,
                                         const DataType box_inflation,
                                         const DataType user_tol,
                                         ContactBoundingBox &box_c,
                                         ContactBoundingBox &box_p,
                                         ContactBoundingBox &box_s)
{
  ContactBoundingBox node_box;
  int num_face_nodes = Nodes_Per_Face();

  box_s.Reset();
  if (num_configs>1) {
    ContactNode<DataType>* node = node_list[0];
    box_c.set_point(node->Variable(NODE_COORD_START));
    box_p.set_point(node->Variable(NODE_COORD_END));
    for(int inode=1 ; inode<num_face_nodes ; ++inode ){
      node = node_list[inode];
      box_c.add_point(node->Variable(NODE_COORD_START));
      box_p.add_point(node->Variable(NODE_COORD_END));
    }
    box_s = box_c+box_p;
  } else {
    ContactNode<DataType>* node = node_list[0];
    box_c.set_point(node->Variable(NODE_COORD_START));
    for(int inode=1 ; inode<num_face_nodes ; ++inode ){
      node = node_list[inode];
      box_c.add_point(node->Variable(NODE_COORD_START));
    }
    box_s = box_c;
  }
  if (auto_tol) {
    DataType max_box_dimension = box_s.max_dimension();
    DataType box_tol = std::max(box_inflation*max_box_dimension,user_tol);
    box_s.add_tolerance(box_tol);
  } else {
    //DataType max_box_dimension = box.max_dimension();
    //DataType box_tol = std::max(max_box_dimension,user_tol);
    //box_s.add_tolerance(box_tol);
    box_s.add_tolerance(user_tol);
  }
}

template<typename DataType>
void 
ContactFace<DataType>::ComputeBoundingBoxForSearch(const int num_configs,
                                         const VariableHandle &NODE_COORD_START,
                                         const VariableHandle &NODE_COORD_END,
                                         const int  auto_tol,
                                         const DataType box_inflation,
					 const DataType* max_node_motion,
                                         const DataType max_remaining_gap_mag,
                                         const DataType user_search_tol,
                                         ContactBoundingBox &box_c,
                                         ContactBoundingBox &box_p,
                                         ContactBoundingBox &box_s)
{
  ContactBoundingBox node_box;
  int num_face_nodes = Nodes_Per_Face();

  box_s.Reset();
  if (num_configs>1) {
    ContactNode<DataType>* node = node_list[0];
    box_c.set_point(node->Variable(NODE_COORD_START));
    box_p.set_point(node->Variable(NODE_COORD_END));
    for(int inode=1 ; inode<num_face_nodes ; ++inode ){
      node = node_list[inode];
      box_c.add_point(node->Variable(NODE_COORD_START));
      box_p.add_point(node->Variable(NODE_COORD_END));
    }
    box_s = box_c+box_p;
  } else {
    ContactNode<DataType>* node = node_list[0];
    box_c.set_point(node->Variable(NODE_COORD_START));
    for(int inode=1 ; inode<num_face_nodes ; ++inode ){
      node = node_list[inode];
      box_c.add_point(node->Variable(NODE_COORD_START));
    }
    box_s = box_c;
  }
  if (auto_tol) {
    DataType max_box_dimension = box_s.max_dimension();
    DataType box_tol = std::max(box_inflation*max_box_dimension,user_search_tol);
    box_s.add_tolerance(box_tol);
  } else {
    //DataType max_box_dimension = box.max_dimension();
    //DataType box_tol = std::max(max_box_dimension,user_tol);
    //box.add_tolerance(box_tol);
    DataType user_tol[3];
    for (int i=0; i<3; ++i) {
      user_tol[i] = max_node_motion[i]+max_remaining_gap_mag+user_search_tol;
    }
    box_s.add_tolerance(user_tol);
  }
}

#endif  // #ifdef ContactFace_C_
