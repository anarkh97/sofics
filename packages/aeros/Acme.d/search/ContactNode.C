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


#ifndef ContactNode_C_
#define ContactNode_C_

#include "allocators.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactFace.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include <cstddef>
#include <new>

using namespace std;

template<typename DataType>
ContactNode<DataType>::ContactNode(ContactFixedSizeAllocator* alloc,
                         ContactSearch::ContactNode_Type Type, 
                         int Block_Index, 
                         int Host_Index_in_Block,
			 ContactType cetype ) 
  : ContactTopologyEntity<DataType>( Block_Index, Host_Index_in_Block, DataArray, cetype), 
    allocators(alloc) {
  node_type = Type;
  shell_node_base_id = -1;
  entity_key = -1;
  physical_type = CONTINUUM_NODE;
}

template<typename DataType>
ContactNode<DataType>* ContactNode<DataType>::new_ContactNode( ContactFixedSizeAllocator* alloc,
					   ContactSearch::ContactNode_Type Type,
					   int Block_Index, 
					   int Host_Index_in_Block,
					   ContactType cetype )
{
  return new (alloc[ContactSearch::ALLOC_ContactNode].New_Frag()) 
    ContactNode(alloc, Type, Block_Index, Host_Index_in_Block, cetype );
}

template<typename DataType>
void ContactNode_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactNode<DataType>),
                100,  // block size
                0,    // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactNode allocator" );
}

template<typename DataType>
ContactNode<DataType>::~ContactNode()
{
  for(int i = 0; i < NodeNodeInteractions.size(); ++i) {
    ContactInteractionEntity<Real>* link=NULL;
    NodeNodeInteractions[i].IteratorStart();
    while ((link=NodeNodeInteractions[i].IteratorForward())) {
      ContactNodeNodeInteraction* cnni = static_cast<ContactNodeNodeInteraction*>(link);
      cnni->~ContactNodeNodeInteraction();
      allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction].Delete_Frag(cnni);
    }
  }
  NodeNodeInteractions.clear();


  for(int i = 0; i < NodeEntityInteractions.size(); ++i) {
    Delete_NodeEntity_Interactions(i);
  }
  NodeEntityInteractions.clear();
}

template<typename DataType>
void ContactNode<DataType>::Delete_Face_Connections( )
{
  faces.clear();
}

template<typename DataType>
void ContactNode<DataType>::Connect_Face(ContactFace<DataType>* Face )
{
  faces.push_back(pair<ContactFace<DataType>*,int>(Face,-1));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
template<typename DataType>
int ContactNode<DataType>::Size_Interactions(int state)
{
  int size = 2*sizeof(int);
 
  ContactInteractionDLL<Real>* interactions = Get_NodeNode_Interactions(state);
  if(interactions != NULL) {
    interactions->IteratorStart();
    ContactInteractionEntity<Real>* entity;
    while ((entity=interactions->IteratorForward())) {
      ContactNodeNodeInteraction* i = 
        static_cast<ContactNodeNodeInteraction*>(entity);
      size += i->Size();
    }
  }
  if(NodeEntityInteractions.size() > state) {
    for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
      ContactNodeEntityInteraction* cnei = NodeEntityInteractions[state][i];
      size += cnei->Size();            
    }
  }
  return size;
}

template<typename DataType>
void ContactNode<DataType>::Pack_Interactions( char* Buffer, int state )
{
  int cnt      = 0;
  char* buffer = Buffer;
  int* i_buf   = reinterpret_cast<int*>(buffer);
  i_buf[cnt++] = owner_proc_array_index;
  i_buf[cnt++] = Number_Interactions(state);
  buffer      += cnt*sizeof(int);
  
  if(state < NodeNodeInteractions.size()) {
    ContactInteractionEntity<Real>* entity;
    NodeNodeInteractions[state].IteratorStart();
    while ((entity=NodeNodeInteractions[state].IteratorForward())) {
      ContactNodeNodeInteraction* i = 
        static_cast<ContactNodeNodeInteraction*>(entity);
      buffer = i->Pack(buffer);
    }    
  }
  if(NodeEntityInteractions.size() > state) {
    for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
      ContactNodeEntityInteraction*    cnei = NodeEntityInteractions[state][i];
      buffer = cnei->Pack(buffer);
    }
  }
}

template<>
inline void ContactNode<Real>::Unpack_Interactions( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  
  //status_flag = (ContactTopologyEntity<Real>::ActiveStatus)i_buf[2];
  int num_interactions = i_buf[1];
  char* buf = buffer+2*sizeof(int);
  for (int n=0; n<num_interactions; ++n) {
    int* ibuf = reinterpret_cast<int*>(buf);
    int  interaction_type = ibuf[0];
    switch (interaction_type) {
    case CT_NNI:
      {
	ContactNodeNodeInteraction* i = 
          ContactNodeNodeInteraction::new_ContactNodeNodeInteraction( 
              allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction]);
	buf = i->Unpack( buf );
	i->Connect_SlaveNode( this );
	Add_NodeNode_Interaction( i, state );
      }
      break;
    case CT_NFI:
      {
	ContactNodeFaceInteraction* i = 
          ContactNodeFaceInteraction::new_ContactNodeFaceInteraction( 
              allocators[ContactSearch::ALLOC_ContactNodeFaceInteraction]);
	buf = i->Unpack( buf );
	i->Connect_Node( this );
	Add_NodeEntity_Interaction( i, state );
      }
      break;
    case CT_NSI:
      {
	ContactNodeSurfaceInteraction* i =
          ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction( 
              allocators[ContactSearch::ALLOC_ContactNodeSurfaceInteraction]);
	buf = i->Unpack( buf );
	i->Connect_Node( this );
	Add_NodeEntity_Interaction( i, state );
      }
      break;
    }
  }
}

template<typename DataType>
void ContactNode<DataType>::Unpack_Interactions( char* buffer, int state )
{
  std::cerr << "ContactNode<DataType>::Unpack_Interactions is not implemented\n";
}

template<>
inline void ContactNode<Real>::Copy_Interactions( ContactNode* src, int state )
{
  if(state < src->NodeNodeInteractions.size()) {
    ContactInteractionEntity<Real>* entity;
    src->NodeNodeInteractions[state].IteratorStart();
    while ((entity=src->NodeNodeInteractions[state].IteratorForward())) {
      ContactNodeNodeInteraction* old_nni = 
        static_cast<ContactNodeNodeInteraction*>(entity);
      ContactNodeNodeInteraction* new_nni = 
        ContactNodeNodeInteraction::new_ContactNodeNodeInteraction( 
            allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction]);
      new_nni->Copy( old_nni );
      new_nni->Connect_SlaveNode( this );
      Add_NodeNode_Interaction( new_nni, state );
    }
  }


  if(src->NodeEntityInteractions.size() > state) {
    for(int i = 0; i < src->NodeEntityInteractions[state].size(); ++i) {
      ContactNodeEntityInteraction* new_cnei = src->NodeEntityInteractions[state][i]->New_Instance(allocators);
      new_cnei->Connect_Node( this );
      Add_NodeEntity_Interaction( new_cnei, state );
    }
  }
}

template<typename DataType>
void ContactNode<DataType>::Copy_Interactions( ContactNode* src, int state )
{
  std::cerr << "ContactNode<DataType>::Copy_Interactions is not implemented\n";
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
template<typename DataType>
int ContactNode<DataType>::Size_Interactions_ForSecondary(int state)
{
  int size = sizeof(int);

  if(state < NodeNodeInteractions.size()) {
    ContactNodeNodeInteraction* nni = 
        static_cast<ContactNodeNodeInteraction*>
        (NodeNodeInteractions[state].HeadEntity());
    size += NodeNodeInteractions[state].NumEntities()*nni->Size();
  } 
  if(NodeEntityInteractions.size() > state) {
    for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
      ContactNodeEntityInteraction* cnei = NodeEntityInteractions[state][i];
      size += cnei->Size_ForSecondary();
    }
  }
  return size;
}

template<typename DataType>
void ContactNode<DataType>::Pack_Interactions_ForSecondary( char* Buffer, int state )
{
  char* buffer = Buffer;
  int* i_buf   = reinterpret_cast<int*> (buffer);
  *i_buf       = Number_Interactions(state);
  buffer      += sizeof(int);
  

  if(state < NodeNodeInteractions.size()) {
    ContactInteractionEntity<Real>* entity;
    NodeNodeInteractions[state].IteratorStart();
    while ((entity=NodeNodeInteractions[state].IteratorForward())) {
      ContactNodeNodeInteraction* i = 
        static_cast<ContactNodeNodeInteraction*>(entity);
      buffer = i->Pack(buffer);
    }
  }
  if(NodeEntityInteractions.size() > state) {
    for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
      ContactNodeEntityInteraction*    cnei = NodeEntityInteractions[state][i];
      buffer = cnei->Pack_ForSecondary(buffer);
    }
  }
}

template<>
inline void ContactNode<Real>::Unpack_Interactions_ForSecondary( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  int num_interactions = *i_buf;
  char* buf = buffer+sizeof(int);
  for (int n=0; n<num_interactions; ++n) {
    int* ibuf = reinterpret_cast<int*> (buf);
    int  interaction_type = ibuf[0];
    switch (interaction_type) {
    case CT_NNI:
      {
        ContactNodeNodeInteraction* i = 
          ContactNodeNodeInteraction::new_ContactNodeNodeInteraction( 
              allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction]);
        buf = i->Unpack( buf );
        i->Connect_SlaveNode( this );
        Add_NodeNode_Interaction( i, state );
      }
      break;
    case CT_NFI:
      {
        ContactNodeFaceInteraction* i = 
          ContactNodeFaceInteraction::new_ContactNodeFaceInteraction( 
              allocators[ContactSearch::ALLOC_ContactNodeFaceInteraction]);
        buf = i->Unpack_ForSecondary( buf );
        i->Connect_Node( this );
        Add_NodeEntity_Interaction( i, state );
      }
      break;
    case CT_NSI:
      {
        ContactNodeSurfaceInteraction* i =
          ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction( 
              allocators[ContactSearch::ALLOC_ContactNodeSurfaceInteraction]);
        buf = i->Unpack_ForSecondary( buf );
        i->Connect_Node( this );
        Add_NodeEntity_Interaction( i, state );
      }
      break;
    }
  }
}

template<typename DataType>
void ContactNode<DataType>::Unpack_Interactions_ForSecondary( char* buffer, int state )
{
  std::cerr << "ContactNode<DataType>::Unpack_Interactions_ForSecondary is not implemented\n";
}

template<>
inline void ContactNode<Real>::Copy_Interactions_ForSecondary( ContactNode* src, int state )
{
  if(state < src->NodeNodeInteractions.size()) {
    ContactInteractionEntity<Real>* entity;
    src->NodeNodeInteractions[state].IteratorStart();
    while ((entity=src->NodeNodeInteractions[state].IteratorForward())) {
      ContactNodeNodeInteraction* old_nni = 
        static_cast<ContactNodeNodeInteraction*>(entity);
      ContactNodeNodeInteraction* new_nni = 
        ContactNodeNodeInteraction::new_ContactNodeNodeInteraction( 
            allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction]);
      new_nni->Copy( old_nni );
      new_nni->Connect_SlaveNode( this );
      Add_NodeNode_Interaction( new_nni, state );
    }
  }

  if(src->NodeEntityInteractions.size() > state) {
    for(int i = 0; i < src->NodeEntityInteractions[state].size(); ++i) {
      ContactNodeEntityInteraction* new_cnei =
        src->NodeEntityInteractions[state][i]->New_Instance(allocators);
      new_cnei->Connect_Node( this );
      Add_NodeEntity_Interaction( new_cnei, state );
    }
  }
}

template<typename DataType>
void ContactNode<DataType>::Copy_Interactions_ForSecondary( ContactNode* src, int state )
{
  std::cerr << "ContactNode<DataType>::Copy_Interactions_ForSecondary is not implemented\n";
}

template<typename DataType>
ContactNodeNodeInteraction* 
ContactNode<DataType>::Get_NodeNode_Interaction(int interaction_number, int state )
{
  ContactInteractionEntity<Real>* entity=NULL;

  if(state >= NodeNodeInteractions.size()) return NULL;

  NodeNodeInteractions[state].IteratorStart();
  while ((entity=NodeNodeInteractions[state].IteratorForward())) {
    if (entity->Index()==interaction_number) break;
  }
  ContactNodeNodeInteraction* nni=NULL;
  if (entity) nni = static_cast<ContactNodeNodeInteraction*>(entity);
  return nni;
}

template<typename DataType>
void 
ContactNode<DataType>::Add_NodeNode_Interaction( 
                ContactNodeNodeInteraction* nni, int state )
{
  PRECONDITION( nni );

  if(NodeNodeInteractions.size() <= state) NodeNodeInteractions.resize(state+1);
  int interaction_number = NodeNodeInteractions[state].NumEntities();
  nni->Index(interaction_number);
  NodeNodeInteractions[state].Append(nni);
}

template<typename DataType>
void 
ContactNode<DataType>::Delete_NodeNode_Interaction( 
                ContactNodeNodeInteraction* nni, int state )
{
  ContactInteractionEntity<Real>* link=NULL;

  if(state >= NodeNodeInteractions.size()) return;

  NodeNodeInteractions[state].IteratorStart();
  while ((link=NodeNodeInteractions[state].IteratorForward())) {
    if (link->Index()==nni->Index()) {
      ContactNodeNodeInteraction* cnni = static_cast<ContactNodeNodeInteraction*>(link);
      cnni->~ContactNodeNodeInteraction();
      allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction].Delete_Frag(cnni);
      NodeNodeInteractions[state].DeletePrev();
      break;
    }
  }
}

template<typename DataType>
void 
ContactNode<DataType>::Delete_NodeNode_Interactions( int state )
{
  ContactInteractionEntity<Real>* link=NULL;

  if(state >= NodeNodeInteractions.size()) return;


  NodeNodeInteractions[state].IteratorStart();
  while ((link=NodeNodeInteractions[state].IteratorForward())) {
    ContactNodeNodeInteraction* cnni = static_cast<ContactNodeNodeInteraction*>(link);
    cnni->~ContactNodeNodeInteraction();
    allocators[ContactSearch::ALLOC_ContactNodeNodeInteraction].Delete_Frag(cnni);
  }
  NodeNodeInteractions[state].Clear();
}

template<typename DataType>
void 
ContactNode<DataType>::Display_NodeNode_Interactions( ContactParOStream& postream,
					    int state )
{

  if(state >= NodeNodeInteractions.size()) return;
  ContactInteractionEntity<Real>* link=NULL;
  NodeNodeInteractions[state].IteratorStart();
  while ((link=NodeNodeInteractions[state].IteratorForward())) {
    ContactNodeNodeInteraction* cnni = static_cast<ContactNodeNodeInteraction*>(link);
    ContactHostGlobalID Node_GID( cnni->MasterNodeEntityData()->host_gid[0], 
                                  cnni->MasterNodeEntityData()->host_gid[1] );
    postream << "Constraint for Node " << cnni->SlaveNode()->Global_ID()
             << " and Node " << Node_GID << "\n";
    postream << "\t\t Slave  node owner = "
             <<cnni->SlaveNode()->Owner()<<"\n";
    postream << "\t\t Master node owner = "
             <<cnni->MasterNodeEntityData()->owner<<"\n";
    postream << "\t\t Distance = " 
             << cnni->Scalar_Var(ContactNodeNodeInteraction::DISTANCE) 
             << "\n";
  }
}

template<typename DataType>
void
ContactNode<DataType>::Add_NodeEntity_Interaction(ContactNodeEntityInteraction* nei,
				        int state )
{
  PRECONDITION( nei );
  PRECONDITION( nei->Node()==this );
  if(state >= NodeEntityInteractions.size()) NodeEntityInteractions.resize(state+1);
  int i = NodeEntityInteractions[state].size();
  nei->Index(i);
  NodeEntityInteractions[state].push_back(nei);
}

template<typename DataType>
void
ContactNode<DataType>::Store_NodeEntity_Interaction(int interaction_number,
					  ContactNodeEntityInteraction* nei,
					  int state )
{
  PRECONDITION( nei );
  PRECONDITION( nei->Node()==this );
  PRECONDITION( interaction_number >= 0  );
  PRECONDITION( state >= 0);

  nei->Index(interaction_number);

  if(state >= NodeEntityInteractions.size()) NodeEntityInteractions.resize(state+1);
  if(interaction_number >= NodeEntityInteractions[state].size()) NodeEntityInteractions[state].resize(interaction_number+1, NULL);

  if( NodeEntityInteractions[state][interaction_number] ){
    NodeEntityInteractions[state][interaction_number]->Delete_Frag(allocators);
  }
  NodeEntityInteractions[state][interaction_number] = nei;
}
  
template<typename DataType>
void 
ContactNode<DataType>::Delete_NodeEntity_Interaction(ContactNodeEntityInteraction* nei, 
                                           int state )
{

  if(state >= NodeEntityInteractions.size()) return;

  std::vector< ContactNodeEntityInteraction* >::iterator iter;

  for(iter = NodeEntityInteractions[state].begin(); iter != NodeEntityInteractions[state].end(); ++iter) {
    ContactNodeEntityInteraction* data = *iter;
    if(data == nei) {
      nei->Delete_Frag(allocators);
      NodeEntityInteractions[state].erase(iter);
      break;
    }
  }
  for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
    NodeEntityInteractions[state][i]->Index(i);
  }
}

template<typename DataType>
void 
ContactNode<DataType>::Delete_NodeEntity_Interactions( int state )
{
  if(state >= NodeEntityInteractions.size()) return;
  for (int i=0; i<NodeEntityInteractions[state].size(); ++i) {
    NodeEntityInteractions[state][i]->Delete_Frag(allocators);
  }
  NodeEntityInteractions[state].clear();
}
  
template<typename DataType>
void 
ContactNode<DataType>::Display_NodeEntity_Interactions( ContactParOStream& postream,
					      int state )
{
  if(state >= NodeEntityInteractions.size()) return;
  for (int i=0; i<NodeEntityInteractions[state].size(); ++i) {
    ContactNodeEntityInteraction* cnei = NodeEntityInteractions[state][i];
    postream << "Constraint for Node " << cnei->Node()->Global_ID();

    ContactNodeFaceInteraction*   cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
    if (cnfi!=NULL) {
      ContactHostGlobalID Face_GID( cnfi->FaceEntityData()->host_gid[0], 
                                    cnfi->FaceEntityData()->host_gid[1] );
      postream << " and Face " << Face_GID << "\n";
    }
    else{
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(cnei);
      postream << " and Analytic surface " << cnsi->Surface()->Global_ID()
               << " id = " << cnsi->SurfaceID() << "\n";
    }

    postream << "\t\t Source  = " << cnei->Source_Name().data() << "\n";
    postream << "\t\t Location  = " << cnei->Location() << "\n";
    postream << "\t\t Gap_Rate*dt = " 
             << cnei->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) 
             << "\n";
    postream << "\t\t Old Gap     = "
             << cnei->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD)
             << "\n";
    DataType* pbdir = cnei->Get_Pushback_Dir();
    postream << "\t\t PB Dir  = " << pbdir[0] << " " << pbdir[1] << " "
             << pbdir[2] << "\n";
    DataType* coords = cnei->Get_Coordinates();
    postream << "\t\t Coords  = " << coords[0] << " "
             << coords[1] << " " << coords[2] << "\n";
    DataType* norm = cnei->Get_Normal();
    postream << "\t\t Normal  = " << norm[0] << " "
             << norm[1] << " " << norm[2] << "\n";
    int node_entity_key = cnei->Get_Node_Key();
    postream << "\t\t Node Entity Key = " 
             << node_entity_key << "\n"; 
    DataType* pnorm = cnei->Get_Physical_Face_Normal();
    postream << "\t\t PF Normal = " << pnorm[0] << " " 
             << pnorm[1] << " " << pnorm[2] << "\n";
    DataType& node_area = cnei->Get_Node_Area();
    postream << "\t\t Node Area = " 
             << node_area << "\n"; 
    DataType& time_to_contact = cnei->Get_Time_To_Contact();
    postream << "\t\t Time to Contact = " 
             << time_to_contact << "\n"; 

    if(cnfi != NULL){
      for (int k=0; k<cnfi->NumSharedFaces(); ++k) {
        typename ContactTopologyEntity<DataType>::connection_data *face_info = cnfi->SharedFaceData(k);
        postream<<"\t\t neighbor face "<<k<<" = ("
                  <<face_info->host_gid[0]<<", "
                  <<face_info->host_gid[1]<<") on proc "
                  <<face_info->owner<<"\n";
      }   
    }
  }
}

template<typename DataType>
void 
ContactNode<DataType>::Print_NodeEntity_Interactions( int state )
{

  if(state >= NodeEntityInteractions.size()) return;
  for (int i=0; i<NodeEntityInteractions[state].size(); ++i) {
    ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(NodeEntityInteractions[state][i]);
    if (cnfi!=NULL) {
      ContactHostGlobalID Face_GID( cnfi->FaceEntityData()->host_gid[0], 
                                    cnfi->FaceEntityData()->host_gid[1] );
      std::cout << "Constraint for Node " << cnfi->Node()->Global_ID()
           << " and Face " << Face_GID << "\n";
      std::cout << "\t\t Source  = " << cnfi->Source_Name().data() << "\n";
      std::cout << "\t\t Gap_Rate*dt = " 
           << cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) 
           << "\n";
      std::cout << "\t\t Old Gap     = "
           << cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD)
           << "\n";
      DataType* pbdir = cnfi->Get_Pushback_Dir();
      std::cout << "\t\t PB Dir  = " << pbdir[0] << " " << pbdir[1] << " "
           << pbdir[2] << "\n";
      DataType* coords = 
        cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
      std::cout << "\t\t Coords  = " << coords[0] << " "
           << coords[1] << " " << coords[2] << "\n";
      DataType* norm =
        cnfi->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
      std::cout << "\t\t Normal  = " << norm[0] << " "
           << norm[1] << " " << norm[2] << "\n";
      DataType& node_entity_key = 
        cnfi->Scalar_Var(ContactNodeFaceInteraction::NODE_ENTITY_KEY);
      std::cout << "\t\t Node Entity Key = " 
           << node_entity_key << "\n"; 
      DataType* pnorm =
        cnfi->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
      std::cout << "\t\t PF Normal = " << pnorm[0] << " " 
           << pnorm[1] << " " << pnorm[2] << "\n";
    }
  }
  for (int i=0; i<NodeEntityInteractions[state].size(); ++i) {
    ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(NodeEntityInteractions[state][i]);
    if (cnsi!=NULL) {
      std::cout << "Node " << cnsi->Node()->Global_ID() << " interacts with "
           << "analytic surface " << cnsi->Surface()->Global_ID() 
           << " id = "<<cnsi->SurfaceID()<< "\n";
      std::cout << "\tGap = " 
           << cnsi->Scalar_Var(ContactNodeSurfaceInteraction::GAP_CUR) << "\n";
      DataType* cpoint = 
        cnsi->Vector_Var(ContactNodeSurfaceInteraction::CONTACT_POINT);
      std::cout << "\tContact Point: " << cpoint[0] << " " << cpoint[1] 
           << " " << cpoint[2] << "\n";
      DataType* snorm = 
        cnsi->Get_Normal();
      std::cout << "\tSurface Normal: " << snorm[0] << " " << snorm[1]
           << " " << snorm[2] << "\n";
      DataType* pfnorm = 
        cnsi->Vector_Var(ContactNodeSurfaceInteraction::PHYSICAL_FACE_NORMAL);
      std::cout << "\tPF Normal: " << pfnorm[0] << " " << pfnorm[1]
           << " " << pfnorm[2] << "\n";
    }
  }
}

template<typename DataType>
void
ContactNode<DataType>::Update_Interactions( ) {
  int tot_num_interactions = 0;
  for(int istate = 0; istate < NodeEntityInteractions.size(); ++istate) {
    tot_num_interactions += NodeEntityInteractions[istate].size();
  }
  //
  // State 0 is the "new" state, State 1 is the "n-1" state, etc.
  //
  if(tot_num_interactions == 0) return;
  //
  //  If node has any interactions defined, ensure it has all valid state for those interactions
  //
  if(NodeEntityInteractions.size() < NumberOfStates()) NodeEntityInteractions.resize( NumberOfStates());
  vector<ContactNodeEntityInteraction*> ne_temp = NodeEntityInteractions[NumberOfStates()-1];
  for(int i=NumberOfStates()-1 ; i>0 ; --i ) {
    NodeEntityInteractions[i] = NodeEntityInteractions[i-1];
  }
  NodeEntityInteractions[0] = ne_temp;
  for(int i=0 ; i < NodeEntityInteractions[0].size() ; ++i ){
    ContactNodeEntityInteraction* cnei = NodeEntityInteractions[0][i];
    if( cnei ){
      cnei->Delete_Frag(allocators);
      tot_num_interactions--;
    }
  }
  NodeEntityInteractions[0].clear();

  if(tot_num_interactions == 0) NodeEntityInteractions.clear();
}

template<typename DataType>
int ContactNode<DataType>::Number_NodeFace_Interactions(const int state) const{
  if(state >= NodeEntityInteractions.size()) return 0;
  int count = 0;
  for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
    ContactNodeEntityInteraction *cnei = NodeEntityInteractions[state][i];
    ContactNodeFaceInteraction *cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
    if(cnfi) count++;
  }
  return count;
}

template<typename DataType>
int ContactNode<DataType>::Number_NodeSurface_Interactions(const int state) const{
  if(state >= NodeEntityInteractions.size()) return 0;
  int count = 0;
  for(int i = 0; i < NodeEntityInteractions[state].size(); ++i) {
    ContactNodeEntityInteraction *cnei = NodeEntityInteractions[state][i];
    ContactNodeSurfaceInteraction *cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(cnei);
    if(cnsi) count++;
  }
  return count;
}

template<typename DataType>
int ContactNode<DataType>::Get_Owning_Entity() {
  //
  //  Check if the node has a defined entity key
  //
  if(entity_key >= 0) return entity_key;
  //
  //  If the node is owned by one and only one face entity, return that key
  //  otherwise return a -1
  //
  int num_face_connect = Number_Face_Connections();
  if(num_face_connect == 0) return -1;
  int face_key = GetFace(0)->Entity_Key();
  for(int iface = 1; iface < num_face_connect; ++iface) {
    if(GetFace(iface)->Entity_Key() == face_key) {
      continue;
    } else {
      return -1;
    }
  }  
  return face_key;
}



template<typename DataType>
int ContactNode<DataType>::Num_Tied_Interactions(int istate) {
  if(istate >= NodeEntityInteractions.size()) return 0;
  int count = 0;
  for(int i = 0; i < NodeEntityInteractions[istate].size(); ++i) {
    ContactNodeEntityInteraction *cnei = NodeEntityInteractions[istate][i];
    if(cnei->Is_Tied()) count++;
  }
  return count;
}

template<typename DataType>
int ContactNode<DataType>::Num_Tracked_Interactions(int istate) {
  if(istate >= NodeEntityInteractions.size()) return 0;
  int count = 0;
  for(int i = 0; i < NodeEntityInteractions[istate].size(); ++i) {
    ContactNodeEntityInteraction *cnei = NodeEntityInteractions[istate][i];
    if(cnei->Is_Tracked()) count++;
  }
  return count;
}

//
//  Return the entity key that the faces of a given physical face belong to.
//  Note, could be ambigous, i.e., faces in multiple keys in that case return the first entity key found
//
template<typename DataType>
int ContactNode<DataType>::GetFacePFEntityKey(int physical_face_num) {
  int num_face = faces.size();
  for(int iface = 0; iface < num_face; ++iface) {
    if(faces[iface].second != physical_face_num) continue;
    return faces[iface].first->Entity_Key();
  }
  return -1;
}

template<typename DataType>
bool ContactNode<DataType>::ConnectedToFace(const ContactHostGlobalID &id) {
  for(int iface = 0; iface < faces.size(); ++iface) {
    if(GetFace(iface)->Global_ID() == id) return true;
  }
  return false;
}

template<typename DataType>
bool sort_face_pair_by_id(const pair<ContactFace<DataType>*,int> &d1, const pair<ContactFace<DataType>*,int> &d2) {
  return(d1.first->Global_ID() < d2.first->Global_ID());
}

template<typename DataType>
void ContactNode<DataType>::SortConnectedFaces() {
  sort(faces.begin(), faces.end(), sort_face_pair_by_id<DataType>);
}

#endif  // #ifdef ContactNode_C_
