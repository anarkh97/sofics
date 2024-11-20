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


#include "ContactSearch.h"
#include "ContactElement.h"
#include "ContactFace.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactElementElementInteraction.h"
#include "contact_assert.h"
#include <cstddef>
#include <iostream>
#include <cstring>
#include <cmath>

ContactElement::ContactElement( ContactFixedSizeAllocator* alloc,
                                ContactSearch::ContactElement_Type Type, 
				int Block_Index, int Host_Index_in_Block, 
				int key)
  : ContactTopologyEntity<Real>( Block_Index, Host_Index_in_Block, DataArray, CT_ELEMENT), allocators(alloc)
{
  type = Type;
  entity_key = key;
  number_of_states = 2;
  ElementElementInteractions = new ContactInteractionDLL<Real>*[number_of_states];
  for( int i=0 ; i<number_of_states ; ++i ){
    ElementElementInteractions[i] = new ContactInteractionDLL<Real>();
  }
}


ContactElement::~ContactElement()
{
  for( int i=0 ; i<number_of_states ; ++i ){
    ContactInteractionEntity<Real>* link=NULL;
    ElementElementInteractions[i]->IteratorStart();
    while ((link=ElementElementInteractions[i]->IteratorForward())) {
      ContactElementElementInteraction* ceei = 
        static_cast<ContactElementElementInteraction*>(link);
      ceei->~ContactElementElementInteraction();
      allocators[ContactSearch::ALLOC_ContactElementElementInteraction].Delete_Frag(ceei);
    }
    delete ElementElementInteractions[i];
  }
  delete [] ElementElementInteractions;
  ElementElementInteractions = NULL;
}

#ifndef CONTACT_NO_MPI
Real ContactElement::MaxSize( VariableHandle POSITION )
{
  int  i, j;
  Real size=0.0;
  for (i=0; i<Nodes_Per_Element()-1; ++i) {
    Real* node_position0 = Node(i)->Variable(POSITION);
    for (j=i+1; j<Nodes_Per_Element(); ++j) {
      Real* node_position1 = Node(j)->Variable(POSITION);
      Real dx = node_position0[0]-node_position1[0];
      Real dy = node_position0[1]-node_position1[1];
      Real dz = node_position0[2]-node_position1[2];
      size    = std::max(size,std::sqrt(dx*dx+dy*dy+dz*dz));
    }
  }
  return (size);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
int ContactElement::Size_Interactions(int state)
{
  int size = sizeof(int);
  
  ContactInteractionDLL<Real>* interactions;
  
  interactions = Get_ElementElement_Interactions(state);
  interactions->IteratorStart();
  while (ContactInteractionEntity<Real>* entity=interactions->IteratorForward()) {
    ContactElementElementInteraction* i = 
      static_cast<ContactElementElementInteraction*>(entity);
    size += i->Size()+sizeof(int);
  }
  
  return size;
}

void ContactElement::Pack_Interactions( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  int num_interactions = Number_Interactions(state);
  
  *i_buf = num_interactions;
  
  buffer += sizeof(int);
  
  ContactInteractionDLL<Real>* interactions = 
    Get_ElementElement_Interactions(state);
  interactions->IteratorStart();
  while (ContactInteractionEntity<Real>* entity=interactions->IteratorForward()) {
    ContactElementElementInteraction* i = 
      static_cast<ContactElementElementInteraction*>(entity);
    int* ibuf = reinterpret_cast<int*>(buffer);
    *ibuf = i->Size();
    i->Pack(buffer+sizeof(int));
    buffer += i->Size()+sizeof(int);
  }
}

void ContactElement::Unpack_Interactions( char* buffer, int state )
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
    case CT_EEI:
      {
	ContactElementElementInteraction* i = 
          ContactElementElementInteraction::new_ContactElementElementInteraction( 
              allocators[ContactSearch::ALLOC_ContactElementElementInteraction]);
	i->Unpack( buf );
	i->Connect_SlaveElement( this );
	Store_ElementElement_Interaction( i, state );
      }
      break;
    }
    buf += size;
  }
}

void ContactElement::Copy_Interactions( ContactElement* src, int state )
{
  ContactInteractionEntity<Real>* entity;
  ContactInteractionDLL<Real>* interactions;
  if (src->ElementElementInteractions[state]->NumEntities()>0) {
    interactions = src->ElementElementInteractions[state];
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactElementElementInteraction* new_eei = 
        ContactElementElementInteraction::new_ContactElementElementInteraction( 
            allocators[ContactSearch::ALLOC_ContactElementElementInteraction]);
      ContactElementElementInteraction* old_eei = 
        static_cast<ContactElementElementInteraction*>(entity);
      new_eei->Copy( old_eei );
      new_eei->Connect_SlaveElement( this );
      Store_ElementElement_Interaction( new_eei, state );
    }
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
int ContactElement::Size_Interactions_ForSecondary(int state)
{
  int size = sizeof(int);
  
  ContactInteractionDLL<Real>* interactions;
  
  interactions = Get_ElementElement_Interactions(state);
  interactions->IteratorStart();
  while (ContactInteractionEntity<Real>* entity=interactions->IteratorForward()) {
    ContactElementElementInteraction* i = 
      static_cast<ContactElementElementInteraction*>(entity);
    size += i->Size()+sizeof(int);
  }
  
  return size;
}

void ContactElement::Pack_Interactions_ForSecondary( char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  int num_interactions = Number_Interactions(state);
  
  *i_buf = num_interactions;
  
  buffer += sizeof(int);
  
  ContactInteractionDLL<Real>* interactions = 
    Get_ElementElement_Interactions(state);
  interactions->IteratorStart();
  while (ContactInteractionEntity<Real>* entity=interactions->IteratorForward()) {
    ContactElementElementInteraction* i = 
      static_cast<ContactElementElementInteraction*>(entity);
    int* ibuf = reinterpret_cast<int*> (buffer);
    *ibuf = i->Size();
    i->Pack(buffer+sizeof(int));
    buffer += i->Size()+sizeof(int);
  }
}

void ContactElement::Unpack_Interactions_ForSecondary( char* buffer, int state )
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
    case CT_EEI:
      {
	ContactElementElementInteraction* i = 
          ContactElementElementInteraction::new_ContactElementElementInteraction( 
              allocators[ContactSearch::ALLOC_ContactElementElementInteraction]);
	i->Unpack( buf );
	i->Connect_SlaveElement( this );
	Store_ElementElement_Interaction( i, state );
      }
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    buf += size;
  }
}

void ContactElement::Copy_Interactions_ForSecondary( ContactElement* src, int state )
{
  ContactInteractionEntity<Real>* entity;
  ContactInteractionDLL<Real>* interactions;
  if (src->ElementElementInteractions[state]->NumEntities()>0) {
    interactions = src->ElementElementInteractions[state];
    interactions->IteratorStart();
    while ((entity=interactions->IteratorForward())) {
      ContactElementElementInteraction* new_eei = 
        ContactElementElementInteraction::new_ContactElementElementInteraction( 
            allocators[ContactSearch::ALLOC_ContactElementElementInteraction]);
      ContactElementElementInteraction* old_eei = 
        static_cast<ContactElementElementInteraction*>(entity);
      new_eei->Copy( old_eei );
      new_eei->Connect_SlaveElement( this );
      Store_ElementElement_Interaction( new_eei, state );
    }
  }
}
#endif

ContactElementElementInteraction* 
ContactElement::Get_ElementElement_Interaction(int interaction_number, int state )
{
  ContactInteractionEntity<Real>* entity=NULL;
  ElementElementInteractions[state]->IteratorStart();
  while ((entity=ElementElementInteractions[state]->IteratorForward())) {
    if (entity->Index()==interaction_number) break;
  }
  ContactElementElementInteraction* eei=NULL;
  if (entity) eei = static_cast<ContactElementElementInteraction*>(entity);
  return eei;
}

void 
ContactElement::Store_ElementElement_Interaction( 
                ContactElementElementInteraction* eei, int state )
{
  PRECONDITION( eei );
  PRECONDITION( state>=0 && state<number_of_states );
  int interaction_number = ElementElementInteractions[state]->NumEntities();
  eei->Index(interaction_number);
  ElementElementInteractions[state]->Append(eei);
}

void 
ContactElement::Delete_ElementElement_Interaction( 
                ContactElementElementInteraction* eei, int state )
{
  ContactInteractionEntity<Real>* link=NULL;
  ElementElementInteractions[state]->IteratorStart();
  while ((link=ElementElementInteractions[state]->IteratorForward())) {
    if (link->Index()==eei->Index()) {
      ContactElementElementInteraction* ceei = static_cast<ContactElementElementInteraction*>(link);
      ceei->~ContactElementElementInteraction();
      allocators[ContactSearch::ALLOC_ContactElementElementInteraction].Delete_Frag(ceei);
      ElementElementInteractions[state]->DeletePrev();
      break;
    }
  }
}

void 
ContactElement::Display_ElementElement_Interactions( ContactParOStream& postream, int state )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=4
  ElementElementInteractions[state]->IteratorStart();
  while (ContactInteractionEntity<Real>* entity=ElementElementInteractions[state]->IteratorForward()) {
    ContactElementElementInteraction* ceei =
       static_cast<ContactElementElementInteraction*> (entity);
   ContactHostGlobalID Slave_GID ( ceei->SlaveElementEntityData()->host_gid[0], 
                                   ceei->SlaveElementEntityData()->host_gid[1] );
   ContactHostGlobalID Master_GID( ceei->MasterElementEntityData()->host_gid[0], 
                                   ceei->MasterElementEntityData()->host_gid[1] );
    postream << "Slave Element: " << Slave_GID << "  proc_index = " 
             << ceei->SlaveElementEntityData()->index_in_owner_proc_array<< "\n";
    postream << "Master Element:" << Master_GID << "  proc_index = " 
             << ceei->MasterElementEntityData()->index_in_owner_proc_array << "\n";
    postream << "\tVolume = " << ceei->Scalar_Var(ContactElementElementInteraction::VOLUME) << "\n";
  }
#endif
}

void
ContactElement::Update_Interactions( )
{
  int i;
  ContactInteractionEntity<Real>* link;
  
  // State 0 is the "new" state, State 1 is the "n-1" state, etc.
  
  ContactInteractionDLL<Real>* temp = ElementElementInteractions[number_of_states-1];
  for( i=number_of_states-1 ; i>0 ; --i ) {
    ElementElementInteractions[i] = ElementElementInteractions[i-1];
  }
  ElementElementInteractions[0] = temp;

  ElementElementInteractions[0]->IteratorStart();
  while ((link=ElementElementInteractions[0]->IteratorForward())) {
    ContactElementElementInteraction* cffi = 
      static_cast<ContactElementElementInteraction*>(link);
    cffi->~ContactElementElementInteraction();
    allocators[ContactSearch::ALLOC_ContactElementElementInteraction].Delete_Frag(cffi);
  }
  ElementElementInteractions[0]->Clear();
  
}

void 
ContactElement::ComputeBoundingBoxForSearch(const int num_configs,
                                            const VariableHandle &NODE_COORD_START,
                                            const VariableHandle &NODE_COORD_END,
                                            const int  auto_tol,
                                            const Real box_inflation,
                                            const Real user_tol,
                                            ContactBoundingBox &box)
{
  ContactBoundingBox node_box;
  int num_element_nodes = Nodes_Per_Element();

  box.Reset();
  if (num_configs>1) {
    for(int inode=0 ; inode<num_element_nodes ; ++inode ){
      ContactNode<Real>* node = Nodes()[inode];
      node_box.set_point(node->Variable(NODE_COORD_START));
      node_box.add_point(node->Variable(NODE_COORD_END));
      box.add_box(node_box);
    }
  } else {
    for(int inode=0 ; inode<num_element_nodes ; ++inode ){
      ContactNode<Real>* node = Nodes()[inode];
      box.add_point(node->Variable(NODE_COORD_START));
    }
  }
  if (auto_tol) {
    Real max_box_dimension = box.max_dimension();
    Real box_tol = std::max(box_inflation*max_box_dimension,user_tol);
    box.add_tolerance(box_tol);
  } else {
    Real max_box_dimension = box.max_dimension();
    Real box_tol = std::max(box_inflation*max_box_dimension,user_tol);
    box.add_tolerance(box_tol);
    //box.add_tolerance(user_tol);
  }
}
