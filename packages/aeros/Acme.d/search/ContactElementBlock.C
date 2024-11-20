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


#include "ContactElementBlock.h"
#include "ContactHexElementL8.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactFixedSizeAllocator.h"

#include <iostream>

ContactElementBlock::ContactElementBlock(
				    ContactSearch::ContactElement_Type Type,
				    int block_index,
				    int Number_of_Elements,
				    int& Next_ID,
                                    int* host_ids,
				    ContactTopology* Topology )
  : topology(Topology)
{
  PRECONDITION( Number_of_Elements >= 0 );
  PRECONDITION( block_index >= 0 );

  ContactSearch* search = topology->Search();
  
  id        = block_index;
  master    = false;
  slave     = false;
  elem_list = new ContactBlockEntityList(Topology->Search()->Get_Allocators(), Number_of_Elements);
  type = Type;
  number_of_elements = Number_of_Elements;
  num_elements_added = Number_of_Elements;
  entity_key = block_index;
  
  int hostID = 0;

  switch( type ){
  case ContactSearch::CARTESIANHEXELEMENTL8:
    for( int j=0 ; j<number_of_elements ; ++j ) {
      ContactElement* element = 
	ContactCartesianHexElementL8::new_ContactCartesianHexElementL8( 
           search->Get_Allocators(),
           block_index, hostID++, entity_key );
      element->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      elem_list->Insert(element);
      ++Next_ID;
    }
    break;
  case ContactSearch::HEXELEMENTL8:
    for( int j=0 ; j<number_of_elements ; ++j ) {
      ContactElement* element = ContactHexElementL8::new_ContactHexElementL8( 
              search->Get_Allocators(),
	      block_index, hostID++, entity_key );
      element->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      elem_list->Insert(element);
      ++Next_ID;
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown element type in ContactElementBlock" << std::endl;
#endif
    POSTCONDITION( 0 );
  }  
}

ContactElementBlock::ContactElementBlock(
			     ContactSearch::ContactElement_Type Type, 
			     int block_index,
			     ContactTopology* Topology )
  : topology(Topology), number_of_elements(0),
    type(Type), entity_key(block_index)
{
  PRECONDITION( block_index >= 0 );
  id        = block_index;
  master    = false;
  slave     = false;
  number_of_elements = 0;
  num_elements_added = 0;
  elem_list = new ContactBlockEntityList(Topology->Search()->Get_Allocators());
}


ContactElementBlock::~ContactElementBlock()
{
  Delete_Elements();
  delete elem_list;
}

void
ContactElementBlock::Add_Elements(int Number_of_Elements, int* host_ids )
{
  PRECONDITION( Number_of_Elements >= 0 );
  
  int j;

  ContactSearch* search = topology->Search();

  switch( type ){
  case ContactSearch::CARTESIANHEXELEMENTL8:
    for( j=0 ; j<Number_of_Elements ; ++j ) {
      ContactElement* element = 
	ContactCartesianHexElementL8::new_ContactCartesianHexElementL8( 
           search->Get_Allocators(), id, -1, entity_key );
      element->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      elem_list->Insert(element);
    }
    break;
  case ContactSearch::HEXELEMENTL8:
    for( j=0 ; j<Number_of_Elements ; ++j ) {
      ContactElement* element = ContactHexElementL8::new_ContactHexElementL8( 
              search->Get_Allocators(), id, -1, entity_key );
      element->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      elem_list->Insert(element);
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown element type in ContactElementBlock" << std::endl;
#endif
    POSTCONDITION( 0 );
  }  
}

void ContactElementBlock::Delete_Element_List()
{
  Delete_Elements();
  elem_list->CleanUp();
  number_of_elements = 0;
}

void ContactElementBlock::Delete_Elements()
{
  if( elem_list ){
    ContactFixedSizeAllocator* alloc=NULL;
    switch (type) {
    case ContactSearch::CARTESIANHEXELEMENTL8: 
      alloc = &topology->Search()->
      Get_Allocator(ContactSearch::ALLOC_ContactCartesianHexElementL8);
      break;
    case ContactSearch::HEXELEMENTL8: 
      alloc = &topology->Search()->
      Get_Allocator(ContactSearch::ALLOC_ContactHexElementL8);
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    POSTCONDITION(alloc!=NULL);
    elem_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity = elem_list->IteratorForward()) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      element->~ContactElement();
      alloc->Delete_Frag(element);
    }
  }
}

void ContactElementBlock::Insert_Element( ContactElement* element )
{
  PRECONDITION( element );
  elem_list->Insert(element);
  number_of_elements = elem_list->NumEntities();
}

#ifndef CONTACT_NO_MPI
void ContactElementBlock::Insert_Element( char* buffer )
{
  PRECONDITION( buffer );
  elem_list->Insert(buffer);
  number_of_elements = elem_list->NumEntities();
}
void ContactElementBlock::Insert_Element_ForSecondary( char* buffer )
{
  PRECONDITION( buffer );
  elem_list->Insert_ForSecondary(buffer);
  number_of_elements = elem_list->NumEntities();
}
#endif

void ContactElementBlock::Delete_Element( ContactElement* element )
{
  PRECONDITION( element );
  elem_list->Delete(element);
  number_of_elements = elem_list->NumEntities();
}

void ContactElementBlock::ComputeBoundingBox(int nconfigs,
                                             VariableHandle POSITION1,
                                             VariableHandle POSITION2,
                                             MPI_Comm& SearchComm)
{
  if( elem_list ){
    local_bounding_box.Reset();
    elem_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity = elem_list->IteratorForward()) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      ContactBoundingBox object_box;
      for (int k=0; k<element->Nodes_Per_Element(); ++k) {
        ContactNode<Real>* node = element->Node(k);
        Real* position = node->Variable(POSITION1);
        object_box.add_point(position);
        if (nconfigs>1) {
          position = node->Variable(POSITION2);
          object_box.add_point(position);
        }
      }
      local_bounding_box.add_box(object_box);
    }
  }
  global_bounding_box.set_box(local_bounding_box);
  contact_global_boundingbox(global_bounding_box, SearchComm);
}
