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


#include "ContactEdgeBlock.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactSearch.h"
#include "ContactTopology.h" 
#include "ContactFixedSizeAllocator.h"
#include <iostream>

ContactEdgeBlock::ContactEdgeBlock( ContactSearch::ContactEdge_Type Type,
				    int block_index,
				    int Entity_Key,
				    int Number_of_Edges,
				    int& Next_ID,
				    ContactTopology* Topology)
  : topology(Topology)
{
  PRECONDITION( Number_of_Edges >= 0 );
  PRECONDITION( Entity_Key > 0 );
  
  ContactSearch* search=topology->Search();

  edge_list = new ContactBlockEntityList(search->Get_Allocators(), Number_of_Edges);
  type = Type;
  number_of_edges = Number_of_Edges;
  num_edges_added = Number_of_Edges;
  entity_key = Entity_Key;
  id = block_index;
  
  MPI_Comm comm = search->Get_Comm();
  int my_proc = contact_processor_number(comm);

  switch( type ){
  case ContactSearch::LINEEDGEL2:
    for( int j=0 ; j<number_of_edges ; ++j ) {
      ContactEdge<Real>* edge = ContactLineEdgeL2<Real>::new_ContactLineEdgeL2(
	       search->Get_Allocator(ContactSearch::ALLOC_ContactLineEdgeL2),
	       block_index,Next_ID);
      edge->Global_ID(my_proc,Next_ID);
      edge_list->Insert(edge);
      ++Next_ID;
    }
    break;
  case ContactSearch::LINEEDGEQ3:
    for( int j=0 ; j<number_of_edges ; ++j ) {
      ContactEdge<Real>* edge = ContactLineEdgeQ3::new_ContactLineEdgeQ3(
	       search->Get_Allocator(ContactSearch::ALLOC_ContactLineEdgeQ3),
	       block_index,Next_ID);
      edge->Global_ID(my_proc,Next_ID);
      edge_list->Insert(edge);
      ++Next_ID;
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown edge type in ContactEdgeBlock" << std::endl;
#endif
    POSTCONDITION( 0 );
  }  
}

ContactEdgeBlock::ContactEdgeBlock( ContactSearch::ContactEdge_Type Type,
				    int block_index, ContactTopology* Topology )
  : topology(Topology)
{
  edge_list = new ContactBlockEntityList(Topology->Search()->Get_Allocators());
  type = Type;
  number_of_edges = 0;
  num_edges_added = 0;
  entity_key = block_index;
  id = block_index;
}

ContactEdgeBlock::~ContactEdgeBlock()
{
  Delete_Edges();
  delete edge_list;
}

void ContactEdgeBlock::Delete_Edge_List()
{
  Delete_Edges();
  edge_list->CleanUp();
  number_of_edges = 0;
  num_edges_added = 0;
}

void ContactEdgeBlock::Delete_Edges()
{
  if( edge_list ){
    ContactFixedSizeAllocator* alloc=NULL;
    switch( type ){
    case ContactSearch::LINEEDGEL2:
      alloc = &topology->Search()->
      Get_Allocator(ContactSearch::ALLOC_ContactLineEdgeL2);
      break;
    case ContactSearch::LINEEDGEQ3:
      alloc = &topology->Search()->
      Get_Allocator(ContactSearch::ALLOC_ContactLineEdgeQ3);
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    POSTCONDITION(alloc!=NULL);
    edge_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity = edge_list->IteratorForward()) {
      ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity);
      edge->~ContactEdge<Real>();
      alloc->Delete_Frag(edge);
    }
  }
}

void ContactEdgeBlock::Insert_Edge( ContactEdge<Real>* edge )
{
  PRECONDITION( edge );
  edge_list->Insert(edge);
  number_of_edges = edge_list->NumEntities();
}

#ifndef CONTACT_NO_MPI
void ContactEdgeBlock::Insert_Edge( char* buffer )
{
  PRECONDITION( buffer );
  edge_list->Insert(buffer);
  number_of_edges = edge_list->NumEntities();
}
#endif

void ContactEdgeBlock::Delete_Edge( ContactEdge<Real>* edge )
{
  PRECONDITION( edge );
  edge_list->Delete(edge);
  number_of_edges = edge_list->NumEntities();
}


