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


#include "ContactNodeBlock.h"
#include "ContactNode.h"
#include "ContactShellNode.h"
#include "ContactTopology.h"
#include "ContactFixedSizeAllocator.h"

#include <iostream>

ContactNodeBlock::ContactNodeBlock( ContactSearch::ContactNode_Type Type,
				    int block_index,
				    int num_nodes,
				    int& Next_ID,
                                    int* exo_ids,
                                    int* host_ids,
				    ContactType* entity_type,
				    ContactTopology* Topology)
  : topology(Topology)
{
  PRECONDITION( num_nodes >= 0 );
  PRECONDITION( block_index >= 0 );

  int offset = (topology->Number_of_Element_Blocks()+
                topology->Number_of_Face_Blocks())>0?1:0;
  
  master          = false;
  slave           = false;
  node_list       = new ContactBlockEntityList(Topology->Search()->Get_Allocators(), num_nodes);
  type            = Type;
  id              = block_index;
  entity_key      = topology->Number_of_Element_Blocks()+
                    topology->Number_of_Face_Blocks() + block_index - offset;
  has_attributes  = 0;
  has_radius_attributes = 0;
  has_normal_attributes = 0;
  number_of_nodes = num_nodes;
  num_nodes_added = number_of_nodes;
  rmax = 0.0;
  
  int hostID = 0;
  switch( type ){
  case ContactSearch::NODE:
    for( int j=0 ; j<number_of_nodes ; ++j) {
      ContactNode<Real>* node;
      switch( entity_type[j] ){
      case( CT_SHELL_NODE ) :
        node = ContactShellNode::new_ContactShellNode( 
	  topology->Search()->Get_Allocators(),
	  ContactSearch::NODE,block_index,hostID++);
        //FIX THIS - what are the exoids and host ids for shell nodes?
        node->Exodus_ID(exo_ids[j]);
        node->Global_ID(host_ids[2*j],host_ids[2*j+1]);
        node_list->Insert(node);
        Next_ID++;
	break;
      case( CT_NODE ) :
        node = ContactNode<Real>::new_ContactNode(topology->Search()->Get_Allocators(),
	  ContactSearch::NODE,block_index,hostID++);
        node->Exodus_ID(exo_ids[j]);
        node->Global_ID(host_ids[2*j],host_ids[2*j+1]);
        node_list->Insert(node);
        Next_ID++;
	break;
      default:
        POSTCONDITION(0);
        break;
      }
    }
    break;
  case ContactSearch::POINT:
    for( int j=0 ; j<number_of_nodes ; ++j) {
      ContactNode<Real>* node = ContactNode<Real>::new_ContactNode(topology->Search()->Get_Allocators(),
	 ContactSearch::POINT,block_index,hostID++);
      node->Exodus_ID(exo_ids[j]);
      node->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      node_list->Insert(node);
      Next_ID++;
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown node type in ContactNodeBlock: " << type << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

ContactNodeBlock::ContactNodeBlock( ContactSearch::ContactNode_Type Type,
				    int block_index, ContactTopology* Topology )
  : topology(Topology)
{
  master    = false;
  slave     = false;
  node_list = new ContactBlockEntityList(Topology->Search()->Get_Allocators());
  has_attributes = 0;
  has_radius_attributes = 0;
  has_normal_attributes = 0;
  type = Type;
  id = block_index;
  number_of_nodes = 0;
  num_nodes_added = 0;
  int offset = (topology->Number_of_Element_Blocks()+
                topology->Number_of_Face_Blocks())>0?1:0;
  entity_key      = topology->Number_of_Element_Blocks()+
                    topology->Number_of_Face_Blocks() + block_index - offset;
}


ContactNodeBlock::~ContactNodeBlock()
{
  Delete_Nodes();
  delete node_list;
}

void
ContactNodeBlock::Add_Nodes( int num_nodes, int* host_ids, int* exo_ids,
                             ContactType* entity_type)
{
  PRECONDITION( num_nodes >= 0 );

  int j;
  
  number_of_nodes += num_nodes;
  
  switch( type ){
  case ContactSearch::NODE:
    for( j=0 ; j<num_nodes ; ++j) {
      ContactNode<Real>* node;
      switch( entity_type[j] ){
      case( CT_SHELL_NODE ) :
        node = ContactShellNode::new_ContactShellNode( 
	  topology->Search()->Get_Allocators(),
	  ContactSearch::NODE,id,-1);
        node->Exodus_ID(exo_ids[j]);
        node->Global_ID(host_ids[2*j],host_ids[2*j+1]);
        node_list->Insert(node);
	break;
      case( CT_NODE ) :
        node = ContactNode<Real>::new_ContactNode(topology->Search()->Get_Allocators(),
	  ContactSearch::NODE,id,-1);
        node->Exodus_ID(exo_ids[j]);
        node->Global_ID(host_ids[2*j],host_ids[2*j+1]);
        node_list->Insert(node);
	break;
      default:
        POSTCONDITION(0);
        break;
      }
    }
    break;
  case ContactSearch::POINT:
    for( j=0 ; j<num_nodes ; ++j) {
      ContactNode<Real>* node = ContactNode<Real>::new_ContactNode(topology->Search()->Get_Allocators(),
	 ContactSearch::POINT,id,-1);
      node->Exodus_ID(exo_ids[j]);
      node->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      node_list->Insert(node);
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown node type in ContactNodeBlock: " << type << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

void 
ContactNodeBlock::Delete_Node_List()
{
  Delete_Nodes();
  node_list->CleanUp();
  number_of_nodes = 0;
  num_nodes_added = 0;
}

void
ContactNodeBlock::UpdateNodeBlock( int block_index,
                                   int num_nodes,
				   int& Next_ID,
				   ContactType* entity_type)
{

  // KHB 2/25/03 : This will have to be updated to deal with shell
  //               nodes for DLB & Adaptivity
#ifdef CONTACT_UPDATE_TOPOLOGY
  PRECONDITION( num_nodes >= 0 );

  int j;

  number_of_nodes = num_nodes;
  Allocate_Node_Array( number_of_nodes );
  switch( type ){
  case ContactSearch::NODE:
    for( j=0 ; j<number_of_nodes ; ++j)
      switch( entity_type[j] ){
      case( CT_SHELL_NODE ) :
	nodes[j] = ContactShellNode::new_ContactShellNode( 
	  topology->Search()->
	  Get_Allocator(ContactSearch::ALLOC_ContactShellNode),
	  ContactSearch::NODE,++Next_ID,block_index,j);
	break;
      case( CT_NODE ) :
	nodes[j] = ContactNode<Real>::new_ContactNode(topology->Search()->Get_Allocator(ContactSearch::ALLOC_ContactNode),
	  ContactSearch::NODE,++Next_ID,block_index,j);
	break;
      default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << "Unknown node type in ContactNodeBlock: " << entity_type[j] << std::endl;
#endif
	POSTCONDITION( 0 );
      }
    num_nodes_added = number_of_nodes;
    break;
  case ContactSearch::POINT:
   for( j=0 ; j<number_of_nodes ; ++j)
      nodes[j] = ContactNode<Real>::new_ContactNode(topology->Search()->Get_Allocator(ContactSearch::ALLOC_ContactNode),
	 ContactSearch::POINT,++Next_ID,block_index,j);
    num_nodes_added = number_of_nodes;
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown node type in ContactNodeBlock: " << type << std::endl;
#endif
    POSTCONDITION( 0 );
  }  
#endif
}

void ContactNodeBlock::Delete_Nodes()
{
  if( node_list ){
    ContactFixedSizeAllocator* alloc;
    switch( type ){
    case ContactSearch::NODE:
      node_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity = node_list->IteratorForward()) {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
	switch( node->Base_Type() ){
	case( CT_NODE ) :
	  {
            alloc = &topology->Search()->
	      Get_Allocator(ContactSearch::ALLOC_ContactNode);
	    node->~ContactNode<Real>();
	    alloc->Delete_Frag( node );
	    break;
	  }
	case( CT_SHELL_NODE ) :
	  { 
            alloc = &topology->Search()->
	      Get_Allocator(ContactSearch::ALLOC_ContactShellNode);
	    (static_cast<ContactShellNode*> (node))->~ContactShellNode();
	    alloc->Delete_Frag( node );
	    break;
	  }
        default:
          POSTCONDITION(0);
          break;
        }
      }
      break;
    case ContactSearch::POINT:
      alloc = &topology->Search()->
	Get_Allocator(ContactSearch::ALLOC_ContactNode);
      node_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity = node_list->IteratorForward()) {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
        node->~ContactNode<Real>();
        alloc->Delete_Frag(node);
      }
      break;
    default:
      POSTCONDITION(0);
      break;
    }
  }
}

void ContactNodeBlock::Insert_Node( ContactNode<Real>* node )
{
  PRECONDITION( node );
  node_list->Insert(node);
  number_of_nodes = node_list->NumEntities();
}

#ifndef CONTACT_NO_MPI
void ContactNodeBlock::Insert_Node( char* buffer )
{
  PRECONDITION( buffer );
  node_list->Insert(buffer);
  number_of_nodes = node_list->NumEntities();
}
void ContactNodeBlock::Insert_Node_ForSecondary( char* buffer )
{
  PRECONDITION( buffer );
  node_list->Insert_ForSecondary(buffer);
  number_of_nodes = node_list->NumEntities();
}
#endif

void ContactNodeBlock::Delete_Node( ContactNode<Real>* node )
{
  PRECONDITION( node );
  node_list->Delete(node);
  number_of_nodes = node_list->NumEntities();
}

void ContactNodeBlock::ComputeBoundingBox(int nconfigs,
                                          VariableHandle POSITION1,
                                          VariableHandle POSITION2,
                                          VariableHandle RADIUS,
                                          MPI_Comm& SearchComm)
{
  if( node_list ){
    local_bounding_box.Reset();
    switch( type ){
    case ContactSearch::NODE:
      node_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity = node_list->IteratorForward()) {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
        if (node->temp_tag==1) {
          Real* position = node->Variable(POSITION1);
          local_bounding_box.add_point(position);
          if (nconfigs>1) {
            position = node->Variable(POSITION2);
            local_bounding_box.add_point(position);
          }
        }
      }
      break;
    case ContactSearch::POINT:
      node_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity = node_list->IteratorForward()) {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
        if (node->temp_tag==1) {
          Real  radius   = *node->Variable(RADIUS);
          Real* position = node->Variable(POSITION1);
          local_bounding_box.add_sphere(position, radius);
          if (nconfigs>1) {
            position = node->Variable(POSITION2);
            local_bounding_box.add_sphere(position, radius);
          }
        }
      }
      break;
    default:
      POSTCONDITION(0);
      break;
    }
  }
  global_bounding_box.set_box(local_bounding_box);
  contact_global_boundingbox(global_bounding_box, SearchComm);
}

void ContactNodeBlock::Display(ContactParOStream& os)
{
  node_list->Display(os);
}

