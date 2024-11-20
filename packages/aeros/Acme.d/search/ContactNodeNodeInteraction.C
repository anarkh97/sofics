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


#include "ContactTopologyEntityList.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include <cstring>
#include <new>

ContactNodeNodeInteraction::ContactNodeNodeInteraction( 
   ContactNode<Real>* SlaveNode, ContactNode<Real>* MasterNode, 
   int node_entity_key, Real distance)
  : ContactInteractionEntity<Real>(DataArray, CT_NNI)
{
  PRECONDITION( SlaveNode && MasterNode);
  slave_node = SlaveNode;
  master_node = MasterNode;
  Set_SlaveNodeEntityData();
  Set_MasterNodeEntityData();
  Scalar_Var(NODE_ENTITY_KEY) = node_entity_key;
  Scalar_Var(DISTANCE) = distance;
}


ContactNodeNodeInteraction::ContactNodeNodeInteraction()
  : ContactInteractionEntity<Real>(DataArray, CT_NNI)
{
  slave_node  = NULL;
  master_node = NULL;
}

ContactNodeNodeInteraction::ContactNodeNodeInteraction( 
                                       ContactNodeNodeInteraction& nni )
: ContactInteractionEntity<Real>(DataArray, CT_NNI)
{
  std::memcpy( &(DataArray[0]), &(nni.DataArray[0]), 
	  (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real) );
  slave_node               = nni.slave_node;
  master_node              = nni.master_node;
  slave_node_entity_data   = nni.slave_node_entity_data;
  master_node_entity_data  = nni.master_node_entity_data;
}

ContactNodeNodeInteraction*
ContactNodeNodeInteraction::new_ContactNodeNodeInteraction(
                                     ContactFixedSizeAllocator& alloc,
                                     ContactNode<Real>* SlaveNode,
                                     ContactNode<Real>* MasterNode,
				     int node_entity_key,
                                     Real distance )
{
  return new (alloc.New_Frag())
    ContactNodeNodeInteraction( SlaveNode, MasterNode, 
                                node_entity_key, distance );
}                                                             

ContactNodeNodeInteraction*
ContactNodeNodeInteraction::new_ContactNodeNodeInteraction(
				        ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactNodeNodeInteraction( );
}
 
 
void 
ContactNodeNodeInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactNodeNodeInteraction),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQNodeNodeInteraction allocator" );
}


ContactNodeNodeInteraction::~ContactNodeNodeInteraction() {}

int ContactNodeNodeInteraction::Restart_Size()
{
  return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+
         sizeof(entity_data)/sizeof(int);
}

void ContactNodeNodeInteraction::Restart_Pack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += PackEntityData(&slave_node_entity_data, &buf_loc[cnt]);
  cnt += PackEntityData(&master_node_entity_data, &buf_loc[cnt]);

  // Pack the data
  std::memcpy( &buf_loc[cnt], &DataArray, 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );

}

void ContactNodeNodeInteraction::Restart_Unpack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += UnPackEntityData(&slave_node_entity_data, &buf_loc[cnt]);
  cnt += UnPackEntityData(&master_node_entity_data, &buf_loc[cnt]);

  // Unpack the Data Array
  std::memcpy( &DataArray, &buf_loc[cnt], 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );
}

void ContactNodeNodeInteraction::Connect_SlaveNode(ContactTopologyEntityList& hash_table)
{
  slave_node = static_cast<ContactNode<Real> *>(hash_table.Find(&slave_node_entity_data));
  POSTCONDITION( slave_node );
}

void ContactNodeNodeInteraction::Connect_SlaveNode(ContactNode<Real>* Node)
{
  slave_node = Node;
  POSTCONDITION( slave_node );
}

void ContactNodeNodeInteraction::Connect_MasterNode(ContactTopologyEntityList& hash_table)
{
  master_node = static_cast<ContactNode<Real> *>(hash_table.Find(&master_node_entity_data));
  //POSTCONDITION( master_node );
}

void ContactNodeNodeInteraction::Connect_MasterNode(ContactNode<Real>* Node)
{
  master_node = Node;
}

int 
ContactNodeNodeInteraction::Set_SlaveNodeEntityData() 
{
  if (slave_node) {
    SetEntityData(&slave_node_entity_data, slave_node);
    return 1;
  }else{      
    return 0;
  }
}

int 
ContactNodeNodeInteraction::Set_MasterNodeEntityData() 
{
  if (master_node) {
    SetEntityData(&master_node_entity_data, master_node);
    return 1;
  }else{      
    return 0;
  }
}

