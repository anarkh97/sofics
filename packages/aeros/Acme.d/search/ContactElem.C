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


#ifndef ContactElem_C_
#define ContactElem_C_

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

template<typename DataType>
ContactElem<DataType>::ContactElem(ContactSearch::ContactElem_Type Type, 
			 int Block_Index, int Host_Index_in_Block, int key)
  : ContactTopologyEntity<DataType>( Block_Index, Host_Index_in_Block, DataArray, CT_ELEM)
{
  type       = Type;
  entity_key = key;
}

template<typename DataType>
ContactElem<DataType>::~ContactElem()
{}

template<typename DataType>
void ContactElem<DataType>::ConnectNode( int num, ContactNode<DataType>* node )
{
  PRECONDITION( num < Nodes_Per_Element() );
  PRECONDITION( Nodes() );
  Nodes()[num] = node;
}



template<typename DataType>
void ContactElem<DataType>::ConnectEdge( int num, ContactEdge<DataType>* edge )
{
  PRECONDITION( num < Edges_Per_Element() );
  PRECONDITION( Edges() );
  Edges()[num] = edge;
}

template<typename DataType>
void ContactElem<DataType>::ConnectFace( int num, ContactFace<DataType>* face )
{
  PRECONDITION( num < Faces_Per_Element() );
  PRECONDITION( Faces() );
  Faces()[num] = face;
}


template<typename DataType>
int ContactElem<DataType>::Size()
{
  return( ContactTopologyEntity<DataType>::Size(DataArray_Length()) + 
	  2*(Nodes_Per_Element()+
             Edges_Per_Element()+
             Faces_Per_Element())*sizeof(int) );
}

template<typename DataType>
void ContactElem<DataType>::Pack( char* buffer )
{
  int i;
  int* i_buf = reinterpret_cast<int*>(buffer);
  // ContactTopologyEntity<DataType> packs in location 0 as ContactFace<DataType> and here we pack
  // in the derived type in location 1.
  i_buf[1] = type;
  ContactTopologyEntity<DataType>::Pack( buffer, DataArray_Length() );
  // Add the global ids of the nodes
  i_buf = reinterpret_cast<int*>(buffer+ContactTopologyEntity<DataType>::Size(DataArray_Length()));
  int index = 0;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ){
    i_buf[index++] = Node(i)->Global_ID().HiInt();
    i_buf[index++] = Node(i)->Global_ID().LoInt();
  }
  for( i=0 ; i<Edges_Per_Element() ; ++i ){
    i_buf[index++] = Edge(i)->Global_ID().HiInt();
    i_buf[index++] = Edge(i)->Global_ID().LoInt();
  }
  for( i=0 ; i<Faces_Per_Element() ; ++i ){
    i_buf[index++] = Face(i)->Global_ID().HiInt();
    i_buf[index++] = Face(i)->Global_ID().LoInt();
  }
  POSTCONDITION( index == 2*(Nodes_Per_Element()+
                             Edges_Per_Element()+
                             Faces_Per_Element()) );
}

template<typename DataType>
void ContactElem<DataType>::Unpack( char* buffer )
{
  ContactTopologyEntity<DataType>::Unpack( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == type );
  // Store off the global ids of the nodes
  int* i_buf = reinterpret_cast<int*>( buffer + ContactTopologyEntity<DataType>::Size(DataArray_Length()) );
  std::memcpy( Node_Ids(), i_buf, 2*Nodes_Per_Element()*sizeof(int) );
  i_buf += 2*Nodes_Per_Element();
  std::memcpy( Edge_Ids(), i_buf, 2*Edges_Per_Element()*sizeof(int) );
  i_buf += 2*Edges_Per_Element();
  std::memcpy( Face_Ids(), i_buf, 2*Faces_Per_Element()*sizeof(int) );
}

#ifndef CONTACT_NO_MPI
template<typename DataType>
DataType ContactElem<DataType>::MaxSize( VariableHandle POSITION )
{
  int  i, j;
  DataType size=0.0, node_positions[8][3];
  
  for (i=0; i<Nodes_Per_Element(); ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  for (i=0,j=1; i<Edges_Per_Element(); ++i,j=(i+1)%Edges_Per_Element()) {
    DataType dx = node_positions[i][0]-node_positions[j][0];
    DataType dy = node_positions[i][1]-node_positions[j][1];
    DataType dz = node_positions[i][2]-node_positions[j][2];
    size = std::max(size,std::sqrt(dx*dx+dy*dy+dz*dz));
  }
  for (i=0,j=1; i<Faces_Per_Element(); ++i,j=(i+1)%Faces_Per_Element()) {
    DataType dx = node_positions[i][0]-node_positions[j][0];
    DataType dy = node_positions[i][1]-node_positions[j][1];
    DataType dz = node_positions[i][2]-node_positions[j][2];
    size = std::max(size,std::sqrt(dx*dx+dy*dy+dz*dz));
  }
  return (size);
}
#endif

#endif  // #define ContactElem_C_
