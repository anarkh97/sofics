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

#ifndef ContactWedgeElementL6_C_
#define ContactWedgeElementL6_C_

#include <algorithm>

#include "allocators.h"
#include "ContactNode.h"
#include "ContactLineEdgeL2.h"
#include "ContactTriFaceL3.h"
#include "ContactQuadFaceL4.h"
#include "ContactWedgeElementL6.h"
#include "ContactFixedSizeAllocator.h"
#include "contact_tolerances.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

template<typename DataType>
ContactWedgeElemL6<DataType>::ContactWedgeElemL6( int Block_Index, 
				        int Host_Index_in_Block, int key ) 
  : ContactElem<DataType>( ContactSearch::WEDGEELEML6,Block_Index,Host_Index_in_Block,key) 
{
  for (int i=0; i<Nodes_Per_Element(); ++i) {
    nodes[i]    = NULL;
    node_ids[i] = -1;
  }
  for (int i=0; i<Edges_Per_Element(); ++i) {
    edges[i]    = NULL;
    edge_ids[i] = -1;
  }
  for (int i=0; i<Faces_Per_Element(); ++i) {
    faces[i]    = NULL;
    face_ids[i] = -1;
  }
}

template<typename DataType>
ContactWedgeElemL6<DataType>* ContactWedgeElemL6<DataType>::new_ContactWedgeElemL6(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key)
{
  return new (alloc.New_Frag())
             ContactWedgeElemL6<DataType>(Block_Index, Host_Index_in_Block, key);
}

template<typename DataType>
void ContactWedgeElemL6_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactWedgeElemL6<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactWedgeElemL6<DataType> allocator" );
}

template<typename DataType>
ContactWedgeElemL6<DataType>::~ContactWedgeElemL6() {}

template<typename DataType>
void ContactWedgeElemL6<DataType>::BuildTopology(int nID, int eID, int fID,
				       ContactFixedSizeAllocator* allocators)
{
  int i;
  ContactNode<DataType>* node;
  ContactEdge<DataType>* edge;
  ContactFace<DataType>* face;
  
  int NextID = nID;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ) {
    node = ContactNode<DataType>::new_ContactNode(allocators,
                                        ContactSearch::NODE, 
                                        ++NextID );
    this->ConnectNode(i, node);
  }
  NextID = eID;
  for( i=0 ; i<Edges_Per_Element() ; ++i ) {
    edge = ContactLineEdgeL2<DataType>::new_ContactLineEdgeL2( 
                        allocators[ContactSearch::ALLOC_ContactLineEdgeL2],
                        ContactSearch::LINEEDGEL2, ++NextID);
    this->ConnectEdge(i, edge);
  }
  edge = Edge(0);
  edge->ConnectNode(0, Node(0));
  edge->ConnectNode(1, Node(1));
  edge = Edge(1);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(4));
  edge = Edge(2);
  edge->ConnectNode(0, Node(4));
  edge->ConnectNode(1, Node(3));
  edge = Edge(3);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(0));
  edge = Edge(4);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(2));
  edge = Edge(5);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(5));
  edge = Edge(6);
  edge->ConnectNode(0, Node(5));
  edge->ConnectNode(1, Node(4));
  edge = Edge(7);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(0));
  edge = Edge(8);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(5));
  
  NextID = fID;
  for( i=0 ; i<3 ; ++i ) {
    face = ContactQuadFaceL4<DataType>::new_ContactQuadFaceL4(allocators, ++NextID );
    this->ConnectFace(i, face);
  }
  for( i=3 ; i<5 ; ++i ) {
    face = ContactTriFaceL3<DataType>::new_ContactTriFaceL3(allocators, ++NextID );
    this->ConnectFace(i, face);
  }
  face = Face(0);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(1));
  face->ConnectNode(2, Node(4));
  face->ConnectNode(3, Node(3));
  face = Face(1);
  face->ConnectNode(0, Node(1));
  face->ConnectNode(1, Node(2));
  face->ConnectNode(2, Node(5));
  face->ConnectNode(3, Node(4));
  face = Face(2);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(3));
  face->ConnectNode(2, Node(5));
  face->ConnectNode(3, Node(2));
  face = Face(3);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(2));
  face->ConnectNode(2, Node(1));
  face = Face(4);
  face->ConnectNode(0, Node(3));
  face->ConnectNode(1, Node(4));
  face->ConnectNode(2, Node(5));
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::DeleteTopology(ContactFixedSizeAllocator* allocators)
{
  int i;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ) {
    ContactNode<DataType>* node = Node(i);
    node->~ContactNode<DataType>();
    allocators[ContactSearch::ALLOC_ContactNode].Delete_Frag(node);
  }
  for( i=0 ; i<Edges_Per_Element() ; ++i ) {
    ContactEdge<DataType>* edge = Edge(i);
    edge->~ContactEdge<DataType>();
    allocators[ContactSearch::ALLOC_ContactLineEdgeL2].Delete_Frag(edge);
  }
  for( i=0 ; i<3 ; ++i ) {
    ContactFace<DataType>* face = Face(i);
    face->~ContactFace<DataType>();
    allocators[ContactSearch::ALLOC_ContactQuadFaceL4].Delete_Frag(face);
  }
  for( i=3 ; i<Faces_Per_Element() ; ++i ) {
    ContactFace<DataType>* face = Face(i);
    face->~ContactFace<DataType>();
    allocators[ContactSearch::ALLOC_ContactTriFaceL3].Delete_Frag(face);
  }
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::UpdateTopology(ContactFace<DataType>* face, 
					VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					VariableHandle NODE_NORMAL,
					Real tol, bool use_node_normals)
{
  int i;
  int num_nodes = face->Nodes_Per_Face();
  for( i=0 ; i<num_nodes ; ++i ){
    DataType* projection;
    ContactNode<DataType>* face_node    = face->Node(i);
    ContactNode<DataType>* elem_node1   = Node(i);
    ContactNode<DataType>* elem_node2   = Node(i+num_nodes);
    DataType* face_node_position  = face_node->Variable(POSITION);
    DataType* elem_node_position1 = elem_node1->Variable(POSITION);
    DataType* elem_node_position2 = elem_node2->Variable(POSITION);
    if (use_node_normals) {
      projection = face_node->Variable(NODE_NORMAL);
    } else {
      projection = face->Variable(FACE_NORMAL);
    }
    for( int k=0 ; k<3 ; ++k ){
      elem_node_position1[k] = face_node_position[k]-projection[k]*tol;
      elem_node_position2[k] = face_node_position[k]+projection[k]*tol;
    }
  }
  for( i=0 ; i<Faces_Per_Element() ; ++i ) {
    Face(i)->Compute_Normal(POSITION,FACE_NORMAL);
  }
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Partial_Face_Normal( int i, VariableHandle CURRENT_POSITION,
                                                                VariableHandle FACE_NORMAL,
                                                                DataType (*face_dface_normal)[3], Real tol,
                                                                DataType (*elem_dface_normal)[3], DataType *dd )
{
  switch(i) {
    case 0 : {
      // compute the derivatives of face's face normal w.r.t face's nodal coordinates
      DataType dface_normal[12][3];
      this->Face(i)->Compute_Partial_Face_Normal(CURRENT_POSITION, dface_normal);

      // compute the derivatives of face's nodal coordinates w.r.t. the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][12] = { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
                             { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
                             { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
                             { 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) {
          dp[k][  j] -= tol*face_dface_normal[k][j];
          dp[k][3+j] -= tol*face_dface_normal[k][j];
          dp[k][6+j] += tol*face_dface_normal[k][j];
          dp[k][9+j] += tol*face_dface_normal[k][j];
        }
      }

      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          elem_dface_normal[j][k] = 0;
          for(int l=0; l<12; ++l)
            elem_dface_normal[j][k] += dp[j][l]*dface_normal[l][k];
        }
      }

      // compute the derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        dd[j] = -(elem_dface_normal[j][0]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[0] + this->Face(i)->Variable(FACE_NORMAL)[0]*dp[j][0]
                 +elem_dface_normal[j][1]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[1] + this->Face(i)->Variable(FACE_NORMAL)[1]*dp[j][1]
                 +elem_dface_normal[j][2]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[2] + this->Face(i)->Variable(FACE_NORMAL)[2]*dp[j][2]);
      }
    } break;
    case 1 : {
      // compute the derivatives of face's face normal w.r.t face's nodal coordinates
      DataType dface_normal[12][3];
      this->Face(i)->Compute_Partial_Face_Normal(CURRENT_POSITION, dface_normal);

      // compute the derivatives of face's nodal coordinates w.r.t. the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][12] = { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
                             { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
                             { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
                             { 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 } };
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) {
          dp[k][  j] -= tol*face_dface_normal[k][j];
          dp[k][3+j] -= tol*face_dface_normal[k][j];
          dp[k][6+j] += tol*face_dface_normal[k][j];
          dp[k][9+j] += tol*face_dface_normal[k][j];
        }
      }

      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          elem_dface_normal[j][k] = 0;
          for(int l=0; l<12; ++l)
            elem_dface_normal[j][k] += dp[j][l]*dface_normal[l][k];
        }
      }

      // compute the derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        dd[j] = -(elem_dface_normal[j][0]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[0] + this->Face(i)->Variable(FACE_NORMAL)[0]*dp[j][0]
                 +elem_dface_normal[j][1]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[1] + this->Face(i)->Variable(FACE_NORMAL)[1]*dp[j][1]
                 +elem_dface_normal[j][2]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[2] + this->Face(i)->Variable(FACE_NORMAL)[2]*dp[j][2]);
      }
    } break;
    case 2 : {
      // compute the derivatives of face's face normal w.r.t face's nodal coordinates
      DataType dface_normal[12][3];
      this->Face(i)->Compute_Partial_Face_Normal(CURRENT_POSITION, dface_normal);

      // compute the derivatives of face's nodal coordinates w.r.t. the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][12] = { { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 } };
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) {
          dp[k][  j] -= tol*face_dface_normal[k][j];
          dp[k][3+j] += tol*face_dface_normal[k][j];
          dp[k][6+j] += tol*face_dface_normal[k][j];
          dp[k][9+j] -= tol*face_dface_normal[k][j];
        }
      }

      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          elem_dface_normal[j][k] = 0;
          for(int l=0; l<12; ++l)
            elem_dface_normal[j][k] += dp[j][l]*dface_normal[l][k];
        }
      }

      // compute the derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        dd[j] = -(elem_dface_normal[j][0]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[0] + this->Face(i)->Variable(FACE_NORMAL)[0]*dp[j][0]
                 +elem_dface_normal[j][1]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[1] + this->Face(i)->Variable(FACE_NORMAL)[1]*dp[j][1]
                 +elem_dface_normal[j][2]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[2] + this->Face(i)->Variable(FACE_NORMAL)[2]*dp[j][2]);
      }
    } break;
    case 3 : {
      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        elem_dface_normal[j][0] = -face_dface_normal[j][0];
        elem_dface_normal[j][1] = -face_dface_normal[j][1];
        elem_dface_normal[j][2] = -face_dface_normal[j][2];
      }

      // compute the derivatives of face's first node's coordinates w.r.t. the nodal coordinates of the master
      // contact face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][3] = { { 1, 0, 0 },
                            { 0, 1, 0 },
                            { 0, 0, 1 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 } };
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) {
          dp[k][j] -= tol*face_dface_normal[k][j];
        }
      }

      // compute the derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        dd[j] = -(elem_dface_normal[j][0]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[0] + this->Face(i)->Variable(FACE_NORMAL)[0]*dp[j][0]
                 +elem_dface_normal[j][1]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[1] + this->Face(i)->Variable(FACE_NORMAL)[1]*dp[j][1]
                 +elem_dface_normal[j][2]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[2] + this->Face(i)->Variable(FACE_NORMAL)[2]*dp[j][2]);
      }

    } break;
    case 4 : {
      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        elem_dface_normal[j][0] = face_dface_normal[j][0];
        elem_dface_normal[j][1] = face_dface_normal[j][1];
        elem_dface_normal[j][2] = face_dface_normal[j][2];
      }

      // compute the derivatives of face's first node's coordinates w.r.t. the nodal coordinates of the master
      // contact face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][3] = { { 1, 0, 0 },
                            { 0, 1, 0 },
                            { 0, 0, 1 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 },
                            { 0, 0, 0 } };
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) {
          dp[k][j] += tol*face_dface_normal[k][j];
        }
      }

      // compute the derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        dd[j] = -(elem_dface_normal[j][0]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[0] + this->Face(i)->Variable(FACE_NORMAL)[0]*dp[j][0]
                 +elem_dface_normal[j][1]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[1] + this->Face(i)->Variable(FACE_NORMAL)[1]*dp[j][1]
                 +elem_dface_normal[j][2]*this->Face(i)->Node(0)->Variable(CURRENT_POSITION)[2] + this->Face(i)->Variable(FACE_NORMAL)[2]*dp[j][2]);
      }

    } break;
  }
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Second_Partial_Face_Normal( int i, VariableHandle CURRENT_POSITION,
                                                                       VariableHandle FACE_NORMAL,
                                                                       DataType (*face_dface_normal)[3],
                                                                       DataType (*face_d2face_normal)[3], Real tol,
                                                                       DataType (*elem_d2face_normal)[3], DataType *d2d )
{
  DataType *node0_position = this->Face(i)->Node(0)->Variable(CURRENT_POSITION);
  DataType *face_normal = this->Face(i)->Variable(FACE_NORMAL);
  switch(i) {
    case 0 : {
      // compute the 1st & 2nd derivatives of face's face normal w.r.t face's nodal coordinates
      DataType dface_normal[12][3];
      this->Face(i)->Compute_Partial_Face_Normal(CURRENT_POSITION, dface_normal);
      DataType d2face_normal[12*12][3];
      this->Face(i)->Compute_Second_Partial_Face_Normal(CURRENT_POSITION, d2face_normal);

      // compute the 1st & 2nd derivatives of face's nodal coordinates w.r.t. the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      // Note: the 2nd derivatives are equal to ±tol*face_d2face_normal
      // XXX the following could be further optimized by taking avantage of the special structure of dp
      DataType dp[9][12] = { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
                             { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
                             { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
                             { 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          DataType t = tol*face_dface_normal[j][k]; 
          dp[j][  k] -= t;
          dp[j][3+k] -= t;
          dp[j][6+k] += t;
          dp[j][9+k] += t;
        }
      }

      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType elem_dface_normal[9][3];
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          elem_dface_normal[j][k] = 0;
          for(int l=0; l<12; ++l)
            elem_dface_normal[j][k] += dp[j][l]*dface_normal[l][k];
        }
      }

      for(int j=0; j<9; ++j)
        for(int k=j; k<9; ++k) {
          elem_d2face_normal[9*j+k][0] = 0;
          elem_d2face_normal[9*j+k][1] = 0;
          elem_d2face_normal[9*j+k][2] = 0;
          for(int l=0; l<3; ++l) {
            DataType t = tol*face_d2face_normal[9*j+k][l];
            elem_d2face_normal[9*j+k][0] += t*(-dface_normal[l][0] - dface_normal[3+l][0] + dface_normal[6+l][0] + dface_normal[9+l][0]);
            elem_d2face_normal[9*j+k][1] += t*(-dface_normal[l][1] - dface_normal[3+l][1] + dface_normal[6+l][1] + dface_normal[9+l][1]); 
            elem_d2face_normal[9*j+k][2] += t*(-dface_normal[l][2] - dface_normal[3+l][2] + dface_normal[6+l][2] + dface_normal[9+l][2]);
          }
        }
  
      for(int j=0; j<9; ++j)
        for(int m=0; m<12; ++m) {
          DataType t[3] = { 0, 0, 0 };
          for(int l=0; l<12; ++l) {
            t[0] += dp[j][l]*d2face_normal[12*l+m][0];
            t[1] += dp[j][l]*d2face_normal[12*l+m][1];
            t[2] += dp[j][l]*d2face_normal[12*l+m][2];
          }
          for(int k=j; k<9; ++k) {
            elem_d2face_normal[9*j+k][0] += t[0]*dp[k][m];
            elem_d2face_normal[9*j+k][1] += t[1]*dp[k][m];
            elem_d2face_normal[9*j+k][2] += t[2]*dp[k][m];
          }
        }

      // compute the second derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        for(int k=j; k<9; ++k) {
          d2d[9*j+k] = -(elem_d2face_normal[9*j+k][0]*node0_position[0] + elem_dface_normal[j][0]*dp[k][0] +
                         elem_dface_normal[k][0]*dp[j][0] - tol*face_normal[0]*face_d2face_normal[9*j+k][0] +
                         elem_d2face_normal[9*j+k][1]*node0_position[1] + elem_dface_normal[j][1]*dp[k][1] +
                         elem_dface_normal[k][1]*dp[j][1] - tol*face_normal[1]*face_d2face_normal[9*j+k][1] +
                         elem_d2face_normal[9*j+k][2]*node0_position[2] + elem_dface_normal[j][2]*dp[k][2] +
                         elem_dface_normal[k][2]*dp[j][2] - tol*face_normal[2]*face_d2face_normal[9*j+k][2]);
        }
      }

    } break;
    case 1 : {
      // compute the derivatives of face's face normal w.r.t face's nodal coordinates
      DataType dface_normal[12][3];
      this->Face(i)->Compute_Partial_Face_Normal(CURRENT_POSITION, dface_normal);
      DataType d2face_normal[12*12][3];
      this->Face(i)->Compute_Second_Partial_Face_Normal(CURRENT_POSITION, d2face_normal);

      // compute the derivatives of face's nodal coordinates w.r.t. the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      // Note: the 2nd derivatives are equal to ±tol*face_d2face_normal
      DataType dp[9][12] = { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
                             { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
                             { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
                             { 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 } };
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          DataType t = tol*face_dface_normal[j][k];
          dp[j][  k] -= t;
          dp[j][3+k] -= t;
          dp[j][6+k] += t;
          dp[j][9+k] += t;
        }
      }

      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType elem_dface_normal[9][3];
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          elem_dface_normal[j][k] = 0;
          for(int l=0; l<12; ++l)
            elem_dface_normal[j][k] += dp[j][l]*dface_normal[l][k];
        }
      }

      for(int j=0; j<9; ++j)
        for(int k=j; k<9; ++k) {
          elem_d2face_normal[9*j+k][0] = 0;
          elem_d2face_normal[9*j+k][1] = 0;
          elem_d2face_normal[9*j+k][2] = 0;
          for(int l=0; l<3; ++l) {
            DataType t = tol*face_d2face_normal[9*j+k][l];
            elem_d2face_normal[9*j+k][0] += t*(-dface_normal[l][0] - dface_normal[3+l][0] + dface_normal[6+l][0] + dface_normal[9+l][0]);
            elem_d2face_normal[9*j+k][1] += t*(-dface_normal[l][1] - dface_normal[3+l][1] + dface_normal[6+l][1] + dface_normal[9+l][1]);
            elem_d2face_normal[9*j+k][2] += t*(-dface_normal[l][2] - dface_normal[3+l][2] + dface_normal[6+l][2] + dface_normal[9+l][2]);
          }
        }

      for(int j=0; j<9; ++j)
        for(int m=0; m<12; ++m) {
          DataType t[3] = { 0, 0, 0 };
          for(int l=0; l<12; ++l) {
            t[0] += dp[j][l]*d2face_normal[12*l+m][0];
            t[1] += dp[j][l]*d2face_normal[12*l+m][1];
            t[2] += dp[j][l]*d2face_normal[12*l+m][2];
          }
          for(int k=j; k<9; ++k) {
            elem_d2face_normal[9*j+k][0] += t[0]*dp[k][m];
            elem_d2face_normal[9*j+k][1] += t[1]*dp[k][m];
            elem_d2face_normal[9*j+k][2] += t[2]*dp[k][m];
          }
        }

      // compute the second derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        for(int k=j; k<9; ++k) {
          d2d[9*j+k] = -(elem_d2face_normal[9*j+k][0]*node0_position[0] + elem_dface_normal[j][0]*dp[k][0] +
                         elem_dface_normal[k][0]*dp[j][0] - tol*face_normal[0]*face_d2face_normal[9*j+k][0] +
                         elem_d2face_normal[9*j+k][1]*node0_position[1] + elem_dface_normal[j][1]*dp[k][1] +
                         elem_dface_normal[k][1]*dp[j][1] - tol*face_normal[1]*face_d2face_normal[9*j+k][1] +
                         elem_d2face_normal[9*j+k][2]*node0_position[2] + elem_dface_normal[j][2]*dp[k][2] +
                         elem_dface_normal[k][2]*dp[j][2] - tol*face_normal[2]*face_d2face_normal[9*j+k][2]);
        }
      }

    } break;
    case 2 : {
      // compute the derivatives of face's face normal w.r.t face's nodal coordinates
      DataType dface_normal[12][3];
      this->Face(i)->Compute_Partial_Face_Normal(CURRENT_POSITION, dface_normal);
      DataType d2face_normal[12*12][3];
      this->Face(i)->Compute_Second_Partial_Face_Normal(CURRENT_POSITION, d2face_normal);

      // compute the derivatives of face's nodal coordinates w.r.t. the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      // Note: the 2nd derivatives are equal to ±tol*face_d2face_normal
      DataType dp[9][12] = { { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 },
                             { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 } };
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          DataType t = tol*face_dface_normal[j][k];
          dp[j][  k] -= t;
          dp[j][3+k] += t;
          dp[j][6+k] += t;
          dp[j][9+k] -= t;
        }
      }

      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      DataType elem_dface_normal[9][3];
      for(int j=0; j<9; ++j) {
        for(int k=0; k<3; ++k) {
          elem_dface_normal[j][k] = 0;
          for(int l=0; l<12; ++l)
            elem_dface_normal[j][k] += dp[j][l]*dface_normal[l][k];
        }
      }

      for(int j=0; j<9; ++j)
        for(int k=j; k<9; ++k) {
          elem_d2face_normal[9*j+k][0] = 0;
          elem_d2face_normal[9*j+k][1] = 0;
          elem_d2face_normal[9*j+k][2] = 0;
          for(int l=0; l<3; ++l) {
            DataType t = tol*face_d2face_normal[9*j+k][l];
            elem_d2face_normal[9*j+k][0] += t*(-dface_normal[l][0] + dface_normal[3+l][0] + dface_normal[6+l][0] - dface_normal[9+l][0]);
            elem_d2face_normal[9*j+k][1] += t*(-dface_normal[l][1] + dface_normal[3+l][1] + dface_normal[6+l][1] - dface_normal[9+l][1]);
            elem_d2face_normal[9*j+k][2] += t*(-dface_normal[l][2] + dface_normal[3+l][2] + dface_normal[6+l][2] - dface_normal[9+l][2]);
          }
        }

      for(int j=0; j<9; ++j)
        for(int m=0; m<12; ++m) {
          DataType t[3] = { 0, 0, 0 };
          for(int l=0; l<12; ++l) {
            t[0] += dp[j][l]*d2face_normal[12*l+m][0];
            t[1] += dp[j][l]*d2face_normal[12*l+m][1];
            t[2] += dp[j][l]*d2face_normal[12*l+m][2];
          }
          for(int k=j; k<9; ++k) {
            elem_d2face_normal[9*j+k][0] += t[0]*dp[k][m];
            elem_d2face_normal[9*j+k][1] += t[1]*dp[k][m];
            elem_d2face_normal[9*j+k][2] += t[2]*dp[k][m];
          }
        }

      // compute the second derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<9; ++j) {
        for(int k=j; k<9; ++k) {
          d2d[9*j+k] = -(elem_d2face_normal[9*j+k][0]*node0_position[0] + elem_dface_normal[j][0]*dp[k][0] +
                         elem_dface_normal[k][0]*dp[j][0] - tol*face_normal[0]*face_d2face_normal[9*j+k][0] +
                         elem_d2face_normal[9*j+k][1]*node0_position[1] + elem_dface_normal[j][1]*dp[k][1] +
                         elem_dface_normal[k][1]*dp[j][1] - tol*face_normal[1]*face_d2face_normal[9*j+k][1] +
                         elem_d2face_normal[9*j+k][2]*node0_position[2] + elem_dface_normal[j][2]*dp[k][2] +
                         elem_dface_normal[k][2]*dp[j][2] - tol*face_normal[2]*face_d2face_normal[9*j+k][2]);
        }
      }

    } break;
    case 3 : {
      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<81; ++j) {
        elem_d2face_normal[j][0] = -face_d2face_normal[j][0];
        elem_d2face_normal[j][1] = -face_d2face_normal[j][1];
        elem_d2face_normal[j][2] = -face_d2face_normal[j][2];
      }

      // compute the second derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][3];
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) dp[k][j] = -tol*face_dface_normal[k][j];
      }
      dp[0][0] += 1; dp[1][1] += 1; dp[2][2] += 1;

      DataType toto[3] = { node0_position[0]+tol*face_normal[0], node0_position[1]+tol*face_normal[1], node0_position[2]+tol*face_normal[2] };
      for(int j=0; j<9; ++j) {
        for(int k=j; k<9; ++k) {
          d2d[9*j+k] = face_d2face_normal[9*j+k][0]*toto[0] + face_dface_normal[j][0]*dp[k][0] + face_dface_normal[k][0]*dp[j][0] +
                       face_d2face_normal[9*j+k][1]*toto[1] + face_dface_normal[j][1]*dp[k][1] + face_dface_normal[k][1]*dp[j][1] +
                       face_d2face_normal[9*j+k][2]*toto[2] + face_dface_normal[j][2]*dp[k][2] + face_dface_normal[k][2]*dp[j][2];
        }
      }

    } break;
    case 4 : {
      // compute the derivatives of face's face normal w.r.t the nodal coordinates of the master contact
      // face from which this wedge element was extruded in UpdateTopology, assuming use_node_normals was false
      for(int j=0; j<81; ++j) {
        elem_d2face_normal[j][0] = face_d2face_normal[j][0];
        elem_d2face_normal[j][1] = face_d2face_normal[j][1];
        elem_d2face_normal[j][2] = face_d2face_normal[j][2];
      }

      // compute the second derivatives of the dot product of face's face normal and the position vector of face's
      // first node w.r.t the nodal coordinates of the master contact face from which this wedge element was
      // extruded in UpdateTopology, assuming use_node_normals was false
      DataType dp[9][3];
      for(int k=0; k<9; ++k) {
        for(int j=0; j<3; ++j) dp[k][j] = tol*face_dface_normal[k][j];
      }
      dp[0][0] += 1; dp[1][1] += 1; dp[2][2] += 1;

      DataType toto[3] = { node0_position[0]+tol*face_normal[0], node0_position[1]+tol*face_normal[1], node0_position[2]+tol*face_normal[2] };
      for(int j=0; j<9; ++j) {
        for(int k=j; k<9; ++k) {
          d2d[9*j+k] = -(face_d2face_normal[9*j+k][0]*toto[0] + face_dface_normal[j][0]*dp[k][0] + face_dface_normal[k][0]*dp[j][0] +
                         face_d2face_normal[9*j+k][1]*toto[1] + face_dface_normal[j][1]*dp[k][1] + face_dface_normal[k][1]*dp[j][1] +
                         face_d2face_normal[9*j+k][2]*toto[2] + face_dface_normal[j][2]*dp[k][2] + face_dface_normal[k][2]*dp[j][2]);
        }
      }

    } break;
  }
}

template<typename DataType>
bool
ContactWedgeElemL6<DataType>::Is_Local_Coordinates_Inside_Element( DataType* local_coords )
{
  if( local_coords[0] >=  0.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >=  0.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >=  0.0 && local_coords[2] <= 1.0 &&
      local_coords[3] >= -1.0 && local_coords[3] <= 1.0 )
    return true;
  return false;
}

template<typename DataType>
bool
ContactWedgeElemL6<DataType>::Is_Local_Coordinates_Near_Element( DataType* local_coords, DataType tolerance )
{
  DataType low_coord  = -(1.+tolerance);
  DataType high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Evaluate_Shape_Functions( DataType* local_coords,
						   DataType* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Global_Coordinates( VariableHandle POSITION,
						     DataType* local_coords,
						     DataType* global_coords )
{
  DataType node_positions[6][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    DataType* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Local_Coordinates( DataType Config_Param,
						    VariableHandle POSITION0, 
						    VariableHandle POSITION1, 
						    VariableHandle FACE_NORMAL,
						    DataType* global_coords,
						    DataType* local_coords )
{
  int i, j;
  DataType node_positions[6][3];
  if (Config_Param == 0.0) {
    for (i=0; i<Nodes_Per_Element(); ++i) {
      DataType* node_position = Node(i)->Variable(POSITION0);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (i=0; i<Nodes_Per_Element(); ++i) {
      DataType* node_position = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    DataType alpha = 1.0 - Config_Param, beta = Config_Param;
    for (i=0; i<Nodes_Per_Element(); ++i) {
      DataType* node_position0 = Node(i)->Variable(POSITION0);
      DataType* node_position1 = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = alpha*node_position0[j]+beta*node_position1[j];
      }
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

template<typename DataType>
bool ContactWedgeElemL6<DataType>::Compute_Local_Coordinates( VariableHandle POSITION,
						    DataType* global_coords,
						    DataType* local_coords )
{
  int i, j;
  DataType node_positions[6][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  return Compute_Local_Coords(node_positions, global_coords, local_coords);
}

/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/* The following functions are supplied for doing generic computations   */
/* on a linear W6 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary W6 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Shape_Functions( DataType* local_coords,
						  DataType* shape_functions )
{
  shape_functions[0] = 0.50*local_coords[0]*(1.0-local_coords[3]);
  shape_functions[1] = 0.50*local_coords[1]*(1.0-local_coords[3]);
  shape_functions[2] = 0.50*local_coords[2]*(1.0-local_coords[3]);
  shape_functions[3] = 0.50*local_coords[0]*(1.0+local_coords[3]);
  shape_functions[4] = 0.50*local_coords[1]*(1.0+local_coords[3]);
  shape_functions[5] = 0.50*local_coords[2]*(1.0+local_coords[3]);
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Shape_Derivatives( DataType* local_coords,
						    DataType  shape_derivs[3][6] )
{
  shape_derivs[0][0] =  0.50*(1.0-local_coords[3]);
  shape_derivs[0][1] =  0.0;
  shape_derivs[0][2] = -0.50*(1.0-local_coords[3]);
  shape_derivs[0][3] =  0.50*(1.0+local_coords[3]);
  shape_derivs[0][4] =  0.0;
  shape_derivs[0][5] = -0.50*(1.0+local_coords[3]);

  shape_derivs[1][0] =  0.0;
  shape_derivs[1][1] =  0.50*(1.0-local_coords[3]);
  shape_derivs[1][2] = -0.50*(1.0-local_coords[3]);
  shape_derivs[1][3] =  0.0;
  shape_derivs[1][4] =  0.50*(1.0+local_coords[3]);
  shape_derivs[1][5] = -0.50*(1.0+local_coords[3]);

  shape_derivs[2][0] = -0.50*local_coords[0];
  shape_derivs[2][1] = -0.50*local_coords[1];
  shape_derivs[2][2] = -0.50*local_coords[2];
  shape_derivs[2][3] =  0.50*local_coords[0];
  shape_derivs[2][4] =  0.50*local_coords[1];
  shape_derivs[2][5] =  0.50*local_coords[2];
}

template<>
inline bool ContactWedgeElemL6<Real>::Compute_Local_Coords( Real node_positions[6][3], 
					       Real* global_coords,
					       Real* local_coords )
{
  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i, j;
  int  nnodes=6;
 if(spatial_tolerance_pre > 0) {
  for (i=0; i<nnodes; ++i) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real dd = sqrt(dx*dx+dy*dy+dz*dz);
    if (dd == 0 || sqrt(dd) < spatial_tolerance_pre) break;
  }
  switch (i) {
  case 0:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    local_coords[3] = -1.0;
    break;
  case 1:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    local_coords[2] =  0.0;
    local_coords[3] = -1.0;
    break;
  case 2:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  1.0;
    local_coords[3] = -1.0;
    break;
  case 3:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    local_coords[3] =  1.0;
    break;
  case 4:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    local_coords[2] =  0.0;
    local_coords[3] =  1.0;
    break;
  case 5:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  1.0;
    local_coords[3] =  1.0;
    break;
  }
  if (i<nnodes) return true;
 }
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=200;
  bool converged = false;
  Real u, u0=0.0, u1, du;
  Real v, v0=0.0, v1, dv;
  Real w, w0=0.0, w1, dw;
  Real f[3], J[3][3], invJ[3][3];
  Real shape_derivatives[3][6], coords[4];
  Real residualNorm2, initialResidualNorm2;
  while (!converged && iterations<max_iterations) {
    coords[0] = u0;
    coords[1] = v0;
    coords[2] = 1.0-u0-v0;
    coords[3] = w0;
    Compute_Global_Coords( node_positions, coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( coords, shape_derivatives );
    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<6; ++j) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
    
    Real detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
                       J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])+
                       J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]));

    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])*detJ;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2])*detJ;
    invJ[0][2] =  (J[1][2]*J[0][1]-J[0][2]*J[1][1])*detJ;
    invJ[1][0] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2])*detJ;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])*detJ;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2])*detJ;
    invJ[2][0] =  (J[1][0]*J[2][1]-J[2][0]*J[1][1])*detJ;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[2][0]*J[0][1])*detJ;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])*detJ;

    // APPLY NEWTON ALGORITHM

    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    residualNorm2 = u*u + v*v + w*w;
    if(iterations == 0) initialResidualNorm2 = residualNorm2;
    if(residualNorm2 < newton_tolerance*initialResidualNorm2 || residualNorm2 < abs_newton_tolerance) converged = true;
    else {
      u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
      v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
      w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
      du = abs(u1-u0);
      dv = abs(v1-v0);
      dw = abs(w1-w0);
      u0 = u1;
      v0 = v1;
      w0 = w1;
      //if (du<newton_tolerance && dv<newton_tolerance && dw<newton_tolerance) converged = true;
      ++iterations;
    }
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if (!converged) {
    std::cerr << "ContactWedgeElemL6::Compute_Local_Coordinates() did not converge\n"
              << std::endl;
  }
#endif
  POSTCONDITION(converged);
  /*if(!converged) {
    std::cerr << "ContactWedgeElemL6<Real>::Compute_Local_Coords did not converge: initialResidualNorm2 = "
              << initialResidualNorm2 << ", residualNorm2 = " << residualNorm2 << std::endl;
  }*/
  if(spatial_tolerance_post > 0) {
    // If it's close to any of the edges, snap to it
    if (u0<1.0+spatial_tolerance_post) {
      u0 = min(u0, 1.0);
    }
    if (u0>-spatial_tolerance_post) {
      u0 = max(u0, 0.0);
    }
    if (v0<1.0+spatial_tolerance_post) {
      v0 = min(v0, 1.0);
    }
    if (v0>-spatial_tolerance_post) {
      v0 = max(v0, 0.0);
    }
    if (abs(w0)<1.0+spatial_tolerance_post) {
      w0 = min(w0, 1.0);
      w0 = max(w0,-1.0);
    }
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = 1.0-u0-v0;
  local_coords[3] = w0;
  return converged;
}

template<typename ActiveScalar>
bool ContactWedgeElemL6<ActiveScalar>::Compute_Local_Coords( ActiveScalar active_node_positions[6][3], 
					       ActiveScalar* active_global_coords,
					       ActiveScalar* active_local_coords )
{
  int  i, j;
  const int  nnodes=6;

  double node_positions[nnodes][3], global_coords[3], local_coords[4];
  for(i=0; i<3; ++i) {
    global_coords[i] = GetActiveScalarValue(active_global_coords[i]);
    for(j=0; j<nnodes; ++j)
      node_positions[j][i] = GetActiveScalarValue(active_node_positions[j][i]);
  }

  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;
  //
  // 1st check for coincidence with one of the face nodes
  //
 if(spatial_tolerance_pre > 0) {
  for (i=0; i<nnodes; ++i) {
    double dx = node_positions[i][0]-global_coords[0];
    double dy = node_positions[i][1]-global_coords[1];
    double dz = node_positions[i][2]-global_coords[2];
    double dd  = dx*dx+dy*dy+dz*dz;
    if (dd == 0 || sqrt(dd) < spatial_tolerance_pre) break;
  }
  switch (i) {
  case 0:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    local_coords[3] = -1.0;
    break;
  case 1:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    local_coords[2] =  0.0;
    local_coords[3] = -1.0;
    break;
  case 2:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  1.0;
    local_coords[3] = -1.0;
    break;
  case 3:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    local_coords[3] =  1.0;
    break;
  case 4:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    local_coords[2] =  0.0;
    local_coords[3] =  1.0;
    break;
  case 5:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  1.0;
    local_coords[3] =  1.0;
    break;
  }
  if (i<nnodes) {
    active_local_coords[0] = local_coords[0];
    active_local_coords[1] = local_coords[1];
    active_local_coords[2] = local_coords[2];
    active_local_coords[3] = local_coords[3];
    return true;
  }
 }
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=200;
  bool converged = false;
  double u, u0=0.0, u1, du;
  double v, v0=0.0, v1, dv;
  double w, w0=0.0, w1, dw;
  double f[3], J[3][3], invJ[3][3];
  double shape_derivatives[3][6], shape_functions[6];
  double residualNorm2, initialResidualNorm2;
  double u0_copy=0.0, v0_copy=0.0, w0_copy=0.0;
  while (!converged && iterations<max_iterations) {
    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = 1.0-u0-v0;
    local_coords[3] = w0;
    // BUILD JACOBIAN AND INVERT
    shape_derivatives[0][0] =  0.50*(1.0-local_coords[3]);
    shape_derivatives[0][1] =  0.0;
    shape_derivatives[0][2] = -0.50*(1.0-local_coords[3]);
    shape_derivatives[0][3] =  0.50*(1.0+local_coords[3]);
    shape_derivatives[0][4] =  0.0;
    shape_derivatives[0][5] = -0.50*(1.0+local_coords[3]);

    shape_derivatives[1][0] =  0.0;
    shape_derivatives[1][1] =  0.50*(1.0-local_coords[3]);
    shape_derivatives[1][2] = -0.50*(1.0-local_coords[3]);
    shape_derivatives[1][3] =  0.0;
    shape_derivatives[1][4] =  0.50*(1.0+local_coords[3]);
    shape_derivatives[1][5] = -0.50*(1.0+local_coords[3]);

    shape_derivatives[2][0] = -0.50*local_coords[0];
    shape_derivatives[2][1] = -0.50*local_coords[1];
    shape_derivatives[2][2] = -0.50*local_coords[2];
    shape_derivatives[2][3] =  0.50*local_coords[0];
    shape_derivatives[2][4] =  0.50*local_coords[1];
    shape_derivatives[2][5] =  0.50*local_coords[2];

    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<6; ++j) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
  
    double detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
                         J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])+
                         J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]));

    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])*detJ;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2])*detJ;
    invJ[0][2] =  (J[1][2]*J[0][1]-J[0][2]*J[1][1])*detJ;
    invJ[1][0] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2])*detJ;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])*detJ;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2])*detJ;
    invJ[2][0] =  (J[1][0]*J[2][1]-J[2][0]*J[1][1])*detJ;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[2][0]*J[0][1])*detJ;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])*detJ;

    // APPLY NEWTON ALGORITHM
    shape_functions[0] = 0.50*local_coords[0]*(1.0-local_coords[3]);
    shape_functions[1] = 0.50*local_coords[1]*(1.0-local_coords[3]);
    shape_functions[2] = 0.50*local_coords[2]*(1.0-local_coords[3]);
    shape_functions[3] = 0.50*local_coords[0]*(1.0+local_coords[3]);
    shape_functions[4] = 0.50*local_coords[1]*(1.0+local_coords[3]);
    shape_functions[5] = 0.50*local_coords[2]*(1.0+local_coords[3]);
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    for( int i=0 ; i<nnodes ; ++i ){
      for (int j=0; j<3; ++j) {
        f[j] += shape_functions[i]*node_positions[i][j];
      }
    }
    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    residualNorm2 = u*u + v*v + w*w;
    if(iterations == 0) initialResidualNorm2 = residualNorm2;
    if(residualNorm2 < newton_tolerance*initialResidualNorm2 || residualNorm2 < abs_newton_tolerance) converged = true;
    else {
      u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
      v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
      w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
      du = abs(u1-u0);
      dv = abs(v1-v0);
      dw = abs(w1-w0);
      u0_copy=u0; v0_copy=v0; w0_copy=w0;
      u0 = u1;
      v0 = v1;
      w0 = w1;
      //if (du<newton_tolerance && dv<newton_tolerance && dw<newton_tolerance) converged = true;
      ++iterations;
    }
  }
  /*if(!converged) {
    std::cerr << " *** WARNING: ContactWedgeElemL6<ActiveScalar>::Compute_Local_Coords did not converge: initialResidualNorm2 = "
              << initialResidualNorm2 << ", residualNorm2 = " << residualNorm2 << std::endl;
  }*/
  {
    //
    // repeat last newton iteration to get the derivatives
    //
    int  iterations=0;
    int  max_iterations=1;
    bool converged = false;
    ActiveScalar u, u0=local_coords[0], u1, du;
    ActiveScalar v, v0=local_coords[1], v1, dv;
    ActiveScalar w, w0=local_coords[3], w1, dw;
    ActiveScalar f[3], J[3][3], invJ[3][3];
    ActiveScalar shape_derivatives[3][6];
    while (!converged && iterations<max_iterations) {
      active_local_coords[0] = u0;
      active_local_coords[1] = v0;
      active_local_coords[2] = 1.0-u0-v0;
      active_local_coords[3] = w0;
      // BUILD JACOBIAN AND INVERT
      Compute_Shape_Derivatives( active_local_coords, shape_derivatives );
      for (i=0; i<3; ++i) {
        J[0][i] = 0.0;
        J[1][i] = 0.0;
        J[2][i] = 0.0;
        for (j=0; j<6; ++j) {
          J[0][i] += shape_derivatives[i][j]*active_node_positions[j][0];
          J[1][i] += shape_derivatives[i][j]*active_node_positions[j][1];
          J[2][i] += shape_derivatives[i][j]*active_node_positions[j][2];
        }
      }
    
      ActiveScalar detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
                                 J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])+
                                 J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]));

      invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])*detJ;
      invJ[0][1] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2])*detJ;
      invJ[0][2] =  (J[1][2]*J[0][1]-J[0][2]*J[1][1])*detJ;
      invJ[1][0] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2])*detJ;
      invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])*detJ;
      invJ[1][2] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2])*detJ;
      invJ[2][0] =  (J[1][0]*J[2][1]-J[2][0]*J[1][1])*detJ;
      invJ[2][1] = -(J[0][0]*J[2][1]-J[2][0]*J[0][1])*detJ;
      invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])*detJ;

      // APPLY NEWTON ALGORITHM
      Compute_Global_Coords( active_node_positions, active_local_coords, f );
      u  = f[0]-active_global_coords[0];
      v  = f[1]-active_global_coords[1];
      w  = f[2]-active_global_coords[2];
      u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
      v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
      w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
      du = abs(u1-u0);
      dv = abs(v1-v0);
      dw = abs(w1-w0);
      u0 = u1;
      v0 = v1;
      w0 = w1;
      //if (du<newton_tolerance && dv<newton_tolerance && dw<newton_tolerance) converged = true;
      ++iterations;
    }
    if(spatial_tolerance_post > 0) {
      // If it's close to any of the edges, snap to it
      if (u0<1.0+spatial_tolerance_post) {
        u0 = min(u0, 1.0);
      }
      if (u0>-spatial_tolerance_post) {
        u0 = max(u0, 0.0);
      }
      if (v0<1.0+spatial_tolerance_post) {
        v0 = min(v0, 1.0);
      }
      if (v0>-spatial_tolerance_post) {
        v0 = max(v0, 0.0);
      }
      if (abs(w0)<1.0+spatial_tolerance_post) {
        w0 = min(w0, 1.0);
        w0 = max(w0,-1.0);
      }
    }
    active_local_coords[0] = u0;
    active_local_coords[1] = v0;
    active_local_coords[2] = 1.0-u0-v0;
    active_local_coords[3] = w0;
  }
  return converged;
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Compute_Global_Coords( DataType node_positions[6][3],
						DataType local_coords[4],
						DataType global_coords[3] )
{
  DataType N[6];
  int  nnodes=6;
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      global_coords[j] += N[i]*node_positions[i][j];
    }
  }
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Interpolate_Scalar( DataType  local_coords[4],
					     DataType  node_scalars[6],
					     DataType& interpolated_scalar )
{
  DataType N[6];
  int  nnodes=6;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

template<typename DataType>
void ContactWedgeElemL6<DataType>::Interpolate_Vector( DataType local_coords[4],
					     DataType node_vectors[6][3],
					     DataType interpolated_vector[3] )
{
  DataType N[6];
  int  nnodes=6;
  interpolated_vector[0] = 0.0;
  interpolated_vector[1] = 0.0;
  interpolated_vector[2] = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      interpolated_vector[j] += N[i]*node_vectors[i][j];
    }
  }
}

#endif  // #define ContactWedgeElementL6_C_
