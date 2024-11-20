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


#include <algorithm>
#include "ContactUtilities.h"

#include "allocators.h"
#include "ContactQuadFaceQ9.h"
#include "ContactFixedSizeAllocator.h"
#include <iostream>
#include "ContactNode.h"
#include <cstdlib>
#include <cmath>
#include <new>

ContactQuadFaceQ9::ContactQuadFaceQ9( ContactFixedSizeAllocator* alloc,
                                      int Block_Index, 
				      int Index_in_Block, int key ) 
  : ContactFace<Real>( alloc, ContactSearch::QUADFACEQ9,
                 Block_Index, Index_in_Block, key, 
                 nodes, edges, Node_Info, Edge_Info)
{}

ContactQuadFaceQ9* ContactQuadFaceQ9::new_ContactQuadFaceQ9(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactQuadFaceQ9].New_Frag())
             ContactQuadFaceQ9(alloc, Block_Index, Index_in_Block, key);
}

void ContactQuadFaceQ9_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactQuadFaceQ9),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQuadFaceQ9 allocator" );
}

ContactQuadFaceQ9::~ContactQuadFaceQ9() {}

void ContactQuadFaceQ9::Compute_Normal(VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL )
{
  int   i;
  Real  shape_derivatives[2][9];
  Real  e1[3], e2[3];
  Real  local_coords[2] = {0.0, 0.0};
  Real  a=0.0, b=0.0, c=0.0;
  Real* node_position;
  Real* face_normal = Variable(FACE_NORMAL);
  
  e1[0] = e1[1] = e1[2] = 0.0;
  e2[0] = e2[1] = e2[2] = 0.0;
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (i=0; i<9; ++i) {
    node_position = Node(i)->Variable(CURRENT_POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  for (i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ      = 1.0/std::sqrt(a*b-c*c);
  face_normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])*detJ;
  face_normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])*detJ;
  face_normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])*detJ;
}

void ContactQuadFaceQ9::Compute_Normal(VariableHandle POSITION,
                                       Real* normal, Real* local_coords )
{
  Real shape_derivatives[2][9];
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<9; ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  Real  a=0.0, b=0.0, c=0.0;
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ = std::sqrt(a*b-c*c);
  normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])/detJ;
  normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])/detJ;
  normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])/detJ;
}

void ContactQuadFaceQ9::Compute_Normal(Real** nodal_positions,
                                       Real* local_coords, Real* normal )
{
  Real shape_derivatives[2][9];
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<9; ++i) {
    Real* node_position = nodal_positions[i];
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  Real  a=0.0, b=0.0, c=0.0;
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ = std::sqrt(a*b-c*c);
  normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])/detJ;
  normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])/detJ;
  normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])/detJ;
}

void ContactQuadFaceQ9::Compute_CharacteristicLength(VariableHandle CURRENT_POSITION,
				                     VariableHandle CHARACTERISTIC_LENGTH )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  ContactNode<Real>* node2 = Node(2);
  ContactNode<Real>* node3 = Node(3);
  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  Real* Position2 = node2->Variable(CURRENT_POSITION);
  Real* Position3 = node3->Variable(CURRENT_POSITION);
  // Compute the vector from node 0 to node 2 & from node 1 to node 3
  Real Vec02[3],Vec13[3];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];
  Vec13[0] = Position3[0] - Position1[0];
  Vec13[1] = Position3[1] - Position1[1];
  Vec13[2] = Position3[2] - Position1[2];
  // Compute the face normal as the cross product of the two vectors
  Real face_normal[3];
  face_normal[0] = Vec02[1]*Vec13[2] - Vec02[2]*Vec13[1];
  face_normal[1] = Vec02[2]*Vec13[0] - Vec02[0]*Vec13[2];
  face_normal[2] = Vec02[0]*Vec13[1] - Vec02[1]*Vec13[0];

  Real* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  *characteristiclength = std::sqrt( face_normal[0]*face_normal[0] + 
		                     face_normal[1]*face_normal[1] +
		                     face_normal[2]*face_normal[2] ) / 2.0;
}

void ContactQuadFaceQ9::Compute_Centroid( VariableHandle CURRENT_POSITION,
					  VariableHandle CENTROID )
{
  // For the time being, compute the centroid as the average of the four nodes
  Real* centroid = Variable(CENTROID);
  Real  local_coords[2] = {0.0, 0.0};

  Compute_Global_Coordinates(CURRENT_POSITION, local_coords, centroid);
}

void ContactQuadFaceQ9::Get_Edge_Nodes( int i, ContactNode<Real>** node )
{
  PRECONDITION( i>=0 && i<4 );
  switch( i ){
  case 0:
    node[0] = Node(0);
    node[1] = Node(1);
    node[2] = Node(4);
    break;
  case 1:
    node[0] = Node(1);
    node[1] = Node(2);
    node[2] = Node(5);
    break;
  case 2:
    node[0] = Node(2);
    node[1] = Node(3);
    node[2] = Node(6);
    break;
  case 3:
    node[0] = Node(3);
    node[1] = Node(0);
    node[2] = Node(7);
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Invalid Edge request in ContactQuadFaceQ9::Get_Face_Nodes" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

int ContactQuadFaceQ9::Get_Edge_Number( ContactNode<Real>** edge_nodes )
{
  PRECONDITION( edge_nodes[0] && edge_nodes[1] && edge_nodes[2] );
  PRECONDITION( edge_nodes[0] != edge_nodes[2] );
  int i1=-1,i2=-1,i3=-1;
  for( int i=0 ; i<Nodes_Per_Face()-1 ; ++i ){
    ContactNode<Real>* node = Node(i);
    if( edge_nodes[0] == node ) i1 = i;
    if( edge_nodes[1] == node ) i2 = i;
    if( edge_nodes[2] == node ) i3 = i;
  }
  PRECONDITION( 0<=i1 && i1<4 );
  PRECONDITION( 0<=i2 && i2<4 );
  PRECONDITION( 4<=i3 && i3<8 );
  switch( i1 ){
  case 0:
    if( i2 == 1 && i3==4 ) return(0);
    if( i2 == 3 && i3==7 ) return(3);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 1:
    if( i2 == 0 && i3==4 ) return(0);
    if( i2 == 2 && i3==5 ) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 2:
    if( i2 == 1 && i3==5 ) return(1);
    if( i2 == 3 && i3==6 ) return(2);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 3:
    if( i2 == 2 && i3==6 ) return(2);
    if( i2 == 0 && i3==7 ) return(3);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  }
  POSTCONDITION( 0 );
  return( -1 );
}

int ContactQuadFaceQ9::Get_Edge_Number( Real* local_coords )
{
  if (local_coords[1]==-1.0 && local_coords[3]==-1.0) return(0);
  if (local_coords[0]== 1.0 && local_coords[2]== 1.0) return(1);
  if (local_coords[1]== 1.0 && local_coords[3]== 1.0) return(2);
  if (local_coords[0]==-1.0 && local_coords[2]==-1.0) return(3);
  return( -1 );
}

void 
ContactQuadFaceQ9::Compute_Edge_Normal( VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					int Edge,
					Real* edge_normal)
{
  int    i;
  Real   shape_derivatives[2][9];
  Real   e1[3], e2[3];
  Real   local_coords[2];
  Real*  node_position;
  Real   a=0.0, b=0.0, c=0.0;
  Real   f_normal[3];
  
  switch (Edge) {
  case 0:
    local_coords[0] = -1.0;
    local_coords[1] =  0.0;
    break;
  case 1:
    local_coords[0] =  0.0;
    local_coords[1] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    break;
  case 3:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    break;
  }
  e1[0] = e1[1] = e1[2] = 0.0;
  e2[0] = e2[1] = e2[2] = 0.0;
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (i=0; i<9; ++i) {
    node_position = Node(i)->Variable(POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  for (i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ   = 1.0/std::sqrt(a*b-c*c);
  f_normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])*detJ;
  f_normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])*detJ;
  f_normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])*detJ;
  Real edge_tangent[3];

  for (i=0; i<3; ++i) {
    switch (Edge) {
    case 0:
      edge_tangent[i] =  e1[i];
      break;
    case 1:
      edge_tangent[i] = e2[i];
      break;
    case 2:
      edge_tangent[i] = -e1[i];
      break;
    case 3:
      edge_tangent[i] = -e2[i];
      break;
    }
  }
  // Compute edge_normal
  acme::Cross(edge_tangent, f_normal, edge_normal);

  // Normalize the edge_normal
  acme::Normalize(edge_normal);
}

void ContactQuadFaceQ9::Get_Close_Edges( Real* local_coords, int& number,
					 int& edge_1_id, int& edge_2_id ){
  Real  xi = local_coords[0];
  Real eta = local_coords[1];

  int edge_case = 0;

  if( eta <= -0.9 ) edge_case += 1;  // edge 0
  if( xi  >=  0.9 ) edge_case += 2;  // edge 1
  if( eta >=  0.9 ) edge_case += 4;  // edge 2
  if( xi  <= -0.9 ) edge_case += 8;  // edge 3

  switch( edge_case ){
  case(0):
    // not near any edge
    number = 0;
    break;
  case(1):
    // near edge 0 only
    number = 1;
    edge_1_id = 0;
    break;
  case(2):
    // near edge 1 only
    number = 1;
    edge_1_id = 1;
    break;
  case(4):
    // Near Edge 2 only
    number = 1; 
    edge_1_id = 2; 
    break;
  case(8):
    // Near Edge 3 only
    number = 1; 
    edge_1_id = 3;
    break;
  case(3):
    // Near Edges 0 and 1
    number = 2; 
    edge_1_id = 0; 
    edge_2_id = 1; 
    break;
  case(6):
    // Near Edges 1 and 2
    number = 2; 
    edge_1_id = 1; 
    edge_2_id = 2; 
    break;
  case(12):
    // Near Edges 2 and 3
    number = 2; 
    edge_1_id = 2; 
    edge_2_id = 3; 
    break;
  case(9):
    // Near Edges 3 and 0
    number = 2; 
    edge_1_id = 3; 
    edge_2_id = 0; 
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: Can't be near more than two edges" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  }
}

bool ContactQuadFaceQ9::Is_Inside_Face( Real* local_coords )
{
  PRECONDITION( 0 );
  return false;
}

ContactFace<Real>* ContactQuadFaceQ9::Neighbor( Real* local_coords )
{
  PRECONDITION( 0 );
  return (ContactFace<Real>*) NULL;
}

void ContactQuadFaceQ9::FacetDecomposition(int& nfacets,
         Real coordinates0[], Real* normals0, VariableHandle POSITION0,
         Real* coordinates1,  Real* normals1, VariableHandle POSITION1,
         Real* coordinates2,  Real* normals2, VariableHandle POSITION2)
{
  int   i, j, n, ii, jj, nn;
  Real  local_coords0[2] = {0.0, 0.0};
  Real  local_coords1[2] = {0.0, 0.0};
  Real  node0[3], node1[3];
  Real  offset0[8] = {-0.5,  0.0,  0.5,  0.5,  0.5,  0.0, -0.5, -0.5};
  Real  offset1[8] = {-0.5, -0.5, -0.5,  0.0,  0.5,  0.5,  0.5,  0.0};
  
  nfacets = 32;
  for (i=0; i<2; ++i) {
    local_coords0[0] = -0.5+1.0*(Real)i;
    for (j=0; j<2; ++j) {
      Real center[3];
      local_coords0[1] = -0.5+1.0*(Real)j;
      Compute_Global_Coordinates( POSITION0, local_coords0, center );
      local_coords1[0] = local_coords0[0]+offset0[0];
      local_coords1[1] = local_coords0[1]+offset1[0];
      Compute_Global_Coordinates( POSITION0, local_coords1, node0 );
      for (n=0; n<8; ++n) {
        for (ii=0, jj=1; ii<8; ++ii, jj=(ii+1)%8) {
          nn = (i*2+j)*8+ii;
          local_coords1[0] = local_coords0[0]+offset0[jj];
          local_coords1[1] = local_coords0[1]+offset1[jj];
          Compute_Global_Coordinates( POSITION0, local_coords1, node1 );
 	  coordinates0[0+9*nn] = node0[0];
	  coordinates0[1+9*nn] = node0[1];
	  coordinates0[2+9*nn] = node0[2];
	  coordinates0[3+9*nn] = node1[0];
	  coordinates0[4+9*nn] = node1[1];
	  coordinates0[5+9*nn] = node1[2];
 	  coordinates0[6+9*nn] = center[0];
	  coordinates0[7+9*nn] = center[1];
	  coordinates0[8+9*nn] = center[2];
          node0[0] = node1[0];
          node0[1] = node1[1];
          node0[2] = node1[2];
	}
      }
    }
  }
  for (i=0; i<nfacets; ++i) {
    // compute vectors from centroid to other  
    // two nodes defining master surface triangle
    Real vx0   = coordinates0[6+9*i] - coordinates0[0+9*i];
    Real vy0   = coordinates0[7+9*i] - coordinates0[1+9*i];
    Real vz0   = coordinates0[8+9*i] - coordinates0[2+9*i];
    Real vx00  = coordinates0[3+9*i] - coordinates0[0+9*i];
    Real vy00  = coordinates0[4+9*i] - coordinates0[1+9*i];
    Real vz00  = coordinates0[5+9*i] - coordinates0[2+9*i];
    // i,j,k components of the normal to the plane
    Real a4i        = -vy0*vz00+vy00*vz0;
    Real a4j        =  vx0*vz00-vx00*vz0;
    Real a4k        = -vx0*vy00+vx00*vy0;
    Real snmag      = 1.0/std::sqrt(a4i*a4i+a4j*a4j+a4k*a4k);
    normals0[0+3*i] = a4i*snmag;
    normals0[1+3*i] = a4j*snmag;
    normals0[2+3*i] = a4k*snmag;
  }
  if (coordinates1!=NULL) {
    for (i=0; i<2; ++i) {
      local_coords0[0] = -0.5+1.0*(Real)i;
      for (j=0; j<2; ++j) {
        Real center[3];
        local_coords0[1] = -0.5+1.0*(Real)j;
        Compute_Global_Coordinates( POSITION1, local_coords0, center );
        local_coords1[0] = local_coords0[0]+offset0[0];
        local_coords1[1] = local_coords0[1]+offset1[0];
        Compute_Global_Coordinates( POSITION1, local_coords1, node0 );
        for (n=0; n<8; ++n) {
          for (ii=0, jj=1; ii<8; ++ii, jj=(ii+1)%8) {
            nn = (i*2+j)*8+ii;
            local_coords1[0] = local_coords0[0]+offset0[jj];
            local_coords1[1] = local_coords0[1]+offset1[jj];
            Compute_Global_Coordinates( POSITION1, local_coords1, node1 );
            coordinates1[0+9*nn] = node0[0];
            coordinates1[1+9*nn] = node0[1];
            coordinates1[2+9*nn] = node0[2];
            coordinates1[3+9*nn] = node1[0];
            coordinates1[4+9*nn] = node1[1];
            coordinates1[5+9*nn] = node1[2];
            coordinates1[6+9*nn] = center[0];
            coordinates1[7+9*nn] = center[1];
            coordinates1[8+9*nn] = center[2];
            node0[0] = node1[0];
            node0[1] = node1[1];
            node0[2] = node1[2];
          }
        }
      }
    }
    for (i=0; i<nfacets; ++i) {
      // compute vectors from centroid to other
      // two nodes defining master surface triangle
      Real vx0   = coordinates1[6+9*i] - coordinates1[0+9*i];
      Real vy0   = coordinates1[7+9*i] - coordinates1[1+9*i];
      Real vz0   = coordinates1[8+9*i] - coordinates1[2+9*i];
      Real vx00  = coordinates1[3+9*i] - coordinates1[0+9*i];
      Real vy00  = coordinates1[4+9*i] - coordinates1[1+9*i];
      Real vz00  = coordinates1[5+9*i] - coordinates1[2+9*i];
      // i,j,k components of the normal to the plane
      Real a4i        = -vy0*vz00+vy00*vz0;
      Real a4j        =  vx0*vz00-vx00*vz0;
      Real a4k        = -vx0*vy00+vx00*vy0;
      Real snmag      = 1.0/std::sqrt(a4i*a4i+a4j*a4j+a4k*a4k);
      normals1[0+3*i] = a4i*snmag;
      normals1[1+3*i] = a4j*snmag;
      normals1[2+3*i] = a4k*snmag;
    }
  }
  if (coordinates2!=NULL) {
    for (i=0; i<2; ++i) {
      local_coords0[0] = -0.5+1.0*(Real)i;
      for (j=0; j<2; ++j) {
        Real center[3];
        local_coords0[1] = -0.5+1.0*(Real)j;
        Compute_Global_Coordinates( POSITION2, local_coords0, center );
        local_coords1[0] = local_coords0[0]+offset0[0];
        local_coords1[1] = local_coords0[1]+offset1[0];
        Compute_Global_Coordinates( POSITION2, local_coords1, node0 );
        for (n=0; n<8; ++n) {
          for (ii=0, jj=1; ii<8; ++ii, jj=(ii+1)%8) {
            nn = (i*2+j)*8+ii;
            local_coords1[0] = local_coords0[0]+offset0[jj];
            local_coords1[1] = local_coords0[1]+offset1[jj];
            Compute_Global_Coordinates( POSITION2, local_coords1, node1 );
            coordinates2[0+9*nn] = node0[0];
            coordinates2[1+9*nn] = node0[1];
            coordinates2[2+9*nn] = node0[2];
            coordinates2[3+9*nn] = node1[0];
            coordinates2[4+9*nn] = node1[1];
            coordinates2[5+9*nn] = node1[2];
            coordinates2[6+9*nn] = center[0];
            coordinates2[7+9*nn] = center[1];
            coordinates2[8+9*nn] = center[2];
            node0[0] = node1[0];
            node0[1] = node1[1];
            node0[2] = node1[2];
          }
        }
      }
    }
    for (i=0; i<nfacets; ++i) {
      // compute vectors from centroid to other
      // two nodes defining master surface triangle
      Real vx0   = coordinates2[6+9*i] - coordinates2[0+9*i];
      Real vy0   = coordinates2[7+9*i] - coordinates2[1+9*i];
      Real vz0   = coordinates2[8+9*i] - coordinates2[2+9*i];
      Real vx00  = coordinates2[3+9*i] - coordinates2[0+9*i];
      Real vy00  = coordinates2[4+9*i] - coordinates2[1+9*i];
      Real vz00  = coordinates2[5+9*i] - coordinates2[2+9*i];
      // i,j,k components of the normal to the plane
      Real a4i        = -vy0*vz00+vy00*vz0;
      Real a4j        =  vx0*vz00-vx00*vz0;
      Real a4k        = -vx0*vy00+vx00*vy0;
      Real snmag      = 1.0/std::sqrt(a4i*a4i+a4j*a4j+a4k*a4k);
      normals2[0+3*i] = a4i*snmag;
      normals2[1+3*i] = a4j*snmag;
      normals2[2+3*i] = a4k*snmag;
    }
  }
}

void ContactQuadFaceQ9::FacetStaticRestriction(int nfacets, 
                                               Real* coordinates, 
                                               Real* normals, 
                                               Real* ctrcl_facets, 
                                               Real* ctrcl)
{
  int ii1=0,ilocc=0,ilocs=0;
  int iistored=0,iconcave=0,iinside=1,iout=2;
  Real projcv,projmv,one_third=1.0/3.0;
  Real dctds,xmc,ymc,zmc,xms,yms,zms,vecmix,vecmiy,vecmiz;
  // Contract the closest point projection of the Node-with-N_Triangles
  // into a single contact with the bilinear QUAD master surface.
  for (int ipass=0; ipass<nfacets; ++ipass) {
    int ii = ipass;
    if (ctrcl_facets[ii*LENGTH+MSPARAM] != 0 ){
      if ( ii == ii1 || (int)(ctrcl_facets[ii1*LENGTH+MSPARAM]+0.5) == 0 ) {
	for (int k=0; k<LENGTH; ++k) {
          ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
	}
        iistored = ii;
      } else {
        // a contact is already stored in location ctrcl_facets(*,ii1),
        // choose the best contact
        //     indicator for sign on current and stored gaps
        //     (+ means both are penetrating or not penetrating,
        //      - means that one is penetrating and one is not penetrating)
        dctds = ctrcl_facets[ii*LENGTH+IPENMAG]*ctrcl_facets[ii1*LENGTH+IPENMAG];
        // location of contact in-the-plane
        ilocc = (int)(ctrcl_facets[ii*LENGTH+ILOCATION]+0.5);
        ilocs = (int)(ctrcl_facets[ii1*LENGTH+ILOCATION]+0.5);
        // Determine if the two tringular surfaces are concave or convex
        // relative to each other
        // centroid of candidate master surface triangle
        xmc = one_third*(coordinates[0+9*ii] + 
                         coordinates[3+9*ii] + 
                         coordinates[6+9*ii]);
        ymc = one_third*(coordinates[1+9*ii] + 
                         coordinates[4+9*ii] + 
                         coordinates[7+9*ii]);
        zmc = one_third*(coordinates[2+9*ii] + 
                         coordinates[5+9*ii] + 
                         coordinates[8+9*ii]);
        // centroid of stored master surface triangle
        xms = one_third*(coordinates[0+9*iistored] + 
                         coordinates[3+9*iistored] + 
                         coordinates[6+9*iistored]);
        yms = one_third*(coordinates[1+9*iistored] + 
                         coordinates[4+9*iistored] + 
                         coordinates[7+9*iistored]);
        zms = one_third*(coordinates[2+9*iistored] + 
                         coordinates[5+9*iistored] + 
                         coordinates[8+9*iistored]);
        // vector from centroid of previous (stored) m.s. to current m.s.
        vecmix = xmc - xms;
        vecmiy = ymc - yms;
        vecmiz = zmc - zms;
        // projection of previous (stored) normal onto vecmi
        projmv = normals[0+3*iistored]*vecmix + 
	         normals[1+3*iistored]*vecmiy + 
	         normals[2+3*iistored]*vecmiz;
        // projection of current normal onto vecmi
        projcv = normals[0+3*ii]*vecmix + 
	         normals[1+3*ii]*vecmiy + 
                 normals[2+3*ii]*vecmiz;
	if(projcv > projmv) {
	  // surfaces are concave
          iconcave = 1;
        } else {
          // surfaces are convex
          iconcave = -1;
	}
        // decision time
	if( std::abs(ilocc) == iinside ) {
          // INSIDE the candidate m.s. triangle
          //==============================
          if( std::abs(ilocs) == iinside ) {
            // INSIDE the stored m.s. triangle
            // ->  a) choose closest
            //     b) choose the m.s. whose normal most opposes s.n. normal
            // =>  c) choose the m.s. whose normal most opposes s.n. contact force
            if( std::fabs(ctrcl_facets[ii*LENGTH+IPENMAG]) < 
		std::fabs(ctrcl_facets[ii1*LENGTH+IPENMAG])) {
	      for (int k=0; k<LENGTH; ++k) {
                ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
	      }
              iistored = ii;
	    }
          } else {
            // OUTSIDE the stored m.s. triangle
            // -=> a) choose the INSIDE contact (i.e. replace stored m.s.)
            for(int k=0; k<LENGTH; ++k) {
              ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
	    }
            iistored = ii;
          }
        }  else if( std::abs(ilocc) == iout ) {
          // OUTSIDE the candidate m.s. triangle
          //====================================
	  if( std::abs(ilocs) == iinside ) {
            // INSIDE the stored m.s. triangle
            //   -=> a) choose the INSIDE contact (i.e. keep stored m.s.)
	  } else if( std::abs(ilocs) == iout ) {
            // OUTSIDE the stored m.s. triangle
            if( dctds == 0.0 && iconcave <= 0 ) {
              // -=> a)  for a convex surface with one zero penetration contact
              if( ctrcl_facets[ii*LENGTH+IPENMAG] < ctrcl_facets[ii1*LENGTH+IPENMAG] ){
		for (int k=0; k<LENGTH; ++k) {
                  ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
		}
                iistored = ii;
              }
            } else if( dctds == 0.0 && iconcave > 0 ) {
              // b)  for a concave surface with one zero penetration contact
              if( ctrcl_facets[ii*LENGTH+IPENMAG] > ctrcl_facets[ii1*LENGTH+IPENMAG] ){
		for (int k=0; k<LENGTH; ++k) {
                  ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
		}
                iistored = ii;
              }
            } else if( dctds > 0 && iconcave <= 0 ) {
              if( std::fabs(ctrcl_facets[ii*LENGTH+IPENMAG]) < 
                  std::fabs(ctrcl_facets[ii1*LENGTH+IPENMAG])){ 
		for (int k=0; k<LENGTH; ++k) {
                  ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
		}
                iistored = ii;
              }
            } else if( dctds > 0.0 && iconcave > 0 ) {
              // b)  for a concave surface, with both contacts either
              //     penetrating or not penetrating, choose:
              // ->  b.1) choose closest
              //     b.2) choose the m.s. whose normal most opposes s.n. normal
              // =>  b.3) choose the m.s. whose normal most opposes s.n.
              //          contact force
              if(std::fabs(ctrcl_facets[ii*LENGTH+IPENMAG])<
		 std::fabs(ctrcl_facets[ii1*LENGTH+IPENMAG])){ 
		for (int k=0; k<LENGTH; ++k) {
                  ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
		}
                iistored = ii;
              }
            } else if( dctds < 0.0 && iconcave <= 0 ) {
              // -=> c)  for a convex surface, with one contact penetrating 
              //         and the other contact not penetrating, choose the 
              //         contact that is penetrating
              if( ctrcl_facets[ii*LENGTH+IPENMAG] < ctrcl_facets[ii1*LENGTH+IPENMAG] ){ 
		for (int k=0; k<LENGTH; ++k) {
                  ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
		}
                iistored = ii;
              }
            } else if( dctds < 0.0 && iconcave > 0 ) {
              //-=> d)  for a concave surface, with one contact penetrating 
              //        and the other contact not penetrating, choose the 
              //        contact that is not penetrating
              if( ctrcl_facets[ii*LENGTH+IPENMAG] > ctrcl_facets[ii1*LENGTH+IPENMAG] ){ 
		for (int k=0; k<LENGTH; ++k) {
                  ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
		}
                iistored = ii;
              }
            }
	  }
        }
      }
    }
  }
  for (int i=0; i<LENGTH; ++i) {
    ctrcl[i] = ctrcl_facets[i];
  }
}

void ContactQuadFaceQ9::FacetDynamicRestriction(int nfacets, 
                                                Real* ctrcl_facets, 
                                                Real* ctrcl)
{
  int i;
  // There are three possibilities with each triangle
  //   1) Accepted =  1
  //   2) Rejected =  0
  //   3) Deferred = -1   (send to Closest Point Proejection)
  //
  // Our logic here is as follows
  //   1) Always keep an accepted subtriangle over rejected or deferred
  //   2) Keep deferred over rejected and send all to CPP
  //   3) Only reject if all subtriangles were rejected

  // Get maximum over all subtriangles
  int max_status = -2;
  int min_status =  2;
  for( i=0 ; i<nfacets ; ++i ){
    int current_status = (int) ctrcl_facets[i*LENGTH+MSPARAM];
    max_status = std::max(max_status,current_status);
    min_status = std::min(min_status,current_status);
  }

  if( min_status == 0  && max_status == 0 ){
    // all rejected
    ctrcl[MSPARAM] = 0;
    return;
  } else if( min_status == -1 && max_status == 0 ){
    // some rejected, some deferred
    ctrcl[MSPARAM] = -1;
    return;
  }
  
  // Contract the moving search of the Node-with-N_Triangles into a single
  // contact with a bilinear QUAD master surface
  //
  // Notes:  The memory is laid out so that we don't need to look at the first
  //         triangle.  If its not in contact it will be ignored above.  If it
  //         is in contact it should be std::left alone and will correctly competed.
  for (int ii=1; ii<nfacets; ++ii) {
    int index;
    if( ctrcl_facets[ii*LENGTH+MSPARAM] == 1 ) {
      if( ctrcl_facets[MSPARAM] != 1 ) {
        // No constraint is yet defined for this quad so store this one
	for (index=0; index<LENGTH; ++index) {
          ctrcl_facets[index] = ctrcl_facets[ii*LENGTH+index];
	}
      } else {
        // A constraint is already defined for this quad.  Take the one
        // with the minimum contact time.
        if( ctrcl_facets[ii*LENGTH+ICTIMC] < ctrcl_facets[ICTIMC] ) {
	  if( ctrcl_facets[ii*LENGTH+ICTIMC] > 0.0 ) {
	    for (index=0; index<LENGTH; ++index) {
              ctrcl_facets[index] = ctrcl_facets[ii*LENGTH+index];
	    }
	  }
	}
      }
    }
  }
  for (i=0; i<LENGTH; ++i) {
    ctrcl[i] = ctrcl_facets[i];
  }
}


void ContactQuadFaceQ9::Smooth_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle NODE_NORMAL,
				       VariableHandle FACE_NORMAL,
				       VariableHandle CURVATURE,
			ContactSearch::Smoothing_Resolution resolution,
				       Real percentage,
				       Real* coordinates,
				       Real* smooth_normal,
                                       Real critical_curvature )
{
  //
  //  3-y-z--------6---------Z-Y-2
  //  |   |                  |   |
  //  u I v        D         w H x
  //  |   |                  |   |
  //  o-p-q------------------r-s-t
  //  |   |                  |   |
  //  |   |                  |   |
  //  7 E |        A         | C 5
  //  |   |                  |   |
  //  |   |                  |   |
  //  i-j-k------------------l-m-n
  //  |   |                  |   |
  //  e F f        B         g G h
  //  |   |                  |   |
  //  0-a-b--------4---------c-d-1
  //
  // Region     A: Use surface normal
  //          B-E: Smooth along one edge
  //          F-I: Smooth along two edges
  //
  // For region A, we always take the face normal.
  //
  // For regions B-E, there are two possibilites to consider
  //   1) No smoothing along the edge in which case just take the face normal.
  //   2) Smoothing is requested along the edge.  The face normal is used as
  //      the values at nodes 0 & 1, and the normal at nodes 2 & 3 is computed
  //      by the edge.
  //
  // For regions F-I, there are three possibilities to consider
  //   1) No smoothing is needed on either edge, in which case just take the
  //      face normal.
  //   2) Smoothing is indicated along both edges.  For this case, we use the
  //      face normal at node 0, the edge normals at nodes 1 & 3, and the node
  //      normal at node 2.
  //   3) Smoothing in indicated along only one of the edges but not the other.
  //      For this case, we use the face normal for node 0, the edge normal at 
  //      node 1 or 3 depending on which edge is smoothed, the face normal at
  //      node 3 or 1 depending on which edge is not smoothed, and we compute
  //      the normal at node 2 is computed by the value along the smoothed 
  //      edge.

  Real upper_bound =  1.0 - percentage;
  Real lower_bound = -1.0 + percentage;
  Real* face_normal = Variable(FACE_NORMAL);
  Real node_position[9][3];
  Real node_normal[9][3];
  Real contact_point[3];
  Real coords[2];
  Real Mag;
  bool smooth_edge0,smooth_edge1;
  Real curvature0,curvature1;

  coords[0] = coordinates[0];
  coords[1] = coordinates[1];
  Compute_Global_Coordinates( CURRENT_POSITION, coords, &(contact_point[0]) );

  if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
      coordinates[1] >= lower_bound && coordinates[1] <= upper_bound ){
    // Region A
    smooth_normal[0] = face_normal[0];
    smooth_normal[1] = face_normal[1];
    smooth_normal[2] = face_normal[2];
    return;
  }
  else if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
	   coordinates[1] < lower_bound ){
    // Region B
    //
    //  1--------4----------0
    //  |                   |
    //  5    B   8          7
    //  |                   |
    //  2--------6----------3
    //
    if( fabs(GetEdgeCurvature(0)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = lower_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = upper_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
    // Compute data for node 4
    coords[0] = 0.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute data for node 5
    coords[0] = lower_bound;
    coords[1] = (-1.0 + lower_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = 0.0;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    GetEdgeSmoothedNormal(0, &(node_normal[6][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute data for node 7
    coords[0] = upper_bound;
    coords[1] = (-1.0 + lower_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[3][0]+node_normal[0][0];
    node_normal[7][1] = node_normal[3][1]+node_normal[0][1];
    node_normal[7][2] = node_normal[3][2]+node_normal[0][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    } 
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }   
  }
  else if( coordinates[0] > upper_bound &&
	   coordinates[1] <= upper_bound && coordinates[1] >= lower_bound ){
    // Region C
    //
    //  0-7-3
    //  |   |
    //  |   |
    //  4 C 6
    //  |   |
    //  |   |
    //  |   |
    //  1-5-2
    //
    if( fabs(GetEdgeCurvature(1)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = upper_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = 1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = 1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
    // Compute data for node 4
    coords[0] = upper_bound;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute data for node 5
    coords[0] = (1.0+upper_bound)/2.0;
    coords[1] = (lower_bound);
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = 1.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    GetEdgeSmoothedNormal(1, &(node_normal[6][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute data for node 7
    coords[0] = (1.0+upper_bound)/2.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[3][0]+node_normal[0][0];
    node_normal[7][1] = node_normal[3][1]+node_normal[0][1];
    node_normal[7][2] = node_normal[3][2]+node_normal[0][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }       
  }
  else if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
	   coordinates[1] > upper_bound ){
    // Region D
    //
    //  3--------6----------2
    //  |                   |
    //  7        D          5
    //  |                   |
    //  0--------4----------1
    //
    if( fabs(GetEdgeCurvature(2)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = upper_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = upper_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = lower_bound;
    coords[1] = 1.0;
    
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
    // Compute data for node 4
    coords[0] = 0.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute data for node 5
    coords[0] = upper_bound;
    coords[1] = (upper_bound+1.0)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = 0.0;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    GetEdgeSmoothedNormal(2, &(node_normal[6][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute data for node 7
    coords[0] = lower_bound;
    coords[1] = (upper_bound+1.0)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[3][0]+node_normal[0][0];
    node_normal[7][1] = node_normal[3][1]+node_normal[0][1];
    node_normal[7][2] = node_normal[3][2]+node_normal[0][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }   
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }    
  }
  else if( coordinates[0] < lower_bound &&
	   coordinates[1] <= upper_bound && coordinates[1] >= lower_bound ){
    // Region E
    //
    //  2-5-1
    //  |   |
    //  |   |
    //  6 E 4
    //  |   |
    //  |   |
    //  |   |
    //  3-7-0
    //
    if( fabs(GetEdgeCurvature(3)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = -1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    GetEdgeSmoothedNormal(3, &(node_normal[2][0]));
    //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = -1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    GetEdgeSmoothedNormal(3, &(node_normal[3][0]));
    //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
    // Compute data for node 4
    coords[0] = lower_bound;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute data for node 5
    coords[0] = (lower_bound-1.0)/2.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = -1.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    GetEdgeSmoothedNormal(3, &(node_normal[6][0]));
    //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute data for node 7
    coords[0] = (lower_bound-1.0)/2.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[3][0]+node_normal[0][0];
    node_normal[7][1] = node_normal[3][1]+node_normal[0][1];
    node_normal[7][2] = node_normal[3][2]+node_normal[0][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }       
  }
  else if( coordinates[0] < lower_bound && coordinates[1] < lower_bound ){
    // Region F
    //
    //  1-4-0
    //  |   |
    //  5 F 7
    //  |   |
    //  2-6-3
    //
    curvature0 = GetEdgeCurvature(3);
    curvature1 = GetEdgeCurvature(0);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = -1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(3, &(node_normal[1][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }    
    // Compute data for node 2
    coords[0] = -1.0;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    node_normal[2][0] = Node(0)->Variable(NODE_NORMAL)[0];
    node_normal[2][1] = Node(0)->Variable(NODE_NORMAL)[1];
    node_normal[2][2] = Node(0)->Variable(NODE_NORMAL)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0]   = Node(0)->Variable(NODE_NORMAL)[0];
      node_normal[2][1]   = Node(0)->Variable(NODE_NORMAL)[1];
      node_normal[2][2]   = Node(0)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(3, &(node_normal[2][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute data for node 3
    coords[0] = lower_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
    // Compute data for node 4
    coords[0] = (-1.0+lower_bound)/2.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = node_normal[0][0]+node_normal[1][0];
    node_normal[4][1] = node_normal[0][1]+node_normal[1][1];
    node_normal[4][2] = node_normal[0][2]+node_normal[1][2];
    Mag = node_normal[4][0]*node_normal[4][0] + 
          node_normal[4][1]*node_normal[4][1] +
          node_normal[4][2]*node_normal[4][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[4][0] *= Mag;
      node_normal[4][1] *= Mag;
      node_normal[4][2] *= Mag;
    }
    // Compute data for node 5
    coords[0] = -1.0;
    coords[1] = (-1.0+lower_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = (-1+lower_bound)/2.0;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    node_normal[6][0] = node_normal[2][0]+node_normal[3][0];
    node_normal[6][1] = node_normal[2][1]+node_normal[3][1];
    node_normal[6][2] = node_normal[2][2]+node_normal[3][2];
    Mag = node_normal[6][0]*node_normal[6][0] + 
          node_normal[6][1]*node_normal[6][1] +
          node_normal[6][2]*node_normal[6][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[6][0] *= Mag;
      node_normal[6][1] *= Mag;
      node_normal[6][2] *= Mag;
    }
    // Compute data for node 7
    coords[0] = lower_bound;
    coords[1] = (-1.+lower_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[0][0]+node_normal[3][0];
    node_normal[7][1] = node_normal[0][1]+node_normal[3][1];
    node_normal[7][2] = node_normal[0][2]+node_normal[3][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }   
  }
  else if( coordinates[0] > upper_bound && coordinates[1] < lower_bound ){
    // Region G
    //
    //  0-7-3
    //  |   |
    //  4 G 6
    //  |   |
    //  1-5-2
    //
    curvature0 = GetEdgeCurvature(0);
    curvature1 = GetEdgeCurvature(1);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = upper_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(0, &(node_normal[1][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }    
    // Compute data for node 2
    coords[0] = 1.0;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(1)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(1)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(1)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }      
    // Compute data for node 3
    coords[0] = 1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
    // Compute data for node 4
    coords[0] = upper_bound;
    coords[1] = (-1.0+lower_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = node_normal[0][0]+node_normal[1][0];
    node_normal[4][1] = node_normal[0][1]+node_normal[1][1];
    node_normal[4][2] = node_normal[0][2]+node_normal[1][2];
    Mag = node_normal[4][0]*node_normal[4][0] + 
          node_normal[4][1]*node_normal[4][1] +
          node_normal[4][2]*node_normal[4][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[4][0] *= Mag;
      node_normal[4][1] *= Mag;
      node_normal[4][2] *= Mag;
    }
    // Compute data for node 5
    coords[0] = (1.0+upper_bound)/2.0;
    coords[1] = -1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = 1.0;
    coords[1] = (-1.0+lower_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    node_normal[6][0] = node_normal[2][0]+node_normal[3][0];
    node_normal[6][1] = node_normal[2][1]+node_normal[3][1];
    node_normal[6][2] = node_normal[2][2]+node_normal[3][2];
    Mag = node_normal[6][0]*node_normal[6][0] + 
          node_normal[6][1]*node_normal[6][1] +
          node_normal[6][2]*node_normal[6][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[6][0] *= Mag;
      node_normal[6][1] *= Mag;
      node_normal[6][2] *= Mag;
    }
    // Compute data for node 7
    coords[0] = (1.0+upper_bound)/2.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[0][0]+node_normal[3][0];
    node_normal[7][1] = node_normal[0][1]+node_normal[3][1];
    node_normal[7][2] = node_normal[0][2]+node_normal[3][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }   
  }
  else if( coordinates[0] > upper_bound && coordinates[1] > upper_bound ){
    // Region H
    //
    //  3-6-2
    //  |   |
    //  7 H 5
    //  |   |
    //  0-4-1
    //
    curvature0 = GetEdgeCurvature(1);
    curvature1 = GetEdgeCurvature(2);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = upper_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = 1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(1, &(node_normal[1][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute data for node 2
    coords[0] = 1.0;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(2)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(2)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(2)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }

    // Compute data for node 3
    coords[0] = upper_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }    
    // Compute data for node 4
    coords[0] = (1.0+upper_bound)/2.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = node_normal[0][0]+node_normal[1][0];
    node_normal[4][1] = node_normal[0][1]+node_normal[1][1];
    node_normal[4][2] = node_normal[0][2]+node_normal[1][2];
    Mag = node_normal[4][0]*node_normal[4][0] + 
          node_normal[4][1]*node_normal[4][1] +
          node_normal[4][2]*node_normal[4][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[4][0] *= Mag;
      node_normal[4][1] *= Mag;
      node_normal[4][2] *= Mag;
    }
    // Compute data for node 5
    coords[0] = 1.0;
    coords[1] = (1.0+upper_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = (1.0+upper_bound)/2.0;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    node_normal[6][0] = node_normal[2][0]+node_normal[3][0];
    node_normal[6][1] = node_normal[2][1]+node_normal[3][1];
    node_normal[6][2] = node_normal[2][2]+node_normal[3][2];
    Mag = node_normal[6][0]*node_normal[6][0] + 
          node_normal[6][1]*node_normal[6][1] +
          node_normal[6][2]*node_normal[6][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[6][0] *= Mag;
      node_normal[6][1] *= Mag;
      node_normal[6][2] *= Mag;
    }
    // Compute data for node 7
    coords[0] = upper_bound;
    coords[1] = (1.0+upper_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[0][0]+node_normal[3][0];
    node_normal[7][1] = node_normal[0][1]+node_normal[3][1];
    node_normal[7][2] = node_normal[0][2]+node_normal[3][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }   
  }
  else if( coordinates[0] < lower_bound && coordinates[1] > upper_bound ){
    // Region I
    //
    //  2-5-1
    //  |   |
    //  6 I 4
    //  |   |
    //  3-7-0
    //
    curvature0 = GetEdgeCurvature(2);
    curvature1 = GetEdgeCurvature(3);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[0][0]));
    // Compute data for node 1
    coords[0] = lower_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[1][0]));
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(2, &(node_normal[1][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }    
    // Compute data for node 2
    coords[0] = -1.0;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[2][0]));
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(3)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(3)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(3)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(3, &(node_normal[2][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute data for node 3
    coords[0] = -1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[3][0]));
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(3, &(node_normal[3][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
    // Compute data for node 4
    coords[0] = lower_bound;
    coords[1] = (1.0+upper_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[4][0]));
    node_normal[4][0] = node_normal[0][0]+node_normal[1][0];
    node_normal[4][1] = node_normal[0][1]+node_normal[1][1];
    node_normal[4][2] = node_normal[0][2]+node_normal[1][2];
    Mag = node_normal[4][0]*node_normal[4][0] + 
          node_normal[4][1]*node_normal[4][1] +
          node_normal[4][2]*node_normal[4][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[4][0] *= Mag;
      node_normal[4][1] *= Mag;
      node_normal[4][2] *= Mag;
    }
    // Compute data for node 5
    coords[0] = (-1.0+lower_bound)/2.0;
    coords[1] = 1.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[5][0]));
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] + 
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute data for node 6
    coords[0] = -1.0;
    coords[1] = (1.0+upper_bound)/2.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[6][0]));
    node_normal[6][0] = node_normal[2][0]+node_normal[3][0];
    node_normal[6][1] = node_normal[2][1]+node_normal[3][1];
    node_normal[6][2] = node_normal[2][2]+node_normal[3][2];
    Mag = node_normal[6][0]*node_normal[6][0] + 
          node_normal[6][1]*node_normal[6][1] +
          node_normal[6][2]*node_normal[6][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[6][0] *= Mag;
      node_normal[6][1] *= Mag;
      node_normal[6][2] *= Mag;
    }
    // Compute data for node 7
    coords[0] = (-1.0+lower_bound)/2.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[7][0]));
    node_normal[7][0] = node_normal[0][0]+node_normal[3][0];
    node_normal[7][1] = node_normal[0][1]+node_normal[3][1];
    node_normal[7][2] = node_normal[0][2]+node_normal[3][2];
    Mag = node_normal[7][0]*node_normal[7][0] + 
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
    // Compute data for node 8
    coords[0] = 0.0;
    coords[1] = 0.0;
    Compute_Global_Coordinates(CURRENT_POSITION,coords,&(node_position[8][0]));
    node_normal[8][0] = 0.0;
    node_normal[8][1] = 0.0;
    node_normal[8][2] = 0.0;
    for (int ii=0; ii<8; ++ii) {
      node_normal[8][0] += node_normal[ii][0];
      node_normal[8][1] += node_normal[ii][1];
      node_normal[8][2] += node_normal[ii][2];
    }
    Mag = node_normal[8][0]*node_normal[8][0] + 
          node_normal[8][1]*node_normal[8][1] +
          node_normal[8][2]*node_normal[8][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0 ){
      Mag = 1.0/Mag;
      node_normal[8][0] *= Mag;
      node_normal[8][1] *= Mag;
      node_normal[8][2] *= Mag;
    }   
  }
  Real local_coords[3];
  Compute_Local_Coords( node_position, contact_point,
			&(local_coords[0]));

  // Now interpolate to get the normal
  Interpolate_Vector( &(local_coords[0]), node_normal, 
				  smooth_normal );

  // Now make smooth_normal a unit vector
  Mag = smooth_normal[0]*smooth_normal[0] + 
        smooth_normal[1]*smooth_normal[1] +
        smooth_normal[2]*smooth_normal[2];
  Mag = std::sqrt(Mag);
  if( Mag > 0.0 ){
    Mag = 1.0/Mag;
    smooth_normal[0] *= Mag;
    smooth_normal[1] *= Mag;
    smooth_normal[2] *= Mag;
  }
}


void ContactQuadFaceQ9::Compute_Node_Areas( VariableHandle POSITION,
                                            VariableHandle FACE_NORMAL,
                                            Real* node_areas )
{
  // In the short term, I'm going to compute the node areas as if this were
  // a quad4 and then distribute the areas to the nine nodes

  // Get the Nodal Positions
  Real* pos_0 = Node(0)->Variable( POSITION );
  Real* pos_1 = Node(1)->Variable( POSITION );
  Real* pos_2 = Node(2)->Variable( POSITION );
  Real* pos_3 = Node(3)->Variable( POSITION );

  // Get the face normal
  Real* fnorm = Variable( FACE_NORMAL );

  // Compute vectors between nodes
  Real vec_0_1[3], vec_2_1[3], vec_2_3[3], vec_3_0[3], vec_0_2[3], vec_3_1[3];
  vec_0_1[0] = pos_1[0] - pos_0[0];
  vec_0_1[1] = pos_1[1] - pos_0[1];
  vec_0_1[2] = pos_1[2] - pos_0[2];
  vec_2_1[0] = pos_1[0] - pos_2[0];
  vec_2_1[1] = pos_1[1] - pos_2[1];
  vec_2_1[2] = pos_1[2] - pos_2[2];
  vec_2_3[0] = pos_3[0] - pos_2[0];
  vec_2_3[1] = pos_3[1] - pos_2[1];
  vec_2_3[2] = pos_3[2] - pos_2[2];
  vec_3_0[0] = pos_0[0] - pos_3[0];
  vec_3_0[1] = pos_0[1] - pos_3[1];
  vec_3_0[2] = pos_0[2] - pos_3[2];
  vec_0_2[0] = pos_2[0] - pos_0[0];
  vec_0_2[1] = pos_2[1] - pos_0[1];
  vec_0_2[2] = pos_2[2] - pos_0[2];
  vec_3_1[0] = pos_1[0] - pos_3[0];
  vec_3_1[1] = pos_1[1] - pos_3[1];
  vec_3_1[2] = pos_1[2] - pos_3[2];

  // Compute three inward normals
  Real s1[3], s2[3], s3[3];
  s1[0] = 3.0 * ( vec_0_2[2]*vec_3_1[1] - vec_0_2[1]*vec_3_1[2] );
  s1[1] = 3.0 * ( vec_0_2[0]*vec_3_1[2] - vec_0_2[2]*vec_3_1[0] );
  s1[2] = 3.0 * ( vec_0_2[1]*vec_3_1[0] - vec_0_2[0]*vec_3_1[1] );
  s2[0] = vec_0_1[2]*vec_2_3[1] - vec_0_1[1]*vec_2_3[2];
  s2[1] = vec_0_1[0]*vec_2_3[2] - vec_0_1[2]*vec_2_3[0];
  s2[2] = vec_0_1[1]*vec_2_3[0] - vec_0_1[0]*vec_2_3[1];
  s3[0] = vec_3_0[2]*vec_2_1[1] - vec_3_0[1]*vec_2_1[2];
  s3[1] = vec_3_0[0]*vec_2_1[2] - vec_3_0[2]*vec_2_1[0];
  s3[2] = vec_3_0[1]*vec_2_1[0] - vec_3_0[0]*vec_2_1[1];

  // Compute the nodal areas
  Real v0[3], v1[3], v2[3], v3[3];
  v0[0] =   s1[0] - s2[0] - s3[0]; 
  v0[1] =   s1[1] - s2[1] - s3[1]; 
  v0[2] =   s1[2] - s2[2] - s3[2]; 
  v1[0] =   s1[0] + s2[0] - s3[0];
  v1[1] =   s1[1] + s2[1] - s3[1];
  v1[2] =   s1[2] + s2[2] - s3[2];
  v2[0] =   s1[0] + s2[0] + s3[0];
  v2[1] =   s1[1] + s2[1] + s3[1];
  v2[2] =   s1[2] + s2[2] + s3[2];
  v3[0] =   s1[0] - s2[0] + s3[0];
  v3[1] =   s1[1] - s2[1] + s3[1];
  v3[2] =   s1[2] - s2[2] + s3[2];

  node_areas[0] = (v0[0]*fnorm[0] + v0[1]*fnorm[1] + v0[2]*fnorm[2] ) / 24.0;
  node_areas[1] = (v1[0]*fnorm[0] + v1[1]*fnorm[1] + v1[2]*fnorm[2] ) / 24.0;
  node_areas[2] = (v2[0]*fnorm[0] + v2[1]*fnorm[1] + v2[2]*fnorm[2] ) / 24.0;
  node_areas[3] = (v3[0]*fnorm[0] + v3[1]*fnorm[1] + v3[2]*fnorm[2] ) / 24.0;

  // Now compute the areas for all nine nodes
  node_areas[4] = (2.0/9.0)*( node_areas[0] + node_areas[1] );
  node_areas[5] = (2.0/9.0)*( node_areas[1] + node_areas[2] );
  node_areas[6] = (2.0/9.0)*( node_areas[2] + node_areas[3] );
  node_areas[7] = (2.0/9.0)*( node_areas[3] + node_areas[0] );
  node_areas[8] = (1.0/9.0)*( node_areas[0] + node_areas[1] +
                              node_areas[2] + node_areas[3]);
  node_areas[0] *= (4.0/9.0);
  node_areas[1] *= (4.0/9.0);
  node_areas[2] *= (4.0/9.0);
  node_areas[3] *= (4.0/9.0);
}

int ContactQuadFaceQ9::FaceEdge_Intersection(VariableHandle POSITION,
					     ContactEdge<Real>* edge, Real* coords)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactQuadFaceQ9::FaceEdge_Intersection not yet implemented\n";
#endif
  POSTCONDITION( 0 );
  return 0;
}

void ContactQuadFaceQ9::Evaluate_Shape_Functions( Real* local_coords,
						  Real* shape_functions )
{
  Compute_Shape_Functions( local_coords, shape_functions );
}

void ContactQuadFaceQ9::Compute_Global_Coordinates( VariableHandle POSITION,
						    Real* local_coords,
						    Real* global_coords )
{
  Real node_positions[9][3];
  for(int i=0; i<9; ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactQuadFaceQ9::Compute_Local_Coordinates( Real Config_Param,
						   VariableHandle POSITION0, 
						   VariableHandle POSITION1, 
						   VariableHandle FACE_NORMAL,
						   Real* global_coords,
						   Real* local_coords )
{
  int  i, j;
  Real node_positions[9][3];
  if (Config_Param == 0.0) {
    for (i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION0);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    Real alpha = 1.0 - Config_Param, beta = Config_Param;
    for (i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position0 = Node(i)->Variable(POSITION0);
      Real* node_position1 = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = alpha*node_position0[j]+beta*node_position1[j];
      }
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactQuadFaceQ9::Compute_Local_Coordinates( VariableHandle POSITION,
						   Real* global_coords,
						   Real* local_coords )
{
  int  i, j;
  Real node_positions[9][3];
  for (i=0; i<Nodes_Per_Face(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/* The following functions are supplied for doing generic computations   */
/* on a quadratic Q9 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary Q9 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/


void ContactQuadFaceQ9::Compute_Shape_Functions( Real local_coords[2],
						 Real shape_functions[9] )
{
  Real xi  = local_coords[0];
  Real eta = local_coords[1];
  shape_functions[0] = 0.25*xi*(xi-1.0)*eta*(eta-1.0);
  shape_functions[1] = 0.25*xi*(xi+1.0)*eta*(eta-1.0);
  shape_functions[2] = 0.25*xi*(xi+1.0)*eta*(eta+1.0);
  shape_functions[3] = 0.25*xi*(xi-1.0)*eta*(eta+1.0);
  shape_functions[4] = 0.50*(1.0-xi*xi)*eta*(eta-1.0);
  shape_functions[5] = 0.50*xi*(xi+1.0)*(1.0-eta*eta);
  shape_functions[6] = 0.50*(1.0-xi*xi)*eta*(eta+1.0);
  shape_functions[7] = 0.50*xi*(xi-1.0)*(1.0-eta*eta);
  shape_functions[8] =      (1.0-xi*xi)*(1.0-eta*eta);
}

void ContactQuadFaceQ9::Compute_Shape_Derivatives( Real local_coords[2],
						   Real shape_derivatives[2][9] )
{
  Real xi  = local_coords[0];
  Real eta = local_coords[1];
  shape_derivatives[0][0] =  0.25*(2.0*xi-1.0)*eta*(eta-1.0);
  shape_derivatives[0][1] =  0.25*(2.0*xi+1.0)*eta*(eta-1.0);
  shape_derivatives[0][2] =  0.25*(2.0*xi+1.0)*eta*(eta+1.0);
  shape_derivatives[0][3] =  0.25*(2.0*xi-1.0)*eta*(eta+1.0);
  shape_derivatives[0][4] = -0.50*(2.0*xi    )*eta*(eta-1.0);
  shape_derivatives[0][5] =  0.50*(2.0*xi+1.0)*(1.0-eta*eta);
  shape_derivatives[0][6] = -0.50*(2.0*xi    )*eta*(eta+1.0);
  shape_derivatives[0][7] =  0.50*(2.0*xi-1.0)*(1.0-eta*eta);
  shape_derivatives[0][8] =      -(2.0*xi    )*(1.0-eta*eta);
  
  shape_derivatives[1][0] =  0.25*xi*(xi-1.0)*(2.0*eta-1.0);
  shape_derivatives[1][1] =  0.25*xi*(xi+1.0)*(2.0*eta-1.0);
  shape_derivatives[1][2] =  0.25*xi*(xi+1.0)*(2.0*eta+1.0);
  shape_derivatives[1][3] =  0.25*xi*(xi-1.0)*(2.0*eta+1.0);
  shape_derivatives[1][4] =  0.50*(1.0-xi*xi)*(2.0*eta-1.0);
  shape_derivatives[1][5] = -0.50*xi*(xi+1.0)*(2.0*eta    );
  shape_derivatives[1][6] =  0.50*(1.0-xi*xi)*(2.0*eta+1.0);
  shape_derivatives[1][7] = -0.50*xi*(xi-1.0)*(2.0*eta    );
  shape_derivatives[1][8] =      -(1.0-xi*xi)*(2.0*eta    );
}

void ContactQuadFaceQ9::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
					      Real global_coords[3],
					      Real local_coords[3] )
{
  int  i, j;
  int  nnodes=9;
  Real spatial_tolerance = 1.0e-10;
  //
  // check for coincidence with one of the face nodes
  //
  for (i=0; i<nnodes; ++i) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real d  = std::sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    break;
  case 3:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    break;
  case 4:
    local_coords[0] =  0.0;
    local_coords[1] = -1.0;
    break;
  case 5:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    break;
  case 6:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    break;
  case 7:
    local_coords[0] = -1.0;
    local_coords[1] =  0.0;
    break;
  case 8:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    break;
  }
  if (i<nnodes) {
    local_coords[2] = 0.0;
    return;
  }
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=200;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real s, s0=0.0, s1, ds=0.0;
  Real t, t0=0.0, t1, dt=0.0;
  Real coords[3];
  Real J[3][2], f[3];
  Real shape_derivatives[2][9];
  while (!converged && iterations<max_iterations) {
    coords[0] = s0;
    coords[1] = t0;
    Compute_Global_Coords(node_positions , coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( coords, shape_derivatives );
    for (i=0; i<2; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<nnodes; ++j) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
    Real JT[2][3];
    JT[0][0] = J[0][0];
    JT[0][1] = J[1][0];
    JT[0][2] = J[2][0];
    JT[1][0] = J[0][1];
    JT[1][1] = J[1][1];
    JT[1][2] = J[2][1];
    
    Real JTJ[2][2];
    JTJ[0][0] = JT[0][0]*J[0][0] + JT[0][1]*J[1][0] + JT[0][2]*J[2][0];
    JTJ[0][1] = JT[0][0]*J[0][1] + JT[0][1]*J[1][1] + JT[0][2]*J[2][1];
    JTJ[1][0] = JT[1][0]*J[0][0] + JT[1][1]*J[1][0] + JT[1][2]*J[2][0];
    JTJ[1][1] = JT[1][0]*J[0][1] + JT[1][1]*J[1][1] + JT[1][2]*J[2][1];
    
    Real invJTJ[2][2];
    Real detJTJ  = 1.0/(JTJ[0][0]*JTJ[1][1]-JTJ[0][1]*JTJ[1][0]);
    invJTJ[0][0] =  JTJ[1][1]*detJTJ;
    invJTJ[0][1] = -JTJ[0][1]*detJTJ;
    invJTJ[1][0] = -JTJ[1][0]*detJTJ;
    invJTJ[1][1] =  JTJ[0][0]*detJTJ;

    // APPLY NEWTON ALGORITHM

    Real dx = f[0]-global_coords[0];
    Real dy = f[1]-global_coords[1];
    Real dz = f[2]-global_coords[2];
    s = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
    t = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;
    
    s1 = s0-(invJTJ[0][0]*s+invJTJ[0][1]*t);
    t1 = t0-(invJTJ[1][0]*s+invJTJ[1][1]*t);
    ds = std::fabs(s1-s0);
    dt = std::fabs(t1-t0);
    s0 = s1;
    t0 = t1;
    if (ds<tolerance && dt<tolerance) converged = true;
    ++iterations;
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  if (!converged) {
    std::cerr << "ContactQuadFaceQ9::Compute_Local_Coords() did not converge" 
	 << std::endl;
    std::cerr << "  Computing Coordinates for point (" << global_coords[0]
	 << "," << global_coords[1] << "," << global_coords[2] << ")"
	 << std::endl;
    std::cerr << "  Face Nodal Coordinates:   (" << node_positions[0][0] 
	 << "," << node_positions[0][1] << "," << node_positions[0][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[1][0] 
	 << "," << node_positions[1][1] << "," << node_positions[1][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[2][0] 
	 << "," << node_positions[2][1] << "," << node_positions[2][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[3][0] 
	 << "," << node_positions[3][1] << "," << node_positions[3][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[4][0] 
	 << "," << node_positions[4][1] << "," << node_positions[4][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[5][0] 
	 << "," << node_positions[5][1] << "," << node_positions[5][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[6][0] 
	 << "," << node_positions[6][1] << "," << node_positions[6][2]
	 << ")" << std::endl;
    std::cerr << "                            (" << node_positions[7][0] 
	 << "," << node_positions[7][1] << "," << node_positions[7][2]
	 << ")" << std::endl;
    std::cerr << "  After " << iterations << "iterations, local_coords = ("
	 << s0 << "," << t0 << ")" << std::endl;
    std::cerr << "  Going to continuing processing anyway!!!" << std::endl;
  }
#endif
  // If it's close to any of the edges, snap to it
  if (std::fabs(s0)<1.0+spatial_tolerance) {
    s0 = std::min(s0, 1.0);
    s0 = std::max(s0,-1.0);
  }
  if (std::fabs(t0)<1.0+spatial_tolerance) {
    t0 = std::min(t0, 1.0);
    t0 = std::max(t0,-1.0);
  }
  local_coords[0] = s0;
  local_coords[1] = t0;
  local_coords[2] = 0.0;
}

void ContactQuadFaceQ9::Compute_Global_Coords( Real node_positions[9][3],
					       Real local_coords[2],
					       Real global_coords[3] )
{
  Real N[9];
  int  nnodes=9;
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      global_coords[j] += N[i]*node_positions[i][j];
    }
  }
}

void ContactQuadFaceQ9::Interpolate_Scalar( Real  local_coords[2],
					    Real  node_scalars[9],
					    Real& interpolated_scalar )
{
  Real N[9];
  int  nnodes=9;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void ContactQuadFaceQ9::Interpolate_Vector( Real local_coords[2],
					    Real node_vectors[9][3],
					    Real interpolated_vector[3] )
{
  Real N[9];
  int  nnodes=9;
  interpolated_vector[0] = 0.0;
  interpolated_vector[1] = 0.0;
  interpolated_vector[2] = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      interpolated_vector[j] += N[i]*node_vectors[i][j];
    }
  }
}
