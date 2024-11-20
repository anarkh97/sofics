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

#include "allocators.h"
#include "ContactUtilities.h"
#include "ContactTriFaceQ6.h"
#include "ContactQuadFaceQ8.h"
#include "ContactFixedSizeAllocator.h"
#include <iostream>
#include "ContactNode.h"
#include <cstdlib>
#include <cmath>
#include <new>

using acme::Dot;
using acme::Cross;
using acme::Normalize;
using acme::Magnitude;
using acme::ScalarTripleProduct;
using acme::Scale;

ContactTriFaceQ6::ContactTriFaceQ6( ContactFixedSizeAllocator* alloc,
                                    int Block_Index, 
				    int Index_in_Block,int key ) 
  : ContactFace<Real>( alloc, ContactSearch::TRIFACEQ6,
                 Block_Index, Index_in_Block, key, 
                 nodes, edges, Node_Info, Edge_Info) 
{}

ContactTriFaceQ6* ContactTriFaceQ6::new_ContactTriFaceQ6(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactTriFaceQ6].New_Frag())
             ContactTriFaceQ6(alloc, Block_Index, Index_in_Block, key);
}

void ContactTriFaceQ6_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactTriFaceQ6),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactTriFaceQ6 allocator" );
}

ContactTriFaceQ6::~ContactTriFaceQ6() {}

void ContactTriFaceQ6::Compute_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL )
{
  int   i;
  Real  shape_derivatives[2][6];
  Real  local_coords[3] = {1.0/3.0, 1.0/3.0, 1.0 - 1.0/3.0 - 1.0/3.0};
  Real* face_normal = Variable(FACE_NORMAL);
  
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (i=0; i<6; ++i) {
    Real* node_position = Node(i)->Variable(CURRENT_POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  Real a = Dot(e1,e1);
  Real b = Dot(e2,e2);
  Real c = Dot(e1,e2);
  Real invDetJ = 1.0 / std::sqrt(a*b-c*c);
  Cross(e1,e2,face_normal);
  Scale(face_normal, invDetJ);
}

void ContactTriFaceQ6::Compute_Normal(VariableHandle POSITION,
				      Real* normal,
				      Real* local_coords )
{
  Real shape_derivatives[2][6];
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<6; ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  Real a = Dot(e1,e1);
  Real b = Dot(e2,e2);
  Real c = Dot(e1,e2);
  Real invDetJ = 1.0 / std::sqrt(a*b-c*c);
  Cross(e1,e2,normal);
  Scale(normal,invDetJ);
}

void ContactTriFaceQ6::Compute_Normal(Real** nodal_positions,
				      Real* local_coords,
				      Real* normal )
{
  Real shape_derivatives[2][6];
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<6; ++i) {
    Real* node_position = nodal_positions[i];
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  Real a = Dot(e1,e1);
  Real b = Dot(e2,e2);
  Real c = Dot(e1,e2);
  Real invDetJ = 1.0 / std::sqrt(a*b-c*c);
  Cross(e1,e2,normal);
  Scale(normal, invDetJ);
}

void ContactTriFaceQ6::Compute_CharacteristicLength( VariableHandle CURRENT_POSITION,
						     VariableHandle CHARACTERISTIC_LENGTH )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  ContactNode<Real>* node2 = Node(2);
  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  Real* Position2 = node2->Variable(CURRENT_POSITION);
  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  Real Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  Real face_normal[3];
  Cross(Vec01,Vec02,face_normal);

  Real* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  *characteristiclength = Magnitude(face_normal) * 0.5;
}

void ContactTriFaceQ6::Compute_Centroid( VariableHandle CURRENT_POSITION,
					 VariableHandle CENTROID )
{
// Assumes straight-sided triangle
  Real* centroid = Variable(CENTROID);
  Real* node_position = Node(0)->Variable(CURRENT_POSITION);
  centroid[0] = node_position[0];
  centroid[1] = node_position[1];
  centroid[2] = node_position[2];
  for( int i=1 ; i<3 ; ++i ){
    node_position = Node(i)->Variable(CURRENT_POSITION);
    centroid[0] += node_position[0];
    centroid[1] += node_position[1];
    centroid[2] += node_position[2];
  }
  Real scale   = 1.0/3.0;
  centroid[0] *= scale;
  centroid[1] *= scale;
  centroid[2] *= scale;
}

void ContactTriFaceQ6::Get_Edge_Nodes( int i, ContactNode<Real>** node )
{
  PRECONDITION( i>=0 && i<3 );
  switch( i ){
  case 0:
    node[0] = Node(0);
    node[1] = Node(1);
    node[2] = Node(3);
    break;
  case 1:
    node[0] = Node(1);
    node[1] = Node(2);
    node[2] = Node(4);
    break;
  case 2:
    node[0] = Node(2);
    node[1] = Node(0);
    node[2] = Node(5);
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Invalid Edge request in ContactTriFaceQ6::Get_Face_Nodes" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

int ContactTriFaceQ6::Get_Edge_Number( ContactNode<Real>** edge_nodes )
{
  PRECONDITION( edge_nodes[0] && edge_nodes[1] && edge_nodes[2] );
  PRECONDITION( edge_nodes[0] != edge_nodes[2] );
  int i1=-1,i2=-1,i3=-1;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    ContactNode<Real>* node = Node(i);
    if( edge_nodes[0] == node ) i1=i;
    if( edge_nodes[1] == node ) i2=i;
    if( edge_nodes[2] == node ) i3=i;
  }
  PRECONDITION( 0<=i1 && i1<3 );
  PRECONDITION( 0<=i2 && i2<3 );
  PRECONDITION( 3<=i3 && i3<6 );
  switch( i1 ){
  case 0:
    if( i2 == 1 && i3==3) return(0);
    if( i2 == 2 && i3==5) return(2);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 1:
    if( i2 == 0 && i3==3) return(0);
    if( i2 == 2 && i3==4) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 2:
    if( i2 == 0 && i3==5) return(2);
    if( i2 == 1 && i3==4) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  }
  POSTCONDITION( 0 );
  return(-1);
}

int ContactTriFaceQ6::Get_Edge_Number( Real* local_coords )
{
  if (local_coords[1]==-1.0 && local_coords[3]==-1.0) return(0);
  if (local_coords[0]==-1.0 && local_coords[2]==-1.0) return(1);
  if ((1.0-local_coords[0]-local_coords[1])==-1 &&
      (1.0-local_coords[2]-local_coords[3])==-1) return(2);
  return( -1 );
}

void
ContactTriFaceQ6::Compute_Edge_Normal( VariableHandle POSITION,
				       VariableHandle FACE_NORMAL,
				       int Edge,
				       Real* edge_normal)
{
  // Edge tangent assumes straight-sided triangle
  ContactNode<Real> *edge_node[3];
  Get_Edge_Nodes( Edge, edge_node );

  // Compute the tangent direction as a vector from node 1 to node 2
  Real* position1 = edge_node[0]->Variable(POSITION);
  Real* position2 = edge_node[1]->Variable(POSITION);

  Real edge_tangent[3];
  edge_tangent[0] = position2[0] - position1[0];
  edge_tangent[1] = position2[1] - position1[1];
  edge_tangent[2] = position2[2] - position1[2];

  // Compute the edge normal direction as the cross product of the edge tangent
  // and the face normal.
  Real* f_normal = Variable(FACE_NORMAL);
  Cross(edge_tangent, f_normal, edge_normal);

  // Normalize the edge_normal
  Normalize(edge_normal);
}

void ContactTriFaceQ6::Get_Close_Edges( Real* local_coords, int& number,
					int& edge_1_id, int& edge_2_id ){
  Real a0 = local_coords[0];
  Real a1 = local_coords[1];
  Real a2 = local_coords[2];

  int edge_case = 0;

  if( a0 <= 0.1 ) edge_case += 1;  // Near Edge 1
  if( a1 <= 0.1 ) edge_case += 2;  // Near Edge 2
  if( a2 <= 0.1 ) edge_case += 4;  // Near Edge 0

  switch( edge_case ){
  case(0):
    // not near any edge
    number = 0;
    break;
  case(1):
    // near edge 1 only
    number = 1;
    edge_1_id = 1;
    break;
  case(2):
    // near edge 2 only
    number = 1;
    edge_1_id = 2;
    break;
  case(4):
    // near edge 0 only
    number = 1;
    edge_1_id = 0;
    break;
  case(3):
    // near edges 1 & 2
    number = 2;
    edge_1_id = 1;
    edge_2_id = 2;
    break;
  case(5): 
    // near edge 0 & 1
    number = 2;
    edge_1_id = 0;
    edge_2_id = 1;
    break;
  case(6):
    // near edges 0 & 2
    number = 2;
    edge_1_id = 0;
    edge_2_id = 2;
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "error can only be near at most 2 edges" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

bool ContactTriFaceQ6::Is_Inside_Face( Real* local_coords )
{
  if( local_coords[0] >= 0.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= 0.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= 0.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

ContactFace<Real>* ContactTriFaceQ6::Neighbor( Real* local_coords )
{
  PRECONDITION( 0 );
  return (ContactFace<Real>*) NULL;
}

void ContactTriFaceQ6::FacetDecomposition(int& nfacets,
          Real* coordinates0, Real* normals0, VariableHandle POSITION0,
          Real* coordinates1, Real* normals1, VariableHandle POSITION1,
          Real* coordinates2, Real* normals2, VariableHandle POSITION2)
{
  Real* node_positiona=NULL;
  Real* node_positionb=NULL;
  Real* node_positionc=NULL;
  nfacets = 4;
  for (int j=0; j<nfacets; ++j) {
    switch(j) {
    case 0:
      node_positiona = Node(0)->Variable(POSITION0);
      node_positionb = Node(3)->Variable(POSITION0);
      node_positionc = Node(5)->Variable(POSITION0);
      break;
    case 1:
      node_positiona = Node(1)->Variable(POSITION0);
      node_positionb = Node(4)->Variable(POSITION0);
      node_positionc = Node(3)->Variable(POSITION0);
      break;
    case 2:
      node_positiona = Node(2)->Variable(POSITION0);
      node_positionb = Node(5)->Variable(POSITION0);
      node_positionc = Node(4)->Variable(POSITION0);
      break;
    case 3:
      node_positiona = Node(3)->Variable(POSITION0);
      node_positionb = Node(4)->Variable(POSITION0);
      node_positionc = Node(5)->Variable(POSITION0);
      break;
    default:
      node_positiona = NULL;
      node_positionb = NULL;
      node_positionc = NULL;
      break;
    }
    POSTCONDITION(node_positiona!=NULL);
    POSTCONDITION(node_positionb!=NULL);
    POSTCONDITION(node_positionc!=NULL);
    coordinates0[0+9*j]  = node_positiona[0];
    coordinates0[1+9*j]  = node_positiona[1];
    coordinates0[2+9*j]  = node_positiona[2];
    coordinates0[3+9*j]  = node_positionb[0];
    coordinates0[4+9*j]  = node_positionb[1];
    coordinates0[5+9*j]  = node_positionb[2];
    coordinates0[6+9*j]  = node_positionc[0];
    coordinates0[7+9*j]  = node_positionc[1];
    coordinates0[8+9*j]  = node_positionc[2];
    // Compute the vector from node a to node b & from node a to node c
    Real vec_ab[3],vec_ac[3];
    vec_ab[0] = node_positionb[0] - node_positiona[0];
    vec_ab[1] = node_positionb[1] - node_positiona[1];
    vec_ab[2] = node_positionb[2] - node_positiona[2];
    vec_ac[0] = node_positionc[0] - node_positiona[0];
    vec_ac[1] = node_positionc[1] - node_positiona[1];
    vec_ac[2] = node_positionc[2] - node_positiona[2];
    // Compute the face normal as the cross product of the two vectors
    Real a4i = vec_ab[1]*vec_ac[2] - vec_ab[2]*vec_ac[1];
    Real a4j = vec_ab[2]*vec_ac[0] - vec_ab[0]*vec_ac[2];
    Real a4k = vec_ab[0]*vec_ac[1] - vec_ab[1]*vec_ac[0];
    Real mag = std::sqrt( a4i*a4i + a4j*a4j + a4k*a4k );
    if( mag > 0.0){
      mag = 1.0/mag;
      normals0[0+3*j] = a4i * mag;
      normals0[1+3*j] = a4j * mag;
      normals0[2+3*j] = a4k * mag;
    }
  }
  if (coordinates1!=NULL) {
    for (int j=0; j<nfacets; ++j) {
      switch(j) {
      case 0:
        node_positiona = Node(0)->Variable(POSITION1);
        node_positionb = Node(3)->Variable(POSITION1);
        node_positionc = Node(5)->Variable(POSITION1);
        break;
      case 1:
        node_positiona = Node(1)->Variable(POSITION1);
        node_positionb = Node(4)->Variable(POSITION1);
        node_positionc = Node(3)->Variable(POSITION1);
        break;
      case 2:
        node_positiona = Node(2)->Variable(POSITION1);
        node_positionb = Node(5)->Variable(POSITION1);
        node_positionc = Node(4)->Variable(POSITION1);
        break;
      case 3:
        node_positiona = Node(3)->Variable(POSITION1);
        node_positionb = Node(4)->Variable(POSITION1);
        node_positionc = Node(5)->Variable(POSITION1);
        break;
      default:
        node_positiona = NULL;
        node_positionb = NULL;
        node_positionc = NULL;
        break;
      }
      POSTCONDITION(node_positiona!=NULL);
      POSTCONDITION(node_positionb!=NULL);
      POSTCONDITION(node_positionc!=NULL);
      coordinates1[0+9*j]  = node_positiona[0];
      coordinates1[1+9*j]  = node_positiona[1];
      coordinates1[2+9*j]  = node_positiona[2];
      coordinates1[3+9*j]  = node_positionb[0];
      coordinates1[4+9*j]  = node_positionb[1];
      coordinates1[5+9*j]  = node_positionb[2];
      coordinates1[6+9*j]  = node_positionc[0];
      coordinates1[7+9*j]  = node_positionc[1];
      coordinates1[8+9*j]  = node_positionc[2];
      // Compute the vector from node a to node b & from node a to node c
      Real vec_ab[3],vec_ac[3];
      vec_ab[0] = node_positionb[0] - node_positiona[0];
      vec_ab[1] = node_positionb[1] - node_positiona[1];
      vec_ab[2] = node_positionb[2] - node_positiona[2];
      vec_ac[0] = node_positionc[0] - node_positiona[0];
      vec_ac[1] = node_positionc[1] - node_positiona[1];
      vec_ac[2] = node_positionc[2] - node_positiona[2];
      // Compute the face normal as the cross product of the two vectors
      Real a4i = vec_ab[1]*vec_ac[2] - vec_ab[2]*vec_ac[1];
      Real a4j = vec_ab[2]*vec_ac[0] - vec_ab[0]*vec_ac[2];
      Real a4k = vec_ab[0]*vec_ac[1] - vec_ab[1]*vec_ac[0];
      Real mag = std::sqrt( a4i*a4i + a4j*a4j + a4k*a4k );
      if( mag > 0.0){
        mag = 1.0/mag;
        normals1[0+3*j] = a4i * mag;
        normals1[1+3*j] = a4j * mag;
        normals1[2+3*j] = a4k * mag;
      }
    }
  }
  if (coordinates2!=NULL) {
    for (int j=0; j<nfacets; ++j) {
      switch(j) {
      case 0:
        node_positiona = Node(0)->Variable(POSITION2);
        node_positionb = Node(3)->Variable(POSITION2);
        node_positionc = Node(5)->Variable(POSITION2);
        break;
      case 1:
        node_positiona = Node(1)->Variable(POSITION2);
        node_positionb = Node(4)->Variable(POSITION2);
        node_positionc = Node(3)->Variable(POSITION2);
        break;
      case 2:
        node_positiona = Node(2)->Variable(POSITION2);
        node_positionb = Node(5)->Variable(POSITION2);
        node_positionc = Node(4)->Variable(POSITION2);
        break;
      case 3:
        node_positiona = Node(3)->Variable(POSITION2);
        node_positionb = Node(4)->Variable(POSITION2);
        node_positionc = Node(5)->Variable(POSITION2);
        break;
      default:
        node_positiona = NULL;
        node_positionb = NULL;
        node_positionc = NULL;
        break;
      }
      POSTCONDITION(node_positiona!=NULL);
      POSTCONDITION(node_positionb!=NULL);
      POSTCONDITION(node_positionc!=NULL);
      coordinates2[0+9*j]  = node_positiona[0];
      coordinates2[1+9*j]  = node_positiona[1];
      coordinates2[2+9*j]  = node_positiona[2];
      coordinates2[3+9*j]  = node_positionb[0];
      coordinates2[4+9*j]  = node_positionb[1];
      coordinates2[5+9*j]  = node_positionb[2];
      coordinates2[6+9*j]  = node_positionc[0];
      coordinates2[7+9*j]  = node_positionc[1];
      coordinates2[8+9*j]  = node_positionc[2];
      // Compute the vector from node a to node b & from node a to node c
      Real vec_ab[3],vec_ac[3];
      vec_ab[0] = node_positionb[0] - node_positiona[0];
      vec_ab[1] = node_positionb[1] - node_positiona[1];
      vec_ab[2] = node_positionb[2] - node_positiona[2];
      vec_ac[0] = node_positionc[0] - node_positiona[0];
      vec_ac[1] = node_positionc[1] - node_positiona[1];
      vec_ac[2] = node_positionc[2] - node_positiona[2];
      // Compute the face normal as the cross product of the two vectors
      Real a4i = vec_ab[1]*vec_ac[2] - vec_ab[2]*vec_ac[1];
      Real a4j = vec_ab[2]*vec_ac[0] - vec_ab[0]*vec_ac[2];
      Real a4k = vec_ab[0]*vec_ac[1] - vec_ab[1]*vec_ac[0];
      Real mag = std::sqrt( a4i*a4i + a4j*a4j + a4k*a4k );
      if( mag > 0.0){
        mag = 1.0/mag;
        normals2[0+3*j] = a4i * mag;
        normals2[1+3*j] = a4j * mag;
        normals2[2+3*j] = a4k * mag;
      }
    }
  }
}

void ContactTriFaceQ6::FacetStaticRestriction(int nfacets, Real* coordinates, 
					      Real* normals, Real* ctrcl_facets,
					      Real* ctrcl)
{
  int ii1=0,ilocc=0,ilocs=0;
  int iistored=0,iconcave=0,iinside=1,iout=2;
  Real projcv,projmv,one_third=1.0/3.0;
  Real dctds,xmc,ymc,zmc,xms,yms,zms,vecmix,vecmiy,vecmiz;
  // Contract the closest point projection of the Node-with-N_Triangles
  // into a single contact with the master surface.
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
        // Determine if the two triangular surfaces are concave or convex
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

void ContactTriFaceQ6::FacetDynamicRestriction(int nfacets, 
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
  // contact with a master surface  
  //
  // Notes:  The memory is laid out so that we don't need to look at the first
  //         triangle.  If its not in contact it will be ignored above.  If it
  //         is in contact it should be left alone and will correctly competed.
  for (int ii=1; ii<nfacets; ++ii) {
    int index;    
    if( ctrcl_facets[ii*LENGTH+MSPARAM] == 1 ) { 
      if( ctrcl_facets[MSPARAM] != 1 ) {
        // No constraint is yet defined for this face so store this one
        for (index=0; index<LENGTH; ++index) {  
          ctrcl_facets[index] = ctrcl_facets[ii*LENGTH+index];
        }      
      } else {   
        // A constraint is already defined for this face.  Take the one
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

void ContactTriFaceQ6::Smooth_Normal( VariableHandle POSITION,
				      VariableHandle NODE_NORMAL,
				      VariableHandle FACE_NORMAL,
				      VariableHandle CURVATURE,
    		       ContactSearch::Smoothing_Resolution resolution,
				      Real percentage,
				      Real* coordinates,
				      Real* smooth_normal,
                                      Real critical_curvature )
{
  /*                       2
                          / \
                         x   x
                        /     \
                       x   G   x
                      / \     / \ 
                     /   x   x   \
                    /     \ /     \
                   /       c       \
                  /       / \       \
                 5       /   \       4
                /   D   /     \   C   \
               /       x       x       \
              /       /    A    \       \
             /       /           \       \
            /       /             \       \
           x---x---x-------w-------x---x---x
          /       /                 \       \ 
         x   E   x         B         x   F   x
        /       /                     \       \
       0---x---x-----------3-----------x---x---1
    
      Region    A: Use surface normal
              B-D: Smooth along one edge
              E-G: Smooth along two edges
  */
  Real upper_bound = 1.0 - percentage;
  Real lower_bound = percentage/2.0;
  Real* face_normal = Variable(FACE_NORMAL);
  Real node_position[8][3];
  Real node_normal[8][3];
  Real contact_point[3];
  Real coords[3];
  Real Mag;
  bool smooth_edge0,smooth_edge1;
  Real curvature0,curvature1;

  if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
      coordinates[1] >= lower_bound && coordinates[1] <= upper_bound &&
      coordinates[2] >= lower_bound && coordinates[2] <= upper_bound ){

    /* Region A
                         x
                        / \
                       /   \
                      /     \
                     x       x
                    /    A    \
                   /           \
                  /             \
                 x-------x-------x
    */

    smooth_normal[0] = face_normal[0];
    smooth_normal[1] = face_normal[1];
    smooth_normal[2] = face_normal[2];
    return;
  } 
  else if( coordinates[0] >= lower_bound && coordinates[1] >= lower_bound &&
	   coordinates[2] <  lower_bound){

    /* Region B
                 1-------4-------0
                /                 \
               5         B         7
              /                     \
             2-----------6-----------3
    */
    if( fabs(*Edge(0)->Variable(CURVATURE)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    } 
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = 1.0-lower_bound;
    coords[1] = lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[2][0]) );
    GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = lower_bound;
    coords[1] = 1.0-lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[3][0]) );
    GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
    // Compute the data for node 4
    coords[0] = 0.5-lower_bound/2.0;
    coords[1] = 0.5-lower_bound/2.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[4][0]) );
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute the data for node 5
    coords[0] = 1.0-1.5*lower_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[5][0]) );
    node_normal[5][0] = node_normal[1][0]+node_normal[2][0];
    node_normal[5][1] = node_normal[1][1]+node_normal[2][1];
    node_normal[5][2] = node_normal[1][2]+node_normal[2][2];
    Mag = node_normal[5][0]*node_normal[5][0] +
          node_normal[5][1]*node_normal[5][1] +
          node_normal[5][2]*node_normal[5][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0){
      Mag = 1.0/Mag;
      node_normal[5][0] *= Mag;
      node_normal[5][1] *= Mag;
      node_normal[5][2] *= Mag;
    }
    // Compute the data for node 6
    node_position[6][0] = Node(3)->Variable(POSITION)[0];
    node_position[6][1] = Node(3)->Variable(POSITION)[1];
    node_position[6][2] = Node(3)->Variable(POSITION)[2];
    GetEdgeSmoothedNormal(0, &(node_normal[6][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute the data for node 7
    coords[0] = lower_bound;
    coords[1] = 1.0-1.5*lower_bound;
    coords[2] = lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[7][0]) );
    node_normal[7][0] = node_normal[0][0]+node_normal[3][0];
    node_normal[7][1] = node_normal[0][1]+node_normal[3][1];
    node_normal[7][2] = node_normal[0][2]+node_normal[3][2];
    Mag = node_normal[7][0]*node_normal[7][0] +
          node_normal[7][1]*node_normal[7][1] +
          node_normal[7][2]*node_normal[7][2];
    Mag = std::sqrt(Mag);
    if( Mag > 0.0){
      Mag = 1.0/Mag;
      node_normal[7][0] *= Mag;
      node_normal[7][1] *= Mag;
      node_normal[7][2] *= Mag;
    }
  }
  else if( coordinates[0] <  lower_bound && coordinates[1] >= lower_bound &&
	   coordinates[2] >= lower_bound ){

    /* Region C
                             3 
                            / \ 
                             7   \
                            /     \
                           0       \
                            \       \
                             \       6
                              \   C   \
                               4       \
                                \       \
                                 \       \
                                  \       \
                                   1---5---2
    */
    if( fabs(*Edge(1)->Variable(CURVATURE)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    } 
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[2][0]) );
    GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = 0.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[3][0]) );
    GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );

    // Compute the data for node 4
    coords[0] = lower_bound;
    coords[1] = 0.5-lower_bound/2.0;
    coords[2] = 0.5-lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[4][0]) );
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute the data for node 5
    coords[0] = lower_bound/2.0;
    coords[1] = 1.0-1.5*lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[5][0]) );
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
    // Compute the data for node 6
    node_position[6][0] = Node(4)->Variable(POSITION)[0];
    node_position[6][1] = Node(4)->Variable(POSITION)[1];
    node_position[6][2] = Node(4)->Variable(POSITION)[2];
    GetEdgeSmoothedNormal(1, &(node_normal[6][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
     //                       &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute the data for node 7
    coords[0] = lower_bound/2.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-1.5*lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[7][0]) );
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
  }
  else if( coordinates[0] >= lower_bound && coordinates[1] <  lower_bound &&
	   coordinates[2] >= lower_bound ){

    /* Region D
                       2
                      / \
                     /   5
                    /     \
                   /       1
                  /       /
                 6       /
                /   D   /
               /       4
              /       /
             /       /
            /       /
           3---7---0
    */
    if( fabs(*Edge(2)->Variable(CURVATURE)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    } 
    // Compute the data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = lower_bound;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[2][0]) );
    GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = 1.0-lower_bound;
    coords[1] = 0.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[3][0]) );
    GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
    // Compute the data for node 4
    coords[0] = 0.5-lower_bound/2.0;
    coords[1] = lower_bound;
    coords[2] = 0.5-lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[4][0]) );
    node_normal[4][0] = face_normal[0];
    node_normal[4][1] = face_normal[1];
    node_normal[4][2] = face_normal[2];
    // Compute the data for node 5
    coords[0] = lower_bound;
    coords[1] = lower_bound/2.0;
    coords[2] = 1.0-1.5*lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[5][0]) );
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
    // Compute the data for node 6
    node_position[6][0] = Node(5)->Variable(POSITION)[0];
    node_position[6][1] = Node(5)->Variable(POSITION)[1];
    node_position[6][2] = Node(5)->Variable(POSITION)[2];
    GetEdgeSmoothedNormal(2, &(node_normal[6][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
     //                       &(node_position[6][0]), &(node_normal[6][0]) );
    // Compute the data for node 7
    coords[0] = 1.0-1.5*lower_bound;
    coords[1] = lower_bound/2.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[7][0]) );
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
  }
  else if( coordinates[0] > upper_bound && coordinates[1] < lower_bound &&
	   coordinates[2] < lower_bound ){

    /* Region E
             1---4---0
            /       /
           5   E   7
          /       /
         2---6---3
    */

    curvature0 = GetEdgeCurvature(2);
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
    // Compute the data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = 1.0-lower_bound;
    coords[1] = 0.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(2, &(node_normal[1][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute the data for node 2
    node_position[2][0] = Node(0)->Variable(POSITION)[0];
    node_position[2][1] = Node(0)->Variable(POSITION)[1];
    node_position[2][2] = Node(0)->Variable(POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(0)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(0)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(0)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } 
    // Compute the data for node 3
    coords[0] = 1.0-lower_bound;
    coords[1] = lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
    // Compute the data for node 4
    coords[0] = 1.0-1.5*lower_bound;
    coords[1] = lower_bound/2.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[4][0]) );
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
    // Compute the data for node 5
    coords[0] = 1.0-lower_bound/2.0;
    coords[1] = 0.0;
    coords[2] = lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[5][0]) );
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
    // Compute the data for node 6
    coords[0] = 1.0-lower_bound/2.0;
    coords[1] = lower_bound/2.0;
    coords[2] = 0.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[6][0]) );
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
    // Compute the data for node 7
    coords[0] = 1.0-1.5*lower_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[7][0]) );
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
  }
  else if( coordinates[0] < lower_bound && coordinates[1] > upper_bound &&
	   coordinates[2] < lower_bound ){

    /* Region F
                                 0---7---3
                                  \       \
                                   4   F   6
                                    \       \
                                     1---5---2
    */

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
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = 1.0-lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(0, &(node_normal[1][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute the data for node 2
    node_position[2][0] = Node(1)->Variable(POSITION)[0];
    node_position[2][1] = Node(1)->Variable(POSITION)[1];
    node_position[2][2] = Node(1)->Variable(POSITION)[2];
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
    // Compute the data for node 3
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
    // Compute the data for node 4
    coords[0] = lower_bound;
    coords[1] = 1.0-1.5*lower_bound;
    coords[2] = lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[4][0]) );
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
    // Compute the data for node 5
    coords[0] = lower_bound/2.0;
    coords[1] = 1.0-lower_bound/2.0;
    coords[2] = 0.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[5][0]) );
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
    // Compute the data for node 6
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound/2.0;
    coords[2] = lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[6][0]) );
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
    // Compute the data for node 7
    coords[0] = lower_bound/2.0;
    coords[1] = 1.0-1.5*lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[7][0]) );
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
  }
  else if( coordinates[0] < lower_bound && coordinates[1] < lower_bound &&
	   coordinates[2] > upper_bound ){

    /* Region G
                           2
                          / \
                         6   5
                        /     \
                       3   G   1 
                        \     /
                         7   4
                          \ /
                           0
    */

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
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = 0.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(1, &(node_normal[1][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute the data for node 2
    node_position[2][0] = Node(2)->Variable(POSITION)[0];
    node_position[2][1] = Node(2)->Variable(POSITION)[1];
    node_position[2][2] = Node(2)->Variable(POSITION)[2];
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
    // Compute the data for node 3
    coords[0] = lower_bound;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
    // Compute the data for node 4
    coords[0] = lower_bound/2.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-1.5*lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[4][0]) );
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
    // Compute the data for node 5
    coords[0] = 0.0;
    coords[1] = lower_bound/2.0;
    coords[2] = 1.0-lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[5][0]) );
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
    // Compute the data for node 6
    coords[0] = lower_bound/2.0;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound/2.0;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[6][0]) );
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
    // Compute the data for node 7
    coords[0] = lower_bound;
    coords[1] = lower_bound/2.0;
    coords[2] = 1.0-1.5*lower_bound;
    Compute_Global_Coordinates( POSITION, coords, &(node_position[7][0]) );
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
  }

  // Get the global coordinates of the contact point
  Compute_Global_Coordinates( POSITION, coordinates, contact_point );

  // Compute Contact Point in sub area
  Real local_coords[3];
  ContactQuadFaceQ8::Compute_Quad_Local_Coords( node_position, contact_point, local_coords);

  // Now interpolate to get the normal
  ContactQuadFaceQ8::Interpolate_Vector( local_coords, node_normal, smooth_normal );

  // Now make smooth_normal a unit vector
  Normalize(smooth_normal);
}

void ContactTriFaceQ6::Compute_Node_Areas( VariableHandle POSITION,
                                           VariableHandle /*FACE_NORMAL*/,
                                           Real* node_areas )
{
  // For the time being, I'm going to compute the area for the 3 node tri
  // and then distribute the area to all 6 nodes.

  // Get the nodal positions
  Real* pos_0 = Node(0)->Variable(POSITION);
  Real* pos_1 = Node(1)->Variable(POSITION);
  Real* pos_2 = Node(2)->Variable(POSITION);

  // Construct vectors between the nodes
  Real vec_0_1[3], vec_0_2[3];
  vec_0_1[0] = pos_1[0] - pos_0[0];
  vec_0_1[1] = pos_1[1] - pos_0[1];
  vec_0_1[2] = pos_1[2] - pos_0[2];
  vec_0_2[0] = pos_2[0] - pos_0[0];
  vec_0_2[1] = pos_2[1] - pos_0[1];
  vec_0_2[2] = pos_2[2] - pos_0[2];

  // Take the cross product
  Real cross[3];
  Cross(vec_0_1, vec_0_2, cross);

  // The total area is then the 1/2 the magnitude of the vector.
  // The node area is then 1/3 of the total area.
  Real node_area = Magnitude(cross) / 6.0;

  node_areas[0] = node_area;
  node_areas[1] = node_area;
  node_areas[2] = node_area;

  // Now distribute the area to all six nodes
  node_areas[3]  = 0.25*( node_areas[0] + node_areas[1] );
  node_areas[4]  = 0.25*( node_areas[1] + node_areas[2] );
  node_areas[5]  = 0.25*( node_areas[2] + node_areas[0] );
  node_areas[0] *= 0.5;
  node_areas[1] *= 0.5;
  node_areas[2] *= 0.5;
}

int ContactTriFaceQ6::FaceEdge_Intersection(VariableHandle POSITION,
					    ContactEdge<Real>* edge, Real* coords)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactTriFaceQ6::FaceEdge_Intersection not yet implemented\n";
#endif
  POSTCONDITION( 0 );
  return 0;
}

void ContactTriFaceQ6::Evaluate_Shape_Functions( Real* local_coords,
						 Real* shape_functions )
{
  Compute_Shape_Functions( local_coords, shape_functions );
}

void ContactTriFaceQ6::Compute_Global_Coordinates( VariableHandle POSITION,
                                                   Real* local_coords,
                                                   Real* global_coords )
{
  Real node_positions[6][3];
  for(int i=0; i<6; ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactTriFaceQ6::Compute_Local_Coordinates( Real Config_Param,
						  VariableHandle POSITION0,
						  VariableHandle POSITION1,
						  VariableHandle FACE_NORMAL,
						  Real* global_coords,
						  Real* local_coords )
{
  int  i, j;
  Real node_positions[6][3];
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

void ContactTriFaceQ6::Compute_Local_Coordinates( VariableHandle POSITION,
						  Real* global_coords,
						  Real* local_coords )
{
  int  i, j;
  Real node_positions[6][3];
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
/* on a linear T3 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary T3 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void ContactTriFaceQ6::Compute_Shape_Functions( Real local_coords[3],
						Real shape_functions[6] )
{
  PRECONDITION( std::fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Real a0 = local_coords[0];
  Real a1 = local_coords[1];
  Real a2 = local_coords[2];
  shape_functions[0] = a0*(2.0*a0-1.0);
  shape_functions[1] = a1*(2.0*a1-1.0);
  shape_functions[2] = a2*(2.0*a2-1.0);
  shape_functions[3] = 4.0*a0*a1;
  shape_functions[4] = 4.0*a1*a2;
  shape_functions[5] = 4.0*a2*a0;
}

void ContactTriFaceQ6::Compute_Shape_Derivatives( Real local_coords[3],
                                                  Real shape_derivatives[2][6])
{
  PRECONDITION( std::fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Real a0 = local_coords[0];
  Real a1 = local_coords[1];
  Real a2 = local_coords[2];
  shape_derivatives[0][0] =  4.0*a0 - 1.0;
  shape_derivatives[0][1] =  0.0;
  shape_derivatives[0][2] = -4.0*a2 + 1.0;
  shape_derivatives[0][3] =  4.0*a1;
  shape_derivatives[0][4] = -4.0*a1;
  shape_derivatives[0][5] =  4.0*(a2-a0);

  shape_derivatives[1][0] =  0.0;
  shape_derivatives[1][1] =  4.0*a1 - 1.0;
  shape_derivatives[1][2] = -4.0*a2 + 1.0;
  shape_derivatives[1][3] =  4.0*a0;
  shape_derivatives[1][4] =  4.0*(a2-a1);
  shape_derivatives[1][5] = -4.0*a0;
}

void ContactTriFaceQ6::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
                                             Real global_coords[3],
                                             Real local_coords[3] )
{
  int  i, j;
  int  nnodes=6;
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
    local_coords[0] = 1.0;
    local_coords[1] = 0.0;
    local_coords[2] = 0.0;
    break;
  case 1:
    local_coords[0] = 0.0;
    local_coords[1] = 1.0;
    local_coords[2] = 0.0;
    break;
  case 2:
    local_coords[0] = 0.0;
    local_coords[1] = 0.0;
    local_coords[2] = 1.0;
    break;
  case 3:
    local_coords[0] = 0.5;
    local_coords[1] = 0.5;
    local_coords[2] = 0.0;
    break;
  case 4:
    local_coords[0] = 0.0;
    local_coords[1] = 0.5;
    local_coords[2] = 0.5;
    break;
  case 5:
    local_coords[0] = 0.5;
    local_coords[1] = 0.0;
    local_coords[2] = 0.5;
    break;
  }
  if (i<nnodes) return;
  /*
     else use newton's method to iterate
    
     Assume straight-sided triangle for first guess
    
                      2
                     /|\
                    / | \
                   /a1|a0\
                  /  /c\  \
                 / /     \ \
                //    a2   \\ 
               0-------------1
  */

  // Compute vectors from the contact point to each node
  Real v_c_0[3] = { node_positions[0][0] - global_coords[0],
                    node_positions[0][1] - global_coords[1],
                    node_positions[0][2] - global_coords[2] };

  Real v_c_1[3] = { node_positions[1][0] - global_coords[0],
                    node_positions[1][1] - global_coords[1],
                    node_positions[1][2] - global_coords[2] };

  Real v_c_2[3] = { node_positions[2][0] - global_coords[0],
                    node_positions[2][1] - global_coords[1],
                    node_positions[2][2] - global_coords[2] };

  // Compute the areas
  Real shape_derivatives[2][6];
  Real midpnt[3] = {1.0/3.0, 1.0/3.0, 1.0 - 1.0/3.0 - 1.0/3.0};
  Compute_Shape_Derivatives( midpnt, shape_derivatives );
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  for (i=0; i<nnodes; ++i) {
    e1[0] += shape_derivatives[0][i]*node_positions[i][0];
    e1[1] += shape_derivatives[0][i]*node_positions[i][1];
    e1[2] += shape_derivatives[0][i]*node_positions[i][2];
    e2[0] += shape_derivatives[1][i]*node_positions[i][0];
    e2[1] += shape_derivatives[1][i]*node_positions[i][1];
    e2[2] += shape_derivatives[1][i]*node_positions[i][2];
  }
  Real a = Dot(e1,e1);
  Real b = Dot(e2,e2);
  Real c = Dot(e1,e2);
  Real face_normal[3];
  Real invDetJ = 1.0/std::sqrt(a*b-c*c);
  Cross(e1,e2,face_normal);
  Scale(face_normal,invDetJ);

  Real a0 = ScalarTripleProduct(v_c_1, v_c_2, face_normal);
  Real a1 = ScalarTripleProduct(v_c_2, v_c_0, face_normal);
  Real a2 = ScalarTripleProduct(v_c_0, v_c_1, face_normal);

  // Compute the local coordinates
  Real area = a0 + a1 + a2;
  a0 /= area;
  a1 /= area;
  a2 /= area;

  int  iterations=0;
  int  max_iterations=200;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real s, s0=a0, s1, ds=0.0;
  Real t, t0=a1, t1, dt=0.0;
  Real coords[3];
  Real J[3][2], f[3];
  while (!converged && iterations<max_iterations) {
    coords[0] = s0;
    coords[1] = t0;
    coords[2] = 1.0-s0-t0;
    Compute_Global_Coords(node_positions, coords, f );
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
    std::cerr << "ContactTriFaceQ6::Compute_Local_Coords() did not converge"
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
    std::cerr << "  After " << iterations << "iterations, local_coords = ("
	 << s0 << "," << t0 << ")" << std::endl;
    std::cerr << "  Going to continue processing anyway!!!" << std::endl;
  }
#endif
  // If it's close to any of the edges, snap to it
  if (s0<1.0+spatial_tolerance) {
    s0 = std::min(s0, 1.0);
  }
  if (s0>-spatial_tolerance) {
    s0 = std::max(s0, 0.0);
  }
  if (t0<1.0+spatial_tolerance) {
    t0 = std::min(t0, 1.0);
  }
  if (t0>-spatial_tolerance) {
    t0 = std::max(t0, 0.0);
  }
  local_coords[0] = s0;
  local_coords[1] = t0;
  local_coords[2] = 1.0-s0-t0;
}

void ContactTriFaceQ6::Compute_Global_Coords( Real node_positions[6][3],
					      Real local_coords[3],
					      Real global_coords[3] )
{
  Real N[6];
  int  nnodes=6;
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

void ContactTriFaceQ6::Interpolate_Scalar( Real  local_coords[3],
					   Real  node_scalars[6],
					   Real& interpolated_scalar )
{
  Real N[6];
  int  nnodes=6;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void ContactTriFaceQ6::Interpolate_Vector( Real local_coords[3],
					   Real node_vectors[6][3],
					   Real interpolated_vector[3] )
{
  Real N[6];
  int  nnodes=6;
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

