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


#ifndef ContactQuadFaceL4_C_
#define ContactQuadFaceL4_C_

#include "allocators.h"
#include "ContactNode.h"
#include "ContactQuadFaceL4.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactUtilities.h"
#include "contact_tolerances.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

template<typename DataType>
ContactQuadFaceL4<DataType>::ContactQuadFaceL4( ContactFixedSizeAllocator* alloc,
                                      int Block_Index, 
				      int Index_in_Block, int key ) 
  : ContactFace<DataType>( alloc, ContactSearch::QUADFACEL4,
                 Block_Index, Index_in_Block, key, 
                 nodes, edges, Node_Info, Edge_Info)
{}

template<typename DataType>
ContactQuadFaceL4<DataType>* ContactQuadFaceL4<DataType>::new_ContactQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactQuadFaceL4].New_Frag())
             ContactQuadFaceL4<DataType>(alloc, Block_Index, Index_in_Block, key);
}

template<typename DataType>
void ContactQuadFaceL4_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactQuadFaceL4<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactQuadFaceL4<DataType> allocator" );
}

template<typename DataType>
ContactQuadFaceL4<DataType>::~ContactQuadFaceL4() {}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Normal(VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL )
{
  using std::sqrt;
  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);
  ContactNode<DataType>* node3 = Node(3);
  DataType* Position0 = node0->Variable(CURRENT_POSITION);
  DataType* Position1 = node1->Variable(CURRENT_POSITION);
  DataType* Position2 = node2->Variable(CURRENT_POSITION);
  DataType* Position3 = node3->Variable(CURRENT_POSITION);
  // Compute the vector from node 0 to node 2 & from node 1 to node 3
  DataType Vec02[3],Vec13[3];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];
  Vec13[0] = Position3[0] - Position1[0];
  Vec13[1] = Position3[1] - Position1[1];
  Vec13[2] = Position3[2] - Position1[2];
  // Compute the face normal as the cross product of the two vectors
  DataType* face_normal = Variable(FACE_NORMAL);
  face_normal[0] = Vec02[1]*Vec13[2] - Vec02[2]*Vec13[1];
  face_normal[1] = Vec02[2]*Vec13[0] - Vec02[0]*Vec13[2];
  face_normal[2] = Vec02[0]*Vec13[1] - Vec02[1]*Vec13[0];
  DataType Mag = sqrt( face_normal[0]*face_normal[0] + 
		       face_normal[1]*face_normal[1] +
		       face_normal[2]*face_normal[2] );
  if( Mag > 0.0){
    Mag = 1.0/Mag;
    face_normal[0] *= Mag;
    face_normal[1] *= Mag;
    face_normal[2] *= Mag;
  }
}

#include <Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.h>
template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Partial_Face_Normal(VariableHandle CURRENT_POSITION, DataType (*dface_normal)[3])
{
  // XXX note: the definition of the "face normal" in this function is the normal to the face at the origin of the local frame
  // which is slightly different from the Compute_Normal function defined above.
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(CURRENT_POSITION)[0],
                                             Node(i)->Variable(CURRENT_POSITION)[1],
                                             Node(i)->Variable(CURRENT_POSITION)[2]);
  }

  int nodes[4] = { 0, 1, 2, 3 };
  DataType local_coords[2] = { 0., 0. };
  FaceQuad4 face(nodes);
  face.GetdUnitNormal<DataType,CoordSetType>(dface_normal, local_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Second_Partial_Face_Normal(VariableHandle CURRENT_POSITION, DataType (*d2face_normal)[3] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(CURRENT_POSITION)[0],
                                             Node(i)->Variable(CURRENT_POSITION)[1],
                                             Node(i)->Variable(CURRENT_POSITION)[2]);
  }

  int nodes[4] = { 0, 1, 2, 3 };
  DataType local_coords[2] = { 0., 0. };
  FaceQuad4 face(nodes);
  face.Getd2UnitNormal<DataType,CoordSetType>(d2face_normal, local_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Normal(VariableHandle POSITION,
                                       DataType* normal, DataType* local_coords )
{
  using std::sqrt;
  DataType shape_derivatives[2][4];
  DataType e1[3] = {0.0, 0.0, 0.0};
  DataType e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<4; ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  DataType  a=0.0, b=0.0, c=0.0;
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  DataType detJ = sqrt(a*b-c*c);
  normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])/detJ;
  normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])/detJ;
  normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])/detJ;
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Normal(DataType** nodal_positions,
                                       DataType* local_coords, DataType* normal )
{
  using std::sqrt;
  DataType shape_derivatives[2][4];
  DataType e1[3] = {0.0, 0.0, 0.0};
  DataType e2[3] = {0.0, 0.0, 0.0};
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<4; ++i) {
    DataType* node_position = nodal_positions[i];
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  DataType  a=0.0, b=0.0, c=0.0;
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  DataType detJ = sqrt(a*b-c*c);
  normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])/detJ;
  normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])/detJ;
  normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])/detJ;
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_CharacteristicLength(VariableHandle CURRENT_POSITION,
				                     VariableHandle CHARACTERISTIC_LENGTH )
{
  using std::sqrt;
  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);
  ContactNode<DataType>* node3 = Node(3);
  DataType* Position0 = node0->Variable(CURRENT_POSITION);
  DataType* Position1 = node1->Variable(CURRENT_POSITION);
  DataType* Position2 = node2->Variable(CURRENT_POSITION);
  DataType* Position3 = node3->Variable(CURRENT_POSITION);
  // Compute the vector from node 0 to node 2 & from node 1 to node 3
  DataType Vec02[3],Vec13[3];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];
  Vec13[0] = Position3[0] - Position1[0];
  Vec13[1] = Position3[1] - Position1[1];
  Vec13[2] = Position3[2] - Position1[2];
  // Compute the face normal as the cross product of the two vectors
  DataType face_normal[3];
  face_normal[0] = Vec02[1]*Vec13[2] - Vec02[2]*Vec13[1];
  face_normal[1] = Vec02[2]*Vec13[0] - Vec02[0]*Vec13[2];
  face_normal[2] = Vec02[0]*Vec13[1] - Vec02[1]*Vec13[0];

  DataType* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  *characteristiclength = sqrt( face_normal[0]*face_normal[0] + 
		                face_normal[1]*face_normal[1] +
		                face_normal[2]*face_normal[2] ) / 2.0;
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Centroid( VariableHandle CURRENT_POSITION,
					  VariableHandle CENTROID )
{
  // For the time being, compute the centroid as the average of the four nodes
  DataType* centroid = Variable(CENTROID);
  DataType* node_position = Node(0)->Variable(CURRENT_POSITION);
  centroid[0] = node_position[0];
  centroid[1] = node_position[1];
  centroid[2] = node_position[2];
  for( int i=1 ; i<4 ; ++i ){
    node_position = Node(i)->Variable(CURRENT_POSITION);
    centroid[0] += node_position[0];
    centroid[1] += node_position[1];
    centroid[2] += node_position[2];
  }
  centroid[0] *= 0.25;
  centroid[1] *= 0.25;
  centroid[2] *= 0.25;
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Get_Edge_Nodes( int i, ContactNode<DataType>** node )
{
  PRECONDITION( i>=0 && i<4 );
  switch( i ){
  case 0:
    node[0] = Node(0);
    node[1] = Node(1);
    break;
  case 1:
    node[0] = Node(1);
    node[1] = Node(2);
    break;
  case 2:
    node[0] = Node(2);
    node[1] = Node(3);
    break;
  case 3:
    node[0] = Node(3);
    node[1] = Node(0);
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Invalid Edge request in ContactQuadFaceL4<DataType>::Get_Edge_Nodes" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

template<typename DataType>
int ContactQuadFaceL4<DataType>::Get_Edge_Number( ContactNode<DataType>** edge_nodes )
{
  PRECONDITION( edge_nodes[0] && edge_nodes[1] );
  PRECONDITION( edge_nodes[0] != edge_nodes[1] );
  int i1=-1,i2=-1;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    ContactNode<DataType>* node = Node(i);
    if( edge_nodes[0] == node ) i1 = i;
    if( edge_nodes[1] == node ) i2 = i;
  }
  PRECONDITION( 0<=i1 && i1<4 );
  PRECONDITION( 0<=i2 && i2<4 );
  switch( i1 ){
  case 0:
    if( i2 == 1 ) return(0);
    if( i2 == 3 ) return(3);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 1:
    if( i2 == 0 ) return(0);
    if( i2 == 2 ) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 2:
    if( i2 == 1 ) return(1);
    if( i2 == 3 ) return(2);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 3:
    if( i2 == 2 ) return(2);
    if( i2 == 0 ) return(3);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  }
  POSTCONDITION( 0 );
  return( -1 );
}

template<typename DataType>
int ContactQuadFaceL4<DataType>::Get_Edge_Number( DataType* local_coords )
{
  using std::abs;
  DataType tol(1.0e-8);
  // PJSA: the explicit typecast of abs argument is a workaround of flaw in Eigen's AutoDiffScalar
  if (abs(DataType(local_coords[1] - -1.0))<tol && 
      abs(DataType(local_coords[3] - -1.0))<tol) return(0);
  if (abs(DataType(local_coords[0] -  1.0))<tol && 
      abs(DataType(local_coords[2] -  1.0))<tol) return(1);
  if (abs(DataType(local_coords[1] -  1.0))<tol && 
      abs(DataType(local_coords[3] -  1.0))<tol) return(2);
  if (abs(DataType(local_coords[0] - -1.0))<tol && 
      abs(DataType(local_coords[2] - -1.0))<tol) return(3);
  return( -1 );
}

template<typename DataType>
void 
ContactQuadFaceL4<DataType>::Compute_Edge_Normal( VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					int Edge,
					DataType* edge_normal)
{
  ContactNode<DataType> *edge_node[2];
  Get_Edge_Nodes( Edge, edge_node );

  // Compute the tangent direction as a vector from node 1 to node 2
  DataType* position1 = edge_node[0]->Variable(POSITION);
  DataType* position2 = edge_node[1]->Variable(POSITION);
  DataType edge_tangent[3];
  edge_tangent[0] = position2[0] - position1[0];
  edge_tangent[1] = position2[1] - position1[1];
  edge_tangent[2] = position2[2] - position1[2];

  // Compute the normal direction as the cross product of the edge tangent
  // and the face normal.
  DataType* f_normal = Variable(FACE_NORMAL);
  edge_normal[0] = edge_tangent[1]*f_normal[2] - edge_tangent[2]*f_normal[1];
  edge_normal[1] = edge_tangent[2]*f_normal[0] - edge_tangent[0]*f_normal[2];
  edge_normal[2] = edge_tangent[0]*f_normal[1] - edge_tangent[1]*f_normal[0];

  acme::Normalize(edge_normal);
}



template<typename DataType>
bool ContactQuadFaceL4<DataType>::Is_Inside_Face( DataType* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 )
    return true;
  return false;
}

template<typename DataType>
ContactFace<DataType>* ContactQuadFaceL4<DataType>::Neighbor( DataType* local_coords )
{
  bool off_edge_0 = local_coords[1] < -1 ? true : false;
  bool off_edge_1 = local_coords[0] >  1 ? true : false;
  bool off_edge_2 = local_coords[1] >  1 ? true : false;
  bool off_edge_3 = local_coords[0] < -1 ? true : false;
  //
  // Handle the cases with no amibiguity (i.e., only off one edge )
  if( off_edge_0 && !off_edge_1 && !off_edge_3 ) 
    return edges[0]->Neighbor_Face( this ); 
  if( off_edge_1 && !off_edge_0 && !off_edge_2 ) 
    return edges[1]->Neighbor_Face( this );
  if( off_edge_2 && !off_edge_1 && !off_edge_3 )
    return edges[2]->Neighbor_Face( this );
  if( off_edge_3 && !off_edge_2 && !off_edge_0 )
    return edges[3]->Neighbor_Face( this );
  // Now handle the case of two edges
  //    Take the neighbor we are furthest within 
  if( off_edge_0 && off_edge_3 ){
    if( local_coords[0] < local_coords[1] ){
      // take the face off edge 3
      return edges[3]->Neighbor_Face( this );
    } else {
      // take the face off edge 0
      return edges[0]->Neighbor_Face( this );
    }
  }
  if( off_edge_0 && off_edge_1 ){
    if( local_coords[0] < -1.0*local_coords[1] ){
      // take the face off edge 0
      return edges[0]->Neighbor_Face( this );
    } else {
      // take the face off edge 1
      return edges[1]->Neighbor_Face( this );
    }
  }
  if( off_edge_1 && off_edge_2 ){
    if( local_coords[0] > local_coords[1] ){
      // take the face off edge 1
      return edges[1]->Neighbor_Face( this );
    } else {
      // take the face off edge 2
      return edges[2]->Neighbor_Face( this );
    }
  }
  if( off_edge_2 && off_edge_3 ){
    if( local_coords[0] < -1.0*local_coords[1] ){
      // take the face off edge 3
      return edges[3]->Neighbor_Face( this );
    } else {
      // take the face off edge 2
      return edges[2]->Neighbor_Face( this );
    }
  }
  // we should never get here.
  POSTCONDITION( 0 );
  return (ContactFace<DataType>*) NULL;
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Get_Close_Edges( DataType* local_coords, int& number,
					 int& edge_1_id, int& edge_2_id ){
  DataType  xi = local_coords[0];
  DataType eta = local_coords[1];

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

template<typename DataType>
bool ContactQuadFaceL4<DataType>::IsPlanar(VariableHandle POSITION)
{
  using std::sqrt;
  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);
  ContactNode<DataType>* node3 = Node(3);
  DataType* Position0 = node0->Variable(POSITION);
  DataType* Position1 = node1->Variable(POSITION);
  DataType* Position2 = node2->Variable(POSITION);
  DataType* Position3 = node3->Variable(POSITION);
  // Compute the vector from node 0 to node 1, from node 0 to node 3,
  // from node 2 to node 1, and node 2 to node 3
  DataType Vec01[3],Vec03[3],Vec23[3],Vec21[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec03[0] = Position3[0] - Position0[0];
  Vec03[1] = Position3[1] - Position0[1];
  Vec03[2] = Position3[2] - Position0[2];
  Vec23[0] = Position3[0] - Position2[0];
  Vec23[1] = Position3[1] - Position2[1];
  Vec23[2] = Position3[2] - Position2[2];
  Vec21[0] = Position1[0] - Position2[0];
  Vec21[1] = Position1[1] - Position2[1];
  Vec21[2] = Position1[2] - Position2[2];
  // Compute the face normal as the cross product of the two vectors
  DataType normal0[3], normal1[3], Mag;
  normal0[0] = Vec01[1]*Vec03[2] - Vec01[2]*Vec03[1];
  normal0[1] = Vec01[2]*Vec03[0] - Vec01[0]*Vec03[2];
  normal0[2] = Vec01[0]*Vec03[1] - Vec01[1]*Vec03[0];
  normal1[0] = Vec23[1]*Vec21[2] - Vec23[2]*Vec21[1];
  normal1[1] = Vec23[2]*Vec21[0] - Vec23[0]*Vec21[2];
  normal1[2] = Vec23[0]*Vec21[1] - Vec23[1]*Vec21[0];
  Mag = sqrt( normal0[0]*normal0[0] + 
	      normal0[1]*normal0[1] +
	      normal0[2]*normal0[2] );
  if( Mag > 0.0){
    Mag = 1.0/Mag;
    normal0[0] *= Mag;
    normal0[1] *= Mag;
    normal0[2] *= Mag;
  }
  Mag = sqrt( normal1[0]*normal1[0] + 
	      normal1[1]*normal1[1] +
	      normal1[2]*normal1[2] );
  if( Mag > 0.0){
    Mag = 1.0/Mag;
    normal1[0] *= Mag;
    normal1[1] *= Mag;
    normal1[2] *= Mag;
  }
  DataType dot = normal0[0]*normal1[0] + 
             normal0[1]*normal1[1] +
             normal0[2]*normal1[2];
  if (dot<=is_planar_tol) 
    return false;
  else 
    return true;
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::FacetDecomposition(int& nfacets,
          DataType* coordinates0, DataType* normals0, VariableHandle POSITION0,
          DataType* coordinates1, DataType* normals1, VariableHandle POSITION1,
          DataType* coordinates2, DataType* normals2, VariableHandle POSITION2) 
{
  using std::sqrt;
  int i, j;
  bool is_planar = IsPlanar(POSITION0);
  if (coordinates1!=NULL && is_planar) {
    is_planar = IsPlanar(POSITION1);
  }
  // if quad is planar then only decompose it into two triangles
  if (is_planar) {
    nfacets             = 2;
    DataType* Position0     = Node(0)->Variable(POSITION0);
    DataType* Position1     = Node(1)->Variable(POSITION0);
    DataType* Position2     = Node(2)->Variable(POSITION0);
    DataType* Position3     = Node(3)->Variable(POSITION0);
    coordinates0[0+3*0] = Position0[0];
    coordinates0[1+3*0] = Position0[1];
    coordinates0[2+3*0] = Position0[2];
    coordinates0[0+3*1] = Position1[0];
    coordinates0[1+3*1] = Position1[1];
    coordinates0[2+3*1] = Position1[2];
    coordinates0[0+3*2] = Position2[0];
    coordinates0[1+3*2] = Position2[1];
    coordinates0[2+3*2] = Position2[2];
    coordinates0[0+3*3] = Position2[0];
    coordinates0[1+3*3] = Position2[1];
    coordinates0[2+3*3] = Position2[2];
    coordinates0[0+3*4] = Position3[0];
    coordinates0[1+3*4] = Position3[1];
    coordinates0[2+3*4] = Position3[2];
    coordinates0[0+3*5] = Position0[0];
    coordinates0[1+3*5] = Position0[1];
    coordinates0[2+3*5] = Position0[2];
    if (coordinates1!=NULL) {
      Position0           = Node(0)->Variable(POSITION1);
      Position1           = Node(1)->Variable(POSITION1);
      Position2           = Node(2)->Variable(POSITION1);
      Position3           = Node(3)->Variable(POSITION1);
      coordinates1[0+3*0] = Position0[0];
      coordinates1[1+3*0] = Position0[1];
      coordinates1[2+3*0] = Position0[2];
      coordinates1[0+3*1] = Position1[0];
      coordinates1[1+3*1] = Position1[1];
      coordinates1[2+3*1] = Position1[2];
      coordinates1[0+3*2] = Position2[0];
      coordinates1[1+3*2] = Position2[1];
      coordinates1[2+3*2] = Position2[2];
      coordinates1[0+3*3] = Position2[0];
      coordinates1[1+3*3] = Position2[1];
      coordinates1[2+3*3] = Position2[2];
      coordinates1[0+3*4] = Position3[0];
      coordinates1[1+3*4] = Position3[1];
      coordinates1[2+3*4] = Position3[2];
      coordinates1[0+3*5] = Position0[0];
      coordinates1[1+3*5] = Position0[1];
      coordinates1[2+3*5] = Position0[2];
    }
    if (coordinates2!=NULL) {
      Position0           = Node(0)->Variable(POSITION2);
      Position1           = Node(1)->Variable(POSITION2);
      Position2           = Node(2)->Variable(POSITION2);
      Position3           = Node(3)->Variable(POSITION2);
      coordinates2[0+3*0] = Position0[0];
      coordinates2[1+3*0] = Position0[1];
      coordinates2[2+3*0] = Position0[2];
      coordinates2[0+3*1] = Position1[0];
      coordinates2[1+3*1] = Position1[1];
      coordinates2[2+3*1] = Position1[2];
      coordinates2[0+3*2] = Position2[0];
      coordinates2[1+3*2] = Position2[1];
      coordinates2[2+3*2] = Position2[2];
      coordinates2[0+3*3] = Position2[0];
      coordinates2[1+3*3] = Position2[1];
      coordinates2[2+3*3] = Position2[2];
      coordinates2[0+3*4] = Position3[0];
      coordinates2[1+3*4] = Position3[1];
      coordinates2[2+3*4] = Position3[2];
      coordinates2[0+3*5] = Position0[0];
      coordinates2[1+3*5] = Position0[1];
      coordinates2[2+3*5] = Position0[2];
    }
  } else {
    nfacets = 4;
	// need to calculate the centroid locally because stored value is
	// only at CURRENT_POSITION and we don't know if CURRENT or PREDICTED
	// position is being asked for.
    DataType centroid[3];
    DataType* node_position = Node(0)->Variable(POSITION0);
    centroid[0] = node_position[0];
    centroid[1] = node_position[1];
    centroid[2] = node_position[2];
    for( i=1 ; i<4 ; ++i ){
      node_position = Node(i)->Variable(POSITION0);
      centroid[0] += node_position[0];
      centroid[1] += node_position[1];
      centroid[2] += node_position[2];
    }
    centroid[0] /= 4;
    centroid[1] /= 4;
    centroid[2] /= 4;
    for (i=0, j=1; i<4; ++i, j=(i+1)%4) {
      DataType* node_position1 = Node(i)->Variable(POSITION0);
      DataType* node_position2 = Node(j)->Variable(POSITION0);
      coordinates0[0+9*i]   = node_position1[0];
      coordinates0[1+9*i]   = node_position1[1];
      coordinates0[2+9*i]   = node_position1[2];
      coordinates0[3+9*i]   = node_position2[0];
      coordinates0[4+9*i]   = node_position2[1];
      coordinates0[5+9*i]   = node_position2[2];
      coordinates0[6+9*i]   = centroid[0];
      coordinates0[7+9*i]   = centroid[1];
      coordinates0[8+9*i]   = centroid[2];
    }
    if (coordinates1!=NULL) {
      node_position = Node(0)->Variable(POSITION1);
      centroid[0]   = node_position[0];
      centroid[1]   = node_position[1];
      centroid[2]   = node_position[2];
      for( i=1 ; i<4 ; ++i ){
        node_position = Node(i)->Variable(POSITION1);
        centroid[0]  += node_position[0];
        centroid[1]  += node_position[1];
        centroid[2]  += node_position[2];
      }
      centroid[0] /= 4;
      centroid[1] /= 4;
      centroid[2] /= 4;
      for (i=0, j=1; i<4; ++i, j=(i+1)%4) {
	DataType* node_position1  = Node(i)->Variable(POSITION1);
	DataType* node_position2  = Node(j)->Variable(POSITION1);
	coordinates1[0+9*i]   = node_position1[0];
	coordinates1[1+9*i]   = node_position1[1];
	coordinates1[2+9*i]   = node_position1[2];
	coordinates1[3+9*i]   = node_position2[0];
	coordinates1[4+9*i]   = node_position2[1];
	coordinates1[5+9*i]   = node_position2[2];
	coordinates1[6+9*i]   = centroid[0];
	coordinates1[7+9*i]   = centroid[1];
	coordinates1[8+9*i]   = centroid[2];
      }
    }
    if (coordinates2!=NULL) {
      node_position = Node(0)->Variable(POSITION2);
      centroid[0]   = node_position[0];
      centroid[1]   = node_position[1];
      centroid[2]   = node_position[2];
      for( i=1 ; i<4 ; ++i ){
        node_position = Node(i)->Variable(POSITION2);
        centroid[0]  += node_position[0];
        centroid[1]  += node_position[1];
        centroid[2]  += node_position[2];
      }
      centroid[0] /= 4;
      centroid[1] /= 4;
      centroid[2] /= 4;
      for (i=0, j=1; i<4; ++i, j=(i+1)%4) {
	DataType* node_position1 = Node(i)->Variable(POSITION2);
	DataType* node_position2 = Node(j)->Variable(POSITION2);
	coordinates2[0+9*i]  = node_position1[0];
	coordinates2[1+9*i]  = node_position1[1];
	coordinates2[2+9*i]  = node_position1[2];
	coordinates2[3+9*i]  = node_position2[0];
	coordinates2[4+9*i]  = node_position2[1];
	coordinates2[5+9*i]  = node_position2[2];
	coordinates2[6+9*i]  = centroid[0];
	coordinates2[7+9*i]  = centroid[1];
	coordinates2[8+9*i]  = centroid[2];
      }
    }
  }
  for (i=0; i<nfacets; ++i) {
    // compute vectors from centroid to other  
    // two nodes defining master surface triangle
    DataType vx0  = coordinates0[6+9*i] - coordinates0[0+9*i];
    DataType vy0  = coordinates0[7+9*i] - coordinates0[1+9*i];
    DataType vz0  = coordinates0[8+9*i] - coordinates0[2+9*i];
    DataType vx00 = coordinates0[3+9*i] - coordinates0[0+9*i];
    DataType vy00 = coordinates0[4+9*i] - coordinates0[1+9*i];
    DataType vz00 = coordinates0[5+9*i] - coordinates0[2+9*i];
    // i,j,k components of the normal to the plane
    DataType a4i        = -vy0*vz00+vy00*vz0;
    DataType a4j        =  vx0*vz00-vx00*vz0;
    DataType a4k        = -vx0*vy00+vx00*vy0;
    DataType snmag      = 1.0/sqrt(a4i*a4i+a4j*a4j+a4k*a4k);
    normals0[0+3*i] = a4i*snmag;
    normals0[1+3*i] = a4j*snmag;
    normals0[2+3*i] = a4k*snmag;
  }
  if (coordinates1!=NULL) {
    for (i=0; i<nfacets; ++i) {
      // compute vectors from centroid to other  
      // two nodes defining master surface triangle
      DataType vx0  = coordinates1[6+9*i] - coordinates1[0+9*i];
      DataType vy0  = coordinates1[7+9*i] - coordinates1[1+9*i];
      DataType vz0  = coordinates1[8+9*i] - coordinates1[2+9*i];
      DataType vx00 = coordinates1[3+9*i] - coordinates1[0+9*i];
      DataType vy00 = coordinates1[4+9*i] - coordinates1[1+9*i];
      DataType vz00 = coordinates1[5+9*i] - coordinates1[2+9*i];
      // i,j,k components of the normal to the plane
      DataType a4i        = -vy0*vz00+vy00*vz0;
      DataType a4j        =  vx0*vz00-vx00*vz0;
      DataType a4k        = -vx0*vy00+vx00*vy0;
      DataType snmag      = 1.0/sqrt(a4i*a4i+a4j*a4j+a4k*a4k);
      normals1[0+3*i] = a4i*snmag;
      normals1[1+3*i] = a4j*snmag;
      normals1[2+3*i] = a4k*snmag;
    }
  }
  if (coordinates2!=NULL) {
    for (i=0; i<nfacets; ++i) {
      // compute vectors from centroid to other  
      // two nodes defining master surface triangle
      DataType vx0  = coordinates2[6+9*i] - coordinates2[0+9*i];
      DataType vy0  = coordinates2[7+9*i] - coordinates2[1+9*i];
      DataType vz0  = coordinates2[8+9*i] - coordinates2[2+9*i];
      DataType vx00 = coordinates2[3+9*i] - coordinates2[0+9*i];
      DataType vy00 = coordinates2[4+9*i] - coordinates2[1+9*i];
      DataType vz00 = coordinates2[5+9*i] - coordinates2[2+9*i];
      // i,j,k components of the normal to the plane
      DataType a4i        = -vy0*vz00+vy00*vz0;
      DataType a4j        =  vx0*vz00-vx00*vz0;
      DataType a4k        = -vx0*vy00+vx00*vy0;
      DataType snmag      = 1.0/sqrt(a4i*a4i+a4j*a4j+a4k*a4k);
      normals2[0+3*i] = a4i*snmag;
      normals2[1+3*i] = a4j*snmag;
      normals2[2+3*i] = a4k*snmag;
    }
  }
}

template<>
inline void ContactQuadFaceL4<Real>::FacetStaticRestriction(int nfacets, Real* coordinates, 
                                         Real* normals, Real* ctrcl_facets, 
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
          }  else if( std::abs(ilocs) == iout ) {
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

template<typename DataType>
void ContactQuadFaceL4<DataType>::FacetStaticRestriction(int nfacets, DataType* coordinates,
                                         DataType* normals, DataType* ctrcl_facets,
                                         DataType* ctrcl)
{
  std::cerr << "ContactQuadFaceL4<DataType>::FacetStaticRestriction is not implemented\n";
  exit(-1);
}

template<>
inline void ContactQuadFaceL4<Real>::FacetDynamicRestriction(int nfacets, 
                                                Real* ctrcl_facets, 
                                                Real* ctrcl)
{
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
  for( int i=0 ; i<nfacets ; ++i ){
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
  //         is in contact it should be left alone and will correctly competed.
  for (int ii=1; ii<nfacets; ++ii) {
    if( ctrcl_facets[ii*LENGTH+MSPARAM] == 1 ) {
      if( ctrcl_facets[MSPARAM] != 1 ) {
        // No constraint is yet defined for this quad so store this one
	for (int index=0; index<LENGTH; ++index) {
          ctrcl_facets[index] = ctrcl_facets[ii*LENGTH+index];
	}
      } else {
        // A constraint is already defined for this quad.  Take the one
        // with the minimum contact time.
        if( ctrcl_facets[ii*LENGTH+ICTIMC] < ctrcl_facets[ICTIMC] ) {
	  if( ctrcl_facets[ii*LENGTH+ICTIMC] > 0.0 ) {
	    for (int index=0; index<LENGTH; ++index) {
              ctrcl_facets[index] = ctrcl_facets[ii*LENGTH+index];
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

template<typename DataType>
void
ContactQuadFaceL4<DataType>::FacetDynamicRestriction(int nfacets,
                                                DataType* ctrcl_facets,
                                                DataType* ctrcl)
{
  std::cerr << "ContactQuadFaceL4<DataType>::FacetDynamicRestriction is not implemented\n";
  exit(-1);
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Smooth_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle NODE_NORMAL,
				       VariableHandle FACE_NORMAL,
				       VariableHandle CURVATURE,
			ContactSearch::Smoothing_Resolution resolution,
				       DataType percentage,
				       DataType* coordinates,
				       DataType* smooth_normal,
                                       DataType critical_curvature )
{
  //
  //  3--k-------------------l--2
  //  | I|        D          |H |
  //  g--h-------------------i--j
  //  |  |                   |  |
  //  |  |                   |  |
  //  | E|        A          |C |
  //  |  |                   |  |
  //  |  |                   |  |
  //  c--d-------------------e--f
  //  | F|        B          |G |
  //  0--a-------------------b--1
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
  
  using std::abs;
  using std::sqrt;

  DataType upper_bound =  1.0 - percentage;
  DataType lower_bound = -1.0 + percentage;
  DataType* face_normal = Variable(FACE_NORMAL);
  DataType node_position[4][3];
  DataType node_normal[4][3];
  DataType contact_point[3];
  DataType coords[2];
  bool smooth_edge0,smooth_edge1;
  DataType curvature0,curvature1;

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
    //  1-------------------0
    //  |        B          |
    //  2-------------------3
    //

    if( abs(GetEdgeCurvature(0)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = lower_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = upper_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
    //			    &(node_position[3][0]), &(node_normal[3][0]) );
  }
  else if( coordinates[0] > upper_bound &&
	   coordinates[1] <= upper_bound && coordinates[1] >= lower_bound ){
    // Region C
    //
    //  0--3
    //  |  |
    //  |  |
    //  |C |
    //  |  |
    //  |  |
    //  |  |
    //  1--2
    //
    if(abs(GetEdgeCurvature(1)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }

    // Compute data for node 0
    coords[0] = upper_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = 1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
    //			    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = 1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
    //                        &(node_position[3][0]), &(node_normal[3][0]) );
  }
  else if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
	   coordinates[1] > upper_bound ){
    // Region D
    //
    //  3-------------------2
    //  |        D          |
    //  0-------------------1
    //
    if(abs(GetEdgeCurvature(2)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = upper_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = upper_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
    //			    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = lower_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
    //			    &(node_position[3][0]), &(node_normal[3][0]) );
  }
  else if( coordinates[0] < lower_bound &&
	   coordinates[1] <= upper_bound && coordinates[1] >= lower_bound ){
    // Region E
    //
    //  2--1
    //  |  |
    //  |  |
    //  | E|
    //  |  |
    //  |  |
    //  |  |
    //  3--0
    //
    if(abs(GetEdgeCurvature(3)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute data for node 2
    coords[0] = -1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(3, &(node_normal[2][0]));
    //Edge(3)->Smooth_Normal( FACE_NORMAL,
    //			    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute data for node 3
    coords[0] = -1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(3, &(node_normal[3][0]));
    //Edge(3)->Smooth_Normal( FACE_NORMAL,
    //			    &(node_position[3][0]), &(node_normal[3][0]) );
  }
  else if( coordinates[0] < lower_bound && coordinates[1] < lower_bound ){
    // Region F
    //
    //  1--0
    //  | F|
    //  2--3
    //
    curvature0 = GetEdgeCurvature(3);
    curvature1 = GetEdgeCurvature(0);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0)>critical_curvature) smooth_edge0 = false;
    if( abs(curvature1)>critical_curvature) smooth_edge1 = false;
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = -1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      this->GetEdgeSmoothedNormal(3, &(node_normal[1][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
      //		        &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute data for node 2
    node_position[2][0] = Node(0)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(0)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(0)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL ){
      node_normal[2][0]   = Node(0)->Variable(NODE_NORMAL)[0];
      node_normal[2][1]   = Node(0)->Variable(NODE_NORMAL)[1];
      node_normal[2][2]   = Node(0)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(3, &(node_normal[2][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } 
    // Compute data for node 3
    coords[0] = lower_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      this->GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
  }
  else if( coordinates[0] > upper_bound && coordinates[1] < lower_bound ){
    // Region G
    //
    //  0--3
    //  | G|
    //  1--2
    //
    curvature0 = GetEdgeCurvature(0);
    curvature1 = GetEdgeCurvature(1);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0)>critical_curvature) smooth_edge0 = false;
    if( abs(curvature1)>critical_curvature) smooth_edge1 = false;
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = upper_bound;
    coords[1] = -1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      this->GetEdgeSmoothedNormal(0, &(node_normal[1][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute data for node 2
    node_position[2][0] = Node(1)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(1)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(1)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1)  ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0]   = Node(1)->Variable(NODE_NORMAL)[0];
      node_normal[2][1]   = Node(1)->Variable(NODE_NORMAL)[1];
      node_normal[2][2]   = Node(1)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute data for node 3
    coords[0] = 1.0;
    coords[1] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      this->GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
  }
  else if( coordinates[0] > upper_bound && coordinates[1] > upper_bound ){
    // Region H
    //
    //  3--2
    //  | H|
    //  0--1
    //
    curvature0 = GetEdgeCurvature(1);
    curvature1 = GetEdgeCurvature(2);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0)>critical_curvature) smooth_edge0 = false;
    if( abs(curvature1)>critical_curvature) smooth_edge1 = false;
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = 1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      this->GetEdgeSmoothedNormal(1, &(node_normal[1][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute data for node 2
    node_position[2][0] = Node(2)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(2)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(2)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0]   = Node(2)->Variable(NODE_NORMAL)[0];
      node_normal[2][1]   = Node(2)->Variable(NODE_NORMAL)[1];
      node_normal[2][2]   = Node(2)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute data for node 3
    coords[0] = upper_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      this->GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }    
  }
  else if( coordinates[0] < lower_bound && coordinates[1] > upper_bound ){
    // Region I
    //
    //  2--1
    //  | I|
    //  3--0
    //
    curvature0 = GetEdgeCurvature(2);
    curvature1 = GetEdgeCurvature(3);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0)>critical_curvature) smooth_edge0 = false;
    if( abs(curvature1)>critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute data for node 1
    coords[0] = lower_bound;
    coords[1] = 1.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      this->GetEdgeSmoothedNormal(2, &(node_normal[1][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }    
    // Compute data for node 2
    node_position[2][0] = Node(3)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(3)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(3)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL ){
      node_normal[2][0]   = Node(3)->Variable(NODE_NORMAL)[0];
      node_normal[2][1]   = Node(3)->Variable(NODE_NORMAL)[1];
      node_normal[2][2]   = Node(3)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(3, &(node_normal[2][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute data for node 3
    coords[0] = -1.0;
    coords[1] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      this->GetEdgeSmoothedNormal(3, &(node_normal[3][0]));
      //Edge(3)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
  }

  // Compute Contact Point in sub area
  DataType local_coords[3];
  Compute_Local_Coords( node_position, contact_point,
			&(local_coords[0]));

  // Now interpolate to get the normal
  Interpolate_Vector( &(local_coords[0]), node_normal, 
				  smooth_normal );

  // Now make smooth_normal a unit vector
  DataType Mag = smooth_normal[0]*smooth_normal[0] + 
             smooth_normal[1]*smooth_normal[1] +
             smooth_normal[2]*smooth_normal[2];
  Mag = sqrt(Mag);
  if( Mag > 0.0){
    Mag = 1.0/Mag;
    smooth_normal[0] *= Mag;
    smooth_normal[1] *= Mag;
    smooth_normal[2] *= Mag;
  }
}

template<typename DataType>
void ContactQuadFaceL4<DataType>:: Compute_Node_Areas( VariableHandle POSITION, 
                                             VariableHandle FACE_NORMAL,
                                             DataType* node_areas )
{
  // Get the Nodal Positions
  DataType* pos_0 = Node(0)->Variable( POSITION );
  DataType* pos_1 = Node(1)->Variable( POSITION );
  DataType* pos_2 = Node(2)->Variable( POSITION );
  DataType* pos_3 = Node(3)->Variable( POSITION );

  // Get the face normal
  DataType* fnorm = Variable( FACE_NORMAL );

  // Compute vectors between nodes
  DataType vec_0_1[3], vec_2_1[3], vec_2_3[3], vec_3_0[3], vec_0_2[3], vec_3_1[3];
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
  DataType s1[3], s2[3], s3[3];
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
  DataType v0[3], v1[3], v2[3], v3[3];
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
}

template<typename DataType>
int ContactQuadFaceL4<DataType>::FaceEdge_Intersection(VariableHandle POSITION,
                                             ContactEdge<DataType>* edge, DataType* coords)
{
  using std::sqrt;
  using std::abs;
  int intersection=0;
#if 0
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactQuadFaceL4<DataType>::FaceEdge_Intersection not yet implemented\n";
#endif
  POSTCONDITION( 0 );
  return intersection;
#else
  // first do a bounding box check
  int i, j;
  DataType edge_min[3], edge_max[3];
  DataType face_min[3], face_max[3];
  DataType* node_position = Node(0)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    face_min[j] = node_position[j];
    face_max[j] = node_position[j];
  }
  for (i=1; i<Nodes_Per_Face(); ++i) {
    node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      face_min[j] = std::min(face_min[j],node_position[j]);
      face_max[j] = std::max(face_max[j],node_position[j]);
    }
  }
  node_position = edge->Node(0)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_min[j] = node_position[j];
    edge_max[j] = node_position[j];
  }
  node_position = edge->Node(1)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_min[j] = std::min(edge_min[j],node_position[j]);
    edge_max[j] = std::max(edge_max[j],node_position[j]);
  }
  if ( edge_min[0]>face_max[0] || 
       edge_min[1]>face_max[1] ||  
       edge_min[2]>face_max[2] ||
       edge_max[0]<face_min[0] || 
       edge_max[1]<face_min[1] || 
       edge_max[2]<face_min[2] )  return intersection;
  
  // Compute edge ray
  DataType dir_mag=0.0;
  DataType edge_pnt[3];
  DataType edge_dir[3];
  DataType* edge_node_position0 = edge->Node(0)->Variable(POSITION);
  DataType* edge_node_position1 = edge->Node(1)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_pnt[j] = edge_node_position0[j];
    edge_dir[j] = edge_node_position1[j]-edge_node_position0[j];
    dir_mag += edge_dir[j]*edge_dir[j];
  }
  dir_mag = 1.0/sqrt(dir_mag);
  edge_dir[0] *= dir_mag;
  edge_dir[1] *= dir_mag;
  edge_dir[2] *= dir_mag;
  DataType tmax;
  if (abs(edge_dir[0])>=abs(edge_dir[1]) && abs(edge_dir[0])>=abs(edge_dir[2])) {
    tmax = (edge_node_position1[0]-edge_node_position0[0])/edge_dir[0];
  } else
  if (abs(edge_dir[1])>=abs(edge_dir[0]) && abs(edge_dir[1])>=abs(edge_dir[2])) {
    tmax = (edge_node_position1[1]-edge_node_position0[1])/edge_dir[1];
  } else
  if (abs(edge_dir[2])>=abs(edge_dir[0]) && abs(edge_dir[2])>=abs(edge_dir[1])) {
    tmax = (edge_node_position1[2]-edge_node_position0[2])/edge_dir[2];
  }
/*
  tmax *= 1.05;
*/
  
  // subdivide the q4 into four tri3's
  DataType center[3];
  center[0] = center[1] = center[2] = 0.0;
  for (i=0; i<Nodes_Per_Face(); ++i) {
    node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      center[j] += node_position[j];
    }
  }
  center[0] *= 0.25;
  center[1] *= 0.25;
  center[2] *= 0.25;
  
  for (i=0; i<4; ++i) {
    int p0 = i;
    int p1 = (i+1)%4;
    DataType* node_position0 = center;
    DataType* node_position1 = Node(p0)->Variable(POSITION);
    DataType* node_position2 = Node(p1)->Variable(POSITION);
    // compute normal of tri3
    DataType normal[3], n_mag;
    DataType dx1  = node_position1[0] - node_position0[0];
    DataType dy1  = node_position1[1] - node_position0[1];
    DataType dz1  = node_position1[2] - node_position0[2];
    DataType dx2  = node_position2[0] - node_position0[0];
    DataType dy2  = node_position2[1] - node_position0[1];
    DataType dz2  = node_position2[2] - node_position0[2];
    normal[0] = dy1*dz2-dz1*dy2;
    normal[1] = dz1*dx2-dx1*dz2;
    normal[2] = dx1*dy2-dy1*dx2;
    n_mag     = 1.0/sqrt(normal[0]*normal[0]+
                         normal[1]*normal[1]+
                         normal[2]*normal[2]);
    normal[0] *= n_mag;
    normal[1] *= n_mag;
    normal[2] *= n_mag;
  
    // determine (edge ray)/(face plane) intersection
    DataType n_dot_d = normal[0]*edge_dir[0] +
                   normal[1]*edge_dir[1] +
                   normal[2]*edge_dir[2];
    // if (edge ray) and (tri3 plane) are parallel => no intersection
    if (abs(n_dot_d)<is_parallel_tol) continue;
    DataType q[3], P[3], t;
    q[0] = node_position0[0]-edge_pnt[0];
    q[1] = node_position0[1]-edge_pnt[1];
    q[2] = node_position0[2]-edge_pnt[2];
    t    = (normal[0]*q[0] + normal[1]*q[1] + normal[2]*q[2])/n_dot_d;
    if (t<0.0 || t>tmax) continue;
    P[0] = edge_pnt[0] + edge_dir[0]*t;
    P[1] = edge_pnt[1] + edge_dir[1]*t;
    P[2] = edge_pnt[2] + edge_dir[2]*t;
 
    // now determine if this point is inside the face
    DataType u0, u1, u2, v0, v1, v2;
    DataType nx = normal[0]*normal[0];
    DataType ny = normal[1]*normal[1];
    DataType nz = normal[2]*normal[2];
    if (nx>=ny && nx>=nz) {
      u0 = P[1] - node_position0[1];
      v0 = P[2] - node_position0[2];
      u1 = node_position1[1] - node_position0[1];
      u2 = node_position2[1] - node_position0[1];
      v1 = node_position1[2] - node_position0[2];
      v2 = node_position2[2] - node_position0[2];
    }
    else if (ny>=nx && ny>=nz) {
      u0 = P[2] - node_position0[2];
      v0 = P[0] - node_position0[0];
      u1 = node_position1[2] - node_position0[2];
      u2 = node_position2[2] - node_position0[2];
      v1 = node_position1[0] - node_position0[0];
      v2 = node_position2[0] - node_position0[0];
    }
    else if (nz>=nx && nz>=ny) {
      u0 = P[0] - node_position0[0];
      v0 = P[1] - node_position0[1];
      u1 = node_position1[0] - node_position0[0];
      u2 = node_position2[0] - node_position0[0];
      v1 = node_position1[1] - node_position0[1];
      v2 = node_position2[1] - node_position0[1];
    }
    /*==========================================*/
    /* BEGIN BARYCENTRIC INTERSECTION ALGORITHM */
    /*==========================================*/
    DataType alpha, beta;
    if (u1!=0.0)    {  /* common case */
        beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
        if ((beta>=0.0) && (beta<=1.0)) {
            alpha        = (u0 - beta*u2)/u1;
            intersection = ((alpha>=0.0) && ((alpha+beta)<=1.0));
        }
    } else {           /* uncommon case */
        beta = u0/u2;
        if ((beta>=0.0) && (beta<=1.0)) {
            alpha        = (v0 - beta*v2)/v1;
            intersection = ((alpha>=0.0) && ((alpha+beta)<=1.0));
        }
    }
    if (intersection) {
      coords[0] = P[0];
      coords[1] = P[1];
      coords[2] = P[2];
      break;
    }
  }
/*.....
c
c now improve on the estimate with a tangent constructed at the 
c estimated intersection
  int max_iterations = 200;
  DataType  J[3][3], invJ[3][3];
  DataType* node_position0 = Node(0)->Variable(POSITION);
  DataType* node_position1 = Node(1)->Variable(POSITION);
  DataType* node_position2 = Node(2)->Variable(POSITION);
  DataType* node_position3 = Node(3)->Variable(POSITION);

  J[0][0] = edge_node_position0[0]-edge_node_position1[0];
  J[1][0] = edge_node_position0[1]-edge_node_position1[1];
  J[2][0] = edge_node_position0[2]-edge_node_position1[2];

  J[0][1] = 0.25*( node_position0
  J[1][1] = 
  J[2][1] = 
  
  J[0][2] = 
  J[1][2] = 
  J[2][2] = 
  
  A(1,2) = PFOURTH*(-QUAD(1,1)+QUAD(2,1)+QUAD(3,1)-QUAD(4,1) ) +
    XVQ(3)*PFOURTH*( QUAD(1,1)-QUAD(2,1)+QUAD(3,1)-QUAD(4,1) )

  A(2,2) = PFOURTH*(-QUAD(1,2)+QUAD(2,2)+QUAD(3,2)-QUAD(4,2) ) +
    XVQ(3)*PFOURTH*( QUAD(1,2)-QUAD(2,2)+QUAD(3,2)-QUAD(4,2) )

  A(3,2) = PFOURTH*(-QUAD(1,3)+QUAD(2,3)+QUAD(3,3)-QUAD(4,3) ) +
    XVQ(3)*PFOURTH*( QUAD(1,3)-QUAD(2,3)+QUAD(3,3)-QUAD(4,3) )

  
  A(1,3) = PFOURTH*(-QUAD(1,1)-QUAD(2,1)+QUAD(3,1)+QUAD(4,1) ) +
    XVQ(2)*PFOURTH*( QUAD(1,1)-QUAD(2,1)+QUAD(3,1)-QUAD(4,1) )

  A(2,3) = PFOURTH*(-QUAD(1,2)-QUAD(2,2)+QUAD(3,2)+QUAD(4,2) ) +
    XVQ(2)*PFOURTH*( QUAD(1,2)-QUAD(2,2)+QUAD(3,2)-QUAD(4,2) )

  A(3,3) = PFOURTH*(-QUAD(1,3)-QUAD(2,3)+QUAD(3,3)+QUAD(4,3) ) +
    XVQ(2)*PFOURTH*( QUAD(1,3)-QUAD(2,3)+QUAD(3,3)-QUAD(4,3) )


  
  DataType detJ  =  J[0][0]*J[1][1]*J[2][2]-J[1][1]*J[2][1]-
                J[0][1]*J[1][0]*J[2][2]-J[2][0]*J[1][2]+
                J[0][2]*J[1][0]*J[2][1]-J[2][0]*J[1][1];
  DataType detJI = 1.0/detJ;
                
  invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])*detJI;
  invJ[0][1] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2])*detJI;
  invJ[0][2] = -(J[1][2]*J[0][1]-J[0][2]*J[1][1])*detJI;
  invJ[1][0] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2])*detJI;
  invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])*detJI;
  invJ[1][2] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2])*detJI;
  invJ[2][0] = -(J[1][0]*J[2][1]-J[2][0]*J[1][1])*detJI;
  invJ[2][1] =  (J[0][0]*J[2][1]-J[2][0]*J[0][1])*detJI;
  invJ[2][2] = -(J[0][0]*J[1][1]-J[0][1]*J[1][0])*detJI;

  if (std::fabs(detJ)>1.0e-30) {
    XVQ(1) = 0
    XM  = PFOURTH*( QUAD(1,1) + QUAD(2,1) + QUAD(3,1) + QUAD(4,1))
    YM  = PFOURTH*( QUAD(1,2) + QUAD(2,2) + QUAD(3,2) + QUAD(4,2))
    ZM  = PFOURTH*( QUAD(1,3) + QUAD(2,3) + QUAD(3,3) + QUAD(4,3))
    X2  = PFOURTH*(-QUAD(1,1) + QUAD(2,1) + QUAD(3,1) - QUAD(4,1))
    Y2  = PFOURTH*(-QUAD(1,2) + QUAD(2,2) + QUAD(3,2) - QUAD(4,2))
    Z2  = PFOURTH*(-QUAD(1,3) + QUAD(2,3) + QUAD(3,3) - QUAD(4,3))
    X3  = PFOURTH*(-QUAD(1,1) - QUAD(2,1) + QUAD(3,1) + QUAD(4,1))
    Y3  = PFOURTH*(-QUAD(1,2) - QUAD(2,2) + QUAD(3,2) + QUAD(4,2))
    Z3  = PFOURTH*(-QUAD(1,3) - QUAD(2,3) + QUAD(3,3) + QUAD(4,3))
    X23 = PFOURTH*( QUAD(1,1) - QUAD(2,1) + QUAD(3,1) - QUAD(4,1))
    Y23 = PFOURTH*( QUAD(1,2) - QUAD(2,2) + QUAD(3,2) - QUAD(4,2))
    Z23 = PFOURTH*( QUAD(1,3) - QUAD(2,3) + QUAD(3,3) - QUAD(4,3))
c using modified Newton method
    RESIDM = 1.E20
     20 K = 1,MAXNUMITER

    RHS(1) = VECTOR(1,1) + XVQ(1)*( VECTOR(2,1)-VECTOR(1,1) )  
       - ( XM + X2*XVQ(2) + X3*XVQ(3) + X23*XVQ(2)*XVQ(3) )

    RHS(2) = VECTOR(1,2) + XVQ(1)*( VECTOR(2,2)-VECTOR(1,2) )  
       - ( YM + Y2*XVQ(2) + Y3*XVQ(3) + Y23*XVQ(2)*XVQ(3) )

    RHS(3) = VECTOR(1,3) + XVQ(1)*( VECTOR(2,3)-VECTOR(1,3) )  
       - ( ZM + Z2*XVQ(2) + Z3*XVQ(3) + Z23*XVQ(2)*XVQ(3) )

    RESIDS = RHS(1)**2 + RHS(2)**2 + RHS(3)**2
    IF( RESIDS.LT.1.0D-012 )THEN
      RETURN
    END IF

    XVQ(1) = XVQ(1)+AI(1,1)*RHS(1)+AI(1,2)*RHS(2)+AI(1,3)*RHS(3)
    XVQ(2) = XVQ(2)+AI(2,1)*RHS(1)+AI(2,2)*RHS(2)+AI(2,3)*RHS(3)
    XVQ(3) = XVQ(3)+AI(3,1)*RHS(1)+AI(3,2)*RHS(2)+AI(3,3)*RHS(3)

    IF( RESIDS.LT.RESIDM )THEN
       RESIDM = RESIDS
       XVQ1M = XVQ(1)
       XVQ2M = XVQ(2)
       XVQ3M = XVQ(3)
    END IF

    IF( abs(XVQ(1)).GT.1000. .OR. 
        ABS(XVQ(2)).GT.1000. .OR.
        ABS(XVQ(3)).GT.1000. )THEN
       // vector and plane have no valid intersection
       IERROR = 1
       RETURN
    END IF

 20 CONTINUE
  } else {
    // vector and plane are co-planar
    IERROR = 1
  }

  IF( RESIDS.GT.1.0D-6 )THEN
     print*
     print*,'inititial vec_quad intersection did not converge'
     print*,'   best configuration produced R = ',residm
     print*,'   xvq =',xvq1m,xvq2m,xvq3m
     print*,'   end  configuration produced R = ',resids
     print*,'   xvq =',xvq(1),xvq(2),xvq(3)
  end if

c now improve on the estimate with a tangent constructed at
c every estimated intersection
c
c   restore configuration to best saved state
         XVQ(1) = XVQ1M
         XVQ(2) = XVQ2M
         XVQ(3) = XVQ3M
         RESIDMS = RESIDM
c
         MAXNUMITER = 20
c using Newton method
         A(1,1) = VECTOR(1,1) - VECTOR(2,1)
         A(2,1) = VECTOR(1,2) - VECTOR(2,2)
         A(3,1) = VECTOR(1,3) - VECTOR(2,3)

          XM  = PFOURTH*( QUAD(1,1) + QUAD(2,1) + QUAD(3,1) + QUAD(4,1))
          YM  = PFOURTH*( QUAD(1,2) + QUAD(2,2) + QUAD(3,2) + QUAD(4,2))
          ZM  = PFOURTH*( QUAD(1,3) + QUAD(2,3) + QUAD(3,3) + QUAD(4,3))
          X2  = PFOURTH*(-QUAD(1,1) + QUAD(2,1) + QUAD(3,1) - QUAD(4,1))
          Y2  = PFOURTH*(-QUAD(1,2) + QUAD(2,2) + QUAD(3,2) - QUAD(4,2))
          Z2  = PFOURTH*(-QUAD(1,3) + QUAD(2,3) + QUAD(3,3) - QUAD(4,3))
          X3  = PFOURTH*(-QUAD(1,1) - QUAD(2,1) + QUAD(3,1) + QUAD(4,1))
          Y3  = PFOURTH*(-QUAD(1,2) - QUAD(2,2) + QUAD(3,2) + QUAD(4,2))
          Z3  = PFOURTH*(-QUAD(1,3) - QUAD(2,3) + QUAD(3,3) + QUAD(4,3))
          X23 = PFOURTH*( QUAD(1,1) - QUAD(2,1) + QUAD(3,1) - QUAD(4,1))
          Y23 = PFOURTH*( QUAD(1,2) - QUAD(2,2) + QUAD(3,2) - QUAD(4,2))
          Z23 = PFOURTH*( QUAD(1,3) - QUAD(2,3) + QUAD(3,3) - QUAD(4,3))

         DO 30 K = 1,MAXNUMITER
            
         A(1,2) = PFOURTH*(-QUAD(1,1)+QUAD(2,1)+QUAD(3,1)-QUAD(4,1) ) +
     *     XVQ(3)*PFOURTH*( QUAD(1,1)-QUAD(2,1)+QUAD(3,1)-QUAD(4,1) )

         A(2,2) = PFOURTH*(-QUAD(1,2)+QUAD(2,2)+QUAD(3,2)-QUAD(4,2) ) +
     *     XVQ(3)*PFOURTH*( QUAD(1,2)-QUAD(2,2)+QUAD(3,2)-QUAD(4,2) )

         A(3,2) = PFOURTH*(-QUAD(1,3)+QUAD(2,3)+QUAD(3,3)-QUAD(4,3) ) +
     *     XVQ(3)*PFOURTH*( QUAD(1,3)-QUAD(2,3)+QUAD(3,3)-QUAD(4,3) )

         
         A(1,3) = PFOURTH*(-QUAD(1,1)-QUAD(2,1)+QUAD(3,1)+QUAD(4,1) ) +
     *     XVQ(2)*PFOURTH*( QUAD(1,1)-QUAD(2,1)+QUAD(3,1)-QUAD(4,1) )

         A(2,3) = PFOURTH*(-QUAD(1,2)-QUAD(2,2)+QUAD(3,2)+QUAD(4,2) ) +
     *     XVQ(2)*PFOURTH*( QUAD(1,2)-QUAD(2,2)+QUAD(3,2)-QUAD(4,2) )

         A(3,3) = PFOURTH*(-QUAD(1,3)-QUAD(2,3)+QUAD(3,3)+QUAD(4,3) ) +
     *     XVQ(2)*PFOURTH*( QUAD(1,3)-QUAD(2,3)+QUAD(3,3)-QUAD(4,3) )


         CALL CINV3(A,AI,DET)

         IF( ABS(DET) .GT. 1.0D-30 )THEN
            XVQ(1) = 0

            RHS(1) = VECTOR(1,1) + XVQ(1)*( VECTOR(2,1)-VECTOR(1,1) )  
     *         - ( XM + X2*XVQ(2) + X3*XVQ(3) + X23*XVQ(2)*XVQ(3) )

            RHS(2) = VECTOR(1,2) + XVQ(1)*( VECTOR(2,2)-VECTOR(1,2) )  
     *         - ( YM + Y2*XVQ(2) + Y3*XVQ(3) + Y23*XVQ(2)*XVQ(3) )

            RHS(3) = VECTOR(1,3) + XVQ(1)*( VECTOR(2,3)-VECTOR(1,3) )  
     *         - ( ZM + Z2*XVQ(2) + Z3*XVQ(3) + Z23*XVQ(2)*XVQ(3) )

            RESIDS = RHS(1)**2 + RHS(2)**2 + RHS(3)**2
            IF( RESIDS.LT.1.0D-12 )RETURN

            XVQ(1) = XVQ(1)+AI(1,1)*RHS(1)+AI(1,2)*RHS(2)+AI(1,3)*RHS(3)
            XVQ(2) = XVQ(2)+AI(2,1)*RHS(1)+AI(2,2)*RHS(2)+AI(2,3)*RHS(3)
            XVQ(3) = XVQ(3)+AI(3,1)*RHS(1)+AI(3,2)*RHS(2)+AI(3,3)*RHS(3)

            IF( RESIDS.LT.RESIDM )THEN
               RESIDM = RESIDS
               XVQ1M = XVQ(1)
               XVQ2M = XVQ(2)
               XVQ3M = XVQ(3)
            END IF

            IF( ABS(XVQ(1)).GT.1000. .OR. 
     *          ABS(XVQ(2)).GT.1000. .OR.
     *          ABS(XVQ(3)).GT.1000. )THEN
c vector and plane have no valid intersection
               IERROR = 1
               RETURN
            END IF
         ELSE
c vector and plane are co-planar
            IERROR = 1
         END IF

 30      CONTINUE
C
         IF( RESIDM.LT.RESIDMS .AND. RESIDS.GT.1.0D-6 )THEN
            print*,'vec_quad intersection did not converge'
            print*,'   best configuration produced R = ',residm
            print*,'   xvq =',xvq1m,xvq2m,xvq3m
            print*,'   end  configuration produced R = ',resids
            print*,'   xvq =',xvq(1),xvq(2),xvq(3)
         ELSE IF( RESIDM.GE.RESIDMS .AND. RESIDS.GT.1.0D-6 )THEN
            print*,'Newtons method could not improve intersection'
            print*,'   best configuration produced R = ',residm
            print*,'   xvq =',xvq1m,xvq2m,xvq3m
            print*,'   end  configuration produced R = ',resids
            print*,'   xvq =',xvq(1),xvq(2),xvq(3)
         END IF
         XVQ(1) = XVQ1M
         XVQ(2) = XVQ2M
         XVQ(3) = XVQ3M

      ELSE
         IERROR = 1
      END IF
c
      RETURN
      END
.....*/
  return intersection;
#endif
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Evaluate_Shape_Functions( DataType* local_coords,
						  DataType* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Global_Coordinates( VariableHandle POSITION,
						    DataType* local_coords,
						    DataType* global_coords )
{
  DataType node_positions[4][3];
  for(int i=0; i<4; ++i ){
    DataType* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Local_Coordinates( DataType Config_Param,
					           VariableHandle POSITION0, 
					           VariableHandle POSITION1, 
						   VariableHandle FACE_NORMAL,
						   DataType* global_coords,
						   DataType* local_coords )
{
  int i, j;
  DataType node_positions[4][3];
  if (Config_Param == 0.0) {
    for (i=0; i<Nodes_Per_Face(); ++i) {
      DataType* node_position = Node(i)->Variable(POSITION0);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (i=0; i<Nodes_Per_Face(); ++i) {
      DataType* node_position = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    DataType alpha = 1.0 - Config_Param, beta = Config_Param;
    for (i=0; i<Nodes_Per_Face(); ++i) {
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
void ContactQuadFaceL4<DataType>::Compute_Local_Coordinates( VariableHandle POSITION,
						   DataType* global_coords,
						   DataType* local_coords )
{
  int i, j;
  DataType node_positions[4][3];
  for (i=0; i<Nodes_Per_Face(); ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
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
/* on a linear Q4 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary Q4 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Shape_Functions( DataType* local_coords,
						 DataType* shape_functions )
{
  shape_functions[0] = 0.25*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_functions[1] = 0.25*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_functions[2] = 0.25*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_functions[3] = 0.25*(1.0-local_coords[0])*(1.0+local_coords[1]);
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Shape_Derivatives( DataType* local_coords,
						   DataType shape_derivs[2][4] )
{
  shape_derivs[0][0] = -0.25*(1.0-local_coords[1]);
  shape_derivs[0][1] =  0.25*(1.0-local_coords[1]);
  shape_derivs[0][2] =  0.25*(1.0+local_coords[1]);
  shape_derivs[0][3] = -0.25*(1.0+local_coords[1]);
  
  shape_derivs[1][0] = -0.25*(1.0-local_coords[0]);
  shape_derivs[1][1] = -0.25*(1.0+local_coords[0]);
  shape_derivs[1][2] =  0.25*(1.0+local_coords[0]);
  shape_derivs[1][3] =  0.25*(1.0-local_coords[0]);
}

template<typename DataType>
void 
ContactQuadFaceL4<DataType>::Compute_Local_Coords( DataType node_positions[MAX_NODES_PER_FACE][3], 
					 DataType global_coords[3],
					 DataType local_coords[3] )
{
  Compute_Quad_Local_Coords(node_positions, global_coords, local_coords);
}

template<>
void 
ContactQuadFaceL4<Real>::Compute_Quad_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3], 
					 Real global_coords[3],
					 Real local_coords[3] )
{
  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;
  int  i, j;
  int  nnodes=4;
  if(spatial_tolerance_pre > 0) {
    //
    // check for coincidence with one of the face nodes
    //
    for (i=0; i<nnodes; ++i) {
      Real dx = node_positions[i][0]-global_coords[0];
      Real dy = node_positions[i][1]-global_coords[1];
      Real dz = node_positions[i][2]-global_coords[2];
      Real dd  = dx*dx+dy*dy+dz*dz;
      if (dd == 0 || sqrt(dd) < spatial_tolerance_pre) break;
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
    }
    if (i<nnodes) {
      local_coords[2] = 0.0;
      return;
    }
  }
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=500;
  bool converged = false;
  Real s, s0=0.0, s1, ds=0.0; 
  Real t, t0=0.0, t1, dt=0.0;
  Real coords[3];
  Real J[3][2], f[3], shape_derivatives[2][4];
  Real detJTJ;
  Real residualNorm2, initialResidualNorm2;

  while (!converged && iterations<max_iterations) {
    coords[0] = s0;
    coords[1] = t0;
    Compute_Global_Coords( node_positions, coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives(coords, shape_derivatives);
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
    detJTJ = 1.0/(JTJ[0][0]*JTJ[1][1]-JTJ[0][1]*JTJ[1][0]);
    invJTJ[0][0] =  JTJ[1][1]*detJTJ;
    invJTJ[0][1] = -JTJ[0][1]*detJTJ;
    invJTJ[1][0] = -JTJ[1][0]*detJTJ;
    invJTJ[1][1] =  JTJ[0][0]*detJTJ;

    // APPLY NEWTON ALGORITHM

    Real dx = f[0]-global_coords[0];
    Real dy = f[1]-global_coords[1];
    Real dz = f[2]-global_coords[2];
    Real rx = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
    Real ry = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;
    residualNorm2 = rx*rx + ry*ry;
    if(iterations == 0) initialResidualNorm2 = residualNorm2;
    if(residualNorm2 < newton_tolerance*initialResidualNorm2 || residualNorm2 < abs_newton_tolerance) converged = true;
    else {
      s = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
      t = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;
    
      s1 = s0-(invJTJ[0][0]*s+invJTJ[0][1]*t);
      t1 = t0-(invJTJ[1][0]*s+invJTJ[1][1]*t);
      ds = abs(s1-s0);
      dt = abs(t1-t0);
      s0 = s1;
      t0 = t1;
      //if (ds<newton_tolerance && dt<newton_tolerance) converged = true;
      ++iterations;
    }
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  if (!converged) {
    std::cerr << "ContactQuadFaceL4<Real>::Compute_Local_Coords() did not converge" 
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
    std::cerr << "  After " << iterations << " iterations, local_coords = ("
	      << s0 << "," << t0 << ")" << std::endl;
    std::cerr << "  Going to continuing processing anyway!!!" << std::endl;
  }
#endif
  POSTCONDITION(converged);
  /*if(!converged) {
    std::cerr << "ContactQuadFaceL4<Real>::Compute_Local_Coords() did not converge, residualNorm2 = " << residualNorm2 
              << ", initialResidualNorm2 = " << initialResidualNorm2 << std::endl;
  }*/
  if(spatial_tolerance_post > 0) {
    // If it's close to any of the edges, snap to it
    if (abs(s0)<1.0+spatial_tolerance_post) {
      s0 = min(s0, 1.0);
      s0 = max(s0,-1.0);
    }
    if (abs(t0)<1.0+spatial_tolerance_post) {
      t0 = min(t0, 1.0);
      t0 = max(t0,-1.0);
    }
  }
  local_coords[0] = s0;
  local_coords[1] = t0;
  local_coords[2] = 0.0;
}

template<typename ActiveScalar>
void 
ContactQuadFaceL4<ActiveScalar>::Compute_Quad_Local_Coords( ActiveScalar active_node_positions[MAX_NODES_PER_FACE][3], 
			 		 ActiveScalar active_global_coords[3],
					 ActiveScalar active_local_coords[3] )
{
  int  i, j;
  const int  nnodes=4;

  double node_positions[nnodes][3], global_coords[3], local_coords[3];
  for(i=0; i<3; ++i) {
    global_coords[i] = GetActiveScalarValue(active_global_coords[i]);
    for(j=0; j<nnodes; ++j)
      node_positions[j][i] = GetActiveScalarValue(active_node_positions[j][i]);
  }

  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;

  if(spatial_tolerance_pre > 0) {
    //
    // check for coincidence with one of the face nodes
    //
    for (i=0; i<nnodes; ++i) {
      double dx = node_positions[i][0]-global_coords[0];
      double dy = node_positions[i][1]-global_coords[1];
      double dz = node_positions[i][2]-global_coords[2];
      double dd = dx*dx+dy*dy+dz*dz;
      if (dd == 0 || sqrt(dd) < spatial_tolerance_pre) break;
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
    }
    if (i<nnodes) {
      active_local_coords[0] = local_coords[0];
      active_local_coords[1] = local_coords[1];
      active_local_coords[2] = 0.0;
      return;
    }
  }
  //
  // else use newton's method to iterate (values only)
  //
  int  iterations=0;
  int  max_iterations=500;
  bool converged = false;
  double s, s0=0.0, s1, ds=0.0; 
  double t, t0=0.0, t1, dt=0.0;
  double J[3][2], f[3], shape_derivatives[2][4], shape_functions[4];
  double JT[2][3], JTJ[2][2], invJTJ[2][2];
  double detJTJ;
  double residualNorm2, initialResidualNorm2;
  double s0_copy=0.0, t0_copy=0.0;
  
  while (!converged && iterations<max_iterations) {
    local_coords[0] = s0;
    local_coords[1] = t0;

    // BUILD JACOBIAN AND INVERT
    shape_derivatives[0][0] = -0.25*(1.0-local_coords[1]);
    shape_derivatives[0][1] =  0.25*(1.0-local_coords[1]);
    shape_derivatives[0][2] =  0.25*(1.0+local_coords[1]);
    shape_derivatives[0][3] = -0.25*(1.0+local_coords[1]);
  
    shape_derivatives[1][0] = -0.25*(1.0-local_coords[0]);
    shape_derivatives[1][1] = -0.25*(1.0+local_coords[0]);
    shape_derivatives[1][2] =  0.25*(1.0+local_coords[0]);
    shape_derivatives[1][3] =  0.25*(1.0-local_coords[0]);

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
    JT[0][0] = J[0][0];
    JT[0][1] = J[1][0];
    JT[0][2] = J[2][0];
    JT[1][0] = J[0][1];
    JT[1][1] = J[1][1];
    JT[1][2] = J[2][1];
    
    JTJ[0][0] = JT[0][0]*J[0][0] + JT[0][1]*J[1][0] + JT[0][2]*J[2][0];
    JTJ[0][1] = JT[0][0]*J[0][1] + JT[0][1]*J[1][1] + JT[0][2]*J[2][1];
    JTJ[1][0] = JT[1][0]*J[0][0] + JT[1][1]*J[1][0] + JT[1][2]*J[2][0];
    JTJ[1][1] = JT[1][0]*J[0][1] + JT[1][1]*J[1][1] + JT[1][2]*J[2][1];
    
    detJTJ  = 1.0/(JTJ[0][0]*JTJ[1][1]-JTJ[0][1]*JTJ[1][0]);
    invJTJ[0][0] =  JTJ[1][1]*detJTJ;
    invJTJ[0][1] = -JTJ[0][1]*detJTJ;
    invJTJ[1][0] = -JTJ[1][0]*detJTJ;
    invJTJ[1][1] =  JTJ[0][0]*detJTJ;

    // APPLY NEWTON ALGORITHM
    shape_functions[0] = 0.25*(1.0-local_coords[0])*(1.0-local_coords[1]);
    shape_functions[1] = 0.25*(1.0+local_coords[0])*(1.0-local_coords[1]);
    shape_functions[2] = 0.25*(1.0+local_coords[0])*(1.0+local_coords[1]);
    shape_functions[3] = 0.25*(1.0-local_coords[0])*(1.0+local_coords[1]);
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    for( int i=0 ; i<nnodes ; ++i ){
      for (int j=0; j<3; ++j) {
        f[j] += shape_functions[i]*node_positions[i][j];
      }
    }

    double dx = f[0]-global_coords[0];
    double dy = f[1]-global_coords[1];
    double dz = f[2]-global_coords[2];
    double rx = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
    double ry = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;
    residualNorm2 = rx*rx + ry*ry;
    if(iterations == 0) initialResidualNorm2 = residualNorm2;
    if(residualNorm2 < newton_tolerance*initialResidualNorm2 || residualNorm2 < abs_newton_tolerance) converged = true;
    else {
      s = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
      t = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;

      s1 = s0-(invJTJ[0][0]*s+invJTJ[0][1]*t);
      t1 = t0-(invJTJ[1][0]*s+invJTJ[1][1]*t);
      ds = abs(s1-s0);
      dt = abs(t1-t0);
      s0_copy=s0; t0_copy=t0;
      s0 = s1;
      t0 = t1;
      //if (ds<newton_tolerance && dt<newton_tolerance) converged = true;
      ++iterations;
    }
  }
  /*if(!converged) {
    std::cerr << "ContactQuadFaceL4<ActiveScalar>::Compute_Local_Coords() did not converge, residualNorm2 = " << residualNorm2 
              << ", initialResidualNorm2 = " << initialResidualNorm2 << std::endl;
  }*/
  {
    //
    // repeat last newton iteration to get derivatives
    //
    int  iterations=0;
    int  max_iterations=1;
    bool converged = false;
    ActiveScalar s, s0=s0_copy, s1, ds=0.0; 
    ActiveScalar t, t0=t0_copy, t1, dt=0.0;
    ActiveScalar J[3][2], f[3], shape_derivatives[2][4];
    ActiveScalar JT[2][3], JTJ[2][2], invJTJ[2][2];
  
    while (!converged && iterations<max_iterations) {
      active_local_coords[0] = s0;
      active_local_coords[1] = t0;

      // BUILD JACOBIAN AND INVERT
      Compute_Shape_Derivatives(active_local_coords, shape_derivatives);
      for (i=0; i<2; ++i) {
        J[0][i] = 0.0;
        J[1][i] = 0.0;
        J[2][i] = 0.0;
        for (j=0; j<nnodes; ++j) {
          J[0][i] += shape_derivatives[i][j]*active_node_positions[j][0];
          J[1][i] += shape_derivatives[i][j]*active_node_positions[j][1];
          J[2][i] += shape_derivatives[i][j]*active_node_positions[j][2];
        }
      }
      JT[0][0] = J[0][0];
      JT[0][1] = J[1][0];
      JT[0][2] = J[2][0];
      JT[1][0] = J[0][1];
      JT[1][1] = J[1][1];
      JT[1][2] = J[2][1];
    
      JTJ[0][0] = JT[0][0]*J[0][0] + JT[0][1]*J[1][0] + JT[0][2]*J[2][0];
      JTJ[0][1] = JT[0][0]*J[0][1] + JT[0][1]*J[1][1] + JT[0][2]*J[2][1];
      JTJ[1][0] = JT[1][0]*J[0][0] + JT[1][1]*J[1][0] + JT[1][2]*J[2][0];
      JTJ[1][1] = JT[1][0]*J[0][1] + JT[1][1]*J[1][1] + JT[1][2]*J[2][1];
    
      ActiveScalar detJTJ  = 1.0/(JTJ[0][0]*JTJ[1][1]-JTJ[0][1]*JTJ[1][0]);
      invJTJ[0][0] =  JTJ[1][1]*detJTJ;
      invJTJ[0][1] = -JTJ[0][1]*detJTJ;
      invJTJ[1][0] = -JTJ[1][0]*detJTJ;
      invJTJ[1][1] =  JTJ[0][0]*detJTJ;

      // APPLY NEWTON ALGORITHM
      Compute_Global_Coords( active_node_positions, active_local_coords, f );
      ActiveScalar dx = f[0]-active_global_coords[0];
      ActiveScalar dy = f[1]-active_global_coords[1];
      ActiveScalar dz = f[2]-active_global_coords[2];
      s = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
      t = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;
    
      s1 = s0-(invJTJ[0][0]*s+invJTJ[0][1]*t);
      t1 = t0-(invJTJ[1][0]*s+invJTJ[1][1]*t);
      ds = abs(s1-s0);
      dt = abs(t1-t0);
      s0 = s1;
      t0 = t1;
      //if (ds<newton_tolerance && dt<newton_tolerance) converged = true;
      ++iterations;
    }
    if(spatial_tolerance_post > 0) {
      // If it's close to any of the edges, snap to it
      if (abs(s0)<1.0+spatial_tolerance_post) {
        s0 = min(s0, 1.0);
        s0 = max(s0,-1.0);
      }
      if (abs(t0)<1.0+spatial_tolerance_post) {
        t0 = min(t0, 1.0);
        t0 = max(t0,-1.0);
      }
    }
    active_local_coords[0] = s0;
    active_local_coords[1] = t0;
    active_local_coords[2] = 0.0;
  }
}

template<typename DataType>
void ContactQuadFaceL4<DataType>::Compute_Global_Coords( DataType node_positions[4][3],
					       DataType local_coords[2],
					       DataType global_coords[3] )
{
  DataType N[4];
  int  nnodes=4;
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
void  ContactQuadFaceL4<DataType>::Interpolate_Scalar( DataType  local_coords[2],
					     DataType  node_scalars[4],
					     DataType& interpolated_scalar )
{
  DataType N[4];
  int  nnodes=4;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

template<typename DataType>
void  ContactQuadFaceL4<DataType>::Interpolate_Vector( DataType local_coords[2],
					     DataType node_vectors[4][3],
					     DataType interpolated_vector[3] )
{
  DataType N[4];
  int  nnodes=4;
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

#endif  // #define ContactQuadFaceL4_C_
