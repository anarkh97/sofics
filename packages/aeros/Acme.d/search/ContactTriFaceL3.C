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


#ifndef ContactTriFaceL3_C_
#define ContactTriFaceL3_C_

#include <algorithm>

#include "allocators.h"
#include "ContactUtilities.h"
#include "ContactTriFaceL3.h"
#include "ContactQuadFaceL4.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include "contact_tolerances.h"
#include <iostream>
#include <cmath>
#include <new>

using acme::Dot;
using acme::Normalize;
using acme::Magnitude;
using acme::ScalarTripleProduct;

template<typename DataType>
ContactTriFaceL3<DataType>::ContactTriFaceL3( ContactFixedSizeAllocator* alloc,
                                    int Block_Index, 
				    int Index_in_Block,int key ) 
  : ContactFace<DataType>( alloc, ContactSearch::TRIFACEL3,
                 Block_Index, Index_in_Block, key, 
                 nodes, edges, Node_Info, Edge_Info) 
{}

template<typename DataType>
ContactTriFaceL3<DataType>* ContactTriFaceL3<DataType>::new_ContactTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactTriFaceL3].New_Frag())
             ContactTriFaceL3<DataType>(alloc, Block_Index, Index_in_Block, key);
}

template<typename DataType>
void ContactTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactTriFaceL3<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactTriFaceL3<DataType> allocator" );
}

template<typename DataType>
ContactTriFaceL3<DataType>::~ContactTriFaceL3() {}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL )
{
  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);

  DataType* Position0 = node0->Variable(CURRENT_POSITION);
  DataType* Position1 = node1->Variable(CURRENT_POSITION);
  DataType* Position2 = node2->Variable(CURRENT_POSITION);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  DataType Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  DataType* face_normal = Variable(FACE_NORMAL);

  acme::Cross(Vec01, Vec02, face_normal);

  Normalize(face_normal);

}

#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>
template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Partial_Face_Normal(VariableHandle CURRENT_POSITION, DataType (*dface_normal)[3] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(CURRENT_POSITION)[0],
                                             Node(i)->Variable(CURRENT_POSITION)[1],
                                             Node(i)->Variable(CURRENT_POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes); 
  face.GetdUnitNormal<DataType,CoordSetType>(dface_normal, (DataType*)NULL, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Second_Partial_Face_Normal(VariableHandle CURRENT_POSITION, DataType (*d2face_normal)[3] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(CURRENT_POSITION)[0],
                                             Node(i)->Variable(CURRENT_POSITION)[1],
                                             Node(i)->Variable(CURRENT_POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes);
  face.Getd2UnitNormal<DataType,CoordSetType>(d2face_normal, (DataType*)NULL, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Normal(VariableHandle POSITION,
				       DataType* normal, DataType* local_coords )
{
  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);

  DataType* Position0    = node0->Variable(POSITION);
  DataType* Position1    = node1->Variable(POSITION);
  DataType* Position2    = node2->Variable(POSITION);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  DataType Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  acme::Cross(Vec01,Vec02,normal);
  Normalize(normal);
}

//
// KHP:Q2: The normal vector variable is the third argument in this
// KHP:Q2: method while it is the second argument in the previous Compute_Normal
// KHP:Q2: method. It would be nice to have the arguments be in a consistent order.
//
template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Normal(DataType** nodal_positions,
				      DataType* local_coords, 
				      DataType* normal )
{
  DataType* Position0    = nodal_positions[0];
  DataType* Position1    = nodal_positions[1];
  DataType* Position2    = nodal_positions[2];

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  DataType Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  acme::Cross(Vec01,Vec02,normal);
  Normalize(normal);
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_CharacteristicLength( VariableHandle CURRENT_POSITION,
						     VariableHandle CHARACTERISTIC_LENGTH )
{
  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);

  DataType* Position0 = node0->Variable(CURRENT_POSITION);
  DataType* Position1 = node1->Variable(CURRENT_POSITION);
  DataType* Position2 = node2->Variable(CURRENT_POSITION);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  DataType Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  DataType face_normal[3];
  acme::Cross(Vec01,Vec02,face_normal);

  DataType* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  //
  // KHP:Q3: Couldn't we compute the square of the characteristic length
  // KHP:Q3: and avoid the sqrt operation?
  //
  *characteristiclength = Magnitude(face_normal) * 0.5;
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Centroid( VariableHandle CURRENT_POSITION,
					 VariableHandle CENTROID )
{
  DataType* centroid = Variable(CENTROID);
  DataType* node_position = Node(0)->Variable(CURRENT_POSITION);
  centroid[0] = node_position[0];
  centroid[1] = node_position[1];
  centroid[2] = node_position[2];
  for( int i=1 ; i<3 ; ++i ){
    node_position = Node(i)->Variable(CURRENT_POSITION);
    centroid[0] += node_position[0];
    centroid[1] += node_position[1];
    centroid[2] += node_position[2];
  }
  DataType scale   = 1.0/3.0;
  centroid[0] *= scale;
  centroid[1] *= scale;
  centroid[2] *= scale;
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Get_Edge_Nodes( int i, ContactNode<DataType>** node )
{
  PRECONDITION( i>=0 && i<3 );
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
    node[1] = Node(0);
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Invalid Edge request in ContactTriFaceL3<DataType>::Get_Face_Nodes" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

template<typename DataType>
int ContactTriFaceL3<DataType>::Get_Edge_Number( ContactNode<DataType>** edge_nodes )
{
  PRECONDITION( edge_nodes[0] && edge_nodes[1] );
  PRECONDITION( edge_nodes[0] != edge_nodes[1] );
  int i1=-1,i2=-1;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    ContactNode<DataType>* node = Node(i);
    if( edge_nodes[0] == node ) i1=i;
    if( edge_nodes[1] == node ) i2=i;
  }
  PRECONDITION( 0<=i1 && i1<3 );
  PRECONDITION( 0<=i2 && i2<3 );
  switch( i1 ){
  case 0:
    if( i2 == 1 ) return(0);
    if( i2 == 2 ) return(2);
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
    if( i2 == 0 ) return(2);
    if( i2 == 1 ) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  }
  POSTCONDITION( 0 );
  return(-1);
}

template<typename DataType>
int ContactTriFaceL3<DataType>::Get_Edge_Number( DataType* local_coords )
{
  using std::abs;
  DataType tol(1.0e-8);
  if (abs(local_coords[1])<tol && abs(local_coords[3])<tol) return(0);
  if (abs(local_coords[0])<tol && abs(local_coords[2])<tol) return(1);
  if (abs((1.0-local_coords[0]-local_coords[1]))<tol &&
      abs((1.0-local_coords[2]-local_coords[3]))<tol) return(2);
  return( -1 );
}

template<typename DataType>
void
ContactTriFaceL3<DataType>::Compute_Edge_Normal( VariableHandle POSITION,
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
  DataType* face_normal = Variable(FACE_NORMAL);
  acme::Cross(edge_tangent, face_normal, edge_normal);
  Normalize(edge_normal);
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Get_Close_Edges( DataType* local_coords, int& number,
					int& edge_1_id, int& edge_2_id ){
  DataType a0 = local_coords[0];
  DataType a1 = local_coords[1];
  DataType a2 = local_coords[2];

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

template<typename DataType>
bool ContactTriFaceL3<DataType>::Is_Inside_Face( DataType* local_coords )
{
  if( local_coords[0] >= 0.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= 0.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= 0.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

template<typename DataType>
ContactFace<DataType>* ContactTriFaceL3<DataType>::Neighbor( DataType* local_coords )
{
  PRECONDITION( 0 );
  return (ContactFace<DataType>*) NULL;
}

template<typename DataType>
void ContactTriFaceL3<DataType>::FacetDecomposition(int& nfacets,
          DataType* coordinates0, DataType* normals0, VariableHandle POSITION0,
          DataType* coordinates1, DataType* normals1, VariableHandle POSITION1,
          DataType* coordinates2, DataType* normals2, VariableHandle POSITION2) 
{
  nfacets = 1;
  for(int j=0 ; j<3 ; ++j ){
	DataType* node_position  = Node(j)->Variable(POSITION0);
	coordinates0[0+3*j]  = node_position[0];
	coordinates0[1+3*j]  = node_position[1];
	coordinates0[2+3*j]  = node_position[2];
  }

  ContactNode<DataType>* node0 = Node(0);
  ContactNode<DataType>* node1 = Node(1);
  ContactNode<DataType>* node2 = Node(2);

  DataType* Position0 = node0->Variable(POSITION0);
  DataType* Position1 = node1->Variable(POSITION0);
  DataType* Position2 = node2->Variable(POSITION0);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  DataType Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  acme::Cross(Vec01,Vec02,normals0);
  Normalize(normals0);
 
  if (coordinates1!=NULL) {
    for(int j=0 ; j<3 ; ++j ){
      DataType* node_position = Node(j)->Variable(POSITION1);
      coordinates1[0+3*j]  = node_position[0];
      coordinates1[1+3*j]  = node_position[1];
      coordinates1[2+3*j]  = node_position[2];
    }

    node0     = Node(0);
    node1     = Node(1);
    node2     = Node(2);

    Position0 = node0->Variable(POSITION1);
    Position1 = node1->Variable(POSITION1);
    Position2 = node2->Variable(POSITION1);

    // Compute the vector from node 0 to node 1 & from node 0 to node 2
    Vec01[0] = Position1[0] - Position0[0];
    Vec01[1] = Position1[1] - Position0[1];
    Vec01[2] = Position1[2] - Position0[2];
    Vec02[0] = Position2[0] - Position0[0];
    Vec02[1] = Position2[1] - Position0[1];
    Vec02[2] = Position2[2] - Position0[2];

    // Compute the face normal as the cross product of the two vectors
    acme::Cross(Vec01,Vec02,normals1);
    Normalize(normals1);
  }

  if (coordinates2!=NULL) {
    for(int j=0 ; j<3 ; ++j ){
      DataType* node_position = Node(j)->Variable(POSITION2);
      coordinates2[0+3*j]  = node_position[0];
      coordinates2[1+3*j]  = node_position[1];
      coordinates2[2+3*j]  = node_position[2];
    }

    node0     = Node(0);
    node1     = Node(1);
    node2     = Node(2);

    Position0 = node0->Variable(POSITION2);
    Position1 = node1->Variable(POSITION2);
    Position2 = node2->Variable(POSITION2);

    // Compute the vector from node 0 to node 1 & from node 0 to node 2
    Vec01[0] = Position1[0] - Position0[0];
    Vec01[1] = Position1[1] - Position0[1];
    Vec01[2] = Position1[2] - Position0[2];
    Vec02[0] = Position2[0] - Position0[0];
    Vec02[1] = Position2[1] - Position0[1];
    Vec02[2] = Position2[2] - Position0[2];

    // Compute the face normal as the cross product of the two vectors
    acme::Cross(Vec01,Vec02,normals2);
    Normalize(normals2);
  }
}

template<typename DataType>
void ContactTriFaceL3<DataType>::FacetStaticRestriction(int nfacets, DataType* coordinates, 
					      DataType* normals, DataType* ctrcl_facets,
					      DataType* ctrcl)
{
  for (int i=0; i<LENGTH; ++i) {
	ctrcl[i] = ctrcl_facets[i];
  }
}
template<typename DataType>
void ContactTriFaceL3<DataType>::FacetDynamicRestriction(int nfacets, 
                                               DataType* ctrcl_facets, 
                                               DataType* ctrcl)
{
  for (int i=0; i<LENGTH; ++i) {
	ctrcl[i] = ctrcl_facets[i];
  }
}



template<typename DataType>
void ContactTriFaceL3<DataType>::Smooth_Normal( VariableHandle CURRENT_POSITION,
				      VariableHandle NODE_NORMAL,
				      VariableHandle FACE_NORMAL,
				      VariableHandle CURVATURE,
    		       ContactSearch::Smoothing_Resolution resolution,
				      DataType percentage,
				      DataType* coordinates,
				      DataType* smooth_normal,
                                      DataType critical_curvature )
{
  /*                   2
                     h/ \g
                     /\G/\
                    /  c  \
                   /  / \  \
                  /  /   \  \
                 /  /     \  \
                / D/   A   \C \
               /  /         \  \
              /  /           \  \
             i--a-------------b--f
            / E/       B       \F \
           0--d-----------------e--1
              
      Region    A: Use surface normal
              B-D: Smooth along one edge
              E-G: Smooth along two edges
  */
  using std::abs;

  DataType upper_bound = 1.0 - percentage;
  DataType lower_bound = percentage/2.0;
  DataType* face_normal = Variable(FACE_NORMAL);
  DataType node_position[4][3];
  DataType node_normal[4][3];
  DataType contact_point[3];
  DataType coords[3];
  bool smooth_edge0,smooth_edge1;
  DataType curvature0,curvature1;

  if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
      coordinates[1] >= lower_bound && coordinates[1] <= upper_bound &&
      coordinates[2] >= lower_bound && coordinates[2] <= upper_bound ){
    // Region A
    smooth_normal[0] = face_normal[0];
    smooth_normal[1] = face_normal[1];
    smooth_normal[2] = face_normal[2];
    return;
  } 
  else if( coordinates[0] >= lower_bound && coordinates[1] >= lower_bound &&
	   coordinates[2] <  lower_bound){
    /* Region B
                  1-------------0
                 /       B       \
                2-----------------3
    */     
    if( abs(GetEdgeCurvature(0)) > critical_curvature) {
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = 1.0-lower_bound;
    coords[1] = lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = lower_bound;
    coords[1] = 1.0-lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
  }
  else if( coordinates[0] <  lower_bound && coordinates[1] >= lower_bound &&
	   coordinates[2] >= lower_bound ){
    /* Region C
                        3            
                       /\
                      0  \
                       \  \
                        \  \
                         \  \
                          \C \
                           \  \
                            \  \
                             1--2
    */
    if( abs(GetEdgeCurvature(1)) > critical_curvature) {
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = 0.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );

  }
  else if( coordinates[0] >= lower_bound && coordinates[1] <  lower_bound &&
	   coordinates[2] >= lower_bound ){
    /* Region D
                       2
                       /\
                      /  1
                     /  /
                    /  /
                   /  /
                  / D/
                 /  /
                /  /
               3--0
    */
    if( abs(GetEdgeCurvature(2)) > critical_curvature) {
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = lower_bound;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    this->GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = 1.0-lower_bound;
    coords[1] = 0.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    this->GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );

  }
  else if( coordinates[0] > upper_bound && coordinates[1] < lower_bound &&
	   coordinates[2] < lower_bound ){
    /* Region E
      
               1--0
              / E/
             2--3
    */
    curvature0 = GetEdgeCurvature(2);
    curvature1 = GetEdgeCurvature(0);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( abs(curvature1) > critical_curvature) smooth_edge1 = false;
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = 1.0-lower_bound;
    coords[1] = 0.0;
    coords[2] = lower_bound;
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
    // Compute the data for node 2
    node_position[2][0] = Node(0)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(0)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(0)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(0)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(0)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(0)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } 
    // Compute the data for node 3
    coords[0] = 1.0-lower_bound;
    coords[1] = lower_bound;
    coords[2] = 0.0;
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
  else if( coordinates[0] < lower_bound && coordinates[1] > upper_bound &&
	   coordinates[2] < lower_bound ){
    /* Region F
       
             0--3
              \F \
               1--2
    */
    curvature0 = GetEdgeCurvature(0);
    curvature1 = GetEdgeCurvature(1);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( abs(curvature1) > critical_curvature) smooth_edge1 = false;
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
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = 1.0-lower_bound;
    coords[2] = 0.0;
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
    // Compute the data for node 2
    node_position[2][0] = Node(1)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(1)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(1)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(1)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(1)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(1)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute the data for node 3
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound;
    coords[2] = lower_bound;
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
  else if( coordinates[0] < lower_bound && coordinates[1] < lower_bound &&
	   coordinates[2] > upper_bound ){
    /* Region G
       
            2
          3/ \1
           \G/
            0
    */
    curvature0 = GetEdgeCurvature(1);
    curvature1 = GetEdgeCurvature(2);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( abs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( abs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Get the global coordinates of the contact point
    Compute_Global_Coordinates( CURRENT_POSITION, coordinates, contact_point );
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = 0.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-lower_bound;
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
    // Compute the data for node 2
    node_position[2][0] = Node(2)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(2)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(2)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(2)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(2)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(2)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      this->GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      this->GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute the data for node 3
    coords[0] = lower_bound;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound;
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
  // Get the global coordinates of the contact point
  Compute_Global_Coordinates( CURRENT_POSITION, coordinates, contact_point );

  // Compute Contact Point in sub area
  DataType local_coords[3];
  ContactQuadFaceL4<DataType>::Compute_Quad_Local_Coords( node_position,
					       contact_point,
					       local_coords);

  // Now interpolate to get the normal
  ContactQuadFaceL4<DataType>::Interpolate_Vector( local_coords, 
					 node_normal, 
					 smooth_normal );

  // Now make smooth_normal a unit vector
  Normalize(smooth_normal);
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Node_Areas( VariableHandle POSITION,
	                                   VariableHandle /*FACE_NORMAL*/,
                                           DataType* node_areas )
{
  // Get the nodal positions
  DataType* pos_0 = Node(0)->Variable(POSITION);
  DataType* pos_1 = Node(1)->Variable(POSITION);
  DataType* pos_2 = Node(2)->Variable(POSITION);

  // Construct vectors between the nodes
  DataType vec_0_1[3], vec_0_2[3];
  vec_0_1[0] = pos_1[0] - pos_0[0];
  vec_0_1[1] = pos_1[1] - pos_0[1];
  vec_0_1[2] = pos_1[2] - pos_0[2];
  vec_0_2[0] = pos_2[0] - pos_0[0];
  vec_0_2[1] = pos_2[1] - pos_0[1];
  vec_0_2[2] = pos_2[2] - pos_0[2];

  // Take the cross product
  DataType cross[3];
  acme::Cross(vec_0_1, vec_0_2, cross);

  // The total area is then the 1/2 the magnitude of the vector.
  // The node area is then 1/3 of the total area.
  DataType node_area = Magnitude(cross) / 6.0;

  node_areas[0] = node_area;
  node_areas[1] = node_area;
  node_areas[2] = node_area;
}


template<typename DataType>
int ContactTriFaceL3<DataType>::FaceEdge_Intersection(VariableHandle POSITION,
                                            ContactEdge<DataType>* edge,DataType* coords)
{
  using std::abs;
  int intersection=0;
#if 0
//#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactTriFaceL3<DataType>::FaceEdge_Intersection not yet implemented\n";
//#endif
  POSTCONDITION( 0 );
  return intersection;
#else
  // Do bounding box check
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
       edge_max[2]<face_min[2] ) return 0;
       
  // Compute edge ray
  DataType  edge_pnt[3];
  DataType  edge_dir[3];
  DataType* edge_node_position0 = edge->Node(0)->Variable(POSITION);
  DataType* edge_node_position1 = edge->Node(1)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_pnt[j] = edge_node_position0[j];
    edge_dir[j] = edge_node_position1[j]-edge_node_position0[j];
  }
  Normalize(edge_dir);
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
  
  DataType* node_position0 = Node(0)->Variable(POSITION);
  DataType* node_position1 = Node(1)->Variable(POSITION);
  DataType* node_position2 = Node(2)->Variable(POSITION);

  // Compute the face normal
  DataType dx1  = node_position1[0] - node_position0[0];
  DataType dy1  = node_position1[1] - node_position0[1];
  DataType dz1  = node_position1[2] - node_position0[2];
  DataType dx2  = node_position2[0] - node_position0[0];
  DataType dy2  = node_position2[1] - node_position0[1];
  DataType dz2  = node_position2[2] - node_position0[2];

  DataType normal[3];
  acme::Cross(dx1,dy1,dz1, // vector a
              dx2,dy2,dz2, // vector b
	      normal);     // vector c = a X b;

  Normalize(normal);
  
  // determine (edge ray)/(face plane) intersection
  DataType n_dot_d = Dot(normal, edge_dir);

  // if (edge ray) and (tri3 plane) are parallel => no intersection
  if (abs(n_dot_d)<is_parallel_tol) {
    return 0;
  }

  DataType q[3] = { node_position0[0] - edge_pnt[0],
                    node_position0[1] - edge_pnt[1],
                    node_position0[2] - edge_pnt[2]};

  DataType t = Dot(normal,q)/n_dot_d;

  if (t<0.0 || t>tmax) {
    return 0;
  }

  DataType P[3] = { edge_pnt[0] + edge_dir[0]*t,
                    edge_pnt[1] + edge_dir[1]*t,
                    edge_pnt[2] + edge_dir[2]*t};

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
  }
  return (intersection);
#endif
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Evaluate_Shape_Functions( DataType* local_coords,
						 DataType* shape_functions )
{
  using std::fabs;
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Compute_Shape_Functions( local_coords, shape_functions );
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Global_Coordinates( VariableHandle POSITION,
						   DataType* local_coords,
						   DataType* global_coords )
{
  DataType node_positions[3][3];
  for(int i=0; i<3; ++i ){
    DataType* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Local_Coordinates( DataType Config_Param,
						  VariableHandle POSITION0,
						  VariableHandle POSITION1,
						  VariableHandle FACE_NORMAL,
						  DataType* global_coords,
						  DataType* local_coords )
{
  int  i, j;
  DataType node_positions[3][3];
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
void ContactTriFaceL3<DataType>::Compute_Local_Coordinates( VariableHandle POSITION,
						  DataType* global_coords,
						  DataType* local_coords )
{
  int  i, j;
  DataType node_positions[3][3];
  for (i=0; i<Nodes_Per_Face(); ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Partial_Local_Coordinates_1( VariableHandle POSITION,
                                                                      DataType* global_coords,
                                                                      DataType (*dlocal_coords)[2] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(POSITION)[0],
                                             Node(i)->Variable(POSITION)[1],
                                             Node(i)->Variable(POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes);
  face.GetdLocalCoords<DataType,CoordSetType>(dlocal_coords, global_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Partial_Local_Coordinates_2( VariableHandle POSITION,
                                                                      DataType* global_coords,
                                                                      DataType dmdX[2], DataType dmdY[2], DataType dmdZ[2] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(POSITION)[0],
                                             Node(i)->Variable(POSITION)[1],
                                             Node(i)->Variable(POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes);
  face.ComputedmdXdmdYAnddmdZ<DataType,CoordSetType>(dmdX, dmdY, dmdZ, global_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Second_Partial_Local_Coordinates_1( VariableHandle POSITION,
                                                                             DataType* global_coords,
                                                                             DataType (*d2local_coords)[2] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(POSITION)[0],
                                             Node(i)->Variable(POSITION)[1],
                                             Node(i)->Variable(POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes);
  face.Getd2LocalCoords<DataType,CoordSetType>(d2local_coords, global_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Second_Partial_Local_Coordinates_2( VariableHandle POSITION,
                                                                             DataType* global_coords,
                                                                             DataType d2mdX2[2], DataType d2mdY2[2], DataType d2mdZ2[2],
                                                                             DataType d2mdXdY[2], DataType d2mdYdZ[2], DataType d2mdXdZ[2] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(POSITION)[0],
                                             Node(i)->Variable(POSITION)[1],
                                             Node(i)->Variable(POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes);
  face.Computed2mdX2d2mdY2Etc<DataType,CoordSetType>(d2mdX2, d2mdY2, d2mdZ2, d2mdXdY, d2mdYdZ, d2mdXdZ, global_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Second_Partial_Local_Coordinates_12( VariableHandle POSITION,
                                                                              DataType* global_coords,
                                                                              DataType (*ddLocalCoordsdX)[2],
                                                                              DataType (*ddLocalCoordsdY)[2],
                                                                              DataType (*ddLocalCoordsdZ)[2] )
{
  typedef std::vector<NodeTemplate<DataType>*> CoordSetType;
  CoordSetType CoordSet(this->Nodes_Per_Face());
  for (int i=0; i<this->Nodes_Per_Face(); ++i) {
    CoordSet[i] = new NodeTemplate<DataType>(Node(i)->Variable(POSITION)[0],
                                             Node(i)->Variable(POSITION)[1],
                                             Node(i)->Variable(POSITION)[2]);
  }

  int nodes[3] = { 0, 1, 2 };
  FaceTri3 face(nodes);
  face.GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ<DataType,CoordSetType>(ddLocalCoordsdX, ddLocalCoordsdY, ddLocalCoordsdZ, global_coords, CoordSet);

  for (int i=0; i<this->Nodes_Per_Face(); ++i) delete CoordSet[i];
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

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Shape_Functions( DataType local_coords[3],
						DataType shape_functions[3] )
{
  using std::fabs; 
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  shape_functions[0] = local_coords[0];
  shape_functions[1] = local_coords[1];
  shape_functions[2] = local_coords[2];
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Shape_Derivatives( DataType local_coords[3],
					          DataType shape_derivs[2][3] )
{
  using std::fabs;
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  shape_derivs[0][0] =  1.0;
  shape_derivs[0][1] =  0.0;
  shape_derivs[0][2] = -1.0;
  shape_derivs[1][0] =  0.0;
  shape_derivs[1][1] =  1.0;
  shape_derivs[1][2] = -1.0;
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Local_Coords( DataType node_positions[MAX_NODES_PER_FACE][3],
					     DataType global_coords[3],
					     DataType local_coords[3] )
{
  using std::sqrt;
  int  i;
  //
  // check for coincidence with one of the face nodes
  //
  // KHP: Couldn't we avoid the sqrt here by squaring both
  // sides of the spatial_tolerance comparison?
  //
  if(spatial_tolerance_pre > 0) {
    for (i=0; i<3; ++i) {
      DataType dx = node_positions[i][0]-global_coords[0];
      DataType dy = node_positions[i][1]-global_coords[1];
      DataType dz = node_positions[i][2]-global_coords[2];
      DataType dd  = dx*dx+dy*dy+dz*dz;
      if (dd == 0 || sqrt(dd) < spatial_tolerance_pre) break;
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
    }
    if (i<3) return;
  }

  /*
                      2
                     /|\
                    / | \
                   /a1|a0\
                  /  /c\  \
                 / /     \ \
                //    a2   \\ 
               0-------------1
  */

  // Compute vectors from the point to each node
  DataType v_c_0[3] = { node_positions[0][0] - global_coords[0],
                        node_positions[0][1] - global_coords[1],
                        node_positions[0][2] - global_coords[2] };

  DataType v_c_1[3] = { node_positions[1][0] - global_coords[0],
                        node_positions[1][1] - global_coords[1],
                        node_positions[1][2] - global_coords[2] };

  DataType v_c_2[3] = { node_positions[2][0] - global_coords[0],
                        node_positions[2][1] - global_coords[1],
                        node_positions[2][2] - global_coords[2] };

  // Compute the normal (as v_0_1 x v_0_2)
  DataType v_0_1[3] = { node_positions[1][0] - node_positions[0][0],
                        node_positions[1][1] - node_positions[0][1],
                        node_positions[1][2] - node_positions[0][2] };
  DataType v_0_2[3] = { node_positions[2][0] - node_positions[0][0],
                        node_positions[2][1] - node_positions[0][1],
                        node_positions[2][2] - node_positions[0][2] };

  DataType Normal[3];
  acme::Cross(v_0_1, v_0_2, Normal);

  // Compute the areas
  DataType a0 = ScalarTripleProduct(v_c_1, v_c_2, Normal);
  DataType a1 = ScalarTripleProduct(v_c_2, v_c_0, Normal);
  DataType a2 = ScalarTripleProduct(v_c_0, v_c_1, Normal);

  // Compute the local coordinates
  DataType area = a0 + a1 + a2;
  DataType L1   = a0/area;
  DataType L2   = a1/area;
  if(spatial_tolerance_post > 0) {
    // If it's close to any of the edges, snap to it
    using std::min; using std::max;
    if (L1<1.0+spatial_tolerance_post) {
      L1 = min(L1, 1.0);
    }
    if (L1>-spatial_tolerance_post) {
      L1 = max(L1, 0.0);
    }
    if (L2<1.0+spatial_tolerance_post) {
      L2 = min(L2, 1.0);
    }
    if (L2>-spatial_tolerance_post) {
      L2 = max(L2, 0.0);
    }
  }

  local_coords[0] = L1;
  local_coords[1] = L2;
  local_coords[2] = 1.0-L1-L2;
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Compute_Global_Coords( DataType node_positions[3][3],
					      DataType local_coords[3],
					      DataType global_coords[3] )
{
  DataType N[3];
  int  nnodes=3;
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

template<typename DataType>
void ContactTriFaceL3<DataType>::Interpolate_Scalar( DataType  local_coords[3],
					   DataType  node_scalars[3],
					   DataType& interpolated_scalar )
{
  DataType N[3];
  int  nnodes=3;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

template<typename DataType>
void ContactTriFaceL3<DataType>::Interpolate_Vector( DataType local_coords[3],
					   DataType node_vectors[3][3],
					   DataType interpolated_vector[3] )
{
  DataType N[3];
  int  nnodes=3;
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

#endif  // #define ContactTriFaceL3_C_
