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


#include "allocators.h"
#include "ContactUtilities.h"
#include "ContactNode.h"
#include "ContactLineFaceL2.h"
#include "ContactFixedSizeAllocator.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

ContactLineFaceL2::ContactLineFaceL2( ContactFixedSizeAllocator* alloc,
                                      int Block_Index, 
				      int Index_in_Block, int key ) 
  : ContactFace<Real>( alloc, ContactSearch::LINEFACEL2,
                 Block_Index, Index_in_Block, key, 
                 nodes, edges, Node_Info, Edge_Info)
{}

ContactLineFaceL2* ContactLineFaceL2::new_ContactLineFaceL2(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactLineFaceL2].New_Frag())
             ContactLineFaceL2(alloc, Block_Index, Index_in_Block, key);
}

void ContactLineFaceL2_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  
  int n2 = sizeof(ContactLineFaceL2);
  alloc.Resize( n2,
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactLineFaceL2 allocator" );
}

ContactLineFaceL2::~ContactLineFaceL2() {}

void ContactLineFaceL2::Evaluate_Shape_Functions( Real* local_coords,
						  Real* shape_functions )
{
  shape_functions[0] = 0.5*(1-local_coords[0]);
  shape_functions[1] = 0.5*(1+local_coords[0]);
}

void ContactLineFaceL2::Compute_Normal( VariableHandle CURRENT_POSITION,
					VariableHandle FACE_NORMAL )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  Real* Position0    = node0->Variable(CURRENT_POSITION);
  Real* Position1    = node1->Variable(CURRENT_POSITION);
  Real* face_normal  = Variable(FACE_NORMAL);
  // The vector from node 0 to 1 can be written as (a,b).
  // The normal to this face is then (b,-a). 
  face_normal[0] = Position1[1] - Position0[1];
  face_normal[1] = Position0[0] - Position1[0];
  face_normal[2] = 0.0;

  acme::Normalize(face_normal);
}

void ContactLineFaceL2::Compute_Normal( VariableHandle POSITION,
					Real* normal, Real* local_coords )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  Real* Position0    = node0->Variable(POSITION);
  Real* Position1    = node1->Variable(POSITION);
  // The vector from node 0 to 1 can be written as (a,b).
  // The normal to this face is then (b,-a). 
  normal[0] = Position1[1] - Position0[1];
  normal[1] = Position0[0] - Position1[0];
  normal[2] = 0.0;

  acme::Normalize(normal);
}

void ContactLineFaceL2::Compute_Normal( Real** nodal_positions, 
					Real* local_coords, Real* normal )
{
  Real* Position0    = nodal_positions[0];
  Real* Position1    = nodal_positions[1];
  // The vector from node 0 to 1 can be written as (a,b).
  // The normal to this face is then (b,-a). 
  normal[0] = Position1[1] - Position0[1];
  normal[1] = Position0[0] - Position1[0];
  normal[2] = 0.0;

  acme::Normalize(normal);
}

void ContactLineFaceL2::Compute_CharacteristicLength( VariableHandle CURRENT_POSITION,
				                      VariableHandle CHARACTERISTIC_LENGTH )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  
  Real* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  Real  dx = Position1[0]-Position0[0];
  Real  dy = Position1[1]-Position0[1];
  Real  dz = Position1[2]-Position0[2];
  *characteristiclength = std::sqrt(dx*dx + dy*dy + dz*dz);
}

void ContactLineFaceL2::Compute_Centroid( VariableHandle CURRENT_POSITION,
					  VariableHandle CENTROID )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  
  Real* centroid = Variable(CENTROID);
  centroid[0] = (Position0[0]+Position1[0])/2;
  centroid[1] = (Position0[1]+Position1[1])/2;
}

void ContactLineFaceL2::Compute_Edge_Normal( VariableHandle, 
					     VariableHandle,
				      	     int, Real*)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "Error: Edges don't exist in 2D" << std::endl;
#endif
  POSTCONDITION( 0 );
}

void ContactLineFaceL2::Compute_Local_Coordinates( Real Config_Param, 
						   VariableHandle POSITION0,
						   VariableHandle POSITION1, 
						   VariableHandle FACE_NORMAL,
						   Real* global_coords, 
						   Real* local_coords )
{
  Real node_positions[2][3];
  if (Config_Param == 0.0) {
    for (int i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION0);
      for (int j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (int i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION1);
      for (int j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    Real alpha = 1.0 - Config_Param, beta = Config_Param;
    for (int i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position0 = Node(i)->Variable(POSITION0);
      Real* node_position1 = Node(i)->Variable(POSITION1);
      for (int j=0; j<3; ++j) {
        node_positions[i][j] = alpha*node_position0[j]+beta*node_position1[j];
      }
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactLineFaceL2::Compute_Local_Coordinates( VariableHandle POSITION, 
						   Real* global_coords, 
						   Real* local_coords )
{
  Real node_positions[2][3];
  for (int i=0; i<Nodes_Per_Face(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactLineFaceL2::Compute_Global_Coordinates( VariableHandle POSITION,
						    Real* local_coords,
						    Real* global_coords )
{
  Real N[2];
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Evaluate_Shape_Functions( local_coords, N );
  for( int i=0 ; i<2 ; ++i){
    Real* node_position = Node(i)->Variable(POSITION);
    global_coords[0]   += N[i]*node_position[0];
    global_coords[1]   += N[i]*node_position[1];
    global_coords[2]   += N[i]*node_position[2];
  }
}

void ContactLineFaceL2::Get_Close_Edges( Real* local_coords, int& number,
					 int& edge_1_id, int& edge_2_id )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactLineFaceL2::Get_Close_Edges not yet implemented" << std::endl;
#endif
  POSTCONDITION( 0 );
}

bool ContactLineFaceL2::Is_Inside_Face( Real* local_coords )
{
  PRECONDITION( 0 );
  return false;
}

ContactFace<Real>* ContactLineFaceL2::Neighbor( Real* local_coords )
{
  PRECONDITION( 0 );
  return (ContactFace<Real>*) NULL;
}

void ContactLineFaceL2::FacetDecomposition(int& nfacets, 
                                           Real* coordinates0,
					   Real* normals0, 
					   VariableHandle POSITION0,
                                           Real* coordinates1, 
					   Real* normals1, 
					   VariableHandle POSITION1,
                                           Real* coordinates2, 
					   Real* normals2, 
					   VariableHandle POSITION2)
{
  nfacets = 1;
  for(int j=0 ; j<2 ; ++j){
    Real* node_position = Node(j)->Variable(POSITION0);
    coordinates0[0+3*j] = node_position[0];
    coordinates0[1+3*j] = node_position[1];
    coordinates0[2+3*j] = node_position[2];
  }
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  Real* Position0 = node0->Variable(POSITION0);
  Real* Position1 = node1->Variable(POSITION0);
  
  // The vector from node 0 to 1 can be written as (a,b).
  // The normal to this face is then (b,-a). 
  normals0[0] = Position1[1] - Position0[1];
  normals0[1] = Position0[0] - Position1[0];
  normals0[2] = 0.0;

  acme::Normalize(normals0);

  if (coordinates1!=NULL) {
    for(int j=0 ; j<2 ; ++j){
      Real* node_position = Node(j)->Variable(POSITION1);
      coordinates1[0+3*j] = node_position[0];
      coordinates1[1+3*j] = node_position[1];
      coordinates1[2+3*j] = node_position[2];
    }
    node0 = Node(0);
    node1 = Node(1);
    Position0 = node0->Variable(POSITION1);
    Position1 = node1->Variable(POSITION1);
    
    // The vector from node 0 to 1 can be written as (a,b). 
    // The normal to this face is then (b,-a). 
    normals1[0] = Position1[1] - Position0[1];
    normals1[1] = Position0[0] - Position1[0];
    normals1[2] = 0.0;

    acme::Normalize(normals1);
  }
}

void 
ContactLineFaceL2::FacetStaticRestriction(int nfacets, Real* coordinates, 
                                          Real* normals, Real* ctrcl_facets, 
                                          Real* ctrcl)
{
  for (int i=0; i<12*nfacets; ++i) {
	ctrcl[i] = ctrcl_facets[i];
  }
}

void 
ContactLineFaceL2::FacetDynamicRestriction(int nfacets,
                                           Real* ctrcl_facets, 
                                           Real* ctrcl)
{
  for (int i=0; i<LENGTH*nfacets; ++i) {
	ctrcl[i] = ctrcl_facets[i];
  }
}


void ContactLineFaceL2::Smooth_Normal( VariableHandle, VariableHandle, 
				       VariableHandle, VariableHandle,
			ContactSearch::Smoothing_Resolution resolution,
				       Real, Real*, Real*,Real )
{
  POSTCONDITION( 0 );
}

int ContactLineFaceL2::FaceEdge_Intersection(VariableHandle POSITION,
                                             ContactEdge<Real>* edge, Real* coords)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactLineFaceL2::FaceEdge_Intersection not yet implemented\n";
#endif
  POSTCONDITION( 0 );
  return 0;
}
  
void 
ContactLineFaceL2::Compute_Node_Areas(VariableHandle, VariableHandle, Real*)
{
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

void ContactLineFaceL2::Compute_Shape_Functions( Real* local_coords,
						 Real* shape_functions )
{
  shape_functions[0] = 0.5*(1.0-local_coords[0]);
  shape_functions[1] = 0.5*(1.0+local_coords[0]);
}

void ContactLineFaceL2::Compute_Shape_Derivatives( Real* local_coords,
						   Real shape_derivs[2][2] )
{
  shape_derivs[0][0] = -0.5;
  shape_derivs[0][1] =  0.5;
  shape_derivs[1][0] =  0.5;
  shape_derivs[1][1] = -0.5;
}

void 
ContactLineFaceL2::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3], 
					 Real global_coords[3],
					 Real local_coords[3] )
{
  int  i;
  Real spatial_tolerance = 1.0e-10;
  //
  // check for coincidence with one of the face nodes
  //
  for (i=0; i<2; ++i) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real d  = std::sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] = 0.0;
    local_coords[1] = 0.0;
    local_coords[2] = 0.0;
    break;
  case 1:
    local_coords[0] = 1.0;
    local_coords[1] = 0.0;
    local_coords[2] = 0.0;
    break;
  }
  if (i<2) return;

  Real dx = node_positions[0][0]-global_coords[0];
  Real dy = node_positions[0][1]-global_coords[1];
  Real dz = node_positions[0][2]-global_coords[2];
  Real Dx = node_positions[0][0]-node_positions[1][0];
  Real Dy = node_positions[0][1]-node_positions[1][1];
  Real Dz = node_positions[0][2]-node_positions[1][2];
  Real d  = std::sqrt(dx*dx+dy*dy+dz*dz);
  Real D  = std::sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
  local_coords[0] = d/D;
  local_coords[1] = 1.0-local_coords[0];
  local_coords[2] = 0.0;
}

void ContactLineFaceL2::Compute_Global_Coords( Real node_positions[2][3],
					       Real local_coords[2],
					       Real global_coords[3] )
{
  Real N[2];
  int  nnodes=2;
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

void  ContactLineFaceL2::Interpolate_Scalar( Real  local_coords[2],
					     Real  node_scalars[2],
					     Real& interpolated_scalar )
{
  Real N[2];
  int  nnodes=2;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void  ContactLineFaceL2::Interpolate_Vector( Real local_coords[2],
					     Real node_vectors[2][3],
					     Real interpolated_vector[3] )
{
  Real N[2];
  int  nnodes=2;
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

