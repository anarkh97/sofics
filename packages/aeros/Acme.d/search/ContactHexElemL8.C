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


#ifndef ContactHexElemL8_C_
#define ContactHexElemL8_C_

#include "allocators.h"
#include "ContactNode.h"
#include "ContactLineEdgeL2.h"
#include "ContactQuadFaceL4.h"
#include "ContactHexElementL8.h"
#include "ContactFixedSizeAllocator.h"
#include "contact_tolerances.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

template<typename DataType>
ContactHexElemL8<DataType>::ContactHexElemL8( int Block_Index, 
				    int Host_Index_in_Block, int key ) 
  : ContactElem<DataType>( ContactSearch::HEXELEML8,Block_Index,Host_Index_in_Block,key) 
{}

template<typename DataType>
ContactHexElemL8<DataType>* ContactHexElemL8<DataType>::new_ContactHexElemL8(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key)
{
  return new (alloc.New_Frag())
             ContactHexElemL8<DataType>(Block_Index, Host_Index_in_Block, key);
}

template<typename DataType>
void ContactHexElemL8_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactHexElemL8<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactHexElemL8<DataType> allocator" );
}

template<typename DataType>
ContactHexElemL8<DataType>::~ContactHexElemL8() {}

template<typename DataType>
void ContactHexElemL8<DataType>::BuildTopology(int nID, int eID, int fID,
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
  edge->ConnectNode(1, Node(5));
  edge = Edge(2);
  edge->ConnectNode(0, Node(5));
  edge->ConnectNode(1, Node(4));
  edge = Edge(3);
  edge->ConnectNode(0, Node(4));
  edge->ConnectNode(1, Node(0));
  edge = Edge(4);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(2));
  edge = Edge(5);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(6));
  edge = Edge(6);
  edge->ConnectNode(0, Node(6));
  edge->ConnectNode(1, Node(5));
  edge = Edge(7);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(3));
  edge = Edge(8);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(7));
  edge = Edge(9);
  edge->ConnectNode(0, Node(7));
  edge->ConnectNode(1, Node(6));
  edge = Edge(10);
  edge->ConnectNode(0, Node(4));
  edge->ConnectNode(1, Node(7));
  edge = Edge(11);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(0));
  
  NextID = fID;
  for( i=0 ; i<Faces_Per_Element() ; ++i ) {
    face = ContactQuadFaceL4<DataType>::new_ContactQuadFaceL4(allocators, ++NextID );
    this->ConnectFace(i, face);
  }
  face = Face(0);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(1));
  face->ConnectNode(2, Node(5));
  face->ConnectNode(3, Node(4));
  face = Face(1);
  face->ConnectNode(0, Node(1));
  face->ConnectNode(1, Node(2));
  face->ConnectNode(2, Node(6));
  face->ConnectNode(3, Node(5));
  face = Face(2);
  face->ConnectNode(0, Node(2));
  face->ConnectNode(1, Node(3));
  face->ConnectNode(2, Node(7));
  face->ConnectNode(3, Node(6));
  face = Face(3);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(4));
  face->ConnectNode(2, Node(7));
  face->ConnectNode(3, Node(3));
  face = Face(4);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(3));
  face->ConnectNode(2, Node(2));
  face->ConnectNode(3, Node(1));
  face = Face(5);
  face->ConnectNode(0, Node(4));
  face->ConnectNode(1, Node(5));
  face->ConnectNode(2, Node(6));
  face->ConnectNode(3, Node(7));
}

template<typename DataType>
void ContactHexElemL8<DataType>::DeleteTopology(ContactFixedSizeAllocator* allocators)
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
  for( i=0 ; i<Faces_Per_Element() ; ++i ) {
    ContactFace<DataType>* face = Face(i);
    face->~ContactFace<DataType>();
    allocators[ContactSearch::ALLOC_ContactQuadFaceL4].Delete_Frag(face);
  }
}

template<typename DataType>
void ContactHexElemL8<DataType>::UpdateTopology(ContactFace<DataType>* face, 
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
bool
ContactHexElemL8<DataType>::Is_Local_Coordinates_Inside_Element( DataType* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= -1.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

template<typename DataType>
bool
ContactHexElemL8<DataType>::Is_Local_Coordinates_Near_Element( DataType* local_coords, DataType tolerance )
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
void ContactHexElemL8<DataType>::Evaluate_Shape_Functions( DataType* local_coords,
						  DataType* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Global_Coordinates( VariableHandle POSITION,
						   DataType* local_coords,
						   DataType* global_coords )
{
  DataType node_positions[8][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    DataType* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Local_Coordinates( DataType Config_Param,
						  VariableHandle POSITION0, 
						  VariableHandle POSITION1, 
						  VariableHandle FACE_NORMAL,
						  DataType* global_coords,
						  DataType* local_coords )
{
  int i, j;
  DataType node_positions[8][3];
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
bool ContactHexElemL8<DataType>::Compute_Local_Coordinates( VariableHandle POSITION,
						  DataType* global_coords,
						  DataType* local_coords )
{
  int i, j;
  DataType node_positions[8][3];
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
/* on a linear H8 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary H8 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Shape_Functions( DataType* local_coords,
						          DataType* shape_functions )
{
  shape_functions[0] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[1] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[2] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[3] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[4] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[5] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[6] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[7] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Shape_Derivatives( DataType* local_coords,
						  DataType shape_derivs[3][8] )
{
  shape_derivs[0][0] = -0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][1] =  0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][2] =  0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][3] = -0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][4] = -0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][5] =  0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][6] =  0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][7] = -0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);

  shape_derivs[1][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][2] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][3] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][4] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][5] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);

  shape_derivs[2][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][2] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][3] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][4] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][5] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
}

template<>
inline bool 
ContactHexElemL8<Real>::Compute_Local_Coords( Real node_positions[8][3], 
					Real global_coords[3],
					Real local_coords[3] )
{
  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i, j;
  int  nnodes=8;
  if(spatial_tolerance_pre > 0) {
    // are we on a node?
    for (i=0; i<nnodes; ++i) {
      Real dx = node_positions[i][0]-global_coords[0];
      Real dy = node_positions[i][1]-global_coords[1];
      Real dz = node_positions[i][2]-global_coords[2];
      Real dd = dx*dx+dy*dy+dz*dz;
      if (dd == 0 || sqrt(dd) < spatial_tolerance_pre) break;
    }
    switch (i) {
    case 0:
      local_coords[0] = -1.0;
      local_coords[1] = -1.0;
      local_coords[2] = -1.0;
      break;
    case 1:
      local_coords[0] =  1.0;
      local_coords[1] = -1.0;
      local_coords[2] = -1.0;
      break;
    case 2:
      local_coords[0] =  1.0;
      local_coords[1] =  1.0;
      local_coords[2] = -1.0;
      break;
    case 3:
      local_coords[0] = -1.0;
      local_coords[1] =  1.0;
      local_coords[2] = -1.0;
      break;
    case 4:
      local_coords[0] = -1.0;
      local_coords[1] = -1.0;
      local_coords[2] =  1.0;
      break;
    case 5:
      local_coords[0] =  1.0;
      local_coords[1] = -1.0;
      local_coords[2] =  1.0;
      break;
    case 6:
      local_coords[0] =  1.0;
      local_coords[1] =  1.0;
      local_coords[2] =  1.0;
      break;
    case 7:
      local_coords[0] = -1.0;
      local_coords[1] =  1.0;
      local_coords[2] =  1.0;
      break;
    }
    if (i<nnodes) return true;
  }
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=400;
  bool converged = false;
  Real u, u0=0.0, u1, du;
  Real v, v0=0.0, v1, dv;
  Real w, w0=0.0, w1, dw;
  Real f[3], J[3][3], invJ[3][3];
  Real shape_derivatives[3][8];
  Real residualNorm2, initialResidualNorm2;
  while (!converged && iterations<max_iterations) {
    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = w0;
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( local_coords, shape_derivatives );
    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<8; ++j) {
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

    Compute_Global_Coords( node_positions, local_coords, f );
    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    residualNorm2 = u*u+v*v+w*w;
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
    std::cerr << "ContactHexElemL8::Compute_Local_Coordinates() did not converge" 
	 << std::endl;
    std::cerr << "                     after "<<max_iterations
         << " iterations:  du = "<<du
         <<";  dv = "<<dv
         <<";  dw = "<<dw<<std::endl;
  }
#endif
  POSTCONDITION(converged);
  /*if(!converged) {
    std::cerr << " *** WARNING: ContactHexElemL8<Real>::Compute_Local_Coordinates() did not converge: initialResidualNorm2 = "
              << initialResidualNorm2 << ", residualNorm2 = " << residualNorm2 << std::endl;
  }*/
  if(spatial_tolerance_post > 0) {
    // If it's close to any of the edges, snap to it
    if (abs(u0)<1.0+spatial_tolerance_post) {
      u0 = min(u0, 1.0);
      u0 = max(u0,-1.0);
    }
    if (abs(v0)<1.0+spatial_tolerance_post) {
      v0 = min(v0, 1.0);
      v0 = max(v0,-1.0);
    }
    if (abs(w0)<1.0+spatial_tolerance_post) {
      w0 = min(w0, 1.0);
      w0 = max(w0,-1.0);
    }
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = w0;
  return converged;
}

template<typename ActiveScalar>
bool
ContactHexElemL8<ActiveScalar>::Compute_Local_Coords( ActiveScalar active_node_positions[8][3],
                                                      ActiveScalar active_global_coords[3],
                                                      ActiveScalar active_local_coords[3] )
{
  int  i, j;
  const int  nnodes=8;

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
    // are we on a node?
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
      local_coords[2] = -1.0;
      break;
    case 1:
      local_coords[0] =  1.0;
      local_coords[1] = -1.0;
      local_coords[2] = -1.0;
      break;
    case 2:
      local_coords[0] =  1.0;
      local_coords[1] =  1.0;
      local_coords[2] = -1.0;
      break;
    case 3:
      local_coords[0] = -1.0;
      local_coords[1] =  1.0;
      local_coords[2] = -1.0;
      break;
    case 4:
      local_coords[0] = -1.0;
      local_coords[1] = -1.0;
      local_coords[2] =  1.0;
      break;
    case 5:
      local_coords[0] =  1.0;
      local_coords[1] = -1.0;
      local_coords[2] =  1.0;
      break;
    case 6:
      local_coords[0] =  1.0;
      local_coords[1] =  1.0;
      local_coords[2] =  1.0;
      break;
    case 7:
      local_coords[0] = -1.0;
      local_coords[1] =  1.0;
      local_coords[2] =  1.0;
      break;
    }
    if (i<nnodes) {
      active_local_coords[0] = local_coords[0];
      active_local_coords[1] = local_coords[1];
      active_local_coords[2] = local_coords[2];
      return true;
    }
  }
  //
  // else use newton's method to iterate (values only)
  //
  int  iterations=0;
  int  max_iterations=400;
  bool converged = false;
  double u, u0=0.0, u1, du;
  double v, v0=0.0, v1, dv;
  double w, w0=0.0, w1, dw;
  double f[3], J[3][3], invJ[3][3];
  double shape_derivatives[3][8], shape_functions[8];
  double residualNorm2, initialResidualNorm2;
  double u0_copy = 0.0, v0_copy = 0.0, w0_copy = 0.0;
  while (!converged && iterations<max_iterations) {
    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = w0;

    // BUILD JACOBIAN AND INVERT
    shape_derivatives[0][0] = -0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
    shape_derivatives[0][1] =  0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
    shape_derivatives[0][2] =  0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
    shape_derivatives[0][3] = -0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
    shape_derivatives[0][4] = -0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
    shape_derivatives[0][5] =  0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
    shape_derivatives[0][6] =  0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);
    shape_derivatives[0][7] = -0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);

    shape_derivatives[1][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
    shape_derivatives[1][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
    shape_derivatives[1][2] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
    shape_derivatives[1][3] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
    shape_derivatives[1][4] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);
    shape_derivatives[1][5] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
    shape_derivatives[1][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
    shape_derivatives[1][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);

    shape_derivatives[2][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
    shape_derivatives[2][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
    shape_derivatives[2][2] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
    shape_derivatives[2][3] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
    shape_derivatives[2][4] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
    shape_derivatives[2][5] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
    shape_derivatives[2][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
    shape_derivatives[2][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);

    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<8; ++j) {
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
    shape_functions[0] = 0.125*(1.0-local_coords[0])*(1.0-local_coords[1])*(1.0-local_coords[2]);
    shape_functions[1] = 0.125*(1.0+local_coords[0])*(1.0-local_coords[1])*(1.0-local_coords[2]);
    shape_functions[2] = 0.125*(1.0+local_coords[0])*(1.0+local_coords[1])*(1.0-local_coords[2]);
    shape_functions[3] = 0.125*(1.0-local_coords[0])*(1.0+local_coords[1])*(1.0-local_coords[2]);
    shape_functions[4] = 0.125*(1.0-local_coords[0])*(1.0-local_coords[1])*(1.0+local_coords[2]);
    shape_functions[5] = 0.125*(1.0+local_coords[0])*(1.0-local_coords[1])*(1.0+local_coords[2]);
    shape_functions[6] = 0.125*(1.0+local_coords[0])*(1.0+local_coords[1])*(1.0+local_coords[2]);
    shape_functions[7] = 0.125*(1.0-local_coords[0])*(1.0+local_coords[1])*(1.0+local_coords[2]);
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
    residualNorm2 = u*u+v*v+w*w;
    if(iterations == 0) initialResidualNorm2 = residualNorm2;
    if(residualNorm2 < newton_tolerance*initialResidualNorm2 || residualNorm2 < abs_newton_tolerance) converged = true;
    else {
      u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
      v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
      w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
      du = abs(u1-u0);
      dv = abs(v1-v0);
      dw = abs(w1-w0);
      u0_copy = u0; v0_copy = v0; w0_copy = w0;
      u0 = u1;
      v0 = v1;
      w0 = w1;
      //if (du<newton_tolerance && dv<newton_tolerance && dw<newton_tolerance) converged = true;
      ++iterations;
    }
  }
  /*if (!converged) {
    std::cerr << " *** WARNING: ContactHexElemL8<ActiveScalar>::Compute_Local_Coords did not converge: initialResidualNorm2 = "
              << initialResidualNorm2 << ", residualNorm2 = " << residualNorm2 << std::endl;
  }*/
  {
    //
    // repeat last newton iteration to get derivatives
    //
    int  iterations=0;
    int  max_iterations=1;
    bool converged = false;
    ActiveScalar u, u0=u0_copy, u1, du;
    ActiveScalar v, v0=v0_copy, v1, dv;
    ActiveScalar w, w0=w0_copy, w1, dw;
    ActiveScalar f[3], J[3][3], invJ[3][3];
    ActiveScalar shape_derivatives[3][8];
    while (!converged && iterations<max_iterations) {
      active_local_coords[0] = u0;
      active_local_coords[1] = v0;
      active_local_coords[2] = w0;

      // BUILD JACOBIAN AND INVERT
      Compute_Shape_Derivatives( active_local_coords, shape_derivatives );
      for (i=0; i<3; ++i) {
        J[0][i] = 0.0;
        J[1][i] = 0.0;
        J[2][i] = 0.0;
        for (j=0; j<8; ++j) {
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
      if (abs(u0)<1.0+spatial_tolerance_post) {
        u0 = min(u0, 1.0);
        u0 = max(u0,-1.0);
      }
      if (abs(v0)<1.0+spatial_tolerance_post) {
        v0 = min(v0, 1.0);
        v0 = max(v0,-1.0);
      }
      if (abs(w0)<1.0+spatial_tolerance_post) {
        w0 = min(w0, 1.0);
        w0 = max(w0,-1.0);
      }
    }
    active_local_coords[0] = u0;
    active_local_coords[1] = v0;
    active_local_coords[2] = w0;
  }
  return converged;
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Global_Coords( DataType node_positions[8][3],
					      DataType local_coords[3],
					      DataType global_coords[3] )
{
  DataType N[8];
  int  nnodes=8;
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
void ContactHexElemL8<DataType>::Interpolate_Scalar( DataType  local_coords[3],
					    DataType  node_scalars[8],
					    DataType& interpolated_scalar )
{
  DataType N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

template<typename DataType>
void ContactHexElemL8<DataType>::Interpolate_Vector( DataType local_coords[3],
					    DataType node_vectors[8][3],
					    DataType interpolated_vector[3] )
{
  DataType N[8];
  int  nnodes=8;
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

#endif  // #define ContactHexElementL8_C_
