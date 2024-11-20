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
#include "ContactNode.h"
#include "ContactLineEdgeL2.h"
#include "ContactQuadFaceL4.h"
#include "ContactHexElementL8.h"
#include "ContactFixedSizeAllocator.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

/*************************************************************************
 
         C A R T E S I A N    H E X   E L E M E N T

*************************************************************************/

ContactCartesianHexElementL8::ContactCartesianHexElementL8( 
                                     ContactFixedSizeAllocator* alloc,
                                     int Block_Index, 
				     int Index_in_Block, int key ) 
  : ContactElement( alloc, ContactSearch::CARTESIANHEXELEMENTL8,
		    Block_Index, Index_in_Block, key) 
{}

ContactCartesianHexElementL8* 
ContactCartesianHexElementL8::new_ContactCartesianHexElementL8(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactCartesianHexElementL8].New_Frag())
             ContactCartesianHexElementL8(alloc, Block_Index, 
                                          Index_in_Block, key);
}

void ContactCartesianHexElementL8_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactCartesianHexElementL8),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactCartesianHexElementL8 allocator" );
}

ContactCartesianHexElementL8::~ContactCartesianHexElementL8() {}

void
ContactCartesianHexElementL8::TetDice( int &ntets, Real thex[][4][3], 
                                       VariableHandle POSITION ) 
{
  switch(ntets) {
  case 5:
    {
      // works fine for planar faces
      ntets = 5;
      int node_index[5][4] = {{0,1,2,5},
			      {2,3,0,7},
			      {7,6,5,2},
			      {5,4,7,0},
			      {0,5,2,7}};
      
      for (int i=0; i<5; ++i) {
	for (int j=0; j<4; ++j) {
	  Real* position = nodes[node_index[i][j]]->Variable(POSITION);
	  for (int k=0; k<3; ++k) {
	    thex[i][j][k] = position[k];
	  }
	}
      }        
      break;
    }
  case 24:
    {
      // Use 24 tets to dice hex
      Real*     node_position;
      Real      centroid3[3];
      Real      centroid4[3];
      int       tetcount = 0;
      int node_index[6][4] = {{0,1,5,4},
			      {1,2,6,5},
			      {2,3,7,6},
			      {0,4,7,3},
			      {0,3,2,1},
			      {4,5,6,7}};
      
      // construct hex-centroid -- this is used for every face decomposition
      centroid4[0] = 0.;
      centroid4[1] = 0.;
      centroid4[2] = 0.;
      for( int i=0; i<8; ++i ) { // 8 nodes per hex
	node_position = nodes[i]->Variable(POSITION);
	centroid4[0] += node_position[0];
	centroid4[1] += node_position[1];
	centroid4[2] += node_position[2];
      }
      centroid4[0] = 0.125*centroid4[0];
      centroid4[1] = 0.125*centroid4[1];
      centroid4[2] = 0.125*centroid4[2];
      
      // now compute decomposition
      tetcount = 0;
      for( int j=0; j<6; ++j ) {// 6 faces on a hex
	centroid3[0] = 0.;
	centroid3[1] = 0.;
	centroid3[2] = 0.;
	for( int i=0; i<4; ++i ) { // 4 nodes per quad face
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  centroid3[0] += node_position[0];
	  centroid3[1] += node_position[1];
	  centroid3[2] += node_position[2];
	}
	centroid3[0] = 0.25*centroid3[0];
	centroid3[1] = 0.25*centroid3[1];
	centroid3[2] = 0.25*centroid3[2];
	
	for(int i = 0; i < 4; i ++){
	  int iplus = (i+1)%4;
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  thex[tetcount][0][0] = node_position[0];
	  thex[tetcount][0][1] = node_position[1];
	  thex[tetcount][0][2] = node_position[2];
	  
	  node_position = nodes[node_index[j][iplus]]->Variable(POSITION);
	  thex[tetcount][1][0] = node_position[0];
	  thex[tetcount][1][1] = node_position[1];
	  thex[tetcount][1][2] = node_position[2];
	  
	  thex[tetcount][2][0] = centroid4[0];
	  thex[tetcount][2][1] = centroid4[1];
	  thex[tetcount][2][2] = centroid4[2];
	  
	  thex[tetcount][3][0] = centroid3[0];
	  thex[tetcount][3][1] = centroid3[1];
	  thex[tetcount][3][2] = centroid3[2];
	  
	  ++tetcount;
	}
      } // end of loop on hex faces
      break;
    }
  default:
    {
      //should never get here and don't want to test for every interaction
      //so just do hard abort
      std::exit(1);
    }
  }
}

void
ContactCartesianHexElementL8::Compute_Volume( VariableHandle NODE_POSITION,
					      VariableHandle ELEMENT_VOLUME )
{
  using std::abs;
  Real* volume = Variable(ELEMENT_VOLUME);
  Real* n0_pos = Node(0)->Variable(NODE_POSITION);
  Real* n6_pos = Node(6)->Variable(NODE_POSITION);
  *volume = abs( (n0_pos[0]-n6_pos[0]) * (n0_pos[1]-n6_pos[1]) * 
		 (n0_pos[2]-n6_pos[2]) );
}

void
ContactCartesianHexElementL8::Compute_Centroid( VariableHandle NODE_POSITION,
					        VariableHandle CENTROID )
{
  Real* centroid = Variable(CENTROID);
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  for (int i=0; i<8; ++i) {
    Real* position = Node(i)->Variable(NODE_POSITION);
    centroid[0] += position[0];
    centroid[1] += position[1];
    centroid[2] += position[2];
  }
  centroid[0] /= 8.0;
  centroid[1] /= 8.0;
  centroid[2] /= 8.0;
}


bool
ContactCartesianHexElementL8::Is_Local_Coordinates_Inside_Element( Real* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= -1.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

bool
ContactCartesianHexElementL8::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance )
{
  Real low_coord  = -(1.+tolerance);
  Real high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

void ContactCartesianHexElementL8::Evaluate_Shape_Functions( Real* local_coords,
						    Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void ContactCartesianHexElementL8::Compute_Global_Coordinates( VariableHandle POSITION,
							       Real* local_coords,
							       Real* global_coords )
{
  Real node_positions[8][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactCartesianHexElementL8::Compute_Local_Coordinates( VariableHandle POSITION,
						     Real* global_coords,
						     Real* local_coords )
{
  int i, j;
  Real node_positions[8][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactCartesianHexElementL8::Interpolate_Scalar_Value( Real* local_coords,
							     Real* node_scalars,
							     Real& interpolated_scalar )
{
  Interpolate_Scalar( local_coords, node_scalars, interpolated_scalar );
  return;
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

void ContactCartesianHexElementL8::Compute_Shape_Functions( Real* local_coords,
							    Real* shape_functions )
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

void 
ContactCartesianHexElementL8::Compute_Local_Coords( Real node_positions[8][3], 
						    Real global_coords[3],
						    Real local_coords[3] )
{
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i;
  int  nnodes=8;


  Real spatial_tolerance = 1.0e-10;
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
  if (i<nnodes) return;

  //
  local_coords[0] = 2.*(global_coords[0] - node_positions[0][0])/
    (node_positions[1][0] - node_positions[0][0]) - 1.;
  local_coords[1] = 2.*(global_coords[1] - node_positions[0][1])/
    (node_positions[3][1] - node_positions[0][1]) - 1.;
  local_coords[2] = 2.*(global_coords[2] - node_positions[0][2])/
    (node_positions[4][2] - node_positions[0][2]) - 1.;
}

void ContactCartesianHexElementL8::Compute_Global_Coords( Real node_positions[8][3],
							  Real local_coords[3],
							  Real global_coords[3] )
{
  Real N[8];
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

void  ContactCartesianHexElementL8::Interpolate_Scalar( Real  local_coords[3],
							Real  node_scalars[8],
							Real& interpolated_scalar )
{
  Real N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}


/*************************************************************************
 
         T R I L I N E A R    H E X   E L E M E N T

*************************************************************************/

ContactHexElementL8::ContactHexElementL8( ContactFixedSizeAllocator* alloc, 
                                          int Block_Index, 
				          int Index_in_Block, int key ) 
  : ContactElement( alloc, ContactSearch::HEXELEMENTL8, 
                    Block_Index, Index_in_Block,
		    key) 
{}

ContactHexElementL8* 
ContactHexElementL8::new_ContactHexElementL8(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactHexElementL8].New_Frag())
             ContactHexElementL8(alloc, Block_Index, Index_in_Block, key);
}

void ContactHexElementL8_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactHexElementL8),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactHexElementL8 allocator" );
}

ContactHexElementL8::~ContactHexElementL8() {}

void
ContactHexElementL8::TetDice( int &ntets, Real thex[5][4][3], 
                              VariableHandle POSITION) 
{
  switch(ntets) {
  case 5:
    {
      // works fine for planar faces
      int node_index[5][4] = {{0,1,2,5},
			      {2,3,0,7},
			      {7,6,5,2},
			      {5,4,7,0},
			      {0,5,2,7}};
      
      for (int i=0; i<5; ++i) {
	for (int j=0; j<4; ++j) {
	  Real* position = nodes[node_index[i][j]]->Variable(POSITION);
	  for (int k=0; k<3; ++k) {
	    thex[i][j][k] = position[k];
	  }
	}
      }        
      break;
    }
  case 24:
    {
      // Use 24 tets to dice hex
      Real*     node_position;
      Real      centroid3[3];
      Real      centroid4[3];
      int       tetcount = 0;
      int node_index[6][4] = {{0,1,5,4},
			      {1,2,6,5},
			      {2,3,7,6},
			      {0,4,7,3},
			      {0,3,2,1},
			      {4,5,6,7}};
      
      // construct hex-centroid -- this is used for every face decomposition
      centroid4[0] = 0.;
      centroid4[1] = 0.;
      centroid4[2] = 0.;
      for( int i=0; i<8; ++i ) { // 8 nodes per hex
	node_position = nodes[i]->Variable(POSITION);
	centroid4[0] += node_position[0];
	centroid4[1] += node_position[1];
	centroid4[2] += node_position[2];
      }
      centroid4[0] = 0.125*centroid4[0];
      centroid4[1] = 0.125*centroid4[1];
      centroid4[2] = 0.125*centroid4[2];
      
      // now compute decomposition
      tetcount = 0;
      for( int j=0; j<6; ++j ) {// 6 faces on a hex
	centroid3[0] = 0.;
	centroid3[1] = 0.;
	centroid3[2] = 0.;
	for( int i=0; i<4; ++i ) { // 4 nodes per quad face
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  centroid3[0] += node_position[0];
	  centroid3[1] += node_position[1];
	  centroid3[2] += node_position[2];
	}
	centroid3[0] = 0.25*centroid3[0];
	centroid3[1] = 0.25*centroid3[1];
	centroid3[2] = 0.25*centroid3[2];
	
	for(int i = 0; i < 4; i ++){
	  int iplus = (i+1)%4;
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  thex[tetcount][0][0] = node_position[0];
	  thex[tetcount][0][1] = node_position[1];
	  thex[tetcount][0][2] = node_position[2];
	  
	  node_position = nodes[node_index[j][iplus]]->Variable(POSITION);
	  thex[tetcount][1][0] = node_position[0];
	  thex[tetcount][1][1] = node_position[1];
	  thex[tetcount][1][2] = node_position[2];
	  
	  thex[tetcount][2][0] = centroid4[0];
	  thex[tetcount][2][1] = centroid4[1];
	  thex[tetcount][2][2] = centroid4[2];
	  
	  thex[tetcount][3][0] = centroid3[0];
	  thex[tetcount][3][1] = centroid3[1];
	  thex[tetcount][3][2] = centroid3[2];
	  
	  ++tetcount;
	}
      } // end of loop on hex faces
      break;
    }
  default:
    {
      //should never get here and don't want to test for every interaction
      //so just do hard abort
      std::exit(1);
    }

  }
}

void ContactHexElementL8::Compute_Volume( VariableHandle NODE_POSITION,
					  VariableHandle ELEMENT_VOLUME )
{
  Real* volume = Variable(ELEMENT_VOLUME);

  Real* n0 = Node(0)->Variable(NODE_POSITION);
  Real* n1 = Node(1)->Variable(NODE_POSITION);
  Real* n2 = Node(2)->Variable(NODE_POSITION);
  Real* n3 = Node(3)->Variable(NODE_POSITION);
  Real* n4 = Node(4)->Variable(NODE_POSITION);
  Real* n5 = Node(5)->Variable(NODE_POSITION);
  Real* n6 = Node(6)->Variable(NODE_POSITION);
  Real* n7 = Node(7)->Variable(NODE_POSITION);
  
  Real x1 = n0[0];
  Real x2 = n1[0];
  Real x3 = n2[0];
  Real x4 = n3[0];
  Real x5 = n4[0];
  Real x6 = n5[0];
  Real x7 = n6[0];
  Real x8 = n7[0];

  Real y1 = n0[1];
  Real y2 = n1[1];
  Real y3 = n2[1];
  Real y4 = n3[1];
  Real y5 = n4[1];
  Real y6 = n5[1];
  Real y7 = n6[1];
  Real y8 = n7[1];

  Real z1 = n0[2];
  Real z2 = n1[2];
  Real z3 = n2[2];
  Real z4 = n3[2];
  Real z5 = n4[2];
  Real z6 = n5[2];
  Real z7 = n6[2];
  Real z8 = n7[2];

  Real rx0 = (y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-
             (z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+
             y8*(z4-z5));
  Real rx1 = (y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-
             (z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+
             y5*(z1-z6));
  Real rx2 = (y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-
             (z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+
             y6*(z2-z7));
  Real rx3 = (y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-
             (z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+
             y7*(z3-z8));
  Real rx4 = (y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-
             (z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+
             y2*(z6-z1));
  Real rx5 = (y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-
             (z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+
             y3*(z7-z2));
  Real rx6 = (y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-
             (z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+
             y4*(z8-z3));
  Real rx7 = (y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-
             (z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+
             y1*(z5-z4));

  *volume = (x1 * rx0 +
	     x2 * rx1 +
	     x3 * rx2 +
	     x4 * rx3 +
	     x5 * rx4 +
	     x6 * rx5 +
	     x7 * rx6 +
	     x8 * rx7) / 12.0;
}

void
ContactHexElementL8::Compute_Centroid( VariableHandle NODE_POSITION,
				       VariableHandle CENTROID )
{
  Real* centroid = Variable(CENTROID);
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  for (int i=0; i<8; ++i) {
    Real* position = Node(i)->Variable(NODE_POSITION);
    centroid[0] += position[0];
    centroid[1] += position[1];
    centroid[2] += position[2];
  }
  centroid[0] /= 8.0;
  centroid[1] /= 8.0;
  centroid[2] /= 8.0;
}

bool
ContactHexElementL8::Is_Local_Coordinates_Inside_Element( Real* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= -1.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

bool
ContactHexElementL8::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance )
{
  Real low_coord  = -(1.+tolerance);
  Real high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

void ContactHexElementL8::Evaluate_Shape_Functions( Real* local_coords,
						    Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void ContactHexElementL8::Compute_Global_Coordinates( VariableHandle POSITION,
						      Real* local_coords,
						      Real* global_coords )
{
  Real node_positions[8][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactHexElementL8::Compute_Local_Coordinates( VariableHandle POSITION,
						     Real* global_coords,
						     Real* local_coords )
{
  int i, j;
  Real node_positions[8][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}


void ContactHexElementL8::Interpolate_Scalar_Value( Real* local_coords,
						    Real* node_scalars,
						    Real& interpolated_scalar )
{
  Interpolate_Scalar( local_coords, node_scalars, interpolated_scalar );
  return;
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

void ContactHexElementL8::Compute_Shape_Functions( Real* local_coords,
						   Real* shape_functions )
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

void ContactHexElementL8::Compute_Shape_Derivatives( Real* local_coords,
						     Real shape_derivs[3][8] )
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

void 
ContactHexElementL8::Compute_Local_Coords( Real node_positions[8][3], 
					   Real global_coords[3],
					   Real local_coords[3] )
{
  using std::abs;
  using std::min;
  using std::max;

  int  i, j;
  int  nnodes=8;
  Real spatial_tolerance = 1.0e-10;

  // are we on a node?
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
  if (i<nnodes) return;
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=400;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real u, u0=0.0, u1, du;
  Real v, v0=0.0, v1, dv;
  Real w, w0=0.0, w1, dw;
  Real f[3], J[3][3], invJ[3][3];
  Real shape_derivatives[3][8];
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
    u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
    v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
    w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
    du = abs(u1-u0);
    dv = abs(v1-v0);
    dw = abs(w1-w0);
    u0 = u1;
    v0 = v1;
    w0 = w1;
    if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
    ++iterations;
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if (!converged) {
    std::cerr << "ContactHexElementL8::Compute_Local_Coordinates() did not converge" 
	 << std::endl;
    std::cerr << "                     after "<<max_iterations
         << " iterations:  du = "<<du
         <<";  dv = "<<dv
         <<";  dw = "<<dw<<std::endl;
  }
#endif
  POSTCONDITION(converged);
  // If it's close to any of the edges, snap to it
  if (abs(u0)<1.0+spatial_tolerance) {
    u0 = min(u0, 1.0);
    u0 = max(u0,-1.0);
  }
  if (abs(v0)<1.0+spatial_tolerance) {
    v0 = min(v0, 1.0);
    v0 = max(v0,-1.0);
  }
  if (abs(w0)<1.0+spatial_tolerance) {
    w0 = min(w0, 1.0);
    w0 = max(w0,-1.0);
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = w0;
}

void ContactHexElementL8::Compute_Global_Coords( Real node_positions[8][3],
						 Real local_coords[3],
						 Real global_coords[3] )
{
  Real N[8];
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

void  ContactHexElementL8::Interpolate_Scalar( Real  local_coords[3],
					       Real  node_scalars[8],
					       Real& interpolated_scalar )
{
  Real N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

