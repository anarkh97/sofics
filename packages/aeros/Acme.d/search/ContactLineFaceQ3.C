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
#include "ContactLineFaceQ3.h"
#include "ContactFixedSizeAllocator.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

ContactLineFaceQ3::ContactLineFaceQ3( ContactFixedSizeAllocator* alloc,
                                      int Block_Index, 
				      int Index_in_Block, int key ) 
  : ContactFace<Real>( alloc, ContactSearch::LINEFACEQ3,
                 Block_Index, Index_in_Block ,key, 
                 nodes, edges, Node_Info, Edge_Info)
{}

ContactLineFaceQ3* ContactLineFaceQ3::new_ContactLineFaceQ3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactLineFaceQ3].New_Frag())
             ContactLineFaceQ3(alloc, Block_Index, Index_in_Block, key);
}

void ContactLineFaceQ3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactLineFaceQ3),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactLineFaceQ3 allocator" );
}


ContactLineFaceQ3::~ContactLineFaceQ3( ) {}


void ContactLineFaceQ3::Evaluate_Shape_Functions( Real* local_coords,
					          Real* shape_functions )
{
  shape_functions[0] = 0.5*(local_coords[0]-1.0)*local_coords[0];
  shape_functions[1] = 0.5*(local_coords[0]+1.0)*local_coords[0];
  shape_functions[3] =     (local_coords[0]*local_coords[0]+1.0);
}


void ContactLineFaceQ3::Compute_Normal( VariableHandle CURRENT_POSITION,
					VariableHandle FACE_NORMAL )
{
  Real  shape_derivatives[2][3];
  Real  e1[3] = {0.0, 0.0, 0.0};
  Real  e2[3] = {0.0, 0.0, 0.0};
  Real  local_coords[2] = {0.5, 0.5};
  Real  a=0.0, b=0.0, c=0.0;
  Real* face_normal = Variable(FACE_NORMAL);
  
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<3; ++i) {
    Real* node_position = Node(i)->Variable(CURRENT_POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ      = 1.0/std::sqrt(a*b-c*c);
  face_normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])*detJ;
  face_normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])*detJ;
  face_normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])*detJ;

  acme::Normalize(face_normal);
}

void ContactLineFaceQ3::Compute_Normal( VariableHandle CURRENT_POSITION,
					Real* normal, Real* local_coords )
{
  Real  shape_derivatives[2][3];
  Real  e1[3] = {0.0, 0.0, 0.0};
  Real  e2[3] = {0.0, 0.0, 0.0};
  Real  a=0.0, b=0.0, c=0.0;
  
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<3; ++i) {
    Real* node_position = Node(i)->Variable(CURRENT_POSITION);
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ      = 1.0/std::sqrt(a*b-c*c);
  normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])*detJ;
  normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])*detJ;
  normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])*detJ;

  acme::Normalize(normal);
}

void ContactLineFaceQ3::Compute_Normal( Real** nodal_positions, 
                                        Real* local_coords, Real* normal)
{
  Real  shape_derivatives[2][3];
  Real  e1[3] = {0.0, 0.0, 0.0};
  Real  e2[3] = {0.0, 0.0, 0.0};
  Real  a=0.0, b=0.0, c=0.0;
  
  Compute_Shape_Derivatives( local_coords, shape_derivatives );
  for (int i=0; i<3; ++i) {
    Real* node_position = nodal_positions[i];
    e1[0] += shape_derivatives[0][i]*node_position[0];
    e1[1] += shape_derivatives[0][i]*node_position[1];
    e1[2] += shape_derivatives[0][i]*node_position[2];
    e2[0] += shape_derivatives[1][i]*node_position[0];
    e2[1] += shape_derivatives[1][i]*node_position[1];
    e2[2] += shape_derivatives[1][i]*node_position[2];
  }
  for (int i=0; i<3; ++i) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real detJ      = 1.0/std::sqrt(a*b-c*c);
  normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])*detJ;
  normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])*detJ;
  normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])*detJ;

  acme::Normalize(normal);
}

void ContactLineFaceQ3::Compute_CharacteristicLength( VariableHandle CURRENT_POSITION,
				                      VariableHandle CHARACTERISTIC_LENGTH )
{
  ContactNode<Real>* node0 = Node(0);
  ContactNode<Real>* node1 = Node(1);
  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  Real* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  *characteristiclength = std::sqrt((Position1[0]-Position0[0])*(Position1[0]-Position0[0]) +
                                    (Position1[1]-Position0[1])*(Position1[1]-Position0[1]) +
                                    (Position1[2]-Position0[2])*(Position1[2]-Position0[2]));
}

void ContactLineFaceQ3::Compute_Centroid( VariableHandle CURRENT_POSITION,
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

void ContactLineFaceQ3::Compute_Edge_Normal( VariableHandle, 
					     VariableHandle,
					     int, Real*)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "Error: Edges don't exist in 2D" << std::endl;
#endif
  POSTCONDITION( 0 );
}

void ContactLineFaceQ3::Compute_Local_Coordinates( Real Config_Param, 
						   VariableHandle POSITION0,
						   VariableHandle POSITION1, 
						   VariableHandle FACE_NORMAL,
						   Real* global_coords, 
						   Real* local_coords )
{
  Real node_positions[3][3];
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

void ContactLineFaceQ3::Compute_Local_Coordinates( VariableHandle POSITION, 
						   Real* global_coords, 
						   Real* local_coords )
{
  Real node_positions[3][3];
  for (int i=0; i<Nodes_Per_Face(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactLineFaceQ3::Compute_Global_Coordinates( VariableHandle POSITION,
						    Real* local_coords,
						    Real* global_coords )
{
  Real N[3];
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Evaluate_Shape_Functions( local_coords, N );
  for( int i=0 ; i<3 ; ++i){
    Real* node_position = Node(i)->Variable(POSITION);
    global_coords[0]   += N[i]*node_position[0];
    global_coords[1]   += N[i]*node_position[1];
    global_coords[2]   += N[i]*node_position[2];
  }
}

void ContactLineFaceQ3::Get_Close_Edges( Real* , int& , int&, int& )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactLineFaceQ3::Get_Close_Edges not yet implemented" << std::endl;
#endif
  POSTCONDITION( 0 );
}

bool ContactLineFaceQ3::Is_Inside_Face( Real* local_coords )
{
  PRECONDITION( 0 );
  return false;
}

ContactFace<Real>* ContactLineFaceQ3::Neighbor( Real* local_coords )
{
  PRECONDITION( 0 );
  return (ContactFace<Real>*) NULL;
}

void ContactLineFaceQ3::FacetDecomposition(int& nfacets, 
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
  nfacets = 5;
  Real* node0_position = Node(0)->Variable(POSITION0);
  Real* node1_position = Node(1)->Variable(POSITION0);
  Real  dx = (node1_position[0]-node0_position[0])/(Real)nfacets;
  Real  dy = (node1_position[1]-node0_position[1])/(Real)nfacets;
  Real  dz = (node1_position[2]-node0_position[2])/(Real)nfacets;
  for (int i=0; i<nfacets; ++i) {
    Real node0[3], node1[3];
    if (i==0) {
      node0[0] = node0_position[0];
      node0[1] = node0_position[1];
      node0[2] = node0_position[2];
    } else {
      node0[0] = node0_position[0]+i*dx;
      node0[1] = node0_position[1]+i*dy;
      node0[2] = node0_position[2]+i*dz;
    }
    if (i==nfacets-1) {
      node1[0] = node1_position[0];
      node1[1] = node1_position[1];
      node1[2] = node1_position[2];
    } else {
      node1[0] = node0_position[0]+(i+1)*dx;
      node1[1] = node0_position[1]+(i+1)*dy;
      node1[2] = node0_position[2]+(i+1)*dz;
    }
    coordinates0[0+6*i] = node0[0];
    coordinates0[1+6*i] = node0[1];
    coordinates0[2+6*i] = node0[2];
    coordinates0[3+6*i] = node1[0];
    coordinates0[4+6*i] = node1[1];
    coordinates0[5+6*i] = node1[2];
    normals0[0+3*i] = node1[1] - node0[1];
    normals0[1+3*i] = node0[0] - node1[0];
    normals0[2+3*i] = 0.0;
    Real Mag = std::sqrt( normals0[0+3*i]*normals0[0+3*i] +
		          normals0[1+3*i]*normals0[1+3*i] +
		          normals0[2+3*i]*normals0[2+3*i]);
    if( Mag ){
      Real invMag = 1.0 / Mag;
      normals0[0+3*i] *= invMag;
      normals0[1+3*i] *= invMag;
      normals0[2+3*i] *= invMag;
    }
  } 
  if (coordinates1!=NULL) {
    node0_position = Node(0)->Variable(POSITION1);
    node1_position = Node(1)->Variable(POSITION1);
    dx = (node1_position[0]-node0_position[0])/(Real)nfacets;
    dy = (node1_position[1]-node0_position[1])/(Real)nfacets;
    dz = (node1_position[2]-node0_position[2])/(Real)nfacets;
    for (int i=0; i<nfacets; ++i) {
      Real node0[3], node1[3];
      if (i==0) {
        node0[0] = node0_position[0];
        node0[1] = node0_position[1];
        node0[2] = node0_position[2];
      } else {
        node0[0] = node0_position[0]+i*dx;
        node0[1] = node0_position[1]+i*dy;
        node0[2] = node0_position[2]+i*dz;
      }
      if (i==nfacets-1) {
        node1[0] = node1_position[0];
        node1[1] = node1_position[1];
        node1[2] = node1_position[2];
      } else {
        node1[0] = node0_position[0]+(i+1)*dx;
        node1[1] = node0_position[1]+(i+1)*dy;
        node1[2] = node0_position[2]+(i+1)*dz;
      }
      coordinates0[0+6*i] = node0[0];
      coordinates0[1+6*i] = node0[1];
      coordinates0[2+6*i] = node0[2];
      coordinates0[3+6*i] = node1[0];
      coordinates0[4+6*i] = node1[1];
      coordinates0[5+6*i] = node1[2];
      normals0[0+3*i] = node1[1] - node0[1];
      normals0[1+3*i] = node0[0] - node1[0];
      normals0[2+3*i] = 0.0;
      Real Mag = std::sqrt( normals0[0+3*i]*normals0[0+3*i] +
                            normals0[1+3*i]*normals0[1+3*i] +
                            normals0[2+3*i]*normals0[2+3*i]);
      if( Mag ){
	Real invMag = 1.0 / Mag;
	normals0[0+3*i] *= invMag;
	normals0[1+3*i] *= invMag;
	normals0[2+3*i] *= invMag;
      }
    } 
  }
  if (coordinates2!=NULL) {
    node0_position = Node(0)->Variable(POSITION2);
    node1_position = Node(1)->Variable(POSITION2);
    dx = (node1_position[0]-node0_position[0])/(Real)nfacets;
    dy = (node1_position[1]-node0_position[1])/(Real)nfacets;
    dz = (node1_position[2]-node0_position[2])/(Real)nfacets;
    for (int i=0; i<nfacets; ++i) {
      Real node0[3], node1[3];
      if (i==0) {
        node0[0] = node0_position[0];
        node0[1] = node0_position[1];
        node0[2] = node0_position[2];
      } else {
        node0[0] = node0_position[0]+i*dx;
        node0[1] = node0_position[1]+i*dy;
        node0[2] = node0_position[2]+i*dz;
      }
      if (i==nfacets-1) {
        node1[0] = node1_position[0];
        node1[1] = node1_position[1];
        node1[2] = node1_position[2];
      } else {
        node1[0] = node0_position[0]+(i+1)*dx;
        node1[1] = node0_position[1]+(i+1)*dy;
        node1[2] = node0_position[2]+(i+1)*dz;
      }
      coordinates0[0+6*i] = node0[0];
      coordinates0[1+6*i] = node0[1];
      coordinates0[2+6*i] = node0[2];
      coordinates0[3+6*i] = node1[0];
      coordinates0[4+6*i] = node1[1];
      coordinates0[5+6*i] = node1[2];
      normals0[0+3*i] = node1[1] - node0[1];
      normals0[1+3*i] = node0[0] - node1[0];
      normals0[2+3*i] = 0.0;
      Real Mag = std::sqrt( normals0[0+3*i]*normals0[0+3*i] +
                            normals0[1+3*i]*normals0[1+3*i] +
                            normals0[2+3*i]*normals0[2+3*i]);
      if( Mag ){
	Real invMag = 1.0 / Mag;
        normals0[0+3*i] *= invMag;
        normals0[1+3*i] *= invMag;
        normals0[2+3*i] *= invMag;
      }
    }
  }
}

void ContactLineFaceQ3::FacetStaticRestriction(int nfacets, Real* coordinates, 
                                               Real* normals, Real* ctrcl_facets, 
                                               Real* ctrcl)
{
  int ii1=0,ilocc=0,ilocs=0;
  int iistored=0,iconcave=0,iinside=1,iout=2;
  Real projcv,projmv;
  Real dctds,xmc,ymc,zmc,xms,yms,zms,vecmix,vecmiy,vecmiz;
  // Contract the closest point projection of the Node-with-N_Lines
  // into a single contact with the bilinear LINE master surface.
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
        // Determine if the two line segments are concave or convex
        // relative to each other
        // centroid of candidate master surface line
        xmc = 0.5*(coordinates[0+6*ii] + 
                   coordinates[3+6*ii]);
        ymc = 0.5*(coordinates[1+6*ii] + 
                   coordinates[4+6*ii]);
        zmc = 0.5*(coordinates[2+6*ii] + 
                   coordinates[5+6*ii]);
        // centroid of stored master surface line
        xms = 0.5*(coordinates[0+6*iistored] + 
                   coordinates[3+6*iistored]);
        yms = 0.5*(coordinates[1+6*iistored] + 
                   coordinates[4+6*iistored]);
        zms = 0.5*(coordinates[2+6*iistored] + 
                   coordinates[5+6*iistored]);
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
          // INSIDE the candidate m.s. line
          //==============================
          if( std::abs(ilocs) == iinside ) {
            // INSIDE the stored m.s. line
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
            // OUTSIDE the stored m.s. line
            // -=> a) choose the INSIDE contact (i.e. replace stored m.s.)
            for(int k=0; k<LENGTH; ++k) {
              ctrcl_facets[ii1*LENGTH+k] = ctrcl_facets[ii*LENGTH+k];
	    }
            iistored = ii;
          }
        }  else if( std::abs(ilocc) == iout ) {
          // OUTSIDE the candidate m.s. line
          //====================================
	  if( std::abs(ilocs) == iinside ) {
            // INSIDE the stored m.s. triangle
            //   -=> a) choose the INSIDE contact (i.e. keep stored m.s.)
	  } else if( std::abs(ilocs) == iout ) {
            // OUTSIDE the stored m.s. line
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

void ContactLineFaceQ3::FacetDynamicRestriction(int nfacets,
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
  
  // Contract the moving search of the Node-with-N_Lines into a single
  // contact with a bilinear LINE master surface
  //
  // Notes:  The memory is laid out so that we don't need to look at the first
  //         line.  If its not in contact it will be ignored above.  If it
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
  for (int i=0; i<LENGTH; ++i) {
    ctrcl[i] = ctrcl_facets[i];
  }
}


void ContactLineFaceQ3::Smooth_Normal( VariableHandle, VariableHandle, 
				       VariableHandle, VariableHandle,
			ContactSearch::Smoothing_Resolution resolution,
				       Real, Real*, Real*,Real )
{
  POSTCONDITION( 0 );
}

int ContactLineFaceQ3::FaceEdge_Intersection(VariableHandle POSITION,
                                             ContactEdge<Real>* edge, Real* coords)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactLineFaceQ3::FaceEdge_Intersection not yet implemented\n";
#endif
  POSTCONDITION( 0 );
  return 0;
}
  
void 
ContactLineFaceQ3::Compute_Node_Areas(VariableHandle, VariableHandle, Real*)
{
}
/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/* The following functions are supplied for doing generic computations   */
/* on a quadratic Q8 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary Q8 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/


void ContactLineFaceQ3::Compute_Shape_Functions( Real local_coords[2],
						 Real shape_functions[3] )
{
  Real xi = local_coords[0];
  shape_functions[0] = -0.5*xi*(1.0-xi);
  shape_functions[1] =  0.5*xi*(1.0+xi);
  shape_functions[2] =  (1.0-xi)*(1.0-xi);
}

void ContactLineFaceQ3::Compute_Shape_Derivatives( Real local_coords[2],
						   Real shape_derivatives[1][3] )
{
  Real xi = local_coords[0];
  shape_derivatives[0][0] = -0.5*(1.0-2.0*xi);
  shape_derivatives[0][1] =  0.5*(1.0+2.0*xi);
  shape_derivatives[0][2] = -2.0*xi;
}

void ContactLineFaceQ3::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
					      Real global_coords[3],
					      Real local_coords[3] )
{
  int  i, j;
  int  nnodes=3;
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
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    break;
  case 2:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    break;
  }
  if (i<nnodes) return;
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
  Real shape_derivatives[1][3];
  while (!converged && iterations<max_iterations) {
    coords[0] = s0;
    coords[1] = t0;
    Compute_Global_Coords(node_positions , coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( coords, shape_derivatives );
    for (i=0; i<1; ++i) {
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
    std::cerr << "ContactLineFaceQ3::Compute_Local_Coords() did not converge" 
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

void ContactLineFaceQ3::Compute_Global_Coords( Real node_positions[3][3],
					       Real local_coords[2],
					       Real global_coords[3] )
{
  Real N[3];
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

void ContactLineFaceQ3::Interpolate_Scalar( Real  local_coords[2],
					    Real  node_scalars[3],
					    Real& interpolated_scalar )
{
  Real N[3];
  int  nnodes=3;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void ContactLineFaceQ3::Interpolate_Vector( Real local_coords[2],
					    Real node_vectors[3][3],
					    Real interpolated_vector[3] )
{
  Real N[3];
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

