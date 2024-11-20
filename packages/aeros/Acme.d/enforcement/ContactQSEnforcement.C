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


#include "ContactTDEnfModel.h"
#include "ContactErrors.h"
#include "ContactQSEnforcement.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "ContactShellNode.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "cstring"
#include <cmath>
#include <cstdio>
#include <unistd.h>

// Internal constants
#undef  NDIM
#define NDIM 3
#undef  ZERO_TOL
#define ZERO_TOL 1.0E-14
#undef  SMALL
#define SMALL 1.0e-6
#undef  CONVERGENCE_TOL
#define CONVERGENCE_TOL 1.0e-10

/* SYMBOLS ______________________________________________________
 * x = position
 * u = displacement
 * v = velocity
 * a = acceleration
 * f = force
 * m = mass
 * nm = outward unit normal vector ( f. nm < 0 for compression )
 *      (nm is typically outward for the master face and the
 *       force is typically stored with the slave node)
 * uj = jump in displacement i.e. relative displacement
 * vj = jump in velocity i.e. relative velocity
 * gap = uj . nm  , where  gap > 0 for separation
 * slip = vj * dt
 * _____________________________________________________________*/


ContactQSEnforcement::ContactQSEnforcement( const Real* enf_data, 
					    ContactSearch* Search,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::TDEnforcement,
			NSIZED, enf_data, true, FRICTION_MODEL_ID )
{  
  if (error != ContactSearch::NO_ERROR) return;
  num_iterations = 1;
  convergence_tolerance = CONVERGENCE_TOL;
  error = ContactSearch::NO_ERROR;
}

ContactQSEnforcement::ContactQSEnforcement( ContactSearch* Search,
					    const Real* restart_data,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::TDEnforcement,
			restart_data )
{  
  if (error != ContactSearch::NO_ERROR) return;
  num_iterations = 1;
  max_interactions = search->Max_Interactions_Per_Node(); 
  error = ContactSearch::NO_ERROR;
}

ContactQSEnforcement::~ContactQSEnforcement()
{
}

void ContactQSEnforcement::Remove_Gaps(Real* Position)
{
  //
  //  Stick slave nodes on to the master surfaces, just move the node from whatever point it is currently at and 
  //  put it directly on the alread computed contact intersection point
  //
  Set_Up();
  //
  // compute contact force from residual
  //
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    ContactNodeEntityInteraction* cnei  = cnei_group.Get_Interaction(0);
    ContactNode<Real>* node = cnei->Node();
    const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
    //
    //  Compute the position of the contact point
    //
    Real contact_point_position[3] = {0.0, 0.0, 0.0};
    for( int k=0 ; k < num_face_nodes ; ++k ){
      ContactNode<Real> *face_node = enforcement_node_list[cnei_group.Get_Face_Node_Index(0, k)];
      Real* face_node_pos = Position + NDIM*(face_node->HostGlobalArrayIndex());
      Real face_shape_function = cnei_group.Get_Face_Shape_Function(0,k);
      for(int j=0 ; j<NDIM ; ++j ){ contact_point_position[j] += face_shape_function*face_node_pos[j];}
    }
    //
    //  Explicity put the slave node at the contact point
    //  NKC, for now just trying to do tied contact, always put node back at point regardless of positive/negative gap...
    //
    Real* slave_node_pos = Position + NDIM*(node->HostGlobalArrayIndex());
    for(int j=0 ; j<NDIM ; ++j ){ slave_node_pos[j] = contact_point_position[j];}
  }
  Clean_Up(); 
}


void ContactQSEnforcement::Modify_Predictor_Velocity(Real* PredictedVelocity, const Real delta_time) 
{
  //
  //  Set up the contact arrays for use in enforcement
  //
  Set_Up();
  //
  //  Modify the slave node prediction model (VEL) to account for anticipated contact.
  //  NKC, assuming tied contact, thus set slave node moition to the face motion
  //
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    ContactNodeEntityInteraction* cnei  = cnei_group.Get_Interaction(0);
    ContactNode<Real>* node = cnei->Node();

    Real contact_point_velocity[3] = {0.0, 0.0, 0.0};

    const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
    for( int k=0 ; k < num_face_nodes ; ++k ){
      ContactNode<Real> *face_node = enforcement_node_list[cnei_group.Get_Face_Node_Index(0, k)];
      Real face_shape_function = cnei_group.Get_Face_Shape_Function(0,k);
      Real* face_node_vel = PredictedVelocity + NDIM*(face_node->HostGlobalArrayIndex());
      for(int j=0 ; j<NDIM ; ++j ){ contact_point_velocity[j] += face_shape_function*face_node_vel[j];}
    }

    Real* slave_node_velocity = PredictedVelocity + NDIM*(node->HostGlobalArrayIndex());
    for(int j=0 ; j<NDIM ; ++j ){ slave_node_velocity[j] = contact_point_velocity[j];}    
  }
  Clean_Up();
}

void ContactQSEnforcement::Modify_Preconditioner(Real* Preconditioner) 
{
  //
  //  Modify the preconditioner, account add stiffness of slave node to the master node.
  //  NKC, tied contact, directly add the stiffness based on the face shape functions
  //
  Set_Up();
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    ContactNodeEntityInteraction* cnei  = cnei_group.Get_Interaction(0);
    ContactNode<Real>* node = cnei->Node();
    Real* slave_node_preconditioner = Preconditioner + NDIM*(node->HostGlobalArrayIndex());

    const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
    for( int k=0 ; k < num_face_nodes ; ++k ){
      ContactNode<Real> *face_node = enforcement_node_list[cnei_group.Get_Face_Node_Index(0, k)];
      Real face_shape_function = cnei_group.Get_Face_Shape_Function(0,k);
      Real* face_node_preconditioner = Preconditioner + NDIM*(face_node->HostGlobalArrayIndex());
      for(int j=0 ; j<NDIM ; ++j ){ face_node_preconditioner[j] += face_shape_function*slave_node_preconditioner[j];}
    }
  }
  Clean_Up();
}

void ContactQSEnforcement::Compute_Contact_Force(Real* Residual) {
  Set_Up();
  //
  // compute contact force from residual
  //
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    ContactNodeEntityInteraction* cnei  = cnei_group.Get_Interaction(0);
    ContactNode<Real>* node = cnei->Node();
    Real * slave_node_residual = Residual + NDIM*(node->HostGlobalArrayIndex());
    const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
    for( int k=0 ; k < num_face_nodes ; ++k ){
      ContactNode<Real> *face_node = enforcement_node_list[cnei_group.Get_Face_Node_Index(0, k)];
      Real face_shape_function = cnei_group.Get_Face_Shape_Function(0,k);
      Real* face_node_residual = Residual + NDIM*(face_node->HostGlobalArrayIndex());
      for(int j=0 ; j<NDIM ; ++j ){ face_node_residual[j] -= face_shape_function*slave_node_residual[j];}
    }
    //
    // slave node is assumed to be in equilibrium, zero out the residual
    //
    for(int j=0 ; j<NDIM ; ++j ){ slave_node_residual[j] = 0.0;}
  }
  Clean_Up();
}

void ContactQSEnforcement::Update_Search_Direction(Real* SearchDirection)
{
  Set_Up();
  //
  //  NKC, again tied contact.
  //  modify the slave node search direction for the face nodes
  //
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    ContactNodeEntityInteraction* cnei  = cnei_group.Get_Interaction(0);
    ContactNode<Real>* node = cnei->Node();
    Real* slave_node_search_direction = SearchDirection + NDIM*(node->HostGlobalArrayIndex());
    const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
    for( int k=0 ; k < num_face_nodes ; ++k ){
      ContactNode<Real> *face_node = enforcement_node_list[cnei_group.Get_Face_Node_Index(0, k)];
      Real face_shape_function = cnei_group.Get_Face_Shape_Function(0,k);
      Real* face_node_search_direction = SearchDirection + NDIM*(face_node->HostGlobalArrayIndex());
      for(int j=0 ; j<NDIM ; ++j ){ slave_node_search_direction[j] -= face_shape_function*face_node_search_direction[j];}
    }
  }
  Clean_Up();
}
