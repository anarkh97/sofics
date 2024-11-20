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


#include "ContactTopologyEntityList.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactTopologyEntityHash.h"
#include "ContactTopology.h"
#include <cstring>
#include <new>

ContactNodeEntityInteraction::ContactNodeEntityInteraction(bool is_tied_,
                                                           bool is_infSlip_,
                                                           bool is_glued_,
                                                           bool is_tracked_,
                                                           ContactNode<Real>* node_,
                                                           ContactTopologyEntity<Real>* entity_,
                                                           int Location,
                                                           InteractionType type_,
                                                           ContactType base_type_) :
  ContactInteractionEntity<Real>((Real*) &DataArray, base_type_),
  is_tied(is_tied_),
  is_infSlip(is_infSlip_),
  is_glued(is_glued_),
  is_tracked(is_tracked_),
  location(Location),
  node(node_),
  entity(entity_),
  type(type_)
{
  PRECONDITION(!(is_tied&&is_infSlip));
  PRECONDITION(!((is_tied||is_infSlip)&&is_glued));
  std::memset(DataArray, 0, (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real)); 
}

ContactNodeEntityInteraction::~ContactNodeEntityInteraction() {}

void  ContactNodeEntityInteraction::Modify_for_Kinematic_Constraints(
					     VariableHandle NUM_KIN_CONSTR,
					     VariableHandle KIN_CONSTR_VECTOR,
					     bool& Would_Violate_Constraints )
{
  // The return value denotes if the interaction would violate the kinematic
  // constraints.  If the node has no degrees of freedom or the constrained
  // direction is parallel to the interaction direction, this function returns
  // true.  If the above conditions are not true, the interaction is modified
  // (if needed) to be consistent with the kinematic constraints.

  const Real tol_1 = 0.996  ; // tolerance for 1 kin. constraint, std::cos(5 deg)
  const Real tol_2 = 0.087  ; // tolerance for 2 kin. constraint, std::sin(5 deg)

  int num_constraints = (int) *node->Variable( NUM_KIN_CONSTR );
  PRECONDITION( num_constraints >= 0 && num_constraints<= 3);
  Real * constraint_vector = node->Variable( KIN_CONSTR_VECTOR );
  Real * pushback_dir = Vector_Var(PUSHBACK_DIR);
  Real& gap = Scalar_Var(GAP_CUR);
  Real dot, mag;
  switch (num_constraints) {
  case 0 :
    // no kinematic constraints so there is nothing to do
    Would_Violate_Constraints = false;
    break;
  case 1:
    // for this case, the constraint_vector is the constrained direction
    dot = pushback_dir[0]*constraint_vector[0] +
          pushback_dir[1]*constraint_vector[1] +
          pushback_dir[2]*constraint_vector[2] ;
    if( std::fabs(dot) < tol_1 ){
      pushback_dir[0] -= dot*constraint_vector[0];
      pushback_dir[1] -= dot*constraint_vector[1];
      pushback_dir[2] -= dot*constraint_vector[2];
      mag = std::sqrt( pushback_dir[0]*pushback_dir[0] +
		       pushback_dir[1]*pushback_dir[1] +
		       pushback_dir[2]*pushback_dir[2] );
      gap *= mag;
      mag  = 1.0/mag;
      pushback_dir[0] *= mag;
      pushback_dir[1] *= mag;
      pushback_dir[2] *= mag;
      Would_Violate_Constraints = false;
    } else {
      Would_Violate_Constraints = true;
    }
    break;
  case 2:
    // for this case, the constraint_vector is the constrained direction
    dot = pushback_dir[0]*constraint_vector[0] +
          pushback_dir[1]*constraint_vector[1] +
          pushback_dir[2]*constraint_vector[2] ;
    if( std::fabs(dot) > tol_2 ){
      mag = dot < 0.0 ? -1.0 : 1.0;
      pushback_dir[0] = mag*constraint_vector[0];
      pushback_dir[1] = mag*constraint_vector[1];
      pushback_dir[2] = mag*constraint_vector[2];
      gap /= std::fabs(dot);
      Would_Violate_Constraints = false;
    } else {
      Would_Violate_Constraints = true;
    }
    break;
  case 3:
    // All DOF are specified so the interaction can not be enforced
    Would_Violate_Constraints = true;
    break;
  }
}

void ContactNodeEntityInteraction::Connect_Node( ContactTopologyEntityList& hash_table )
{   
  node = static_cast<ContactNode<Real> *>(hash_table.Find( &node_entity_data ));
  POSTCONDITION( node );
  Set_NodeEntityData();
}

void ContactNodeEntityInteraction::Connect_Node( ContactTopologyEntityHash& hash_table )
{
  ContactHostGlobalID gid( node_entity_data.host_gid[0], 
                           node_entity_data.host_gid[1] );
  node = static_cast<ContactNode<Real> *>(hash_table.find( gid ));
  POSTCONDITION( node );
  Set_NodeEntityData();
}

void ContactNodeEntityInteraction::Connect_Node( ContactTopology* topology )
{ 
  int block = node_entity_data.block_id;
  //int index = node_entity_data.index_in_block;
  //node = static_cast<ContactNode<Real> *>(topology->Node_Block(block)->NodeList()->Find( index ));
  node = static_cast<ContactNode<Real> *>
           (topology->Node_Block(block)->NodeList()->Find( &node_entity_data ));
  POSTCONDITION( node );
  Set_NodeEntityData();
}

void ContactNodeEntityInteraction::Connect_Node( ContactNode<Real>* Node )
{   
  node = Node;
  POSTCONDITION( node );
  Set_NodeEntityData();
}

CString ContactNodeEntityInteraction::Source_Name()
{
  switch( (ContactNodeEntityInteraction::InteractionSource) Scalar_Var(SOURCE) ){
  case( CLOSEST_POINT_PROJECTION_1 ):
    return CString( "Closest Point Projection, 1 configuration" );
  case( CLOSEST_POINT_PROJECTION_2 ):
    return CString( "Closest Point Projection, 2 configuration" );
  case( CLOSEST_POINT_PROJECTION_2_FROM_MOVING ):
    return CString( "Closest Point Projection, 2 configuration (From Moving)" );
  case( MOVING_INTERSECTION ):
    return CString( "Moving Intersection" );
  case( RETRIEVED_TIED ):
    return CString( "Retrieved Tied" );
  case( RETRIEVED_GLUED ):
    return CString( "Retrieved Glued" );
  case( UNKNOWN_SOURCE ):
    return CString( "Unknown" );
  default:
    POSTCONDITION(0);
    break;
  }
  POSTCONDITION(0);
  return CString( "Source Error" );
}
