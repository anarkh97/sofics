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


#include "ContactEnfModel.h"

ContactEnfModel::
ContactEnfModel( int ID, ContactEnforcement::Enforcement_Model_Types Type,
		 ContactTopology* Topology ) 
  : errors(nullptr), topology( Topology ), id(ID), type(Type), 
    number_of_nodes(-1), 
    physical_face_normals_0(nullptr), physical_face_normals_1(nullptr),
    node_state_data_0( nullptr ), node_state_data_1(nullptr),
    interaction_state_data_0(nullptr), interaction_state_data_1( nullptr )
{
  std::memset(message,0,81);
}

ContactEnfModel::
ContactEnfModel( ContactEnforcement::Enforcement_Model_Types Type,
		 ContactTopology* Topology ) 
  : errors(nullptr), topology( Topology ), id(-1), type(Type), 
    number_of_nodes(Topology->Number_of_Nodes()), 
    physical_face_normals_0( nullptr ), physical_face_normals_1( nullptr ),
    node_state_data_0( nullptr ), node_state_data_1( nullptr ),
    interaction_state_data_0( nullptr ), interaction_state_data_1( nullptr )
{
  std::memset(message,0,81);
}

ContactEnfModel::~ContactEnfModel()
{
  delete [] physical_face_normals_0;
  delete [] physical_face_normals_1;
  delete [] node_state_data_0;
  delete [] node_state_data_1;
  delete [] interaction_state_data_0;
  delete [] interaction_state_data_1;
}

ContactSearch::ContactErrorCode
ContactEnfModel::Initialize_State_Data()
{
  int state = 0;
  number_of_nodes = topology->Number_of_Nodes();
  if( Num_Node_State_Variables() && number_of_nodes ){
    node_state_data_0 = new Real[number_of_nodes*Num_Node_State_Variables()];
    node_state_data_1 = new Real[number_of_nodes*Num_Node_State_Variables()];
    for( int i=0 ; i<number_of_nodes ; ++i ){
      Initialize_Node_State_Data( Node_State_Data( i,state ) );
    }
  }
  /*
  if( Num_Interaction_State_Variables() && number_of_nodes ){
    interaction_state_data_0 = new 
      Real[number_of_nodes*Num_Interaction_State_Variables()];
    interaction_state_data_1 = new 
      Real[number_of_nodes*Num_Interaction_State_Variables()];
    for( int i=0 ; i<number_of_nodes ; i++ ){
      Initialize_Interaction_State_Data( NFI_State_Data( i,0,state ) );
      Initialize_Interaction_State_Data( NFI_State_Data( i,1,state ) );
      Initialize_Interaction_State_Data( NFI_State_Data( i,2,state ) );
    }
  }
  */
  if( Num_Interaction_State_Variables() ){
    return ContactSearch::INTERNAL_ERROR;
  }

  /*
  if( number_of_nodes ){
    physical_face_normals_0 = new Real[number_of_nodes*3*3];
    std::memset( physical_face_normals_0, 0, number_of_nodes*3*3*sizeof(Real) );
    physical_face_normals_1 = new Real[number_of_nodes*3*3];
  }
  */

  return ContactSearch::NO_ERROR;
}

void ContactEnfModel::Update_State()
{
  /*
  Real* tmp = physical_face_normals_0;
  physical_face_normals_0 = physical_face_normals_1;
  physical_face_normals_1 = physical_face_normals_0;
  */
  Real* tmp = node_state_data_0;
  node_state_data_0 = node_state_data_1;
  node_state_data_1 = tmp;
  /*
  tmp = interaction_state_data_0;
  interaction_state_data_0 = interaction_state_data_1;
  interaction_state_data_1 = tmp;
  */
}

void ContactEnfModel::Match_Interactions()
{
  // This function will rearrange, the old interaction information to
  // be consistent with the ordering this time step by using the physical
  // face normals.
  /*
  for( i=0 ; i<number_node_face_interactions ; i++ ){
    ContactNodeEntityInteraction* cnfi = nfi_interactions[i];
    int node_index = cnfi->Node()->Contact_Index();
    Real* new_pfnorm = ;
  }
  */
}


ContactSearch::ContactErrorCode
ContactEnfModel::Update_For_Topology_Change( int number_of_new_nodes,
					     int* old_to_new_map )
{
  int i,j;

  Real* old_ns_data_0 = node_state_data_0;
  Real* old_ns_data_1 = node_state_data_1;
  
  int num_node_state_var = Num_Node_State_Variables();
  if( num_node_state_var && number_of_new_nodes ){
    // Allocate new node state data and initialize as new 
    node_state_data_0 = new Real[number_of_nodes*Num_Node_State_Variables()];
    node_state_data_1 = new Real[number_of_nodes*Num_Node_State_Variables()];
    for( i=0 ; i<number_of_nodes ; ++i ){
      Initialize_Node_State_Data( Node_State_Data( i,0 ) );
      Initialize_Node_State_Data( Node_State_Data( i,1 ) );
    }
    for( i=0 ; i<number_of_nodes ; ++i ){
      int new_index = old_to_new_map[i];
      if( new_index >= 0 ){
	int old_offset = i*num_node_state_var;
	int new_offset = new_index*num_node_state_var;
	Real* old_state_data_0 = old_ns_data_0 + old_offset;
	Real* old_state_data_1 = old_ns_data_1 + old_offset;
	Real* new_state_data_0 = node_state_data_0 + new_offset;
	Real* new_state_data_1 = node_state_data_1 + new_offset;
	for( j=0 ; j<num_node_state_var ; ++j ){
	  new_state_data_0[j] = old_state_data_0[j];
	  new_state_data_1[j] = old_state_data_1[j];
	}
      }
    }
    delete [] old_ns_data_0;
    delete [] old_ns_data_1;
  } 
  number_of_nodes = number_of_new_nodes;
  return ContactSearch::NO_ERROR;
}
