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


#ifndef ContactEnfModel_h_
#define ContactEnfModel_h_

#include "Contact_Defines.h"
#include "ContactEnforcement.h"
#include "ContactSearch.h"
#include <cstring>

#undef OMIT_NS
#define OMIT_NS 1 

class ContactTable;

class ContactEnfModel {

 public:
  ContactEnfModel( int ID, ContactEnforcement::Enforcement_Model_Types Type,
		   ContactTopology* Topology );
  ContactEnfModel(ContactEnforcement::Enforcement_Model_Types Type,
		  ContactTopology* Topology );
  virtual ~ContactEnfModel();

  int ID() { return id; };
  ContactEnforcement::Enforcement_Model_Types Type() {return type;};
  ContactSearch::ContactErrorCode Initialize_State_Data();
  void Update_State();
  void Match_Interactions();

  virtual int Num_Interaction_State_Variables() = 0;
  virtual int Num_Node_State_Variables() = 0;
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables ) = 0;
  // virtual void Initialize_Interaction_State_Data( Real* data ) = 0;
  virtual void Initialize_Node_State_Data( Real* data ) = 0;
  virtual int Restart_Size() = 0;
  virtual int Extract_Restart_Data( Real* restart_data ) = 0;
  virtual int Implant_Restart_Data( const Real* restart_data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, const Real* data ) = 0;

  ContactSearch::ContactErrorCode
    Update_For_Topology_Change( int number_of_new_nodes, int* old_to_new_map );

  ContactErrors* errors;
  char message[81];

 protected:

  ContactTopology* topology;
  int id;
  ContactEnforcement::Enforcement_Model_Types type;
  int number_of_nodes;
  Real* physical_face_normals_0;
  Real* physical_face_normals_1;
  Real* node_state_data_0;
  Real* node_state_data_1;
  Real* interaction_state_data_0;
  Real* interaction_state_data_1;
  Real* Node_State_Data( int node_index, int state ){
    int offset = node_index*Num_Node_State_Variables();
    if( state == 0 )
      return node_state_data_0+offset;
    else
      return node_state_data_1+offset;
  };
  int Restart_Size_State_Data()
    { return number_of_nodes*Num_Node_State_Variables(); };
  int Extract_Restart_State_Data( Real* buf )
    { 
      std::memcpy( buf, node_state_data_0, 
	      number_of_nodes*Num_Node_State_Variables()*sizeof(Real) ); 
      return number_of_nodes*Num_Node_State_Variables();
    };
  int Implant_Restart_State_Data( const Real* buf )
    { 
      node_state_data_0 = new Real[number_of_nodes*Num_Node_State_Variables()];
      std::memcpy( node_state_data_0, buf,
	      number_of_nodes*Num_Node_State_Variables()*sizeof(Real) ); 
      return number_of_nodes*Num_Node_State_Variables();
    };
  /*
  Real* NFI_State_Data( int node_index, int physical_face_index,
			int state ){
    int offset = node_index*3*Num_Interaction_State_Variables() + 
                 physical_face_index*Num_Interaction_State_Variables();
    if( state == 0 )
      return interaction_state_data_0+offset;
    else
      return interaction_state_data_1+offset;
  };
  */
  
 private:
  
  ContactEnfModel(ContactEnfModel&);
  ContactEnfModel& operator=(ContactEnfModel&);
  
};

#endif
