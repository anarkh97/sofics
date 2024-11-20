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


#include "ContactErrors.h"
#include "ContactTopology.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactTiedKinematics.h"
#include <cstring>
#include <cstdio>
#include <cmath>

ContactTiedKinematics::ContactTiedKinematics( const Real* enf_data,
					      ContactSearch* Search,
			       ContactSearch::ContactErrorCode& error )
  : ContactEnforcement( error, Search, ContactEnforcement::TiedKinematics,
			NSIZED, enf_data, false )
{
  if (error != ContactSearch::NO_ERROR) return;
  // plot vars
  number_global_plot_vars = 0;
  number_element_plot_vars = 0;
  number_nodal_plot_vars = 3;
#ifndef CONTACT_NO_EXODUS_OUTPUT
  final_position = new Real[number_of_nodes*3];
#else
  final_position = NULL;
#endif
}

ContactTiedKinematics::ContactTiedKinematics( ContactSearch* Search,
					      const Real* restart_data,
                               ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::TiedKinematics,
                        restart_data )
{
  if (error != ContactSearch::NO_ERROR) return;
  number_global_plot_vars = 0;
  number_nodal_plot_vars = 3;
  number_element_plot_vars = 0;

#ifndef CONTACT_NO_EXODUS_OUTPUT
  final_position = new Real[number_of_nodes*3];
#else
  final_position = NULL;
#endif
  error = ContactSearch::NO_ERROR;
}

ContactTiedKinematics::~ContactTiedKinematics()
{
  if( final_position ) delete [] final_position;
}

void ContactTiedKinematics::Get_Nodal_Plot_Variable( int var_num, Real* data )
{
#ifndef CONTACT_NO_EXODUS_OUTPUT
  PRECONDITION( var_num < 3 );
  std::memcpy( data, final_position+var_num*number_of_nodes,
          number_of_nodes*sizeof(Real) );
#endif
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Extract_General_Restart_Variable( Real* data )
{
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Implant_General_Restart_Variable( Real* data )
{
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Implant_Nodal_Restart_Variable( int n, Real* data )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Extract_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Implant_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Extract_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Implant_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Extract_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Implant_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTiedKinematics::Compute_Position( Real* position )
{
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;
#ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(communicator)==0) {
    std::cout << "  Tied Kinematics\n"<<std::flush;
  }
#endif

  VariableHandle CURRENT_POSITION =
      topology->Variable_Handle( ContactTopology::Current_Position );

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
  for( int i=0 ; i<Number_Debug_Nodes() ; ++i){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      Real* ipos = dn->Variable(CURRENT_POSITION);
      postream << "TiedKinematics Initial Position for Node ("
	       << dn->Exodus_ID() << ") = " 
	       << ipos[0] << " " << ipos[1] << " " << ipos[2] << "\n";
    }
  }
  postream.flush();
#endif

  // Call base class Set_Up
  error_code = Set_Up();
  if( error_code ) return error_code;

  // Compute the final position, initially the current position
  for( int i=0 ; i<number_of_nodes ; ++i){
    Real* pos = position+3*i;
    if( enforcement_node_list[i]->Ownership() == ContactTopologyEntity<Real>::OWNED ){
      Real* ipos = enforcement_node_list[i]->Variable(CURRENT_POSITION);
      pos[0] = ipos[0];
      pos[1] = ipos[1];
      pos[2] = ipos[2];
    } else {
      pos[0] = 0.0;
      pos[1] = 0.0;
      pos[2] = 0.0;
    }
  }

  // Now compute the new positions for the slave nodes
  for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
       cnfi != NULL;
       cnfi = node_entity_list.Node_Face_Iterator_Next()) {
    int contact_index = cnfi->Node()->ProcArrayIndex();
    Real* pos = position+3*contact_index;
    cnfi->Face()->Compute_Global_Coordinates( CURRENT_POSITION,
	      cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES),
					      pos );
  }

  // Now we must swapadd the position array to get the new coordinates to
  // all processors (not just the owning processor).
  swapadd_data_array( position, 3 );

#ifndef CONTACT_NO_EXODUS_OUTPUT
  for( int i=0 ; i<number_of_nodes; ++i){
    for( int j=0 ; j<dimensionality ; ++j)
      final_position[j*number_of_nodes+i] = position[i*dimensionality+j];
  }  
#endif

  // map back to host-code ordering
  Real* tmp = new Real[3*number_of_nodes];
  std::memcpy(tmp, position, 3*number_of_nodes*sizeof(Real));
  for( int i=0 ; i<number_of_nodes; ++i){
    int index  = enforcement_node_list[i]->HostGlobalArrayIndex();
    Real* pos  = position+3*index;
    Real* ipos = tmp+3*i;
    pos[0] = ipos[0];
    pos[1] = ipos[1];
    pos[2] = ipos[2];
  } 
  delete [] tmp;

  // Call the base class clean up functions
  Clean_Up();

#ifdef CONTACT_DEBUG_NODE
  for( int i=0 ; i<Number_Debug_Nodes() ; ++i){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      int index = dn->ProcArrayIndex();
      postream << "TiedKinematics  Final  Position for Node ("
	       << dn->Exodus_ID() << ") = " 
	       << position[index*dimensionality+0] << " " 
	       << position[index*dimensionality+1] << " " 
	       << position[index*dimensionality+2] << "\n";
    }
  }
  postream.flush();
#endif

  return (ContactSearch::ContactErrorCode)
    contact_global_error_check( error_code, communicator );
}


ContactSearch::ContactErrorCode
ContactTiedKinematics::Update_For_Topology_Change( int new_number_of_nodes,
						   int* old_to_new_map )
{
  // We don't need to map final_position to the new topology
  if( final_position ) delete [] final_position;
  if( new_number_of_nodes ) 
    final_position = new Real[3*new_number_of_nodes];
  else
    final_position = NULL;
  number_of_nodes = new_number_of_nodes;
  return ContactSearch::NO_ERROR;
}
