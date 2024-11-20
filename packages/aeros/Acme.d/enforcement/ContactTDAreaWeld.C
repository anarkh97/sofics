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


// this is a near copy of ContactTDSpotWeld

#include "ContactTDAreaWeld.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactErrors.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "contact_assert.h"
#include <cstdio>
#undef LOCAL_PRINT_FLAG
#define LOCAL_PRINT_FLAG 0
#if (CONTACT_DEBUG_PRINT_LEVEL>=2) || (LOCAL_PRINT_FLAG>0)
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif

ContactTDAreaWeld::ContactTDAreaWeld( int ID, int* int_data, 
					Real* real_data, 
					ContactTopology* topo )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_SPOT_WELD, topo )
{
  C_n = real_data[0]; // normal_capacity
  C_t = real_data[1]; // tangential_capacity
  S_i = int_data[0];  // failure_steps  
  failed_model_id = int_data[1];
  failed_model = NULL;
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding area weld enforcement model, ID ="<<ID<<std::endl;
    std::cout<<"  failure model id = "<<failed_model_id<<std::endl;
    std::cout<<"  C_n = "<<C_n<<std::endl;
    std::cout<<"  C_t = "<<C_t<<std::endl;
    std::cout<<"  S_i = "<<S_i<<std::endl;
    std::cout<<std::flush;
  }
#endif
}


ContactTDAreaWeld::ContactTDAreaWeld( ContactTopology* topo )
  : ContactTDEnfModel(ContactEnforcement::TD_SPOT_WELD, topo)
{
  C_n = 0.0;
  C_t = 0.0;
  S_i = 0;
  failed_model = NULL;
}

ContactTDAreaWeld::~ContactTDAreaWeld()
{
}

ContactSearch::ContactErrorCode
ContactTDAreaWeld::Initialize_Model( int num_models, ContactEnfModel** models,
                                      int /* num_tables */,
                                      ContactTable** /* tables */)
{
  for( int i=0 ; i<num_models ; ++i){
    if( failed_model_id == models[i]->ID() ){
      failed_model = (ContactTDEnfModel*) models[i];
      return ContactSearch::NO_ERROR;
    }
  }
  return ContactSearch::INVALID_ID;
}


void ContactTDAreaWeld::Initialize_for_Time_Step()
{
  for( int i=0 ; i<number_of_nodes ; ++i)
    Node_State_Data(i,0)[DECREASED_STRENGTH_THIS_STEP] = 0.0;
}


int ContactTDAreaWeld::Num_Node_State_Variables()
{
  // Comment out the failed mode because of restart for now.  We know its 0
  return SIZE_NODE_STATE_DATA;//+failed_model->Num_Node_State_Variables();
}

int ContactTDAreaWeld::Num_Interaction_State_Variables()
{
  // Comment out the failed mode because of restart for now.  We know its 0
  return 0;//+failed_model->Num_Interaction_State_Variables();
}

int ContactTDAreaWeld::Interaction_Type( ContactNodeEntityInteraction* cnfi )
{
  if( Node_State_Data( cnfi->Node()->fcs_index, 0 )[STATUS] == 0 )
    return TDEM_TIED;
  else
    return failed_model->Interaction_Type( cnfi );
}

int ContactTDAreaWeld::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  words_added += Extract_Restart_State_Data( restart_data );
  restart_data[words_added++] = id;
  restart_data[words_added++] = C_n;
  restart_data[words_added++] = C_t;
  restart_data[words_added++] = S_i;
  restart_data[words_added++] = failed_model_id;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDAreaWeld::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  words_read += Implant_Restart_State_Data( restart_data );
  id = (int) restart_data[words_read++];
  C_n = restart_data[words_read++];
  C_t = restart_data[words_read++];
  S_i = (int) restart_data[words_read++];
  failed_model_id = (int) restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDAreaWeld::Extract_Nodal_Restart_Variable(int n, Real* data )
{
  if ( topology->Have_Shells()) {
    ContactShellHandler * shell_handler = topology->Shell_Handler();
    
    if( n>Num_Node_State_Variables()*
	shell_handler->Max_Num_Acme_Nodes_for_Host_Node() 
	|| n<=0 ){
      std::sprintf( message, "Variable index %d is unreasonable", n );
      errors->Add_Error_Message(message);
      return ContactSearch::INVALID_ID;
    }
    
    int host_to_acme_node = 0;
    int var_num = n;
    while (var_num > Num_Node_State_Variables()){
      var_num -= Num_Node_State_Variables();
      host_to_acme_node++;
    }
    
    for ( int i = 0; i < shell_handler->Number_Host_Code_Nodes(); ++i) {
      if (host_to_acme_node < shell_handler->Num_Acme_Nodes_for_Host_Node(i)){
	int node_idx = 
	  shell_handler->Acme_Node_for_Host_Node(i,host_to_acme_node);
	data[i] = 
	  node_state_data_0[node_idx*Num_Node_State_Variables()+var_num-1];
      } else {
        data[i] = 0.0;
      }
    }

  } else {
    
    if( n>Num_Node_State_Variables() || n<=0 ){
      std::sprintf( message, "Variable index %d is unreasonable", n );
      errors->Add_Error_Message(message);
      return ContactSearch::INVALID_ID;
    }
    
    for ( int i = 0; i < number_of_nodes; ++i)
      data[i] = node_state_data_0[i*Num_Node_State_Variables()+n-1];
  }

  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactTDAreaWeld::Implant_Nodal_Restart_Variable(int n, const Real* data )
{
  if ( topology->Have_Shells()) {
    ContactShellHandler * shell_handler = topology->Shell_Handler();

    if( n>Num_Node_State_Variables()*
	shell_handler->Max_Num_Acme_Nodes_for_Host_Node() 
	|| n<=0 ){
      std::sprintf( message, "Variable index %d is unreasonable", n );
      errors->Add_Error_Message(message);
      return ContactSearch::INVALID_ID;
    }
    
    int host_to_acme_node = 0;
    int var_num = n;
    while (var_num > Num_Node_State_Variables()){
      var_num -= Num_Node_State_Variables();
      host_to_acme_node++;
    }
    
    for ( int i = 0; i < shell_handler->Number_Host_Code_Nodes(); ++i) {
      if (host_to_acme_node < shell_handler->Num_Acme_Nodes_for_Host_Node(i)){
	int node_idx = 
	  shell_handler->Acme_Node_for_Host_Node(i,host_to_acme_node);
	node_state_data_0[node_idx*Num_Node_State_Variables()+var_num-1] =  
	  data[i];
      }
      else continue;
    }
    
  } else {

    if( n>Num_Node_State_Variables() || n<=0 ){
      std::sprintf( message, "Variable index %d is unreasonable", n );
      errors->Add_Error_Message(message);
      return ContactSearch::INVALID_ID;
    }
    
    for ( int i = 0; i < number_of_nodes; ++i)
      node_state_data_0[i*Num_Node_State_Variables()+n-1] =  data[i];

  }
  return ContactSearch::NO_ERROR;
}


bool ContactTDAreaWeld::Active_Interaction
(ContactNodeEntityInteraction* cnfi, Real gap)
{   
  Real& status = Node_State_Data( cnfi->Node()->fcs_index, 0 )[STATUS];
  if( status == 0 ){
    return true;
  } else {
    return failed_model->Active_Interaction( cnfi,gap);
  }
} 
  
bool ContactTDAreaWeld::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{
#if LOCAL_PRINT_FLAG>=2
  Real force_n = Dot(force,normal);
  Real force_t[3];
  Remove_Component(force,normal,force_t);
  Real s_n = Dot(slip,normal);
  Real s_t[3];
  Remove_Component(slip,normal,s_t);
  std::cout << " gap " << gap << ", s_n " << s_n 
	  <<  " inc_s_n " << Dot(rel_disp,normal) << "\n";
  std::cout << " f_n " << force_n << "," << " s_n " << s_n << "\n";
#endif


  Scale(force,1/area,force);
  Real& status = Node_State_Data( cnfi->Node()->fcs_index, 0 )[STATUS];
  if( status == 0 ){
    Real& strength=Node_State_Data( cnfi->Node()->fcs_index,0)[STRENGTH];
    Real& decremented = Node_State_Data( cnfi->Node()->fcs_index,0)
      [DECREASED_STRENGTH_THIS_STEP];
    Real f_n = Dot(force,normal);
    // The weld will support any amount of normal compressive load
    // Note: force*normal < 0.0 implies compressive
    Real norm_ratio = std::max(-f_n,0.0)/C_n; // I think this is wrong
    Real f_t[3];
    Remove_Component(force,normal,f_t);
    Real mag_f_t = Magnitude(f_t);
    Real tang_ratio = mag_f_t/C_t;
    Real trial = norm_ratio*norm_ratio + tang_ratio*tang_ratio;
#if LOCAL_PRINT_FLAG>=2
    std::cout << " strength: " << strength << ", n_ratio: " << norm_ratio
	     << ", t_ratio: " << tang_ratio << ", trial: " << trial <<"\n";
#endif
    if(strength == S_i) { // still welded
      if ( trial > 1.0) { // beginning failure
	if( decremented == 0.0 ){
	  strength -= 1.0;
	  decremented = 1.0;
	}
	Real scale = (strength+1)/S_i; // NOTE should be : (strength+1)/S_i;
#if LOCAL_PRINT_FLAG>=2
    std::cout << " scale: " << scale << "f_n: " << f_n << "\n";
#endif
	if( f_n < -scale*C_n ) f_n = -scale*C_n;
#if LOCAL_PRINT_FLAG>=2
    std::cout << " f_n+: " << f_n << "\n";
#endif
	Linear_Combination(f_n,normal,scale,f_t,force); 
      }
    }
    else { // failing
      if( decremented == 0.0 ){
	strength -= 1.0;
	decremented = 1.0;
      }
      Real scale = (strength+1)/S_i;
      if( f_n < -scale*C_n ) f_n = -scale*C_n;
#if LOCAL_PRINT_FLAG>=2
    std::cout << " f_n+: " << f_n << " scale:"<< scale <<" C_n "<< C_n << "\n";
#endif
      Linear_Combination(f_n,normal,scale,f_t,force); 
      if (strength <= 0.0) {
	Node_State_Data( cnfi->Node()->fcs_index, 0 )[STATUS] = 1.0;
	if (Needs_Glued_Search(failed_model->Interaction_Type(cnfi))) {
          cnfi->Is_Glued( true );
        } else {
          cnfi->Is_Glued( false );
        }
      }
    }
#if LOCAL_PRINT_FLAG>=2
      if (Dot(force,normal)*f_n < -ZERO_TOL) std::cout << "WARNING : normal force flipped\n";
        if (Dot(force,f_t) < -ZERO_TOL) std::cout << "WARNING : tangential force flipped\n";
#endif
  } else {
    failed_model->Limit_Force( cnfi,gap, rel_disp, slip, normal, dt, area, force );
  }
  Scale(force,area,force);
  return false;
}
