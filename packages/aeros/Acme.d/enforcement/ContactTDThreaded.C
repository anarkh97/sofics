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


#include "ContactTDThreaded.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactErrors.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "contact_assert.h"
#include <cstdio>
#undef  LOCAL_PRINT_FLAG
#define LOCAL_PRINT_FLAG 0
#if (CONTACT_DEBUG_PRINT_LEVEL>=2) || (LOCAL_PRINT_FLAG>0)
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif

ContactTDThreaded::ContactTDThreaded( int ID, int* int_data, 
					Real* real_data, 
					ContactTopology* topo )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_THREADED, topo )
{
  Ntable_id = int_data[0];
  Ttable_id = int_data[1];
  TNtable_id = int_data[2];
  S_i = int_data[3]; // failure_steps
  failed_model_id = int_data[4];
  C_n = real_data[0]; // normal_capacity
  C_t = real_data[1]; // tangential_capacity
  fail_exp = real_data[2]; // failure criterion exponent
  F_n = NULL;
  F_t = NULL;
  F_tn = NULL;
  failed_model = NULL;
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding threaded enforcement model, ID ="<<ID<<std::endl;
    std::cout<<"  failure model id = "<<failed_model_id<<std::endl;
    std::cout<<"  Ntable_id  = "<<Ntable_id<<std::endl;
    std::cout<<"  Ttable_id  = "<<Ttable_id<<std::endl;
    std::cout<<"  TNtable_id = "<<TNtable_id<<std::endl;
    std::cout<<"  fail_exp   = "<<fail_exp<<std::endl;
    std::cout<<"  C_n        = "<<C_n<<std::endl;
    std::cout<<"  C_t        = "<<C_t<<std::endl;
    std::cout<<"  S_i        = "<<S_i<<std::endl;
    std::cout<<std::flush;
  }
#endif
}


ContactTDThreaded::ContactTDThreaded( ContactTopology* topo )
  : ContactTDEnfModel(ContactEnforcement::TD_THREADED, topo)
{
  Ntable_id = 0;
  Ttable_id = 0;
  TNtable_id = 0;
  S_i = 0;
  failed_model_id = 0;
  C_n = 0.0;
  C_t = 0.0;
  fail_exp = 2;
  F_n = NULL;
  F_t = NULL;
  F_tn = NULL;
  failed_model = NULL;
}

ContactTDThreaded::~ContactTDThreaded()
{
}

ContactSearch::ContactErrorCode
ContactTDThreaded::Initialize_Model( int num_models, ContactEnfModel** models,
				      int num_tables, 
				      ContactTable** tables )
{
  // search for tables 
  for( int i=0 ; i<num_tables ; ++i){
	if( Ntable_id == tables[i]->ID() ){
	  F_n = tables[i];
	}
	if( Ttable_id == tables[i]->ID() ){
	  F_t = tables[i];
	}
	if( TNtable_id == tables[i]->ID() ){
	  F_tn = tables[i];
	}
  }
  if( !F_n || !F_t || !F_tn )
	return ContactSearch::INVALID_ID;

  // search for failed model
  for( int i=0 ; i<num_models ; ++i){
    if( failed_model_id == models[i]->ID() ){
      failed_model = (ContactTDEnfModel*) models[i];
      return ContactSearch::NO_ERROR;
    }
  }
  return ContactSearch::INVALID_ID;
}

void ContactTDThreaded::Initialize_for_Time_Step()
{
  for( int i=0 ; i<number_of_nodes ; ++i)
    Node_State_Data(i,0)[DECREASED_STRENGTH_THIS_STEP] = 0.0;
}


int ContactTDThreaded::Num_Node_State_Variables()
{
  // Comment out the failed mode because of restart for now.  We know its 0
  return SIZE_NODE_STATE_DATA;//+failed_model->Num_Node_State_Variables();
}

int ContactTDThreaded::Num_Interaction_State_Variables()
{
  // Comment out the failed mode because of restart for now.  We know its 0
  return 0;//+failed_model->Num_Interaction_State_Variables();
}

int ContactTDThreaded::Interaction_Type( ContactNodeEntityInteraction* cnfi )
{
  if( Node_State_Data( cnfi->Node()->fcs_index, 0 )[STRENGTH] > 0 ) {
    if (cnfi->Get_Type() == ContactNodeEntityInteraction::NODE_SURFACE_INTERACTION) {
      return true;
    } else {
      return TDEM_TIED;
    }
  }
  else {
    return failed_model->Interaction_Type( cnfi );
  }
}

int ContactTDThreaded::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  words_added += Extract_Restart_State_Data( restart_data );
  restart_data[words_added++] = id;
  restart_data[words_added++] = Ntable_id;
  restart_data[words_added++] = Ttable_id;
  restart_data[words_added++] = TNtable_id;
  restart_data[words_added++] = S_i; 
  restart_data[words_added++] = failed_model_id;
  restart_data[words_added++] = C_n;
  restart_data[words_added++] = C_t;
  restart_data[words_added++] = fail_exp;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDThreaded::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  words_read += Implant_Restart_State_Data( restart_data );
  id = (int) restart_data[words_read++];
  Ntable_id = (int) restart_data[words_read++];
  Ttable_id = (int) restart_data[words_read++];
  TNtable_id = (int) restart_data[words_read++];
  S_i = (int) restart_data[words_read++];
  failed_model_id  = (int) restart_data[words_read++];
  C_n = restart_data[words_read++];
  C_t = restart_data[words_read++];
  fail_exp = restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDThreaded::Extract_Nodal_Restart_Variable(int n, Real* data )
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
      }
      else data[i] = 0.0;
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
ContactTDThreaded::Implant_Nodal_Restart_Variable(int n, const Real* data )
{
  if (  topology->Have_Shells()) {
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

bool ContactTDThreaded::Active_Interaction(ContactNodeEntityInteraction* cnfi, Real gap) {   
  Real& strength=Node_State_Data( cnfi->Node()->fcs_index, 0 )[STRENGTH];
  if (strength > 0) {
    return true;
  } else {
    return failed_model->Active_Interaction( cnfi,gap);
  }
} 
  
bool ContactTDThreaded::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{
#if LOCAL_PRINT_FLAG>=2
  Real f_n = Dot(force,normal);
  Real f_t[3];
  Remove_Component(force,normal,f_t);
  Real s_n = Dot(rel_disp,normal);
  Real s_t[3];
  Remove_Component(rel_disp,normal,s_t);
  std::cout << " f_n " << f_n << "," << " s_n " << s_n << "\n";
#endif
  Scale(force,1/area,force);
  Real& strength=Node_State_Data( cnfi->Node()->fcs_index, 0 )[STRENGTH];
  if (strength > 0) {
    Real& already_decremented = Node_State_Data(cnfi->Node()->fcs_index,0)
      [DECREASED_STRENGTH_THIS_STEP];
    Real s_n = Dot(rel_disp,normal);
    Real s_t[3];
    Remove_Component(rel_disp,normal,s_t);
    Real mag_s_t = Magnitude(s_t);
    if(strength == S_i) { // still welded
#if LOCAL_PRINT_FLAG>1
        std::cout << " THREADED : WELDED \n";
#endif
      Real f_n = Dot(force,normal);
      //assume table has positive (x,y)
      if (f_n < 0) f_n = -F_n->Interpolate_Value(s_n); 
      Real norm_ratio = std::fabs(std::max(-f_n,0.0)/C_n); // no limit for compression
      Real mag_f_t = -(F_t->Interpolate_Value(mag_s_t))
                     *(F_tn->Interpolate_Value(s_n)); // Table lookup
      Real tang_ratio = std::fabs(mag_f_t/C_t);
      Real failure_crit = std::pow(norm_ratio,fail_exp) + std::pow(tang_ratio,fail_exp);
#if LOCAL_PRINT_FLAG>=2
      std::cout << "Strength: " << strength 
	      << ", n_r: " << norm_ratio << ", t_r: "<< tang_ratio 
	      << ", fail_crit: " << failure_crit << "\n";
#endif
      if ( failure_crit >= 1.0) { // beginning failure
#if LOCAL_PRINT_FLAG>1
        std::cout << " BEGINNING FAILURE \n";
#endif
        if( already_decremented == 0.0 ){
          strength -= 1.0;
          already_decremented = 1.0;
        }
        // store failure values 
        Real& f_n_fail =Node_State_Data( cnfi->Node()->fcs_index, 0 )
           [NORMAL_FORCE_AT_FAILURE];
        Real& f_t_fail =Node_State_Data( cnfi->Node()->fcs_index, 0 )
           [TANGENTIAL_FORCE_AT_FAILURE];
        f_n_fail = f_n;
        f_t_fail = mag_f_t; // this is negative
      }
      if (mag_s_t > ZERO_TOL) {
        Linear_Combination(f_n,normal,mag_f_t/mag_s_t,s_t,force);
      } else {
        Scale(normal,f_n,force);
      }

    }
    else { // failing, scale (stored) failure values
#if LOCAL_PRINT_FLAG>1
        std::cout << " FAILING \n";
#endif
      Real scale = strength/S_i;
      Real& f_n_fail =Node_State_Data( cnfi->Node()->fcs_index, 0 )
           [NORMAL_FORCE_AT_FAILURE];
      Real f_n = (s_n > 0) 
               ? scale*f_n_fail:    //failure
                 Dot(force,normal); //impenetrability
#if LOCAL_PRINT_FLAG>=2
      std::cout << "Strength: " << strength 
	      << " scale: " << scale << ", f_n_fail: " << f_n_fail <<'\n';
#endif
      if (mag_s_t > ZERO_TOL) { 
        Real& f_t_fail =Node_State_Data( cnfi->Node()->fcs_index, 0 )
           [TANGENTIAL_FORCE_AT_FAILURE];
#if LOCAL_PRINT_FLAG>=2
      std::cout << "f_t_fail: " << f_t_fail << '\n'; 
#endif
        Linear_Combination(f_n,normal,scale*f_t_fail/mag_s_t,s_t,force); 
      } else {
        Scale(normal,f_n,force); 
      }
      if( already_decremented == 0.0 ){
        strength -= 1.0;
        already_decremented = 1.0;
      }
      if (strength <= 0.0) {
        if (Needs_Glued_Search(failed_model->Interaction_Type(cnfi))) {
          cnfi->Is_Glued( true );
        } else {
          cnfi->Is_Glued( false );
        }
      }
    }
  }
  else { // completely failed
#if LOCAL_PRINT_FLAG>1
        std::cout << " COMPLETELY FAILED \n";
#endif
    failed_model->Limit_Force( cnfi,gap, rel_disp, slip, normal, dt, area, force );
  }
#if LOCAL_PRINT_FLAG>=2
  if (Dot(force,normal)*f_n < -ZERO_TOL) { 
	  std::cout << "WARNING : normal force flipped\n";
	  std::cout << "  normal force : "  << Dot(force,normal)
	       << ", original : " << f_n << "\n";
  }
  if (Dot(force,f_t) < -ZERO_TOL && Magnitude(f_t) > ZERO_TOL) {
	  std::cout << "WARNING : tangential force flipped\n";
	  std::cout << "  force : "  << force[0] <<", "<<force[1]<<", "<<force[2]
	       << "\n";
	  std::cout << "  original f_t : " << f_t[0] <<", "<<f_t[1]<<", "<<f_t[2] 
	       << "\n";
  }
  if (Dot(force,rel_disp) > ZERO_TOL && Magnitude(rel_disp) > ZERO_TOL) {
	  std::cout << "WARNING : force does not oppose constraint\n";
	  std::cout << "  force : "  << force[0] <<", "<<force[1]<<", "<<force[2]
	       << "\n";
	  std::cout << "  rel_disp : " << rel_disp[0] <<", "<<rel_disp[1]<<", "<<rel_disp[2] 
	       << "\n";
  }
#endif
  Scale(force,area,force);
  return false;
}

