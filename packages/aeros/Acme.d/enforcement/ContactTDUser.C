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


#include "ContactTDUser.h"
#include "ContactTDUserStubs.h"
#include "ContactTDEnforcement.h"
#include "ContactErrors.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "contact_assert.h"
#include "ContactNodeEntityInteraction.h"
#include <cstdio>

ContactTDUser::ContactTDUser( int ID, int* int_data, 
                              Real* real_data,
                              ContactTopology* topology,
                              ContactEnforcement* enf )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_USER, topology )
{
  //  int_data[0]      number of integer user data entries, NI
  //  int_data[1]      number of real user data entries, NR
  //  int_data[2]      number of nodal state variables
  //  int_data[3]      failure model id
  //  int_data[4]      user integer data entry 1
  //  int_data[5]      user integer data entry 2
  //  .
  //  .
  //  .
  //  int_data[NI+3]   user integer data entry NI
  //
  //  real_data[0]     user real data entry 1
  //  real_data[1]     user real data entry 2
  //  .
  //  .
  //  .
  //  real_data[NR-1]  user real data entry NR
  //
  enforcement     = enf;
  num_i           = int_data[0];
  num_r           = int_data[1];
  num_nstate_vars = int_data[2];
  failed_model_id = int_data[3];
  if (num_i>0) {
    idata = new int[num_i];
    std::memcpy( idata, &int_data[4], num_i*sizeof(int)); 
  } else {
    idata = NULL;
  }
  if (num_r>0) {
    rdata = new Real[num_r];
    std::memcpy( rdata, real_data, num_r*sizeof(Real)); 
  } else {
    rdata = NULL;
  }
  failed_model                    = NULL;
  InitializeModel                 = FORTRAN(user_initialize_model);
  InitializeTimeStep              = FORTRAN(user_initialize_time_step);
  InitializeNodeStateData         = FORTRAN(user_init_node_state_data);
  Is_Active                       = FORTRAN(user_is_active);
  LimitForce                      = FORTRAN(user_limit_force);
  InteractionType                 = FORTRAN(user_interaction_type);
}

ContactTDUser::ContactTDUser(ContactTopology* topology)
  : ContactTDEnfModel(ContactEnforcement::TD_USER, topology)
{
  num_i                           = 0;
  num_r                           = 0;
  num_nstate_vars                 = 0;
  idata                           = NULL;
  rdata                           = NULL;
  failed_model                    = NULL;
  InitializeModel                 = FORTRAN(user_initialize_model);
  InitializeTimeStep              = FORTRAN(user_initialize_time_step);
  InitializeNodeStateData         = FORTRAN(user_init_node_state_data);
  Is_Active                       = FORTRAN(user_is_active);
  LimitForce                      = FORTRAN(user_limit_force);
  InteractionType                 = FORTRAN(user_interaction_type);
}

ContactTDUser::~ContactTDUser()
{
  if (idata) delete [] idata;
  if (rdata) delete [] rdata;
}

int ContactTDUser::Extract_Restart_Data( Real* restart_data )
{
  int i;
  int words_added = 0;
  words_added += Extract_Restart_State_Data( restart_data );
  restart_data[words_added++] = failed_model_id;
  restart_data[words_added++] = num_i;
  restart_data[words_added++] = num_r;
  for (i=0; i<num_i; ++i)
    restart_data[words_added++] = idata[i];
  for (i=0; i<num_r; ++i)
    restart_data[words_added++] = rdata[i];
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDUser::Implant_Restart_Data( const Real* restart_data )
{
  int i;
  int words_read = 0;
  words_read += Implant_Restart_State_Data( restart_data );
  failed_model_id = (int) restart_data[words_read++];
  num_i = (int) restart_data[words_read++];
  num_r = (int) restart_data[words_read++];
  if (num_i>0) {
    idata = new int[num_i];
    for (i=0; i<num_i; ++i)
      idata[i] = (int) restart_data[words_read++];
  }
  if (num_r>0) {
    rdata = new Real[num_r];
    for (i=0; i<num_r; ++i)
      rdata[i] = restart_data[words_read++];
  }
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDUser::Extract_Nodal_Restart_Variable(int n, Real* data )
{ 
  if ( topology->Have_Shells()) {
    ContactShellHandler * shell_handler = topology->Shell_Handler();

    if( n>Num_Node_State_Variables()*shell_handler->
	Max_Num_Acme_Nodes_for_Host_Node() 
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
ContactTDUser::Implant_Nodal_Restart_Variable(int n, const Real* data )
{
  if (  topology->Have_Shells()) {
    ContactShellHandler * shell_handler = topology->Shell_Handler();

    if( n>Num_Node_State_Variables()*shell_handler-> 
	Max_Num_Acme_Nodes_for_Host_Node() 
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

ContactSearch::ContactErrorCode
ContactTDUser::Initialize_Model( int num_models, 
                                 ContactEnfModel** models,
                                 int  num_tables ,
                                 ContactTable**  tables )
{
  int status;
  // search for failed model
  if (failed_model_id>=0) {
    for( int i=0 ; i<num_models ; ++i){
      if( failed_model_id == models[i]->ID() ){
        failed_model = (ContactTDEnfModel*) models[i];
        break;
      }
    }
    if (failed_model==NULL) {
      return ContactSearch::INVALID_ID;
    }
  }
  InitializeModel(reinterpret_cast<void*>(enforcement), 
                  id, rdata, idata, status);
  if( !status ) {
    return ContactSearch::INTERNAL_ERROR;
  }
  return ContactSearch::NO_ERROR;
}

void ContactTDUser::Initialize_for_Time_Step() 
{
  int status;
  InitializeTimeStep(reinterpret_cast<void*>(enforcement), 
                     id, rdata, idata, status);
};

void ContactTDUser::Initialize_Node_State_Data(Real* state_data) 
{
  int status;
  InitializeNodeStateData(reinterpret_cast<void*>(enforcement), 
                          id, rdata, idata, state_data, status);
};

int ContactTDUser::Num_Node_State_Variables()
{
  return num_nstate_vars;
}

int ContactTDUser::Num_Interaction_State_Variables()
{
  return 0;  // this is an unused function at the moment
}

bool ContactTDUser::Limit_Force(ContactNodeEntityInteraction* nfi,
                                Real gap, Real* rel_disp, Real* slip, 
                                Real* normal, Real dt, Real area, Real* force) 
{ 
  int status;
  int itype = nfi->Get_Type();
  int index = nfi->Node()->fcs_index+1;
  LimitForce(reinterpret_cast<void*>(enforcement),
             reinterpret_cast<void*>(nfi),  
             id, rdata, idata, itype, index, gap, rel_disp, 
             slip, normal, dt, area, force, status);
  if( status == 1 ) {
    return true;
  } else if (status==0) {
    return false;
  }
  if (failed_model) {
    return failed_model->Limit_Force( nfi, gap, rel_disp, 
                                      slip, normal, dt, 
                                      area, force );
  }
  return false;
}

bool
ContactTDUser::Active_Interaction(ContactNodeEntityInteraction* nfi, Real gap)
{
  int status;
  int itype = nfi->Get_Type();
  int index = nfi->Node()->fcs_index+1;
  Is_Active(reinterpret_cast<void*>(enforcement), 
            reinterpret_cast<void*>(nfi), 
            id, rdata, idata, itype, index, gap, status);
  if( status == 1 ) {
    return true;
  } else if (status==0) {
    return false;
  }
  if (failed_model) {
    return failed_model->Active_Interaction( nfi, gap );
  }
  return false;
}

int
ContactTDUser::Interaction_Type( ContactNodeEntityInteraction* nei )
{
  int status;
  int itype = nei->Get_Type();
  int index = nei->Node()->fcs_index+1;
  InteractionType(reinterpret_cast<void*>(enforcement), 
             reinterpret_cast<void*>(nei),   
             id, rdata, idata, itype, index, status);
  if( status > -1 ) {
     if ( status < ContactTDEnfModel::TDEM_NUM_MODEL_TYPES) return status;
  } else {
     if (failed_model) return failed_model->Interaction_Type( nei );
  }
  return 0;
}


void
ContactTDUser::Set_Initialize_Model_Fn(CONTACT_INIT_MODEL_FN *fn)
{
  InitializeModel = fn;
}

void
ContactTDUser::Set_Initialize_Time_Step_Fn(CONTACT_INIT_TIME_STEP_FN *fn)
{
  InitializeTimeStep = fn;
}

void
ContactTDUser::Set_Initialize_Node_State_Data_Fn(CONTACT_INIT_NODE_STATE_DATA_FN *fn)
{
  InitializeNodeStateData = fn;
}

void
ContactTDUser::Set_Limit_Force_Fn(CONTACT_LIMIT_FORCE_FN *fn)
{
  LimitForce = fn;
}

void
ContactTDUser::Set_Active_Fn(CONTACT_INTERACTION_ACTIVE_FN *fn)
{
  Is_Active = fn;
}

void
ContactTDUser::Set_Interaction_Type_Fn(CONTACT_INTERACTION_TYPE_FN *fn)
{
  InteractionType = fn;
}

