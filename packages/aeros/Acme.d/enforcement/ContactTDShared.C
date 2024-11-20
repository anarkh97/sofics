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


/*
ContactTDShared: provides an interface to a class of shared models -- Simod 
(i.e., models that can be used by ACME and application codes that have different
interface descriptors such as boundary conditions and interface elements).
Reference the header file for additional information.
   
Search for "unresolved" to find unresolved implementation issues!!!!!!!!!!!!!!!
*/
   
#include "ContactTDShared.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactErrors.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "contact_assert.h"
#include <cstdio>

#undef LOCAL_PRINT_FLAG
#define LOCAL_PRINT_FLAG 0

#ifdef CONTACT_SIMOD_SUPPORT
#include "simod_time_backdoor.h"
  using simod_util_spc::extracted_time;

/* constructors: =========================================================== */

ContactTDShared::ContactTDShared(int ID, int* int_data, 
					Real* real_data, 
					ContactTopology* topology )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_SHARED, topology ),
    failed_model(NULL),          /* set within Initialize_Model */
    failed_model_id(int_data[0])

/* assumption: the shared model has already been created by a host code. */
{
}

/* ------------------------------------------------------------------------- */

// if this for restart how do I get the function pointer
ContactTDShared::ContactTDShared( ContactTopology* topology )
  : ContactTDEnfModel(ContactEnforcement::TD_SHARED, topology),
    failed_model(NULL),          /* set within Initialize_Model */
    failed_model_id(0),
    shared_model(NULL)
{
}

/* destructor: ============================================================= */

ContactTDShared::~ContactTDShared()
{
}

/* accessors & mutators: initialization ==================================== */

ContactSearch::ContactErrorCode ContactTDShared::Initialize_Model
  (int num_models, 
   ContactEnfModel** models,
   int            /* num_tables */,
   ContactTable** /* tables */)
{
  // check if the host code has set the embedded model
  if (shared_model == NULL) {
	return  ContactSearch::ID_NOT_FOUND;
  }
  // match id with an instantiated model
  for( int i=0 ; i<num_models ; ++i){
    if(failed_model_id == models[i]->ID() ){
      failed_model = (ContactTDEnfModel*) models[i];
      return ContactSearch::NO_ERROR;
    }
  }
  return ContactSearch::INVALID_ID;
}

/* ------------------------------------------------------------------------- */

void ContactTDShared::Initialize_Node_State_Data(Real* state_data)
{
  shared_model->initializeInternalVariables(state_data);
}

/* ------------------------------------------------------------------------- */

void ContactTDShared::Initialize_for_Time_Step()
  // unresolved: do any of the Simod models need to do something here?
{
}

/* accessors & mutators: state variable queries ============================ */

/* unresolved: what is the distinction between interaction state variables and
   nodel s.v.?   Should the former always be zero for Simod models as used
   below? */
int ContactTDShared::Num_Interaction_State_Variables()
 // except for the name, this function is unchanged from ContactTDAreaWeld  
{
  // rj -- Comment out the failed mode because of restart for now.  
  return 0; //+failed_model->Num_Interaction_State_Variables();
}

/* ------------------------------------------------------------------------- */

int ContactTDShared::Num_Node_State_Variables()
{
  // rj -- Comment out the failed mode because of restart for now. 
  return shared_model->num_internal_var();//+failed_model->Num_Node_State_Variables();
}

/* accessors: model type queries =========================================== */

int ContactTDShared::Interaction_Type( ContactNodeEntityInteraction* cnfi )
{
  // a pointer to the nodal state data 
  Real* state_data = Node_State_Data( cnfi->Node()->fcs_index, 0 );
  if( shared_model->active(state_data) )
    return shared_model->Interaction_Type();
  else
    return failed_model->Interaction_Type(cnfi);
}

/* accessors: active model queries ========================================= */

bool ContactTDShared::Active_Interaction
(ContactNodeEntityInteraction* cnfi, Real gap)
{
 /* define a pointer to the nodal state data that can be passed
     to the shared model */
  Real* state_data = Node_State_Data( cnfi->Node()->fcs_index, 0 );
  if(shared_model->active(state_data))
    return true;
  else
    return failed_model->Active_Interaction(cnfi,gap);
} 
  
/* accessors & mutators: restart data ====================================== */

int ContactTDShared::Extract_Restart_Data( Real* restart_data )
{
  // extract all state data for nodes associated with this model
  int words_added = Extract_Restart_State_Data(restart_data);
  restart_data[words_added++] = id;
  restart_data[words_added++] = failed_model_id;
  // save the pointer to the current, soon to be old, model object
  long int temp = reinterpret_cast<long int>(shared_model);
  restart_data[words_added++] = static_cast<double>(temp);
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

/* ------------------------------------------------------------------------- */

int ContactTDShared::Implant_Restart_Data( const Real* restart_data )
{
  // implant all state data for nodes associated with this model
  int words_read = Implant_Restart_State_Data( restart_data );
  id = (int) restart_data[words_read++];
  failed_model_id = (int) restart_data[words_read++];
  // obtain the previous pointer to the model, i.e., before restart
  long int temp = static_cast<long int>(restart_data[words_read++]);
  SharedInterfaceModel* old_model = reinterpret_cast<SharedInterfaceModel*>(temp);
  // give simod the old model pointer to identify the new model pointer
  bool simod_ready_for_implant = false;
  PRECONDITION(simod_ready_for_implant);
  // shared_model = SharedInterfaceModel::new_model(old_model);
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

/* ------------------------------------------------------------------------- */

ContactSearch::ContactErrorCode
ContactTDShared::Extract_Nodal_Restart_Variable(int n, Real* data )
  // except for the name, this function is unchanged from ContactTDAreaWeld
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

/* ------------------------------------------------------------------------- */

ContactSearch::ContactErrorCode
ContactTDShared::Implant_Nodal_Restart_Variable(int n, const Real* data )
// except for the name, this function is unchanged from ContactTDAreaWeld
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

/* accessors & mutators: force calculations ============================= */

bool ContactTDShared::Limit_Force
  (ContactNodeEntityInteraction* cnfi,
   Real gap, 
   Real* rel_disp,
   Real* slip,
   Real* normal,
   Real dt,
   Real area,
   Real* force) 
{
 /* define a pointer to the nodal state data that can be passed
     to the shared model */
  Real* state_data = Node_State_Data( cnfi->Node()->fcs_index, 0 );
  Real time = extracted_time();
  bool activeModel = shared_model->calcForces3d
                       (state_data,gap,rel_disp,slip,normal,time,dt,area,force);
  if (!activeModel) // use the linked model
    activeModel = failed_model->Limit_Force(cnfi,gap,rel_disp,slip,normal,dt,area,force);
  return activeModel;
  // unresolved: why does the AreaWeld model not use the "activeModel result" from
  //             the failed_model and always returns false for this function?
}

/* ------------------------------------------------------------------------- */

#else
// dummy model
ContactTDShared::ContactTDShared(int ID, int* int_data,
                                        Real* real_data,
                                        ContactTopology* topology )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_SHARED, topology )
{
	std::cerr << "ERROR : dummy model, no SIMOD support in ACME\n";
}

ContactTDShared::ContactTDShared(ContactTopology* topology )
  : ContactTDEnfModel(ContactEnforcement::TD_SHARED, topology)
{
	std::cerr << "ERROR : dummy model, no SIMOD support in ACME\n";
}

ContactTDShared::~ContactTDShared() {}

ContactSearch::ContactErrorCode 
ContactTDShared::Initialize_Model(int , ContactEnfModel** , int , ContactTable** ) 
{
    return ContactSearch::NO_ERROR;
}

void ContactTDShared::Initialize_Node_State_Data(Real* ) {}

void ContactTDShared::Initialize_for_Time_Step() {;}
int ContactTDShared::Num_Interaction_State_Variables(){return 0;}
int ContactTDShared::Num_Node_State_Variables() {return 0;}
int  ContactTDShared::Interaction_Type(ContactNodeEntityInteraction* ) {return 0;}
bool ContactTDShared::Active_Interaction( ContactNodeEntityInteraction* , Real ) {return false;}
int ContactTDShared::Restart_Size() { return 0; }
int ContactTDShared::Extract_Restart_Data(Real* ) { return 0; }
int ContactTDShared::Implant_Restart_Data(const Real* ) { return 0; }
ContactSearch::ContactErrorCode ContactTDShared::Extract_Nodal_Restart_Variable (int n, Real* data ) { return ContactSearch::INVALID_ID;}
ContactSearch::ContactErrorCode ContactTDShared::Implant_Nodal_Restart_Variable(int , const Real* ) { return ContactSearch::INVALID_ID;}
bool ContactTDShared::Limit_Force (ContactNodeEntityInteraction* ,Real  , Real* , Real* , Real* , Real  , Real  ,  Real* ) {return false;}

#endif 
