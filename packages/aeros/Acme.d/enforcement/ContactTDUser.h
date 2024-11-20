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


#ifndef ContactTDUser_h_
#define ContactTDUser_h_

#include "ContactTDEnfModel.h"
#include "ContactTDUserSubTypes.h"
#include "ContactTable.h"
#include "ContactSearch.h"

class ContactTDUser : public ContactTDEnfModel {

 public:
  ContactTDUser( int ID, int* int_data, Real* real_data,
			 ContactTopology* topology,
                         ContactEnforcement* enf );
  ContactTDUser( ContactTopology* topology );
  ~ContactTDUser();
  
  int Restart_Size() { return 4+num_i+num_r+Restart_Size_State_Data(); };
  int Extract_Restart_Data(double* restart_data);
  int Implant_Restart_Data(const double* restart_data);
  ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable(int n, Real* data );
  ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable(int n, const Real* data );

  ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables );
  void Initialize_for_Time_Step();
  void Initialize_Node_State_Data(Real* state_data);
  
  int Num_Node_State_Variables();
  int Num_Interaction_State_Variables();
  
  int  Interaction_Type(ContactNodeEntityInteraction*);

  bool Limit_Force(ContactNodeEntityInteraction*,
        Real, Real*, Real*,Real*,Real, Real, Real*);
        
  bool Active_Interaction(ContactNodeEntityInteraction*,Real);
  
  void Set_Initialize_Model_Fn(CONTACT_INIT_MODEL_FN *fn);
  void Set_Initialize_Time_Step_Fn(CONTACT_INIT_TIME_STEP_FN *fn);
  void Set_Initialize_Node_State_Data_Fn(CONTACT_INIT_NODE_STATE_DATA_FN *fn);
  void Set_Interaction_Type_Fn(CONTACT_INTERACTION_TYPE_FN *fn);
  void Set_Active_Fn(CONTACT_INTERACTION_ACTIVE_FN *fn);
  void Set_Limit_Force_Fn(CONTACT_LIMIT_FORCE_FN *fn);
  
  ContactTDEnfModel* Failure_Model() {return failed_model;};

 private:

  ContactEnforcement* enforcement;
  ContactTDEnfModel*  failed_model;
  int   failed_model_id;
  
  int   num_i;
  int   num_r;
  int   num_nstate_vars;
  int*  idata;
  Real* rdata;
  
  CONTACT_INIT_MODEL_FN *                 InitializeModel;
  CONTACT_INIT_TIME_STEP_FN *             InitializeTimeStep;
  CONTACT_INIT_NODE_STATE_DATA_FN *       InitializeNodeStateData;
  CONTACT_INTERACTION_ACTIVE_FN *         Is_Active;
  CONTACT_LIMIT_FORCE_FN *                LimitForce;
  CONTACT_INTERACTION_TYPE_FN *           InteractionType;

};
#endif
