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


#ifndef ContactTDSpotWeld_h_
#define ContactTDSpotWeld_h_

#include "ContactTDEnfModel.h"
#include "ContactSearch.h"
#include "contact_assert.h"

class ContactTDSpotWeld : public ContactTDEnfModel {

 public:
  ContactTDSpotWeld( int ID, int* int_data, Real* real_data, 
		      ContactTopology* topology );
  ContactTDSpotWeld(ContactTopology* topology);
  ~ContactTDSpotWeld();

  virtual int Num_Interaction_State_Variables();
  virtual int Num_Node_State_Variables();
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
                      int num_tables, ContactTable** tables  );
  // virtual void Initialize_Interaction_State_Data(Real* state_data);
  virtual void Initialize_Node_State_Data(Real* state_data)
    { state_data[STATUS]   = 0.0;
      state_data[STRENGTH] = S_i;
      state_data[DECREASED_STRENGTH_THIS_STEP] = 0.0;};
  virtual int Interaction_Type( ContactNodeEntityInteraction* );
  virtual bool Active_Interaction (ContactNodeEntityInteraction* ,Real );
  virtual bool Limit_Force(ContactNodeEntityInteraction*, Real, Real*, Real*,Real*,Real, Real, Real*);

  virtual int Restart_Size() { return 5+Restart_Size_State_Data(); };
  virtual int Extract_Restart_Data(Real* restart_data);
  virtual int Implant_Restart_Data(const Real* restart_data);
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable(int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable(int n, const Real* data );
  void Initialize_for_Time_Step();

 private:
  int failed_model_id;
  ContactTDEnfModel* failed_model;
  Real C_n; // normal capacity
  Real C_t; // tangential capacity
  int  S_i; // initial strength (i.e. number of timesteps to complete failure)

  enum NODE_STATE_DATA_HANDLE { STATUS=0, STRENGTH, 
				DECREASED_STRENGTH_THIS_STEP,
				SIZE_NODE_STATE_DATA };
};

#endif
