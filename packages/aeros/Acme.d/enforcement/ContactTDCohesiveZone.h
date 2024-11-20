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


#ifndef ContactTDCohesiveZone_h_
#define ContactTDCohesiveZone_h_

#include "ContactTDEnfModel.h"
#include "ContactTable.h"
#include "ContactSearch.h"

class ContactTDCohesiveZone : public ContactTDEnfModel {

 public:
  ContactTDCohesiveZone( int ID, int* int_data, Real* real_data,
			 ContactTopology* topology );
  ContactTDCohesiveZone( ContactTopology* topology );
  ~ContactTDCohesiveZone();

  virtual int Num_Interaction_State_Variables() { return 0; };
  virtual int Num_Node_State_Variables() { return 0; };
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables );
  // virtual void Initialize_Interaction_State_Data(Real*) {};
  virtual void Initialize_Node_State_Data(Real* state_data) {};

  virtual int Restart_Size() { return 4; };
  virtual int Extract_Restart_Data(double* restart_data);
  virtual int Implant_Restart_Data(const double* restart_data);
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable(int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable(int n, const Real* data );

  virtual int  Interaction_Type(ContactNodeEntityInteraction*) 
  {return TDEM_SPRINGY;};
  virtual bool Active_Interaction(ContactNodeEntityInteraction*,Real);
  virtual bool Limit_Force(ContactNodeEntityInteraction*,
        Real, Real*, Real*,Real*,Real, Real, Real*);

 private:
  int table_id;
  ContactTable* F; 
  Real g_n_crit;
  Real g_t_crit;
};
#endif
