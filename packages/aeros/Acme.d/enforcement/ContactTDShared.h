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
   (i.e., models that can be used by ACME and application codes that have 
   different interface descriptors such as boundary conditions and interface 
   elements).
   This interface can be used with models that have the most general form
   supported by ACME TDEnforcement models (i.e., they receive all of the 
   kinematic data and return interfacial forces).  
   
   ContactTDAreaWeld was used as a template for this model with reference to 
   other enforcement models as well.
   
   Initial version: Feb 03
   Argument requirements: (1) int_data[0]: linked/failed model ID number
*/

#ifndef ContactTDShared_h_
#define ContactTDShared_h_

#include "ContactTDEnfModel.h"
#include "ContactSearch.h"
#include "contact_assert.h"
#ifdef CONTACT_SIMOD_SUPPORT
#include "SharedInterfaceModel.h"
  using SharedInterfaceModel_spc::SharedInterfaceModel;
  
  
class ContactTDShared : public ContactTDEnfModel {

  public:
  
  /* constructors & destructor */
  ContactTDShared(int ID,
                    int* int_data,
                    Real* real_data, 
	            ContactTopology* topology);
 
  ContactTDShared(ContactTopology* topology);
  /* can't define shared model type, so this must be used in conjunction with
         a function that initializes the model's data.  The ContactEnforcement
	 constructor uses this with Implant_Restart_Data. */
                    
  virtual ~ContactTDShared();

  /* model initializers */
  void Link_External_Function(void* function_pointer)
    {shared_model = static_cast<SharedInterfaceModel*>(function_pointer);}

  virtual ContactSearch::ContactErrorCode Initialize_Model
      ( int num_models, 
        ContactEnfModel** models,
        int num_tables, 
        ContactTable** tables  );
  virtual void Initialize_Node_State_Data(Real* state_data);
  virtual void Initialize_for_Time_Step();
  // virtual void Initialize_Interaction_State_Data(Real* state_data);
    
  /* state variable queries */
  virtual int Num_Interaction_State_Variables();
  virtual int Num_Node_State_Variables();    
    
  /* model "type" queries */
  virtual int  Interaction_Type(ContactNodeEntityInteraction* cnfi);
  
  /* active model queries */
  virtual bool Active_Interaction( ContactNodeEntityInteraction* cnfi, Real gap);
 
  /* restart functions */
  /* There are two independent approaches to restart implemented in ACME.
  As examples of host codes that use each, the first one is used by Alegra,
  and the second one is used by Sierra applications.  The first approach
  requires the parameter data (for the ACME model -- not the Simod model)
  to be written to the restart file, while for the second approach the paremeter
  data is reinitiated upon restart.  In the first case the pointer to the new
  version of the Simod model will be known, while in the second case at best
  the old pointer is known.  Simod must be modified to include a function
  for providing the new address of a model given the old address, i.e., to
  map between the old object's address and the new object's address.  */

  /* The following three restart functions implement the first approach.
  The parameter data is extracted from the model and then later implanted
  into the model. */
  virtual int Restart_Size() { return 3+Restart_Size_State_Data(); };
  virtual int Extract_Restart_Data(Real* restart_data);
  virtual int Implant_Restart_Data(const Real* restart_data);
    
  /* These are the restart functions needed by Sierra -- only nodal state 
  variables are saved in the host code. */
  virtual ContactSearch::ContactErrorCode Extract_Nodal_Restart_Variable(int n, Real* data);
  virtual ContactSearch::ContactErrorCode Implant_Nodal_Restart_Variable(int n, const Real* data);
   
  /* surface force calculations */
  // Reference ContactTDEnforcement.C for some variable definitions
  virtual bool Limit_Force
      (ContactNodeEntityInteraction* cnei,
       Real  gap,      // <rel_disp,normal>, >0 for separation
       Real* rel_disp, // relative displacement between two points
       Real* slip,     // velocity_jump*dt
       Real* normal,   // ourward unit normal vector to the face
       Real  dt,       // change in time to current time
       Real  area,     // tributary area associated with the node
       Real* force);   // total force vector, typically on the slave node

 private:
   SharedInterfaceModel* shared_model; /* pointer to the initial shared model */
   int failed_model_id;                /* id of model used after this fails */
   ContactTDEnfModel* failed_model;    /* pointer to the "failed model." */  
/* Note that by linking to this class multiple links can occur -- at a cost */
};
#else
// dummy model 
class ContactTDShared : public ContactTDEnfModel {
  public:
  ContactTDShared(int , int* , Real* , ContactTopology* );
  ContactTDShared(ContactTopology* );
  virtual ~ContactTDShared();
  virtual ContactSearch::ContactErrorCode Initialize_Model
      ( int , ContactEnfModel** , int , ContactTable** );
  virtual void Initialize_Node_State_Data(Real* );
  virtual void Initialize_for_Time_Step();
  virtual int Num_Interaction_State_Variables();
  virtual int Num_Node_State_Variables();
  virtual int  Interaction_Type(ContactNodeEntityInteraction* );
  virtual bool Active_Interaction( ContactNodeEntityInteraction* , Real );
  virtual int Restart_Size();
  virtual int Extract_Restart_Data(Real* );
  virtual int Implant_Restart_Data(const Real*);
  virtual ContactSearch::ContactErrorCode Extract_Nodal_Restart_Variable
       (int n, Real* data );
  virtual ContactSearch::ContactErrorCode Implant_Nodal_Restart_Variable
       (int , const Real* );
  virtual bool Limit_Force (ContactNodeEntityInteraction* ,
     Real  , Real* , Real* , Real* , Real  , Real  ,  Real* );
};
#endif
#endif
