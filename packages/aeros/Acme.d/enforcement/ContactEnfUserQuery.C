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


#include "ContactEnfUserQuery.h"
#include "ContactEnforcement.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactEnfModel.h"
#include "ContactTDUser.h"
#include "ContactTable.h"
#include "ContactNodeEntityInteraction.h"
#include "contact_assert.h"

#include <iostream>
#include <cstring>

void FORTRAN(userquery_table_last_abscissa)(void *obj, int *id, Real* value)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)obj;
  ContactSearch*      search      = enforcement->Search();
  for( int i=0 ; i<search->Num_Tables() ; ++i ){
    if( *id == search->Tables()[i]->ID() ){
      *value = search->Tables()[i]->Last_Abscissa();
      break;
    }
  }
}  

void FORTRAN(userquery_table_interpolate)(void *obj, int *id, 
                                          Real* abscissa, Real* ordinate)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)obj;
  ContactSearch*      search      = enforcement->Search();
  for( int i=0 ; i<search->Num_Tables() ; ++i ){
    if( *id == search->Tables()[i]->ID() ){
      *ordinate = search->Tables()[i]->Interpolate_Value(*abscissa);
      break;
    }
  }
}  

void FORTRAN(userquery_number_of_nodes)(void *obj, int *nnodes)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)obj;
  ContactTopology*    topology    = enforcement->Topology();
  *nnodes = topology->Number_of_Nodes();
}
  
void FORTRAN(userquery_node_state_data)(void *obj, int* id, int* node, 
                                        int* offset, Real* value)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)obj;
  int nmodels = enforcement->Num_Models();
  ContactEnfModel** models = enforcement->Models();
  for( int i=0 ; i<nmodels ; ++i ){
    if( *id == models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (models[i]);
      int c_node   = *node-1;
      int c_offset = *offset-1;
      *value = model->Get_Node_State_Data(c_node, c_offset);
      break;
    }
  }
}

void FORTRAN(userset_node_state_data)(void *obj, int* id, int* node,
                                      int* offset, Real* value)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)obj;
  int nmodels = enforcement->Num_Models();
  ContactEnfModel** models = enforcement->Models();
  for( int i=0 ; i<nmodels ; ++i ){
    if( *id == models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (models[i]);
      int c_node   = *node-1;
      int c_offset = *offset-1;
      model->Set_Node_State_Data(*value, c_node, c_offset);
      break;
    }
  }
}

void FORTRAN(userset_nfi_failure_tied)(void *enf_obj, int* model_id, 
                                       void *nfi_obj)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)enf_obj;
  ContactNodeEntityInteraction* nfi = (ContactNodeEntityInteraction *)nfi_obj;
  int nmodels = enforcement->Num_Models();
  ContactEnfModel** models = enforcement->Models();
  for( int i=0 ; i<nmodels ; ++i ){
    if( *model_id == models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>(models[i]);
      nfi->Is_Glued( model->Failure_Model()->Needs_Glued_Search( model->Failure_Model()->Interaction_Type(nfi) ) );
      break;
    }
  }
}
