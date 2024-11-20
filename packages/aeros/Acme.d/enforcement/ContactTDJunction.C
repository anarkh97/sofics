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


#include "ContactTDJunction.h"
#include "ContactTDEnforcement.h"
#include "ContactNodeEntityInteraction.h"
#include "contact_assert.h"
#undef  LOCAL_PRINT_FLAG
#define LOCAL_PRINT_FLAG 0
#if (CONTACT_DEBUG_PRINT_LEVEL>=2) || (LOCAL_PRINT_FLAG>0)
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif

ContactTDJunction::ContactTDJunction( int ID, int* int_data, 
					      Real* real_data,
					      ContactTopology* topo )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_JUNCTION, topo )
{
   Ntable_id = int_data[0];
   Ttable_id = int_data[1];
   F_n = NULL; // force-displ law
   F_t = NULL; // force-displ law
   NCutoff = 0.0;
   TCutoff = real_data[0]; // normal distance cutoff for tangential force
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding junction enforcement model, ID ="<<ID<<std::endl;
    std::cout<<"  Ntable_id = "<<Ntable_id<<std::endl;
    std::cout<<"  Ttable_id = "<<Ttable_id<<std::endl;
    std::cout<<"  NCutoff   = "<<NCutoff<<std::endl;
    std::cout<<"  TCutoff   = "<<TCutoff<<std::endl;
    std::cout<<std::flush;
  }
#endif
}

ContactTDJunction::ContactTDJunction(ContactTopology* topo)
  : ContactTDEnfModel(ContactEnforcement::TD_JUNCTION, topo)
{
   Ntable_id = 0;
   Ttable_id = 0;
   F_n = NULL; // force-displ law
   F_t = NULL; // force-displ law
   NCutoff = 0.0; // normal distance cutoff for tangential force
   TCutoff = 0.0; 
}

ContactSearch::ContactErrorCode
ContactTDJunction::Initialize_Model( int /* num_models */, 
                                    ContactEnfModel** /* models */,
                                      int  num_tables ,
                                      ContactTable**  tables )
{
  // initialize normal force-displ. law and cutoff from the table
  for( int i=0 ; i<num_tables ; ++i){
    if( Ntable_id == tables[i]->ID() ){
      F_n = tables[i]; 
    }
    if( Ttable_id == tables[i]->ID() ){
      F_t = tables[i];
    }
  }
  if( !F_n || !F_t )
    return ContactSearch::INVALID_ID;
  else
   NCutoff = F_n->Last_Abscissa();
  return ContactSearch::NO_ERROR;
}


int ContactTDJunction::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = Ntable_id;
  restart_data[words_added++] = Ttable_id;
  restart_data[words_added++] = TCutoff;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDJunction::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  Ntable_id = (int) restart_data[words_read++];
  Ttable_id = (int) restart_data[words_read++];
  TCutoff   = (int) restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDJunction::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model 
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDJunction::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model 
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactTDJunction::~ContactTDJunction(){}

bool ContactTDJunction::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{ 
  bool modified = false;
  Real g = Dot(rel_disp,normal); // since the gap read in has been manipulated
  Real f_n = Dot(force,normal);
  Real f_t[3];
  Add(force,-f_n,normal,f_t);
  Real mag_t = Magnitude(f_t); // magnitude of f_t (not altered)
  Real mag_f_t = mag_t;  // tangential force magnitude (pot. altered)
  // determine normal force 
  if (f_n < 0)  { // tension
    if (g < NCutoff)  {
      f_n = -(area)*( F_n->Interpolate_Value(g) );
      modified = true;
    } else {
      f_n = 0.0;
      modified = true;
    }
  }  // else compression

  // determine tangential force (magnitude)
  if ( g < TCutoff ) {
      Real h[3];
      Scale(slip,1.0,h);  // incremental slip
      Real h_n = Dot(h,normal);
      Add(h,-h_n,normal,h); // incremental tangential slip
      Real mag_h = Magnitude(h);
      Real limit = (area)*F_t->Interpolate_Value(mag_h/dt);
      if ( mag_f_t > limit) {
        mag_f_t = limit;
        modified = true;
      } // else stick
  } else {
      mag_f_t = 0.0;
      modified = true;
  }

  if ( mag_t > ZERO_TOL) {
     Linear_Combination(f_n,normal,mag_f_t/mag_t,f_t,force);
  } else {
     Scale(normal,f_n,force);   
  }
  return modified;
}

bool
ContactTDJunction::Active_Interaction(ContactNodeEntityInteraction*,Real gap)
{
   return (gap < NCutoff) ? true : false;
}


