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


#include "ContactEnforcement.h"
#include "ContactTDCohesiveZone.h"
#include "contact_assert.h"
#undef LOCAL_PRINT_FLAG 
#define LOCAL_PRINT_FLAG 0
#if (CONTACT_DEBUG_PRINT_LEVEL>=2) || (LOCAL_PRINT_FLAG>0)
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif

ContactTDCohesiveZone::ContactTDCohesiveZone( int ID, int* int_data, 
					      Real* real_data,
					      ContactTopology* topo )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_COHESIVE_ZONE, topo )
{
   table_id = int_data[0];
   F = NULL; // force-displ law
   g_n_crit = real_data[0]; // normal cutoff distance
   g_t_crit = real_data[1]; // tangential cutoff distance
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding cohesive zone enforcement model, ID ="<<ID<<std::endl;
    std::cout<<"  table_id = "<<table_id<<std::endl;
    std::cout<<"  g_n_crit = "<<g_n_crit<<std::endl;
    std::cout<<"  g_t_crit = "<<g_t_crit<<std::endl;
    std::cout<<std::flush;
  }
#endif
}

ContactTDCohesiveZone::ContactTDCohesiveZone(ContactTopology* topo)
  : ContactTDEnfModel(ContactEnforcement::TD_COHESIVE_ZONE, topo)
{
   table_id = 0;
   F = NULL; // force-displ law
   g_n_crit = 1.0; // normal cutoff distance
   g_t_crit = 1.0; // tangential cutoff distance
}

ContactSearch::ContactErrorCode
ContactTDCohesiveZone::Initialize_Model( int /* num_models */,
                                    ContactEnfModel** /* models */,
                                      int  num_tables ,
                                      ContactTable**  tables )
{
  for( int i=0 ; i<num_tables ; ++i){
    if( table_id == tables[i]->ID() ){
      F = tables[i];
      return ContactSearch::NO_ERROR;
    }
  }
  return ContactSearch::INVALID_ID;
}


int ContactTDCohesiveZone::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = table_id;
  restart_data[words_added++] = g_n_crit;
  restart_data[words_added++] = g_t_crit;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDCohesiveZone::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  table_id = (int) restart_data[words_read++];
  g_n_crit = restart_data[words_read++];
  g_t_crit = restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDCohesiveZone::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model 
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDCohesiveZone::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model 
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactTDCohesiveZone::~ContactTDCohesiveZone(){}

bool ContactTDCohesiveZone::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{ 
// the force computed by the tied assumption is not a good predictor
// for separation
#if LOCAL_PRINT_FLAG>=2
  Real force_n = Dot(force,normal);
  Real force_t[3];
  Remove_Component(force,normal,force_t);
  Real mag_force_t = Magnitude(force_t);
#endif
  Scale(force,1/area,force);
  Real f_n = Dot(force,normal);
  Real s_n = Dot(rel_disp,normal);
  if (s_n < g_n_crit) { //using norm rel_disp instead of gap
    POSTCONDITION(g_n_crit > 0);
    Real d_n = s_n/g_n_crit;  // non-dimensionalize
    Real s_t[3];
    Remove_Component(rel_disp,normal,s_t);
    Real mag_s_t = Magnitude(s_t);
    Real d_t = ( g_t_crit > ZERO_TOL ) ? mag_s_t/g_t_crit : 0.0; // non-dim
    Real lambda = std::sqrt(d_n*d_n + d_t*d_t);
    Real sigma = std::fabs(F->Interpolate_Value(lambda)); // assuming positive values
    // take constraint value instead of a L'Hopital analysis at lambda = 0
    if  ((s_n > 0.0) && (lambda > ZERO_TOL))
      f_n = -(sigma/lambda * d_n); // replace tension value
    Real mag_f_t = 0.0;
    if (g_t_crit > ZERO_TOL && mag_s_t > ZERO_TOL*g_t_crit ) {
      mag_f_t = -(sigma/lambda * g_n_crit/g_t_crit * d_t);
      Linear_Combination(f_n,normal,mag_f_t/mag_s_t,s_t,force);
    } else {
      Scale(normal,f_n,force);
    }
    Scale(force,area,force);
#if LOCAL_PRINT_FLAG>=2
  std::cout << "f_n- " << force_n << ", f_n+ " << f_n << "\n";
  std::cout << "s_n "<<s_n << ", sigma " << sigma << ", lambda " << lambda << "\n";
  std::cout << "|s_t| "<< mag_s_t 
	  << ", |f_t|- " << mag_force_t << ", |f_t| " << mag_f_t << "\n";
  if (f_n*s_n > ZERO_TOL && std::fabs(s_n) > ZERO_TOL) 
	  std::cout << "normal force and disp not opposed \n";
#endif
    return true;
  }
  else {
    Zero(force);
#if LOCAL_PRINT_FLAG>=2
    std::cout << " Force ZEROd " << force[0] << force[1] << force[2] << "\n";
#endif
    return true;
  }
}

bool 
ContactTDCohesiveZone::Active_Interaction(ContactNodeEntityInteraction*,Real gap)
{
   return (gap < g_n_crit ? true : false);
}
