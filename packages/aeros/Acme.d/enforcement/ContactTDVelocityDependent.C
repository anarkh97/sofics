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


#include "ContactTDVelocityDependent.h"
#include "contact_assert.h"

#include "ContactTDEnforcement.h"
#include "ContactNodeEntityInteraction.h"

ContactTDVelocityDependent::ContactTDVelocityDependent( int ID, int* /*int_data*/,
						      Real* real_data,
						      ContactTopology* topology)
  : ContactTDEnfModel( ID, ContactEnforcement::TD_VELOCITY_DEPENDENT, topology )
{
  mu_static  = real_data[0]; 
  mu_dynamic = real_data[1];
  vel_decay  = real_data[2];
}


ContactTDVelocityDependent::ContactTDVelocityDependent(ContactTopology* topology)
  : ContactTDEnfModel(ContactEnforcement::TD_VELOCITY_DEPENDENT, topology)
{
  mu_static  = 1.0; 
  mu_dynamic = 1.0;
  vel_decay  = 1.0;
}

ContactTDVelocityDependent::~ContactTDVelocityDependent(){}

int ContactTDVelocityDependent::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = mu_static;
  restart_data[words_added++] = mu_dynamic;
  restart_data[words_added++] = vel_decay;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDVelocityDependent::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id  = (int) restart_data[words_read++];
  mu_static  = restart_data[words_read++];
  mu_dynamic = restart_data[words_read++];
  vel_decay  = restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDVelocityDependent::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDVelocityDependent::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

bool ContactTDVelocityDependent::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{ 
  bool modified = false;
  Scale(force,1/area,force);
  Real f_n = Dot(force,normal);
  if (f_n > 0.0) {
    Real f_t[3];
    Add(force,-f_n,normal,f_t);
    Real mag_t = Magnitude(f_t);
    Real slip_rate = Magnitude(slip)/dt;
    Real mu = (std::fabs(-vel_decay*slip_rate) > ZERO_TOL) ?
       (mu_static - mu_dynamic) *std::exp(-vel_decay*slip_rate) + mu_dynamic:
       mu_static;
    mu = std::max(mu,mu_static); 
    Real limit = std::max ( 0.0, mu*f_n   ) ; // limit >= 0.0
    if ( limit < mag_t && mag_t > ZERO_TOL ) {
       Linear_Combination(f_n,normal,limit/mag_t,f_t,force);
       modified = true;
    }
    Scale(force,area,force);
    return modified;
  }
  else {
    Zero(force);
    return true;
  }
}



