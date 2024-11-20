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


#include "ContactTDConstantFriction.h"
#include "contact_assert.h"

#include "ContactTDEnforcement.h"
#include "ContactNodeEntityInteraction.h"

#define LOCAL_PRINT_FLAG 0
#if (CONTACT_DEBUG_PRINT_LEVEL>=2) || (LOCAL_PRINT_FLAG>0)
#include <algorithm>
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif


ContactTDConstantFriction::ContactTDConstantFriction( int ID, int* /*int_data*/,
						      Real* real_data,
						      ContactTopology* topo)
  : ContactTDEnfModel( ID, ContactEnforcement::TD_CONSTANT_FRICTION, topo )
{
  mu = real_data[0];
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding contant friction enforcement model, ID ="<<ID<<std::endl;
    std::cout<<"  mu = "<<mu<<std::endl<<std::flush;
  }
#endif
}


ContactTDConstantFriction::ContactTDConstantFriction(ContactTopology* topo)
  : ContactTDEnfModel(ContactEnforcement::TD_CONSTANT_FRICTION, topo)
{
  mu = -1.0;
}

ContactTDConstantFriction::~ContactTDConstantFriction(){}

int ContactTDConstantFriction::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = mu;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDConstantFriction::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  mu = restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDConstantFriction::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDConstantFriction::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

bool ContactTDConstantFriction::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{ 
  bool modified = false;
  Scale(force,1/area,force);
  Real f_n = Dot(force,normal);
  if (f_n > 0.0) {
    Real f_t[3];
    Add(force,-f_n,normal,f_t);
    Real mag_t = Magnitude(f_t);
    Real limit = std::max( 0.0, mu*f_n   ) ; // limit >= 0.0
    if ( mag_t > limit && mag_t > ZERO_TOL ) {
       Linear_Combination(f_n,normal,limit/mag_t,f_t,force);
       modified = true;
#if LOCAL_PRINT_FLAG > 0
       if (Dot(f_t,slip) > 0) std::cout << " WARNING : force and slip not opposed "
	       << Dot(f_t,slip) << "\n"
	       << " f_t " << f_t[0] << ", " << f_t[1] << ", " << f_t[2] << "\n" 
	       << " slp " << slip[0] << ", " << slip[1] << ", " << slip[2] << "\n";
#endif
    }
    Scale(force,area,force);
    return modified;
  }
  else {
    Zero(force);
    return true;
  }
}
