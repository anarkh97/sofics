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
#include "ContactTDFrictionless.h"
#include "contact_assert.h"
#if (CONTACT_DEBUG_PRINT_LEVEL>=2)
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif

ContactTDFrictionless::ContactTDFrictionless( int ID, int* /*int_data*/, 
					      Real* real_data,
					      ContactTopology* topo )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_FRICTIONLESS, topo )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding frictionless enforcement model, ID ="<<ID<<std::endl<<std::flush;
  }
#endif
}

ContactTDFrictionless::ContactTDFrictionless(ContactTopology* topo)
  : ContactTDEnfModel(ContactEnforcement::TD_FRICTIONLESS, topo)
{
}

int ContactTDFrictionless::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDFrictionless::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDFrictionless::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDFrictionless::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactTDFrictionless::~ContactTDFrictionless(){}

bool ContactTDFrictionless::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{ 
  Real f_n = Dot(force,normal);
  if (f_n > 0) {   // NOTE: THIS IS BACKWARD (REJ 4/11)
    Scale(normal,f_n,force);
    return false;
  }
  else {
    Zero(force);
    return true;
  }
}
