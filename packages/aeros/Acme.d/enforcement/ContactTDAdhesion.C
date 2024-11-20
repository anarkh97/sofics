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


#include "ContactTDAdhesion.h"
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

ContactTDAdhesion::ContactTDAdhesion( int ID, int* int_data, 
					      Real* real_data,
					      ContactTopology* topo )
  : ContactTDEnfModel( ID, ContactEnforcement::TD_ADHESION, topo )
{
   table_id = int_data[0];
   F_n = NULL; // force-displ law
   Cutoff = 0.0; // adhesion cutoff 
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding adhesion enforcement model, ID ="<<ID<<std::endl;
    std::cout<<"  table id = "<<table_id<<std::endl;
    std::cout<<"  cutoff   = "<<Cutoff<<std::endl<<std::flush;
  }
#endif
}

ContactTDAdhesion::ContactTDAdhesion(ContactTopology* topo)
  : ContactTDEnfModel(ContactEnforcement::TD_ADHESION, topo)
{
   table_id = 0;
   F_n = NULL; // force-displ law
   Cutoff = 0.0; // adhesion cutoff
}

ContactSearch::ContactErrorCode
ContactTDAdhesion::Initialize_Model( int /* num_models */, 
                                    ContactEnfModel** /* models */,
                                      int  num_tables ,
                                      ContactTable**  tables )
{
  // initialize normal force-displ. law and cutoff from the table
  for( int i=0 ; i<num_tables ; ++i){
    if( table_id == tables[i]->ID() ){
      F_n = tables[i];
      break;
    }
  }
  if( !F_n )
    return ContactSearch::INVALID_ID;
  else
    Cutoff = F_n->Last_Abscissa();
#if LOCAL_PRINT_FLAG > 0
    std::cout << " TDAdhesion, Cutoff " << Cutoff << "\n";
#endif
  return ContactSearch::NO_ERROR;
}


int ContactTDAdhesion::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = table_id;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDAdhesion::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  table_id = (int) restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDAdhesion::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model 
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDAdhesion::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model 
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactTDAdhesion::~ContactTDAdhesion(){}

bool ContactTDAdhesion::Limit_Force(ContactNodeEntityInteraction* cnfi,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{ 
#if LOCAL_PRINT_FLAG>=2
std::cout << "\nArea:"<<area<<"\n";
//  std::cout << "gap " << gap << ", rel_disp.n " << Dot(rel_disp,normal) << ", slip.n " << Dot(slip,normal) << "\n";

  Real force_n = Dot(force,normal);
  Real force_t[3];
  Remove_Component(force,normal,force_t);
#endif

  Real g = Dot(rel_disp,normal); // since the gap read in has been manipulated
  Real f_n = Dot(force,normal);
  if (f_n >= 0)  {  // compression
#if LOCAL_PRINT_FLAG>=2
std::cout << "AD in Compression " << f_n << " gap " <<  gap << ", Cutoff " << Cutoff << "\n" ;
#endif
    return false;
  }
  else { // tension
    if (g < Cutoff)  {
      f_n = -F_n->Interpolate_Value(g);
      Scale(normal,(f_n*area),force);
#if LOCAL_PRINT_FLAG>=2
std::cout << "AD in Tension " << force_n << " -> " << -F_n->Interpolate_Value(g) << " with gap " << g << " < " << Cutoff << "\n";
if (Dot(force,normal)*force_n < 0) std::cout << "WARNING : normal force flipped\n";
#endif
      return true;
    } else {
#if LOCAL_PRINT_FLAG>=2
    std::cout << "AD ZEROd " << f_n << " gap " <<  g << " > " << Cutoff << "\n" ;
#endif
      Zero(force);
      return true;
    }
  }
}

bool
ContactTDAdhesion::Active_Interaction(ContactNodeEntityInteraction*,Real gap)
{
   return (gap < Cutoff) ? true : false;
}

