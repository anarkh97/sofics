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


#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactEnforcement.h"
#include "ContactTable.h"

void FORTRAN(UserQuery_Table_Last_Abscissa)(void *data, int *id, Real *value)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  ContactSearch*      search      = enforcement->Search();
  int num_tables = search->Num_Tables();
  ContactTable** tables = search->Tables();
  for( int i=0 ; i<num_tables ; ++i ){
    ContactTable* table = tables[i];
    if( *id == table->ID() ){
      *value = table->Last_Abscissa();
      break;
    }
  }
}
