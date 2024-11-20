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


#include "ContactEnforcementData.h"
#include "contact_assert.h"
#include <cstring>

ContactEnforcementData::ContactEnforcementData( int Number_Entity_Keys,
						int Size_data_per_pair,
						const Real* Data )
  : number_entity_keys( Number_Entity_Keys ),
    size_data_per_pair( Size_data_per_pair )
{
  data = new Real[size_data_per_pair*number_entity_keys*number_entity_keys];
  if (Data) {
     std::memcpy( data, Data, size_data_per_pair*number_entity_keys*number_entity_keys*
   	     sizeof(Real) );
  }
}

ContactEnforcementData::ContactEnforcementData()
  : number_entity_keys( -1 ),
    size_data_per_pair( -1 )
{
  data = NULL;
}

int ContactEnforcementData::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = number_entity_keys;
  restart_data[words_added++] = size_data_per_pair;
  int size_data = number_entity_keys*number_entity_keys*size_data_per_pair;
  std::memcpy( &(restart_data[words_added]), data, size_data*sizeof(Real) );
  words_added += size_data;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactEnforcementData::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;

  // Get the size data
  number_entity_keys = (int) restart_data[0];
  size_data_per_pair = (int) restart_data[1];
  words_read +=2;

  int size_data = number_entity_keys*number_entity_keys*size_data_per_pair;
  data = new Real[size_data];
  std::memcpy( data, &(restart_data[words_read]), size_data*sizeof(Real) );
  words_read += size_data;
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactEnforcementData::~ContactEnforcementData()
{
  if( data ) delete [] data;
}
