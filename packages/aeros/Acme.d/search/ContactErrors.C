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


#include "ContactErrors.h"
#include <cstring>

ContactErrors::ContactErrors()
{
  number_of_errors = 0;
  errors_array_size = 10;
  errors = new CString*[errors_array_size];
  std::memset( errors, 0, errors_array_size*sizeof(CString*) );
}

ContactErrors::~ContactErrors()
{
  for( int i=0 ; i<number_of_errors ; ++i )
    delete errors[i];
  delete [] errors;
}


const char* ContactErrors::Error_Message( int i )
{
  return errors[i]->data();
}

void ContactErrors::Add_Error_Message( const char* Message )
{
  if( number_of_errors == errors_array_size ) Reallocate();
  // The last error message wasn't deleted in Error_Message so delete it
  // here to avoid a memory lead (if needed)
  errors[number_of_errors++] = new CString(Message);
}

void ContactErrors::Reallocate()
{
  CString** temp = errors;
  errors_array_size += 10;
  errors = new CString*[errors_array_size];
  std::memset( errors, 0, errors_array_size*sizeof(CString*) );
  for( int i=0 ; i<number_of_errors ; ++i )
    errors[i] = temp[i];
  delete [] temp;
}
