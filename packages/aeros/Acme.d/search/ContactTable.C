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


#include "ContactTable.h"
#include "contact_assert.h"
#include <cstring>

ContactTable::ContactTable( int ID, int Num_Points, Real* Abscissa, 
			    Real* Ordinate )
{
  PRECONDITION( ID > 0 );
  PRECONDITION( Num_Points > 0 );
  PRECONDITION( Abscissa );
  PRECONDITION( Ordinate );
  id = ID;
  num_points = Num_Points;
  abscissa = new Real[num_points];
  ordinate = new Real[num_points];
  std::memcpy( abscissa, Abscissa, num_points*sizeof(Real) );
  std::memcpy( ordinate, Ordinate, num_points*sizeof(Real) );
}


ContactTable::~ContactTable()
{
  delete [] abscissa;
  delete [] ordinate;
}

Real ContactTable::Interpolate_Value( Real Abscissa )
{
  for( int i=0 ; i<num_points-1 ; ++i){
    if( Abscissa >= abscissa[i] && Abscissa <= abscissa[i+1] ){
      // linearly interpolate
      return( ordinate[i] +
	      ((Abscissa      - abscissa[i] )/
	       (abscissa[i+1] - abscissa[i] )) *
	       (ordinate[i+1] - ordinate[i]) );
    }
  }

  // Out of the bounds.  Return the value at the bound
  if( Abscissa < abscissa[0] )
    return ordinate[0];
  else
    return ordinate[num_points-1];
}


int ContactTable::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = num_points;
  std::memcpy( &restart_data[words_added], abscissa, num_points*sizeof(Real) );
  words_added += num_points;
  std::memcpy( &restart_data[words_added], ordinate, num_points*sizeof(Real) );
  words_added += num_points;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTable::Implant_Restart_Data( Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  num_points = (int) restart_data[words_read++];
  abscissa = new Real[num_points];
  ordinate = new Real[num_points];
  std::memcpy( abscissa, &restart_data[words_read], num_points*sizeof(Real) );
  words_read += num_points;
  std::memcpy( ordinate, &restart_data[words_read], num_points*sizeof(Real) );
  words_read += num_points;
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}
