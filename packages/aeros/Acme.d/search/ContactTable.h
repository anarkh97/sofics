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


#ifndef ContactTable_h_
#define ContactTable_h_

#include "Contact_Defines.h"

class ContactTable {
  
 public:
  
  ContactTable( int ID, int Num_Points, Real* abscissa, Real* ordinate );
  ContactTable() : id(-1),num_points(0),abscissa(NULL),ordinate(NULL) {};
  ~ContactTable();
  
  Real Interpolate_Value( Real abscissa );
  int ID() { return id; };
  int Restart_Size() { return 2+2*num_points; };
  int Extract_Restart_Data( Real* restart_data );
  int Implant_Restart_Data( Real* restart_data );
  inline Real Last_Ordinate()  {return ordinate[num_points-1];}
  inline Real First_Ordinate() {return ordinate[0];}
  inline Real Last_Abscissa()  {return abscissa[num_points-1];}
  inline Real First_Abscissa() {return abscissa[0];}
  inline int Num_Points() {return num_points;};
  inline Real* Abscissas() {return abscissa;};
  inline Real* Ordinates() {return ordinate;};

 private:
  
  ContactTable(ContactTable&);
  ContactTable& operator=(ContactTable&);
  
  int id;
  int num_points;
  Real* abscissa;
  Real* ordinate;

};


#endif
