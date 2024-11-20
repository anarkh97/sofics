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


#ifndef ContactErrors_h_
#define ContactErrors_h_

#include "Contact_Defines.h"
#include "CString.h"

class ContactErrors{

 public:
  ContactErrors();
  ~ContactErrors();

  void Add_Error_Message(const char* message);

  inline int Number_of_Errors() {return number_of_errors;};
  const char* Error_Message( int i );

 private:
  
  ContactErrors(ContactErrors&);
  ContactErrors& operator=(ContactErrors&);
  
  int errors_array_size;
  int number_of_errors;
  CString** errors;

  void Reallocate();

};

#endif
