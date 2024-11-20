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


#ifndef ContactEnforcementData_h_
#define ContactEnforcementData_h_

#include "Contact_Defines.h"
#include "contact_assert.h"

class ContactTopology;

class ContactEnforcementData {

 public:
  ContactEnforcementData( int number_entity_keys, int size_data_per_pair,
			  const Real* Data );
  ContactEnforcementData();
  ~ContactEnforcementData();

  Real* Data_Array() { return data; };
  Real Get_Data( int index, int f_key, int n_key  )
    { PRECONDITION (f_key >= 0);
      PRECONDITION (n_key >= 0);
      PRECONDITION (f_key < number_entity_keys);
      PRECONDITION (n_key < number_entity_keys);
      PRECONDITION (index < size_data_per_pair );
      return data[(f_key*number_entity_keys+n_key)*size_data_per_pair+index]; 
    };
  void Set_Data( int index, int f_key, int n_key , Real value )
    { PRECONDITION (f_key >= 0);
      PRECONDITION (n_key >= 0);
      PRECONDITION (f_key < number_entity_keys);
      PRECONDITION (n_key < number_entity_keys);
      PRECONDITION (index < size_data_per_pair );
      data[(f_key*number_entity_keys+n_key)*size_data_per_pair+index]=value; 
    };
  int Number_Entity_Keys() { return number_entity_keys; };
  int Size_Data_Per_Pair() { return size_data_per_pair; };
  int Restart_Size() { return 2+number_entity_keys*number_entity_keys*
			      size_data_per_pair; };
  int Extract_Restart_Data( Real* restart_data );
  int Implant_Restart_Data( const Real* restart_data );

 private:
  
  ContactEnforcementData(ContactEnforcementData&);
  ContactEnforcementData& operator=(ContactEnforcementData&);
  
  int number_entity_keys;
  int size_data_per_pair;
  Real* data;
  
};

#endif
