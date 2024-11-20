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


#ifndef ContactInteractionEntity_C_
#define ContactInteractionEntity_C_

#include "ContactInteractionEntity.h"
#include "ContactParOStream.h"
#include "contact_assert.h"
#include <cstddef>
#include <cstring>

template<typename DataType>
ContactInteractionEntity<DataType>::ContactInteractionEntity(DataType *data_array_, ContactType base_type_)
  : data_array(data_array_), base_type(base_type_), index(-1), proc_index(-1)
{}

template<typename DataType>
ContactInteractionEntity<DataType>::~ContactInteractionEntity()
{}

template<typename DataType>
void 
ContactInteractionEntity<DataType>::Display(ContactParOStream& postream)
{
  postream<<"ContactEntity: \n";
  postream<<"               entity type:       ";
  switch (Base_Type()) {
  case CT_NNI:
    postream<<"CT_NNI\n";
    break;
  case CT_NFI:
    postream<<"CT_NFI\n";
    break;
  case CT_NSI:
    postream<<"CT_NSI\n";
    break;
  case CT_FFI:
    postream<<"CT_FFI\n";
    break;
  case CT_FCI:
    postream<<"CT_FCI\n";
    break;
  case CT_EEI:
    postream<<"CT_EEI\n";
    break;
  case CT_UNKNOWN:
    postream<<"CT_UNKNOWN\n";
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  postream<<"               index:             "<<index<<"\n";
  postream<<"               proc_index:        "<<proc_index<<"\n";
}

#endif  // ifdef ContactInteractionEntity_C_
