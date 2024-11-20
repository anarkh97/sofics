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


#ifndef ContactShellTriFaceL3_C_
#define ContactShellTriFaceL3_C_

#include "ContactShellTriFaceL3.h"
#include "ContactFixedSizeAllocator.h"
#include <new>

template<typename DataType>
ContactShellTriFaceL3<DataType>::ContactShellTriFaceL3(ContactFixedSizeAllocator* alloc,
                                             int Block_Index, 
					     int Index_in_Block, int key ) 
  : ContactTriFaceL3<DataType>( alloc, Block_Index, Index_in_Block, key )
{
  // Reset the type (because TriFaceL3 set it to its type)
  this->face_type = ContactSearch::SHELLTRIFACEL3;
}

template<typename DataType>
ContactShellTriFaceL3<DataType>* ContactShellTriFaceL3<DataType>::new_ContactShellTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactShellTriFaceL3].New_Frag())
             ContactShellTriFaceL3<DataType>(alloc, Block_Index, Index_in_Block, key);
}

template<typename DataType>
void ContactShellTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellTriFaceL3<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactShellTriFaceL3 allocator" );
}

template<typename DataType>
ContactShellTriFaceL3<DataType>::~ContactShellTriFaceL3() 
{
}

#endif // ContactShellTriFaceL3_C_
