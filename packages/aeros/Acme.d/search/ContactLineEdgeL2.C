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


#ifndef ContactLineEdgeL2_C_
#define ContactLineEdgeL2_C_

#include "ContactLineEdgeL2.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include <new>

template<typename DataType>
ContactLineEdgeL2<DataType>::ContactLineEdgeL2( int Blk_Index, int Host_Index_in_Blk ) 
  : ContactEdge<DataType>( ContactSearch::LINEEDGEL2, Blk_Index, Host_Index_in_Blk, nodes )
{}


template<typename DataType>
ContactLineEdgeL2<DataType>* ContactLineEdgeL2<DataType>::new_ContactLineEdgeL2(
		    ContactFixedSizeAllocator& alloc,
		    int Block_Index, int Host_Index_in_Block )
{
  return new (alloc.New_Frag())
    ContactLineEdgeL2<DataType>( Block_Index, Host_Index_in_Block );
}

template<typename DataType>
void ContactLineEdgeL2_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactLineEdgeL2<DataType>),
		100,    // block size
		0,      // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactLineEdgeL2<DataType> allocator" );
}

template<typename DataType>
ContactLineEdgeL2<DataType>::~ContactLineEdgeL2() {}

#endif  // #define ContactLineEdgeL2_C_
