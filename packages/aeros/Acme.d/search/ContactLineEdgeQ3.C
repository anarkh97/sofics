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


#include "ContactLineEdgeQ3.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include <new>

ContactLineEdgeQ3::ContactLineEdgeQ3( int Blk_Index, int Index_in_Blk ) 
  : ContactEdge<Real>( ContactSearch::LINEEDGEQ3, Blk_Index, Index_in_Blk, nodes )
{}

ContactLineEdgeQ3* ContactLineEdgeQ3::new_ContactLineEdgeQ3(
		    ContactFixedSizeAllocator& alloc,
		    int Block_Index, int Index_in_Block )
{
  return new (alloc.New_Frag())
    ContactLineEdgeQ3( Block_Index, Index_in_Block );
}

void ContactLineEdgeQ3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactLineEdgeQ3),
		100,    // block size
		0 );  // initial block size
  alloc.Set_Name( "ContactLineEdgeQ3 allocator" );
}

ContactLineEdgeQ3::~ContactLineEdgeQ3() {}
