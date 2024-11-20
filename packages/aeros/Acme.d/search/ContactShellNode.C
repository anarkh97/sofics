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


#include "ContactShellNode.h"
#include "ContactFixedSizeAllocator.h"
#include <cstddef>
#include <new>

ContactShellNode::ContactShellNode( ContactFixedSizeAllocator* alloc,
                                    ContactSearch::ContactNode_Type Type, 
				    int Block_Index, 
				    int Index_in_Block )
  : ContactNode<Real>( alloc, Type, Block_Index, Index_in_Block, 
		 CT_SHELL_NODE )
{
  previous_lofting[0] = 0.0;
  previous_lofting[1] = 0.0;
  previous_lofting[2] = 0.0;
}

ContactShellNode::~ContactShellNode()
{
}

bool ContactShellNode::ConnectedToAllShellFaces() {
  int num_faces = Number_Face_Connections();
  for(int iface = 0; iface < num_faces; ++iface) {
    if(!ContactSearch::Is_a_Shell_Face(GetFace(iface)->FaceType())) return false;
  }
  return true;
}


ContactShellNode* 
ContactShellNode::new_ContactShellNode( ContactFixedSizeAllocator* alloc,
					ContactSearch::ContactNode_Type Type,
					int Block_Index, 
					int Index_in_Block )
{
  return new (alloc[ContactSearch::ALLOC_ContactShellNode].New_Frag()) 
    ContactShellNode( alloc, Type, Block_Index, Index_in_Block );
}

void ContactShellNode_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellNode),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactShellNode allocator" );
}
