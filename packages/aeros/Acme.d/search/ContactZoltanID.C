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


#ifndef CONTACT_NO_MPI

#include "ContactZoltanID.h"
#include "ContactHostGlobalID.h"
#include "Contact_Defines.h"

// Classes to implement Zoltan local and global IDs.  
//
// An array of two unsigned ints is used to describe the local ID. 
// 
//    lid[0] = entity_type
//    lid[1] = index on processer
//
//An array of three unsigned ints is used to describe the global ID.
// 
//    gid[0] = entity_type
//    gid[1] = host_id[0]
//    gid[1] = host_id[1]


ContactZoltanLID::ContactZoltanLID() : 
  type(-1), index(-1)
{}

ContactZoltanLID::~ContactZoltanLID(){}

//=======================================================================

ContactZoltanGID::ContactZoltanGID() : 
  type(-1), hi(-1), lo(-1)
{}

ContactZoltanGID::~ContactZoltanGID(){}

void 
ContactZoltanGID::ZoltanGID(int id_type, ContactHostGlobalID* gid, 
                            LB_ID_PTR id)
{
  id[0] = id_type;
  id[1] = gid->HiInt();
  id[2] = gid->LoInt();
}

#endif
