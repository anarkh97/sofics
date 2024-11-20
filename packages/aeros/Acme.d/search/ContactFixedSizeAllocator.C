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


#include <cstring>

#include "ContactFixedSizeAllocator.h"
#include <iostream>

using namespace std;

ContactFixedSizeAllocator::ContactFixedSizeAllocator()
  : object_name  (0),
    storage_bytes(0),
    block_size   (0),
    initial_size (0),
    num_obj      (0),
    num_free     (0)
{ }

ContactFixedSizeAllocator::ContactFixedSizeAllocator(
         int ssize, int blen, int ilen, const char* name, size_t rlen)
  : object_name  (0),
    storage_bytes(0),
    block_size   (0),
    initial_size (0),
    num_obj      (0),
    num_free     (0)
{
  PRECONDITION(ssize >= 0);
  PRECONDITION(blen > 1);
  PRECONDITION(ilen >= 0 && ilen != 1);
  
  if (ilen == 0) ilen = blen;  // Initial length defaults to block size.
  if (rlen == 0) rlen = sizeof(Real);

  // FSM links unused frags by using the memory allocated to the frags, so each
  // frag must be large enough to hold the link pointer when the frag is free.
  storage_bytes = ssize >= sizeof(char*) ? ssize : sizeof(char*);
  
  // In order to align the free frag link pointers on 4-byte boundaries,
  // the storage bytes must be a multiple of 4.

  int max_alignment = std::max(sizeof(char*), rlen);
  while ( storage_bytes % max_alignment ) ++storage_bytes;
  
  num_obj = blen;

  block_size   = blen * storage_bytes;
  initial_size = ilen * storage_bytes;
  
  Set_Name(name);
}

ContactFixedSizeAllocator::~ContactFixedSizeAllocator()
{
  PRECONDITION(Sanity_Check());
  
  Purge();
  
  if (object_name != 0) delete [] object_name;
}

void ContactFixedSizeAllocator::Resize(int ssize, int blen, int ilen, size_t rlen)
{
  PRECONDITION(Sanity_Check());
  PRECONDITION(ssize >= 0);
  PRECONDITION(blen > 1);
  PRECONDITION(ilen >= 0 && ilen != 1);
  
  if (ilen == 0) ilen = blen;  // Initial length defaults to block size.
  if (rlen == 0) rlen = sizeof(Real);
  
  num_obj = blen;

  Purge();
 
  storage_bytes = ssize >= sizeof(char*) ? ssize : sizeof(char*);
  int max_alignment = std::max(sizeof(char*), rlen);
  while ( storage_bytes % max_alignment ) ++storage_bytes;
  block_size   = blen * storage_bytes;
  initial_size = ilen * storage_bytes;
}

void ContactFixedSizeAllocator::Set_Name(const char* name)
{
  PRECONDITION(name == 0 || std::strlen(name) > 0);
  
  if (object_name != 0) {
    delete [] object_name;
    object_name = 0;
  }
  
  if (name != 0) {
    object_name = new char[std::strlen(name) + 1];
    std::strcpy(object_name, name);
  }
}

void ContactFixedSizeAllocator::Purge()
{
  PRECONDITION(Sanity_Check());
  
// Can't purge if memory was given out by global new/delete operators.
#ifndef CONTACTALLOCATOR_BYPASS
  
  for(int i = 0; i < mem_blocks.size(); ++i) {
    delete [] mem_blocks[i];
  }
  mem_blocks.clear();
  free_spots.clear();
  num_free = 0;
#endif
}


int ContactFixedSizeAllocator::Sanity_Check() const
{
  // Sanity checks.
  PRECONDITION(storage_bytes >= 0);
  PRECONDITION(block_size   >= 0);
  PRECONDITION(initial_size >= 0);
  
  PRECONDITION( !( storage_bytes > 0 && !block_size ) );
  PRECONDITION( !( storage_bytes > 0 && !initial_size ) );

  PRECONDITION( !( storage_bytes > 0 && storage_bytes < sizeof(char*)) );
  PRECONDITION( !( storage_bytes % sizeof(char*)) );
  PRECONDITION( !( storage_bytes > 0 && block_size % storage_bytes) );
  PRECONDITION( !( storage_bytes > 0 && initial_size % storage_bytes) );
  
  PRECONDITION( !( storage_bytes > 0 && block_size == storage_bytes) );
  PRECONDITION( !( storage_bytes > 0 && initial_size == storage_bytes) );
  
  return 1;
}
