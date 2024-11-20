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


#include "ContactSequentialAllocator.h"

const int ContactSequentialAllocator::MEMORY_BOUNDARY_SIZE = sizeof(double);

ContactSequentialAllocator::ContactSequentialAllocator()
  : num_frags(0), block_size(0), initial_multiplier(0),
    free_ptr(0), bytes_free(0), unused_bytes(0), first_block(0)
{ }

ContactSequentialAllocator::ContactSequentialAllocator(int bsize, int mult)
  : num_frags(0), block_size(0), initial_multiplier(0),
    free_ptr(0), bytes_free(0), unused_bytes(0), first_block(0)
{
  PRECONDITION(bsize > 0);
  int m = bsize % MEMORY_BOUNDARY_SIZE;
  if (m != 0)
    block_size = bsize + MEMORY_BOUNDARY_SIZE - m;
  else
    block_size = bsize;
  
  PRECONDITION(mult > 0);
  initial_multiplier = mult;
}

void ContactSequentialAllocator::Resize(int bsize, int mult)
{
  // Check state.
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size >= 0 && !(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier >= 0);
  PRECONDITION(bytes_free >= 0);

  Purge_Mem();
  
  PRECONDITION(bsize > 0);
  int m = bsize % MEMORY_BOUNDARY_SIZE;
  if (m != 0)
    block_size = bsize + MEMORY_BOUNDARY_SIZE - m;
  else
    block_size = bsize;
  
  PRECONDITION(mult > 0);
  initial_multiplier = mult;

  free_ptr     = 0;
  first_block  = 0;
  bytes_free   = 0;
  num_frags    = 0;
  unused_bytes = 0;
}

void ContactSequentialAllocator::New_Block()
{
  // Check state.
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size > 0 && !(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier > 0);
  PRECONDITION(bytes_free >= 0);

  if (free_ptr)  // Is there at least one block already allocated?
  {
    char* next_block = *( (char**)( &free_ptr[ bytes_free ] ) );
    if (!next_block)
    {
      next_block = new char[ block_size + sizeof(char*) ];  // Get new block.
      PRECONDITION(next_block);
      
      *( (char**)( &free_ptr[ bytes_free ] ) ) = next_block;  // Link new block.
      *( (char**)( &next_block[ block_size ] ) ) = 0;  // Set next pntr to null.
    }
    free_ptr   = next_block;
    bytes_free = block_size;
  }
  else  // No blocks currently in memory.
  {
    bytes_free  = block_size * initial_multiplier;  // First block has multiple.
    first_block = free_ptr = new char[ bytes_free + sizeof(char*) ];
    PRECONDITION(first_block);
    *( (char**)( &first_block[ bytes_free ] ) ) = 0;  // Set next pntr to null.
  }
}

void ContactSequentialAllocator::Reset()
{
  // Check state.
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size >= 0 && !(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier >= 0);
  PRECONDITION(bytes_free >= 0);
  PRECONDITION( !( first_block && (block_size == 0 || initial_multiplier == 0) ) );

  if (first_block) {  // Has an initial block been allocated?
    bytes_free = block_size * initial_multiplier;
    free_ptr   = first_block;
  }
  else {
    bytes_free = 0;
    free_ptr   = 0;
  }
  
  num_frags    = 0;
  unused_bytes = 0;
}

void ContactSequentialAllocator::Purge_Mem()
{
  // Check state.
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size >= 0 && !(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier >= 0);
  PRECONDITION(bytes_free >= 0);
  
  char* tmp = first_block;
  if (tmp)
  {
    // Set first_block to next block pointer.
    first_block = *( (char**)( &tmp[ block_size * initial_multiplier ] ) );
    delete [] tmp;
    tmp = first_block;
    
    while (tmp) {
      first_block = *( (char**)( &tmp[ block_size ] ) );  // Next block pointer.
      delete [] tmp;
      tmp = first_block;
    }
  }
  free_ptr     = 0;
  num_frags    = 0;
  bytes_free   = 0;
  unused_bytes = 0;
}

int ContactSequentialAllocator::Num_Blocks() const
{
  // Check state.
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size >= 0 && !(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier >= 0);
  PRECONDITION(bytes_free >= 0);
  
  int blocks = 0;
  char* tmp = first_block;
  if (tmp)
  {
    ++blocks;
    tmp = *( (char**)( &tmp[ block_size * initial_multiplier ] ) );
    
    while (tmp) {
      ++blocks;
      tmp = *( (char**)( &tmp[ block_size ] ) );
    }
  }
  PRECONDITION(blocks >= 0);
  return blocks;
}

int ContactSequentialAllocator::Bytes_Used() const
{
  // Check state.
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size >= 0 && !(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier >= 0);
  PRECONDITION(bytes_free >= 0);
  
  int blocks = Num_Blocks();
  if (blocks)
    return block_size * initial_multiplier + (blocks - 1) * block_size;
  else
    return 0;
}
