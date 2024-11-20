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


#ifndef ContactSequentialAllocator_h_
#define ContactSequentialAllocator_h_

#include "contact_assert.h"
#include "Contact_Defines.h"

class ContactSequentialAllocator {
public:
  
  static const int MEMORY_BOUNDARY_SIZE;
  
  ContactSequentialAllocator();
  ContactSequentialAllocator(int bytes_per_block,
                             int initial_multiplier);  // multiplies # of bytes
 ~ContactSequentialAllocator() { Purge_Mem(); }
  
  void  Resize(int bytes_per_block, int initial_multiplier);
  
  inline Real* New_Real_Frag(int num_words);
  inline void* New_Frag(int num_bytes);
  void  Reset();                  // Set pointer back to beginning of memory.
  void  Purge_Mem();              // Free up all memory to std::system.
  
  int   Num_Frags()          const { return num_frags;          }
  int   Block_Size()         const { return block_size;         }
  int   Initial_Multiplier() const { return initial_multiplier; }
  
  int   Num_Blocks() const;  // computed
  int   Bytes_Used() const;  // computed
  
 private:
  
  int num_frags;                // Current number of memory frags given out.
  int block_size;               // Number of bytes per allocated block.
  int initial_multiplier;       // Very first allocated block has size
                                //   block_size times this value.
  char* free_ptr;               // Next std::free memory location.
  int   bytes_free;             // Number of bytes std::free in current block.
  int   unused_bytes;           // Cummulative waste on ends of blocks.
  
  char* first_block;            // Pointer to very first block.
  
  
  // Internal methods.
  
  ContactSequentialAllocator(const ContactSequentialAllocator&);
  ContactSequentialAllocator& operator=(const ContactSequentialAllocator&);
  
  void  New_Block();
};

inline Real* ContactSequentialAllocator::New_Real_Frag( int num_words )
{ 
  return (Real*) New_Frag( num_words*sizeof(Real) );
}

inline void* ContactSequentialAllocator::New_Frag(int bytes_requested)
{
  PRECONDITION(num_frags >= 0);
  PRECONDITION(block_size > 0);
  PRECONDITION(!(block_size % MEMORY_BOUNDARY_SIZE));
  PRECONDITION(initial_multiplier > 0);
  PRECONDITION(bytes_free >= 0);
  PRECONDITION( !( Num_Blocks() == 0 &&
                   bytes_requested > block_size * initial_multiplier) );
  PRECONDITION( !( Num_Blocks() == 1 && bytes_requested > bytes_free
                                     && bytes_requested > block_size) );
  PRECONDITION( !( Num_Blocks() > 1 && bytes_requested > block_size) );
  
  if (bytes_requested <= 0) return 0;
  
  // This bumps the bytes requested to a multiple of MEMORY_BOUNDARY_SIZE.
  int m = bytes_requested % MEMORY_BOUNDARY_SIZE;
  if (m != 0) bytes_requested += MEMORY_BOUNDARY_SIZE - m;
  
  if (bytes_free < bytes_requested)
  {
    unused_bytes += bytes_free;
    New_Block();
  }
  
  char* p = free_ptr;
  free_ptr = &free_ptr[ bytes_requested ]; // Shift std::free pointer down.
  bytes_free -= bytes_requested;    // Reduce number of std::free bytes in the block.
  ++num_frags;
  
  return (void*)p;
}

#endif
