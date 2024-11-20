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


#ifndef ContactEntityDataHash_
#define ContactEntityDataHash_

#include "ContactInteractionEntity.h"
#include "ContactHostGlobalID.h"
#include <iostream>
#include <vector>

class ContactEntityDataHash {
  
 public:
  ContactEntityDataHash( std::vector<ContactInteractionEntity<Real>::entity_data*> *list );
  ~ContactEntityDataHash();

  struct hash {
    ContactInteractionEntity<Real>::entity_data* entity;
    struct hash *next;
  };

  struct hash_list {
    struct hash* ptr;
    struct hash_list *next;
  };

  ContactInteractionEntity<Real>::entity_data* find( ContactInteractionEntity<Real>::entity_data*, int, hash*, 
                                                     ContactInteractionEntity<Real>::entity_data* );
  friend std::ostream& operator<<( std::ostream& os, const ContactEntityDataHash& hash );

 private:
  
  ContactEntityDataHash(ContactEntityDataHash&);
  ContactEntityDataHash& operator=(ContactEntityDataHash&);
  
  int nbins;
  int number_of_entities;

  hash* hash_space;
  hash** bins;
  hash_list* hash_space_list;

  void create_space();
  int hash_func(int);  
  int CompareEntityData(ContactInteractionEntity<Real>::entity_data*,
                        ContactInteractionEntity<Real>::entity_data*);
};

inline int
ContactEntityDataHash::hash_func(int val)
{
#if defined(JENKINS96_HASH_FUNC)
  unsigned int mask = nbins-1;
  unsigned int c    = (unsigned int)val;
  unsigned int a    = 0x9e3779b9;
  unsigned int b    = 0x7b6552e3;
  a=a-b;  a=a-c;  a=a^(c>>13);
  b=b-c;  b=b-a;  b=b^(a<<8);
  c=c-a;  c=c-b;  c=c^(b>>13);
  a=a-b;  a=a-c;  a=a^(c>>12);
  b=b-c;  b=b-a;  b=b^(a<<16);
  c=c-a;  c=c-b;  c=c^(b>>5);
  a=a-b;  a=a-c;  a=a^(c>>3);
  b=b-c;  b=b-a;  b=b^(a<<10);
  c=c-a;  c=c-b;  c=c^(b>>15);
  c=(c & mask);
  return (int)c;
#elif defined(JENKINS32_HASH_FUNC)
  unsigned int mask = nbins-1;
  unsigned int key  = (unsigned int)val;
  key += (key << 12);
  key ^= (key >> 22);
  key += (key <<  4);
  key ^= (key >>  9);
  key += (key << 10);
  key ^= (key >>  2);
  key += (key <<  7);
  key ^= (key >> 12);
  key  = (key & mask);
  return (int)key;
#elif defined(WANG32_HASH_FUNC)
  unsigned int mask = nbins-1;
  unsigned int key  = (unsigned int)val;
  key += ~(key << 15);
  key ^=  (key >> 10);
  key +=  (key <<  3);
  key ^=  (key >>  6);
  key += ~(key << 11);
  key ^=  (key >> 16);
  key  =  (key & mask);
  return (int)key;
#elif defined(SEDGEWICK_HASH_FUNC)
  return ((16161*val)%nbins);
#elif defined(KNUTH_HASH_FUNC)
  return ((2654435761U*(unsigned)val)%nbins);
#else
  return (val % nbins);
#endif
}

#endif // ContactEntityDataHash_
