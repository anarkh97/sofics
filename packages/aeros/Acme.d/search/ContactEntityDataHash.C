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


#include "ContactEntityDataHash.h"

ContactEntityDataHash::ContactEntityDataHash( std::vector<ContactInteractionEntity<Real>::entity_data*> *list)
{
  number_of_entities = list->size();
  hash_space = NULL;
  hash_space_list = NULL;

  // Create space for bins
  int nbins_orig = std::max(BIN_MINIMUM,(int)(BIN_FRACTION*number_of_entities));
  // Round up nbins to make it more prime-like
#if defined(PRIME_HASH_BINS)
  int s[] = {1,5,3,3,9,3,1,3,19,15,1,5,1,3,9,3,15,3,39,5,39,57,3,35,1};
  int cnt = 0;
  for (nbins=BIN_MINIMUM; nbins<nbins_orig; nbins *= 2, ++cnt);
  nbins -= s[cnt];
#elif defined(POWER2_HASH_BINS)
  // Round up nbins to make it a power of 2
  for (nbins=2; nbins<nbins_orig; nbins *= 2);
#else
  // Round up nbins to make it more prime-like
  nbins = nbins_orig;
  if ( !(nbins % 2) ) nbins += 1;
  if ( !(nbins % 3) ) nbins += 2;
  if ( !(nbins % 6) ) nbins += 6;
#endif

  // Allocate the bins
  bins = new hash*[nbins];

  // Initialize the bins to empty
  for( int i=0 ; i<nbins ; ++i )
    bins[i] = NULL;

  // Create hash space
  create_space();

  // Add each entity to the hash table
  for(int i = 0; i < list->size(); ++i) {
    ContactInteractionEntity<Real>::entity_data *entity = (*list)[i];
    find(entity, 1, NULL, entity);
  }
}

ContactEntityDataHash::~ContactEntityDataHash()
{
  hash_list* next_list;
  while( hash_space_list != NULL ){
    delete [] hash_space_list->ptr;
    next_list = hash_space_list->next;
    delete hash_space_list;
    hash_space_list = next_list;
  }
  delete [] bins;
}



void ContactEntityDataHash::create_space()
{
  int n = number_of_entities;
  if( n<=0 ) n=1;

  hash_space = new hash[n];

  hash_list* new_list = new hash_list;
  new_list->next = hash_space_list;
  new_list->ptr = hash_space;
  hash_space_list = new_list;
}


ContactInteractionEntity<Real>::entity_data* 
ContactEntityDataHash::find( ContactInteractionEntity<Real>::entity_data* data, 
                             int flag, hash* newpt, 
                             ContactInteractionEntity<Real>::entity_data* Entity )
{
  hash* ptr;
  hash** previous;
  int ibin;

  ibin = hash_func(data->index_in_owner_proc_array);
  
  previous = bins + ibin;
  ptr = bins[ibin];
  
  while( ptr != NULL && CompareEntityData(data, ptr->entity)<0 ){
    previous = &(ptr->next);
    ptr = ptr->next;
  }
  
  if( ptr == NULL || CompareEntityData(data, ptr->entity)>0 ){
    if( flag ){
      if( newpt == NULL ) newpt = hash_space++;
      *previous = newpt;
      newpt->entity = Entity;
      newpt->next = ptr;
    }
    return(NULL);
  }
  else if( flag == 2 ){
    ptr->entity = Entity;
  }
  return ptr->entity;
}

int
ContactEntityDataHash::CompareEntityData(ContactInteractionEntity<Real>::entity_data* e1,
                                         ContactInteractionEntity<Real>::entity_data* e2)
{
  if (e1->type < e2->type) {
    return -1;
  } else if (e1->type > e2->type) {
    return 1;
  }
  if (e1->owner < e2->owner) {
    return -1;
  } else if (e1->owner > e2->owner) {
    return 1;
  }
  if (e1->index_in_owner_proc_array < e2->index_in_owner_proc_array) {
    return -1;
  } else if (e1->index_in_owner_proc_array > e2->index_in_owner_proc_array) {
    return 1;
  }
  return 0;
}

std::ostream& operator<<( std::ostream& os, const ContactEntityDataHash& HASH )
{
  os << "Number of Entries = " << HASH.number_of_entities << std::endl;
  for( int i=0 ; i<HASH.nbins ; ++i ){
    os << "Outputting information for bin " << i << std::endl;
    struct ContactEntityDataHash::hash* ptr = HASH.bins[i];
    while( ptr != NULL ){
      ContactHostGlobalID global_id(ptr->entity->host_gid[0],
                                    ptr->entity->host_gid[1]);
      os << "Global ID = " << global_id << std::endl;
      ptr = ptr->next;
    }
  }
  return os;
}
