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


#include <algorithm>

#include "ContactTopologyEntityHash.h"

using namespace std;

ContactTopologyEntityHash::ContactTopologyEntityHash() {
}

ContactTopologyEntityHash::ContactTopologyEntityHash(int n, ContactTopologyEntity<Real>** e) : implementation(n, e)
{
}

ContactTopologyEntityHash::ContactTopologyEntityHash(std::vector<ContactTopologyEntity<Real>*> *v) : implementation(v)
{
}

ContactTopologyEntityHash::~ContactTopologyEntityHash() {
}

std::ostream& operator<<( std::ostream& os, const ContactTopologyEntityHash& HASH ) {
  return HASH.stream_data(os);
}

ContactParOStream& operator<<( ContactParOStream& os, const ContactTopologyEntityHash& HASH ) {
  return HASH.stream_data(os);
}

ContactParOStream &ContactTopologyEntityHash::stream_data(ContactParOStream &os) const {
  return stream_data(os);
}

std::ostream &ContactTopologyEntityHash::stream_data(std::ostream &os) const {
  return stream_data(os);
}


namespace topology_hash_0 {

implementation::implementation() : nbins(0), nbins_orig(0), number_of_entities(0),
				   hash_space(0), bins(0), next_free(0)
#ifdef CONTACT_ANALYZE_HASH
      , hash_collisions(0)
#endif
{}

implementation::implementation(int n, ContactTopologyEntity<Real> **entities) :
  nbins(0), nbins_orig(0), number_of_entities(n),
  hash_space(0), bins(0), next_free(0)
#ifdef CONTACT_ANALYZE_HASH
      , hash_collisions(0)
#endif
{
  ComputeNbins(n);
  create_space();
  // Add each entity to the hash table
  for(int i=n-1 ;i>=0 ; --i) {
    insert( entities[i]->Global_ID(), entities[i] );
  }
}

implementation::implementation(std::vector<ContactTopologyEntity<Real>*> *v) :
  nbins(0), nbins_orig(0), number_of_entities(v->size()),
  hash_space(0), bins(0), next_free(0)
#ifdef CONTACT_ANALYZE_HASH
      , hash_collisions(0)
#endif
{
  // Create space for bins
  ComputeNbins(v->size());

  // Create hash space
  create_space();

  // Add each entity to the hash table
  for(int i = 0; i < v->size(); ++i) {
    ContactTopologyEntity<Real> *entity = (*v)[i];
    insert(entity->Global_ID(), entity);
  }
}

implementation::~implementation()
{
  ClearHash();
}

void implementation::ReHash(int number_entities, ContactTopologyEntity<Real> **entities)
{
  ClearHash();
  SetupHash(number_entities, entities);
}

void implementation::ClearHash()
{
  delete [] hash_space;
  delete [] bins;
  nbins              = 0;
  number_of_entities = 0;
  bins               = NULL;
  hash_space         = NULL;
  next_free          = NULL;
}

void implementation::insert( ContactHostGlobalID &Global_ID, ContactTopologyEntity<Real>* Entity) {
  hash* ptr;
  hash** previous;

  int ibin = hash_func(Global_ID.LoInt());

  previous = bins + ibin;
  ptr      = bins[ibin];
#ifdef CONTACT_ANALYZE_HASH  
  if (ptr!=NULL) ++hash_collisions;
#endif
  while( ptr != NULL && Global_ID < ptr->global_id ){
    previous = &(ptr->next);
    ptr = ptr->next;
  }
  if( ptr == NULL || Global_ID > ptr->global_id ){
    hash *newpt = next_free++;
    *previous = newpt;
    newpt->global_id = Global_ID;
    newpt->entity = Entity;
    newpt->next = ptr;
  }
}

ContactTopologyEntity<Real> *implementation::find(ContactHostGlobalID &Global_ID)
{
  hash* ptr;
  int ibin = hash_func(Global_ID.LoInt());
  ptr      = bins[ibin];

  while( ptr != NULL && Global_ID < ptr->global_id ){
    ptr = ptr->next;
  }
  if( ptr == NULL || Global_ID > ptr->global_id ){
    return(NULL);
  }
  return ptr->entity;
}

ContactTopologyEntity<Real> *implementation::find(ContactInteractionEntity<Real>::entity_data *data)
{
  ContactHostGlobalID id(data->host_gid[0], data->host_gid[1]);
  return find(id);
}

void implementation::SetupHash(int number_entities, ContactTopologyEntity<Real> **entities)
{
  number_of_entities = number_entities;

  // Create space for bins
  ComputeNbins(number_entities);

  // Create hash space
  create_space();

  // Add each entity to the hash table
  for( int i=number_entities-1 ; i>=0 ; --i )
    insert( entities[i]->Global_ID(), entities[i] );
}

void implementation::create_space()
{
  int n = number_of_entities;
  if( n<=0 ) n=1;
  hash_space          = new hash[n];
  next_free = hash_space;
}

void implementation::ComputeNbins(int n)
{
  nbins_orig = std::max(BIN_MINIMUM,(int)(BIN_FRACTION*n));
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
  if ( !(nbins % 6) ) nbins += 7;
#endif

  // Allocate the bins
  bins = new hash*[nbins];

  // Initialize the bins to empty
  for( int i=0 ; i<nbins ; ++i ) bins[i] = NULL;
}

int implementation::hash_func(int lo_int)
{
#if defined(JENKINS96_HASH_FUNC)
  unsigned int mask = nbins-1;
  unsigned int c    = (unsigned int)lo_int;
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
  unsigned int key  = (unsigned int)lo_int;
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
  unsigned int key  = (unsigned int)lo_int;
  key += ~(key << 15);
  key ^=  (key >> 10);
  key +=  (key <<  3);
  key ^=  (key >>  6);
  key += ~(key << 11);
  key ^=  (key >> 16);
  key  =  (key & mask);
  return (int)key;
#elif defined(SEDGEWICK_HASH_FUNC)
  return ((16161*lo_int)%nbins);
#elif defined(KNUTH_HASH_FUNC)
  return ((2654435761U*(unsigned)lo_int)%nbins);
#else
  return (lo_int % nbins);
#endif
}

std::ostream& implementation::stream_data(std::ostream& os) const
{
  Real avg=0.0;
  int  min=100000, max=0;
  os << "Number of Entries = " << number_of_entities << std::endl;
  for( int i=0 ; i<nbins ; ++i ){
    //os << "Outputting information for bin " << i << std::endl;
    struct hash* ptr = bins[i];
    int cnt=0;
    while( ptr != NULL ){
      //os << "    Global ID = " << ptr->entity->Global_ID() 
      //   << "  proc index = " << ptr->entity->ProcArrayIndex() 
      //   << "  block index = " << ptr->entity->BlockArrayIndex() << std::endl;
      ptr = ptr->next;
      ++cnt;
    }
    min  = std::min(min,cnt);
    max  = std::max(max,cnt);
    avg += cnt;
  }
  os << "Min = " << min << std::endl;
  os << "Max = " << max << std::endl;
  os << "Avg = " << avg/(Real)nbins << std::endl;
#ifdef CONTACT_ANALYZE_HASH
  os << "Num Collisions = " << hash_collisions << std::endl;
#endif
  return os;
}

ContactParOStream& implementation::stream_data(ContactParOStream& os) const
{
  int nused=0;
  int nempty=0;
  Real avg=0.0;
  int  min=100000, max=0;
  for( int i=0 ; i<nbins ; ++i ){
    struct hash* ptr = bins[i];
    if (ptr==NULL) {
      ++nempty;
    } else {
      ++nused;
      int cnt=0;
      while( ptr != NULL ){
        ptr = ptr->next;
        ++cnt;
      }
      min  = std::min(min,cnt);
      max  = std::max(max,cnt);
      avg += cnt;
    }
  }
#if defined(JENKINS96_HASH_FUNC)
  os << "  Jenkins 96 bit hash function\n";
#elif defined(JENKINS32_HASH_FUNC)
  os << "  Jenkins 32 bit hash function\n";
#elif defined(WANG32_HASH_FUNC)
  os << "  Wang 32 bit hash function\n";
#elif defined(SEDGEWICK_HASH_FUNC)
  os << "  Sedgewick modulus hash function\n";
#elif defined(KNUTH_HASH_FUNC)
  os << "  Knuth modulus hash function\n";
#else
  os << "  Plain modulus hash function\n";
#endif
#if defined(POWER2_HASH_BINS)
  os << "  Power of 2 number of hash bins\n";
#elif defined(PRIME_HASH_BINS)
  os << "  Prime number of hash bins\n";
#else
  os << "  Prime-like number of hash bins\n";
#endif
  os << "  Number of Entries     = " << number_of_entities << "\n";
  os << "  Number of Bins (orig) = " << nbins_orig << "\n";
  os << "  Number of Bins        = " << nbins << "\n";
  os << "  Number of Empty Bins  = " << nempty << "\n";
  os << "  Number of Used Bins   = " << nused << "\n";
  os << "    Min = " << min << "\n";
  os << "    Max = " << max << "\n";
  if (nused==0) {
    os << "      Avg = 0\n";
  } else {
    os << "    Avg = " << avg/(Real)nused << "\n";
  }
#ifdef CONTACT_ANALYZE_HASH
  os << "  Number of Collisions  = " << hash_collisions << "\n";
#endif
  return os;
}

} // namespace topology_hash_0

#ifdef CONTACT_HAVE_COMPILER_HASH
namespace topology_hash_1 {

  implementation::implementation(int n, ContactTopologyEntity<Real> **e)
  {
    SetupHash(n, e);
  }

  implementation::implementation(std::vector<ContactTopologyEntity<Real>*> *v)
  {
    SetupHash(v->size(), &(*v)[0]);
  }

  std::ostream& implementation::stream_data(std::ostream& os) const {
    return stream_data(os);
  }

  ContactParOStream& implementation::stream_data(ContactParOStream& os) const {
    return stream_data(os);
  }

  void implementation::ReHash(int n, ContactTopologyEntity<Real> **e) {
    ClearHash();
    SetupHash(n, e);
  }

  void implementation::AddEntities(ContactBlockEntityList *link_list) {
    link_list->IteratorStart();
    while(ContactTopologyEntity<Real> *entity = link_list->IteratorForward()) {
      hash_map<int, ContactTopologyEntity<Real> *>::value_type vt(entity->Global_ID().LoInt(), entity);
      container.insert(vt);
    }
  }

  void implementation::SetupHash(int n, ContactTopologyEntity<Real> **entities) {
    ContactTopologyEntity<Real> *e = 0;
    for (int i=0; i<n; ++i) {
      e = entities[i];
      hash_map<int, ContactTopologyEntity<Real> *>::value_type vt(e->Global_ID().LoInt(), e);
      container.insert(vt);
    }
  }

  void implementation::ClearHash() {
    // erase contents
    container.erase(container.begin(), container.end());
  }

  ContactTopologyEntity<Real> *implementation::find(ContactHostGlobalID &id) {
    // no check on whether container must contain id or not
    hash_map<int, ContactTopologyEntity<Real> *>::iterator it = container.find(id.LoInt());
    if (it != container.end()) {
      return (*it).second;
    }
    else {
      return 0;
    }
  }

  ContactTopologyEntity<Real> *implementation::find(ContactInteractionEntity<Real>::entity_data *data)
  {
    ContactHostGlobalID id(data->host_gid[0], data->host_gid[1]);
    return find(id);
  }

  void implementation::insert(ContactHostGlobalID &id, ContactTopologyEntity<Real> *e) {
    // no check on whether insert succeeded (id may already have been present,
    // but with a different entity
    hash_map<int, ContactTopologyEntity<Real> *>::value_type vt(id.LoInt(), e);
    container.insert(vt);
  }

} // namespace topology_hash_1
#endif
