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

#include "ContactBlockEntityList.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSearch.h"
#include "ContactNode.h"
#include "ContactShellNode.h"
#include "ContactEdge.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactFace.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineFaceQ3.h"
#include "ContactElement.h"
#include "ContactHexElementL8.h"

using namespace std;

#define HashID_FindPtr(a,b)     while((a)!=NULL && b>(a)->id) (a)=(a)->next

#define HashID_InvalidPtr(a,b)     ((a)==NULL)||((b)<(a)->id)

#define Hash_AtTheEnd(a)           (a)==NULL
#define HashID_InTheMiddle(a,b)    (b)<(a)->id

ContactBlockEntityList::ContactBlockEntityList(ContactFixedSizeAllocator* alloc, int n)
  : ContactTopologyDLL(), do_hash(1),
    hash_size(0), hash_size_orig(0), 
    hash_nbins(0),hash_nbins_orig(0),
    hash_bins(NULL), 
#ifdef CONTACT_ANALYZE_HASH
    hash_collisions(0),
#endif
    allocators(alloc)
{ 
  SetupHash(n);
}

ContactBlockEntityList::ContactBlockEntityList(ContactFixedSizeAllocator* alloc)
  : ContactTopologyDLL(), do_hash(1),
    hash_size(0), hash_size_orig(0), 
    hash_nbins(0),hash_nbins_orig(0),
    hash_bins(NULL), 
#ifdef CONTACT_ANALYZE_HASH
    hash_collisions(0),
#endif
    allocators(alloc)
{ 
}

ContactBlockEntityList::~ContactBlockEntityList()
{
  CleanUp();
}

void
ContactBlockEntityList::CleanUp()
{
  int i;
  if (hash_nbins>0) {
    for (i=0; i<hash_nbins; ++i) {
      hash* hlink = hash_bins[i].head;
      while (hlink) {
        hash* ptr = hlink;
        hlink = hlink->next;
        delete ptr;
      }
    }
    delete [] hash_bins;
    hash_nbins = 0;
    hash_size  = 0;
    hash_size_orig  = 0;
  }
  Clear();
}

void ContactBlockEntityList::Insert( ContactTopologyEntity<Real>* entity )
{
  if (do_hash) {
  ContactHostGlobalID ID = entity->Global_ID();

  int       ibin = hash_func(ID.LoInt());
  hash_bin* bin  = &hash_bins[ibin];
  hash*     ptr  = bin->head;
  
  if (Hash_AtTheEnd(ptr)) {
    //=============================================
    // HASH BIN IS EMPTY SO WE KNOW THIS ENTITY 
    // IS NOT YET IN THE HASH TABLE SO JUST ADD IT
    //=============================================

    hash* hlink   = new hash;    
    hlink->list_index = Append(entity);
    hlink->id     = ID;
    hlink->entity = entity;
    PrependHash(bin, hlink);
  } else {
    //==================================================
    // HASH BIN IS NOT EMPTY SO SEE IF IT'S IN THIS BIN
    //==================================================
#ifdef CONTACT_ANALYZE_HASH
    ++hash_collisions;
#endif
    HashID_FindPtr(ptr,ID);
    if( Hash_AtTheEnd(ptr) ) {
      //===========================================
      // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
      //===========================================
      hash* hlink   = new hash;
      hlink->list_index = Append(entity);
      hlink->id     = ID;
      hlink->entity = entity;
      AppendHash(bin, hlink);
    } else if (HashID_InTheMiddle(ptr,ID)) {
      //===========================================
      // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
      //===========================================
      hash* hlink   = new hash;
      hlink->list_index = Append(entity);
      hlink->id     = ID;
      hlink->entity = entity;
      InsertBeforeHash(bin, hlink, ptr);
    } else {
      //=========================================
      // ENTITY IS ALREADY IN THE HASH TABLE BUT 
      // CHECK THE OWNERSHIP OF THE TWO ENTITIES
      //=========================================
      PRECONDITION (ID == ptr->id);
      ContactSearch::Delete_ContactTopologyEntity(entity, allocators);
    }
  }
  } else {
    Append(entity);
  }
}

#ifndef CONTACT_NO_MPI
void ContactBlockEntityList::Insert( char* buffer )
{
  int*  ibuffer = reinterpret_cast<int*>(buffer);
  
  if (do_hash) {
    ContactHostGlobalID ID(ibuffer[ContactTopologyEntity<Real>::GID_HI],
                           ibuffer[ContactTopologyEntity<Real>::GID_LO]);
                           
    int       ibin = hash_func(ID.LoInt());
    hash_bin* bin  = &hash_bins[ibin];
    hash*     ptr  = bin->head;
    
    if (Hash_AtTheEnd(ptr)) {
      //=============================================
      // HASH BIN IS EMPTY SO WE KNOW THIS ENTITY 
      // IS NOT YET IN THE HASH TABLE SO JUST ADD IT
      //=============================================
      ContactTopologyEntity<Real> *entity = CreateEntity(buffer);
      hash* hlink   = new hash;
      hlink->list_index = Append(entity);
      hlink->id     = ID;
      hlink->entity = entity;
      PrependHash(bin, hlink);
    } else {
      //==================================================
      // HASH BIN IS NOT EMPTY SO SEE IF IT'S IN THIS BIN
      //==================================================
#ifdef CONTACT_ANALYZE_HASH
      ++hash_collisions;
#endif
      HashID_FindPtr(ptr,ID);
      if (Hash_AtTheEnd(ptr)) {
        //===========================================
        // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
        //===========================================
        ContactTopologyEntity<Real> *entity = CreateEntity(buffer);
        hash* hlink   = new hash;
        hlink->list_index = Append(entity);
        hlink->id     = ID;
        hlink->entity = entity;
        AppendHash(bin, hlink);
      } else if (HashID_InTheMiddle(ptr,ID)) {
        //===========================================
        // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
        //===========================================
        ContactTopologyEntity<Real> *entity = CreateEntity(buffer);
        hash* hlink   = new hash;
        hlink->list_index = Append(entity);
        hlink->id     = ID;
        hlink->entity = entity;
        InsertBeforeHash(bin, hlink, ptr);
      } else {
        //=========================================
        // ENTITY IS ALREADY IN THE HASH TABLE
        //=========================================
        PRECONDITION (ID == ptr->id);
        if (ibuffer[ContactTopologyEntity<Real>::OWNERSHIP]==ContactTopologyEntity<Real>::OWNED &&
            ptr->entity->Ownership()!=ContactTopologyEntity<Real>::OWNED) {
          //===============================================
          // IF THE NEW ENTITY IS OWNED AND THE OLD ENTITY
          // IS NOT_OWNED THEN REPLACE IT WITH THE NEW ONE
          //===============================================
          // delete old entity and replace it with this one
          ContactTopologyEntity<Real>* old_entity = ptr->entity;
          ContactSearch::Delete_ContactTopologyEntity(old_entity, allocators);
          ContactTopologyEntity<Real>* entity = CreateEntity(buffer);
          SetEntity(ptr->list_index, entity);
        }
      }
    }
  } else {
    ContactTopologyEntity<Real> *entity = CreateEntity(buffer);
    Append(entity);
  }
}

void ContactBlockEntityList::Insert_ForSecondary( char* buffer )
{
  int*  ibuffer = reinterpret_cast<int*> (buffer);
  
  if (do_hash) {
    ContactHostGlobalID ID(ibuffer[ContactTopologyEntity<Real>::GID_HI],
                           ibuffer[ContactTopologyEntity<Real>::GID_LO]);
                           
    int       ibin = hash_func(ID.LoInt());
    hash_bin* bin  = &hash_bins[ibin];
    hash*     ptr  = bin->head;
    
    if (Hash_AtTheEnd(ptr)) {
      //=============================================
      // HASH BIN IS EMPTY SO WE KNOW THIS ENTITY 
      // IS NOT YET IN THE HASH TABLE SO JUST ADD IT
      //=============================================
      ContactTopologyEntity<Real> *entity = CreateEntity_ForSecondary(buffer);
      hash* hlink   = new hash;
      hlink->list_index = Append(entity);
      hlink->id     = ID;
      hlink->entity = entity;
      PrependHash(bin, hlink);
    } else {
      //==================================================
      // HASH BIN IS NOT EMPTY SO SEE IF IT'S IN THIS BIN
      //==================================================
#ifdef CONTACT_ANALYZE_HASH
      ++hash_collisions;
#endif
      HashID_FindPtr(ptr,ID);
      if (Hash_AtTheEnd(ptr)) {
        //===========================================
        // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
        //===========================================
        ContactTopologyEntity<Real> *entity = CreateEntity_ForSecondary(buffer);
        hash* hlink   = new hash;
        hlink->list_index = Append(entity);
        hlink->id     = ID;
        hlink->entity = entity;
        AppendHash(bin, hlink);
      } else if (HashID_InTheMiddle(ptr,ID)) {
        //===========================================
        // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
        //===========================================
        ContactTopologyEntity<Real> *entity = CreateEntity_ForSecondary(buffer);
        hash* hlink   = new hash;
        hlink->list_index = Append(entity);
        hlink->id     = ID;
        hlink->entity = entity;
        InsertBeforeHash(bin, hlink, ptr);
      } else {
        //=========================================
        // ENTITY IS ALREADY IN THE HASH TABLE
        //=========================================
        PRECONDITION (ID == ptr->id);
        if (ibuffer[ContactTopologyEntity<Real>::OWNERSHIP]==ContactTopologyEntity<Real>::OWNED &&
            ptr->entity->Ownership()!=ContactTopologyEntity<Real>::OWNED) {
          //===============================================
          // IF THE NEW ENTITY IS OWNED AND THE OLD ENTITY
          // IS NOT_OWNED THEN REPLACE IT WITH THE NEW ONE
          //===============================================
          // delete old entity and replace it with this one
          ContactTopologyEntity<Real>* old_entity = ptr->entity;
          ContactSearch::Delete_ContactTopologyEntity(old_entity, allocators);
          ContactTopologyEntity<Real>* entity = CreateEntity_ForSecondary(buffer);
          SetEntity(ptr->list_index, entity);
        }
      }
    }
  } else {
    ContactTopologyEntity<Real> *entity = CreateEntity_ForSecondary(buffer);
    Append(entity);
  }
}
#endif

void ContactBlockEntityList::Delete( ContactTopologyEntity<Real>* entity )
{
  if (!do_hash) {
    PRECONDITION(0);
    return;
  }
  ContactHostGlobalID ID = entity->Global_ID();

  int       ibin = hash_func(ID.LoInt());
  hash_bin* bin  = &hash_bins[ibin];
  hash*     ptr  = bin->head;
  
  HashID_FindPtr(ptr,ID);
  if (HashID_InvalidPtr(ptr,ID)) return;
  
  RemoveEntity(ptr->list_index);

  // delete the hash table link
  if (ptr == bin->head && ptr == bin->tail) {
    // only one entry in the hash bin
    bin->head = NULL;
    bin->tail = NULL;
  } else if (ptr == bin->head) {
    // entry is at the beginning of the hash bin
    hash* next = ptr->next;
    next->prev = NULL;
    bin->head  = next;
  } else if (ptr == bin->tail) {
    // entry is at the end of the hash bin
    hash* prev = ptr->prev;
    prev->next = NULL;
    bin->tail  = prev;
  } else {
    // entry is in the middle of the hash bin
    hash* prev = ptr->prev;
    hash* next = ptr->next;
    prev->next = next;
    next->prev = prev;
  }
  delete ptr;
  --hash_size;

  ContactSearch::Delete_ContactTopologyEntity(entity, allocators);
}

void
ContactBlockEntityList::SetupHash(int n)
{
  do_hash = 1;
  if (hash_nbins>0) {
    for (int i=0; i<hash_nbins; ++i) {
      hash* hlink = hash_bins[i].head;
      while (hlink) {
        hash* prev = hlink;
        hlink = hlink->next;
        delete prev;
      }
    }
    delete [] hash_bins;
  }
#ifdef CONTACT_ANALYZE_HASH
  hash_collisions = 0;
#endif
  
  // find the number of bins
  ComputeNbins(n);
  
  // Allocate the bins
  hash_bins = new hash_bin[hash_nbins];

  // Initialize the bins to empty
  for( int i=0 ; i<hash_nbins ; ++i ) {
    hash_bins[i].head = NULL;
    hash_bins[i].tail = NULL;
  }
}

void
ContactBlockEntityList::Rehash()
{
  if (!do_hash) {
    PRECONDITION(0);
    return;
  }
  SetupHash(NumEntities());
  
  int index = 0;
  std::vector< ContactTopologyEntity<Real>* > *entity_data = EntityList();
  std::vector< ContactTopologyEntity<Real>* >::iterator current_entity = entity_data->begin();
  while(current_entity < entity_data->end()) {
    ContactTopologyEntity<Real> *entity = *current_entity;
    if (entity!=NULL) {
      ContactHostGlobalID ID = entity->Global_ID();

      int       ibin = hash_func(ID.LoInt());
      hash_bin* bin  = &hash_bins[ibin];
      hash*     ptr  = bin->head;
      
      if (Hash_AtTheEnd(ptr)) {
        hash* hlink   = new hash;
        hlink->list_index = index;
        hlink->id     = ID;
        hlink->entity = entity;
        PrependHash(bin, hlink);
      } else {
        HashID_FindPtr(ptr,ID);
        if (Hash_AtTheEnd(ptr)) {
          // at the end of the bin, insert after ptr
          hash* hlink   = new hash;
          hlink->list_index = index;
          hlink->id     = ID;
          hlink->entity = entity;
          AppendHash(bin, hlink);
        } else if (HashID_InTheMiddle(ptr,ID)) {
          // middle of the bin, insert before ptr
          hash* hlink   = new hash;
          hlink->list_index = index;
          hlink->id     = ID;
          hlink->entity = entity;
          InsertBeforeHash(bin, hlink, ptr);
        }
      }
    }
    current_entity++;
    index++;
  }
}

void
ContactBlockEntityList::PrependHash( hash_bin* bin, hash* hlink )
{
  if (bin->head == NULL) {
    // HASH BIN IS EMPTY
    hlink->prev = NULL;
    hlink->next = NULL;
    bin->head   = hlink;
    bin->tail   = hlink;
  } else {
    hash* next  = bin->head;
    next->prev  = hlink;
    hlink->prev = NULL;
    hlink->next = next;
    bin->head   = hlink;
  }
  ++hash_size;
}

void
ContactBlockEntityList::AppendHash( hash_bin* bin, hash* hlink )
{
  if (bin->tail == NULL) {
    // HASH BIN IS EMPTY
    hlink->prev = NULL;
    hlink->next = NULL;
    bin->head   = hlink;
    bin->tail   = hlink;
  } else {
    hash* prev  = bin->tail;
    prev->next  = hlink;
    hlink->prev = prev;
    hlink->next = NULL;
    bin->tail   = hlink;
  }
  ++hash_size;
}

void
ContactBlockEntityList::InsertBeforeHash( hash_bin* bin, hash* hlink, hash* ptr )
{
  PRECONDITION (ptr != NULL);
  if (ptr == bin->head) {
    PrependHash(bin, hlink);
  } else {
    hash* prev  = ptr->prev;
    hash* next  = ptr;
    hlink->prev = prev;
    hlink->next = next;
    prev->next  = hlink;
    next->prev  = hlink;
    ++hash_size;
  }
}

ContactTopologyEntity<Real>* 
ContactBlockEntityList::Find( ContactHostGlobalID& ID )
{
  if (do_hash) {
    if (hash_nbins==0) return(NULL);
    int   ibin = hash_func(ID.LoInt());
    hash* ptr  = hash_bins[ibin].head;
    
    HashID_FindPtr(ptr,ID);
    if( HashID_InvalidPtr(ptr,ID) ) return(NULL);
    POSTCONDITION(ID == ptr->id);
    return ptr->entity;
  } else {
    POSTCONDITION(0);
    return NULL;
  }
}

ContactTopologyEntity<Real>* 
ContactBlockEntityList::Find( ContactTopologyEntity<Real>::connection_data* cdata )
{
  if (do_hash) {
    if (hash_nbins==0) return(NULL);
    ContactHostGlobalID GID( cdata->host_gid[0], cdata->host_gid[1] );
    return Find(GID);
  } else {
    POSTCONDITION(0);
    return NULL;
  }
}

ContactTopologyEntity<Real>* 
ContactBlockEntityList::Find( ContactInteractionEntity<Real>::entity_data* edata )
{
  if (do_hash) {
    if (hash_nbins==0) return(NULL);
    ContactHostGlobalID GID( edata->host_gid[0], edata->host_gid[1] );
    return Find(GID);
  } else {
    POSTCONDITION(0);
    return NULL;
  }
}

#ifndef CONTACT_NO_MPI

ContactTopologyEntity<Real>* 
ContactBlockEntityList::CreateEntity(char* buffer)
{
  ContactTopologyEntity<Real>* entity = NULL;
  int*  ibuffer     = reinterpret_cast<int*>(buffer);
  int   entity_type = ibuffer[ContactTopologyEntity<Real>::BASE_TYPE];
  switch (entity_type) {
  case CT_NODE:
    {
      ContactNode<Real>* node = ContactNode<Real>::new_ContactNode(allocators, ContactSearch::NODE );
      POSTCONDITION(node);
      node->Unpack(buffer);
      entity = node;
    }
    break;
  case CT_SHELL_NODE:
    {
      ContactShellNode* node = ContactShellNode::new_ContactShellNode( 
                allocators, ContactSearch::NODE );
      POSTCONDITION(node);
      node->Unpack(buffer); 
      entity = node;
    }
    break;
  case CT_EDGE:
    {
      ContactFixedSizeAllocator* alloc;
      ContactEdge<Real>* edge = NULL;
      int edge_type = ibuffer[ContactTopologyEntity<Real>::DERIVED_TYPE];
      switch( edge_type ){
      case ContactSearch::LINEEDGEL2 :
        alloc = &allocators[ContactSearch::ALLOC_ContactLineEdgeL2];
        edge  = ContactLineEdgeL2<Real>::new_ContactLineEdgeL2(*alloc);
        break;
      case ContactSearch::LINEEDGEQ3 :
        alloc = &allocators[ContactSearch::ALLOC_ContactLineEdgeQ3];
        edge  = ContactLineEdgeQ3::new_ContactLineEdgeQ3(*alloc);
        break;
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(edge);
      edge->Unpack(buffer);
      entity = edge;
    }
    break;
  case CT_FACE:
    {
      ContactSearch::ContactFace_Type face_type = (ContactSearch::ContactFace_Type)ibuffer[ContactTopologyEntity<Real>::DERIVED_TYPE];
      ContactFace<Real>* face = ContactSearch::New_ContactFace(face_type, allocators);
      POSTCONDITION(face);
      face->Unpack(buffer);
      entity = face;
    }
    break;
  case CT_ELEMENT:
    {
      ContactElement* element = NULL;
      int element_type = ibuffer[ContactTopologyEntity<Real>::DERIVED_TYPE];
      switch( element_type ){
      case ContactSearch::HEXELEMENTL8 :
        element = ContactHexElementL8::new_ContactHexElementL8(allocators);
        break;
      case ContactSearch::CARTESIANHEXELEMENTL8 :
        element = ContactCartesianHexElementL8::new_ContactCartesianHexElementL8(allocators);
        break;
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(element);
      element->Unpack(buffer);
      entity = element;
    }
    break;
  default:
    POSTCONDITION(false);
    break;
  }
  POSTCONDITION(entity);
  return entity;
}

ContactTopologyEntity<Real>* 
ContactBlockEntityList::CreateEntity_ForSecondary(char* buffer)
{
  ContactTopologyEntity<Real>* entity = NULL;
  int*  ibuffer     = reinterpret_cast<int*> (buffer);
  int   entity_type = ibuffer[ContactTopologyEntity<Real>::BASE_TYPE];
  switch (entity_type) {
  case CT_NODE:
    {
      ContactNode<Real>* node = ContactNode<Real>::new_ContactNode(allocators, ContactSearch::NODE );
      POSTCONDITION(node);
      node->Unpack_ForSecondary(buffer);
      entity = node;
    }
    break;
  case CT_SHELL_NODE:
    {
      ContactShellNode* node = ContactShellNode::new_ContactShellNode( 
                allocators, ContactSearch::NODE );
      POSTCONDITION(node);
      node->Unpack_ForSecondary(buffer); 
      entity = node;
    }
    break;
  case CT_EDGE:
    POSTCONDITION(false);
    break;
  case CT_FACE:
    {
      ContactSearch::ContactFace_Type face_type = (ContactSearch::ContactFace_Type)ibuffer[ContactTopologyEntity<Real>::DERIVED_TYPE];
      ContactFace<Real>* face = ContactSearch::New_ContactFace(face_type, allocators);
      POSTCONDITION(face);
      face->Unpack_ForSecondary(buffer);
      entity = face;
    }
    break;
  case CT_ELEMENT:
    {
      ContactElement* element = NULL;
      int element_type = ibuffer[ContactTopologyEntity<Real>::DERIVED_TYPE];
      switch( element_type ){
      case ContactSearch::HEXELEMENTL8 :
        element = ContactHexElementL8::new_ContactHexElementL8(allocators);
        break;
      case ContactSearch::CARTESIANHEXELEMENTL8 :
        element = ContactCartesianHexElementL8::new_ContactCartesianHexElementL8(allocators);
        break;
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(element);
      element->Unpack_ForSecondary(buffer);
      entity = element;
    }
    break;
  default:
    POSTCONDITION(false);
    break;
  }
  POSTCONDITION(entity);
  return entity;
}
#endif

void 
ContactBlockEntityList::ComputeNbins(int n)
{
  hash_nbins_orig = std::max(BIN_MINIMUM,(int)(BIN_FRACTION*n));
#if defined(PRIME_HASH_BINS)
  // Round up nbins to make it prime
  int s[] = {1,5,3,3,9,3,1,3,19,15,1,5,1,3,9,3,15,3,39,5,39,57,3,35,1};
  int cnt = 0;
  for (hash_nbins=BIN_MINIMUM; hash_nbins<hash_nbins_orig; hash_nbins*=2, ++cnt);
  hash_nbins -= s[cnt];
#elif defined(POWER2_HASH_BINS)
  // Round up nbins to make it a power of 2
  for (hash_nbins=2; hash_nbins<hash_nbins_orig; hash_nbins *= 2);
#else
  // Round up nbins to make it more prime-like
  hash_nbins = hash_nbins_orig;
  if ( !(hash_nbins % 2) ) hash_nbins += 1;
  if ( !(hash_nbins % 3) ) hash_nbins += 2;
  if ( !(hash_nbins % 6) ) hash_nbins += 6;
#endif
}

std::ostream& operator<<( std::ostream& os, const ContactBlockEntityList& LIST )
{
  //os << "Linked List -----------------------------------" << std::endl << std::endl;
  //os << "Number of Entries = " << LIST.NumEntities() << std::endl << std::endl;
  //LIST.DisplayList();
  
  
  os << "Hash Table ------------------------------------" << std::endl << std::endl;
  os << "Number of Bins    = " << LIST.hash_nbins << std::endl;
  os << "Number of Entries = " << LIST.hash_size << std::endl;
  os << std::endl;
  Real avg=0.0;
  int  min=100000, max=0;
  for( int i=0 ; i<LIST.hash_nbins ; ++i ){
    struct ContactBlockEntityList::hash* ptr = LIST.hash_bins[i].head;
    if (ptr==NULL) continue;
    //os << "  Bin " << i 
    //   << "  head = "<<LIST.hash_bins[i].head
    //   << "  tail = "<<LIST.hash_bins[i].tail
    //   << std::endl << std::endl;
    int cnt=0;
    while( ptr != NULL ){
    //  os << "    Address   = " << ptr <<std::endl;
    //  os << "    Global ID = " << ptr->id  << std::endl;
    //  os << "    prev      = " << ptr->prev << std::endl;
    //  os << "    next      = " << ptr->next << std::endl<<std::endl;
      ptr = ptr->next;
      ++cnt;
    }
    min  = std::min(min,cnt);
    max  = std::max(max,cnt);
    avg += cnt;
  }
  os << "Min = " << min << std::endl;
  os << "Max = " << max << std::endl;
  os << "Avg = " << avg/(Real)LIST.hash_nbins << std::endl;
#ifdef CONTACT_ANALYZE_HASH
  os << "Number of Collisions = " << LIST.hash_collisions << std::endl;
#endif
  return os;
}

ContactParOStream& operator<<( ContactParOStream& os, 
			       const ContactBlockEntityList& LIST )
{
  if (LIST.do_hash) {
  int nempty=0;
  int nused=0;
  Real avg=0.0;
  int  min=100000, max=0;
  for( int i=0 ; i<LIST.hash_nbins ; ++i ){
    struct ContactBlockEntityList::hash* ptr = LIST.hash_bins[i].head;
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
  os << "    Jenkins 96 bit hash function\n";
#elif defined(JENKINS32_HASH_FUNC)
  os << "    Jenkins 32 bit hash function\n";
#elif defined(WANG32_HASH_FUNC)
  os << "    Wang 32 bit hash function\n";
#elif defined(SEDGEWICK_HASH_FUNC)
  os << "    Sedgewick modulus hash function\n";
#elif defined(KNUTH_HASH_FUNC)
  os << "    Knuth modulus hash function\n";
#else
  os << "    Plain modulus hash function\n";
#endif
#if defined(POWER2_HASH_BINS)
  os << "    Power of 2 number of hash bins\n";
#elif defined(PRIME_HASH_BINS)
  os << "    Prime number of hash bins\n";
#else
  os << "    Prime-like number of hash bins\n";
#endif
  os << "    Initial Guess for N   = " << LIST.hash_size_orig << "\n";
  os << "    Number of Entries     = " << LIST.hash_size << "\n";
  os << "    Number of Bins (orig) = " << LIST.hash_nbins_orig << "\n";
  os << "    Number of Bins        = " << LIST.hash_nbins << "\n";
  os << "    Number of Empty Bins  = " << nempty << "\n";
  os << "    Number of Used Bins   = " << nused << "\n";
  os << "      Min = " << min << "\n";
  os << "      Max = " << max << "\n";
  if (nused==0) {
    os << "      Avg = 0\n";
  } else {
    os << "      Avg = " << avg/(Real)nused << "\n";
  }
  if (LIST.hash_size==0) {
    os << "    Actual Scale Factor   = n/a\n";
  } else {
    os << "    Actual Scale Factor   = " << (Real)LIST.hash_nbins/(Real)LIST.hash_size << "\n";
  }
#ifdef CONTACT_ANALYZE_HASH
  os << "    Number of Collisions  = " << LIST.hash_collisions << "\n";
#endif
  } else {
  os << "    No hash table present\n";
  }
  return os;
}
