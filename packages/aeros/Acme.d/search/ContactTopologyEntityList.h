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


#ifndef ContactTopologyEntityList_
#define ContactTopologyEntityList_

#include "ContactTopologyEntity.h"
#include "ContactInteractionEntity.h"
#include "ContactTopologyEntityHash.h"
#include "ContactHostGlobalID.h"
#include <iostream>

class ContactNodeBlock;
class ContactEdgeBlock;
class ContactFaceBlock;
class ContactElementBlock;

class ContactTopologyEntityList : public ContactTopologyEntityHash {
  
 public:
  ContactTopologyEntityList( );
  ~ContactTopologyEntityList( );

  void CleanUp();
  void Rehash();
  void BuildList( ContactNodeBlock** , int, int );
  void BuildList( ContactEdgeBlock** , int, int );
  void BuildList( ContactFaceBlock** , int, int );
  void BuildList( ContactElementBlock** , int, int );
  void BuildList( ContactNodeBlock** , ContactNodeBlock** , int, int );
  void BuildList( ContactEdgeBlock** , ContactEdgeBlock** , int, int );
  void BuildList( ContactFaceBlock** , ContactFaceBlock** , int, int );
  void BuildList( ContactElementBlock** , ContactElementBlock** , int, int );
  
  void SortByNodeGID();
  
  inline ContactTopologyEntity<Real>* Find( ContactTopologyEntity<Real>::connection_data* data )
    { ContactHostGlobalID GID( data->host_gid[0], data->host_gid[1] );
      return find(GID); };
  inline ContactTopologyEntity<Real>* Find( ContactInteractionEntity<Real>::entity_data* data )
    { ContactHostGlobalID GID( data->host_gid[0], data->host_gid[1] );
      return find(GID); };
  inline ContactTopologyEntity<Real>* Find( ContactHostGlobalID& GID ) 
    {return find(GID);};
  inline ContactTopologyEntity<Real>* Find( int index) {return entity_list[index];};
  inline void IteratorStart() {current=0;};
  inline void IteratorEnd  () {current=num_entities-1;};
  inline ContactTopologyEntity<Real>* IteratorForward() 
    { if (current==num_entities) return NULL; 
      else return entity_list[current++]; };
  inline ContactTopologyEntity<Real>* IteratorBackward()
    { if (current<0) return NULL; 
      else return entity_list[current--]; };
  inline int NumEntities() {return num_entities;};
  inline ContactTopologyEntity<Real>** EntityList() {return entity_list;};
  inline ContactTopologyEntity<Real>** BlockEntityList(int i) 
    {return num_entities>0?&entity_list[block_offset[i]]:NULL;};
  inline int Num_Blocks() {return nblocks;};
  inline int BlockNumEntities(int i) {return num_entities>0?block_cnt[i]:0;};
  
  
  void Display(ContactParOStream&);

 private:
  int  current;
  int  num_entities;
  int  nblocks;
  int* block_cnt;
  int* block_offset;
  ContactTopologyEntity<Real>** entity_list;
};

#endif // ContactTopologyEntityList_
