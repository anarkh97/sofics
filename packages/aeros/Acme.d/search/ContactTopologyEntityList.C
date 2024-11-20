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


#include "ContactTopologyEntityList.h"
#include "ContactNodeBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactElementBlock.h"
#include "ContactTopology.h"


ContactTopologyEntityList::ContactTopologyEntityList( )
  : ContactTopologyEntityHash( )
{ 
  current      = 0;
  num_entities = 0;
  nblocks      = 0;
  block_cnt    = NULL;
  block_offset = NULL;
  entity_list  = NULL;
}

ContactTopologyEntityList::~ContactTopologyEntityList()
{
  CleanUp();
}

void
ContactTopologyEntityList::CleanUp()
{
  ClearHash();
  if (block_cnt)    delete [] block_cnt;
  if (block_offset) delete [] block_offset;
  if (entity_list)  delete [] entity_list;
  block_cnt    = NULL;
  block_offset = NULL;
  entity_list  = NULL;
  num_entities = 0;
  nblocks      = 0;
}

void
ContactTopologyEntityList::BuildList( ContactNodeBlock** blocks, int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks]; 
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks[i]->NodeList()->NumEntities() 
                 == blocks[i]->Number_of_Nodes());
    block_cnt[i]    = blocks[i]->Number_of_Nodes();
    block_offset[i] = cnt;
    cnt += blocks[i]->Number_of_Nodes();
  }
  entity_list = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    cnt = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* node_list = blocks[i]->NodeList();
    node_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=node_list->IteratorForward()) {
      list[cnt] = entity;
      ++cnt;
    }
    if (sort) ContactTopology::SortEntityList(cnt, list);
    for (int j=0; j<cnt; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::BuildList( ContactEdgeBlock** blocks, int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks[i]->EdgeList()->NumEntities() 
                 == blocks[i]->Number_of_Edges());
    block_cnt[i]    = blocks[i]->Number_of_Edges();
    block_offset[i] = cnt;
    cnt += blocks[i]->Number_of_Edges();
  }
  entity_list  = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    cnt  = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* edge_list = blocks[i]->EdgeList();
    edge_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=edge_list->IteratorForward()) {
      list[cnt] = entity;
      ++cnt;
    }
    if (sort) ContactTopology::SortEntityList(cnt, list);
    for (int j=0; j<cnt; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::BuildList( ContactFaceBlock** blocks, int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks[i]->FaceList()->NumEntities() 
                 == blocks[i]->Number_of_Faces());
    block_cnt[i]    = blocks[i]->Number_of_Faces();
    block_offset[i] = cnt;
    cnt += blocks[i]->Number_of_Faces();
  }
  entity_list  = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    cnt  = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* face_list = blocks[i]->FaceList();
    face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=face_list->IteratorForward()) {
      list[cnt] = entity;
      ++cnt;
    }
    if (sort) ContactTopology::SortEntityList(cnt, list);
    for (int j=0; j<cnt; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::BuildList( ContactElementBlock** blocks, int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks[i]->ElemList()->NumEntities() 
                 == blocks[i]->Number_of_Elements());
    block_cnt[i]    = blocks[i]->Number_of_Elements();
    block_offset[i] = cnt;
    cnt += blocks[i]->Number_of_Elements();
  }
  entity_list  = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    cnt  = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* elem_list = blocks[i]->ElemList();
    elem_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=elem_list->IteratorForward()) {
      list[cnt] = entity;
      ++cnt;
    }
    if (sort) ContactTopology::SortEntityList(cnt, list);
    for (int j=0; j<cnt; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void
ContactTopologyEntityList::BuildList( ContactNodeBlock** blocks1, 
                                      ContactNodeBlock** blocks2, 
                                      int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks1[i]->NodeList()->NumEntities() 
                 == blocks1[i]->Number_of_Nodes());
    PRECONDITION(blocks2[i]->NodeList()->NumEntities() 
                 == blocks2[i]->Number_of_Nodes());
    block_offset[i] = cnt;
    block_cnt[i]    = blocks1[i]->Number_of_Nodes()+
                      blocks2[i]->Number_of_Nodes();
    cnt            += block_cnt[i];
  }
  entity_list = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    int k  = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* node_list = blocks1[i]->NodeList();
    node_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=node_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    node_list = blocks2[i]->NodeList();
    node_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=node_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    if (sort) ContactTopology::SortEntityList(k, list);
    for (int j=0; j<k; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::BuildList( ContactEdgeBlock** blocks1, 
                                      ContactEdgeBlock** blocks2,
                                      int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks1[i]->EdgeList()->NumEntities() 
                 == blocks1[i]->Number_of_Edges());
    PRECONDITION(blocks2[i]->EdgeList()->NumEntities() 
                 == blocks2[i]->Number_of_Edges());
    block_offset[i] = cnt;
    block_cnt[i]    = blocks1[i]->Number_of_Edges()+
                      blocks2[i]->Number_of_Edges();
    cnt            += block_cnt[i];
  }
  entity_list  = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    int k = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* edge_list = blocks1[i]->EdgeList();
    edge_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=edge_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    edge_list = blocks2[i]->EdgeList();
    edge_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=edge_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    if (sort) ContactTopology::SortEntityList(k, list);
    for (int j=0; j<k; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::BuildList( ContactFaceBlock** blocks1, 
                                      ContactFaceBlock** blocks2, 
                                      int num_blocks, int sort)
{
  CleanUp();
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks1[i]->FaceList()->NumEntities() 
                 == blocks1[i]->Number_of_Faces());
    PRECONDITION(blocks2[i]->FaceList()->NumEntities() 
                 == blocks2[i]->Number_of_Faces());
    block_offset[i] = cnt;
    block_cnt[i]    = blocks1[i]->Number_of_Faces()+
                      blocks2[i]->Number_of_Faces();
    cnt            += block_cnt[i];
  }
  entity_list  = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    int k = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* face_list = blocks1[i]->FaceList();
    face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=face_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    face_list = blocks2[i]->FaceList();
    face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=face_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    if (sort) ContactTopology::SortEntityList(k, list);
    for (int j=0; j<k; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::BuildList( ContactElementBlock** blocks1, 
                                      ContactElementBlock** blocks2, 
                                      int num_blocks, int sort)
{
  CleanUp();
  if (num_blocks==0) return;
  
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; ++i) {
    PRECONDITION(blocks1[i]->ElemList()->NumEntities() 
                 == blocks1[i]->Number_of_Elements());
    PRECONDITION(blocks2[i]->ElemList()->NumEntities() 
                 == blocks2[i]->Number_of_Elements());
    block_offset[i] = cnt;
    block_cnt[i]    = blocks1[i]->Number_of_Elements()+
                      blocks2[i]->Number_of_Elements();
    cnt            += block_cnt[i];
  }
  entity_list  = new ContactTopologyEntity<Real>*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; ++i) {
    int k = 0;
    ContactTopologyEntity<Real>** list = &entity_list[num_entities];
    ContactBlockEntityList* elem_list = blocks1[i]->ElemList();
    elem_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=elem_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    elem_list = blocks2[i]->ElemList();
    elem_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity=elem_list->IteratorForward()) {
      PRECONDITION(k<cnt);
      list[k] = entity;
      ++k;
    }
    if (sort) ContactTopology::SortEntityList(k, list);
    for (int j=0; j<k; ++j) {
      list[j]->ProcArrayIndex(num_entities++);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::Rehash( )
{
  ReHash(num_entities, entity_list);
}

void 
ContactTopologyEntityList::Display(ContactParOStream& postream)
{
  for (int i=0; i<num_entities; ++i) {
    entity_list[i]->Display(postream);
  }
}

void 
ContactTopologyEntityList::SortByNodeGID( )
{
  if (num_entities>1) {
    int n = 0;
    int cnt=1;
    int old_block = entity_list[0]->BlockID();
    ContactTopologyEntity<Real>** list = &entity_list[0];
    for (int i=1; i<num_entities; ++i) {
      int new_block = entity_list[i]->BlockID();
      if (new_block != old_block || i==num_entities-1) {
        ContactTopology::SortEntityList1(cnt, (ContactEdge<Real>**)list);
        for (int j=0; j<cnt; ++j) {
          list[j]->ProcArrayIndex(n++);
        }
        list = &entity_list[i];
	cnt  = 1;
      } else {
        ++cnt;
      }
    }
  }
}
