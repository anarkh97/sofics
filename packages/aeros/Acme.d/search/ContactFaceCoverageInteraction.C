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


#include "allocators.h"
#include "ContactFaceCoverageInteraction.h"
#include "CString.h"
#include "ContactTopologyEntityList.h"
#include "ContactFixedSizeAllocator.h"
#include <cstddef>
#include <cstring>
#include <new>

ContactFaceCoverageInteraction::ContactFaceCoverageInteraction( )
   : ContactInteractionEntity<Real>(DataArray, CT_FCI)
{
  num_vertices = 0;
  slave_face   = NULL;
  head         = NULL;
  tail         = NULL;
  slave_face_entity_data.type                      = -1;
  slave_face_entity_data.owner                     = -1;
  slave_face_entity_data.block_id                  = -1;
  slave_face_entity_data.index_in_host_array       = -1;
  slave_face_entity_data.index_in_owner_proc_array = -1;
  slave_face_entity_data.host_gid[0]               = -1;
  slave_face_entity_data.host_gid[1]               = -1;
}                                                                        

ContactFaceCoverageInteraction::ContactFaceCoverageInteraction( 
           ContactFace<Real>* face ) 
	: ContactInteractionEntity<Real>(DataArray, CT_FCI)
{
  PRECONDITION( face );
  slave_face = face;
  slave_face_global_id = face->Global_ID();
  num_vertices = 0;
  head         = NULL;
  tail         = NULL;
  slave_face_entity_data.type                      = -1;
  slave_face_entity_data.owner                     = -1;
  slave_face_entity_data.block_id                  = -1;
  slave_face_entity_data.index_in_host_array       = -1;
  slave_face_entity_data.index_in_owner_proc_array = -1;
  slave_face_entity_data.host_gid[0]               = -1;
  slave_face_entity_data.host_gid[1]               = -1;
}

ContactFaceCoverageInteraction::ContactFaceCoverageInteraction( 
                                           ContactFaceCoverageInteraction& fci )
: ContactInteractionEntity<Real>(DataArray, CT_FCI)
{
  num_vertices = fci.num_vertices;
  slave_face   = fci.slave_face;
  slave_face_global_id = fci.slave_face_global_id;
  if (num_vertices>0) {
    ContactFaceCoverageVertex* llnode;
    for( llnode=fci.Head(); llnode; llnode=llnode->next ){
      AddVertex(llnode);
    }
  }
}

ContactFaceCoverageInteraction* 
ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFace<Real>* face )
{
  return new (alloc.New_Frag())
    ContactFaceCoverageInteraction( face );
}

ContactFaceCoverageInteraction* 
ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
				     ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactFaceCoverageInteraction( );
}


ContactFaceCoverageInteraction* 
ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFaceCoverageInteraction& cfci )
{
  return new (alloc.New_Frag())
    ContactFaceCoverageInteraction( cfci );
}


void ContactFaceCoverageInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactFaceCoverageInteraction),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQFaceCoverageInteraction allocator" );
}


ContactFaceCoverageInteraction::~ContactFaceCoverageInteraction()
{
  ContactFaceCoverageVertex* llnode = head;
  while( llnode ){
    ContactFaceCoverageVertex* next = llnode->next;
    delete llnode;
    llnode = next;
  }
  head = NULL;
  tail = NULL;
  num_vertices = 0;
}

void ContactFaceCoverageInteraction::AddVertex(ContactFaceCoverageVertex* v) 
{ 
  if( head ){
    tail->next = v;
    tail = v;
  } else {
    head = v;
    tail = v;
  }
  num_vertices++;
}

void ContactFaceCoverageInteraction::AddVertex(Real x, Real y) 
{
  ContactFaceCoverageVertex* vv = new struct ContactFaceCoverageVertex;
  vv->slave_x = x;
  vv->slave_y = y;
  vv->next    = NULL;
  if( head ){
    tail->next = vv;
    tail = vv;
  } else {
    head = vv;
    tail = vv;
  }
  num_vertices++;
}

#ifndef CONTACT_NO_MPI
int ContactFaceCoverageInteraction::Size()
{
  return(ContactInteractionEntity<Real>::Size() + 
         sizeof(entity_data) +
         1*sizeof(int) + 
         DataArray_Length()*sizeof(Real) +
         (num_vertices)*sizeof(struct ContactFaceCoverageVertex));
}
#endif

#ifndef CONTACT_NO_MPI
void ContactFaceCoverageInteraction::Pack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity<Real>::Pack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity<Real>::Size());
  cnt += PackEntityData(&slave_face_entity_data, &i_buf[cnt]);
  i_buf[cnt++] = num_vertices;
  
  char* buf = buffer+ContactInteractionEntity<Real>::Size()+cnt*sizeof(int);
  std::memcpy( buf, DataArray, DataArray_Length()*sizeof(Real));
}
#endif

#ifndef CONTACT_NO_MPI
void ContactFaceCoverageInteraction::Unpack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity<Real>::Unpack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity<Real>::Size());
  cnt += UnPackEntityData(&slave_face_entity_data, &i_buf[cnt]);
  num_vertices = i_buf[cnt++];
  
  char* buf = buffer+ContactInteractionEntity<Real>::Size()+cnt*sizeof(int);
  std::memcpy( DataArray, buf, DataArray_Length()*sizeof(Real));
}
#endif

#ifndef CONTACT_NO_MPI
void ContactFaceCoverageInteraction::Copy( ContactFaceCoverageInteraction* src )
{
  ContactInteractionEntity<Real>::Copy( src );
  slave_face_entity_data  = src->slave_face_entity_data;
  num_vertices            = src->num_vertices;
  std::memcpy( DataArray,   src->DataArray, DataArray_Length()*sizeof(Real));
}
#endif

void ContactFaceCoverageInteraction::Connect_SlaveFace( ContactTopologyEntityList& hash_table )
{   
  slave_face = static_cast<ContactFace<Real> *>(hash_table.Find( slave_face_global_id ));
  POSTCONDITION( slave_face );
}

int ContactFaceCoverageInteraction::Data_Size()
{
  return 1+2*num_vertices;
}

int ContactFaceCoverageInteraction::Restart_Size()
{
  return sizeof(entity_data)/sizeof(int)+1+2*num_vertices;
}

void ContactFaceCoverageInteraction::Restart_Pack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += PackEntityData(&slave_face_entity_data, &buf_loc[cnt]);
  buf_loc[cnt++] = num_vertices;

  // Pack the data
  for( ContactFaceCoverageVertex* ll_node=head; 
       ll_node; ll_node=ll_node->next ){
    buf_loc[cnt++] = ll_node->slave_x;
    buf_loc[cnt++] = ll_node->slave_y;
  }
}

void ContactFaceCoverageInteraction::Restart_Unpack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += UnPackEntityData(&slave_face_entity_data, &buf_loc[cnt]);
  num_vertices = (int) buf_loc[cnt++];
  
  // Unpack the Data Array
  for (int i=0; i<num_vertices; ++i) {
    Real x = buf_loc[cnt++];
    Real y = buf_loc[cnt++];
    AddVertex(x, y);
  }
}



int ContactFaceCoverageInteraction::Set_SlaveFaceEntityData() 
{
  if (slave_face) {
    SetEntityData(&slave_face_entity_data, slave_face);
    return 1;
  }else{      
    return 0;
  }
}
