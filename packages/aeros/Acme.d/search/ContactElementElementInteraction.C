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
#include "ContactElementElementInteraction.h"
#include "CString.h"
#include "ContactTopologyEntityHash.h"
#include "ContactTopologyEntityList.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include "ContactElementBlock.h"
#include <cstddef>
#include <cstring>
#include <new>

ContactElementElementInteraction::ContactElementElementInteraction( )
   : ContactInteractionEntity<Real>(DataArray, CT_EEI)
{
  slave_element  = NULL;
  master_element = NULL;
  slave_element_entity_data.type                       = -1;
  slave_element_entity_data.owner                      = -1;
  slave_element_entity_data.block_id                   = -1;
  slave_element_entity_data.index_in_host_array        = -1;
  slave_element_entity_data.index_in_owner_proc_array  = -1;
  slave_element_entity_data.host_gid[0]                = -1;
  slave_element_entity_data.host_gid[1]                = -1;
  master_element_entity_data.type                      = -1;
  master_element_entity_data.owner                     = -1;
  master_element_entity_data.block_id                  = -1;
  master_element_entity_data.index_in_host_array       = -1;
  master_element_entity_data.index_in_owner_proc_array = -1;
  master_element_entity_data.host_gid[0]               = -1;
  master_element_entity_data.host_gid[1]               = -1;
}                                                                               

ContactElementElementInteraction::ContactElementElementInteraction( ContactElement* Selement,
							ContactElement* Melement,
                                                        Real Volume )
: ContactInteractionEntity<Real>(DataArray, CT_EEI)
{
  PRECONDITION( Selement && Melement );
  slave_element  = Selement;
  master_element = Melement;
  Set_SlaveElementEntityData();
  Set_MasterElementEntityData();
  Scalar_Var(VOLUME) = Volume;
}

ContactElementElementInteraction::ContactElementElementInteraction( 
                                       ContactElementElementInteraction& eei )
: ContactInteractionEntity<Real>(DataArray, CT_EEI)
{
  slave_element              = eei.slave_element;
  master_element             = eei.master_element;
  slave_element_entity_data  = eei.slave_element_entity_data;
  master_element_entity_data = eei.master_element_entity_data;
  std::memcpy( &(DataArray[0]), &(eei.DataArray[0]), 
	  (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real) );
}

ContactElementElementInteraction* 
ContactElementElementInteraction::new_ContactElementElementInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactElement* Selement,
				     ContactElement* Melement,
                                     Real Volume )
{
  return new (alloc.New_Frag())
    ContactElementElementInteraction( Selement, Melement, Volume );
}

ContactElementElementInteraction* 
ContactElementElementInteraction::new_ContactElementElementInteraction(
				     ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactElementElementInteraction( );
}


ContactElementElementInteraction* 
ContactElementElementInteraction::new_ContactElementElementInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactElementElementInteraction& ceei )
{
  return new (alloc.New_Frag())
    ContactElementElementInteraction( ceei );
}


void ContactElementElementInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactElementElementInteraction),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQElementElementInteraction allocator" );
}


ContactElementElementInteraction::~ContactElementElementInteraction()
{
}

int ContactElementElementInteraction::Size()
{
  return(ContactInteractionEntity<Real>::Size() +
         2*sizeof(entity_data) +
         sizeof(Real) +
         DataArray_Length()*sizeof(Real));
}

void ContactElementElementInteraction::Pack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity<Real>::Pack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity<Real>::Size());
  cnt += PackEntityData(&slave_element_entity_data, &i_buf[cnt]);
  cnt += PackEntityData(&master_element_entity_data, &i_buf[cnt]);
  
  char* buf = buffer+ContactInteractionEntity<Real>::Size()+cnt*sizeof(int);
  
  // Pack the data
  std::memcpy( buf, DataArray, DataArray_Length()*sizeof(Real) );
}

void ContactElementElementInteraction::Unpack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity<Real>::Unpack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity<Real>::Size());
  cnt += UnPackEntityData(&slave_element_entity_data, &i_buf[cnt]);
  cnt += UnPackEntityData(&master_element_entity_data, &i_buf[cnt]);
  
  char* buf = buffer+ContactInteractionEntity<Real>::Size()+cnt*sizeof(int);
  std::memcpy( DataArray, buf, DataArray_Length()*sizeof(Real) );
}

void ContactElementElementInteraction::Copy( ContactElementElementInteraction* src )
{
  ContactInteractionEntity<Real>::Copy( src );
  slave_element_entity_data  = src->slave_element_entity_data;
  master_element_entity_data = src->master_element_entity_data;
  std::memcpy( DataArray,      src->DataArray_Buffer(), DataArray_Length()*sizeof(Real) );
}

void ContactElementElementInteraction::Connect_SlaveElement( ContactTopologyEntityList& hash_table )
{
  slave_element = static_cast<ContactElement *>(hash_table.Find( &slave_element_entity_data ));
  POSTCONDITION( slave_element );
}

void ContactElementElementInteraction::Connect_MasterElement( ContactTopologyEntityList& hash_table )
{
  master_element = static_cast<ContactElement *>(hash_table.Find( &master_element_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_element );
}

void ContactElementElementInteraction::Connect_SlaveElement( ContactTopologyEntityHash& hash_table )
{
  ContactHostGlobalID gid( slave_element_entity_data.host_gid[0], 
                           slave_element_entity_data.host_gid[1] );
  slave_element = static_cast<ContactElement *>(hash_table.find( gid ));
  POSTCONDITION( slave_element );
}

void ContactElementElementInteraction::Connect_MasterElement( ContactTopologyEntityHash& hash_table )
{
  ContactHostGlobalID gid( master_element_entity_data.host_gid[0], 
                           master_element_entity_data.host_gid[1] );
  master_element = static_cast<ContactElement *>(hash_table.find( gid ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_element );
}

void ContactElementElementInteraction::Connect_SlaveElement( ContactTopology* topology )
{   
  int block = slave_element_entity_data.block_id;
  slave_element = static_cast<ContactElement *>
    (topology->Element_Block(block)->ElemList()->Find( &slave_element_entity_data ));
  POSTCONDITION( slave_element );
}
void ContactElementElementInteraction::Connect_MasterElement( ContactTopology* topology )
{
  int block = master_element_entity_data.block_id;
  master_element = static_cast<ContactElement *>
    (topology->Element_Block(block)->ElemList()->Find( &master_element_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_element );
}

void ContactElementElementInteraction::Connect_SlaveElement( ContactElement* Element )
{   
  slave_element = Element;
  POSTCONDITION( slave_element );
}
void ContactElementElementInteraction::Connect_MasterElement( ContactElement* Element )
{
  master_element = Element;
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_element );
}

int ContactElementElementInteraction::Data_Size()
{
  return 3;
}

int ContactElementElementInteraction::Restart_Size()
{
  return 2*sizeof(entity_data)/sizeof(int)+
         NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;
}

void ContactElementElementInteraction::Restart_Pack( Real* buffer )
{
  int cnt=0;
  Real* buf_loc = buffer;

  cnt += PackEntityData(&slave_element_entity_data, &buf_loc[cnt]);
  cnt += PackEntityData(&master_element_entity_data, &buf_loc[cnt]);
  
  // Pack the data
  std::memcpy( &buf_loc[cnt], &DataArray, 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );
}

void ContactElementElementInteraction::Restart_Unpack( Real* buffer )
{
  int cnt=0;
  Real* buf_loc = buffer;

  cnt += UnPackEntityData(&slave_element_entity_data, &buf_loc[cnt]);
  cnt += UnPackEntityData(&master_element_entity_data, &buf_loc[cnt]);
  
  // Unpack the Data Array
  std::memcpy( &DataArray, &buf_loc[cnt], 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );
}

int ContactElementElementInteraction::Set_SlaveElementEntityData() 
{
  if (slave_element) {
    SetEntityData(&slave_element_entity_data, slave_element);
    return 1;
  }else{      
    return 0;
  }
}

int ContactElementElementInteraction::Set_MasterElementEntityData() 
{
  if (master_element) {
    SetEntityData(&master_element_entity_data, master_element);
    return 1;
  }else{      
    return 0;
  }
}
