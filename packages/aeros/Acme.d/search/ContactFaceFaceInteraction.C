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


#ifndef ContactFaceFaceInteraction_C_
#define ContactFaceFaceInteraction_C_

#include "allocators.h"
#include "ContactFaceFaceInteraction.h"
#include "CString.h"
#include "ContactTopologyEntityList.h"
#include "ContactTopologyEntityHash.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include "ContactFaceBlock.h"
#include <cstddef>
#include <cstring>
#include <new>

template <typename DataType>
ContactFaceFaceInteraction<DataType>::ContactFaceFaceInteraction( )
   : ContactInteractionEntity<DataType>(DataArray, CT_FFI)
{
  slave_face  = NULL;
  master_face = NULL;
  num_edges   = 0;
  num_derivatives = 0;
  num_second_derivatives = 0;
  vertices    = NULL;
}

template <typename DataType>
ContactFaceFaceInteraction<DataType>::ContactFaceFaceInteraction( ContactFace<DataType>* Sface,
							ContactFace<DataType>* Mface,
							int Nedges, 
                                                        int* FaceEdge,
                                                        int* EdgeMaster,
                                                        DataType* Sarea, 
                                                        DataType* Marea,
                                                        DataType (*Sarea_derivatives)[42],
                                                        DataType (*Marea_derivatives)[42],
                                                        DataType (*Sarea_second_derivatives)[42],
                                                        DataType (*Marea_second_derivatives)[42] )
: ContactInteractionEntity<DataType>(DataArray, CT_FFI)
{
  PRECONDITION( Sface && Mface );
  num_edges   = Nedges;
  slave_face  = Sface;
  master_face = Mface;
  Set_SlaveFaceEntityData();
  Set_MasterFaceEntityData();
  num_derivatives = (Sarea_derivatives && Marea_derivatives) ?
                    3*(slave_face->Nodes_Per_Face() + master_face->Nodes_Per_Face()) : 0;
  num_second_derivatives = (Sarea_second_derivatives && Marea_second_derivatives) ?
                           num_derivatives*num_derivatives : 0;
  if (num_edges>0) {
    int i, j;
    vertices = new ContactFaceFaceVertex<DataType>[num_edges+1];
    for (i=0; i<num_edges; ++i) {
      vertices[i].slave_x          = Sarea[2*i];
      vertices[i].slave_y          = Sarea[2*i+1];
      vertices[i].master_x         = Marea[2*i];
      vertices[i].master_y         = Marea[2*i+1];
      vertices[i].slave_edge_id    = FaceEdge[i];
      vertices[i].master_edge_flag = EdgeMaster[i];
      if(num_derivatives > 0) {
        vertices[i].slave_x_derivatives = new DataType[num_derivatives];
        vertices[i].slave_y_derivatives = new DataType[num_derivatives];
        vertices[i].master_x_derivatives = new DataType[num_derivatives];
        vertices[i].master_y_derivatives = new DataType[num_derivatives];
        for(j = 0; j<num_derivatives; ++j) {
          vertices[i].slave_x_derivatives[j] = Sarea_derivatives[j][2*i];
          vertices[i].slave_y_derivatives[j] = Sarea_derivatives[j][2*i+1];
          vertices[i].master_x_derivatives[j] = Marea_derivatives[j][2*i];
          vertices[i].master_y_derivatives[j] = Marea_derivatives[j][2*i+1];
        }
      }
      if(num_second_derivatives > 0) {
        vertices[i].slave_x_second_derivatives = new DataType[num_second_derivatives];
        vertices[i].slave_y_second_derivatives = new DataType[num_second_derivatives];
        vertices[i].master_x_second_derivatives = new DataType[num_second_derivatives];
        vertices[i].master_y_second_derivatives = new DataType[num_second_derivatives];
        for(j = 0; j<num_second_derivatives; ++j) {
          vertices[i].slave_x_second_derivatives[j] = Sarea_second_derivatives[j][2*i];
          vertices[i].slave_y_second_derivatives[j] = Sarea_second_derivatives[j][2*i+1];
          vertices[i].master_x_second_derivatives[j] = Marea_second_derivatives[j][2*i];
          vertices[i].master_y_second_derivatives[j] = Marea_second_derivatives[j][2*i+1];
        }
      }
    }
    i = num_edges;
    vertices[i].slave_x          = Sarea[2*i];
    vertices[i].slave_y          = Sarea[2*i+1];
    vertices[i].master_x         = Marea[2*i];
    vertices[i].master_y         = Marea[2*i+1];
    vertices[i].slave_edge_id    = 0;
    vertices[i].master_edge_flag = 0;
  }
}

template <typename DataType>
ContactFaceFaceInteraction<DataType>::ContactFaceFaceInteraction( 
                                       ContactFaceFaceInteraction& ffi )
: ContactInteractionEntity<DataType>(DataArray, CT_FFI)
{
  slave_face              = ffi.slave_face;
  master_face             = ffi.master_face;
  slave_face_entity_data  = ffi.slave_face_entity_data;
  master_face_entity_data = ffi.master_face_entity_data;
  num_edges               = ffi.num_edges;
  if (num_edges>0) {
    vertices = new ContactFaceFaceVertex<DataType>[num_edges+1];
    for (int i=0; i<num_edges+1; ++i) {
      vertices[i] = ffi.vertices[i];
    }
  }
}

template <typename DataType>
ContactFaceFaceInteraction<DataType>* 
ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFace<DataType>* Sface,
				     ContactFace<DataType>* Mface,
				     int Nedges, 
                                     int* FaceEdge,
                                     int* EdgeMaster,
                                     DataType* Sarea, DataType* Marea,
                                     DataType (*Sarea_derivatives)[42], DataType (*Marea_derivatives)[42],
                                     DataType (*Sarea_second_derivatives)[42], DataType (*Marea_second_derivatives)[42] )
{
  return new (alloc.New_Frag())
    ContactFaceFaceInteraction( Sface, Mface, Nedges, 
                                FaceEdge, EdgeMaster, Sarea, Marea, Sarea_derivatives, Marea_derivatives,
                                Sarea_second_derivatives, Marea_second_derivatives );
}

template <typename DataType>
ContactFaceFaceInteraction<DataType>* 
ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactFaceFaceInteraction( );
}

template <typename DataType>
ContactFaceFaceInteraction<DataType>* 
ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFaceFaceInteraction& cffi )
{
  return new (alloc.New_Frag())
    ContactFaceFaceInteraction( cffi );
}

template <typename DataType>
void ContactFaceFaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactFaceFaceInteraction<DataType>),
                100,  // block size
                0,    // initial block size
                sizeof(DataType)); 
  alloc.Set_Name( "ContactQFaceFaceInteraction allocator" );
}

template <typename DataType>
ContactFaceFaceInteraction<DataType>::~ContactFaceFaceInteraction()
{
  if (num_derivatives > 0) {
    for (int i=0; i<num_edges; ++i) {
      delete [] vertices[i].slave_x_derivatives;
      delete [] vertices[i].slave_y_derivatives;
      delete [] vertices[i].master_x_derivatives;
      delete [] vertices[i].master_y_derivatives;
    }
  }
  if (num_second_derivatives > 0) {
    for (int i=0; i<num_edges; ++i) {
      delete [] vertices[i].slave_x_second_derivatives;
      delete [] vertices[i].slave_y_second_derivatives;
      delete [] vertices[i].master_x_second_derivatives;
      delete [] vertices[i].master_y_second_derivatives;
    }
  }
  if (vertices) delete [] vertices;
}

template <typename DataType>
int ContactFaceFaceInteraction<DataType>::Size()
{
  return(ContactInteractionEntity<DataType>::Size() + 
         2*sizeof(typename ContactInteractionEntity<DataType>::entity_data)+
         1*sizeof(int) + 
         DataArray_Length()*sizeof(DataType) +
         (num_edges+1)*sizeof(ContactFaceFaceVertex<DataType>));
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Pack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity<DataType>::Pack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity<DataType>::Size());
  cnt += this->PackEntityData(&slave_face_entity_data, &i_buf[cnt]);
  cnt += this->PackEntityData(&master_face_entity_data, &i_buf[cnt]);
  i_buf[cnt++] = num_edges;
  
  char* buf = buffer+ContactInteractionEntity<DataType>::Size()+cnt*sizeof(int);
  std::memcpy( buf, DataArray, DataArray_Length()*sizeof(DataType));

  buf += DataArray_Length()*sizeof(DataType);
  std::memcpy( buf, vertices,(num_edges+1)*sizeof(ContactFaceFaceVertex<DataType>));
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Unpack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity<DataType>::Unpack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity<DataType>::Size());
  cnt += this->UnPackEntityData(&slave_face_entity_data, &i_buf[cnt]);
  cnt += this->UnPackEntityData(&master_face_entity_data, &i_buf[cnt]);
  num_edges = i_buf[cnt++];
  
  char* buf = buffer+ContactInteractionEntity<DataType>::Size()+cnt*sizeof(int);
  std::memcpy( DataArray, buf, DataArray_Length()*sizeof(DataType));

  buf += DataArray_Length()*sizeof(DataType);
  vertices  = new ContactFaceFaceVertex<DataType>[num_edges+1];
  std::memcpy( vertices,buf,(num_edges+1)*sizeof(ContactFaceFaceVertex<DataType>));
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Copy( ContactFaceFaceInteraction* src )
{
  ContactInteractionEntity<DataType>::Copy( src );
  slave_face_entity_data  = src->slave_face_entity_data;
  master_face_entity_data = src->master_face_entity_data;
  num_edges               = src->num_edges;
  std::memcpy( DataArray,   src->DataArray, DataArray_Length()*sizeof(DataType));
  vertices                = new ContactFaceFaceVertex<DataType>[num_edges+1];
  std::memcpy( vertices,    src->vertices,(num_edges+1)*sizeof(ContactFaceFaceVertex<DataType>));
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_SlaveFace( ContactTopologyEntityList& hash_table )
{   
  slave_face = static_cast<ContactFace<DataType> *>(hash_table.Find( &slave_face_entity_data ));
  POSTCONDITION( slave_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_MasterFace( ContactTopologyEntityList& hash_table )
{
  master_face = static_cast<ContactFace<DataType> *>(hash_table.Find( &master_face_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_SlaveFace( ContactTopologyEntityHash& hash_table )
{   
  slave_face = static_cast<ContactFace<DataType> *>(hash_table.find( &slave_face_entity_data ));
  POSTCONDITION( slave_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_MasterFace( ContactTopologyEntityHash& hash_table )
{
  master_face = static_cast<ContactFace<DataType> *>(hash_table.find( &master_face_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_SlaveFace( ContactTopology* topology )
{   
  int block = slave_face_entity_data.block_id;
  slave_face = static_cast<ContactFace<DataType> *>
    (topology->Face_Block(block)->FaceList()->Find( &slave_face_entity_data ));
  POSTCONDITION( slave_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_MasterFace( ContactTopology* topology )
{
  int block = master_face_entity_data.block_id;
  master_face = static_cast<ContactFace<DataType> *>
    (topology->Face_Block(block)->FaceList()->Find( &master_face_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( slave_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_SlaveFace( ContactFace<DataType>* Face )
{   
  slave_face = Face;
  POSTCONDITION( slave_face );
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Connect_MasterFace( ContactFace<DataType>* Face )
{
  master_face = Face;
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_face );
}

template <typename DataType>
int ContactFaceFaceInteraction<DataType>::Data_Size()
{
  return 1+num_edges+num_edges+4*num_edges*(1+num_derivatives+num_second_derivatives);
}

template <typename DataType>
int ContactFaceFaceInteraction<DataType>::Restart_Size()
{
  return 2*sizeof(typename ContactInteractionEntity<DataType>::entity_data)/sizeof(int)+1+6*(num_edges+1);
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Restart_Pack( DataType* buffer )
{
  int i;
  int cnt=0;
  DataType* buf_loc = buffer;

  cnt += this->PackEntityData(&slave_face_entity_data, &buf_loc[cnt]);
  cnt += this->PackEntityData(&master_face_entity_data, &buf_loc[cnt]);

  // Pack the data
  buf_loc[cnt++] = num_edges;
  for (i=0; i<num_edges+1; ++i) {
    buf_loc[cnt++] = vertices[i].slave_x;
    buf_loc[cnt++] = vertices[i].slave_y;
    buf_loc[cnt++] = vertices[i].master_x;
    buf_loc[cnt++] = vertices[i].master_y;
    buf_loc[cnt++] = vertices[i].slave_edge_id;
    buf_loc[cnt++] = vertices[i].master_edge_flag;
  }
}

template <typename DataType>
void ContactFaceFaceInteraction<DataType>::Restart_Unpack( DataType* buffer )
{
  int i;
  int cnt=0;
  DataType* buf_loc = buffer;

  cnt += this->UnPackEntityData(&slave_face_entity_data, &buf_loc[cnt]);
  cnt += this->UnPackEntityData(&master_face_entity_data, &buf_loc[cnt]);
  
  // Unpack the Data Array
  num_edges = (int) *buf_loc++;
  vertices = new ContactFaceFaceVertex<DataType>[num_edges+1];
  for (i=0; i<num_edges+1; ++i) {
    vertices[i].slave_x          = buf_loc[cnt++];
    vertices[i].slave_y          = buf_loc[cnt++];
    vertices[i].master_x         = buf_loc[cnt++];
    vertices[i].master_y         = buf_loc[cnt++];
    vertices[i].slave_edge_id    = (int) buf_loc[cnt++];
    vertices[i].master_edge_flag = (int) buf_loc[cnt++];
  }
}

template <typename DataType>
int ContactFaceFaceInteraction<DataType>::Set_SlaveFaceEntityData() 
{
  if (slave_face) {
    this->SetEntityData(&slave_face_entity_data, slave_face);
    return 1;
  }else{
    return 0;
  }
}

template <typename DataType>
int ContactFaceFaceInteraction<DataType>::Set_MasterFaceEntityData() 
{
  if (master_face) {
    this->SetEntityData(&master_face_entity_data, master_face);
    return 1;
  }else{
    return 0;
  }
}

#if (MAX_FFI_DERIVATIVES > 0)
template <>
void
ContactFaceFaceInteraction<Real>::Set_Derivatives( ContactFaceFaceInteraction<ActiveScalar> *active_cffi)
{
  num_derivatives = (active_cffi) ?
                    3*(slave_face->Nodes_Per_Face() + master_face->Nodes_Per_Face()) : 0;
  if (num_edges>0) {
    int i, j, k, l;
    for (i=0; i<num_edges; ++i) {
      if(num_derivatives > 0) {
        vertices[i].slave_x_derivatives = new Real[num_derivatives];
        vertices[i].slave_y_derivatives = new Real[num_derivatives];
        vertices[i].master_x_derivatives = new Real[num_derivatives];
        vertices[i].master_y_derivatives = new Real[num_derivatives];
        for(j = 0; j<3*slave_face->Nodes_Per_Face(); ++j) {
          vertices[i].slave_x_derivatives[j] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->slave_x, j);
          vertices[i].slave_y_derivatives[j] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->slave_y, j);
          vertices[i].master_x_derivatives[j] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->master_x, j);
          vertices[i].master_y_derivatives[j] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->master_y, j);
        }
        for(j = 0, k = 3*slave_face->Nodes_Per_Face(), l=MAX_FFI_DERIVATIVES/2; j<3*master_face->Nodes_Per_Face(); ++j, ++k, ++l) {
          vertices[i].slave_x_derivatives[k] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->slave_x, l);
          vertices[i].slave_y_derivatives[k] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->slave_y, l);
          vertices[i].master_x_derivatives[k] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->master_x, l);
          vertices[i].master_y_derivatives[k] = GetActiveScalarDerivative(active_cffi->Get_Vertex(i)->master_y, l);
        }
      }
    }
  }
}
#endif

#endif // ContactFaceFaceInteraction_C_
