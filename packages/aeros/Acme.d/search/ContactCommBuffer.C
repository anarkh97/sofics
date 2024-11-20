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


#include "ContactCommBuffer.h"
#include "ContactSymComm.h"
#include "ContactAsymComm.h"

#include <cstddef>

ContactCommBuffer::ContactCommBuffer()
{
  buffer_size = 0;
  buffer = NULL;
  handles_size = 0;
  handles = NULL;
}

ContactCommBuffer::~ContactCommBuffer()
{
  if( buffer ) delete [] buffer;
  if( handles ) delete [] handles;
}

#ifndef CONTACT_NO_MPI
void ContactCommBuffer::Buffers( ContactSymComm& Comm, int size_per_object,
				 char** send_buffer, char** recv_buffer,
				 RequestHandle** recv_handles, RequestHandle** send_handles )
{
  int bsize = 2*size_per_object*(Comm.Size());
  if( bsize > buffer_size ){
    if( buffer ) delete [] buffer;
    buffer_size = (int) 1.25*bsize;
    buffer = new char[buffer_size];
  }
  *send_buffer = buffer;
  *recv_buffer = buffer+size_per_object*(Comm.Size());
  int hsize = Comm.Num_Comm_Partners() * 2;
  if( hsize > handles_size ){
    if( handles ) delete [] handles;
    handles_size = (int) (1.25*hsize+1);
    handles = new RequestHandle[handles_size];
  }
  *recv_handles = handles;
  *send_handles = handles + Comm.Num_Comm_Partners();
}

void ContactCommBuffer::Buffers( ContactAsymComm& Comm, int size_per_object,
				 char** send_buffer, char** recv_buffer,
				 RequestHandle** recv_handles )
{
  int bsize = size_per_object*(Comm.Size_Import() + Comm.Size_Export());
  if( bsize > buffer_size ){
    if( buffer ) delete [] buffer;
    buffer_size = (int) 1.25*bsize;
    buffer = new char[buffer_size];
  }
  *send_buffer = buffer;
  *recv_buffer = buffer+size_per_object*Comm.Size_Export();
  int hsize = Comm.Num_Import_Comm_Partners();
  if( hsize > handles_size ){
    if( handles ) delete [] handles;
    handles_size = (int) (1.25*hsize+1);
    handles = new RequestHandle[handles_size];
  }
  *recv_handles = handles;
}

void ContactCommBuffer::Buffers_Import( ContactAsymComm& Comm, int size_per_object,
				 char** send_buffer, char** recv_buffer,
				 RequestHandle** recv_handles, RequestHandle** send_handles )
{
  int bsize = size_per_object*(Comm.Size_Import() + Comm.Size_Export());
  if( bsize > buffer_size ){
    if( buffer ) delete [] buffer;
    buffer_size = (int) 1.25*bsize;
    buffer = new char[buffer_size];
  }
  *send_buffer = buffer;
  *recv_buffer = buffer+size_per_object*Comm.Size_Export();
  int hsize = Comm.Num_Import_Comm_Partners() + Comm.Num_Export_Comm_Partners();
  if( hsize > handles_size ){
    if( handles ) delete [] handles;
    handles_size = (int) (1.25*hsize+1);
    handles = new RequestHandle[handles_size];
  }
  *recv_handles = handles;
  *send_handles = handles + Comm.Num_Import_Comm_Partners();
}

void ContactCommBuffer::Buffers_Export( ContactAsymComm& Comm, int size_per_object,
				 char** send_buffer, char** recv_buffer,
				 RequestHandle** recv_handles, RequestHandle** send_handles )
{
  int bsize = size_per_object*(Comm.Size_Import() + Comm.Size_Export());
  if( bsize > buffer_size ){
    if( buffer ) delete [] buffer;
    buffer_size = (int) 1.25*bsize;
    buffer = new char[buffer_size];
  }
  *send_buffer = buffer;
  *recv_buffer = buffer+size_per_object*Comm.Size_Import();
  int hsize = Comm.Num_Export_Comm_Partners() + Comm.Num_Import_Comm_Partners();
  if( hsize > handles_size ){
    if( handles ) delete [] handles;
    handles_size = (int) (1.25*hsize+1);
    handles = new RequestHandle[handles_size];
  }
  *recv_handles = handles;
  *send_handles = handles + Comm.Num_Export_Comm_Partners();
}

void ContactCommBuffer::Buffers( int size_send_buf, int size_recv_buf,
				 int size_request_handle,
				 char** send_buf, char** recv_buf,
				 RequestHandle** request_handles )
{
  int bsize = size_send_buf + size_recv_buf;
  if( bsize > buffer_size ){
    if( buffer ) delete [] buffer;
    buffer_size = (int) 1.25*bsize;
    buffer = new char[buffer_size];
  }
  *send_buf = buffer;
  *recv_buf = buffer + size_send_buf;
  int hsize = size_request_handle;
  if( hsize > handles_size ){
    if( handles ) delete [] handles;
    handles_size = (int) (1.25*hsize+1);
    handles = new RequestHandle[handles_size];
  }
  *request_handles = handles;
}
#endif
