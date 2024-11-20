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


#ifndef ContactCommBuffer_h_
#define ContactCommBuffer_h_

#include "Contact_Defines.h"
#include "Contact_Communication.h"

#ifndef CONTACT_NO_MPI
class ContactSymComm;
class ContactAsymComm;
#endif

class ContactCommBuffer {
  
 public:

  ContactCommBuffer(); 
  ~ContactCommBuffer();

#ifndef CONTACT_NO_MPI
  void Buffers( ContactSymComm&,    int, char**, char**, 
                RequestHandle**, RequestHandle**);
  void Buffers( ContactAsymComm&,   int, char**, char**, RequestHandle** );
  void Buffers_Export( ContactAsymComm&,   int, char**, char**, 
                       RequestHandle**, RequestHandle** );
  void Buffers_Import( ContactAsymComm&,   int, char**, char**, 
                       RequestHandle**, RequestHandle** );
  void Buffers( int size_send_buf, int size_recv_buff, int size_request_handle,
		char** send_buf, char** recv_buf, 
		RequestHandle** request_handles );
#endif

 private:
  
  ContactCommBuffer(ContactCommBuffer&);
  ContactCommBuffer& operator=(ContactCommBuffer&);

  int buffer_size;
  char* buffer;
  int handles_size;
  RequestHandle* handles;
};

#endif
