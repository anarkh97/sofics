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


#ifndef ContactParOStream_h_
#define ContactParOStream_h_

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#else
typedef int MPI_Comm;
#endif

class ContactParOStream {
public:
  
  ContactParOStream( MPI_Comm& communicator,
		     int orank = 0,
		     char* fname = 0 );
 ~ContactParOStream();
  
  ContactParOStream& operator<<(const char*);  // null terminated string
  ContactParOStream& operator<<(int);
  ContactParOStream& operator<<(long);
  ContactParOStream& operator<<(unsigned int);
  ContactParOStream& operator<<(unsigned long);
  ContactParOStream& operator<<(float);
  ContactParOStream& operator<<(double);
  ContactParOStream& operator<<(void*);        // prints an address
  
  void flush();  // implies a parallel synchronization among the communicator
  
private:
  
  ContactParOStream(ContactParOStream&);
  ContactParOStream& operator=(ContactParOStream&);
  
  enum { MAX_FILE_NAME_SIZE = 512 };
  char file_name[MAX_FILE_NAME_SIZE];
  
  int out_rank;
  MPI_Comm com;
  
  int   buf_len;         // number of chars of client data
  int   buf_size;
  int   buf_block_size;  // right now determined in constructor
  char* buf;
  
  enum { TMP_BUF_SIZE = 1024, TMP_BUF_SIZE_WITH_NULL};
  char tmp_buf[TMP_BUF_SIZE_WITH_NULL];
  
  void append_to_buffer(const char* s);
  void output_buffer(char* str, int* seg_sizes);
};

#endif
