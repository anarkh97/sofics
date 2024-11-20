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


#include <iostream>

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "ContactParOStream.h"
#include "Contact_Communication.h"

#include <unistd.h>


ContactParOStream::ContactParOStream( MPI_Comm& communicator,
				      int orank,
				      char* fname )
  : out_rank(orank),
    com(communicator),
    buf_len(0),
    buf_size(0),
    buf_block_size(32768),
    buf(0)
{
  if (fname != 0)
  {
    std::strncpy(file_name, fname, MAX_FILE_NAME_SIZE - 1);
    file_name[MAX_FILE_NAME_SIZE - 1] = '\0';
  }
  else
    std::strcpy(file_name, "std::cout");
}

ContactParOStream::~ContactParOStream()
{
  if (buf != 0) delete [] buf;  buf = 0;
}


void ContactParOStream::append_to_buffer(const char* s)
{
  int size = std::strlen(s);
  if (size==0) return;
  assert(buf_len <= buf_size);
  
  if (buf_size - buf_len < size)
  {
    int factor    = ( (buf_len + size) / buf_block_size ) + 1;
    buf_size      = factor * buf_block_size;
    char* new_buf = new char[buf_size];
    if (buf != 0) {
      std::memcpy(new_buf, buf, buf_len);
      delete [] buf;
    }
    buf = new_buf;
  }
  
  std::memcpy(buf + buf_len, s, size);
  buf_len += size;
}

ContactParOStream& ContactParOStream::operator<<(const char* str)
{
  append_to_buffer(str);
  return *this;
}

ContactParOStream& ContactParOStream::operator<<(int i)
{
  std::sprintf(tmp_buf, "%d", i);
  append_to_buffer(tmp_buf);
  return *this;
}

ContactParOStream& ContactParOStream::operator<<(long i)
{
  std::sprintf(tmp_buf, "%ld", i);
  append_to_buffer(tmp_buf);
  return *this;
}

ContactParOStream& ContactParOStream::operator<<(unsigned int i)
{
  std::sprintf(tmp_buf, "%u", i);
  append_to_buffer(tmp_buf);
  return *this;
}

ContactParOStream& ContactParOStream::operator<<(unsigned long i)
{
  std::sprintf(tmp_buf, "%lu", i);
  append_to_buffer(tmp_buf);
  return *this;
}

// ContactParOStream& ContactParOStream::operator<<(long long i)
// {
//   std::sprintf(tmp_buf, "%lu", i);
//   append_to_buffer(tmp_buf);
//   return *this;
// }

ContactParOStream& ContactParOStream::operator<<(float f)
{
  std::sprintf(tmp_buf, "%g", f);
  append_to_buffer(tmp_buf);
  return *this;
}

ContactParOStream& ContactParOStream::operator<<(double d)
{
  std::sprintf(tmp_buf, "%g", d);
  append_to_buffer(tmp_buf);
  return *this;
}

ContactParOStream& ContactParOStream::operator<<(void* p)
{
  std::sprintf(tmp_buf, "%p", p);
  append_to_buffer(tmp_buf);
  return *this;
}

void ContactParOStream::output_buffer(char* str, int* sizes)
{
  std::cout.setf(std::ios::fixed);
  std::cout.precision(16);
  std::cout.flush();
  
  int offset    = 0;
  int num_procs = contact_number_of_processors(com);
  for (int p=0 ; p<num_procs ; ++p)
  {
    if (sizes[p] == 0) continue;
    
    char tmp_buf_tmp[TMP_BUF_SIZE_WITH_NULL];
    
    int b      = 0;
    tmp_buf_tmp[0] = '\0';
    char prefix[64];
    int  prefix_len = 0;
    
    if (num_procs>1) {
      std::sprintf(prefix, "P%d: ", p);
      prefix_len = std::strlen(prefix);
      std::strcat(tmp_buf_tmp, prefix);  
      b += prefix_len;
    } else {
      prefix_len = 0;
      prefix[0]  = '\0';
    }
    
    for (int i = 0; i < sizes[p]; ++i)
    {
      char c = str[offset + i];
      if (c == '\n' && i < sizes[p]-1)
      {
        assert( b <= TMP_BUF_SIZE );
        if (b + prefix_len + 1 >= TMP_BUF_SIZE) 
        { 
          std::cout << tmp_buf_tmp; 
          std::cout.flush(); 
          b = 0; 
        }
        tmp_buf_tmp[b++] = '\n';  
        tmp_buf_tmp[b  ] = '\0';
        std::strcat(tmp_buf_tmp, prefix);  
        b += prefix_len;
        assert( b <= TMP_BUF_SIZE );
      }
      else {
        assert( b <= TMP_BUF_SIZE );
        if (b == TMP_BUF_SIZE) 
        { 
          std::cout << tmp_buf_tmp; 
          std::cout.flush(); 
          b = 0; 
        }
        tmp_buf_tmp[b++] = c;  
        tmp_buf_tmp[b  ] = '\0';
        assert( b <= TMP_BUF_SIZE );
      }
    }
    
    // if the last char in a processor's buffer was not a
    // EOL, then append one here before the final output
    // of the processor's buffer
    if (str[offset + sizes[p] - 1] != '\n') {
      assert( b <= TMP_BUF_SIZE );
      if (b == TMP_BUF_SIZE) 
      { 
        std::cout << tmp_buf_tmp; 
        std::cout.flush(); 
        b = 0; 
      }
      tmp_buf_tmp[b++] = '\n';  
      tmp_buf_tmp[b  ] = '\0';
      assert( b <= TMP_BUF_SIZE );
    }
    if (b > 0) std::cout << tmp_buf_tmp;
    std::cout.flush();
    
    offset += sizes[p];
    std::cout.flush();
    //sleep(1); // boy this sure makes it slow!
  }
}

void ContactParOStream::flush()
{
  if( contact_number_of_processors(com) == 1 ){
    output_buffer( buf, &buf_len );
  } else {
    if( contact_processor_number(com) == out_rank )
    {
      int  com_size   = contact_number_of_processors(com);
      int* recv_sizes = new int[com_size];
      contact_gather( com, buf_len, recv_sizes, out_rank );
      
      int recv_size = 0;
      for (int i = 0; i < com_size; ++i) {
	recv_size += recv_sizes[i];
      }
      char* recv_buf = new char[recv_size+256];
      
      contact_hetro_gather( com, buf, buf_len, recv_buf, recv_sizes, out_rank);
      
      output_buffer(recv_buf, recv_sizes);
      
      delete [] recv_sizes;
      delete [] recv_buf;
    }
    else
    {
      contact_gather( com, buf_len, NULL, out_rank );
      contact_hetro_gather( com, buf, buf_len, NULL, NULL, out_rank );
    }
    contact_global_sync( com );
  }
  buf_len = 0;
  std::memset( buf, 0, buf_size );
}

