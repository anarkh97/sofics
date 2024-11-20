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


#include "Contact_Communication.h"
#include "contact_assert.h"
#include "ContactCommBuffer.h"
#include "ContactAsymComm.h"
#include "ContactSymComm.h"
#include "ContactParOStream.h"
#include "ContactSearch.h"
#include "ContactTopologyEntity.h"
#include "ContactEdge.h"
#include "ContactFace.h"
#include "ContactScratchManager.h"

#include <cstdlib>
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <vector>

using namespace std;


int contact_processor_number( MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int myrank;
  MPI_Comm_rank( communicator, &myrank );
  return( myrank );
#else
  return( 0 );
#endif
}

int contact_any_processor_number()
{
  return -1;
}

int contact_any_message_number()
{
  return -1;
}

int contact_number_of_processors( MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int size;
  MPI_Comm_size( communicator, &size );
  return size;
#else
  return 1;
#endif
}


void contact_global_sync( MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Barrier( communicator );
  if( ierr ) std::cerr << "Error in contact_global_sync" << std::endl;
#else
  MPI_Barrier( communicator );
#endif
#endif
}

void contact_print_sync_start( MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int rank = contact_processor_number( communicator );
  if (rank!=0) {
    int flag=1, mesg=10000;
    MPI_Status status;
    int src  = rank-1;
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    int ierr = MPI_Recv(&flag, 1, MPI_INT, src, mesg, communicator, &status);
    if( ierr ) std::cerr << "Error in contact_print_sync_start" << std::endl;
#else
    MPI_Recv(&flag, 1, MPI_INT, src, mesg, communicator, &status);
#endif
  }
#endif
}

void contact_print_sync_end( MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int flag=1, mesg=10000, src, dest;
  MPI_Status status;
  int nprocs = contact_number_of_processors( communicator );
  int rank   = contact_processor_number( communicator );

  std::cout<<std::flush;sleep(2);
  if (rank<nprocs-1) {
    dest = rank+1;
  } else {
    dest = 0;
  }
  int ierr = MPI_Send(&flag, 1, MPI_INT, dest, mesg, communicator);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if( ierr ) std::cerr << "Error in contact_print_sync_end" << std::endl;
#endif
  if (rank==0) {
    src  = nprocs-1;
    ierr = MPI_Recv(&flag, 1, MPI_INT, src, 
                    mesg, communicator, &status);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    if( ierr ) std::cerr << "Error in contact_print_sync_end" << std::endl;
#endif
  }
  contact_global_sync(communicator);
#endif
}

int contact_global_sum( int value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int result;
  MPI_Allreduce( &value, &result, 1, MPI_INT, MPI_SUM, communicator );
  return result;
#else
  return value;
#endif
}

int contact_global_error_check( int value, MPI_Comm& communicator )
{
  //
  //  NKC 2/23/05
  //
  //  The all reduce error check takes a non-trivial ammount of time.  Would
  //  like to optinally replace this with something that does not require 
  //  communication.
  //  Below is the code from sierra, need to tweak a bit more to ensure
  //  that it functions correctly.  Possibly make this available as an
  //  option?  Need to do some more profiling as well to determine if these
  //  operations really do take significant quantities of time.  Entirely
  //  possible that in the 1000 processor range all-reduce operations could
  //  be quite slow.
  //
  //#define NON_SYNCED_ERRORS 
  //#ifdef NON_SYNCED_ERRORS
  //if(value != ContactSearch::NO_ERROR) {
  //  std::cout<<"exiting ACME with error value "<<value<<std::endl;
  //  MPI_Abort( MPI_COMM_WORLD , MPI_ERR_OTHER ); // First  try to die
  //  std::exit( EXIT_FAILURE );                        // Second try to die
  //  std::abort();                                     // Final  try to die
  //  return value;
  //}
  //#else
#ifndef CONTACT_NO_MPI
  int result;
  MPI_Allreduce( &value, &result, 1, MPI_INT, MPI_MAX, communicator );
  return result;
#else
  return value;
#endif
  //#endif
}

Real contact_global_sum( Real value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
  Real result;
  MPI_Allreduce( &value, &result, 1, MPI_DOUBLE, MPI_SUM, communicator );
  return result;
#else
  return value;
#endif
}

void contact_global_sum( Real* ivalues, Real* ovalues, int nvalues, 
			 MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Allreduce( ivalues, ovalues, nvalues, MPI_DOUBLE, MPI_SUM, communicator);
#else
  std::memcpy( ovalues, ivalues, nvalues*sizeof(Real) );
#endif
}

void contact_global_sum( int* ivalues, int* ovalues, int nvalues, 
			 MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Allreduce( ivalues, ovalues, nvalues, MPI_INT, MPI_SUM, communicator);
#else
  std::memcpy( ovalues, ivalues, nvalues*sizeof(int) );
#endif
}

int contact_global_minimum( int value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int result;
  MPI_Allreduce( &value, &result, 1, MPI_INT, MPI_MIN, communicator );
  return result;
#else
  return value;
#endif
}

Real contact_global_minimum( Real value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
  Real result;
  MPI_Allreduce( &value, &result, 1, MPI_DOUBLE, MPI_MIN, communicator );
  return result;
#else
  return value;
#endif
}

void contact_global_minimum( int* ivalues, int* ovalues, int nvalues,
			     MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Allreduce( ivalues, ovalues, nvalues, MPI_INT, MPI_MIN, communicator );
#else
  for (int i=0; i < nvalues; ++i)
    ovalues[i] = ivalues[i];
#endif
}

void contact_global_minimum( Real* ivalues, Real* ovalues, int nvalues,
			     MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
  MPI_Allreduce( ivalues, ovalues, nvalues, MPI_DOUBLE, MPI_MIN, communicator );
#else
  for (int i=0; i < nvalues; ++i)
    ovalues[i] = ivalues[i];
#endif
}


int contact_global_maximum( int value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int result;
  MPI_Allreduce( &value, &result, 1, MPI_INT, MPI_MAX, communicator );
  return result;
#else
  return value;
#endif
}

Real contact_global_maximum( Real value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
  Real result;
  MPI_Allreduce( &value, &result, 1, MPI_DOUBLE, MPI_MAX, communicator );
  return result;
#else
  return value;
#endif
}


bool contact_global_or( bool value, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int argument = value;
  int result;
  MPI_Allreduce( &argument, &result, 1, MPI_INT, MPI_MAX, communicator );
  return result;
#else
  return value;
#endif
}


int contact_broadcast( int * array, int size, int sender, 
		       MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int ierr;
  ierr = MPI_Bcast( array, size, MPI_INT, sender, communicator );
  return ierr;
#else
  return 0;
#endif
}

int contact_broadcast( Real * array, int size, int sender, 
		       MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int ierr;
  ierr = MPI_Bcast( array, size, MPI_DOUBLE, sender, communicator );
  return ierr;
#else
  return 0;
#endif
}

int contact_broadcast( char * array, int size, int sender, 
		       MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  int ierr;
  ierr = MPI_Bcast( array, size, MPI_CHAR, sender, communicator );
  return ierr;
#else
  return 0;
#endif
}

void contact_global_maximum( int* ivalues, int* ovalues, int nvalues,
			     MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Allreduce( ivalues, ovalues, nvalues, MPI_INT, MPI_MAX, communicator );
#else
  for (int i=0; i < nvalues; ++i)
    ovalues[i] = ivalues[i];
#endif
}

void contact_global_maximum( Real* ivalues, Real* ovalues, int nvalues,
			     MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
  MPI_Allreduce( ivalues, ovalues, nvalues, MPI_DOUBLE, MPI_MAX, communicator );
#else
  for (int i=0; i < nvalues; ++i)
    ovalues[i] = ivalues[i];
#endif
}

RequestHandle contact_nonblocking_receive( int mesg, int* data, int len, 
					   int proc, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Request request = MPI_REQUEST_NULL;
  int isource = proc;
  if( proc == -1 ) isource = MPI_ANY_SOURCE;
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Irecv( data, len, MPI_INT, isource, mesg, communicator, &request );
  if( ierr ) std::cerr << "Error in contact_nonblocking_receive" << std::endl;
#else
  MPI_Irecv( data, len, MPI_INT, isource, mesg, communicator, &request );
#endif
  return request;
#else
  return -1;
#endif
}

RequestHandle contact_nonblocking_receive( int mesg, Real* data, int len, 
					   int proc, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
  MPI_Request request = MPI_REQUEST_NULL;
  int isource = proc;
  if( proc == -1 ) isource = MPI_ANY_SOURCE;
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Irecv( data, len, MPI_DOUBLE, isource, mesg, communicator, &request );
  if( ierr ) std::cerr << "Error in contact_nonblocking_receive" << std::endl;
#else
  MPI_Irecv( data, len, MPI_DOUBLE, isource, mesg, communicator, &request );
#endif
  return request;
#else
  return -1;
#endif
}

RequestHandle contact_nonblocking_receive( int mesg, char* data, int len, 
					   int proc, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Request request = MPI_REQUEST_NULL;
  int isource = proc;
  if( proc == -1 ) isource = MPI_ANY_SOURCE;
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Irecv( data, len, MPI_CHAR, isource, mesg, communicator, &request );
  if( ierr ) std::cerr << "Error in contact_nonblocking_receive" << std::endl;
#else
  MPI_Irecv( data, len, MPI_CHAR, isource, mesg, communicator, &request );
#endif
  return request;
#else
  return -1;
#endif
}


RequestHandle contact_nonblocking_send( int mesg, Real* data, int len, 
					int proc, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Request request(MPI_REQUEST_NULL);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Isend( data, len, MPI_DOUBLE, proc, mesg, communicator, &request );
  if( ierr ) std::cerr << "Error in contact_nonblocking_send" << std::endl;
#else
  MPI_Isend( data, len, MPI_DOUBLE, proc, mesg, communicator, &request );
#endif
  return request;
#else
  return -1;
#endif
}

RequestHandle contact_nonblocking_send( int mesg, int* data, int len, 
					int proc, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Request request(MPI_REQUEST_NULL);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Isend( data, len, MPI_INT, proc, mesg, communicator, &request );
  if( ierr ) std::cerr << "Error in contact_nonblocking_send" << std::endl;
#else
  MPI_Isend( data, len, MPI_INT, proc, mesg, communicator, &request );
#endif
  return request;
#else
  return -1;
#endif
}

RequestHandle contact_nonblocking_send( int mesg, char* data, int len, 
					int proc, MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  MPI_Request request(MPI_REQUEST_NULL);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Isend( data, len, MPI_CHAR, proc, mesg, communicator, &request );
  if( ierr ) std::cerr << "Error in contact_nonblocking_send" << std::endl;
#else
  MPI_Isend( data, len, MPI_CHAR, proc, mesg, communicator, &request );
#endif
  return request;
#else
  return -1;
#endif
}

void contact_blocking_send( int mesg, int* data, int len, int proc,
			    MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Send( data, len, MPI_INT, proc, mesg, communicator );
  if( ierr ) std::cerr << "Error in contact_blocking_send" << std::endl;
#else
  MPI_Send( data, len, MPI_INT, proc, mesg, communicator );
#endif
#endif
}


void contact_blocking_send( int mesg, Real* data, int len, int proc,
			    MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
  PRECONDITION( sizeof(Real) == sizeof(double) );
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Send( data, len, MPI_DOUBLE, proc, mesg, communicator );
  if( ierr ) std::cerr << "Error in contact_blocking_send" << std::endl;
#else
  MPI_Send( data, len, MPI_DOUBLE, proc, mesg, communicator );
#endif
#endif
}

void contact_blocking_send( int mesg, char* data, int len, int proc,
			    MPI_Comm& communicator )
{
#ifndef CONTACT_NO_MPI
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Send( data, len, MPI_CHAR, proc, mesg, communicator );
  if( ierr ) std::cerr << "Error in contact_blocking_send" << std::endl;
#else
  MPI_Send( data, len, MPI_CHAR, proc, mesg, communicator );
#endif
#endif
}


int contact_wait_msg_done( RequestHandle request )
{
#ifndef CONTACT_NO_MPI
  MPI_Status status;
  status.MPI_SOURCE = 0;
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Wait( &request, &status );
  if( ierr ) std::cerr << "Error in contact_wait_msg_done" << std::endl;
#else
  MPI_Wait( &request, &status );
#endif
  return status.MPI_SOURCE;
#else
  return -1;
#endif
}

int contact_test_msg_done( RequestHandle request )
{
#ifndef CONTACT_NO_MPI
  MPI_Status status;
  int flag;
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int ierr = MPI_Test( &request, &flag, &status );
  if( ierr ) std::cerr << "Error in contact_test_msg_done" << std::endl;
#else
  MPI_Test( &request, &flag, &status );
#endif
  return flag;
#else
  return -1;
#endif
}

void contact_swapadd_reg_data( MPI_Comm& search_comm,  
			       ContactSymComm& comm,
			       ContactCommBuffer& comm_buffer,
			       VariableHandle var,
			       int size_of_var )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10002;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = size_of_var*sizeof(Real);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, 
                      &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy(sb, entity_list[j]->Variable(var), bytes_per_var );
      sb += size_of_var;
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = size_of_var*comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
				   	           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_of_var*comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = size_of_var*comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
#ifndef CONTACT_USE_BLOCKING_SEND
    contact_wait_msg_done(send_handles[i]);
#endif
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy(sb, entity_list[j]->Variable(var), bytes_per_var );
      std::memset(entity_list[j]->Variable(var), 0, bytes_per_var );
      sb += size_of_var;
    }
  }  

  // Assemble the data back to contact entity in ascending processor order
  bool added_mine = false;
  rb = recv_buf;
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    //
    // Wait till current message has arrived
    //
    contact_wait_msg_done( recv_handles[i] );
    if( !added_mine && contact_processor_number(search_comm)<comm.Comm_Proc_ID( i ) ){
      // Add in my data
      added_mine = true;
      for(int j=0 ; j<num_procs ; ++j ){
	int num_to_proc = comm.Num_to_Proc( j );
	ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
	for(int k=0 ; k<num_to_proc ; ++k ){
	  Real *vardata = entity_list[k]->Variable(var);
	  for(int l=0 ; l<size_of_var ; ++l )
	    vardata[l] += *sb++;
	}
      }
    }
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      Real *vardata = entity_list[j]->Variable(var);
      for(int k=0 ; k<size_of_var ; ++k )
	vardata[k] += *rb++;
    }
  }

  // If my data hasn't been added, add it
  if( !added_mine ){
    for(int i=0 ; i<num_procs ; ++i ){
      int num_to_proc = comm.Num_to_Proc( i );
      ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
      for(int j=0 ; j<num_to_proc ; ++j ){
	Real *vardata = entity_list[j]->Variable(var);
	for(int k=0 ; k<size_of_var ; ++k )
	  vardata[k] += *sb++;
      }
    }
  }
#endif
}

void contact_swapadd_temp_tag( MPI_Comm& search_comm,  
			       ContactSymComm& comm,
			       ContactCommBuffer& comm_buffer )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10002;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = sizeof(int);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, 
                      &crcv_buf, &recv_handles, &send_handles);
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<int*>(csnd_buf);
  rb = recv_buf = reinterpret_cast<int*>(crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      *sb++ = entity_list[j]->temp_tag;
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
				   	           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
#ifndef CONTACT_USE_BLOCKING_SEND
    contact_wait_msg_done(send_handles[i]);
#endif
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      *sb++ = entity_list[j]->temp_tag;
      entity_list[j]->temp_tag = 0;
    }
  }  

  // Assemble the data back to contact entity in ascending processor order
  bool added_mine = false;
  rb = recv_buf;
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    //
    // Wait till current message has arrived
    //
    contact_wait_msg_done( recv_handles[i] );
    if( !added_mine && contact_processor_number(search_comm)<comm.Comm_Proc_ID( i ) ){
      // Add in my data
      added_mine = true;
      for(int j=0 ; j<num_procs ; ++j ){
	int num_to_proc = comm.Num_to_Proc( j );
	ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
	for(int k=0 ; k<num_to_proc ; ++k ){
          entity_list[k]->temp_tag += *sb++;
	}
      }
    }
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      entity_list[j]->temp_tag += *rb++;
    }
  }

  // If my data hasn't been added, add it
  if( !added_mine ){
    for(int i=0 ; i<num_procs ; ++i ){
      int num_to_proc = comm.Num_to_Proc( i );
      ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
      for(int j=0 ; j<num_to_proc ; ++j ){
	entity_list[j]->temp_tag += *sb++;
      }
    }
  }
#endif
}

void contact_swapadd_context( MPI_Comm& search_comm,  
			      ContactSymComm& comm,
			      ContactCommBuffer& comm_buffer )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10002;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = sizeof(int);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, 
                      &crcv_buf, &recv_handles, &send_handles);
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<int*>(csnd_buf);
  rb = recv_buf = reinterpret_cast<int*>(crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      *sb++ = entity_list[j]->GetContext();
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
				   	           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

#ifndef CONTACT_USE_BLOCKING_SEND
  for(int i=0 ; i<num_procs ; ++i ){
    contact_wait_msg_done(send_handles[i]);
  }  
#endif

  // Assemble the data back to contact entity in ascending processor order
  rb = recv_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    //
    // Wait till current message has arrived
    //
    contact_wait_msg_done( recv_handles[i] );
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      unsigned int new_context = (*rb++) & (ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE|ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
      unsigned int old_context = entity_list[j]->GetContext();
      unsigned int context     = old_context|new_context;
      entity_list[j]->SetContext(context);
    }
  }

#endif
}

void contact_swap_edge_faces( MPI_Comm& search_comm,  
			      ContactSymComm& comm,
			      ContactCommBuffer& comm_buffer )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10003;
  int len,num_to_proc;
  int var_len = 4;
  ContactTopologyEntity<Real>** entity_list;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = var_len*sizeof(int);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, &crcv_buf, 
                      &recv_handles, &send_handles);
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<int*>(csnd_buf);
  rb = recv_buf = reinterpret_cast<int*>(crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    num_to_proc = comm.Num_to_Proc( i );
    entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity_list[j]);
      PRECONDITION(edge);
      *sb++ = edge->Face(0)->Owner();
      *sb++ = edge->Face(0)->BlockID();
      *sb++ = edge->Face(0)->Global_ID().HiInt();
      *sb++ = edge->Face(0)->Global_ID().LoInt();
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    len = var_len*comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += var_len*comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    len = var_len*comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Assemble the data back to contact entity in ascending processor order
  rb = recv_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    contact_wait_msg_done( recv_handles[i] );
    num_to_proc = comm.Num_to_Proc( i );
    entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity_list[j]);
      PRECONDITION(edge);
      edge->FaceInfo()->owner                  = *rb++;
      edge->FaceInfo()->owner_proc_array_index = *rb++;
      edge->FaceInfo()->host_gid[0]            = *rb++;
      edge->FaceInfo()->host_gid[1]            = *rb++;
    }
  }

#ifndef CONTACT_USE_BLOCKING_SEND
  // Wait till all of my messages have arrived
  for(int i=0 ; i<num_procs ; ++i ){
    contact_wait_msg_done( send_handles[i] );
  }
#endif

#endif
}

void contact_reduce_add_scratch_var( MPI_Comm& search_comm,  
				     ContactAsymComm& comm,
				     ContactCommBuffer& comm_buffer,
				     ScratchVariable &var)
{

#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10005;
  int my_proc = contact_processor_number(search_comm);
  //
  //  Step one, copy and add all ghost copies of the data onto the localy owned master copy
  //
  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  //
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = var.Get_Size()*sizeof(Real);

  comm_buffer.Buffers_Export(comm, bytes_per_var, &csnd_buf, &crcv_buf, 
                             &recv_handles, &send_handles);

  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);
  //
  // Gather my data into send_buf
  //
  int *data_send_count = comm.Num_Import_from_Procs();
  int *data_receive_count = comm.Num_Export_to_Procs();
  int *import_index_list = comm.Import_Index_List();
  int total_import = comm.Size_Import();

  for(int i = 0; i < total_import; ++i) {
    std::memcpy(sb, var.Get_Scratch(import_index_list[i]), bytes_per_var);
    sb += var.Get_Size();
  }

  // Post Receives
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    int len = var.Get_Size() * data_receive_count[i];
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
				  	           comm.Export_Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < comm.Num_Import_Comm_Partners() && comm.Import_Comm_Proc_ID(proc_index) < my_proc) {
    sb += var.Get_Size()*data_send_count[proc_index];
    proc_index++;
  }
  if (proc_index == comm.Num_Import_Comm_Partners()) {
    sb = send_buf;
    proc_index = 0;
  }
  for(int i = proc_index, ipart = 0 ; ipart<comm.Num_Import_Comm_Partners() ; ++ipart ){
    int len = var.Get_Size() * data_send_count[i];
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Import_Comm_Proc_ID( i ),
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, 
                           comm.Import_Comm_Proc_ID( i ), 
                           search_comm );
#endif
    sb += len;
    if (++i == comm.Num_Import_Comm_Partners()) {
      i  = 0;
      sb = send_buf;
    }
  }

  rb = recv_buf;

  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    // Wait till current message has arrived
    contact_wait_msg_done( recv_handles[i] );
    //
    //  Add the data onto the master copy
    //
    int num_from_proc = data_receive_count[i];

    int* export_index_list = comm.Export_Index_List( i );

    for(int j=0 ; j<num_from_proc ; ++j ){
      Real *vardata = var.Get_Scratch(export_index_list[j]);
      for(int k = 0; k < var.Get_Size(); ++k) {    
        vardata[k] += rb[k];
      }
      rb += var.Get_Size();
    }
  }
  
#ifndef CONTACT_USE_BLOCKING_SEND
  //
  //  Assure that all sends have completed before reusing the buffers
  //
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    contact_wait_msg_done( send_handles[i] );
  }
#endif

  //
  //  Step 2, copy the reduced variables to all ghost node copies
  //
  contact_import_scratch_var(search_comm,
                             comm,
                             comm_buffer,
                             var);
#endif
}



void contact_reduce_add_scratch_vars( MPI_Comm& search_comm,  
				      ContactAsymComm& comm,
				      ContactCommBuffer& comm_buffer,
				      int num_vars,
				      ScratchVariable** vars)
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10006;
  int my_proc = contact_processor_number(search_comm);

  //
  //  Step one, copy and add all ghost copies of
  //  the data onto the localy owned master copy
  //

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;

  int *var_size = new int[num_vars];
  int *var_bytes = new int[num_vars];

  int size_of_all_vars = 0;
  for(int i = 0; i < num_vars; ++i) {
    size_of_all_vars += vars[i]->Get_Size();
    var_size[i] = vars[i]->Get_Size();
    var_bytes[i] = var_size[i] * sizeof(Real);
  }
  int bytes_per_entity = size_of_all_vars*sizeof(Real);
  comm_buffer.Buffers_Export(comm, bytes_per_entity, &csnd_buf, 
                             &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);
  // Gather my data into send_buf

  int *import_index_list = comm.Import_Index_List();
  int total_import = comm.Size_Import();

  for(int i=0 ; i<total_import ; ++i ){
    for(int j=0 ; j<num_vars ; ++j ){
      std::memcpy(sb, vars[j]->Get_Scratch(import_index_list[i]), var_bytes[j] );
      sb += var_size[j];
    }
  }


  // Post Receives
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    int len = size_of_all_vars*comm.Num_Export_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
					           comm.Export_Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < comm.Num_Import_Comm_Partners() && comm.Import_Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_of_all_vars*comm.Num_Import_from_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == comm.Num_Import_Comm_Partners()) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<comm.Num_Import_Comm_Partners() ; ++j ){
    int len = size_of_all_vars*comm.Num_Import_from_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Import_Comm_Proc_ID( i ),
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, 
                           comm.Import_Comm_Proc_ID( i ), 
                           search_comm );
#endif
    sb += len;
    if (++i == comm.Num_Import_Comm_Partners()) {
      i  = 0;
      sb = send_buf;
    }
  }

  rb = recv_buf;

  int* data_receive_count = comm.Num_Export_to_Procs();
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    // Wait till current message has arrived
    contact_wait_msg_done( recv_handles[i] );
    //
    //  Add the data onto the master copy
    //
    int num_from_proc = data_receive_count[i];
    int* export_index_list = comm.Export_Index_List( i );
    for(int j=0 ; j<num_from_proc ; ++j ){
      for(int k=0 ; k<num_vars ; ++k ){
        Real *vardata = vars[k]->Get_Scratch(export_index_list[j]);
	for(int l=0 ; l<var_size[k] ; ++l ) {
	  vardata[l] += rb[l];
	}
        rb += var_size[k];
      }
    }
  }
  
#ifndef CONTACT_USE_BLOCKING_SEND
  //
  //  Assure that all sends have completed before reusing the buffers
  //
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    contact_wait_msg_done( send_handles[i] );
  }
#endif

  delete [] var_size;
  delete [] var_bytes;

  //
  //  Step 2, copy the reduced variables to all ghost node copies
  //
  contact_import_scratch_vars(search_comm,
                              comm,
                              comm_buffer,
			      num_vars,
                              vars);
#endif
}


void contact_swapadd_data_array( MPI_Comm& search_comm,  
				 ContactSymComm& comm,
				 ContactCommBuffer& comm_buffer,
				 int* data,
				 int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10007;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = size_per_entity*sizeof(int);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, &crcv_buf, 
                      &recv_handles, &send_handles);
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<int*>(csnd_buf);
  rb = recv_buf = reinterpret_cast<int*>(crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy(sb, &data[entity_list[j]->ProcArrayIndex()*size_per_entity],
	     bytes_per_var );
      sb += size_per_entity;
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = size_per_entity*comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_per_entity*comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = size_per_entity*comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
#ifndef CONTACT_USE_BLOCKING_SEND
    contact_wait_msg_done(send_handles[i]);
#endif
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      int *idata = &data[entity_list[j]->ProcArrayIndex()*size_per_entity];
      std::memcpy(sb, idata, bytes_per_var );
      std::memset(idata, 0, bytes_per_var );
      sb += size_per_entity;
    }
  }  

  // Assemble the data back to contact entity in ascending processor order
  bool added_mine = false;
  rb = recv_buf;
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    // Wait untill the current message has arrived
    contact_wait_msg_done( recv_handles[i] );

    if( !added_mine && 
	contact_processor_number(search_comm)<comm.Comm_Proc_ID( i ) ){
      // Add in my data
      added_mine = true;
      for(int j=0 ; j<num_procs ; ++j ){
	int num_to_proc = comm.Num_to_Proc( j );
	ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
	for(int k=0 ; k<num_to_proc ; ++k ){
	  int *idata = &data[entity_list[k]->ProcArrayIndex()*size_per_entity];
	  for(int l=0 ; l<size_per_entity ; ++l )
	    idata[l] += *sb++;
	}
      }
    }
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      int *idata = &data[entity_list[j]->ProcArrayIndex()*size_per_entity];
      for(int k=0 ; k<size_per_entity ; ++k )
	idata[k] += *rb++;
    }
  }

  // If my data hasn't been added, add it
  if( !added_mine ){
    for(int j=0 ; j<num_procs ; ++j ){
      int num_to_proc = comm.Num_to_Proc( j );
      ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
      for(int k=0 ; k<num_to_proc ; ++k ){
	int *idata = &data[entity_list[k]->ProcArrayIndex()*size_per_entity];
	for(int l=0 ; l<size_per_entity ; ++l )
	  idata[l] += *sb++;
      }
    }
  }

#endif
}


void contact_swapadd_data_array_fcs( MPI_Comm& search_comm,  
				     ContactSymComm& comm,
				     ContactCommBuffer& comm_buffer,
				     int* data,
				     int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10007;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = size_per_entity*sizeof(int);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, &crcv_buf, 
                      &recv_handles, &send_handles);
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<int*>(csnd_buf);
  rb = recv_buf = reinterpret_cast<int*>(crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy(sb, &data[entity_list[j]->fcs_index*size_per_entity],
	     bytes_per_var );
      sb += size_per_entity;
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = size_per_entity*comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_per_entity*comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = size_per_entity*comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
#ifndef CONTACT_USE_BLOCKING_SEND
    contact_wait_msg_done(send_handles[i]);
#endif
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      int *idata = &data[entity_list[j]->fcs_index*size_per_entity];
      std::memcpy(sb, idata, bytes_per_var );
      std::memset(idata, 0, bytes_per_var );
      sb += size_per_entity;
    }
  }  

  // Assemble the data back to contact entity in ascending processor order
  bool added_mine = false;
  rb = recv_buf;
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    // Wait untill the current message has arrived
    contact_wait_msg_done( recv_handles[i] );

    if( !added_mine && 
	contact_processor_number(search_comm)<comm.Comm_Proc_ID( i ) ){
      // Add in my data
      added_mine = true;
      for(int j=0 ; j<num_procs ; ++j ){
	int num_to_proc = comm.Num_to_Proc( j );
	ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
	for(int k=0 ; k<num_to_proc ; ++k ){
	  int *idata = &data[entity_list[k]->fcs_index*size_per_entity];
	  for(int l=0 ; l<size_per_entity ; ++l )
	    idata[l] += *sb++;
	}
      }
    }
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      int *idata = &data[entity_list[j]->fcs_index*size_per_entity];
      for(int k=0 ; k<size_per_entity ; ++k )
	idata[k] += *rb++;
    }
  }

  // If my data hasn't been added, add it
  if( !added_mine ){
    for(int j=0 ; j<num_procs ; ++j ){
      int num_to_proc = comm.Num_to_Proc( j );
      ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
      for(int k=0 ; k<num_to_proc ; ++k ){
	int *idata = &data[entity_list[k]->fcs_index*size_per_entity];
	for(int l=0 ; l<size_per_entity ; ++l )
	  idata[l] += *sb++;
      }
    }
  }

#endif
}



void contact_swapadd_data_array( MPI_Comm& search_comm,  
				 ContactSymComm& comm,
				 ContactCommBuffer& comm_buffer,
				 Real* data,
				 int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10008;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = size_per_entity*sizeof(Real);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy(sb, &data[entity_list[j]->ProcArrayIndex()*size_per_entity],
	     bytes_per_var );
      sb += size_per_entity;
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = size_per_entity*comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_per_entity*comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = size_per_entity*comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
#ifndef CONTACT_USE_BLOCKING_SEND
    contact_wait_msg_done(send_handles[i]);
#endif
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      Real *rdata = &data[entity_list[j]->ProcArrayIndex()*size_per_entity];
      std::memcpy(sb, rdata, bytes_per_var );
      std::memset(rdata, 0, bytes_per_var );
      sb += size_per_entity;
    }
  }  

  // Assemble the data back to contact entity in ascending processor order
  bool added_mine = false;
  rb = recv_buf;
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    contact_wait_msg_done( recv_handles[i] );
    if( !added_mine && 
	contact_processor_number(search_comm)<comm.Comm_Proc_ID( i ) ){
      // Add in my data
      added_mine = true;
      for(int j=0 ; j<num_procs ; ++j ){
	int num_to_proc = comm.Num_to_Proc( j );
	ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
	for(int k=0 ; k<num_to_proc ; ++k ){
	  Real *rdata = &data[entity_list[k]->ProcArrayIndex()*size_per_entity];
	  for(int l=0 ; l<size_per_entity ; ++l )
	    rdata[l] += *sb++;
	}
      }
    }
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      Real *rdata = &data[entity_list[j]->ProcArrayIndex()*size_per_entity];
      for(int k=0 ; k<size_per_entity ; ++k )
	rdata[k] += *rb++;
    }
  }

  // If my data hasn't been added, add it
  if( !added_mine ){
    for(int j=0 ; j<num_procs ; ++j ){
      int num_to_proc = comm.Num_to_Proc( j );
      ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
      for(int k=0 ; k<num_to_proc ; ++k ){
	Real *rdata = &data[entity_list[k]->ProcArrayIndex()*size_per_entity];
	for(int l=0 ; l<size_per_entity ; ++l )
	  rdata[l] += *sb++;
      }
    }
  }
#endif
}



void contact_swapadd_data_array_fcs( MPI_Comm& search_comm,  
				     ContactSymComm& comm,
				     ContactCommBuffer& comm_buffer,
				     Real* data,
				     int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10008;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = size_per_entity*sizeof(Real);
  comm_buffer.Buffers(comm, bytes_per_var, &csnd_buf, &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);

  // Gather my data into send_buf
  int num_procs = comm.Num_Comm_Partners();
  for(int i=0 ; i<num_procs ; ++i ){
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy(sb, &data[entity_list[j]->fcs_index*size_per_entity],
	     bytes_per_var );
      sb += size_per_entity;
    }
  }
  
  // Post Receives
  for(int i=0 ; i<num_procs ; ++i ){
    int len = size_per_entity*comm.Num_to_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
                                                   comm.Comm_Proc_ID( i ), 
					           search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < num_procs && comm.Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_per_entity*comm.Num_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == num_procs) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<num_procs ; ++j ){
    int len = size_per_entity*comm.Num_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Comm_Proc_ID( i ), 
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == num_procs) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
#ifndef CONTACT_USE_BLOCKING_SEND
    contact_wait_msg_done(send_handles[i]);
#endif
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      Real *rdata = &data[entity_list[j]->fcs_index*size_per_entity];
      std::memcpy(sb, rdata, bytes_per_var );
      std::memset(rdata, 0, bytes_per_var );
      sb += size_per_entity;
    }
  }  

  // Assemble the data back to contact entity in ascending processor order
  bool added_mine = false;
  rb = recv_buf;
  sb = send_buf;
  for(int i=0 ; i<num_procs ; ++i ){
    contact_wait_msg_done( recv_handles[i] );
    if( !added_mine && 
	contact_processor_number(search_comm)<comm.Comm_Proc_ID( i ) ){
      // Add in my data
      added_mine = true;
      for(int j=0 ; j<num_procs ; ++j ){
	int num_to_proc = comm.Num_to_Proc( j );
	ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
	for(int k=0 ; k<num_to_proc ; ++k ){
	  Real *rdata = &data[entity_list[k]->fcs_index*size_per_entity];
	  for(int l=0 ; l<size_per_entity ; ++l )
	    rdata[l] += *sb++;
	}
      }
    }
    int num_to_proc = comm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      Real *rdata = &data[entity_list[j]->fcs_index*size_per_entity];
      for(int k=0 ; k<size_per_entity ; ++k )
	rdata[k] += *rb++;
    }
  }

  // If my data hasn't been added, add it
  if( !added_mine ){
    for(int j=0 ; j<num_procs ; ++j ){
      int num_to_proc = comm.Num_to_Proc( j );
      ContactTopologyEntity<Real>** entity_list = comm.Entity_List( j );
      for(int k=0 ; k<num_to_proc ; ++k ){
	Real *rdata = &data[entity_list[k]->fcs_index*size_per_entity];
	for(int l=0 ; l<size_per_entity ; ++l )
	  rdata[l] += *sb++;
      }
    }
  }
#endif
}

void contact_import_scratch_var( MPI_Comm& search_comm,  
				 ContactAsymComm& comm,
				 ContactCommBuffer& comm_buffer,
				 ScratchVariable &var)
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10009;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = var.Get_Size()*sizeof(Real);
  comm_buffer.Buffers_Import(comm, bytes_per_var, &csnd_buf, &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);
 
  int total_num_export = comm.Size_Export();
  int *export_index_list = comm.Export_Index_List();
 
  // Gather my data into send_buf
  for(int i=0 ; i < total_num_export ; ++i ){
    std::memcpy( sb, var.Get_Scratch(export_index_list[i]), bytes_per_var );
    sb += var.Get_Size();
  }

  // Post Receives
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    int len = var.Get_Size()*comm.Num_Import_from_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
					 comm.Import_Comm_Proc_ID( i ), 
					 search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < comm.Num_Export_Comm_Partners() && comm.Export_Comm_Proc_ID(proc_index) < my_proc) {
      sb += var.Get_Size()*comm.Num_Export_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == comm.Num_Export_Comm_Partners()) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<comm.Num_Export_Comm_Partners() ; ++j ){
    int len = var.Get_Size()*comm.Num_Export_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Export_Comm_Proc_ID( i ),
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, 
                           comm.Export_Comm_Proc_ID( i ), 
                           search_comm );
#endif
    sb += len;
    if (++i == comm.Num_Export_Comm_Partners()) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Copy the data from recv_buf into data
  rb = recv_buf;
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    //
    //  Wait until the current message is received
    //
    contact_wait_msg_done( recv_handles[i] );
    int num_from_proc = comm.Num_Import_from_Proc( i );
    int* import_index_list = comm.Import_Index_List( i );
    for(int j=0 ; j<num_from_proc ; ++j ){
      std::memcpy( var.Get_Scratch(import_index_list[j]), rb, bytes_per_var );
      rb += var.Get_Size();
    }
  }

#ifndef CONTACT_USE_BLOCKING_SEND
  for(int i = 0; i < comm.Num_Export_Comm_Partners(); ++i) {
    contact_wait_msg_done( send_handles[i] );
  }
#endif

#endif
}

void contact_import_scratch_vars( MPI_Comm& search_comm,  
			          ContactAsymComm& comm,
			          ContactCommBuffer& comm_buffer,
			          int num_vars,
			          ScratchVariable **vars)
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10010;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;

  int *var_bytes = new int[num_vars];
  int *var_size = new int[num_vars];
  int sizeof_Real = sizeof(Real);
  int size_of_all_vars = 0;
  for(int i=0 ; i<num_vars ; ++i ) {
    int cur_var_size = vars[i]->Get_Size();
    size_of_all_vars += cur_var_size;
    var_size[i] = cur_var_size;
    var_bytes[i] = cur_var_size * sizeof_Real;
  }
  int bytes_per_entity = size_of_all_vars*sizeof(Real);

  comm_buffer.Buffers_Import(comm, bytes_per_entity, &csnd_buf,
                             &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);
  
  int total_export       = comm.Size_Export();
  int *export_index_list = comm.Export_Index_List(); 


  // Gather my data into send_buf
  for(int i=0 ; i<total_export ; ++i ){
    for(int j = 0; j < num_vars; ++j) {
      std::memcpy( sb, vars[j]->Get_Scratch(export_index_list[i]), var_bytes[j] );
      sb += var_size[j];
    }
  }

  // Post Receives
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    int len = size_of_all_vars*comm.Num_Import_from_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
					 comm.Import_Comm_Proc_ID( i ), 
					 search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < comm.Num_Export_Comm_Partners() && comm.Export_Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_of_all_vars*comm.Num_Export_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == comm.Num_Export_Comm_Partners()) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<comm.Num_Export_Comm_Partners() ; ++j ){
    int len = size_of_all_vars*comm.Num_Export_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Export_Comm_Proc_ID( i ),
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, 
                           comm.Export_Comm_Proc_ID( i ), 
                           search_comm );
#endif
    sb += len;
    if (++i == comm.Num_Export_Comm_Partners()) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Copy the data from recv_buf into data
  rb = recv_buf;
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    //
    //  Wait until the current message is received
    //
    contact_wait_msg_done( recv_handles[i] );
    int num_from_proc = comm.Num_Import_from_Proc( i );
    int* import_index_list = comm.Import_Index_List( i );
    for(int j=0; j<num_from_proc; ++j ){
      for(int k=0; k < num_vars; ++k) {
        std::memcpy( vars[k]->Get_Scratch(import_index_list[j]), rb, var_bytes[k] );
        rb += var_size[k];
      }
    }
  }

#ifndef CONTACT_USE_BLOCKING_SEND
  for(int i = 0; i < comm.Num_Export_Comm_Partners(); ++i) {
    contact_wait_msg_done( send_handles[i] );
  }
#endif

  delete [] var_bytes;
  delete [] var_size;

#endif
}



void contact_import_reg_data( MPI_Comm& search_comm,  
			      ContactAsymComm& comm,
			      ContactCommBuffer& comm_buffer,
			      VariableHandle vh,
			      int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10011;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want Real* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  int bytes_per_var = size_per_entity*sizeof(Real);
  comm_buffer.Buffers_Import(comm, bytes_per_var, &csnd_buf, 
                             &crcv_buf, &recv_handles, &send_handles);
  Real *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<Real*> (csnd_buf);
  rb = recv_buf = reinterpret_cast<Real*> (crcv_buf);
  
  // Gather my data into send_buf
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    int num_to_proc = comm.Num_Export_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Export_Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      std::memcpy( sb, entity_list[j]->Variable(vh), bytes_per_var );
      sb += size_per_entity;
    }
  }

  // Post Receives
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    int len = size_per_entity*comm.Num_Import_from_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
					 comm.Import_Comm_Proc_ID( i ), 
					 search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < comm.Num_Export_Comm_Partners() && comm.Export_Comm_Proc_ID(proc_index) < my_proc) {
      sb += size_per_entity*comm.Num_Export_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == comm.Num_Export_Comm_Partners()) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<comm.Num_Export_Comm_Partners() ; ++j ){
    int len = size_per_entity*comm.Num_Export_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                comm.Export_Comm_Proc_ID( i ),
    			                        search_comm );
#else
    contact_blocking_send( mesg, sb, len, 
                           comm.Export_Comm_Proc_ID( i ), 
                           search_comm );
#endif
    sb += len;
    if (++i == comm.Num_Export_Comm_Partners()) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Copy the data from recv_buf into data
  rb = recv_buf;
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    contact_wait_msg_done( recv_handles[i] );
    int num_from_proc = comm.Num_Import_from_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Import_Entity_List( i );
    for(int j=0 ; j<num_from_proc ; ++j ){
      std::memcpy( entity_list[j]->Variable(vh), rb, bytes_per_var );
      rb += size_per_entity;
    }
  }

#ifndef CONTACT_USE_BLOCKING_SEND
  //
  //  Assure that all sends have completed before reusing the buffers
  //
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    contact_wait_msg_done( send_handles[i] );
  }
#endif

#endif
}



void contact_import_node_status( MPI_Comm& search_comm,  
			         ContactAsymComm& comm,
			         ContactCommBuffer& comm_buffer )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( search_comm ) == 1 ) return;
  int mesg = 10012;
  int my_proc = contact_processor_number(search_comm);

  // comm_buffer treates everything as char*, here we want int* buffers
  // so we do an explicit cast to that type.
  char* csnd_buf;
  char* crcv_buf;
  RequestHandle* recv_handles;
  RequestHandle* send_handles;
  comm_buffer.Buffers_Import(comm, sizeof(int), &csnd_buf, 
                             &crcv_buf, &recv_handles, &send_handles);
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = reinterpret_cast<int*>(csnd_buf);
  rb = recv_buf = reinterpret_cast<int*>(crcv_buf);
  
  // Gather my data into send_buf
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    int num_to_proc = comm.Num_Export_to_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Export_Entity_List( i );
    for(int j=0 ; j<num_to_proc ; ++j ){
      *sb++ = int(entity_list[j]->GetContext() & 
                  (ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE | ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE));
    }
  }

  // Post Receives
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    int len = comm.Num_Import_from_Proc( i );
    recv_handles[i] = contact_nonblocking_receive( mesg, rb, len, 
					 comm.Import_Comm_Proc_ID( i ), 
					 search_comm );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( search_comm );

  // Send my data
  sb = send_buf;
  int proc_index = 0;
  while (proc_index < comm.Num_Export_Comm_Partners() && comm.Export_Comm_Proc_ID(proc_index) < my_proc) {
      sb += comm.Num_Export_to_Proc( proc_index );
      proc_index++;
  }
  if (proc_index == comm.Num_Export_Comm_Partners()) {
      sb = send_buf;
      proc_index = 0;
  }
  for(int i = proc_index, j=0 ; j<comm.Num_Export_Comm_Partners() ; ++j ){
    int len = comm.Num_Export_to_Proc( i );
#ifndef CONTACT_USE_BLOCKING_SEND
    send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                         comm.Export_Comm_Proc_ID( i ),
    			                 search_comm );
#else
    contact_blocking_send( mesg, sb, len, comm.Export_Comm_Proc_ID( i ), search_comm );
#endif
    sb += len;
    if (++i == comm.Num_Export_Comm_Partners()) {
      i  = 0;
      sb = send_buf;
    }
  }

  // Copy the data from recv_buf into data
  rb = recv_buf;
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    contact_wait_msg_done( recv_handles[i] );
    int num_from_proc = comm.Num_Import_from_Proc( i );
    ContactTopologyEntity<Real>** entity_list = comm.Import_Entity_List( i );
    for(int j=0 ; j<num_from_proc ; ++j ){
      unsigned int status = (unsigned int)(*rb++);
      entity_list[j]->ClearContextBit(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE | ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
      entity_list[j]->SetContextBit(status);
    }
  }

#ifndef CONTACT_USE_BLOCKING_SEND
  //
  //  Assure that all sends have completed before reusing the buffers
  //
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    contact_wait_msg_done( send_handles[i] );
  }
#endif

#endif
}


void contact_gather( MPI_Comm& communicator, 
		     char* send_values, int num_send_values, 
		     char* recv_values, int root_rank)
{
#ifndef CONTACT_NO_MPI
  MPI_Gather(send_values, num_send_values, MPI_CHAR,
             recv_values, num_send_values, MPI_CHAR,
             root_rank, communicator );
#endif
}


void contact_gather( MPI_Comm& communicator,
		     int* send_values, int num_send_values, 
		     int* recv_values, int root_rank)
{
#ifndef CONTACT_NO_MPI
  MPI_Gather(send_values, num_send_values, MPI_INT,
             recv_values, num_send_values, MPI_INT,
             root_rank, communicator);
#endif
}

void contact_gather( MPI_Comm& communicator,
		     int send_value, int* recv_values, int root_rank )
{
#ifndef CONTACT_NO_MPI
  MPI_Gather(&send_value, 1, MPI_INT,
             recv_values, 1, MPI_INT,
             root_rank, communicator);
#endif
}


void contact_hetro_gather( MPI_Comm& communicator,
			   char* send_values, int  num_send_values,
			   char* recv_values, int* recv_sizes,
                           int root_rank )
{
#ifndef CONTACT_NO_MPI
  int comm_size = contact_number_of_processors( communicator );
  if (comm_size==1) {
    std::memcpy(recv_values,send_values,num_send_values);
  } else {
    int* offsets = new int[comm_size];
    if( contact_processor_number( communicator ) == root_rank )
    {
      offsets[0] = 0;
      for (int i = 1; i < comm_size; ++i) {
        offsets[i] = offsets[i-1]+recv_sizes[i-1];
      }
    }
    MPI_Gatherv(send_values, num_send_values, MPI_CHAR,
                recv_values, recv_sizes, offsets, MPI_CHAR,
                root_rank, communicator);
    delete [] offsets;
  }
#endif
}

void contact_hetro_gather( MPI_Comm& communicator,
			   int* send_values, int num_send_values,
			   int* recv_values, int* recv_sizes,
                           int root_rank )
{
#ifndef CONTACT_NO_MPI
  int comm_size = contact_number_of_processors( communicator );
  if (comm_size==1) {
    std::memcpy(recv_values,send_values,num_send_values*sizeof(int));
  } else {
    int* offsets = new int[comm_size];
    if( contact_processor_number( communicator ) == root_rank )
    {
      offsets[0] = 0;
      for (int i = 1; i < comm_size; ++i) {
        offsets[i] = offsets[i-1]+recv_sizes[i-1];
      }
    }
    MPI_Gatherv(send_values, num_send_values, MPI_INT,
                recv_values, recv_sizes, offsets, MPI_INT,
                root_rank, communicator);
    delete [] offsets;
  }
#endif
}


void contact_communicate_packed_buffers( MPI_Comm& communicator,
					 ContactAsymComm& comm,
					 int size_per_entity,
					 char* send_buf,
					 char* recv_buf,
					 RequestHandle* request_handles )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( communicator ) == 1 ) return;
  int mesg = 10013;
  int len;

  // Post Receives
  char* rb = recv_buf;
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i ){
    len = size_per_entity*comm.Num_Import_from_Proc( i );
    request_handles[i] = 
      contact_nonblocking_receive( mesg, rb, len, comm.Import_Comm_Proc_ID(i),
				   communicator );
    rb += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( communicator );

  // Send my data
  char* sb = send_buf;
  for(int i=0 ; i<comm.Num_Export_Comm_Partners() ; ++i ){
    len = size_per_entity*comm.Num_Export_to_Proc( i );
    contact_blocking_send( mesg, sb, len, comm.Export_Comm_Proc_ID( i ), communicator );
    sb += len;
  }

  // Wait till all of my messages have arrived
  for(int i=0 ; i<comm.Num_Import_Comm_Partners() ; ++i )
    contact_wait_msg_done( request_handles[i] );

#endif 
}

void contact_global_boundingbox( ContactBoundingBox& box, MPI_Comm& comm)
{
#ifndef CONTACT_NO_MPI
  Real l_min[3], l_max[3];
  Real g_min[3], g_max[3];
  box.get_min_max_points(l_min, l_max);
  contact_global_minimum(l_min, g_min, 3, comm);
  contact_global_maximum(l_max, g_max, 3, comm);
  box.set_box(g_min, g_max);
#endif
}



#ifndef CONTACT_NO_MPI
namespace ACME {
  //
  //  Parallel_Data_Exchange:  Exchange real data between a set of send lists and a set of 
  //  receive lists
  //
  void Parallel_Data_Exchange(Real_Vector* &send_lists,
                              Real_Vector* &recv_lists,
                              Int_Vector   &comm_partners,
                              MPI_Comm &mpi_communicator) {
    //
    //  Send the message sizes and allocate the receive lists
    //
    int num_procs = comm_partners.size();
    Int_Vector send_msg_sizes(num_procs);
    Int_Vector recv_msg_sizes(num_procs);
    MPI_Request* recv_handles = new MPI_Request[num_procs];

    for(int iproc = 0; iproc < num_procs; ++iproc) {
      send_msg_sizes[iproc] = send_lists[iproc].size();
    }    
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        MPI_Irecv(&recv_msg_sizes[iproc], 1, MPI_INT, iproc, 2900, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        MPI_Send(&send_msg_sizes[iproc], 1, MPI_INT, iproc, 2900, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
        recv_lists[iproc].resize(recv_msg_sizes[iproc]);
      }
    }
    //
    //  Send the messages
    //
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        MPI_Irecv(recv_lists[iproc].get_buffer(), recv_lists[iproc].size(), MPI_DOUBLE, iproc, 2800, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(send_lists[iproc].size() > 0) {
        MPI_Send(send_lists[iproc].get_buffer(), send_lists[iproc].size(), MPI_DOUBLE, iproc, 2800, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
      }
    }
    delete [] recv_handles;
  }

  //
  //  Parallel_Data_Exchange:  Exchange sets of object bounding boxes.
  //
  void Parallel_Data_Exchange(std::vector< std::vector<ObjectBoundingBox> > &send_lists,
                              std::vector< std::vector<ObjectBoundingBox> > &recv_lists,
                              std::vector< int > &comm_partners_send,
                              std::vector< int > &comm_partners_recv,
                              MPI_Comm &mpi_communicator ) {
    //
    //  Determine the number of processors involved in this communication
    //
    int num_procs = comm_partners_send.size();
    PRECONDITION(num_procs == send_lists.size() && num_procs == recv_lists.size() &&
                 num_procs == comm_partners_recv.size());
    //
    //  Send the message sizes and allocate the receive lists
    //
    vector<int> send_msg_sizes(num_procs);
    vector<int> recv_msg_sizes(num_procs);
    vector<MPI_Request> recv_handles(num_procs);

    for(int iproc = 0; iproc < num_procs; ++iproc) {
      send_msg_sizes[iproc] = send_lists[iproc].size();
      PRECONDITION(send_msg_sizes[iproc] == 0 || comm_partners_send[iproc] != 0);
    }    
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners_recv[iproc]) {
        MPI_Irecv(&recv_msg_sizes[iproc], 1, MPI_INT, iproc, 2900, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners_send[iproc]) {
        MPI_Send(&send_msg_sizes[iproc], 1, MPI_INT, iproc, 2900, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners_recv[iproc]) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
        recv_lists[iproc].resize(recv_msg_sizes[iproc]);
      }
    }
    //
    //  Send the actual messages as raw byte streams.
    //
    int class_size = sizeof(ObjectBoundingBox);

    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        char* recv_buffer = (char*)(&(recv_lists[iproc][0]));
        MPI_Irecv(recv_buffer, recv_lists[iproc].size() * class_size, MPI_CHAR, 
                  iproc, 2800, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(send_lists[iproc].size() > 0) {
        char* send_buffer = (char*)(&(send_lists[iproc][0]));

        MPI_Send(send_buffer, send_lists[iproc].size() * class_size, MPI_CHAR, 
                 iproc, 2800, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
      }
    }
  }



  //
  //  Parallel_Data_Exchange:  Exchange boolean lists
  //
  void Parallel_Data_Exchange(std::vector< vector<char> > &send_lists,
                              std::vector< vector<char> > &recv_lists,
                              std::vector< int > &comm_partners,
                              MPI_Comm &mpi_communicator ) {
    //
    //  Determine the number of processors involved in this communication
    //
    int num_procs = comm_partners.size();
    PRECONDITION(num_procs == send_lists.size() && num_procs == recv_lists.size());
    //
    //  Send the message sizes and allocate the receive lists
    //
    vector<int> send_msg_sizes(num_procs);
    vector<int> recv_msg_sizes(num_procs);
    vector<MPI_Request> recv_handles(num_procs);

    for(int iproc = 0; iproc < num_procs; ++iproc) {
      send_msg_sizes[iproc] = send_lists[iproc].size();
      PRECONDITION(send_msg_sizes[iproc] == 0 || comm_partners[iproc] != 0);
    }    
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        MPI_Irecv(&recv_msg_sizes[iproc], 1, MPI_INT, iproc, 2900, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        MPI_Send(&send_msg_sizes[iproc], 1, MPI_INT, iproc, 2900, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
        recv_lists[iproc].resize(recv_msg_sizes[iproc]);
      }
    }
    //
    //  Send the actual messages as raw byte streams.
    //
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        char* recv_buffer = (char*)&recv_lists[iproc][0];
        int recv_size = recv_lists[iproc].size();
        MPI_Irecv(recv_buffer, recv_size, MPI_CHAR, 
                  iproc, 2800, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(send_lists[iproc].size() > 0) {
        char* send_buffer = (char*)&send_lists[iproc][0];
        int send_size = send_lists[iproc].size();
        MPI_Send(send_buffer, send_size, MPI_CHAR, 
                 iproc, 2800, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
      }
    }
  }



  //
  //  Parallel_Data_Exchange:  Exchange integer lists, determine comm plan on the fly
  //
  void Parallel_Data_Exchange(std::vector< vector<int> > &send_lists,
                              std::vector< vector<int> > &recv_lists,
                              MPI_Comm &mpi_communicator ) {
    //
    //  Determine the number of processors involved in this communication
    //
    const int msg_tag = 10232;
    int num_procs;
    MPI_Comm_size(mpi_communicator, &num_procs);
    int my_proc;
    MPI_Comm_rank(mpi_communicator, &my_proc);
    PRECONDITION(num_procs == send_lists.size() && num_procs == recv_lists.size());
    //
    //  Determine the number of messages each processor will receive
    //
    vector<int> local_number_to_receive(num_procs, 0);
    vector<int> global_number_to_receive(num_procs, 0);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(send_lists[iproc].size() > 0) local_number_to_receive[iproc] = 1;
    }
    MPI_Reduce(&local_number_to_receive[0], 
               &global_number_to_receive[0], 
               num_procs, 
               MPI_INT,
               MPI_SUM, 
               0,
               mpi_communicator);
    int num_to_recv;
    MPI_Scatter(&global_number_to_receive[0], 1, MPI_INT, &num_to_recv, 1, MPI_INT, 0, mpi_communicator);
    //
    //  Now each processor knows how many messages it will recive, but does not know the message lenghts or where
    //  the messages will be recived from.  Need to extract this information.
    //
    vector<MPI_Request> recv_handles(num_procs);
    //
    //  Post a recieve for each expected message.
    //
    vector<int> recv_size_buffers(num_to_recv);
    for(int imsg = 0; imsg < num_to_recv; ++imsg) {
      int *recv_buffer = &(recv_size_buffers[imsg]);
      MPI_Irecv(recv_buffer, 1, MPI_INT, MPI_ANY_SOURCE,
                msg_tag, mpi_communicator, &recv_handles[imsg]);
    }
    MPI_Barrier(mpi_communicator);
    //
    //  Send message lengths
    //
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(send_lists[iproc].size() > 0) {
        int send_length = send_lists[iproc].size();
        MPI_Send(&send_length, 1, MPI_INT, iproc, msg_tag, mpi_communicator);
      }
    }
    //
    //  Get each message and place the length in the proper place in the length array
    //
    for(int imsg = 0; imsg < num_to_recv; ++imsg) {
      MPI_Status status;
      MPI_Wait(&recv_handles[imsg], &status);
      recv_lists[status.MPI_SOURCE].resize(recv_size_buffers[imsg]);
    }
    //
    //  Send the actual messages as raw byte streams.
    //
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        char* recv_buffer = (char*)&recv_lists[iproc][0];
        int recv_size = recv_lists[iproc].size();
        MPI_Irecv(recv_buffer, recv_size, MPI_INT, 
                  iproc, 2800, mpi_communicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(send_lists[iproc].size() > 0) {
        char* send_buffer = (char*)&send_lists[iproc][0];
        int send_size = send_lists[iproc].size();
        MPI_Send(send_buffer, send_size, MPI_INT, 
                 iproc, 2800, mpi_communicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recv_lists[iproc].size() > 0) {
        MPI_Status status;
        MPI_Wait( &recv_handles[iproc], &status );
      }
    }
  }
} // end namespace ACME
#endif
