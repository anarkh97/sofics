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


// This file represents ALL of the communication calls that the contact library
// makes.  All communication should occur through one of these wrapper 
// functions.
//
// Currently we only support MPI with a host code supplied communicator.

#ifndef Contact_Communication_h_
#define Contact_Communication_h_

class ContactAsymComm;
class ContactSymComm;
class ContactCommBuffer;
class ContactBoundingBox;
class ScratchVariable;

#include "Contact_Defines.h"
#include "ContactVector.h"
#include <vector>
#include "ContactBoundingBox.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
typedef MPI_Request RequestHandle;
#else
typedef int RequestHandle;
#endif

int contact_processor_number( MPI_Comm& communicator );
int contact_any_processor_number();
int contact_any_message_number();
int contact_number_of_processors( MPI_Comm& communicator );

void contact_global_sync( MPI_Comm& );
void contact_print_sync_start( MPI_Comm& );
void contact_print_sync_end( MPI_Comm& );

// Global sum of a single value
int  contact_global_sum        ( int,  MPI_Comm& );
int  contact_global_error_check( int,  MPI_Comm& );
Real contact_global_sum        ( Real, MPI_Comm& );

// Global sum of an array of values from ivalues (which is untouched) and
// the result is put in ovalues
void contact_global_sum(     int*  ivalues, int*  ovalues, int nvalues,
			     MPI_Comm& );
void contact_global_sum(     Real* ivalues, Real* ovalues, int nvalues,
			     MPI_Comm& );

// Global min/max of a single value of array of values
int  contact_global_minimum( int,  MPI_Comm& );
Real contact_global_minimum( Real, MPI_Comm& );

int  contact_global_maximum( int,  MPI_Comm& );
Real contact_global_maximum( Real, MPI_Comm& );

bool contact_global_or( bool, MPI_Comm& communicator );

int contact_broadcast( int * array, int size, int sender, 
		       MPI_Comm& communicator );
int contact_broadcast( Real * array, int size, int sender, 
		       MPI_Comm& communicator );
int contact_broadcast( char * array, int size, int sender, 
		       MPI_Comm& communicator );

// Global min/max of an array of values from ivalues (which is untouched) and
// the result is put in ovalues
void contact_global_minimum(int*  ivalues,int*  ovalues,int nvalues,MPI_Comm&);
void contact_global_minimum(Real* ivalues,Real* ovalues,int nvalues,MPI_Comm&);
void contact_global_maximum(int*  ivalues,int*  ovalues,int nvalues,MPI_Comm&);
void contact_global_maximum(Real* ivalues,Real* ovalues,int nvalues,MPI_Comm&);

void contact_blocking_send( int mesg, int* data, int len, int proc, 
			    MPI_Comm& communicator );
void contact_blocking_send( int mesg, Real* data, int len, int proc, 
			    MPI_Comm& communicator );
void contact_blocking_send( int mesg, char* data, int len, int proc, 
			    MPI_Comm& communicator );

RequestHandle contact_nonblocking_send(int mesg, Real* data, int len, int proc,
                                       MPI_Comm& communicator);
RequestHandle contact_nonblocking_send(int mesg,int* data, int len, int proc,
                                       MPI_Comm& communicator);
RequestHandle contact_nonblocking_send(int mesg,char* data, int len, int proc,
                                       MPI_Comm& communicator);


RequestHandle contact_nonblocking_receive( int mesg, int* data, int len,
					   int proc, MPI_Comm& communicator );
RequestHandle contact_nonblocking_receive( int mesg, Real* data, int len,
					   int proc, MPI_Comm& communicator );
RequestHandle contact_nonblocking_receive( int mesg, char* data, int len,
					   int proc, MPI_Comm& communicator );

int contact_wait_msg_done( RequestHandle request );
int contact_test_msg_done( RequestHandle request );

void contact_swapadd_reg_data( MPI_Comm& communicator,
			       ContactSymComm& comm,
			       ContactCommBuffer& comm_buffer,
			       VariableHandle var,
			       int size_of_var );

void contact_swapadd_temp_tag( MPI_Comm& communicator,
			       ContactSymComm& comm,
			       ContactCommBuffer& comm_buffer );

void contact_swapadd_context( MPI_Comm& communicator,
			       ContactSymComm& comm,
			       ContactCommBuffer& comm_buffer );

void contact_reduce_add_scratch_var( MPI_Comm& communicator,
				     ContactAsymComm& comm,
			             ContactCommBuffer& comm_buffer,
			             ScratchVariable &var);

void contact_reduce_add_scratch_vars( MPI_Comm& communicator,
				      ContactAsymComm& comm,
				      ContactCommBuffer& comm_buffer,
				      int num_vars,
				      ScratchVariable** vars);

void contact_swap_edge_faces( MPI_Comm& communicator,
		              ContactSymComm& comm,
		              ContactCommBuffer& comm_buffer );

void contact_swapadd_scratch_var( MPI_Comm& communicator,
				  ContactSymComm& comm,
			          ContactCommBuffer& comm_buffer,
			          ScratchVarHandle var,
			          int size_of_var );

void contact_swapadd_scratch_vars( MPI_Comm& communicator,
		      	           ContactSymComm& comm,
				   ContactCommBuffer& comm_buffer,
			           int num_vars,
			           ScratchVarHandle* vars,
			           int* sizes_of_var );

void contact_swapadd_data_array_fcs( MPI_Comm& communicator,
				     ContactSymComm& comm,
				     ContactCommBuffer& comm_buffer,
				     Real* data,
				     int size_per_entity );

void contact_swapadd_data_array_fcs( MPI_Comm& communicator,
				     ContactSymComm& comm,
				     ContactCommBuffer& comm_buffer,
				     int* data,
				     int size_per_entity );

void contact_swapadd_data_array( MPI_Comm& communicator,
				 ContactSymComm& comm,
				 ContactCommBuffer& comm_buffer,
				 Real* data,
				 int size_per_entity );

void contact_swapadd_data_array( MPI_Comm& communicator,
				 ContactSymComm& comm,
				 ContactCommBuffer& comm_buffer,
				 int* data,
				 int size_per_entity );

void contact_import_scratch_var( MPI_Comm& communicator,
			         ContactAsymComm& comm,
			         ContactCommBuffer& comm_buffer,
			         ScratchVariable &vh);

void contact_import_scratch_vars( MPI_Comm& communicator,
			          ContactAsymComm& comm,
			          ContactCommBuffer& comm_buffer,
                                  int num_vars,
			          ScratchVariable** vars);

void contact_import_reg_data( MPI_Comm& communicator,
			      ContactAsymComm& comm,
			      ContactCommBuffer& comm_buffer,
			      VariableHandle var,
			      int size_of_var );
                              
void contact_import_node_status( MPI_Comm& communicator,
			         ContactAsymComm& comm,
			         ContactCommBuffer& comm_buffer );

void contact_gather( MPI_Comm& communicator, 
		     char* send_values, int num_send_values, 
		     char* recv_values, int root_rank);

void contact_gather( MPI_Comm& communicator,
		     int* send_values, int num_send_values, 
		     int* recv_values, int root_rank);

void contact_gather( MPI_Comm& communicator,
		     int send_value, int* recv_values, int root_rank );

void contact_hetro_gather( MPI_Comm& communicator,
			   char* send_values, int num_send_values,
			   char* recv_values, int* recv_sizes,
                           int root_rank );

void contact_communicate_packed_buffers( MPI_Comm& communicator,
					 ContactAsymComm& comm,
					 int size_per_entity,
					 char* send_buf,
					 char* recv_buf,
					 RequestHandle* request_handles );
                                         
void contact_global_boundingbox( ContactBoundingBox&,  MPI_Comm& );

namespace ACME {
#ifndef CONTACT_NO_MPI
  //
  //  Parallel_Data_Exchange:  Exchange real data between a set of send lists and a set of 
  //  receive lists
  //
  void Parallel_Data_Exchange(Real_Vector* &send_lists,
                              Real_Vector* &recv_lists,
                              Int_Vector   &comm_partners,
                              MPI_Comm &mpi_communicator);
  //
  //  Parallel_Data_Exchange:  Exchange sets of object bounding boxes.
  //
  void Parallel_Data_Exchange(std::vector< std::vector<ObjectBoundingBox> > &boxes_to_send,
                              std::vector< std::vector<ObjectBoundingBox> > &boxes_to_recv,
                              std::vector< int > &comm_partners_send,
                              std::vector< int > &comm_partners_recv,
                              MPI_Comm &mpi_communicator );

  //
  //  Parallel_Data_Exchange:  Exchange boolean lists
  //
  void Parallel_Data_Exchange(std::vector< std::vector<char> > &send_list,
                              std::vector< std::vector<char> > &recv_list,
                              std::vector< int > &comm_partners,
                              MPI_Comm &mpi_communicator );

  //
  //  Parallel_Data_Exchange:  Exchange integer lists with an unknown comm plan
  //
  void Parallel_Data_Exchange(std::vector< std::vector<int> > &send_list,
                              std::vector< std::vector<int> > &recv_list,
                              MPI_Comm &mpi_communicator );

#endif
}

#endif
