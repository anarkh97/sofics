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


#ifndef contact_timer_h_
#define contact_timer_h_

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif

class CString;

class ContactTimer {

 public:

  ContactTimer(
#ifndef CONTACT_NO_MPI
	       MPI_Comm& communicator
#else
	       int communicator
#endif
	       );

  ~ContactTimer();

  int Register_Timer( const CString );

  void Start_Timer( int timer_handle );
  void Stop_Timer( int timer_handle );

  void Output_Processor_Timings();
  void Output_Timings();

 private:

  // Don't allow copying 
  ContactTimer( const ContactTimer& );
  ContactTimer operator=(const ContactTimer& );

#ifndef CONTACT_NO_MPI
  MPI_Comm communicator;
#else
  int communicator;
#endif

  int allocated_size;
  int number_timers;
  double* times;
  double* start_times;
  CString** names;
  int* timing_calls;

  void Increase_Size();
  
};

#endif
