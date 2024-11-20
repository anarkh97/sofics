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

 
#include "ContactTimer.h"
#include "CString.h"
#include <iostream>
#include <cstring>

ContactTimer::ContactTimer(
#ifndef CONTACT_NO_MPI
			   MPI_Comm& Comm
#else
			   int Comm
#endif
			   )
  : communicator( Comm )
{
  allocated_size = 0;
  number_timers = 0;
  times = NULL;
  names = NULL;
  start_times = NULL;
  timing_calls = NULL;
}

ContactTimer::~ContactTimer()
{
  delete [] times;
  for( int i=0 ; i<number_timers ; ++i) delete names[i];
  delete [] names;
  delete [] start_times;
  delete [] timing_calls;
}

int ContactTimer::Register_Timer( const CString Name )
{
#ifdef CONTACT_TIMINGS
  int i;
  if( number_timers == allocated_size ) Increase_Size();
  
  // Truncate the name to 34 characters
  char name[35] = {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};

  int len = Name.size();
  len = len > 34 ? 34 : len;
  for( i=0 ; i<len ; ++i) name[i] = Name[i];
  for( ; i<34 ; ++i) name[i] = ' ';
  times[number_timers] = 0.0;
  names[number_timers] = new CString( &(name[0]) );
  timing_calls[number_timers] = 0;
  number_timers++;
  return number_timers-1;
#else
  return 0;
#endif
}


void ContactTimer::Increase_Size()
{
#ifdef CONTACT_TIMINGS
  int old_allocated_size = allocated_size;
  allocated_size += 10;

  double* old_times = times;
  times = new double[allocated_size];
  std::memcpy( times, old_times, old_allocated_size*sizeof(double) );
  delete [] old_times;

  CString** old_names = names;
  names = new CString*[allocated_size];
  std::memcpy( names, old_names, old_allocated_size*sizeof(CString*) );
  delete [] old_names;

  double* old_start_time = start_times;
  start_times = new double[allocated_size];
  std::memcpy( start_times, old_start_time, old_allocated_size*sizeof(double) );
  delete [] old_start_time;

  int* old_timing_calls = timing_calls;
  timing_calls = new int[allocated_size];
  std::memcpy( timing_calls, old_timing_calls, old_allocated_size*sizeof(int) );
  delete [] old_timing_calls;
#endif
}

void ContactTimer::Start_Timer( int timer_handle )
{
#ifdef CONTACT_TIMINGS
  timing_calls[timer_handle] += 1;
#ifndef CONTACT_NO_MPI
  start_times[timer_handle] = MPI_Wtime();
#else
  start_times[timer_handle] = 0.0;
#endif
#endif
}

void ContactTimer::Stop_Timer( int timer_handle )
{
#ifdef CONTACT_TIMINGS
#ifndef CONTACT_NO_MPI
  times[timer_handle] += (MPI_Wtime() - start_times[timer_handle]);
#else
  times[timer_handle] = 0.0;
#endif
#endif
}

void ContactTimer::Output_Timings()
{
#ifdef CONTACT_TIMINGS
  if( number_timers == 0 ) return;

  int my_proc = 0;
#ifndef CONTACT_NO_MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &my_proc );
#endif
  int num_procs = 1;
#ifndef CONTACT_NO_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
#endif
  std::cout.setf(std::ios::showpoint);
  std::cout.setf(std::ios::fixed);
  std::cout.precision(4);
  if( num_procs > 1 ){
    if( my_proc == 0 ){
      std::cout << "\n\n";
      std::cout << "Timer #   # Calls    Name                                     Min          Max           Ave\n";
      std::cout << "-------  ---------  ------------------------------------   ----------   ----------   ----------\n";
    }
  } else {
      std::cout << "\n\n";
      std::cout << "Timer #   # Calls    Name                                     Time   \n";
      std::cout << "-------  ---------  ------------------------------------   ----------\n";
  }
#ifndef CONTACT_NO_MPI
  int n = number_timers;
#endif
  double* min_val=NULL;
  double* max_val=NULL;
  double* ave_val=NULL;
  if( num_procs > 1 ){
    min_val = new double[number_timers];
    max_val = new double[number_timers];
    ave_val = new double[number_timers];
#ifndef CONTACT_NO_MPI
    MPI_Allreduce( times, min_val, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce( times, max_val, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce( times, ave_val, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  }
  if( my_proc == 0 ){
    for( int i=0 ; i<number_timers ; ++i){
      if( num_procs > 1 ) ave_val[i] /= num_procs;
      std::cout.width(4);
      std::cout << i;
      std::cout.width(12);
      std::cout << timing_calls[i];
      std::cout << "     " << *(names[i]);
      std::cout << "    ";
      if( num_procs > 1 ){
	std::cout.setf(std::ios::scientific,std::ios::floatfield);
	std::cout << min_val[i];
	std::cout << "   ";
	std::cout.setf(std::ios::scientific,std::ios::floatfield);
	std::cout << max_val[i]; 
	std::cout << "   ";
	std::cout.setf(std::ios::scientific,std::ios::floatfield);
	std::cout << ave_val[i];
	std::cout << "\n";
      } else {
	std::cout.setf(std::ios::scientific,std::ios::floatfield);
	std::cout << times[i] << "\n";
      }
    }
    std::cout.flush();
  }
  if( num_procs > 1 ){
    delete [] min_val;
    delete [] max_val;
    delete [] ave_val;
  }
#endif
}
