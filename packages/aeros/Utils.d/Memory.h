#ifndef _MEMORY_H_
#define _MEMORY_H_

// function to return memory used
#include <Threads.d/Paral.h>
#include <Utils.d/Memory.h>

extern long currentSizes;

#if defined(sgi) && ! defined(_OPENMP)
extern pid_t * childProc;
#endif


inline
long
memoryUsed()
{
#if defined(sgi) && ! defined(_OPENMP)
 if(childProc == 0)
   return currentSizes;
 else
   return threadManager->memoryUsed();
#else
 return currentSizes;
#endif
}

#ifdef TFLOP
extern "C" int heap_info(int*, int*, int*, int*);
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#include <Sandia.d/SInterface.h>
extern Feti98Params* parameters;
#endif

inline
double *
checkAndAllocateMemory(int size, const char *message)
{
  double *pointer=0;
#ifdef TFLOP
  int fragments, total_free, largest_free, total_used;
  heap_info(&fragments, &total_free, &largest_free, &total_used);

//  for debugging
 if((parameters->verbose_flag & (1 << 3)) != 0)
   fprintf(stdout,"Processor %d is about to allocate %d bytes\n",
           scom->myID(), sizeof(double)*size);

  if( sizeof(double)*size < largest_free)
    pointer = new double[size];
  else {
    fprintf(stdout,"%s",message);
    fprintf(stdout,"Requested %d bytes when %d bytes were available\n",
                   sizeof(double)*size, largest_free);
  }
#else
  pointer = new double[size];
#endif
  return pointer;
}

#endif
