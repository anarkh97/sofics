#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <sys/types.h>


#include <unistd.h>

//--- UH ---
#ifndef NO_MALLOC_DOT_H
#include <malloc.h>
#endif
//--- UH ---

#ifndef WINDOWS
#include <sys/mman.h>
#endif
#include <fcntl.h>

#include <Threads.d/Paral.h>
#include <Utils.d/DistHelper.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>

long currentSizes = 0;


bool loud = false;

void * operator new(size_t size) //deprecated in c++11: throw(std::bad_alloc)
{
  //currentSizes += size;
  currentSizes = (long(size)/long(1024)) + currentSizes; // PJSA convert to KB since currentSizes is long (max 2147483647 = 2 GB, not enough)
  return malloc(size);
}

long
ThreadManager::getLocalMem()
{
 return currentSizes;
}

void operator delete(void *p) throw()
{
 free(p);
}

typedef void (*P)(void *, size_t);
ThreadManager::ThreadManager(int nThr)
{
#if defined(_OPENMP)
 numThreads = nThr;
// fprintf(stderr, "Forcing %d threads\n",numThreads);
 omp_set_dynamic(0);
 omp_set_num_threads(numThreads);
#else
 nThr = 1;
 numThreads = nThr;
#endif
 //fprintf(stderr, " In ThreadManager::ThreadManager(...), Creating %d threads\n", numThreads);
}

long
ThreadManager::memoryUsed()
{
 return currentSizes;
}

void
ThreadManager::execParal(int ntasks, TaskDescr **td)
{
  timer = 0;
  int i;
#ifdef _OPENMP
  std::runtime_error *e = 0;
  #pragma omp parallel for schedule(static,1) shared(e)
  for(i = 0; i < ntasks; ++i) {
    try {
      td[i]->run();
    }
    catch(std::runtime_error &_e) {
      #pragma omp critical
      e = new std::runtime_error(_e);
    }
  }
  if(e) throw(*e);
#else
  for(i=0; i < ntasks; ++i)
    td[i]->run();
#endif
}

void
ThreadManager::execParal(int ntasks, TaskDescr *td)
{
  timer = 0;
  int i;
#ifdef _OPENMP
  std::runtime_error *e = 0;
  #pragma omp parallel for schedule(static,1) shared(e)
  for(i = 0; i < ntasks; ++i) {
    try {
      td->runFor(i);
    }
    catch(std::runtime_error &_e) {
      #pragma omp critical
      e = new std::runtime_error(_e);
    }
  }
  if(e) throw(*e);
#else
  for(i = 0; i < ntasks; ++i)
    td->runFor(i);
#endif
}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr **td)
{
  timer = &thisTimer;
  double initTime = getTime();
  long initMem = threadManager->getLocalMem();
  int i;
#ifdef _OPENMP
  std::runtime_error *e = 0;
  #pragma omp parallel for schedule(static,1) shared(e)
  for(i = 0; i < ntasks; ++i) {
    try {
      td[i]->run();
    }
    catch(std::runtime_error &_e) {
      #pragma omp critical
      e = new std::runtime_error(_e);
    }
  }
  if(e) throw(*e);
#else
  for(i = 0; i < ntasks; ++i)
    td[i]->run();
#endif
  long finalMem = threadManager->getLocalMem();
  timer->addTo(0, finalMem-initMem, getTime()-initTime);
}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr *td)
{
  timer = &thisTimer;

  double initTime = getTime();
  long initMem = threadManager->getLocalMem();

  int i;
#ifdef _OPENMP
  std::runtime_error *e = 0;
  #pragma omp parallel for schedule(static,1) shared(e)
  for(i = 0; i < ntasks; ++i) {
    try {
      td->runFor(i);
    }
    catch(std::runtime_error &_e) {
      #pragma omp critical
      e = new std::runtime_error(_e);
    }
  }
  if(e) throw(*e);
#else
  for(i = 0; i < ntasks; ++i)
    td->runFor(i);
#endif

  long finalMem = threadManager->getLocalMem();
  timer->addTo(0, finalMem-initMem, getTime()-initTime);
}

#if defined(sgi) &&  !defined(_OPENMP)
barrier_t *
ThreadManager::getBarrier()
{
 return new_barrier(usPtr);
}
#endif

