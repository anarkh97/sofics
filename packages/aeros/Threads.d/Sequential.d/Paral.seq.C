#include <cstdlib>
#include <sys/types.h>
#include <unistd.h>
#include <malloc.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <Threads.d/Paral.h>
#include <Timers.d/GetTime.h>
#include <Timers.d/DistTimer.h>

// isParal = 0 sequential
// isParal = 1 parallel

// NOTE: if you change this variable in a 
//       sub-directory, you can force 
//       sequential execution for that part of the code.

int isParal = 1;

// int x = 0;
int oneMb    = 1024*1024;
int initSize = 2*oneMb;
int numMb    = 32;
int growth   = numMb*oneMb;

pid_t pproc = 0;
int numProc = 0;
pid_t *childProc = NULL;
void **arenas;
int *maxSizes;
long *curSizes;
long currentSizes = 0;
void *pparena = 0;
int zeroFd = -1;

#ifdef TFLOP
extern "C" int heap_info(int*, int*, int*, int*);
#endif

void * operator new(size_t size)
{
#ifdef TFLOP
  int fragments, total_free, largest_free, total_used;
  heap_info(&fragments, &total_free, &largest_free, &total_used);
#endif
  void *newMem = malloc(size);
  if(newMem == 0 && size != 0) {
    fprintf(stderr,"MEMORY ERROR: Trying to allocate %d bytes\n",size);
    fprintf(stderr,"            : Memory allocated by new so far %d bytes\n",
                   currentSizes);
#ifdef TFLOP
    fprintf(stderr,"Before trying to allocated:\n");
    fprintf(stderr,"memory usage: total free %d, largest free %d, total_used %d\n",
                total_free, largest_free, total_used);
#endif
  }
  currentSizes += size;
  return newMem;
}

void operator delete(void *p)
{
 free(p);
}



long
ThreadManager::getLocalMem()
{
 return currentSizes;
}

void
ThreadLock::lock()
{
// ussetlock(lockV);
}

void
ThreadLock::unlock()
{
// usunsetlock(lockV);
}


ThreadLock::ThreadLock()
{
// lockV = usnewlock(usPtr);
}

typedef void (*P)(void *, size_t);
ThreadManager::ThreadManager(int nThr)
{
 nThr = 1;
 numThreads = nThr;
 
 firstProc = 0;
}

ThreadManager::~ThreadManager()
{
}

long
ThreadManager::memoryUsed()
{
 return currentSizes;
}

void
ThreadManager::memUsage()
{
}

void
ThreadManager::execTasks(int ntasks, TaskDescr **td)
{
 int i;
 for(i = 0; i < ntasks; ++i) 
   td[i]->run();
}

void
ThreadManager::execTasks(int ntasks, TaskDescr *td)
{
 int i;

 for(i = 0; i < ntasks; ++i) {
   td->runFor(i);
 }
}

void
ThreadManager::execParal(int ntasks, TaskDescr **td)
{
 timer = NULL;
 int i;
 for(i=0; i < ntasks; ++i)
    td[i]->run();
}

void
ThreadManager::execParal(int ntasks, TaskDescr *td)
{
 timer = NULL;
 int i;
 for(i = 0; i < ntasks; ++i) {
   td->runFor(i);
 }
}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr **td)
{
 timer = &thisTimer;
 double initTime = getTime();
 long initMem  = threadManager->getLocalMem();
 int i;
 for(i=0; i < ntasks; ++i)
    td[i]->run();

 long finalMem = threadManager->getLocalMem();
 timer->addTo(0, finalMem-initMem, getTime()-initTime);

}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr *td)
{
 timer = &thisTimer;
 double initTime = getTime();
 long initMem  = threadManager->getLocalMem();

 int i;
 for(i = 0; i < ntasks; ++i)
   td->runFor(i);

 long finalMem = threadManager->getLocalMem();
 timer->addTo(0, finalMem-initMem, getTime()-initTime);
}

