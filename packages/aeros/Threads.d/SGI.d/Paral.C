#include <cstdlib>
#include <sys/types.h>
#include <sys/prctl.h>
#include <ulocks.h>
#include <unistd.h>
#include <malloc.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <Threads.d/Paral.h>
#include <Timers.d/Timing.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>

usptr_t * usPtr;

// isParal = 0 sequential
// isParal = 1 parallel

// NOTE: if you change this variable in a 
//       sub-directory, you can force 
//       sequential execution for that part of the code.

int isParal = 1;

int x = 0;
int initSize=2*1024*1024;
//int growth = 16*1024*1024;
int growth = 32*1024*1024; 
//int growth =  64*1024*1024;

ulock_t allocLock = NULL;

pid_t pproc = 0;
int numProc = 0;
pid_t *childProc = NULL;
void **arenas;
int *maxSizes;
long *curSizes;
long currentSizes = 0;
void *pparena = 0;
int zeroFd = -1;


void *
arenaGrow(size_t size, void *)
{
 // sbrk must be sequentialized
 /*double t0 = -getTime();
 ussetlock(allocLock);
 t0 += getTime();
 double t1 = -getTime();
 void *p = sbrk(size);
 t1 += getTime();
 usunsetlock(allocLock);
 fprintf(stderr,"Growth timings: %f %f\n",t0,t1);
 return p;
 return malloc(size);
*/

 void *gr =  mmap(0, size, PROT_READ|PROT_WRITE , MAP_PRIVATE, zeroFd, 0);
 if(gr == 0) { fprintf(stderr,"Out of memory\n"); exit(-1); }

// Test on touching initial memory allocation
/*
 char *mem = (char *)gr;
 int i;
 for(i=0; i < size; i += 4096)
   mem[i] =  0;
*/
 return gr;
}

void * operator new(size_t size)
{

// fprintf(stderr,"allocated %20d bytes  Running Total %14.3f Mb \n",
//                 size, memoryUsed()/(1024.0*1024.0));

 // get my process ID
 pid_t myPid = getpid();

 if(pproc == 0 || myPid == pproc) {
   if(pparena == 0) {
      pparena = malloc(initSize);
      zeroFd = open("/dev/zero",O_RDWR);
      acreate(pparena,  initSize, 0, NULL, &arenaGrow);
      amallopt(M_BLKSZ, growth, pparena);
   }
   currentSizes += size;
   return amalloc(size,pparena);
 }

 // Find my process number
 int iProc; 
 for(iProc = 0; iProc < numProc && myPid != childProc[iProc]; ++iProc)
  ;

 if(iProc == numProc) {
   fprintf(stderr,"Error\n");
   return malloc(size);
 }

 curSizes[iProc] += size;
 void *p = amalloc(size,arenas[iProc]); 
 return p;
}

long
ThreadManager::getLocalMem()
{
 // get my process ID
 pid_t myPid = getpid();
 if(pproc == 0 || myPid == pproc)
   return currentSizes;

 // Find my process number
 int iProc;
 for(iProc = 0; iProc < numProc && myPid != childProc[iProc]; ++iProc)
  ;
 if(iProc == numProc) {
   fprintf(stderr,"Error\n");
   return 0;
 }
 return curSizes[iProc];
}

void operator delete(void *p)
{
// fprintf(stderr,"deleting %x\n",p);
 pid_t myPid = getpid();
 if(myPid == pproc || pproc==0) {
    afree(p,pparena);
    return;
 }
 myPid = getpid();
 // Find my process number
 int iProc; 
 for(iProc = 0; iProc < numProc && myPid != childProc[iProc]; ++iProc)
  ;
 if(iProc == numProc) {
   fprintf(stderr,"Error\n");
   free(p);
   return;
 }
 // Free using my arena...
 afree(p,arenas[iProc]);
 return; 
}




void
ThreadLock::lock()
{
 ussetlock(lockV);
}

void
ThreadLock::unlock()
{
 usunsetlock(lockV);
}


ThreadLock::ThreadLock()
{
 lockV = usnewlock(usPtr);
}

typedef void (*P)(void *, size_t);
ThreadManager::ThreadManager(int nThr)
{
 numThreads = nThr;
 
 usconfig(CONF_INITUSERS, numThreads);
 usPtr = usinit("/dev/zero");

 firstProc = 0;
 allDone = usnewsema(usPtr , 0);
 readyProc =  usnewsema(usPtr , 0);
 sprocListLock = usnewlock(usPtr);
 allocLock = usnewlock(usPtr);

 allProc = new OneSproc[numThreads-1];

 // Setup arenas
 numProc = numThreads-1;
 arenas = new void *[numThreads-1];
 curSizes = new long[numThreads-1];
 maxSizes = new int[numThreads-1];
 childProc = new pid_t[numThreads-1];
 int i;
 for(i = 0; i < numThreads-1; ++i) {
    OneSproc *thisProc = allProc + i;
    thisProc->done = allDone;
    thisProc->wait = usnewsema(usPtr , 0);
    if(thisProc->wait ==0) perror("Could not get a semaphore:");
    thisProc->tman = this;
    thisProc->myNum = i;
    thisProc->step = numThreads;
//    pid_t pid = sproc(OneSproc::run, PR_SALL, thisProc);
    pid_t pid = sprocsp((P)OneSproc::run, PR_SALL, thisProc,NULL,10000000);
    if(pid < 0) perror("Thread did not start:");

    childProc[i] = pid;
    curSizes[i] = maxSizes[i] = 0;
 }
 for(i = 0; i < numThreads-1; ++i) {
//  set up a memory arena;
    arenas[i] = malloc(initSize);
    acreate(arenas[i], initSize, 0, NULL, &arenaGrow);
    amallopt(M_BLKSZ, growth, arenas[i]);
    curSizes[i] = maxSizes[i] = 0;
 }
 pproc = getpid();

 fprintf(stderr," ... Setting Memory Arenas          ...\n");
}

ThreadManager::~ThreadManager()
{
 int i;
 for(i = 0; i < numThreads-1; ++i)
   allProc[i].allTasks = 0;
 for(i = 0; i < numThreads-1; ++i)
   usvsema(allProc[i].wait);
 long totSizes = 0;
 for(i = 0; i < numThreads-1; ++i)
  totSizes += curSizes[i];
 // fprintf(stderr,
 //  "Total memory consumption in megabytes = %10.3f\n",totSizes/(1024.*1024.));
}

long
ThreadManager::memoryUsed()
{
 int i;
 long totSizes = currentSizes;
 for(i = 0; i < numThreads-1; ++i)
    totSizes += curSizes[i];

 return totSizes;
}

long
memoryUsed()
{
 if(childProc == NULL)
   return currentSizes;
 else
   return threadManager->memoryUsed();
}

void
ThreadManager::memUsage()
{
 long totSizes = 0;

 long minMem = currentSizes;
 long avgMem = 0;
 long maxMem = currentSizes;

 int i;
 for(i = 0; i < numThreads-1; ++i) {
    minMem = min(minMem, curSizes[i]);
    maxMem = max(maxMem, curSizes[i]);
    totSizes += curSizes[i];
 }

 minMem = min( minMem, currentSizes );
 maxMem = max( maxMem, currentSizes );

 totSizes += currentSizes;

 avgMem = totSizes/numThreads;

//fprintf(stderr,"%d         min             avg             max             total\n",numThreads);
 //fprintf(stderr,"%12.3f Mb %12.3f Mb %12.3f Mb %12.3f Mb\n",
 //        minMem/(1024.0*1024.0), avgMem/(1024.0*1024.0), 
 //        maxMem/(1024.0*1024.0), totSizes/(1024.0*1024.0));
 fprintf(stderr,"min %12.3f Mb avg %12.3f Mb max %12.3f Mb tot %12.3f Mb\n",
         minMem/(1024.0*1024.0), avgMem/(1024.0*1024.0), 
         maxMem/(1024.0*1024.0), totSizes/(1024.0*1024.0));
}

void
OneSproc::run(void *p)
{
 OneSproc *sp = (OneSproc *) p;
 ThreadManager *tman = sp->tman;

 while(1) {
     uspsema(sp->wait);
     if(sp->allTasks == 0) return;
     long initMem;
     double initTime;
     if(tman->timer) {
       initTime = getTime();
       initMem  = threadManager->getLocalMem();
     }
     int index;
     for(index = sp->myNum; index < sp->numTasks; index+= sp->step) {
       if(tman->single) sp->allTasks[0]->runFor(index);
       else sp->allTasks[index]->run();
     }
     if(tman->timer) {
       long finalMem = threadManager->getLocalMem();
       tman->timer->addTo(sp->myNum, finalMem-initMem, getTime()-initTime);
     }
     usvsema(sp->done);
 }
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
  if(isParal == 0)
    for(i=0; i < ntasks; ++i)
      td[i]->run();
  else {
    for(i = 0; i < numThreads-1; ++i) {
       allProc[i].allTasks = td;
       allProc[i].numTasks = ntasks;
    }
    single = 0;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);
 
    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td[index]->run();
 
    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
  }
}

void
ThreadManager::execParal(int ntasks, TaskDescr *td)
{
  timer = NULL;
  int i;
  if(isParal == 0)
    for(i = 0; i < ntasks; ++i) {
      td->runFor(i);
    }
  else {
    for(i = 0; i < numThreads-1; ++i) {
      allProc[i].allTasks = &td;
      allProc[i].numTasks = ntasks;
    }
    single = 1;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);

    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td->runFor(index);

    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
  }
}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr **td)
{

    timer = &thisTimer;
    int i;
    for(i = 0; i < numThreads-1; ++i) {
      allProc[i].allTasks = td;
      allProc[i].numTasks = ntasks;
    }
    single = 0;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);

    double initTime = getTime();
    long initMem  = threadManager->getLocalMem();
    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td[index]->run();
    long finalMem = threadManager->getLocalMem();
    timer->addTo(0, finalMem-initMem, getTime()-initTime);

    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);

}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr *td)
{
    timer = &thisTimer;
    int i;
    for(i = 0; i < numThreads-1; ++i) {
       allProc[i].allTasks = &td;
       allProc[i].numTasks = ntasks;
    }
    single = 0;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);

    double initTime = getTime();
    long initMem  = threadManager->getLocalMem();

    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td->runFor(index);

    long finalMem = threadManager->getLocalMem();
    timer->addTo(0, finalMem-initMem, getTime()-initTime);

    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
}


barrier_t *
ThreadManager::getBarrier()
{
 return new_barrier(usPtr);
}
