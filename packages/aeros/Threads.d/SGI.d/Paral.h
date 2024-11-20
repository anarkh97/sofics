#ifndef _PARAL_H_
#define _PARAL_H_

#include <ulocks.h>

// function to return memory used whether
// threadmanager has been created or not
long memoryUsed();

class ThreadLock {
      ulock_t lockV;
   public:
      ThreadLock();
      void lock();
      void unlock();
};

/***************************************************************************
 * TaskDescr is a purely abstract class. It should be subclassed
 * to specific task objects that contain sufficient data to refer to the
 * work to be done, and overload the run() routine to perform the work
 ***************************************************************************/

class TaskDescr {
   public:
	virtual void run() {};
	virtual void runFor(int) = 0;
};

class ThreadManager;
class DistTimer;

class OneSproc {
   OneSproc *next;
   TaskDescr **allTasks;
   int numTasks;
   int step;
   usema_t *wait;
   usema_t *done;
   ThreadManager *tman;
   int myNum;

   static void run(void *);
   friend class ThreadManager;
};

/***************************************************************************
 * The ThreadManager class manages the creation, and deletion of a given
 * number of threads. Once it has been created, it can be given a list
 * of TaskDescr's that it will run on the available threads. In the default
 * case of a sequential machine, it will loop over the tasks to execute them
 *
 * NB: Because of this paradigm in which we do not know how many tasks are
 *     really running in parallel, a Task should not intend to communicate
 *     with another. It is only allowed to use Locks to insure proper
 *     update of shared data between tasks
 ***************************************************************************/
class ThreadManager {
	int numThreads;		// Number of threads
        int single;
	OneSproc *firstProc; 	// List of threads
        OneSproc *allProc;
        ulock_t sprocListLock;
        usema_t *readyProc;
        usema_t *allDone;
        DistTimer *timer;
   protected:
	static void threadStart(void *);
        void dispatch(TaskDescr *task);
   public:
	ThreadManager(int);	// Create a Manager with n threads
	~ThreadManager();
	void execTasks(int, TaskDescr **);
				// run the n given Tasks
        void execTasks(int, TaskDescr *);
        void execParal(int, TaskDescr **);
        void execParal(int, TaskDescr *);
        void execTimedParal(DistTimer &, int, TaskDescr **);
        void execTimedParal(DistTimer &, int, TaskDescr *);
        void memUsage();
        void memRequests();
        long getLocalMem();
	long memoryUsed();
        int numThr() { return numThreads; }
        barrier_t *getBarrier();

        friend class OneSproc;
};

extern ThreadManager *threadManager;

#endif
