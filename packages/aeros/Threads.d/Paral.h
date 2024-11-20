#ifndef _PARAL_H_
#define _PARAL_H_

#include <mutex>
#include <functional>

#if defined(_OPENMP)
#include <omp.h>
#endif

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
#if defined(sgi) &&  !defined(_OPENMP)
	usema_t *wait;
   usema_t *done;
#endif
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
	DistTimer *timer;
public:
	ThreadManager(int);	// Create a Manager with n threads
	~ThreadManager() = default;

	void execParal(int, TaskDescr **td);
	void execParal(int, TaskDescr *td);
	void execTimedParal(DistTimer &, int, TaskDescr **);
	void execTimedParal(DistTimer &, int, TaskDescr *);

	/** \brief Call in parallel a function with one integer argument. */
	template <typename Ftor>
	void callParal(int nObjects, const Ftor &f) {
		struct TD : public TaskDescr {
			const Ftor &f;
			TD(const Ftor &f) : f(f) {}
			void run(){}
			void runFor(int i) { f(i); }
		};
		TD td{f};
		execParal(nObjects, &td);
	}

	long getLocalMem();
	long memoryUsed();
	int numThr() { return numThreads; }

	friend class OneSproc;
};

extern ThreadManager *threadManager;

#endif
