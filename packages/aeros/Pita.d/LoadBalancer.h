#ifndef PITA_LOADBALANCER_H
#define PITA_LOADBALANCER_H

#include "Fwk.h"

class Connectivity;

namespace Pita {

typedef int TaskCount;
typedef int TaskRank;
typedef int WorkerCount;
typedef int WorkerRank;

class LoadBalancer : public Fwk::PtrInterface<LoadBalancer> {
public:
  EXPORT_PTRINTERFACE_TYPES(LoadBalancer);
  
  class TaskIterator;

  TaskCount totalTasks() const { return totalTasks_; }
  TaskCount maxWorkload() const { return maxWorkload_; }
  WorkerCount availableWorkers() const { return availableWorkers_; }

  TaskCount completedTasks() const { return completedTasks_; }
  void completedTasksInc(TaskCount increment = TaskCount(1));

  TaskCount totalWorkload(WorkerRank pr) const;
  TaskCount currentGlobalWorkload() const; 
  TaskCount currentWorkload(WorkerRank pr) const; 

  TaskIterator tasks(WorkerRank cr) const;
  WorkerRank worker(TaskRank tr) const;

  TaskRank firstCurrentTask() const;
  TaskRank firstWaitingTask() const { return firstWaitingTask_; }

  static Ptr New(TaskCount totalTasks, WorkerCount availableWorkers, TaskCount maxWorkload) {
    return new LoadBalancer(totalTasks, availableWorkers, maxWorkload);
  }

protected:
  LoadBalancer(TaskCount totalTasks, WorkerCount availableWorkers, TaskCount maxWorkload);
  virtual ~LoadBalancer();

  void updateFirstWaitingTask();
  
  friend class TaskIterator;

private:
  TaskCount totalTasks_;
  TaskCount maxWorkload_;
  WorkerCount availableWorkers_;
  TaskCount completedTasks_;

  Connectivity * taskToWorker_;
  Connectivity * workerToTasks_;

  TaskRank firstWaitingTask_;

  DISALLOW_COPY_AND_ASSIGN(LoadBalancer);
};

class LoadBalancer::TaskIterator {
public:
  TaskIterator & operator++();
  TaskIterator operator++(int);
  TaskRank operator*() const;
  operator bool() const;

  // Use default copy constructor, assignement operator

protected:
  TaskIterator(const LoadBalancer * parent, WorkerRank worker);

  friend class LoadBalancer;

private:
  LoadBalancer::PtrConst parent_;
  WorkerRank worker_;
  TaskCount offset_;
};

OStream & operator<<(OStream & out, const LoadBalancer & tmgr);

} // end namespace Pita

#endif /* PITA_LOADBALANCER_H */
