#include "LoadBalancer.h"

#include <Utils.d/Connectivity.h>

namespace Pita {

LoadBalancer::LoadBalancer(TaskCount totTasks, WorkerCount availWorkers, TaskCount maxLoad) :
  totalTasks_(totTasks),
  maxWorkload_(maxLoad),
  availableWorkers_(availWorkers),
  completedTasks_(0),
  taskToWorker_(NULL),
  workerToTasks_(NULL)
{
  int * arrayTasks = new int[totalTasks() + 1];
  for (int i = 0; i <= totalTasks(); ++i)
    arrayTasks[i] = i;

  int * targetWorker = new int[totalTasks()];
  
  int ratio, remain;
  if (totalTasks() < availableWorkers() * maxWorkload()) {
    ratio = totalTasks() / availableWorkers();
    remain = totalTasks() - ratio * availableWorkers();
  } else {
    ratio = maxWorkload();
    remain = 0;
  }

  int currentWorker = 0;
  int remainingSpotsForCurrentWorker = ratio + (remain > currentWorker ? 1 : 0);
  for (int currentTask = 0; currentTask < totalTasks(); ++currentTask) {
    if (remainingSpotsForCurrentWorker <= 0) {
      currentWorker = (currentWorker + 1) % availableWorkers();
      remainingSpotsForCurrentWorker = ratio + (remain > currentWorker ? 1 : 0);
    }
    targetWorker[currentTask] = currentWorker;
    --remainingSpotsForCurrentWorker;
  }

  taskToWorker_ = new Connectivity(totalTasks(), arrayTasks, targetWorker, 1);
  workerToTasks_ = taskToWorker_->reverse();
  workerToTasks_->sortTargets(); 
  updateFirstWaitingTask();
}
  
LoadBalancer::~LoadBalancer() {
  delete taskToWorker_;
  delete workerToTasks_;
}

void
LoadBalancer::completedTasksInc(TaskCount increment) {
  completedTasks_ += increment;
  updateFirstWaitingTask();
}

TaskCount
LoadBalancer::totalWorkload(WorkerRank pr) const {
  return TaskCount(workerToTasks_->num(pr));
}

TaskCount
LoadBalancer::currentGlobalWorkload() const {
  return firstWaitingTask() - completedTasks();
}

TaskCount
LoadBalancer::currentWorkload(WorkerRank pr) const {
  TaskCount result(0);
  for (TaskIterator it = tasks(pr); it && (*it < firstWaitingTask()); ++it) {
    if (*it >= completedTasks()) {
      ++result;
    }
  }
  return result;
}

LoadBalancer::TaskIterator
LoadBalancer::tasks(WorkerRank pr) const {
  return TaskIterator(this, pr); 
}

WorkerRank
LoadBalancer::worker(TaskRank tr) const {
  return (tr >= 0 && tr < totalTasks())    ?
         taskToWorker_->getTargetValue(tr) :
         -1;
}

TaskRank
LoadBalancer::firstCurrentTask() const { 
  return TaskRank(completedTasks());
}

void
LoadBalancer::updateFirstWaitingTask() {
  firstWaitingTask_ = std::min(totalTasks(), maxWorkload() *  availableWorkers() + completedTasks());
}

// LoadBalancer::TaskIterator

LoadBalancer::TaskIterator::TaskIterator(const LoadBalancer * parent, WorkerRank worker) :
  parent_(parent),
  worker_(worker),
  offset_(0)
{
  if (parent_->totalWorkload(worker_) == TaskCount(0))
    parent_ = NULL;
}

LoadBalancer::TaskIterator &
LoadBalancer::TaskIterator::operator++() {
  if (parent_) {
    if ((++offset_) >= parent_->totalWorkload(worker_))
      parent_ = NULL;
  }
  return *this;
}

LoadBalancer::TaskIterator
LoadBalancer::TaskIterator::operator++(int) {
  TaskIterator temp(*this);
  ++(*this);
  return temp;
}

TaskRank
LoadBalancer::TaskIterator::operator*() const {
  TaskRank result = parent_ ?
                    *(parent_->workerToTasks_->operator[](worker_) + offset_) :
                    -1;
  return result;
}

LoadBalancer::TaskIterator::operator bool() const {
  return (parent_ != NULL);
}

// Output

OStream & operator<<(OStream & out, const LoadBalancer & tmgr) {
  for (int worker = 0; worker < tmgr.availableWorkers(); ++worker) {
    out << "# " << worker << " ->";
    for (LoadBalancer::TaskIterator it = tmgr.tasks(worker); it; ++it) {
      out << ' ' << *it;
    }
    out << '\n';
  }
 return out; 
}

} // end namespace Pita
