#include "TimedExecution.h"

namespace Pita {

TimedExecution::TimedExecution(TaskManager * taskMgr) :
  taskMgr_(taskMgr),
  currentIteration_(taskMgr->iteration()),
  timer_(new HierarchicalTimer())
{}

TimedExecution::~TimedExecution() {
  log() << "\n" << *timer_;
}

void
TimedExecution::targetIterationIs(IterationRank lastIteration) {
  timer_->startAll();
  timer_->start("Main iterations");
  while (currentIteration() <= lastIteration) {
#ifndef NDEBUG
    log() << "-> Iteration: " << currentIteration() << "\n";
#endif /* NDEBUG */
    timer_->newIteration();

    while (TaskManager::Phase::Ptr currentPhase = taskMgr_->phase()) {
#ifndef NDEBUG
      log() << "  -> Phase: " << currentPhase->name() << "\n";
#endif /* NDEBUG */
      timer_->start(currentPhase->name());
      for (TaskManager::Phase::TaskIterator task_it = currentPhase->task(); task_it; ++task_it) {
        NamedTask::Ptr currentTask = *task_it;
#ifndef NDEBUG
        log() << "    -> Task: " << currentTask->name() << "\n";
#endif /* NDEBUG */
        timer_->start(currentTask->name());
        currentTask->iterationIs(currentIteration());
        timer_->stop();
      } 
      timer_->stop();
      taskMgr_->phaseInc();
    }

    taskMgr_->iterationInc();
    currentIteration_ = taskMgr_->iteration();
  }
  timer_->stop();
  timer_->stopAll();
}

} /* end namespace Pita */
