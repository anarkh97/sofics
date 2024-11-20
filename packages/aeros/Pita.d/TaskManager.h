#ifndef PITA_TASKMANAGER_H
#define PITA_TASKMANAGER_H

#include "Fwk.h"
#include "Types.h"

#include "NamedTask.h"

#include <list>

namespace Pita {

class TaskManager : public Fwk::PtrInterface<TaskManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(TaskManager);

  // Iteration control
  IterationRank iteration() const { return iteration_; }
  virtual void iterationInc() = 0;

  // Phases
  class Phase;
  virtual Phase * phase() = 0;
  virtual void phaseInc() = 0;

protected:
  // Auxiliary implementation type
  template <typename E>
  class Iterator {
  public:
    E * operator*() { return it_->ptr(); }
    E * operator->() { return it_->ptr(); }
    Iterator<E> & operator++() { ++it_; return *this; }
    Iterator<E> operator++(int) { Iterator<E> temp(*this); this->operator++(); return temp; }
    operator bool() const { return it_ != it_end_; } 

    explicit Iterator(std::list<Fwk::Ptr<E> > & container) :
      it_(container.begin()),
      it_end_(container.end())
    {}

  private:
    typename std::list<Fwk::Ptr<E> >::iterator it_;
    typename std::list<Fwk::Ptr<E> >::iterator it_end_;
  };

public:
  class Phase : public NamedInterface {
  public:
    EXPORT_PTRINTERFACE_TYPES(Phase);

    typedef Iterator<NamedTask> TaskIterator;

    TaskIterator task() { return TaskIterator(task_); }

  protected:
    Phase(const String & name, const std::list<NamedTask::Ptr> & taskList) :
      NamedInterface(name),
      task_(taskList)
    {}
    
    // Pilfers the contents of taskList
    Phase(const String & name, std::list<NamedTask::Ptr> & taskList) :
      NamedInterface(name),
      task_()
    {
      task_.swap(taskList);
    }

    friend class TaskManager;

  private:
    typedef std::list<NamedTask::Ptr> TaskList;
    TaskList task_; 

    DISALLOW_COPY_AND_ASSIGN(Phase);
  };

protected:
  explicit TaskManager(IterationRank initialIter) :
    iteration_(initialIter)
  {}

  void setIteration(IterationRank iter) { iteration_ = iter; }

  typedef Phase::TaskList TaskList;

  static Phase * phaseNew(const String & name, TaskList & taskList) {
    return new Phase(name, taskList);
  }

  static Phase * phaseNew(const String & name, NamedTask::Ptr task) {
    TaskList taskList(1, task);
    return phaseNew(name, taskList);
  }

private:
  IterationRank iteration_;

  DISALLOW_COPY_AND_ASSIGN(TaskManager);
};

} /* end namespace Pita */

#endif /* PITA_TASKMANAGER_H */
