#ifndef PITA_HTS_LINEARTASKMANAGER_H
#define PITA_HTS_LINEARTASKMANAGER_H

#include "../TaskManager.h"

#include "LinearLocalNetwork.h"

#include "JumpConvergenceEvaluator.h"
#include "LinearProjectionNetwork.h"
#include "../RemoteStateMpiImpl.h"

#include <memory>

namespace Pita { namespace Hts {

class LinearTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearTaskManager);

  virtual Phase * phase() { return phaseIt_ ? *phaseIt_ : NULL; }
  virtual void phaseInc() { ++phaseIt_; }

protected:
  LinearTaskManager(IterationRank initialIteration,
                    SliceMapping * mapping,
                    AffinePropagatorManager * propMgr,
                    CorrectionPropagator<DynamState>::Manager * fullCorrMgr,
                    JumpConvergenceEvaluator * jumpCvgMgr,
                    LinSeedDifferenceEvaluator::Manager * jumpErrorMgr,
                    LinearProjectionNetwork * correctionMgr,
                    RemoteState::MpiManager * commMgr);

  LinearLocalNetwork * network() { return network_.ptr(); }
  LinearProjectionNetwork * correctionMgr() { return correctionMgr_.ptr(); }
  RemoteState::MpiManager * commMgr() { return commMgr_.ptr(); }
  JumpConvergenceEvaluator * jumpCvgMgr() { return jumpCvgMgr_.ptr(); }

  typedef std::vector<Fwk::Ptr<Phase> > PhaseList;
  PhaseList & phases() { return phase_; } // TODO remove ?

  class PhaseIteratorImpl {
  public:
    virtual Phase * get() = 0;
    virtual void next() = 0;
    virtual bool hasNext() const = 0;
    virtual PhaseIteratorImpl * clone() const = 0;
  };

  class PhaseIterator {
  public:
    Phase * operator*() { return pimpl_->get(); }
    Phase * operator->() { return pimpl_->get(); }
    PhaseIterator & operator++() { pimpl_->next(); return *this; }
    PhaseIterator operator++(int) { PhaseIterator temp(*this); pimpl_->next(); return temp; }
    operator bool() const { return pimpl_->hasNext(); }

    explicit PhaseIterator(PhaseIteratorImpl * pimpl) : pimpl_(pimpl) {}
    PhaseIterator(const PhaseIterator & other) : pimpl_(other.pimpl_->clone()) {}
    PhaseIterator & operator=(const PhaseIterator & other) { pimpl_.reset(other.pimpl_->clone()); return *this; }

  private:
    std::unique_ptr<PhaseIteratorImpl> pimpl_;
  };

  class HtsPhaseIteratorImpl : public PhaseIteratorImpl {
  public:
    virtual Phase * get() { return it_->ptr(); }
    virtual void next() { ++it_; }
    virtual bool hasNext() const { return it_ != it_end_; }
    virtual HtsPhaseIteratorImpl * clone() const { return new HtsPhaseIteratorImpl(*this); }

    explicit HtsPhaseIteratorImpl(PhaseList & container) :
      it_(container.begin()),
      it_end_(container.end())
    {}

  private:
    PhaseList::iterator it_;
    PhaseList::iterator it_end_;
  };

  void updatePhaseIt() { phaseIt_ = PhaseIterator(new HtsPhaseIteratorImpl(phase_)); }

  void scheduleNormalIteration();
  void scheduleFinePropagation();
  void scheduleCorrection();

  void schedulePhase(const String & phaseName, const LinearLocalNetwork::TaskList & networkTaskList);

private:
  LinearLocalNetwork::Ptr network_;
  JumpConvergenceEvaluator::Ptr jumpCvgMgr_;
  LinearProjectionNetwork::Ptr correctionMgr_;
  RemoteState::MpiManager::Ptr commMgr_;

  PhaseList phase_;
  PhaseIterator phaseIt_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARTASKMANAGER_H */
