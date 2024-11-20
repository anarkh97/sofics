#ifndef PITA_HTS_NLLOCALNETWORK_H
#define PITA_HTS_NLLOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "SliceMapping.h"

#include "../Seed.h"
#include "../RemoteStateMpiImpl.h"

#include "NlPropagatorManager.h"
#include "../JumpBuilder.h"
#include "CorrectionReductor.h"
#include "CorrectionReconstructor.h"
#include "../SeedUpdater.h"

#include "NlProjectionNetwork.h"

#include "JumpConvergenceEvaluator.h"
#include "../SeedDifferenceEvaluator.h"

#include <map>
#include <list>
#include <functional>
#include <algorithm>

namespace Pita { namespace Hts {

// Half-open interval
class ActivationRange {
public:
  HalfSliceRank begin() const { return begin_; }
  HalfSliceRank end() const { return end_; }

  HalfSliceCount size() const { return end_ - begin_; }
  bool isEmpty() const { return begin_ == end_; }

  bool contains(const ActivationRange &other) const;
  bool isContainedIn(const ActivationRange &other) const;

  ActivationRange();
  ActivationRange(HalfSliceRank b, HalfSliceRank e);

private:
  HalfSliceRank begin_, end_;
};

inline
ActivationRange::ActivationRange() :
  begin_(0), end_(0)
{}

inline
ActivationRange::ActivationRange(HalfSliceRank b, HalfSliceRank e) :
  begin_(b), end_(e)
{
  if (end_ < begin_) {
    end_ = begin_;
  }
}

inline
bool
ActivationRange::contains(const ActivationRange &other) const {
  return other.isEmpty() || (begin() <= other.begin() && end() >= other.end());
}

inline
bool
ActivationRange::isContainedIn(const ActivationRange &other) const {
  return other.contains(*this);
}

// Activation rank and dependency interval for a task
class ActivationInfo {
public:
  HalfSliceRank rankInPhase() const { return rankInPhase_; }
  HalfSliceRank deactivation() const { return deactivation_; }
  HalfSliceRank activation() const { return activation_; }

  // Implicit conversion from HalfSliceRank
  ActivationInfo(HalfSliceRank rankInPhase);
  ActivationInfo(HalfSliceRank rankInPhase, HalfSliceCount dependency);
  ActivationInfo(HalfSliceRank rankInPhase, HalfSliceCount deactivationOffset, HalfSliceCount activationOffset);
 
  bool isActiveIn(const ActivationRange &) const;

  struct Comparator;
private:
  HalfSliceRank rankInPhase_;
  HalfSliceRank deactivation_;
  HalfSliceRank activation_;
};

inline
ActivationInfo::ActivationInfo(HalfSliceRank rankInPhase) :
  rankInPhase_(rankInPhase), activation_(rankInPhase), deactivation_(rankInPhase)
{}

inline
ActivationInfo::ActivationInfo(HalfSliceRank rankInPhase, HalfSliceCount dependency) :
  rankInPhase_(rankInPhase), activation_(rankInPhase), deactivation_(rankInPhase)
{
  HalfSliceRank bound = rankInPhase + dependency;
  (bound > rankInPhase_ ? activation_ : deactivation_) = bound;
}

inline
ActivationInfo::ActivationInfo(HalfSliceRank rankInPhase, HalfSliceCount deactivationOffset, HalfSliceCount activationOffset) :
  rankInPhase_(rankInPhase), deactivation_(rankInPhase + deactivationOffset), activation_(rankInPhase + activationOffset)
{}

inline
bool
ActivationInfo::isActiveIn(const ActivationRange &range) const {
  return deactivation() >= range.begin() && activation() < range.end();
}

// Total ordering
struct ActivationInfo::Comparator : public std::binary_function<ActivationInfo, ActivationInfo, bool> {
  bool operator()(const ActivationInfo & a, const ActivationInfo & b) const {
    if (a.deactivation() == b.deactivation()) {
      if (a.rankInPhase() == b.rankInPhase()) {
        return a.activation() < b.activation();
      }
      return a.rankInPhase() < b.rankInPhase();
    }
    return a.deactivation() < b.deactivation();
  }
};

inline
bool
includedIn(const ActivationInfo & a, const ActivationInfo & b) {
  return a.activation() >= b.activation() && a.deactivation() <= b.deactivation();
}


class NlLocalNetwork : public Fwk::PtrInterface<NlLocalNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlLocalNetwork);

  Seed::Manager * seedManager() { return seedMgr_.ptr(); }
  
  // Tasks
  typedef std::list<NamedTask::Ptr> TaskList;
  
  TaskList activeFinePropagators() const { return toTaskList(finePropagators_[activeParity()]); }
  TaskList activePropagatedSeedSyncs() const { return toTaskList(propagatedSeedSyncs_[activeParity()]); }
  TaskList activeJumpBuilders() const { return toTaskList(jumpBuilders_[activeParity()]); }
  TaskList activeCorrectionPropagators() const { return toTaskList(correctionPropagators_[activeParity()]); }
  TaskList activeSeedUpdaters() const { return toTaskList(seedUpdaters_[activeParity()]); }

  TaskList activeCondensations() const { return toTaskList(condensations_[activeParity()]); }
  TaskList activeProjectionBuilders() const { return toTaskList(projectionBuilders_[activeParity()]); }
  
  typedef std::map<SliceRank, Seed::Ptr> MainSeedMap;
  MainSeedMap activeSeeds() const { return seeds_[activeParity()]; }
  
  void applyConvergenceStatus();

  NlLocalNetwork(SliceMapping * mapping,
                 RemoteState::MpiManager * commMgr,
                 NlPropagatorManager * propMgr,
                 CorrectionReductor::Manager * corrRedMgr,
                 CorrectionReconstructor::Manager * corrReconMgr,
                 BasisCondensationManager * condensMgr,
                 ProjectionBuildingFactory * projBuildMgr,
                 JumpConvergenceEvaluator * jumpCvgMgr,
                 NonLinSeedDifferenceEvaluator::Manager * jumpEvalMgr);

private:
  // Build recurrent tasks
  void init();

  void addForwardSlice(HalfSliceRank);
  void addBackwardSlice(HalfSliceRank);

  void addMainSeed(HalfSliceRank seedRank);

  void addForwardPropagation(HalfSliceRank mainSeedRank);
  void addBackwardPropagation(HalfSliceRank mainSeedRank);

  void addPropagatedSeedSend(HalfSliceRank seedRank);
  void addPropagatedSeedRecv(HalfSliceRank seedRank);

  void addJumpBuilder(HalfSliceRank seedRank);
 
  void addProjectionBuilding(HalfSliceRank seedRank);

  void addCorrectionReductor(HalfSliceRank seedRank);
  void addCorrectionReconstructor(HalfSliceRank seedRank);
  void addCorrectionSend(HalfSliceRank seedRank);
  void addCorrectionRecv(HalfSliceRank seedRank);
  void addReducedCorrectionSend(HalfSliceRank seedRank);
  void addReducedCorrectionRecv(HalfSliceRank seedRank);

  void addSeedUpdater(HalfSliceRank seedRank);

  // Create new task
  NamedTask::Ptr forwardHalfSliceNew(HalfSliceRank sliceRank);
  NamedTask::Ptr backwardHalfSliceNew(HalfSliceRank sliceRank);
  
  NamedTask::Ptr propagatedSeedSendNew(HalfSliceRank seedRank);
  NamedTask::Ptr propagatedSeedRecvNew(HalfSliceRank seedRank);
  NamedTask::Ptr jumpBuilderNew(HalfSliceRank seedRank);

  NamedTask::Ptr correctionReductorNew(HalfSliceRank seedRank);
  NamedTask::Ptr correctionReconstructorNew(HalfSliceRank seedRank);
  NamedTask::Ptr correctionSendNew(HalfSliceRank seedRank);
  NamedTask::Ptr correctionRecvNew(HalfSliceRank seedRank);
  NamedTask::Ptr reducedCorrectionSendNew(HalfSliceRank seedRank);
  NamedTask::Ptr reducedCorrectionRecvNew(HalfSliceRank seedRank);

  NamedTask::Ptr seedUpdater(HalfSliceRank seedRank);
  NamedTask::Ptr seedUpdaterNew(HalfSliceRank seedRank);
  
  // Cpu map
  CpuRank localCpu() const { return commMgr_->localCpu(); }
  CpuRank hostCpu(HalfSliceRank slice) const { return mapping_->hostCpu(slice); }

  // Current computational state
  HalfSliceRank firstActiveSlice() const { return mapping_->firstActiveSlice(); }
  HalfSliceRank firstInactiveSlice() const { return mapping_->firstInactiveSlice(); }
  ActivationRange activeRange() const { return ActivationRange(mapping_->firstActiveSlice(), mapping_->firstInactiveSlice()); }
  bool isCurrentlyActive(const ActivationInfo &info) const { return info.isActiveIn(activeRange()); }

  static int parity(HalfSliceRank sliceRank) { return sliceRank.value() % 2; } 
  int activeParity() const { return parity(firstActiveSlice()); } 
 
  // Seed instances 
  SharedState<DynamState> * fullSeedGet(const SeedId & id);
  SharedState<Vector> * reducedSeedGet(const SeedId & id);

  // Task managers
  NlPropagatorManager * propMgr() { return propMgr_.ptr(); }
  JumpBuilder::Manager * jumpMgr() { return jumpMgr_.ptr(); }
  NonLinSeedDifferenceEvaluator::Manager * jumpEvalMgr() { return jumpEvalMgr_.ptr(); }
  SeedUpdater::Manager * seedUpMgr() { return seedUpMgr_.ptr(); }
  CorrectionReductor::Manager * corrRedMgr() { return corrRedMgr_.ptr(); }
  CorrectionReconstructor::Manager * corrReconMgr() { return corrReconMgr_.ptr(); }
  BasisCondensationManager * condensMgr() { return condensMgr_.ptr(); }
  ProjectionBuildingFactory * projBuildMgr() { return projBuildMgr_.ptr(); }


private:
  SliceMapping::Ptr mapping_;
  RemoteState::Manager::Ptr commMgr_;
 
  Seed::Manager::Ptr seedMgr_;
  ReducedSeed::Manager::Ptr reducedSeedMgr_;

  NlPropagatorManager::Ptr propMgr_;
  
  JumpBuilder::Manager::Ptr jumpMgr_;
  SeedUpdater::Manager::Ptr seedUpMgr_;
  CorrectionReductor::Manager::Ptr corrRedMgr_;
  CorrectionReconstructor::Manager::Ptr corrReconMgr_;

  BasisCondensationManager::Ptr condensMgr_;
  ProjectionBuildingFactory::Ptr projBuildMgr_;

  JumpConvergenceEvaluator::Ptr jumpCvgMgr_;

  NonLinSeedDifferenceEvaluator::Manager::Ptr jumpEvalMgr_;
  
  MainSeedMap seeds_[2];
  
  typedef std::map<HalfSliceRank, Seed::Ptr> SeedMap;
  SeedMap seedCorrection_;

  typedef std::map<ActivationInfo, NamedTask::Ptr, ActivationInfo::Comparator> TaskMap;
  TaskMap finePropagators_[2];
  TaskMap propagatedSeedSyncs_[2];
  TaskMap jumpBuilders_[2];
  TaskMap correctionPropagators_[2];
  TaskMap seedUpdaters_[2];
  TaskMap condensations_[2];
  TaskMap projectionBuilders_[2];

  static TaskList toTaskList(const TaskMap &);

  DISALLOW_COPY_AND_ASSIGN(NlLocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLLOCALNETWORK_H */
