#ifndef PITA_HTS_NLTASKMANAGER_H
#define PITA_HTS_NLTASKMANAGER_H

#include "../TaskManager.h"

#include "SliceMapping.h"
#include "NlPropagatorManager.h"
#include "../SeedInitializer.h"
#include "../PostProcessingManager.h"
#include "../SeedDifferenceEvaluator.h"
#include "../RemoteStateMpiImpl.h"
#include "NlProjectionNetwork.h"
#include "JumpConvergenceEvaluator.h"
#include "GlobalStateSharing.h"

namespace Pita { namespace Hts {

class NlLocalNetwork;

class NlTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NlTaskManager);

  virtual void iterationInc(); // overriden
  virtual Phase * phase();     // overriden
  virtual void phaseInc();     // overriden

  CpuRank localCpu() const { return commMgr_->localCpu(); }

  NlTaskManager(SliceMapping * mapping, RemoteState::MpiManager * commMgr,
                NlPropagatorManager * propagatorMgr,
                SeedInitializer * seedInitializer,
                GlobalStateSharing * basisUpdateMgr,
                PostProcessing::Manager * postProcessingMgr,
                JumpConvergenceEvaluator * jumpCvgMgr, NonLinSeedDifferenceEvaluator::Manager * jumpEvaluatorMgr,
                double projectorTolerance, IterationRank lastIteration);

protected:
  void initialize();

  // Execution
  void setPhase(Phase * p) { phase_ = p; }

  typedef void (NlTaskManager::*Continuation)();
  void setContinuation(Continuation c) { continuation_ = c; }

  void noop() {}
  void scheduleNothing();

  void scheduleInitialSeed();
  void fillFirstProjectionBasis(); 
  void scheduleProjectionBasisCondensation();
  void scheduleFinePropagation();
  void schedulePropagatedSeedSynchronization();
  void scheduleJumpEvaluation();
  void scheduleDataSharing();
  void scheduleConvergence();
  void applyConvergence();
  void scheduleProjectionBuilding();
  void scheduleCorrectionPropagation();
  void scheduleSeedUpdate();
  void enrichProjectionBasis();

  static Phase * phaseNew(const String & name, const TaskList & tmpList) {
    TaskList newList(tmpList);
    return TaskManager::phaseNew(name, newList);
  }
  using TaskManager::phaseNew; 

private:
  SliceMapping::Ptr mapping_;
  RemoteState::MpiManager::Ptr commMgr_;

  NlPropagatorManager::Ptr propagatorMgr_;
  SeedInitializer::Ptr seedInitializer_;
  PostProcessing::Manager::Ptr postProcessingMgr_;
  JumpConvergenceEvaluator::Ptr jumpCvgMgr_;
  NonLinSeedDifferenceEvaluator::Manager::Ptr jumpEvaluatorMgr_;

  double projectorTolerance_;

  Phase::Ptr phase_;
  Continuation continuation_;

  Fwk::Ptr<NlLocalNetwork> localNetwork_;
  
  GlobalStateSharing::Ptr basisUpdateMgr_;
  NlProjectionNetwork::Ptr projectionNetwork_;
  IterationRank lastIteration_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLTASKMANAGER_H */
