#include "NlTaskManager.h"

#include "../JumpBuilder.h"
#include "../ConcurrentBasisManager.h"
#include "../NearSymmetricSolver.h"
#include "../DynamStateReductor.h"
#include "../DynamStateReconstructor.h"

#include "../InitialSeedTask.h"

#include "NlLocalNetwork.h"

namespace Pita { namespace Hts {

NlTaskManager::NlTaskManager(SliceMapping * mapping, RemoteState::MpiManager * commMgr,
                             NlPropagatorManager * propagatorMgr,
                             SeedInitializer * seedInitializer,
                             GlobalStateSharing * basisUpdateMgr,
                             PostProcessing::Manager * postProcessingMgr,
                             JumpConvergenceEvaluator * jumpCvgMgr, NonLinSeedDifferenceEvaluator::Manager * jumpEvaluatorMgr,
                             double projectorTolerance, IterationRank lastIteration) :
  TaskManager(IterationRank(0)),
  mapping_(mapping),
  commMgr_(commMgr),
  propagatorMgr_(propagatorMgr),
  seedInitializer_(seedInitializer),
  postProcessingMgr_(postProcessingMgr),
  jumpCvgMgr_(jumpCvgMgr),
  jumpEvaluatorMgr_(jumpEvaluatorMgr),
  projectorTolerance_(projectorTolerance),
  phase_(NULL),
  continuation_(&NlTaskManager::noop),
  localNetwork_(NULL),
  basisUpdateMgr_(basisUpdateMgr),
  projectionNetwork_(NULL),
  lastIteration_(lastIteration)
{
  initialize();
  scheduleInitialSeed();
}

void
NlTaskManager::iterationInc() {
  setIteration(iteration().next());
  if (iteration() < lastIteration_) {
    basisUpdateMgr_->mappingIs(*mapping_);
  }
  scheduleConvergence();
}

NlTaskManager::Phase *
NlTaskManager::phase() {
  return phase_.ptr();
}

void
NlTaskManager::phaseInc() {
  (this->*continuation_)();
}

void
NlTaskManager::initialize() {
  // Projection (global data sharing)
  NlDynamOps::Ptr dynOps = propagatorMgr_->dynamOpsNew();
  projectionNetwork_ = new NlProjectionNetwork(basisUpdateMgr_.ptr(), dynOps.ptr(), commMgr_->vectorSize(), projectorTolerance_);

  // Parallel time-integration
  propagatorMgr_->postProcessingManagerIs(postProcessingMgr_.ptr());
  propagatorMgr_->concurrentBasisManagerIs(projectionNetwork_->concurrentMgr());
  propagatorMgr_->propagatedBasisManagerIs(projectionNetwork_->endBasisMgr());

  // Recurring tasks 
  localNetwork_ = new NlLocalNetwork(mapping_.ptr(),
                                     commMgr_.ptr(),
                                     propagatorMgr_.ptr(),
                                     projectionNetwork_->corrRedMgr(),
                                     projectionNetwork_->corrReconMgr(),
                                     projectionNetwork_->condensMgr(),
                                     projectionNetwork_->projBuildMgr(),
                                     jumpCvgMgr_.ptr(),
                                     jumpEvaluatorMgr_.ptr());

  basisUpdateMgr_->seedMgrIs(localNetwork_->seedManager());
}

void
NlTaskManager::scheduleNothing() {
  setPhase(NULL);
  setContinuation(&NlTaskManager::noop);
}

void
NlTaskManager::scheduleInitialSeed() {
  TaskList initialSeedInformation;

  typedef NlLocalNetwork::MainSeedMap MainSeedMap;
  MainSeedMap primalSeeds = localNetwork_->activeSeeds();
  for (MainSeedMap::iterator it = primalSeeds.begin(); it != primalSeeds.end(); ++it) {
    Seed::Ptr targetSeed = it->second;
    Seed::Status initialSeedStatus = (it->first == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
    NamedTask::Ptr task = new InitialSeedTask(targetSeed.ptr(), seedInitializer_.ptr(), it->first, initialSeedStatus);
    initialSeedInformation.push_back(task);
  }
  
  setPhase(phaseNew("Initial Seed Information", initialSeedInformation));
  setContinuation(&NlTaskManager::fillFirstProjectionBasis); 
}

void
NlTaskManager::fillFirstProjectionBasis() {
  int initSeedCount = (mapping_->activeSlices().value() / 2) + 1;
  for (int i = 0; i < initSeedCount; ++i) {
    projectionNetwork_->sharedProjectionBasis()->lastStateIs(seedInitializer_->initialSeed(SliceRank(i)));
  }

  scheduleProjectionBasisCondensation();
}

void
NlTaskManager::scheduleProjectionBasisCondensation() {
  setPhase(phaseNew("Projection Basis Condensation", localNetwork_->activeCondensations()));
  setContinuation(&NlTaskManager::scheduleFinePropagation);
}

void
NlTaskManager::scheduleFinePropagation() {
  setPhase(phaseNew("Parallel Fine Propagation", localNetwork_->activeFinePropagators()));
  setContinuation(&NlTaskManager::schedulePropagatedSeedSynchronization);
}

void
NlTaskManager::scheduleDataSharing() {
  if (iteration() < lastIteration_) {
    setPhase(phaseNew("Data Sharing", basisUpdateMgr_));
    setContinuation(&NlTaskManager::enrichProjectionBasis);
  } else {
    projectionNetwork_->sharedProjectionBasis()->stateBasisDel();
    scheduleProjectionBasisCondensation();
  }
}

void
NlTaskManager::scheduleConvergence() {
  setPhase(phaseNew("Convergence", jumpCvgMgr_));
  setContinuation(&NlTaskManager::applyConvergence);
}

void
NlTaskManager::applyConvergence() {
  localNetwork_->applyConvergenceStatus();
  scheduleProjectionBuilding();
}

void
NlTaskManager::schedulePropagatedSeedSynchronization() {
  setPhase(phaseNew("Propagated Seed Synchronization", localNetwork_->activePropagatedSeedSyncs()));
  setContinuation(&NlTaskManager::scheduleJumpEvaluation);
}

void
NlTaskManager::scheduleJumpEvaluation() {
  setPhase(phaseNew("Jump Evaluation", localNetwork_->activeJumpBuilders()));
  setContinuation(&NlTaskManager::scheduleNothing);
}

void
NlTaskManager::scheduleProjectionBuilding() {
  setPhase(phaseNew("Projection Building", localNetwork_->activeProjectionBuilders()));
  setContinuation(&NlTaskManager::scheduleCorrectionPropagation);
}

void
NlTaskManager::scheduleCorrectionPropagation() {
  setPhase(phaseNew("Correction Propagation", localNetwork_->activeCorrectionPropagators()));
  setContinuation(&NlTaskManager::scheduleSeedUpdate);
}

void
NlTaskManager::scheduleSeedUpdate() {
  setPhase(phaseNew("Seed Update", localNetwork_->activeSeedUpdaters()));
  setContinuation(&NlTaskManager::scheduleDataSharing);
}

void
NlTaskManager::enrichProjectionBasis() {
  projectionNetwork_->sharedProjectionBasis()->firstStateBasisIs(basisUpdateMgr_->consolidatedBasis());
  scheduleProjectionBasisCondensation();
}

} /* end namespace Hts */ } /* end namespace Pita */
