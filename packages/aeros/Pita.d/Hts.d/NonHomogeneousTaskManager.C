#include "NonHomogeneousTaskManager.h"

#include "../InitialSeedTask.h"

namespace Pita { namespace Hts {

NonHomogeneousTaskManager::NonHomogeneousTaskManager(SliceMapping * mapping,
                                                     RemoteState::MpiManager * commMgr,
                                                     AffinePropagatorManager * propMgr,
                                                     LinearProjectionNetwork * correctionMgr,
                                                     JumpConvergenceEvaluator * jumpCvgMgr,
                                                     LinSeedDifferenceEvaluator::Manager * jumpErrorMgr,
                                                     SeedInitializer * initializer,
                                                     CorrectionPropagator<DynamState>::Manager * fullCorrMgr) :
  LinearTaskManager(IterationRank(-1), mapping, propMgr, fullCorrMgr, jumpCvgMgr, jumpErrorMgr, correctionMgr, commMgr),
  seedInit_(initializer)
{
  schedulePreIteration();
  updatePhaseIt();
}

void
NonHomogeneousTaskManager::iterationInc() {
  phases().clear();
  
  IterationRank nextIteration = iteration().next();

  if (nextIteration == IterationRank(0)) {
    scheduleIterationZero();
  } else {
    correctionMgr()->prepareProjection();
    scheduleNormalIteration();
  }

  updatePhaseIt();  
  setIteration(nextIteration);
}

void
NonHomogeneousTaskManager::schedulePreIteration() {
  scheduleBasicSeedInitialization();
  schedulePhase("Affine term precomputation", network()->halfTimeSlices());
  schedulePhase("Propagated seed synchronization", network()->activeLeftSeedSyncs());
  schedulePhase("Trivial jump evaluation", network()->activeJumpAssemblers());
}

void
NonHomogeneousTaskManager::scheduleBasicSeedInitialization() {
  TaskList taskList;

  DynamState zeroState(seedInit_->vectorSize(), 0.0);

  LinearLocalNetwork::SeedMap mainSeeds = network()->mainSeeds();
  for (LinearLocalNetwork::SeedMap::iterator it = mainSeeds.begin();
      it != mainSeeds.end();
      ++it) {
  
    Seed * target = it->second.ptr();
    if (it->first.value() % 2 == 0) {
      SliceRank rank(it->first.value() / 2);
      Seed::Status status = (rank == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
      taskList.push_back(new InitialSeedTask(target, seedInit_.ptr(), rank, status));
    } else {
      target->stateIs(zeroState);
      target->statusIs(Seed::SPECIAL);
      target->iterationIs(iteration());
    }
  }

  Phase::Ptr zeroSeedInitialization = phaseNew("Seed initialization", taskList);
  phases().push_back(zeroSeedInitialization);
}

void
NonHomogeneousTaskManager::scheduleIterationZero() {
  network()->convergedSlicesInc(HalfSliceCount(1));
  network()->applyConvergenceStatus();

  schedulePhase("Coarse propagation", network()->activeCoarseTimeSlices());
  schedulePhase("Correction synchronization", network()->activeFullCorrectionSyncs());
  schedulePhase("Seed update", network()->activeSeedAssemblers());
  scheduleFinePropagation(); 
}

} /* end namespace Hts */ } /* end namespace Pta */
