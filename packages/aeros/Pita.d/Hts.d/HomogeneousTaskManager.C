#include "HomogeneousTaskManager.h"

#include "../InitialSeedTask.h"

namespace Pita { namespace Hts {

HomogeneousTaskManager::HomogeneousTaskManager(SliceMapping * mapping,
                                               RemoteState::MpiManager * commMgr,
                                               AffinePropagatorManager * propMgr,
                                               LinearProjectionNetwork * correctionMgr,
                                               JumpConvergenceEvaluator * jumpCvgMgr,
                                               LinSeedDifferenceEvaluator::Manager * jumpErrorMgr,
                                               SeedInitializer * initializer) :
  LinearTaskManager(IterationRank(0), mapping, propMgr, NULL, jumpCvgMgr, jumpErrorMgr, correctionMgr, commMgr),
  initializer_(initializer)
{
  scheduleIterationZero();
  updatePhaseIt();
}

void
HomogeneousTaskManager::iterationInc() {
  phases().clear();
  
  IterationRank nextIteration = iteration().next();
  
  correctionMgr()->prepareProjection();
  scheduleNormalIteration();
  updatePhaseIt();

  setIteration(nextIteration);
}

void
HomogeneousTaskManager::scheduleIterationZero() {
  scheduleSeedInitialization();
  scheduleFinePropagation();
}


void
HomogeneousTaskManager::scheduleSeedInitialization() {
  TaskList taskList;
  
  LinearLocalNetwork::MainSeedMap mainSeeds = network()->activeMainSeeds();
  for (LinearLocalNetwork::MainSeedMap::iterator it = mainSeeds.begin();
      it != mainSeeds.end();
      ++it) {
    SliceRank slice = it->first;
    Seed::Status initialSeedStatus = (slice == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
    taskList.push_back(new InitialSeedTask(it->second.ptr(), initializer_.ptr(), slice, initialSeedStatus));
  }
  
  Phase::Ptr finePropagation = phaseNew("Seed initialization", taskList);
  phases().push_back(finePropagation);
}

} /* end namespace Hts */ } /* end namespace Pita */
