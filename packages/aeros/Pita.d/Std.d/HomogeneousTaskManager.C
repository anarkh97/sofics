#include "HomogeneousTaskManager.h"

namespace Pita { namespace Std {

HomogeneousTaskManager::HomogeneousTaskManager(SliceMapping * mapping,
                                               RemoteState::Manager * commMgr,
                                               LinearPropagatorManager * propagatorMgr,
                                               LinearProjectionNetwork * projectionMgr,
                                               JumpConvergenceEvaluator * jumpCvgEval,
                                               LinSeedDifferenceEvaluator::Manager * jumpOutMgr,
                                               SeedInitializer * initializer) :
  LinearTaskManager<HomogeneousTaskManager>(IterationRank(0), mapping, commMgr, propagatorMgr, projectionMgr, jumpCvgEval, jumpOutMgr, initializer)
{
  setContinuation(&HomogeneousTaskManager::scheduleFinePropagation);
}

} /* end namespace Std */ } /* end namespace Pita */
