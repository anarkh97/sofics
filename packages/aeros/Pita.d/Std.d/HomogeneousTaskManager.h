#ifndef PITA_STD_HOMOGENEOUSTASKMANAGER_H
#define PITA_STD_HOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"

namespace Pita { namespace Std {

class HomogeneousTaskManager : public LinearTaskManager<HomogeneousTaskManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(HomogeneousTaskManager);

  HomogeneousTaskManager(SliceMapping * mapping,
                         RemoteState::Manager * commMgr,
                         LinearPropagatorManager * propagatorMgr,
                         LinearProjectionNetwork * projectionMgr,
                         JumpConvergenceEvaluator * jumpCvgEval,
                         LinSeedDifferenceEvaluator::Manager * jumpOutMgr,
                         SeedInitializer * initializer);
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_HOMOGENEOUSTASKMANAGER_H */
