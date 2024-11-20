#ifndef PITA_STD_NONHOMOGENEOUSTASKMANAGER_H
#define PITA_STD_NONHOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"

namespace Pita { namespace Std {

class NonHomogeneousTaskManager : public LinearTaskManager<NonHomogeneousTaskManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousTaskManager);

  virtual void iterationInc(); // overriden

  NonHomogeneousTaskManager(SliceMapping * mapping, RemoteState::Manager * commMgr,
                            LinearPropagatorManager * propagatorMgr, LinearProjectionNetwork * projectionMgr,
                            JumpConvergenceEvaluator * jumpCvgEval,
                            LinSeedDifferenceEvaluator::Manager * jumpOutMgr,
                            SeedInitializer * seedInit,
                            CorrectionPropagator<DynamState>::Manager * fullCorrPropMgr);

protected:
  void schedulePrecomputation();
  void scheduleCoarseCorrectionPropagation();

private:
  CorrectionPropagator<DynamState>::Manager::Ptr fullCorrPropMgr_;
  
  DISALLOW_COPY_AND_ASSIGN(NonHomogeneousTaskManager);
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_NONHOMOGENEOUSTASKMANAGER_H */
