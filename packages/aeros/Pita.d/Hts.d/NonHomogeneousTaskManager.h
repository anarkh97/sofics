#ifndef PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H
#define PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"
#include "../SeedInitializer.h"

namespace Pita { namespace Hts {

class NonHomogeneousTaskManager : public LinearTaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousTaskManager);

  virtual void iterationInc(); //overriden

  NonHomogeneousTaskManager(SliceMapping * mapping,
                            RemoteState::MpiManager * commMgr,
                            AffinePropagatorManager * propMgr,
                            LinearProjectionNetwork * correctionMgr,
                            JumpConvergenceEvaluator * jumpCvgMgr,
                            LinSeedDifferenceEvaluator::Manager * jumpErrorMgr,
                            SeedInitializer * initializer,
                            CorrectionPropagator<DynamState>::Manager * fullCorrMgr);

protected:
  void schedulePreIteration();
  void scheduleIterationZero();
  void scheduleBasicSeedInitialization();

private:
  SeedInitializer::Ptr seedInit_;
};

} /* end namespace Hts */ } /* end namespace Pta */

#endif /* PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H */
