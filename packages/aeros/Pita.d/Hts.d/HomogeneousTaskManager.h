#ifndef PITA_HTS_HOMOGENEOUSTASKMANAGER_H
#define PITA_HTS_HOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"

#include "../SeedInitializer.h"

namespace Pita { namespace Hts {

class HomogeneousTaskManager : public LinearTaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(HomogeneousTaskManager);

  // overriden
  virtual void iterationInc();

  HomogeneousTaskManager(SliceMapping * mapping,
                         RemoteState::MpiManager * commMgr,
                         AffinePropagatorManager * propMgr,
                         LinearProjectionNetwork * correctionMgr,
                         JumpConvergenceEvaluator * jumpCvgMgr,
                         LinSeedDifferenceEvaluator::Manager * jumpErrorMgr,
                         SeedInitializer * initializer);

protected:
  void scheduleIterationZero();
  void scheduleSeedInitialization();

private:
  SeedInitializer::Ptr initializer_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HOMOGENEOUSTASKMANAGER_H */
