#ifndef PITA_HTS_REDUCEDCORRECTIONMANAGER_H
#define PITA_HTS_REDUCEDCORRECTIONMANAGER_H

#include "Fwk.h"

#include "../JumpBuilder.h"
#include "../JumpProjection.h"
#include "../UpdatedSeedAssembler.h"
#include "../CorrectionPropagator.h"


namespace Pita { namespace Hts {

class LinearProjectionNetwork;

class ReducedCorrectionManager : public Fwk::PtrInterface<ReducedCorrectionManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedCorrectionManager);

  JumpProjection::Manager * jumpProjMgr() const { return jumpProjMgr_.ptr(); }
  CorrectionPropagator<Vector>::Manager * rcpMgr() const { return rcpMgr_.ptr(); }
  CorrectionPropagator<DynamState>::Manager * fcpMgr() const { return fcpMgr_.ptr(); }
  UpdatedSeedAssembler::Manager * usaMgr() const { return usaMgr_.ptr(); }
  
  ReducedCorrectionManager(LinearProjectionNetwork * correctionMgr, CorrectionPropagator<DynamState>::Manager * fcpMgr);

private:
  JumpProjection::Manager::Ptr jumpProjMgr_;
  CorrectionPropagator<Vector>::Manager::Ptr rcpMgr_;
  CorrectionPropagator<DynamState>::Manager::Ptr fcpMgr_;
  UpdatedSeedAssembler::Manager::Ptr usaMgr_;

  DISALLOW_COPY_AND_ASSIGN(ReducedCorrectionManager);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_REDUCEDCORRECTIONMANAGER_H */
