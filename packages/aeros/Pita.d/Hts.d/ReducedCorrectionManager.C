#include "ReducedCorrectionManager.h"

#include "../ReducedCorrectionPropagatorImpl.h"
#include "../JumpProjection.h"
#include "../UpdatedSeedAssemblerImpl.h"

#include "LinearProjectionNetwork.h"

namespace Pita { namespace Hts {

ReducedCorrectionManager::ReducedCorrectionManager(LinearProjectionNetwork * correctionMgr,
                                                   CorrectionPropagator<DynamState>::Manager * fcpMgr) :
  jumpProjMgr_(JumpProjection::Manager::New(correctionMgr->projectionBasis())),
  rcpMgr_(ReducedCorrectionPropagatorImpl::Manager::New(correctionMgr->reprojectionMatrix(), correctionMgr->normalMatrixSolver())),
  fcpMgr_(fcpMgr),
  usaMgr_(UpdatedSeedAssemblerImpl::Manager::New(correctionMgr->propagatedBasis()))
{}

} /* end namespace Hts */ } /* end namespace Pita */
