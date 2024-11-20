#include "NlPropagatorManager.h"

namespace Pita { namespace Hts {

NlPropagatorManager::NlPropagatorManager(NlDynamTimeIntegrator * integrator,
                                         Seconds fineTimeStep,
                                         TimeStepCount halfSliceRatio,
                                         Seconds initialTime) :
  integrator_(integrator),
  fineTimeStep_(fineTimeStep),
  halfSliceRatio_(halfSliceRatio),
  initialTime_(initialTime),
  postProcessingMgr_(NULL),
  concurrentBasisMgr_(NULL),
  propagatedBasisMgr_(NULL)
{}

Fwk::Ptr<GenIntegratorPropagator<NlDynamTimeIntegrator> >
NlPropagatorManager::createNewInstance(const HalfSliceId & id) {
  GenIntegratorPropagator<NlDynamTimeIntegrator>::Ptr result = GenIntegratorPropagator<NlDynamTimeIntegrator>::New(integrator_.ptr());

  HalfSliceRank initialSeedRank = id.rank();
  Seconds timeStep = fineTimeStep_;

  if (id.direction() == BACKWARD) {
    initialSeedRank = initialSeedRank.next();
    timeStep = -timeStep;
  }

  Seconds sliceLength = fineTimeStep_ * halfSliceRatio_.value();
  Seconds sliceInitialTime = initialTime_ + sliceLength * initialSeedRank.value();
  
  result->initialTimeIs(sliceInitialTime);
  result->timeStepCountIs(halfSliceRatio_);
  result->timeStepSizeIs(timeStep);

  // Attach post-processor
  if (postProcessingMgr_) {
    postProcessingMgr_->outputFileSetIs(result.ptr(), PostProcessor::FileSetId(id.rank().value()));
  }

  // Prepare concurrent propagation of projection basis
  if (propagatedBasisMgr_) {
    // Request projection basis
    DynamStatePlainBasis::Ptr concurrentBasis = propagatedBasisMgr_->instance(id);
    if (!concurrentBasis) {
      concurrentBasis = propagatedBasisMgr_->instanceNew(id);
    }

    // Attach concurrent basis propagator
    if (concurrentBasisMgr_) {
      concurrentBasisMgr_->referencePropagatorIs(concurrentBasis.ptr(), result.ptr());
    }
  }

  return result;
}

} /* end namespace Hts */ } /* end namespace Pita */
