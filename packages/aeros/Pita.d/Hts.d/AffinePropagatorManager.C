#include "AffinePropagatorManager.h"

#include "../AffineIntegratorPropagator.h"

#include <cassert>

namespace Pita { namespace Hts {

AffinePropagatorManager::AffinePropagatorManager(AffineBasisCollector * collector,
                                                 GenFineIntegratorManager<AffineGenAlphaIntegrator> * integratorMgr,
                                                 PostProcessing::Manager * postProcessingMgr,
                                                 TimeStepCount halfSliceRatio,
                                                 Seconds initialTime,
                                                 AffineDynamPropagator::ConstantTerm defaultMode) :
  collector_(collector),
  integratorMgr_(integratorMgr),
  postProcessingMgr_(postProcessingMgr),
  fineTimeStep_(integratorMgr->fineTimeStepSize()),
  halfSliceRatio_(halfSliceRatio),
  initialTime_(initialTime),
  defaultMode_(defaultMode)
{}

AffineDynamPropagator *
AffinePropagatorManager::instance(const HalfSliceId & id) const {
  return const_cast<AffineDynamPropagator *>(collector_->source(id));
}

size_t
AffinePropagatorManager::instanceCount() const { 
  return collector_->sourceCount();
}

AffineDynamPropagator *
AffinePropagatorManager::instanceNew(const HalfSliceId & id) {
  // Create propagator 
  AffineGenAlphaIntegrator::Ptr integrator = integratorMgr_->fineIntegrator(id.direction());
  AffineIntegratorPropagator::Ptr newPropagator = AffineIntegratorPropagator::New(integrator.ptr());

  // Set up propagator
  Seconds halfCoarseTimeStep = fineTimeStep_ * halfSliceRatio_.value();
  HalfSliceRank initialSeedRank = (id.direction() == FORWARD) ? id.rank() : id.rank() + HalfSliceCount(1);
  Seconds sliceInitialTime = initialTime_ + halfCoarseTimeStep * initialSeedRank.value();

  newPropagator->initialTimeIs(sliceInitialTime);
  newPropagator->timeStepCountIs(halfSliceRatio_);
  newPropagator->constantTermIs(defaultMode_);
  
  // Attach AffineBasisCollector Reactor 
  collector_->sourceIs(id, newPropagator.ptr());
 
  // Attach PostProcessing Reactor
  if (postProcessingMgr_) {
    this->postProcessingMgr_->outputFileSetIs(newPropagator.ptr(), PostProcessor::FileSetId(id.rank().value()));
  }

  return newPropagator.ptr();
}

void
AffinePropagatorManager::instanceDel(const HalfSliceId & id) {
  collector_->sourceIs(id, NULL);
}

} /* end namespace Hts */ } /* end namespace Pita */
