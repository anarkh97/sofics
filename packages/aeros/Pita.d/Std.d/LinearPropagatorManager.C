#include "LinearPropagatorManager.h"

namespace Pita { namespace Std {

LinearPropagatorManager::LinearPropagatorManager(AffineDynamTimeIntegrator * sharedIntegrator,
                                                 PostProcessing::Manager * postProcessingMgr,
                                                 AffineBasisCollector * collector,
                                                 TimeStepCount timeSliceRatio,
                                                 Seconds initialTime,
                                                 AffineDynamPropagator::ConstantTerm defaultMode) :
  sharedIntegrator_(sharedIntegrator),
  postProcessingMgr_(postProcessingMgr),
  collector_(collector),
  fineTimeStep_(sharedIntegrator->timeStepSize()),
  timeSliceRatio_(timeSliceRatio),
  initialTime_(initialTime),
  defaultMode_(defaultMode),
  instance_()
{}

AffineIntegratorPropagator *
LinearPropagatorManager::instance(const SliceRank & id) const {
  InstanceMap::const_iterator it = instance_.find(id);
  return (it != instance_.end()) ? it->second.ptr() : NULL;
}

size_t
LinearPropagatorManager::instanceCount() const {
  return instance_.size();
}

AffineIntegratorPropagator *
LinearPropagatorManager::instanceNew(const SliceRank & id) {
  // Find insertion point
  InstanceMap::iterator it = instance_.lower_bound(id);
  if (it != instance_.end() && it->first == id) {
    throw NameInUseException();
  }

  // Build and add a new instance
  AffineIntegratorPropagator::Ptr newInstance = createNewInstance(id);
  instance_.insert(it, std::make_pair(id, newInstance));

  return newInstance.ptr();
}

void
LinearPropagatorManager::instanceDel(const SliceRank & id) {
  instance_.erase(id);
}

AffineIntegratorPropagator::Ptr
LinearPropagatorManager::createNewInstance(SliceRank rank) {
  // Create propagator
  AffineIntegratorPropagator::Ptr newInstance = AffineIntegratorPropagator::New(sharedIntegrator_.ptr());

  // Set-up propagator
  Seconds sliceInitialTime = initialTime_ + fineTimeStep_ * (rank.value() * timeSliceRatio_.value());
  newInstance->initialTimeIs(sliceInitialTime);
  newInstance->timeStepCountIs(timeSliceRatio_);
  newInstance->constantTermIs(defaultMode_);

  // Data collection for correction
  if (collector_) {
    collector_->sourceIs(rank, newInstance.ptr());
  }
  
  // Set-up post-processing
  if (postProcessingMgr_) {
    postProcessingMgr_->outputFileSetIs(newInstance.ptr(), PostProcessor::FileSetId(rank.value()));
  }

  return newInstance;
}

} /* end namespace Std */ } /* end namespace Pita */
