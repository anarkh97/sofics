#include "IntegratorSeedInitializer.h"

#include <cassert>

namespace Pita {

IntegratorSeedInitializer::IntegratorSeedInitializer(DynamTimeIntegrator * i,
                                                     Seconds tss,
                                                     TimeStepCount tsbs,
                                                     DynamState is,
                                                     Seconds it) :
  SeedInitializer(i->vectorSize(), SliceRank(0)),
  integrator_(i),
  timeStepSize_(tss),
  timeStepsBetweenSeeds_(tsbs),
  lastStateTime_(it),
  seed_(DynamStatePlainBasis::New(i->vectorSize()))
{
  if (is.vectorSize() == 0) {
    is = i->currentState();
  }
  assert(is.vectorSize() == i->vectorSize());
  seed_->lastStateIs(is);
}

DynamState
IntegratorSeedInitializer::initialSeed(SliceRank rank) const {
  if (rank > lastSlice()) { // State is not in cache
    // Resetting the initial condition and/or the time-step would be too inefficient:
    // Therefore, it is assumed that the integrator has not been changed externally

    // Fails if the integrator has obviously been modified externally
    assert(integrator_->timeStepSize() == timeStepSize_);
    assert(integrator_->currentTime() == lastStateTime_);

    // Compute and put in cache target value and intermediate results
    int targetIndex = rank.value();
    for (int index = lastSlice().value(); index < targetIndex; ++index) {
      integrator_->timeStepCountInc(timeStepsBetweenSeeds());
      seed_->lastStateIs(integrator_->currentState()); 
    }

    // Preserve time
    lastStateTime_ = integrator_->currentTime();

    // Update number of available slices
    const_cast<IntegratorSeedInitializer *>(this)->setSlices(rank);
  }

  return seed_->state(rank.value());
}

} // end namespace Pita
