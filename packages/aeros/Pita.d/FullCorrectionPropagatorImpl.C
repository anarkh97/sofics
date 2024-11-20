#include "FullCorrectionPropagatorImpl.h"

#include <cassert>

namespace Pita {

/* FullCorrectionPropagatorImpl implementation */
  
FullCorrectionPropagatorImpl::FullCorrectionPropagatorImpl(const String & name, DynamPropagator * propagator) :
  CorrectionPropagator<DynamState>(name),
  propagator_(propagator)
{}

void
FullCorrectionPropagatorImpl::iterationIs(IterationRank ir) {
  assert(correction()->iteration() == jump()->iteration() || correction()->status() == Seed::INACTIVE);
  assert(jump()->status() == Seed::CONVERGED || correction()->status() != Seed::INACTIVE);
  assert(jump()->iteration().next() == ir);

  DynamState initialState = jump()->state();

  if (correction()->status() != Seed::INACTIVE) {
    initialState += correction()->state();
  }

  propagator_->initialStateIs(initialState);
  
  nextCorrection()->statusIs(correction()->status() == Seed::ACTIVE ? Seed::ACTIVE : jump()->status());
  nextCorrection()->stateIs(propagator_->finalState());
  nextCorrection()->iterationIs(jump()->iteration());
 
  setIteration(jump()->iteration());
}

/* Manager implementation */

FullCorrectionPropagatorImpl::Manager::Manager(DynamPropagator * sharedPropagator) :
  sharedPropagator_(sharedPropagator)
{}

FullCorrectionPropagatorImpl *
FullCorrectionPropagatorImpl::Manager::createNewInstance(const String & key) {
  return new FullCorrectionPropagatorImpl(String("Propagate Correction ") + key, sharedPropagator()); 
}

} /* end namespace Pita */
