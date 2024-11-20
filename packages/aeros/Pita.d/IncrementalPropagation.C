#include "IncrementalPropagation.h"

namespace Pita {

IncrementalPropagation::IncrementalPropagation(const Fwk::String & name,
                                               AffineDynamPropagator * propagator) :
  NamedTask(name),
  propagator_(propagator),
  previousSeedState_()
{}

void
IncrementalPropagation::seedIs(const Seed * s) {
  previousSeedState_ = (s->status() == Seed::ACTIVE) ? s->state() : DynamState();
  seed_ = s;
}

void
IncrementalPropagation::propagatedSeedIs(Seed * ps) {
  propagatedSeed_ = ps;
}

void
IncrementalPropagation::iterationIs(IterationRank i) {
  // Sanity checks
  assert(seed()->iteration() == i);
  assert(seed()->status() != Seed::INACTIVE);

  // Perform propagation
  if (previousSeedState_.vectorSize() != 0) {
    propagator()->constantTermIs(AffineDynamPropagator::HOMOGENEOUS);
    propagator()->initialStateIs(seed()->state() - previousSeedState_);
  } else {
    propagator()->initialStateIs(seed()->state());
  }
 
  propagatedSeed()->statusIs(seed()->status());
  if (propagatedSeed()->state().vectorSize() != 0) {
    propagatedSeed()->stateIs(propagatedSeed()->state() + propagator()->finalState());
  } else {
    propagatedSeed()->stateIs(propagator()->finalState());
  }

  previousSeedState_ = seed()->state();

  propagatedSeed()->iterationIs(seed()->iteration());
}

} /* end namespace Pita */
