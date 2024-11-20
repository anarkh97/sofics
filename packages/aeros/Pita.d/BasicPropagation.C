#include "BasicPropagation.h"

namespace Pita {

BasicPropagation::BasicPropagation(const Fwk::String & name,
                                   DynamPropagator * propagator) :
  NamedTask(name),
  propagator_(propagator)
{}

void
BasicPropagation::seedIs(const Seed * s) {
  seed_ = s;
}

void
BasicPropagation::propagatedSeedIs(Seed * ps) {
  propagatedSeed_ = ps;
}

void
BasicPropagation::iterationIs(IterationRank i) {
  // Sanity checks
  assert(seed()->iteration() == i);
  assert(seed()->status() != Seed::INACTIVE);

  propagator()->initialStateIs(seed()->state());
 
  propagatedSeed()->statusIs(seed()->status());
  propagatedSeed()->stateIs(propagator()->finalState());
  propagatedSeed()->iterationIs(seed()->iteration());
}

} /* end namespace Pita */
