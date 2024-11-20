#include "JumpProjection.h"

#include <cassert>

namespace Pita {

JumpProjection::JumpProjection(const String & name,
                               const JumpProjection::Manager * manager) :
  NamedTask(name),
  manager_(manager)
{}

void
JumpProjection::iterationIs(IterationRank ir) {
  assert(seedJump()->iteration().next() == ir); // Happens at the next iteration

  Vector projectionResult(reducedBasisSize());
  for (size_t index = 0; index < reducedBasisSize(); ++index) {
    projectionResult[index] = seedJump()->state() * reducedBasis()->state(index);
  }

  reducedSeedJump()->stateIs(projectionResult);
  reducedSeedJump()->iterationIs(seedJump()->iteration());
  reducedSeedJump()->statusIs(seedJump()->status());

  setIteration(ir);
}

JumpProjection::Manager::Manager(const DynamStateBasis * drb) :
  defaultReducedBasis_(drb)
{}

JumpProjection *
JumpProjection::Manager::createNewInstance(const String & key) {
  String instanceName = String("Project Jump ") + key;
  return new JumpProjection(instanceName, this);
}

} // end namespace Pita
