#include "JumpBuilder.h"

namespace Pita {

void
JumpBuilder::iterationIs(IterationRank ir) {
  assert(actualSeed()->iteration() == ir);
  assert(actualSeed()->status() != Seed::INACTIVE);
  assert(actualSeed()->status() != Seed::SPECIAL);
  assert(predictedSeed()->iteration() == ir);

  DynamState newJumpSeed = actualSeed()->state() - predictedSeed()->state();
  seedJump()->statusIs(actualSeed()->status());
  seedJump()->stateIs(newJumpSeed);
  seedJump()->iterationIs(actualSeed()->iteration());

  setIteration(ir);
} 

JumpBuilder *
JumpBuilder::ManagerImpl::createNewInstance(const String & key) {
  String taskName = "Evaluate Jump " + toString(key);
  return new JumpBuilder(taskName);
}

} /* end namespace */
