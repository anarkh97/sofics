#include "InitialSeedTask.h"

namespace Pita { 

Fwk::String
InitialSeedTask::buildName(SliceRank seedRank) {
  return String("Initial Seed ") + toString(seedRank);
}

InitialSeedTask::InitialSeedTask(Seed * targetSeed,
                                 SeedInitializer * initializer,
                                 SliceRank seedRank,
                                 Seed::Status status) :
  NamedTask(buildName(seedRank)),
  initializer_(initializer),
  seedRank_(seedRank),
  targetSeed_(targetSeed),
  status_(status)
{}

void
InitialSeedTask::iterationIs(IterationRank i) {
  targetSeed_->stateIs(initializer_->initialSeed(seedRank_));
  targetSeed_->iterationIs(i);
  targetSeed_->statusIs(status_);
}

} /* end namespace Pita */
