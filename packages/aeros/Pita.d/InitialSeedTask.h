#ifndef PITA_INITIAL_SEEDTASK_H
#define PITA_INITIAL_SEEDTASK_H

#include "NamedTask.h"
#include "SeedInitializer.h"
#include "Seed.h"

namespace Pita { 

class InitialSeedTask : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(InitialSeedTask);

  void iterationIs(IterationRank i);

  InitialSeedTask(Seed * targetSeed, SeedInitializer * initializer,
                  SliceRank seedRank, Seed::Status status);

protected:
  static String buildName(SliceRank seedRank);

private:
  Seed::Ptr targetSeed_;
  SeedInitializer::Ptr initializer_;
  SliceRank seedRank_;
  Seed::Status status_;
};

} /* end namespace Pita */

#endif /* PITA_INITIAL_SEEDTASK_H */
