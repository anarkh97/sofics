#ifndef PITA_INTEGRATORSEEDINITIALIZER_H
#define PITA_INTEGRATORSEEDINITIALIZER_H

#include "Fwk.h"
#include "Types.h"

#include "SeedInitializer.h"

#include "DynamTimeIntegrator.h"
#include "DynamState.h"
#include "DynamStatePlainBasis.h"

namespace Pita {

class IntegratorSeedInitializer : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(IntegratorSeedInitializer);

  virtual DynamState initialSeed(SliceRank rank) const; // Overriden
  
  const DynamTimeIntegrator * integrator() const { return integrator_.ptr(); }
  Seconds timeStepSize() const { return timeStepSize_; }
  TimeStepCount timeStepsBetweenSeeds() const { return timeStepsBetweenSeeds_; }

  // integrator must have initialCondition and timeStepSize correctly initialized must not be modified externally
  static Ptr New(DynamTimeIntegrator * integrator,
                 TimeStepCount timeStepsBetweenSeeds) {
    return new IntegratorSeedInitializer(integrator, integrator->timeStepSize(), timeStepsBetweenSeeds, integrator->currentState(), integrator->currentTime());
  }


protected:
  IntegratorSeedInitializer(DynamTimeIntegrator * i,
                            Seconds tss,
                            TimeStepCount tsbs,
                            DynamState is,
                            Seconds it);

private:
  mutable DynamTimeIntegrator::Ptr integrator_;

  Seconds timeStepSize_;
  TimeStepCount timeStepsBetweenSeeds_;

  mutable Seconds lastStateTime_;
  mutable DynamStatePlainBasis::Ptr seed_;
};

} // end namespace Pita

#endif /* PITA_INTEGRATORSEEDINITIALIZER_H */
