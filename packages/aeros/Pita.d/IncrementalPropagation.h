#ifndef PITA_INCREMENTALPROPAGATION_H
#define PITA_INCREMENTALPROPAGATION_H

#include "Fwk.h"
#include "Types.h"

#include "NamedTask.h"

#include "Seed.h"
#include "AffineDynamPropagator.h" 

namespace Pita {

class IncrementalPropagation : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(IncrementalPropagation);

  const Seed * seed() const { return seed_.ptr(); }
  Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }

  void seedIs(const Seed * s);
  void propagatedSeedIs(Seed * ps);

  AffineDynamPropagator * propagator() const { return propagator_.ptr(); }
  
  virtual void iterationIs(IterationRank i); // Overriden

  static Ptr New(const Fwk::String & name, AffineDynamPropagator * propagator) {
    return new IncrementalPropagation(name, propagator);
  }

protected:
  IncrementalPropagation(const Fwk::String & name, AffineDynamPropagator * propagator);

private:
  Seed::PtrConst seed_;
  Seed::Ptr propagatedSeed_;

  AffineDynamPropagator::Ptr propagator_;
  DynamState previousSeedState_;

  DISALLOW_COPY_AND_ASSIGN(IncrementalPropagation);
};

} /* end namespace Pita */

#endif /* PITA_INCREMENTALPROPAGATION_H */
