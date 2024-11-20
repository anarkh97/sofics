#ifndef PITA_BASICPROPAGATION_H
#define PITA_BASICPROPAGATION_H

#include "Fwk.h"
#include "Types.h"

#include "NamedTask.h"

#include "Seed.h"
#include "DynamPropagator.h" 

namespace Pita {

class BasicPropagation : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(BasicPropagation);

  const Seed * seed() const { return seed_.ptr(); }
  Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }

  void seedIs(const Seed * s);
  void propagatedSeedIs(Seed * ps);

  DynamPropagator * propagator() const { return propagator_.ptr(); }
  
  virtual void iterationIs(IterationRank i); // Overriden

  static Ptr New(const Fwk::String & name, DynamPropagator * propagator) {
    return new BasicPropagation(name, propagator);
  }

protected:
  BasicPropagation(const Fwk::String & name, DynamPropagator * propagator);

private:
  Seed::PtrConst seed_;
  Seed::Ptr propagatedSeed_;

  DynamPropagator::Ptr propagator_;

  DISALLOW_COPY_AND_ASSIGN(BasicPropagation);
};

} /* end namespace Pita */

#endif /* PITA_BASICPROPAGATION_H */
