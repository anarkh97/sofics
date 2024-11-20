#ifndef PITA_LINEARIZEDPROPAGATOR_H
#define PITA_LINEARIZEDPROPAGATOR_H

#include "Fwk.h"
#include "DynamState.h"

namespace Pita {

class PitaNonLinDynamic;
class NlDynamTimeIntegrator;

class LinearizedPropagator : public Fwk::PtrInterface<LinearizedPropagator> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearizedPropagator);
  
  size_t vectorSize() const;
  const PitaNonLinDynamic * probDesc() const;

  DynamState & finalState(DynamState & modifiedInitialState) const;

  static LinearizedPropagator::Ptr New(const NlDynamTimeIntegrator * baseIntegrator) {
    return new LinearizedPropagator(baseIntegrator);
  }

protected:
  explicit LinearizedPropagator(const NlDynamTimeIntegrator * baseIntegrator);

private:
  Fwk::Ptr<const NlDynamTimeIntegrator> baseIntegrator_;

  // Internal implementation
  mutable GenVector<double> temp1_;
  mutable GenVector<double> temp2_;
};  
  
} /* end namespace Pita */

#endif /* PITA_LINEARIZEDPROPAGATOR_H */
