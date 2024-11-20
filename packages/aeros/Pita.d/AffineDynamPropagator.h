#ifndef PITA_AFFINEDYNAMPROPAGATOR_H
#define PITA_AFFINEDYNAMPROPAGATOR_H

#include "DynamPropagator.h"

namespace Pita {

class AffineDynamPropagator : public DynamPropagator {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineDynamPropagator);

  enum ConstantTerm {
    NONHOMOGENEOUS = 0,
    HOMOGENEOUS
  };

  ConstantTerm constantTerm() const { return constantTerm_; }
  virtual void constantTermIs(ConstantTerm ct) = 0;

protected:
  explicit AffineDynamPropagator(size_t vectorSize, ConstantTerm constantTerm) :
    DynamPropagator(vectorSize),
    constantTerm_(constantTerm)
  {}

  void setConstantTerm(ConstantTerm at) { constantTerm_ = at; }

private:
  ConstantTerm constantTerm_;
};

} /* end namespace Pita */

#endif /* PITA_AFFINEDYNAMPROPAGATOR_H */
