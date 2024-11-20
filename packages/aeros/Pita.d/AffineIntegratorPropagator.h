#ifndef PITA_AFFINEINTEGRATORPROPAGATOR_H
#define PITA_AFFINEINTEGRATORPROPAGATOR_H

#include "AffineDynamPropagator.h"
#include "AffineDynamTimeIntegrator.h"

namespace Pita {

// TODO: Allow time-integrator sharing by adding attribute (Seconds timeStepSize)
class AffineIntegratorPropagator : public AffineDynamPropagator {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineIntegratorPropagator);

  AffineDynamTimeIntegrator * integrator() const { return integrator_.ptr(); }
  Seconds initialTime() const { return initialTime_; }
  TimeStepCount timeStepCount() const { return timeStepCount_; }
  
  virtual void initialStateIs(const DynamState & initialState);
  virtual void constantTermIs(ConstantTerm ct);
  void integratorIs(AffineDynamTimeIntegrator * integrator) { integrator_ = integrator; };
  void initialTimeIs(Seconds t0) { initialTime_ = t0; } 
  void timeStepCountIs(TimeStepCount timeStepCount) { timeStepCount_ = timeStepCount; }

  static AffineIntegratorPropagator::Ptr New(AffineDynamTimeIntegrator * integrator) {
    return new AffineIntegratorPropagator(integrator);
  }
  
protected:
  explicit AffineIntegratorPropagator(AffineDynamTimeIntegrator * integrator);

private:
  AffineDynamTimeIntegrator::Ptr integrator_;
  Seconds initialTime_;
  TimeStepCount timeStepCount_;
  SliceRank sliceRank_;
};

}

#endif /* PITA_AFFINEINTEGRATORPROPAGATOR_H */
