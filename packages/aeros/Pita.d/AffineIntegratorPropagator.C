#include "AffineIntegratorPropagator.h"

namespace Pita {
  
AffineIntegratorPropagator::AffineIntegratorPropagator(AffineDynamTimeIntegrator * integrator) :
  AffineDynamPropagator(integrator ? integrator->vectorSize() : 0, NONHOMOGENEOUS),
  integrator_(integrator),
  initialTime_(Seconds(0.0)),
  timeStepCount_(TimeStepCount(1))
{}

void
AffineIntegratorPropagator::initialStateIs(const DynamState & initialState) {
  setInitialState(initialState);
  initialStateNotify();

  if (integrator()) {
    integrator()->initialConditionIs(initialState, initialTime());
    integrator()->timeStepCountInc(this->timeStepCount());
    setFinalState(integrator()->currentState());
  } else {
    setFinalState(initialState);
  }

  finalStateNotify();
}

void
AffineIntegratorPropagator::constantTermIs(AffineIntegratorPropagator::ConstantTerm ct) {
  if (ct == HOMOGENEOUS) {
    integrator()->externalForceStatusIs(AffineDynamTimeIntegrator::HOMOGENEOUS);
  } else {
    integrator()->externalForceStatusIs(AffineDynamTimeIntegrator::NONHOMOGENEOUS);
  }

  setConstantTerm(ct);
}

} // end namespace Pita
