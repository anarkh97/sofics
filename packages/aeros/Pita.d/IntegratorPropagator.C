#include "IntegratorPropagator.h"
#include "DynamTimeIntegrator.h"

namespace Pita {
  
IntegratorPropagatorRoot::IntegratorPropagatorRoot(DynamTimeIntegrator * integrator) :
  DynamPropagator(integrator ? integrator->vectorSize() : 0),
  integrator_(integrator),
  initialTime_(integrator ? integrator->initialTime() : Seconds(0.0)),
  timeStepCount_(TimeStepCount(1)),
  timeStepSize_(integrator ? integrator->timeStepSize() : Seconds(0.0))
{}

void
IntegratorPropagatorRoot::initialStateIs(const DynamState & initialState) {
  setInitialState(initialState);
  initialStateNotify();

  if (integrator()) {
    integrator()->timeStepSizeIs(timeStepSize());
    integrator()->initialConditionIs(initialState, initialTime());
    integrator()->timeStepCountInc(timeStepCount());
    setFinalState(integrator()->currentState());
  } else {
    setFinalState(initialState);
  }

  finalStateNotify();
}

// Obsolete
IntegratorPropagator::IntegratorPropagator(DynamTimeIntegrator * integrator) :
  DynamPropagator(integrator ? integrator->vectorSize() : 0),
  integrator_(integrator),
  initialTime_(integrator ? integrator->initialTime() : Seconds(0.0)),
  timeStepCount_(TimeStepCount(1)),
  timeStepSize_(integrator ? integrator->timeStepSize() : Seconds(0.0))
{}

void
IntegratorPropagator::initialStateIs(const DynamState & initialState) {
  setInitialState(initialState);
  initialStateNotify();

  if (integrator()) {
    integrator()->timeStepSizeIs(timeStepSize());
    integrator()->initialConditionIs(initialState, initialTime());
    integrator()->timeStepCountInc(timeStepCount());
    setFinalState(integrator()->currentState());
  } else {
    setFinalState(initialState);
  }

  finalStateNotify();
}

} // end namespace Pita
