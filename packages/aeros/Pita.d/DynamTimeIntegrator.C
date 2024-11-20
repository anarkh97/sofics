#include "DynamTimeIntegrator.h"

namespace Pita {

DynamTimeIntegrator::DynamTimeIntegrator(size_t vectorSize) :
  vectorSize_(vectorSize),
  initialState_(vectorSize),
  currentState_(vectorSize),
  initialTime_(0.0), 
  currentTime_(0.0),
  timeStepSize_(0),
  timeStepCount_(0)
{}
  
void
DynamTimeIntegrator::initialConditionIs(const DynamState & initialState, Seconds t0) {
  setInitialTime(t0);
  setInitialState(initialState);
  setTimeStepCount(TimeStepCount(0));
  setCurrentTime(t0);
  setCurrentState(initialState);
}

} // end namespace Pita
