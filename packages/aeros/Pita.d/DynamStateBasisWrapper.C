#include "DynamStateBasisWrapper.h"

namespace Pita {

DynamStateBasisWrapper::DynamStateBasisWrapper(size_t vectorSize, size_t stateCount, double * data) :
  DynamStateBasis(vectorSize),
  stateCount_(stateCount),
  data_(data)
{}

  
} // end namespace Pita
