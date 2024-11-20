#include "DynamStatePlainBasis.h"

namespace Pita {

DynamState &
DynamStatePlainBasis::internalState(size_t index) {
  DynamState & result = state_[index];
  unshare(result);
  return result;
}

void
DynamStatePlainBasis::firstStateIs(const DynamState & ds) {
  if (ds.vectorSize() != this->vectorSize())
    throw Fwk::RangeException(); 
  prependState(ds);
}

void
DynamStatePlainBasis::firstStateBasisIs(const DynamStateBasis * dsb) {
  if (dsb->vectorSize() != this->vectorSize())
    throw Fwk::RangeException();
  size_t statesToAdd = dsb->stateCount();
  for (int i = statesToAdd - 1; i >= 0; --i) {
    prependState(dsb->state(i));
  }
}

void
DynamStatePlainBasis::firstStateBasisIs(const DynamStatePlainBasis * dsb) {
  if (dsb->vectorSize() != this->vectorSize())
    throw Fwk::RangeException();
  prependStateBasis(dsb);
}

void
DynamStatePlainBasis::lastStateIs(const DynamState & ds) {
  if (ds.vectorSize() != this->vectorSize())
    throw Fwk::RangeException(); 
  appendState(ds);
}

void
DynamStatePlainBasis::lastStateBasisIs(const DynamStateBasis * dsb) {
  if (dsb->vectorSize() != this->vectorSize())
    throw Fwk::RangeException();
  size_t statesToAdd = dsb->stateCount();
  for (size_t i = 0; i < statesToAdd; ++i) {
    appendState(dsb->state(i));
  }
}

void
DynamStatePlainBasis::lastStateBasisIs(const DynamStatePlainBasis * dsb) {
  if (dsb->vectorSize() != this->vectorSize())
    throw Fwk::RangeException();
  appendStateBasis(dsb);
}


} // end namespace Pita
