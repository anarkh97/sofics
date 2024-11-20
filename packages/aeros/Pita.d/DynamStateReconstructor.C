#include "DynamStateReconstructor.h"

namespace Pita {

DynamStateReconstructor::DynamStateReconstructor(const DynamStateBasis * rb) :
  reconstructionBasis_(rb)
{}

DynamState
DynamStateReconstructor::result(const Vector & c) const {
  int rbs = static_cast<int>(reducedBasisSize());
  if (c.size() != rbs) {
    throw Fwk::RangeException("Dimension mismatch");
  }

  DynamState answer(vectorSize(), 0.0);
  for (int i = 0; i < rbs; ++i) {
    if (c[i] != 0.0) {
      answer.linAdd(c[i], reconstructionBasis_->state(i));
    }
  }
  return answer;
}

} // end namespace Pita
