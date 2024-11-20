#include "DynamStateOps.h"

#include <Math.d/SparseMatrix.h>

namespace Pita {

DynamState
mult(const DynamOps * metric, const DynamState & state) {
  DynamState temp(state.vectorSize());
  const_cast<SparseMatrix*>(metric->stiffnessMatrix())->mult(state.displacement(), temp.displacement());
  const_cast<SparseMatrix*>(metric->massMatrix())->mult(state.velocity(), temp.velocity());

  return temp;
}

void
mult(const DynamOps * metric, const DynamState & state, double * resultBuffer) {
  const_cast<SparseMatrix*>(metric->stiffnessMatrix())->mult(state.displacement(), resultBuffer);
  const_cast<SparseMatrix*>(metric->massMatrix())->mult(state.velocity(), resultBuffer + state.vectorSize());
}

} /* end namespace Pita */
