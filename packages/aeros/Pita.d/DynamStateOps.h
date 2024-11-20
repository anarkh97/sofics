#ifndef PITA_DYNAMSTATEOPS_H
#define PITA_DYNAMSTATEOPS_H

#include "DynamState.h"
#include "DynamOps.h"

namespace Pita {

// Computes [metric->K * state.displacement, metric->M * state.velocity]
DynamState mult(const DynamOps * metric, const DynamState & state);

// Put [metric->K * state.displacement, metric->M * state.velocity] in resultBuffer 
// (resultBuffer must point to a buffer of size at least 2 * state.vectorSize())
void mult(const DynamOps * metric, const DynamState & state, double * resultBuffer);

inline
double
dotProduct(const DynamOps * metric, const DynamState & state1, const DynamState & state2) {
  return mult(metric, state1) * state2;
}

inline
double
energy(const DynamOps * metric, const DynamState & state) {
  return 0.5 * dotProduct(metric, state, state);
}

} /* end namespace Pita */

#endif /* PITA_DYNAMSTATEOPS_H */
