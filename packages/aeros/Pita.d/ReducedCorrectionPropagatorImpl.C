#include "ReducedCorrectionPropagatorImpl.h"

#include <Math.d/FullSquareMatrix.h>

#include <cassert>

namespace Pita {

// ReducedCorrectionPropagatorImpl implementation

ReducedCorrectionPropagatorImpl::ReducedCorrectionPropagatorImpl(const String & name,
                                                                 const FullSquareMatrix * reprojectionMatrix,
                                                                 const RankDeficientSolver * projectionSolver):
  CorrectionPropagator<Vector>(name),
  reprojectionMatrix_(reprojectionMatrix),
  projectionSolver_(projectionSolver)
{}

void
ReducedCorrectionPropagatorImpl::iterationIs(IterationRank ir) {
  assert(correction()->iteration() == jump()->iteration() || correction()->status() == Seed::INACTIVE);
  assert(jump()->status() == Seed::CONVERGED || correction()->status() != Seed::INACTIVE);
  assert(jump()->iteration().next() == ir);

  Vector result = jump()->state();
  
  // result == reprojectionMatrix * correction + jump
  if (correction()->status() != Seed::INACTIVE) {
    const_cast<FullSquareMatrix *>(reprojectionMatrix())->multiply(
        const_cast<Vector &>(correction()->state()),
        result,
        1.0,
        FullSquareMatrix::TRANSPOSED);
  }
  
  // result == normalMatrix^{-1} * (reprojectionMatrix * correction + jump)
  projectionSolver()->solution(result);

  nextCorrection()->stateIs(result);
  nextCorrection()->statusIs(correction()->status() == Seed::ACTIVE ? Seed::ACTIVE : jump()->status());
  nextCorrection()->iterationIs(jump()->iteration());

  setIteration(ir);
}

// ReducedCorrectionPropagatorImpl::Manager implementation

ReducedCorrectionPropagatorImpl::Manager::Manager(const FullSquareMatrix * defaultReprojectionMatrix,
                                           const RankDeficientSolver * defaultProjectionSolver) :
  defaultReprojectionMatrix_(defaultReprojectionMatrix),
  defaultProjectionSolver_(defaultProjectionSolver)
{}

ReducedCorrectionPropagatorImpl *
ReducedCorrectionPropagatorImpl::Manager::createNewInstance(const String & r) {
  String instanceName(String("Propagate Reduced Correction ") + r);
  return new ReducedCorrectionPropagatorImpl(instanceName, defaultReprojectionMatrix(), defaultProjectionSolver());
}

} /* end namespace Pita */
