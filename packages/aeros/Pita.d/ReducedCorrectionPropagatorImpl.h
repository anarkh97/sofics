#ifndef PITA_REDUCEDCORRECTIONPROPAGATORIMPL_H
#define PITA_REDUCEDCORRECTIONPROPAGATORIMPL_H

#include "CorrectionPropagator.h"
#include "Seed.h"

#include "RankDeficientSolver.h"

template <typename Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

namespace Pita {

class ReducedCorrectionPropagatorImpl : public CorrectionPropagator<Vector> {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedCorrectionPropagatorImpl);
  class Manager;

  // Overriden mutators
  virtual void iterationIs(IterationRank i);

  // Reduced-space operators
  const FullSquareMatrix * reprojectionMatrix() const { return reprojectionMatrix_; }
  const RankDeficientSolver * projectionSolver() const { return projectionSolver_; }

  void reprojectionMatrixIs(const FullSquareMatrix * rm) { reprojectionMatrix_ = rm; }
  void projectionSolverIs(const RankDeficientSolver * ps) { projectionSolver_ = ps; }

protected:
  ReducedCorrectionPropagatorImpl(const String & name,
                           const FullSquareMatrix * reprojectionMatrix,
                           const RankDeficientSolver * projectionSolver);
  
private:
  const FullSquareMatrix * reprojectionMatrix_;
  const RankDeficientSolver * projectionSolver_;
};


class ReducedCorrectionPropagatorImpl::Manager : public CorrectionPropagator<Vector>::Manager, private Fwk::GenManagerImpl<ReducedCorrectionPropagatorImpl, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  ReducedCorrectionPropagatorImpl * instance(const String & r) const { return Fwk::GenManagerImpl<ReducedCorrectionPropagatorImpl, String>::instance(r); }
  InstanceCount instanceCount() const { return Fwk::GenManagerImpl<ReducedCorrectionPropagatorImpl, String>::instanceCount(); }

  ReducedCorrectionPropagatorImpl * instanceNew(const String & r) { return Fwk::GenManagerImpl<ReducedCorrectionPropagatorImpl, String>::instanceNew(r); } 
  void instanceDel(const String & r) { Fwk::GenManagerImpl<ReducedCorrectionPropagatorImpl, String>::instanceDel(r); }

  const FullSquareMatrix * defaultReprojectionMatrix() const { return defaultReprojectionMatrix_; }
  const RankDeficientSolver * defaultProjectionSolver() const { return defaultProjectionSolver_.ptr(); }
  
  void defaultReprojectionMatrixIs(const FullSquareMatrix * drm) { defaultReprojectionMatrix_ = drm; }
  void defaultPojectionSolverIs(const RankDeficientSolver * ds) { defaultProjectionSolver_ = ds; }

  static Ptr New(const FullSquareMatrix * defaultReprojectionMatrix, const RankDeficientSolver * defaultProjectionSolver) {
    return new Manager(defaultReprojectionMatrix, defaultProjectionSolver);
  }
  
protected:
  Manager(const FullSquareMatrix * defaultReprojectionMatrix, const RankDeficientSolver * defaultProjectionSolver);

  virtual ReducedCorrectionPropagatorImpl * createNewInstance(const String & r);

private:
  const FullSquareMatrix * defaultReprojectionMatrix_;
  RankDeficientSolver::PtrConst defaultProjectionSolver_;
};

} /* end namespace Pita */

#endif /* PITA_REDUCEDCORRECTIONPROPAGATORIMPL_H */
