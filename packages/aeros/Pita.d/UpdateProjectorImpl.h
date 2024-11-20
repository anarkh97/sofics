#ifndef PITA_UPDATEPROJECTORIMPL_H
#define PITA_UPDATEPROJECTORIMPL_H

#include "UpdateProjector.h"

#include "DynamStateBasis.h"
#include "RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

namespace Pita {

class UpdateProjectorImpl : public UpdateProjector {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdateProjectorImpl);
  class Manager;

  virtual size_t reducedBasisSize() const;

  /* Results */  
  using UpdateProjector::reducedJump;
  Vector & reducedJump() { return const_cast<Vector &>(const_cast<const UpdateProjectorImpl *>(this)->reducedJump()); }
 
  using UpdateProjector::updateComponents;
  Vector & updateComponents() { return const_cast<Vector &>(const_cast<const UpdateProjectorImpl *>(this)->updateComponents()); }

  /* Sources */
  virtual void jumpIs(const Seed * jb);
  virtual void reducedCorrectionIs(const ReducedSeed * rc);

  /* Reduced basis and operators */
  const DynamStateBasis * reductionBasis() const;
  const FullSquareMatrix * reprojectionMatrix() const;
  const RankDeficientSolver * solver() const;

  // HACK
  void doProjection();

protected:
  class JumpReactor;
  class CorrectionReactor;

  explicit UpdateProjectorImpl(const Manager * manager);

  friend class Manager;

private:
  const Manager * manager_;

  DynamStateBasis::PtrConst reductionBasis_;
  const FullSquareMatrix * reprojectionMatrix_;
  RankDeficientSolver::PtrConst solver_;

  Fwk::Ptr<JumpReactor> jumpReactor_;
  Fwk::Ptr<CorrectionReactor> correctionReactor_;
};

/* UpdateProjectorImpl::Manager definition */

class UpdateProjectorImpl::Manager : public UpdateProjector::Manager, private Fwk::GenManagerImpl<UpdateProjectorImpl, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  // Overriden members
  virtual UpdateProjectorImpl * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual UpdateProjectorImpl * instanceNew(const String & key) { return Impl::instanceNew(key); } 
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  // Added members
  const DynamStateBasis * defaultReductionBasis() const { return defaultReductionBasis_.ptr(); }
  const FullSquareMatrix * defaultReprojectionMatrix() const { return defaultReprojectionMatrix_; }
  const RankDeficientSolver * defaultSolver() const { return defaultSolver_.ptr(); }

  void defaultReductionBasisIs(const DynamStateBasis * drb) { defaultReductionBasis_ = drb; }
  void defaultReprojectionMatrixIs(const FullSquareMatrix * drm) { defaultReprojectionMatrix_ = drm; }
  void defaultSolverIs(const RankDeficientSolver * ds) { defaultSolver_ = ds; }

  static Ptr New(const DynamStateBasis * defaultReductionBasis,
                 const FullSquareMatrix * defaultReprojectionMatrix,
                 const RankDeficientSolver * defaultSolver) {
    return new Manager(defaultReductionBasis, defaultReprojectionMatrix, defaultSolver);
  }

protected:
  Manager(const DynamStateBasis * drb, const FullSquareMatrix * drm, const RankDeficientSolver * ds);

  virtual UpdateProjectorImpl * createNewInstance(const String & key);

private:
  typedef Fwk::GenManagerImpl<UpdateProjectorImpl, String> Impl;
  
  DynamStateBasis::PtrConst defaultReductionBasis_;
  const FullSquareMatrix * defaultReprojectionMatrix_;
  RankDeficientSolver::PtrConst defaultSolver_;
};


inline
const DynamStateBasis *
UpdateProjectorImpl::reductionBasis() const {
  return manager_->defaultReductionBasis();
}

inline
const FullSquareMatrix *
UpdateProjectorImpl::reprojectionMatrix() const {
  return manager_->defaultReprojectionMatrix();
}

inline
const RankDeficientSolver *
UpdateProjectorImpl::solver() const {
  return manager_->defaultSolver();
}

} /* end namespace Pita */

#endif /* PITA_UPDATEPROJECTORIMPL_H */
