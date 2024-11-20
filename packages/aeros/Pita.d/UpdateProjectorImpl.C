#include "UpdateProjectorImpl.h"

namespace Pita {

/* UpdateProjectorImpl::JumpReactor definition */

class UpdateProjectorImpl::JumpReactor : public Seed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpReactor);

  virtual void onState(); // overriden

  JumpReactor(const Seed * notifier, UpdateProjectorImpl * parent);

private:
  UpdateProjectorImpl * parent_;
}; 

/* UpdateProjectorImpl::JumpReactor implementation */

UpdateProjectorImpl::JumpReactor::JumpReactor(const Seed * notifier, UpdateProjectorImpl * parent) :
  Seed::NotifieeConst(notifier),
  parent_(parent)
{}

void
UpdateProjectorImpl::JumpReactor::onState() {
  // Aliases
  const DynamState & jump = notifier()->state();
  const DynamStateBasis * basis = parent_->reductionBasis();
  const RankDeficientSolver * solver = parent_->solver();
  
  // reducedJump <- 0
  Vector & reducedJump = parent_->reducedJump();
  reducedJump.reset(basis->stateCount());
 
  // Fill-in relevant part of reducedJump 
  int factorRank = solver->factorRank(); 
  for (int i = 0; i < factorRank; ++i) {
    int index = solver->factorPermutation(i);
    reducedJump[index] = jump * basis->state(index);
  }
  
  log() << "Reduced jump of size " << reducedJump.size() << "\n";
}

/* UpdateProjectorImpl::CorrectionReactor definition */

class UpdateProjectorImpl::CorrectionReactor : public ReducedSeed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReactor);

  virtual void onState(); // overriden

  CorrectionReactor(const ReducedSeed * notifier, UpdateProjectorImpl * parent);

private:
  UpdateProjectorImpl * parent_;
};

/* UpdateProjectorImpl::CorrectionReactor implementation */

UpdateProjectorImpl::CorrectionReactor::CorrectionReactor(const ReducedSeed * notifier, UpdateProjectorImpl * parent) :
  ReducedSeed::NotifieeConst(notifier),
  parent_(parent)
{}

void
UpdateProjectorImpl::CorrectionReactor::onState() {
  // Aliases
  const Vector & reducedCorrection = notifier()->state();
  Vector & updateComponents = parent_->updateComponents();
  
  // updateComponents <- reducedJump
  updateComponents = parent_->reducedJump();

  log() << "updateComponents.size = " << updateComponents.size() << "\n";
  
  // updateComponents <- reprojectionMatrix * reducedCorrection + updateComponents
  const_cast<FullSquareMatrix *>(parent_->reprojectionMatrix())->multiply(const_cast<Vector &>(reducedCorrection), updateComponents, 1.0, FullSquareMatrix::TRANSPOSED);

  // updateComponents <- normalMatrix^{-1} * updateComponents 
  parent_->solver()->solution(updateComponents);
};

/* UpdateProjectorImpl implementation */

UpdateProjectorImpl::UpdateProjectorImpl(const UpdateProjectorImpl::Manager * manager) :
  UpdateProjector(),
  manager_(manager),
  jumpReactor_(new JumpReactor(NULL, this)),
  correctionReactor_(new CorrectionReactor(NULL, this))
{}

void
UpdateProjectorImpl::doProjection() {
  // reducedJump <- 0
  reducedJump().reset(reductionBasis()->stateCount());
 
  // Fill-in relevant part of reducedJump 
  int factorRank = solver()->factorRank(); 
  for (int i = 0; i < factorRank; ++i) {
    int index = solver()->factorPermutation(i);
    reducedJump()[index] = jump()->state() * reductionBasis()->state(index);
  }
  
  log() << "Reduced jump of size " << reducedJump().size() << "\n";
}

size_t
UpdateProjectorImpl::reducedBasisSize() const {
  return reductionBasis()->stateCount(); 
}

void
UpdateProjectorImpl::jumpIs(const Seed * j) {
  setJump(j);
}

void
UpdateProjectorImpl::reducedCorrectionIs(const ReducedSeed * rc) {
  correctionReactor_->notifierIs(rc);
  setReducedCorrection(rc);
}

/* UpdateProjectorImpl::Manager implementation */

UpdateProjectorImpl::Manager::Manager(const DynamStateBasis * drb,
                                      const FullSquareMatrix * drm,
                                      const RankDeficientSolver * ds) :
  defaultReductionBasis_(drb),
  defaultReprojectionMatrix_(drm),
  defaultSolver_(ds)
{}

UpdateProjectorImpl *
UpdateProjectorImpl::Manager::createNewInstance(const String & key) {
  return new UpdateProjectorImpl(this);
}

} /* namespace Pita */
