#include "NlProjectionNetwork.h"
#include "../NearSymmetricSolver.h"
#include "../DynamStateOps.h"

#include <cassert>

namespace Pita { namespace Hts {

NlProjectionNetwork::NlProjectionNetwork(GlobalStateSharing * sharing,
                                         NlDynamOps * dynamOps,
                                         size_t vectorSize, double projectionTolerance) :
  vectorSize_(vectorSize),
  sharing_(sharing),
  sharedProjectionBasis_(DynamStatePlainBasis::New(vectorSize)),
  midBasisMgr_(new MidBasisManager(sharedProjectionBasis_.ptr())),
  concurrentMgr_(ConcurrentBasis::Manager::New()),
  endBasisMgr_(SimpleBasisManager::New(vectorSize)),
  projBasisMgr_(SimpleBasisManager::New(vectorSize)),
  solverMgr_(SimpleSolverManager<NearSymmetricSolver>::New(projectionTolerance)),
  condensMgr_(NULL),
  projBuildMgr_(NULL),
  corrRedMgr_(NULL),
  corrReconMgr_(NULL)
{
  double condensationTolerance = projectionTolerance * 0.1; // TODO as parameter
  condensMgr_ = new BasisCondensationManager(vectorSize, condensationTolerance, endBasisMgr_.ptr(), midBasisMgr_.ptr());
  projBuildMgr_ = new ProjectionBuildingFactory(dynamOps, endBasisMgr_.ptr(), projBasisMgr_.ptr(), solverMgr_.ptr()); 

  CorrectionReductor::Manager::OperatorManager::Ptr redMgr = ReductorManager::New(projBasisMgr_.ptr(), solverMgr_.ptr());
  corrRedMgr_ = CorrectionReductor::Manager::New(redMgr.ptr());
  
  CorrectionReconstructor::Manager::OperatorManager::Ptr reconMgr = ReconstructorManager::New(endBasisMgr_.ptr());
  corrReconMgr_ = CorrectionReconstructor::Manager::New(reconMgr.ptr());
}


BasisCondensation::BasisCondensation(const String & name,
                                     const DynamStateBasis * origin,
                                     DynamStatePlainBasis * target,
                                     BasisTransform * operation) :
  NamedTask(name),
  origin_(origin),
  target_(target),
  operation_(operation)
{}

void
BasisCondensation::iterationIs(IterationRank ir) {
  operation_->inputValueIs(*origin_);

  target_->stateBasisDel();
  target_->lastStateBasisIs(operation_->output());

  log() << "*** Condensation: " << origin_->stateCount() << " to " << target_->stateCount() << "\n";

  setIteration(ir);
}


DynamStateReductor *
ReductorManager::instance(const HalfSliceRank & key) const {
  DynamStatePlainBasis::Ptr basis = basisMgr_->instance(keyToKey(key));
  RankDeficientSolver::Ptr solver = solverMgr_->instance(keyToKey(key));
  
  assert(basis);
  assert(solver);

  return new DynamStateReductor(basis.ptr(), solver.ptr());
}

size_t
ReductorManager::instanceCount() const {
  return basisMgr_->instanceCount();
}

DynamStateReductor *
ReductorManager::instanceNew(const HalfSliceRank & key) {
  HalfSliceId id = keyToKey(key);

  DynamStatePlainBasis::Ptr basis = basisMgr_->instance(id);
  if (!basis) {
    basis = basisMgr_->instanceNew(id);
  }

  RankDeficientSolver::Ptr solver = solverMgr_->instance(id);
  if (!solver) {
    solver = solverMgr_->instanceNew(id);
  }

  return new DynamStateReductor(basis.ptr(), solver.ptr());
}

void
ReductorManager::instanceDel(const HalfSliceRank & key) {
  basisMgr_->instanceDel(keyToKey(key));
  solverMgr_->instanceDel(keyToKey(key));
}

DynamStateReconstructor *
ReconstructorManager::instance(const HalfSliceRank & key) const {
  DynamStatePlainBasis::Ptr basis = basisMgr_->instance(keyToKey(key));
  return new DynamStateReconstructor(basis.ptr());
}

size_t
ReconstructorManager::instanceCount() const {
  return basisMgr_->instanceCount();
}

DynamStateReconstructor *
ReconstructorManager::instanceNew(const HalfSliceRank & key) {
  HalfSliceId id = keyToKey(key);

  DynamStatePlainBasis::Ptr basis = basisMgr_->instance(id);
  if (!basis) {
    basis =  basisMgr_->instanceNew(id);
  }

  return new DynamStateReconstructor(basis.ptr());
}

void
ReconstructorManager::instanceDel(const HalfSliceRank & key) {
  basisMgr_->instanceDel(keyToKey(key));
}

ProjectionBuilding::ProjectionBuilding(const String & name,
                                       NlDynamOps * dynOps,
                                       const DynamStateBasis * projectionBasis,
                                       DynamStatePlainBasis * dualBasis,
                                       RankDeficientSolver * solver) :
  NamedTask(name),
  dynamOps_(dynOps),
  refDisp_(projectionBasis->vectorSize(), 0.0),
  refTime_(0.0),
  projectionBasis_(projectionBasis),
  dualBasis_(dualBasis),
  solver_(solver)
{
  assert(projectionBasis->vectorSize() == dualBasis->vectorSize());
}

void
ProjectionBuilding::iterationIs(IterationRank ir) {
  // Rebuild metric Q
  dynamOps_->displacementIs(refDisp_, refTime_);

  // Compute dual basis D = QV where V = projectionBasis
  dualBasis_->stateBasisDel();
  for (DynamStateBasis::IteratorConst it = projectionBasis()->state(); it; ++it) {
    dualBasis_->lastStateIs(mult(dynamOps(), *it));
  }

  // Form transposed normal matrix V^T D == (D^T V)^T
  FullSquareMatrix normalMatrix(projectionBasis()->stateCount());
  double * m = normalMatrix.data();

  for (DynamStateBasis::IteratorConst row_it = projectionBasis()->state(); row_it; ++row_it) {
    for (DynamStateBasis::IteratorConst col_it = dualBasis()->state(); col_it; ++col_it) {
      *m = (*col_it) * (*row_it);
      ++m;
    }
  }

  // Perform factorization, leave original ordering
  solver_->transposedMatrixIs(normalMatrix);

  log() << "*** Numerical rank = " << solver_->factorRank() << " / " << solver_->matrixSize()  << "\n";

  // Publish result
  setIteration(ir);
}

void
ProjectionBuilding::referenceDisplacementIs(const GenVector<double> & disp, Seconds time) {
  refDisp_ = disp;
  refTime_ = time;
}

void
ProjectionBuilding::toleranceIs(double tol) {
  solver_->toleranceIs(tol);
}

ProjectionBuildingFactory::ProjectionBuildingFactory(NlDynamOps * dynamOps,
                                                     SimpleBasisManager * endBasisMgr,
                                                     SimpleBasisManager * projBasisMgr,
                                                     SimpleSolverManager<NearSymmetricSolver> * solverMgr) :
  dynamOps_(dynamOps),
  endBasisMgr_(endBasisMgr),
  projBasisMgr_(projBasisMgr),
  solverMgr_(solverMgr)
{}

ProjectionBuilding *
ProjectionBuildingFactory::instanceNew(HalfSliceRank key) {
  String taskName = "Projection Building " + toString(key);
  HalfSliceId sliceId(key, BACKWARD);

  DynamStateBasis::PtrConst propagatedBasis = endBasisMgr_->instance(sliceId);
  if (!propagatedBasis) {
    propagatedBasis = endBasisMgr_->instanceNew(sliceId);
  }
  
  DynamStatePlainBasis::Ptr dualBasis = projBasisMgr_->instance(sliceId);
  if (!dualBasis) {
    dualBasis = projBasisMgr_->instanceNew(sliceId);
  }

  RankDeficientSolver::Ptr solver = solverMgr_->instance(sliceId);
  if (!solver) {
    solver = solverMgr_->instanceNew(sliceId);
  }

  return new ProjectionBuilding(taskName, dynamOps_.ptr(), propagatedBasis.ptr(), dualBasis.ptr(), solver.ptr());
}

BasisCondensationManager::BasisCondensationManager(size_t vectorSize, double tolerance,
                                                   SimpleBasisManager * endBasisMgr, MidBasisManager * midBasisMgr) :
  vectorSize_(vectorSize),
  tolerance_(tolerance),
  endBasisMgr_(endBasisMgr),
  midBasisMgr_(midBasisMgr)
{}

BasisCondensation *
BasisCondensationManager::createNewInstance(const HalfSliceId & sliceId) {
  HalfSliceRank basisRank = (sliceId.direction() == FORWARD) ? sliceId.rank() : sliceId.rank().next();
  DynamStateBasis::PtrConst origin = midBasisMgr_->instance(basisRank);
  if (!origin) {
    origin = midBasisMgr_->instanceNew(basisRank);
  }
  
  DynamStatePlainBasis::Ptr target = endBasisMgr_->instance(sliceId);
  if (!target) {
    target = endBasisMgr_->instanceNew(sliceId);
  }

  BasisTransform::Ptr operation = new CholeskyOrtho(vectorSize_, tolerance_);

  return new BasisCondensation("Condense Basis " + toString(sliceId), origin.ptr(), target.ptr(), operation.ptr());
}

} /* end namespace Hts */ } /* end namespace Pita */
