#include "LinearTaskManager.h"

#include "ReducedCorrectionManager.h"

namespace Pita { namespace Hts {

class ProjectionBasis : public NamedTask {
public:
  void iterationIs(IterationRank i) {
    correctionMgr_->buildProjection();
    size_t reducedBasisSize = correctionMgr_->reducedBasisSize();
    commMgr_->reducedStateSizeIs(reducedBasisSize);
  }

  ProjectionBasis(LinearProjectionNetwork * correctionMgr, RemoteState::MpiManager * commMgr) :
    NamedTask("Projection building"),
    correctionMgr_(correctionMgr),
    commMgr_(commMgr)
  {}

private:
  LinearProjectionNetwork::Ptr correctionMgr_;
  RemoteState::MpiManager::Ptr commMgr_;
};

LinearTaskManager::LinearTaskManager(IterationRank initialIteration,
                                     SliceMapping * mapping,
                                     AffinePropagatorManager * propMgr,
                                     CorrectionPropagator<DynamState>::Manager * fullCorrMgr,
                                     JumpConvergenceEvaluator * jumpCvgMgr,
                                     LinSeedDifferenceEvaluator::Manager * jumpErrorMgr,
                                     LinearProjectionNetwork * correctionMgr,
                                     RemoteState::MpiManager * commMgr) :
  TaskManager(initialIteration),
  network_(new LinearLocalNetwork(mapping, propMgr,
           new ReducedCorrectionManager(correctionMgr, fullCorrMgr),
           jumpCvgMgr, commMgr, jumpErrorMgr)),
  jumpCvgMgr_(jumpCvgMgr),
  correctionMgr_(correctionMgr),
  commMgr_(commMgr),
  phase_(),

  phaseIt_(new HtsPhaseIteratorImpl(phase_))
{}

void
LinearTaskManager::scheduleNormalIteration() {
  jumpCvgMgr()->iterationIs(iteration().next()); // TODO Better
  network()->applyConvergenceStatus();

  scheduleCorrection();
  scheduleFinePropagation();
}

void
LinearTaskManager::scheduleFinePropagation() {
  schedulePhase("Fine propagation", network()->activeHalfTimeSlices());
  schedulePhase("Propagated seed synchronization", network()->activeLeftSeedSyncs());
  schedulePhase("Jump evaluation", network()->activeJumpAssemblers());
}

void
LinearTaskManager::scheduleCorrection() {
  // Projection Basis
  {
    TaskList projectionBuilding;
    projectionBuilding.push_back(new ProjectionBasis(correctionMgr_.ptr(), commMgr_.ptr()));
    Phase::Ptr projectionBasis = phaseNew("Projection basis", projectionBuilding);
    phases().push_back(projectionBasis);
  }
  schedulePhase("Jump Projection", network()->activeJumpProjectors());
  schedulePhase("Correction", network()->activeFullTimeSlices());
  schedulePhase("Correction synchronization", network()->activeCorrectionSyncs());
  schedulePhase("Seed update", network()->activeSeedAssemblers());
}

void
LinearTaskManager::schedulePhase(const String & phaseName, const LinearLocalNetwork::TaskList & networkTaskList) {
  TaskList taskList(networkTaskList.begin(), networkTaskList.end());

  Phase::Ptr projectionBasis = phaseNew(phaseName, taskList);
  phases().push_back(projectionBasis);
}

} /* end namespace Hts */ } /* end namespace Pita */
