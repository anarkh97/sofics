#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "../DynamState.h"
#include "../LinearDynamOps.h"
#include "../DynamStateOps.h"

#include <Problems.d/DynamDescr.h>
#include <Utils.d/SolverInfo.h>
#include <Driver.d/Domain.h>

#include "SliceMapping.h"
#include "../RemoteStateMpiImpl.h"

#include "../IncrementalPropagation.h"
#include "../PostProcessingManager.h"
#include "../IncrementalPostProcessor.h"

#include "../LinearGenAlphaIntegrator.h"
#include "LinearPropagatorManager.h"

#include "LinearProjectionNetwork.h"
#include "../NearSymmetricSolver.h"

#include "JumpConvergenceEvaluator.h"
#include "../SeedDifferenceEvaluator.h"

#include "../FullCorrectionPropagatorImpl.h"

#include "../IntegratorSeedInitializer.h"
#include "RemoteSeedInitializerServer.h"
#include "../RemoteSeedInitializerProxy.h"
#include "../UserProvidedSeedInitializer.h"
#include "../SimpleSeedInitializer.h"

#include "HomogeneousTaskManager.h"
#include "NonHomogeneousTaskManager.h"
#include "../TimedExecution.h"

#include "../CommSplitter.h"
#include "../IntegratorPropagator.h"
#include "../RemoteDynamPropagatorProxy.h"
#include "RemoteCoarseCorrectionServer.h"

#include <Timers.d/GetTime.h>

namespace Pita { namespace Std {

class ReducedLinearDriverImpl::BasisSizeReactor : public LinearProjectionNetwork::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisSizeReactor);

  virtual void onProjectionOperators();

  BasisSizeReactor(const LinearProjectionNetwork * notifier, RemoteState::MpiManager * target);

private:
  RemoteState::MpiManager::Ptr target_;
};

ReducedLinearDriverImpl::BasisSizeReactor::BasisSizeReactor(const LinearProjectionNetwork * notifier,
                                                            RemoteState::MpiManager * target) :
  LinearProjectionNetwork::NotifieeConst(notifier),
  target_(target)
{}

void
ReducedLinearDriverImpl::BasisSizeReactor::onProjectionOperators() {
  target_->reducedStateSizeIs(notifier()->reducedBasisSize());
}

ReducedLinearDriverImpl::ReducedLinearDriverImpl(SingleDomainDynamic * pbDesc,
                                                 GeoSource * geoSource,
                                                 Domain * domain,
                                                 SolverInfo * solverInfo,
                                                 Communicator * baseComm) :
  LinearDriverImpl(pbDesc, geoSource, domain, solverInfo, baseComm)
{}

// Main routine
void
ReducedLinearDriverImpl::solve() {
  log() << "Begin Linear Pita\n";
  double tic = getTime(); // Total time

  preprocess();
  log() << "\n";

  printSummary();
  log() << "\n";

  // Run process-specific tasks
  if (remoteCoarse_) {
    CommSplitter::Ptr commSplitter = CommSplitter::New(baseComm(), CpuCount(1)); // Specializes one cpu
    if (commSplitter->localGroup() == CommSplitter::STANDARD) {
      solveParallel(commSplitter->splitComm(), commSplitter->interComm());
    } else {
      solveCoarse(commSplitter->interComm());
    }
  } else {
    solveParallel(baseComm(), NULL);
  }

  double toc = getTime();
  log() << "\n" << "Total time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\n" << "End Linear Pita\n";
}

void
ReducedLinearDriverImpl::preprocess() {
  double tic = getTime();
  
  probDesc()->preProcess();
  
  // Space-domain 
  vectorSize_ = probDesc()->solVecInfo();
  dynamOpsMgr_= LinearDynamOps::Manager::New(probDesc());

  // Time-domain 
  fineTimeStep_ = Seconds(solverInfo()->getTimeStep());
  sliceRatio_ = TimeStepCount(solverInfo()->pitaTimeGridRatio);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(solverInfo()->initialTime);
  finalTime_ = Seconds(solverInfo()->tmax);
 
  SliceCount numSlices(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (sliceRatio_.value() * fineTimeStep_.value()))));
  finalTime_ = fineTimeStep_.operator*(Seconds(numSlices.value() * sliceRatio_.value())); // Must have only complete time-slices
  
  // Main options 
  noForce_ = solverInfo()->pitaNoForce;
  userProvidedSeeds_ = solverInfo()->pitaReadInitSeed;
  remoteCoarse_ = solverInfo()->pitaRemoteCoarse && (baseComm()->numCPUs() > 1);

  // Load balancing 
  CpuCount numCpus(baseComm()->numCPUs() - (remoteCoarse_ ? 1 : 0));
  SliceCount maxActive(solverInfo()->pitaProcessWorkloadMax);
  mapping_ = SliceMapping::New(numSlices, numCpus, maxActive);

  // Other parameters 
  lastIteration_ = IterationRank(solverInfo()->pitaMainIterMax);
  jumpCvgRatio_ = solverInfo()->pitaJumpCvgRatio;
  if (jumpCvgRatio_ == 0.0) {
    // Use default value
    const int schemeOrder = 2;
    jumpCvgRatio_ = std::pow(static_cast<double>(sliceRatio_.value()), schemeOrder);
  }
  projectorTolerance_ = solverInfo()->pitaProjTol;
  coarseRhoInfinity_ = 1.0; // TODO Could be set in input file

  // PITA-specific output
  jumpMagnOutput_ = solverInfo()->pitaJumpMagnOutput;

  double toc = getTime();
  log() << "\n";
  log() << "Total preprocessing time = " << (toc - tic) / 1000.0 << " s\n";
}

void
ReducedLinearDriverImpl::printSummary() {
  log() << "Slices = " << mapping_->totalSlices() << ", MaxActive = " << mapping_->maxWorkload() << ", Cpus = " << mapping_->availableCpus() << "\n";
  log() << "dt = " << fineTimeStep_ << ", J = " << sliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", Tf = Slices*Dt = " << finalTime_ << "\n";
  log() << "Iteration count = " << lastIteration_ << "\n"; 
  if (jumpCvgRatio_ >= 0.0) { log() << "Jump-based convergence ratio = " << jumpCvgRatio_ << "\n"; }
  if (noForce_) { log() << "No external force\n"; }
  if (remoteCoarse_) { log() << "Remote coarse time-integration\n"; }
  if (userProvidedSeeds_) { log() << "Reading user-provided initial seed information\n"; }
  log() << "VectorSize = " << vectorSize_ << " dofs\n";
  log() << "Projector tol = " << projectorTolerance_ << "\n";
}

void
ReducedLinearDriverImpl::solveParallel(Communicator * timeComm, Communicator * coarseComm) {
  log() << "\n";
  log() << "Parallel solver initialization\n";

  double tic = getTime();
   
  // Linear operators
  const double fineRhoInfinity = 1.0; // No numerical damping
  GeneralizedAlphaParameter fineIntegrationParam(fineTimeStep_, fineRhoInfinity);
  LinearDynamOps::Ptr dynOps = dynamOpsMgr_->dynOpsNew(fineIntegrationParam); // Computationally expensive step

  // Projection-based correction
  RankDeficientSolver::Ptr normalMatrixSolver = NearSymmetricSolver::New(projectorTolerance_);
  LinearProjectionNetwork::Ptr projectionMgr = LinearProjectionNetwork::New(mapping_.ptr(),
                                                                            timeComm,
                                                                            dynOps.ptr(),
                                                                            vectorSize_,
                                                                            normalMatrixSolver.ptr());
  
  // Fine-grid time integration and feedback into correction
  LinearGenAlphaIntegrator::Ptr fineIntegrator;
  AffineDynamPropagator::ConstantTerm constTermStatus;
  if (noForce_) {
    fineIntegrator = new HomogeneousGenAlphaIntegrator(dynamOpsMgr_.ptr(), fineIntegrationParam);
    constTermStatus = AffineDynamPropagator::HOMOGENEOUS;
  } else {
    fineIntegrator = new AffineGenAlphaIntegrator(dynamOpsMgr_.ptr(), fineIntegrationParam);
    constTermStatus = AffineDynamPropagator::NONHOMOGENEOUS;
  }
  PostProcessing::Manager::Ptr postProcessingMgr = buildPostProcessor(CpuRank(timeComm->myID()));
  AffineBasisCollector::Ptr collector = projectionMgr->collector();
  LinearPropagatorManager::Ptr finePropagatorManager = LinearPropagatorManager::New(
      fineIntegrator.ptr(), postProcessingMgr.ptr(), collector.ptr(),
      sliceRatio_, initialTime_, constTermStatus);

  // MPI-based point-to-point communication
  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(timeComm, vectorSize_);
  BasisSizeReactor::PtrConst basisSizeReactor = new BasisSizeReactor(projectionMgr.ptr(), commMgr.ptr());

  // Jump-based convergence policy
  JumpConvergenceEvaluator::Ptr jumpCvgEval;
  if (jumpCvgRatio_ >= 0.0) {
    jumpCvgEval = AccumulatedJumpConvergenceEvaluator::New(jumpCvgRatio_, dynOps.ptr(), mapping_.ptr(), timeComm);
  } else {
    jumpCvgEval = TrivialConvergenceEvaluator::New(mapping_.ptr());
  }
 
  // Jump output
  LinSeedDifferenceEvaluator::Manager::Ptr jumpOutMgr;
  if (jumpMagnOutput_) {
    jumpOutMgr = LinSeedDifferenceEvaluator::Manager::New(dynOps.ptr());
  }

  // Initial seed information on coarse time-grid
  SeedInitializer::Ptr seedInitializer = buildSeedInitializer(coarseComm);
  
  // Generate tasks
  TaskManager::Ptr taskMgr;
  if (noForce_) {
    taskMgr = new HomogeneousTaskManager(mapping_.ptr(),
                                         commMgr.ptr(),
                                         finePropagatorManager.ptr(),
                                         projectionMgr.ptr(),
                                         jumpCvgEval.ptr(),
                                         jumpOutMgr.ptr(),
                                         seedInitializer.ptr());
  } else {
    // Coarse time-grid propagator
    CorrectionPropagator<DynamState>::Manager::Ptr fullCorrPropMgr = buildCoarseCorrection(coarseComm);
    taskMgr = new NonHomogeneousTaskManager(mapping_.ptr(),
                                            commMgr.ptr(),
                                            finePropagatorManager.ptr(),
                                            projectionMgr.ptr(),
                                            jumpCvgEval.ptr(),
                                            jumpOutMgr.ptr(),
                                            seedInitializer.ptr(),
                                            fullCorrPropMgr.ptr());
  }

  TimedExecution::Ptr execution = TimedExecution::New(taskMgr.ptr());
  
  double toc = getTime();
  log() << "\n";
  log() << "Total initialization time = " << (toc - tic) / 1000.0 << " s\n";
  tic = toc;
 
  // Perform the algorithm 
  execution->targetIterationIs(lastIteration_);

  toc = getTime();
  log() << "\n";
  log() << "Total solve time = " << (toc - tic) / 1000.0 << " s\n";
}

void
ReducedLinearDriverImpl::solveCoarse(Communicator * clientComm) {
  double tic = getTime();
  RemoteCoarseServer::Ptr coarseServer;

  if (noForce_) {
    SeedInitializer::Ptr seedInitializer = buildSeedInitializer();
    coarseServer = RemoteSeedInitializerServer::New(clientComm, seedInitializer.ptr(), mapping_.ptr());

    log() << "\n";
    log() << "Remote coarse initialization\n";
  } else {
    DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator();

    RemoteDynamPropagatorServer::Ptr propagatorServer = new RemoteDynamPropagatorServer(coarsePropagator.ptr(), clientComm);
    coarseServer = RemoteCoarseCorrectionServer::New(propagatorServer.ptr(), mapping_.ptr());

    mapping_->convergedSlicesInc(); // Since the coarse correction occurs at the second iteration (iteration 0 after iteration -1)
    
    log() << "\n";
    log() << "Remote coarse propagation\n";
  }
  double toc = getTime();
  log() << "\n";
  log() << "Total initialization time = " << (toc - tic) / 1000.0 << " s\n";
  tic = getTime();
  
  coarseServer->statusIs(RemoteCoarseServer::BUSY);
  
  toc = getTime();
  log() << "\n";
  log() << "Total solve time = " << (toc - tic) / 1000.0 << " s\n";
}


PostProcessing::Manager::Ptr
ReducedLinearDriverImpl::buildPostProcessor(CpuRank localCpu) const {
  std::vector<int> localFileId;
  for (SliceMapping::SliceIterator it = mapping_->hostedSlice(localCpu); it; ++it) {
    localFileId.push_back((*it).value());
  }
  IncrementalPostProcessor::Ptr pitaPostProcessor = IncrementalPostProcessor::New(geoSource(), localFileId.size(), &localFileId[0], probDesc()->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<IncrementalPostProcessor> LinearIntegratorReactor;
  return PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());
}

SeedInitializer::Ptr
ReducedLinearDriverImpl::buildSeedInitializer(Communicator * clientComm) const {
  if (userProvidedSeeds_) {
    return UserProvidedSeedInitializer::New(vectorSize_, geoSource(), domain());
  }

  if (!noForce_) {
    return SimpleSeedInitializer::New(initialSeed());
  }

  if (clientComm == NULL) {
    // Local time-integration
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    coarseIntegrator->initialConditionIs(initialSeed(), initialTime_);
    return IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1));
  }

  return RemoteSeedInitializerProxy::New(clientComm, vectorSize_);
}

LinearGenAlphaIntegrator::Ptr
ReducedLinearDriverImpl::buildCoarseIntegrator() const {
  return new HomogeneousGenAlphaIntegrator(dynamOpsMgr_.ptr(), GeneralizedAlphaParameter(coarseTimeStep_, coarseRhoInfinity_));
}


CorrectionPropagator<DynamState>::Manager::Ptr
ReducedLinearDriverImpl::buildCoarseCorrection(Communicator * coarseComm) const {
  DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator(coarseComm);
  return FullCorrectionPropagatorImpl::Manager::New(coarsePropagator.ptr());
}

DynamPropagator::Ptr
ReducedLinearDriverImpl::buildCoarsePropagator(Communicator * coarseComm) const {
  if (!coarseComm) {
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    IntegratorPropagator::Ptr localPropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
    localPropagator->timeStepCountIs(TimeStepCount(1));
    return localPropagator;
  }

  return RemoteDynamPropagatorProxy::New(vectorSize_, coarseComm, CpuRank(0));
}

} /* end namespace Std */ } /* end namespace Pita */

// Global state
extern GeoSource * geoSource;
extern Communicator * structCom;
extern Domain * domain;

// Entrypoint
Pita::LinearDriver::Ptr
linearPitaDriverNew(SingleDomainDynamic * pbDesc) {
  return Pita::Std::ReducedLinearDriverImpl::New(pbDesc, geoSource, domain, &domain->solInfo(), structCom);
}
