#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "../DynamState.h"
#include "../LinearDynamOps.h"
#include "../DynamStateOps.h"

#include <Problems.d/DynamDescr.h>
#include <Utils.d/SolverInfo.h>
#include <Driver.d/Domain.h>

#include "../LinearGenAlphaIntegrator.h"

#include "../IntegratorPropagator.h"
#include "../AffineIntegratorPropagator.h"

#include "SliceMapping.h"

#include "../RankDeficientSolver.h"
#include "../NearSymmetricSolver.h"
#include "../PivotedCholeskySolver.h"
#include "../LeastSquareSolver.h"

#include "LinearProjectionNetwork.h"

#include "../Seed.h"

#include "../RemoteStateMpiImpl.h"

#include "../IncrementalPostProcessor.h"
#include "LinearFineIntegratorManager.h"
#include "AffinePropagatorManager.h"

#include "../FullCorrectionPropagatorImpl.h"
#include "JumpConvergenceEvaluator.h"

#include "../IntegratorSeedInitializer.h"
#include "RemoteSeedInitializerServer.h"
#include "../RemoteSeedInitializerProxy.h"
#include "../UserProvidedSeedInitializer.h"
#include "../SimpleSeedInitializer.h"

#include "../RemoteDynamPropagatorProxy.h"
#include "RemoteCoarseCorrectionServer.h"

#include "../CommSplitter.h"
#include <Comm.d/Communicator.h>

#include <Driver.d/Domain.h>

#include "../TimedExecution.h"
#include "HomogeneousTaskManager.h"
#include "NonHomogeneousTaskManager.h"

#include "../SeedDifferenceEvaluator.h"

namespace Pita { namespace Hts {

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
  log() << "Begin Reversible Linear Pita\n";
  double tic = getTime(); // Total time

  preprocess();

  /* Summarize problem and parameters */
  log() << "\n"; 
  log() << "Slices = " << mapping_->totalSlices() << ", MaxActive = " << mapping_->maxWorkload() << ", Cpus = " << mapping_->availableCpus() << "\n";
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", Tf = Slices*(J/2)*dt = " << finalTime_ << "\n";
  log() << "Iteration count = " << lastIteration_ << "\n"; 
  if (jumpCvgRatio_ >= 0.0) { log() << "Jump-based convergence ratio = " << jumpCvgRatio_ << "\n"; }
  if (noForce_) { log() << "No external force\n"; }
  if (remoteCoarse_) { log() << "Remote coarse time-integration\n"; }
  if (userProvidedSeeds_) { log() << "Reading user-provided initial seed information\n"; }
  log() << "VectorSize = " << vectorSize_ << " dofs\n";
  log() << "Projector tol = " << projectorTolerance_ << "\n";

  /* Determine process task */
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
  log() << "\n" << "End Reversible Linear Pita\n";
}

void
ReducedLinearDriverImpl::preprocess() {
  double tic = getTime();
  
  probDesc()->preProcess();
  
  /* Space-domain */
  vectorSize_ = probDesc()->solVecInfo();
  dynamOpsMgr_= LinearDynamOps::Manager::New(probDesc());

  /* Time-domain */
  fineTimeStep_ = Seconds(solverInfo()->getTimeStep());
  halfSliceRatio_ = TimeStepCount(solverInfo()->pitaTimeGridRatio / 2);
  sliceRatio_ = TimeStepCount(halfSliceRatio_.value() * 2);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(solverInfo()->initialTime);
  finalTime_ = Seconds(solverInfo()->tmax);
 
  HalfSliceCount numSlices(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (halfSliceRatio_.value() * fineTimeStep_.value()))));
  FullSliceCount fullTimeSlices((numSlices.value() / 2) + (numSlices.value() % 2));
  numSlices = HalfSliceCount(fullTimeSlices.value() * 2);
  finalTime_ = fineTimeStep_.operator*(Seconds(numSlices.value() * halfSliceRatio_.value())); // To have a whole number of full time-slices 
  
  /* Main options */
  noForce_ = solverInfo()->pitaNoForce;
  userProvidedSeeds_ = solverInfo()->pitaReadInitSeed;
  remoteCoarse_ = solverInfo()->pitaRemoteCoarse && (baseComm()->numCPUs() > 1);

  /* Load balancing */ 
  CpuCount numCpus(baseComm()->numCPUs() - (remoteCoarse_ ? 1 : 0));
  HalfSliceCount maxActive(solverInfo()->pitaProcessWorkloadMax);
  mapping_ = SliceMapping::New(fullTimeSlices, numCpus, maxActive);

  /* Other parameters */ 
  lastIteration_ = IterationRank(solverInfo()->pitaMainIterMax);
  jumpCvgRatio_ = solverInfo()->pitaJumpCvgRatio;
  if (jumpCvgRatio_ == 0.0) {
    // Use default value
    const int schemeOrder = 2;
    jumpCvgRatio_ = std::pow(static_cast<double>(sliceRatio_.value()), schemeOrder);
  }
  projectorTolerance_ = solverInfo()->pitaProjTol;
  coarseRhoInfinity_ = 1.0; // TODO Could be set in input file

  /* PITA-specific output */
  jumpMagnOutput_ = solverInfo()->pitaJumpMagnOutput;

  double toc = getTime();
  log() << "\n";
  log() << "Total preprocessing time = " << (toc - tic) / 1000.0 << " s\n";
}


void
ReducedLinearDriverImpl::solveParallel(Communicator * timeComm, Communicator * coarseComm) {
  log() << "\n";
  log() << "Parallel solver initialization\n";

  double tic = getTime();
  
  /* Local process identity */
  CpuRank myCpu(timeComm->myID());

  /* Fine time integration operators */
  const double noDampingRhoInfinity = 1.0; // No damping allowed for time-reversible
  GeneralizedAlphaParameter fineIntegrationParam(fineTimeStep_, noDampingRhoInfinity);
  LinearDynamOps::Ptr dynamOps = dynamOpsMgr_->dynOpsNew(fineIntegrationParam); // Expensive operation 
  
  /* Post processing */ 
  PostProcessing::Manager::Ptr postProcessingMgr = buildPostProcessor(myCpu);

  /* Correction */
  RankDeficientSolver::Ptr normalMatrixSolver = NearSymmetricSolver::New(projectorTolerance_); // TODO Other solver implementations ?
  
  LinearProjectionNetwork::Ptr correctionMgr = LinearProjectionNetwork::New(
      vectorSize_,
      timeComm,
      mapping_.ptr(),
      dynamOps.ptr(),
      normalMatrixSolver.ptr());

  /* Local fine propagation */
  GenFineIntegratorManager<AffineGenAlphaIntegrator>::Ptr fineIntegratorMgr = LinearFineIntegratorManager<AffineGenAlphaIntegrator>::New(dynamOpsMgr_.ptr(), fineIntegrationParam);

  AffinePropagatorManager::Ptr propagatorMgr = AffinePropagatorManager::New(
      correctionMgr->collector(),
      fineIntegratorMgr.ptr(), 
      postProcessingMgr.ptr(),
      halfSliceRatio_,
      initialTime_,
      noForce_ ? AffineDynamPropagator::HOMOGENEOUS : AffineDynamPropagator::NONHOMOGENEOUS);

  /* Local communication */
  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(timeComm, vectorSize_);

  // Convergence criterion
  JumpConvergenceEvaluator::Ptr jumpCvgMgr;
  if (jumpCvgRatio_ >= 0.0) {
    jumpCvgMgr = AccumulatedJumpConvergenceEvaluator::New(jumpCvgRatio_, dynamOps.ptr(), mapping_.ptr(), timeComm);
  } else {
    jumpCvgMgr = TrivialConvergenceEvaluator::New(mapping_.ptr());
  }

  /* Jump error */
  LinSeedDifferenceEvaluator::Manager::Ptr jumpErrorMgr;
  if (jumpMagnOutput_) {
    jumpErrorMgr = LinSeedDifferenceEvaluator::Manager::New(dynamOps.ptr());
  }

  /* Seed initialization */
  SeedInitializer::Ptr seedInitializer = buildSeedInitializer(coarseComm); 
  
  double toc = getTime();
  log() << "\n";
  log() << "Total initialization time = " << (toc - tic) / 1000.0 << " s\n";
  tic = toc;
  
  LinearTaskManager::Ptr taskManager;
  if (noForce_) {
    taskManager = new HomogeneousTaskManager(
        mapping_.ptr(),
        commMgr.ptr(),
        propagatorMgr.ptr(),
        correctionMgr.ptr(),
        jumpCvgMgr.ptr(),
        jumpErrorMgr.ptr(),
        seedInitializer.ptr());
  } else {
    CorrectionPropagator<DynamState>::Manager::Ptr fullCorrPropMgr = buildCoarseCorrection(coarseComm);
    taskManager = new NonHomogeneousTaskManager(
        mapping_.ptr(),
        commMgr.ptr(),
        propagatorMgr.ptr(),
        correctionMgr.ptr(),
        jumpCvgMgr.ptr(),
        jumpErrorMgr.ptr(),
        seedInitializer.ptr(),
        fullCorrPropMgr.ptr());
  }

  TimedExecution::Ptr timedExecution = TimedExecution::New(taskManager.ptr()); 
  timedExecution->targetIterationIs(lastIteration_);

  toc = getTime();
  log() << "\n";
  log() << "Total solve time = " << (toc - tic) / 1000.0 << " s\n";
}

void
ReducedLinearDriverImpl::solveCoarse(Communicator * timeComm) {
  double tic = getTime();
  RemoteCoarseServer::Ptr coarseServer;

  if (noForce_) {
    SeedInitializer::Ptr seedInitializer = buildSeedInitializer();
    coarseServer = RemoteSeedInitializerServer::New(timeComm, seedInitializer.ptr(), mapping_.ptr());

    log() << "\n";
    log() << "Remote coarse initialization\n";
  } else {
    DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator();

    RemoteDynamPropagatorServer::Ptr propagatorServer = new RemoteDynamPropagatorServer(coarsePropagator.ptr(), timeComm);
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

CorrectionPropagator<DynamState>::Manager::Ptr
ReducedLinearDriverImpl::buildCoarseCorrection(Communicator * coarseComm) const {
  DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator(coarseComm);
  return FullCorrectionPropagatorImpl::Manager::New(coarsePropagator.ptr());
}

LinearGenAlphaIntegrator::Ptr
ReducedLinearDriverImpl::buildCoarseIntegrator() const {
  return new HomogeneousGenAlphaIntegrator(dynamOpsMgr_.ptr(), GeneralizedAlphaParameter(coarseTimeStep_, coarseRhoInfinity_));
}

DynamPropagator::Ptr
ReducedLinearDriverImpl::buildCoarsePropagator(Communicator * coarseComm) const {
  if (coarseComm == NULL) {
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    IntegratorPropagator::Ptr localPropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
    localPropagator->timeStepCountIs(TimeStepCount(1));
    return localPropagator;
  }

  return RemoteDynamPropagatorProxy::New(vectorSize_, coarseComm, CpuRank(0));
}


SeedInitializer::Ptr
ReducedLinearDriverImpl::buildSeedInitializer(Communicator * timeComm) const {
  if (userProvidedSeeds_) {
    return UserProvidedSeedInitializer::New(vectorSize_, geoSource(), domain());
  }

  if (!noForce_) {
    return SimpleSeedInitializer::New(initialSeed());
  }

  if (timeComm == NULL) {
    // Local time-integration
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    coarseIntegrator->initialConditionIs(initialSeed(), initialTime_);
    return IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1));
  }

  return RemoteSeedInitializerProxy::New(timeComm, vectorSize_);
}

} /* end namespace Hts */ } /* end namespace Pita */

// Global state
extern GeoSource * geoSource;
extern Communicator * structCom;
extern Domain * domain;

/* Entrypoint */
Pita::LinearDriver::Ptr
linearReversiblePitaDriverNew(SingleDomainDynamic * pbDesc) {
  return Pita::Hts::ReducedLinearDriverImpl::New(pbDesc, geoSource, domain, &domain->solInfo(), structCom);
}
