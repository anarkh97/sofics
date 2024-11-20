#include "NlDriver.h"

#include "../NlDynamTimeIntegrator.h"

#include "../NlPostProcessor.h"
#include "../PostProcessingManager.h"

#include "../IntegratorPropagator.h"
#include "NlPropagatorManager.h"

#include "../IntegratorSeedInitializer.h"
#include "../UserProvidedSeedInitializer.h"

#include "../RemoteStateMpiImpl.h"
#include "GlobalStateSharing.h"

#include "JumpConvergenceEvaluator.h"
#include "../SeedDifferenceEvaluator.h"

#include "NlTaskManager.h"
#include "../TimedExecution.h"

#include <Timers.d/GetTime.h>
#include <Driver.d/GeoSource.h>

#include "../NlDynamOps.h"

namespace Pita { namespace Hts {

NlDriver::NlDriver(PitaNonLinDynamic * probDesc,
                   GeoSource * geoSource,
                   SolverInfo * solverInfo,
                   Communicator * baseComm) :
  NlDriverImpl(probDesc),
  geoSource_(geoSource),
  solverInfo_(solverInfo),
  baseComm_(baseComm)
{}

void
NlDriver::solve() {
  log() << "Begin Reversible NonLinear Pita\n";
  double tic = getTime();

  preprocess();
  summarizeParameters();
  solveParallel();

  double toc = getTime();
  log() << "\n" << "Total time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\n" << "End Reversible NonLinear Pita\n";
}

void
NlDriver::preprocess() {
  double tic = getTime();
  
  probDesc()->preProcess();
  
  // Space-domain 
  vectorSize_ = probDesc()->solVecInfo();

  // Time-domain 
  fineTimeStep_ = Seconds(solverInfo()->getTimeStep());
  halfSliceRatio_ = TimeStepCount(solverInfo()->pitaTimeGridRatio / 2);
  sliceRatio_ = TimeStepCount(halfSliceRatio_.value() * 2);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(solverInfo()->initialTime);
  finalTime_ = Seconds(solverInfo()->tmax);
 
  HalfSliceCount numSlices(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (halfSliceRatio_.value() * fineTimeStep_.value()))));
  FullSliceCount fullTimeSlices((numSlices.value() / 2) + (numSlices.value() % 2));
  numSlices = HalfSliceCount(fullTimeSlices.value() * 2);
  finalTime_ = fineTimeStep_.operator*(Seconds(numSlices.value() * halfSliceRatio_.value())); // To have a whole number of primal full time-slices 

  // Load balancing 
  CpuCount numCpus(baseComm()->numCPUs());
  localCpu_ = CpuRank(baseComm()->myID());
  HalfSliceCount maxActive(solverInfo()->pitaProcessWorkloadMax);
  mapping_ = SliceMapping::New(fullTimeSlices, numCpus, maxActive);

  // Other parameters
  lastIteration_ = IterationRank(solverInfo()->pitaMainIterMax);
  jumpCvgRatio_ = solverInfo()->pitaJumpCvgRatio;
  if (jumpCvgRatio_ == 0.0) {
    // Use default value
    const int schemeOrder = 2;
    jumpCvgRatio_ = std::pow(static_cast<double>(sliceRatio_.value()), schemeOrder);
  }
  projectorTolerance_ = solverInfo()->pitaProjTol;
  userProvidedSeeds_ = solverInfo()->pitaReadInitSeed;
  globalBasisEnrichment_ = solverInfo()->pitaGlobalBasisImprovement; // 1 => {M}, 2 => {L, R}, 3 => {M, L, R}

  // PITA-specific output
  jumpMagnOutput_ = solverInfo()->pitaJumpMagnOutput;

  double toc = getTime();
  log() << "\n";
  log() << "Total preprocessing time = " << (toc - tic) / 1000.0 << " s\n";
}

void
NlDriver::summarizeParameters() const {
  log() << "\n"; 
  log() << "Slices = " << mapping_->totalSlices() << ", MaxActive = " << mapping_->maxWorkload() << ", Cpus = " << mapping_->availableCpus() << "\n";
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", Tf = Slices*(J/2)*dt = " << finalTime_ << "\n";
  log() << "Iteration count = " << lastIteration_ << "\n"; 
  if (jumpCvgRatio_ >= 0.0) { log() << "Jump-based convergence ratio = " << jumpCvgRatio_ << "\n"; }
  if (userProvidedSeeds_) { log() << "Reading user-provided initial seed information\n"; }
  log() << "VectorSize = " << vectorSize_ << " dofs\n";
  log() << "Projector tol = " << projectorTolerance_ << "\n";
}

void
NlDriver::solveParallel() {
  log() << "\n";
  log() << "Parallel solver initialization\n";
  double tic = getTime();

  // Initial conditions
  Seconds initTime(0.0);
  DynamState initState = initialState();

  // Time-integration
  NlDynamTimeIntegrator::Ptr integrator = NlDynamTimeIntegrator::New(probDesc());
  NlPropagatorManager::Ptr propagatorMgr = NlPropagatorManager::New(
      integrator.ptr(),
      fineTimeStep_,
      halfSliceRatio_,
      initTime);

  // Initial Seeds
  SeedInitializer::Ptr seedInitializer;
  if (userProvidedSeeds_) {
    seedInitializer = UserProvidedSeedInitializer::New(vectorSize_, geoSource(), domain());
  } else {
    integrator->timeStepSizeIs(coarseTimeStep_);
    integrator->initialConditionIs(initState, initTime);
    seedInitializer = IntegratorSeedInitializer::New(integrator.ptr(), TimeStepCount(1));
  }
  
  // Post-processing
  PostProcessing::Manager::Ptr postProcessingMgr = buildPostProcessor();

  // Convergence criterion
  JumpConvergenceEvaluator::Ptr jumpCvgMgr;
  if (jumpCvgRatio_ >= 0.0) {
    NlDynamOps::Ptr dirtyOps = NlDynamOps::New(probDesc());
    jumpCvgMgr = AccumulatedJumpConvergenceEvaluator::New(jumpCvgRatio_, dirtyOps.ptr(), mapping_.ptr(), baseComm());
  } else {
    jumpCvgMgr = TrivialConvergenceEvaluator::New(mapping_.ptr());
  }

  // Jump evaluation (Output only)
  NonLinSeedDifferenceEvaluator::Manager::Ptr jumpEvalMgr;
  if (jumpMagnOutput_) {
    jumpEvalMgr = NonLinSeedDifferenceEvaluator::Manager::New(probDesc());
  }

  // Communications
  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(baseComm(), vectorSize_);

  // Basis update
  const SeedType enrichmentTypes[3] = { MAIN_SEED, LEFT_SEED, RIGHT_SEED };
  const GlobalStateSharing::Strategy globalUpdateStrategy(
      enrichmentTypes + (globalBasisEnrichment_ == 2),
      enrichmentTypes + 1 + 2 * (globalBasisEnrichment_ >= 2)); 
  GlobalStateSharing::Ptr basisUpdateMgr = new GlobalStateSharing(baseComm(), vectorSize_, globalUpdateStrategy);

  // Execute the algorithm 
  NlTaskManager::Ptr taskManager = new NlTaskManager(mapping_.ptr(), commMgr.ptr(),
                                                     propagatorMgr.ptr(),
                                                     seedInitializer.ptr(),
                                                     basisUpdateMgr.ptr(),
                                                     postProcessingMgr.ptr(),
                                                     jumpCvgMgr.ptr(), jumpEvalMgr.ptr(),
                                                     projectorTolerance_, lastIteration_);
  TimedExecution::Ptr execution = TimedExecution::New(taskManager.ptr()); 
  execution->targetIterationIs(lastIteration_);

  double toc = getTime();
  log() << "\n";
  log() << "Total solve time = " << (toc - tic) / 1000.0 << " s\n";
}

PostProcessing::Manager::Ptr
NlDriver::buildPostProcessor() { 
  std::vector<int> ts;
  for (SliceMapping::SliceIterator s = mapping_->hostedSlice(localCpu_); s; ++s) {
    ts.push_back((*s).value()); 
  }
  NlPostProcessor::Ptr pitaPostProcessor = NlPostProcessor::New(geoSource(), ts.size(), &ts[0], probDesc());
  typedef PostProcessing::IntegratorReactorImpl<NlPostProcessor> NlIntegratorReactor;
  NlIntegratorReactor::Builder::Ptr ppBuilder = NlIntegratorReactor::Builder::New(pitaPostProcessor.ptr());
  return PostProcessing::Manager::New(ppBuilder.ptr());
}

} /* end namespace Hts */ } /* end namespace Pita */

#include <Driver.d/Domain.h>

extern GeoSource * geoSource;
extern Domain * domain;
extern Communicator * structCom;

Pita::NlDriver::Ptr
nlReversiblePitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor) {
  return Pita::Hts::NlDriver::New(problemDescriptor, geoSource, domain, &domain->solInfo(), structCom);
}
