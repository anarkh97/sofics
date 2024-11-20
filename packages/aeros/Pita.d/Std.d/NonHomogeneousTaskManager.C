#include "NonHomogeneousTaskManager.h"

#include "../SeedInitializer.h"

namespace Pita { namespace Std {

NonHomogeneousTaskManager::NonHomogeneousTaskManager(SliceMapping * mapping,
                                                     RemoteState::Manager * commMgr,
                                                     LinearPropagatorManager * propagatorMgr,
                                                     LinearProjectionNetwork * projectionMgr,
                                                     JumpConvergenceEvaluator * jumpCvgEval,
                                                     LinSeedDifferenceEvaluator::Manager * jumpOutMgr,
                                                     SeedInitializer * seedInit,
                                                     CorrectionPropagator<DynamState>::Manager * fullCorrPropMgr) :
  LinearTaskManager<NonHomogeneousTaskManager>(IterationRank(-1), mapping, commMgr, propagatorMgr, projectionMgr, jumpCvgEval, jumpOutMgr, seedInit),
  fullCorrPropMgr_(fullCorrPropMgr)
{
  setContinuation(&NonHomogeneousTaskManager::schedulePrecomputation);
}

void
NonHomogeneousTaskManager::iterationInc() {
  if (iteration() == IterationRank(-1)) {
    mapping()->convergedSlicesInc();
    applyConvergence();
  
    scheduleCoarseCorrectionPropagation();
    setIteration(iteration().next());
  } else {
    LinearTaskManager<NonHomogeneousTaskManager>::iterationInc();
  }
}

void
NonHomogeneousTaskManager::schedulePrecomputation() {
  setRecurrentPhase("Affine Term Precomputation", &SliceTasks::finePropagation);
  setContinuation(&NonHomogeneousTaskManager::schedulePropagatedSeedSynchronization);
}

void
NonHomogeneousTaskManager::scheduleCoarseCorrectionPropagation() {
  TaskList correctionPropagation;
 
  for (SliceMapping::SliceIterator sl_it = mapping()->hostedSlice(localCpu()); sl_it; ++sl_it) {
    SliceRank slice = *sl_it;
    SliceRank previousSlice = slice.previous();
    SliceRank nextSlice = slice.next();

    if (slice < mapping()->firstActiveSlice()) {
      continue;
    }

    if (slice >= mapping()->firstInactiveSlice()) {
      break;
    }

    // Receive previous  
    if (previousSlice >= mapping()->firstActiveSlice() && mapping()->hostCpu(previousSlice) != localCpu()) {
      CpuRank peer = mapping()->hostCpu(previousSlice);
      
      Seed::Ptr correction = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));
      RemoteState::Writer<DynamState>::Ptr activity = commMgr_->writerNew(correction.ptr(), peer);
      RemoteStateTask::Ptr task = RemoteStateTask::New("Receive " + correction->name(), activity.ptr());
      
      correctionPropagation.push_back(task);
    } else {
      seedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)))->statusIs(Seed::INACTIVE);
    }

    // Local propagation
    {
      Seed::Ptr jump = seedMgr_->instance(toString(SeedId(SEED_JUMP, slice)));
      Seed::Ptr correction = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));
      Seed::Ptr nextCorrection = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, nextSlice)));
      if (!nextCorrection) {
        nextCorrection = seedMgr_->instanceNew(toString(SeedId(SEED_CORRECTION, nextSlice)));
      }

      CorrectionPropagator<DynamState>::Ptr task = fullCorrPropMgr_->instanceNew(toString(slice));

      task->jumpIs(jump.ptr());
      task->correctionIs(correction.ptr());
      task->nextCorrectionIs(nextCorrection.ptr());

      correctionPropagation.push_back(task);
    }

    // Send next
    if (nextSlice < mapping()->firstInactiveSlice() && mapping()->hostCpu(nextSlice) != localCpu()) {
      CpuRank peer = mapping()->hostCpu(nextSlice);

      Seed::Ptr correction = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, nextSlice)));
      RemoteState::Reader<DynamState>::Ptr activity = commMgr_->readerNew(correction.ptr(), peer);
      RemoteStateTask::Ptr task = RemoteStateTask::New("Send " + correction->name(), activity.ptr());
      
      correctionPropagation.push_back(task);
    }
  }
  
  setPhase(phaseNew("Coarse Correction Propagation", correctionPropagation));
  
  setContinuation(&NonHomogeneousTaskManager::scheduleSeedUpdate);
}

} /* end namespace Std */ } /* end namespace Pita */
