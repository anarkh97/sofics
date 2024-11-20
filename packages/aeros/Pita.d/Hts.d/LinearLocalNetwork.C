#include "LinearLocalNetwork.h"

#include "../IncrementalPropagation.h"

#include <algorithm>
#include <iterator>

#include <cassert>

namespace Pita { namespace Hts {

/* LinearLocalNetwork Constructor */

LinearLocalNetwork::LinearLocalNetwork(SliceMapping * mapping,
                                       AffinePropagatorManager * propMgr,
                                       ReducedCorrectionManager * redCorrMgr,
                                       JumpConvergenceEvaluator * jumpCvgMgr,
                                       RemoteState::Manager * commMgr,
                                       LinSeedDifferenceEvaluator::Manager * jumpErrorMgr) :
  LocalNetwork(mapping, commMgr),
  jumpBuildMgr_(JumpBuilder::ManagerImpl::New()),
  jumpCvgMgr_(jumpCvgMgr),
  propMgr_(propMgr),
  redCorrMgr_(redCorrMgr),
  jumpErrorMgr_(jumpErrorMgr),
  fullCorrectionBuilder_(mapping, redCorrMgr->fcpMgr(), commMgr, fullSeedGetter()),
  reducedCorrectionBuilder_(mapping, redCorrMgr->rcpMgr(), commMgr, reducedSeedGetter())
{
  init();
  setStatus(ACTIVE);
}

/* Task building - Helper classes */

void
LinearLocalNetwork::init() {
  for (SliceMapping::SliceIterator it = hostedSlice(localCpu()); it; ++it) {
    HalfSliceRank primalSliceRank = *it;
    HalfSliceRank dualSliceRank = primalSliceRank.next();
  
    buildForwardPropagation(primalSliceRank);
    if (primalSliceRank != HalfSliceRank(0)) {
      buildBackwardPropagation(dualSliceRank);
      buildPrimalCorrectionNetwork(primalSliceRank);
    }

    buildDualCorrectionNetwork(dualSliceRank);
  }
}

void
LinearLocalNetwork::buildForwardPropagation(HalfSliceRank sliceRank) {
  const HalfSliceRank nextSliceRank = sliceRank.next();
  const HalfSliceId id(sliceRank, FORWARD);

  AffineDynamPropagator::Ptr propagator = propMgr_->instanceNew(id);
  IncrementalPropagation::Ptr forwardSlice = IncrementalPropagation::New(Fwk::String("Propagate ") + toString(id), propagator.ptr());
  
  Seed::Ptr mainSeed = fullSeedGet(SeedId(MAIN_SEED, sliceRank));
  forwardSlice->seedIs(mainSeed.ptr());
  forwardSlice->propagatedSeedIs(fullSeedGet(SeedId(LEFT_SEED, nextSliceRank)));
  halfTimeSlice_[sliceRank.value() % 2].insert(std::make_pair(sliceRank, forwardSlice));

  mainSeed_[sliceRank] = mainSeed;
}

void
LinearLocalNetwork::buildBackwardPropagation(HalfSliceRank sliceRank) {
  const HalfSliceRank previousSliceRank = sliceRank.previous();
  const HalfSliceId id(previousSliceRank, BACKWARD);
  
  AffineDynamPropagator::Ptr propagator = propMgr_->instanceNew(id);
  IncrementalPropagation::Ptr backwardSlice = IncrementalPropagation::New(Fwk::String("Propagate ") + toString(id), propagator.ptr());
  
  Seed::Ptr mainSeed = fullSeedGet(SeedId(MAIN_SEED, sliceRank));
  backwardSlice->seedIs(mainSeed.ptr());
  backwardSlice->propagatedSeedIs(fullSeedGet(SeedId(RIGHT_SEED, previousSliceRank)));
  halfTimeSlice_[sliceRank.value() % 2].insert(std::make_pair(sliceRank.previous(), backwardSlice));

  mainSeed_[sliceRank] = mainSeed;
}

void
LinearLocalNetwork::buildPrimalCorrectionNetwork(HalfSliceRank sliceRank) {
  // Parallel
  buildPropagatedSeedRecv(sliceRank);
  buildJumpAssembly(sliceRank);
  if (jumpErrorMgr_) {
    buildJumpEstimator(sliceRank);
  }
  buildJumpProjection(sliceRank);

  // Sequential
  buildCorrectionPropagator(sliceRank);

  // Parallel
  buildCorrectionSynchronizationSend(sliceRank + HalfSliceCount(2));
  buildSeedUpdater(sliceRank);
}

void
LinearLocalNetwork::buildDualCorrectionNetwork(HalfSliceRank sliceRank) {
  // Entirely parallel (provides support for the primal network)
  if (sliceRank < firstInactiveSlice()) {
    buildPropagatedSeedSend(sliceRank);
  }

  buildCorrectionSynchronizationRecv(sliceRank);

  buildSeedUpdater(sliceRank);
}

void
LinearLocalNetwork::buildPropagatedSeedRecv(HalfSliceRank sliceRank) {
  HalfSliceRank previousSliceRank = sliceRank.previous();
  CpuRank previousCpu = hostCpu(previousSliceRank);

  assert(previousCpu != CpuRank(-1));

  if (previousCpu != localCpu()) {
    RemoteState::SeedWriter::Ptr writer = commMgr()->writerNew(fullSeedGet(SeedId(LEFT_SEED, sliceRank)), previousCpu);
    if (writer) {
      String taskName = String("Receive Propagated Seed ") + toString(sliceRank);
      RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, writer.ptr());
      leftSeedSync_[sliceRank.previous().value() % 2].insert(std::make_pair(sliceRank, task));
    }
  }
}

void
LinearLocalNetwork::buildPropagatedSeedSend(HalfSliceRank sliceRank) {
  if (sliceRank < firstInactiveSlice()) {
    assert(hostCpu(sliceRank.previous()) == localCpu());

    CpuRank targetCpu(hostCpu(sliceRank));

    assert(targetCpu != CpuRank(-1));
    
    if (targetCpu != localCpu()) {
      RemoteState::SeedReader::Ptr reader = commMgr()->readerNew(fullSeedGet(SeedId(LEFT_SEED, sliceRank)), targetCpu);
      if (reader) {
        String taskName = String("Send Propagated Seed ") + toString(sliceRank);
        RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, reader.ptr());
        leftSeedSync_[sliceRank.previous().value() % 2].insert(std::make_pair(sliceRank, task));
      }
    }
  }
}

void
LinearLocalNetwork::buildJumpAssembly(HalfSliceRank sliceRank) {
  String taskName = toString(sliceRank);
  JumpBuilder::Ptr task = jumpBuildMgr_->instanceNew(taskName);

  task->predictedSeedIs(fullSeedGet(SeedId(RIGHT_SEED, sliceRank)));
  task->actualSeedIs(fullSeedGet(SeedId(LEFT_SEED, sliceRank)));
  task->seedJumpIs(fullSeedGet(SeedId(SEED_JUMP, sliceRank)));

  jumpCvgMgr_->localJumpIs(sliceRank, fullSeedGet(SeedId(SEED_JUMP, sliceRank)));

  jumpAssembler_[sliceRank.previous().value() % 2].insert(std::make_pair(sliceRank, task));
}

void
LinearLocalNetwork::buildJumpProjection(HalfSliceRank sliceRank) {
  String taskName = toString(sliceRank);
  JumpProjection::Ptr task = redCorrMgr_->jumpProjMgr()->instanceNew(taskName);

  task->seedJumpIs(fullSeedGet(SeedId(SEED_JUMP, sliceRank)));
  task->reducedSeedJumpIs(reducedSeedGet(SeedId(SEED_JUMP, sliceRank)));

  jumpProjector_[sliceRank.value() % 2].insert(std::make_pair(sliceRank, task));
}

void
LinearLocalNetwork::buildJumpEstimator(HalfSliceRank sliceRank) {
  LinSeedDifferenceEvaluator::Ptr jumpErrorEvaluator = jumpErrorMgr_->instanceNew(fullSeedGet(SeedId(SEED_JUMP, sliceRank)));
  jumpErrorEvaluator->referenceSeedIs(fullSeedGet(SeedId(LEFT_SEED, sliceRank)));
}

void
LinearLocalNetwork::buildSeedUpdater(HalfSliceRank sliceRank) {
  String taskName = toString(sliceRank);
  UpdatedSeedAssembler::Ptr task = redCorrMgr_->usaMgr()->instance(taskName);

  if (!task) {
    task = redCorrMgr_->usaMgr()->instanceNew(taskName);

    task->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, sliceRank)));
    task->propagatedSeedIs(fullSeedGet(SeedId(LEFT_SEED, sliceRank)));
    task->updatedSeedIs(fullSeedGet(SeedId(MAIN_SEED, sliceRank)));
    task->correctionComponentsIs(reducedSeedGet(SeedId(SEED_CORRECTION, sliceRank)));
   
    reducedSeedCorrection_[sliceRank] = reducedSeedGet(SeedId(SEED_CORRECTION, sliceRank));
    seedCorrection_[sliceRank] = fullSeedGet(SeedId(SEED_CORRECTION, sliceRank));

    seedAssembler_[sliceRank.value() % 2].insert(std::make_pair(sliceRank, task));
  }
}

void 
LinearLocalNetwork::buildCorrectionPropagator(HalfSliceRank sliceRank) {
  buildReducedCorrectionPropagator(sliceRank);
  if (redCorrMgr_->fcpMgr()) {
    buildFullCorrectionPropagator(sliceRank);
  }
}

void
LinearLocalNetwork::buildCorrectionSynchronizationRecv(HalfSliceRank sliceRank) {
  reducedCorrectionBuilder_.buildCorrectionSyncRecv(sliceRank, correctionSync_[sliceRank.value() % 2]);
  if (redCorrMgr_->fcpMgr()) {
    fullCorrectionBuilder_.buildCorrectionSyncRecv(sliceRank, alternateCorrectionSync_[sliceRank.value() % 2]);
  }
}

void
LinearLocalNetwork::buildCorrectionSynchronizationSend(HalfSliceRank sliceRank) {
  reducedCorrectionBuilder_.buildCorrectionSyncSend(sliceRank, correctionSync_[sliceRank.value() % 2]);
  if (redCorrMgr_->fcpMgr()) {
    fullCorrectionBuilder_.buildCorrectionSyncSend(sliceRank, alternateCorrectionSync_[sliceRank.value() % 2]);
  }
}

void
LinearLocalNetwork::buildReducedCorrectionPropagator(HalfSliceRank sliceRank) {
  TaskMap & fullTimeSlice = correctionPropagator_[sliceRank.value() % 2];
  reducedCorrectionBuilder_.buildCorrectionRecv(sliceRank, fullTimeSlice);
  reducedCorrectionBuilder_.buildLocalPropagator(sliceRank, fullTimeSlice);
  reducedCorrectionBuilder_.buildCorrectionSend(sliceRank + HalfSliceCount(2), fullTimeSlice);
}

void
LinearLocalNetwork::buildFullCorrectionPropagator(HalfSliceRank sliceRank) {
  TaskMap & fullTimeSlice = alternateCorrectionPropagator_[sliceRank.value() % 2];
  fullCorrectionBuilder_.buildCorrectionRecv(sliceRank, fullTimeSlice);
  fullCorrectionBuilder_.buildLocalPropagator(sliceRank, fullTimeSlice);
  fullCorrectionBuilder_.buildCorrectionSend(sliceRank + HalfSliceCount(2), fullTimeSlice);
}

/* Active tasks: LinearLocalNetwork */

void
LinearLocalNetwork::applyConvergenceStatus() {
  mainSeed_.erase(mainSeed_.begin(), mainSeed_.lower_bound(firstActiveSlice()));

  // Deactivate full correction
  {
    SeedMap::iterator firstActiveCorrection = seedCorrection_.upper_bound(firstActiveSlice());
    for (SeedMap::iterator it = seedCorrection_.begin(); it != firstActiveCorrection; ++it) {
      log() << "Inactivate full " << it->second->name() << "\n";
      it->second->statusIs(Seed::INACTIVE);
    }
    seedCorrection_.erase(seedCorrection_.begin(), firstActiveCorrection);
    
    if (redCorrMgr_->usaMgr()->instance(toString(firstActiveSlice()))) {
      Seed::PtrConst leadMainSeed = fullSeedGet(SeedId(MAIN_SEED, firstActiveSlice()));
      Seed::Ptr leadLeftSeed = fullSeedGet(SeedId(LEFT_SEED, firstActiveSlice()));
      if ((leadMainSeed->status() == Seed::ACTIVE || leadMainSeed->status() == Seed::CONVERGED) && leadMainSeed->iteration() > leadLeftSeed->iteration()) {
        // Convergence front has moved in a non-staggered fashion
        // It should only happen when the complete active time-domain has converged
        // The leading seed is used to create a 'fake' propagated seed to account for this exception
        leadLeftSeed->stateIs(leadMainSeed->state());
        leadLeftSeed->iterationIs(leadMainSeed->iteration());
        leadLeftSeed->statusIs(Seed::CONVERGED);
      }
    }
  }
  
  // Deactivate reduced correction
  {
    ReducedSeedMap::iterator firstActiveCorrection = reducedSeedCorrection_.upper_bound(firstActiveSlice());
    for (ReducedSeedMap::iterator it = reducedSeedCorrection_.begin(); it != firstActiveCorrection; ++it) {
      log() << "Inactivate reduced " << it->second->name() << "\n";
      it->second->statusIs(Seed::INACTIVE);
    }
    reducedSeedCorrection_.erase(reducedSeedCorrection_.begin(), firstActiveCorrection);
  }

  eraseInactive(halfTimeSlice_[0]);
  eraseInactive(halfTimeSlice_[1]);
  
  eraseInactive(leftSeedSync_[0]);
  eraseInactive(leftSeedSync_[1]);
  eraseInactive(jumpAssembler_[0]);
  eraseInactive(jumpAssembler_[1]);
  eraseInactive(jumpProjector_[0]);
  eraseInactive(jumpProjector_[1]);
  
  eraseInactive(correctionPropagator_[0]);
  eraseInactive(correctionPropagator_[1]);
  eraseInactive(correctionSync_[0]);
  eraseInactive(correctionSync_[1]);
  
  eraseInactive(alternateCorrectionPropagator_[0]);
  eraseInactive(alternateCorrectionPropagator_[1]);
  eraseInactive(alternateCorrectionSync_[0]);
  eraseInactive(alternateCorrectionSync_[1]);

  eraseInactive(seedAssembler_[0]);
  eraseInactive(seedAssembler_[1]);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::halfTimeSlices() const {
  TaskList result = getAll(halfTimeSlice_[0]);
  TaskList other = getAll(halfTimeSlice_[1]);
  result.insert(result.end(), other.begin(), other.end());
  return result;
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeHalfTimeSlices() const {
  return getActive(halfTimeSlice_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeFullTimeSlices() const {
  return getActive(correctionPropagator_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeCoarseTimeSlices() const {
  return getActive(alternateCorrectionPropagator_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeJumpAssemblers() const {
  return getActive(jumpAssembler_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeJumpProjectors() const {
  return getActive(jumpProjector_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeLeftSeedSyncs() const {
  return getActive(leftSeedSync_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeCorrectionSyncs() const {
  return getActive(correctionSync_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeFullCorrectionSyncs() const {
  return getActive(alternateCorrectionSync_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::activeSeedAssemblers() const {
  return getActive(seedAssembler_);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::getActive(const TaskMap task[2]) const {
  return getAll(task[firstActiveSlice().value() % 2]);
}

LinearLocalNetwork::TaskList
LinearLocalNetwork::getAll(const TaskMap & task) {
  TaskList result;

  TaskMap::const_iterator it_end = task.end();
  for (TaskMap::const_iterator it = task.begin(); it != it_end; ++it) {
    result.push_back(it->second);
  }

  return result;
}

LinearLocalNetwork::SeedMap
LinearLocalNetwork::mainSeeds() const {
  return SeedMap(mainSeed_);
}

LinearLocalNetwork::MainSeedMap
LinearLocalNetwork::activeMainSeeds() const {
  MainSeedMap msp;

  SeedMap::const_iterator it_end = mainSeed_.end();
  for (SeedMap::const_iterator it = mainSeed_.begin(); it != it_end; ++it) {
    if ((it->first.value() - firstActiveSlice().value()) % 2 == 0) {
      msp.insert(msp.end(), std::make_pair(SliceRank(it->first.value() / 2), it->second));
    }
  }

  return msp;
}

} /* end namespace Hts */ } /* end namespace Pita */
