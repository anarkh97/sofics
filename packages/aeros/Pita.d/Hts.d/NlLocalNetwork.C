#include "NlLocalNetwork.h"

#include "../BasicPropagation.h"
#include "../ReducedSeedWriterTask.h"

namespace Pita { namespace Hts {

NlLocalNetwork::NlLocalNetwork(SliceMapping * mapping,
                               RemoteState::MpiManager * commMgr,
                               NlPropagatorManager * propMgr,
                               CorrectionReductor::Manager * corrRedMgr,
                               CorrectionReconstructor::Manager * corrReconMgr,
                               BasisCondensationManager * condensMgr,
                               ProjectionBuildingFactory * projBuildMgr,
                               JumpConvergenceEvaluator * jumpCvgMgr,
                               NonLinSeedDifferenceEvaluator::Manager * jumpEvalMgr) :
  mapping_(mapping),
  commMgr_(commMgr),
  seedMgr_(Seed::Manager::New()),
  reducedSeedMgr_(ReducedSeed::Manager::New()),
  propMgr_(propMgr),
  jumpMgr_(JumpBuilder::ManagerImpl::New()),
  seedUpMgr_(SeedUpdater::ManagerImpl::New()),
  corrRedMgr_(corrRedMgr),
  corrReconMgr_(corrReconMgr),
  condensMgr_(condensMgr),
  projBuildMgr_(projBuildMgr),
  jumpCvgMgr_(jumpCvgMgr),
  jumpEvalMgr_(jumpEvalMgr)
{
  init();
}

void
NlLocalNetwork::init() {
  for (SliceMapping::SliceIterator it = mapping_->hostedSlice(localCpu()); it; ++it) {
    HalfSliceRank sliceRank = *it;
    addForwardSlice(sliceRank);
    addBackwardSlice(sliceRank);
  }
}

void
NlLocalNetwork::addForwardSlice(HalfSliceRank beginSeedRank) {
  // Tasks executed before convergence check
  {
    addForwardPropagation(beginSeedRank);
    addPropagatedSeedSend(beginSeedRank + HalfSliceCount(1));
  }

  // Preparatory tasks executed after convergence check
  {
    addProjectionBuilding(beginSeedRank);

    addCorrectionReductor(beginSeedRank);
    addCorrectionRecv(beginSeedRank);
    addReducedCorrectionSend(beginSeedRank + FullSliceCount(1));
  
    addSeedUpdater(beginSeedRank);
  }
  
  addMainSeed(beginSeedRank);
}

void
NlLocalNetwork::addBackwardSlice(HalfSliceRank beginSeedRank) {
  // Tasks executed before convergence check
  {
    addBackwardPropagation(beginSeedRank);
    addPropagatedSeedRecv(beginSeedRank);
    addJumpBuilder(beginSeedRank);
  }
  
  // Preparatory tasks executed after convergence check
  HalfSliceRank endSeedRank = beginSeedRank.next();
  
  {
    addCorrectionReconstructor(endSeedRank);
    addReducedCorrectionRecv(endSeedRank);
    addCorrectionSend(endSeedRank);
  
    addSeedUpdater(endSeedRank);
  }
  
  addMainSeed(endSeedRank);
}

void
NlLocalNetwork::addForwardPropagation(HalfSliceRank sliceRank) {
  ActivationInfo info(sliceRank);
  int iterParity = parity(sliceRank);

  finePropagators_[iterParity][info] = forwardHalfSliceNew(sliceRank);
  condensations_[iterParity][info] = condensMgr_->instanceNew(HalfSliceId(sliceRank, FORWARD));
}

void
NlLocalNetwork::addPropagatedSeedSend(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = propagatedSeedSendNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank.previous());
  propagatedSeedSyncs_[iterParity][info] = task;
}

void
NlLocalNetwork::addReducedCorrectionRecv(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank - FullSliceCount(1), HalfSliceCount(1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = reducedCorrectionRecvNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][info] = task;
}

void
NlLocalNetwork::addCorrectionReconstructor(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank.previous(), HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = correctionReconstructorNew(seedRank);

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][info] = task;
}

void
NlLocalNetwork::addCorrectionSend(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = correctionSendNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][info] = task;
}


void
NlLocalNetwork::addBackwardPropagation(HalfSliceRank sliceRank) {
  ActivationInfo info(sliceRank, HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  int iterParity = parity(sliceRank.next());
  finePropagators_[iterParity][info] = backwardHalfSliceNew(sliceRank);
  condensations_[iterParity][info] = condensMgr_->instanceNew(HalfSliceId(sliceRank, BACKWARD));
}


void
NlLocalNetwork::addPropagatedSeedRecv(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = propagatedSeedRecvNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank.previous());
  propagatedSeedSyncs_[iterParity][info] = task;
}

void
NlLocalNetwork::addJumpBuilder(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = jumpBuilderNew(seedRank);

  int iterParity = parity(seedRank.previous());
  jumpBuilders_[iterParity][info] = task;

  jumpCvgMgr_->localJumpIs(seedRank, fullSeedGet(SeedId(SEED_JUMP, seedRank)));
  
  if (jumpEvalMgr()) {
    NonLinSeedDifferenceEvaluator::Ptr jumpEvaluator = jumpEvalMgr()->instanceNew(fullSeedGet(SeedId(SEED_JUMP, seedRank)));
    jumpEvaluator->referenceSeedIs(fullSeedGet(SeedId(LEFT_SEED, seedRank)));
  }
}

void
NlLocalNetwork::addProjectionBuilding(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(1));
  if (!isCurrentlyActive(info)) return;
  
  NamedTask::Ptr task = projBuildMgr_->instanceNew(seedRank);

  int iterParity = parity(seedRank);
  projectionBuilders_[iterParity][info] = task;
}

void
NlLocalNetwork::addCorrectionReductor(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = correctionReductorNew(seedRank);

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][info] = task;
}

void
NlLocalNetwork::addCorrectionRecv(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank.previous(), HalfSliceCount(-1), HalfSliceCount(1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = correctionRecvNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][info] = task;
}

void
NlLocalNetwork::addReducedCorrectionSend(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank.previous(), HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = reducedCorrectionSendNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][info] = task;
}

void
NlLocalNetwork::addSeedUpdater(HalfSliceRank seedRank) {
  ActivationInfo info(seedRank, HalfSliceCount(1), HalfSliceCount(-1));
  if (!isCurrentlyActive(info)) return;

  NamedTask::Ptr task = seedUpdater(seedRank);
  if (!task) {
    task = seedUpdaterNew(seedRank);
  }

  int iterParity = parity(seedRank);
  seedUpdaters_[iterParity][info] = task;
}

void
NlLocalNetwork::addMainSeed(HalfSliceRank seedRank) {
  int iterParity = parity(seedRank);
  SliceRank fullSeedRank(seedRank.value() / 2); 
  seeds_[iterParity][fullSeedRank] = fullSeedGet(SeedId(MAIN_SEED, seedRank));
}


NamedTask::Ptr
NlLocalNetwork::forwardHalfSliceNew(HalfSliceRank sliceRank) {
  const HalfSliceId id(sliceRank, FORWARD);

  DynamPropagator::Ptr propagator = propMgr()->instanceNew(id);
  const Fwk::String name = Fwk::String("Propagate ") + toString(id);
  BasicPropagation::Ptr result = BasicPropagation::New(name, propagator.ptr());

  result->seedIs(fullSeedGet(SeedId(MAIN_SEED, sliceRank)));
  result->propagatedSeedIs(fullSeedGet(SeedId(LEFT_SEED, sliceRank.next())));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::backwardHalfSliceNew(HalfSliceRank sliceRank) {
  const HalfSliceId id(sliceRank, BACKWARD);

  DynamPropagator::Ptr propagator = propMgr()->instanceNew(id);
  const Fwk::String name = Fwk::String("Propagate ") + toString(id);
  BasicPropagation::Ptr result = BasicPropagation::New(name, propagator.ptr());
  
  result->seedIs(fullSeedGet(SeedId(MAIN_SEED, sliceRank.next())));
  result->propagatedSeedIs(fullSeedGet(SeedId(RIGHT_SEED, sliceRank)));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::propagatedSeedSendNew(HalfSliceRank seedRank) {
  assert(hostCpu(seedRank.previous()) == localCpu());
  CpuRank targetCpu(hostCpu(seedRank));
  
  RemoteStateTask::Ptr result = NULL;
  if (targetCpu != localCpu()) {
    RemoteState::SeedReader::Ptr reader = commMgr_->readerNew(fullSeedGet(SeedId(LEFT_SEED, seedRank)), targetCpu);
    String taskName = String("Send Propagated Seed ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, reader.ptr());
  }
  return result;
}

NamedTask::Ptr
NlLocalNetwork::propagatedSeedRecvNew(HalfSliceRank seedRank) {
  assert(hostCpu(seedRank) == localCpu());
  CpuRank previousCpu = hostCpu(seedRank.previous());

  RemoteStateTask::Ptr result = NULL;
  if (previousCpu != localCpu()) {
    RemoteState::SeedWriter::Ptr writer = commMgr_->writerNew(fullSeedGet(SeedId(LEFT_SEED, seedRank)), previousCpu);
    String taskName = String("Receive Propagated Seed ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, writer.ptr());
  }
  return result;
}

NamedTask::Ptr
NlLocalNetwork::jumpBuilderNew(HalfSliceRank seedRank) {
  JumpBuilder::Ptr result = jumpMgr()->instanceNew(toString(seedRank));

  result->predictedSeedIs(fullSeedGet(SeedId(RIGHT_SEED, seedRank)));
  result->actualSeedIs(fullSeedGet(SeedId(LEFT_SEED, seedRank)));
  result->seedJumpIs(fullSeedGet(SeedId(SEED_JUMP, seedRank)));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionReductorNew(HalfSliceRank seedRank) {
  CorrectionReductor::Ptr result = corrRedMgr()->instanceNew(seedRank);

  result->jumpIs(fullSeedGet(SeedId(SEED_JUMP, seedRank)));
  result->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)));
  result->nextCorrectionIs(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank + FullSliceCount(1))));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionReconstructorNew(HalfSliceRank seedRank) {
  CorrectionReconstructor::Ptr result = corrReconMgr()->instanceNew(seedRank);

  result->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)));
  result->correctionComponentsIs(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank)));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionSendNew(HalfSliceRank seedRank) {
  CpuRank nextHeadCpu = hostCpu(seedRank);

  RemoteStateTask::Ptr result = NULL;
  if (nextHeadCpu != localCpu()) {
    RemoteState::SeedReader::Ptr reader = commMgr_->readerNew(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)), nextHeadCpu);
    String taskName = String("Send Correction ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, reader.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionRecvNew(HalfSliceRank seedRank) {
  CpuRank previousTailCpu = hostCpu(seedRank.previous());

  RemoteStateTask::Ptr result = NULL;
  if (previousTailCpu != localCpu()) {
    RemoteState::SeedWriter::Ptr writer = commMgr_->writerNew(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)), previousTailCpu);
    String taskName = String("Recv Correction ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, writer.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::reducedCorrectionSendNew(HalfSliceRank seedRank) {
  CpuRank tailCpu = hostCpu(seedRank.previous());

  RemoteStateTask::Ptr result = NULL;
  if (tailCpu != localCpu()) {
    RemoteState::ReducedSeedReader::Ptr reader = commMgr_->readerNew(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank)), tailCpu);
    String taskName = String("Send Reduced Correction ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, reader.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::reducedCorrectionRecvNew(HalfSliceRank seedRank) {
  CpuRank headCpu = hostCpu(seedRank - FullSliceCount(1));

  RemoteStateTask::Ptr result = NULL;
  if (headCpu != localCpu()) {
    RemoteState::ReducedSeedWriter::Ptr writer = commMgr_->writerNew(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank)), headCpu);
    RemoteState::MpiReducedSeedWriter::Ptr mpiWriter = ptr_cast<RemoteState::MpiReducedSeedWriter>(writer); // Safe cast
    String taskName = String("Recv Reduced Correction ") + toString(seedRank);
    CorrectionReconstructor::Ptr reconstructor = corrReconMgr()->instance(seedRank);
    assert(reconstructor);
    result = ReducedSeedWriterTask<CorrectionReconstructor>::New(taskName, mpiWriter.ptr(), reconstructor.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::seedUpdater(HalfSliceRank seedRank) {
  return seedUpMgr()->instance(toString(seedRank));
}

NamedTask::Ptr
NlLocalNetwork::seedUpdaterNew(HalfSliceRank seedRank) {
  SeedUpdater::Ptr result = seedUpMgr()->instanceNew(toString(seedRank));

  result->propagatedSeedIs(fullSeedGet(SeedId(LEFT_SEED, seedRank)));
  result->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)));
  result->updatedSeedIs(fullSeedGet(SeedId(MAIN_SEED, seedRank)));
    
  seedCorrection_[seedRank] = fullSeedGet(SeedId(SEED_CORRECTION, seedRank));

  return result;
}

Seed *
NlLocalNetwork::fullSeedGet(const SeedId & id) {
  const String name = toString(id);
  Seed::Ptr result = seedMgr_->instance(name);
  if (!result) {
    result = seedMgr_->instanceNew(name);
  }
  return result.ptr();
}

ReducedSeed *
NlLocalNetwork::reducedSeedGet(const SeedId & id) {
  const String name = toString(id);
  ReducedSeed::Ptr result = reducedSeedMgr_->instance(name);
  if (!result) {
    result = reducedSeedMgr_->instanceNew(name);
  }
  return result.ptr();
}

void
NlLocalNetwork::applyConvergenceStatus() {
  // Deactivate correction
  {
    SeedMap::iterator firstActiveCorrection = seedCorrection_.upper_bound(firstActiveSlice());
    for (SeedMap::iterator it = seedCorrection_.begin(); it != firstActiveCorrection; ++it) {
      log() << "Inactivate " << it->second->name() << "\n";
      it->second->statusIs(Seed::INACTIVE);
    }
    seedCorrection_.erase(seedCorrection_.begin(), firstActiveCorrection);

    if (seedUpdater(firstActiveSlice())) {
      Seed::PtrConst leadMainSeed = fullSeedGet(SeedId(MAIN_SEED, firstActiveSlice()));
      Seed::Ptr leadLeftSeed = fullSeedGet(SeedId(LEFT_SEED, firstActiveSlice()));
      if ((leadMainSeed->status() == Seed::ACTIVE || leadMainSeed->status() == Seed::CONVERGED) && leadMainSeed->iteration() > leadLeftSeed->iteration()) {
        // The seed to be updated is more recent than the forward-propagated state.
        // It may happen when the convergence front has moved in a non-staggered fashion,
        // but only happen when the complete active time-domain has converged.
        // The leading seed is used as a 'fake' propagated seed to account for this exception.
        leadLeftSeed->stateIs(leadMainSeed->state());
        leadLeftSeed->iterationIs(leadMainSeed->iteration());
        leadLeftSeed->statusIs(Seed::CONVERGED);
      }
    }
  }

  for (int parity = 0; parity < 2; ++parity) {
    HalfSliceRank start = firstActiveSlice();
    finePropagators_[parity].erase(finePropagators_[parity].begin(), finePropagators_[parity].lower_bound(start));
    propagatedSeedSyncs_[parity].erase(propagatedSeedSyncs_[parity].begin(), propagatedSeedSyncs_[parity].lower_bound(start));
    jumpBuilders_[parity].erase(jumpBuilders_[parity].begin(), jumpBuilders_[parity].lower_bound(start));
    correctionPropagators_[parity].erase(correctionPropagators_[parity].begin(), correctionPropagators_[parity].lower_bound(start));
    seedUpdaters_[parity].erase(seedUpdaters_[parity].begin(), seedUpdaters_[parity].lower_bound(start));

    condensations_[parity].erase(condensations_[parity].begin(), condensations_[parity].lower_bound(start));
    projectionBuilders_[parity].erase(projectionBuilders_[parity].begin(), projectionBuilders_[parity].lower_bound(start));

    SliceRank seedStart((start.value() + parity) / 2);
    seeds_[parity].erase(seeds_[parity].begin(), seeds_[parity].lower_bound(seedStart));
  }
}

NlLocalNetwork::TaskList
NlLocalNetwork::toTaskList(const TaskMap &tm) {
  TaskList result;

  TaskMap::const_iterator it_end = tm.end();
  for (TaskMap::const_iterator it = tm.begin(); it != it_end; ++it) {
    result.push_back(it->second);
  }

  return result;
}

} /* end namespace Hts */ } /* end namespace Pita */
