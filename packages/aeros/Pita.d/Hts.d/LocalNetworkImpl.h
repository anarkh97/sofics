#ifndef PITA_HTS_LOCALNETWORKIMPL_H
#define PITA_HTS_LOCALNETWORKIMPL_H

#include "Fwk.h"
#include "Types.h"
#include "../SharedState.h"

#include "../NamedTask.h"

#include "../JumpBuilder.h"

#include "../SeedUpdater.h"
#include "../UpdatedSeedAssembler.h"

#include "../CorrectionPropagator.h"

#include "../RemoteStateTask.h"

#include "SliceMapping.h"

#include <map>
#include <deque>
#include <cassert>

namespace Pita { namespace Hts {

namespace LocalNetworkImpl {

/* Seed getter functors */

template <typename S>
class SeedGetter {
public:
  typedef S StateType;
  typedef typename SharedState<S>::Manager StateManager;

  explicit SeedGetter(const Fwk::Ptr<StateManager> & stateMgr) :
    stateMgr_(stateMgr)
  {}

  SharedState<S> * operator()(const SeedId & id);
  SharedState<S> * operator()(const SeedId & id) const;

  StateManager * stateMgr() { return stateMgr_.ptr(); }
  const StateManager * stateMgr() const { return stateMgr_.ptr(); }

private:
  Fwk::Ptr<StateManager> stateMgr_;
};

template <typename S>
inline
SharedState<S> *
SeedGetter<S>::operator()(const SeedId & id) const {
  String seedName = toString(id);
  return stateMgr()->instance(seedName);
}

template <typename S>
inline
SharedState<S> *
SeedGetter<S>::operator()(const SeedId & id) {
  String seedName = toString(id);
  SharedState<S> * result = stateMgr()->instance(seedName);
  return (result != NULL) ? result : stateMgr()->instanceNew(seedName);
}

typedef SeedGetter<DynamState> FullSeedGetter;
typedef SeedGetter<Vector> ReducedSeedGetter;

/* Base factory interface */

class NamedTaskFactory : public Fwk::PtrInterface<NamedTaskFactory> {
public:
  EXPORT_PTRINTERFACE_TYPES(NamedTaskFactory);

  virtual NamedTask * taskGet(HalfSliceRank sliceRank) = 0;
};

/* Manager-based factory implementation */

template <typename T>
class Decorator {
public:
  void operator()(HalfSliceRank sliceRank, T * task) {}
};

template <typename T>
class GenNamedTaskFactory : public NamedTaskFactory {
public:
  EXPORT_PTRINTERFACE_TYPES(GenNamedTaskFactory);
  typedef typename T::Manager TaskManager;

  virtual T * taskGet(HalfSliceRank sliceRank);

  GenNamedTaskFactory(TaskManager * taskManager, Decorator<T> decorator) :
    taskManager_(taskManager),
    decorator_(decorator)
  {}

protected:
  void decorate(HalfSliceRank sliceRank, T * task) { decorator_(sliceRank, task); }

private:
  Fwk::Ptr<TaskManager> taskManager_;
  Decorator<T> decorator_;
};

template <typename T>
T *
GenNamedTaskFactory<T>::taskGet(HalfSliceRank sliceRank) {
  String taskName = toString(sliceRank);
  T * task = taskManager_->instance(taskName);
  if (!task) {
    task = taskManager_->instanceNew(taskName);
    decorate(sliceRank, task);
  }
  return task;
}

/* Specializations */

// JumpBuilder
template <>
class Decorator<JumpBuilder> {
public:
  explicit Decorator(FullSeedGetter fsg) :
    fullSeed_(fsg)
  {}

  void operator()(HalfSliceRank sliceRank, JumpBuilder * task) {
    task->predictedSeedIs(fullSeed_(SeedId(RIGHT_SEED, sliceRank)));
    task->actualSeedIs(fullSeed_(SeedId(LEFT_SEED, sliceRank)));
    task->seedJumpIs(fullSeed_(SeedId(SEED_JUMP, sliceRank)));
  }

private:
  FullSeedGetter fullSeed_;
};

typedef GenNamedTaskFactory<JumpBuilder> JumpBuilderFactory;

inline
JumpBuilderFactory::Ptr
JumpBuilderFactoryNew(JumpBuilder::Manager * mgr, FullSeedGetter fsg) {
  return new JumpBuilderFactory(mgr, Decorator<JumpBuilder>(fsg));
}

// SeedUpdater
template <>
class Decorator<SeedUpdater> {
public:
  explicit Decorator(FullSeedGetter fsg) :
    fullSeed_(fsg)
  {}

  void operator()(HalfSliceRank sliceRank, SeedUpdater * task) {
    task->correctionIs(fullSeed_(SeedId(SEED_CORRECTION, sliceRank)));
    task->propagatedSeedIs(fullSeed_(SeedId(LEFT_SEED, sliceRank)));
    task->updatedSeedIs(fullSeed_(SeedId(MAIN_SEED, sliceRank)));
  }

private:
  FullSeedGetter fullSeed_;
};

typedef GenNamedTaskFactory<SeedUpdater> SeedUpdaterFactory;

inline
SeedUpdaterFactory::Ptr
SeedUpdaterFactoryNew(SeedUpdater::Manager * mgr, FullSeedGetter fsg) {
  return new SeedUpdaterFactory(mgr, Decorator<SeedUpdater>(fsg));
}


// CorrectionPropagator<S>
template <typename S>
class Decorator<CorrectionPropagator<S> > {
public:
  explicit Decorator(SeedGetter<S> sg) :
    seed_(sg)
  {}

  void operator()(HalfSliceRank sliceRank, CorrectionPropagator<S> * task) {
    task->jumpIs(seed_(SeedId(SEED_JUMP, sliceRank)));
    task->correctionIs(seed_(SeedId(SEED_CORRECTION, sliceRank)));
    task->nextCorrectionIs(seed_(SeedId(SEED_CORRECTION, sliceRank + HalfSliceCount(2))));
  }

private:
    SeedGetter<S> seed_;
};

template <typename S>
inline
typename GenNamedTaskFactory<CorrectionPropagator<S> >::Ptr
CorrectionPropagatorFactoryNew(typename CorrectionPropagator<S>::Manager * mgr, SeedGetter<S> sg) {
  return new GenNamedTaskFactory<CorrectionPropagator<S> >(mgr, sg);
}

// UpdatedSeedAssembler
template <>
class Decorator<UpdatedSeedAssembler> {
public:
  Decorator(FullSeedGetter fsg, ReducedSeedGetter rsg) :
    delegate_(fsg), reducedSeed_(rsg)
  {}

  void operator()(HalfSliceRank sliceRank, UpdatedSeedAssembler * task) {
    delegate_(sliceRank, task),
    task->correctionComponentsIs(reducedSeed_(SeedId(SEED_CORRECTION, sliceRank)));
  }

private:
  Decorator<SeedUpdater> delegate_;
  ReducedSeedGetter reducedSeed_;
};

typedef GenNamedTaskFactory<UpdatedSeedAssembler> SeedAssemblerFactory;

inline
SeedAssemblerFactory::Ptr
SeedAssemblerFactoryNew(UpdatedSeedAssembler::Manager * mgr, FullSeedGetter fsg, ReducedSeedGetter rsg) {
  return new SeedAssemblerFactory(mgr, Decorator<UpdatedSeedAssembler>(fsg, rsg));
}

// Helper class
template <typename S>
class CorrectionPropagatorBuilder {
public:
  CorrectionPropagatorBuilder(const SliceMapping * mapping,
                              typename CorrectionPropagator<S>::Manager * cpMgr,
                              RemoteState::Manager * commMgr,
                              LocalNetworkImpl::SeedGetter<S> seedGetter) :
    mapping_(mapping),
    cpMgr_(cpMgr),
    commMgr_(commMgr),
    seedGetter_(seedGetter)
  {}

  typedef std::map<HalfSliceRank, NamedTask::Ptr> TaskMap;

  void buildLocalPropagator(HalfSliceRank sliceRank, TaskMap & propagation);
  void buildCorrectionRecv(HalfSliceRank sliceRank, TaskMap & propagation);
  void buildCorrectionSend(HalfSliceRank sliceRank, TaskMap & propagation);
  void buildCorrectionSyncRecv(HalfSliceRank sliceRank, TaskMap & propagation);
  void buildCorrectionSyncSend(HalfSliceRank sliceRank, TaskMap & propagation);

  const SliceMapping * mapping() const { return mapping_; }
  CpuRank localCpu() const { return commMgr_->localCpu(); }
  typename CorrectionPropagator<S>::Manager * cpMgr() const { return cpMgr_; }
  RemoteState::Manager * commMgr() const { return commMgr_; }

protected:
  SharedState<S> * getSeed(const SeedId & seedId) {
    return seedGetter_(seedId);
  }

private:
  const SliceMapping * mapping_;
  typename CorrectionPropagator<S>::Manager * cpMgr_;
  RemoteState::Manager * commMgr_;
  LocalNetworkImpl::SeedGetter<S> seedGetter_;
};

template <typename S>
void
CorrectionPropagatorBuilder<S>::buildLocalPropagator(
    HalfSliceRank sliceRank,
    CorrectionPropagatorBuilder<S>::TaskMap & propagation) {
  
  HalfSliceRank nextSliceRank = sliceRank + HalfSliceCount(1);
  HalfSliceRank nextFullSliceRank = sliceRank + HalfSliceCount(2);
  
  if (nextSliceRank < mapping()->firstInactiveSlice()) {
    typename CorrectionPropagator<S>::Ptr fullSlice = cpMgr()->instanceNew(toString(sliceRank));
    fullSlice->jumpIs(getSeed(SeedId(SEED_JUMP, sliceRank)));
    fullSlice->correctionIs(getSeed(SeedId(SEED_CORRECTION, sliceRank)));
    fullSlice->nextCorrectionIs(getSeed(SeedId(SEED_CORRECTION, nextFullSliceRank)));
    propagation.insert(std::make_pair(sliceRank, fullSlice));
  }
}

template <typename S>
void
CorrectionPropagatorBuilder<S>::buildCorrectionRecv(
    HalfSliceRank sliceRank,
    CorrectionPropagatorBuilder<S>::TaskMap & propagation) {
  
  HalfSliceRank previousFullSliceRank = sliceRank - HalfSliceCount(2);
  CpuRank previousFullSliceCpu(mapping()->hostCpu(previousFullSliceRank));
  
  if (previousFullSliceCpu != CpuRank(-1) && previousFullSliceCpu != localCpu()) {
    typename RemoteState::Writer<S>::Ptr writer = commMgr()->writerNew(getSeed(SeedId(SEED_CORRECTION, sliceRank)), previousFullSliceCpu);
    if (writer) {
      String taskName = String("Receive Correction ") + toString(sliceRank);
      RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, writer.ptr());
      propagation.insert(std::make_pair(sliceRank.previous(), task));
    }
  }
}

template <typename S>
void
CorrectionPropagatorBuilder<S>::buildCorrectionSend(
    HalfSliceRank sliceRank,
    CorrectionPropagatorBuilder<S>::TaskMap & propagation) {
 
  CpuRank sliceCpu = mapping()->hostCpu(sliceRank);
  
  if (sliceCpu != CpuRank(-1) && sliceCpu != localCpu()) {
    typename RemoteState::Reader<S>::Ptr reader = commMgr()->readerNew(getSeed(SeedId(SEED_CORRECTION, sliceRank)), sliceCpu);
    if (reader) {
      String taskName = String("Send Correction ") + toString(sliceRank);
      RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, reader.ptr());
      propagation.insert(std::make_pair(sliceRank.previous(), task));
    }
  }
}

template <typename S>
void
CorrectionPropagatorBuilder<S>::buildCorrectionSyncRecv(
    HalfSliceRank sliceRank,
    CorrectionPropagatorBuilder<S>::TaskMap & synchronization) {

  assert(localCpu() == mapping()->hostCpu(sliceRank.previous()));

  CpuRank otherTargetCpu(mapping()->hostCpu(sliceRank));
  if (localCpu() == otherTargetCpu) // No need for optimization
    return;

  HalfSliceRank originSliceRank = sliceRank - HalfSliceCount(2);
  CpuRank originCpu(mapping()->hostCpu(originSliceRank));

  if (originCpu != CpuRank(-1) && originCpu != localCpu()) {
    typename RemoteState::Writer<S>::Ptr correctionWriter = commMgr()->writerNew(getSeed(SeedId(SEED_CORRECTION, sliceRank)), originCpu);
    if (correctionWriter) {
      String taskName = String("Receive Correction ") + toString(sliceRank);
      RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, correctionWriter.ptr());
      synchronization.insert(synchronization.end(), std::make_pair(sliceRank - HalfSliceCount(2), task));
    }
  }
}

template <typename S>
void
CorrectionPropagatorBuilder<S>::buildCorrectionSyncSend(
    HalfSliceRank sliceRank,
    CorrectionPropagatorBuilder<S>::TaskMap & synchronization) {

  assert(localCpu() == mapping()->hostCpu(sliceRank - HalfSliceCount(2)));

  CpuRank targetCpu(mapping()->hostCpu(sliceRank.previous()));
  
  CpuRank otherTargetCpu(mapping()->hostCpu(sliceRank));
  if (targetCpu == otherTargetCpu) // No need for optimization 
    return;

  if (targetCpu != CpuRank(-1) && targetCpu != localCpu()) {
    typename RemoteState::Reader<S>::Ptr correctionReader = commMgr()->readerNew(getSeed(SeedId(SEED_CORRECTION, sliceRank)), targetCpu);
    if (correctionReader) {
      String taskName = String("Send Correction ") + toString(sliceRank);
      RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, correctionReader.ptr());
      synchronization.insert(synchronization.end(), std::make_pair(sliceRank + HalfSliceCount(2), task));
    }
  }
}

} /* end namespace LocalNetworkImpl */

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LOCALNETWORKIMPL_H */
