#ifndef PITA_HTS_LOCALNETWORK_H
#define PITA_HTS_LOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "SliceMapping.h"

#include "../Seed.h"

#include <map>
#include <deque>

#include "LocalNetworkImpl.h"

namespace Pita { namespace Hts {

class LocalNetwork : public Fwk::PtrInterface<LocalNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(LocalNetwork);

  enum Status {
    EMPTY = 0,
    ACTIVE
  };

  Status status() const { return status_; }
  virtual void statusIs(Status s) = 0;

  typedef std::map<HalfSliceRank, NamedTask::Ptr> TaskMap;
  typedef std::deque<NamedTask::Ptr> TaskList;

  typedef std::map<HalfSliceRank, Seed::Ptr> SeedMap;
  typedef std::map<HalfSliceRank, ReducedSeed::Ptr> ReducedSeedMap;
  typedef std::map<SliceRank, Seed::Ptr> MainSeedMap;
 
  // Basic info
  FullSliceCount totalFullSlices() const { return mapping_->totalFullSlices(); }
  CpuCount availableCpus() const { return mapping_->availableCpus(); }
  HalfSliceCount maxWorkload() const { return mapping_->maxWorkload(); }
  CpuRank localCpu() const { return commMgr()->localCpu(); }

  // Current computational state
  HalfSliceRank firstActiveSlice() const { return mapping_->firstActiveSlice(); }
  HalfSliceRank firstInactiveSlice() const { return mapping_->firstInactiveSlice(); }
  HalfSliceCount convergedSlices() const { return mapping_->convergedSlices(); }

  // Change computational state
  void convergedSlicesInc(HalfSliceCount inc) { mapping_->convergedSlicesInc(inc); }

public: 
  Seed::Manager * seedManager() { return fullSeedGetter_.stateMgr(); }
  ReducedSeed::Manager * reducedSeedManager() { return reducedSeedGetter_.stateMgr(); }

protected:
  LocalNetwork(SliceMapping * mapping, RemoteState::Manager * commMgr) :
    status_(EMPTY), mapping_(mapping), commMgr_(commMgr),
    fullSeedGetter_(Seed::Manager::New()), reducedSeedGetter_(ReducedSeed::Manager::New())
  {}
  
  LocalNetworkImpl::SeedGetter<DynamState> fullSeedGetter() { return fullSeedGetter_; }
  LocalNetworkImpl::SeedGetter<Vector> reducedSeedGetter() { return reducedSeedGetter_; }

  void setStatus(Status s) { status_ = s; }

  SliceMapping::SliceIterator hostedSlice(CpuRank cpu) const { return mapping_->hostedSlice(cpu); }
  CpuRank hostCpu(HalfSliceRank slice) const { return mapping_->hostCpu(slice); }

  RemoteState::Manager * commMgr() { return commMgr_.ptr(); }
  const RemoteState::Manager * commMgr() const { return commMgr_.ptr(); }

  SharedState<DynamState> * fullSeedGet(const SeedId & id) { return fullSeedGetter_.operator()(id); }
  SharedState<Vector> * reducedSeedGet(const SeedId & id) { return reducedSeedGetter_.operator()(id); }

  static int parity(HalfSliceRank sliceRank) { return sliceRank.value() % 2; } 
  int activeParity() const { return parity(firstActiveSlice()); } 

private:
  Status status_;

  SliceMapping::Ptr mapping_;
  RemoteState::Manager::Ptr commMgr_;
  
  LocalNetworkImpl::SeedGetter<DynamState> fullSeedGetter_;
  LocalNetworkImpl::SeedGetter<Vector> reducedSeedGetter_;

  DISALLOW_COPY_AND_ASSIGN(LocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARLOCALNETWORK_H */
