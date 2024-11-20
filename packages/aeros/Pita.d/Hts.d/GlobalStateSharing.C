#include "GlobalStateSharing.h"

#include "../DynamStatePlainBasis.h"
#include "../DynamStateBasisWrapper.h"
#include <Comm.d/Communicator.h>

#include <numeric>

namespace Pita { namespace Hts {

GlobalStateSharing::GlobalStateSharing(Communicator * timeComm, size_t vectorSize, Strategy strategy) :
  NamedTask("Exchange Global Projection Basis"),
  timeComm_(timeComm),
  vectorSize_(vectorSize),
  strategy_(strategy),
  seedMgr_(NULL),
  stateCount_(0),
  localStates_(),
  buffer_(),
  bufferCounts_(),
  bufferStrides_(),
  consolidatedBasis_(DynamStatePlainBasis::New(vectorSize))
{}

inline
CpuRank
GlobalStateSharing::localCpu() const {
  return CpuRank(timeComm_->myID());
}

void
GlobalStateSharing::mappingIs(const SliceMapping & mapping) {
  // Partition buffer
  const int cpuCount = mapping.availableCpus().value();
  bufferCounts_.sizeIs(cpuCount);
  bufferStrides_.sizeIs(cpuCount);

  stateCount_ = 0;
  const int stateSize = 2 * vectorSize();
  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    int statesOnCpu = 0;
    for (SliceMapping::SliceIterator it = mapping.hostedSlice(CpuRank(cpu)); it; ++it) {
      const HalfSliceRank rank = *it;
      
      if (rank < mapping.firstActiveSlice()) {
        continue;
      }
      if (rank >= mapping.firstInactiveSlice()) {
        break;
      }

      const HalfSliceCount distance = rank - mapping.firstActiveSlice();
      statesOnCpu += strategy_.stateCount(distance);
    }

    bufferCounts_[cpu] = statesOnCpu * stateSize;
    stateCount_ += statesOnCpu;
  }

  bufferStrides_[0] = 0;
  std::partial_sum(bufferCounts_.array(),
                   bufferCounts_.array() + cpuCount - 1,
                   bufferStrides_.array() + 1);
  
  // Allocate buffer
  const size_t newBufferSize = stateCount_ * stateSize;
  if (buffer_.size() < newBufferSize) {
    buffer_.sizeIs(newBufferSize);
  }

  // Enqueue local seeds to be exchanged
  assert(localStates_.empty());
  
  for (SliceMapping::SliceIterator it = mapping.hostedSlice(localCpu()); it; ++it) {
    const HalfSliceRank rank = *it;

    if (rank < mapping.firstActiveSlice()) {
      continue;
    }
    if (rank >= mapping.firstInactiveSlice()) {
      break;
    }

    const HalfSliceCount distance = rank - mapping.firstActiveSlice();

    for (Strategy::ConstIterator jt = strategy_.sliceSeedTypeBegin(distance);
         jt != strategy_.sliceSeedTypeEnd(distance);
         ++jt) {
      const SeedType type = *jt;
      const HalfSliceRank seedRank = (type == RIGHT_SEED) ? rank : rank.next();
      Seed::Ptr seed = seedMgr_->instance(toString(SeedId(type, seedRank)));
      if (seed) {
        localStates_.push(seed);
      }
    }
  }
}

void
GlobalStateSharing::iterationIs(IterationRank iter) {
  // Fill local buffer
  const int myCpu = localCpu().value();
  double * bufferBegin = buffer_.array() + bufferStrides_[myCpu];
  const int stateSize = 2 * vectorSize();

  assert(localStates_.size() * stateSize == bufferCounts_[myCpu]);

  while (!localStates_.empty()) {
    const DynamState s = localStates_.front()->state();
    bufferStateCopy(s, bufferBegin);
    bufferBegin += stateSize;
    localStates_.pop();
  }

  // Exchange data
  timeComm_->allGatherv(buffer_.array(),
                        bufferCounts_.array(),
                        bufferStrides_.array());

  // Expose data
  consolidatedBasis_ = DynamStateBasisWrapper::New(vectorSize(), stateCount_, buffer_.array());
  setIteration(iter);
}

void
GlobalStateSharing::Strategy::init() {
  if (seedTypes_.find(MAIN_SEED) != seedTypes_.end()) {
    sliceSeedTypes_[0].insert(MAIN_SEED);
  }
  if (seedTypes_.find(LEFT_SEED) != seedTypes_.end()) {
    sliceSeedTypes_[0].insert(LEFT_SEED);
  }
  if (seedTypes_.find(RIGHT_SEED) != seedTypes_.end()) {
    sliceSeedTypes_[1].insert(RIGHT_SEED);
  }
}

} /* end namespace Hts */ } /* end namespace Pita */
