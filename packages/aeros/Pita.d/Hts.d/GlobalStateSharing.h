#ifndef PITA_HTS_GLOBALSTATESHARING_H
#define PITA_HTS_GLOBALSTATESHARING_H

#include "Fwk.h"
#include "Types.h"

#include "../NamedTask.h"

#include "SliceMapping.h"

#include "../DynamStateBasis.h"
#include "../DynamState.h"
#include "../Seed.h"

class Communicator;

#include "../SimpleBuffer.h"
#include <map>
#include <queue>
#include <set>

namespace Pita { namespace Hts {

class GlobalStateSharing : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(GlobalStateSharing);

  // Parameters
  CpuRank localCpu() const;
  size_t vectorSize() const { return vectorSize_; }

  // Strategy
  class Strategy {
  public:
    typedef std::set<SeedType>::const_iterator ConstIterator;

    ConstIterator allSeedTypesBegin() const { return seedTypes_.begin(); }
    ConstIterator allSeedTypesEnd()   const { return seedTypes_.end();   }

    ConstIterator sliceSeedTypeBegin(HalfSliceCount distanceFromFirst) const {
      return sliceSeedTypes_[parity(distanceFromFirst)].begin();
    }
    ConstIterator sliceSeedTypeEnd(HalfSliceCount distanceFromFirst) const {
      return sliceSeedTypes_[parity(distanceFromFirst)].end();
    }

    int stateTypeCount() const { return seedTypes_.size(); }
    int stateCount(HalfSliceCount distanceFromFirst) const { return sliceSeedTypes_[parity(distanceFromFirst)].size(); }

    template<typename InputIterator>
    Strategy(InputIterator begin, InputIterator end) :
      seedTypes_(begin, end)
    {
      init();
    }

  private:
    std::set<SeedType> seedTypes_;
    std::set<SeedType> sliceSeedTypes_[2];

    static int parity(HalfSliceCount distance) { return distance.value() % 2; }
    void init();
  };

  // Input
  const Seed::Manager * seedMgr() const { return seedMgr_.ptr(); }
  void seedMgrIs(Seed::Manager * sg) { seedMgr_ = sg; }

  // Output
  const DynamStateBasis * consolidatedBasis() const { return consolidatedBasis_.ptr(); }
  
  // Setup exchange parameters
  virtual void mappingIs(const SliceMapping & m);
  
  // Execution
  virtual void iterationIs(IterationRank iter); // overriden

  GlobalStateSharing(Communicator * timeComm, size_t vectorSize, Strategy strategy);

private:
  Communicator * timeComm_;
  size_t vectorSize_;
  Strategy strategy_;

  Seed::Manager::Ptr seedMgr_;

  std::queue<Seed::PtrConst> localStates_;
  int stateCount_;
  SimpleBuffer<double> buffer_;
  SimpleBuffer<int> bufferCounts_;
  SimpleBuffer<int> bufferStrides_; 

  DynamStateBasis::Ptr consolidatedBasis_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_GLOBALSTATESHARING_H */
