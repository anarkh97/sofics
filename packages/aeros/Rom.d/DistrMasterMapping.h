#ifndef ROM_DISTRMASTERMAPPING_H
#define ROM_DISTRMASTERMAPPING_H

#include "MasterMapping.h"
#include "DistrDomainUtils.h"
#include "BlockCyclicMap.h"

#include <Driver.d/SubDomain.h>

#include <vector>

namespace Rom {

class DistrMasterMapping {
public:
  int subDomainCount() const { return subMapping_.size(); }
  int localNodeCount() const { return localNodes_.size(); }
  int masterNodeCount() const { return masterNodes_.size(); }

  typedef std::vector<MasterMapping>::const_iterator SubMasterMappingIt;
  SubMasterMappingIt begin() const { return subMapping_.begin(); }
  SubMasterMappingIt end()   const { return subMapping_.end();   }

  typedef std::vector<int>::const_iterator NodeIndexIt;
  NodeIndexIt localNodeBegin() const { return localNodes_.begin(); }
  NodeIndexIt localNodeEnd()   const { return localNodes_.end();   }
  NodeIndexIt masterNodeBegin() const { return masterNodes_.begin(); }
  NodeIndexIt masterNodeEnd()   const { return masterNodes_.end();   }

  DistrMasterMapping() {}
  template <typename SubDomFwdIt>
  DistrMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast, bool internalNodes=false);

protected:
  std::vector<MasterMapping> subMapping_;
  std::vector<int> localNodes_;
  std::vector<int> masterNodes_;
};

template <typename SubDomFwdIt>
DistrMasterMapping::DistrMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast, bool internalNodes) {
  for (SubDomFwdIt subDomIt = subDomFirst; subDomIt != subDomLast; ++subDomIt) {
    SubDomain &sd = *subDomIt;
    auto &globalNodes = sd.getGlNodes();

    localNodes_.insert(localNodes_.end(), globalNodes.begin(), globalNodes.end());

    std::vector<bool> masterNodes;
    if(internalNodes) internal_node_flags(sd, std::back_inserter(masterNodes));
    else master_node_flags(sd, std::back_inserter(masterNodes));
    subMapping_.push_back(MasterMapping(globalNodes.begin(), globalNodes.end(), masterNodes.begin()));

    std::vector<bool>::const_iterator flagIt = masterNodes.begin();
    for (auto nd : globalNodes) {
      if (*flagIt++) {
        masterNodes_.push_back(nd);
      }
    }
  }
}

class DistrMpcMasterMapping : public DistrMasterMapping {
public:
  template <typename SubDomFwdIt>
  DistrMpcMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast);
};

template <typename SubDomFwdIt>
DistrMpcMasterMapping::DistrMpcMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast) {
  for (SubDomFwdIt subDomIt = subDomFirst; subDomIt != subDomLast; ++subDomIt) {
    SubDomain &sd = *subDomIt;
    const int * const globalBegin = sd.getGlMPCs();
    const int * const globalEnd = globalBegin + sd.getNumMpc();

    localNodes_.insert(localNodes_.end(), globalBegin, globalEnd);

    std::vector<bool> masterNodes(sd.getMpcMaster(), sd.getMpcMaster()+sd.getNumMpc());
    subMapping_.push_back(MasterMapping(globalBegin, globalEnd, masterNodes.begin()));

    std::vector<bool>::const_iterator flagIt = masterNodes.begin();
    for (const int *globalIt = globalBegin; globalIt != globalEnd; ++globalIt) {
      if (*flagIt++) {
        masterNodes_.push_back(*globalIt);
      }
    }
  }
}

class DistrTrivialMasterMapping : public DistrMasterMapping {
public:
  template <typename SubDomFwdIt>
  DistrTrivialMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast, int globalLen,
                            int blockSize, Communicator *com, int numLocalSub);
private:
  BlockCyclicMap bcMap_;
};

template <typename SubDomFwdIt>
DistrTrivialMasterMapping::DistrTrivialMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast,
                                                     int globalLen, int blockSize, Communicator *com,
                                                     int numLocalSub)
  : bcMap_(globalLen, blockSize, com->numCPUs(), numLocalSub)
{
  for (SubDomFwdIt subDomIt = subDomFirst; subDomIt != subDomLast; ++subDomIt) {
    SubDomain &sd = *subDomIt;
    const int locLen = bcMap_.subLen(com->myID(), sd.localSubNum());
    for(int i=0; i<locLen; ++i) {
      int glIndx = bcMap_.localToGlobal(com->myID(), sd.localSubNum(), i);
      localNodes_.push_back(glIndx); 
      masterNodes_.push_back(glIndx);
    }
  }
}

} // end namespace Rom


#endif /* ROM_DISTRMASTERMAPPING_H */
