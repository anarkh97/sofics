#ifndef ROM_DISTRVECNODEDOF6CONVERSION_H
#define ROM_DISTRVECNODEDOF6CONVERSION_H

#include "RestrictedVecNodeDof6Conversion.h"
#include "VecNodeDof6Conversion.h"
#include "DistrDomainUtils.h"
#include "BlockCyclicMap.h"

#include <Driver.d/SubDomain.h> 

#include <vector>
#include <cassert>

namespace Rom {

template<int DOFS_PER_NODE=6>
class DistrVecNodeDofConversion {
public:
  int subDomainCount() const { return subDomains_.size(); }

  template <typename NodeDofType, typename VecType>
  const NodeDofType &paddedNodeDof6(const VecType &origin, NodeDofType &target) const;

  template <typename NodeDofType, typename VecType>
  const NodeDofType &unpaddedNodeDof6(const VecType &origin, NodeDofType &target) const;

  template <typename NodeDofType, typename VecType>
  const VecType &vector(const NodeDofType &origin, VecType &target) const;
  
  template <typename NodeDofType, typename VecType>
  const VecType &paddedMasterVector(const NodeDofType &origin, VecType &target) const;

  template <typename NodeDofType, typename VecType>
  const VecType &unpaddedMasterVector(const NodeDofType &origin, VecType &target) const;

  template <typename SubDomPtrFwdIt>
  DistrVecNodeDofConversion(SubDomPtrFwdIt first, SubDomPtrFwdIt last, bool = true);
  template <typename SubDomPtrFwdIt>
  DistrVecNodeDofConversion(SubDomPtrFwdIt first, SubDomPtrFwdIt last, int globalLen, int blockSize,
                            Communicator *com, int numLocalSub);

  ~DistrVecNodeDofConversion();

private:
  typedef std::vector<const SubDomain *> SubDomContainer;
  typedef std::vector<const VecNodeDofConversion<DOFS_PER_NODE> *> ConversionContainer;
  typedef std::vector<const RestrictedVecNodeDofConversion<DOFS_PER_NODE> *> RestrictedConversionContainer;
  SubDomContainer subDomains_;
  BlockCyclicMap bcMap_;
  Communicator *com_;
  ConversionContainer subConversions_;
  RestrictedConversionContainer subRestrictedConversions_;

  // Disallow copy and assignment
  DistrVecNodeDofConversion(const DistrVecNodeDofConversion &);
  DistrVecNodeDofConversion &operator=(const DistrVecNodeDofConversion &);
};

typedef DistrVecNodeDofConversion<6> DistrVecNodeDof6Conversion;
typedef DistrVecNodeDofConversion<1> DistrVecNodeDof1Conversion;

template <int DOFS_PER_NODE>
template <typename SubDomPtrFwdIt>
DistrVecNodeDofConversion<DOFS_PER_NODE>::DistrVecNodeDofConversion(SubDomPtrFwdIt first, SubDomPtrFwdIt last, bool flag) :
  subDomains_(first, last), com_(NULL)
{
  for (SubDomPtrFwdIt it = first; it != last; ++it) {
    SubDomain * const s = *it;
    if(flag) {
      subConversions_.push_back(new VecNodeDofConversion<DOFS_PER_NODE>(*s->getCDSA()));

      std::vector<bool> masterFlags;
      master_node_flags(*s, std::back_inserter(masterFlags));
      subRestrictedConversions_.push_back(new RestrictedVecNodeDofConversion<DOFS_PER_NODE>(*s->getCDSA(),
                                                                                            masterFlags.begin(),
                                                                                            masterFlags.end()));
    }
    else {
      subConversions_.push_back(new VecNodeDofConversion<DOFS_PER_NODE>(s->getCDSA()->numNodes()));
    }
  }
}

template <int DOFS_PER_NODE>
template <typename SubDomPtrFwdIt>
DistrVecNodeDofConversion<DOFS_PER_NODE>::DistrVecNodeDofConversion(SubDomPtrFwdIt first, SubDomPtrFwdIt last, int globalLen,
                                                                    int blockSize, Communicator *com, int numLocalSub) :
  subDomains_(first, last), bcMap_(globalLen, blockSize, com->numCPUs(), numLocalSub), com_(com)
{ 
  for (SubDomPtrFwdIt it = first; it != last; ++it) {
    SubDomain * const s = *it;
    int localLen = bcMap_.subLen(com_->myID(), s->localSubNum());
    subConversions_.push_back(new VecNodeDofConversion<DOFS_PER_NODE>(localLen));
  }
}

template <typename NodeDofType>
class SubNodeDofAdapter {
public:
  typedef double Scalar; 

  explicit SubNodeDofAdapter(const NodeDofType &nodeDof, const SubDomain &subDom, const BlockCyclicMap &bcMap, Communicator *com) :
    nodeDof_(nodeDof), subDomain_(subDom), bcMap_(bcMap), com_(com)
  {}

  const Scalar *operator[](int locIdx) const {
    SubDomain &s = const_cast<SubDomain &>(subDomain_);
    const int globIdx = (com_) ? bcMap_.localToGlobal(com_->myID(), s.localSubNum(), locIdx) : s.localToGlobal(locIdx);
    const Scalar *result = nodeDof_[globIdx];
    assert(result);
    return result;
  }

  Scalar *operator[](int locIdx) {
    const SubNodeDofAdapter &self = *this;
    return const_cast<Scalar *>(self[locIdx]);
  }

private:
  const NodeDofType &nodeDof_;
  const SubDomain &subDomain_;
  const BlockCyclicMap &bcMap_;
  Communicator *com_;
};

template <int DOFS_PER_NODE>
template <typename NodeDofType, typename VecType>
const NodeDofType &
DistrVecNodeDofConversion<DOFS_PER_NODE>::paddedNodeDof6(const VecType &origin, NodeDofType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    SubNodeDofAdapter<NodeDofType> subNodeDof(target, *subDomains_[iSub], bcMap_, com_);
    const GenStackVector<double> subVector(const_cast<VecType &>(origin).subData(iSub),
                                           const_cast<VecType &>(origin).subLen(iSub));
    if(subRestrictedConversions_.empty()) subConversions_[iSub]->paddedNodeDof6(subVector, subNodeDof);
    else subRestrictedConversions_[iSub]->paddedNodeDof6(subVector, subNodeDof);
  }

  return target;
}

template <int DOFS_PER_NODE>
template <typename NodeDofType, typename VecType>
const NodeDofType &
DistrVecNodeDofConversion<DOFS_PER_NODE>::unpaddedNodeDof6(const VecType &origin, NodeDofType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    SubNodeDofAdapter<NodeDofType> subNodeDof(target, *subDomains_[iSub], bcMap_, com_);
    const GenStackVector<double> subVector(const_cast<VecType &>(origin).subData(iSub),
                                           const_cast<VecType &>(origin).subLen(iSub));
    subRestrictedConversions_[iSub]->unpaddedNodeDof6(subVector, subNodeDof);
  }

  return target;
}

template <int DOFS_PER_NODE>
template <typename NodeDofType, typename VecType>
const VecType &
DistrVecNodeDofConversion<DOFS_PER_NODE>::vector(const NodeDofType &origin, VecType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    const SubNodeDofAdapter<NodeDofType> subNodeDof(origin, *subDomains_[iSub], bcMap_, com_);
    GenStackVector<double> subVector(target.subData(iSub), target.subLen(iSub));
    subConversions_[iSub]->vector(subNodeDof, subVector);
  }

  return target;
};

template <int DOFS_PER_NODE>
template <typename NodeDofType, typename VecType>
const VecType &
DistrVecNodeDofConversion<DOFS_PER_NODE>::paddedMasterVector(const NodeDofType &origin, VecType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    const SubNodeDofAdapter<NodeDofType> subNodeDof(origin, *subDomains_[iSub], bcMap_, com_);
    GenStackVector<double> subVector(target.subData(iSub), target.subLen(iSub));
    subRestrictedConversions_[iSub]->paddedVector(subNodeDof, subVector);
  }

  return target;
};

template <int DOFS_PER_NODE>
template <typename NodeDofType, typename VecType>
const VecType &
DistrVecNodeDofConversion<DOFS_PER_NODE>::unpaddedMasterVector(const NodeDofType &origin, VecType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    const SubNodeDofAdapter<NodeDofType> subNodeDof(origin, *subDomains_[iSub], bcMap_, com_);
    GenStackVector<double> subVector(target.subData(iSub), target.subLen(iSub));
    subRestrictedConversions_[iSub]->unpaddedVector(subNodeDof, subVector);
  }

  return target;
};

} // end namespace Rom

#endif /* ROM_DISTRVECNODEDOF6CONVERSION_H */
