#ifndef ROM_RESTRICTEDVECNODEDOF6CONVERSION_H
#define ROM_RESTRICTEDVECNODEDOF6CONVERSION_H

#include "DofSetUtils.h"
#include "SimpleBuffer.h"

#include <Utils.d/dofset.h>

#include <vector>
#include <algorithm>

#include <cassert>

namespace Rom {

template<int DOFS_PER_NODE=6>
class RestrictedVecNodeDofConversion {
public:
  int dofSetNodeCount() const { return dofSetNodeCount_; }
  int nodeCount() const { return nodeCount_; }
  int vectorSize() const { return vectorSize_; }

  template <typename NodeDofsType, typename VecType>
  const NodeDofsType &paddedNodeDof6(const VecType &origin, NodeDofsType &target) const;

  template <typename NodeDofsType, typename VecType>
  const NodeDofsType &unpaddedNodeDof6(const VecType &origin, NodeDofsType &target) const;

  template <typename NodeDofsType, typename VecType>
  const VecType &paddedVector(const NodeDofsType &origin, VecType &target) const;

  template <typename NodeDofsType, typename VecType>
  const VecType &unpaddedVector(const NodeDofsType &origin, VecType &target) const;

  template <typename BoolFwdIt>
  RestrictedVecNodeDofConversion(const DofSetArray &dsa, BoolFwdIt nodeMaskBegin,
                                                         BoolFwdIt nodeMaskEnd);
private:
  int dofSetNodeCount_;
  int vectorSize_;

  std::vector<NodeDof> locationId_;
  
  typedef SimpleBuffer<int[DOFS_PER_NODE]> DofLocation;
  DofLocation dofLocation_;

  std::vector<bool> nodeMask_;
  int nodeCount_;

  void initialize(const DofSetArray &dsa);

  // Disallow copy and assignment
  RestrictedVecNodeDofConversion(const RestrictedVecNodeDofConversion &);
  RestrictedVecNodeDofConversion &operator=(const RestrictedVecNodeDofConversion &);
};

typedef RestrictedVecNodeDofConversion<6> RestrictedVecNodeDof6Conversion;
typedef RestrictedVecNodeDofConversion<1> RestrictedVecNodeDof1Conversion;

template <int DOFS_PER_NODE>
template <typename NodeDofsType, typename VecType>
const NodeDofsType &
RestrictedVecNodeDofConversion<DOFS_PER_NODE>::paddedNodeDof6(const VecType &origin, NodeDofsType &target) const {
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
    if (nodeMask_[iNode]) {
      for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        target[iNode][iDof] = (loc >= 0) ? origin[loc] : 0.0;
      }
    }
  }

  return target;
}

template <int DOFS_PER_NODE>
template <typename NodeDofsType, typename VecType>
const NodeDofsType &
RestrictedVecNodeDofConversion<DOFS_PER_NODE>::unpaddedNodeDof6(const VecType &origin, NodeDofsType &target) const {
  int pos = 0;
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
    if (nodeMask_[iNode]) {
      for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        target[iNode][iDof] = (loc >= 0) ? origin[pos++] : 0.0;
      }
    }
  }

  return target;
}

template <int DOFS_PER_NODE>
template <typename NodeDofsType, typename VecType>
const VecType &
RestrictedVecNodeDofConversion<DOFS_PER_NODE>::paddedVector(const NodeDofsType &origin, VecType &target) const {
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
    if (nodeMask_[iNode]) {
      for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        if (loc >= 0) {
          target[loc] = origin[iNode][iDof];
        }
      }
    } else {
      for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        if (loc >= 0) {
          target[loc] = 0.0;
        }
      }
    }
  }
  return target;
}

template <int DOFS_PER_NODE>
template <typename NodeDofsType, typename VecType>
const VecType &
RestrictedVecNodeDofConversion<DOFS_PER_NODE>::unpaddedVector(const NodeDofsType &origin, VecType &target) const {
  int pos = 0;
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
    if (nodeMask_[iNode]) {
      for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        if (loc >= 0) {
          target[pos++] = origin[iNode][iDof];
        }
      }
    }
  }
  return target;
}

template <int DOFS_PER_NODE>
template <typename BoolFwdIt>
RestrictedVecNodeDofConversion<DOFS_PER_NODE>::RestrictedVecNodeDofConversion(const DofSetArray &dsa,
                                                                              BoolFwdIt nodeMaskBegin,
                                                                              BoolFwdIt nodeMaskEnd) :
  dofSetNodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  locationId_(vectorSize()),
  dofLocation_(dofSetNodeCount()),
  nodeMask_(nodeMaskBegin, nodeMaskEnd),
  nodeCount_(std::count(nodeMask_.begin(), nodeMask_.end(), true))
{
  initialize(dsa);
}

} /* end namespace Rom */

#endif /* ROM_RESTRICTEDVECNODEDOF6CONVERSION_H */
