#ifndef ROM_VECNODEDOF6CONVERSION_H
#define ROM_VECNODEDOF6CONVERSION_H

#include "SimpleBuffer.h"

#include <cassert>

class DofSetArray;

namespace Rom {

template<int DOFS_PER_NODE=6>
class VecNodeDofConversion {
public:
  int dofSetNodeCount() const { return dofSetNodeCount_; }
  int vectorSize() const { return vectorSize_; }

  template <typename NodeDofsType, typename VecType>
  const NodeDofsType &paddedNodeDof6(const VecType &origin, NodeDofsType &target) const;

  template <typename NodeDofsType, typename VecType>
  const VecType &vector(const NodeDofsType &origin, VecType &target) const;

  explicit VecNodeDofConversion(const DofSetArray &);
  explicit VecNodeDofConversion(const int numNodes);

private:
  int dofSetNodeCount_;
  int vectorSize_;

  typedef SimpleBuffer<int[DOFS_PER_NODE]> DofLocation;
  DofLocation dofLocation_;

  // Disallow copy and assignment
  VecNodeDofConversion(const VecNodeDofConversion &);
  VecNodeDofConversion &operator=(const VecNodeDofConversion &);
};

typedef VecNodeDofConversion<6> VecNodeDof6Conversion;
typedef VecNodeDofConversion<1> VecNodeDof1Conversion;

template <int DOFS_PER_NODE>
template <typename NodeDofsType, typename VecType>
const NodeDofsType &
VecNodeDofConversion<DOFS_PER_NODE>::paddedNodeDof6(const VecType &origin, NodeDofsType &target) const {
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
    for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      target[iNode][iDof] = (loc >= 0) ? origin[loc] : 0.0;
    }
  }

  return target;
}

template <int DOFS_PER_NODE>
template <typename NodeDofsType, typename VecType>
const VecType &
VecNodeDofConversion<DOFS_PER_NODE>::vector(const NodeDofsType &origin, VecType &target) const {
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
    for (int iDof = 0; iDof < DOFS_PER_NODE; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      if (loc >= 0) {
        target[loc] = origin[iNode][iDof];
      }
    }
  }

  return target;
}

} /* end namespace Rom */

#endif /* ROM_VECNODEDOF6CONVERSION_H */
