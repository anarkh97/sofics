#ifndef ROM_VECNODEDOF6MAP_H
#define ROM_VECNODEDOF6MAP_H

#include "DofSetUtils.h"
#include "SimpleBuffer.h"

#include <vector>
#include <cassert>

class DofSetArray;

namespace Rom {

class VecNodeDof6Map {
public:
  int nodeCount() const { return nodeCount_; }
  int vectorSize() const { return vectorSize_; }

  NodeDof nodeDof(int vecLoc) const;

  template <typename IndexOut>
  IndexOut locations(int nodeRank, IndexOut result) const;

  explicit VecNodeDof6Map(const DofSetArray &);

private:
  int nodeCount_;
  int vectorSize_;

  std::vector<NodeDof> locationId_;
  
  typedef SimpleBuffer<int[6]> DofLocation;
  DofLocation dofLocation_;

  // Disallow copy and assignment
  VecNodeDof6Map(const VecNodeDof6Map &);
  VecNodeDof6Map &operator=(const VecNodeDof6Map &);
};

template <typename IndexOut>
IndexOut
VecNodeDof6Map::locations(int nodeRank, IndexOut result) const {
  assert(nodeRank >= 0 && nodeRank < nodeCount());
  const int *nodeLoc = dofLocation_[nodeRank];

  for (int iDof = 0; iDof < 6; ++iDof) {
    const int loc = nodeLoc[iDof];
    if (loc >= 0) {
      *result++ = loc;
    }
  }

  return result;
}

} /* end namespace Rom */

#endif /* ROM_VECNODEDOF6MAP_H */
