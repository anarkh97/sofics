#include "VecNodeDof6Map.h"

#include <Utils.d/dofset.h>

#include <deque>

namespace Rom {

VecNodeDof6Map::VecNodeDof6Map(const DofSetArray &dsa) :
  nodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  locationId_(vectorSize()),
  dofLocation_(nodeCount())
{
  for (int iNode = 0, iNodeEnd = nodeCount(); iNode != iNodeEnd; ++iNode) {
    for (int iDof = 0; iDof < DOF_ID_COUNT; ++iDof) {
      const NodeDof::DofType dofId = DOF_ID(iDof);
      const int loc = const_cast<DofSetArray &>(dsa).locate(iNode, dofId);

      dofLocation_[iNode][iDof] = loc;
      if (loc >= 0) {
        locationId_[loc] = NodeDof(iNode, dofId);
      }
    }
  }
}

NodeDof
VecNodeDof6Map::nodeDof(int vecLoc) const {
  assert(vecLoc >= 0 && vecLoc < vectorSize());
  return locationId_[vecLoc]; 
}

} /* end namespace Rom */
