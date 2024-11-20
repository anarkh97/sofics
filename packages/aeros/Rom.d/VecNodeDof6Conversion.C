#include "VecNodeDof6Conversion.h"
#include "DofSetUtils.h"

#include <Utils.d/dofset.h>

namespace Rom {

template<int DOFS_PER_NODE>
VecNodeDofConversion<DOFS_PER_NODE>::VecNodeDofConversion(const DofSetArray &dsa) :
  dofSetNodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  dofLocation_(dofSetNodeCount())
{
  for (int iNode = 0, iNodeEnd = dofSetNodeCount(); iNode < iNodeEnd; ++iNode) {
    for (int iDof = 0; iDof != DOFS_PER_NODE; ++iDof) {
      const int loc = const_cast<DofSetArray &>(dsa).locate(iNode, DOF_ID(iDof));
      dofLocation_[iNode][iDof] = loc;
    }
  }
}

template<int DOFS_PER_NODE>
VecNodeDofConversion<DOFS_PER_NODE>::VecNodeDofConversion(const int numNodes) :
  dofSetNodeCount_(numNodes),
  vectorSize_(numNodes*DOFS_PER_NODE),
  dofLocation_(dofSetNodeCount())
{
  for (int iNode = 0, loc = 0, iNodeEnd = dofSetNodeCount(); iNode < iNodeEnd; ++iNode) {
    for (int iDof = 0; iDof != DOFS_PER_NODE; ++iDof, ++loc) {
      dofLocation_[iNode][iDof] = loc;
    }
  }
}

template
VecNodeDofConversion<6>::VecNodeDofConversion(const DofSetArray &dsa);

template
VecNodeDofConversion<6>::VecNodeDofConversion(const int numNodes);

template
VecNodeDofConversion<1>::VecNodeDofConversion(const DofSetArray &dsa);

template
VecNodeDofConversion<1>::VecNodeDofConversion(const int numNodes);

} /* end namespace Rom */
