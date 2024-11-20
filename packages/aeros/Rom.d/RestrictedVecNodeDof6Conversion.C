#include "RestrictedVecNodeDof6Conversion.h"

namespace Rom {

template<int DOFS_PER_NODE>
void
RestrictedVecNodeDofConversion<DOFS_PER_NODE>::initialize(const DofSetArray &dsa) {
  for (int iNode = 0; iNode < dofSetNodeCount(); ++iNode) {
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

template
void RestrictedVecNodeDofConversion<6>::initialize(const DofSetArray &);

template
void RestrictedVecNodeDofConversion<1>::initialize(const DofSetArray &);

} /* end namespace Rom */
