#include "NodalRestrictionMapping.h"
#include "DofSetUtils.h"

#include <Utils.d/dofset.h>

namespace Rom {

NodalRestrictionMapping::InfoType
NodalRestrictionMapping::extractOriginalInfo(const DofSetArray &dsa) {
  return const_cast<DofSetArray &>(dsa).size();
}

void
NodalRestrictionMapping::addSampleNode(int iNode, const DofSetArray &dsa) {
  for (int iDof = 0; iDof < DOF_ID_COUNT; ++iDof) {
    const NodeDof::DofType dofId = DOF_ID(iDof);
    const int originLoc = const_cast<DofSetArray &>(dsa).locate(iNode, dofId);
    if (originLoc >= 0) {
      originIndex_.push_back(originLoc);
    }
  }
}

std::ostream &
operator<<(std::ostream &out, const NodalRestrictionMapping &source) {
  out << "NodalRestrictionMapping: " << source.originInfo() << "->" << source.restrictedInfo() << " :";
  int index = 0;
  typedef std::vector<NodalRestrictionMapping::IndexType>::const_iterator idx_it;
  for (idx_it it = source.originIndex_.begin(); it != source.originIndex_.end(); ++it) {
    out << " " << *it << "->" << index++;
  }

  return out;
}

} /* end namespace Rom */
