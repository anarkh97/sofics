#include "DistrNodeDof6Buffer.h"

#include "RenumberingUtils.h"

#include <algorithm>
#include <iterator>
#include <cstddef>

namespace Rom {

template<int DOFS_PER_NODE>
const double *
DistrNodeDofBuffer<DOFS_PER_NODE>::operator[](int globalNodeIdx) const {
  std::map<int, int>::const_iterator it = localNodeIndices_.find(globalNodeIdx);
  return (it != localNodeIndices_.end()) ? buffer_[it->second] : NULL;
}
  
template<int DOFS_PER_NODE>
void
DistrNodeDofBuffer<DOFS_PER_NODE>::initialize() {
  // Sort and make sure every index is unique (local to global)
  std::sort(globalNodeIndices_.begin(), globalNodeIndices_.end());
  globalNodeIndices_.erase(std::unique(globalNodeIndices_.begin(), globalNodeIndices_.end()), globalNodeIndices_.end());

  // Compute the inverse mapping (global to local) 
  inverse_numbering(globalNodeIndices_.begin(), globalNodeIndices_.end(),
                    std::inserter(localNodeIndices_, localNodeIndices_.end()));
  
  // Resize internal buffer (local indexing)
  buffer_.sizeIs(localNodeCount());
}

template
const double * DistrNodeDofBuffer<6>::operator[](int globalNodeIdx) const;

template
void DistrNodeDofBuffer<6>::initialize();

template
const double * DistrNodeDofBuffer<1>::operator[](int globalNodeIdx) const;

template
void DistrNodeDofBuffer<1>::initialize();

} /* end namespace Rom */
