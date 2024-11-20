#include "MappedNodeDof6Buffer.h"

#include "RenumberingUtils.h"

#include <iterator>
#include <cstddef>

namespace Rom {

template<int DOFS_PER_NODE>
const double *
MappedNodeDofBuffer<DOFS_PER_NODE>::operator[](int iNode) const {
  std::map<int, int>::const_iterator it = underlyingNodeIndices_.find(iNode);
  return (it != underlyingNodeIndices_.end()) ? buffer_[it->second] : NULL;
}
 
template<int DOFS_PER_NODE>
void
MappedNodeDofBuffer<DOFS_PER_NODE>::initialize() {
  // Compute the inverse mapping (usual to underlying) 
  inverse_numbering(nodeIndices_.begin(), nodeIndices_.end(),
                    std::inserter(underlyingNodeIndices_, underlyingNodeIndices_.end()));
  
  // Resize internal buffer (underlying indexing)
  buffer_.sizeIs(nodeIndices_.size());
}

template
const double * MappedNodeDofBuffer<6>::operator[](int iNode) const;

template
void MappedNodeDofBuffer<6>::initialize();

template
const double * MappedNodeDofBuffer<1>::operator[](int iNode) const;

template
void MappedNodeDofBuffer<1>::initialize();

} /* end namespace Rom */
