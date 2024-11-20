#ifndef ROM_MAPPEDNODEDOF6BUFFER_H
#define ROM_MAPPEDNODEDOF6BUFFER_H

#include "NodeDof6Buffer.h"

#include <vector>
#include <map>

namespace Rom {

template<int DOFS_PER_NODE>
class MappedNodeDofBuffer {
public:
  int nodeCount() const { return nodeIndices_.size(); }

  typedef std::vector<int>::const_iterator NodeItConst;
  NodeItConst nodeIndexBegin() const { return nodeIndices_.begin(); }
  NodeItConst nodeIndexEnd()   const { return nodeIndices_.end(); }

  const double *operator[](int iNode) const;
  double *operator[](int iNode) {
    const MappedNodeDofBuffer &self = *this;
    return const_cast<double *>(self[iNode]);
  }

  double *array() { return buffer_.array(); }
  const double *array() const { return buffer_.array(); }

  const NodeDofBuffer<DOFS_PER_NODE> &underlyingBuffer() const { return buffer_; }
  NodeDofBuffer<DOFS_PER_NODE> &underlyingBuffer() { return buffer_; }

  // Range [first, last) should not have duplicated elements
  template <typename IdxInIt>
  MappedNodeDofBuffer(IdxInIt first, IdxInIt last);

private:
  void initialize();

  std::map<int, int> underlyingNodeIndices_;
  std::vector<int> nodeIndices_;
  NodeDofBuffer<DOFS_PER_NODE> buffer_;
};

typedef MappedNodeDofBuffer<6> MappedNodeDof6Buffer;

template <int DOFS_PER_NODE>
template <typename IdxInIt>
MappedNodeDofBuffer<DOFS_PER_NODE>::MappedNodeDofBuffer(IdxInIt first, IdxInIt last) :
  underlyingNodeIndices_(),
  nodeIndices_(first, last),
  buffer_()
{
  initialize();
}

} /* end namespace Rom */

#endif /* ROM_MAPPEDNODEDOF6BUFFER_H */
