#ifndef ROM_NODEDOF6BUFFER_H
#define ROM_NODEDOF6BUFFER_H

#include "SimpleBuffer.h"

#include <cstddef>

namespace Rom {

template<int DOFS_PER_NODE>
class NodeDofBuffer {
public:
  explicit NodeDofBuffer(size_t nodeCount = 0) :
    nodeCount_(nodeCount),
    buffer_(DOFS_PER_NODE * nodeCount)
  {}

  size_t size() const { return nodeCount_; }
  void sizeIs(size_t nodeCount) { buffer_.sizeIs(DOFS_PER_NODE * nodeCount); nodeCount_ = nodeCount; }
  void zero() { for(int i = 0; i < buffer_.size(); ++i) buffer_[i] = 0.0; }; 

  int dofsPerNode() const { return DOFS_PER_NODE; }

  const double *operator[](size_t iNode) const { return buffer_.array() + (DOFS_PER_NODE * iNode); }
  double *operator[](size_t iNode) { return buffer_.array() + (DOFS_PER_NODE * iNode); }

  double *array() { return buffer_.array(); }
  const double *array() const { return buffer_.array(); }

private:
  size_t nodeCount_;
  SimpleBuffer<double> buffer_;

  // Disallow copy and assignment
  NodeDofBuffer(const NodeDofBuffer &);
  NodeDofBuffer &operator=(const NodeDofBuffer &);
};

typedef NodeDofBuffer<1> NodeDof1Buffer;
typedef NodeDofBuffer<6> NodeDof6Buffer;

} /* end namespace Rom */

#endif /* ROM_NODEDOF6BUFFER_H */
