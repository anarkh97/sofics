#ifndef ROM_DISTRNODEDOF6BUFFER_H
#define ROM_DISTRNODEDOF6BUFFER_H

#include "NodeDof6Buffer.h"

#include <vector>
#include <map>

namespace Rom {

template<int DOFS_PER_NODE>
class DistrNodeDofBuffer {
public:
  int localNodeCount() const { return globalNodeIndices_.size(); }
  int dofsPerNode() const { return DOFS_PER_NODE; }

  typedef std::vector<int>::const_iterator NodeItConst;
  NodeItConst globalNodeIndexBegin() const { return globalNodeIndices_.begin(); }
  NodeItConst globalNodeIndexEnd()   const { return globalNodeIndices_.end(); }

  const double *operator[](int globalNodeIdx) const;
  double *operator[](int globalNodeIdx) {
    const DistrNodeDofBuffer &self = *this;
    return const_cast<double *>(self[globalNodeIdx]);
  }

  void zero() { buffer_.zero(); }

  template <typename IdxInIt>
  DistrNodeDofBuffer(IdxInIt first, IdxInIt last, int offset = 0);

private:
  void initialize();

  std::map<int, int> localNodeIndices_;
  std::vector<int> globalNodeIndices_;
  NodeDofBuffer<DOFS_PER_NODE> buffer_;
};

typedef DistrNodeDofBuffer<6> DistrNodeDof6Buffer;
typedef DistrNodeDofBuffer<1> DistrNodeDof1Buffer;

template <int DOFS_PER_NODE>
template <typename IdxInIt>
DistrNodeDofBuffer<DOFS_PER_NODE>::DistrNodeDofBuffer(IdxInIt first, IdxInIt last, int offset) :
  localNodeIndices_(),
  globalNodeIndices_(first, last),
  buffer_()
{
  if(offset != 0) {
    for(std::vector<int>::iterator it = globalNodeIndices_.begin(); it != globalNodeIndices_.end(); ++it) (*it) -= offset;
  }
  initialize();
}

} /* end namespace Rom */

#endif /* ROM_DISTRNODEDOF6BUFFER_H */
