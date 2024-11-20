#ifndef ROM_BASISFILEITERATOR_H
#define ROM_BASISFILEITERATOR_H

#include "BasisBinaryFile.h"
#include "VecNodeDof6Conversion.h"

#include "NodeDof6Buffer.h"
#include "MappedNodeDof6Buffer.h"

#include <string>
#include <utility>
#include <cstddef>

namespace Rom {

// Input Stream

template<int DOFS_PER_NODE>
class BasisInputStream {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
  BasisBinaryInputFile& file() { return file_; }

  // Check the readiness of the stream
  int currentVectorRank() const { return file_.currentStateIndex(); }
  operator const void*() const;

  BasisInputStream(const std::string &fileName, const VecNodeDofConversion<DOFS_PER_NODE> &converter);

  template <typename VectorBufferType, int dofs_per_node>
  friend BasisInputStream<dofs_per_node> &operator>>(BasisInputStream<dofs_per_node> &, VectorBufferType &);
  
  // Convenience overload (to bind with non-const rvalues)
  template<int dofs_per_node>
  friend BasisInputStream<dofs_per_node> &operator>>(BasisInputStream<dofs_per_node> &, double *);
  
  template <typename VectorBufferType, int dofs_per_node>
  friend BasisInputStream<dofs_per_node> &operator>>(BasisInputStream<dofs_per_node> &, std::pair<double, VectorBufferType> &);

private:
  template <typename VectorBufferType> void performInput(VectorBufferType &);
  template <typename VectorBufferType> void performUncheckedInput(VectorBufferType &);
 
  void checkInput();

  template <typename VectorBufferType>
  void performInput(std::pair<double, VectorBufferType> &);

  BasisBinaryInputFile file_;
  bool isValid_;
  const VecNodeDofConversion<DOFS_PER_NODE> &converter_;
  MappedNodeDofBuffer<DOFS_PER_NODE> buffer_;
};

template <int DOFS_PER_NODE>
inline
BasisInputStream<DOFS_PER_NODE>::operator const void*() const {
  return isValid_ ? this : NULL;
}

template <int DOFS_PER_NODE>
inline
void
BasisInputStream<DOFS_PER_NODE>::checkInput() {
  isValid_ = file_.validCurrentState();
}

template <int DOFS_PER_NODE>
template <typename VectorBufferType>
inline
void
BasisInputStream<DOFS_PER_NODE>::performUncheckedInput(VectorBufferType &target) {
  file_.currentStateBuffer(buffer_.underlyingBuffer()); // Fill buffer in one pass, ignoring the file node mapping
  converter_.vector(buffer_, target); // Fill vector according to the file node mapping 
  file_.currentStateIndexInc();
}

template <int DOFS_PER_NODE>
template <typename VectorBufferType>
inline
void
BasisInputStream<DOFS_PER_NODE>::performInput(VectorBufferType &target) {
  checkInput();
  if (isValid_) {
    performUncheckedInput(target);
  }
}
 
template <int DOFS_PER_NODE>
template <typename VectorBufferType>
inline
void
BasisInputStream<DOFS_PER_NODE>::performInput(std::pair<double, VectorBufferType> &target) {
  checkInput();
  if (isValid_) {
    target.first = file_.currentStateHeaderValue();
    performUncheckedInput(target.second);
  }
}

// Input operations

template <typename VectorBufferType, int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
operator>>(BasisInputStream<DOFS_PER_NODE> &in, VectorBufferType &target) {
  in.performInput(target);
  return in;
}

template <int DOFS_PER_NODE>
inline
BasisInputStream<DOFS_PER_NODE> &
operator>>(BasisInputStream<DOFS_PER_NODE> &in, double *target) {
  in.performInput(target);
  return in;
}

template <typename VectorBufferType, int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
operator>>(BasisInputStream<DOFS_PER_NODE> &in, std::pair<double, VectorBufferType> &target) {
  in.performInput(target);
  return in;
}

template <typename FwdIter, int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readVectors(BasisInputStream<DOFS_PER_NODE> &in, FwdIter first, FwdIter last) {
  FwdIter it = first;
  while (it != last && in >> *it++) {
    // Nothing to do
  }
  return in;
}

// Output Stream

template<int DOFS_PER_NODE>
class BasisOutputStream {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
  
  BasisOutputStream(const std::string &fileName, const VecNodeDofConversion<DOFS_PER_NODE> &converter, bool restart);

  template <typename VectorBufferType, int dofs_per_node>
  friend BasisOutputStream<dofs_per_node> &operator<<(BasisOutputStream<dofs_per_node> &, const VectorBufferType &source);
  
  template <typename VectorBufferType, int dofs_per_node>
  friend BasisOutputStream<dofs_per_node> &operator<<(BasisOutputStream<dofs_per_node> &, const std::pair<double, VectorBufferType> &source);
  
private:
  template <typename VectorBufferType>
  const NodeDofBuffer<DOFS_PER_NODE> &convert(const VectorBufferType &);

  BasisBinaryOutputFile file_;
  const VecNodeDofConversion<DOFS_PER_NODE> &converter_;

  NodeDofBuffer<DOFS_PER_NODE> buffer_;
};

// Output operations

template <int DOFS_PER_NODE>
template <typename VectorBufferType>
inline
const NodeDofBuffer<DOFS_PER_NODE> &
BasisOutputStream<DOFS_PER_NODE>::convert(const VectorBufferType &source) {
  return converter_.paddedNodeDof6(source, buffer_);
}

template <typename VectorBufferType, int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
operator<<(BasisOutputStream<DOFS_PER_NODE> &out, const std::pair<double, VectorBufferType> &source) {
  out.file_.stateAdd(out.convert(source.second), source.first);
  return out;
}

template <typename VectorBufferType, int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
operator<<(BasisOutputStream<DOFS_PER_NODE> &out, const VectorBufferType &source) {
  out.file_.stateAdd(out.convert(source));
  return out;
}

template <typename InputIter, int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeVectors(BasisOutputStream<DOFS_PER_NODE> &out, InputIter first, InputIter last) {
  InputIter it = first;
  while (it != last) {
    out << *it++;
  }
  return out;
}

} /* end namespace Rom */

#endif /* ROM_BASISFILEITERATOR_H */
