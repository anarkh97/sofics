#include "VecBasisFile.h"

#include "BasisFileStream.h"
#include "NodalRestrictionMapping.h"

#include <Math.d/Vector.h>

#include <cassert>
#include <algorithm>
#include <new>

namespace Rom {

// Output
template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
operator<<(BasisOutputStream<DOFS_PER_NODE> &out, const VecBasis &source) {
  assert(out.vectorSize() == source.vectorSize());
  return writeVectors(out, source.begin(), source.end());
}

template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeVectors(BasisOutputStream<DOFS_PER_NODE> &out, const VecBasis &source, int countMax) {
  assert(out.vectorSize() == source.vectorSize());
  VecBasis::const_iterator last = std::min(source.end(), source.begin() + countMax);
  return writeVectors(out, source.begin(), last);
}

template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeExtendedVectors(BasisOutputStream<DOFS_PER_NODE> &out, const VecBasis &source, const NodalRestrictionMapping &mapping) {
  assert(out.vectorSize() == mapping.originInfo());
  assert(source.vectorSize() == mapping.restrictedInfo());
  
  Vector buffer(mapping.originInfo());
  for (VecBasis::const_iterator vecIt = source.begin(); vecIt != source.end(); ++vecIt) {
    mapping.extension(*vecIt, buffer);
    out << buffer;
  }
  return out;
}

template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeRestrictedVectors(BasisOutputStream<DOFS_PER_NODE> &out, const VecBasis &source, const NodalRestrictionMapping &mapping) {
  assert(out.vectorSize() == mapping.restrictedInfo());
  assert(source.vectorSize() == mapping.originInfo());
  
  Vector buffer(mapping.restrictedInfo());
  for (VecBasis::const_iterator vecIt = source.begin(); vecIt != source.end(); ++vecIt) {
    mapping.restriction(*vecIt, buffer);
    out << buffer;
  }
  return out;
}

// Input
template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
operator>>(BasisInputStream<DOFS_PER_NODE> &in, VecBasis &sink) {
  sink.dimensionIs(in.size() - in.currentVectorRank(), in.vectorSize());
  return readVectors(in, sink.begin(), sink.end());
}

template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readVectors(BasisInputStream<DOFS_PER_NODE> &in, VecBasis &sink, int countMax) {
  const int count = std::max(std::min(in.size() - in.currentVectorRank(), countMax), 0);
  sink.dimensionIs(count, in.vectorSize());
  return readVectors(in, sink.begin(), sink.end());
}

template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readVectors(BasisInputStream<DOFS_PER_NODE> &in, VecBasis &sink, int countMax, int localSize, int offset) {
  sink.dimensionIs(countMax, in.vectorSize());
  return readVectors(in, sink.begin()+offset, sink.begin()+offset+localSize);
}

template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readRestrictedVectors(BasisInputStream<DOFS_PER_NODE> &in, VecBasis &sink, const NodalRestrictionMapping &mapping) {
  assert(in.vectorSize() == mapping.originInfo());
  
  sink.dimensionIs(in.size() - in.currentVectorRank(), mapping.restrictedInfo());
  Vector buffer(mapping.originInfo());
  for (VecBasis::iterator vecIt = sink.begin(); vecIt != sink.end(); ++vecIt) {
    in >> buffer;
    mapping.restriction(buffer, *vecIt);
  }
  return in;
}

template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readExtendedVectors(BasisInputStream<DOFS_PER_NODE> &in, VecBasis &sink, const NodalRestrictionMapping &mapping) {
  assert(in.vectorSize() == mapping.restrictedInfo());
  
  sink.dimensionIs(in.size() - in.currentVectorRank(), mapping.originInfo());
  Vector buffer(mapping.restrictedInfo());
  for (VecBasis::iterator vecIt = sink.begin(); vecIt != sink.end(); ++vecIt) {
    in >> buffer;
    mapping.extension(buffer, *vecIt);
  }
  return in;
}

template
BasisOutputStream<6> &
operator<<(BasisOutputStream<6> &out, const VecBasis &source);

template
BasisOutputStream<6> &
writeVectors(BasisOutputStream<6> &out, const VecBasis &source, int countMax);

template
BasisOutputStream<6> &
writeExtendedVectors(BasisOutputStream<6> &out, const VecBasis &source, const NodalRestrictionMapping &mapping);

template
BasisOutputStream<6> &
writeRestrictedVectors(BasisOutputStream<6> &out, const VecBasis &source, const NodalRestrictionMapping &mapping);

template
BasisOutputStream<1> &
operator<<(BasisOutputStream<1> &out, const VecBasis &source);

template
BasisOutputStream<1> &
writeVectors(BasisOutputStream<1> &out, const VecBasis &source, int countMax);

template
BasisOutputStream<1> &
writeExtendedVectors(BasisOutputStream<1> &out, const VecBasis &source, const NodalRestrictionMapping &mapping);

template
BasisOutputStream<1> &
writeRestrictedVectors(BasisOutputStream<1> &out, const VecBasis &source, const NodalRestrictionMapping &mapping);


template
BasisInputStream<6> & operator>>(BasisInputStream<6> &in, VecBasis &sink);

template
BasisInputStream<6> & readVectors(BasisInputStream<6> &in, VecBasis &sink, int countMax);

template
BasisInputStream<6> & readVectors(BasisInputStream<6> &in, VecBasis &sink, int countMax, int localSize, int offset);

template
BasisInputStream<6> & readRestrictedVectors(BasisInputStream<6> &in, VecBasis &sink, const NodalRestrictionMapping &mapping);

template
BasisInputStream<6> & readExtendedVectors(BasisInputStream<6> &in, VecBasis &sink, const NodalRestrictionMapping &mapping);

template
BasisInputStream<1> & operator>>(BasisInputStream<1> &in, VecBasis &sink);

template
BasisInputStream<1> & readVectors(BasisInputStream<1> &in, VecBasis &sink, int countMax);

template
BasisInputStream<1> & readVectors(BasisInputStream<1> &in, VecBasis &sink, int countMax, int localSize, int offset);

template
BasisInputStream<1> & readRestrictedVectors(BasisInputStream<1> &in, VecBasis &sink, const NodalRestrictionMapping &mapping);

template
BasisInputStream<1> & readExtendedVectors(BasisInputStream<1> &in, VecBasis &sink, const NodalRestrictionMapping &mapping);

} /* end namespace Rom */
