#ifndef ROM_VECBASISFILE_H
#define ROM_VECBASISFILE_H

#include "VecBasis.h"

#include <Math.d/Vector.h>

namespace Rom {

template<int DOFS_PER_NODE> class BasisInputStream;
template<int DOFS_PER_NODE> class BasisOutputStream;

class NodalRestrictionMapping;

// Output
// Output the full content of the basis
template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
operator<<(BasisOutputStream<DOFS_PER_NODE> &, const VecBasis &);

template<int DOFS_PER_NODE>
inline
BasisOutputStream<DOFS_PER_NODE> &
writeVectors(BasisOutputStream<DOFS_PER_NODE> &out, const VecBasis &source) {
  return out << source;
}

// Output a partial content of the basis
template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeVectors(BasisOutputStream<DOFS_PER_NODE> &, const VecBasis &, int countMax);

// Output the full extended content of the basis
template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeExtendedVectors(BasisOutputStream<DOFS_PER_NODE> &, const VecBasis &, const NodalRestrictionMapping &);

// Output the full extended content of the basis
template<int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE> &
writeRestrictedVectors(BasisOutputStream<DOFS_PER_NODE> &, const VecBasis &, const NodalRestrictionMapping &);

// Input
// Reset basis with the full content of the stream (inverse of operator<<)
template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
operator>>(BasisInputStream<DOFS_PER_NODE> &, VecBasis &);

template<int DOFS_PER_NODE>
inline
BasisInputStream<DOFS_PER_NODE> &
readVectors(BasisInputStream<DOFS_PER_NODE> &in, VecBasis &target) {
  return in >> target;
}

// Reset basis with a partial content of the stream (inverse of writeVectors)
template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readVectors(BasisInputStream<DOFS_PER_NODE> &, VecBasis &, int countMax);

// Reset basis with a partial content of the stream (inverse of writeVectors)
// or add new vectors to existing basis
template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readVectors(BasisInputStream<DOFS_PER_NODE> &, VecBasis &, int countMax, int localsize, int offset);

// Reset basis with the extended content of the stream (inverse of writeRestrictedVectors)
template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readExtendedVectors(BasisInputStream<DOFS_PER_NODE> &, VecBasis &, const NodalRestrictionMapping &);

// Reset basis with the restricted content of the stream (inverse of writeExtendedVectors)
template<int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE> &
readRestrictedVectors(BasisInputStream<DOFS_PER_NODE> &, VecBasis &, const NodalRestrictionMapping &);

} /* end namespace Rom */

#endif /* ROM_VECBASISFILE_H */
