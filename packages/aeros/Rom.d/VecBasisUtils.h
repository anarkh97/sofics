#ifndef ROM_VECBASISUTILS_H
#define ROM_VECBASISUTILS_H

#include "VecBasis.h"
#include "NodalRestrictionMapping.h"

namespace Rom {

template <typename Scalar>
void
restrict_basis(const NodalRestrictionMapping &mapping,
               const GenVecBasis<Scalar> &input, GenVecBasis<Scalar> &result) {
  result.dimensionIs(input.vectorCount(), mapping.restrictedInfo());
  
  typename GenVecBasis<Scalar>::iterator jt = result.begin();
  const typename GenVecBasis<Scalar>::const_iterator itEnd = input.end();
  for (typename GenVecBasis<Scalar>::const_iterator it = input.begin(); it != itEnd; ++it) {
    mapping.restriction(*it, *jt++);
  }
}

// Computes result^{T} := targetPod^{T} * originPod * originProjection^{T} 
template <typename Scalar>
const GenVecBasis<Scalar> &
combine_projections(const GenVecBasis<Scalar> &targetPod,
                    const GenVecBasis<Scalar> &originPod,
                    const GenVecBasis<Scalar> &originProjection,
                    GenVecBasis<Scalar> &result);

} /* end namespace Rom */

#endif /* ROM_VECBASISUTILS_H */
