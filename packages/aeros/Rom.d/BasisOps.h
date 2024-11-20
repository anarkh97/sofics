#ifndef ROM_BASISOPS_H
#define ROM_BASISOPS_H

#include <cassert>

namespace Rom {

template <typename VecSetType, typename VecType, typename OtherVecType>
inline
void
checkDimensions(const VecSetType &basis, const VecType &xFull, const OtherVecType &xReduced) {
  assert(basis.size() == xFull.size());
  assert(basis.numVec() == xReduced.size());
}

template <typename VecSetType, typename VecType, typename OtherVecType>
const OtherVecType &
reduce(const VecSetType &basis, const VecType &xFull, OtherVecType &xReduced) {
  checkDimensions(basis, xFull, xReduced);

  for (int i = 0; i < basis.numVec(); ++i) {
    //xReduced[i] = basis[i]*xFull;
    //workaround lack of const-correctness in GenDistrVector class
    xReduced[i] = const_cast<typename VecSetType::VecType &>(basis[i]) * const_cast<VecType &>(xFull);
  }

  return xReduced;
}

template <typename VecSetType, typename VecType, typename OtherVecType>
const OtherVecType &
expand(const VecSetType &basis, const VecType &xReduced, OtherVecType &xFull) {
  checkDimensions(basis, xFull, xReduced);

  xFull.zero();
  for (int i = 0; i < basis.numVec(); ++i) {
    xFull += xReduced[i]*basis[i];
  }

  return xFull; 
}

template <typename VecSetType, typename VecType, typename OtherVecType>
const OtherVecType &
reduceAdd(const VecSetType &basis, const VecType &xFull, OtherVecType &xReduced) {
  checkDimensions(basis, xFull, xReduced);

  for (int i = 0; i < basis.numVec(); ++i) {
    xReduced[i] += basis[i] * xFull;
  }

  return xReduced;
}

template <typename VecSetType, typename VecType, typename OtherVecType>
const OtherVecType &
expandAdd(const VecSetType &basis, const VecType &xReduced, OtherVecType &xFull) {
  checkDimensions(basis, xFull, xReduced);

  for (int i = 0; i < basis.numVec(); ++i) {
    xFull += xReduced[i]*basis[i];
  }

  return xFull; 
}

} /* end namespace Rom */

#endif /* ROM_BASISOPS_H */
