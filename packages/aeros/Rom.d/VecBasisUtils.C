#include "VecBasisUtils.h"
#include "SimpleBuffer.h"

#include <Math.d/Vector.h>

#include <Utils.d/linkfc.h>

#include <stdexcept>
#include <complex>
#include <cassert>

extern "C" {
  // Matrix-matrix multiplication: C := alpha * op(A) * op(B) + beta * C
  void _FORTRAN(dgemm)(const char *trans, const char *transb,
                       const int *m, const int *n, const int *k,
                       const double *alpha, const double *a, const int *lda,
                       const double *b, const int *ldb,
                       const double *beta, double *c, const int *ldc);
}

namespace Rom {

namespace { // anonymous

  template <typename Scalar>
  inline
  void
  check_dimensions(const GenVecBasis<Scalar> &targetPod,
                   const GenVecBasis<Scalar> &originPod,
                   const GenVecBasis<Scalar> &originProjection,
                   GenVecBasis<Scalar> &result) {
    assert(targetPod.vectorSize()        == originPod.vectorSize());
    assert(originProjection.vectorSize() == result.vectorSize());
    assert(targetPod.vectorCount()       == result.vectorCount());
    assert(originPod.vectorCount()       == originProjection.vectorCount());
  }

  template <typename Scalar>
  inline
  Scalar *
  get_buffer(GenVecBasis<Scalar> &basis) {
    return &basis[0][0];
  }
  
  template <typename Scalar>
  inline
  const Scalar *
  get_buffer(const GenVecBasis<Scalar> &basis) {
    return &basis[0][0];
  }

} // end anonymous namespace

template <>
const GenVecBasis<double> &
combine_projections<double>(const GenVecBasis<double> &targetPod,
                            const GenVecBasis<double> &originPod,
                            const GenVecBasis<double> &originProjection,
                            GenVecBasis<double> &result) {
  check_dimensions(targetPod, originPod, originProjection, result);

  const int fullVecSize    = targetPod.vectorSize();
  const int reducedVecSize = originProjection.vectorSize();
  const int tempRows       = originPod.vectorCount();
  const int tempCols       = targetPod.vectorCount();

  SimpleBuffer<double> tempBuffer(tempRows * tempCols);

  const double zero = 0.0;
  const double one  = 1.0;

  const char notrans = 'N';
  const char trans   = 'T';

  // temp := originPod^{T} * targetPod
  _FORTRAN(dgemm)(&trans, &notrans,
                  &tempRows, &tempCols, &fullVecSize,
                  &one, get_buffer(originPod), &fullVecSize,
                  get_buffer(targetPod), &fullVecSize,
                  &zero, tempBuffer.array(), &tempRows);

  // result := originProjection * temp
  _FORTRAN(dgemm)(&notrans, &notrans,
                  &reducedVecSize, &tempCols, &tempRows,
                  &one, get_buffer(originProjection), &reducedVecSize,
                  tempBuffer.array(), &tempRows,
                  &zero, get_buffer(result), &reducedVecSize);

  return result;
}

template <>
const GenVecBasis<complex<double> > &
combine_projections<complex<double> >(const GenVecBasis<complex<double> > &targetPod,
                                           const GenVecBasis<complex<double> > &originPod,
                                           const GenVecBasis<complex<double> > &originProjection,
                                           GenVecBasis<complex<double> > &result) {
  throw std::logic_error("Not implemented");
}

} /* end namespace Rom */
