#include "NearSymmetricSolver.h"

#include <algorithm>

// Kernel routines
#include "Lapack32.d/dsyequb.h"

extern "C" {
  // Blas: Rank-1 update A += alpha*x*y'
  void _FORTRAN(dger)(const int* m, const int* n, const double* alpha,
                      const double* x, const int* incx, const double* y, const int* incy,
                      double* a, const int* lda);

  // Blas: Triangular system resolution
  void _FORTRAN(dtrsv)(const char * uplo, const char * trans, const char * diag, const int * n,
                       const double * a, const int * lda, double * x, const int * incx);

  // Lapack: Row pivoting
  void _FORTRAN(dlaswp)(const int* n, double* a, const int* lda, const int* k1, const int* k2,
                        const int* ipiv, const int* incx);

  // Blas: Column (vector) interchange
  void _FORTRAN(dswap)(const int* n, double* x, const int* incx, double* y, const int* incy);

  // Blas: Vector scaling
  void _FORTRAN(dscal)(const int* n, const double* da, double* dx, const int* incx);

  // Lapack: Perform scaling
  void _FORTRAN(dlaqge)(const int* m, const int* n, double* a, const int* lda,
                        const double* r, const double* c, const double* rowcnd, const double* colcnd,
                        const double* amax, char* equed);
}

namespace Pita {

// Constructor

NearSymmetricSolver::NearSymmetricSolver(double tol) :
  RankDeficientSolver(tol),
  transposedMatrix_(),
  pivots_(),
  scaling_()
{}

// Factor

void
NearSymmetricSolver::transposedMatrixIs(const FullSquareMatrix & tm) {
  setTransposedMatrix(tm);
  setMatrixSize(transposedMatrix_.dim());
    
  // Rescale matrix
  scaling_.sizeIs(matrixSize());

  double scond, amax;

  // Hand-made routine to replace dsyequb
  equilibrateSym(getMatrixSize(), transposedMatrix_.data(), tolerance(), scaling_.array(), &scond, &amax);

  char equed;
  _FORTRAN(dlaqge)(&getMatrixSize(), &getMatrixSize(), transposedMatrix_.data(), &getMatrixSize(),
                   scaling_.array(), scaling_.array(), &scond, &scond, &amax, &equed); 

  switch (equed) {
    case 'N':
      setRescalingStatus(NO_RESCALING);
      break;
    case 'R':
      setRescalingStatus(ROW_RESCALING);
      break;
    case 'B':
      setRescalingStatus(SYMMETRIC_RESCALING);
      break;
    default:
      throw Fwk::InternalException("Invalid rescaling status");
  }

  // Initialize permutation
  setVectorSize(matrixSize());
  setOrdering(PERMUTED);
  getFactorPermutation().sizeIs(vectorSize());
  for (int i = 0; i < matrixSize(); ++i) {
    getFactorPermutation()[i] = i + 1;
  }
  pivots_.sizeIs(matrixSize());

  // Numerical constants
  const int int_one = 1;
  const double minus_one = -1.0;

  double absTol = 0.0;

  int k;
  for (k = 0; k < matrixSize(); ++k) {
    // Find largest pivot on diagonal
    int p;
    double piv_val;
    {
      p = k;
      const double * piv_addr = &transposedMatrix()[p][p];
      piv_val = *piv_addr;

      const int piv_stride = matrixSize() + 1;
      for (int i = k + 1; i < matrixSize(); ++i) {
        piv_addr += piv_stride;
        if (*piv_addr > piv_val) {
          p = i;
          piv_val = *piv_addr;
        }
      }
    }

    pivots_[k] = p + 1;

    // Perform symmetric permutation
    std::swap(getFactorPermutation()[k], getFactorPermutation()[p]);
    std::swap(scaling_.array()[k], scaling_.array()[p]);

    _FORTRAN(dswap)(&getMatrixSize(), &transposedMatrix_[k][0], &int_one, &transposedMatrix_[p][0], &int_one);
    _FORTRAN(dswap)(&getMatrixSize(), &transposedMatrix_[0][k], &getMatrixSize(), &transposedMatrix_[0][p], &getMatrixSize());

    // Determine absolute tolerance (First iteration only)
    if (k == 0) {
      absTol = piv_val * tolerance();
    }
    
    // Check for singularity
    if (piv_val <= absTol) {
      break;
    }

    // Compute column of L
    const double pivot_inverse = 1.0 / piv_val;
    const int remainder_size = matrixSize() - (k + 1);
    _FORTRAN(dscal)(&remainder_size, &pivot_inverse, &transposedMatrix_[k][k+1], &int_one);

    // Update remaining block of A
    _FORTRAN(dger)(&remainder_size, &remainder_size, &minus_one,
                   &transposedMatrix_[k][k+1], &int_one, &transposedMatrix_[k+1][k], &getMatrixSize(),
                   &transposedMatrix_[k+1][k+1], &getMatrixSize());
  }

  setFactorRank(k);
}

// Solve

const Vector &
NearSymmetricSolver::solution(Vector & rhs) const {
  if (rhs.size() != vectorSize()) {
    throw Fwk::RangeException("Size mismatch"); 
  }

  // 1) rhs <- P * rhs
  const int int_one = 1;
  if (ordering() == PERMUTED) {
    _FORTRAN(dlaswp)(&int_one, rhs.data(), &getVectorSize(), &int_one, &getFactorRank(),
                     pivots_.array(), &int_one);
  }

  // 2) rhs <- S^{-1} * rhs
  if (rescalingStatus() != NO_RESCALING) {
    for (int i = 0; i < factorRank(); ++i) {
      rhs[i] *= scaling_[i];
    }
  }

  // 3) rhs <- L^{-1} * rhs
  const char lower = 'L';
  const char non_trans = 'N';
  const char unit_diag = 'U';
  _FORTRAN(dtrsv)(&lower, &non_trans, &unit_diag, &getFactorRank(), transposedMatrix().data(),
                  &getMatrixSize(), rhs.data(), &int_one);

  // 4) rhs <- U^{-1} * rhs
  const char upper = 'U';
  const char non_unit_diag = 'N';
  _FORTRAN(dtrsv)(&upper, &non_trans, &non_unit_diag, &getFactorRank(), transposedMatrix().data(),
                  &getMatrixSize(), rhs.data(), &int_one);
  
  // 5) rhs <- S^{-1} * rhs
  if (rescalingStatus() == SYMMETRIC_RESCALING) {
    for (int i = 0; i < factorRank(); ++i) {
      rhs[i] *= scaling_[i];
    }
  }

  // 6) Pad with zeros
  std::fill(rhs.data() + getFactorRank(), rhs.data() + getVectorSize(), 0.0);

  // 7) rhs <- P^T * rhs
  const int int_minus_one = -1;
  if (ordering() == PERMUTED) {
    _FORTRAN(dlaswp)(&int_one, rhs.data(), &getVectorSize(), &int_one, &getFactorRank(),
                     pivots_.array(), &int_minus_one);
  }

  return rhs;
}

// Reorder

void
NearSymmetricSolver::orderingIs(Ordering o) {
  if (ordering() == o)
    return;

  if (o == COMPACT) {
    setVectorSize(factorRank());
    getFactorPermutation().sizeIs(vectorSize());

    for (int i = 0; i < factorRank(); ++i) {
      getFactorPermutation()[i] = i + 1; // Fortran numbering
    }
  } else {
    throw Fwk::RangeException("Forbidden transition");
  }

  setOrdering(o);
}

void
NearSymmetricSolver::toleranceIs(double tol) {
  // TODO Extend/shrink factorization ?
  setTolerance(tol);
}

} /* end namespace Pita */
