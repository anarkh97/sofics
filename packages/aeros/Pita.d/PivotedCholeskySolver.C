#include "PivotedCholeskySolver.h"
#include <algorithm>

// Kernel routines
#include "Lapack32.d/dpstrf.h"
#include "Lapack32.d/dsyequb.h"

extern "C" {
  // Blas: Backward/Forward substitution
  void _FORTRAN(dtrsv)(const char * uplo, const char * trans, const char * diag, const int * n,
                       const double * a, const int * lda, double * x, const int * incx);
  
  // Blas: Vector scaling
  void _FORTRAN(dscal)(const int* n, const double* da, double* dx, const int* incx);

  // Lapack: Perform symmetric scaling
  void _FORTRAN(dlaqsy)(const char* uplo, const int* n, double* a, const int* lda,
                        const double* s, const double* scond, const double* amax, char* equed);
}

namespace Pita {

// Constructor

PivotedCholeskySolver::PivotedCholeskySolver(double tolerance) :
  RankDeficientSolver(tolerance), 
  choleskyFactor_(),
  scaling_()
{}

// Factor

void
PivotedCholeskySolver::matrixIs(const SymFullMatrix & matrix) {
  int newSize = matrix.dim();
  
  // Fill in lower triangular part of choleskyFactor
  choleskyFactor_.setSize(newSize);
  const double * sourceData = const_cast<SymFullMatrix &>(matrix).data();
  for (int i = 0; i < newSize; ++i) {
    const double * sourceDataEnd = (sourceData + i) + 1;
    std::copy(sourceData, sourceDataEnd, choleskyFactor_[i]);
    sourceData = sourceDataEnd;
  }

  performFactorization();
}

void
PivotedCholeskySolver::transposedMatrixIs(const FullSquareMatrix & matrix) {
  choleskyFactor_.copy(matrix);
  performFactorization();
}

void
PivotedCholeskySolver::performFactorization() {
  setMatrixSize(choleskyFactor_.dim());
  
  const char upper = 'U';   // Lower triangular in C indexing == upper triangular in Fortran indexing 
  SimpleBuffer<double> workspace(3 * matrixSize());
  
  // Rescale matrix
  scaling_.sizeIs(matrixSize());
  int info;
  double scond, amax;

  // Hand-made routine to replace dsyequb
  equilibrateSym(getMatrixSize(), choleskyFactor_.data(), tolerance(), scaling_.array(), &scond, &amax);
  
  char equed;
  _FORTRAN(dlaqsy)(&upper, &getMatrixSize(), choleskyFactor_.data(), &getMatrixSize(),
                   scaling_.array(), &scond, &amax, &equed);

  switch (equed) {
    case 'N':
      setRescalingStatus(NO_RESCALING);
      break;
    case 'Y':
      setRescalingStatus(SYMMETRIC_RESCALING);
      break;
    default:
      throw Fwk::InternalException("Invalid rescaling status");
  }
  
  // Initialize permutation
  setVectorSize(matrixSize());
  setOrdering(PERMUTED);
  getFactorPermutation().sizeIs(vectorSize());

  setFactorRank(0);
  if (matrixSize() == 0) {
    return;
  }

  // Determine the lower bound for the pivot
  double absTol;
  if (tolerance() >= 0) {
    double maxDiagVal;
    {
      const double * diag_addr = choleskyFactor_.data();
      maxDiagVal = *diag_addr;

      const int diag_stride = matrixSize() + 1;
      for (int i = 1; i < matrixSize(); ++i) {
        diag_addr += diag_stride;
        if (*diag_addr > maxDiagVal) {
          maxDiagVal = *diag_addr;
        }
      }
    }
    absTol = maxDiagVal * tolerance();
  } else {
    absTol = tolerance();
  }

  // Perform factorization
  int numericalRank;

  _FORTRAN(dpstrf)(&upper, &getMatrixSize(), choleskyFactor_.data(), &getMatrixSize(),
                   getFactorPermutation().array(), &numericalRank, &absTol, workspace.array(), &info);
  assert(info >= 0);

  // Reorder scaling coefficients
  std::copy(scaling_.array(), scaling_.array() + matrixSize(), workspace.array());
  for (int i = 0; i < matrixSize(); ++i) {
    scaling_[i] = workspace[factorPermutation(i)];
  }

  setFactorRank(numericalRank);
}

// Reordering

void
PivotedCholeskySolver::orderingIs(Ordering o) {
  if (ordering() == 0)
    return;

  if (o == COMPACT) {
    setVectorSize(factorRank());
    getFactorPermutation().sizeIs(vectorSize());

    for (int i = 0; i < factorRank(); ++i) {
      getFactorPermutation()[i] = i + 1; // Fortran numbering
    }
  }

  setOrdering(o);
}

// Solve

const Vector &
PivotedCholeskySolver::solution(Vector & rhs) const {
  if (rhs.size() != vectorSize()) {
    throw Fwk::RangeException("Size mismatch"); 
  }

  // 1) rhs <- P * rhs
  SimpleBuffer<double> perm_vec(factorRank());
  
  double * rhs_data;
  if (ordering() == PERMUTED) {
    // Pointers at the beginning / end of permutation array
    const int * fp_ptr = getFactorPermutation().array();
    const int * fp_ptr_end = getFactorPermutation().array() + factorRank();

    for (double * perm_vec_ptr = perm_vec.array(); fp_ptr != fp_ptr_end; ++perm_vec_ptr, ++fp_ptr) {
      *perm_vec_ptr = rhs[*fp_ptr - 1]; // Offset to get C indexing from Fortran indexing
    }
    rhs_data = perm_vec.array();
  } else {
    rhs_data = rhs.data();
  }
  
  // 2) rhs <- S^{-1} * rhs
  if (rescalingStatus() != NO_RESCALING) {
    for (int i = 0; i < factorRank(); ++i) {
      rhs_data[i] *= scaling_[i];
    }
  }

  // 3) rhs <- U^{-1} * U^{-T} * rhs
  const char upper = 'U';             // Lower triangular in C indexing == upper triangular in Fortran indexing
  const char forwardSubTrans = 'T';   // To solve U^{-T}
  const char backwardSubTrans = 'N';  // To solve U^{-1}
  const char diag = 'N';              // CholeskyFactor has non-unit diagonal
  const int incx = 1;                 // Rhs vector elements are contiguous

  // Forward substitution: rhs <- U^{-T} * rhs
  _FORTRAN(dtrsv)(&upper, &forwardSubTrans, &diag, &getFactorRank(), choleskyFactor().data(),
                  &getMatrixSize(), rhs_data, &incx);
  // Backward substitution: rhs <- U^{-1} * rhs
  _FORTRAN(dtrsv)(&upper, &backwardSubTrans, &diag, &getFactorRank(), choleskyFactor().data(),
                  &getMatrixSize(), rhs_data, &incx);

  // 4) rhs <- S^{-1} * rhs
  if (rescalingStatus() == SYMMETRIC_RESCALING) {
    for (int i = 0; i < factorRank(); ++i) {
      rhs_data[i] *= scaling_[i];
    }
  }
  
  // 5) rhs <- P^T * rhs 
  if (ordering() == PERMUTED) {
    // Replace rhs by solution
    rhs.zero();
    
    // Pointers at the beginning / end of permutation array
    const int * fp_ptr = getFactorPermutation().array();
    const int * fp_ptr_end = getFactorPermutation().array() + factorRank();
    for (const double * perm_vec_ptr = perm_vec.array(); fp_ptr != fp_ptr_end; ++perm_vec_ptr, ++fp_ptr) {
      rhs[*fp_ptr - 1] = *perm_vec_ptr; // Offset to get C indexing from Fortran indexing
    }
  }
  
  return rhs;
}

} /* end namespace Pita */
