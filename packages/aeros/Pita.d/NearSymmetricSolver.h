#ifndef PITA_NEARSYMMETRICSOLVER_H
#define PITA_NEARSYMMETRICSOLVER_H

#include "Fwk.h"
#include "RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

namespace Pita {

// Solves by a direct method the linear algebra problem of the form A * x = b
// where A is a near-symmetric positive semi-definite matrix
// (Such as the one obtained by forming explicitely A = V^T * Q * V with Q symmetric positive definite)
// with diagonal rescaling and using Gaussian elimination with complete symmetric pivoting.
// The factorization has the form: P^T * S * A * S * P = L * U 
// where L is unit lower triangular, U is upper triangular, P is a permutation and S is diagonal.

class NearSymmetricSolver : public RankDeficientSolver {
public:
  EXPORT_PTRINTERFACE_TYPES(NearSymmetricSolver);
  
  const FullSquareMatrix & transposedMatrix() const { return transposedMatrix_; } // The whole matrix is used
  virtual const Vector & solution(Vector & rhs) const; // In-place solution: rhs modified

  virtual void toleranceIs(double tol);
  virtual void transposedMatrixIs(const FullSquareMatrix & tm); 
  virtual void orderingIs(Ordering o);

  static Ptr New(double tolerance) {
    return new NearSymmetricSolver(tolerance);
  }

protected:
  explicit NearSymmetricSolver(double tol);

  void setTransposedMatrix(const FullSquareMatrix & tm) { transposedMatrix_.copy(tm); } 

private:
  FullSquareMatrix transposedMatrix_;

  SimpleBuffer<int> pivots_;

  SimpleBuffer<double> scaling_;
};

} /* end namespace Pita */

#endif /* PITA_NEARSYMMETRICSOLVER_H */
