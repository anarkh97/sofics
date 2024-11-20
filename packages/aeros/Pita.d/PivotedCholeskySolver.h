#ifndef PITA_PIVOTEDCHOLESKYSOLVER_H
#define PITA_PIVOTEDCHOLESKYSOLVER_H

#include "Fwk.h"
#include "RankDeficientSolver.h"
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SymFullMatrix.h>
#include "SimpleBuffer.h"

namespace Pita {

// Solves by a direct method the linear algebra problem of the form A * x = b
// where A is a symmetric positive semi-definite matrix
// with diagonal rescaling and using a Cholesky algorithm with complete pivoting.
// Only the lower triangular part (in C indexing) of A is accessed.
// The factorization has the form: P^T * S * A * S * P = U^T * U 
// where U is upper triangular, P is a permutation and S is diagonal.

class PivotedCholeskySolver : public RankDeficientSolver {
public:
  EXPORT_PTRINTERFACE_TYPES(PivotedCholeskySolver);

  const FullSquareMatrix & choleskyFactor() const { return choleskyFactor_; } // Only lower triangular part is relevant
  virtual const Vector & solution(Vector & rhs) const; // In-place solution: rhs modified

  void matrixIs(const SymFullMatrix & matrix);
  virtual void transposedMatrixIs(const FullSquareMatrix & matrix); // Use only lower triangular part of transposed matrix
  virtual void orderingIs(Ordering o);
  
  static Ptr New(double tolerance = -1.0) {
    return new PivotedCholeskySolver(tolerance);
  }

protected:
  explicit PivotedCholeskySolver(double tolerance);

  void performFactorization();

private:
  FullSquareMatrix choleskyFactor_;

  SimpleBuffer<double> scaling_;
}; 

} // end namespace Pita

#endif /* PITA_PIVOTEDCHOLESKYSOLVER_H */
