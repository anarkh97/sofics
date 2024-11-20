#ifndef PITA_LEASTSQUARESOLVER_H
#define PITA_LEASTSQUARESOLVER_H

#include "Fwk.h"

#include "RankDeficientSolver.h"

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>

#include "SimpleBuffer.h"

namespace Pita {

class LeastSquareSolver : public RankDeficientSolver {
public:
  EXPORT_PTRINTERFACE_TYPES(LeastSquareSolver);
  
  const FullSquareMatrix & transposedMatrix() const { return transposedMatrix_; }
  
  virtual void transposedMatrixIs(const FullSquareMatrix & transposedMatrix); // The matrix is copied
  virtual void toleranceIs(double tol);
  virtual void orderingIs(Ordering o); // COMPACT reordering not implemented

  // In place solution, argument rhs is modified
  virtual const Vector & solution(Vector & rhs) const;

  static Ptr New(double tolerance) {
    return new LeastSquareSolver(tolerance);
  }

protected:
  explicit LeastSquareSolver(double tol);
  
  void updateFactorRank();

private:
  FullSquareMatrix transposedMatrix_; // To have column-major ordering

  SimpleBuffer<double> tau_;
  mutable SimpleBuffer<double> workspace_;

  DISALLOW_COPY_AND_ASSIGN(LeastSquareSolver);
};

} /* end namespace Pita */

#endif /* PITA_LEASTSQUARESOLVER_H */
