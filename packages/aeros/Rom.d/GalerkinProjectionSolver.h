#ifndef ROM_GALERKINPROJECTIONSOLVER_H
#define ROM_GALERKINPROJECTIONSOLVER_H

#include "PodProjectionSolver.h"

#include "CholeskyUtils.h"

#include <Math.d/FullSquareMatrix.h>

#include <cstddef>
#include <cassert>

namespace Rom {

template <typename Scalar>
class GenGalerkinProjectionSolver : public GenDBSparsePodProjectionSolver<Scalar> {
public:
  GenGalerkinProjectionSolver(const Connectivity *cn, const DofSetArray *dsa, const ConstrainedDSA *c_dsa,
      bool pivotFlag=false);
  ~GenGalerkinProjectionSolver();
private:
  GenFullSquareMatrix<Scalar> reducedMatrix_;
  bool pivotFlag_;
  int *ipiv_;
 
  // Overriden 
  virtual void resetSolver(int vCount, int vSize) override;
  virtual void assembleAndFactorReducedSystem(double Mcoef) override;
  virtual double getReducedRhsNorm() const override;
  virtual void solveReducedSystem(GenVector<Scalar> &) override;
 
  // Implementation 
  void performFactor();
  void performSolve();
};

template <typename Scalar>
GenGalerkinProjectionSolver<Scalar>::GenGalerkinProjectionSolver(const Connectivity *cn,
                                                                 const DofSetArray *dsa,
                                                                 const ConstrainedDSA *c_dsa,
                                                                 bool pivotFlag):
  GenDBSparsePodProjectionSolver<Scalar>(cn, dsa, c_dsa),
  reducedMatrix_(),
  pivotFlag_(pivotFlag),
  ipiv_(NULL)
{}

template <typename Scalar>
GenGalerkinProjectionSolver<Scalar>::~GenGalerkinProjectionSolver()
{
  if(ipiv_) delete [] ipiv_;
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::resetSolver(int vCount, int) {
  reducedMatrix_.setSize(vCount);
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::assembleAndFactorReducedSystem(double Mcoef) {
  for (int row = 0; row < this->basisSize(); ++row) {
    const GenVector<Scalar> &action = this->lastReducedMatrixAction()[row];
    for (int col = row; col < this->basisSize(); ++col) {
      reducedMatrix_[row][col] = action * this->projectionBasis()[col];
    }
    reducedMatrix_[row][row] += Mcoef;
  }

  performFactor();
}

template <typename Scalar>
double
GenGalerkinProjectionSolver<Scalar>::getReducedRhsNorm() const {
  return this->lastReducedSolution().norm();
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::solveReducedSystem(GenVector<Scalar> &rhs) {
  if(pivotFlag_) {
    if(!ipiv_)
      throw std::runtime_error("You cannot have pivoting if you did not compute it.");
    ldlt_solve_upper(reducedMatrix_, rhs.data(), ipiv_);
  }
  else
    cholesky_solve_upper(reducedMatrix_, rhs.data());
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::performFactor() {

  if(pivotFlag_) {
    if(!ipiv_)
      throw std::runtime_error("You cannot have pivoting if you did not compute it.");
    ldlt_factor_upper(reducedMatrix_, ipiv_);
  }
  else {
    cholesky_factor_upper(reducedMatrix_);
  }
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::performSolve() {

  if(pivotFlag_) {
    ldlt_solve_upper(reducedMatrix_, this->getReducedSolution().data(), ipiv_);
  }
  else {
    cholesky_solve_upper(reducedMatrix_, this->getReducedSolution().data());
  }
}

typedef GenGalerkinProjectionSolver<double> GalerkinProjectionSolver;

} /* end namespace Rom */

#endif /* ROM_GALERKINPROJECTIONSOLVER_H */
