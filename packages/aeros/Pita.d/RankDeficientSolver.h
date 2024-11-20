#ifndef PITA_RANKDEFICIENTSOLVER_H
#define PITA_RANKDEFICIENTSOLVER_H

#include "Fwk.h"
#include "SimpleBuffer.h"
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>

namespace Pita {

class RankDeficientSolver : public Fwk::PtrInterface<RankDeficientSolver> {
public:
  EXPORT_PTRINTERFACE_TYPES(RankDeficientSolver);

  double tolerance() const { return tolerance_; } // Relative tolerance
  int matrixSize() const { return matrixSize_; }
  int factorRank() const { return factorRank_; }
  
  int vectorSize() const { return vectorSize_; }
  
  enum Ordering {
    COMPACT = 0,
    PERMUTED
  };
  
  enum RescalingStatus {
    NO_RESCALING,
    ROW_RESCALING,
    SYMMETRIC_RESCALING
  };

  RescalingStatus rescalingStatus() const { return rescalingStatus_; }
  
  Ordering ordering() const { return ordering_; }
  int factorPermutation(int index) const { return factorPermutation_[index] - 1; } // 0 <= index <= factorRank

  virtual const Vector & solution(Vector & rhs) const = 0; // In-place solution: rhs modified

  virtual void toleranceIs(double tol) { setTolerance(tol); }
  virtual void transposedMatrixIs(const FullSquareMatrix & tm) = 0;
  virtual void orderingIs(Ordering o) = 0;

protected:
  explicit RankDeficientSolver(double tol) :
    tolerance_(tol),
    matrixSize_(0),
    factorRank_(0),
    ordering_(COMPACT),
    factorPermutation_(),
    rescalingStatus_(NO_RESCALING)
  {}

  const double & getTolerance() const { return tolerance_; }
  const int & getMatrixSize() const { return matrixSize_; }
  const int & getFactorRank() const { return factorRank_; }
  const int & getVectorSize() const { return vectorSize_; }
  
  const SimpleBuffer<int> & getFactorPermutation() const { return factorPermutation_; }
  SimpleBuffer<int> & getFactorPermutation() { return factorPermutation_; }

  void setMatrixSize(int ms) { matrixSize_ = ms; }
  void setTolerance(double tol) { tolerance_ = tol; }
  void setFactorRank(int fr) { factorRank_ = fr; }
  void setVectorSize(int vs) { vectorSize_ = vs; }
  void setOrdering(Ordering ps) { ordering_ = ps; }
  void setRescalingStatus(RescalingStatus rs) { rescalingStatus_ = rs; }

private:
  double tolerance_;
  int matrixSize_;
  int factorRank_;
  
  int vectorSize_;
  Ordering ordering_;
  SimpleBuffer<int> factorPermutation_; // Indices starting at 1 (Fortran convention)

  RescalingStatus rescalingStatus_;
  
  DISALLOW_COPY_AND_ASSIGN(RankDeficientSolver);
};

} /* end namespace Pita */

#endif /* PITA_RANKDEFICIENTSOLVER_H */
