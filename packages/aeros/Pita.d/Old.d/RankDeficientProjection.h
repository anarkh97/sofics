#ifndef RANKDEFICIENTPROJECTION_H
#define RANKDEFICIENTPROJECTION_H

#include "DynamStateSet.h"

#include <Pita.d/NearSymmetricSolver.h>
#include <Math.d/FullSquareMatrix.h>

namespace Pita { namespace Old {

class RankDeficientProjection {
public:
  typedef double Scalar;
  typedef DynamState<Scalar> State;
  typedef DynamStateSet<Scalar> StateSet;

  //** Raw basis characteristics
  int vectorSize() const { return rawBasis_.vectorSize(); }
  int stateCount() const { return rawBasis_.numStates(); }
  const State & state(int i) const { return rawBasis_[i]; }

  void vectorSizeIs(int vectorSize);
  void lastStateIs(const State &);
  void lastStateSetIs(const StateSet &);
  void lastStateSetIs(int stateCount, Scalar * data);
  void stateDel();

  //** Orthogonalization
  double relativeTolerance() const { return solver_->tolerance(); }
  void relativeToleranceIs(double relTol) { solver_->toleranceIs(relTol); }

  enum OGStatus {
    TO_DO,
    DONE
  };
  OGStatus orthoStatus() const { return orthoStatus_; }

  // Orthogonalize with respect to provided metric
  // to determine an independent basis (independentBasis is modified for efficiency)
  void metricIs(const GenSparseMatrix<Scalar> * K, const GenSparseMatrix<Scalar> * M,
                StateSet & independentBasis);

  //**Projection
  // Projection of inputState on finalSet with respect to the metric originally provided
  // (outputState is modified for efficiency)
  const State & projection(const StateSet & finalSet, const State & inputState, State & outputState) const;

  //** Constructors
  explicit RankDeficientProjection(int vectorSize, double relTol);
  
  // Copy and assignment to simulate value semantics
  RankDeficientProjection(const RankDeficientProjection & other);
  RankDeficientProjection & operator=(const RankDeficientProjection & other);

private:
  OGStatus orthoStatus_;

  StateSet rawBasis_;
  StateSet metricBasis_;
  
  FullSquareMatrix normalMatrix_;
  Pita::NearSymmetricSolver::Ptr solver_;
};

} /* end namespace Old */ } /* end namespace Pita */

#endif /* RANKDEFICIENTPROJECTION_H */
