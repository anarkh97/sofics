#include "RankDeficientProjection.h"

#include <Pita.d/NearSymmetricSolver.h>

namespace Pita { namespace Old {

RankDeficientProjection::RankDeficientProjection(int vectorSize, double relTol) :
  orthoStatus_(TO_DO),
  rawBasis_(vectorSize),
  metricBasis_(vectorSize),
  normalMatrix_(),
  solver_(Pita::NearSymmetricSolver::New(relTol))
{}

RankDeficientProjection::RankDeficientProjection(const RankDeficientProjection & other) :
  orthoStatus_(other.orthoStatus_),
  rawBasis_(other.rawBasis_),
  metricBasis_(other.metricBasis_),
  normalMatrix_(other.normalMatrix_),
  solver_(Pita::NearSymmetricSolver::New(other.solver_->tolerance()))
{}

RankDeficientProjection &
RankDeficientProjection::operator=(const RankDeficientProjection & other) {
  orthoStatus_ = other.orthoStatus_;
  rawBasis_ = other.rawBasis_;
  metricBasis_ = other.metricBasis_;
  normalMatrix_ = other.normalMatrix_;
  solver_ = Pita::NearSymmetricSolver::New(other.solver_->tolerance());

  return *this;
}

void
RankDeficientProjection::vectorSizeIs(int vSize) {
  orthoStatus_ = TO_DO;
  rawBasis_.reset(vSize, 0);
  metricBasis_.reset(vSize, 0);
  normalMatrix_.reSize(0);
}

void
RankDeficientProjection::lastStateIs(const State & state) {
  rawBasis_.addState(state);
  orthoStatus_ = TO_DO;
}

void
RankDeficientProjection::lastStateSetIs(const StateSet & stateSet) {
  rawBasis_.resize(rawBasis_.numStates() + stateSet.numStates());
  for (int i = 0; i < stateSet.numStates(); ++i) {
    rawBasis_.addState(stateSet[i]);
  }
  assert(rawBasis_.numStates() == rawBasis_.maxNumStates());
  if (stateSet.numStates() > 0) {
    orthoStatus_ = TO_DO;
  }
}

void
RankDeficientProjection::lastStateSetIs(int stateCount, Scalar * data) {
  rawBasis_.resize(rawBasis_.numStates() + stateCount);
  for (int i = 0; i < stateCount; ++i) {
    State tempState(vectorSize(), data + 2 * vectorSize() * i);
    rawBasis_.addState(tempState);
  }
  if (stateCount > 0) {
    orthoStatus_ = TO_DO;
  }
}

void
RankDeficientProjection::stateDel() {
  orthoStatus_ = TO_DO;
  rawBasis_.reset(rawBasis_.vectorSize(), 0);
  metricBasis_.reset(metricBasis_.vectorSize(), 0);
  normalMatrix_.reSize(0);
}

void
RankDeficientProjection::metricIs(const GenSparseMatrix<Scalar> * K, const GenSparseMatrix<Scalar> * M,
                                  StateSet & independentBasis) {
  // Compute dual states
  int previousStateCount = metricBasis_.numStates();
  int newStateCount = rawBasis_.numStates();
  metricBasis_.resize(newStateCount);
  for (int i = previousStateCount; i < newStateCount; ++i) {
    metricBasis_.addState(State(rawBasis_[i], K, M));
  }

  assert(previousStateCount == normalMatrix_.dim());
  assert(metricBasis_.numStates() == rawBasis_.numStates());

  // Assemble normal matrix
  normalMatrix_.reSize(newStateCount);
  for (int r = 0; r < previousStateCount; ++r) {
    for (int c = previousStateCount; c < newStateCount; ++c) {
      normalMatrix_[r][c] = rawBasis_[r] * metricBasis_[c]; 
    }
  }
  for (int r = previousStateCount; r < newStateCount; ++r) {
    for (int c = 0; c < newStateCount; ++c) {
      normalMatrix_[r][c] = rawBasis_[r] * metricBasis_[c]; 
    }
  }

  // Factor normal matrix
  solver_->transposedMatrixIs(normalMatrix_);
  
  int numericalRank = solver_->factorRank();
  independentBasis.reset(vectorSize(), numericalRank);
  for (int i = 0; i < numericalRank; ++i) {
    int index = solver_->factorPermutation(i);
    independentBasis.addState(rawBasis_[index]); // Simple copy
  }

  orthoStatus_ = DONE;
}

const RankDeficientProjection::State &
RankDeficientProjection::projection(const StateSet & finalSet, const State & inputState, State & outputState) const {
  assert(orthoStatus() == DONE);
  assert(solver_->factorRank() == finalSet.numStates());

  Vector components(metricBasis_.numStates());
  for (int i = 0; i < finalSet.numStates(); ++i) {
    components[solver_->factorPermutation(i)] = metricBasis_[solver_->factorPermutation(i)] * inputState;
  }

  solver_->solution(components);

  outputState = 0.0; 
  for (int i = 0; i < finalSet.numStates(); ++i) {
    outputState.linAdd(components[solver_->factorPermutation(i)], finalSet[i]);
  }

  return outputState;
}

} /* end namespace Old */ } /* end namespace Pita */
