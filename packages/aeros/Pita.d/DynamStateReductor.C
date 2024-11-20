#include "DynamStateReductor.h"

namespace Pita {

DynamStateReductor::DynamStateReductor(const DynamStateBasis * reductionBasis,
                                       const RankDeficientSolver * solver) :
  reductionBasis_(reductionBasis),
  solver_(solver)
{}

Vector
DynamStateReductor::result(const DynamState & is) const {
  Vector answer(reducedBasisSize());
  return result(is, answer);
}

const Vector &
 DynamStateReductor::result(const DynamState & is, Vector & answer) const {
  assert(int(reducedBasisSize()) == solver()->matrixSize());

  if (is.vectorSize() != vectorSize()) {
    throw Fwk::RangeException("Dimension mismatch");
  }

  if (answer.size() != reducedBasisSize()) {
    answer.initialize(reducedBasisSize());
  }

  if (reducedBasisSize() != 0) {
    // Assemble relevant part of rhs
    for (int i = 0; i < solver_->factorRank(); ++i) {
      int index = solver()->factorPermutation(i);
      answer[index] = is * reductionBasis()->state(index);
    }

    /*log() << "Rhs = ";
    for (int i = 0; i < reducedBasisSize(); ++i) {
      log() << answer[i] << " ";
    }
    log() << "\n";*/
    
    // Perform in place resolution
    solver()->solution(answer);

    /*log() << "Solution = ";
    for (int i = 0; i < reducedBasisSize(); ++i) {
      log() << answer[i] << " ";
    }
    log() << "\n";*/
  }

  return answer;
}

} /* end namespace Pita */
