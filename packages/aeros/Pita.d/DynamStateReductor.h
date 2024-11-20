#ifndef PITA_DYNAMSTATEREDUCTOR_H
#define PITA_DYNAMSTATEREDUCTOR_H

#include "Fwk.h"

#include "DynamState.h"
#include <Math.d/Vector.h>

#include "DynamStateBasis.h"
#include "RankDeficientSolver.h"

namespace Pita {

class DynamStateReductor : public Fwk::PtrInterface<DynamStateReductor> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReductor);

  // Main function
  Vector result(const DynamState & is) const;
  // Optimized main function
  const Vector & result(const DynamState & is, Vector & answer) const;

  // Projection elements
  // Should have: reductionBasis->vectorSize() == solver->matrixSize()
  const DynamStateBasis * reductionBasis() const { return reductionBasis_.ptr(); }
  const RankDeficientSolver * solver() const { return solver_.ptr(); }

  // Sizes
  size_t vectorSize() const { return reductionBasis()->vectorSize(); }
  size_t reducedBasisSize() const { return reductionBasis()->stateCount(); }
 
  DynamStateReductor(const DynamStateBasis * reductionBasis,
                     const RankDeficientSolver * solver);

private:
  DynamStateBasis::PtrConst reductionBasis_;
  RankDeficientSolver::PtrConst solver_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReductor);
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATEREDUCTOR_H */
