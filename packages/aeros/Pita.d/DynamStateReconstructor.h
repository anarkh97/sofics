#ifndef PITA_DYNAMSTATERECONSTRUCTOR_H
#define PITA_DYNAMSTATERECONSTRUCTOR_H

#include "Fwk.h"

#include "DynamState.h"
#include <Math.d/Vector.h>

#include "DynamStateBasis.h"

namespace Pita {

class DynamStateReconstructor : public Fwk::PtrInterface<DynamStateReconstructor> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReconstructor);

  // Main function
  DynamState result(const Vector & components) const;

  // Operator properties
  size_t vectorSize() const { return reconstructionBasis_->vectorSize(); }
  size_t reducedBasisSize() const { return reconstructionBasis_->stateCount(); }
 
  // Basis 
  const DynamStateBasis * reconstructionBasis() const { return reconstructionBasis_.ptr(); }
  void reconstructionBasisIs(const DynamStateBasis * rb) { reconstructionBasis_ = rb; }
  
  explicit DynamStateReconstructor(const DynamStateBasis * rb);
 
private:
  DynamStateBasis::PtrConst reconstructionBasis_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReconstructor);
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATERECONSTRUCTOR_H */
