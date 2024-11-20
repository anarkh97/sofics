#ifndef PITA_DYNAMOPS_H
#define PITA_DYNAMOPS_H

#include "Fwk.h"

template <typename> class GenSparseMatrix;

namespace Pita {

class DynamOps : public Fwk::PtrInterface<DynamOps> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamOps);

  virtual const GenSparseMatrix<double> * massMatrix() const = 0;
  virtual const GenSparseMatrix<double> * stiffnessMatrix() const = 0;

protected:
  DynamOps() {}

private:
  DISALLOW_COPY_AND_ASSIGN(DynamOps);
};
  
} // end namespace Pita

#endif /* PITA_DYNAMOPS_H */
