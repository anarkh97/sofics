#ifndef PITA_NLDYNAMOPS_H
#define PITA_NLDYNAMOPS_H

#include "DynamOps.h"
#include <Math.d/Vector.h>
#include "Types.h"

class GeomState;

namespace Pita {

class PitaNonLinDynamic;
  
class NlDynamOps : public DynamOps {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDynamOps);

  virtual const GenSparseMatrix<double> * massMatrix() const;
  virtual const GenSparseMatrix<double> * stiffnessMatrix() const;

  void displacementIs(const GenVector<double> & disp, Seconds time = Seconds(0.0));
  
  static Ptr New(PitaNonLinDynamic * probDesc) {
    return new NlDynamOps(probDesc);
  }
  
protected:
  explicit NlDynamOps(PitaNonLinDynamic * probDesc);
  
  ~NlDynamOps();
  
private:
  PitaNonLinDynamic * probDesc_;
  GeomState * refState_;
  GeomState * geomState_;
  GenVector<double> internalForceDummy_;
  GenVector<double> residualDummy_;
};
  
} // end namespace Pita

#endif /* PITA_NLDYNAMOPS_H */
