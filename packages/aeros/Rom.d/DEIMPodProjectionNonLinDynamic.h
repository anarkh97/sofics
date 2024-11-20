#ifndef ROM_DEIMPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DEIMPODPROJECTIONNONLINDYNAMIC_H

#include "LumpedPodProjectionNonLinDynamic.h"

#include <map>

namespace Rom {

class DEIMPodProjectionNonLinDynamic : public LumpedPodProjectionNonLinDynamic {
public:
  explicit DEIMPodProjectionNonLinDynamic(Domain *);

  virtual void preProcess();

  virtual void reBuild(ModalGeomState &, int, double, double);
  virtual double getStiffAndForce(ModalGeomState &, Vector &, Vector &, double = -1, ModalGeomState * = NULL, bool = false);


private:

  virtual void getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                          Corotator **allCorot, FullSquareMatrix *kelArray,
                                          Vector &residual, double lambda, double time, GeomState *refState,
                                          FullSquareMatrix *melArray, bool forceOnly);
  void computeRedKDyn(double);
  void readInterpolationBasis();
  void buildSparseInterpolationBasis();
  void buildReducedLinearOperator();

  GenFullSquareMatrix<double> *kelArrayCopy;

  VecBasis deimBasis_;
  VecBasis LinearStiffness_; //V^TKV, use expand to multiply K_red by modal coordinates

};

} /* end namespace Rom */

#endif /* ROM_DEIMPODPROJECTIONNONLINDYNAMIC_H */
