#ifndef ROM_UDEIMPODPROJECTIONNONLINDYNAMIC_H
#define ROM_UDEIMPODPROJECTIONNONLINDYNAMIC_H

#include "LumpedPodProjectionNonLinDynamic.h"

#include <map>

namespace Rom {

class UDEIMPodProjectionNonLinDynamic : public LumpedPodProjectionNonLinDynamic {
public:
  explicit UDEIMPodProjectionNonLinDynamic(Domain *);

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

  VecBasis udeimBasis_;
  VecBasis LinearStiffness_; //V^TKV, use expand to multiply K_red by modal coordinates

  int unassembledInfo;
  int numOfIndices;

  std::map<int, std::vector<int> > unassembledElemDOFMask;

};

} /* end namespace Rom */

#endif /* ROM_UDEIMPODPROJECTIONNONLINDYNAMIC_H */
