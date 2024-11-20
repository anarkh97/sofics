#ifndef ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H
#define ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H

#include "PodProjectionNonLinDynamic.h"

#include <map>

namespace Rom {

class LumpedPodProjectionNonLinDynamic : public PodProjectionNonLinDynamic {
public:
  explicit LumpedPodProjectionNonLinDynamic(Domain *);

  virtual void preProcess();
  virtual void updateStates(ModalGeomState *refState, ModalGeomState& geomState, double time);
  virtual void setLocalReducedMesh(int j);

private:
  virtual void getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                          Corotator **allCorot, FullSquareMatrix *kelArray,
                                          Vector &residual, double lambda, double time, GeomState *refState,
                                          FullSquareMatrix *melArray, bool forceOnly);

protected:
  std::vector<std::map<int, double> > packedElementWeights_;
  std::vector<int> packedWeightedNodes_;
  std::set<int> packedWeightedElems_;
  void buildPackedElementWeights();
  int localReducedMeshId_;
  std::vector<std::vector<int> > localPackedWeightedNodes_;
};

} /* end namespace Rom */

#endif /* ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H */
