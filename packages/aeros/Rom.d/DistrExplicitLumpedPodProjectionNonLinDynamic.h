#ifndef ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include <vector>
#include <map>

namespace Rom {

class DistrExplicitLumpedPodProjectionNonLinDynamic : public DistrExplicitPodProjectionNonLinDynamicBase {
public:
  explicit DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  void updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n);
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex); // Alternate internal force computation
  MDDynamMat * buildOps(double, double, double);

  virtual void setLocalReducedMesh(int j);

private:
  void buildPackedElementWeights();
  void subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex);
  void subUpdateWeightedNodesOnly(int iSub, DistrVector &v);
  void subSetVelocityWeightedNodesOnly(int iSub, DistrVector &v);
  void subInitWeightedStiffOnly(int iSub);
  void subBuildPackedElementWeights(int iSub);

  std::vector<std::vector<int> > packedWeightedNodes_;
  std::vector<std::vector<std::vector<int> > > localPackedWeightedNodes_;
  
protected:
  void subTransformWeightedNodesOnly(int iSub, DistrVector &v, int type);

  std::vector<std::vector<std::map<int, double> > > packedElementWeights_;

  int localReducedMeshId_;

  SubDOp *K;
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H */
