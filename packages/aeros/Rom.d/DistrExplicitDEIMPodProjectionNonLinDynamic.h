#ifndef ROM_DISTREXPLICITDEIMPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITDEIMPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include <vector>
#include <map>
#include <utility>

namespace Rom {

class DistrExplicitDEIMPodProjectionNonLinDynamic : public DistrExplicitLumpedPodProjectionNonLinDynamic {
public:
  explicit DistrExplicitDEIMPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex); // Alternate internal force computation

private:
  void buildInterpolationBasis();
  void buildReducedLinearOperator();
  void subBuildInterpolationBasis(int iSub, std::vector< std::vector<std::pair<int,int> > > &maskedIndicesBuf);
  void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);
  void subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex);
  void subGetFollowerForceOnly(int iSub, DistrVector &f, double &t, int &tIndex);

  DistrInfo unassembledInfo;
  DistrVecBasis deimBasis_;
  DistrVecBasis ReducedStiffness;
  GenFullSquareMatrix<double> **kelArrayCopy;
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITDEIMPODPROJECTIONNONLINDYNAMIC_H */
