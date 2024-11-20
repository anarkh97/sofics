#ifndef ROM_DISTREXPLICITUDEIMPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITUDEIMPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include <vector>
#include <map>
#include <utility>

namespace Rom {

class DistrExplicitUDEIMPodProjectionNonLinDynamic : public DistrExplicitLumpedPodProjectionNonLinDynamic {
public:
  explicit DistrExplicitUDEIMPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex); // Alternate internal force computation

private:
  void buildInterpolationBasis();
  void buildReducedLinearOperator();
  void subBuildUnassembledMask(int iSub);
  void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);
  void subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex);
  void subGetFollowerForceOnly(int iSub, DistrVector &f, double &t, int &tIndex);

  int numOfIndices;

  DistrInfo unassembledInfo;
  DistrVecBasis udeimBasis_;
  DistrVecBasis ReducedStiffness;
  GenFullSquareMatrix<double> **kelArrayCopy;
  std::vector<std::map<int, std::vector<int> > > unassembledElemDOFMask;
  std::vector<std::map<int,std::vector<int> > > columnKey; //element->DOF->columnMap
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITUDEIMPODPROJECTIONNONLINDYNAMIC_H */
