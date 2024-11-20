#ifndef ROM_DEIM_CONSTRAINTSAMPLINGDRIVER_H
#define ROM_DEIM_CONSTRAINTSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"
#include "VecBasisOps.h"
#include "FileNameInfo.h"
#include "VecBasisFile.h"
#include "VecNodeDof6Map.h"

#include "VecNodeDof6Conversion.h"

#include <set>

#ifdef USE_EIGEN3
#include <Eigen/LU>
#endif

namespace Rom {

class DEIMConstraintSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit DEIMConstraintSamplingDriver(Domain *);

private:

  void readInBasis(VecBasis &podBasis_, BasisId::Type type, BasisId::Level level, bool vectorQuant = true, int i = 0, int podSizeMax = 0, bool normalized = false);
  template <typename Scalar> void writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type,BasisId::Level, bool vectorQuant = true);

  int  snapSize(BasisId::Type type, std::vector<int> &snapshotCounts);
  void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis_,
                          const VecNodeDof6Conversion *vecDofConversion,
                          std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config);
  void buildConstraintArray(const VecBasis &displac);
  void OrthoConstraintSnap(VecBasis &constraintBasis,std::vector<double> &SVs);  

  void writeProjConstraintSnap(); 

  void filterRows(VecBasis &constraintBasis);
  void computeInterpIndices(VecBasis &constraintBasis, std::vector<int> &maskIndices); 
  void computeAndWriteDEIMBasis(VecBasis &constraintBasis, std::vector<int> &maskIndices);  
  void writeSampledMesh(std::vector<int> &maskIndices);

  void computeErrorBound();

  VecNodeDof6Conversion *converter6;
  VecNodeDof1Conversion *converter1;

//  VecBasis podBasis;
  VecBasis dualBasis;
  VecBasis deimBasis;

  double MPOSingularValue; //m+1 singular value of basis V
};

} /* end namespace Rom */

Rom::DriverInterface *deimConstraintSamplingDriverNew(Domain *);
 
#endif /* ROM_BASIS_ORTHODRIVER_H */
