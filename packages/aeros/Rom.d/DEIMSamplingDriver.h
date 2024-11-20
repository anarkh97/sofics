#ifndef ROM_DEIM_SAMPLINGDRIVER_H
#define ROM_DEIM_SAMPLINGDRIVER_H

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

class DEIMSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit DEIMSamplingDriver(Domain *);

private:

  void readInBasis(VecBasis &podBasis, BasisId::Type type, BasisId::Level level, int podSizeMax, bool normalized = false);
  template <typename Scalar> void writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type,BasisId::Level);

  int  snapSize(BasisId::Type type, std::vector<int> &snapshotCounts);
  void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                          const VecNodeDof6Conversion *vecDofConversion,
                          std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config);
  void buildForceArray(VecBasis &forceBasis, const VecBasis &displac,
                       const VecBasis *veloc, const VecBasis *accel,std::vector<double> timeStamps_,
                       std::vector<int> snapshotCounts_);
  void OrthoForceSnap(VecBasis &forceBasis,std::vector<double> &SVs);  

  int  elementCount() const; 

  void writeProjForceSnap(VecBasis &forceBasis,std::vector<double> &timeStamps); 

#ifdef USE_EIGEN3
  void getFullNodeIndices(Eigen::Matrix<double,Eigen::Dynamic,1> res,int MaxCoeff,std::vector<int> &container, std::set<int> &auxilaryIndices);
#endif

  void computeInterpIndices(VecBasis &forceBasis, std::vector<int> &maskIndices); 
  void computeAndWriteDEIMBasis(VecBasis &forceBasis, std::vector<int> &maskIndices);  
  void writeSampledMesh(std::vector<int> &maskIndices);

  void computeErrorBound();

  VecNodeDof6Conversion *converter;
  VecNodeDof6Map *nodeDofMap;

  VecBasis podBasis_;
  VecBasis deimBasis;

  FullSquareMatrix *kelArrayCopy;

  double MPOSingularValue; //m+1 singular value of basis V
};

} /* end namespace Rom */

Rom::DriverInterface *deimSamplingDriverNew(Domain *);

#endif /* ROM_BASIS_ORTHODRIVER_H */
