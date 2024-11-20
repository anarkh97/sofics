#ifndef ROM_UDEIM_SAMPLINGDRIVER_H
#define ROM_UDEIM_SAMPLINGDRIVER_H

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

class UDEIMSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit UDEIMSamplingDriver(Domain *);

private:

  void readInBasis(VecBasis &podBasis, BasisId::Type type, BasisId::Level level, int podSizeMax, bool normalized = false);
  template <typename Scalar> void writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type,BasisId::Level);

  int  snapSize(BasisId::Type type, std::vector<int> &snapshotCounts);
  void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                          const VecNodeDof6Conversion *vecDofConversion,
                          std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config);
  void buildForceArray(VecBasis &unassembledForceBasis,const VecBasis &displac,
                       const VecBasis *veloc,const VecBasis *accel,std::vector<double> timeStamps_,
                       std::vector<int> snapshotCounts_);
  void OrthoForceSnap(VecBasis &forceBasis,std::vector<double> &SVs);  
  void computeAssembledIndices(std::vector<int> &umaskIndices, std::vector<int> &amaskIndices, std::set<int> &selectedElems, std::vector<std::pair<int,int> > &elemRankDOFContainer); 

  int  elementCount() const; 
  int  unassembledVecInfo();

  void writeUnassembledForceSnap(VecBasis &unassembledForceBasis); 
  void assembleBasisVectors(VecBasis &assembledForceBasis, VecBasis &unassembledForceBasis);
  void readUnassembledForceSnap(VecBasis &unassembledForceBasis, std::vector<double> &SVs);

  void computeInterpIndices(VecBasis &forceBasis, std::vector<int> &maskIndices); 
  void computeAndWriteUDEIMBasis(VecBasis &unassembledForceBuf,VecBasis &assembledForceBuf,std::vector<int> &umaskIndices,std::vector<int> &amaskIndices, std::vector<double> &singularVals);  
  void writeSampledMesh(std::vector<int> &maskIndices, std::set<int> &selectedElems, std::vector<std::pair<int,int> > &elemRankDOFContainer);
  void getFullElemIndices(int ,std::set<int> &,std::set<int> &);
  void getFullNodeIndices(int ,int ,std::set<int> &,std::set<int> &);
  void computeErrorBound(std::vector<int> &umaskIndices);

  VecNodeDof6Conversion *converter;
  VecNodeDof6Map *nodeDofMap;

#ifdef USE_EIGEN3
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> compressedDBTranspose;
#endif

  VecBasis podBasis_;
  std::map<int, std::pair<int,int> > uDOFaDOFmap;// key: unassembled index, MapValue: paired element and element DOF

  FullSquareMatrix *kelArrayCopy; 
};

} /* end namespace Rom */

Rom::DriverInterface *udeimSamplingDriverNew(Domain *);

#endif /* ROM_BASIS_ORTHODRIVER_H */
