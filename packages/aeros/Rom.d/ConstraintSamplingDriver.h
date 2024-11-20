#ifndef ROM_CONSTRAINTSAMPLINGDRIVER_H
#define ROM_CONSTRAINTSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"
#include "VecBasisOps.h"
#include "FileNameInfo.h"
#include "VecBasisFile.h"
#include "VecNodeDof6Map.h"

#include "VecNodeDof6Conversion.h"

#include "SparseNonNegativeLeastSquaresSolver.h"

#include <set>

#ifdef USE_EIGEN3
#include <Eigen/LU>
#endif

namespace Rom {

class ConstraintSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit ConstraintSamplingDriver(Domain *);

  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>& solver() { return solver_; }

private:

  void readInBasis(VecBasis &podBasis_, BasisId::Type type, BasisId::Level level, bool vectorQuant = true, int i = 0, int podSizeMax = 0, bool normalized = false);
  template <typename Scalar> void writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type,BasisId::Level, bool vectorQuant = true);

  int  snapSize(BasisId::Type type, std::vector<int> &snapshotCounts);
  void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis_,
                          const VecNodeDof6Conversion *vecDofConversion,
                          std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config);
  void buildConstraintArray(const VecBasis &dualBasis, const VecBasis &displac);

  void setSolverFlags(); 
  int NNZRows(const VecBasis &basis);
  void writeSampledMesh(std::vector<double> &solution);
  void expandSolution(std::vector<double> &solution, std::vector<double> &expandedSolution);
  void writeWeightedDualBasis(std::vector<double> &solution, std::vector<double> &expandedSolution);

  VecNodeDof6Conversion *converter6;
  VecNodeDof1Conversion *converter1;

  VecBasis dualBasis;

  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> solver_;

};

} /* end namespace Rom */

Rom::DriverInterface *constraintSamplingDriverNew(Domain *);
 
#endif /* ROM_BASIS_ORTHODRIVER_H */
