#ifndef ROM_ELEMENTSAMPLINGDRIVER_H
#define ROM_ELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"

#include "SparseNonNegativeLeastSquaresSolver.h"
#include <Math.d/Vector.h>
#include "VecBasis.h"
#include <vector>
#include "MeshDesc.h"
#include "FileNameInfo.h"

class Domain;
class Corotator;
class GeomState; 
class StaticTimers;
class FileNameInfo;
class MeshDesc;

template <typename Scalar> class GenFullSquareMatrix; 
template <typename Scalar> class GenVector;
typedef GenVector<double> Vector;

namespace Rom {
void outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const int podVectorCount);
void outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const std::vector<int> &localBasisSize);
template<typename WeightsVecType, typename ElemIdsVecType>
void outputFullWeights(const WeightsVecType &weights, const ElemIdsVecType &elemIds, int j=-1);
std::string getMeshFilename(const FileNameInfo &fileInfo);

template<typename MatrixBufferType = std::vector<double>, typename SizeType = size_t>
class ElementSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve(); // overriden
  void computeSolution(Vector &solution, double relativeTolerance, bool verboseFlag = true);
  void postProcess(Vector &solution, bool verboseFlag = true);
  VecBasis& podBasis() { return podBasis_; }
  VecBasis& displac() { return displac_; }
  VecBasis* veloc() { if(!veloc_) veloc_ = new VecBasis; return veloc_; }
  VecBasis* accel() { if(!accel_) accel_ = new VecBasis; return accel_; }
  int vectorSize() const;
  void timeStampsIs(const std::vector<double> &tst) { timeStamps_ = tst; }
  void snapshotCountsIs(const std::vector<int> &scs) { snapshotCounts_ = scs; }

  explicit ElementSamplingDriver(Domain *);
  ~ElementSamplingDriver();

  virtual void preProcess();
  template<typename VecBasisType>
  void assembleTrainingData(VecBasisType &podBasis, int podVectorCount, VecBasisType &displac,
                            VecBasisType *veloc, VecBasisType *accel, VecBasisType *ndscfgCoords,
                            VecBasisType *ndscfgMassOrthogonalBases, int j=-1);

  void clean();
  SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>& solver() { return solver_; }
  int elementCount() const;

protected:
  void preProcessGlobal(AllOps<double>& allOps);
  void preProcessLocal(AllOps<double>& allOps, int j);
  void buildDomainCdsa();
  void postProcessLocal(Vector &solution, std::vector<int> packedToInput, int j,
                        std::vector<int> &sampleElemIds, std::map<int, double> &weights,
                        bool verboseFlag=true);
  void postProcessGlobal(std::vector<int> &sampleElemIds, std::vector<std::map<int, double> > &weights, bool verboseFlag=true);
  void makePackedToInput(std::vector<int> &packedToInput);
  void addContactElems(std::vector<int> &sampleElemIds, std::vector<std::map<int, double> > &weights);
  Domain *domain_;

  Corotator **corotators_;
  GeomState *geomState_;
  GenFullSquareMatrix<double> *kelArray_, *melArray_;

  VecBasis podBasis_;
  VecBasis displac_;
  VecBasis *veloc_;
  VecBasis *accel_;
  VecBasis *ndscfgCoords_;
  VecBasis *ndscfgMassOrthogonalBases_;
  std::vector<double> timeStamps_;
  std::vector<int> snapshotCounts_;
  SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType> solver_;
};

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *);

#endif /* ROM_ELEMENTSAMPLINGDRIVER_H */
