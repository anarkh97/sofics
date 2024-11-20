#ifndef ROM_DISTRSNAPSHOTCLUSTERINGSOLVER_H
#define ROM_DISTRSNAPSHOTCLUSTERINGSOLVER_H

#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
class Communicator;
class SCDoubleMatrix;

#include <Eigen/Core>
#include <vector>

namespace Rom {

class DistrSnapshotClusteringSolver {
public:
  // Local data distribution
  int localRows() const { return localRows_; }

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by colCount]
  double *matrixColBuffer(int col);
  const double *clusterColBuffer(int i, int col) const;
  double getClusterColTimeStamp(int i, int col) const; 
  // Local cluster centroid buffer: [localRows by numClusters]
  const double *clusterCentroidBuffer(int i) const;

  int clusterCol(int i, int col) const;
  int clusterColCount(int i) const;

  int getNumClusters() {return numClusters_;}
  int solverType() const { return solverType_; }
  void solverTypeIs(int solTyp) { solverType_ = solTyp; }

  int kmMaxIter() const { return kmMaxIter_; }
  void kmMaxIterIs(int kmMaxIter) { kmMaxIter_ = kmMaxIter; }

  int kmSeed() const { return kmSeed_; }
  void kmSeedIs(int kmSeed) { kmSeed_ = kmSeed; }
  void setNNLSTolerance(double _tol) { nnlsTol = _tol; }
  void setColTimeStamp(int col, double timeStamp) { timeStamps_[col] = timeStamp; }
  
  void recomputeCentroids();

  void solve();

  DistrSnapshotClusteringSolver(Communicator * comm, int rowCount, int colCount, int localRows, int numClusters,
                                int blockSize);

private:
  // Disallow copy & assignment
  DistrSnapshotClusteringSolver(const DistrSnapshotClusteringSolver &);
  DistrSnapshotClusteringSolver & operator=(const DistrSnapshotClusteringSolver &);

  Communicator * communicator_;
  int rowCount_, colCount_, localRows_, numClusters_, blockSize_;

  int solverType_;
  int kmMaxIter_;
  int kmSeed_;
  double nnlsTol;

  Eigen::MatrixXd matrixBuffer_;
  Eigen::MatrixXd centroidBuffer_;
  std::vector<std::vector<int> > clusterCols_;
  std::vector<double> timeStamps_; 
};

/* Helper functions for buffer access */
inline
double *
DistrSnapshotClusteringSolver::matrixColBuffer(int col) {
  return matrixBuffer_.data() + col*localRows_;
}

inline
const double *
DistrSnapshotClusteringSolver::clusterColBuffer(int i, int col) const {
  return matrixBuffer_.data() + clusterCols_[i][col]*localRows_;
}

inline
double
DistrSnapshotClusteringSolver::getClusterColTimeStamp(int i, int col) const {
  return timeStamps_[clusterCols_[i][col]];
}

inline
int
DistrSnapshotClusteringSolver::clusterCol(int i, int col) const {
  return clusterCols_[i][col];
}

inline
int
DistrSnapshotClusteringSolver::clusterColCount(int i) const {
  return clusterCols_[i].size();
}

inline
const double *
DistrSnapshotClusteringSolver::clusterCentroidBuffer(int i) const {
  return centroidBuffer_.data() + i*localRows_;
}

} // end namespace Rom
#endif

#endif /* ROM_DISTRSNAPSHOTCLUSTERINGSOLVER_H */
