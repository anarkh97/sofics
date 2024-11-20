#ifndef ROM_DISTRSNAPSHOTROWCLUSTERINGSOLVER_H
#define ROM_DISTRSNAPSHOTROWCLUSTERINGSOLVER_H

#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
class Communicator;
class SCDoubleMatrix;

#include <Eigen/Core>
#include <vector>

namespace Rom {

class DistrSnapshotRowClusteringSolver {
public:
  // Local data distribution
  int colCount() const { return colCount_; }
  int localRows() const { return localRows_; }

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by colCount]
  double *matrixColBuffer(int col);
  int clusterLocalRow(int i, int row) const;
  int clusterLocalRowCount(int i) const;

  int getNumClusters() {return numClusters_;}
  int solverType() const { return solverType_; }
  void solverTypeIs(int solTyp) { solverType_ = solTyp; }

  int kmMaxIter() const { return kmMaxIter_; }
  void kmMaxIterIs(int kmMaxIter) { kmMaxIter_ = kmMaxIter; }

  int kmSeed() const { return kmSeed_; }
  void kmSeedIs(int kmSeed) { kmSeed_ = kmSeed; }
  void setNNLSTolerance(double _tol) { nnlsTol = _tol; }
  void setColTimeStamp(int col, double timeStamp) { timeStamps_[col] = timeStamp; }
  double getColTimeStamp(int col) { return timeStamps_[col]; }

  void solve();

  DistrSnapshotRowClusteringSolver(Communicator * comm, int rowCount, int colCount, int localRows, int numClusters,
                                   int blockSize);

private:
  // Disallow copy & assignment
  DistrSnapshotRowClusteringSolver(const DistrSnapshotRowClusteringSolver &);
  DistrSnapshotRowClusteringSolver & operator=(const DistrSnapshotRowClusteringSolver &);

  Communicator * communicator_;
  int rowCount_, colCount_, localRows_, numClusters_, blockSize_;

  int solverType_;
  int kmMaxIter_;
  int kmSeed_;
  double nnlsTol;

  Eigen::MatrixXd matrixBuffer_;
  std::vector<std::vector<int> > clusterLocalRows_;
  std::vector<double> timeStamps_;
};

/* Helper functions for buffer access */
inline
double *
DistrSnapshotRowClusteringSolver::matrixColBuffer(int col) {
  return matrixBuffer_.data() + col*localRows_;
}

inline
int
DistrSnapshotRowClusteringSolver::clusterLocalRow(int i, int row) const {
  return clusterLocalRows_[i][row];
}

inline
int
DistrSnapshotRowClusteringSolver::clusterLocalRowCount(int i) const {
  return clusterLocalRows_[i].size();
}

} // end namespace Rom
#endif

#endif /* ROM_DISTRSNAPSHOTROWCLUSTERINGSOLVER_H */
