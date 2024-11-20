#ifndef ROM_DISTRNONNEGATIVEMATRIXFACTORIZATION_H
#define ROM_DISTRNONNEGATIVEMATRIXFACTORIZATION_H

#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
class Communicator;
class SCDoubleMatrix;

#include <Eigen/Core>
#include <vector> 
#include <Math.d/SCMatrix.d/SCIntMatrix.h>

namespace Rom {

class DistrNonnegativeMatrixFactorization {
public:
  // Local data distribution
  int localRows() const { return localRows_; }
  int basisSize() const { return basisDimension_; }

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by colCount]
  double *matrixColBuffer(int col);
  // Local basis buffer: [localRows by basisDimension]
  const double *basisColBuffer(int col) const;

  void setEnergy(double energy) {energy_ = energy;}
  int energySVD(double energy, std::vector<int> rows, std::vector<int> cols); // svd dual snapshots to find singular value cutoff

  void solve();

  DistrNonnegativeMatrixFactorization(Communicator * comm, int rowCount, int colCount, int localRows, int basisDimension,
                                      int blockSize, int maxIter, double tol, int method, int nsub, int pqnNumInnerIter, double pqnAlpha);

  void solveNNLS_MRHS(SCDoubleMatrix &A, SCDoubleMatrix &B, SCDoubleMatrix &X, int flag, int SSCflag = 0);

private:
  // Disallow copy & assignment
  DistrNonnegativeMatrixFactorization(const DistrNonnegativeMatrixFactorization &);
  DistrNonnegativeMatrixFactorization & operator=(const DistrNonnegativeMatrixFactorization &);

  Communicator * communicator_;
  int rowCount_, colCount_, localRows_, basisDimension_, blockSize_, maxIter_, method_, nsub_, pqnNumInnerIter_;
  double tol_, pqnAlpha_, energy_;

  Eigen::MatrixXd matrixBuffer_;
  Eigen::MatrixXd basisBuffer_;
};

/* Helper functions for buffer access */
inline
double *
DistrNonnegativeMatrixFactorization::matrixColBuffer(int col) {
  return matrixBuffer_.data() + col*localRows_;
}

inline
const double *
DistrNonnegativeMatrixFactorization::basisColBuffer(int col) const {
  return basisBuffer_.data() + col*localRows_;
}

} // end namespace Rom
#endif

#endif /* ROM_DISTRNONNEGATIVEMATRIXFACTORIZATION_H */
