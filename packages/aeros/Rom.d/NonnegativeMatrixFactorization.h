#ifndef ROM_NONNEGATIVEMATRIXFACTORIZATION_H
#define ROM_NONNEGATIVEMATRIXFACTORIZATION_H

#include "SimpleBuffer.h"

#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

namespace Rom {

class NonnegativeMatrixFactorization {
public:
  int rowCount()           const { return rowCount_; }
  int colCount()           const { return colCount_; }

  void robSizeIs(int rows, int cols);
  void matrixSizeIs(int rows, int cols);
  void maxIterIs(int maxIter);
  void toleranceIs(double tol);
  void nmfcAlphaIs(double alp);
  void nmfcBetaIs(double bet);
  void nmfcGammaIs(double gam);
  void basisDimensionIs(int basisDimension);
  void maxBasisDimensionIs(int maxBasisDimension);
  void numRandInitIs(int numRandInit);
  const double *robBuffer() const;
  const double *matrixBuffer() const;
  const double *matrixCol(int c) const;
  const double *robCol(int c) const;
  const double &matrixEntry(int r, int c) const;
  
  double *matrixBuffer();
  double *robBuffer();
  double *matrixCol(int c);
  double *robCol(int c);
  double &matrixEntry(int r, int c);

  void solve(int basisDimensionSmall);

  NonnegativeMatrixFactorization(int basisDimension, int method);
  ~NonnegativeMatrixFactorization() { fclose(bFile); }

private:
  int rowCount_, colCount_, maxBasisDimension_, basisDimension_, numRandInit_, method_, maxIter_;
  double tol_, nmfcAlpha, nmfcBeta, nmfcGamma;

  FILE *bFile; 

  SimpleBuffer<double> matrixBuffer_;
  SimpleBuffer<double> robBuffer_;

  // Disallow copy and assignment
  NonnegativeMatrixFactorization(const NonnegativeMatrixFactorization &);
  NonnegativeMatrixFactorization &operator=(const NonnegativeMatrixFactorization &);

#ifdef USE_EIGEN3
  Eigen::VectorXd solveNNLS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::VectorXd> &b);
  Eigen::MatrixXd solveNNLS_MRHS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::MatrixXd> &B);
  int findColumnWithLargestMagnitude(const Eigen::Ref<const Eigen::MatrixXd> &_X);
#endif
};

inline
const double *
NonnegativeMatrixFactorization::matrixBuffer() const {
  return matrixBuffer_.array();
}

inline
const double *
NonnegativeMatrixFactorization::matrixCol(int c) const {
  return matrixBuffer() + (c * rowCount_);
}

inline
const double &
NonnegativeMatrixFactorization::matrixEntry(int r, int c) const {
  return *(matrixCol(c) + r);
}

inline
double *
NonnegativeMatrixFactorization::matrixBuffer() {
  return const_cast<double *>(const_cast<const NonnegativeMatrixFactorization *>(this)->matrixBuffer());
}

inline
double *
NonnegativeMatrixFactorization::matrixCol(int c) {
  return const_cast<double *>(const_cast<const NonnegativeMatrixFactorization *>(this)->matrixCol(c));
}

inline
double &
NonnegativeMatrixFactorization::matrixEntry(int r, int c) {
  return const_cast<double &>(const_cast<const NonnegativeMatrixFactorization *>(this)->matrixEntry(r, c));
}

inline
const double *
NonnegativeMatrixFactorization::robBuffer() const {
  return robBuffer_.array();
}

inline
const double *
NonnegativeMatrixFactorization::robCol(int c) const {
  return robBuffer() + (c * rowCount_);
}

inline
double *
NonnegativeMatrixFactorization::robBuffer() {
  return const_cast<double *>(const_cast<const NonnegativeMatrixFactorization *>(this)->robBuffer());
}

inline
double *
NonnegativeMatrixFactorization::robCol(int c) {
  return const_cast<double *>(const_cast<const NonnegativeMatrixFactorization *>(this)->robCol(c));
}



} /* end namespace Rom */

#endif /* ROM_NONNEGATIVEMATRIXFACTORIZATION_H */
