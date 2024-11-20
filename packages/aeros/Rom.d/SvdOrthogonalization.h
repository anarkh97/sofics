#ifndef ROM_SVDORTHOGONALIZATION_H
#define ROM_SVDORTHOGONALIZATION_H

#include "SimpleBuffer.h"

namespace Rom {

class SvdOrthogonalization {
public:
  int rowCount()           const { return rowCount_; }
  int colCount()           const { return colCount_; }
  int singularValueCount() const { return smallestDimension_; }

  void matrixSizeIs(int rows, int cols);

  const double *matrixBuffer() const;
  const double *matrixCol(int c) const;
  const double &matrixEntry(int r, int c) const;
  
  double *matrixBuffer();
  double *matrixCol(int c);
  double &matrixEntry(int r, int c);

  const double *singularValueBuffer() const;
  const double &singularValue(int rank) const;

  void solve();

  SvdOrthogonalization();

private:
  int rowCount_, colCount_, smallestDimension_;

  SimpleBuffer<double> matrixBuffer_;
  SimpleBuffer<double> sigmaBuffer_;

  // Disallow copy and assignment
  SvdOrthogonalization(const SvdOrthogonalization &);
  SvdOrthogonalization &operator=(const SvdOrthogonalization &);
};

inline
const double *
SvdOrthogonalization::matrixBuffer() const {
  return matrixBuffer_.array();
}

inline
const double *
SvdOrthogonalization::matrixCol(int c) const {
  return matrixBuffer() + (c * rowCount_);
}

inline
const double &
SvdOrthogonalization::matrixEntry(int r, int c) const {
  return *(matrixCol(c) + r);
}

inline
double *
SvdOrthogonalization::matrixBuffer() {
  return const_cast<double *>(const_cast<const SvdOrthogonalization *>(this)->matrixBuffer());
}

inline
double *
SvdOrthogonalization::matrixCol(int c) {
  return const_cast<double *>(const_cast<const SvdOrthogonalization *>(this)->matrixCol(c));
}

inline
double &
SvdOrthogonalization::matrixEntry(int r, int c) {
  return const_cast<double &>(const_cast<const SvdOrthogonalization *>(this)->matrixEntry(r, c));
}

inline
const double *
SvdOrthogonalization::singularValueBuffer() const {
  return sigmaBuffer_.array();
}

inline
const double &
SvdOrthogonalization::singularValue(int rank) const {
  return sigmaBuffer_[rank];
}

} /* end namespace Rom */

#endif /* ROM_SVDORTHOGONALIZATION_H */
