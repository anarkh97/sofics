#include "SvdOrthogonalization.h"

#include <Utils.d/linkfc.h>

#include <algorithm>
#include <cstddef>
#include <cassert>

extern "C" {
  void _FORTRAN(dgesvd)(const char *jobu, const char *jobvt,
                        const int *m, const int *n, double *a, const int *lda,
                        double *s, double *u, const int *ldu, double *vt, const int *ldvt,
                        double *work, const int *lwork, int *info);
}

namespace Rom {

SvdOrthogonalization::SvdOrthogonalization() :
  rowCount_(0),
  colCount_(0),
  smallestDimension_(0),
  matrixBuffer_(0),
  sigmaBuffer_(0)
{}

void
SvdOrthogonalization::matrixSizeIs(int rows, int cols) {
  const int minDim = std::min(rows, cols);

  matrixBuffer_.sizeIs(rows * cols);
  sigmaBuffer_.sizeIs(minDim);

  rowCount_ = rows;
  colCount_ = cols;
  smallestDimension_ = minDim;
}

void
SvdOrthogonalization::solve() {
  const char overwrite = 'O';
  const char ignore = 'N';

  const int dummy_ld = 1; 
  const int lda = (rowCount_ > 0) ? rowCount_ : dummy_ld;

  int info;

  const int query = -1;
  double queryAns;
  _FORTRAN(dgesvd)(&overwrite, &ignore,
                   &rowCount_, &colCount_, matrixBuffer_.array(), &lda,
                   sigmaBuffer_.array(), NULL, &dummy_ld, NULL, &dummy_ld,
                   &queryAns, &query, &info);
  assert(info == 0);

  const int lwork = static_cast<int>(queryAns);
  SimpleBuffer<double> work(lwork); 
  _FORTRAN(dgesvd)(&overwrite, &ignore,
                   &rowCount_, &colCount_, matrixBuffer_.array(), &lda,
                   sigmaBuffer_.array(), NULL, &dummy_ld, NULL, &dummy_ld,
                   work.array(), &lwork, &info);
  assert(info == 0);
}

} /* end namespace Rom */
