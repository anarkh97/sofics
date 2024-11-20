#include "SparseNonNegativeLeastSquaresSolver.C"
#include <vector>

namespace Rom {

template
SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>
::SparseNonNegativeLeastSquaresSolver();

template
void
SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>
::problemSizeIs(long eqnCount, long unkCount);

template
void
SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>
::solve();

} // end namespace Rom
