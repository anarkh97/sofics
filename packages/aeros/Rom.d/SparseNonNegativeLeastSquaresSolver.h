#ifndef ROM_SPARSENONNEGATIVELEASTSQUARESSOLVER_H
#define ROM_SPARSENONNEGATIVELEASTSQUARESSOLVER_H

#include "SimpleBuffer.h"
#include <vector>

namespace Rom {

template<typename MatrixBufferType, typename SizeType>
class SparseNonNegativeLeastSquaresSolver {
public:
  typedef typename MatrixBufferType::value_type Scalar;

  // Problem size
  long equationCount() const { return equationCount_; }
  long unknownCount() const { return unknownCount_; }

  void problemSizeIs(long eqnCount, long unkCount);

  // Options
  double relativeTolerance() const { return relativeTolerance_; }
  void relativeToleranceIs(double relTol) { relativeTolerance_ = relTol; }

  bool verboseFlag() const { return verboseFlag_; }
  void verboseFlagIs(bool verFlg) { verboseFlag_ = verFlg; }

  bool scalingFlag() const { return scalingFlag_; }
  void scalingFlagIs(bool scaFlg) { scalingFlag_ = scaFlg; }

  bool centerFlag() const { return centerFlag_; }
  void centerFlagIs(bool cenFlg) { centerFlag_ = cenFlg; }

  bool reverseFlag() const { return reverseFlag_; }
  void reverseFlagIs(bool revFlg) { reverseFlag_ = revFlg; }

  bool projectFlag() const { return projectFlag_; }
  void projectFlagIs(bool posFlg) { projectFlag_ = posFlg; }

  bool positivityFlag() const { return positivity_; }
  void positivityIs(bool posFlg) { positivity_ = posFlg; }

  int solverType() const { return solverType_; }
  void solverTypeIs(int solTyp) { solverType_ = solTyp; }

  double maxSizeRatio() const { return maxSizeRatio_; }
  void maxSizeRatioIs(double maxSze) { maxSizeRatio_ = maxSze; }

  double maxIterRatio() const { return maxIterRatio_; }
  void maxIterRatioIs(double maxIte) { maxIterRatio_ = maxIte; }

  int maxNumElems() const { return maxNumElems_; }
  void maxNumElemsIs(int maxElem) { maxNumElems_ = maxElem; }

  void useHotStart(bool hs) { hotStart_ = hs; }

  // Buffers: Internal column-major ordering, zero-based indexing
  // Matrix buffer: [equationCount by unknownCount]
  Scalar matrixEntry(int row, int col) const;
  typename MatrixBufferType::const_iterator matrixColBuffer(int col) const;
  typename MatrixBufferType::const_iterator matrixBuffer() const;

  Scalar & matrixEntry(int row, int col);
  typename MatrixBufferType::iterator matrixColBuffer(int col);
  typename MatrixBufferType::iterator matrixBuffer();

  // Rhs buffer: [equationCount]
  Scalar rhsEntry(int row) const;
  const Scalar * rhsBuffer() const;
  
  Scalar & rhsEntry(int row);
  Scalar * rhsBuffer();

  // Primal solution buffer: [unknownCount]
  Scalar solutionEntry(int row) const;
  const Scalar * solutionBuffer() const;

  // Dual solution buffer: [unknownCount]
  Scalar dualSolutionEntry(int row) const;
  const Scalar * dualSolutionBuffer() const;

  // Error magnitude
  double errorMagnitude() const { return errorMagnitude_; }

  // Perform solve
  void solve();

  // Constructor
  SparseNonNegativeLeastSquaresSolver();

private:

  long int equationCount_;
  long int unknownCount_;

  double relativeTolerance_;
  MatrixBufferType matrixBuffer_;

  SimpleBuffer<Scalar> rhsBuffer_;
  SimpleBuffer<Scalar> solutionBuffer_;
  SimpleBuffer<Scalar> dualSolutionBuffer_;

  std::vector<long int> indices; 

  double errorMagnitude_;
  bool verboseFlag_;
  bool scalingFlag_;
  bool centerFlag_;
  bool reverseFlag_;
  bool projectFlag_;
  bool positivity_;
  bool hotStart_;
  int solverType_;
  double maxSizeRatio_;
  double maxIterRatio_;
  int maxNumElems_;

  // Disallow copy & assignment
  SparseNonNegativeLeastSquaresSolver(const SparseNonNegativeLeastSquaresSolver &);
  SparseNonNegativeLeastSquaresSolver &operator=(const SparseNonNegativeLeastSquaresSolver &);
};

// Inline member functions
template<typename MatrixBufferType, typename SizeType>
typename MatrixBufferType::const_iterator
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::matrixBuffer() const {
  return matrixBuffer_.begin();
} 

template<typename MatrixBufferType, typename SizeType>
typename MatrixBufferType::const_iterator
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::matrixColBuffer(int col) const {
  return matrixBuffer() + (col * equationCount_);
}

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::matrixEntry(int row, int col) const {
  return matrixColBuffer(col)[row];
}

template<typename MatrixBufferType, typename SizeType>
typename MatrixBufferType::iterator
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::matrixBuffer() {
  return matrixBuffer_.begin();
} 

template<typename MatrixBufferType, typename SizeType>
typename MatrixBufferType::iterator
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::matrixColBuffer(int col) {
  return matrixBuffer() + (col * equationCount_);
}

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar &
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::matrixEntry(int row, int col) {
  return matrixColBuffer(col)[row];
}

template<typename MatrixBufferType, typename SizeType>
const typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar *
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::rhsBuffer() const {
  return rhsBuffer_.array();
} 

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::rhsEntry(int row) const {
  return rhsBuffer()[row];
}

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar *
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::rhsBuffer() {
  return rhsBuffer_.array();
} 

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar &
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::rhsEntry(int row) {
  return rhsBuffer()[row];
}

template<typename MatrixBufferType, typename SizeType>
const typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar *
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::solutionBuffer() const {
  return solutionBuffer_.array();
} 

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::solutionEntry(int row) const {
  return solutionBuffer()[row];
}

template<typename MatrixBufferType, typename SizeType>
const typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar *
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::dualSolutionBuffer() const {
  return dualSolutionBuffer_.array();
} 

template<typename MatrixBufferType, typename SizeType>
typename SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::Scalar
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::dualSolutionEntry(int row) const {
  return dualSolutionBuffer()[row];
}

} // end namespace Rom

#endif /* ROM_SPARSENONNEGATIVELEASTSQUARESSOLVER_H */
