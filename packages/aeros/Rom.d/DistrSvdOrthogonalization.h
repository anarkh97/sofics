#ifndef ROM_DISTRSVDORTHOGONALIZATION_H
#define ROM_DISTRSVDORTHOGONALIZATION_H

class Communicator;

#include "SimpleBuffer.h"

namespace Rom {

// ScaLAPACK-based distributed SVD orthogonalization
class DistrSvdOrthogonalization {
public:
  typedef double Scalar;
 
  // Underlying communicator
  const Communicator *communicator() const { return communicator_; }

  // Cpu topology
  int rowCpus() const { return rowCpus_; }
  int colCpus() const { return colCpus_; }
  int localCpuRow() const { return localCpuRow_; }
  int localCpuCol() const { return localCpuCol_; }

  // Blocking factor (common to all arrays and dimensions)
  int blockSize() const { return blockSize_; }
  void blockSizeIs(int size);

  // Problem size
  int vectorSize()  const { return vectorSize_;  }
  int vectorCount() const { return vectorCount_; }
  int singularValueCount() const { return singularValueCount_; }

  void problemSizeIs(int vSize, int vCount);

  // Local data distribution
  int localRows() const { return localRows_; }
  int localCols() const { return localCols_; }

  // Global/Local (zero-based) index mapping
  // Local to global
  int globalRowIdx(int localIdx) const;
  int globalColIdx(int localIdx) const;

  // Global to local: host cpu
  int rowHostCpu(int globalIdx) const;
  int colHostCpu(int globalIdx) const;

  // Global to local: local index on host cpu
  int localRowIdx(int globalIdx) const;
  int localColIdx(int localIdx) const;

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by localCols]
  Scalar matrixEntry(int row, int col) const;
  const Scalar *matrixColBuffer(int col) const;
  const Scalar *matrixBuffer() const;
  
  Scalar &matrixEntry(int row, int col);
  Scalar *matrixColBuffer(int col);
  Scalar *matrixBuffer();
  
  int localMatrixLeadDim() const { return localMatrixLeadDim_; }
  
  // Local basis buffer: [localRows by localCols]
  Scalar basisEntry(int row, int col) const;
  const Scalar *basisColBuffer(int col) const;
  const Scalar *basisBuffer() const;
  
  int localBasisLeadDim() const { return localMatrixLeadDim_; }
 
  // Singular values: [vectorCount]
  Scalar singularValue(int rank) const;
  const Scalar *singularValueBuffer() const;

  // Cosntructor and destructor
  DistrSvdOrthogonalization(Communicator *comm, int rowCpus, int colCpus);
  ~DistrSvdOrthogonalization();

  // Solve
  void solve();

private:
  Communicator * communicator_;
  int blacsHandle_;

  typedef int Context;
  
  int rowCpus_, colCpus_;
  int localCpuRow_, localCpuCol_;
  Context context_;

  int blockSize_;
  
  typedef int ArrayDesc[9];

  int vectorSize_, vectorCount_, singularValueCount_;
  int localRows_, localCols_, localBasisCols_;
  int localMatrixLeadDim_;
  ArrayDesc matrixDesc_, basisDesc_;

  SimpleBuffer<Scalar> matrixBuffer_;
  SimpleBuffer<Scalar> basisBuffer_;
  SimpleBuffer<Scalar> sigmaBuffer_;

  // Private functions
  void reset();

  static const int DEFAULT_BLOCK_SIZE;
  
  // Addressable constants
  static const int INT_ZERO;
  static const int INT_ONE;
  static const int INT_MINUS_ONE;

  // Disallow copy & assignment
  DistrSvdOrthogonalization(const DistrSvdOrthogonalization &);
  DistrSvdOrthogonalization & operator=(const DistrSvdOrthogonalization &);
};


/* Helper functions for buffer access */

inline
const DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::matrixBuffer() const {
  return matrixBuffer_.array();
} 

inline
const DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::matrixColBuffer(int col) const {
  return matrixBuffer() + (col * localMatrixLeadDim());
}

inline
DistrSvdOrthogonalization::Scalar
DistrSvdOrthogonalization::matrixEntry(int row, int col) const {
  return matrixColBuffer(col)[row];
}

inline
DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::matrixBuffer() {
  return matrixBuffer_.array();
} 

inline
DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::matrixColBuffer(int col) {
  return matrixBuffer() + (col * localMatrixLeadDim());
}

inline
DistrSvdOrthogonalization::Scalar &
DistrSvdOrthogonalization::matrixEntry(int row, int col) {
  return matrixColBuffer(col)[row];
}

inline
const DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::basisBuffer() const {
  return basisBuffer_.array();
} 

inline
const DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::basisColBuffer(int col) const {
  return basisBuffer() + (col * localBasisLeadDim());
}

inline
DistrSvdOrthogonalization::Scalar
DistrSvdOrthogonalization::basisEntry(int row, int col) const {
  return basisColBuffer(col)[row];
}

inline
const DistrSvdOrthogonalization::Scalar *
DistrSvdOrthogonalization::singularValueBuffer() const {
  return sigmaBuffer_.array();
}

inline
DistrSvdOrthogonalization::Scalar
DistrSvdOrthogonalization::singularValue(int rank) const {
  return sigmaBuffer_[rank];
}

} // end namespace Rom

#endif /* ROM_DISTRSVDORTHOGONALIZATION_H */
