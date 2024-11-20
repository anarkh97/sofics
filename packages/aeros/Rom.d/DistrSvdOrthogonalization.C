#include "DistrSvdOrthogonalization.h"

#include <Comm.d/Communicator.h>

#include <Utils.d/linkfc.h>

#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <cassert>

#ifdef USE_SCALAPACK

#include <mpi.h>

extern "C" {
  // Context & cpu topology management
  int Csys2blacs_handle(MPI_Comm comm);
  void Cfree_blacs_system_handle(int handle);
  void Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);
  void Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myrow, int *mycol);
  void Cblacs_gridexit(int ictxt);

  // Block-cycling
  int _FORTRAN(numroc)(const int *n, const int *nb, const int *iproc, const int *isrcproc, const int *nprocs);
  void _FORTRAN(descinit)(int *desc, const int *m, const int *n, const int *mb, const int *nb,
                          const int *irsrc, const int *icsrc, const int *ictxt, const int *lld, int *info);

  // Index mapping 
  int _FORTRAN(indxl2g)(const int *indxloc, const int *nb, const int *iproc, const int *isrcproc, const int *nprocs);
  int _FORTRAN(indxg2l)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc_dummy, const int *nprocs);
  int _FORTRAN(indxg2p)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc, const int *nprocs);

  // Computations
  void _FORTRAN(pdgesvd)(const char *jobu, const char *jobvt, const int *m, const int *n,
                         double *a, const int *ia, const int *ja, const int desca[9],
                         double *s,
                         double *u, const int *iu, const int *ju, const int descu[9],
                         double *vt, const int *ivt, const int *jvt, const int descvt[9],
                         double *work, const int *lwork, int *info);
}

#endif /* USE_SCALAPACK */

namespace Rom {

const int DistrSvdOrthogonalization::DEFAULT_BLOCK_SIZE = 32;

const int DistrSvdOrthogonalization::INT_ZERO = 0;
const int DistrSvdOrthogonalization::INT_ONE = 1;
const int DistrSvdOrthogonalization::INT_MINUS_ONE = -1;

DistrSvdOrthogonalization::DistrSvdOrthogonalization(Communicator * comm, int rowCpus, int colCpus) : 
  communicator_(comm)
{
#ifdef USE_SCALAPACK
  if (rowCpus * colCpus > communicator_->numCPUs()) {
    throw std::invalid_argument("Not enough cpus"); 
  }
  
  if (rowCpus * colCpus < communicator_->numCPUs()) {
    throw std::invalid_argument("Too many cpus"); 
  }

  blacsHandle_ = Csys2blacs_handle(*communicator_->getCommunicator());
  context_ = blacsHandle_;
  char order[] = "R";
  Cblacs_gridinit(&context_, order, rowCpus, colCpus);
  Cblacs_gridinfo(context_, &rowCpus_, &colCpus_, &localCpuRow_, &localCpuCol_);
 
  assert(rowCpus_ == rowCpus);
  assert(colCpus_ == colCpus);

  // Empty problem
  vectorSize_ = 0;
  vectorCount_ = 0;
  singularValueCount_ = 0;

  blockSizeIs(DEFAULT_BLOCK_SIZE);

#else /* USE_SCALAPACK */
  throw std::runtime_error("ScaLAPACK not available");
#endif /* USE_SCALAPACK */
}

DistrSvdOrthogonalization::~DistrSvdOrthogonalization() {
#ifdef USE_SCALAPACK
  Cblacs_gridexit(context_);
  Cfree_blacs_system_handle(blacsHandle_);
#endif /* USE_SCALAPACK */
}

void
DistrSvdOrthogonalization::blockSizeIs(int size) {
  blockSize_ = size;

  reset(); 
}

void
DistrSvdOrthogonalization::problemSizeIs(int vSize, int vCount) {
  vectorSize_ = vSize;
  vectorCount_ = vCount;
  singularValueCount_ = std::min(vSize, vCount);

  reset();
}

void
DistrSvdOrthogonalization::reset() {
#ifdef USE_SCALAPACK
  localRows_ = _FORTRAN(numroc)(&vectorSize_, &blockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_);
  localCols_ = _FORTRAN(numroc)(&vectorCount_, &blockSize_, &localCpuCol_, &INT_ZERO, &colCpus_);
  localBasisCols_ = _FORTRAN(numroc)(&singularValueCount_, &blockSize_, &localCpuCol_, &INT_ZERO, &colCpus_);

  const int reqMatrixBufferSize = localCols_ * localRows_;
  const int reqBasisBufferSize = localBasisCols_ * localRows_;
  const int reqSigmaBufferSize = singularValueCount_;

  matrixBuffer_.sizeIs(reqMatrixBufferSize);
  basisBuffer_.sizeIs(reqBasisBufferSize);
  sigmaBuffer_.sizeIs(reqSigmaBufferSize);

  localMatrixLeadDim_ = std::max(localRows_, 1);

  int info;
  _FORTRAN(descinit)(matrixDesc_, &vectorSize_, &vectorCount_, &blockSize_, &blockSize_,
                     &INT_ZERO, &INT_ZERO, &context_, &localMatrixLeadDim_, &info);
  assert(info == 0);
  _FORTRAN(descinit)(basisDesc_, &vectorSize_, &singularValueCount_, &blockSize_, &blockSize_,
                     &INT_ZERO, &INT_ZERO, &context_, &localMatrixLeadDim_, &info);
  assert(info == 0);
#endif /* USE_SCALAPACK */
}

int
DistrSvdOrthogonalization::globalRowIdx(int localRowIdx) const {
#ifdef USE_SCALAPACK
  const int localRowIdx_f = localRowIdx + 1;
  return _FORTRAN(indxl2g)(&localRowIdx_f, &blockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_) - 1;
#else /* USE_SCALAPACK */
  return 0;
#endif /* USE_SCALAPACK */
}

int
DistrSvdOrthogonalization::globalColIdx(int localColIdx) const {
#ifdef USE_SCALAPACK
  const int localColIdx_f = localColIdx + 1;
  return _FORTRAN(indxl2g)(&localColIdx_f, &blockSize_, &localCpuCol_, &INT_ZERO, &colCpus_) - 1;
#else /* USE_SCALAPACK */
  return 0;
#endif /* USE_SCALAPACK */
}

int
DistrSvdOrthogonalization::rowHostCpu(int globalRowIdx) const {
#ifdef USE_SCALAPACK
  const int globalRowIdx_f = globalRowIdx + 1;
  return _FORTRAN(indxg2p)(&globalRowIdx_f, &blockSize_, NULL, &INT_ZERO, &rowCpus_);
#else /* USE_SCALAPACK */
  return 0;
#endif /* USE_SCALAPACK */
}

int
DistrSvdOrthogonalization::colHostCpu(int globalColIdx) const {
#ifdef USE_SCALAPACK
  const int globalColIdx_f = globalColIdx + 1;
  return _FORTRAN(indxg2p)(&globalColIdx_f, &blockSize_, NULL, &INT_ZERO, &colCpus_);
#else /* USE_SCALAPACK */
  return 0;
#endif /* USE_SCALAPACK */
}

int
DistrSvdOrthogonalization::localRowIdx(int globalRowIdx) const {
#ifdef USE_SCALAPACK
  const int globalRowIdx_f = globalRowIdx + 1;
  return _FORTRAN(indxg2l)(&globalRowIdx_f, &blockSize_, NULL, NULL, &rowCpus_) - 1;
#else /* USE_SCALAPACK */
  return 0;
#endif /* USE_SCALAPACK */
}

int
DistrSvdOrthogonalization::localColIdx(int globalColIdx) const {
#ifdef USE_SCALAPACK
  const int globalColIdx_f = globalColIdx + 1;
  return _FORTRAN(indxg2l)(&globalColIdx_f, &blockSize_, NULL, NULL, &colCpus_) - 1;
#else /* USE_SCALAPACK */
  return 0;
#endif /* USE_SCALAPACK */
}

void
DistrSvdOrthogonalization::solve() {
#ifdef USE_SCALAPACK
  const char *compute_u = "V";
  const char *no_vt = "N";
  int info;

  Scalar workspaceQuery;
  _FORTRAN(pdgesvd)(compute_u, no_vt, &vectorSize_, &vectorCount_,
                    matrixBuffer_.array(), &INT_ONE, &INT_ONE, matrixDesc_,
                    sigmaBuffer_.array(),
                    basisBuffer_.array(), &INT_ONE, &INT_ONE, basisDesc_,
                    NULL, NULL, NULL, NULL,
                    &workspaceQuery, &INT_MINUS_ONE, &info);
  assert(info == 0);

  const int lwork = static_cast<int>(workspaceQuery);
  SimpleBuffer<Scalar> workspace(lwork);
  _FORTRAN(pdgesvd)(compute_u, no_vt, &vectorSize_, &vectorCount_,
                    matrixBuffer_.array(), &INT_ONE, &INT_ONE, matrixDesc_,
                    sigmaBuffer_.array(),
                    basisBuffer_.array(), &INT_ONE, &INT_ONE, basisDesc_,
                    NULL, NULL, NULL, NULL,
                    workspace.array(), &lwork, &info);

  assert(info == 0);
#endif /* USE_SCALAPACK */
}
 
} // end namespace Rom
