#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
#include "DistrSnapshotRowClusteringSolver.h"

#include <Comm.d/Communicator.h>
#include <Math.d/SCMatrix.d/SCDoubleMatrix.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/linkfc.h>
#include <Rom.d/ClusterSolvers.d/Kmeans.d/Kmeans.h>
#include <Rom.d/ClusterSolvers.d/SSC.d/SSC.h>

#include <mpi.h>

extern "C" {
  // Context & cpu topology management
  int Csys2blacs_handle(MPI_Comm comm);

  void Cfree_blacs_system_handle(int handle);
  void Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);

  void Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myrow, int *mycol);
  void Cblacs_gridexit(int ictxt);

  // Index mapping 
  int _FORTRAN(indxg2l)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc_dummy, const int *nprocs);
  int _FORTRAN(indxg2p)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc, const int *nprocs);
}

namespace Rom {

DistrSnapshotRowClusteringSolver
::DistrSnapshotRowClusteringSolver(Communicator * comm, int rowCount, int colCount, int localRows, int numClusters,
                                int blockSize) : 
  communicator_(comm),
  rowCount_(rowCount),
  colCount_(colCount),
  localRows_(localRows),
  numClusters_(numClusters),
  blockSize_(blockSize),
  matrixBuffer_(localRows,colCount),
  clusterLocalRows_(numClusters),
  solverType_(1),
  kmMaxIter_(1000),
  kmSeed_(1)
{
  timeStamps_.resize(colCount);
}

void
DistrSnapshotRowClusteringSolver::solve()
{
  // block-cyclic redistribution of snapshot matrix
  int blacsHandle = Csys2blacs_handle(*communicator_->getCommunicator());
  int context = blacsHandle;
  char order[] = "R";
  int rowCpus = 1; // XXX
  int colCpus = communicator_->numCPUs();
  Cblacs_gridinit(&context, order, rowCpus, colCpus);

  SCDoubleMatrix Xtranspose(context, colCount_, rowCount_, blockSize_, blockSize_, *communicator_->getCommunicator()); // global matrix

  int localBlacsHandle = Csys2blacs_handle(MPI_COMM_SELF);
  int localContext = localBlacsHandle;
  Cblacs_gridinit(&localContext, order, 1, 1);

  SCDoubleMatrix Xtransposei(localContext, colCount_, localRows_, blockSize_, blockSize_, MPI_COMM_SELF); // local matrix

  int *localRows = new int[communicator_->numCPUs()];
  communicator_->allGather(&localRows_, 1, localRows, 1);

  int offset;
  for(int i=0,jb=1; i<communicator_->numCPUs(); i++) {
    if(i == communicator_->myID()) {
      for(int k=0; k<colCount_; ++k) {
        Xtransposei.setMatrixRow(k+1, matrixColBuffer(k));
      }
      Xtransposei.copyRedist(colCount_, localRows_, 1, 1, Xtranspose, 1, jb, context);
      offset = jb-1;
    }
    else {
      SCDoubleMatrix::copyRedist(colCount_, localRows[i], Xtranspose, 1, jb, context);
    }
    jb += localRows[i]; //ia += localRows[i];
  }

  /*double Xnorm = X.froNorm();
  for(int i=0; i<communicator_->numCPUs(); ++i) { if(i == communicator_->myID()) {
  std::cerr << "matrixBuffer_.norm() = " << matrixBuffer_.norm() << std::endl;
  std::cerr << "X.froNorm() = " << Xnorm << std::endl;
  } communicator_->sync(); }*/

  Cblacs_gridexit(localContext);
  Cfree_blacs_system_handle(localBlacsHandle);

  std::vector<std::vector<int> > clusterRows(numClusters_);
  switch(solverType_) {
    default :
    case 0 : { // Random assignment 
      filePrint(stderr, " ... Using Random Clustering        ...\n");
      std::vector<int> clusterAssignment(rowCount_);
      for(int k=0; k<rowCount_; ++k) clusterAssignment[k] = rand()%numClusters_;
      // make a list of the rows assigned to each cluster
      for(int k=0; k<rowCount_; ++k) {
        int i = clusterAssignment[k];
        clusterRows[i].push_back(k);
      }
    } break;

    case 1: { // K-means
      filePrint(stderr, " ... Using K-means Clustering       ...\n");
      Kmeans solver = Kmeans();
      solver.setNumClusters(numClusters_);
      solver.setMaxIter(kmMaxIter_);
      solver.setSeed(kmSeed_);
      int status = solver.cluster(Xtranspose);
      solver.printTimes();
      solver.getClusterColumns(clusterRows);
    } break;

    case 2: { // sparse subspace clustering
      filePrint(stderr, " ... Using Sparse Subspace Clustering ...\n");
      SparseSubspaceClustering solver = SparseSubspaceClustering();
      solver.setNumClusters(numClusters_);
      solver.setMaxIter(kmMaxIter_);
      solver.setCommunicator(communicator_);
      solver.setSparseTolerance(nnlsTol);
      int status = solver.cluster(Xtranspose);
      solver.printTimes();
      int numClust = solver.getNumClusters();
      numClusters_ = numClust;
      // resize containers 
      clusterRows.resize(numClust);
      solver.getClusterColumns(clusterRows);
    } break;
  }

  for(int i = 0; i < numClusters_; ++i) {
    for(int j = 0; j < clusterRows[i].size(); ++j) {
      int k = clusterRows[i][j];
      if(offset <= k && k < offset+localRows[communicator_->myID()])
        clusterLocalRows_[i].push_back(k-offset);
    }
  }

  delete [] localRows;

  Cblacs_gridexit(context);
  Cfree_blacs_system_handle(blacsHandle);
}

} // end namespace Rom
#endif
