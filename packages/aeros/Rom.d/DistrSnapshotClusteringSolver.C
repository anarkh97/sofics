#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
#include "DistrSnapshotClusteringSolver.h"

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

DistrSnapshotClusteringSolver
::DistrSnapshotClusteringSolver(Communicator * comm, int rowCount, int colCount, int localRows, int numClusters,
                                int blockSize) : 
  communicator_(comm),
  rowCount_(rowCount),
  colCount_(colCount),
  localRows_(localRows),
  numClusters_(numClusters),
  blockSize_(blockSize),
  matrixBuffer_(localRows,colCount),
  centroidBuffer_(localRows,numClusters),
  clusterCols_(numClusters),
  solverType_(1),
  kmMaxIter_(1000),
  kmSeed_(1)
{
  timeStamps_.resize(colCount); 
}

void
DistrSnapshotClusteringSolver::solve()
{
  // block-cyclic redistribution of snapshot matrix
  int blacsHandle = Csys2blacs_handle(*communicator_->getCommunicator());
  int context = blacsHandle;
  char order[] = "R";
  int rowCpus = communicator_->numCPUs();
  int colCpus = 1; // XXX
  Cblacs_gridinit(&context, order, rowCpus, colCpus);

  SCDoubleMatrix X(context, rowCount_, colCount_, blockSize_, blockSize_, *communicator_->getCommunicator()); // global matrix

  int localBlacsHandle = Csys2blacs_handle(MPI_COMM_SELF);
  int localContext = localBlacsHandle;
  Cblacs_gridinit(&localContext, order, 1, 1);

  SCDoubleMatrix Xi(localContext, localRows_, colCount_, blockSize_, blockSize_, MPI_COMM_SELF); // local matrix

  int *localRows = new int[communicator_->numCPUs()];
  communicator_->allGather(&localRows_, 1, localRows, 1);

  for(int i=0,ia=1; i<communicator_->numCPUs(); i++) {
    if(i == communicator_->myID()) {
      for(int k=0; k<colCount_; ++k) { // XXX this copy could be avoided
        Xi.setMatrixColumn(k+1, matrixColBuffer(k));
      }
      Xi.copyRedist(localRows_, colCount_, 1, 1, X, ia, 1, context); 
    }
    else {
      SCDoubleMatrix::copyRedist(localRows[i], colCount_, X, ia, 1, context);
    }
    ia += localRows[i];
  }
  delete [] localRows;

  /*double Xnorm = X.froNorm();
  for(int i=0; i<communicator_->numCPUs(); ++i) { if(i == communicator_->myID()) {
  std::cerr << "matrixBuffer_.norm() = " << matrixBuffer_.norm() << std::endl;
  std::cerr << "X.froNorm() = " << Xnorm << std::endl;
  } communicator_->sync(); }*/

  Cblacs_gridexit(localContext);
  Cfree_blacs_system_handle(localBlacsHandle);

  switch(solverType_) {
    default :
    case 0 : { // Random assignment 
      filePrint(stderr, " ... Using Random Clustering        ...\n");
      std::vector<int> clusterAssignment(colCount_);
      for(int k=0; k<colCount_; ++k) clusterAssignment[k] = rand()%numClusters_;
      // make a list of the columns assigned to each cluster
      // and compute the centroids
      centroidBuffer_.setZero();
      for(int k=0; k<colCount_; ++k) {
        int i = clusterAssignment[k];
        clusterCols_[i].push_back(k);
        centroidBuffer_.col(i) += matrixBuffer_.col(k);
      }
      for(int i=0; i<numClusters_; ++i) {
        centroidBuffer_.col(i) /= clusterCols_[i].size();
      }
    } break;

    case 1: { // K-means
      filePrint(stderr, " ... Using K-means Clustering       ...\n");
      Kmeans solver = Kmeans();
      solver.setNumClusters(numClusters_);
      solver.setMaxIter(kmMaxIter_);
      solver.setSeed(kmSeed_);
      int status = solver.cluster(X);
      solver.printTimes();
      solver.getClusterColumns(clusterCols_);
      // Recomputing the centroids is easier/faster than mapping from SCDoubleMatrix to Eigen.
      centroidBuffer_.setZero();
      for(int i=0; i<numClusters_; ++i) {
        for (int j=0; j<clusterCols_[i].size(); j++) {
           centroidBuffer_.col(i) += matrixBuffer_.col(clusterCols_[i][j]);
        }
      }
      for(int i=0; i<numClusters_; ++i) {
        centroidBuffer_.col(i) /= clusterCols_[i].size();
      }
    } break;

    case 2: { // sparse subspace clustering
      filePrint(stderr, " ... Using Sparse Subspace Clustering ...\n");
      SparseSubspaceClustering solver = SparseSubspaceClustering();
      solver.setNumClusters(numClusters_);
      solver.setMaxIter(kmMaxIter_);
      solver.setCommunicator(communicator_);
      solver.setSparseTolerance(nnlsTol);
      int status = solver.cluster(X);
      solver.printTimes();
      int numClust = solver.getNumClusters();
      numClusters_ = numClust;
      // resize containers 
      clusterCols_.resize(numClust);
      centroidBuffer_.resize(localRows_,numClust);
      solver.getClusterColumns(clusterCols_);
      // Recomputing the centroids is easier/faster than mapping from SCDoubleMatrix to Eigen.
      centroidBuffer_.setZero();
      for(int i=0; i<numClusters_; ++i) {
        for (int j=0; j<clusterCols_[i].size(); j++) {
           centroidBuffer_.col(i) += matrixBuffer_.col(clusterCols_[i][j]);
        }
      }
      for(int i=0; i<numClusters_; ++i) {
        centroidBuffer_.col(i) /= clusterCols_[i].size();
      }
    } break;
  }

  Cblacs_gridexit(context);
  Cfree_blacs_system_handle(blacsHandle);
}

void
DistrSnapshotClusteringSolver::recomputeCentroids(){

  centroidBuffer_.setZero();
  for(int i=0; i<numClusters_; ++i) {
    for (int j=0; j<clusterCols_[i].size(); j++) {
       centroidBuffer_.col(i) += matrixBuffer_.col(clusterCols_[i][j]);
    }
  }
  for(int i=0; i<numClusters_; ++i) {
    centroidBuffer_.col(i) /= clusterCols_[i].size();
  }

}

} // end namespace Rom
#endif
