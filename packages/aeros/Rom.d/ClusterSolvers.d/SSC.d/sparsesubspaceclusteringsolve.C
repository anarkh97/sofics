#if defined(USE_MPI) && defined(USE_SCALAPACK) && defined(USE_EIGEN3)
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <cassert>

#include "SSC.h"
#ifdef CLUSTER_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif

#include "Rom.d/DistrNonnegativeMatrixFactorization.h"

extern "C" {

   // this function computes all eigenvalues and eigenvectors 
   void _FORTRAN(pdsyev)(const char *jobz, const char *uplo, const int *n, const double *a, 
                         const int *ia, const int *ja, const int desca[9], double *w, double *z, 
                         const int *iz, const int *jz, const int descz[9], double *work, 
                         const int *lwork, int *info);
}


int
SparseSubspaceClustering::sparsesubspacecluster(SCDoubleMatrix & snapshots) {

    // take the snapshots and compute connectivity
    // this class is defined in the Rom namespace, for you noobs, that means you must use namespace delimiter
    Rom::DistrNonnegativeMatrixFactorization solver(_otherComm, snapshots.getNumberOfRows(), snapshots.getNumberOfCols(), 
                                                    snapshots.getNumberOfRowsLocal(), snapshots.getNumberOfRows(), 
                                                    snapshots.getRowBlockingFactor(), snapshots.getNumberOfRows()*10, 
                                                    sparseTol, 1, 0, 0, 0.0);

    SCDoubleMatrix ConnectivityGraph(_context, snapshots.getNumberOfCols(), snapshots.getNumberOfCols(), _mb, _nb, _comm); // allocate space for snapshot connectivity
    ConnectivityGraph.zero();

    //snapshots.normalizeColumns();     

    if (_mypid == 0) {
      std::cout << " Building Connectivity Graph" << std::endl;
    }

    startTime(TIME_SPARSESUBSPACECLUSTERING_MAIN_LOOP);
    {
      if(_mypid == 0) std::cout << " Making copy of Data" << std::endl; 
      SCDoubleMatrix CopyMat(snapshots,true); // copy for targets 
      if(_mypid == 0) std::cout << " Done, entering solver" << std::endl;
      solver.solveNNLS_MRHS(snapshots, CopyMat, ConnectivityGraph, 0, 1); // call NNLS for multiple right hand sides to get connectivity
    }

    if (_mypid == 0) {
      std::cout << " Symmetrizing Snapshot Connectivity" << std::endl;
    }

    // write connectivity graph for visualization
    ConnectivityGraph.write("graphfile.txt");

    // symmetrize connectivity graph
    if(_mypid == 0) std::cout << " Making copy of graph" << std::endl;
    SCDoubleMatrix CopyGraph(ConnectivityGraph,true); 
    if(_mypid == 0) std::cout << " Done. Symmetrizing" << std::endl;
    CopyGraph.add(ConnectivityGraph,'T',ConnectivityGraph.getNumberOfRows(), ConnectivityGraph.getNumberOfCols(), 1.0, 1.0); // Symmetrize graph

    ConnectivityGraph.write("symmetrizedGraphfile.txt");

    // compute Laplacian matrix L = D^(-1/2)*A*D^(-1/2) -> required for eigevectors to cluster on surface of hypersphere
    if (_mypid == 0) {
      std::cout << " Forming Laplacian Matrix" << std::endl;
    }
    ConnectivityGraph.computeLaplacian('1','Y');
    ConnectivityGraph.write("Laplacian.txt");

    // Now compute eigenvectors of Laplacian matrix
    if (_mypid == 0) {
      std::cout << " Spectral clustering of Laplacian Matrix" << std::endl;
    }
    computeEigenvectors(ConnectivityGraph);

    int status = sparsesubspaceclusterInit(*_eigVectors);// initialize clusters for eigenvector rows computed in last function call
    if (status != 0) {
        if (_mypid == 0) {
            std::cout << " SparseSubspaceClustering initialization failed." << std::endl;
        }
        return 1;
    }

    // perform kmeans clustering on normalized rows of eigenvectors
    if (_mypid == 0) {
      std::cout << " Kmeans clustering of Largest Laplacian Eigenvectors" << std::endl;
    }
    bool done = false;
    _snapshotCentroids->zero();
    SCIntMatrix snapshotCentroidsOld = SCIntMatrix(*_snapshotCentroids);
    _iter = 0;
    while (!done && _iter < _max_iter) {
        tagSnapshot(*_eigVectors);
        getCentroids(*_eigVectors, *_evCentroids);
        if (_snapshotCentroids->isEqual(snapshotCentroidsOld) == 0) {
            done = true;
        } else {
            _snapshotCentroids->copy(snapshotCentroidsOld);
        }
        output();
        _iter++;
    }

    // take clustering and compute centroids of high dimensional vectors
    getCentroids(snapshots, *_centroids);
    estimateClusterDimensionality(CopyGraph);

    stopTime(TIME_SPARSESUBSPACECLUSTERING_MAIN_LOOP);
    _wallclock_total[TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE_COPY] = snapshots.getTime(SCDBL_TIME_GETL2COLDIST_COPY);
    _wallclock_total[TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE_ADD] = snapshots.getTime(SCDBL_TIME_GETL2COLDIST_ADD);
    // write partitioning vector for visualization
    _snapshotCentroids->write("snapshotCentroids.txt");
    return 0;
}

void
SparseSubspaceClustering::tagSnapshot(SCDoubleMatrix & snapshots) {
    startTime(TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS);
    for (int is=1; is<=_n; is++) {
        double distmin = SPARSESUBSPACECLUSTERING_LARGE;
        int icmin = -1;
        startTime(TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE);
        for (int ic=1; ic<=_numClusters; ic++) {
            double dist = snapshots.getL2ColDistance(is, *_evCentroids, ic); // Fortran indices required, compute distance to centroids
            if (dist < distmin) { // find closest centroid
                distmin = dist;
                icmin = ic; 
            }
        }
        stopTime(TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE);
        _snapshotCentroids->setElement(is,1,icmin); // Fortran indices required, set the indice for closest centroid
    }
    stopTime(TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS);
}

void
SparseSubspaceClustering::computeEigenvectors(SCDoubleMatrix & graph) {

    const char *compute_vectors = "V"; // V tells scala to compute eigenvalues and eigenvectors
    const char *compute_values  = "N"; // N tells scala to compute eigenvalues only
    const char *use_lowerupper  = "L"; // L tells scala to use lower part, U for upper

    //double * W = new double[graph.getNumberOfRows()]; // container for eigenvalues
    std::vector<double> W(graph.getNumberOfRows(),0.0);

    int info;
    double workspaceQuery;
    const int MINUS_ONE = -1;
    const int N         = graph.getNumberOfRows();     // size of matrix whose eigenvalues we compute
    const int IA        = 1; // global row index for beggining of submatrix, use 1 for whole matrix
    const int JA        = 1; // global row index for beggining of submatrix, use 1 for whole matrix

    SCDoubleMatrix EVContainer(graph); EVContainer.zero();

    // the first call only computes the amount of workspace needed. 
    _FORTRAN(pdsyev)(compute_vectors, use_lowerupper, &N, graph.getMatrix(), &IA, &JA, 
                     graph.getDesc(), &W[0], EVContainer.getMatrix(), &IA, &JA, EVContainer.getDesc(), 
                     &workspaceQuery, &MINUS_ONE, &info);

    assert(info == 0);

    {
      if (_mypid == 0) {
        std::cout << " Computing Eigenvalue range" << std::endl;
      }
      // the second call actually computes the eigenvalues, upper(lower) part of input matrix is destroyed on exit (see scalapack documentation)
      const int lwork = static_cast<int>(workspaceQuery);
      double * workspace = new double[lwork];
      _FORTRAN(pdsyev)(compute_vectors, use_lowerupper, &N, graph.getMatrix(), &IA, &JA, 
                       graph.getDesc(), &W[0], EVContainer.getMatrix(), &IA, &JA, EVContainer.getDesc(), 
                       workspace, &lwork, &info);
      assert(info == 0);
  
      delete[] workspace;
    }

    double dWold = W[N-1] - W[N-2];
    double dWnew, dWmax = dWold; 
    if(_numClusters < 2) { // simple logic to determine optimial number of clusters
      for(int row = 0; row < N; row++){
        dWnew = W[N-1-row] - W[N-2-row];
        if(dWnew > dWmax) dWmax= dWnew;
        if(dWnew < 0.1*dWmax){
          _numClusters = row + 1;
          if(_mypid == 0) std::cout << "Dividing data into " << _numClusters << " clusters" << std::endl;
          break;
        } else {
          dWold = dWnew;
        }
      }
    }

    // print eigenvalues to screen
    if(_mypid == 0)
      for(int asdf = 0; asdf < _numClusters+5; asdf++) std::cout << " eigenvalue " << asdf << " is " << W[N-1-asdf] << std::endl;

    // store eigenvectors in transposed format
    _eigVectors = new SCDoubleMatrix(_context, _numClusters, graph.getNumberOfRows(), _mb, _nb, _comm); // number of clusters x number of snapshots
    _eigVectors->zero();
    // extract last _numClusters numbero of eigenvectors (associated with largest eigenvalues)
    EVContainer.add(*_eigVectors, 'T', _numClusters, EVContainer.getNumberOfRows(), 1.0, 1.0,1,N-_numClusters+1,1,1);
    _eigVectors->normalizeColumns(); // eigenvector rows need to be normalized   

    // write normalized eigenvectors to visualize on hyper-sphere
    _eigVectors->write("normalizedEigenVectors.txt");
 
}

void
SparseSubspaceClustering::getCentroids(SCDoubleMatrix & snapshots, SCDoubleMatrix & _centroidContainer) {
    
    startTime(TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS);
    int nss;
    for (int ic=1; ic<=_numClusters; ic++) { // Fortran Index, loop over clusters
        nss = _snapshotCentroids->countValue(ic); // count how many times the values ic shows up in this vector
        double fac = 1.0/((float) nss);
        startTime(TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS_SUM);
        snapshots.sumOfColumns(_centroidContainer, ic, *_snapshotCentroids, ic, fac); // add all columns associated with ic and average
        stopTime(TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS_SUM);
    }
    stopTime(TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS);
}

void
SparseSubspaceClustering::estimateClusterDimensionality(SCDoubleMatrix &graph){

    for (int ic=1; ic<=_numClusters; ic++) { // loop over clusters
      int nnzmax = 0; 
      for (int row=1; row<=graph.getNumberOfCols(); row++){ // loop over rows
        int rowmask = _snapshotCentroids->getElement(row,1);
        _snapshotCentroids->SCBaseMatrix::distributeVector(&rowmask,1);
        int nnzrow = 0; 
        if(rowmask == ic){ // check which cluster this row is in
          for(int col=1; col<=graph.getNumberOfRows(); col++){
            int colmask = _snapshotCentroids->getElement(col,1);
            _snapshotCentroids->SCBaseMatrix::distributeVector(&colmask,1);
            if(colmask == ic) { // check which cluster this column is in
              double elem = graph.getElement(row,col); 
              graph.SCBaseMatrix::distributeVector(&elem,1);
              if(elem > 1e-16) nnzrow++; 
            }
          }
        }
        if(nnzrow>nnzmax) nnzmax = nnzrow; 
      }
      if(_mypid == 0) std::cout << "Maximum connectivity for cluster " << ic << " is " << nnzmax << std::endl;
    } 

}

int
SparseSubspaceClustering::sparsesubspaceclusterInit(SCDoubleMatrix & snapshots) {
    //std::cout << "Entering sparsesubspaceclusternit" << std::endl;
    startTime(TIME_SPARSESUBSPACECLUSTERING_INIT); 
    
    // check that the number of clusters asked for is smaller than number of snapshots
    if (_numClusters > _n) {
        std::cout << "Not enough snapshots for " << _numClusters;
        std::cout << " clusters. Maximun is " << _n << "." << std::endl;
        return 1;
    }
 
    // initialize random number generator
    srandom(_seed); 

    // allocate space for the centroid labels of the snapshots
    _snapshotCentroids = new SCIntMatrix(_context, _n, 1, _mb, _nb, _comm);  // number of snapshots x 1
    int nb = _n / _npcol;
    if (nb == 0) nb = 1;
    _centroids         = new SCDoubleMatrix(_context, _m, _numClusters, _mb, _nb, _comm);           // snapshot length x number of clusters
    _evCentroids       = new SCDoubleMatrix(_context, _numClusters, _numClusters, _mb, _nb, _comm); // number of clusters x number of clusters


    // Get initial centroids
    int one = 1, jx, ic; 
    int * centroidIndices = new int[_numClusters];
    int j = random() % _n;
    for (int i=0; i<_numClusters; i++) {
        bool unique = false;
        while (i>0 && !unique) {
            j = random() % _n;
            unique = true;
            for (int k=0; k<i; k++) {
                if (centroidIndices[k] == j) {
                    unique = false;
                    break;
                }
            }
        }
        centroidIndices[i] = j;
        jx = j+1;
        ic = i+1;
        snapshots.pdcopy(snapshots.getNumberOfRows(), one, jx, one, *_evCentroids, one, ic, one);
    }

    printCentroidIndices(centroidIndices);
    delete[] centroidIndices;
    stopTime(TIME_SPARSESUBSPACECLUSTERING_INIT);
    return 0;
}


void
SparseSubspaceClustering::printCentroidIndices(int * centroidIndices) {
    if (_mypid == 0) {
        std::cout << " Initial cluster centroid indices are: ";
        for (int i=0; i<_numClusters; i++) {
            std::cout << centroidIndices[i];
            if (i < _numClusters-1) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl << std::endl;
    }
}

#endif
