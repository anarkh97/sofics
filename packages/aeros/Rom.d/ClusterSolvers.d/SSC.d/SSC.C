#if defined(USE_MPI)
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <algorithm>

#include "SSC.h"
#ifdef CLUSTER_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif


// constructor     
SparseSubspaceClustering::SparseSubspaceClustering() {
    setDefaults();
}


// destructor
SparseSubspaceClustering::~SparseSubspaceClustering() {
    delete _centroids;
    delete _evCentroids;
    delete _snapshotCentroids;
    delete _eigVectors;
}


void
SparseSubspaceClustering::setDefaults() {
    _seed = 1;
    _numClusters = 10;
    _max_iter = SPARSESUBSPACECLUSTERING_MAX_ITER_DEFAULT;
    sparseTol = 0.01;
}


int
SparseSubspaceClustering::close() {
    MPI_Barrier(MPI_COMM_WORLD);
    int one=1;
    _FORTRAN(blacs_exit)(&one);
    return 0;
}


double
SparseSubspaceClustering::getWallTime(){
    struct timeval time;
    gettimeofday(&time,NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}



void
SparseSubspaceClustering::output() {
    if (_mypid == 0) {
        if (_iter%SPARSESUBSPACECLUSTERING_HEADER_INCR == 0) {
            header();
        }
        printf("   %6d\n", _iter);
    }
}


void
SparseSubspaceClustering::header() {
    std::cout << "   ";
    std::cout << "iter   ";
    std::cout << std::endl;

    std::cout << "-- ";
    std::cout << "------ ";
    std::cout << std::endl;
}


void
SparseSubspaceClustering::printTimes() {
    double times_max[N_TIMES_SPARSESUBSPACECLUSTERING];
    MPI_Allreduce(_wallclock_total, times_max, N_TIMES_SPARSESUBSPACECLUSTERING, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int i=0; i<N_TIMES_SPARSESUBSPACECLUSTERING; i++) {
        _wallclock_total[i] = times_max[i];
    }   

    if (_mypid == 0) {
        std::cout << std::endl;
        std::cout << "Wallclock Times (seconds):"                                               << std::endl;
        std::cout << "    Solver                  : " << times_max[TIME_SPARSESUBSPACECLUSTERING_MAIN_LOOP]      << std::endl;
        std::cout << "        kmeans init         : " << times_max[TIME_SPARSESUBSPACECLUSTERING_INIT]    << std::endl;
        std::cout << "        tag snapshots       : " << times_max[TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS]  << std::endl;
        std::cout << "           compute l2 norms : " << times_max[TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE]  << std::endl;
        std::cout << "           copy             : " << times_max[TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE_COPY]  << std::endl;
        std::cout << "           add              : " << times_max[TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE_ADD]  << std::endl;
        std::cout << "        get centroids       : " << times_max[TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS]  << std::endl;
        std::cout << "           sum              : " << times_max[TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS_SUM]  << std::endl;
    }
}


void
SparseSubspaceClustering::summary(SCDoubleMatrix& A) {
    if (_mypid == 0) {
        std::cout << std::endl;
        std::cout << "################### Solver Summary ########################" << std::endl;
        std::cout << "Solver: Sparse Subspace Clustering with " << _numClusters << " clusters." << std::endl;
        std::cout << "Matrix size: " << _m  << " x " << _n << std::endl;
        std::cout << "Blacs context:"                    << std::endl;
        std::cout << "    mprow        = " << _mprow      << std::endl;
        std::cout << "    npcol        = " << _npcol      << std::endl;
        std::cout << "Scalapack blocking factors:"       << std::endl;
        std::cout << "    mb           = " << _mb         << std::endl;
        std::cout << "    nb           = " << _nb         << std::endl;
        std::cout << "############################################################" << std::endl << std::endl;
    }
}

int
SparseSubspaceClustering::init(SCDoubleMatrix& A) {
    // Define MPI Process Grid if not defined
    _FORTRAN(blacs_pinfo)(&_mypid, &_nprocs);

    // This is just for convenience
    _m  = A.getNumberOfRows();
    _n  = A.getNumberOfCols();
    _mb = A.getRowBlockingFactor();
    _nb = A.getColBlockingFactor();
    _context = A.getContext();
    _comm = A.getMpiComm();
    _FORTRAN(blacs_gridinfo)(&_context, &_mprow, &_npcol, &_myrow, &_mycol);

    for (int i=0; i<N_TIMES_SPARSESUBSPACECLUSTERING; i++) {
        _wallclock[i] = 0.0;
        _wallclock_total[i] = 0.0;
    }
    return 0;
}


int
SparseSubspaceClustering::cluster(SCDoubleMatrix& A) {
    init(A);
    summary(A);
    sparsesubspacecluster(A);
    return 0;
}


void
SparseSubspaceClustering::getClusterColumns(std::vector<std::vector<int> >& clusterCols) {
    _snapshotCentroids->setScope('A');
    for (int i=0; i<_n; i++) {
        int icluster = _snapshotCentroids->getElement(i+1,1)-1;  // _snapshotCentroids holds Fortran indices
        clusterCols[icluster].push_back(i);
    }
    _snapshotCentroids->setScope();
}

#endif
