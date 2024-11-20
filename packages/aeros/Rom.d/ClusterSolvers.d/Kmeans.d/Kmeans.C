#if defined(USE_MPI)
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <algorithm>

#include "Kmeans.h"
#ifdef CLUSTER_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif



Kmeans::Kmeans() {
    setDefaults();
}


Kmeans::~Kmeans() {
    delete _centroids;
    delete _snapshotCentroids;
}


void
Kmeans::setDefaults() {
    _seed = 1;
    _numClusters = 10;
    _max_iter = KMEANS_MAX_ITER_DEFAULT;
}


int
Kmeans::close() {
    MPI_Barrier(MPI_COMM_WORLD);
    int one=1;
    _FORTRAN(blacs_exit)(&one);
    return 0;
}


double
Kmeans::getWallTime(){
    struct timeval time;
    gettimeofday(&time,NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


std::string
Kmeans::intToString(int i) {
    std::ostringstream ss;
    ss << i;
    return ss.str();
}



void
Kmeans::output() {
    if (_mypid == 0) {
        if (_iter%KMEANS_HEADER_INCR == 0) {
            header();
        }
        printf("   %6d\n", _iter);
    }
}


void
Kmeans::header() {
    std::cout << "   ";
    std::cout << "iter   ";
    std::cout << std::endl;

    std::cout << "-- ";
    std::cout << "------ ";
    std::cout << std::endl;
}


void
Kmeans::printTimes() {
    double times_max[N_TIMES_KMEANS];
    MPI_Allreduce(_wallclock_total, times_max, N_TIMES_KMEANS, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int i=0; i<N_TIMES_KMEANS; i++) {
        _wallclock_total[i] = times_max[i];
    }   

    if (_mypid == 0) {
        std::cout << std::endl;
        std::cout << "Wallclock Times (seconds):"                                               << std::endl;
        std::cout << "    Solver                  : " << times_max[TIME_KMEANS_MAIN_LOOP]      << std::endl;
        std::cout << "        kmeans init         : " << times_max[TIME_KMEANS_KMEANS_INIT]    << std::endl;
        std::cout << "        tag snapshots       : " << times_max[TIME_KMEANS_TAG_SNAPSHOTS]  << std::endl;
        std::cout << "           compute l2 norms : " << times_max[TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE]  << std::endl;
        std::cout << "           copy             : " << times_max[TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE_COPY]  << std::endl;
        std::cout << "           add              : " << times_max[TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE_ADD]  << std::endl;
        std::cout << "        get centroids       : " << times_max[TIME_KMEANS_GET_CENTROIDS]  << std::endl;
        std::cout << "           sum              : " << times_max[TIME_KMEANS_GET_CENTROIDS_SUM]  << std::endl;
    }
}


void
Kmeans::summary(SCDoubleMatrix& A) {
    if (_mypid == 0) {
        std::cout << std::endl;
        std::cout << "################### Solver Summary ########################" << std::endl;
        std::cout << "Solver: k-means with " << _numClusters << " clusters." << std::endl;
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


void
Kmeans::writeSolution(bool compact) {
}


int
Kmeans::init(SCDoubleMatrix& A) {
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

    for (int i=0; i<N_TIMES_KMEANS; i++) {
        _wallclock[i] = 0.0;
        _wallclock_total[i] = 0.0;
    }
    return 0;
}


int
Kmeans::cluster(SCDoubleMatrix& A) {
    init(A);
    summary(A);
    kmeanscluster(A);
    return 0;
}


void
Kmeans::getClusterColumns(std::vector<std::vector<int> >& clusterCols) {
    _snapshotCentroids->setScope('A');
    for (int i=0; i<_n; i++) {
        int icluster = _snapshotCentroids->getElement(i+1,1)-1;  // _snapshotCentroids holds Fortran indices
        clusterCols[icluster].push_back(i);
    }
    _snapshotCentroids->setScope();
}

/*
void
Kmeans::getCluster(std::vector<std::vector<int> > clusterCols, Eigen::MatrixXd centroidBuffer) {
    _snapshotCentroids.setScope('A');
    for (int i=0; i<_n; i++) {
        int icluster = _snapshotCentroids.getElement(i+1,1)-1;  // _snapshotCentroids holds Fortran indices
        clusterCols[icluster].push_back(i-1);
    }
    _snapshotCentroids.setScope();

   nrows_local = centroidBuffer.rows();
   int row_lengths = new int[_mprow];
   char scope = 'C';
   char top = ' ';
   if (_mycol == 0) {
       for (int i=0; i<_mprow; i++) {
            if (_myrow == i) {
                row_lengths[i] = nrows_local;
                _FORTRAN(dgebs2d)(&_context, &scope, &top, &one, &one, &nrows_local, &one);
            } else {
                _FORTRAN(dgebr2d)(&_context, &scope, &top, &one, &one, &(row_lengths[i]), &one, &i, &zero);
            }
        }
   }
}
*/


#endif
