#if defined(USE_MPI)
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <algorithm>

#include "Kmeans.h"
#ifdef CLUSTER_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif

int
Kmeans::kmeanscluster(SCDoubleMatrix & snapshots) {
    startTime(TIME_KMEANS_MAIN_LOOP);
    int status = kmeansInit(snapshots);
    if (status != 0) {
        if (_mypid == 0) {
            std::cout << "Kmeans initialization failed." << std::endl;
        }
        return 1;
    }
    bool done = false;
    _snapshotCentroids->zero();
    SCIntMatrix snapshotCentroidsOld = SCIntMatrix(*_snapshotCentroids);
    _iter = 0;
    while (!done && _iter < _max_iter) {
        tagSnapshot(snapshots);
        getCentroids(snapshots);
        if (_snapshotCentroids->isEqual(snapshotCentroidsOld) == 0) {
            done = true;
        } else {
            _snapshotCentroids->copy(snapshotCentroidsOld);
        }
        output();
        _iter++;
    }
    stopTime(TIME_KMEANS_MAIN_LOOP);
    _wallclock_total[TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE_COPY] = snapshots.getTime(SCDBL_TIME_GETL2COLDIST_COPY);
    _wallclock_total[TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE_ADD] = snapshots.getTime(SCDBL_TIME_GETL2COLDIST_ADD);
     _centroids->write("centroids.txt");
    _snapshotCentroids->write("snapshotCentroids.txt");
    return 0;
}


void
Kmeans::tagSnapshot(SCDoubleMatrix & snapshots) {
    //std::cout << "Entering tagSnapshot" << std::endl;
    startTime(TIME_KMEANS_TAG_SNAPSHOTS);
    for (int is=1; is<=_n; is++) {
        double distmin = KMEANS_LARGE;
        int icmin = -1;
        startTime(TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE);
        for (int ic=1; ic<=_numClusters; ic++) {
            double dist = snapshots.getL2ColDistance(is, *_centroids, ic); // Fortran indices required
            if (dist < distmin) {
                distmin = dist;
                icmin = ic;
            }
        }
        stopTime(TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE);
        //_snapshotCentroids->setElement(1,is,icmin); // Fortran indices required
        _snapshotCentroids->setElement(is,1,icmin); // Fortran indices required
    }
    //_snapshotCentroids->write("snapshotCentroids.txt");
    stopTime(TIME_KMEANS_TAG_SNAPSHOTS);
    //std::cout << "Exiting tagSnapshot" << std::endl;
}


void
Kmeans::getCentroids(SCDoubleMatrix & snapshots) {
    //std::cout << "Entering getCentroids" << std::endl;
    startTime(TIME_KMEANS_GET_CENTROIDS);
    int nss;
    
for (int ic=1; ic<=_numClusters; ic++) { // Fortran Index
        nss = _snapshotCentroids->countValue(ic);
        double fac = 1.0/((float) nss);
        startTime(TIME_KMEANS_GET_CENTROIDS_SUM);
        snapshots.sumOfColumns(*_centroids, ic, *_snapshotCentroids, ic, fac);
        //snapshots.sumOfColumns(*_centroids, ic, *_snapshotCentroids, ic);
        stopTime(TIME_KMEANS_GET_CENTROIDS_SUM);
        //double fac = 1.0/((float) nss);
        //_centroids->add(*_centroids, 'N', _m, 1, fac, 0.0, 1, ic, 1, ic);
    }
    stopTime(TIME_KMEANS_GET_CENTROIDS);
    //std::cout << "Exiting getCentroids" << std::endl;
}


int
Kmeans::kmeansInit(SCDoubleMatrix & snapshots) {
    //std::cout << "Entering kmeansInit" << std::endl;
    startTime(TIME_KMEANS_KMEANS_INIT);
    if (_numClusters > _n) {
        std::cout << "Not enough snapshots for " << _numClusters;
        std::cout << " clusters. Maximun is " << _n << "." << std::endl;
        return 1;
    }
    srandom(_seed); 
    //_snapshotCentroids = new SCIntMatrix(_context, 1, _n, _mb, _nb, _comm);
    _snapshotCentroids = new SCIntMatrix(_context, _n, 1, _mb, _nb, _comm);
    int nb = _n / _npcol;
    if (nb == 0) nb = 1;
    _centroids = new SCDoubleMatrix(_context, _m, _numClusters, _mb, _nb, _comm);

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
        snapshots.pdcopy(_m, one, jx, one, *_centroids, one, ic, one);
    }
    printCentroidIndices(centroidIndices);
    //_centroids->write("centroids.txt");
    delete[] centroidIndices;
    stopTime(TIME_KMEANS_KMEANS_INIT);
    //std::cout << "Exiting kmeansInit" << std::endl;
    return 0;
}


void
Kmeans::printCentroidIndices(int * centroidIndices) {
    if (_mypid == 0) {
        std::cout << "Initial cluster centroid indices are: ";
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
