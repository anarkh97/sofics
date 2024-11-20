#ifndef KMEANS_H_
#define KMEANS_H_

#include <mpi.h>
#include <cstring>
#include <vector>


#ifdef CLUSTER_DEV
#include "SCDoubleMatrix.h"
#include "SCIntMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCDoubleMatrix.h"
#include "Math.d/SCMatrix.d/SCIntMatrix.h"
#endif


#define KMEANS_HEADER_INCR 30
#define KMEANS_RESIDUAL_NO_DECREASE_REL_TOL 1.0E-14
#define KMEANS_RESIDUAL_NO_DECREASE_COUNT 5

#define TIME_KMEANS_MAIN_LOOP 0
#define TIME_KMEANS_TAG_SNAPSHOTS 1
#define TIME_KMEANS_GET_CENTROIDS 2
#define TIME_KMEANS_KMEANS_INIT 3
#define TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE 4
#define TIME_KMEANS_GET_CENTROIDS_SUM 5
#define TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE_COPY 6
#define TIME_KMEANS_TAG_SNAPSHOTS_DISTANCE_ADD 7
#define N_TIMES_KMEANS 8

#define KMEANS_LARGE 1.0e+300
#define KMEANS_MAX_ITER_DEFAULT 1000


class Kmeans {

    public:
        Kmeans();
        ~Kmeans();

        int setMatrixRow(int j, double *row);
        int setMatrixColumn(int j, double *col);
        int cluster(SCDoubleMatrix& snapshots);
        int close();

        //int getContext() {return _context;};
        double getTime(int i) { return _wallclock_total[i]; }
        double getComputeTime() {return getTime(TIME_KMEANS_MAIN_LOOP);}
        void summary(SCDoubleMatrix& A);
        void printTimes();
        int getStatus() {return _status;}
        void setNumClusters( double numClusters ) {_numClusters = numClusters;}

        void setMaxIter(int max_iter) {_max_iter = max_iter;}
        void writeSolution(bool compact=false);
        double getWallTime();
        void startTime(int i) {_wallclock[i] = -getWallTime();}
        void stopTime(int i)  {_wallclock[i] += getWallTime(); _wallclock_total[i] += _wallclock[i];}
        void printCentroidIndices(int * centroidIndices);
        void setSeed(unsigned seed) {_seed = seed;}
        void getClusterColumns(std::vector<std::vector<int> >& clusterCols);

    private:
        std::string _stopString;

        int _numClusters;

        // Timings
        double _wallclock[N_TIMES_KMEANS];
        double _wallclock_total[N_TIMES_KMEANS];

        int _max_iter;

        int _m;
        int _n;
        int _mb;
        int _nb;
        int _mprow;
        int _npcol;
        int _myrow;
        int _mycol;
        int _context;
        MPI_Comm _comm;

        int _nprocs;
        int _mypid;
        int _iter;
        int _status;
        unsigned _seed;
        SCDoubleMatrix * _centroids;
        SCIntMatrix * _snapshotCentroids;

        int init(SCDoubleMatrix & A);
        int kmeanscluster(SCDoubleMatrix & snapshots);
        void setDefaults();
        std::string intToString(int);
        void output();
        void header();
        int kmeansInit(SCDoubleMatrix & snapshots);
        void tagSnapshot(SCDoubleMatrix & snapshots);
        void getCentroids(SCDoubleMatrix & snapshots);
};

#endif // KMEANS_H_
