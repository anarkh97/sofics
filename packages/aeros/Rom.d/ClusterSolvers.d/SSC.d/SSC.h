#ifndef SPARSESUBSPACECLUSTERING_H_
#define SPARSESUBSPACECLUSTERING_H_

#include <mpi.h>
#include <cstring>
#include <vector>

#include "Rom.d/DistrNonnegativeMatrixFactorization.h"
#ifdef CLUSTER_DEV
#include "SCDoubleMatrix.h"
#include "SCIntMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCDoubleMatrix.h"
#include "Math.d/SCMatrix.d/SCIntMatrix.h"
#endif


#define SPARSESUBSPACECLUSTERING_HEADER_INCR 30
#define SPARSESUBSPACECLUSTERING_RESIDUAL_NO_DECREASE_REL_TOL 1.0E-14
#define SPARSESUBSPACECLUSTERING_RESIDUAL_NO_DECREASE_COUNT 5

#define TIME_SPARSESUBSPACECLUSTERING_MAIN_LOOP 0
#define TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS 1
#define TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS 2
#define TIME_SPARSESUBSPACECLUSTERING_INIT 3
#define TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE 4
#define TIME_SPARSESUBSPACECLUSTERING_GET_CENTROIDS_SUM 5
#define TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE_COPY 6
#define TIME_SPARSESUBSPACECLUSTERING_TAG_SNAPSHOTS_DISTANCE_ADD 7
#define N_TIMES_SPARSESUBSPACECLUSTERING 8

#define SPARSESUBSPACECLUSTERING_LARGE 1.0e+300
#define SPARSESUBSPACECLUSTERING_MAX_ITER_DEFAULT 1000

class Communicator;

class SparseSubspaceClustering {

    public:
        SparseSubspaceClustering();
        ~SparseSubspaceClustering();

        int setMatrixRow(int j, double *row);
        int setMatrixColumn(int j, double *col);
        int cluster(SCDoubleMatrix& snapshots);
        int close();
        void setCommunicator(Communicator * otherComm) { _otherComm = otherComm; }

        //int getContext() {return _context;};
        double getTime(int i) { return _wallclock_total[i]; }
        double getComputeTime() {return getTime(TIME_SPARSESUBSPACECLUSTERING_MAIN_LOOP);}
        void summary(SCDoubleMatrix& A);
        void printTimes();
        int getStatus() {return _status;}
        int getNumClusters() {return _numClusters;}
        void setNumClusters( double numClusters ) {_numClusters = numClusters;}
        void setSparseTolerance(double tol) { sparseTol = tol; }

        void setMaxIter(int max_iter) {_max_iter = max_iter;}
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
        double _wallclock[N_TIMES_SPARSESUBSPACECLUSTERING];
        double _wallclock_total[N_TIMES_SPARSESUBSPACECLUSTERING];

        double sparseTol;

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
        Communicator * _otherComm;

        int _nprocs;
        int _mypid;
        int _iter;
        int _status;
        unsigned _seed;
        SCDoubleMatrix * _centroids;
        SCDoubleMatrix * _evCentroids;
        SCDoubleMatrix * _eigVectors;
        SCIntMatrix    * _snapshotCentroids; // container for which centroid a snapshot belongs to

        int init(SCDoubleMatrix & A);
        int sparsesubspacecluster(SCDoubleMatrix & snapshots);
        void computeEigenvectors(SCDoubleMatrix & graph);
        void setDefaults();
        void output();
        void header();
        int  sparsesubspaceclusterInit(SCDoubleMatrix & snapshots);
        void tagSnapshot(SCDoubleMatrix & snapshots);
        void getCentroids(SCDoubleMatrix & snapshots, SCDoubleMatrix & _centroidContainer);
        void estimateClusterDimensionality(SCDoubleMatrix &graph);
};

#endif // SPARSESUBSPACECLUSTERING_H_
