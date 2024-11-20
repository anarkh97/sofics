#ifndef NMF_H_
#define NMF_H_

#include <cstring>
#include <vector>


#ifdef NMF_DEV
#include "SCDoubleMatrix.h"
#include "SCIntMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCDoubleMatrix.h"
#include "Math.d/SCMatrix.d/SCIntMatrix.h"
#endif


#define HEADER_INCR_NMF_PQN 30

#define TIME_NMF_PQN_MAIN_LOOP 0
#define TIME_NMF_PQN_GRADF 1
#define TIME_NMF_PQN_UPDATE 2
#define TIME_NMF_PQN_HESSIAN 3
#define TIME_NMF_PQN_FIXEDSET 4
#define TIME_NMF_PQN_RESIDUAL 5
#define N_NMF_PQN_TIMES 6


class Nmf {

    public:
        Nmf(SCDoubleMatrix& A, SCDoubleMatrix& W, SCDoubleMatrix& H);
        ~Nmf();

        int solve();

        int writeMatrix(std::string filename = std::string("A.txt"), bool compact=false);
        int getContext() {return _context;};
        double getTime(int i) { return _wallclock_total[i]; }
        double getComputeTime() {return getTime(TIME_NMF_PQN_MAIN_LOOP);}
        void setMaxIter(int max_iter) {_max_iter = max_iter;}
        void setAlpha(double alpha) {_alpha = alpha;}
        void setTol(double tol) {_tol = tol;}
        void setNumInnerIter(int num_inner_iter) {_num_inner_iter = num_inner_iter;}
        void summary();
        void printTimes();
        int getStatus() {return _status;}
        void setMaxIterRatio( double maxite ) {_max_iter=maxite*_n;}

        int init(SCDoubleMatrix& A, SCDoubleMatrix& W, SCDoubleMatrix& H);
        void writeSolution(bool compact=false);
        double conditionNumber();
        double minSingularValue();
        double maxSingularValue();
        void writeSingularValues(std::string filename = std::string("svd.txt"));
        void loadMaxTimes();
        void printTimes(bool debug=true);


    private:

        SCDoubleMatrix * _A;
        SCDoubleMatrix * _W;
        SCDoubleMatrix * _Ht;
        SCDoubleMatrix * _gradf_W;
        SCDoubleMatrix * _gradf_Ht;

        SCIntMatrix * _fixedset_W;
        SCIntMatrix * _fixedset_Ht;

        // Timings
        double _wallclock[N_NMF_PQN_TIMES];
        double _wallclock_total[N_NMF_PQN_TIMES];

        double _alpha;
        double _rnorm2;
        double _froNorm;
        double _tol;
        double _Wincrement;

        int _iter;
        int _num_inner_iter;
        int _max_iter;
        int _m;
        int _n;
        int _k;
        int _mb;
        int _nb;
        int _mprow;
        int _npcol;
        int _myrow;
        int _mycol;
        int _context;
        int _nprocs;
        int _mypid;
        int _method;
        int _status;

        int initnmf();
        void fixedSet(SCDoubleMatrix &A, SCDoubleMatrix &gradf, SCIntMatrix &fixedset);
        int update(SCDoubleMatrix& C, SCDoubleMatrix& grad_C, SCIntMatrix& fixed_C, SCDoubleMatrix& B, SCDoubleMatrix& hessian_B);
        void output();
        void header();
        std::string intToString(int i);
        double getWallTime();
        void startTime(int i) {_wallclock[i] = -getWallTime();}
        void stopTime(int i) {_wallclock[i] += getWallTime(); _wallclock_total[i] += _wallclock[i];}
        int distributeVector(char scope, char top, double *vec, int n);
        int distributeVector(char scope, char top, int *ivec, int n);
        void getHessian(SCDoubleMatrix D, SCDoubleMatrix& hessian, SCDoubleMatrix& hessianInverse);
        void gradf(SCDoubleMatrix& A, SCDoubleMatrix& B, SCDoubleMatrix& C, SCDoubleMatrix& hessian, char transpose, SCDoubleMatrix& gradf_B);
        void finalOutput();
        void residual();
};

#endif // NMF_H_
