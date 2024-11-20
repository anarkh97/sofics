#if defined(USE_MPI)
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <algorithm>

#include "Nmf.h"
#ifdef NNLS_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif


Nmf::Nmf(SCDoubleMatrix& A, SCDoubleMatrix& W, SCDoubleMatrix& Ht) {
    init(A, W, Ht);
}

Nmf::~Nmf() {
    delete _gradf_W;
    delete _gradf_Ht;
    delete _fixedset_W;
    delete _fixedset_Ht;
}


int
Nmf::init(SCDoubleMatrix& A, SCDoubleMatrix& W, SCDoubleMatrix& Ht) {
    // Define MPI Process Grid if not defined
    _FORTRAN(blacs_pinfo)(&_mypid, &_nprocs);

    _A  = &A;
    _W  = &W;
    _Ht = &Ht;

    _context = A.getContext();
    _m = A.getNumberOfRows();
    _n = A.getNumberOfCols();
    _k = W.getNumberOfCols();
    _mb = A.getRowBlockingFactor();
    _nb = A.getColBlockingFactor();

    _FORTRAN(blacs_gridinfo)(&_context, &_mprow, &_npcol, &_myrow, &_mycol);

    // Initialize vectors and matrices
    _gradf_W  = new SCDoubleMatrix(W, false);
    _gradf_Ht = new SCDoubleMatrix(Ht, false);
    //_hessian  = new SCDoubleMatrix(_context, _k, _k, _mb, _nb);
    //_hessianInverse  = new SCDoubleMatrix(*_hessian, false);

    _fixedset_W  = new SCIntMatrix(_context, _m, _k, _mb, _nb, MPI_COMM_WORLD);
    _fixedset_Ht = new SCIntMatrix(_context, _n, _k, _mb, _nb, MPI_COMM_WORLD);

    for (int i=0; i<N_NMF_PQN_TIMES; i++) {
        _wallclock[i] = 0.0;
        _wallclock_total[i] = 0.0;
    }

    // Test
    //SCDoubleMatrix colnorms = SCDoubleMatrix(_context, 1, _n, _mb, _nb, MPI_COMM_WORLD);
    //_A->scaleColumnsByL2Norm(colnorms);

    _alpha = 0.4;
    _tol = 1.0e-10;
    return 0;
}


void
Nmf::summary() {
    if (_mypid == 0) {
        std::cout << "################### Solver Summary ########################" << std::endl;
        std::cout << "Nonnegative Matrix Factorization based on the PQN algorithm." << std::endl;
        std::cout << "   Matrix size: " << std::endl;
        std::cout << "      M = " << _m << std::endl;
        std::cout << "      N = " << _n << std::endl;
        std::cout << "      K = " << _k << std::endl;
        std::cout << "   NMF Parameters:" << std::endl;
        std::cout << "      Maximun outer iterations   = " << _max_iter << std::endl;
        std::cout << "      Number of inner iterations = " << _num_inner_iter << std::endl;
        std::cout << "      Tolerance on delta W       = " << _tol << std::endl;
        std::cout << "   Scalapack blocking factors:" << std::endl;
        std::cout << "      MB = " << _mb << std::endl;
        std::cout << "      NB = " << _mb << std::endl;
        std::cout << "   Line search step size (alpha) = " << _alpha << std::endl;
        std::cout << "Blacs context:"                    << std::endl;
        std::cout << "   mprow = " << _mprow << std::endl;
        std::cout << "   npcol = " << _npcol << std::endl;
        std::cout << "############################################################" << std::endl << std::endl;
    }
}


int
Nmf::writeMatrix(std::string filename, bool compact) {
    _A->write(filename, compact);
    return 0;
}


double
Nmf::getWallTime(){
    struct timeval time;
    gettimeofday(&time,NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


std::string
Nmf::intToString(int i) {
    std::ostringstream ss;
    ss << i;
    return ss.str();
}


int
Nmf::distributeVector(char scope, char top, double *vec, int n) {
    int zero = 0;
    int one = 1;
    std::string col = std::string("C");
    std::string s(1,scope);
    if (col.compare(s) == 0) {
        if (_myrow == 0) {
            _FORTRAN(dgebs2d)(&_context, &scope, &top, &one, &n, vec, &one);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &top, &one, &n, vec, &one, &zero, &_mycol);
        }
    } else {
        if (_mycol == 0) {
            _FORTRAN(dgebs2d)(&_context, &scope, &top, &n, &one, vec, &one);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &top, &n, &one, vec, &n, &_myrow, &zero);
        }
    }
    return 0;
}


int
Nmf::distributeVector(char scope, char top, int *ivec, int n) {
    int zero = 0;
    int one = 1;
    std::string col = std::string("C");
    std::string s(1,scope);
    if (col.compare(s) == 0) {
        if (_myrow == 0) {
            _FORTRAN(igebs2d)(&_context, &scope, &top, &one, &n, ivec, &one);
        } else {
            _FORTRAN(igebr2d)(&_context, &scope, &top, &one, &n, ivec, &one, &zero, &_mycol);
        }
    } else {
        if (_mycol == 0) {
            _FORTRAN(igebs2d)(&_context, &scope, &top, &n, &one, ivec, &one);
        } else {
            _FORTRAN(igebr2d)(&_context, &scope, &top, &n, &one, ivec, &n, &_myrow, &zero);
        }
    }
    return 0;
}


void
Nmf::writeSolution(bool compact) {
    _W->write("W.nmf", compact);
    _Ht->write("Ht.nmf", compact);
}


double
Nmf::conditionNumber() {
    return _A->conditionNumber();
}


double
Nmf::minSingularValue() {
    return _A->minSingularValue();
}


double
Nmf::maxSingularValue() {
    return _A->maxSingularValue();
}


void
Nmf::writeSingularValues(std::string filename) {
    return _A->writeSingularValues(filename);
}


void
Nmf::header() {
    if (_mypid == 0) {
        std::cout << "   ";
        std::cout << "iter   ";
        std::cout << "||A-WH||/||A|| ";
        std::cout << "||W-W_old||_F  ";
        std::cout << std::endl;

        std::cout << "-- ";
        std::cout << "------ ";
        std::cout << "-------------- ";
        std::cout << "-------------- ";
        std::cout << std::endl;
    }
}


void
Nmf::output() {
    if (_mypid == 0) {
        if ((_iter-1)%HEADER_INCR_NMF_PQN == 0) {
            header();
        }   
        printf("   %6d %14.5e %14.5e\n", _iter, _froNorm, _Wincrement);
    }
}


void
Nmf::finalOutput() {
    if (_mypid == 0) {
        std::cout << "Frobenius Norm of A - WH: " << _froNorm << std::endl;
    }
    loadMaxTimes();
    //_A->write("A.nmf");
    //_W->write("W.nmf");
    //_Ht->write("Ht.nmf");
}


void
Nmf::loadMaxTimes() {
    double times_max[N_NMF_PQN_TIMES];
    MPI_Allreduce(_wallclock_total, times_max, N_NMF_PQN_TIMES, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int i=0; i<N_NMF_PQN_TIMES; i++) {
        _wallclock_total[i] = times_max[i];
    }   
}


void
Nmf::printTimes(bool debug) {
    if (debug) {
        if (_mypid == 0) {
            std::cout << std::endl;
            std::cout << "Wallclock Times (seconds):"                                                     << std::endl;
            std::cout << "    Solver            : " << _wallclock_total[TIME_NMF_PQN_MAIN_LOOP]           << std::endl;
            std::cout << "        gradf         : " << _wallclock_total[TIME_NMF_PQN_GRADF]               << std::endl;
            std::cout << "        update        : " << _wallclock_total[TIME_NMF_PQN_UPDATE]              << std::endl;
            std::cout << "        freeset       : " << _wallclock_total[TIME_NMF_PQN_FIXEDSET]            << std::endl;
            std::cout << "        hessian       : " << _wallclock_total[TIME_NMF_PQN_HESSIAN]             << std::endl;
            std::cout << "        residual      : " << _wallclock_total[TIME_NMF_PQN_RESIDUAL]            << std::endl;
        }
    } else {
        if (_mypid == 0) {
            std::cout << std::endl;
            std::cout << "Wallclock Times (seconds):"                                                     << std::endl;
            std::cout << "    Solver            : " << _wallclock_total[TIME_NMF_PQN_MAIN_LOOP]           << std::endl;
            std::cout << "        gradf         : " << _wallclock_total[TIME_NMF_PQN_GRADF]               << std::endl;
            std::cout << "        update        : " << _wallclock_total[TIME_NMF_PQN_UPDATE]              << std::endl;
            std::cout << "        freeset       : " << _wallclock_total[TIME_NMF_PQN_FIXEDSET]            << std::endl;
            std::cout << "        hessian       : " << _wallclock_total[TIME_NMF_PQN_HESSIAN]             << std::endl;
            std::cout << "        residual      : " << _wallclock_total[TIME_NMF_PQN_RESIDUAL]            << std::endl;
        }
    }
}


#endif
