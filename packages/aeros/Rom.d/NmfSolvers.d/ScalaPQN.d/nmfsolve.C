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

#include "Nmf.h"
#ifdef NNLS_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif


int
Nmf::solve() {
    MPI_Barrier(MPI_COMM_WORLD);
    startTime(TIME_NMF_PQN_MAIN_LOOP);
    SCDoubleMatrix hessian = SCDoubleMatrix(_context, _k, _k, _mb, _nb, MPI_COMM_WORLD);
    SCDoubleMatrix hessianInverse = SCDoubleMatrix(hessian, false);
    double A_froNorm = _A->froNorm();
    bool done = false;
    _iter = 0;
    while (!done) {
        // Make a copy for stopping criteria
        SCDoubleMatrix W_copy(*_W);
        getHessian(*_W, hessian, hessianInverse);
        //hessian.write("hessianW_" + intToString(_iter) + ".txt");
        //hessianInverse.write("hessianWInverse_" + intToString(_iter) + ".txt");
        //std::cout << "Updating H" << std::endl;
        for (int iter=0; iter<_num_inner_iter; iter++) {
            gradf(*_A, *_Ht, *_W, hessian, 'T', *_gradf_Ht);
            fixedSet(*_Ht, *_gradf_Ht, *_fixedset_Ht);
            update(*_Ht, *_gradf_Ht, *_fixedset_Ht, *_W, hessianInverse);
        }
        getHessian(*_Ht, hessian, hessianInverse);
        //hessian.write("hessianHt_" + intToString(_iter) + ".txt");
        //hessianInverse.write("hessianHtInverse_" + intToString(_iter) + ".txt");
        //std::cout << "Updating W" << std::endl;
        for (int iter=0; iter<_num_inner_iter; iter++) {
            gradf(*_A, *_W, *_Ht, hessian, 'N', *_gradf_W);
            fixedSet(*_W, *_gradf_W, *_fixedset_W);
            update(*_W, *_gradf_W, *_fixedset_W, *_Ht, hessianInverse);
        }
        SCDoubleMatrix Err(*_A);
        _W->multiply(*_Ht, Err, 'N', 'T', -1.0, 1.0); // Err = A-W*H;
        _froNorm = Err.froNorm()/A_froNorm;
        _W->add(W_copy, 'N', _m, _k, 1.0, -1.0); // W_copy = W-W_copy
        _Wincrement = W_copy.froNorm();
        output();
        _iter++;
        if (_Wincrement < _tol || _iter >= _max_iter) {
            done = true;
        }
    }
    finalOutput();
    stopTime(TIME_NMF_PQN_MAIN_LOOP);

    return 0;
}


void
Nmf::residual() {
    startTime(TIME_NMF_PQN_RESIDUAL);
    SCDoubleMatrix R = SCDoubleMatrix(*_A);
    _W->multiply(*_Ht, R, 'N', 'T', -1.0, 1.0);
    _froNorm = R.froNorm();
    //_froNorm = _A->froNorm();
    //std::cout << "fro norm is " << _froNorm << std::endl;
    stopTime(TIME_NMF_PQN_RESIDUAL);
}


void
Nmf::getHessian(SCDoubleMatrix D, SCDoubleMatrix& hessian, SCDoubleMatrix& hessianInverse) {
    startTime(TIME_NMF_PQN_HESSIAN);
    D.multiply(D, hessian, 'T', 'N');
    SCDoubleMatrix hessianDecomp  = SCDoubleMatrix(hessian);
    hessianInverse.loadIdentityMatrix();
    int decompStatus = hessianDecomp.choldecomp();
    if (decompStatus != 0) {
        std::cout << "Problem in Nmf::getHessian() on proc " << _mypid << ", iter = " << _iter << ": decompStatus = " << decompStatus << std::endl;
    }
    int solveStatus = hessianDecomp.cholsolve(hessianInverse);
    if (solveStatus != 0) {
        std::cout << "Problem in Nmf::getHessian() on proc " << _mypid << ", iter = " << _iter << ": solveStatus = " << solveStatus << std::endl;
    }
    stopTime(TIME_NMF_PQN_HESSIAN);
}


//multiply(SCDoubleMatrix &B, SCDoubleMatrix &C, char transA, char transB, double alpha, double beta, int m, int n, int k) 
void
Nmf::gradf(SCDoubleMatrix& A, SCDoubleMatrix& B, SCDoubleMatrix& C, SCDoubleMatrix& hessian, char transpose, SCDoubleMatrix& gradf_B) {
    startTime(TIME_NMF_PQN_GRADF);
    A.multiply(C, gradf_B, transpose, 'N');
    B.multiply(hessian, gradf_B, 'N', 'N', 1.0, -1.0);
    stopTime(TIME_NMF_PQN_GRADF);
}


int
Nmf::update(SCDoubleMatrix& C, SCDoubleMatrix& grad_C, SCIntMatrix& fixed_C, SCDoubleMatrix& B, SCDoubleMatrix& hessianInverse) {
    startTime(TIME_NMF_PQN_UPDATE);
    SCDoubleMatrix U = SCDoubleMatrix(grad_C, false);
    grad_C.zeroout(fixed_C);
    grad_C.multiply(hessianInverse, U, 'N', 'N');
    //U.write("U.txt");
    U.zeroout(fixed_C);
    double *matrixC = C.getMatrix();
    double *matrixU = U.getMatrix();
    for (int i=0; i<C.getSizelocal(); i++) {
        double a = matrixC[i] - _alpha * matrixU[i];
        if (a < 0.0) {
            matrixC[i] = 0;
        } else {
            matrixC[i] = a;
        }
    }
    stopTime(TIME_NMF_PQN_UPDATE);
    return 0;
}


void
Nmf::fixedSet(SCDoubleMatrix &A, SCDoubleMatrix &gradf, SCIntMatrix &fixedset) {
    startTime(TIME_NMF_PQN_FIXEDSET);
    if (!(A.isSameShape(gradf) && A.isSameShape(fixedset))) {
        std::cout << "Problem in Nmf::fixedSet. Matrices are not the same size" << std::endl;
        exit(0);
    } 
    double *matrix = A.getMatrix();
    double *grad   = gradf.getMatrix();
    int *fixed     = fixedset.getMatrix();
    for (int j=0; j<A.getSizelocal(); j++) {
        if ((matrix[j] == 0.0) && grad[j] > 0.0) {
            fixed[j] = 1;
        } else {
            fixed[j] = 0;
        }
    }
    stopTime(TIME_NMF_PQN_FIXEDSET);
}

#endif
