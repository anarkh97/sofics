#if defined(USE_MPI) && defined(USE_EIGEN3)
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <cstdio>

#include "Plh.h"

#define X_EQUALS_Z 0
#define X_UPDATED_OLD_QR 1
#define X_UPDATED_NEW_QR 2


int
Plh::solve() {
    MPI_Barrier(_comm);
    initplh();
    _iter = 0;
    _iter_total = 0;
    _subiter = 0;
    bool done = false;
    bool reject;
    int iqr = -1;
    int count;

    startTime(TIME_MAIN_LOOP);
    while (_iter < _max_iter && _nP < _maxNP && _rnorm2 > _rtol && !done) {
        startTime(TIME_ITER);
        gradf();
        done = nextVector(); 
        if (!done) {
            iqr = _nP;
            count = 0;
            _subiter = 0;
            reject = true;  // Force through the loop once.
            while (reject && !done) {
                reject = updateQR(iqr);
                while (reject && !done) { // loop to reject linearly dependent columns
                    done = rejectVector();
                    if (_mypid == 0 && _verbose > 0) {
                         std::cout << "Linear dependent. _wmax.i = " << _wmax.i << ", _wmax.x = " << _wmax.x;
                         std::cout << ", _iter = " << _iter << std::endl;
                    }
                    if (done) break;
                    reject = updateQR(iqr);
                }
                if (done) break;
                if (count == 0) {
                    reject = updateQtb(iqr);
                } else {
                    reject = updateQtb();
                }
                if (reject) {
                    if (_mypid == 0 && _verbose > 0) {
                        std::cout << "Positive z. _wmax.i = " << _wmax.i << ", _iter = " << _iter << std::endl;
                    }
                }
                count++;
            }
            if (done) break;
            solveR();
            _zmin = _zQR->getMin(1, _nP);
            while (_zmin.x <= 0.0 && !_OMP) { // loop to reject elements that violate constraints
                // Downdate
                startTime(TIME_DOWNDATE);
                if(_ddmask) _wmask->setElement(1,_zmin.i,0.0); // Don't allow this vector back until mask is reset.
                updateX();
                //iqr = moveFromPToZSwap();
                iqr = moveFromPToZShift();
                updateQR(iqr);
                updateQtb();
                solveR();
                _zmin = _zQR->getMin(1, _nP);
                 sub_iteration_output(iqr);
                _subiter++;
                _iter_total++;
                stopTime(TIME_DOWNDATE);
            }
            _xQR->zero();
            _zQR->copy(*_xQR, _nP);
            copyxQR2x();
        }
//        if(_PFP) updateVertex(); // this could be done way better that solving Q*R^-T*1 each time
        computeResidual();
        iteration_output();
        if (_iter % _residualIncr == 0) writeResidual();
        stopTime(TIME_ITER);
        _iter++;
        _iter_total++;
    }
    //writeSet();
    copyToHotInd();
    // Write the final residual regardless of _residualIncr 
    if ((_iter-1) % _residualIncr != 0) writeResidual();
    if (_col_scaling) _x->hadamardProduct(*_colnorms);
    loadMaxTimes();
    if (_iter < _max_iter) {
        _status = 1;
    } else {
        _status = 3;
    }
    MPI_Barrier(_comm);
    stopTime(TIME_MAIN_LOOP);
    return _nP;
}


int
Plh::initplh() {
    startTime(TIME_INIT_NNLS);
    _x->zero();
    _xQR->zero();
    _nP = 0;
    _nZ = _n; 
    _QtoA->identityPermutation();
    //_QtoA->zero();

    int lw1 = _Q->getLworkPxormqr(*_Q);
    int lw2 = _Q->getLworkPxgeqrf();
    _lwork_qr = std::max(lw1, lw2);
    if(!_work_qr) _work_qr = new double[_lwork_qr];

    //_b->copy(*_Qtb, _m);
    // _bQR is for convenience
    _b->copyRedist(_m, 1, 1, 1, *_bQR, 1, 1, _contextQR);
    _bQR->copy(*_Qtb, _m);
    _bQR->copy(*_rQR, _m);

    _A->multiply(*_b,  *_Atb, 'T', _m, _n,  1.0, 0.0);

    _rnorm2 = _b->norm2();
    if (_rtol > 0.0) {
        _rtol *= _rnorm2;
    }
    if (_mypid == 0 && _verbose > 0) {
        std::cout << "Initial residual      = " << _rnorm2 << std::endl;
        std::cout << "Residual norm target  = " << _rtol << std::endl;
    }

    if (_ddmask) {
        _wmask->set(1.0);
    }

    if (_col_scaling && (_colnorms == NULL)) {
        startTime(TIME_COLUMNSCALING);
        if (_mypid == 0 && _verbose > 0) {
            std::cout << "Scaling the matrix columns." << std::endl;
        }
        _colnorms = new SCDoubleMatrix(_context, 1, _n, _mb, _nb, _comm);
        _A->scaleColumnsByL2Norm(*_colnorms);           
        stopTime(TIME_COLUMNSCALING);
    }

    if(_PFP) { // initialize working arrays for Polytope Faces Pursuit
      _vertex->zero();
    }

    stopTime(TIME_INIT_NNLS);
    return 0;
}


void
Plh::writeResidual() {
    startTime(TIME_WRITE_RESIDUAL);
    if (_mypid == 0 && _residualFilePtr != NULL) {
        fprintf(_residualFilePtr, "%6d %12.4e %12.4e\n", _iter, _wallclock_total[TIME_ITER], _rnorm2);
    }
    stopTime(TIME_WRITE_RESIDUAL);
}


double
Plh::residual2Norm() {
    startTime(TIME_RESIDUAL2NORM);
    SCDoubleMatrix * res = _workm;
    _b->copy(*res, _m); // Result in _res
    _A->multiply(*_x, *res, 'N', _A->getNumberOfRows(), _A->getNumberOfCols(), 1.0, -1.0);
    _rnorm2 = res->norm2();
    stopTime(TIME_RESIDUAL2NORM);
    return _rnorm2;
}


bool
Plh::updateQR(int iqr) {
    startTime(TIME_UPDATEQR);
    bool reject = false;
    char side = 'L';
    char trans = 'T';
    int zero=0, one=1;
    int info;
    int col_proc;
    int dummy;
    double norm, ajj, diff;
    int ia;
    char scope = 'C';
    char top = ' ';
    double factor = 0.01;
    int ncols = _nP - iqr + 1;
    if (ncols > 0) {
        startTime(TIME_COPYREDIST);
        _QtoA->setScope('A');
        if (_context == _contextQR) {
            for (int iq=iqr; iq<=_nP; iq++) {
                ia = _QtoA->getElement(1,iq);
                _A->pdcopy(_A->getNumberOfRows(), 1, ia, 1, *_Q, 1, iq, 1);
            }
        } else {
            for (int iq=iqr; iq<=_nP; iq++) {
                ia = _QtoA->getElement(1,iq);
                _A->copyRedist(_A->getNumberOfRows(), 1, 1, ia, *_Q, 1, iq, _contextQR);
            }
        }
        _QtoA->setScope();
        stopTime(TIME_COPYREDIST);
        startTime(TIME_PDORMQR);
        int k = iqr-1;
        if (k > 0) {
            _FORTRAN(pdormqr)(&side, &trans, &_m, &ncols, &k,
                _Q->getMatrix(), &one, &one, _Q->getDesc(), _Q->getTau(),
                _Q->getMatrix(), &one, &iqr, _Q->getDesc(),
                _work_qr, &_lwork_qr, &info);
        }
        int nrows = _m - iqr + 1;
        _FORTRAN(pdgeqrf)(&nrows, &ncols, _Q->getMatrix(), &iqr, &iqr, _Q->getDesc(), _Q->getTau(),
             _work_qr, &_lwork_qr, &info);
        stopTime(TIME_PDORMQR);

        // Check for near linear dependence
        double diff;
        if (iqr == _nP && _nP > 1) {
            startTime(TIME_LDCHECK);
            int icol = _FORTRAN(indxg2p)(&_nP, &_nbq, &dummy, &zero, &_npcolQR);
            if (_mycolQR == icol) {
                double norm;
                int m = _nP-1;
                _FORTRAN(pdnrm2)(&m, &norm, _Q->getMatrix(), &one, &_nP, _Q->getDesc(), &one);
                _Q->setScope('C');
                double a = _Q->getElement(_nP, _nP);
                diff = norm + fabs(a)*.01 - norm;
            }
            distributeVector(_contextQR, 'R', top, &diff, 1, icol);
            if (diff <= 0.0) {
                reject = true;
            }
            stopTime(TIME_LDCHECK);
        }
    }
    stopTime(TIME_UPDATEQR);

    return reject;
}

void
Plh::updateVertex() { // this seems more numerically accurate than the fast update, but it is much slower

    // first do forward triangular solve: y = R^-T*1
    {
      startTime(TIME_SOLVER);
      char uplo  = 'U';
      char trans = 'T';
      char diag  = 'N';
      int zero=0, one=1, nrhs = 1, info;
      _oneVecQR->set(1.0); // reinitialize vector to all ones
      _trslvQR->zero();
      startTime(TIME_PDTRSV);
      // overwrite first _nP rows of _neVecQR
      _FORTRAN(pdtrsv)(&uplo, &trans, &diag, &_nP, _Q->getMatrix(), &one, &one, _Q->getDesc(),
                       _oneVecQR->getMatrix(), &one, &one, _oneVecQR->getDesc(), &one); 
      _oneVecQR->add(*_trslvQR,trans,_nP,1,1.0,0.0);
      stopTime(TIME_PDTRSV);
      stopTime(TIME_SOLVER);
    }
     
    // then multiply by Q: vertex = Q*y

    {
      char side  = 'L';
      char trans = 'N';
      int j      = 1; 
      int one    = 1;
      int k      = _nP;
      int m      = _m;
      int info;
      // overwrite first _m rows of _trslvQR
      _FORTRAN(pdormqr)(&side, &trans, &m, &one, &k,
              _Q->getMatrix(),  &j, &j, _Q->getDesc(), _Q->getTau(),
              _trslvQR->getMatrix(), &j, &one, _trslvQR->getDesc(),
              _work_qr, &_lwork_qr, &info);   
      // transfer the result (first _m rows of _trslvQR) into vertex container
      mcopyQtoA(_trslvQR, _trslv);
      _trslv->add(*_vertex,trans,_m,1,1.0,0.0); // now add it over
    }

    return;
}

bool
Plh::updateQtb(int iqr) {
    startTime(TIME_UPDATEQTB);
    bool reject = false;
    char side   = 'L';
    char trans  = 'T';
    char top    = ' ';
    int zero=0, one=1;
    int j, k, m, info;
    int dummy;
    if (iqr == _nP) {
        j = _nP;
        k = 1;
        m = _m - _nP + 1;
    } else {
        _bQR->copy(*_Qtb, _m);
        j = 1;
        k = _nP;
        m = _m;
    }
    _FORTRAN(pdormqr)(&side, &trans, &m, &one, &k,
        _Q->getMatrix(),  &j, &j, _Q->getDesc(), _Q->getTau(),
        _Qtb->getMatrix(), &j, &one, _Qtb->getDesc(),
        _work_qr, &_lwork_qr, &info);

    // Check for z at _nP > 0
    double ztest;
    double a, q;
    if (iqr == _nP && _nP > 1) {
        startTime(TIME_PZCHECK);
        int icol = _FORTRAN(indxg2p)(&_nP, &_nbq, &dummy, &zero, &_npcolQR);
        _Q->setScope('C');
        if (_mycolQR == icol) {
            a = _Q->getElement(_nP, _nP);
        }
        distributeVector(_contextQR, 'R', top, &a, 1, icol);
        if (_mycolQR == 0) {
            q = _Qtb->getElement(_nP, 1);
        }
        distributeVector(_contextQR, 'R', top, &q, 1, 0);
        ztest = q/a;
        if (ztest <= 0.0 && !_OMP) {
            reject = true;
        }
        stopTime(TIME_PZCHECK);
    }

    stopTime(TIME_UPDATEQTB);
    return reject;
}


void
Plh::solveR() {
    startTime(TIME_SOLVER);
    char uplo = 'U';
    char trans = 'N';
    char diag = 'N';
    int d = std::min(_m, _n);
    int zero=0, one=1, nrhs = 1, info;
    _zQR->zero();
    _Qtb->copy(*_zQR, _nP); // copy Q^T*b into working array that will be over-written with solution
    startTime(TIME_PDTRSV);
    _FORTRAN(pdtrsv)(&uplo, &trans, &diag, &_nP, _Q->getMatrix(), &one, &one, _Q->getDesc(),
        _zQR->getMatrix(), &one, &one, _zQR->getDesc(), &one);
    stopTime(TIME_PDTRSV);
    stopTime(TIME_SOLVER);
    return;
}


bool
Plh::nextVector() {
    startTime(TIME_NEXT_VECTOR);
    // For convenience of finding the maximum w in Z, set all w in P to -LARGE
    bool retval = true; // Means we are done. Default is true.
    if(_hotStart) { // force add next element in hotstart set
      int j = hotIndices[_hotInd];
      _hotInd++;
      if(_hotInd >= _sizeHS) _hotStart = false;
      if(j == _constraint && _hotStart) { // theres only one constraint, don't violate it
        j = hotIndices[_hotInd];
        _hotInd++; if(_hotInd >= _sizeHS) _hotStart = false;
      } 
      if(j != _constraint) { // protect against case where last hot indice was constrained
        _nP++;
        _nZ--;
        _wmax.x = _w->getElement(1,j);
        _wmax.i = j; 
        _QtoA->setElement(1,_nP, j);
        retval = false; 
        return retval; 
      }
    } 
    startTime(TIME_GETMAX);
    if(_OMP) _w->elementWiseAbsoluteValue();
    _wmax = _w->getMax(); // getMax() distributes to all processors
    stopTime(TIME_GETMAX);
    if (_wmax.x > 0.0 && _nZ > 0) {
        _nP++;
        _nZ--;
        _QtoA->setElement(1,_nP, _wmax.i);
        retval = false;
    }
    stopTime(TIME_NEXT_VECTOR);
    return retval;
}


bool
Plh::rejectVector() {
    startTime(TIME_REJECT_VECTOR);
    bool done = false;
    if(_hotStart) {
       int j = hotIndices[_hotInd]; _hotInd++; 
       if(_hotInd >= _sizeHS) _hotStart = false; 
       if(j == _constraint && _hotStart) { // theres only one constraint, don't violate it
         j = hotIndices[_hotInd];
         _hotInd++; if(_hotInd >= _sizeHS) _hotStart = false;
       }
       if(j != _constraint) {
         _wmax.x = _w->getElement(1,j);
         _wmax.i = j;     
         _QtoA->setElement(1,_nP,_wmax.i);
         return done;
       }
    } 

    _w->setElement(1,_wmax.i,0.0);
    DoubleInt maxval;
    startTime(TIME_GETMAX);
     _wmax = _w->getMax(); // getMax() distributes to all processors
    stopTime(TIME_GETMAX);
    if (_wmax.x <= 0.0) {
        done = true;
        _nP--;
    } else {
        _QtoA->setElement(1,_nP,_wmax.i);
    }
    
    stopTime(TIME_REJECT_VECTOR);
    return done;
}


int
Plh::mcopyQtoA(SCDoubleMatrix * xQR, SCDoubleMatrix * x) {
    startTime(TIME_MCOPYQTOA);
    x->zero();
    int m = xQR->getNumberOfRows();
    int n = xQR->getNumberOfCols();
    if (_context == _contextQR) {
        int p;
        if (m == 1) {
             p = n;
        } else if (n == 1) {
             p = m;
        } else {
             std::cout << "Problem in mcopyQtoA" << std::endl;
             exit(-1);
        }
        xQR->pdcopy(p, 1, 1, 1, *x, 1, 1, 1);
    } else {
        if (m == 1) {
            xQR->copyRedist(1, _n, 1, 1, *x, 1, 1, _contextQR);
        } else if (n == 1) {
            xQR->copyRedist(_m, 1, 1, 1, *x, 1, 1, _contextQR);
        } else {
            std::cout << "Problem in mcopyQtoA()";
            exit(-1);
        }
    }
    stopTime(TIME_MCOPYQTOA);
    return 0;
}


int
Plh::copyxQR2x() {
    // Copy _xQR to _x
    startTime(TIME_COPYXQRTOX);
    _x->zero();
    if (_context == _contextQR) {
        _xQR->pdcopy(_nP, 1, 1, 1, *_x, 1, 1, 1);
    } else {
        _xQR->copyRedist(1, _nP, 1, 1, *_x, 1, 1, _contextQR);
    }
    //_xQR->write("xQR_" + intToString(_iter) + ".txt");
    //_x->write("x_before_reorder" + intToString(_iter) + ".txt");
    _x->reorder(*_QtoA, _nP);
    //_x->write("x_after_reorder_" + intToString(_iter) + ".txt");
    stopTime(TIME_COPYXQRTOX);
    return 0;
}


void
Plh::computeResidual() {
    if(_PFP){
      double lambda = 1.0/_wmax.x;
      mcopyQtoA(_rQR, _workm); // copy to same context
      _workm->add(*_vertex, 'N', _m, 1, lambda, 1.0); // vertex^(k+1/2) = vertex^(k-1) + lambda^(k)*resdidual^(k-1);
    }
    _Qtb->copy(*_rQR, _m);
    for (int i=1; i<=_nP; i++) {
        _rQR->setElement(i,1,0.0);
    }
    char side = 'L';
    char trans = 'N';
    int one = 1;
    int info;
    if (_nP > 0) {
        startTime(TIME_PDORMQR_GRADF);
        _FORTRAN(pdormqr)(&side, &trans, &_m, &one, &_nP,
            _Q->getMatrix(), &one, &one, _Q->getDesc(), _Q->getTau(),
            _rQR->getMatrix(), &one, &one, _rQR->getDesc(),
            _work_qr, &_lwork_qr, &info);
        stopTime(TIME_PDORMQR_GRADF);
    }
    if(_PFP){
      double lambda = -1.0/_wmax.x;
      mcopyQtoA(_rQR, _workm);
      _workm->add(*_vertex, 'N', _m, 1, lambda, 1.0); // vertex^(k+1/2) = -lambda^(k)*resdidual^(k-1);
    }
    _rnorm2 = _rQR->norm2();
}


int
Plh::gradf() {
    startTime(TIME_GRADF);
    int j, k;
    double x;

    mcopyQtoA(_rQR, _workm);                                   // dump contents of _rQR to _workm
    MPI_Barrier(_comm);
    startTime(TIME_MULT_GRADF);
    if(_hotStart && _hotInd < _sizeHS-1) {
      _w->zero();
      if(_PFP) _Atv->set(1.0);
      j = hotIndices[ _hotInd];
      _A->multiply('T', _m, 1, 1.0, 1, j,
                   *_workm, 1, 1, 1, 0.0,
                   *_w, 1, j, 1);
      if(_PFP) {
        _A->multiply('T', _m, 1, 1.0, 1, j, 
                     *_vertex, 1, 1, 1,-1.0,
                     *_Atv, 1, j, 1);
      }
      if(_PFP) _w->invHadamardProduct(*_Atv);
    } else {
      _A->multiply(*_workm, *_w, 'T', _m, _n, 1.0, 0.0);       // this permorfs A^T*residual
      if(_PFP) {                                               // extra step for Polytope faces pursuit
        _Atv->set(1.0);                                        // compute 1 - A^T*_vertex
        _A->multiply(*_vertex, *_Atv, 'T', _m, _n, -1.0, 1.0);
        _w->invHadamardProduct(*_Atv);                         // compute A^T*residual/(1 - A^T*_vertex) | A^T&residual > 0
      }
    }
    stopTime(TIME_MULT_GRADF);
    for (int i=1; i<=_nP; i++) {                             // zero out the elements already in the active set
        j = _QtoA->getElement(1,i);
        _w->setElement(1,j,0.0);
    }
    if (_ddmask) _w->hadamardProduct(*_wmask);               // if downdating, mask stuff
    if (_constraint > -1) _w->setElement(1,_constraint,0.0); // enforce an element not to be selected
    stopTime(TIME_GRADF);
    return 0;
}


int
Plh::updateX() {
    startTime(TIME_UPDATEX);
    getAlpha();
    computeX();
    _xQR->setElement(1, _ialpha, 0.0);
    stopTime(TIME_UPDATEX);
    return 0;
}


int
Plh::computeX() {
    double *x = _xQR->getMatrix();
    double *z = _zQR->getMatrix();
    int n = _xQR->getNumberOfColsLocal();
    int zero=0;
    if (_alpha > 0.0) {
        for (int i=0; i<n; i++) {
            x[i] += _alpha * (z[i] - x[i]);
            if (x[i] < 0.0) {
                x[i] = 0.0;  // Only < 0 due to roundoff
            }
        }
    }
    return 0;
}


int
Plh::getAlpha() {
    _alpha  = 0.0;
    _ialpha = -1;
    _xQR->setScope();
    _zQR->setScope();
    if (_myrowQR == 0) {
        int zero = 0, one = 1, p;
        int iz;
        int iproc = 0;
        double x, z, a;
        int ialpha = -1;
        double alpha = 2.0;
        int count = 0, ialoc;
        for (int ia=1; ia<=_nP; ia++) {
            z = _zQR->getElement(1, ia);
            if (z <= 0.0) {
                x = _xQR->getElement(1, ia);
                a = x/(x-z);
                if (alpha > a) {
                     alpha = a;
                     ialpha = ia;
                }
            }
        }
        
        // Broadcast local alpha
        double * minalpha  = new double[_npcolQR];
        int * iminalpha = new int[_npcolQR];
        char top = ' ';
        char scope = 'R';
        for (int pc=0; pc<_npcolQR; pc++) {
            if (_mycolQR == pc) {
                minalpha[pc]  = alpha;
                iminalpha[pc] = ialpha;
                _FORTRAN(dgebs2d)(&_contextQR, &scope, &top, &one, &one, &(minalpha[pc]),  &one);
                _FORTRAN(igebs2d)(&_contextQR, &scope, &top, &one, &one, &(iminalpha[pc]), &one);
            } else {
                _FORTRAN(dgebr2d)(&_contextQR, &scope, &top, &one, &one, &(minalpha[pc]),  &one, &zero, &pc);
                _FORTRAN(igebr2d)(&_contextQR, &scope, &top, &one, &one, &(iminalpha[pc]), &one, &zero, &pc);
            }
        }
        _alpha  = minalpha[0];
        _ialpha = iminalpha[0];
        for (int pc=1; pc<_npcolQR; pc++) {
            if (iminalpha[pc] != -1) {
                if (_ialpha == -1) {
                    _alpha = minalpha[pc];
                    _ialpha = iminalpha[pc];
                } else if (minalpha[pc] < _alpha) {
                    _alpha = minalpha[pc];
                    _ialpha = iminalpha[pc];
                }
            }
        }
        delete[] minalpha;
        delete[] iminalpha;
    }
    distributeVector(_contextQR, 'C', ' ', &_ialpha, 1);
    distributeVector(_contextQR, 'C', ' ', &_alpha, 1);
    _xQR->setElement(1, _ialpha, 0.0);
    return 0;
}


// contextQR routine
int
Plh::moveFromPToZSwap() {
    startTime(TIME_MOVEFROMPTOZ);
    int qmin = _nP;
    int zero= 0;
    char scope = 'R';
    char top = ' ';
    int ra, ca, ldia=-1, rdest=-1, cdest;
    int nP = _nP;
    int k;
    bool found;

    _xQR->setScope('A');
    _QtoA->setScope('A');
    //if (_myrowQR == 0) {
        double x, y;
        int i = 1;
        int j = _nP;
        while (i <= j) {
            x = _xQR->getElement(1,i);
            if (x <= 0.0) {
                if (i < qmin) {
                    qmin = i;
                }
                found = false;
                while (j >= i && !found) {
                    y = _xQR->getElement(1,j);
                    if (y > 0.0) {
                        found = true;
                        k = _QtoA->getElement(1,j);
                        _QtoA->setElement(1,i,k);
                        _xQR->setElement(1,i,y);
                        _nP--;
                        _nZ++;
                    } else {
                        _nP--;
                        _nZ++;
                    }
                    j--;
                }
            }
            i++;
        }
        //Cigamn2d(_contextQR, &scope, &top, 1, 1, (char *) &qmin, 1, &ra, &ca, ldia, rdest, cdest);
    //}
    _xQR->setScope();
    _QtoA->setScope();
    stopTime(TIME_MOVEFROMPTOZ);
    //if (qmin > _nP) {
    //    std::cout << "qmin > _nP. qmin = " << qmin << ", _nP = " << _nP << std::endl;
    //}

    // std::cout << "qmin = " << qmin << ", _nP = " << _nP << std::endl;
    return qmin;
}


// contextQR routine
int
Plh::moveFromPToZShift() {
    startTime(TIME_MOVEFROMPTOZ);
    int qmin = _nP;
    //int zero= 0;
    //char scope = 'R';
    //char top = ' ';
    int nP = _nP;
    bool found;

    if (_contextQR != _context) {
        _xQR->setScope('A');
        _QtoA->setScope('A');
        double x, y;
        int k, i=1, n=0;
        while (i <= _nP) {
            x = _xQR->getElement(1,i);
            if (x <= 0.0) {
                if (i < qmin) {
                    qmin = i;
                }
                n++;
            } else {
                if (n > 0) {
                    k = _QtoA->getElement(1,i);
                    _QtoA->setElement(1,i-n,k);
                    _xQR->setElement(1,i-n,x);
                }
            }
            i++;
        }
        _nP -= n;
        _nZ += n;
        _xQR->setScope();
        _QtoA->setScope();
    } else {
        _xQR->setScope();
        _QtoA->setScope();
        if (_myrowQR == 0) {
            double x, y;
            int k, i=1, n=0;
            while (i <= _nP) {
                x = _xQR->getElement(1,i);
                if (x <= 0.0) {
                    if (i < qmin) {
                        qmin = i;
                    }
                    n++;
                } else {
                    if (n > 0) {
                        k = _QtoA->getElement(1,i);
                        _QtoA->setElement(1,i-n,k);
                        _xQR->setElement(1,i-n,x);
                    }
                }
                i++;
            }
            _nP -= n;
            _nZ += n;
        }
        int buf[3];
        buf[0] = qmin;
        buf[1] = _nP;
        buf[2] = _nZ;
        distributeVector(_contextQR, 'C', ' ', buf, 3);
        if (_myrowQR != 0) {
            qmin = buf[0];
            _nP  = buf[1];
            _nZ  = buf[2];
        }
    }
    stopTime(TIME_MOVEFROMPTOZ);
    return qmin;
}


void
Plh::writeSet(std::string filename) {
    FILE * fptr = fopen(filename.c_str(), "w");
    int j;
    for (int i=1; i<= _nP; i++) {
        j = _QtoA->getElement(1,i);
        fprintf(fptr, "%d\n", j);
    }
    fclose(fptr);
}

void
Plh::writeSet() {
    int j;
    fprintf(stderr," ... writing _QtoA ... \n");
    for (int i=1; i<= _nP; i++) {
        j = _QtoA->getElement(1,i);
        fprintf(stderr, "%d ", j);
    }
    fprintf(stderr, "\n");
}

void
Plh::writeHotSet() {
    fprintf(stderr," ... writing hotIndices ... \n");
    for (std::vector<int>::iterator it = hotIndices.begin();
         it != hotIndices.end(); it++) {
        fprintf(stderr, "%d ", *it);
    }
    fprintf(stderr, "\n");
}

void
Plh::copyToHotInd(){
    int j;
    hotIndices.clear();
    for(int i=1; i<= _nP; i++){
      j = _QtoA->getElement(1,i,'A');
      hotIndices.push_back(j);
    }
    //fprintf(stderr,"size of set is %d \n", hotIndices.size());
}

#endif
