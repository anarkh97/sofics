#ifdef USE_MPI
#ifdef SCARRAYS_DEV
#include "SCBaseMatrix.h"
#include "scpblas.h"
#include "scblacs.h"
#else
#include "Math.d/SCMatrix.d/SCBaseMatrix.h"
#include "Math.d/SCMatrix.d/scpblas.h"
#include "Math.d/SCMatrix.d/scblacs.h"
#endif

#include <iostream>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <sys/time.h>

int SCBaseMatrix::ZERO=0;

SCBaseMatrix::SCBaseMatrix(std::string filename, int context, int mb, int nb, MPI_Comm comm) :
    _context(context), _mb(mb), _nb(nb), _comm(comm) {
    getMatrixSize(filename);
    init();
}


SCBaseMatrix::SCBaseMatrix(int context, int m, int n, int mb, int nb, MPI_Comm comm) :
    _context(context), _m(m), _n(n), _mb(mb), _nb(nb), _comm(comm) {
    init();
}


SCBaseMatrix::~SCBaseMatrix() {
    if (_row_col_comm_set) {
        MPI_Comm_free(&_row_comm);
        MPI_Comm_free(&_col_comm);
    }
}


int
SCBaseMatrix::getMatrixSize(std::string filename) {
    int count;
    FILE * fptr = fopen(filename.c_str(), "rb");
    count = 0;
    count += fread(&_m, sizeof(int), 1, fptr);
    count += fread(&_n, sizeof(int), 1, fptr);
    fclose(fptr);
    return count;
}


void
SCBaseMatrix::init() {
    int zero = 0;
    _FORTRAN(blacs_pinfo)(&_mypid, &_nprocs);
    _FORTRAN(blacs_gridinfo)(&_context, &_mprow, &_npcol, &_myrow, &_mycol);
    _mlocal = _FORTRAN(numroc)(&_m, &_mb, &_myrow, &zero, &_mprow);
    _nlocal = _FORTRAN(numroc)(&_n, &_nb, &_mycol, &zero, &_npcol);
    _lld = std::max(1, _mlocal);
    this->setTopology();
    this->setScope();
    this->setSrc();
    int info;
    _FORTRAN(descinit)(_desc, &_m, &_n, &_mb, &_nb, &_rsrc, &_csrc, &_context, &_lld, &info);
    _sizelocal = std::max(_mlocal,1)*std::max(_nlocal,1);
    _row_col_comm_set = false;
}


void
SCBaseMatrix::printDescription() {
    std::cout << "Matrix description:" << std::endl;
    for (int i=0; i < DLEN_; i++) {
        std::cout << "   mypid = " << _mypid << ", _desc[" << i << "] = " << _desc[i] << std::endl;
    }
}


void
SCBaseMatrix::printSummary() {
    std::cout << "mypid = " << _mypid << ", _mprow  = " << _mprow << ", _npcol   = " << _npcol  << std::endl;
    std::cout << "mypid = " << _mypid << ", _myrow  = " << _myrow << ", _mycol   = " << _mycol  << std::endl;
    std::cout << "mypid = " << _mypid << ", _mlocal = " << _mlocal << ", _nlocal = " << _nlocal << std::endl;
}


// For distributing quantities other than the class matrix
int
SCBaseMatrix::distributeVector(int *vec, int n) {
    int zero = 0;
    int one = 1;
    if (_m == 1) {
        char scope = 'C';
        if (_myrow == 0) {
            _FORTRAN(igebs2d)(&_context, &scope, &_top, &one, &n, vec, &one);
        } else {
            _FORTRAN(igebr2d)(&_context, &scope, &_top, &one, &n, vec, &one, &zero, &_mycol);
        }
    } else {
        char scope = 'R';
        if (_mycol == 0) {
            _FORTRAN(igebs2d)(&_context, &scope, &_top, &n, &one, vec, &one);
        } else {
            _FORTRAN(igebr2d)(&_context, &scope, &_top, &n, &one, vec, &n, &_myrow, &zero);
        }
    }
    return 0;
}


int
SCBaseMatrix::distributeVector(double *vec, int n) {
    int zero = 0;
    int one = 1;
    if (_m == 1) {
        char scope = 'C';
        if (_myrow == 0) {
            _FORTRAN(dgebs2d)(&_context, &scope, &_top, &one, &n, vec, &one);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &_top, &one, &n, vec, &one, &zero, &_mycol);
        }
    } else {
        char scope = 'R';
        if (_mycol == 0) {
            _FORTRAN(dgebs2d)(&_context, &scope, &_top, &n, &one, vec, &one);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &_top, &n, &one, vec, &n, &_myrow, &zero);
        }
    }
    return 0;
}


void
SCBaseMatrix::setScope(char scope) {
    _scope = scope;
}


void
SCBaseMatrix::setScope() {
    if (_n == 1) {
        _scope = 'C';
    } else if (_m == 1) {
        _scope = 'R';
    } else {
        _scope = 'A';
   }
}



// From Scalapack Documentation:
//
// IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
// IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
// NpA0 = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
// 
// IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
// ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
// ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
// MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
// NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
//
// SIDE = 'L'
// LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) + NB_A * NB_A
//
//  c1 = (NB_A*(NB_A-1))/2
//  c2 = (NqC0 + MpC0)*NB_A 

// SIDE = 'R'
// LWORK >= MAX( (NB_A*(NB_A-1))/2,
//              ( NqC0 + MAX( NpA0 + NUMROC( NUMROC( N+ICOFFC, NB_A, 0, 0, NPCOL ), NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A
//            ) + NB_A * NB_A
//
//  c3 =  ( NqC0 + MAX( NpA0 + NUMROC( NUMROC( N+ICOFFC, NB_A, 0, 0, NPCOL ), NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A
//     =  ( NqC0 + MAX(NpA0 + d2, MpC0))*NB_A
//
//  d1 = NUMROC( N+ICOFFC, NB_A, 0, 0, NPCOL )
//  d2 = NUMROC( d1, NB_A, 0, 0, LCMQ )
//
int
SCBaseMatrix::getLworkPxormqr(SCBaseMatrix& C, int ia, int ja, int rsrca, int csrca, int ic, int jc, int rsrcc, int csrcc) {
    int n, m;
    int zero = 0;
    int mbC = C._mb;
    int nbC = C._nb;

    int iroffa = (ia-1) % _mb;
    int icoffa = (ja-1) % _nb;
    int iarow = _FORTRAN(indxg2p)(&ia, &_mb, &_myrow, &rsrca, &_mprow);
    n = _n + iroffa;
    int NpA0 = _FORTRAN(numroc)(&n, &_mb, &_myrow, &iarow, &_mprow);

    int iroffc = (ic-1) % mbC;
    int icoffc = (jc-1) % nbC;
    int icrow = _FORTRAN(indxg2p)(&ic, &mbC, &_myrow, &rsrcc, &_mprow);
    int iccol = _FORTRAN(indxg2p)(&jc, &nbC, &_mycol, &csrcc, &_npcol);
    m = _m + iroffc;
    n = _n + icoffc;
    int MpC0 = _FORTRAN(numroc)(&m, &mbC, &_myrow, &icrow, &_mprow);
    int NqC0 = _FORTRAN(numroc)(&n, &nbC, &_mycol, &iccol, &_npcol);

    int c1  = (_nb*(_nb-1))/2;
    int c2 = (NqC0 + MpC0)*_nb;
    int lworkL = std::max(c1, c2) + _nb*_nb;

    int lcmq = _FORTRAN(ilcm)(&_mprow, &_npcol) / _npcol;
    n = _n + icoffc;
    int d1 = _FORTRAN(numroc)(&n,  &_nb, &zero, &zero, &_npcol);
    int d2 = _FORTRAN(numroc)(&d1, &_nb, &zero, &zero, &lcmq);
    int c3 = (NqC0 + std::max(NpA0 + d2, MpC0))*_nb;
    int lworkR = std::max(c1, c3) + _nb*_nb;

    int lwork = std::max(lworkL, lworkR);

    return lwork;
}


int
SCBaseMatrix::getLworkPxgeqrf(int ia, int ja, int rsrc, int csrc) {
    int zero = 0;
    int iroff = (ia-1) % _mb;
    int icoff = (ja-1) % _nb;
    int iarow = _FORTRAN(indxg2p)(&ia, &_mb, &_myrow, &rsrc, &_mprow);
    int iacol = _FORTRAN(indxg2p)(&ja, &_nb, &_mycol, &csrc, &_npcol);
    int m = _m + iroff;
    int n = _n + icoff;
    int Mp0 = _FORTRAN(numroc)(&m, &_mb, &_myrow, &iarow, &_mprow);
    int Nq0 = _FORTRAN(numroc)(&n, &_nb, &_mycol, &iacol, &_npcol);
    int lwork = _nb * (Mp0 + Nq0 + _nb);
    return lwork;
}


char
SCBaseMatrix::getPivroc() {
    char pivroc;
    if (_m == 1) {
        pivroc = 'R';
    } else if (_n == 1) {
        pivroc = 'C';
    } else {
        pivroc = 'N';
    }
    return pivroc;
}


// Problem with _mprow == _npcol case
int
SCBaseMatrix::getLworkPdlapiv(char rowcol) {
    char pivroc = this->getPivroc();
    int ldw = 0;
    int zero = 0;
    if (rowcol == pivroc) {
        if (rowcol == 'C') {
            int m = _m + _m; // Should be _m + _mb but doesn't work. This maximizes it.
            int locc = _FORTRAN(numroc)(&m, &_nb, &_mycol, &zero, &_npcol);
            if (_mprow == _npcol) {
                ldw = locc + _mb;
            } else {
                int lcm = _FORTRAN(ilcm)(&_mprow, &_npcol);
                int locr = _FORTRAN(numroc)(&_m, &_mb, &_myrow, &zero, &_mprow);
                double c1 = ceil((double) locr / ((double) (_mb)));
                double c2 = (double) lcm / ((double) _mprow);
                int k = (int) ceil (c1/c2);
                ldw = locc + _mb * k;
            }
        } else if (rowcol == 'R') {
            int n = _n + _n;
            int locr = _FORTRAN(numroc)(&n, &_mb, &_myrow, &zero, &_mprow);
            if (_mprow == _npcol) {
                ldw = locr + _n;
            } else {
                int lcm = _FORTRAN(ilcm)(&_mprow, &_npcol);
                int locc = _FORTRAN(numroc)(&_n, &_nb, &_mycol, &zero, &_npcol);
                double c1 = ceil(((double) locc) / ((double) (_nb)));
                double c2 = ((double) lcm) / ((double) _npcol);
                int k = (int) ceil (c1/c2);
                ldw = locr + _nb * k;
            }
        }
    }
    return ldw;
}


double
SCBaseMatrix::getWallTime(){
    struct timeval time;
    gettimeofday(&time,NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


int
SCBaseMatrix::getProc(int ig, int jg) {
    int irow = _FORTRAN(indxg2p)(&ig, &_mb, &_myrow, &_rsrc, &_mprow);
    int jcol = _FORTRAN(indxg2p)(&jg, &_nb, &_mycol, &_csrc, &_npcol);
    int p = _FORTRAN(blacs_pnum)(&_context, &irow, &jcol);
    return p;
}


int
SCBaseMatrix::getLocalOffset(int ig, int jg) {
    int dummy=0;
    int iloc = _FORTRAN(indxg2l)(&ig, &_mb, &dummy, &dummy, &_mprow);
    int jloc = _FORTRAN(indxg2l)(&jg, &_nb, &dummy, &dummy, &_npcol);
    int offset = _mlocal*iloc + _nlocal*jloc;
    return offset;
}


void
SCBaseMatrix::setRowColComms() {
    MPI_Comm_split(_comm, _myrow, _mycol, &_row_comm);
    MPI_Comm_split(_comm, _mycol, _myrow, &_col_comm);
    _row_col_comm_set = true;
}
#endif
