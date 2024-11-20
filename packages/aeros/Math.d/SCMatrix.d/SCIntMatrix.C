#ifdef USE_MPI
#ifdef SCARRAYS_DEV
#include "SCIntMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCIntMatrix.h"
#endif

#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <algorithm>


SCIntMatrix::SCIntMatrix(int context, int m, int n, int mb, int nb, MPI_Comm comm, bool pvec) :
    _pvec(pvec), SCBaseMatrix(context, m, n, mb, nb, comm) {

    init();
}


SCIntMatrix::SCIntMatrix(const SCIntMatrix& matrix) : 
    _pvec(matrix._pvec), SCBaseMatrix(matrix._context, matrix._m, matrix._n, matrix._mb, matrix._nb, matrix._comm) {
    SCBaseMatrix::init();
    SCIntMatrix::init();
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] = matrix._matrix[i];
    }
}


// This assmues the block size for this array is the same as the matrix it will permute.
void
SCIntMatrix::init() {
    int dim = _sizelocal;
    if (_pvec) {
        int zero=0;
        int m, n, mr, mc, nr, nc;

        m = 2*_m - 1;
        mr = _FORTRAN(numroc)(&m, &_mb, &_myrow, &zero, &_mprow) + _mb;
        m = _n + _mb;
        mc = _FORTRAN(numroc)(&m, &_mb, &_myrow, &zero, &_mprow);
        m = std::max(mr, mc);

        n = 2*_n - 1;
        nr = _FORTRAN(numroc)(&n, &_nb, &_mycol, &zero, &_npcol) + _nb;
        n = _m + _nb;
        nc = _FORTRAN(numroc)(&n, &_nb, &_mycol, &zero, &_npcol);
        n = std::max(nr, nc);

        dim = std::max(n, m);
        if (dim < _sizelocal) {
            std::cout << "Problem in SCIntMatrix::init()." << std::endl;
            exit(-1);
        }
    }
    _matrix = new int[dim];
}


SCIntMatrix::~SCIntMatrix() {
    delete[] _matrix;
}

void
SCIntMatrix::writeLocal(std::string filename) {
    int value;
    FILE *f;
    f = fopen(filename.c_str(), "w");
    int m = std::max(1, _mlocal);
    int n = std::max(1, _nlocal);
    // std::cout << "m = " << m << ", n = " << n << std::endl;
    int ig, jg, il, jl, k;
    int zero = 0;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            il = i + 1;
            jl = j + 1;
            ig = _FORTRAN(indxl2g)(&il, &_mb, &_myrow, &zero, &_mprow);
            jg = _FORTRAN(indxl2g)(&jl, &_nb, &_mycol, &zero, &_npcol);
            k = i + j*m;
            value = _matrix[k];
            fprintf(f, "%d %d (%d,%d)\n", k, value, ig, jg);
        }
    }
    fclose(f);
}


// For gdb
void
SCIntMatrix::write(const char * fname) {
    write(std::string(fname), _m, _n);
}


void
SCIntMatrix::write(std::string filename, int m, int n) {
    char scope = 'A';
    char blank = ' ';
    int alpha;
    int one = 1;
    FILE *f;
    if (m == 0) m = _m;
    if (n == 0) n = _n;
    if (_mypid == 0) {
        f = fopen(filename.c_str(), "w");
    }
    for (int i=1; i<=m; i++) {
        for (int j=1; j<=n; j++) {
            _FORTRAN(pielget)(&scope, &blank, &alpha, _matrix, &i, &j, _desc);
            if (_mypid == 0) {
                fprintf(f, "%d ", alpha);
                if (_m == 1 || _n == 1) { 
                    fprintf(f, "\n");
                }    
            }
        }
    }
    if (_mypid == 0) {
        fclose(f);
    }
}

// j is an element of [1,_n]; Starts at 1 for Fortran
int
SCIntMatrix::setMatrixColumn(int j, int *col) {
    if (j > 0 && j <= _n) {
        for (int i=1; i<=_m; i++) {
            _FORTRAN(pielset)(_matrix, &i, &j, _desc, &(col[i-1]));
        }   
    } else {
        std::cerr << "Problem in SCIntMatrix::setMatrixColumn. j = " << j << " is not a valid column index." << std::endl;
    }
    return 0;
}


// Note, only for vectors. Not matrices!
int
SCIntMatrix::distributeVector() {
    int zero = 0;
    int one = 1;
    if (_m == 1) {
        char scope = 'C';
        if (_myrow == 0) {
            _FORTRAN(igebs2d)(&_context, &scope, &_top, &one, &_nlocal, _matrix, &one);
        } else {
            _FORTRAN(igebr2d)(&_context, &scope, &_top, &one, &_nlocal, _matrix, &one, &zero, &_mycol);
        }
    } else if (_n == 1) {
        char scope = 'R';
        if (_mycol == 0) {
            _FORTRAN(igebs2d)(&_context, &scope, &_top, &_mlocal, &one, _matrix, &one);
        } else {
            _FORTRAN(igebr2d)(&_context, &scope, &_top, &_mlocal, &one, _matrix, &_mlocal, &_myrow, &zero);
        }
    } else {
        return 1;
    }
    return 0;
}


int
SCIntMatrix::zero() {
    std::memset(_matrix, 0, _sizelocal*sizeof(int));
    return 0;
}


int
SCIntMatrix::identityPermutation() {
    // Column permutation must be distributed over all process rowws
    // Construct on process row 0
    int zero = 0;
    if (_m == 1) {
        for (int i=1; i <=_nlocal; i++) {
            _matrix[i-1] = _FORTRAN(indxl2g)(&i, &_nb, &_mycol, &zero, &_npcol);
        }
    } else {
        for (int i=1; i <=_mlocal; i++) {
            _matrix[i-1] = _FORTRAN(indxl2g)(&i, &_mb, &_myrow, &zero, &_mprow);
        }
    }
    return 0;
}


void
SCIntMatrix::setElement(int i, int j, int value) {
    _FORTRAN(pielset)(_matrix, &i, &j, _desc, &value);
}


int
SCIntMatrix::getElement(int i, int j) {
    int value;
    _FORTRAN(pielget)(&_scope, &_top, &value, _matrix, &i, &j, _desc);
    return value;
}

int
SCIntMatrix::getElement(int i, int j, char scope) {
    int value;
    _FORTRAN(pielget)(&scope, &_top, &value, _matrix, &i, &j, _desc);
    return value;
}

int
SCIntMatrix::getLocalElements(int *elems) {
    for (int i=0; i<_sizelocal; i++) {
        elems[i] = _matrix[i];
    }
    return _sizelocal;
}


// New permute
int
SCIntMatrix::permute(char direc, char rowcol, SCIntMatrix &ip, int m, int n) {
    int zero = 0, one = 1;

    if (m == 0) m = _m;
    if (n == 0) n = _n;

    int *iwork = NULL;  // iwork only needed if ipiv is needs to be transposed.
    int lwork = ip.getLworkPdlapiv(rowcol);
    if (lwork > 0) {
        iwork = (int *) malloc(lwork * sizeof(int));
    }
    char pivroc = ip.getPivroc();
    float * matrix = (float *) _matrix;
    _FORTRAN(fpslapiv)(&direc, &rowcol, &pivroc, &m, &n, matrix, &one, &one, _desc,
              ip.getMatrix(), &one, &one, ip.getDesc(), iwork);
    if (iwork != NULL) {
        free(iwork);
    }
    return 0;
}


// Only for vectors.
int
SCIntMatrix::reorder(SCIntMatrix& order) {
    if (_m == 1) {
        if (_myrow == 0) {
            SCIntMatrix *mtmp = new SCIntMatrix(*this);
            int value;
            int k;
            mtmp->setScope();
            for (int j=1; j<=_n; j++) {
                k = order.getElement(1,j);
                if (k != j) {
                    value = mtmp->getElement(1,j);
                    this->setElement(1,k,value);
                }
            }
            delete mtmp;
        }
    } else if (_n == 1) {
        if (_mycol == 0) {
            SCIntMatrix *mtmp = new SCIntMatrix(*this);
            int value;
            int k;
            mtmp->setScope();
            for (int i=1; i<=_m; i++) {
                k = order.getElement(i,1);
                if (k != i) {
                    value = mtmp->getElement(i,1);
                    this->setElement(k,1,value);
                }
            }
        delete mtmp;
        }
    }
    return 0;
}


void
SCIntMatrix::swap(int i, int j) {
    int tmp;
    if (_m == 1) {
        if (_myrow == 0) {
            tmp = this->getElement(1,i);
            this->setElement(1, i, this->getElement(1,j));
            this->setElement(1, j, tmp);
        }
    } else if (_n == 1) {
        if (_mycol == 0) {
            tmp = this->getElement(i,1);
            this->setElement(i, 1, this->getElement(j,1));
            this->setElement(j, 1, tmp);
        }
    } else {
        std::cout << "SCIntMatrix::swap(i,j) only handles vectors" << std::endl;
    }
}
 

int
SCIntMatrix::isEqual(SCIntMatrix& imat) {
    int * matrix = imat.getMatrix();
    int count = 0;
    for (int i=0; i<_sizelocal; i++) {
        if (_matrix[i] != matrix[i]) {
            count++;
        }
    }
    int minusone=-1, zero=0, one=1;
    _FORTRAN(igsum2d)(&_context, &_scope, &_top, &one, &one, &count, &one, &minusone, &zero);
    SCBaseMatrix::distributeVector(&count, 1);
    return count;
}


int
SCIntMatrix::countValue(int value) {
    int localCount = 0;
    for (int i=0; i<_sizelocal; i++) {
        if (_matrix[i] == value) localCount++;
    }
    int count = localCount;
    int minusone=-1, zero=0, one=1;
    _FORTRAN(igsum2d)(&_context, &_scope, &_top, &one, &one, &count, &one, &minusone, &zero);
    //SCBaseMatrix::distributeVector(&count, 1);
    return count;
}


int
SCIntMatrix::copy(SCIntMatrix& A) {
    for (int i=0; i<_sizelocal; i++) {
        A._matrix[i] = _matrix[i];
    }
    return 0;
}
#endif
