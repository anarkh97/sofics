#include <Utils.d/dbg_alloca.h>
#include <cstdio>

#include <Comm.d/Communicator.h>
extern Communicator *structCom;

#include <Utils.d/linkfc.h>

#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, double *e, const int &f,
                  double *g, const int &h, const double &i, double *j, const int &k)
{
 _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}

inline void Tgemv(const char &a, const int &b, const int &c,
                  const complex<double> &d, complex<double> *e, const int &f,
                  complex<double> *g, const int &h, const complex<double> &i, complex<double> *j, const int &k)
{
 _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif

template<class Scalar>
GenDistBLKSparse<Scalar>::GenDistBLKSparse(const Connectivity *cn, const EqNumberer *dsa,
                                           double trbm, const SolverCntl &scntl, int fRow, int nRow) :
 GenBLKSparseMatrix<Scalar>(cn, dsa, trbm, scntl)
{
 firstRow = fRow;
 numRows  = nRow;

 // allocate memory for the number of rows we wish to store
 nlines = new Scalar[numRows*this->numUncon];
}

template<class Scalar>
GenDistBLKSparse<Scalar>::~GenDistBLKSparse()
{
 if(nlines) { delete [] nlines; nlines = 0; }
}

template<class Scalar>
void
GenDistBLKSparse<Scalar>::factor()
{
 GenBLKSparseMatrix<Scalar>::factor();

 // zero the n-lines
 int i;
 for(i=0; i<numRows*this->numUncon; ++i)
   nlines[i]=0.0;   

 Scalar **rows = (Scalar **) dbg_alloca(sizeof(Scalar*)*numRows);

 for(i=0; i<numRows; ++i) {
   nlines[i*this->numUncon+firstRow+i] = 1.0;
   rows[i] = nlines+i*this->numUncon;
 }

 for(i=0; i<numRows; ++i)
   GenBLKSparseMatrix<Scalar>::reSolve(rows[i]);

 // Delete the Memory to store GtG in Sparse
 delete [] this->lnz;      this->lnz      = nullptr;
 delete [] this->lindx;    this->lindx    = nullptr;
 delete [] this->xlindx;   this->xlindx   = nullptr;
 delete [] this->snode;    this->snode    = nullptr;
 delete [] this->xsuper;   this->xsuper   = nullptr;
 delete [] this->xlnz;     this->xlnz     = nullptr;
 delete [] this->perm;     this->perm     = nullptr;
 delete [] this->invp;     this->invp     = nullptr;
 delete [] this->invsuper; this->invsuper = nullptr;
}

template<class Scalar>
void
GenDistBLKSparse<Scalar>::reSolve(Scalar *rhs)
{  
 Scalar *partialSum = (Scalar *)dbg_alloca(sizeof(Scalar)*numRows);

 Tgemv('T', this->numUncon, numRows, 1.0, nlines, this->numUncon,
       rhs, 1, 0.0, partialSum, 1);

 // zero the rhs
 int i;
 for(i=0; i<this->numUncon; ++i)
   rhs[i] = 0.0;

 for(i=0; i<numRows; ++i)
   rhs[firstRow+i] = partialSum[i];

#ifdef DISTRIBUTED
 structCom->globalSum(this->numUncon, rhs);
#endif
}

template<class Scalar>
void
GenDistBLKSparse<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
 if(this->neqs() > 0) reSolve(rhs.data());
}


