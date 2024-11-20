#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <Comm.d/Communicator.h>
#include <Utils.d/linkfc.h>

extern Communicator *structCom;

extern "C" {
  void _FORTRAN(dgemv)(const char &, const int &,const int &,
                       const double &, double *, const int &,
                       double *, const int &, const double &, double *, const int &);
  void _FORTRAN(zgemv)(const char &, const int &,const int &,
                       const complex<double> &, complex<double> *, const int &,
                       complex<double> *, const int &, const complex<double> &, complex<double> *, const int &);
}
                                                                                                                                  
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

template<class Scalar, class AnyVector, class AnyOperator>
GenDistPCGSolver<Scalar,AnyVector,AnyOperator>::GenDistPCGSolver(AnyOperator *K, int precno, int maxiter, 
                                                                 double tolerance, int fRow, int nRow) :
 GenPCGSolver<Scalar,AnyVector,AnyOperator>(K, precno, maxiter, tolerance)
{
 numUncon = K->dim();
 firstRow = fRow;
 numRows  = nRow;

 // allocate memory for the number of rows we wish to store
 nlines = new double[numRows*numUncon];
}

template<class Scalar, class AnyVector, class AnyOperator>
GenDistPCGSolver<Scalar,AnyVector,AnyOperator>::~GenDistPCGSolver()
{
 if(nlines) { delete [] nlines; nlines = 0; }
}

template<class Scalar, class AnyVector, class AnyOperator>
void
GenDistPCGSolver<Scalar,AnyVector,AnyOperator>::factor()
{
 GenPCGSolver<Scalar,AnyVector,AnyOperator>::factor();

 // zero the n-lines
 int i;
 for(i=0; i<numRows*numUncon; ++i)
   nlines[i]=0.0;   

 Scalar **rows = (Scalar **) dbg_alloca(sizeof(Scalar*)*numRows);

 for(i=0; i<numRows; ++i) {
   nlines[i*numUncon+firstRow+i] = 1.0;
   rows[i] = nlines+i*numUncon;
 }

 GenPCGSolver<Scalar,AnyVector,AnyOperator>::reSolve(numRows, rows);

 // Delete the Memory to store GtG in Sparse
 //PCGSolver::getOperator().deleteMemory(); 

}

template<class Scalar, class AnyVector, class AnyOperator>
void
GenDistPCGSolver<Scalar,AnyVector,AnyOperator>::reSolve(Scalar *rhs) const
{
 Scalar *partialSum = (Scalar *)dbg_alloca(sizeof(Scalar)*numRows);

 Tgemv('T',numUncon,numRows,1.0,nlines,numUncon,
       rhs, 1, 0.0, partialSum, 1);

 // zero the rhs
 int i;
 for(i=0; i<numUncon; ++i)
   rhs[i] = 0.0;

 for(i=0; i<numRows; ++i)
   rhs[firstRow+i] = partialSum[i];

#ifdef DISTRIBUTED
 structCom->globalSum(numUncon, rhs);
#endif
}

