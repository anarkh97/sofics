/*****************************************************************************
 *                   Copyright (C) 1999 CMSoft                               *
 *                                                                           *
 *  These lines of code and declarations contain unpublished proprietary     *
 *  information of CMSoft. They may not be copied or duplicated in whole     *
 *  or part without prior authorization from CMSoft.                         *
 *                                                                           *
 *****************************************************************************/

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <Utils.d/linkfc.h>
#include <Math.d/Vector.h>
#include <Math.d/BLAS.h>


inline void Tgesv(const int &a, const int &b, double *c, const int &d, int *e, double *f,
                  const int &g, int &h)
{
 _FORTRAN(dgesv)(a,b,c,d,e,f,g,h);
}

inline void Tgesv(const int &a, const int &b, complex<double> *c, const int &d, int *e, complex<double> *f,
                  const int &g, int &h)
{
 _FORTRAN(zgesv)(a,b,c,d,e,f,g,h);
}

#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, const double *e, const int &f,
                  const double *g, const int &h, const double &i, double *j, const int &k)
{
 _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}

inline void Tgemv(const char &a, const int &b, const int &c,
                  const complex<double> &d, const complex<double> *e, const int &f,
                  const complex<double> *g, const int &h, const complex<double> &i, complex<double> *j, const int &k)
{
 _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif 


template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix() 
{
 nrow    = 0;
 ncolumn = 0;
 ld      = 0;
 v       = 0;
}

template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix(int nr) 
{
 nrow    = nr;
 ncolumn = nr;
 ld      = nr;
 v       = new Scalar[nrow*ncolumn];
}

template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix(int nr, int nc) 
{
 nrow    = nr;
 ncolumn = nc;
 ld      = nc;
 v       = new Scalar[nrow*ncolumn];
}

template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix(int nr, int nc, Scalar initialValue)
{
 nrow    = nr;
 ncolumn = nc;
 ld      = nc;
 v       = new Scalar[nrow*ncolumn];
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = initialValue;
}

template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix(int nr, int nc, Scalar* array)
{
 nrow    = nr;
 ncolumn = nc;
 ld      = nc;
 v       = new Scalar[nrow*ncolumn];
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = array[i];
}

// MLX This routine was wrong June 28th 99
template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix(const GenFSFullMatrix<Scalar> &m) 
{
 nrow    = m.nrow;
 ncolumn = m.ncolumn;
 ld      = ncolumn;
 v       = new Scalar[nrow*ld];
 int i,j;
 for(i=0; i < nrow; ++i)
   for(j = 0; j < ncolumn; ++j)
     (*this)[i][j] = m[i][j];
}

template<class Scalar> 
GenFSFullMatrix<Scalar>::GenFSFullMatrix(const GenFSFullMatrix<Scalar> &m, int nr, int sr,int nc, int sc)
{
 // sr = start row number
 // sc = start column number

 nrow    = nr;
 ncolumn = nc;
 ld      = nc;
 v       = new Scalar[nrow*ncolumn] ;

 int i, j;
 for(i=0; i < nrow; ++i)
   for(j=0; j < ncolumn; ++j)
     (*this)[i][j] = m[i+sr][j+sc] ;
}

template<class Scalar> 
GenFSFullMatrix<Scalar>::~GenFSFullMatrix()
{
  if(v) { 
    delete [] v;
    v=0;
  }
}

template<class Scalar>
void
GenFSFullMatrix<Scalar>::copy(Scalar* array)
{
  for(int i=0; i<nrow*ncolumn; ++i) v[i] = array[i];
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::reSize(int nr, int nc, Scalar initVal)
{
 nrow    = nr;
 ncolumn = nc;
 ld = nc;
 if(v) { 
   delete [] v;
   v = 0;
 }
 v = new Scalar[nrow*ncolumn];
 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = initVal;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::operator=(const GenFSFullMatrix<Scalar> &m)
{
 if(m.v == v) return;

 if(nrow != m.nrow || ncolumn != m.ncolumn) {
   if(v) delete [] v ;
   nrow    = m.nrow; 
   ncolumn = m.ncolumn ;
   ld = ncolumn;
   v       = new Scalar[nrow*ld] ;
 }

 int i,j;
 for(i=0; i < nrow; ++i)
   for(j = 0; j < ncolumn; ++j)
     (*this)[i][j] = m[i][j];
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::operator=(const Scalar c)
{
 int i,j;
 for(i=0; i<nrow; ++i)
   for(j = 0; j < ncolumn; ++j)
     (*this)[i][j] = c;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::operator*=(const Scalar c)
{
 int i,j;
 for(i=0; i<nrow; ++i)
   for(j = 0; j < ncolumn; ++j)
     (*this)[i][j] *= c;
}

template<class Scalar> 
GenFSFullMatrix<Scalar>
GenFSFullMatrix<Scalar>::operator*(GenFSFullMatrix<Scalar> &m)
{
 if(ncolumn != m.nrow) {
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar>-GenFSFullMatrix<Scalar> operator mult\n");
   return GenFSFullMatrix<Scalar>(1,1); 
 }

 GenFSFullMatrix<Scalar> C(nrow,m.ncolumn) ;

 //mult(m,C);

 int i,j,k;
 for(i = 0 ; i < nrow ; ++i)
  for(j=0; j < m.ncolumn ; ++j) {
    C[i][j] = 0.0;
    for(k = 0;  k < ncolumn; ++k)
      C[i][j] += (*this)[i][k] * m[k][j];
   }

 return C;
}

template<class Scalar> 
GenFSFullMatrix<Scalar>
GenFSFullMatrix<Scalar>::operator%(GenFSFullMatrix<Scalar> &B)
{
 if(ncolumn != B.ncolumn) { 
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar>-GenFSFullMatrix<Scalar> operator percent\n");
   return GenFSFullMatrix<Scalar>(1,1) ; 
 }
 GenFSFullMatrix<Scalar> C(nrow,B.nrow) ;

 //multTr(B,C);

 int i,j,k;
 for(i=0; i<nrow; ++i)
  for(j=0; j<B.nrow; ++j) {
    C[i][j] = 0.0;
    for(k=0;  k<ncolumn; ++k)
      C[i][j] += (*this)[i][k] * B[j][k];
   }

 return C;
}

template<class Scalar> 
// C = A^T*B
GenFSFullMatrix<Scalar>
GenFSFullMatrix<Scalar>::operator^(GenFSFullMatrix<Scalar> &B)
{
 if(nrow != B.nrow) { 
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar>-GenFSFullMatrix<Scalar> operator ^\n");
   return GenFSFullMatrix<Scalar>(1,1); 
 }

 GenFSFullMatrix<Scalar> C(ncolumn,B.ncolumn);

 Tgemm('N','T',B.numCol(),numCol(),numRow(),1.0,
                 B.data(),B.ldim(),data(),ldim(),0.0,C.data(),C.ldim());

 return C;
}

template<class Scalar> 
GenFSFullMatrix<Scalar>
GenFSFullMatrix<Scalar>::invert()
{
 if(nrow != ncolumn) {
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar> invert() non-square matrix\n");
   return GenFSFullMatrix<Scalar>(1,1); 
 }
 GenFSFullMatrix<Scalar> res(*this) ;
 GenFSFullMatrix<Scalar> inv(nrow,nrow) ;
 int info;
 int *ipiv = new int[nrow];
 int i,j;
 for(i = 0; i < nrow; ++i) {
   for(j = 0; j < nrow; ++j)
     inv[i][j] = 0.0;
   inv[i][i] = 1.0;
 }

 Tgesv(nrow, nrow, res.data(), nrow, ipiv, inv.data(), nrow, info);
/*
  for(i = 0 ; i < nrow ; ++i)
    for(j=i+1; j < nrow ; ++j) {
       Scalar p = res[j][i] = -res[j][i]/res[i][i] ;
       for(k = i+1; k < ncolumn; ++k)
         res[j][k] += p*res[i][k] ;
    }
 for(i = nrow-1 ; i >=0 ; --i)
  for(j=0; j < ncolumn; ++j)
   {
    Scalar piv = 1.0/res[i][i] ;
    if(j < i)  inv[i][j] = res[i][j] ;
    if(j == i) inv[i][j] = 1.0 ;
    if(j > i)  inv[i][j] = 0.0 ;
    for(k = i+1; k <nrow; ++k)
       inv[i][j] -= res[i][k]*inv[k][j] ;
    inv[i][j] *= piv ;
   }
*/
 return inv ;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::luFactor()
{
 if(nrow != ncolumn) {
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar> luFactor() non-square matrix\n");
   return; 
 }

 int i,j,k;
 for(i=0; i<nrow; ++i) {
   Scalar invD = (*this)[i][i] = 1.0/(*this)[i][i];
   for(j=i+1; j < nrow; ++j) {
      Scalar c = ( (*this)[j][i] *= invD );
      for(k = i+1; k < nrow; ++k)
        (*this)[j][k] -= c*(*this)[i][k];
   }
 }
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::symLuFactor()
{
 // LOWER Factorization!!
 if(nrow != ncolumn) {
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar> luFactor() non-square matrix\n");
   return;
 }

 int i,j,k;
 for(i=0; i<nrow; ++i) {
   Scalar invD = (*this)[i][i] = 1.0/(*this)[i][i];
   for(j=i+1; j < nrow; ++j) {
      Scalar c = (*this)[j][i] * invD ;
      for(k = j; k < nrow; ++k)
        (*this)[k][j] -= c*(*this)[k][i];
   }
 }
}

// Deactivated since buggy
/*template<class Scalar> 
void
GenFSFullMatrix<Scalar>::solve(Scalar *x)
{
 int i,j;
 // Forward elimination
 for(i=0; i<nrow; ++i) {
   for(j=i+1; j < nrow; ++j) 
      x[j] -= (*this)[i][j]*x[i];
 }

 // Backward substitution
 for(i=nrow; i--; ) { // BUG: the index i is out of bounds ! 
   for(j = i+1; j < nrow; ++j)
     x[i] -= (*this)[j][i]*x[j];
   x[i] *= (*this)[i][i];
 }
}*/

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::print(const char *msg, const char *msg2)
{
 if(*msg) std::cerr << msg << "\n";
 int i,j;
 for(i = 0 ; i < nrow ; ++i) {
   for(j=0; j < ncolumn ; ++j)
     std::cerr << msg2 << "(" << i+1 << "," << j+1 << ")" << (*this)[i][j] << std::endl;
   std::cerr << "\n";
 }
}

template<class Scalar> 
double
GenFSFullMatrix<Scalar>::max()
{
 double max = ScalarTypes::norm( (*this)[0][0] );
 int i,j;
 for(i = 0; i<nrow; ++i)
   for(j = 0; j<ncolumn; ++j)
     if( ScalarTypes::norm((*this)[i][j]) > max) 
       max = ScalarTypes::norm( (*this)[i][j] );

 return max;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::zero()
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = 0.0;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::identity()
{
 zero();

 int i;
 for(i=0; i<nrow; ++i)
   (*this)[i][i] = 1.0;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::add(GenFSFullMatrix<Scalar> &mat, int fRow, int fCol)
{
  // fRow = first Row to position submatrix mat
  // fCol = first Column to position submatrix mat

  int mrow = mat.numRow();
  int mcol = mat.numCol();

  int icol,irow;
  for(icol = 0; icol < mcol; ++icol) {
    for(irow = 0; irow < mrow; ++irow) {
      (*this)[fRow+irow][fCol+icol] += mat[irow][icol];
      }
  }
}

// mult = matrix-vector multiplication, y = alpha*A*x + beta*y
template<class Scalar> 
void
GenFSFullMatrix<Scalar>::mult(const GenVector<Scalar> &x, GenVector<Scalar> &y, 
                              Scalar alpha, Scalar beta)
{ 
 if(ldim() == 0) return;
   Tgemv('T',numCol(),numRow(),alpha,data(),ldim(),x.data(),1,
                                       beta,y.data(),1);
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::mult(Scalar *x, Scalar *y,Scalar alpha, Scalar beta)
{
 if(ldim() == 0) return;
 Tgemv('T',numCol(),numRow(),alpha,data(),ldim(),x,1,beta,y,1);
}

// trMult = transpose matrix-vector multiplication, y = A^t*x
template<class Scalar> 
void
GenFSFullMatrix<Scalar>::trMult(const GenVector<Scalar> &x, GenVector<Scalar> &y, 
                                Scalar alpha, Scalar beta)
{ 
 if(ldim() == 0) return;
 Tgemv('N',numCol(),numRow(),alpha,data(),ldim(),x.data(),1,beta,
                 y.data(),1);
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::trMult(const Scalar *x, Scalar *y,Scalar alpha, Scalar beta)
{
 if(ldim() == 0) return;
 Tgemv('N',numCol(),numRow(),alpha,data(),ldim(),x,1,beta,y,1);
}

// C=A^T*B^T
template<class Scalar> 
void
GenFSFullMatrix<Scalar>::trMultTr(const GenFSFullMatrix<Scalar> &B, GenFSFullMatrix<Scalar> &C,
                       Scalar alpha,Scalar beta)
{
 if(ldim() == 0) return;
 Tgemm('T','T',B.numRow(),numCol(),numRow(),alpha,
                 B.data(),B.ldim(),data(),ldim(),beta,C.data(),C.ldim());
}

// C = A*B^T
template<class Scalar> 
void
GenFSFullMatrix<Scalar>::multTr(const GenFSFullMatrix<Scalar> &B,GenFSFullMatrix<Scalar> &C,Scalar alpha,Scalar beta)
{
 if(C.numRow() != numRow() || C.numCol() != B.numRow())
   fprintf(stderr, "Inconsistent sizes\n");
 Tgemm('T','N',B.numRow(),numRow(),numCol(),alpha,
                 B.data(),B.ldim(),data(),ldim(),beta,C.data(),C.ldim());
}

// C = A^T*B 
template<class Scalar> 
void
GenFSFullMatrix<Scalar>::trMult(const GenFSFullMatrix<Scalar> &B,GenFSFullMatrix<Scalar> &C,Scalar alpha,Scalar beta)
{
 if(numRow() != B.numRow()) {
   fprintf(stderr," *** ERROR: incompatible dimensions GenFSFullMatrix<Scalar>::trMult\n");
 }                                                                              
 Tgemm('N','T',B.numCol(),numCol(),numRow(),alpha,
                 B.data(),B.ldim(),data(),ldim(),beta,C.data(),C.ldim());
}

// C = A*B 
template<class Scalar> 
void
GenFSFullMatrix<Scalar>::mult(const GenFSFullMatrix<Scalar> &B,GenFSFullMatrix<Scalar> &C,Scalar alpha, Scalar beta)
{
 Tgemm('N','N',B.numCol(),numRow(),numCol(),alpha,
                 B.data(),B.ldim(),data(),ldim(),beta,C.data(),C.ldim());
}


template<class Scalar> 
GenFSFullMatrix<Scalar>
GenFSFullMatrix<Scalar>::transpose()
{
 GenFSFullMatrix<Scalar> res(ncolumn,nrow);

 int i,j;
 for(i=0; i<nrow; ++i)
   for(j=0; j<ncolumn; ++j)
      res[j][i] = (*this)[i][j];

 return res;
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::Lm1Mult(Scalar *a, int nc, int lda)
{ 
 int i, j, col;
 // Forward elimination (L^-1)
 for(i=0; i<nrow; ++i)
   for(col = 0; col < nc ; ++col) {
     Scalar c = a[i+col*lda]*(*this)[i][i];
     for(j=i+1; j < nrow; ++j) 
       a[j+col*lda] -= (*this)[j][i]*c;
   }
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::Um1Mult(Scalar *a, int nc, int lda)
{ 
 int i, j, col;
 // Backward substitution
 for(i=nrow; i--; )
   for(col = 0; col < nc ; ++col) {
     a[i+col*lda] *= (*this)[i][i];
     for(j=0; j < i; ++j) 
       a[j+col*lda] -= (*this)[i][j]*a[i+col*lda];
   }
}

template<class Scalar> 
void
GenFSFullMatrix<Scalar>::Um1TMult(Scalar *a, int nc, int lda)
{ 
 int i, j, col;
 // Forward elimination (U^T^-1)
 for(i=0; i<nrow; ++i)
   for(col = 0; col < nc ; ++col) {
     Scalar c = (a[i+col*lda] *= (*this)[i][i]);
     for(j=i+1; j < nrow; ++j) 
       a[j+col*lda] -= (*this)[j][i]*c;
   }
}

template<class Scalar> 
int
GenFSFullMatrix<Scalar>::symLuFactor(int *perm, double tol, Scalar *origDiag, int *sing)
{
 // LOWER Factorization!! with permutation and detection of zero pivots
 if(nrow != ncolumn) {
   fprintf(stderr,"*** ERROR: GenFSFullMatrix<Scalar> luFactor() non-square matrix\n");
   return - 1;
 }

 int nzem = 0;
 int i,j,k;
 for(i=0; i<nrow; ++i)
   perm[i] = i;

// Jing Li's problem, check with Michel
 for(i=0; i<nrow; ++i) {
#ifdef USE_PERM
   // First find the largest diagonal term
   Scalar diagMax = (*this)[i][i];
   int maxIndex = i;
   for(k = i; k < nrow; ++k)
     if( ScalarTypes::norm((*this)[k][k]) > ScalarTypes::norm(diagMax) ) {
         diagMax = (*this)[k][k];
         maxIndex = k;
     }
   // Perform and record the permutation
   //if(perm[i] != i) {
   if(maxIndex != i) {
     //fprintf(stderr, "Permuting\n");
     std::swap(perm[i], perm[maxIndex]);
     for(k = 0; k < i; ++k)
       std::swap((*this)[i][k], (*this)[maxIndex][k]);
     std::swap( (*this)[i][i], (*this)[maxIndex][maxIndex]);
     for(k = i+1; k < maxIndex; ++k) 
       std::swap((*this) [k][i], (*this)[maxIndex][k]);
     for(k = maxIndex+1; k < nrow; ++k)
       std::swap((*this)[k][i], (*this)[k][maxIndex]);

   }
#endif
   // Detect singularity
   if( ScalarTypes::norm((*this)[i][i]) <= (std::abs(tol)*ScalarTypes::norm(origDiag[perm[i]])) ) {
     nzem++;
     if(sing != 0) sing[i] = 1;  // PJSA
// For Debugging
/*
     if(origDiag[perm[i]] != 0.0)
       fprintf(stderr, "irow = %d, Found RBM: %e\n", i,
               (*this)[i][i]/(tol*origDiag[perm[i]]));
     else
       fprintf(stderr, "irow = %d, Found RBM: 0.0\n", i);
*/
     (*this)[i][i] = 0.0;
     continue;
   }
   Scalar invD = (*this)[i][i] = 1.0/(*this)[i][i];
   for(j=i+1; j < nrow; ++j) {
      Scalar c = (*this)[j][i] * invD ;
      for(k = j; k < nrow; ++k)
        (*this)[k][j] -= c*(*this)[k][i];
   }
 }
 return nzem;
}

template<class Scalar>
void
GenFSFullMatrix<Scalar>::Lm1Mult(Scalar *a, int nc, int lda, int *perm)
{ 
 int i, j, col;
 // Forward elimination (L^-1)
 for(i=0; i<nrow; ++i)
   for(col = 0; col < nc ; ++col) {
#ifdef USE_PERM
     Scalar c = a[perm[i]+col*lda]*(*this)[i][i];
     for(j=i+1; j < nrow; ++j) 
       a[perm[j]+col*lda] -= (*this)[j][i]*c;
#else
     Scalar c = a[i+col*lda]*(*this)[i][i];
     for(j=i+1; j < nrow; ++j)
       a[j+col*lda] -= (*this)[j][i]*c;
#endif
   }
}

template<class Scalar>
void
GenFSFullMatrix<Scalar>::Um1Mult(Scalar *a, int nc, int lda, int *perm)
{ 
 int i, j, col;
 // Backward substitution
 for(i=nrow; i--; )
   for(col = 0; col < nc ; ++col) {
#ifdef USE_PERM
     a[perm[i]+col*lda] *= (*this)[i][i];
     for(j=0; j < i; ++j) 
       a[perm[j]+col*lda] -= (*this)[i][j]*a[perm[i]+col*lda];
#else
     a[i+col*lda] *= (*this)[i][i];
     for(j=0; j < i; ++j)
       a[j+col*lda] -= (*this)[i][j]*a[i+col*lda];
#endif
   }
}

template<class Scalar>
void
GenFSFullMatrix<Scalar>::Um1TMult(Scalar *a, int nc, int lda, int *perm)
{ 
 int i, j, col;
 // Forward elimination (U^T^-1)
 for(i=0; i<nrow; ++i)
   for(col = 0; col < nc ; ++col) {
#ifdef USE_PERM
     Scalar c = (a[perm[i]+col*lda] *= (*this)[i][i]);
     for(j=i+1; j < nrow; ++j) 
       a[perm[j]+col*lda] -= (*this)[j][i]*c;
#else
     Scalar c = (a[i+col*lda] *= (*this)[i][i]);
     for(j=i+1; j < nrow; ++j)
       a[j+col*lda] -= (*this)[j][i]*c;
#endif
   }
}
