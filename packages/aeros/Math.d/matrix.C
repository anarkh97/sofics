#include <cstdio>
#include <strings.h>
#include <iostream>
#include <complex>

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Threads.d/PHelper.h>
#include <Utils.d/linkfc.h>
#include "BLAS.h"
#include "matrix.h"

#ifdef TFLOP
#include <bsd.h>
#endif

#ifdef WINDOWS
#define bzero(s, n) __builtin_memset (s, 0, n)
#endif

extern "C" {
  // triangular factorization of a real general matrix using Gaussian elimination with complete pivoting  
  void _FORTRAN(dgecp)(const int &, double *, const int &, const double &, const int &, int &, int *, int *, int &);
  void _FORTRAN(zgecp)(const int &, complex<double> *, const int &, const double &, const int &, int &, int *, int *, int &);
  // solves a real general linear system using the triangular factorization computed by dgecp/zgecp
  void _FORTRAN(dgers)(const int &, const int &, double *, const int &, const int &, const int *, 
                       const int *, double *, const int &, int &);
  void _FORTRAN(zgers)(const int &, const int &, complex<double> *, const int &, const int &, const int *, 
                       const int *, complex<double> *, const int &, int &);
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

#ifndef _TGECP__
#define _TGECP__
inline void Tgecp(const int &a, double *b, const int &c, const double &d, const int &e, int &f, int *g, int *h, int &i)
{
 _FORTRAN(dgecp)(a,b,c,d,e,f,g,h,i);
}
inline void Tgecp(const int &a, complex<double> *b, const int &c, const double &d, const int &e, int &f, int *g, int *h, int &i)
{
 _FORTRAN(zgecp)(a,b,c,d,e,f,g,h,i);
}
#endif

#ifndef _TGERS__
#define _TGERS__
inline void Tgers(const int &a, const int &b, double *c, const int &d, const int &e, const int *f,
                       const int *g, double *h, const int &i, int &j)
{
  _FORTRAN(dgers)(a,b,c,d,e,f,g,h,i,j);
}
inline void Tgers(const int &a, const int &b, complex<double> *c, const int &d, const int &e, const int *f,
                       const int *g, complex<double> *h, const int &i, int &j)
{
  _FORTRAN(zgers)(a,b,c,d,e,f,g,h,i,j);
}
#endif


template<class Scalar> 
GenFullM<Scalar>::GenFullM() 
{
 nrow    = 0;
 ncolumn = 0;
 v = 0;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(int nr) 
{
 nrow    = nr;
 ncolumn = nr;
 if(nrow*ncolumn == 0)
   v = 0;
 else
   v = new Scalar[nrow*ncolumn];
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(int nr, int nc) 
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn == 0)
   v = 0;
 else
   v = new Scalar[nrow*ncolumn];
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar>
GenFullM<Scalar>::GenFullM(int nr, int nc, Scalar init_val)
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn == 0)
   v = 0;
 else
   v = new Scalar[nrow*ncolumn];

 for(int i=0; i<nrow*ncolumn; ++i) v[i] = init_val;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(const GenFullM<Scalar> &m) 
{
 nrow    = m.nrow;
 ncolumn = m.ncolumn;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn];
 else v = 0;

 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = m.v[i];
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar>
GenFullM<Scalar>::GenFullM(const GenFullM<Scalar> &m, double scale_factor)
{
 nrow    = m.nrow;
 ncolumn = m.ncolumn;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn];
 else v = 0;

 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = m.v[i]*scale_factor;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(const GenFullM<Scalar> &m, int nr, int sr, int nc, int sc)
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn] ;
 else v = 0;

 int i,j;
 for(i=0; i < nrow; ++i)
  for(j=0; j < ncolumn; ++j)
    (*this)[i][j] = m[i+sr][j+sc] ;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(const GenFullM<Scalar> &m, int nr, int *rows, int nc, int *cols)
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn] ;
 else v = 0;

 int i,j;
 for(i=0; i < nrow; ++i)
  for(j=0; j < ncolumn; ++j)
    (*this)[i][j] = m[rows[i]][cols[j]] ;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(const GenFullM<Scalar> &m, int nr, int *rows, int nc, int sc)
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn];
 else v = 0;

 int i,j;
 for(i=0; i < nrow; ++i)
  for(j=0; j < ncolumn; ++j)
    (*this)[i][j] = m[rows[i]][j+sc] ;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(const GenFullM<Scalar> &m, int nr, int sr, int nc, int *cols)
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn] ;
 else v = 0;

 int i,j;
 for(i=0; i < nrow; ++i)
  for(j=0; j < ncolumn; ++j)
    (*this)[i][j] = m[i+sr][cols[j]] ;
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::GenFullM(Scalar *data, int nr, int nc, int flag)
{
 int i, j;
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn];
 else v = 0;

 if(flag == 0) // data is stored by rows, same as v
   for(i=0; i < nrow*ncolumn; ++i)
     v[i] = data[i];
 else
   for(i=0; i<nrow; ++i) // data is stored by columns, unlike v
     for(j=0; j<ncolumn; ++j)
       v[j + i*ncolumn] = data[i + j*nrow];
 iprow = 0;
 ipcol = 0;
 ndef = 0;
}

template<class Scalar> 
GenFullM<Scalar>::~GenFullM()
{
#if !defined(DEBUG_OPENMP)
 if(v) { delete [] v; v=0; }
#endif
 if(iprow) { delete [] iprow; iprow = 0; }
 if(ipcol) { delete [] ipcol; ipcol = 0; }
}

template<class Scalar>
void
GenFullM<Scalar>::copy(Scalar *array)
{
  for(int i=0; i<nrow*ncolumn; ++i) v[i] = array[i];
}

template<class Scalar> 
void
GenFullM<Scalar>::setNewSize(int nr, int nc, Scalar initVal)
{
 nrow    = nr;
 ncolumn = nc;
 if(v) delete [] v;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn];
 else v = 0; 

 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = initVal;
}

template<class Scalar> 
void
GenFullM<Scalar>::setNewSize(int nr, Scalar initVal)
{
 nrow = nr;
 ncolumn = nr;
 if(v) delete [] v;
 if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn];
 else v = 0;

 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = initVal;
}

template<class Scalar> 
void
GenFullM<Scalar>::operator=(const GenFullM<Scalar> &m)
{
 if(m.v == v) return ;

 if(nrow != m.nrow || ncolumn != m.ncolumn) {
   if(v) delete [] v ;
   nrow = m.nrow ; ncolumn = m.ncolumn ;
   if(nrow*ncolumn > 0) v = new Scalar[nrow*ncolumn] ;
   else  v = 0;
 }

 // copy data
 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = m.v[i] ;

}

template<class Scalar> 
void
GenFullM<Scalar>::operator=(const Scalar c)
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = c;
}

template<class Scalar>
void
GenFullM<Scalar>::operator*=(const Scalar c)
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] *= c;
}

//template<class Scalar>
//GenFullM<Scalar>
//GenFullM<Scalar>::operator+(const GenFullM<Scalar> &m) const
//{
//  GenFullM<Scalar> res(nrow+m.nrow,ncolumn);
//  int i,j;
//  for(i=0; i<nrow; ++i)
//    for(j=0; j<ncolumn; ++j)
//      res[i][j] = (*this)[i][j];
//  for(i=nrow; i<m.nrow; ++i)
//    for(j=0; j<ncolumn; ++j)
//      res[i+nrow][j] + m[i][j];
//  return res;
//}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::operator*(const GenFullM<Scalar> &m) const
{
 if(ncolumn != m.nrow) {
   std::cerr << " *** ERROR in GenFullM<Scalar>::operator*(GenFullM<Scalar> &m), ncolumn != m.nrow \n";
   return GenFullM<Scalar>(1,1) ; //error
 }
 GenFullM<Scalar> res(nrow,m.ncolumn) ;
 int i,j,k;
 for(i = 0 ; i < nrow ; ++i)
  for(j=0; j < m.ncolumn ; ++j) {
    res[i][j] = 0.0;
    for(k = 0;  k < ncolumn; ++k)
      res[i][j] += (*this)[i][k] * m[k][j];
   }
 return res;
}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::operator*=(const GenFullM<Scalar> &m)
{
 if(ncolumn != m.nrow) {
   std::cerr << " *** ERROR in GenFullM<Scalar>::operator*=(GenFullM<Scalar> &m), ncolumn != m.nrow \n";
   return GenFullM<Scalar>(1,1) ; //error
 }
 if(m.ncolumn != ncolumn) {
   std::cerr << " *** ERROR in GenFullM<Scalar>::operator*=(GenFullM<Scalar> &m), m.ncolumn != ncolumn \n";
   return GenFullM<Scalar>(1,1) ; //error
 }
 GenFullM<Scalar> res(nrow,m.ncolumn) ;
 int i,j,k;
 for(i = 0 ; i < nrow ; ++i)
  for(j=0; j < m.ncolumn ; ++j) {
    res[i][j] = 0.0;
    for(k = 0;  k < ncolumn; ++k)
      res[i][j] += (*this)[i][k] * m[k][j];
   }

 *this = res;
 return (*this);
}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::operator^=(const GenFullM<Scalar> &m)
{
 if(nrow != m.nrow) {
   std::cerr << " *** ERROR in GenFullM<Scalar>::operator^=(GenFullM<Scalar> &m), nrow != m.nrow \n";
   return GenFullM<Scalar>(1,1) ; //error
 }
 if(m.ncolumn != nrow) {
   std::cerr << " *** ERROR in GenFullM<Scalar>::operator^=(GenFullM<Scalar> &m), m.ncolumn != nrow \n";
   return GenFullM<Scalar>(1,1) ; //error
 }
 GenFullM<Scalar> res(nrow,ncolumn);
 int i,j,k;
 for(j=0; j < ncolumn ; ++j) 
   for(i = 0 ; i < nrow ; ++i) {
    res[i][j] = 0.0;
    for(k = 0;  k < nrow; ++k)
      res[i][j] += m[i][k] * (*this)[k][j];
   }

 *this = res;
 return (*this);
}

template<class Scalar> 
GenVector<Scalar>
GenFullM<Scalar>::operator*(const GenVector<Scalar> &x) const
{
  if(ncolumn != x.size()) return GenVector<Scalar>(1) ; //error
  GenVector<Scalar> res(nrow);
  int i,j;
  for(i = 0 ; i < nrow ; ++i) {
    res[i] = 0.0;
    for(j=0; j < x.size(); ++j)
      res[i] += (*this)[i][j] * x[j];
  }
 return res;
}

template<class Scalar> 
GenVector<Scalar>
GenFullM<Scalar>::operator^(const GenVector<Scalar> &x) const
{
  if(nrow != x.size()) return GenVector<Scalar>(1) ; //error
  GenVector<Scalar> res(ncolumn,0.0);
  int i,j;
  for(i = 0; i <nrow; ++i) {
    for(j=0; j < ncolumn; ++j)
      res[j] += (*this)[i][j] * x[i];
  }
 return res;
 
}

template<class Scalar>
GenFullM<Scalar>
GenFullM<Scalar>::operator += (const GenFullM<Scalar> &M2)
{
 if((ncolumn!=M2.ncolumn)&(nrow!=M2.nrow)){
  std::cerr <<" !!! In GenFullM::operator += : matrices don't have the same size !!!"<< std::endl;
  return *this;
 }  
 int length = ncolumn*nrow;
 int i;
 for(i = 0; i < length; ++i)
    v[i] += M2.v[i];
 return *this;
}

template<class Scalar>
GenFullM<Scalar>
GenFullM<Scalar>::operator -= (const GenFullM<Scalar> &M2)
{
 if((ncolumn!=M2.ncolumn)&(nrow!=M2.nrow)){
  std::cerr <<" !!! In GenFullM::operator -= : matrices don't have the same size !!!"<< std::endl;
  return *this;
 }
 int length = ncolumn*nrow;
 int i;
 for(i = 0; i < length; ++i)
    v[i] -= M2.v[i];
 return *this;
}

template<class Scalar> 
void
GenFullM<Scalar>::mult(Scalar *x, Scalar *y, Scalar alpha, Scalar beta)
{
 Tgemv('T',ncolumn, nrow, alpha, data(), ncolumn,
       x, 1, beta, y, 1);
}

template<class Scalar> 
void
GenFullM<Scalar>::trMult(Scalar *x, Scalar *y, Scalar alpha, Scalar beta)
{
 Tgemv('N',ncolumn, nrow, alpha, data(), ncolumn,
       x, 1, beta, y, 1);
}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::operator%(const GenFullM<Scalar> &m) const
{
 if(ncolumn != m.ncolumn)
		 throw "Incompatible matrix sizes in operator^.";
 GenFullM<Scalar> res(nrow,m.nrow) ;
 int i,j,k;
 for(i=0; i<nrow; ++i)
  for(j=0; j<m.nrow; ++j) {
    res[i][j] = 0.0;
    for(k=0;  k<ncolumn; ++k)
      res[i][j] += (*this)[i][k] * m[j][k];
   }
 return res;
}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::operator^(const GenFullM<Scalar> &m) const
{
 if(nrow != m.nrow)
	 throw "Incompatible matrix sizes in operator^.";

 GenFullM<Scalar> res(ncolumn,m.ncolumn) ;
 int i,j,k;
 for(i=0; i < ncolumn; ++i)
  for(j=0; j < m.ncolumn; ++j)
   {
    res[i][j] = 0.0 ;
    for(k = 0 ;  k < nrow ; ++k)
      res[i][j] += (*this)[k][i] * m[k][j];
   }
 return res;
}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::invert()
{
 if(nrow != ncolumn) { std::cerr << " *** ERROR: GenFullM<Scalar>::invert(), nrow != ncolumn \n"; return GenFullM<Scalar>(); }
 GenFullM<Scalar> res(*this);
 res.factor();
 GenFullM<Scalar> inv(nrow,nrow,0.0);
 for(int i = nrow-1 ; i >=0 ; --i) {
   inv[i][i] = 1.0;
   res.reSolve(inv[i]);
 }
 return inv.transpose();
}

template<class Scalar>
GenFullM<Scalar>
GenFullM<Scalar>::Invert(double tol)
{
 if(nrow != ncolumn) { std::cerr << " *** ERROR: GenFullM<Scalar>::Invert(double tol), nrow != ncolumn \n"; return GenFullM<Scalar>(); }
 GenFullM<Scalar> res(*this);
 res.Factor(tol);
 GenFullM<Scalar> inv(nrow,nrow,0.0);
 for(int i = nrow-1 ; i >=0 ; --i) {
   inv[i][i] = 1.0;
   res.ReSolve(inv[i]);
 }
 return inv.transpose();
}

template<class Scalar> 
GenFullM<Scalar>
GenFullM<Scalar>::transpose()
{
 GenFullM<Scalar> res(ncolumn,nrow);

 int i,j;
 for(i=0; i<nrow; ++i)
   for(j=0; j<ncolumn; ++j)
      res[j][i] = (*this)[i][j]; 

 return res; 
}

template<class Scalar> 
void
GenFullM<Scalar>::factor()
{
	if(nrow != ncolumn)
		throw std::invalid_argument("Trying to factor a non square matrix.");
 int i,j,k;
 for(i=0; i<nrow; ++i) {
   Scalar invD = (*this)[i][i] = 1.0/(*this)[i][i];
   for(j=i+1; j < nrow; ++j) {
      Scalar c = ( (*this)[j][i] *=  invD );
      for(k = i+1; k < nrow; ++k)
        (*this)[j][k] -= c*(*this)[i][k];
   }
 }
}

template<class Scalar>
void
GenFullM<Scalar>::Factor(double tol, bool print_ndef)
{
  // triangular factorization of a real general matrix using Gaussian elimination with complete pivoting
  if(nrow != ncolumn) std::cerr << " *** WARNING: GenFullM<Scalar>::Factor(double tol), nrow != ncolumn \n";
  
  if(iprow) { delete [] iprow; iprow = 0; }
  if(ipcol) { delete [] ipcol; ipcol = 0; }
  iprow = new int[nrow]; // output: row permutation
  ipcol = new int[ncolumn]; // output: column permutation
  int info; // output: error flag
  
  // change ordering of v to be [column][row] instead of [row][column]
  Scalar *v_copy = new Scalar[nrow*ncolumn];
  for(int i=0; i<nrow*ncolumn; ++i) v_copy[i] = v[i];
  for(int i=0; i<nrow; ++i)
    for(int j=0; j<ncolumn; ++j)
      v[i*nrow+j] = v_copy[j*ncolumn+i];
  delete [] v_copy;

  Tgecp(nrow, v, ncolumn, tol, 0, ndef, iprow, ipcol, info);

  if(info != 0) std::cerr << " *** WARNING: error in Tgecp, info = " << info << std::endl;
  else if(ndef > 0 && print_ndef) std::cerr << "Matrix factored by dgecp/zgecp is rank deficient, ndef = " << ndef << std::endl; 
}

template<class Scalar>
void
GenFullM<Scalar>::ReSolve(Scalar *b)
{
  // solves a real general linear system using the triangular factorization computed by Factor(int)
  int info;
  Tgers(nrow, 1, v, nrow, ndef, iprow, ipcol, b, nrow, info);
                                                                                                                                                          
  if(info != 0) std::cerr << " *** WARNING: error in Tgers, info = " << info << std::endl;
}

template<class Scalar> 
void
GenFullM<Scalar>::reSolve(Scalar *x)
{
 int i,j;
 // Forward elimination
 for(i=0; i<nrow; ++i){
   for(j=i+1; j < nrow; ++j) 
      x[j] -= (*this)[j][i]*x[i];
 }

 // Backward substitution
 for(i=nrow; i--; ) {
   for(j = i+1; j < nrow; ++j)
     x[i] -= (*this)[i][j]*x[j];
   x[i] *= (*this)[i][i];
 }
}

template<class Scalar>
void
GenFullM<Scalar>::symmetrize_from_uptriag()
{
 int i, j;
 for(i=0; i<nrow; ++i)
   for(j=0;  j<i; ++j)
     (*this)[i][j] = (*this)[j][i];
}

template<>
inline void
GenFullM<std::complex<double>>::print(const char *msg, const char *msg2,FILE *f)
{
 if(*msg) fprintf(f,"%s\n",msg);
 int i,j;
 for(i = 0 ; i < nrow ; ++i) {
   for(j=0; j < ncolumn ; ++j) {
     fprintf(f,"%s(%d,%d) = %16.11e,%16.11e;",msg2,i+1,j+1,std::real((*this)[i][j]), std::imag((*this)[i][j]));
     fprintf(f,"\n") ;
   }
 }
}

template<class Scalar> 
void
GenFullM<Scalar>::print(const char *msg, const char *msg2,FILE *f)
{
 if(*msg) fprintf(f,"%s\n",msg);
 int i,j;
 for(i = 0 ; i < nrow ; ++i) {
   for(j=0; j < ncolumn ; ++j) {
     fprintf(f,"%s(%d,%d) = %16.11e;",msg2,i+1,j+1,(*this)[i][j]);
     fprintf(f,"\n") ;
   }
 }
}

template<class Scalar>
void
GenFullM<Scalar>::Print()
{
 for(int i = 0 ; i < nrow ; ++i) 
   for(int j=0; j < ncolumn ; ++j) 
     std::cerr << "(" << i+1 << "," << j+1 << ") = " << (*this)[i][j] << std::endl;
}

template<class Scalar> 
double
GenFullM<Scalar>::max()
{
	using ScalarTypes::sqNorm;
	double max = (v) ? sqNorm((*this)[0][0]) : 0.0;
	int i,j;
	for(i = 0; i<nrow; ++i)
		for(j = 0; j<ncolumn; ++j)
			if( sqNorm( (*this)[i][j] ) > max) max = sqNorm( (*this)[i][j] );

	return max;
}

inline double scalarAbs(double d) { return fabs(d); }
inline double scalarAbs(DComplex d) { return abs(d); }

template<class Scalar> 
double
GenFullM<Scalar>::maxAbs()
{
 double max = 0.0;
 int i,j;
 for(i = 0; i<nrow; ++i)
   for(j = 0; j<ncolumn; ++j) 
     if(scalarAbs((*this)[i][j]) > max) max = scalarAbs((*this)[i][j]);
 return max;
}

template<class Scalar> 
void
GenFullM<Scalar>::zero()
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = 0.0;
}

template<class Scalar> 
void
GenFullM<Scalar>::negate()
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = -v[i];
}

template<class Scalar> 
void
GenFullM<Scalar>::clean_up()
{
 if(v) { delete [] v; v=0; }
}

template<>
void
GenAssembledFullM<complex<double> >::addImaginary(const FullSquareMatrix &mat, const int *dofs);

template<>
void
GenAssembledFullM<double>::addImaginary(const FullSquareMatrix &kel, const int *dofs);

template<>
void
GenAssembledFullM<complex<double> >::add(const FullSquareMatrixC &mat, const int *dofs);

template<>
void
GenAssembledFullM<double>::add(const FullSquareMatrixC &kel, const int *dofs);

template<class Scalar> 
void
GenAssembledFullM<Scalar>::add(const FullSquareMatrix &mat, gsl::span<const int> dofs)
{
 int kndof = mat.dim();
 int i,j,ri,cj;
 for(i=0; i<kndof; ++i) {
   if(dofs[i] == -1) continue;
   if((ri = rowMap[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(dofs[j] == -1) continue;
     if( (cj = colMap[dofs[j]]) == -1) continue;
     (*this)[ri][cj] += mat[i][j];
   }
 }
}

template<class Scalar>
void
GenAssembledFullM<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  if(dof < 0) return;
  int cdof = rowMap[dof];
  if(cdof < 0) return;
  (*this)[cdof][cdof] += dmass;
}

template<class Scalar> 
void
GenAssembledFullM<Scalar>::add(Scalar **matrix)
{
 int i,j;
 for(i=0; i<this->nrow; ++i) {
   for(j=0; j<this->ncolumn; ++j) {
     (*this)[i][j] += matrix[i][j];
   }
 }
}

template<class Scalar> 
void
GenAssembledFullM<Scalar>::addBoeing(int nlines, const int *Kai, const int *Kaj,
      const double *nz, int *map, Scalar multiplier)
{
 if(this->numRow() > 0 && this->numCol() > 0) {
   int i, j;
   for(i = 0; i < nlines; ++i) {
     if(map[i] == -1) continue;
     int dofI = rowMap[map[i]];
     if(dofI < 0) continue;
     for(j = Kai[i]; j < Kai[i+1]; j++) {
       if(map[Kaj[j-1]-1] == -1) continue;
       int dofJ = colMap[map[Kaj[j-1]-1]];
       if(dofJ < 0) continue;
       (*this)[dofI][dofJ] += (nz[j-1]*multiplier);
       if(dofI != dofJ)
         (*this)[dofJ][dofI] += (nz[j-1]*multiplier);
     }
   }
 }
}

template<class Scalar> 
void
GenFullM<Scalar>::add(const GenFullM<Scalar> &mat, int fRow, int fCol)
{
  int mrow = mat.numRow();
  int mcol = mat.numCol();

  int icol,irow;
  for(icol = 0; icol < mcol; ++icol) {
    for(irow = 0; irow < mrow; ++irow) {
      (*this)[fRow+irow][fCol+icol] += mat[irow][icol];
    }
  }
}

template<class Scalar>
void
GenFullM<Scalar>::addrows(GenFullM<Scalar> &mat, int *rows)
{
  int mrow = mat.numRow();
  int mcol = mat.numCol();
  if(mcol > ncolumn) { std::cerr << " *** ERROR in GenFullM<Scalar>::addrows(GenFullM<Scalar> &mat, int *rows) \n"; return; }

  int icol,irow;
  for(icol = 0; icol < mcol; ++icol) {
    for(irow = 0; irow < mrow; ++irow) {
      (*this)[rows[irow]][icol] += mat[irow][icol];
    }
  }
}

template<class Scalar>
void
GenFullM<Scalar>::add(const GenFullM<Scalar> &mat, int *rc)
{
  int mrow = mat.numRow();
  int mcol = mat.numCol();
  if(mrow != mcol) { std::cerr << " *** ERROR in GenFullM<Scalar>::add(const GenFullM<Scalar> &mat, int *rc) \n"; return; }

  int icol,irow;
  for(icol = 0; icol < mcol; ++icol) {
    for(irow = 0; irow < mrow; ++irow) {
      (*this)[rc[irow]][rc[icol]] += mat[irow][icol];
    }
  }
}




template<class Scalar> 
void
GenFullM<Scalar>::add(GenVector<Scalar> &vec, int fRow, int fCol)
{
  int size = vec.size();

  int irow;
  for(irow = 0; irow < size; ++irow) {
    (*this)[fRow+irow][fCol] += vec[irow];
  }
}


template<class Scalar> 
void
GenFullM<Scalar>::subtract(GenFullM<Scalar> &mat, int fRow, int fCol)
{
  int mrow = mat.numRow();
  int mcol = mat.numCol();

  int icol,irow;
  for(icol = 0; icol < mcol; ++icol) {
    for(irow = 0; irow < mrow; ++irow) {
      (*this)[fRow+irow][fCol+icol] -= mat[irow][icol];
      }
  }
}

#ifndef _TGEMM__
#define _TGEMM__
inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const double &f, const double *g, const int &h,
                  double *i, const int &j, const double &k, double *l,
                  const int &m)
{
 _FORTRAN(dgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const complex<double> &f, const complex<double> *g, const int &h,
                  complex<double> *i, const int &j, const complex<double> &k, complex<double> *l,
                  const int &m)
{
 _FORTRAN(zgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}
#endif

template<class Scalar> 
void
GenFullM<Scalar>::paralAddProd(GenFullM<Scalar> &A, GenFullM<Scalar> &B, Scalar coefC, Scalar coefAB)
{
 // Check matrix sizes consistency
 if(A.nrow != nrow || B.ncolumn != ncolumn || A.ncolumn != B.nrow) {
   fprintf(stderr, "Inconsitent sizes\n");
   return;
 }

 int nt = threadManager->numThr();
 execParal(nt, this, &GenFullM<Scalar>::subAddProd, nt, &A, &B, coefC, coefAB);
}

template<class Scalar> 
void
GenFullM<Scalar>::subAddProd(int i, int nThreads, GenFullM<Scalar> *A, GenFullM<Scalar> *B, 
                  Scalar coefC, Scalar coefAB)
{
 int dd = nrow/nThreads;
 int remainder = nrow%nThreads;
 int rowStart = dd*i+( (i < remainder) ? i : remainder );
 int rowStop  = dd*(i+1) +( (i+1 < remainder) ? i+1 : remainder );

 Tgemm('N', 'N', ncolumn, rowStop-rowStart, B->nrow, coefAB, 
           B->v, B->ncolumn,
           (*A)[rowStart] , A->ncolumn,
           coefC, (*this)[rowStart], ncolumn);
}

template<class Scalar> 
void
GenFullM<Scalar>::paralAddUTProd(GenFullM<Scalar> &A, GenFullM<Scalar> &B, Scalar coefC, Scalar coefAB)
{
 // Check matrix sizes consistency
 if(A.nrow != nrow || B.nrow != ncolumn || A.ncolumn != B.ncolumn) {
   fprintf(stderr, "Inconsitent sizes\n");
   return;
 }

 int nt = threadManager->numThr();
 execParal(nt, this, &GenFullM<Scalar>::subAddUTProd, nt, &A, &B, coefC, coefAB);
}

template<class Scalar> 
void
GenFullM<Scalar>::subAddUTProd(int i, int nThreads, GenFullM<Scalar> *A, GenFullM<Scalar> *B,
                  Scalar coefC, Scalar coefAB)
{
 int dd = nrow/nThreads;
 int remainder = nrow%nThreads;
 int rowStart = dd*i+( (i < remainder) ? i : remainder );
 int rowStop  = dd*(i+1) +( (i+1 < remainder) ? i+1 : remainder );

 Tgemm('T', 'N', ncolumn, rowStop-rowStart, B->ncolumn, coefAB,
           B->v, B->ncolumn,
           (*A)[rowStart] , A->ncolumn,
           coefC, (*this)[rowStart], ncolumn);
}

template<class Scalar> 
void
GenFullM<Scalar>::subAddTTProd(int i, int nThreads, GenFullM<Scalar> *A, GenFullM<Scalar> *B,
                  Scalar coefC, Scalar coefAB)
{
 int dd = nrow/nThreads;
 int remainder = nrow%nThreads;
 int rowStart = dd*i+( (i < remainder) ? i : remainder );
 int rowStop  = dd*(i+1) +( (i+1 < remainder) ? i+1 : remainder );

 Tgemm('T', 'T', ncolumn, rowStop-rowStart, B->ncolumn, coefAB,
           B->v, B->ncolumn,
           (*A)[rowStart] , A->ncolumn,
           coefC, (*this)[rowStart], ncolumn);
}

template<class Scalar> 
void
GenFullM<Scalar>::paralAddTTProd(GenFullM<Scalar> &A, GenFullM<Scalar> &B, Scalar coefC, Scalar coefAB)
{
 // Check matrix sizes consistency
 if(A.ncolumn != nrow || B.nrow != ncolumn || A.nrow != B.ncolumn) {
   fprintf(stderr, "Inconsitent sizes\n");
   return;
 }

 int nt = threadManager->numThr();
 execParal(nt, this, &GenFullM<Scalar>::subAddUTProd, nt, &A, &B, coefC, coefAB);
}


template<class Scalar> 
void
GenFullM<Scalar>::paralZero()
{
 int nt = threadManager->numThr();
 execParal(nt, this, &GenFullM<Scalar>::subZero, nt);
}



template<class Scalar> 
void
GenFullM<Scalar>::subZero(int i, int nThreads)
{
 int dd = nrow/nThreads;
 int remainder = nrow%nThreads;
 int rowStart = dd*i+( (i < remainder) ? i : remainder );
 int rowStop  = dd*(i+1) +( (i+1 < remainder) ? i+1 : remainder );

#ifndef sgi
 bzero((char *)(v+rowStart*ncolumn), (rowStop-rowStart)*ncolumn*sizeof(Scalar));
#else
 bzero(v+rowStart*ncolumn, (rowStop-rowStart)*ncolumn*sizeof(Scalar));
#endif

}

template<class Scalar> 
void
GenFullM<Scalar>::transposeAssign(GenFullM<Scalar> &source)
{
 if(source.nrow != ncolumn || source.ncolumn != nrow) return; //error

 // unroll to avoid cache trashing

 int i,j;
 for(i = 0; i + 7 < ncolumn; i+= 8) {
   for(j = 0; j < nrow; ++j)  {
     Scalar *vv = v+j*ncolumn;
     Scalar *ww = source.v+j;
     vv[i] = ww[i*nrow]; 
     vv[i+1] = ww[(i+1)*nrow]; 
     vv[i+2] = ww[(i+2)*nrow]; 
     vv[i+3] = ww[(i+3)*nrow]; 
     vv[i+4] = ww[(i+4)*nrow]; 
     vv[i+5] = ww[(i+5)*nrow]; 
     vv[i+6] = ww[(i+6)*nrow]; 
     vv[i+7] = ww[(i+7)*nrow]; 
   }
 }
}

template<class Scalar> 
void
GenFullM<Scalar>::transposeMult(GenFullM<Scalar> &m, GenFullM<Scalar> &res)
{
 // res = this^T * m
 if(nrow != m.nrow) {
   fprintf(stderr," *** ERROR: incompatible dimensions in GenFullM<Scalar>::tranposeMult \n");
 }

 int i,j,k;
 for(i=0; i < ncolumn; ++i)
  for(j=0; j < m.ncolumn; ++j) {
    res[i][j] = 0.0 ;
    for(k = 0 ;  k < nrow ; ++k)
      res[i][j] += (*this)[k][i] * m[k][j];
   }
}

template<class Scalar>
void
GenFullM<Scalar>::transposeMultD(GenFullM<Scalar> &m, GenVector<Scalar> &d, GenFullM<Scalar> &res)
{
 // res = this^T * D * m
 if((nrow != m.nrow) || (nrow != d.size())) {
   fprintf(stderr," *** ERROR: incompatible dimensions in GenFullM<Scalar>::tranposeMultD \n");
 }

 int i,j,k;
 for(i=0; i < ncolumn; ++i)
  for(j=0; j < m.ncolumn; ++j) {
    res[i][j] = 0.0 ;
    for(k = 0 ;  k < nrow ; ++k)
      res[i][j] += (*this)[k][i] * d[k] * m[k][j];
   }
}

template<class Scalar>
Scalar*
GenFullM<Scalar>::Column(int col) const {

 Scalar *vector = new Scalar[nrow];
 int i;

 for (i=0; i<nrow; ++i)
   vector[i] = (*this)[i][col];

 return vector;

}

template class GenFullM<double>;
template class GenFullM<std::complex<double>>;
template class GenAssembledFullM<double>;
template class GenAssembledFullM<std::complex<double>>;
