//
// Created by Michel Lesoinne on 11/14/17.
//

#ifndef FEM_BLAS_H
#define FEM_BLAS_H
#include <complex>
#include <Utils.d/linkfc.h>
extern "C" {

//void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
//                     const int &, const double &, const double *, const int &,
//                     double *, const int &, const double &, double *, const int &);
void _FORTRAN(dgemm)(const char &trans, const char &transb,
                     const int &m, const int &n, const int &k,
                     const double &alpha, const double *a, const int &lda,
                     const double *b, const int &ldb,
                     const double &beta, double *c, const int &ldc);

void _FORTRAN(dgemv)(const char &, const int &,const int &,
                     const double &, const double *, const int &,
                     const double *, const int &, const double &, double *, const int &);

void _FORTRAN(zgemv)(const char &trans, const int &m, const int &n,
                     const ComplexD &lapha, const ComplexD *a, const int &lda,
                     const ComplexD *x, const int &incx,
                     const ComplexD &beta, ComplexD *y, const int &incy);

void _FORTRAN(zgemm)(const char &transa, const char &transb,
                     const int &m,const int &n, const int &k,
                     const ComplexD &alpha, const ComplexD *a, const int &lda,
                     const ComplexD *b, const int &ldb,
                     const ComplexD &beta, ComplexD *c, const int &ldc);

void _FORTRAN(dgesv)(const int &, const int &, double *, const int &, int *, double *,
                     const int &, int &);

void _FORTRAN(zgesv)(const int &, const int &, std::complex<double> *, const int &, int *, std::complex<double> *,
                     const int &, int &);

}
extern "C" {
void _FORTRAN(dsvdc)(double *, int &, int &, int&, double *,
                     double *, double *, int &, double *, int &,
                     double *, const int &, int &);

void _FORTRAN(zsvdc)(std::complex<double> *, int &, int &, int&, std::complex<double> *,
                     std::complex<double> *, std::complex<double> *, int &, std::complex<double> *, int &,
                     std::complex<double> *, const int &, int &);
}

inline void Tsvdc(double *a, int &b, int &c, int&d, double *e,
                  double *f, double *g, int &h, double *i, int &j,
                  double *k, const int &l, int &m)
{
  _FORTRAN(dsvdc)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

inline void Tsvdc(std::complex<double> *a, int &b, int &c, int&d, std::complex<double> *e,
                  std::complex<double> *f, std::complex<double> *g, int &h, std::complex<double> *i, int &j,
                  std::complex<double> *k, const int &l, int &m)
{
  _FORTRAN(zsvdc)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}
#ifndef _TGEMM__
#define _TGEMM__
inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const double &f, const double *g, const int &h,
                  const double *i, const int &j, const double &k, double *l,
                  const int &m)
{
  _FORTRAN(dgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const std::complex<double> &f, const std::complex<double> *g, const int &h,
                  const std::complex<double> *i, const int &j, const std::complex<double> &k, std::complex<double> *l,
                  const int &m)
{
  _FORTRAN(zgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}
#endif
#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, const double *e, const int &f,
                  const double *g, const int &h, const double &i, double *j, const int &k)
{
  _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
inline void Tgemv(const char &a, const int &b, const int &c,
                  const std::complex<double> &d, const std::complex<double> *e, const int &f,
                  const std::complex<double> *g, const int &h, const std::complex<double> &i, std::complex<double> *j, const int &k)
{
  _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif

#endif //FEM_BLAS_H
