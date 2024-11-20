#ifndef _MYCOMPLEX_H_
#define _MYCOMPLEX_H_
#include "Math.d/ComplexD.h"

namespace ScalarTypes 
{
 inline double conj(double a) { return a; }
 inline complex<double> conj(complex<double> a) { return std::conj(a); }

 inline double sqNorm(double a) { return a*a; }
 inline double sqNorm(complex<double> a)
    { return a.real()*a.real()+a.imag()*a.imag(); }

 inline double norm(double a) { return fabs(a); }
 inline double norm(complex<double> a) { return ::sqrt(sqNorm(a)); }

 // inline double real(double a) { return a; }
 // inline double real(complex<double> a) { return a.real(); }

 inline complex<double> sqrt(complex<double> a) { 
   double r  = ::sqrt(std::abs(a));
   double th = arg(a)/2.;
   return complex<double>(r*cos(th), r*sin(th));
 }
 inline double sqrt(double a) { return ::sqrt(a); }

 inline double Real(double a) { return a; }
 inline double Real(complex<double> a) { return a.real(); }

 inline double Imag(double a) { return 0.0; }
 inline double Imag(complex<double> a) { return a.imag(); }

 inline void initScalar(double &s, const double r, const double i = 0.0) { s = r; } 
 inline void initScalar(complex<double> &s, const double r, const double i = 0.0) 
   { s = complex<double>(r, i); }

 inline void addScalar(double &s, const double r, const double i = 0.0) { s += r; } 
 inline void addScalar(complex<double> &s, const double r, const double i = 0.0)
   { s += complex<double>(r, i); }

 inline void addReal(double &s, double _s) { s += _s; }
 inline void addComplex(double &s, complex<double> _s) { s += _s.real(); }
 inline void addReal(complex<double> &s, double _s) { addScalar(s, _s); }
 inline void addComplex(complex<double> &s, complex<double> _s) { s += _s; }

 inline void copy(double &s, complex<double> _s) { s = _s.real(); }
 inline void copy(complex<double> &s, complex<double> _s) { s = _s; }

 inline double doublify(double d) { return d; }
 inline double doublify(DComplex c) { return ::ScalarTypes::norm(c); }

 inline double d_doublify(double d, double g) { return g; }
 inline DComplex cd_doublify(DComplex c, DComplex g)
   {
     double n = ::ScalarTypes::norm(c);
     return (n==0)? 0 : (::ScalarTypes::conj(c)*g)/n;
   }
 inline double d_doublify(DComplex c, DComplex g) 
   { 
     return ::ScalarTypes::Real(cd_doublify(c,g));
   }

 template<class Scalar> Scalar id();
 template<> inline double id<double>() { return 1.0; }
 template<> inline DComplex id<DComplex>() { return DComplex(1.0, 1.0); }

 template<class Scalar> Scalar unify(Scalar);
 template<> inline double unify<double>(double) { return 1.0; }
 template<> inline DComplex unify<DComplex>(DComplex c) 
   { 
     double n = ::ScalarTypes::norm(c); 
     return n>0? c/n : 0.0;
   }

 inline bool greaterThan(double s1, double s2) { return (s1 > s2); }
 inline bool greaterThan(complex<double> s1, complex<double> s2) 
  { return ((s1.real() > s2.real() && s1.imag() >= s2.imag()) || (s1.real() >= s2.real() && s1.imag() > s2.imag())); }
 inline bool greaterThanEq(double s1, double s2) { return (s1 >= s2); }
 inline bool greaterThanEq(complex<double> s1, complex<double> s2)
  { return (s1.real() >= s2.real() && s1.imag() >= s2.imag()); }

 inline bool lessThan(double s1, double s2) { return (s1 < s2); }
 inline bool lessThan(complex<double> s1, complex<double> s2)
  { return ((s1.real() < s2.real() && s1.imag() <= s2.imag()) || (s1.real() <= s2.real() && s1.imag() < s2.imag())); }
 inline bool lessThanEq(double s1, double s2) { return (s1 <= s2); }
 inline bool lessThanEq(complex<double> s1, complex<double> s2)
  { return (s1.real() <= s2.real() && s1.imag() <= s2.imag()); }
}

#include <Utils.d/linkfc.h>

extern "C"      {
   void _FORTRAN(dsyev)(const char &, const char &, const int &, 
                        double *, const int &, double *, double *, 
                        const int &, int &);

   void _FORTRAN(zgeev)(const char &, const char &, const int &, DComplex *,
                        const int &, DComplex *, DComplex *, const int &,
                        DComplex *, const int &, DComplex *,
                        const int &, double *, int &);
}

inline void Tdsyev(const char &a, const char &b, const int &c,
                   double *d, const int &e, double *f, double *g,
                   const int &h, int &i)
{ _FORTRAN(dsyev)(a,b,c,d,e,f,g,h,i); }

inline void Tdsyev(const char &a, const char &b, const int &c,
                   DComplex *d, const int &e, DComplex *f, DComplex *g,
                   const int &h, int &i)
{ 
  double rwork[6];
  _FORTRAN(zgeev)(a,a,c,d,e,f,NULL,1,NULL,1,g,h,rwork,i);
}
                   
#endif
