#ifndef _COMPLEXD_H_
#define _COMPLEXD_H_
#define _STANDARD_C_PLUS_PLUS
#ifdef _STANDARD_C_PLUS_PLUS
#include <complex>
using std::complex;
#else
#include <complex.h>
#endif
#ifdef COMPLEX_NON_TEMPLATE
typedef complex ComplexD;
typedef complex DComplex;
#else
typedef complex<double> ComplexD;
typedef complex<double> DComplex;
#endif
#endif
