#ifndef _INTEGFUNCTION_H_
#define _INTEGFUNCTION_H_

#include <cmath>
#include <complex>
using std::complex;

class IntegFunctionL2d {
public:
 virtual void evaluate(double *x, double *N, double *cross,
                       double nsign, double w)=0;
};

class IntegFunctionA2d {
public:
 virtual void evaluate(double *x, double *N, double (*dNdx)[2],
                       double w, double det)=0;
};

class IntegFunctionL3d {
public:
 virtual void evaluate(double *x, double *N, double *cross,
                       double nsign, double *tau, double tsign, double w)=0;
};

class IntegFunctionA3d {
public:
 virtual void evaluate(double *x, double *N, double *cross,
                       double nsign, double w)=0;
};


class IntegFunctionAt3d {
public:
 virtual void evaluate(double *x, double *N, double *tau1, double *tau2,
                       double nsign, double w)=0;
};

class IntegFunctionAG3d {
public:
 virtual void evaluate(double *x, double *N, double (*dNdx)[3], double *cross,
                       double nsign, double w)=0;
};

class IntegFunctionAC3d {
public:
 virtual void evaluate(double *x, double *N, double *cross,
                       double nsign, double (*sd)[3], double (*surfgrad)[2],
                       double w)=0;
};

class IntegFunctionV3d {
public:
 virtual void evaluate(double *x, double *N, double (*dNdx)[3],
                       double w, double det)=0;
};

#endif
