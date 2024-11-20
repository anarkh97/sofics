#ifndef _ISOPARAMUTILS2D_H_
#define _ISOPARAMUTILS2D_H_

#include <cmath>
#include <complex>
#include <Element.d/Helm.d/IntegFunction.h>

using std::complex;

class IsoParamUtils2d {
public:
  int order;

  IsoParamUtils2d(int _o) { order = _o; }


  complex<double> expdiff(complex<double> alpha,
            complex<double> eaa, complex<double> eab, double L);
  double sidelen(double x, double y);

  virtual int getordersq();

  template<class Scalar> void zeroOut(int n,Scalar *v);
  template<class Scalar> void symmetrize(int n, Scalar *v);
  template<class Scalar> void copyOut(int m,int n, Scalar *K,
                                      double *Kr, double *Ki);

  virtual void lagGalShapeFunction(int ng, double *xi, double *N, 
                         int secondDerivsFlag = 0);
  virtual void spectralLagGalShapeFunction(int ng, double *xi, double *N, 
                         int secondDerivsFlag = 0);

  virtual void jmatrix(double *N, double *xyz, double (*J)[2]);
  double detj(double (*J)[2]);
  void invj(double (*J)[2], double det, double (*j)[2]);
  virtual void crossj(double (*J)[2], int axis, double *cross);

  virtual void lineInt2d(double *xy, int faceindex, IntegFunctionL2d &f, int gorder = 7);
  virtual void lineLineInt2d(double *xy, IntegFunctionL2d &f, int gorder = 7);

  virtual void areaInt2d(double *xy, IntegFunctionA2d &f, int gorder = 7);
  virtual void spectralAreaInt2d(double *xy, IntegFunctionA2d &f);

  virtual void elementcenter(double *xyz, double *cxyz);
  virtual void sidecenter(double *xyz, int faceindex, double *scxyz);
  virtual void faceindeces(int faceindex, int *fi);
  virtual void facemap(int &faceindex, int *fn, int *n, int *map);

  virtual int isStraight(double *xyz, int faceindex);
};

class IsoParamUtils2dTri: public IsoParamUtils2d {
public:
  IsoParamUtils2dTri(int _o): IsoParamUtils2d(_o) { ; }

  virtual int getordersq();
  virtual void lagGalShapeFunction(int mo, int ng, double *xi,
                              double *N, int secondDerivsFlag=0);

  virtual void jmatrix(double *N, double *xyz, double (*J)[2]);
  virtual void crossj(double (*J)[2], int faceindex, double *cross);

  virtual void lineInt2d(double *xy, int faceindex, IntegFunctionL2d &f, int gorder = 7);
  virtual void areaInt2d(double *xy, IntegFunctionA2d &f, int gorder = 7*7);

  virtual void elementcenter(double *xyz, double *cxyz);
  virtual void sidecenter(double *xyz, int faceindex, double *scxyz);
  virtual void faceindeces(int faceindex, int *fi);
  virtual void facemap(int &faceindex, int *fn, int *n, int *map);

  virtual int isStraight(double *xyz, int faceindex);
};


class GalFunction2d : public IntegFunctionA2d {
 int n;
 double kappa;
 double *K;
public:
 GalFunction2d(int _n, double _kappa, double *_K) {
   n = _n; kappa = _kappa; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {
   int i,j;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += w*det*(-kappa*kappa*N[i]*N[j]
        +dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]);
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


template<class Scalar> class GalStiffFunction2d : public IntegFunctionA2d {
 int n;
 Scalar *K;
public:
 GalStiffFunction2d(int _n, Scalar*_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {
   int i,j;
   double wdet = w*det;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += wdet*(
        dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]);
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


template<class Scalar> class GalMassFunction2d: public IntegFunctionA2d {
 int n;
 Scalar *K;
public:
 GalMassFunction2d(int _n, Scalar*_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {
   int i,j;
   double wdet = w*det;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += wdet*N[i]*N[j];
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


class NeumannBCGalFunction2d : public IntegFunctionL2d {
 double *incdir;
 double kappa;
 int n;
 int k; // k>0 for derivative of function ie: d^kf/dkappa^k
 complex<double>* v;
public:
 NeumannBCGalFunction2d(int _n, double _kappa, double *_incdir,
                        complex<double>* _v , int _k = 0) {
   n = _n; kappa = _kappa; incdir = _incdir; v = _v; k = _k;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   complex<double> e;
   complex<double> a = complex<double>(0.0,
                            (incdir[0]*x[0]+incdir[1]*x[1]));
   complex<double> b = complex<double>(0.0,(nsign*
            (cross[0]*incdir[0]+cross[1]*incdir[1])));
   if (k==0)
     e = b * kappa * exp(a*kappa);
   else 
     e = b * pow(a,k-1) * (double(k)+a*kappa) * exp(a*kappa);
   for(int i=0;i<n;i++)
       v[i] += w* N[i]* e;
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class RobinBCGalFunction2d : public IntegFunctionL2d {
 double *incdir;
 double kappa;
 int n;
 complex<double>* v;
public:
public:
 RobinBCGalFunction2d(int _n, double _kappa, double *_incdir,
                      complex<double>* _v) {
   n = _n; kappa = _kappa; incdir = _incdir; v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w,
                          complex<double> *K) {
   int i,j;
   for(j=0;j<ncolumns();j++)
     for(i=0;i<nrows();i++)
       K[j*nrows()+i] += w*
         N[i]*complex<double>(0.0,kappa)*
         (nsign*(cross[0]*incdir[0]+cross[1]*incdir[1])-
          sqrt(cross[0]*cross[0]+cross[1]*cross[1]))*
         exp(complex<double>(0.0,kappa*
                              (incdir[0]*x[0]+incdir[1]*x[1])));

 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class SomGalFunction2d : public IntegFunctionL2d {
 int n;
 double *K;
public:
 SomGalFunction2d(int _n, double *_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   int i,j;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += w*
         sqrt(cross[0]*cross[0]+cross[1]*cross[1])*N[i]*N[j];
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


class WetInterfaceGalFunction2d : public IntegFunctionL2d {
 int n;
 double *K;
public:
 WetInterfaceGalFunction2d(int _n, double *_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int i,j;
   for(j=0;j<n;j++)
     for(i=0;i<n;i++) {
       K[j*n+i] += w* cross[0]*N[i]*N[j];
       K[n*n+j*n+i] += w* cross[1]*N[i]*N[j];
     }
 }
 int nrows() { return n; }
 int ncolumns() { return 2*n; }
};


class WetInterfaceBCGalFunction2d : public IntegFunctionL2d {
 int n;
 double kappa;
 double *dir;
 int k;
 complex<double> *K;
public:
 WetInterfaceBCGalFunction2d(int _n, double _kappa, double *_dir,
                           complex<double> *_KK, int _k = 0) {
   n = _n; kappa = _kappa; dir = _dir; K = _KK; k = _k;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> idx = complex<double>(0.0,
                              (x[0]*dir[0]+x[1]*dir[1]));
   complex<double> edx = exp(kappa*idx);
   double cdir = cross[0]*dir[0]+cross[1]*dir[1];
   complex<double> wgen =
     w*edx*complex<double>(0,cdir)* pow(idx,k-1)*(double(k)+idx*kappa);
   complex<double> we = w*edx*pow(idx,k);
   for(int i=0;i<n;i++) {
     K[3*i+0] += we*N[i]*cross[0];
     K[3*i+1] += we*N[i]*cross[1];
     K[3*i+2] -= wgen*N[i];
   }
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class LEStiffFunction2d : public IntegFunctionA2d {
 int n;
 double lambda;
 double mu;
 double *K;
public:
 LEStiffFunction2d(int _n, double E, double nu, double *_K) {
   n = _n; K = _K;
   lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
   mu = E/(2.0*(1.0+nu));
 }
 double elasticForm(double lambda, double mu,
          double (*gradu)[2], double (*gradv)[2], double *u, double *v) {
  double epsu[2][2], epsv[2][2];
  epsu[0][0] = gradu[0][0];
  epsu[1][1] = gradu[1][1];
  epsu[0][1] = 0.5*(gradu[0][1]+gradu[1][0]);
  epsu[1][0] = epsu[0][1];
  epsv[0][0] = gradv[0][0];
  epsv[1][1] = gradv[1][1];
  epsv[0][1] = 0.5*(gradv[0][1]+gradv[1][0]);
  epsv[1][0] = epsv[0][1];
 
 
  return
    2.0*mu*(epsu[0][0]*epsv[0][0]+epsu[0][1]*epsv[0][1]+
            epsu[1][0]*epsv[1][0]+epsu[1][1]*epsv[1][1]) +
    lambda*(gradu[0][0]+gradu[1][1])*(gradv[0][0]+gradv[1][1]);
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det){
   int nn = n/2;
   double wdet = w*det;
   for(int j=0;j<nn;j++) for(int i=j;i<nn;i++) {

     double gradu[2][2] = {{dNdx[i][0], dNdx[i][1]},{0.0,0.0}} ;
     double u[2] = {N[i],0.0};

     double gradv[2][2] = {{dNdx[j][0], dNdx[j][1]},{0.0,0.0}} ;
     double v[2] = {N[j],0.0};

     double gradw[2][2] = {{0.0,0.0},{dNdx[i][0], dNdx[i][1]}} ;
     double w[2] = {0.0,N[i]};

     double gradx[2][2] = {{0.0,0.0},{dNdx[j][0], dNdx[j][1]}} ;
     double x[2] = {0.0,N[j]};

     K[2*j*n+2*i] += wdet*elasticForm(lambda, mu, gradu, gradv, u, v);
     K[2*j*n+2*i+1] += wdet*elasticForm(lambda, mu, gradw, gradv, w, v);
     if (i>j)
       K[(2*j+1)*n+2*i] += wdet*elasticForm(lambda, mu, gradu, gradx, u, x);
     K[(2*j+1)*n+2*i+1] += wdet*elasticForm(lambda, mu, gradw, gradx, w, x);
   }
 }
};


class LEMassFunction2d : public IntegFunctionA2d {
 int n;
 double rho;
 double *K;
public:
 LEMassFunction2d(int _n, double _rho, double *_K) {
   n = _n; rho = _rho; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det){
   int nn = n/2;
   double wdetrho = w*det*rho;
   for(int j=0;j<nn;j++) for(int i=j;i<nn;i++) {

     K[2*j*n+2*i] += wdetrho*N[i]*N[j];
     K[(2*j+1)*n+2*i+1] += wdetrho*N[i]*N[j];
   }
 }
};

class AreaFunction2d : public IntegFunctionA2d {
 double *area;
public:
 AreaFunction2d(double *_area) { area = _area; }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) { *area += w*det; }
};


#ifdef _TEMPLATE_FIX_
#include <Element.d/Helm.d/IsoParamUtils2dT.C>
#endif

#endif
