#ifndef _ISOPARAMUTILS_H_
#define _ISOPARAMUTILS_H_

#include <cmath>
#include <complex>
#include <Element.d/Helm.d/IntegFunction.h>

using std::complex;

class IsoParamUtils {
public:
  int order;

  IsoParamUtils(int _o) { order = _o; }

  double sidelen(double x, double y, double z);

  virtual int getorderc();
  virtual int getordersq(int faceindex=0);

  template<class Scalar> void zeroOut(int n,Scalar *v);
  template<class Scalar> void symmetrize(int n, Scalar *v);

  virtual void lagGalShapeFunction(int ng, double *xi, double *N, 
                         int secondDerivsFlag = 0);
  virtual void spectralLagGalShapeFunction(int ng, double *xi, double *N, 
                         int secondDerivsFlag = 0);

  void crossproduct(double *af, double *bf, double *nf);
  virtual void jmatrix(double *N, double *xyz, double (*J)[3]);
  double detj(double (*J)[3]);
  void invj(double (*J)[3], double det, double (*j)[3]);
  virtual void crossj(double (*J)[3], int axis, double *cross);
  virtual void crossj(double (*J)[3], int axis, double *cross, 
                      double *tau1, double *tau2);

  virtual void getSurfDer(double *N, double *N2, double *xyz,
                       int axis, double (*sd)[3], double (*surfgrad)[2]);

  virtual void lineInt3d(double *xyz, int faceindex, int edgeindex,
               double *nf, IntegFunctionL3d &f, int gorder = 7);

  virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f, int gorder = 7);
  virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f, int gorder = 7);
  virtual void surfGInt3d(double *xyz, int faceindex, IntegFunctionAG3d &f, int gorder = 7);
  virtual void surfSurfInt3d(double *xyz, IntegFunctionA3d &f, int gorder =7);
  virtual void spectralSurfSurfInt3d(double *xyz, IntegFunctionA3d &f);

  virtual void surfCurvInt3d(double *xyz, int faceindex, IntegFunctionAC3d &f, int gorder = 7);

  virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f, int gorder = 7);
  virtual void spectralVolumeInt3d(double *xyz, IntegFunctionV3d &f);

  virtual void elementcenter(double *xyz, double *cxyz);
  virtual void cornerindeces(int faceindex, int *ci);
  virtual void sidecenter(double *xyz, int faceindex, double *scxyz);
  virtual void faceindeces(int faceindex, int *fi);
  virtual void facemap(int &faceindex, int *fn, int *n, int *map);
  virtual void adjacentfaces(int faceindex, int *baf);

  virtual int isFlat(double *xyz, int faceindex);
  virtual int isStraight(double *xyz, int faceindex, int edgeindex);

};


class IsoParamUtilsTetra: public IsoParamUtils {
public:
  IsoParamUtilsTetra(int _o): IsoParamUtils(_o) { ; }

  virtual int getorderc();
  virtual int getordersq(int faceindex=0);

  virtual void lagGalShapeFunction(int mo, int ng, double *xi,
                              double *N, int secondDerivsFlag = 0);

  virtual void crossj(double (*J)[3], int faceindex, double *cross);

  virtual void getSurfDer(double *N, double *N2, double *xyz,
                       int faceindex, double (*sd)[3], double (*surfgrad)[2]);

  virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f, int gorder = 7*7);
  virtual void surfGInt3d(double *xyz, int faceindex, IntegFunctionAG3d &f, int gorder = 7*7);
  virtual void surfSurfInt3d(double *xyz, IntegFunctionA3d &f, int gorder = 7*7);
  virtual void surfCurvInt3d(double *xyz, int faceindex, IntegFunctionAC3d &f, int gorder = 7*7);

  virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f, int gorder = 7*7*7);

  virtual void elementcenter(double *xyz, double *cxyz);
  virtual void cornerindeces(int faceindex, int *ci);
  virtual void sidecenter(double *xyz, int faceindex, double *scxyz);
  virtual void faceindeces(int faceindex, int *fi);
  virtual void facemap(int &faceindex, int *fn, int *n, int *map);
};


class IsoParamUtilsPrism: public IsoParamUtils {
public:
  IsoParamUtilsPrism(int _o): IsoParamUtils(_o) { ; }

  virtual int getorderc();
  virtual int getordersq(int faceindex=0);

  virtual void lagGalShapeFunctionTr(int mo, int ng, double *xi,
                              double *N, int secondDerivsFlag = 0);

  virtual void crossj(double (*J)[3], int faceindex, double *cross);

  virtual void getSurfDer(double *N, double *N2, double *xyz,
                       int faceindex, double (*sd)[3], double (*surfgrad)[2]) {
    fprintf(stderr,"IsoParamUtilsPrism::getSurfDer is not implemented.\n");
  }
  virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f, int gorder = 7);
  virtual void surfGInt3d(double *xyz, int faceindex, IntegFunctionAG3d &f,
                          int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPrism::surfGInt3d is not implemented.\n");
  }
  virtual void surfSurfInt3d(double *xyz, IntegFunctionA3d &f,
                             int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPrism::surfSurfInt3d is not implemented.\n");
  }
  virtual void surfCurvInt3d(double *xyz, int faceindex, IntegFunctionAC3d &f,
                             int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPrism::surfCurvInt3d is not implemented.\n");
  }

  virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f, int gorder = 7);

  virtual void elementcenter(double *xyz, double *cxyz);
  virtual void cornerindeces(int faceindex, int *ci);
  virtual void sidecenter(double *xyz, int faceindex, double *scxyz);
  virtual void faceindeces(int faceindex, int *fi);
  virtual void facemap(int &faceindex, int *fn, int *n, int *map) {
    fprintf(stderr,"IsoParamUtilsPrism::facemap is not implemented.\n");
  }
};


class IsoParamUtilsPyramid: public IsoParamUtilsPrism {
public:
  IsoParamUtilsPyramid(int _o): IsoParamUtilsPrism(_o) { ; }

  virtual int getorderc() { return 5; }
  virtual int getordersq(int faceindex=0) { return (faceindex<=1)?4:3; }

  virtual void crossj(double (*J)[3], int faceindex, double *cross) {
    fprintf(stderr,"IsoParamUtilsPyramid::crossj is not implemented.\n");
  }

  virtual void getSurfDer(double *N, double *N2, double *xyz,
                       int faceindex, double (*sd)[3], double (*surfgrad)[2]) {
    fprintf(stderr,"IsoParamUtilsPyramid::getSurfDer is not implemented.\n");
  }
  virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f, int gorder = 7);
  virtual void surfGInt3d(double *xyz, int faceindex, IntegFunctionAG3d &f,
                          int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPyramid::surfGInt3d is not implemented.\n");
  }
  virtual void surfSurfInt3d(double *xyz, IntegFunctionA3d &f,
                             int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPyramid::surfSurfInt3d is not implemented.\n");
  }
  virtual void surfCurvInt3d(double *xyz, int faceindex, IntegFunctionAC3d &f,
                             int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPyramid::surfCurvInt3d is not implemented.\n");
  }
  virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f, int gorder = 7) {
    fprintf(stderr,"IsoParamUtilsPyramid::volumeInt3d is not implemented.\n");
  }
  virtual void elementcenter(double *xyz, double *cxyz);
  virtual void cornerindeces(int faceindex, int *ci);
  virtual void sidecenter(double *xyz, int faceindex, double *scxyz);
  virtual void faceindeces(int faceindex, int *fi);
  virtual void facemap(int &faceindex, int *fn, int *n, int *map) {
    fprintf(stderr,"IsoParamUtilsPyramid::facemap is not implemented.\n");
  }
};


class GalFunction : public IntegFunctionV3d {
 int n;
 double kappa;
 double *K;
public:
 GalFunction(int _n, double _kappa, double *_K) {
   n = _n; kappa = _kappa; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
   for(int j=0;j<ncolumns();j++)
     for(int i=j;i<nrows();i++)
       K[j*nrows()+i] += w*det*(-kappa*kappa*N[i]*N[j]
        +dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]+dNdx[i][2]*dNdx[j][2]);
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


class NeumannBCGalFunction : public IntegFunctionA3d {
 double *incdir;
 double kappa;
 int pflag;
 int n;
 int k;  // PJSA: k>0 for derivative of function ie: d^kf/dkappa^k
 complex<double> *v;
public:
 NeumannBCGalFunction(int _n, double _kappa, double *_incdir, int _pflag,
  complex<double>* _v, int _k = 0) {
   n = _n; kappa = _kappa; incdir = _incdir; pflag = _pflag; k = _k; v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   complex<double> e;
   if (pflag==0) {
     double aa = incdir[0]*x[0]+incdir[1]*x[1]+incdir[2]*x[2];
     complex<double> a = complex<double>(0.0,aa);
     complex<double> b = complex<double>(0.0,(nsign*
              (cross[0]*incdir[0]+cross[1]*incdir[1]+cross[2]*incdir[2])));
     if (k==0) 
       e = b * kappa * exp(a*kappa);
     else if (aa==0) {
       if (k==1) e = b;
       else e = 0;
     }
     else  {
       e = b * pow(a,k-1) * (double(k)+a*kappa) * exp(a*kappa);
     }
   } else {
     double r = sqrt((x[0]-incdir[0])*(x[0]-incdir[0])+
                     (x[1]-incdir[1])*(x[1]-incdir[1])+
                     (x[2]-incdir[2])*(x[2]-incdir[2]));
     double nx = cross[0]*(x[0]-incdir[0])+
                 cross[1]*(x[1]-incdir[1])+
                 cross[2]*(x[2]-incdir[2]);
     complex<double> e1 =
                 nsign*nx*exp(complex<double>(0.0,kappa*r))/(r*r*r*4.0*M_PI);
     complex<double> a = complex<double>(0.0,r);
     if (k==0) e = e1 * (a*kappa-1.0);
     else e = e1 * ( pow(a,k)*(double(k-1) +a*kappa) );
   }
     for(int i=0;i<n;i++)
       v[i] += w* N[i]* e;
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class RobinBCGalFunction {
 double *incdir;
 double kappa;
 int n;
 complex<double>* v;
public:
public:
 RobinBCGalFunction(int _n, double _kappa, double *_incdir, complex<double>* _v) {
   n = _n; kappa = _kappa; incdir = _incdir; v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w,
                          complex<double> *K) {
   int i,j;
   for(j=0;j<ncolumns();j++)
     for(i=0;i<nrows();i++)
       K[j*nrows()+i] += w*
         N[i]*complex<double>(0.0,kappa)*
         (nsign*(cross[0]*incdir[0]+cross[1]*incdir[1]+cross[2]*incdir[2])-
          sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]))*
         exp(complex<double>(0.0,kappa*
                              (incdir[0]*x[0]+incdir[1]*x[1]+incdir[2]*x[2])));

 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class SomGalFunction : public IntegFunctionA3d {
 int n;
 double kappa;
 double *K;
public:
 SomGalFunction(int _n, double *_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   int i,j;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += w*
         sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2])*N[i]*N[j];
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


class WetInterfaceGalFunction : public IntegFunctionA3d {
 int n;
 double kappa;
 double *K;
public:
 WetInterfaceGalFunction(int _n, double *_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   for(int j=0;j<n;j++) 
     for(int i=0;i<n;i++) {
       K[j*n+i] += w* cross[0]*N[i]*N[j];
       K[n*n+j*n+i] += w* cross[1]*N[i]*N[j];
       K[2*n*n+j*n+i] += w* cross[2]*N[i]*N[j]; 
     }
 }
 int nrows() { return n; }
 int ncolumns() { return 3*n; }
};


class WetInterfaceBCGalFunction : public IntegFunctionA3d {
 int n;
 double kappa;
 double *dir;
 int pflag;
 complex<double> *K;
 int k;  // k>0 for derivative of function ie: d^kf/dkappa^k
public:
 WetInterfaceBCGalFunction(int _n, double _kappa, double *_dir, int _pflag,
                           complex<double> *_KK, int _k=0) {
   n = _n; kappa = _kappa; dir = _dir; pflag = _pflag; K = _KK; k = _k; 
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> we, wgen;
   if (pflag == 0) {
     complex<double> idx = complex<double>(0.0,
                                (x[0]*dir[0]+x[1]*dir[1]+x[2]*dir[2])); 
     complex<double> edx = exp(kappa*idx); 
     double cdir = cross[0]*dir[0]+cross[1]*dir[1]+cross[2]*dir[2];
     wgen = w*edx*complex<double>(0,cdir)* pow(idx,k-1)*(double(k)+idx*kappa);
     we = w*edx*pow(idx,k);
   } else {
     double r = sqrt((x[0]-dir[0])*(x[0]-dir[0])+
                     (x[1]-dir[1])*(x[1]-dir[1])+
                     (x[2]-dir[2])*(x[2]-dir[2]));
     double nx = cross[0]*(x[0]-dir[0])+
                 cross[1]*(x[1]-dir[1])+
                 cross[2]*(x[2]-dir[2]);
     complex<double> e1 = nx*exp(complex<double>(0.0,kappa*r))/(r*r*r*4.0*M_PI);
     complex<double> a = complex<double>(0.0,r);
     if (k==0) wgen = w * e1 * (a*kappa-1.0);
     else wgen = w * e1 * ( pow(a,k)*(double(k-1) +a*kappa) );
     we = w*exp(complex<double>(0.0,kappa*r))/(r*4.0*M_PI)*pow(a,k);
   }
   for(int i=0;i<n;i++) {
     K[4*i+0] += we*N[i]*cross[0];
     K[4*i+1] += we*N[i]*cross[1];
     K[4*i+2] += we*N[i]*cross[2];
     K[4*i+3] -= wgen*N[i];
   }
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class WetInterfaceHeteroBCGalFunction : public IntegFunctionA3d {
 int n;
 int imat;
 complex<double> (*diri)[3];
 complex<double> *coefi;
 complex<double> *K;
public:
 WetInterfaceHeteroBCGalFunction(int _n, int _imat,
                       complex<double> (*_diri)[3], complex<double> *_coefi,
                           complex<double> *_KK) {
  n = _n; imat = _imat; diri = _diri; coefi = _coefi; K = _KK;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   complex<double> we, wgen;
   if (imat==1) {
     complex<double> e0 =
        exp(diri[0][0]*x[0]+diri[0][1]*x[1]+diri[0][2]*x[2]);
     complex<double> e1 =
        exp(diri[1][0]*x[0]+diri[1][1]*x[1]+diri[1][2]*x[2]);
     we = w * ( coefi[0] * e0 + coefi[1] * e1);
     wgen = w*(
        (cross[0]*diri[0][0]+cross[1]*diri[0][1]+cross[2]*diri[0][2])*
           coefi[0]*e0 +
        (cross[0]*diri[1][0]+cross[1]*diri[1][1]+cross[2]*diri[1][2])*
           coefi[1]*e1);
   }
   else {
     complex<double> e2 =
        exp(diri[2][0]*x[0]+diri[2][1]*x[1]+diri[2][2]*x[2]);
     we = w * coefi[2]*e2;
     wgen = w * coefi[2]*e2 *
        (cross[0]*diri[2][0]+cross[1]*diri[2][1]+cross[2]*diri[2][2]);
   }
   int i;
   for(i=0;i<n;i++) {
     K[4*i+0] += we*N[i]*cross[0];
     K[4*i+1] += we*N[i]*cross[1];
     K[4*i+2] += we*N[i]*cross[2];
     K[4*i+3] -= wgen*N[i];
//     K[4*i+3] += wgen*N[i];
   }
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class WetInterfaceHeteroBCDerivGalFunction : public IntegFunctionA3d {
 int n;
 int imat;
 complex<double> (*diri)[3];
 complex<double> *coefi;
 complex<double> *kappai;
 int dero;
 complex<double> *K;
public:
 WetInterfaceHeteroBCDerivGalFunction(int _n, int _imat,
                       complex<double> (*_diri)[3], complex<double> *_coefi,
                       complex<double> *_kappa, int _dero,
                       complex<double> *_KK) {
  n = _n; imat = _imat; diri = _diri; coefi = _coefi; kappai = _kappa;
  dero = _dero; K = _KK;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   complex<double> we, wgen;
   if (imat==1) {
// RT: Not implemented yet
     we = 0.0;
     wgen = 0.0;
   }
   else {

// coefi[2] has cscale_factor and multiplies solid dofs only

     complex<double> e2 =
        exp(diri[2][0]*x[0]+diri[2][1]*x[1]+diri[2][2]*x[2]);
     complex<double> idx = (diri[2][0]*x[0]+diri[2][1]*x[1]+diri[2][2]*x[2]) /
                           kappai[2];
     complex<double> idn = (cross[0]*diri[2][0]+
                            cross[1]*diri[2][1]+
                            cross[2]*diri[2][2]) / kappai[2];

     we = w * coefi[2]*pow(idx,dero)*e2;
     wgen = w * e2 * idn *
            (double(dero)*pow(idx,dero-1) + kappai[2]*pow(idx,dero));
   }
   int i;
   for(i=0;i<n;i++) {
     K[4*i+0] += we*N[i]*cross[0];
     K[4*i+1] += we*N[i]*cross[1];
     K[4*i+2] += we*N[i]*cross[2];
     K[4*i+3] -= wgen*N[i];
   }
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class FFPGalFunction : public IntegFunctionAG3d {
 int n;
 double kappa;
 int nffp;
 double *ffpdir;
 complex<double> *sol;
 complex<double>* v;
public:
 FFPGalFunction(int _n, double _kappa, int _nffp, double *_ffpdir,
                complex<double> *_sol, complex<double>* _v) {
   n = _n; kappa = _kappa; nffp = _nffp; ffpdir = _ffpdir; sol = _sol; v = _v;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double *cross,
               double nsign, double w) {
   
   int i;
   complex<double> s(0.0,0.0);
   complex<double> gs[3] = {0.0,0.0,0.0};
   for(i=0;i<n;i++) {
    s += N[i]*sol[i];
    gs[0] +=  dNdx[i][0]*sol[i];
    gs[1] +=  dNdx[i][1]*sol[i];
    gs[2] +=  dNdx[i][2]*sol[i];
   }

   for(i=0;i<nffp;i++) {
     complex<double> we = (nsign*w/(4.0*M_PI)) * exp(complex<double>(0.0,-kappa*
           (ffpdir[3*i+0]*x[0]+ffpdir[3*i+1]*x[1]+ffpdir[3*i+2]*x[2])));

     complex<double> f = 
      (gs[0]+complex<double>(0.0,kappa*ffpdir[3*i+0])*s)*cross[0]+
      (gs[1]+complex<double>(0.0,kappa*ffpdir[3*i+1])*s)*cross[1]+
      (gs[2]+complex<double>(0.0,kappa*ffpdir[3*i+2])*s)*cross[2];
     
     v[i] += we*f;
   }
 }
};

class KirchhoffGalFunction : public IntegFunctionAG3d {
 int n;
 double kappa;
 int nffp;
 double *ffploc;
 complex<double> *sol;
 complex<double>* v;
public:
 KirchhoffGalFunction(int _n, double _kappa, int _nffp, double *_ffploc,
                      complex<double> *_sol, complex<double>* _v) {
   n = _n; kappa = _kappa; nffp = _nffp; ffploc = _ffploc; sol = _sol; v = _v;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double *cross,
               double nsign, double w) {

   int i;
   complex<double> s(0.0,0.0);
   complex<double> gs[3] = {0.0,0.0,0.0};
   for(i=0;i<n;i++) {
    s += N[i]*sol[i];
    gs[0] +=  dNdx[i][0]*sol[i];
    gs[1] +=  dNdx[i][1]*sol[i];
    gs[2] +=  dNdx[i][2]*sol[i];
   }

   for(i=0;i<nffp;i++) {
     double rv[3];
     rv[0] = x[0] - ffploc[3*i+0];
     rv[1] = x[1] - ffploc[3*i+1];
     rv[2] = x[2] - ffploc[3*i+2];

     double r = sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);

     double gr[3];
     gr[0] = rv[0]/r;
     gr[1] = rv[1]/r;
     gr[2] = rv[2]/r;

     complex<double> we = - (nsign*w/(4.0*M_PI*r*r)) * exp(complex<double>(0.0, kappa*r));

     complex<double> sikr = s * complex<double>(-1.0, kappa*r);

     complex<double> f =
      (r*gs[0] - sikr * gr[0])*cross[0]+
      (r*gs[1] - sikr * gr[1])*cross[1]+
      (r*gs[2] - sikr * gr[2])*cross[2];

     v[i] += we*f;
   }

 }
};

class LENeumInterfaceBCGalFunction: public IntegFunctionA3d {
 int n;
 double omega;
 double lambda;
 double mu;
 double rho;
 double *dir;
 complex<double> *K;
public:
 LENeumInterfaceBCGalFunction(int _n, double _omega, double E, double nu,
                            double _rho, double *_dir, complex<double> *_KK) {
   n = _n; omega= _omega; rho = _rho;  dir = _dir; K = _KK;
   lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
   mu = E/(2.0*(1.0+nu));
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int nn = n/3;

   double kp = omega*sqrt(rho)/sqrt(lambda+2.0*mu); 
   complex<double> we = w*
     exp(complex<double>(0.0,kp*(x[0]*dir[0]+x[1]*dir[1]+x[2]*dir[2])));
   complex<double> wike = we*complex<double>(0.0,kp);
   //double l = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
   double sun[3] = {lambda*cross[0], lambda*cross[1], lambda*cross[2]} ;
   int i;
   for(i=0;i<3;i++) {
    sun[i] += 2.0*mu*
      (dir[i]*dir[0]*cross[0]+dir[i]*dir[1]*cross[1]+dir[i]*dir[2]*cross[2]);
   }
   for(i=0;i<nn;i++) {
//     K[3*i+0] -= wike*sun[0]*N[i];
//    K[3*i+1] -= wike*sun[1]*N[i];
//     K[3*i+2] -= wike*sun[2]*N[i];

/*
if (x[2]>0.0) {
     K[3*i+0] += w*cross[0]*N[i];
     K[3*i+1] += w*cross[1]*N[i];
     K[3*i+2] += w*cross[2]*N[i];
}
*/
   double r = sqrt((x[0]-0.1524)*(x[0]-0.1524)+x[1]*x[1]+x[2]*x[2])/0.0015;
   double hat_;
   if (r<0.99)
     hat_ = 1.0/exp(-1.0/(r*r-1.0));
   else 
    hat_ = 0.0; 
   double hat = hat_;
//     K[3*i+0] += hat*w*cross[0]*N[i];
//     K[3*i+1] += hat*w*cross[1]*N[i];
//     K[3*i+2] += hat*w*cross[2]*N[i];
     K[3*i+0] += w*cross[0]*N[i];

   }
 }
 int nrows() { return n; }
 int ncolumns() { return 1; }
};


class Som2GalFunction : public IntegFunctionAC3d {
 int n;
 double kappa;
 complex<double> *Mat;
public:
 Som2GalFunction(int _n, double _kappa, complex<double> *_Mat) {
   n = _n; kappa = _kappa; Mat = _Mat;
 }
void evaluate(double *x, double *SN, double *cross, double nsign,
              double (*sd)[3], double (*surfgrad)[2], double w) {
   double l = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
   double E = sd[0][0]*sd[0][0]+sd[0][1]*sd[0][1]+sd[0][2]*sd[0][2];
   double F = sd[0][0]*sd[1][0]+sd[0][1]*sd[1][1]+sd[0][2]*sd[1][2];
   double G = sd[1][0]*sd[1][0]+sd[1][1]*sd[1][1]+sd[1][2]*sd[1][2];
   double L = nsign*(cross[0]*sd[2][0]+cross[1]*sd[2][1]+cross[2]*sd[2][2])/l;
   double M = nsign*(cross[0]*sd[4][0]+cross[1]*sd[4][1]+cross[2]*sd[4][2])/l;
   double N = nsign*(cross[0]*sd[3][0]+cross[1]*sd[3][1]+cross[2]*sd[3][2])/l;

   double H = fabs((G*L-2.0*F*M+E*N)/(2.0*(E*G-F*F)));
   double K = (L*N-M*M)/(E*G-F*F);

   complex<double> r;
   complex<double> c = complex<double>(H,-kappa) +
                       complex<double>(0.0,(K-H*H)/2.0)/
                       complex<double>(kappa,2.0*H);

   double dN[2][2] = { {(M*F-L*G)/(E*G-F*F),(N*F-M*G)/(E*G-F*F)},
                       {(L*F-M*E)/(E*G-F*F),(M*F-N*E)/(E*G-F*F)}};

   int i,j;
   for(j=0;j<ncolumns();j++) {
     for(i=j;i<nrows();i++) {

       Mat[j*nrows()+i] += w*l*c*SN[i]*SN[j];

       double surfgradcoef1[2] = { (surfgrad[i][0]*G-surfgrad[i][1]*F)/(E*G-F*F),
                                  (surfgrad[i][1]*E-surfgrad[i][0]*F)/(E*G-F*F)};
       double surfgradcoef2[2] = { (surfgrad[j][0]*G-surfgrad[j][1]*F)/(E*G-F*F),
                                  (surfgrad[j][1]*E-surfgrad[j][0]*F)/(E*G-F*F)};


       complex<double> MX[2][2] = {
          { complex<double>(kappa,dN[0][0]), complex<double>(0.0,dN[0][1]) },
          { complex<double>(0.0,dN[1][0]), complex<double>(kappa,dN[1][1]) }};

       complex<double> invMX[2][2] = { {MX[1][1],-MX[0][1]},{-MX[1][0],MX[0][0]}};
       invMX[0][0] /= MX[0][0]*MX[1][1]-MX[0][1]*MX[1][0];
       invMX[1][0] /= MX[0][0]*MX[1][1]-MX[0][1]*MX[1][0];
       invMX[0][1] /= MX[0][0]*MX[1][1]-MX[0][1]*MX[1][0];
       invMX[1][1] /= MX[0][0]*MX[1][1]-MX[0][1]*MX[1][0];

       complex<double> tmp[2] = {
        complex<double>(0.0,0.5) *
          (invMX[0][0]*surfgradcoef1[0]+invMX[0][1]*surfgradcoef1[1]),
        complex<double>(0.0,0.5) *
          (invMX[1][0]*surfgradcoef1[0]+invMX[1][1]*surfgradcoef1[1])};

       double gsurfgrad2[3] = {
         surfgradcoef2[0]*sd[0][0]+surfgradcoef2[1]*sd[1][0],
         surfgradcoef2[0]*sd[0][1]+surfgradcoef2[1]*sd[1][1],
         surfgradcoef2[0]*sd[0][2]+surfgradcoef2[1]*sd[1][2]
       };
       complex<double> gsurftmp[3] = {
         tmp[0]*sd[0][0]+tmp[1]*sd[1][0],
         tmp[0]*sd[0][1]+tmp[1]*sd[1][1],
         tmp[0]*sd[0][2]+tmp[1]*sd[1][2]
       };

       Mat[j*nrows()+i] += w*(gsurftmp[0]*gsurfgrad2[0]+
                       gsurftmp[1]*gsurfgrad2[1]+gsurftmp[2]*gsurfgrad2[2])*l;
     }
   }
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};

//#define LESLOW
#ifdef LESLOW
class LEStiffFunction {
 int n;
 double lambda;
 double mu;
 double *K;
public:
 LEStiffFunction(int _n, double E, double nu, double *_K) {
   n = _n; K = _K;
   lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
   mu = E/(2.0*(1.0+nu));

 }
 double elasticForm(double lambda, double mu,
          double (*gradu)[3], double (*gradv)[3], double *u, double *v) {
  double epsu[3][3], epsv[3][3];
  epsu[0][0] = gradu[0][0];
  epsu[1][1] = gradu[1][1];
  epsu[2][2] = gradu[2][2];
  epsu[0][1] = 0.5*(gradu[0][1]+gradu[1][0]);
  epsu[0][2] = 0.5*(gradu[0][2]+gradu[2][0]);
  epsu[1][2] = 0.5*(gradu[1][2]+gradu[2][1]);
  epsu[1][0] = epsu[0][1];
  epsu[2][0] = epsu[0][2];
  epsu[2][1] = epsu[1][2];
  epsv[0][0] = gradv[0][0];
  epsv[1][1] = gradv[1][1];
  epsv[2][2] = gradv[2][2];
  epsv[0][1] = 0.5*(gradv[0][1]+gradv[1][0]);
  epsv[0][2] = 0.5*(gradv[0][2]+gradv[2][0]);
  epsv[1][2] = 0.5*(gradv[1][2]+gradv[2][1]);
  epsv[1][0] = epsv[0][1];
  epsv[2][0] = epsv[0][2];
  epsv[2][1] = epsv[1][2];


  return
    2.0*mu*(
     epsu[0][0]*epsv[0][0]+epsu[0][1]*epsv[0][1]+epsu[0][2]*epsv[0][2]+
     epsu[1][0]*epsv[1][0]+epsu[1][1]*epsv[1][1]+epsu[1][2]*epsv[1][2]+
     epsu[2][0]*epsv[2][0]+epsu[2][1]*epsv[2][1]+epsu[2][2]*epsv[2][2]) +
    lambda*(gradu[0][0]+gradu[1][1]+gradu[2][2])*
           (gradv[0][0]+gradv[1][1]+gradv[2][2]);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det){
   int nn = n/3;
   double wdet = w*det;
   for(int j=0;j<nn;j++) for(int i=j;i<nn;i++) 
    for(int a=0;a<3;a++) for(int b=0; b<3; b++) {

     double gradu[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
     double u[3]={0.0,0.0,0.0};
     gradu[a][0] = dNdx[i][0];
     gradu[a][1] = dNdx[i][1];
     gradu[a][2] = dNdx[i][2];
     u[a] = N[i];

     double gradv[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
     double v[3]={0.0,0.0,0.0};
     gradv[b][0] = dNdx[j][0];
     gradv[b][1] = dNdx[j][1];
     gradv[b][2] = dNdx[j][2];
     u[b] = N[j];
     if ((i!=j) || ((i==j) && (a>=b))) {
       K[(3*j+b)*n+3*i+a] += wdet*elasticForm(lambda, mu, gradu, gradv, u, v);
    }
 }
};
#else
class LEStiffFunction : public IntegFunctionV3d {
 int n;
 double lambda;
 double mu;
 double *K;
public:
 LEStiffFunction(int _n, double E, double nu, double *_K) {
   n = _n; K = _K;
   lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
   mu = E/(2.0*(1.0+nu));

 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det){
   int nn = n/3;
   double wdet = w*det;
   double wdet2mu = wdet*2.0*mu;
   double wdetlambda = wdet*lambda;
   int i;
   double *e = (double*)alloca(sizeof(double)*3*nn);
   for(i=0;i<nn;i++) { 
     e[3*i+0] = dNdx[i][0];
     e[3*i+1] = dNdx[i][1];
     e[3*i+2] = dNdx[i][2];
   }
   
   int j; 
   for(j=0;j<nn;j++) for(i=j+1;i<nn;i++) {
     int a,b;
     a = 0; b = 0;
     K[(3*j+b)*n+3*i+a] += wdet2mu*(
      e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1]+0.5*e[3*i+2]*e[3*j+2] ) +
                           wdetlambda*e[3*i+0]*e[3*j+0];
     a = 1; b = 1;
     K[(3*j+b)*n+3*i+a] += wdet2mu*(
      e[3*i+1]*e[3*j+1]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+2]*e[3*j+2] ) +
                           wdetlambda*e[3*i+1]*e[3*j+1];
     a = 2; b = 2;
     K[(3*j+b)*n+3*i+a] += wdet2mu*(
      e[3*i+2]*e[3*j+2]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1] ) +
                           wdetlambda*e[3*i+2]*e[3*j+2];
     a = 0; b = 1;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+1]*e[3*j+0] ) +
                           wdetlambda*e[3*i+0]*e[3*j+1];
     a = 1; b = 0;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+0]*e[3*j+1] ) +
                           wdetlambda*e[3*i+1]*e[3*j+0];
     a = 0; b = 2;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+2]*e[3*j+0] ) +
                           wdetlambda*e[3*i+0]*e[3*j+2];
     a = 2; b = 0;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+0]*e[3*j+2] ) +
                           wdetlambda*e[3*i+2]*e[3*j+0];
     a = 1; b = 2;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+2]*e[3*j+1] ) +
                           wdetlambda*e[3*i+1]*e[3*j+2];
     a = 2; b = 1;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+1]*e[3*j+2] ) +
                           wdetlambda*e[3*i+2]*e[3*j+1];
   }
   for(j=0;j<nn;j++) {
     int a,b;
     int i = j;
     a = 0; b = 0;
     K[(3*j+b)*n+3*i+a] += wdet2mu*(
      e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1]+0.5*e[3*i+2]*e[3*j+2] ) +
                           wdetlambda*e[3*i+0]*e[3*j+0];
     a = 1; b = 1;
     K[(3*j+b)*n+3*i+a] += wdet2mu*(
      e[3*i+1]*e[3*j+1]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+2]*e[3*j+2] ) +
                           wdetlambda*e[3*i+1]*e[3*j+1];
     a = 2; b = 2;
     K[(3*j+b)*n+3*i+a] += wdet2mu*(
      e[3*i+2]*e[3*j+2]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1] ) +
                           wdetlambda*e[3*i+2]*e[3*j+2];
     a = 1; b = 0;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+0]*e[3*j+1] ) +
                           wdetlambda*e[3*i+1]*e[3*j+0];
     a = 2; b = 0;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+0]*e[3*j+2] ) +
                           wdetlambda*e[3*i+2]*e[3*j+0];
     a = 2; b = 1;
     K[(3*j+b)*n+3*i+a] += wdet2mu*( 
      0.5*e[3*i+1]*e[3*j+2] ) +
                           wdetlambda*e[3*i+2]*e[3*j+1];
   }
 }
};
#endif


class LEARubberStiffFunction : public IntegFunctionV3d {
 int n;
 complex<double> *K;
public:
 LEARubberStiffFunction(int _n, complex<double> *_K) {
   n = _n; K = _K; 
 }

 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det){
   int nn = n/3;
   double wdet = w*det;

   int i;
   double *e = (double*)alloca(sizeof(double)*3*nn);
   for(i=0;i<nn;i++) { 
     e[3*i+0] = dNdx[i][0];
     e[3*i+1] = dNdx[i][1];
     e[3*i+2] = dNdx[i][2];
   }
   {
     complex<double> wdet2mu = wdet*2.0;
     complex<double> wdetlambda = wdet;
     for(int j=0;j<nn;j++) for(i=j+1;i<nn;i++) {
       int a,b;
       a = 0; b = 0;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*(
        e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1]+0.5*e[3*i+2]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+0]*e[3*j+0];
       a = 1; b = 1;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*(
        e[3*i+1]*e[3*j+1]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+2]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+1]*e[3*j+1];
       a = 2; b = 2;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*(
        e[3*i+2]*e[3*j+2]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+2]*e[3*j+2];
       a = 0; b = 1;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+1]*e[3*j+0] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+0]*e[3*j+1];
       a = 1; b = 0;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+0]*e[3*j+1] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+1]*e[3*j+0];
       a = 0; b = 2;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+2]*e[3*j+0] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+0]*e[3*j+2];
       a = 2; b = 0;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+0]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+2]*e[3*j+0];
       a = 1; b = 2;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+2]*e[3*j+1] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+1]*e[3*j+2];
       a = 2; b = 1;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+1]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+2]*e[3*j+1];
     }
     for(int j=0;j<nn;j++) {
       int a,b;
       int i = j;
       a = 0; b = 0;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*(
        e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1]+0.5*e[3*i+2]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+0]*e[3*j+0];
       a = 1; b = 1;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*(
        e[3*i+1]*e[3*j+1]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+2]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+1]*e[3*j+1];
       a = 2; b = 2;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*(
        e[3*i+2]*e[3*j+2]+0.5*e[3*i+0]*e[3*j+0]+0.5*e[3*i+1]*e[3*j+1] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+2]*e[3*j+2];
       a = 1; b = 0;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+0]*e[3*j+1] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+1]*e[3*j+0];
       a = 2; b = 0;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+0]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+2]*e[3*j+0];
       a = 2; b = 1;
       K[0*n*n+(3*j+b)*n+3*i+a] += wdet2mu*( 
        0.5*e[3*i+1]*e[3*j+2] );
       K[1*n*n+(3*j+b)*n+3*i+a] += wdetlambda*e[3*i+2]*e[3*j+1];
     }
   }

 }
};



class LEMassFunction : public IntegFunctionV3d {
 int n;
 double rho;
 double *K;
public:
 LEMassFunction(int _n, double _rho, double *_K) {
   n = _n; rho = _rho; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det){
   int nn = n/3;
   double wdetrho = w*det*rho;

   for(int j=0;j<nn;j++) for(int i=j;i<nn;i++) {
     double v = wdetrho*N[i]*N[j];
     K[3*j*n+3*i] += v;
     K[(3*j+1)*n+3*i+1] += v;
     K[(3*j+2)*n+3*i+2] += v;
   }
 }
};


class LEGravityFunction : public IntegFunctionV3d {
 int n;
 double rho;
 double *g;
 double *v;
public:
 LEGravityFunction(int _n, double _rho, double *_g, double *_v) {
   n = _n; rho = _rho; g = _g; v = _v;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det){
   int nn = n/3;
   double wdetrho = w*det*rho;

   for(int j=0;j<nn;j++) {
     double d = wdetrho*N[j];
     v[3*j] += d*g[0];
     v[3*j+1] += d*g[1];
     v[3*j+2] += d*g[2];
   }
 }
};


#ifdef _TEMPLATE_FIX_
#include <Element.d/Helm.d/IsoParamUtilsT.C>
#endif

#endif
