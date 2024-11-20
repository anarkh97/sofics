#ifndef _PML_H_
#define _PML_H_

#include <alloca.h>
#include <Element.d/Helm.d/IntegFunction.h>


/*
struct PMLData {
 virtual int pmltype()=0;
};


struct CylindricalPMLData: public PMLData {
 virtual int pmltype() { return 1; }
 double R,S;
 double Rmz, Smz;
 double Rpz, Spz;
};


struct SphericalPMLData: public PMLData {
 virtual int pmltype() { return 2; }
 double R,S;
};


struct CartesianPMLData: public PMLData {
 virtual int pmltype() { return 0; }
 double Rmx, Rmy, Rmz, Smx, Smy, Smz;
 double Rpx, Rpy, Rpz, Spx, Spy, Spz;
};*/


struct SphericalCoord {
// notation as in http://mathworld.wolfram.com/SphericalCoordinates.html
 double r;
 double ct;
 double st;
 double cp;
 double sp;
 SphericalCoord(double *xyz) {
   double z2 = xyz[0]*xyz[0]+xyz[1]*xyz[1];
   double R = sqrt(z2);
   r = sqrt(xyz[2]*xyz[2]+z2);
   if (R == 0.0) {
     ct = 1.0; st = 0.0;
   } else {
     ct = xyz[0]/R;
     st = xyz[1]/R;
   }
   if (r == 0.0) {
     cp = 1.0; sp = 0.0;
   } else {
     cp = xyz[2]/r;
     sp = R/r; 
   }
 }
 void cartGrad2sphGrad(double *gxyz, double *grtp) {
   grtp[0] = sp*(ct*gxyz[0] + st*gxyz[1]) + cp*gxyz[2];
   grtp[1] = cp*(ct*gxyz[0] + st*gxyz[1]) - sp*gxyz[2];
   grtp[2] = -st*gxyz[0] + ct*gxyz[1];
 }
 void sphGrad2cartGrad(double *grtp, double *gxyz) {
   gxyz[0] = ct*(sp*grtp[0] + cp*grtp[1]) - st*grtp[2];
   gxyz[1] = st*(sp*grtp[0] + cp*grtp[1]) + ct*grtp[2];
   gxyz[2] = cp*grtp[0] - sp*grtp[1];
 }
 double radius() { return r; }
};

struct CylindricalCoord {
// notation as in http://mathworld.wolfram.com/CylindricalCoordinates.html
public:
 double r;
 double z;
 double ct;
 double st;
 CylindricalCoord(double *xyz) {
   r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
   if (r == 0.0) { 
     ct = 1.0; st = 0.0;
   } else {
     ct = xyz[0]/r;
     st = xyz[1]/r;
   }
   z = xyz[2];
 }
 void cartGrad2cylGrad(double *gxyz, double *grtz) {
   grtz[0] = ct*gxyz[0] + st*gxyz[1];
   grtz[1] = -st*gxyz[0] + ct*gxyz[1];
   grtz[2] = gxyz[2];
 }
 void cartGrad2cylGrad2d(double *gxyz, double *grtz) {
   grtz[0] = ct*gxyz[0] + st*gxyz[1];
   grtz[1] = -st*gxyz[0] + ct*gxyz[1];
 }
 void cylGrad2cartGrad(double *grtz, double *gxyz) {
   gxyz[0] = ct*grtz[0] - st*grtz[1];
   gxyz[1] = st*grtz[0] + ct*grtz[1];
   gxyz[2] = grtz[2];
 }
 double radius() { return r; }
};


struct PMLFunction {
 double R;
 double S;
 double gamma;
 
 complex<double> alpha(double r) {
  double xi = (r-R)/(S-R);
  if (xi<0.0) return complex<double>(1.0,0.0);
  if (xi>1.0) return complex<double>(1.0,0.0);
  double SR = S-R;
  double rRSR = (r-R)/SR;
  double integral = SR* ( exp(gamma*rRSR)/gamma*(rRSR-1.0/gamma) - 
                     rRSR*rRSR/2.0 + 1.0/(gamma*gamma));
  
  return complex<double>(1.0,integral/r);
 }
 complex<double> beta(double r) {
  double xi = (r-R)/(S-R);
//fprintf(stderr,"hehe %f %f %f %f\n",xi,r, R,S);
  if (xi<0.0 || xi>1.0) return complex<double>(1.0,0.0);
  return complex<double>(1.0,xi*(exp(gamma*xi)-1.0));
 }
 PMLFunction(double _RR, double _SS, double _gamma) {
   R = _RR; S = _SS;
   gamma = _gamma; 
 }
 PMLFunction() {}
};


class CartPMLGalStiffFunction : public IntegFunctionV3d {
 int n;
 double *Rm, *Rp, *Sm, *Sp;
 double gamma;
 PMLFunction fm[3],fp[3];
 complex<double> *K;
public:
 CartPMLGalStiffFunction(int _n, double *_Rm, double *_Sm,
                         double *_Rp, double *_Sp, double _gamma,
                         complex<double> *_K) {
   n = _n;
   Rm = _Rm;
   Rp = _Rp;
   Sm = _Sm;
   Sp = _Sp;
   gamma = _gamma;
   K = _K;

   fm[0] = PMLFunction(Rm[0],Sm[0],gamma);
   fm[1] = PMLFunction(Rm[1],Sm[1],gamma);
   fm[2] = PMLFunction(Rm[2],Sm[2],gamma);
   fp[0] = PMLFunction(Rp[0],Sp[0],gamma);
   fp[1] = PMLFunction(Rp[1],Sp[1],gamma);
   fp[2] = PMLFunction(Rp[2],Sp[2],gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   complex<double> beta[3] = {1.0,1.0,1.0}; 
   
   int i;
   for(i=0;i<3;i++)
     if (x[i]>Rp[i]) beta[i] = fp[i].beta(x[i]);
     else if (x[i]<Rm[i])  beta[i] = fm[i].beta(x[i]);


   complex<double> wdetb = w*det*beta[0]*beta[1]*beta[2];
   complex<double> invb[3] = { 1.0/(beta[0]*beta[0]), 
                               1.0/(beta[1]*beta[1]),
                               1.0/(beta[2]*beta[2])}; 

   int j;
   for(j=0;j<n;j++)
     for(i=j;i<n;i++)
       K[j*n+i] += wdetb*(
        invb[0]*dNdx[i][0]*dNdx[j][0]+
        invb[1]*dNdx[i][1]*dNdx[j][1]+
        invb[2]*dNdx[i][2]*dNdx[j][2]);
 }
};


class CartPMLGalStiffFunction2d: public IntegFunctionA2d {
 int n;
 double *Rm, *Rp, *Sm, *Sp;
 double gamma;
 PMLFunction fm[2],fp[2];
 complex<double> *K;
public:
 CartPMLGalStiffFunction2d(int _n, double *_Rm, double *_Sm,
                         double *_Rp, double *_Sp, double _gamma,
                         complex<double> *_K) {
   n = _n;
   Rm = _Rm;
   Rp = _Rp;
   Sm = _Sm;
   Sp = _Sp;
   gamma = _gamma;
   K = _K;

   fm[0] = PMLFunction(Rm[0],Sm[0],gamma);
   fm[1] = PMLFunction(Rm[1],Sm[1],gamma);
   fp[0] = PMLFunction(Rp[0],Sp[0],gamma);
   fp[1] = PMLFunction(Rp[1],Sp[1],gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {

   complex<double> beta[2] = {1.0,1.0}; 

   for(int i=0;i<2;i++)
     if (x[i]>Rp[i])  beta[i] = fp[i].beta(x[i]);
     else if (x[i]<Rm[i])  beta[i] = fm[i].beta(x[i]);

   complex<double> wdetb = w*det*beta[0]*beta[1];
   complex<double> invb[2] = { 1.0/(beta[0]*beta[0]), 
                               1.0/(beta[1]*beta[1]) };

   int i,j;
   for(j=0;j<n;j++)
     for(i=j;i<n;i++)
       K[j*n+i] += wdetb*(
        invb[0]*dNdx[i][0]*dNdx[j][0]+
        invb[1]*dNdx[i][1]*dNdx[j][1]);
 }
};


class CartPMLGalMassFunction : public IntegFunctionV3d {
 int n;
 double *Rm, *Rp, *Sm, *Sp;
 double gamma;
 PMLFunction fm[3],fp[3];
 complex<double> *K;
public:
 CartPMLGalMassFunction(int _n, double *_Rm, double *_Sm,
                        double *_Rp, double *_Sp, double _gamma,
                        complex<double> *_K) {
   n = _n;
   Rm = _Rm;
   Rp = _Rp;
   Sm = _Sm;
   Sp = _Sp;
   gamma = _gamma;
   K = _K;

   fm[0] = PMLFunction(Rm[0],Sm[0],gamma);
   fm[1] = PMLFunction(Rm[1],Sm[1],gamma);
   fm[2] = PMLFunction(Rm[2],Sm[2],gamma);
   fp[0] = PMLFunction(Rp[0],Sp[0],gamma);
   fp[1] = PMLFunction(Rp[1],Sp[1],gamma);
   fp[2] = PMLFunction(Rp[2],Sp[2],gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   complex<double> beta[3] = {1.0,1.0,1.0}; 

   for(int i=0;i<3;i++)
     if (x[i]>Rp[i]) beta[i] = fp[i].beta(x[i]);
     else if (x[i]<Rm[i])  beta[i] = fm[i].beta(x[i]);

   complex<double> wdetb = w*det*beta[0]*beta[1]*beta[2];

//fprintf(stderr,"hehe %f %f\n",real(beta[0]*beta[1]*beta[2]),imag(beta[0]*beta[1]*beta[2]));

   int i,j;
   for(j=0;j<n;j++)
     for(i=j;i<n;i++)
       K[j*n+i] += wdetb*N[i]*N[j];
 }
};


class CartPMLGalMassFunction2d : public IntegFunctionA2d {
 int n;
 double *Rm, *Rp, *Sm, *Sp;
 double gamma;
 PMLFunction fm[2],fp[2];
 complex<double> *K;
public:
 CartPMLGalMassFunction2d(int _n, double *_Rm, double *_Sm,
                          double *_Rp, double *_Sp, double _gamma,
                          complex<double> *_K) {
   n = _n;
   Rm = _Rm;
   Rp = _Rp;
   Sm = _Sm;
   Sp = _Sp;
   gamma = _gamma;
   K = _K;

   fm[0] = PMLFunction(Rm[0],Sm[0],gamma);
   fm[1] = PMLFunction(Rm[1],Sm[1],gamma);
   fp[0] = PMLFunction(Rp[0],Sp[0],gamma);
   fp[1] = PMLFunction(Rp[1],Sp[1],gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {

   complex<double> beta[2] = {1.0,1.0};

   for(int i=0;i<2;i++)
     if (x[i]>Rp[i]) beta[i] = fp[i].beta(x[i]);
     else if (x[i]<Rm[i])  beta[i] = fm[i].beta(x[i]);

   complex<double> wdetb = w*det*beta[0]*beta[1];

   int i,j;
   for(j=0;j<n;j++)
     for(i=j;i<n;i++)
       K[j*n+i] += wdetb*N[i]*N[j];
 }
};


class SphPMLGalStiffFunction : public IntegFunctionV3d {
 int n;
 double R, S;
 double gamma;
 PMLFunction f;
 complex<double> *K;
public:
 SphPMLGalStiffFunction(int _n, double _RR, double _SS, double _gamma,
                        complex<double> *_KK) {
   n = _n;
   R = _RR;
   S = _SS;
   gamma = _gamma;
   K = _KK;

   f = PMLFunction(R,S,gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   SphericalCoord sc(x); 
   double r = sc.radius();
   complex<double> beta = f.beta(r);
   complex<double> alpha = f.alpha(r);

   complex<double> wdetb = w*det*beta*alpha*alpha;
  
   complex<double> (*grad)[3] = (complex<double>(*)[3])
                                  alloca(sizeof(complex<double>)*3*n);
   for(int i=0; i<n; i++) {
     double grtp[3];
     sc.cartGrad2sphGrad(dNdx[i],grtp);
     grad[i][0] = grtp[0]/beta;
     grad[i][1] = grtp[1]/alpha;
     grad[i][2] = grtp[2]/alpha;
   } 
  
   for(int j=0;j<n;j++)
     for(int i=j;i<n;i++)
       K[j*n+i] += wdetb*(
        grad[i][0]*grad[j][0]+
        grad[i][1]*grad[j][1]+
        grad[i][2]*grad[j][2]);
 }
};


class SphPMLGalMassFunction : public IntegFunctionV3d {
 int n;
 double R, S;
 double gamma;
 PMLFunction f;
 complex<double> *K;
public:
 SphPMLGalMassFunction(int _n, double _RR, double _SS, double _gamma,
                       complex<double> *_KK) {
   n = _n;
   R = _RR;
   S = _SS;
   gamma = _gamma;
   K = _KK;

   f = PMLFunction(R,S,gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   SphericalCoord sc(x); 
   double r = sc.radius();
   complex<double> beta = f.beta(r);
   complex<double> alpha = f.alpha(r);

   complex<double> wdetb = w*det*beta*alpha*alpha;
  
   for(int j=0;j<n;j++)
     for(int i=j;i<n;i++)
       K[j*n+i] += wdetb*N[i]*N[j];
 }
};


class CylPMLGalStiffFunction : public IntegFunctionV3d {
 int n;
 double R, S, Dm, Dp, Lm, Lp;
 double gamma;
 PMLFunction fr,fm,fp;
 complex<double> *K;
public:
 CylPMLGalStiffFunction(int _n, double _RR, double _SS, double _Dm, double _Dp,
                        double _Lm, double _Lp, double _gamma,
                        complex<double> *_KK) {
   n = _n;
   R = _RR; S = _SS;
   Dm = _Dm; Dp = _Dp; Lm = _Lm; Lp = _Lp;
   gamma = _gamma;
   K = _KK;

   fr = PMLFunction(R,S,gamma);
   fm = PMLFunction(Dm,Lm,gamma);
   fp = PMLFunction(Dp,Lp,gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   CylindricalCoord cc(x); 

   double r = cc.radius();
   complex<double> betar = fr.beta(r);
   complex<double> alpha = fr.alpha(r);

   complex<double> betaz = 1.0;
   if (x[2]>Dp) betaz = fp.beta(x[2]);
   else if (x[2]<Dm)  betaz = fm.beta(x[2]);

   complex<double> wdetb = w*det*betar*betaz*alpha;
  
   complex<double> (*grad)[3] = (complex<double>(*)[3])
                                  alloca(sizeof(complex<double>)*3*n);
   for(int i=0; i<n; i++) {
     double grtz[3];
     cc.cartGrad2cylGrad(dNdx[i],grtz);
     grad[i][0] = grtz[0]/betar;
     grad[i][1] = grtz[1]/alpha;
     grad[i][2] = grtz[2]/betaz;
   } 
  
   for(int j=0;j<n;j++)
     for(int i=j;i<n;i++)
       K[j*n+i] += wdetb*(
        grad[i][0]*grad[j][0]+
        grad[i][1]*grad[j][1]+
        grad[i][2]*grad[j][2]);
 }
};


class CylPMLGalStiffFunction2d : public IntegFunctionA2d {
 int n;
 double R, S;
 double gamma;
 PMLFunction fr;
 complex<double> *K;
public:
 CylPMLGalStiffFunction2d(int _n, double _RR, double _SS, double _gamma,
                          complex<double> *_KK) {
   n = _n;
   R = _RR; S = _SS;
   gamma = _gamma;
   K = _KK;

   fr = PMLFunction(R,S,gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {

   CylindricalCoord cc(x); 

   double r = cc.radius();
   complex<double> betar = fr.beta(r);
   complex<double> alpha = fr.alpha(r);

   complex<double> wdetb = w*det*betar*alpha;

   complex<double> (*grad)[2] = (complex<double>(*)[2])
                                  alloca(sizeof(complex<double>)*2*n);
   for(int i=0; i<n; i++) {
     double grtz[2];
     cc.cartGrad2cylGrad2d(dNdx[i],grtz);
     grad[i][0] = grtz[0]/betar;
     grad[i][1] = grtz[1]/alpha;
   }
  
   for(int j=0;j<n;j++)
     for(int i=j;i<n;i++) {
       K[j*n+i] += wdetb*(
        grad[i][0]*grad[j][0]+
        grad[i][1]*grad[j][1]);
}
 }
};


class CylPMLGalMassFunction : public IntegFunctionV3d {
 int n;
 double R, S, Dm, Dp, Lm, Lp;
 double gamma;
 PMLFunction fr,fm,fp;
 complex<double> *K;
public:
 CylPMLGalMassFunction(int _n, double _RR, double _SS, double _Dm, double _Dp,
                       double _Lm, double _Lp, double _gamma,
                       complex<double> *_KK) {
   n = _n;
   R = _RR; S = _SS;
   Dm = _Dm; Dp = _Dp; Lm = _Lm; Lp = _Lp;
   gamma = _gamma;
   K = _KK;

   fr = PMLFunction(R,S,gamma);
   fm = PMLFunction(Dm,Lm,gamma);
   fp = PMLFunction(Dp,Lp,gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   CylindricalCoord cc(x); 

   double r = cc.radius();
   complex<double> betar = fr.beta(r);
   complex<double> alpha = fr.alpha(r);

   complex<double> betaz = 1.0;
   if (x[2]>Dp) betaz = fp.beta(x[2]);
   else if (x[2]<Dm)  betaz = fm.beta(x[2]);

   complex<double> wdetb = w*det*betar*betaz*alpha;
  
   for(int j=0;j<n;j++)
     for(int i=j;i<n;i++)
       K[j*n+i] += wdetb*N[i]*N[j];
 }
};


class CylPMLGalMassFunction2d : public IntegFunctionA2d {
 int n;
 double R, S;
 double gamma;
 PMLFunction fr;
 complex<double> *K;
public:
 CylPMLGalMassFunction2d(int _n, double _RR, double _SS, double _gamma,
                         complex<double> *_KK) {
   n = _n;
   R = _RR; S = _SS;
   gamma = _gamma;
   K = _KK;

   fr = PMLFunction(R,S,gamma);
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {

   CylindricalCoord cc(x); 

   double r = cc.radius();
   complex<double> betar = fr.beta(r);
   complex<double> alpha = fr.alpha(r);

   complex<double> wdetb = w*det*betar*alpha;

   for(int j=0;j<n;j++)
     for(int i=j;i<n;i++)
       K[j*n+i] += wdetb*N[i]*N[j];
 }
};

#endif
