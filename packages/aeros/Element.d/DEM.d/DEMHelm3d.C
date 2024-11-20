#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>

#ifndef SANDIA
#include <Element.d/DEM.d/DEMHelm3d.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/PML.h>
#else
#include "DEMHelm3d.h"
#include "IsoParamUtils.h"
#endif

class HelmDGMEMatricesFunction3d : public IntegFunctionA3d {
 int oc;
 double kappa;
 int ndir;
 complex<double> *dirs;
 int nldir;
 int *nldirs;
 complex<double> *ldirs;
 complex<double> *K;
 complex<double> *L;
 complex<double> *PL;
 complex<double> *PP;
 complex<double> *PE;
 int arbflag;
 int fi;
 double *xsc;
 double *xc;
public:
 HelmDGMEMatricesFunction3d(int _oc, double _kappa,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs,
                      double *_xsc, double *_xc, int _arbflag, int _fi,
                      complex<double> *_K, complex<double> *__L,
                      complex<double> *_PL,complex<double> *_PP,
                      complex<double> *_PE) {
   oc = _oc; kappa = _kappa; ndir = _ndir; dirs = _dirs; nldir = _nldir; 
   ldirs = _ldirs; K = _K; L = __L; PL = _PL; PP = _PP; PE =_PE;
   arbflag = _arbflag; fi = _fi;
   xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   double wc = w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

   for(int j=0;j<nldir;j++) {
     complex<double> l = exp(ldirs[j*3+0]*(x[0]-xsc[0])+
                             ldirs[j*3+1]*(x[1]-xsc[1])+
                             ldirs[j*3+2]*(x[2]-xsc[2]) );
     l *= wc;
     for(int i=0;i<ndir;i++)
       L[j*ndir+i] += e[i]*l;
     if (PL!=0) for(int i=0;i<oc;i++)
       PL[j*oc+i] += N[i]*l;
   }

   complex<double> ikc =
     (arbflag==0)? 0.0: complex<double>(0.0, w*kappa*
                 sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
   for(int j=0;j<ndir;j++) {
     complex<double> ic = w*nsign*
              (cross[0]*dirs[j*3+0]+cross[1]*dirs[j*3+1]+cross[2]*dirs[j*3+2]);
     ic -= ikc;
     for(int i=j;i<ndir;i++)
       K[j*ndir+i] += ic*e[i]*e[j];
   }
// Sommerfeld part of bc
   if (PP!=0)  {
     for(int j=0;j<oc;j++) {
       for(int i=j;i<oc;i++) PP[j*oc+i] -= ikc*N[i]*N[j];
       for(int i=0;i<ndir;i++) PE[i*oc+j] -= ikc*N[j]*e[i];
     }
   }
   
   delete[] e;
 }
};


class HelmDGMELMatrixFunction3d : public IntegFunctionA3d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 int nldir;
 int *nldirs;
 complex<double> *ldirs;
 complex<double> *L;
 int fi;
 double *xsc;
 double *xc;
public:
 HelmDGMELMatrixFunction3d(double _kappa,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs, double *_xsc, double *_xc,
                      complex<double> *__L, int _fi) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; nldir = _nldir; ldirs = _ldirs;
   L = __L; fi = _fi; xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
  
   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   double wc = w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
   for(int j=0;j<nldir;j++) {
     complex<double> l = exp(ldirs[j*3+0]*(x[0]-xsc[0])+
                             ldirs[j*3+1]*(x[1]-xsc[1])+
                             ldirs[j*3+2]*(x[2]-xsc[2]) );
     l *= wc;
     for(int i=0;i<ndir;i++)
       L[j*ndir+i] += e[i]*l;
   }
   delete[] e;
 }
};


class HelmDGMSomEEMatrixFunction3d : public IntegFunctionA3d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *K;
public:
 HelmDGMSomEEMatrixFunction3d(double _kappa,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   complex<double> ikc = complex<double>(0.0, w*kappa*
           sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] -= ikc*e[i]*e[j];
   delete[] e; 
 }
};


class HelmDGMEeMatrixFunction3d : IntegFunctionA3d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *K;
public:
 HelmDGMEeMatrixFunction3d(double _kappa,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[2*ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   complex<double> *ikc = e+ndir;
   for(int i=0;i<ndir;i++) 
     ikc[i] = w*nsign*
              (cross[0]*dirs[i*3+0]+cross[1]*dirs[i*3+1]+cross[2]*dirs[i*3+2]);

   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] += ikc[i]*e[i]*e[j];
   delete[] e;
 }
};

static double sidelen(double x, double y, double z) {
 return sqrt(x*x+y*y+z*z);
}



void DGMHelm3d::HelmDGMEMatricesExactFace3d(double *xyz,
                    int ndir, complex<double> *dirs,
                    int nldir, complex<double> *ldirs,
                    double kappa,  int sflag,
                    double *xsc,  double *xc, int faceindex,
                    complex<double> *K, complex<double> *L) {


 int nfc = nFaceCorners(faceindex);
 if (nfc>4) {
   fprintf(stderr,"Face with more than 4 vertices not supported in DGMHelm3d::HelmDGMEMatricesExactFace3d\n"); 
   exit(-1);
 }

 int *corner = faceCornerI(faceindex);
 double (*polyxyz)[3] = new double[nfc+1][3];
 int orderc = nGeomNodes();

 for(int j=0;j<nfc;j++) {
   polyxyz[j][0] = xyz[0*orderc+corner[j]];
   polyxyz[j][1] = xyz[1*orderc+corner[j]];
   polyxyz[j][2] = xyz[2*orderc+corner[j]];
 }
 polyxyz[nfc][0] = polyxyz[0][0];
 polyxyz[nfc][1] = polyxyz[0][1];
 polyxyz[nfc][2] = polyxyz[0][2];
 delete[] corner;

 double xf[3],af[3],bf[3];
 xf[0] = polyxyz[0][0];
 xf[1] = polyxyz[0][1];
 xf[2] = polyxyz[0][2];
 af[0] = polyxyz[1][0] - xf[0]; 
 af[1] = polyxyz[1][1] - xf[1]; 
 af[2] = polyxyz[1][2] - xf[2]; 
 bf[0] = polyxyz[2][0] - xf[0]; 
 bf[1] = polyxyz[2][1] - xf[1]; 
 bf[2] = polyxyz[2][2] - xf[2];
 double area = 0.5*sidelen(af[1]*bf[2]-af[2]*bf[1],
                     af[2]*bf[0]-af[0]*bf[2],af[0]*bf[1]-af[1]*bf[0]);

 if (nfc==4) {
   double cf[3];
   cf[0] = polyxyz[3][0] - xf[0]; 
   cf[1] = polyxyz[3][1] - xf[1]; 
   cf[2] = polyxyz[3][2] - xf[2];
   area += 0.5*sidelen(cf[1]*bf[2]-cf[2]*bf[1],
                     cf[2]*bf[0]-cf[0]*bf[2],cf[0]*bf[1]-cf[1]*bf[0]);
 }
 double l1 = sidelen(af[0],af[1],af[2]);
 af[0] /= l1; af[1] /= l1; af[2] /= l1;
 double alpha = af[0]*bf[0]+af[1]*bf[1]+af[2]*bf[2];
 bf[0] -= alpha*af[0];
 bf[1] -= alpha*af[1];
 bf[2] -= alpha*af[2];
 double l2 = sidelen(bf[0],bf[1],bf[2]);
 bf[0] /= l2; bf[1] /= l2; bf[2] /= l2;
 double nf[3];
 nf[0] = af[1]*bf[2]-af[2]*bf[1];
 nf[1] = af[2]*bf[0]-af[0]*bf[2];
 nf[2] = af[0]*bf[1]-af[1]*bf[0];

 double tau[4][3];
 double aftau[4];
 double bftau[4];
 double lentau[4];
 for(int e=0;e<nfc;e++) {
   tau[e][0] = polyxyz[e+1][0]-polyxyz[e][0];
   tau[e][1] = polyxyz[e+1][1]-polyxyz[e][1];
   tau[e][2] = polyxyz[e+1][2]-polyxyz[e][2];
   lentau[e] = sidelen(tau[e][0],tau[e][1],tau[e][2]);
   tau[e][0] /= lentau[e]; tau[e][1] /= lentau[e]; tau[e][2] /= lentau[e];
   aftau[e]= tau[e][0]*af[0]+tau[e][1]*af[1]+tau[e][2]*af[2];
   bftau[e]= tau[e][0]*bf[0]+tau[e][1]*bf[1]+tau[e][2]*bf[2];
 }

 complex<double> *edir = new complex<double>[ndir*nfc];
 complex<double> *taudir = new complex<double>[ndir*nfc];
 for(int j=0;j<nfc;j++)
   for(int l=0;l<ndir;l++) {
     edir[l*nfc+j] = exp(dirs[l*3+0]*(polyxyz[j][0]-xc[0])+
		      dirs[l*3+1]*(polyxyz[j][1]-xc[1])+
		      dirs[l*3+2]*(polyxyz[j][2]-xc[2]));
     taudir[l*nfc+j] =
       dirs[l*3+0]*tau[j][0]+dirs[l*3+1]*tau[j][1]+dirs[l*3+2]*tau[j][2];
   }

 complex<double> *diraf = new complex<double>[ndir];
 complex<double> *dirbf = new complex<double>[ndir];
 for(int l=0;l<ndir;l++) {
   diraf[l] = dirs[l*3+0]*af[0]+dirs[l*3+1]*af[1]+dirs[l*3+2]*af[2];
   dirbf[l] = dirs[l*3+0]*bf[0]+dirs[l*3+1]*bf[1]+dirs[l*3+2]*bf[2];
 }

 complex<double> *eldir = new complex<double>[nldir*nfc];
 complex<double> *tauldir = new complex<double>[nldir*nfc];
 for(int j=0;j<nfc;j++)
   for(int l=0;l<nldir;l++) {
     eldir[l*nfc+j] = exp(ldirs[l*3+0]*(polyxyz[j][0]-xsc[0])+
		       ldirs[l*3+1]*(polyxyz[j][1]-xsc[1])+
		       ldirs[l*3+2]*(polyxyz[j][2]-xsc[2]));
     tauldir[l*nfc+j] =
       ldirs[l*3+0]*tau[j][0]+ldirs[l*3+1]*tau[j][1]+ldirs[l*3+2]*tau[j][2];

   }


 complex<double> *ldiraf = new complex<double>[nldir];
 complex<double> *ldirbf = new complex<double>[nldir];
 for(int l=0;l<nldir;l++) {
   ldiraf[l] = ldirs[l*3+0]*af[0]+ldirs[l*3+1]*af[1]+ldirs[l*3+2]*af[2];
   ldirbf[l] = ldirs[l*3+0]*bf[0]+ldirs[l*3+1]*bf[1]+ldirs[l*3+2]*bf[2];
 }


 for(int kk=0;kk<nldir;kk++) {
   for(int l=0;l<ndir;l++) {
     complex<double> bt1 = diraf[l]+ldiraf[kk];
     complex<double> bt2 = dirbf[l]+ldirbf[kk];
     complex<double> betatau = std::norm(bt1)+std::norm(bt2);
     if (abs(betatau)<1e-12) {
       L[(kk)*ndir+l] +=area*edir[l*nfc+0]*eldir[kk*nfc+0];
     } else {

       for(int e=0;e<nfc;e++) {
         complex<double> nbetatau =
           (conj(bt1)* bftau[e] - conj(bt2)* aftau[e]) 
             / betatau;
         complex<double> ebetatau =
            tauldir[kk*nfc+e]+taudir[l*nfc+e];
         if (abs(ebetatau)<1e-6) {
           L[kk*ndir+l] += nbetatau*lentau[e]*edir[l*nfc+e]*eldir[kk*nfc+e];
         } else {
           int ep1 = (e==nfc-1)?0:e+1;
           L[kk*ndir+l] += nbetatau/ebetatau * 
             (edir[l*nfc+ep1]*eldir[kk*nfc+ep1]-edir[l*nfc+e]*eldir[kk*nfc+e]);
         }
       }
     }
   }
 }

 for(int kk=0;kk<ndir;kk++) {
   complex<double> C;
   if (sflag==0)
     C=dirs[kk*3+0]*nf[0]+dirs[kk*3+1]*nf[1]+dirs[kk*3+2]*nf[2];
   else
     C=dirs[kk*3+0]*nf[0]+dirs[kk*3+1]*nf[1]+dirs[kk*3+2]*nf[2] -
       complex<double>(0.0,kappa);
   for(int l=kk;l<ndir;l++) {
     complex<double> bt1 = diraf[l]+diraf[kk];
     complex<double> bt2 = dirbf[l]+dirbf[kk];
     complex<double> betatau = std::norm(bt1)+std::norm(bt2);
     if (abs(betatau)<1e-12) {
       K[(kk)*ndir+l] += area*C*edir[kk*nfc+0]*edir[l*nfc+0];
     } else {
       for(int e=0;e<nfc;e++) {
         complex<double> nbetatau =
           (conj(bt1)* bftau[e] - conj(bt2)* aftau[e]) 
             / betatau;
         complex<double> ebetatau =
            taudir[kk*nfc+e]+taudir[l*nfc+e];
         if (abs(ebetatau)<1e-6)
           K[(kk)*ndir+l] += C*nbetatau*lentau[e]*edir[kk*nfc+e]*edir[l*nfc+e];
         else {
           int ep1 = (e==nfc-1)?0:e+1;
           complex<double> CC = C*nbetatau/ebetatau;
           K [(kk)*ndir+l] += CC*
             (edir[kk*nfc+ep1]*edir[l*nfc+ep1]-edir[kk*nfc+e]*edir[l*nfc+e]);
         }
       }
     }
   }
 }

 delete[] polyxyz;
 delete[] eldir;
 delete[] tauldir;
 delete[] edir;
 delete[] taudir;
 delete[] diraf;
 delete[] dirbf;
 delete[] ldiraf;
 delete[] ldirbf;
}


void DGMHelm3d::HelmDGMEMatrices3d(double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc,
                    complex<double> *K, complex<double> *L) {
 
 int nldir = 0;
 int nf = nFaces();
 for(int i=0;i<nf;i++) nldir += nldirs[i];
 for(int i=0;i<nldir*ndir;i++) L[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) K[i] = 0.0;

 int order = o;
 IsoParamUtils ipu(order);

 int c = 0;
 for(int faceindex=1;faceindex<=nf;faceindex++) {
     int isFlt = isFlatAndStraight(xyz,faceindex);
//     fprintf(stderr," FLAT  %i   \n",isFlt);
     
   if (!isFlt) {
//#define TIME
#ifdef TIME
 struct timespec tp1;
 struct timespec tp2;
 clock_gettime(CLOCK_REALTIME, &tp1);
#endif
     HelmDGMEMatricesFunction3d f(nGeomNodes(),kappa,ndir,dirs,
                  nldirs[faceindex-1],ldirs+c*3,
                  xsc + (faceindex-1)*3, xc, sflags[faceindex-1], faceindex,
                  K,L+c*ndir,0,0,0);
     surfInt3d(xyz, faceindex, f);
#ifdef TIME
 clock_gettime(CLOCK_REALTIME, &tp2);
 fprintf(stderr,"gauss: %ld nanoseconds\n",tp2.tv_nsec-tp1.tv_nsec);
#endif

   } else {
#ifdef TIME
 struct timespec tp1;
 struct timespec tp2;
 clock_gettime(CLOCK_REALTIME, &tp1);
#endif

     HelmDGMEMatricesExactFace3d(xyz, ndir, dirs, nldirs[faceindex-1],ldirs+c*3,
                    kappa,  sflags[faceindex-1],
                    xsc+(faceindex-1)*3,  xc, faceindex,
                    K, L+c*ndir);

#ifdef TIME
 clock_gettime(CLOCK_REALTIME, &tp2);
 fprintf(stderr,"exact: %ld nanoseconds\n",tp2.tv_nsec-tp1.tv_nsec);
#endif

/*
      complex<double>* KK = new complex<double>[ndir*ndir];
      complex<double>* LL = new complex<double>[ndir*nldirs[faceindex-1]];
      for(int i=0;i<nldirs[faceindex-1]*ndir;i++) LL[i] = 0.0;
      for(int i=0;i<ndir*ndir;i++) KK[i] = 0.0;
     HelmDGMEMatricesExactFace3d(xyz, ndir, dirs, nldirs[faceindex-1],ldirs+c*3,
                    kappa,  sflags[faceindex-1],
                    xsc+(faceindex-1)*3,  xc, faceindex,
                    KK, LL);
      complex<double>* KKK = new complex<double>[ndir*ndir];
      complex<double>* LLL = new complex<double>[ndir*nldirs[faceindex-1]];
      for(int i=0;i<nldirs[faceindex-1]*ndir;i++) LLL[i] = 0.0;
      for(int i=0;i<ndir*ndir;i++) KKK[i] = 0.0;
     HelmDGMEMatricesFunction3d f(nGeomNodes(),kappa,ndir,dirs,
                  nldirs[faceindex-1],ldirs+c*3,
                  xsc + (faceindex-1)*3, xc, sflags[faceindex-1], faceindex,
                  KKK,LLL,0,0,0);
     surfInt3d(xyz, faceindex, f);
 for(int i=0;i<ndir;i++)
 for(int j=0;j<ndir;j++)
fprintf(stderr,"%d %d: %e %e  %e %e\n",i+1,j+1,
real(KK[i*ndir+j]),imag(KK[i*ndir+j]),
real(KKK[i*ndir+j]),imag(KKK[i*ndir+j]));
      delete[] KK;
      delete[] KKK;
      delete[] LL;
      delete[] LLL;*/
   }
   c += nldirs[faceindex-1];
 } 

 for(int j=0;j<ndir;j++) {
   for(int i=0;i<j;i++)
     K[j*ndir+i] = K[i*ndir+j];
 }

}


class HelmDGMPMLEEMatrixFunction3d : public IntegFunctionV3d {
 double *pmldata;
 double *cxy;
 int ndir;
 complex<double> *dirs;
 complex<double> *K;
public:
 HelmDGMPMLEEMatrixFunction3d(double *_pmldata,
                      int _ndir, complex<double> *_dirs,
                            double *_cxy, complex<double> *_K) {
   pmldata = _pmldata;
   ndir = _ndir; dirs = _dirs;
   cxy = _cxy;
   K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {

   double kappa = pmldata[0];
   double a = pmldata[1];
   double b = pmldata[2];
   double gamma = pmldata[3];

   PMLFunction fp(a,b,gamma);
   PMLFunction fm(-a,-b,gamma);
   PMLFunction *fx = &fp;
   PMLFunction *fy = &fp;
   PMLFunction *fz = &fp;
   if (x[0]<0.0) fx = &fm;
   if (x[1]<0.0) fy = &fm;
   if (x[2]<0.0) fz = &fm;
   complex<double> wdetb = w*det*fx->beta(x[0])*fy->beta(x[1])*fz->beta(x[2]);

   complex<double> xt[3];
   xt[0] = x[0]*fx->alpha(x[0]);
   xt[1] = x[1]*fy->alpha(x[1]);
   xt[2] = x[2]*fz->alpha(x[2]);

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++)
     e[i] = exp(dirs[i*3+0]*xt[0] + dirs[i*3+1]*xt[1]+
                dirs[i*3+2]*xt[2]);

   for(int j=0;j<ndir;j++)
     for(int i=j;i<ndir;i++)
       K[j*ndir+i] += wdetb*
         (dirs[i*3+0]*dirs[j*3+0] + dirs[i*3+1]*dirs[j*3+1] +
          dirs[i*3+2]*dirs[j*3+2]-
                         complex<double>(kappa*kappa,0.0) ) * e[i] * e[j];
   delete[] e;
 }
};


class HelmDGMPMLELMatrixFunction3d_ : public IntegFunctionAt3d {
 double *pmldata;
 int ndir;
 complex<double> *dirs;
 int nldir;
 complex<double> *ldirs;
 complex<double> *L;
 int fi;
 double *xsc;
 double *xc;
public:
 HelmDGMPMLELMatrixFunction3d_(double *_pmldata,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs,
                      double *_xsc, double *_xc,
                      complex<double> *__L, int _fi) {
   pmldata = _pmldata;
   ndir = _ndir; dirs = _dirs; nldir = _nldir; ldirs = _ldirs;
   L = __L; fi = _fi; xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *tau1, double *tau2,
               double nsign, double w) {

   double kappa = pmldata[0];
   double a = pmldata[1];
   double b = pmldata[2];
   double gamma = pmldata[3];

   PMLFunction fp(a,b,gamma);
   PMLFunction fm(-a,-b,gamma);
   PMLFunction *fx = &fp;
   PMLFunction *fy = &fp;
   PMLFunction *fz = &fp;
   if (x[0]<0.0) fx = &fm;
   if (x[1]<0.0) fy = &fm;
   if (x[2]<0.0) fz = &fm;
   complex<double> cross[3] = {0,0,0};
   cross[0] = fy->beta(x[1])*fz->beta(x[2])*(tau1[1]*tau2[2] - tau1[2]*tau2[1]);
   cross[1] = fz->beta(x[2])*fx->beta(x[0])*(tau1[2]*tau2[0] - tau1[0]*tau2[2]);
   cross[2] = fx->beta(x[0])*fy->beta(x[1])*(tau1[0]*tau2[1] - tau1[1]*tau2[0]);
   complex<double> wcb = w  * sqrt(
                            cross[0]*conj(cross[0])+
                            cross[1]*conj(cross[1]) + cross[2]*conj(cross[2]));

   complex<double> xt[3];
   xt[0] = x[0]*fx->alpha(x[0]);
   xt[1] = x[1]*fy->alpha(x[1]);
   xt[2] = x[2]*fz->alpha(x[2]);

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++)
     e[i] = exp(dirs[i*3+0]*xt[0] + dirs[i*3+1]*xt[1]+
                dirs[i*3+2]*xt[2]);

   for(int j=0;j<nldir;j++)
     for(int i=0;i<ndir;i++) {
       complex<double> t = wcb*e[i]*
                           exp(ldirs[j*3+0]*(xt[0]) +
                               ldirs[j*3+1]*(xt[1]) +
                               ldirs[j*3+2]*(xt[2]) );

       L[j*ndir+i] += t;
     }
   delete[] e;
 }
};


class HelmDGMPMLELMatrixFunction3d : public IntegFunctionA3d {
 double *pmldata;
 int ndir;
 complex<double> *dirs;
 int nldir;
 complex<double> *ldirs;
 complex<double> *L;
 int fi;
 double *xsc;
 double *xc;
public:
 HelmDGMPMLELMatrixFunction3d(double *_pmldata,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs,
                      double *_xsc, double *_xc,
                      complex<double> *__L, int _fi) {
   pmldata = _pmldata;
   ndir = _ndir; dirs = _dirs; nldir = _nldir; ldirs = _ldirs;
   L = __L; fi = _fi; xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross,
               double nsign, double w) {

   double kappa = pmldata[0];
   double a = pmldata[1];
   double b = pmldata[2];
   double gamma = pmldata[3];

   PMLFunction fp(a,b,gamma);
   PMLFunction fm(-a,-b,gamma);
   PMLFunction *fx = &fp;
   PMLFunction *fy = &fp;
   PMLFunction *fz = &fp;
   if (x[0]<0.0) fx = &fm;
   if (x[1]<0.0) fy = &fm;
   if (x[2]<0.0) fz = &fm;
   complex<double> cr[3] = {0,0,0};
   cr[0] = fy->beta(x[1])*fz->beta(x[2])*cross[0];
   cr[1] = fz->beta(x[2])*fx->beta(x[0])*cross[1];
   cr[2] = fx->beta(x[0])*fy->beta(x[1])*cross[2];
   complex<double> wcb = w  * sqrt(
                            cr[0]*conj(cr[0])+
                            cr[1]*conj(cr[1]) + cr[2]*conj(cr[2]));

   complex<double> xt[3];
   xt[0] = x[0]*fx->alpha(x[0]);
   xt[1] = x[1]*fy->alpha(x[1]);
   xt[2] = x[2]*fz->alpha(x[2]);

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++)
     e[i] = exp(dirs[i*3+0]*xt[0] + dirs[i*3+1]*xt[1]+
                dirs[i*3+2]*xt[2]);

   for(int j=0;j<nldir;j++)
     for(int i=0;i<ndir;i++) {
       complex<double> t = wcb*e[i]*
                           exp(ldirs[j*3+0]*(xt[0]) +
                               ldirs[j*3+1]*(xt[1]) +
                               ldirs[j*3+2]*(xt[2]) );

       L[j*ndir+i] += t;
     }
   delete[] e;
 }
};


void DGMHelm3d::HelmDGMPMLEMatrices3d(double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc, PMLProps *pml,
                    complex<double> *K, complex<double> *L) {

 int nldir = 0;
 int nf = nFaces();
 for(int i=0;i<nf;i++) nldir += nldirs[i];
 for(int i=0;i<nldir*ndir;i++) L[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) K[i] = 0.0;

 double pmldata[4] = { kappa, pml->Rx,pml->Sx,pml->gamma};
 HelmDGMPMLEEMatrixFunction3d f(pmldata,ndir,dirs,xc,K);
 volumeInt3d(xyz, f);

 int c = 0;
 for(int faceindex=1;faceindex<=nf;faceindex++) {
/*
     HelmDGMPMLELMatrixFunction3d_ f(pmldata,ndir,dirs,
                nldirs[faceindex-1],ldirs+c*3, xsc + (faceindex-1)*3, xc,
                L+c*ndir,faceindex);
     ipu->surftInt3d(xyz, faceindex, f,gorder);
*/
     HelmDGMPMLELMatrixFunction3d f(pmldata,ndir,dirs,
                nldirs[faceindex-1],ldirs+c*3, xsc + (faceindex-1)*3, xc,
                L+c*ndir,faceindex);
     surfInt3d(xyz, faceindex, f);
     c += nldirs[faceindex-1];
 }
 for(int j=0;j<ndir;j++) {
   for(int i=0;i<j;i++)
     K[j*ndir+i] = K[i*ndir+j];
 }
}


HexDGMElement3d::HexDGMElement3d(int _n, int* nodenums) {
 o = int(pow(double(_n),1.0/3.0)+0.5);
 init(_n,nodenums);
}

int *HexDGMElement3d::faceCornerI(int fi) {
 int *fc = new int[4];
 int osq = o*o;
 int oc = osq*o;
 if (fi==1) {
   fc[0] = 1 - 1;
   fc[1] = osq-o+1 - 1;
   fc[2] = osq - 1;
   fc[3] = o - 1;
 } else if (fi==2) {
   fc[0] = 1 - 1;
   fc[1] = o - 1;
   fc[2] = osq*(o-1)+o - 1;
   fc[3] = osq*(o-1)+1 - 1;
 } else if (fi==3) {
   fc[0] = o - 1;
   fc[1] = osq - 1;
   fc[2] = oc - 1;
   fc[3] = osq*(o-1)+o - 1;
 } else if (fi==4) {
   fc[0] = osq-o+1 - 1;
   fc[1] = oc-o+1 - 1;
   fc[2] = oc - 1;
   fc[3] = osq - 1;
 } else if (fi==5) {
   fc[0] = 1 - 1;
   fc[1] = osq*(o-1)+1 - 1;
   fc[2] = oc-o+1 - 1;
   fc[3] = osq-o+1 - 1;
 } else {
   fc[0] = osq*(o-1)+1 - 1;
   fc[1] = osq*(o-1)+o - 1;
   fc[2] = oc - 1;
   fc[3] = oc-o+1 - 1;
 }
 return fc;
}

void HexDGMElement3d::externalNormal(double *xyz, int faceindex, double *n) {
 double xc[3];
 IsoParamUtils ipu(o);
 ipu.sidecenter(xyz,faceindex,n);
 ipu.elementcenter(xyz,xc);
 n[0] -= xc[0];
 n[1] -= xc[1];
 n[2] -= xc[2];
}

void HexDGMElement3d::volumeInt3d(double *xyz, IntegFunctionV3d &f) {
 IsoParamUtils ipu(o);
 ipu.volumeInt3d(xyz, f);
}

void HexDGMElement3d::surfInt3d(double *xyz, int faceindex,
                        IntegFunctionA3d &f) {
 IsoParamUtils ipu(o);
 int gorder = 7;
 ipu.surfInt3d(xyz, faceindex, f, gorder);
}
void HexDGMElement3d::surftInt3d(double *xyz, int faceindex,
                        IntegFunctionAt3d &f) {
 IsoParamUtils ipu(o);
 int gorder = 7;
 ipu.surftInt3d(xyz, faceindex, f, gorder);
}

int HexDGMElement3d::isFlatAndStraight(double *xyz,int faceindex) {
 IsoParamUtils ipu(o);
 int flag = ipu.isFlat(xyz,faceindex);
 for(int ei=0;ei<4;ei++) flag = flag && ipu.isStraight(xyz,faceindex,ei);
 return flag;
}


TetraDGMElement3d::TetraDGMElement3d(int _n, int* nodenums) {
 init(_n,nodenums);
 if (_n==4) o = 2;
 else if (_n==10) o = 3;
 else if (_n==20) o = 4;
 else if (_n==35) o = 5;
 else if (_n==56) o = 6;
 else {
   fprintf(stderr,"TetraDGMElement3d::TetraDGMElement3d: order too high\n");
   exit(-1);
 }
}

int *TetraDGMElement3d::faceCornerI(int fi) {
 int *fc = new int[3];
 int osq = ((o)*(o+1))/2;
 int oc = ((o)*(o+1)*(o+2))/6;
 if (fi==1) {
   fc[0] = o - 1;
   fc[1] = osq - 1;
   fc[2] = oc - 1;
 } else if (fi==2) {
   fc[0] = 1 - 1;
   fc[1] = oc - 1;
   fc[2] = osq - 1;
 } else if (fi==3) {
   fc[0] = 1 - 1;
   fc[1] = o - 1;
   fc[2] = oc - 1;
 } else {
   fc[0] = 1 - 1;
   fc[1] = osq - 1;
   fc[2] = o - 1;
 }
 return fc;
}

void TetraDGMElement3d::externalNormal(double *xyz, int faceindex, double *n) {
 double xc[3];
 IsoParamUtilsTetra ipu(o);
 ipu.sidecenter(xyz,faceindex,n);
 ipu.elementcenter(xyz,xc);
 n[0] -= xc[0];
 n[1] -= xc[1];
 n[2] -= xc[2];
}
void TetraDGMElement3d::volumeInt3d(double *xyz, IntegFunctionV3d &f) {
 IsoParamUtilsTetra ipu(o);
 int gorder = 13;
 ipu.volumeInt3d(xyz, f);
}

void TetraDGMElement3d::surfInt3d(double *xyz, int faceindex,
                        IntegFunctionA3d &f) {
 IsoParamUtilsTetra ipu(o);
 int gorder = 13;
 ipu.surfInt3d(xyz, faceindex, f, gorder);
}
void TetraDGMElement3d::surftInt3d(double *xyz, int faceindex,
                        IntegFunctionAt3d &f) {
 fprintf(stderr,"TetraDGMElement3d::surftInt3d is not implemented.\n");
}

int TetraDGMElement3d::isFlatAndStraight(double *xyz,int faceindex) {
 IsoParamUtilsTetra ipu(o);
// return ipu.isFlat(xyz,faceindex);
 return 0;
}


PrismDGMElement3d::PrismDGMElement3d(int _n, int* nodenums) {
 init(_n,nodenums);
 if (_n==6) o = 2;
 else if (_n==18) o = 3;
 else if (_n==40) o = 4;
 else if (_n==75) o = 5;
}

int *PrismDGMElement3d::faceCornerI(int fi) {
 int *fc = new int[nFaceCorners(fi)];
 int osq = ((o)*(o+1))/2;
 if (fi==1) {
   fc[0] = 1 - 1;
   fc[1] = osq - 1;
   fc[2] = o - 1;
 } else if (fi==2) {
   fc[0] = (o-1)*osq + 1 - 1;
   fc[1] = (o-1)*osq + o - 1;
   fc[2] = (o-1)*osq + osq - 1;
 } else if (fi==3) {
   fc[0] = o - 1;
   fc[1] = osq - 1;
   fc[2] = (o-1)*osq + osq - 1;
   fc[3] = (o-1)*osq + o - 1;
 } else if (fi==4) {
   fc[0] = osq - 1;
   fc[1] = 1 - 1;
   fc[2] = (o-1)*osq + 1 - 1;
   fc[3] = (o-1)*osq + osq - 1;
 } else {
   fc[0] = 1 - 1;
   fc[1] = o - 1;
   fc[2] = (o-1)*osq + o - 1;
   fc[3] = (o-1)*osq + 1 - 1;
 }
 return fc;
}

void PrismDGMElement3d::externalNormal(double *xyz, int faceindex, double *n) {
 double xc[3];
 IsoParamUtilsPrism ipu(o);
 ipu.sidecenter(xyz,faceindex,n);
 ipu.elementcenter(xyz,xc);
 n[0] -= xc[0];
 n[1] -= xc[1];
 n[2] -= xc[2];
}

void PrismDGMElement3d::volumeInt3d(double *xyz, IntegFunctionV3d &f) {
 IsoParamUtilsPrism ipu(o);
 int gorder = 7;
 ipu.volumeInt3d(xyz, f, gorder);
}

void PrismDGMElement3d::surfInt3d(double *xyz, int faceindex,
                        IntegFunctionA3d &f) {
 IsoParamUtilsPrism ipu(o);
 int gorder = 7;
 ipu.surfInt3d(xyz, faceindex, f, gorder);
}
void PrismDGMElement3d::surftInt3d(double *xyz, int faceindex,
                        IntegFunctionAt3d &f) {
 fprintf(stderr,"PrismDGMElement3d::surftInt3d is not implemented.\n");
}

int PrismDGMElement3d::isFlatAndStraight(double *xyz,int faceindex) {
 IsoParamUtilsPrism ipu(o);
// return ipu.isFlat(xyz,faceindex);
 return 0;
}


PyramidDGMElement3d::PyramidDGMElement3d(int _n, int* nodenums) {
 o = 2;
 init(_n,nodenums);
}

int *PyramidDGMElement3d::faceCornerI(int fi) {
 int *fc = new int[nFaceCorners(fi)];
 if (fi==1)  {
   fc[0] = 0; fc[1] = 1; fc[2] = 3; fc[3] = 2;
 } else if (fi==2) {
   fc[0] = 0; fc[1] = 1; fc[2] = 4;
 } else if (fi==3) {
   fc[0] = 1; fc[1] = 3; fc[2] = 4;
 } else if (fi==4) {
   fc[0] = 3; fc[1] = 2; fc[2] = 4;
 } else {
   fc[0] = 2; fc[1] = 0; fc[2] = 4;
 }
 return fc;
}

void PyramidDGMElement3d::externalNormal(double *xyz, int faceindex, double *n) {
 double xc[3];
 IsoParamUtilsPyramid ipu(o);
 ipu.sidecenter(xyz,faceindex,n);
 ipu.elementcenter(xyz,xc);
 n[0] -= xc[0];
 n[1] -= xc[1];
 n[2] -= xc[2];
}


void PyramidDGMElement3d::volumeInt3d(double *xyz, IntegFunctionV3d &f) {
 fprintf(stderr,"PyramidDGMElement3d::volumeInt3dis not implemented.\n");
}

void PyramidDGMElement3d::surfInt3d(double *xyz, int faceindex,
                        IntegFunctionA3d &f) {
 IsoParamUtilsPyramid ipu(o);
 int gorder = 7;
 ipu.surfInt3d(xyz, faceindex, f, gorder);
}

void PyramidDGMElement3d::surftInt3d(double *xyz, int faceindex,
                        IntegFunctionAt3d &f) {
 fprintf(stderr,"PyramidDGMElement3d::surftInt3d is not implemented.\n");
}

int PyramidDGMElement3d::isFlatAndStraight(double *xyz,int faceindex) {
 IsoParamUtilsPyramid ipu(o);
// return ipu.isFlat(xyz,faceindex);
 return 0;
}


DGMHelm3d_6::DGMHelm3d_6(int _nnodes, int* nodenums) :
  HexDGMElement3d(_nnodes,nodenums) {
 ndir = 6;
}

DGMHelm3d_6t::DGMHelm3d_6t(int _nnodes, int* nodenums) :
  TetraDGMElement3d(_nnodes,nodenums) {
 ndir = 6;
}

DGMHelm3d_6p::DGMHelm3d_6p(int _nnodes, int* nodenums) :
  PrismDGMElement3d(_nnodes,nodenums) {
 ndir = 6;
}

DGMHelm3d_6pd::DGMHelm3d_6pd(int _nnodes, int* nodenums) :
  PyramidDGMElement3d(_nnodes,nodenums) {
 ndir = 6;
}

DGMHelm3d_26::DGMHelm3d_26(int _nnodes, int* nodenums) :
  HexDGMElement3d(_nnodes,nodenums) {
 ndir = 26;
}

DGMHelm3d_26t::DGMHelm3d_26t(int _nnodes, int* nodenums) :
  TetraDGMElement3d(_nnodes,nodenums) {
 ndir = 26;
}

DGMHelm3d_26p::DGMHelm3d_26p(int _nnodes, int* nodenums) :
  PrismDGMElement3d(_nnodes,nodenums) {
 ndir = 26;
}

DGMHelm3d_26pd::DGMHelm3d_26pd(int _nnodes, int* nodenums) :
  PyramidDGMElement3d(_nnodes,nodenums) {
 ndir = 26;
}

DGMHelm3d_56::DGMHelm3d_56(int _nnodes, int* nodenums) :
  HexDGMElement3d(_nnodes,nodenums) {
 ndir = 56;
}

DGMHelm3d_56t::DGMHelm3d_56t(int _nnodes, int* nodenums) :
  TetraDGMElement3d(_nnodes,nodenums) {
 ndir = 56;
}

DGMHelm3d_56p::DGMHelm3d_56p(int _nnodes, int* nodenums) :
  PrismDGMElement3d(_nnodes,nodenums) {
 ndir = 56;
}

DGMHelm3d_56pd::DGMHelm3d_56pd(int _nnodes, int* nodenums) :
  PyramidDGMElement3d(_nnodes,nodenums) {
 ndir = 56;
}

DGMHelm3d_98::DGMHelm3d_98(int _nnodes, int* nodenums) :
  HexDGMElement3d(_nnodes,nodenums) {
 ndir = 98;
}

void DGMHelm3d_6::dir(int n, complex<double> *d) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1}
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

void DGMHelm3d_6t::dir(int n, complex<double> *d) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1}
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

void DGMHelm3d_6p::dir(int n, complex<double> *d) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1} 
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

void DGMHelm3d_6pd::dir(int n, complex<double> *d) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1} 
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

double* DGMHelm3d::getCubeDir(int n) {
 static double *a[3] = {0,0,0}; 
 if (a[n-2]==0) {
   a[n-2] = new double[(6*n*n+2)*3];
   int c = 0;
   for(int kk=0;kk<=n;kk++) for (int jj=0;jj<=n;jj++) for(int ii=0;ii<=n;ii++) {
     double cube[3];
     if (kk==0 || kk==n) {
       cube[0] = 0.5*tan(M_PI/4.0*(-1.0+2.0*ii/double(n)));
       cube[1] = 0.5*tan(M_PI/4.0*(-1.0+2.0*jj/double(n)));
       cube[2] = kk/double(n) -0.5;
     } else if (jj==0 || jj==n) {
       cube[0] = 0.5*tan(M_PI/4.0*(-1.0+2.0*ii/double(n)));
       cube[1] = jj/double(n) -0.5;
       cube[2] = 0.5*tan(M_PI/4.0*(-1.0+2.0*kk/double(n)));
     } else if (ii==0 || ii==n) {
       cube[0] = ii/double(n) -0.5;
       cube[1] = 0.5*tan(M_PI/4.0*(-1.0+2.0*jj/double(n)));
       cube[2] = 0.5*tan(M_PI/4.0*(-1.0+2.0*kk/double(n)));
     } else continue;
     double l = sqrt(cube[0]*cube[0]+cube[1]*cube[1]+cube[2]*cube[2]);
     a[n-2][c*3+0] = cube[0]/l;
     a[n-2][c*3+1] = cube[1]/l;
     a[n-2][c*3+2] = cube[2]/l;
     c++;
   }
 } 
 return a[n-2];
}

void DGMHelm3d_26::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_26t::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_26p::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_26pd::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_56::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_56t::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_56p::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_56pd::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_98::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(4);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}


void DGMHelm3d_1_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = std::max( e1->getOmega()/e1->getSpeedOfSound(),
                      (e2!=0)?e2->getOmega()/e2->getSpeedOfSound():0.0);

 d[0] = 0.0;
 d[1] = 0.0;
 d[2] = 0.0;
}

void DGMHelm3d_4_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = std::max( e1->getOmega()/e1->getSpeedOfSound(),
                      (e2!=0)?e2->getOmega()/e2->getSpeedOfSound():0.0);
 double a[][2] = { 
    { 1.0/sqrt(2.0), 1.0/sqrt(2.0) },
    { -1.0/sqrt(2.0), 1.0/sqrt(2.0) }, 
    { 1.0/sqrt(2.0), -1.0/sqrt(2.0) },
    { -1.0/sqrt(2.0), -1.0/sqrt(2.0) } };
 d[0] = complex<double>(0.0,0.6*lkappa*(a[n][0]*tau1[0]+a[n][1]*tau2[0]));
 d[1] = complex<double>(0.0,0.6*lkappa*(a[n][0]*tau1[1]+a[n][1]*tau2[1]));
 d[2] = complex<double>(0.0,0.6*lkappa*(a[n][0]*tau1[2]+a[n][1]*tau2[2]));
}

void DGMHelm3d_8_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = std::max( e1->getOmega()/e1->getSpeedOfSound(),
                      (e2!=0)?e2->getOmega()/e2->getSpeedOfSound():0.0);
 if (n<4) {
   double a[][2] = { 
      { 1.0/sqrt(2.0), 1.0/sqrt(2.0) },
      { -1.0/sqrt(2.0), 1.0/sqrt(2.0) }, 
      { 1.0/sqrt(2.0), -1.0/sqrt(2.0) },
      { -1.0/sqrt(2.0), -1.0/sqrt(2.0) } };
   d[0] = complex<double>(0.0,0.8*lkappa*(a[n][0]*tau1[0]+a[n][1]*tau2[0]));
   d[1] = complex<double>(0.0,0.8*lkappa*(a[n][0]*tau1[1]+a[n][1]*tau2[1]));
   d[2] = complex<double>(0.0,0.8*lkappa*(a[n][0]*tau1[2]+a[n][1]*tau2[2]));
 } else {
   n -= 4;
   double a[][2] = { 
      { 1,0} , { -1,0}, {0,1}, {0,-1} };
   d[0] = complex<double>(0.0,0.5*lkappa*(a[n][0]*tau1[0]+a[n][1]*tau2[0]));
   d[1] = complex<double>(0.0,0.5*lkappa*(a[n][0]*tau1[1]+a[n][1]*tau2[1]));
   d[2] = complex<double>(0.0,0.5*lkappa*(a[n][0]*tau1[2]+a[n][1]*tau2[2]));
 }
}

void DGMHelm3d_12_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = std::max( e1->getOmega()/e1->getSpeedOfSound(),
                      (e2!=0)?e2->getOmega()/e2->getSpeedOfSound():0.0);
 if (n<6) {
   double c = cos(n*M_PI/6.0+M_PI/12.0);
   double s = sin(n*M_PI/6.0+M_PI/12.0);
   d[0] = complex<double>(0.0,0.8*lkappa*(c*tau1[0]+s*tau2[0]));
   d[1] = complex<double>(0.0,0.8*lkappa*(c*tau1[1]+s*tau2[1]));
   d[2] = complex<double>(0.0,0.8*lkappa*(c*tau1[2]+s*tau2[2]));
 } else {
   double c = cos(n*M_PI/6.0);
   double s = sin(n*M_PI/6.0);
   d[0] = complex<double>(0.0,0.5*lkappa*(c*tau1[0]+s*tau2[0]));
   d[1] = complex<double>(0.0,0.5*lkappa*(c*tau1[1]+s*tau2[1]));
   d[2] = complex<double>(0.0,0.5*lkappa*(c*tau1[2]+s*tau2[2]));
 }
}


void DGMHelm3d::init(int _nnodes, int* nodenums) {
 if (_nnodes<0) {
   _nnodes=-_nnodes;
 }
 nn = new int[_nnodes]; 
 for(int i=0;i<_nnodes;i++) nn[i] = nodenums[i];
 lm = new DEMLM*[nFaces()];
 for(int i=0;i<nFaces();i++) lm[i] = 0;
 bc = new int[nFaces()];
 for(int i=0;i<nFaces();i++) bc[i] = 0;
 ndir = 0;
 storeMatricesF = true;
 condensedF = true;
}


void DGMHelm3d::getRef(double *xyz,double *cxyz) {
/*
 if (o>0) {
   IsoParamUtils ipu(o);
   ipu.elementcenter(xyz,cxyz);  
 } else {
   IsoParamUtilsTetra ipu(-o);
   ipu.elementcenter(xyz,cxyz);  
 }
*/
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
}

void DGMHelm3d::createM(complex<double>*M) {

// IsoParamUtils *ipu = (o>0)? new IsoParamUtils(o):new IsoParamUtilsTetra(-o);
// int os = ipu->getordersq();
 int oc = nGeomNodes();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 int nldir[6] = {0,0,0,0,0,0};
 int nf = nFaces();
 int tnldir = 0;
 for(int fi=0;fi<nf;fi++) if (lm[fi]!=0) {
   nldir[fi] = lm[fi]->nDofs();
   tnldir += nldir[fi];
   DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[fi]);
//   DGMHelm3d_Eva_LM *le = dynamic_cast<DGMHelm3d_Eva_LM*>(lm[fi]);
//   if (l==0 && le==0) {
   if (l==0) {
     fprintf(stderr,"DGMHelm3d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMHelm3d_LM type specific
 double xlref[18] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<nf;fi++) 
// ipu->sidecenter(xyz,fi+1,xlref+fi*3);

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];

 complex<double>* cldir = new complex<double>[3*tnldir];
/*
 int corner[6] = { 1-1, o-1, os-1, os-o+1-1, 1-1, o-1};
 if (o<0) {
   corner[1] = -o-1;
   corner[3] = 1-1;
   corner[4] = -o-1;
 }
 int *fi = new int[os];
 int cc = 0;
 double sign[6] = { 1,1,1,1,1,1};
 for(int i=0;i<nf;i++) {
   ipu->faceindeces(i+1,fi);
   int cmin = 1;
   for(int j=2;j<=nFaceCorners(i+1);j++)
     if ( nn[fi[corner[j]]]<nn[fi[corner[cmin]]]) cmin = j;
   double tau1[3], tau2[3];
   if (nn[fi[corner[cmin-1]]] < nn[fi[corner[cmin+1]]]) {
     tau1[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau2[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   } else {
     tau2[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau1[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   }
   double xsc[3];
   double xc[3];
   ipu->sidecenter(xyz,i+1,xsc);
   ipu->elementcenter(xyz,xc);
   xsc[0] -= xc[0];
   xsc[1] -= xc[1];
   xsc[2] -= xc[2];
   double cross[3] = {
        tau1[1]*tau2[2]-tau1[2]*tau2[1],
        tau1[2]*tau2[0]-tau1[0]*tau2[2],
        tau1[0]*tau2[1]-tau1[1]*tau2[0]};
   if (cross[0]*xsc[0]+cross[1]*xsc[1]+cross[2]*xsc[2]>0.0) sign[i] = -1.0;
   double l = sqrt(tau1[0]*tau1[0]+tau1[1]*tau1[1]+tau1[2]*tau1[2]);
   tau1[0] /= l;
   tau1[1] /= l;
   tau1[2] /= l;
   double p = tau1[0]*tau2[0]+tau1[1]*tau2[1]+tau1[2]*tau2[2];
   tau2[0] -= p*tau1[0];
   tau2[1] -= p*tau1[1];
   tau2[2] -= p*tau1[2];
   l = sqrt(tau2[0]*tau2[0]+tau2[1]*tau2[1]+tau2[2]*tau2[2]);
   tau2[0] /= l;
   tau2[1] /= l;
   tau2[2] /= l;
   if (nldir[i]>0) {
     DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[i]);
     for(int j=0;j<lm[i]->nDofs();j++) {
       l->ldir(j,tau1,tau2,cldir+3*cc);
       cc++;
     }
   }
 }
 delete[] fi;
*/
 int corner[6];
 int cc = 0;
 double sign[6] = { 1,1,1,1,1,1};
 for(int i=0;i<nf;i++) {
   int *fc = faceCornerI(i+1);
   for(int j=0;j<nFaceCorners(i+1);j++) corner[j] = fc[j];
   corner[nFaceCorners(i+1)] = corner[0];
   corner[nFaceCorners(i+1)+1] = corner[1];
   delete[] fc;
   int cmin = 1;
   for(int j=2;j<=nFaceCorners(i+1);j++)
     if ( nn[corner[j]]<nn[corner[cmin]]) cmin = j;
   double tau1[3], tau2[3];
   if (nn[corner[cmin-1]] < nn[corner[cmin+1]]) {
     tau1[0] = xyz[0*oc+corner[cmin-1]]-xyz[0*oc+corner[cmin]];
     tau1[1] = xyz[1*oc+corner[cmin-1]]-xyz[1*oc+corner[cmin]];
     tau1[2] = xyz[2*oc+corner[cmin-1]]-xyz[2*oc+corner[cmin]];
     tau2[0] = xyz[0*oc+corner[cmin+1]]-xyz[0*oc+corner[cmin]];
     tau2[1] = xyz[1*oc+corner[cmin+1]]-xyz[1*oc+corner[cmin]];
     tau2[2] = xyz[2*oc+corner[cmin+1]]-xyz[2*oc+corner[cmin]];
   } else {
     tau2[0] = xyz[0*oc+corner[cmin-1]]-xyz[0*oc+corner[cmin]];
     tau2[1] = xyz[1*oc+corner[cmin-1]]-xyz[1*oc+corner[cmin]];
     tau2[2] = xyz[2*oc+corner[cmin-1]]-xyz[2*oc+corner[cmin]];
     tau1[0] = xyz[0*oc+corner[cmin+1]]-xyz[0*oc+corner[cmin]];
     tau1[1] = xyz[1*oc+corner[cmin+1]]-xyz[1*oc+corner[cmin]];
     tau1[2] = xyz[2*oc+corner[cmin+1]]-xyz[2*oc+corner[cmin]];
   }
   double xsc[3];
//   double xc[3];
//   ipu->sidecenter(xyz,i+1,xsc);
//   ipu->elementcenter(xyz,xc);
   externalNormal(xyz,i+1,xsc);
   double cross[3] = {
        tau1[1]*tau2[2]-tau1[2]*tau2[1],
        tau1[2]*tau2[0]-tau1[0]*tau2[2],
        tau1[0]*tau2[1]-tau1[1]*tau2[0]};
   if (cross[0]*xsc[0]+cross[1]*xsc[1]+cross[2]*xsc[2]>0.0) sign[i] = -1.0;
   double l = sqrt(tau1[0]*tau1[0]+tau1[1]*tau1[1]+tau1[2]*tau1[2]);
   tau1[0] /= l;
   tau1[1] /= l;
   tau1[2] /= l;
   double p = tau1[0]*tau2[0]+tau1[1]*tau2[1]+tau1[2]*tau2[2];
   tau2[0] -= p*tau1[0];
   tau2[1] -= p*tau1[1];
   tau2[2] -= p*tau1[2];
   l = sqrt(tau2[0]*tau2[0]+tau2[1]*tau2[1]+tau2[2]*tau2[2]);
   tau2[0] /= l;
   tau2[1] /= l;
   tau2[2] /= l;
   if (nldir[i]>0) {
     DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[i]);
/*
fprintf(stderr,"kuku %f %f %f  %f %f %f  %f\n",
tau1[0], tau1[1], tau1[2],
tau2[0], tau2[1], tau2[2], sign[i]);
*/
     for(int j=0;j<lm[i]->nDofs();j++) {
       l->ldir(j,tau1,tau2,cldir+3*cc);
       cc++;
     }
   }
 }

 int arbFlag[6] = { 0,0,0,0,0,0};
 for(int i=0;i<nf;i++) if (bc[i]==1) arbFlag[i] = 1;


 PMLProps* pml = getPMLProps();
 if (pml->PMLtype==0)
   HelmDGMEMatrices3d(xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, kee, kel); 
 else
   HelmDGMPMLEMatrices3d(xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, pml, kee, kel); 

 cc = 0;
 for(int i=0;i<nf;i++) {
  for(int j=cc;j<cc+nldir[i];j++) for(int k=0;k<ndir;k++)
    kel[j*ndir+k] *= sign[i];
  cc += nldir[i];
 }
/*
 for(int i=0;i<ndir*ndir;i++) fprintf(stderr,"K %d %d: %e %e\n",i/ndir,i%ndir,real(kee[i]),imag(kee[i]));
 for(int i=0;i<tnldir*ndir;i++) fprintf(stderr,"L %d %d: %e %e\n",i/ndir,i%ndir,real(kel[i]),imag(kel[i]));
*/

 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(tnldir+j)*(ndir+tnldir)+(tnldir+i)] = kee[j*ndir+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(tnldir+i)*(ndir+tnldir)+j] = kel[j*ndir+i]/rho;
   M[j*(ndir+tnldir)+(tnldir+i)] = kel[j*ndir+i]/rho;
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++) 
   M[j*(ndir+tnldir)+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] xyz;
 delete[] cldir;
 delete[] cdir;
// delete ipu;
}


class HelmDGMWetEMatrixFunction3d : public IntegFunctionA3d {
 double omega;
 int ndir;
 complex<double> *dirs;
 double *xc;
 DEMElement *deme2;
 complex<double> *K;
public:
 HelmDGMWetEMatrixFunction3d(double _omega,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      DEMElement *_deme2,
                      complex<double> *_K) {
   omega = _omega; ndir = _ndir; dirs = _dirs; xc = _xc; deme2 = _deme2;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
/*
   for(int i=0;i<ndir;i++) { e[i] = omega*omega*w*
    exp(dirs[i*3+0]*(x[0]-xc[0]) + dirs[i*3+1]*(x[1]-xc[1]) + 
        dirs[i*3+2]*(x[2]-xc[2]) );
   }

   int ndir2 = deme2->nEnrichmentDofs();
   complex<double> *e2 = new complex<double>[ndir2*3];
   deme2->enrichmentF(x,e2);
   for(int j=0;j<ndir2;j++) {
     for(int i=0;i<ndir;i++) {
       complex<double> v = nsign*e[i]*
          (cross[0]*e2[3*j+0]+cross[1]*e2[3*j+1]+cross[2]*e2[3*j+2]);
       K[(j+ndir)*(ndir+ndir2)+i] += v;
       K[i*(ndir+ndir2)+j+ndir] += v;
     }
   }
   delete[] e2;
   delete[] e;
*/

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) {
       e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                  dirs[i*3+1]*(x[1]-xc[1]) +
                  dirs[i*3+2]*(x[2]-xc[2]) );
   }

   int oc2 = deme2->nPolynomialDofs();
   int ndir2 = deme2->nEnrichmentDofs();
   int ndof = ndir + oc2 + ndir2;
   double oown = omega*omega*w*nsign;

   complex<double> *e2 = new complex<double>[ndir2*3];
   deme2->enrichmentF(x,e2); 

   for(int j=0;j<ndir2;j++) {
     complex<double> oownc = oown*
       (cross[0]*e2[j*3+0]+cross[1]*e2[j*3+1]+cross[2]*e2[j*3+2]);
     for(int i=0;i<ndir;i++) {
       K[(ndir+oc2+j)*ndof+i] += oownc * e[i];
       K[(i)*ndof+ndir+oc2+j] += oownc * e[i];
     }
   }
   delete[] e2;

   if (oc2>0) {
     double *u = new double[oc2*3];
    deme2->polynomialF(x,u);

     for(int j=0;j<oc2;j++) {
       double oownc = oown*
          (cross[0]*u[3*j+0]+cross[1]*u[3*j+1]+cross[2]*u[3*j+2]);
       for(int i=0;i<ndir;i++) {
         K[(ndir+j)*ndof+i] += oownc * e[i];
         K[(i)*ndof+ndir+j] += oownc * e[i];
       }
     }
     delete[] u;
   }
   delete[] e;
 }
};



void DGMHelm3d::interfMatrix(int fi, DEMElement* deme2,
                              complex<double> *K) {

 int oc = nGeomNodes();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double omega = getOmega();

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 int ndofs = ndir+deme2->nPolynomialDofs()+deme2->nEnrichmentDofs();
 for(int i=0;i<ndofs*ndofs;i++) K[i] = 0.0;

 HelmDGMWetEMatrixFunction3d f(omega,ndir,cdir,xref,deme2,K);
 surfInt3d(xyz, fi, f);

 delete[] cdir;
}


class HelmDGMENeumVFunction3d : public IntegFunctionA3d {
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
public:
 HelmDGMENeumVFunction3d(
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc,
                      complex<double> *_v) {
   ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int i;
   complex<double> *e = new complex<double>[ndir];
   for(i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                  dirs[i*3+1]*(x[1]-xc[1]) +
                                  dirs[i*3+2]*(x[2]-xc[2]));

   complex<double> ince = w* (nsign*
        (incdir[0]*cross[0]+incdir[1]*cross[1]+incdir[2]*cross[2])) *
        exp(incdir[0]*x[0] + incdir[1]*x[1] + incdir[2]*x[2]);

   for(i=0;i<ndir;i++)
     v[i] += ince*e[i];
   delete[] e; 
 }
};


void DGMHelm3d::HelmDGMENeumV(double *xyz,
                     int ndir, complex<double> *dirs,
                     double kappa, complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {

 HelmDGMENeumVFunction3d f(ndir,dirs,incdir,xc,v);
 surfInt3d(xyz,faceindex, f);
}

class HelmDGMNeumFVFunction3d : public IntegFunctionA3d {
 int order;
 int ndir;
 complex<double> *dirs;
 complex<double> *f;
 double *xc;
 complex<double> *v;
public:
 HelmDGMNeumFVFunction3d(int _order,
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_f, double *_xc,
                      complex<double> *_v) {
   order = _order; ndir = _ndir; dirs = _dirs; f = _f; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   int oc = order*order*order;
  
   double wc = w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]); 
   complex<double> ff = 0.0;
   for(int i=0;i<oc;i++) ff = f[i]*N[i];

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1])+
                                      dirs[i*3+2]*(x[2]-xc[2]));
   for(int i=0;i<ndir;i++)
     v[i] += wc*ff*e[i];
   delete[] e;
 }
};


void HelmDGMNeumFV(int order, double *xy,
                     int ndir, complex<double> *dirs,
                     complex<double> *fv, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils ipu(order);
 int gorder = (order>0)?7:13;

 HelmDGMNeumFVFunction3d f(order,ndir,dirs,fv,xc,v);
 ipu.surfInt3d(xy, faceindex, f, gorder);
}

class HelmDGMERobVFunction3d : public IntegFunctionA3d {
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
 double kappa;
public:
 HelmDGMERobVFunction3d(
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc, double _kappa,
                      complex<double> *_v) {
   ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc; kappa = _kappa;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int i;
   complex<double> *e = new complex<double>[ndir];
   for(i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                  dirs[i*3+1]*(x[1]-xc[1]) +
                                  dirs[i*3+2]*(x[2]-xc[2]));

   complex<double> ince = w* ((nsign*
        (incdir[0]*cross[0]+incdir[1]*cross[1]+incdir[2]*cross[2])) -
complex<double>(0.0,kappa*
 sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2])) ) *
        exp(incdir[0]*x[0] + incdir[1]*x[1] + incdir[2]*x[2]);

   for(i=0;i<ndir;i++)
     v[i] += ince*e[i];
   delete[] e; 
 }
};


void HelmDGMERobV(int order, double *xyz,
                     int ndir, complex<double> *dirs,
                     double kappa, complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils *ipu = (order>0)? new IsoParamUtils(order):new IsoParamUtilsTetra(-order);
 int gorder = (order>0)?7:13;

 HelmDGMERobVFunction3d f(ndir,dirs,incdir,xc,kappa,v);
 ipu->surfInt3d(xyz,faceindex, f, gorder);
 delete ipu;
}



void DGMHelm3d::createRHS(complex<double>*v) {
 int oc = nGeomNodes();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 complex<double> *vv = new complex<double>[ndir];
 for(int i=0;i<ndir;i++) vv[i] = 0.0;

 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2) {
     double *dir = getWaveDirection();
     complex<double> incdir[3] = {complex<double>(0.0,kappa*dir[0]),
                                  complex<double>(0.0,kappa*dir[1]),
                                  complex<double>(0.0,kappa*dir[2])};
     HelmDGMENeumV(xyz, ndir, cdir, kappa, incdir, i+1 , xref, vv);
   } else if (bc[i]==4) {
     complex<double> *f = getField();
     HelmDGMNeumFV(o, xyz, ndir, cdir, f, i+1, xref, vv);
   }
//   else if (bc[i]==1)
//     HelmDGMERobV(o, xyz, ndir, cdir, kappa, incdir, i+1 , xref, vv);
 }

// for(int i=0;i<ndir;i++) 
//   fprintf(stderr,"rhs %d %e %e\n",i,real(vv[i]),imag(vv[i]));


 delete[] cdir;
 delete[] xyz;

 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[i] = 0;
 for(int i=0;i<ndir;i++) 
   v[tnldir+i] = vv[i]/rho;
 delete[] vv;
}


void DGMHelm3d::createSol(double *xyz,
                        complex<double>* sol, complex<double> *nodalSol) {

 for(int i=0;i<8;i++) nodalSol[i] = 0.0;
 for(int i=0;i<ndir;i++) {
   complex<double> d[3];
   dir(i,d);
   nodalSol[0] += sol[i]*exp( d[0]*xyz[0]+d[1]*xyz[1]+d[2]*xyz[2] );
 }
}


DEMHelm3d::DEMHelm3d(int _nnodes, int* nodenums) :
  HexDGMElement3d(_nnodes, nodenums) {}


DEMHelm3d_6::DEMHelm3d_6(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 6;
}

DEMHelm3d_26::DEMHelm3d_26(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 26;
}

DEMHelm3d_56::DEMHelm3d_56(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 56;
}

DEMHelm3d_98::DEMHelm3d_98(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 98;
}

void DEMHelm3d_6::dir(int n, complex<double> *d) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1}
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

void DEMHelm3d_26::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DEMHelm3d_56::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DEMHelm3d_98::dir(int n, complex<double>* d) {
 double kappa = getOmega()/getSpeedOfSound();
 double *a = getCubeDir(4);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}


class HelmPMatricesFunction3d : public IntegFunctionV3d {
 int o;
 int ndir;
 complex<double> *dirs;
 double *xc;
 double kappa;
 complex<double> *PP;
 complex<double> *PE;
public:
 HelmPMatricesFunction3d(int _o,  double _kappa,
                         int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_PE, complex<double> *_PP) {
   o = _o; kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   PP = _PP; PE = _PE;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
    
   int oc = o*o*o;

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1])+
                                      dirs[i*3+2]*(x[2]-xc[2]));
   for(int j=0;j<ndir;j++)
     for(int i=0;i<oc;i++)
       PE[j*oc+i] += w*det*e[j]*(-kappa*kappa*N[i] +
     dNdx[i][0]*dirs[j*3+0]+dNdx[i][1]*dirs[j*3+1]+dNdx[i][2]*dirs[j*3+2]);
   delete[] e;
   for(int j=0;j<oc;j++)
     for(int i=j;i<oc;i++)
       PP[j*oc+i] += w*det*(-kappa*kappa*N[i]*N[j]
        +dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]+dNdx[i][2]*dNdx[j][2]);
 }
};


void HelmDEMMatrices3d(int order, double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc,
                    complex<double> *kee, complex<double> *kel,
                    complex<double> *kpp, complex<double> *kpl,
                    complex<double> *kpe) {

 IsoParamUtils ipu(order);
 int oc = ipu.getorderc();

 int nldir =nldirs[0] + nldirs[1] + nldirs[2] +
            nldirs[3] + nldirs[4] + nldirs[5];
 for(int i=0;i<nldir*ndir;i++) kel[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) kee[i] = 0.0;
 for(int i=0;i<oc*ndir;i++) kpe[i] = 0.0;
 for(int i=0;i<oc*nldir;i++) kpl[i] = 0.0;
 for(int i=0;i<oc*oc;i++) kpp[i] = 0.0;

 int c = 0;
 for(int faceindex=1;faceindex<=6;faceindex++) {
   HelmDGMEMatricesFunction3d f(order,kappa,ndir,dirs,
                nldirs[faceindex-1],ldirs+c*3,
                xsc + (faceindex-1)*3, xc, sflags[faceindex-1], faceindex,
                kee,kel+c*ndir,kpl+c*oc,kpp,kpe);
   ipu.surfInt3d(xyz, faceindex, f);
   c += nldirs[faceindex-1];
 }
 for(int j=0;j<ndir;j++) {
   for(int i=0;i<j;i++)
     kee[j*ndir+i] = kee[i*ndir+j];
 }

 HelmPMatricesFunction3d f(order,kappa,ndir,dirs, xc, kpe,kpp);
 ipu.volumeInt3d(xyz, f);

 for(int j=0;j<oc;j++) {
   for(int i=0;i<j;i++)
     kpp[j*oc+i] = kpp[i*oc+j];
 }
}


void DEMHelm3d::createM(complex<double>*M) {

 IsoParamUtils ipu(o);
 int os = ipu.getordersq();
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 int nldir[6] = {0,0,0,0,0,0};
 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<6;i++) tnldir += nldir[i];

 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) {
   DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[fi]);
   if (l==0) {
     fprintf(stderr,"DGMHelm3d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMHelm3d_LM type specific
 double xlref[18] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<6;fi++) 
// ipu.sidecenter(xyz,fi+1,xlref+fi*3);

 complex<double>* cldir = new complex<double>[3*tnldir];
 int corner[6] = { 1-1, o-1, os-1, os-o+1-1, 1-1, o-1};
 int *fi = new int[os];
 int cc = 0;
 double sign[6] = { 1,1,1,1,1,1};
 for(int i=0;i<6;i++) {
   ipu.faceindeces(i+1,fi);
   int cmin = 1;
   for(int j=2;j<5;j++) if ( nn[fi[corner[j]]]<nn[fi[corner[cmin]]]) cmin = j;
   double tau1[3], tau2[3];
   if (nn[fi[corner[cmin-1]]] < nn[fi[corner[cmin+1]]]) {
     tau1[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau2[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   } else {
     tau2[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau1[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   }
   double xsc[3];
   double xc[3];
   ipu.sidecenter(xyz,i+1,xsc);
   ipu.elementcenter(xyz,xc);
   xsc[0] -= xc[0];
   xsc[1] -= xc[1];
   xsc[2] -= xc[2];
   double cross[3] = {
        tau1[1]*tau2[2]-tau1[2]*tau2[1],
        tau1[2]*tau2[0]-tau1[0]*tau2[2],
        tau1[0]*tau2[1]-tau1[1]*tau2[0]};
   if (cross[0]*xsc[0]+cross[1]*xsc[1]+cross[2]*xsc[2]>0.0) sign[i] = -1.0;
   double l = sqrt(tau1[0]*tau1[0]+tau1[1]*tau1[1]+tau1[2]*tau1[2]);
   tau1[0] /= l;
   tau1[1] /= l;
   tau1[2] /= l;
   double p = tau1[0]*tau2[0]+tau1[1]*tau2[1]+tau1[2]*tau2[2];
   tau2[0] -= p*tau1[0];
   tau2[1] -= p*tau1[1];
   tau2[2] -= p*tau1[2];
   l = sqrt(tau2[0]*tau2[0]+tau2[1]*tau2[1]+tau2[2]*tau2[2]);
   tau2[0] /= l;
   tau2[1] /= l;
   tau2[2] /= l;
   if (nldir[i]>0) {
     DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[i]);
     for(int j=0;j<lm[i]->nDofs();j++) {
       l->ldir(j,tau1,tau2,cldir+3*cc);
       cc++;
     }
   }
 }
 delete[] fi;

 int arbFlag[6] = { 0,0,0,0,0,0};
 for(int i=0;i<6;i++) if (bc[i]==1) arbFlag[i] = 1;

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];
 complex<double>* kpp = new complex<double>[oc*oc];
 complex<double>* kpe = new complex<double>[ndir*oc];
 complex<double>* kpl = new complex<double>[oc*tnldir];
 HelmDEMMatrices3d(o, xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, kee, kel, kpp, kpl, kpe); 

 cc = 0;
 for(int i=0;i<6;i++) {
  for(int j=cc;j<cc+nldir[i];j++) {
    for(int k=0;k<ndir;k++) kel[j*ndir+k] *= sign[i];
    for(int k=0;k<oc;k++) kpl[j*oc+k] *= sign[i];
  }
  cc += nldir[i];
 }

 for(int i=0;i<oc;i++) for(int j=0;j<oc;j++) M[j*(oc+tnldir+ndir)+i] = kpp[j*oc+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(oc+tnldir+j)*(oc+ndir+tnldir)+(oc+tnldir+i)] = kee[j*ndir+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(oc+tnldir+i)*(oc+ndir+tnldir)+oc+j] = kel[j*ndir+i]/rho;
   M[(oc+j)*(oc+ndir+tnldir)+(oc+tnldir+i)] = kel[j*ndir+i]/rho;
 }
 for(int i=0;i<oc;i++) for(int j=0;j<ndir;j++) {
   M[i*(oc+ndir+tnldir)+(oc+tnldir+j)] = kpe[j*oc+i]/rho;
   M[(oc+tnldir+j)*(oc+ndir+tnldir)+i] = kpe[j*oc+i]/rho;
 }
 for(int i=0;i<oc;i++) for(int j=0;j<tnldir;j++) {
   M[i*(oc+ndir+tnldir)+(oc+j)] = kpl[j*oc+i]/rho;
   M[(oc+j)*(oc+ndir+tnldir)+i] = kpl[j*oc+i]/rho;
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++)
   M[(oc+j)*(oc+ndir+tnldir)+oc+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] kpp;
 delete[] kpe;
 delete[] kpl;
 delete[] xyz;
 delete[] cldir;
 delete[] cdir;
}


class HelmDEMWetMatrixFunction3d : public IntegFunctionA3d {
 int o;
 double omega;
 int ndir;
 complex<double> *dirs;
 double *xc;
 DEMElement *deme2;
 complex<double> *K;
public:
 HelmDEMWetMatrixFunction3d(int _o, double _omega,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      DEMElement *_deme2,
                      complex<double> *_K) {
   o= _o; omega = _omega; ndir = _ndir; dirs = _dirs; xc = _xc; deme2 = _deme2;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) {
       e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                  dirs[i*3+1]*(x[1]-xc[1]) +
                  dirs[i*3+2]*(x[2]-xc[2]) );
   }

   int oc = o*o*o;
   int oc2 = deme2->nPolynomialDofs();
   int ndir2 = deme2->nEnrichmentDofs();
   int ndof = oc + ndir + oc2 + ndir2;
   double oown = omega*omega*w*nsign;

   complex<double> *e2 = new complex<double>[ndir2*3];
   deme2->enrichmentF(x,e2); 

   for(int j=0;j<ndir2;j++) {
     complex<double> oownc = oown*
       (cross[0]*e2[j*3+0]+cross[1]*e2[j*3+1]+cross[2]*e2[j*3+2]);
     for(int i=0;i<ndir;i++) {
       K[(oc+ndir+oc2+j)*ndof+oc+i] += oownc * e[i];
       K[(oc+i)*ndof+oc+ndir+oc2+j] += oownc * e[i];
     }
     for(int i=0;i<oc;i++) {
       K[(oc+ndir+oc2+j)*ndof+i] += oownc * N[i];
       K[(i)*ndof+oc+ndir+oc2+j] += oownc * N[i];
     }
   }
   delete[] e2;

   double *u = new double[oc2*3];
   deme2->polynomialF(x,u);

   for(int j=0;j<oc2;j++) {
     double oownc = oown*
        (cross[0]*u[3*j+0]+cross[1]*u[3*j+1]+cross[2]*u[3*j+2]);
     for(int i=0;i<ndir;i++) {
       K[(oc+ndir+j)*ndof+oc+i] += oownc * e[i];
       K[(oc+i)*ndof+oc+ndir+j] += oownc * e[i];
     }
     for(int i=0;i<oc;i++) {
       K[(oc+ndir+j)*ndof+i] += oownc * N[i];
       K[(i)*ndof+oc+ndir+j] += oownc * N[i];
     }
   }
   delete[] u;
   delete[] e;
 }
};

void DEMHelm3d::interfMatrix(int fi, DEMElement* deme2,
                              complex<double> *K) {

 IsoParamUtils ipu(o);
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double omega = getOmega();

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 int ndofs = nPolynomialDofs()+ndir+
             deme2->nPolynomialDofs()+deme2->nEnrichmentDofs();

 for(int i=0;i<ndofs*ndofs;i++) K[i] = 0.0;

 HelmDEMWetMatrixFunction3d f(o,omega,ndir,cdir,xref,deme2,K);
 ipu.surfInt3d(xyz, fi, f);

 delete[] cdir;

}

class HelmDEMNeumVFunction3d : public IntegFunctionA3d {
 int order;
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
public:
 HelmDEMNeumVFunction3d(int _order,
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc,
                      complex<double> *_v) {
   order = _order; ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   int oc = order*order*order;
   
   complex<double> ince = w* (nsign* 
                  (incdir[0]*cross[0]+incdir[1]*cross[1]+incdir[2]*cross[2])) *
                  exp(incdir[0]*x[0] + incdir[1]*x[1] + incdir[2]*x[2]);

   for(int i=0;i<oc;i++) v[i] += ince*N[i];

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1])+
                                      dirs[i*3+2]*(x[2]-xc[2]));
   for(int i=0;i<ndir;i++)
     v[oc+i] += ince*e[i];
   delete[] e;
 }
};


class HelmDEMNeumFVFunction3d : public IntegFunctionA3d {
 int order;
 int ndir;
 complex<double> *dirs;
 complex<double> *f;
 double *xc;
 complex<double> *v;
public:
 HelmDEMNeumFVFunction3d(int _order,
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_f, double *_xc,
                      complex<double> *_v) {
   order = _order; ndir = _ndir; dirs = _dirs; f = _f; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   int oc = order*order*order;
  
   double wc = w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]); 
   complex<double> ff = 0.0;
   for(int i=0;i<oc;i++) ff = f[i]*N[i];

   for(int i=0;i<oc;i++) v[i] += wc*ff*N[i];

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1])+
                                      dirs[i*3+2]*(x[2]-xc[2]));
   for(int i=0;i<ndir;i++)
     v[oc+i] += wc*ff*e[i];
   delete[] e;
 }
};


void HelmDEMNeumV(int order, double *xy,
                     int ndir, complex<double> *dirs,
                     complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils ipu(order);
 int gorder = (order>0)?7:13;

 HelmDEMNeumVFunction3d f(order,ndir,dirs,incdir,xc,v);
 ipu.surfInt3d(xy, faceindex, f, gorder);
}


void HelmDEMNeumFV(int order, double *xy,
                     int ndir, complex<double> *dirs,
                     complex<double> *ff, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils ipu(order);
 int gorder = (order>0)?7:13;

 HelmDEMNeumFVFunction3d f(order,ndir,dirs,ff,xc,v);
 ipu.surfInt3d(xy, faceindex, f, gorder);
}



void DEMHelm3d::createRHS(complex<double>*v) {
 IsoParamUtils ipu(o);
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);


 complex<double> *vv = new complex<double>[nPolynomialDofs()+ndir];
 for(int i=0;i<nPolynomialDofs()+ndir;i++) vv[i] = 0.0;

 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2) {
     double *dr = getWaveDirection();
     complex<double> incdir[3] = {complex<double>(0.0,kappa*dr[0]),
                                  complex<double>(0.0,kappa*dr[1]),
                                  complex<double>(0.0,kappa*dr[2])};
     HelmDEMNeumV(o, xyz, ndir, cdir, incdir, i+1 , xref, vv);
   } else if (bc[i]==4) {
     complex<double> *f = getField();
     HelmDEMNeumFV(o, xyz, ndir, cdir, f, i+1, xref, vv);
   }
 }

 delete[] cdir;
 delete[] xyz;

 for(int i=0;i<nPolynomialDofs();i++)
   v[i] = vv[i]/rho;
 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[nPolynomialDofs()+i] = 0;
 for(int i=0;i<ndir;i++)
   v[nPolynomialDofs()+tnldir+i] = vv[nPolynomialDofs()+i]/rho;
 delete[] vv;
}
