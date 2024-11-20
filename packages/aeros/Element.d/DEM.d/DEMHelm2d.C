#include <cstdio>

#include <Element.d/DEM.d/DEMHelm2d.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>
#include <Element.d/Helm.d/PML.h>

class HelmDGMELMatrixFunction : public IntegFunctionL2d {
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
 HelmDGMELMatrixFunction(double _kappa,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs, double *_xsc, double *_xc,
                      complex<double> *__L, int _fi) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; nldir = _nldir; ldirs = _ldirs;
   L = __L; fi = _fi; xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
  
   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                  dirs[i*2+1]*(x[1]-xc[1]));

   for(int j=0;j<nldir;j++) {
     complex<double> l = exp(ldirs[j*2+0]*(x[0]-xsc[0])+
                             ldirs[j*2+1]*(x[1]-xsc[1]));

     for(int i=0;i<ndir;i++)
       L[j*ndir+i] += e[i]*l*w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
   }
   delete[] e;
 }
};


class HelmDGMSomEEMatrixFunction : public IntegFunctionL2d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *K;
public:
 HelmDGMSomEEMatrixFunction(double _kappa,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                  dirs[i*2+1]*(x[1]-xc[1]));

   complex<double> ikc = complex<double>(0.0,
                           w*kappa*sqrt(cross[0]*cross[0]+cross[1]*cross[1]));
   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] -= ikc*e[i]*e[j];
   delete[] e; 
 }
};


class HelmDGMEeMatrixFunction : public IntegFunctionL2d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *K;
public:
 HelmDGMEeMatrixFunction(double _kappa,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[2*ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                  dirs[i*2+1]*(x[1]-xc[1]));

   complex<double> *ikc = e+ndir;
   for(int i=0;i<ndir;i++) 
     ikc[i] = w*nsign*(cross[0]*dirs[i*2+0]+cross[1]*dirs[i*2+1]);

   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] += ikc[i]*e[i]*e[j];
   delete[] e;
 }
};


void HelmDGMEMatrices2d(int order, double *xy,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc,
                    complex<double> *K, complex<double> *L) {


 IsoParamUtils2d *ipu =
     (order>0)?new IsoParamUtils2d(order):new IsoParamUtils2dTri(-order);

 int nldir = nldirs[0] + nldirs[1] + nldirs[2];
 if (order>0) nldir += nldirs[3];
 for(int i=0;i<nldir*ndir;i++) L[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) K[i] = 0.0;

 int ic[5] = { 0, order-1, order*order-1, order*(order-1), 0 }; 
 int nf = 4;
 if (order<0) {
  ic[0] = -order-1;
  ic[1] = (-order*(-order+1))/2-1;
  ic[2] = 0;
  ic[3] = ic[0];
  nf = 3;
 }
 
 int ordersq = ipu->getordersq();
 double x[5] = {
   xy[0*ordersq+ic[0]],
   xy[0*ordersq+ic[1]],
   xy[0*ordersq+ic[2]],
   xy[0*ordersq+ic[3]],
   xy[0*ordersq+ic[4]]
 };
 double y[5] = {
   xy[1*ordersq+ic[0]],
   xy[1*ordersq+ic[1]],
   xy[1*ordersq+ic[2]],
   xy[1*ordersq+ic[3]],
   xy[1*ordersq+ic[4]]
 };

 complex<double> *e = new complex<double>[ndir*(nf+1)];
 for(int k=0;k<nf+1;k++)
   for(int i=0;i<ndir;i++)
     e[k*ndir+i] = exp(dirs[i*2+0]*(x[k]-xc[0]) +
                                   dirs[i*2+1]*(y[k]-xc[1]));

 int faceindex;
 int c = 0;
 complex<double> *dn = new complex<double>[ndir];
 for(faceindex=1;faceindex<=nf;faceindex++) {

   if (ipu->isStraight(xy,faceindex)) {
     double tau[2] = {
       x[faceindex]-x[faceindex-1],
       y[faceindex]-y[faceindex-1]
     };
     double len = ipu->sidelen(tau[0],tau[1]);
     tau[0] /= len;
     tau[1] /= len;
     double n[2] = { tau[1], -tau[0] };

     if (sflags[faceindex-1]==0)
       for(int i=0;i<ndir;i++) dn[i] = (dirs[i*2+0]*n[0] + dirs[i*2+1]*n[1]);
     else 
       for(int i=0;i<ndir;i++) dn[i] = (dirs[i*2+0]*n[0] + dirs[i*2+1]*n[1]) -
                                       complex<double>(0.0,kappa);

     for(int j=0;j<ndir;j++) {
       for(int i=0;i<ndir;i++) {
         complex<double> expalpha1 = e[ndir*(faceindex-1)+i]*
                                     e[ndir*(faceindex-1)+j];
         complex<double> expalpha2 = e[ndir*(faceindex)+i]*
                                     e[ndir*(faceindex)+j];
         complex<double> alpha = (dirs[i*2+0]+dirs[j*2+0])*tau[0] +
                                 (dirs[i*2+1]+dirs[j*2+1])*tau[1];
         K[j*ndir+i] += dn[i] * ipu->expdiff(alpha,expalpha1,expalpha2,len);
       }
     }

     if (nldirs[faceindex-1]!=0) {
       double *scxy = xsc + (faceindex-1)*2;
       complex<double> *ldrs = ldirs+c*2;
  
       for(int j=0;j<nldirs[faceindex-1];j++) {
  
         complex<double> l1 = exp(ldrs[j*2+0]*(x[faceindex-1]-scxy[0])+
                                  ldrs[j*2+1]*(y[faceindex-1]-scxy[1]));
         complex<double> l2 = exp(ldrs[j*2+0]*(x[faceindex]-scxy[0])+
                                  ldrs[j*2+1]*(y[faceindex]-scxy[1]));
  
         for(int i=0;i<ndir;i++) {
           complex<double> expalpha1 = e[ndir*(faceindex-1)+i]*l1;
           complex<double> expalpha2 = e[ndir*(faceindex)+i]*l2;
           complex<double> alpha = (dirs[i*2+0]+ldrs[j*2+0])*tau[0] +
                                   (dirs[i*2+1]+ldrs[j*2+1])*tau[1];
           L[(c+j)*ndir+i] += ipu->expdiff(alpha,expalpha1,expalpha2,len);
         }
       }
       c += nldirs[faceindex-1];
     }
     
   } else {
     HelmDGMEeMatrixFunction f(kappa,ndir,dirs,xc,K);
     ipu->lineInt2d(xy, faceindex, f);
     if (nldirs[faceindex-1]!=0) {
       HelmDGMELMatrixFunction f(kappa,ndir,dirs,
                  nldirs[faceindex-1],ldirs+c*2, xsc + (faceindex-1)*2, xc,
                  L+c*ndir,faceindex);
       ipu->lineInt2d(xy, faceindex, f);
       c += nldirs[faceindex-1];
     }
     if (sflags[faceindex-1]!=0) {
       HelmDGMSomEEMatrixFunction f(kappa,ndir,dirs,xc,K);
       ipu->lineInt2d(xy, faceindex, f);
     }
   }
 }
 delete[] e;
 delete[] ipu;
 delete[] dn;
}


class HelmDGMPMLEEMatrixFunction : public IntegFunctionA2d {
 double *pmldata;
 double *cxy;
 int ndir;
 complex<double> *dirs;
 complex<double> *K;
public:
 HelmDGMPMLEEMatrixFunction(double *_pmldata, 
                      int _ndir, complex<double> *_dirs,
                            double *_cxy, complex<double> *_K) {
   pmldata = _pmldata; 
   ndir = _ndir; dirs = _dirs;
   cxy = _cxy;
   K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {

   double kappa = pmldata[0];
   double a = pmldata[1];
   double b = pmldata[2];
   double gamma = pmldata[3];

   PMLFunction fp(a,b,gamma);
   PMLFunction fm(-a,-b,gamma);
   PMLFunction *fx = &fp;
   PMLFunction *fy = &fp;
   if (x[0]<0.0) fx = &fm;
   if (x[1]<0.0) fy = &fm;
   complex<double> wdetb = w*det*fx->beta(x[0])*fy->beta(x[1]);

   complex<double> xt[2];
   xt[0] = x[0]*fx->alpha(x[0]);
   xt[1] = x[1]*fy->alpha(x[1]);

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*xt[0] + dirs[i*2+1]*xt[1]);

   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] += wdetb*
                        (dirs[i*2+0]*dirs[j*2+0] + dirs[i*2+1]*dirs[j*2+1]-
                         complex<double>(kappa*kappa,0.0) ) * e[i] * e[j];
   delete[] e;
 }
};


class HelmDGMPMLELMatrixFunction : public IntegFunctionL2d {
 double *pmldata;
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
 HelmDGMPMLELMatrixFunction(double *_pmldata,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs,
                      double *_xsc, double *_xc,
                      complex<double> *__L, int _fi) {
   pmldata = _pmldata;
   ndir = _ndir; dirs = _dirs; nldir = _nldir; ldirs = _ldirs;
   L = __L; fi = _fi; xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   double kappa = pmldata[0];
   double a = pmldata[1];
   double b = pmldata[2];
   double gamma = pmldata[3];


   PMLFunction fp(a,b,gamma);
   PMLFunction fm(-a,-b,gamma);
   PMLFunction *fx = &fp;
   PMLFunction *fy = &fp;
   if (x[0]<0.0) fx = &fm;
   if (x[1]<0.0) fy = &fm;
   complex<double> wcb = w* sqrt(
                         cross[0]*cross[0]*fx->beta(x[1])*conj(fx->beta(x[1]))+
                         cross[1]*cross[1]*fx->beta(x[0])*conj(fx->beta(x[0])));
//RT: mess
/*   complex<double> wcb = w* sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
   double ll = sqrt(cross[0]*cross[0]+ cross[1]*cross[1]);
   double tau[2] = { nsign*cross[1]/ll,-nsign*cross[0]/ll };
   if (fabs(x[0])>=a || fabs(x[1])>=a) {
      if (fabs(tau[0])>fabs(tau[1]))
          wcb *= fx->beta(x[0]);
      else
          wcb *= fy->beta(x[1]);
//                   wcb *= tau[0]*fx->beta(x[0])+tau[1]*fy->beta(x[1]);
   }*/

   complex<double> xt[2];
   xt[0] = x[0]*fx->alpha(x[0]);
   xt[1] = x[1]*fy->alpha(x[1]);

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*xt[0] + dirs[i*2+1]*xt[1]);

   for(int j=0;j<nldir;j++)
     for(int i=0;i<ndir;i++) {
       complex<double> t = wcb*e[i]*
                           exp(ldirs[j*2+0]*(xt[0]) +
                               ldirs[j*2+1]*(xt[1]));

       L[j*ndir+i] += t;
     }
   delete[] e;
 }
};



void HelmDGMPMLEMatrices2d(int order, double *xy,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc, PMLProps *pml,
                    complex<double> *K, complex<double> *L) {


 IsoParamUtils2d *ipu =
     (order>0)?new IsoParamUtils2d(order):new IsoParamUtils2dTri(-order);

 int nldir = nldirs[0] + nldirs[1] + nldirs[2];
 if (order>0) nldir += nldirs[3];
 for(int i=0;i<nldir*ndir;i++) L[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) K[i] = 0.0;

 int nf = 4;
 if (order<0) {
  nf = 3;
 }
 
 int ordersq = ipu->getordersq();

 double pmldata[4] = { kappa, pml->Rx,pml->Sx,pml->gamma};
 HelmDGMPMLEEMatrixFunction f(pmldata,ndir,dirs,xc,K);
 ipu->areaInt2d(xy, f);

 int faceindex;
 int c = 0;

 for(faceindex=1;faceindex<=nf;faceindex++) {
   if (nldirs[faceindex-1]!=0) {
     HelmDGMPMLELMatrixFunction f(pmldata,ndir,dirs,
                nldirs[faceindex-1],ldirs+c*2, xsc + (faceindex-1)*2, xc,
                L+c*ndir,faceindex);
     ipu->lineInt2d(xy, faceindex, f);
     c += nldirs[faceindex-1];
   }
 }
 delete[] ipu;
}


DGMHelm2d_4::DGMHelm2d_4(int _nnodes, int* nodenums) :
  DGMHelm2d(_nnodes,nodenums) {
 ndir = 4;
}

DGMHelm2d_4t::DGMHelm2d_4t(int _nnodes, int* nodenums) :
  DGMHelm2d(-_nnodes,nodenums) {
 ndir = 4;
}

DGMHelm2d_8::DGMHelm2d_8(int _nnodes, int* nodenums) :
  DGMHelm2d(_nnodes,nodenums) {
 ndir = 8;
}

DGMHelm2d_8t::DGMHelm2d_8t(int _nnodes, int* nodenums) :
  DGMHelm2d(-_nnodes,nodenums) {
 ndir = 8;
}

DGMHelm2d_16::DGMHelm2d_16(int _nnodes, int* nodenums) :
  DGMHelm2d(_nnodes,nodenums) {
 ndir = 16;
}

DGMHelm2d_32::DGMHelm2d_32(int _nnodes, int* nodenums) :
  DGMHelm2d(_nnodes,nodenums) {
 ndir = 32;
}

DGMHelm2d_Eva::DGMHelm2d_Eva(int _nnodes, int* nodenums) :
 DGMHelm2d(_nnodes,nodenums) {
}

DGMHelm2d_Eva2_8::DGMHelm2d_Eva2_8(int _nnodes, int* nodenums) :
  DGMHelm2d_Eva(_nnodes,nodenums) {
 ndir = 10;
 gamma = 4.0;
 normal[0] = 0.0; normal[1] = 1.0;
}

complex<double> DGMHelm2d_4::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] =  {
   1,0,-1,0,0,1,0,-1};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DGMHelm2d_4t::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] =  {
   1,0,-1,0,0,1,0,-1};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DGMHelm2d_8::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 1,0,
                1/sqrt(2.0),1/sqrt(2.0),
                0,1,
               -1/sqrt(2.0),1/sqrt(2.0),
               -1,0,
               -1/sqrt(2.0),-1/sqrt(2.0),
                0,-1,
                1/sqrt(2.0),-1/sqrt(2.0)};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DGMHelm2d_8t::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 1,0,
                1/sqrt(2.0),1/sqrt(2.0),
                0,1,
               -1/sqrt(2.0),1/sqrt(2.0),
               -1,0,
               -1/sqrt(2.0),-1/sqrt(2.0),
                0,-1,
                1/sqrt(2.0),-1/sqrt(2.0)};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DGMHelm2d_16::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 
             cos(2*M_PI*double(0)/16), sin(2*M_PI*double(0)/16),
             cos(2*M_PI*double(1)/16), sin(2*M_PI*double(1)/16),
             cos(2*M_PI*double(2)/16), sin(2*M_PI*double(2)/16),
             cos(2*M_PI*double(3)/16), sin(2*M_PI*double(3)/16),
             cos(2*M_PI*double(4)/16), sin(2*M_PI*double(4)/16),
             cos(2*M_PI*double(5)/16), sin(2*M_PI*double(5)/16),
             cos(2*M_PI*double(6)/16), sin(2*M_PI*double(6)/16),
             cos(2*M_PI*double(7)/16), sin(2*M_PI*double(7)/16),
             cos(2*M_PI*double(8)/16), sin(2*M_PI*double(8)/16),
             cos(2*M_PI*double(9)/16), sin(2*M_PI*double(9)/16),
             cos(2*M_PI*double(10)/16), sin(2*M_PI*double(10)/16),
             cos(2*M_PI*double(11)/16), sin(2*M_PI*double(11)/16),
             cos(2*M_PI*double(12)/16), sin(2*M_PI*double(12)/16),
             cos(2*M_PI*double(13)/16), sin(2*M_PI*double(13)/16),
             cos(2*M_PI*double(14)/16), sin(2*M_PI*double(14)/16),
             cos(2*M_PI*double(15)/16), sin(2*M_PI*double(15)/16) };
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DGMHelm2d_32::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 
             cos(2*M_PI*double(0)/32), sin(2*M_PI*double(0)/32),
             cos(2*M_PI*double(1)/32), sin(2*M_PI*double(1)/32),
             cos(2*M_PI*double(2)/32), sin(2*M_PI*double(2)/32),
             cos(2*M_PI*double(3)/32), sin(2*M_PI*double(3)/32),
             cos(2*M_PI*double(4)/32), sin(2*M_PI*double(4)/32),
             cos(2*M_PI*double(5)/32), sin(2*M_PI*double(5)/32),
             cos(2*M_PI*double(6)/32), sin(2*M_PI*double(6)/32),
             cos(2*M_PI*double(7)/32), sin(2*M_PI*double(7)/32),
             cos(2*M_PI*double(8)/32), sin(2*M_PI*double(8)/32),
             cos(2*M_PI*double(9)/32), sin(2*M_PI*double(9)/32),
             cos(2*M_PI*double(10)/32), sin(2*M_PI*double(10)/32),
             cos(2*M_PI*double(11)/32), sin(2*M_PI*double(11)/32),
             cos(2*M_PI*double(12)/32), sin(2*M_PI*double(12)/32),
             cos(2*M_PI*double(13)/32), sin(2*M_PI*double(13)/32),
             cos(2*M_PI*double(14)/32), sin(2*M_PI*double(14)/32),
             cos(2*M_PI*double(15)/32), sin(2*M_PI*double(15)/32),
             cos(2*M_PI*double(16)/32), sin(2*M_PI*double(16)/32),
             cos(2*M_PI*double(17)/32), sin(2*M_PI*double(17)/32),
             cos(2*M_PI*double(18)/32), sin(2*M_PI*double(18)/32),
             cos(2*M_PI*double(19)/32), sin(2*M_PI*double(19)/32),
             cos(2*M_PI*double(20)/32), sin(2*M_PI*double(20)/32),
             cos(2*M_PI*double(21)/32), sin(2*M_PI*double(21)/32),
             cos(2*M_PI*double(22)/32), sin(2*M_PI*double(22)/32),
             cos(2*M_PI*double(23)/32), sin(2*M_PI*double(23)/32),
             cos(2*M_PI*double(24)/32), sin(2*M_PI*double(24)/32),
             cos(2*M_PI*double(25)/32), sin(2*M_PI*double(25)/32),
             cos(2*M_PI*double(26)/32), sin(2*M_PI*double(26)/32),
             cos(2*M_PI*double(27)/32), sin(2*M_PI*double(27)/32),
             cos(2*M_PI*double(28)/32), sin(2*M_PI*double(28)/32),
             cos(2*M_PI*double(29)/32), sin(2*M_PI*double(29)/32),
             cos(2*M_PI*double(30)/32), sin(2*M_PI*double(30)/32),
             cos(2*M_PI*double(31)/32), sin(2*M_PI*double(31)/32) };
 return complex<double>(0.0,kappa*a[n]);
}


complex<double> DGMHelm2d_Eva2_8::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 1,0,
                1/sqrt(2.0),1/sqrt(2.0),
                0,1,
               -1/sqrt(2.0),1/sqrt(2.0),
               -1,0,
               -1/sqrt(2.0),-1/sqrt(2.0),
                0,-1,
                1/sqrt(2.0),-1/sqrt(2.0)};
 if (n<16) return complex<double>(0.0,kappa*a[n]);
 complex<double> b[] = { 
   complex<double>( gamma*kappa*normal[0],
                   sqrt(1+gamma*gamma)*kappa*(-normal[1]) ),
   complex<double>( gamma*kappa*normal[1],
                   sqrt(1+gamma*gamma)*kappa*(normal[0]) ),
   complex<double>( -gamma*kappa*normal[0],
                   sqrt(1+gamma*gamma)*kappa*(-normal[1]) ),
   complex<double>( -gamma*kappa*normal[1],
                   sqrt(1+gamma*gamma)*kappa*(normal[0]) ) };
   return b[n-16];
}


complex<double> DGMHelm2d_1_LM::coef(int n) {
 return 0.0;
}

complex<double> DGMHelm2d_2_LM::coef(int n) {
 if (n==0) return 0.5;
 else return -0.5;
}

complex<double> DGMHelm2d_4_LM::coef(int n) {
 if (n==0) return 0.2;
 else if (n==1) return -0.2;
 else if (n==2) return 0.75;
 else return -0.75;
}

complex<double> DGMHelm2d_8_LM::coef(int n) {
 double a[] = {
   0.1, -0.1, 0.707, -0.707,
   cos(M_PI/8), -cos(M_PI/8), sin(M_PI/8), -sin(M_PI/8) };
 return a[n];
}

void DGMHelm2d_Eva1_2_LM::init() {

 DGMHelm2d *eh1 = dynamic_cast<DGMHelm2d*>(e1);
 int o = eh1->o;
 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 e1->getNodalCoord(os,e1->nn,xyz);

 int fi=-1;
 for(int i=0;i<e1->nFaces();i++) 
   if (e1->lm[i]==this) { fi = i+1; break; }

 int nc[2];
 if (fi==1) { nc[0] = 0; nc[1] = o-1; }
 else if (fi==2) { nc[0] = o-1; nc[1] = o*o-1; }
 else if (fi==3) { nc[0] = o*o-1; nc[1] = o*(o-1); }
 else if (fi==4) { nc[0] = o*(o-1); nc[1] = 0; }

 double tau[2] = { xyz[nc[1]] - xyz[nc[0]], xyz[os+nc[1]] - xyz[os+nc[0]] };
 double l = sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
 tau[0] /= l;
 tau[1] /= l;

 double lkappa = std::max(e1->getProperty()->kappaHelm,
                          e2->getProperty()->kappaHelm);
 

 DGMHelm2d_Eva *ee1 = dynamic_cast<DGMHelm2d_Eva*>(e1);
 DGMHelm2d_Eva *ee2 = dynamic_cast<DGMHelm2d_Eva*>(e2);
 // Average normals and gamma
 gamma = 0.0;
 normal[0] = normal[1] = 0;
 int cn = 0;
 if (ee1) {
  normal[0] += ee1->normal[0];
  normal[1] += ee1->normal[1];
  gamma += ee1->gamma;
  cn++;
 }
 if (ee2) {
  normal[0] += ee2->normal[0];
  normal[1] += ee2->normal[1];
  gamma += ee2->gamma;
  cn++;
 }
 if (cn==0) {
    fprintf(stderr,
     "Undefined normal in DGMHelm2d_Eva1_2_LM:init. Exiting.\n");
    exit(-1);
 }
 l = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
 normal[0] /= l;
 normal[1] /= l;
 gamma /= cn;

 double s = std::abs(normal[0]*tau[0]+normal[1]*tau[1]);
 if (s>0.8) ndofs = 3;
 else ndofs = 4;

 delete[] xyz;
}


complex<double> DGMHelm2d_Eva1_2_LM::ldir(int n, double tau[2]) {
 double lkappa = std::max(e1->getProperty()->kappaHelm,
                          e2->getProperty()->kappaHelm);
 if (n==0) return complex<double>(0.0, 0.5*lkappa*tau[0]);
 else if (n==1) return complex<double>(0.0, 0.5*lkappa*tau[1]);
 else if (n==2) return complex<double>(0.0, -0.5*lkappa*tau[0]);
 else if (n==3) return complex<double>(0.0, -0.5*lkappa*tau[1]);
 else {
   if (ndofs==3) {
     if (n==4) return complex<double>(gamma*lkappa*normal[0],0.0);
     else      return complex<double>(gamma*lkappa*normal[1],0.0);
   } else {
     double g = sqrt(1+gamma*gamma);
     if (n==4) return complex<double>(0.0, g*lkappa*tau[0]);
     else if (n==5) return complex<double>(0.0, g*lkappa*tau[1]);
     else if (n==6) return complex<double>(0.0, -g*lkappa*tau[0]);
     else return complex<double>(0.0, -g*lkappa*tau[1]);
   }
 }
}



DGMHelm2d::DGMHelm2d(int _nnodes, int* nodenums) {
 if (_nnodes<0) {
   _nnodes=-_nnodes;
   o = - int(-1.0+sqrt(1.0+8.0*double(_nnodes)))/2;
 }
 else {
   o = int(sqrt(double(_nnodes))+0.5);
 }
 nn = new int[_nnodes]; 
 for(int i=0;i<_nnodes;i++) nn[i] = nodenums[i];
 lm = new DEMLM*[nFaces()];
 for(int i=0;i<nFaces();i++) lm[i] = 0;
 bc = new int[nFaces()];
 for(int i=0;i<nFaces();i++) bc[i] = 0;
 ndir = 0;
}


complex<double> DGMHelm2d::dir(int n) {
 return 0.0;
}

void DGMHelm2d::getRef(double *xyz,double *xy) {
// IsoParamUtils2d ipu(o);
// ipu.elementcenter(xyz,xy);
 xy[0] = xy[1] = 0.0;
}

void DGMHelm2d::createM(complex<double>*M) {

 IsoParamUtils2d *ipu =
     (o>0)?new IsoParamUtils2d(o):new IsoParamUtils2dTri(-o);
 int os = ipu->getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[2];
 getRef(xyz,xref);


 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 int nf = (o<0)?3:4;
 int nldir[4] = {0,0,0,0};
 for(int fi=0;fi<nf;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<nf;i++) tnldir += nldir[i];

 for(int fi=0;fi<nf;fi++) if (lm[fi]!=0) {
   DGMHelm2d_LM *l = dynamic_cast<DGMHelm2d_LM*>(lm[fi]);
   DGMHelm2d_Eva_LM *le = dynamic_cast<DGMHelm2d_Eva_LM*>(lm[fi]);
   if (l==0 && le==0) {
     fprintf(stderr,"DGMHelm2d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMHelm2d_LM type specific
 double xlref[8] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<4;fi++) 
//   if (lm[fi]!=0) ipu.sidecenter(xyz,fi,xlref+fi*2);

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];

 complex<double>* cldir = new complex<double>[2*tnldir];
 int corner[5] = { 1-1, o-1, os-1, os-o+1-1, 1-1};
 if (o<0) 
  { corner[0] = -o-1; corner[1] = os-1; corner[2] = 0; corner[3] = corner[0]; }
   
 int cc = 0;
 double sign[4] = { 1,1,1,1};
 for(int i=0;i<nf;i++) {
   double tau[2]= { xyz[corner[i+1]]-xyz[corner[i]],
                    xyz[os+corner[i+1]] - xyz[os+corner[i]] };
   if (nn[corner[i]]>nn[corner[i+1]]) {
     tau[0] = -tau[0];
     tau[1] = -tau[1];
     sign[i] = -1.0;
   }
   double l = sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
   tau[0] /= l;
   tau[1] /= l;
   if (nldir[i]>0) {
     DGMHelm2d_LM *l = dynamic_cast<DGMHelm2d_LM*>(lm[i]);
     DGMHelm2d_Eva_LM *le = dynamic_cast<DGMHelm2d_Eva_LM*>(lm[i]);
     if (le!=0) {
       for(int j=0;j<lm[i]->nDofs();j++) {
         cldir[2*cc+0] = le->ldir(2*j,tau);
         cldir[2*cc+1] = le->ldir(2*j+1,tau);
         cc++;
       }
     }
     else {
       // Max kappa
       double lkappa = std::max(l->e1->getProperty()->kappaHelm,
                                l->e2->getProperty()->kappaHelm);
       for(int j=0;j<lm[i]->nDofs();j++) {
         cldir[2*cc+0] = complex<double>(0.0, lkappa)*l->coef(j)*tau[0];
         cldir[2*cc+1] = complex<double>(0.0, lkappa)*l->coef(j)*tau[1];
         cc++;
       }
     }
   }
 }

 int arbFlag[4] = { 0,0,0,0};
 for(int i=0;i<nf;i++) if (bc[i]==1) arbFlag[i] = 1;

 PMLProps* pml = getPMLProps();
 if (pml->PMLtype==0)
   HelmDGMEMatrices2d(o, xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, kee, kel); 
 else 
   HelmDGMPMLEMatrices2d(o, xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, pml, kee, kel); 

 cc = 0;
 for(int i=0;i<nf;i++) {
  for(int j=cc;j<cc+nldir[i];j++) for(int k=0;k<ndir;k++)
    kel[j*ndir+k] *= sign[i];
  cc += nldir[i];
 }
// for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++)
//fprintf(stderr,"zzz %d %d %e %e\n",i+1,j+1,real(kel[j*ndir+i]),imag(kel[j*ndir+i]));

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
 delete[] ipu;
}


class HelmDGMWetEMatrixFunction : public IntegFunctionL2d {
 double omega;
 int ndir;
 complex<double> *dirs;
 double *xc;
 DEMElement *deme2;
 complex<double> *K;
public:
 HelmDGMWetEMatrixFunction(double _omega,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      DEMElement *_deme2,
                      complex<double> *_K) {
   omega = _omega; ndir = _ndir; dirs = _dirs; xc = _xc; deme2 = _deme2;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
/*
   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) { e[i] = omega*omega*w*
           exp(dirs[i*2+0]*(x[0]-xc[0]) + dirs[i*2+1]*(x[1]-xc[1]) );
   }

   int ndir2 = deme2->nEnrichmentDofs();
   complex<double> *e2 = new complex<double>[ndir2*2];
   deme2->enrichmentF(x,e2);
   for(int j=0;j<ndir2;j++) {
// Get deme2 function
     for(int i=0;i<ndir;i++) {
       complex<double> v = nsign*e[i]*(cross[0]*e2[2*j+0]+cross[1]*e2[2*j+1]);
       K[(j+ndir)*(ndir+ndir2)+i] += v;
       K[i*(ndir+ndir2)+j+ndir] += v;
     }
   }
   delete[] e2;
   delete[] e;
*/

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) { e[i] =
           exp(dirs[i*2+0]*(x[0]-xc[0]) + dirs[i*2+1]*(x[1]-xc[1]) );
   }

   int os2 = deme2->nPolynomialDofs();
   int ndir2 = deme2->nEnrichmentDofs();
   int ndof = ndir + os2 + ndir2;
   double oown = omega*omega*w*nsign;

   complex<double> *e2 = new complex<double>[ndir2*2];
   deme2->enrichmentF(x,e2); 

   for(int j=0;j<ndir2;j++) {
     complex<double> oownc = oown*(cross[0]*e2[j*2+0]+cross[1]*e2[j*2+1]);
     for(int i=0;i<ndir;i++) {
       K[(ndir+os2+j)*ndof+i] += oownc * e[i];
       K[(i)*ndof+ndir+os2+j] += oownc * e[i];
     }
   }
   delete[] e2;

   if (os2>0)  {
     double *u = new double[os2*2];
     deme2->polynomialF(x,u);
  
     for(int j=0;j<os2;j++) {
       double oownc = oown*(cross[0]*u[2*j+0]+cross[1]*u[2*j+1]);
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



void DGMHelm2d::interfMatrix(int fi, DEMElement* deme2,
                              complex<double> *K) {

 IsoParamUtils2d *ipu =
     (o>0)?new IsoParamUtils2d(o):new IsoParamUtils2dTri(-o);
 int os = ipu->getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double omega = getOmega();

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 int ndofs = nPolynomialDofs()+ndir+
             deme2->nPolynomialDofs()+deme2->nEnrichmentDofs();
 for(int i=0;i<ndofs*ndofs;i++) K[i] = 0.0;

 HelmDGMWetEMatrixFunction f(omega,ndir,cdir,xref,deme2,K);
 ipu->lineInt2d(xyz, fi, f);

 delete[] cdir;
}


class HelmDGMNeumVFunction : public IntegFunctionL2d {
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
public:
 HelmDGMNeumVFunction(
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc,
                      complex<double> *_v) {
   ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int i;
   complex<double> *e = new complex<double>[ndir];
   for(i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                  dirs[i*2+1]*(x[1]-xc[1]));


   complex<double> ince = w* (nsign*
        (incdir[0]*cross[0]+incdir[1]*cross[1])
         )
         * exp(incdir[0]*x[0] + incdir[1]*x[1]);

   for(i=0;i<ndir;i++)
     v[i] += ince*e[i];
   delete[] e; 
 }
};


void HelmDGMNeumV2d(int order, double *xy,
                     int ndir, complex<double> *dirs,
                     complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils2d *ipu =
     (order>0)?new IsoParamUtils2d(order):new IsoParamUtils2dTri(-order);

 HelmDGMNeumVFunction f(ndir,dirs,incdir,xc,v);
 ipu->lineInt2d(xy, faceindex, f);
 delete ipu;
}



void DGMHelm2d::createRHS(complex<double>*v) {
 int os;
 if (o>0) { 
   IsoParamUtils2d ipu(o);
   os = ipu.getordersq();
 } else {
   IsoParamUtils2dTri ipu(-o);
   os = ipu.getordersq();
 }
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 complex<double> *vv = new complex<double>[ndir];
 for(int i=0;i<ndir;i++) vv[i] = 0.0;

 double *dir = getWaveDirection();
 complex<double> incdir[2] = {complex<double>(0.0,kappa*dir[0]),
                              complex<double>(0.0,kappa*dir[1])};
 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2) {
     HelmDGMNeumV2d(o, xyz, ndir, cdir, incdir, i+1 , xref, vv);
   }
 }
/*
 complex<double> *vv2 = new complex<double>[ndir];
 for(int i=0;i<ndir;i++) vv2[i] = 0.0;
 complex<double> incdir2[2] = {complex<double>(0.0,-kappa*1.0),
                              complex<double>(0.0,kappa*0.0)};
 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     HelmDGMNeumV2d(o, xyz, ndir, cdir, incdir2, i+1 , xref, vv2);
 }
 for(int i=0;i<ndir;i++) vv[i] = vv[i]+
                                 0.93680538552900*vv2[i];
 delete[] vv2;
*/

 delete[] cdir;
 delete[] xyz;

 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[i] = 0;
 for(int i=0;i<ndir;i++) 
   v[tnldir+i] = vv[i]/rho;
 delete[] vv;
}


void DGMHelm2d::createSol(double *xyz,
                        complex<double>* sol, complex<double> *nodalSol) {

 for(int i=0;i<8;i++) nodalSol[i] = 0.0;
 for(int i=0;i<ndir;i++) nodalSol[0] += 
   sol[i]*exp( dir(2*i)*xyz[0]+dir(2*i+1)*xyz[1] );
}

DEMHelm2d::DEMHelm2d(int _nnodes, int* nodenums) : DGMHelm2d(_nnodes, nodenums)  { }


DEMHelm2d_4::DEMHelm2d_4(int _o, int* nodenums) : DEMHelm2d(_o,nodenums) { ndir = 4; }
DEMHelm2d_4t::DEMHelm2d_4t(int _o, int* nodenums) : DEMHelm2d(-_o,nodenums) { ndir = 4; }
DEMHelm2d_8::DEMHelm2d_8(int _o, int* nodenums) : DEMHelm2d(_o,nodenums) { ndir = 8; }
DEMHelm2d_8t::DEMHelm2d_8t(int _o, int* nodenums) : DEMHelm2d(-_o,nodenums) { ndir = 8; }
DEMHelm2d_16::DEMHelm2d_16(int _o, int* nodenums) : DEMHelm2d(_o,nodenums) { ndir = 16; }
DEMHelm2d_32::DEMHelm2d_32(int _o, int* nodenums) : DEMHelm2d(_o,nodenums) { ndir = 32; }

complex<double> DEMHelm2d_4::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] =  {
   1,0,-1,0,0,1,0,-1};
 return complex<double>(0.0,kappa*a[n]);
}


complex<double> DEMHelm2d_4t::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] =  {
   1,0,-1,0,0,1,0,-1};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DEMHelm2d_8::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 1,0,
                1/sqrt(2.0),1/sqrt(2.0),
                0,1,
               -1/sqrt(2.0),1/sqrt(2.0),
               -1,0,
               -1/sqrt(2.0),-1/sqrt(2.0),
                0,-1,
                1/sqrt(2.0),-1/sqrt(2.0)};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DEMHelm2d_8t::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 1,0,
                1/sqrt(2.0),1/sqrt(2.0),
                0,1,
               -1/sqrt(2.0),1/sqrt(2.0),
               -1,0,
               -1/sqrt(2.0),-1/sqrt(2.0),
                0,-1,
                1/sqrt(2.0),-1/sqrt(2.0)};
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DEMHelm2d_16::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 
             cos(2*M_PI*double(0)/16), sin(2*M_PI*double(0)/16),
             cos(2*M_PI*double(1)/16), sin(2*M_PI*double(1)/16),
             cos(2*M_PI*double(2)/16), sin(2*M_PI*double(2)/16),
             cos(2*M_PI*double(3)/16), sin(2*M_PI*double(3)/16),
             cos(2*M_PI*double(4)/16), sin(2*M_PI*double(4)/16),
             cos(2*M_PI*double(5)/16), sin(2*M_PI*double(5)/16),
             cos(2*M_PI*double(6)/16), sin(2*M_PI*double(6)/16),
             cos(2*M_PI*double(7)/16), sin(2*M_PI*double(7)/16),
             cos(2*M_PI*double(8)/16), sin(2*M_PI*double(8)/16),
             cos(2*M_PI*double(9)/16), sin(2*M_PI*double(9)/16),
             cos(2*M_PI*double(10)/16), sin(2*M_PI*double(10)/16),
             cos(2*M_PI*double(11)/16), sin(2*M_PI*double(11)/16),
             cos(2*M_PI*double(12)/16), sin(2*M_PI*double(12)/16),
             cos(2*M_PI*double(13)/16), sin(2*M_PI*double(13)/16),
             cos(2*M_PI*double(14)/16), sin(2*M_PI*double(14)/16),
             cos(2*M_PI*double(15)/16), sin(2*M_PI*double(15)/16) };
 return complex<double>(0.0,kappa*a[n]);
}

complex<double> DEMHelm2d_32::dir(int n) {
 double kappa = getOmega()/getSpeedOfSound();
 double a[] = { 
             cos(2*M_PI*double(0)/32), sin(2*M_PI*double(0)/32),
             cos(2*M_PI*double(1)/32), sin(2*M_PI*double(1)/32),
             cos(2*M_PI*double(2)/32), sin(2*M_PI*double(2)/32),
             cos(2*M_PI*double(3)/32), sin(2*M_PI*double(3)/32),
             cos(2*M_PI*double(4)/32), sin(2*M_PI*double(4)/32),
             cos(2*M_PI*double(5)/32), sin(2*M_PI*double(5)/32),
             cos(2*M_PI*double(6)/32), sin(2*M_PI*double(6)/32),
             cos(2*M_PI*double(7)/32), sin(2*M_PI*double(7)/32),
             cos(2*M_PI*double(8)/32), sin(2*M_PI*double(8)/32),
             cos(2*M_PI*double(9)/32), sin(2*M_PI*double(9)/32),
             cos(2*M_PI*double(10)/32), sin(2*M_PI*double(10)/32),
             cos(2*M_PI*double(11)/32), sin(2*M_PI*double(11)/32),
             cos(2*M_PI*double(12)/32), sin(2*M_PI*double(12)/32),
             cos(2*M_PI*double(13)/32), sin(2*M_PI*double(13)/32),
             cos(2*M_PI*double(14)/32), sin(2*M_PI*double(14)/32),
             cos(2*M_PI*double(15)/32), sin(2*M_PI*double(15)/32),
             cos(2*M_PI*double(16)/32), sin(2*M_PI*double(16)/32),
             cos(2*M_PI*double(17)/32), sin(2*M_PI*double(17)/32),
             cos(2*M_PI*double(18)/32), sin(2*M_PI*double(18)/32),
             cos(2*M_PI*double(19)/32), sin(2*M_PI*double(19)/32),
             cos(2*M_PI*double(20)/32), sin(2*M_PI*double(20)/32),
             cos(2*M_PI*double(21)/32), sin(2*M_PI*double(21)/32),
             cos(2*M_PI*double(22)/32), sin(2*M_PI*double(22)/32),
             cos(2*M_PI*double(23)/32), sin(2*M_PI*double(23)/32),
             cos(2*M_PI*double(24)/32), sin(2*M_PI*double(24)/32),
             cos(2*M_PI*double(25)/32), sin(2*M_PI*double(25)/32),
             cos(2*M_PI*double(26)/32), sin(2*M_PI*double(26)/32),
             cos(2*M_PI*double(27)/32), sin(2*M_PI*double(27)/32),
             cos(2*M_PI*double(28)/32), sin(2*M_PI*double(28)/32),
             cos(2*M_PI*double(29)/32), sin(2*M_PI*double(29)/32),
             cos(2*M_PI*double(30)/32), sin(2*M_PI*double(30)/32),
             cos(2*M_PI*double(31)/32), sin(2*M_PI*double(31)/32) };
 return complex<double>(0.0,kappa*a[n]);
}


class HelmPLMatrixFunction : public IntegFunctionL2d {
 int o;
 double kappa;
 int nldir;
 int *nldirs;
 complex<double> *ldirs;
 complex<double> *L;
 int fi;
 double *xsc;
public:
 HelmPLMatrixFunction(int _o, double _kappa,
                      int _nldir, complex<double> *_ldirs, double *_xsc,
                      complex<double> *__L, int _fi) {
   o= _o; kappa = _kappa; nldir = _nldir; ldirs = _ldirs;
   L = __L; fi = _fi; xsc = _xsc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
  
   for(int j=0;j<nldir;j++) {
     complex<double> l = exp(ldirs[j*2+0]*(x[0]-xsc[0])+
                             ldirs[j*2+1]*(x[1]-xsc[1]));

     for(int i=0;i<o*o;i++)
       L[j*o*o+i] += N[i]*l*w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
   }
 }
};


class HelmSomPPEMatrixFunction : public IntegFunctionL2d {
 int o;
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *PP;
 complex<double> *PE;
public:
 HelmSomPPEMatrixFunction(int _o, double _kappa,
                          int _ndir, complex<double> *_dirs, double *_xc,
                          complex<double> *_PP, complex<double> *_PE) {
   o = _o; kappa = _kappa; 
   ndir = _ndir; dirs = _dirs; xc = _xc;
   PP = _PP; PE = _PE;

 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                      dirs[i*2+1]*(x[1]-xc[1]) );

   complex<double> ikc = complex<double>(0.0,
                           w*kappa*sqrt(cross[0]*cross[0]+cross[1]*cross[1]));
   for(int j=0;j<o*o;j++) {
     for(int i=0;i<o*o;i++)
       PP[j*o*o+i] -= ikc*N[i]*N[j];
     for(int i=0;i<ndir;i++)
       PE[i*o*o+j] -= ikc*N[j]*e[i];
   }
   delete[] e;

 }
};


class HelmPPMatrixFunction : public IntegFunctionA2d {
 int o;
 double kappa;
 complex<double> *K;
public:
 HelmPPMatrixFunction(int _o, double _kappa,
                      complex<double> *_K) {
   o = _o; kappa = _kappa; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {
   for(int j=0;j<o*o;j++)
     for(int i=0;i<o*o;i++)
       K[j*o*o+i] += w*det*(-kappa*kappa*N[i]*N[j]
        +dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]);
 }
};


class HelmPEMatrixFunction : public IntegFunctionA2d {
 int o;
 int ndir;
 complex<double> *dirs;
 double *xc;
 double kappa;
 complex<double> *K;
public:
 HelmPEMatrixFunction(int _o,  double _kappa, int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   o = _o; kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                  dirs[i*2+1]*(x[1]-xc[1]));
   for(int j=0;j<ndir;j++)
     for(int i=0;i<o*o;i++)
       K[j*o*o+i] += w*det*e[j]*(-kappa*kappa*N[i]
        +dNdx[i][0]*dirs[j*2+0]+dNdx[i][1]*dirs[j*2+1]);
   delete[] e;
 }
};


void HelmPMatrices2d(int order, double *xy,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc, double *xc,
                    complex<double> *K, complex<double> *L, complex<double> *KE) {


 IsoParamUtils2d ipu(order);
 int os = order*order;

 int nldir = nldirs[0] + nldirs[1] + nldirs[2] + nldirs[3];
 for(int i=0;i<nldir*os;i++) L[i] = 0.0;
 for(int i=0;i<os*os;i++) K[i] = 0.0;
 for(int i=0;i<ndir*os;i++) KE[i] = 0.0;

 HelmPPMatrixFunction f(order,kappa,K);
 ipu.areaInt2d(xy, f);
 HelmPEMatrixFunction fa(order,kappa,ndir,dirs,xc,KE);
 ipu.areaInt2d(xy, fa);
 int faceindex;
 int c = 0;
 for(faceindex=1;faceindex<=4;faceindex++) {

     if (nldirs[faceindex-1]!=0) {
       HelmPLMatrixFunction f(order,kappa,
                  nldirs[faceindex-1],ldirs+c*2, xsc + (faceindex-1)*2,
                  L+c*os,faceindex);
       ipu.lineInt2d(xy, faceindex, f);
       c += nldirs[faceindex-1];
     }
     if (sflags[faceindex-1]!=0) {
       HelmSomPPEMatrixFunction f(order,kappa,ndir,dirs,xc,K,KE);
       ipu.lineInt2d(xy, faceindex, f);
     }
 }
}

void DEMHelm2d::createM(complex<double>*M) {

 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 int nldir[4] = {0,0,0,0};
 for(int fi=0;fi<4;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<4;i++) tnldir += nldir[i];

 for(int fi=0;fi<4;fi++) if (lm[fi]!=0) {
   DGMHelm2d_LM *l = dynamic_cast<DGMHelm2d_LM*>(lm[fi]);
   DGMHelm2d_Eva_LM *le = dynamic_cast<DGMHelm2d_Eva_LM*>(lm[fi]);
   if (l==0 && le==0) {
     fprintf(stderr,"DGMHelm2d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMHelm2d_LM type specific
 double xlref[8] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<4;fi++) 
//   if (lm[fi]!=0) ipu.sidecenter(xyz,fi,xlref+fi*2);


 complex<double>* cldir = new complex<double>[2*tnldir];
 int corner[5] = { 1-1, o-1, os-1, os-o+1-1, 1-1};
 int cc = 0;
 double sign[4] = { 1,1,1,1};
 for(int i=0;i<4;i++) {
   double tau[2]= { xyz[corner[i+1]]-xyz[corner[i]],
                    xyz[os+corner[i+1]] - xyz[os+corner[i]] };
   if (nn[corner[i]]>nn[corner[i+1]]) {
     tau[0] = -tau[0];
     tau[1] = -tau[1];
     sign[i] = -1.0;
   }
   double l = sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
   tau[0] /= l;
   tau[1] /= l;
   if (nldir[i]>0) {
     DGMHelm2d_LM *l = dynamic_cast<DGMHelm2d_LM*>(lm[i]);
     DGMHelm2d_Eva_LM *le = dynamic_cast<DGMHelm2d_Eva_LM*>(lm[i]);
     if (le!=0) {
       for(int j=0;j<lm[i]->nDofs();j++) {
         cldir[2*cc+0] = le->ldir(2*j,tau);
         cldir[2*cc+1] = le->ldir(2*j+1,tau);
         cc++;
       }
     }
     else {
       // Max kappa
       double lkappa = std::max(l->e1->getProperty()->kappaHelm,
                                l->e2->getProperty()->kappaHelm);
       for(int j=0;j<lm[i]->nDofs();j++) {
         cldir[2*cc+0] = complex<double>(0.0, lkappa)*l->coef(j)*tau[0];
         cldir[2*cc+1] = complex<double>(0.0, lkappa)*l->coef(j)*tau[1];
         cc++;
       }
     }
   }
 }

 int arbFlag[4] = { 0,0,0,0};
 for(int i=0;i<4;i++) if (bc[i]==1) arbFlag[i] = 1;

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];
 HelmDGMEMatrices2d(o, xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, kee, kel); 
 complex<double>* kpp = new complex<double>[o*o*o*o];
 complex<double>* kpe = new complex<double>[ndir*o*o];
 complex<double>* kpl = new complex<double>[o*o*tnldir];
 HelmPMatrices2d(o, xyz, ndir, cdir, nldir, cldir, kappa, arbFlag, xlref, xref,
                    kpp, kpl, kpe); 

 cc = 0;
 for(int i=0;i<4;i++) {
  for(int j=cc;j<cc+nldir[i];j++) {
    for(int k=0;k<ndir;k++) kel[j*ndir+k] *= sign[i];
    for(int k=0;k<o*o;k++) kpl[j*o*o+k] *= sign[i];
  }
  cc += nldir[i];
 }

 for(int i=0;i<o*o;i++) for(int j=0;j<o*o;j++) M[j*(o*o+tnldir+ndir)+i] = kpp[j*o*o+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(o*o+tnldir+j)*(o*o+ndir+tnldir)+(o*o+tnldir+i)] = kee[j*ndir+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(o*o+tnldir+i)*(o*o+ndir+tnldir)+o*o+j] = kel[j*ndir+i]/rho;
   M[(o*o+j)*(o*o+ndir+tnldir)+(o*o+tnldir+i)] = kel[j*ndir+i]/rho;
 }
 for(int i=0;i<o*o;i++) for(int j=0;j<ndir;j++) {
   M[i*(o*o+ndir+tnldir)+(o*o+tnldir+j)] = kpe[j*o*o+i]/rho;
   M[(o*o+tnldir+j)*(o*o+ndir+tnldir)+i] = kpe[j*o*o+i]/rho;
 }
 for(int i=0;i<o*o;i++) for(int j=0;j<tnldir;j++) {
   M[i*(o*o+ndir+tnldir)+(o*o+j)] = kpl[j*o*o+i]/rho;
   M[(o*o+j)*(o*o+ndir+tnldir)+i] = kpl[j*o*o+i]/rho;
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++) 
   M[(o*o+j)*(o*o+ndir+tnldir)+o*o+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] kpp;
 delete[] kpe;
 delete[] kpl;

 delete[] xyz;
 delete[] cldir;
 delete[] cdir;
}


class HelmDEMWetMatrixFunction : public IntegFunctionL2d {
 int o;
 double omega;
 int ndir;
 complex<double> *dirs;
 double *xc;
 DEMElement *deme2;
 complex<double> *K;
public:
 HelmDEMWetMatrixFunction(int _o, double _omega,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      DEMElement *_deme2,
                      complex<double> *_K) {
   o= _o; omega = _omega; ndir = _ndir; dirs = _dirs; xc = _xc; deme2 = _deme2;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) { e[i] =
           exp(dirs[i*2+0]*(x[0]-xc[0]) + dirs[i*2+1]*(x[1]-xc[1]) );
   }

   int os = o*o;
   int os2 = deme2->nPolynomialDofs();
   int ndir2 = deme2->nEnrichmentDofs();
   int ndof = os + ndir + os2 + ndir2;
   double oown = omega*omega*w*nsign;

   complex<double> *e2 = new complex<double>[ndir2*2];
   deme2->enrichmentF(x,e2); 

   for(int j=0;j<ndir2;j++) {
     complex<double> oownc = oown*(cross[0]*e2[j*2+0]+cross[1]*e2[j*2+1]);
     for(int i=0;i<ndir;i++) {
       K[(os+ndir+os2+j)*ndof+os+i] += oownc * e[i];
       K[(os+i)*ndof+os+ndir+os2+j] += oownc * e[i];
     }
     for(int i=0;i<os;i++) {
       K[(os+ndir+os2+j)*ndof+i] += oownc * N[i];
       K[(i)*ndof+os+ndir+os2+j] += oownc * N[i];
     }
   }
   delete[] e2;

   double *u = new double[os2*2];
   deme2->polynomialF(x,u);

   for(int j=0;j<os2;j++) {
     double oownc = oown*(cross[0]*u[2*j+0]+cross[1]*u[2*j+1]);
     for(int i=0;i<ndir;i++) {
       K[(os+ndir+j)*ndof+os+i] += oownc * e[i];
       K[(os+i)*ndof+os+ndir+j] += oownc * e[i];
     }
     for(int i=0;i<os;i++) {
       K[(os+ndir+j)*ndof+i] += oownc * N[i];
       K[(i)*ndof+os+ndir+j] += oownc * N[i];
     }
   }
   delete[] u;
   delete[] e;
 }
};

void DEMHelm2d::interfMatrix(int fi, DEMElement* deme2,
                              complex<double> *K) {

 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double omega = getOmega();

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 int ndofs = nPolynomialDofs()+ndir+
             deme2->nPolynomialDofs()+deme2->nEnrichmentDofs();
 for(int i=0;i<ndofs*ndofs;i++) K[i] = 0.0;

 HelmDEMWetMatrixFunction f(o,omega,ndir,cdir,xref,deme2,K);
 ipu.lineInt2d(xyz, fi, f);

 delete[] cdir;

}


class HelmDEMNeumVFunction : public IntegFunctionL2d {
 int order;
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
public:
 HelmDEMNeumVFunction(int _order,
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc,
                      complex<double> *_v) {
   order = _order; ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*2+0]*(x[0]-xc[0]) +
                                  dirs[i*2+1]*(x[1]-xc[1]));

   complex<double> ince = w* (nsign* (incdir[0]*cross[0]+incdir[1]*cross[1])) *
                          exp(incdir[0]*x[0] + incdir[1]*x[1]);

   for(int i=0;i<order*order;i++) v[i] += ince*N[i];
   for(int i=0;i<ndir;i++)
     v[order*order+i] += ince*e[i];
   delete[] e; 
 }
};


void HelmDEMNeumV2d(int order, double *xy,
                     int ndir, complex<double> *dirs,
                     complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils2d ipu(order);

 HelmDEMNeumVFunction f(order,ndir,dirs,incdir,xc,v);
 ipu.lineInt2d(xy, faceindex, f);
}


void DEMHelm2d::createRHS(complex<double>*v) {
 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double kappa = getOmega()/getSpeedOfSound();
 double rho = getRho();

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 complex<double> *vv = new complex<double>[nPolynomialDofs()+ndir];
 for(int i=0;i<nPolynomialDofs()+ndir;i++) vv[i] = 0.0;

 double *dir = getWaveDirection();
 complex<double> incdir[2] = {complex<double>(0.0,kappa*dir[0]),
                              complex<double>(0.0,kappa*dir[1])};
 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     HelmDEMNeumV2d(o, xyz, ndir, cdir, incdir, i+1 , xref, vv);
 }
/*
 complex<double> *vv2 = new complex<double>[nPolynomialDofs()+ndir];
 for(int i=0;i<nPolynomialDofs()+ndir;i++) vv2[i] = 0.0;
 complex<double> incdir2[2] = {complex<double>(0.0,-kappa*1.0),
                              complex<double>(0.0,kappa*0.0)};
 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     HelmDEMNeumV2d(o, xyz, ndir, cdir, incdir2, i+1 , xref, vv2);
 }
 for(int i=0;i<nPolynomialDofs()+ndir;i++) vv[i] = vv[i]+
                                 0.93680538552900*vv2[i];
 delete[] vv2;
*/


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
