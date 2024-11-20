#include <cstdio>
#include <complex>

#include <Element.d/DEM.d/DEMLE2d.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>

template<class Scalar>
Scalar elasticForm(double lambda, double mu, double omega, double rho,
         Scalar (*gradu)[2], Scalar (*gradv)[2], Scalar *u, Scalar *v) {
 Scalar epsu[2][2], epsv[2][2];
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
   lambda*(gradu[0][0]+gradu[1][1])*(gradv[0][0]+gradv[1][1]) -
   omega*omega*rho*(u[0]*v[0]+u[1]*v[1]);
}


template<class Scalar>
void traction(double lambda, double mu, Scalar grad[2][2], double *n,
                     Scalar *t) {
Scalar e[2][2];
e[0][0] = grad[0][0];
e[1][1] = grad[1][1];
e[0][1] = 0.5*(grad[0][1]+grad[1][0]);
e[1][0] = e[0][1];

t[0] = 2.0*mu*(e[0][0]*n[0]+e[0][1]*n[1]) + lambda*(e[0][0]+e[1][1])*n[0];
t[1] = 2.0*mu*(e[1][0]*n[0]+e[1][1]*n[1]) + lambda*(e[0][0]+e[1][1])*n[1];
}

class SolidDGMEEFunction : public IntegFunctionA2d {
 double *materialc;
 int ndir;
 complex<double> *dirs;
 complex<double> *coef;
 double *xc;
 complex<double> *KEE;
public:
 SolidDGMEEFunction(double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef, 
                    double *_xc, complex<double> *_KEE) {
   materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   xc = _xc; KEE = _KEE;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {
   double &omega = materialc[1];
   double &lambda = materialc[2];
   double &mu = materialc[3];
   double &rho = materialc[4];
   double wdet = w*det;

   complex<double> (*u)[2] = new complex<double>[ndir][2];
   complex<double> (*gradu)[2][2] = new complex<double>[ndir][2][2];

   for(int k=0;k<ndir;k++) {
      complex <double> eu = exp(dirs[k*2+0]*(x[0]-xc[0])+
                                dirs[k*2+1]*(x[1]-xc[1]));
      u[k][0] = coef[k*2+0]*eu;
      u[k][1] = coef[k*2+1]*eu;
      gradu[k][0][0] =  coef[k*2+0]*dirs[k*2+0]*eu;
      gradu[k][0][1] =  coef[k*2+0]*dirs[k*2+1]*eu;
      gradu[k][1][0] =  coef[k*2+1]*dirs[k*2+0]*eu;
      gradu[k][1][1] =  coef[k*2+1]*dirs[k*2+1]*eu;
   }

   for(int l=0;l<ndir;l++)
     for(int k=0;k<ndir;k++) {
       KEE[l*ndir+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradu[k], gradu[l], u[k], u[l]);
   }
   
   delete[] u;
   delete[] gradu;
 }
};
                                                                                                                           

class SolidDGMELFunction : public IntegFunctionL2d {
 double *materialc;
 int ndir;
 complex<double> *dirs;
 complex<double> *coef;
 int nl;
 complex<double> *ldirs;
 complex<double> *lcoef;
 double *xsc;
 double *xc;
 complex<double> *L;
public:
 SolidDGMELFunction(double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef,
                    int _nl, complex<double> *_ldirs, complex<double> *_lcoef,
                    double *_xsc, double *_xc, complex<double> *__L) {
   materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   nl = _nl; ldirs = _ldirs; lcoef = _lcoef;
   xsc = _xsc; xc = _xc; L = __L;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   double &omega = materialc[1];

   complex<double> (*u)[2] = new complex<double>[ndir][2];
   complex<double> (*la)[2] = new complex<double>[nl][2];

   for (int k=0;k<ndir;k++) {
     complex<double> e = exp(dirs[k*2+0]*(x[0]-xc[0]) +
                             dirs[k*2+1]*(x[1]-xc[1]));
     u[k][0] =  coef[k*2+0]*e;
     u[k][1] =  coef[k*2+1]*e;
   };

   for (int l=0;l<nl;l++) {
     complex<double> eu = exp(ldirs[l*2+0]*(x[0]-xsc[0])+
                              ldirs[l*2+1]*(x[1]-xsc[1]));
     la[l][0] = lcoef[l*2+0]*eu;
     la[l][1] = lcoef[l*2+1]*eu;
   };

   for (int l=0;l<nl;l++) 
     for (int k=0;k<ndir;k++) {
      L[l*ndir+k] += w*omega*omega*(u[k][0]*la[l][0]+u[k][1]*la[l][1])*
                     sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
   }
   delete[] u;
   delete[] la;
 }
};


                                                     
void SolidDGMEMatrices2d(int order, double *xy, double *materialc,
                    int ndir, complex<double> *dirs, complex<double> *coef,
                    int *nldirs, complex<double> *ldirs, complex<double> *lcoef,
                    double *xsc, double *xc,
                    complex<double> *K, complex<double> *L) {
                                                                                                                           
 IsoParamUtils2d ipu(order);
 SolidDGMEEFunction f(materialc,ndir,dirs,coef,xc,K);
 ipu.areaInt2d(xy, f);

 int c = 0;
 for(int faceindex=1;faceindex<=4;faceindex++) 
   if (nldirs[faceindex-1]>0) {
     SolidDGMELFunction g(materialc,ndir,dirs,coef,
                          nldirs[faceindex-1],ldirs+c*2, lcoef+c*2,
                          xsc+(faceindex-1)*2,xc,L+c*ndir);
     c += nldirs[faceindex-1];
     ipu.lineInt2d(xy, faceindex, g);
  }
}


DGMLE2d::DGMLE2d(int _nnodes, int* nodenums) {
 nn = new int[_nnodes]; 
 for(int i=0;i<_nnodes;i++) nn[i] = nodenums[i];
 o = int(sqrt(double(_nnodes)));
 lm = new DEMLM*[nFaces()];
 for(int i=0;i<nFaces();i++) lm[i] = 0;
 bc = new int[nFaces()];
 for(int i=0;i<nFaces();i++) bc[i] = 0;
 ndir = 0;
}


void DGMLE2d::getRef(double *xyz,double *xy) {
 IsoParamUtils2d ipu(o);
// ipu.elementcenter(xyz,xy);
 xy[0] = xy[1] = 0.0;
}


void DGMLE2d::createM(complex<double>*M) {

 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 complex<double>* cndir = new complex<double>[ndir*2];
 for(int i=0;i<ndir;i++) {
   complex<double> d[4];
   dir(i,d);
   cdir[2*i+0] = d[0];
   cdir[2*i+1] = d[1];
   cndir[2*i+0] = d[2];
   cndir[2*i+1] = d[3];
 }

 int nldir[4] = {0,0,0,0};
 for(int fi=0;fi<4;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<4;i++) tnldir += nldir[i];

 double xlref[8] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<4;fi++) 
//   if (lm[fi]!=0) ipu.sidecenter(xyz,fi,xlref+fi*2);

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];

 for(int fi=0;fi<4;fi++) if (lm[fi]!=0) {
   DGMLE2d_LM *l = dynamic_cast<DGMLE2d_LM*>(lm[fi]);
   if (l==0) {
     fprintf(stderr,"DGMLE2d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMLE2d_LM type specific
 complex<double>* cldir = new complex<double>[2*tnldir];
 complex<double>* clndir = new complex<double>[2*tnldir];

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
     DGMLE2d_LM *l = dynamic_cast<DGMLE2d_LM*>(lm[i]);
     if (l!=0) {
       for(int j=0;j<lm[i]->nDofs();j++) {
         complex<double> d[4];
         l->ldir(j,tau,d);
         cldir[2*cc+0] = d[0];
         cldir[2*cc+1] = d[1];
         clndir[2*cc+0] = d[2];
         clndir[2*cc+1] = d[3];
         cc++;
       }
     }
   }
 }

 double materialc[5];
 materialc[1] = getOmega();
 double E = getE();
 double nu = getNu();
 materialc[2] = nu*E/((1.0+nu)*(1.0-2.0*nu));
 materialc[3] = E/(2.0*(1.0+nu));
 materialc[4] = getRho();

 SolidDGMEMatrices2d(o, xyz, materialc,
                    ndir, cdir, cndir, nldir, cldir, clndir,
                    xlref, xref, kee, kel);

 cc = 0;
 for(int i=0;i<4;i++) {
  for(int j=cc;j<cc+nldir[i];j++) for(int k=0;k<ndir;k++)
    kel[j*ndir+k] *= sign[i];
  cc += nldir[i];
 }

 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(tnldir+j)*(ndir+tnldir)+(tnldir+i)] = kee[j*ndir+i];
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(tnldir+i)*(ndir+tnldir)+j] = kel[j*ndir+i];
   M[j*(ndir+tnldir)+(tnldir+i)] = kel[j*ndir+i];
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++) 
   M[j*(ndir+tnldir)+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] xyz;
 delete[] cldir;
 delete[] clndir;
 delete[] cdir;
 delete[] cndir;
}


void DGMLE2d::enrichmentF(double *x, complex<double> *f) {

 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 delete[] xyz;

 for(int i=0;i<ndir;i++) {
   complex<double> d[4];
   dir(i,d);
   complex<double> e = exp(d[0]*(x[0]-xref[0]) + d[1]*(x[1]-xref[1]));
   f[i*2+0] = d[2]*e;
   f[i*2+1] = d[3]*e;
 }
}


void DGMLE2d::polynomialF(double *x, double *f) {


//  Find x in the coordinate system of the element using Newton's method
 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double xi[2] = {0,0};
 double *NN = new double[2*2*o];
 double *N = new double[3*os];

 while (1) {
   ipu.lagGalShapeFunction(2,xi,NN);
   for(int mx=0;mx<o;mx++) {
     for(int my=0;my<o;my++) {
       N[my*o+mx] = NN[0*o+mx]*NN[1*o+my];
       N[os+my*o+mx] = NN[2*o+0*o+mx]*NN[1*o+my];
       N[2*os+my*o+mx] = NN[0*o+mx]*NN[2*o+1*o+my];
     }
   }
   double J[2][2];
   ipu.jmatrix(N,xyz,J);
   double detj = ipu.detj(J);
   double invJ[2][2];
   ipu.invj(J,detj,invJ);

   double xx[2];
   xx[0] = -x[0];
   xx[1] = -x[1];
   for(int m=0;m<os;m++) {
     xx[0] += N[m]*xyz[m];
     xx[1] += N[m]*xyz[m+os];
   }

   double delta[2];

   delta[0] = -(invJ[0][0]*xx[0]+invJ[1][0]*xx[1]);
   delta[1] = -(invJ[0][1]*xx[0]+invJ[1][1]*xx[1]);

   xi[0] += delta[0];
   xi[1] += delta[1];

   if (fabs(delta[0])+fabs(delta[1])<1e-6) break;
 }

// Compute the values of shape functions

 ipu.lagGalShapeFunction(2,xi,NN);

 for(int i=0;i<2*os;i++) {
   int ii = i/2;
   double v = NN[0*o+ii%o]*NN[1*o+ii/o];
   f[i*2+0] = f[i*2+1] = 0;
   if (i%2==0) f[i*2+0] = v;
   else f[i*2+1] = v;
 }
 delete[] xyz;
 delete[] N;
 delete[] NN;
}


DGMLE2d_4::DGMLE2d_4(int _nnodes, int* nodenums) :
  DGMLE2d(_nnodes,nodenums) {
 ndir = 8;
}

DGMLE2d_16::DGMLE2d_16(int _nnodes, int* nodenums) :
  DGMLE2d(_nnodes,nodenums) {
 ndir = 32;
}

void DGMLE2d_4::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

// double a[][2] = { {1,0},{-1,0},{0,1},{0,-1} };
// double b[][2] = { {0,1},{0,1},{-1,0},{1,0} };
 double a[][2] = { {1,0},{0,1},{-1,0},{0,-1} };
 double b[][2] = { {0,-1},{1,0},{0,1},{-1,0} };
 if (n<4) {
   d[0] = complex<double>(0.0,kp*a[n][0]);
   d[1] = complex<double>(0.0,kp*a[n][1]);
   d[2] = a[n][0];
   d[3] = a[n][1];
 } else {
   n -= 4;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = b[n][0];
   d[3] = b[n][1];
 }
}

void DGMLE2d_16::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][2] = {
             {cos(2*M_PI*double(0)/16), sin(2*M_PI*double(0)/16)},
             {cos(2*M_PI*double(1)/16), sin(2*M_PI*double(1)/16)},
             {cos(2*M_PI*double(2)/16), sin(2*M_PI*double(2)/16)},
             {cos(2*M_PI*double(3)/16), sin(2*M_PI*double(3)/16)},
             {cos(2*M_PI*double(4)/16), sin(2*M_PI*double(4)/16)},
             {cos(2*M_PI*double(5)/16), sin(2*M_PI*double(5)/16)},
             {cos(2*M_PI*double(6)/16), sin(2*M_PI*double(6)/16)},
             {cos(2*M_PI*double(7)/16), sin(2*M_PI*double(7)/16)},
             {cos(2*M_PI*double(8)/16), sin(2*M_PI*double(8)/16)},
             {cos(2*M_PI*double(9)/16), sin(2*M_PI*double(9)/16)},
             {cos(2*M_PI*double(10)/16), sin(2*M_PI*double(10)/16)},
             {cos(2*M_PI*double(11)/16), sin(2*M_PI*double(11)/16)},
             {cos(2*M_PI*double(12)/16), sin(2*M_PI*double(12)/16)},
             {cos(2*M_PI*double(13)/16), sin(2*M_PI*double(13)/16)},
             {cos(2*M_PI*double(14)/16), sin(2*M_PI*double(14)/16)},
             {cos(2*M_PI*double(15)/16), sin(2*M_PI*double(15)/16)} };
 double b[][2] = {
             {-sin(2*M_PI*double(0)/16), cos(2*M_PI*double(0)/16)},
             {-sin(2*M_PI*double(1)/16), cos(2*M_PI*double(1)/16)},
             {-sin(2*M_PI*double(2)/16), cos(2*M_PI*double(2)/16)},
             {-sin(2*M_PI*double(3)/16), cos(2*M_PI*double(3)/16)},
             {-sin(2*M_PI*double(4)/16), cos(2*M_PI*double(4)/16)},
             {-sin(2*M_PI*double(5)/16), cos(2*M_PI*double(5)/16)},
             {-sin(2*M_PI*double(6)/16), cos(2*M_PI*double(6)/16)},
             {-sin(2*M_PI*double(7)/16), cos(2*M_PI*double(7)/16)},
             {-sin(2*M_PI*double(8)/16), cos(2*M_PI*double(8)/16)},
             {-sin(2*M_PI*double(9)/16), cos(2*M_PI*double(9)/16)},
             {-sin(2*M_PI*double(10)/16), cos(2*M_PI*double(10)/16)},
             {-sin(2*M_PI*double(11)/16), cos(2*M_PI*double(11)/16)},
             {-sin(2*M_PI*double(12)/16), cos(2*M_PI*double(12)/16)},
             {-sin(2*M_PI*double(13)/16), cos(2*M_PI*double(13)/16)},
             {-sin(2*M_PI*double(14)/16), cos(2*M_PI*double(14)/16)},
             {-sin(2*M_PI*double(15)/16), cos(2*M_PI*double(15)/16)} };
 if (n<16) {
   d[0] = complex<double>(0.0,kp*a[n][0]);
   d[1] = complex<double>(0.0,kp*a[n][1]);
   d[2] = a[n][0];
   d[3] = a[n][1];
 } else {
   n -= 16;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = b[n][0];
   d[3] = b[n][1];
 }
}


void DGMLE2d_1_LM::ldir(int n, double tau[2], complex<double> *d) {
 double E = e1->getE();
 if (n==0) {
   d[0] = 0.0;
   d[1] = 0.0;
   d[2] = E*1.0;
   d[3] = E*0.0;
 } else {
   d[0] = 0.0;
   d[1] = 0.0;
   d[2] = E*0.0;
   d[3] = E*1.0;
 }
}

void DGMLE2d_4_LM::ldir(int n, double tau[2], complex<double> *d) {
 double omega = e1->getOmega();
 double E = e1->getE();
 double nu = e1->getNu();
 double rho = e1->getRho();
 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 if (e2) {
// Max kp, ks from two elements
 }
 
 double normal[2] = { -tau[1], tau[0] };
 
}


class SolidNaturalEFunction : public IntegFunctionL2d {
 complex<double> *idir;
 complex<double> *cidir;
 double *materialc;
 int numW;
 complex<double> *dirs;
 complex<double> *coef;
 double *xc;
 complex<double> *v;
public:
 SolidNaturalEFunction(double* _materialc, int _numW, complex<double> *_dirs, complex<double> *_coef,
                  complex<double> *_idir, complex<double> *_cidir, double *_xc,
                  complex<double> *_v) {
   materialc = _materialc; numW = _numW; dirs = _dirs; coef = _coef;
   idir = _idir; cidir = _cidir; xc = _xc; v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   double lambda = materialc[2];
   double mu = materialc[3];
   complex<double> grad[2][2];
   complex<double> t[2];
   complex<double> e = exp(idir[0]*x[0]+idir[1]*x[1]);
   grad[0][0] = cidir[0]*idir[0]*e;
   grad[0][1] = cidir[0]*idir[1]*e;
   grad[1][0] = cidir[1]*idir[0]*e;
   grad[1][1] = cidir[1]*idir[1]*e;

   traction(lambda,mu,grad,cross,t);
   
   int ii;
   for(ii=0;ii<numW;ii++) {
     complex<double> e = exp(dirs[2*ii+0]*(x[0]-xc[0])+dirs[2*ii+1]*(x[1]-xc[1])); 
     complex<double> u[2];
     u[0] = coef[2*ii+0]*e;
     u[1] = coef[2*ii+1]*e;
     v[ii] += nsign*w*(t[0]*u[0]+t[1]*u[1]);
   };                      
 };
};


void SolidDGMNaturalEVector2d(int order, double *xy, double *materialc,
                     int numW, complex<double> *dirs, complex<double> *coef,
                     complex<double> *idir, complex<double> *cidir,
                     int faceindex, double *xc, complex<double> *v) {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 SolidNaturalEFunction f(materialc,numW,dirs,coef,idir,cidir,xc,v);
 ipu.lineInt2d(xy, faceindex, f);
}


void DGMLE2d::createRHS(complex<double>*v) {
 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 complex<double>* cndir = new complex<double>[ndir*2];
 for(int i=0;i<ndir;i++) {
   complex<double> d[4];
   dir(i,d);
   cdir[2*i+0] = d[0];
   cdir[2*i+1] = d[1];
   cndir[2*i+0] = d[2];
   cndir[2*i+1] = d[3];
 }

 complex<double> *vv = new complex<double>[ndir];
 for(int i=0;i<ndir;i++) vv[i] = 0.0;

 double materialc[5];
 materialc[1] = getOmega();
 double E = getE();
 double nu = getNu();
 materialc[2] = nu*E/((1.0+nu)*(1.0-2.0*nu));
 materialc[3] = E/(2.0*(1.0+nu));
 materialc[4] = getRho();

 double omega = materialc[1];
 double lambda = materialc[2];
 double mu = materialc[3];
 double rho = materialc[4];
 double kp = omega*sqrt(rho/(lambda+2*mu));

 complex<double> idir[2] = {complex<double>(0.0,kp*1.0),
                              complex<double>(0.0,0.0)};
// complex<double> cidir[2] = {complex<double>(0.0,1.218403954095460e-11),
//                              complex<double>(0.0,0.0)};
// complex<double> cidir[2] = {complex<double>(0.0,-0.93621651068154e-11),
//                              complex<double>(0.0,0.0)};
 complex<double> cidir[2] = {complex<double>(1.0,0.0),
                              complex<double>(0.0,0.0)};

 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     SolidDGMNaturalEVector2d(o, xyz, materialc, ndir, cdir, cndir, idir, cidir, i+1 , xref, vv);
;
 }

 delete[] cdir;
 delete[] xyz;

 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[i] = 0;
 for(int i=0;i<ndir;i++) 
   v[tnldir+i] = omega*omega*vv[i];
 delete[] vv;
}


void DGMLE2d::createSol(double *xyz,
                        complex<double>* sol, complex<double> *nodalSol) {

 for(int i=0;i<8;i++) nodalSol[i] = 0.0;
 for(int i=0;i<ndir;i++) {
   complex<double> d[4];
   dir(i,d);
   nodalSol[1] += sol[i]*d[2]*exp( d[0]*xyz[0]+d[1]*xyz[1] );
   nodalSol[2] += sol[i]*d[3]*exp( d[0]*xyz[0]+d[1]*xyz[1] );
 }
}


DEMLE2d::DEMLE2d(int _nnodes, int* nodenums) : DGMLE2d(_nnodes, nodenums)  { }

DEMLE2d_4::DEMLE2d_4(int _o, int* nodenums) : DEMLE2d(_o,nodenums) { ndir = 1; }

void DEMLE2d_4::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][2] = { {1,0},{-1,0},{0,1},{0,-1} };
 double b[][2] = { {0,1},{0,1},{-1,0},{1,0} };
 if (n<4) {
   d[0] = complex<double>(0.0,kp*a[n][0]);
   d[1] = complex<double>(0.0,kp*a[n][1]);
   d[2] = a[n][0];
   d[3] = a[n][1];
 } else {
   n -= 4;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = b[n][0];
   d[3] = b[n][1];
 }
}


class LEPELFunction2d : public IntegFunctionL2d {
 int o;
 double *materialc;
 int ndir;
 complex<double> *dirs;
 complex<double> *coef;
 int nl;
 complex<double> *ldirs;
 complex<double> *lcoef;
 double *xsc;
 double *xc;
 complex<double> *kel;
 complex<double> *kpl;
public:
 LEPELFunction2d(int _o, double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef,
                    int _nl, complex<double> *_ldirs, complex<double> *_lcoef,
                    double *_xsc, double *_xc,
                    complex<double> *_kel, complex<double> *_kpl) {
   o = _o; materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   nl = _nl; ldirs = _ldirs; lcoef = _lcoef;
   xsc = _xsc; xc = _xc; kel = _kel; kpl = _kpl;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   double &omega = materialc[1];

   complex<double> (*u)[2] = new complex<double>[ndir][2];
   complex<double> (*la)[2] = new complex<double>[nl][2];

   for (int k=0;k<ndir;k++) {
     complex<double> e = exp(dirs[k*2+0]*(x[0]-xc[0]) +
                             dirs[k*2+1]*(x[1]-xc[1]));
     u[k][0] =  coef[k*2+0]*e;
     u[k][1] =  coef[k*2+1]*e;
   };

   int tos = 2*o*o;
   complex<double> (*v)[2] = new complex<double>[tos][2];

   for(int k=0;k<tos;k+=2) {
     v[k][0] = N[k/2]; v[k][1] = 0.0;
     v[k+1][1] = N[k/2]; v[k+1][0] = 0.0;
   }

   for (int l=0;l<nl;l++) {
     complex<double> eu = exp(ldirs[l*2+0]*(x[0]-xsc[0])+
                              ldirs[l*2+1]*(x[1]-xsc[1]));
     la[l][0] = lcoef[l*2+0]*eu;
     la[l][1] = lcoef[l*2+1]*eu;
   };

   for (int l=0;l<nl;l++) {
     for (int k=0;k<ndir;k++)
      kel[l*ndir+k] += w*omega*omega*(u[k][0]*la[l][0]+u[k][1]*la[l][1])*
                     sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
     for (int k=0;k<tos;k++) 
      kpl[l*tos+k] += w*omega*omega*(v[k][0]*la[l][0]+v[k][1]*la[l][1])*
                     sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
   }
   delete[] v;
   delete[] u;
   delete[] la;
 }
};


class LEPEFunction2d : public IntegFunctionA2d {
 int o;
 double *materialc;
 int ndir;
 complex<double> *dirs;
 complex<double> *coef;
 double *xc;
 complex<double> *kee;
 complex<double> *kpe;
 complex<double> *kpp;
public:
 LEPEFunction2d(int _o, double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef, 
                    double *_xc, complex<double> *_kee,
                    complex<double> *_kpe, complex<double> *_kpp) {
   o = _o; materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   xc = _xc; kee = _kee; kpe = _kpe; kpp = _kpp;
 }
 void evaluate(double *x, double *N, double (*dNdx)[2], double w, double det) {
   double &omega = materialc[1];
   double &lambda = materialc[2];
   double &mu = materialc[3];
   double &rho = materialc[4];
   double wdet = w*det;

   complex<double> (*u)[2] = new complex<double>[ndir][2];
   complex<double> (*gradu)[2][2] = new complex<double>[ndir][2][2];

   for(int k=0;k<ndir;k++) {
      complex <double> eu = exp(dirs[k*2+0]*(x[0]-xc[0])+
                                dirs[k*2+1]*(x[1]-xc[1]));
      u[k][0] = coef[k*2+0]*eu;
      u[k][1] = coef[k*2+1]*eu;
      gradu[k][0][0] =  coef[k*2+0]*dirs[k*2+0]*eu;
      gradu[k][0][1] =  coef[k*2+0]*dirs[k*2+1]*eu;
      gradu[k][1][0] =  coef[k*2+1]*dirs[k*2+0]*eu;
      gradu[k][1][1] =  coef[k*2+1]*dirs[k*2+1]*eu;
   }

   int tos = 2*o*o;
   complex<double> (*v)[2] = new complex<double>[tos][2];
   complex<double> (*gradv)[2][2] = new complex<double>[tos][2][2];

   for(int k=0;k<tos;k+=2) {
     v[k][0] = N[k/2]; v[k][1] = 0.0;
     v[k+1][1] = N[k/2]; v[k+1][0] = 0.0;
     gradv[k][0][0] = dNdx[k/2][0];
     gradv[k][0][1] = dNdx[k/2][1];
     gradv[k][1][0] = gradv[k][1][1] = 0.0;
     gradv[k+1][1][0] = dNdx[k/2][0];
     gradv[k+1][1][1] = dNdx[k/2][1];
     gradv[k+1][0][0] = gradv[k+1][0][1] = 0.0;
   }

   for(int l=0;l<ndir;l++) {
     for(int k=0;k<ndir;k++)
       kee[l*ndir+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradu[k], gradu[l], u[k], u[l]);
     for(int k=0;k<tos;k++) 
       kpe[l*tos+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradv[k], gradu[l], v[k], u[l]);
   }
   for(int l=0;l<tos;l++) 
     for(int k=0;k<tos;k++) 
       kpp[l*tos+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradv[k], gradv[l], v[k], v[l]);
   
   delete[] v;
   delete[] gradv;
   delete[] u;
   delete[] gradu;
 }
};

void LEDEMMatrices2d(int order, double *xy, double *materialc,
                    int ndir, complex<double> *dirs, complex<double> *coef,
                    int *nldirs, complex<double> *ldirs, complex<double> *lcoef,
                    double *xsc, double *xc,
                    complex<double> *kee, complex<double> *kel,
                    complex<double> *kpp, complex<double> *kpl,
                    complex<double> *kpe) {

 IsoParamUtils2d ipu(order);
 int os = ipu.getordersq();
 int tos = 2*os;

 int nldir =nldirs[0] + nldirs[1] + nldirs[2] +
            nldirs[3];
 for(int i=0;i<nldir*ndir;i++) kel[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) kee[i] = 0.0;
 for(int i=0;i<tos*ndir;i++) kpe[i] = 0.0;
 for(int i=0;i<tos*nldir;i++) kpl[i] = 0.0;
 for(int i=0;i<tos*tos;i++) kpp[i] = 0.0;

 int c = 0;
 for(int faceindex=1;faceindex<=4;faceindex++) {
   LEPELFunction2d f(order,materialc,ndir,dirs,coef,
                nldirs[faceindex-1],ldirs+c*2,lcoef+c*2,
                xsc + (faceindex-1)*2, xc,
                kel+c*ndir,kpl+c*tos);
   ipu.lineInt2d(xy, faceindex, f);
   c += nldirs[faceindex-1];
 }

 LEPEFunction2d f(order,materialc,ndir,dirs, coef, xc, kee, kpe,kpp);
 ipu.areaInt2d(xy, f); 
}



void DEMLE2d::createM(complex<double>*M) {

 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 complex<double>* cndir = new complex<double>[ndir*2];
 for(int i=0;i<ndir;i++) {
   complex<double> d[4];
   dir(i,d);
   cdir[2*i+0] = d[0];
   cdir[2*i+1] = d[1];
   cndir[2*i+0] = d[2];
   cndir[2*i+1] = d[3];
 }

 int nldir[4] = {0,0,0,0};
 for(int fi=0;fi<4;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<4;i++) tnldir += nldir[i];

 double xlref[8] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<4;fi++) 
//   if (lm[fi]!=0) ipu.sidecenter(xyz,fi,xlref+fi*2);

 for(int fi=0;fi<4;fi++) if (lm[fi]!=0) {
   DGMLE2d_LM *l = dynamic_cast<DGMLE2d_LM*>(lm[fi]);
   if (l==0) {
     fprintf(stderr,"DGMLE2d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMLE2d_LM type specific
 complex<double>* cldir = new complex<double>[2*tnldir];
 complex<double>* clndir = new complex<double>[2*tnldir];

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
     DGMLE2d_LM *l = dynamic_cast<DGMLE2d_LM*>(lm[i]);
     if (l!=0) {
       for(int j=0;j<lm[i]->nDofs();j++) {
         complex<double> d[4];
         l->ldir(j,tau,d);
         cldir[2*cc+0] = d[0];
         cldir[2*cc+1] = d[1];
         clndir[2*cc+0] = d[2];
         clndir[2*cc+1] = d[3];
         cc++;
       }
     }
   }
 }

 double materialc[5];
 materialc[1] = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();
 materialc[2] = nu*E/((1.0+nu)*(1.0-2.0*nu));
 materialc[3] = E/(2.0*(1.0+nu));

 int tos = 2*os;

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];
 complex<double>* kpp = new complex<double>[tos*tos];
 complex<double>* kpe = new complex<double>[ndir*tos];
 complex<double>* kpl = new complex<double>[tos*tnldir];

 LEDEMMatrices2d(o, xyz, materialc,
                    ndir, cdir, cndir, nldir, cldir, clndir,
                    xlref, xref, kee, kel,kpp,kpl,kpe);

 cc = 0;
 for(int i=0;i<4;i++) {
  for(int j=cc;j<cc+nldir[i];j++) {
    for(int k=0;k<ndir;k++) kel[j*ndir+k] *= sign[i];
    for(int k=0;k<tos;k++) kpl[j*tos+k] *= sign[i];
  }
  cc += nldir[i];
 }

 for(int i=0;i<tos;i++) for(int j=0;j<tos;j++)
   M[j*(tos+tnldir+ndir)+i] = kpp[j*tos+i];
 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(tos+tnldir+j)*(tos+ndir+tnldir)+(tos+tnldir+i)] = kee[j*ndir+i];
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(tos+tnldir+i)*(tos+ndir+tnldir)+tos+j] = kel[j*ndir+i];
   M[(tos+j)*(tos+ndir+tnldir)+(tos+tnldir+i)] = kel[j*ndir+i];
 }
 for(int i=0;i<tos;i++) for(int j=0;j<ndir;j++) {
   M[i*(tos+ndir+tnldir)+(tos+tnldir+j)] = kpe[j*tos+i];
   M[(tos+tnldir+j)*(tos+ndir+tnldir)+i] = kpe[j*tos+i];
 }
 for(int i=0;i<tos;i++) for(int j=0;j<tnldir;j++) {
   M[i*(tos+ndir+tnldir)+(tos+j)] = kpl[j*tos+i];
   M[(tos+j)*(tos+ndir+tnldir)+i] = kpl[j*tos+i];
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++)
   M[(tos+j)*(tos+ndir+tnldir)+tos+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] kpp;
 delete[] kpe;
 delete[] kpl;

 delete[] xyz;
 delete[] cldir;
 delete[] clndir;
 delete[] cdir;
 delete[] cndir;
}



class LENaturalPEFunction2d : public IntegFunctionL2d {
 int o;
 complex<double> *idir;
 complex<double> *cidir;
 double *materialc;
 int ndir;
 complex<double> *dirs;
 complex<double> *coef;
 double *xc;
 complex<double> *v;
public:
 LENaturalPEFunction2d(int _o, double* _materialc,
                     int _ndir, complex<double> *_dirs, complex<double> *_coef,
                  complex<double> *_idir, complex<double> *_cidir, double *_xc,
                  complex<double> *_v) {
   o = _o; materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   idir = _idir; cidir = _cidir; xc = _xc; v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
   double lambda = materialc[2];
   double mu = materialc[3];
   complex<double> grad[2][2];
   complex<double> t[2];
   complex<double> e = exp(idir[0]*x[0]+idir[1]*x[1]);
   grad[0][0] = cidir[0]*idir[0]*e;
   grad[0][1] = cidir[0]*idir[1]*e;
   grad[1][0] = cidir[1]*idir[0]*e;
   grad[1][1] = cidir[1]*idir[1]*e;

   traction(lambda,mu,grad,cross,t);

   int tos = 2*o*o;
  
   for(int ii=0;ii<tos;ii+=2) {
     v[ii] += nsign*w*(t[0]*N[ii/2]);
     v[ii+1] += nsign*w*(t[1]*N[ii/2]);
   }
   
   for(int ii=0;ii<ndir;ii++) {
     complex<double> e = exp(dirs[2*ii+0]*(x[0]-xc[0])+dirs[2*ii+1]*(x[1]-xc[1])); 
     complex<double> u[2];
     u[0] = coef[2*ii+0]*e;
     u[1] = coef[2*ii+1]*e;
     v[tos+ii] += nsign*w*(t[0]*u[0]+t[1]*u[1]);
   };                      
 };
};


void LEDEMNaturalVector2d(int order, double *xy, double *materialc,
                     int numW, complex<double> *dirs, complex<double> *coef,
                     complex<double> *idir, complex<double> *cidir,
                     int faceindex, double *xc, complex<double> *v) {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 LENaturalPEFunction2d f(order,materialc,numW,dirs,coef,idir,cidir,xc,v);
 ipu.lineInt2d(xy, faceindex, f);
}


void DEMLE2d::createRHS(complex<double>*v) {
 IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 getNodalCoord(os,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 complex<double>* cndir = new complex<double>[ndir*2];
 for(int i=0;i<ndir;i++) {
   complex<double> d[4];
   dir(i,d);
   cdir[2*i+0] = d[0];
   cdir[2*i+1] = d[1];
   cndir[2*i+0] = d[2];
   cndir[2*i+1] = d[3];
 }

 double materialc[5];
 materialc[1] = getOmega();
 double E = getE();
 double nu = getNu();
 materialc[2] = nu*E/((1.0+nu)*(1.0-2.0*nu));
 materialc[3] = E/(2.0*(1.0+nu));
 materialc[4] = getRho();

 double omega = materialc[1];
 double lambda = materialc[2];
 double mu = materialc[3];
 double rho = materialc[4];
 double kp = omega*sqrt(rho/(lambda+2*mu));

 complex<double> idir[2] = {complex<double>(0.0,kp*1.0),
                              complex<double>(0.0,0.0)};
// complex<double> cidir[2] = {complex<double>(0.0,1.218403954095460e-11),
//                              complex<double>(0.0,0.0)};
// complex<double> cidir[2] = {complex<double>(0.0,-0.93621651068154e-11),
//                              complex<double>(0.0,0.0)};
 complex<double> cidir[2] = {complex<double>(0.0,1.0),
                              complex<double>(0.0,0.0)};

 complex<double> *vv = new complex<double>[nPolynomialDofs()+ndir];
 for(int i=0;i<nPolynomialDofs()+ndir;i++) vv[i] = 0.0;

 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     LEDEMNaturalVector2d(o, xyz, materialc,
                          ndir, cdir, cndir, idir, cidir, i+1 , xref, vv);
 }

 delete[] cdir;
 delete[] xyz;

 for(int i=0;i<nPolynomialDofs();i++) v[i] = omega*omega*vv[i];
 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[nPolynomialDofs()+i] = 0;
 for(int i=0;i<ndir;i++)
   v[nPolynomialDofs()+tnldir+i] = omega*omega*vv[nPolynomialDofs()+i];
 delete[] vv;
}
