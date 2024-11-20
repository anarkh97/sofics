#include <cstdio>
#include <complex>

#include <Element.d/DEM.d/DEMLE3d.h>
#include <Element.d/Helm.d/IsoParamUtils.h>

template<class Scalar>
Scalar elasticForm(double lambda, double mu, double omega, double rho,
         Scalar (*gradu)[3], Scalar (*gradv)[3], Scalar *u, Scalar *v) {
 Scalar epsu[3][3], epsv[3][3];
 epsu[0][0] = gradu[0][0];
 epsu[1][1] = gradu[1][1];
 epsu[2][2] = gradu[2][2];
 epsu[0][1] = 0.5*(gradu[0][1]+gradu[1][0]);
 epsu[1][0] = epsu[0][1];
 epsu[0][2] = 0.5*(gradu[0][2]+gradu[2][0]);
 epsu[2][0] = epsu[0][2];
 epsu[1][2] = 0.5*(gradu[1][2]+gradu[2][1]);
 epsu[2][1] = epsu[1][2];

 epsv[0][0] = gradv[0][0];
 epsv[1][1] = gradv[1][1];
 epsv[2][2] = gradv[2][2];
 epsv[0][1] = 0.5*(gradv[0][1]+gradv[1][0]);
 epsv[1][0] = epsv[0][1];
 epsv[0][2] = 0.5*(gradv[0][2]+gradv[2][0]);
 epsv[2][0] = epsv[0][2];
 epsv[1][2] = 0.5*(gradv[1][2]+gradv[2][1]);
 epsv[2][1] = epsv[1][2];

 return
   2.0*mu*(epsu[0][0]*epsv[0][0]+epsu[0][1]*epsv[0][1]+epsu[0][2]*epsv[0][2]+
           epsu[1][0]*epsv[1][0]+epsu[1][1]*epsv[1][1]+epsu[1][2]*epsv[1][2]+
           epsu[2][0]*epsv[2][0]+epsu[2][1]*epsv[2][1]+epsu[2][2]*epsv[2][2]) +
   lambda*(gradu[0][0]+gradu[1][1]+gradu[2][2])*
          (gradv[0][0]+gradv[1][1]+gradv[2][2]) -
   omega*omega*rho*(u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
}


void traction(double lambda, double mu,
              complex<double> *d, complex<double> *t, double *n,
              complex<double> *s) {
 complex<double> ltdotd =  lambda*(d[0]*t[0] + d[1]*t[1] + d[2]*t[2]);
 complex<double> td[3][3];
 td[0][0] = mu*(d[0]*t[0]+d[0]*t[0]) + ltdotd;
 td[0][1] = mu*(d[0]*t[1]+d[1]*t[0]);
 td[1][0] = td[0][1];
 td[0][2] = mu*(d[0]*t[2]+d[2]*t[0]);
 td[2][0] = td[0][2];
 td[1][1] = mu*(d[1]*t[1]+d[1]*t[1]) + ltdotd;
 td[1][2] = mu*(d[1]*t[2]+d[2]*t[1]);
 td[2][1] = td[1][2];
 td[2][2] = mu*(d[2]*t[2]+d[2]*t[2]) + ltdotd;
 s[0] = td[0][0]*n[0]+td[0][1]*n[1]+td[0][2]*n[2];
 s[1] = td[1][0]*n[0]+td[1][1]*n[1]+td[1][2]*n[2];
 s[2] = td[2][0]*n[0]+td[2][1]*n[1]+td[2][2]*n[2];
}


class SolidDGMEEFunction3d : public IntegFunctionV3d {
 double *materialc;
 int ndir;
 complex<double> *dirs;
 complex<double> *coef;
 double *xc;
 complex<double> *KEE;
public:
 SolidDGMEEFunction3d(double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef, 
                    double *_xc, complex<double> *_KEE) {
   materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   xc = _xc; KEE = _KEE;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
   double &omega = materialc[1];
   double &lambda = materialc[2];
   double &mu = materialc[3];
   double &rho = materialc[4];
   double wdet = w*det;

   complex<double> (*u)[3] = new complex<double>[ndir][3];
   complex<double> (*gradu)[3][3] = new complex<double>[ndir][3][3];


   for(int k=0;k<ndir;k++) {
      complex <double> eu = exp(dirs[k*3+0]*(x[0]-xc[0])+
                                dirs[k*3+1]*(x[1]-xc[1])+
                                dirs[k*3+2]*(x[2]-xc[2]));
      u[k][0] = coef[k*3+0]*eu;
      u[k][1] = coef[k*3+1]*eu;
      u[k][2] = coef[k*3+2]*eu;
      gradu[k][0][0] =  coef[k*3+0]*dirs[k*3+0]*eu;
      gradu[k][0][1] =  coef[k*3+0]*dirs[k*3+1]*eu;
      gradu[k][0][2] =  coef[k*3+0]*dirs[k*3+2]*eu;
      gradu[k][1][0] =  coef[k*3+1]*dirs[k*3+0]*eu;
      gradu[k][1][1] =  coef[k*3+1]*dirs[k*3+1]*eu;
      gradu[k][1][2] =  coef[k*3+1]*dirs[k*3+2]*eu;
      gradu[k][2][0] =  coef[k*3+2]*dirs[k*3+0]*eu;
      gradu[k][2][1] =  coef[k*3+2]*dirs[k*3+1]*eu;
      gradu[k][2][2] =  coef[k*3+2]*dirs[k*3+2]*eu;
   }

   for(int l=0;l<ndir;l++)
     for(int k=0;k<ndir;k++) {
       KEE[l*ndir+k] += wdet*omega*omega*
           elasticForm(lambda, mu, omega, rho, gradu[k], gradu[l], u[k], u[l]);
   }
   delete[] u;
   delete[] gradu;
 }
};
                                                                                                                           

class SolidDGMELFunction3d : public IntegFunctionA3d {
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
 SolidDGMELFunction3d(double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef,
                    int _nl, complex<double> *_ldirs, complex<double> *_lcoef,
                    double *_xsc, double *_xc, complex<double> *__L) {
   materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   nl = _nl; ldirs = _ldirs; lcoef = _lcoef;
   xsc = _xsc; xc = _xc; L = __L;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   double &omega = materialc[1];

   complex<double> (*u)[3] = new complex<double>[ndir][3];
   complex<double> (*la)[3] = new complex<double>[nl][3];

   for (int k=0;k<ndir;k++) {
     complex<double> e = exp(dirs[k*3+0]*(x[0]-xc[0]) +
                             dirs[k*3+1]*(x[1]-xc[1]) +
                             dirs[k*3+2]*(x[2]-xc[2]));
     u[k][0] =  coef[k*3+0]*e;
     u[k][1] =  coef[k*3+1]*e;
     u[k][2] =  coef[k*3+2]*e;
   };

   for (int l=0;l<nl;l++) {
     complex<double> eu = exp(ldirs[l*3+0]*(x[0]-xsc[0])+
                              ldirs[l*3+1]*(x[1]-xsc[1])+
                              ldirs[l*3+2]*(x[2]-xsc[2]));
     la[l][0] = lcoef[l*3+0]*eu;
     la[l][1] = lcoef[l*3+1]*eu;
     la[l][2] = lcoef[l*3+2]*eu;
   };

   for (int l=0;l<nl;l++) 
     for (int k=0;k<ndir;k++) {
      L[l*ndir+k] += w*omega*omega*
                    (u[k][0]*la[l][0]+u[k][1]*la[l][1]+u[k][2]*la[l][2])*
                sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
   }
   delete[] u;
   delete[] la;
 }
};


                                                     
void SolidDGMEMatrices3d(int order, double *xy, double *materialc,
                    int ndir, complex<double> *dirs, complex<double> *coef,
                    int *nldirs, complex<double> *ldirs, complex<double> *lcoef,
                    double *xsc, double *xc,
                    complex<double> *K, complex<double> *L) {
                                                                                                                           
 IsoParamUtils ipu(order);
 SolidDGMEEFunction3d f(materialc,ndir,dirs,coef,xc,K);
 ipu.volumeInt3d(xy, f);

 int c = 0;
 for(int faceindex=1;faceindex<=6;faceindex++) 
   if (nldirs[faceindex-1]>0) {
     SolidDGMELFunction3d g(materialc,ndir,dirs,coef,
                          nldirs[faceindex-1],ldirs+c*3, lcoef+c*3,
                          xsc+(faceindex-1)*3,xc,L+c*ndir);
     c += nldirs[faceindex-1];
     ipu.surfInt3d(xy, faceindex, g);
  }
}


DGMLE3d::DGMLE3d(int _nnodes, int* nodenums) {
 nn = new int[_nnodes]; 
 for(int i=0;i<_nnodes;i++) nn[i] = nodenums[i];
 o = int(pow(double(_nnodes),1.0/3.0));
 lm = new DEMLM*[nFaces()];
 for(int i=0;i<nFaces();i++) lm[i] = 0;
 bc = new int[nFaces()];
 for(int i=0;i<nFaces();i++) bc[i] = 0;
 ndir = 0;
}


void DGMLE3d::getRef(double *xyz,double *cxyz) {
 IsoParamUtils ipu(o);
// ipu.elementcenter(xyz,xy);
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
}


void DGMLE3d::createM(complex<double>*M) {

 IsoParamUtils ipu(o);
 int os = ipu.getordersq();
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 complex<double>* cndir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) {
   complex<double> d[6];
   dir(i,d);
   cdir[3*i+0] = d[0];
   cdir[3*i+1] = d[1];
   cdir[3*i+2] = d[2];
   cndir[3*i+0] = d[3];
   cndir[3*i+1] = d[4];
   cndir[3*i+2] = d[5];
 }

 int nldir[6] = {0,0,0,0,0,0};
 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<6;i++) tnldir += nldir[i];

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];

 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) {
   DGMLE3d_LM *l = dynamic_cast<DGMLE3d_LM*>(lm[fi]);
   if (l==0) {
     fprintf(stderr,"DGMLE3d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMLE3d_LM type specific
 double xlref[18] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<6;fi++)
// ipu.sidecenter(xyz,fi+1,xlref+fi*3);

 complex<double>* cldir = new complex<double>[3*tnldir];
 complex<double>* clndir = new complex<double>[3*tnldir];

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
     DGMLE3d_LM *l = dynamic_cast<DGMLE3d_LM*>(lm[i]);
     if (l!=0) {
       for(int j=0;j<lm[i]->nDofs();j++) {
         complex<double> d[6];
         l->ldir(j,tau1,tau2,d);
         cldir[2*cc+0] = d[0];
         cldir[2*cc+1] = d[1];
         cldir[2*cc+2] = d[2];
         clndir[3*cc+0] = d[3];
         clndir[3*cc+1] = d[4];
         clndir[3*cc+2] = d[5];
         cc++;
       }
     }
   }
 }

 double materialc[5];
 materialc[1] = getOmega();
 double E = getE();
 double nu = getNu();
 materialc[4] = getRho();
 materialc[2] = nu*E/((1.0+nu)*(1.0-2.0*nu));
 materialc[3] = E/(2.0*(1.0+nu));

 SolidDGMEMatrices3d(o, xyz, materialc,
                    ndir, cdir, cndir, nldir, cldir, clndir,
                    xlref, xref, kee, kel);

 cc = 0;
 for(int i=0;i<6;i++) {
  for(int j=cc;j<cc+nldir[i];j++) for(int k=0;k<ndir;k++)
    kel[j*ndir+k] *= sign[i];
  cc += nldir[i];
 }

 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++) {
   M[(tnldir+j)*(ndir+tnldir)+(tnldir+i)] = kee[j*ndir+i];
 }
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


void DGMLE3d::enrichmentF(double *x, complex<double> *f) {
 IsoParamUtils ipu(o);
 int os = ipu.getordersq();
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double xref[3];
 getRef(xyz,xref);

 delete[] xyz;

 for(int i=0;i<ndir;i++) {
   complex<double> d[6];
   dir(i,d);
   complex<double> e =
      exp(d[0]*(x[0]-xref[0]) + d[1]*(x[1]-xref[1]) +  d[2]*(x[2]-xref[2]));
   f[i*3+0] = d[3]*e;
   f[i*3+1] = d[4]*e;
   f[i*3+2] = d[5]*e;
 }
}


void DGMLE3d::polynomialF(double *x, double *f) {


//  Find x in the coordinate system of the element using Newton's method
 IsoParamUtils ipu(o);
 int os = ipu.getordersq();
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double xi[3] = {0,0,0};
 double *NN = new double[3*2*o];
 double *N = new double[4*oc];

 while (1) {
   ipu.lagGalShapeFunction(3,xi,NN);
   for(int mx=0;mx<o;mx++) {
     for(int my=0;my<o;my++) {
       for(int mz=0;mz<o;mz++) {
         N[mz*o*o+my*o+mx] = NN[0*o+mx]*NN[1*o+my]*NN[2*o+mz];
         N[oc+mz*o*o+my*o+mx] = NN[3*o+0*o+mx]*NN[1*o+my]*NN[2*o+mz];
         N[2*oc+mz*o*o+my*o+mx] = NN[0*o+mx]*NN[3*o+1*o+my]*NN[2*o+mz];
         N[3*oc+mz*o*o+my*o+mx] = NN[0*o+mx]*NN[1*o+my]*NN[3*o+2*o+mz];
       }
     }
   }
   double J[3][3];
   ipu.jmatrix(N,xyz,J);
   double detj = ipu.detj(J);
   double invJ[3][3];
   ipu.invj(J,detj,invJ);

   double xx[3];
   xx[0] = -x[0];
   xx[1] = -x[1];
   xx[2] = -x[2];
   for(int m=0;m<oc;m++) {
     xx[0] += N[m]*xyz[m];
     xx[1] += N[m]*xyz[m+oc];
     xx[2] += N[m]*xyz[m+2*oc];
   }

   double delta[3];

   delta[0] = -(invJ[0][0]*xx[0]+invJ[1][0]*xx[1]+invJ[2][0]*xx[2]);
   delta[1] = -(invJ[0][1]*xx[0]+invJ[1][1]*xx[1]+invJ[2][1]*xx[2]);
   delta[2] = -(invJ[0][2]*xx[0]+invJ[1][2]*xx[1]+invJ[2][2]*xx[2]);

   xi[0] += delta[0];
   xi[1] += delta[1];
   xi[2] += delta[2];

   if (fabs(delta[0])+fabs(delta[1])+fabs(delta[2])<1e-6) break;
 }

// Compute the values of shape functions

 ipu.lagGalShapeFunction(3,xi,NN);

 for(int i=0;i<3*oc;i++) {
   int ii = i/3;
   double v = NN[0*o+(ii%(o*o))%o]*NN[1*o+(ii%(o*o))/o]*N[2*o+ii/(o*o)];

   f[i*3+0] = f[i*3+1] = f[i*3+2] = 0;
   if (i%3==0) f[i*3+0] = v;
   else if (i%3==1) f[i*3+1] = v;
   else f[i*3+2] = v;
 }
 delete[] xyz;
 delete[] N;
 delete[] NN;
}


DGMLE3d_6::DGMLE3d_6(int _nnodes, int* nodenums) :
  DGMLE3d(_nnodes,nodenums) {
 ndir = 18;
}


void DGMLE3d_6::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][3] = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0}, {0,0,1}, {0,0,-1} };
 double b[][3] = { {0,1,0},{0,1,0},{1,0,0},{1,0,0},{1,0,0},{1,0,0} };
 double c[][3] = { {0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,1,0},{0,1,0} };
 if (n<6) {
   d[0] = complex<double>(0.0,kp*a[n][0]);
   d[1] = complex<double>(0.0,kp*a[n][1]);
   d[2] = complex<double>(0.0,kp*a[n][2]);
   d[3] = a[n][0];
   d[4] = a[n][1];
   d[5] = a[n][2];
 } else if (n<12) {
   n -= 6;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = complex<double>(0.0,ks*a[n][2]);
   d[3] = b[n][0];
   d[4] = b[n][1];
   d[5] = b[n][2];
 } else {
   n -= 12;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = complex<double>(0.0,ks*a[n][2]);
   d[3] = c[n][0];
   d[4] = c[n][1];
   d[5] = c[n][2];
 }
}


DGMLE3d_26::DGMLE3d_26(int _nnodes, int* nodenums) :
  DGMLE3d(_nnodes,nodenums) {
 ndir = 26*3;
}


void DGMLE3d_26::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, -1, 0}, {1, 0, -1}, {0, 1, -1} };
 double b[][3] = { {0, 1, 0}, {1, 0, 0}, {1, 0, 0}, {0, -1, 1}, {0, 1, 1}, {0, 1, 1}, {0, 1, -1}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0} };
 double k;
 switch (n%6) {
  case 0: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,kp*a[n][0]/l);
    d[1] = complex<double>(0.0,kp*a[n][1]/l);
    d[2] = complex<double>(0.0,kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 1: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,-kp*a[n][0]/l);
    d[1] = complex<double>(0.0,-kp*a[n][1]/l);
    d[2] = complex<double>(0.0,-kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 2: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 3: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 4: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
  }
  case 5: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
   }
 }
}

DGMLE3d_50::DGMLE3d_50(int _nnodes, int* nodenums) :
  DGMLE3d(_nnodes,nodenums) {
 ndir = 50*3;
}


void DGMLE3d_50::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][3] = {
   {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1},
   {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 0},
   {1, 0, 1}, {0, 1, 1}, {1, -1, 0}, {1, 0, -1}, {0, 1, -1},
   {1, -0.5, -0.5}, {1, -0.5, 0.5}, {1, 0.5, -0.5}, {1, 0.5, 0.5},
   {-0.5, 1, -0.5}, {-0.5, 1, 0.5}, {0.5, 1, -0.5}, {0.5, 1, 0.5},
   {-0.5, -0.5, 1}, {-0.5, 0.5, 1}, {0.5, -0.5, 1}, {0.5, 0.5, 1}
 };
 double b[][3] = {
   {0, 1, 0}, {1, 0, 0}, {1, 0, 0}, {0, -1, 1},
   {0, 1, 1}, {0, 1, 1}, {0, 1, -1}, {0, 0, 1},
   {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0},
   {0, 1, -1}, {0, 1, 1}, {0, 1, 1}, {0, 1, -1},
   {1, 0, -1}, {1, 0, 1}, {1, 0, 1}, {1, 0, -1},
   {1, -1, 0}, {1, 1, 0}, {1, 1, 0}, {1, -1, 0}
 };
 double k;
 double aa[3], bb[3];
 switch (n%6) {
  case 0: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,kp*a[n][0]/l);
    d[1] = complex<double>(0.0,kp*a[n][1]/l);
    d[2] = complex<double>(0.0,kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 1: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,-kp*a[n][0]/l);
    d[1] = complex<double>(0.0,-kp*a[n][1]/l);
    d[2] = complex<double>(0.0,-kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 2: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 3: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 4: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
  }
  case 5: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
   }
 }
}


void DGMLE3d_3_LM::ldir(int n, double *tau1, double *tau2, complex<double> *d) {
 if (n==0) {
   d[0] = 0.0; d[1] = 0.0; d[2] = 0.0; d[3] = 1.0; d[4] = 0.0; d[5] = 0.0;
 } else if (n==1) {
   d[0] = 0.0; d[1] = 0.0; d[2] = 0.0; d[3] = 0.0; d[4] = 1.0; d[5] = 0.0;
 } else {
   d[0] = 0.0; d[1] = 0.0; d[2] = 0.0; d[3] = 0.0; d[4] = 0.0; d[5] = 1.0;
 }
}

void DGMLE3d_15_LM::ldir(int n, double *tau1, double *tau2, complex<double> *d) {
 double omega = e1->getOmega();
 double E = (e1->getProperty()->E+e2->getProperty()->E)/2.0;
 double nu = (e1->getProperty()->nu+e2->getProperty()->nu)/2.0;
 double rho = (e1->getProperty()->rho+e2->getProperty()->rho)/2.0;

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double no[3];
 no[0] = tau1[1]*tau1[2]-tau1[2]*tau1[1];
 no[1] = tau1[2]*tau1[0]-tau1[0]*tau1[2];
 no[2] = tau1[0]*tau1[1]-tau1[1]*tau1[0];

 double a = sqrt(2.0);
 double localdirLM[][3] = {
   {0, 0, 1}, {0, 0, 1}, {0, 0, 1},
   {1, 1, a}, {-1, -1, a}, {1, -1, a}, {-1, 1, a},
   {1, 1, a}, {-1, -1, a}, {1, -1, a}, {-1, 1, a},
   {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0} 
 };
 double localcdirLM[][3] = {
   {0, 0, 1}, {0, 1, 0}, {1, 0, 0},
   {1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {1, 1, 0},
   {a, a, -2}, {a, a, 2}, {a, -a, -2}, {a, -a, 2},
   {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1},
 };
 double l = sqrt(localdirLM[n][0]*localdirLM[n][0]+
                 localdirLM[n][1]*localdirLM[n][1]+
                 localdirLM[n][2]*localdirLM[n][2]);
 double lc = sqrt(localcdirLM[n][0]*localcdirLM[n][0]+
                  localcdirLM[n][1]*localcdirLM[n][1]+
                  localcdirLM[n][2]*localcdirLM[n][2]);
 complex<double> lmd[3]={0,0,0};
 complex<double> lmdc[3]={0,0,0};
 if (n<3) {
   for(int i=0;i<3;i++) d[i] = lmd[i] = complex<double>(0.0,kp)*
     (localdirLM[n][0]*tau1[i]+localdirLM[n][1]*tau2[i]+localdirLM[n][0]*no[i])/l;
   for(int i=0;i<3;i++) lmdc[i] = 
 (localcdirLM[n][0]*tau1[i]+localcdirLM[n][1]*tau2[i]+localcdirLM[n][0]*no[i])/lc;
 } else {
   for(int i=0;i<3;i++) d[i] = lmd[i] = complex<double>(0.0,ks)*
     (localdirLM[n][0]*tau1[i]+localdirLM[n][1]*tau2[i]+localdirLM[n][0]*no[i])/l;
   for(int i=0;i<3;i++) lmdc[i] = 
 (localcdirLM[n][0]*tau1[i]+localcdirLM[n][1]*tau2[i]+localcdirLM[n][0]*no[i])/lc;
 }

 traction(lambda,mu,lmd,lmdc,no,d+3);
}


void DGMLE3d_28_LM::ldir(int n, double *tau1, double *tau2, complex<double> *d) {
 double omega = e1->getOmega();
 double E = (e1->getProperty()->E+e2->getProperty()->E)/2.0;
 double nu = (e1->getProperty()->nu+e2->getProperty()->nu)/2.0;
 double rho = (e1->getProperty()->rho+e2->getProperty()->rho)/2.0;

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double no[3];
 no[0] = tau1[1]*tau1[2]-tau1[2]*tau1[1];
 no[1] = tau1[2]*tau1[0]-tau1[0]*tau1[2];
 no[2] = tau1[0]*tau1[1]-tau1[1]*tau1[0];

 double a= 1;
 double b= 2.5;
 double localdirLM[][3] = {
   {1, 0, b}, 
   {-1, 0, b}, 
   {0, 1, b}, 
   {0, -1, b}, 
   {1, 0, -b}, 
   {-1, 0, -b}, 
   {0, 1, -b}, 
   {0, -1, -b}, 

   {1, 0, a}, 
   {-1, 0, a}, 
   {0, 1, a}, 
   {0, -1, a}, 

   {1, 0, a}, 
   {-1, 0, a}, 
   {0, 1, a}, 
   {0, -1, a}, 

   {1, 1, sqrt(2.0)*a},
   {-1, -1, sqrt(2.0)*a},
   {1, -1, sqrt(2.0)*a},
   {-1, 1, sqrt(2.0)*a},

   {1, 1, sqrt(2.0)*a},
   {-1, -1, sqrt(2.0)*a},
   {1, -1, sqrt(2.0)*a},
   {-1, 1, sqrt(2.0)*a},

   {1, 1, 0},
   {-1, 1, 0},
   {1, -1, 0},
   {-1, -1, 0},
 };
 double localcdirLM[][3] = {
   {1, 0, b}, 
   {-1, 0, b}, 
   {0, 1, b}, 
   {0, -1, b}, 
   {1, 0, -b}, 
   {-1, 0, -b}, 
   {0, 1, -b}, 
   {0, -1, -b}, 

   {0, 1, 0}, 
   {0, 1, 0}, 
   {1, 0, 0}, 
   {1, 0, 0}, 

   {-a, 0, 1}, 
   {a, 0, 1}, 
   {0, -a, 1}, 
   {0, a, 1}, 
   
   {1, -1, 0},
   {1, -1, 0},
   {1, 1, 0},
   {1, 1, 0},

   {sqrt(2.0)*a, sqrt(2.0)*a, -2},
   {sqrt(2.0)*a, sqrt(2.0)*a, 2},
   {sqrt(2.0)*a, -sqrt(2.0)*a, -2},
   {sqrt(2.0)*a, -sqrt(2.0)*a, 2},

   {0, 0, 1},
   {0, 0, 1},
   {0, 0, 1},
   {0, 0, 1},
 };

 double l = sqrt(localdirLM[n][0]*localdirLM[n][0]+
                 localdirLM[n][1]*localdirLM[n][1]+
                 localdirLM[n][2]*localdirLM[n][2]);
 double lc = sqrt(localcdirLM[n][0]*localcdirLM[n][0]+
                  localcdirLM[n][1]*localcdirLM[n][1]+
                  localcdirLM[n][2]*localcdirLM[n][2]);
 complex<double> lmd[3]={0,0,0};
 complex<double> lmdc[3]={0,0,0};
 if (n<8) {
   for(int i=0;i<3;i++) d[i] = lmd[i] = complex<double>(0.0,kp)*
     (localdirLM[n][0]*tau1[i]+localdirLM[n][1]*tau2[i]+localdirLM[n][0]*no[i])/l;
   for(int i=0;i<3;i++) lmdc[i] = 
 (localcdirLM[n][0]*tau1[i]+localcdirLM[n][1]*tau2[i]+localcdirLM[n][0]*no[i])/lc;
 } else {
   for(int i=0;i<3;i++) d[i] = lmd[i] = complex<double>(0.0,ks)*
     (localdirLM[n][0]*tau1[i]+localdirLM[n][1]*tau2[i]+localdirLM[n][0]*no[i])/l;
   for(int i=0;i<3;i++) lmdc[i] = 
 (localcdirLM[n][0]*tau1[i]+localcdirLM[n][1]*tau2[i]+localcdirLM[n][0]*no[i])/lc;
 }

 traction(lambda,mu,lmd,lmdc,no,d+3);
}


void DGMLE3d::createRHS(complex<double>*v) {
/* IsoParamUtils2d ipu(o);
 int os = ipu.getordersq();
 double *xyz= new double[3*os];
 cs.getCoordinates(nn,os,xyz,xyz+os,xyz+2*os);

 double kappa = prop ->kappaHelm;

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*2];
 for(int i=0;i<2*ndir;i++) cdir[i] = dir(i);

 complex<double> incdir[2] = {complex<double>(0.0,kappa*1.0),
                              complex<double>(0.0,kappa*0.0)};

 complex<double> *vv = new complex<double>[ndir];
 for(int i=0;i<ndir;i++) vv[i] = 0.0;

 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     HelmDGMENeumV2d(o, xyz, ndir, cdir, kappa, incdir, i+1 , xref, vv);
;
 }

 delete[] cdir;
 delete[] xyz;

 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[i] = 0;
 for(int i=0;i<ndir;i++) 
   v[tnldir+i] = vv[i];
 delete[] vv;*/
}


void DGMLE3d::createSol(double *xyz,
                        complex<double>* sol, complex<double> *nodalSol) {

 for(int i=0;i<8;i++) nodalSol[i] = 0.0;
 for(int i=0;i<ndir;i++) {
   complex<double> d[6];
   dir(i,d);
   nodalSol[1] += sol[i]*d[3]*exp( d[0]*xyz[0]+d[1]*xyz[1]+d[2]*xyz[2] );
   nodalSol[2] += sol[i]*d[4]*exp( d[0]*xyz[0]+d[1]*xyz[1]+d[2]*xyz[2] );
   nodalSol[3] += sol[i]*d[5]*exp( d[0]*xyz[0]+d[1]*xyz[1]+d[2]*xyz[2] );
 }
}


DEMLE3d::DEMLE3d(int _nnodes, int* nodenums) : DGMLE3d(_nnodes, nodenums)  { }


DEMLE3d_6::DEMLE3d_6(int _nnodes, int* nodenums) :
  DEMLE3d(_nnodes,nodenums) {
 ndir = 18;
}


void DEMLE3d_6::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][3] = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0}, {0,0,1}, {0,0,-1} };
 double b[][3] = { {0,1,0},{0,1,0},{1,0,0},{1,0,0},{1,0,0},{1,0,0} };
 double c[][3] = { {0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,1,0},{0,1,0} };
 if (n<6) {
   d[0] = complex<double>(0.0,kp*a[n][0]);
   d[1] = complex<double>(0.0,kp*a[n][1]);
   d[2] = complex<double>(0.0,kp*a[n][2]);
   d[3] = a[n][0];
   d[4] = a[n][1];
   d[5] = a[n][2];
 } else if (n<12) {
   n -= 6;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = complex<double>(0.0,ks*a[n][2]);
   d[3] = b[n][0];
   d[4] = b[n][1];
   d[5] = b[n][2];
 } else {
   n -= 12;
   d[0] = complex<double>(0.0,ks*a[n][0]);
   d[1] = complex<double>(0.0,ks*a[n][1]);
   d[2] = complex<double>(0.0,ks*a[n][2]);
   d[3] = c[n][0];
   d[4] = c[n][1];
   d[5] = c[n][2];
 }
}


DEMLE3d_26::DEMLE3d_26(int _nnodes, int* nodenums) :
  DEMLE3d(_nnodes,nodenums) {
 ndir = 26*3;
}


void DEMLE3d_26::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, -1, 0}, {1, 0, -1}, {0, 1, -1} };
 double b[][3] = { {0, 1, 0}, {1, 0, 0}, {1, 0, 0}, {0, -1, 1}, {0, 1, 1}, {0, 1, 1}, {0, 1, -1}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0} };
 double k;
 switch (n%6) {
  case 0: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,kp*a[n][0]/l);
    d[1] = complex<double>(0.0,kp*a[n][1]/l);
    d[2] = complex<double>(0.0,kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 1: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,-kp*a[n][0]/l);
    d[1] = complex<double>(0.0,-kp*a[n][1]/l);
    d[2] = complex<double>(0.0,-kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 2: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 3: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 4: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
  }
  case 5: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
   }
 }
}

DEMLE3d_50::DEMLE3d_50(int _nnodes, int* nodenums) :
  DEMLE3d(_nnodes,nodenums) {
 ndir = 50*3;
}


void DEMLE3d_50::dir(int n, complex<double> *d) {
 double omega = getOmega();
 double E = getE();
 double nu = getNu();
 double rho = getRho();

 double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
 double mu = E/(2.0*(1.0+nu));
 double kp = omega*sqrt(rho/(lambda+2*mu));
 double ks = omega*sqrt(rho/mu);

 double a[][3] = {
   {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1},
   {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 0},
   {1, 0, 1}, {0, 1, 1}, {1, -1, 0}, {1, 0, -1}, {0, 1, -1},
   {1, -0.5, -0.5}, {1, -0.5, 0.5}, {1, 0.5, -0.5}, {1, 0.5, 0.5},
   {-0.5, 1, -0.5}, {-0.5, 1, 0.5}, {0.5, 1, -0.5}, {0.5, 1, 0.5},
   {-0.5, -0.5, 1}, {-0.5, 0.5, 1}, {0.5, -0.5, 1}, {0.5, 0.5, 1}
 };
 double b[][3] = {
   {0, 1, 0}, {1, 0, 0}, {1, 0, 0}, {0, -1, 1},
   {0, 1, 1}, {0, 1, 1}, {0, 1, -1}, {0, 0, 1},
   {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0},
   {0, 1, -1}, {0, 1, 1}, {0, 1, 1}, {0, 1, -1},
   {1, 0, -1}, {1, 0, 1}, {1, 0, 1}, {1, 0, -1},
   {1, -1, 0}, {1, 1, 0}, {1, 1, 0}, {1, -1, 0}
 };
 double k;
 double aa[3], bb[3];
 switch (n%6) {
  case 0: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,kp*a[n][0]/l);
    d[1] = complex<double>(0.0,kp*a[n][1]/l);
    d[2] = complex<double>(0.0,kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 1: {
    n /= 6;
    double l = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    d[0] = complex<double>(0.0,-kp*a[n][0]/l);
    d[1] = complex<double>(0.0,-kp*a[n][1]/l);
    d[2] = complex<double>(0.0,-kp*a[n][2]/l);
    d[3] = a[n][0]/l;
    d[4] = a[n][1]/l;
    d[5] = a[n][2]/l;
    break;
  }
  case 2: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 3: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = b[n][0]/lb;
    d[4] = b[n][1]/lb;
    d[5] = b[n][2]/lb;
    break;
  }
  case 4: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,ks*a[n][0]/la);
    d[1] = complex<double>(0.0,ks*a[n][1]/la);
    d[2] = complex<double>(0.0,ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
  }
  case 5: {
    n /= 6;
    double la = sqrt(a[n][0]*a[n][0]+a[n][1]*a[n][1]+a[n][2]*a[n][2]);
    double lb = sqrt(b[n][0]*b[n][0]+b[n][1]*b[n][1]+b[n][2]*b[n][2]);
    d[0] = complex<double>(0.0,-ks*a[n][0]/la);
    d[1] = complex<double>(0.0,-ks*a[n][1]/la);
    d[2] = complex<double>(0.0,-ks*a[n][2]/la);
    d[3] = (a[n][1]*b[n][2]-a[n][2]*b[n][1])/(la*lb);
    d[4] = (a[n][2]*b[n][0]-a[n][0]*b[n][2])/(la*lb);
    d[5] = (a[n][0]*b[n][1]-a[n][1]*b[n][0])/(la*lb);
    break;
   }
 }
}


class LEPELFunction3d : public IntegFunctionA3d {
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
 LEPELFunction3d(int _o, double *_materialc,
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

   complex<double> (*u)[3] = new complex<double>[ndir][3];
   complex<double> (*la)[3] = new complex<double>[nl][3];

   for (int k=0;k<ndir;k++) {
     complex<double> e = exp(dirs[k*3+0]*(x[0]-xc[0]) +
                             dirs[k*3+1]*(x[1]-xc[1]) +
                             dirs[k*3+2]*(x[2]-xc[2]));
     u[k][0] =  coef[k*3+0]*e;
     u[k][1] =  coef[k*3+1]*e;
     u[k][2] =  coef[k*3+2]*e;
   };

   int toc = 3*o*o*o;
   complex<double> (*v)[3] = new complex<double>[toc][3];

   for(int k=0;k<toc;k+=3) {
     v[k][0] = N[k/3]; v[k][1] = 0.0; v[k][2] = 0.0;
     v[k+1][1] = N[k/3]; v[k+1][0] = 0.0; v[k+1][2] = 0.0;
     v[k+2][0] = 0.0; v[k+2][1] = 0.0; v[k+2][2] = N[k/3];
   }

   for (int l=0;l<nl;l++) {
     complex<double> eu = exp(ldirs[l*3+0]*(x[0]-xsc[0])+
                              ldirs[l*3+1]*(x[1]-xsc[1])+
                              ldirs[l*3+2]*(x[2]-xsc[2]));
     la[l][0] = lcoef[l*3+0]*eu;
     la[l][1] = lcoef[l*3+1]*eu;
     la[l][2] = lcoef[l*3+2]*eu;
   };

   for (int l=0;l<nl;l++) {
     for (int k=0;k<ndir;k++) 
       kel[l*ndir+k] += w*omega*omega*
                    (u[k][0]*la[l][0]+u[k][1]*la[l][1]+u[k][2]*la[l][2])*
                sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
     for (int k=0;k<toc;k++) 
       kpl[l*toc+k] += w*omega*omega*
                      (v[k][0]*la[l][0]+v[k][1]*la[l][1]+v[k][2]*la[l][2])*
                sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
   }
    
   delete[] v;
   delete[] u;
   delete[] la;
 }
};


class LEPEFunction3d : public IntegFunctionV3d {
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
 LEPEFunction3d(int _o, double *_materialc,
                    int _ndir, complex<double> *_dirs, complex<double> *_coef, 
                    double *_xc, complex<double> *_kee,
                    complex<double> *_kpe, complex<double> *_kpp) {
   o = _o; materialc = _materialc; ndir = _ndir; dirs = _dirs; coef = _coef;
   xc = _xc; kee = _kee; kpe = _kpe; kpp = _kpp;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
   double &omega = materialc[1];
   double &lambda = materialc[2];
   double &mu = materialc[3];
   double &rho = materialc[4];
   double wdet = w*det;

   complex<double> (*u)[3] = new complex<double>[ndir][3];
   complex<double> (*gradu)[3][3] = new complex<double>[ndir][3][3];

   for(int k=0;k<ndir;k++) {
      complex <double> eu = exp(dirs[k*3+0]*(x[0]-xc[0])+
                                dirs[k*3+1]*(x[1]-xc[1])+
                                dirs[k*3+2]*(x[2]-xc[2]));
      u[k][0] = coef[k*3+0]*eu;
      u[k][1] = coef[k*3+1]*eu;
      u[k][2] = coef[k*3+2]*eu;
      gradu[k][0][0] =  coef[k*3+0]*dirs[k*3+0]*eu;
      gradu[k][0][1] =  coef[k*3+0]*dirs[k*3+1]*eu;
      gradu[k][0][2] =  coef[k*3+0]*dirs[k*3+2]*eu;
      gradu[k][1][0] =  coef[k*3+1]*dirs[k*3+0]*eu;
      gradu[k][1][1] =  coef[k*3+1]*dirs[k*3+1]*eu;
      gradu[k][1][2] =  coef[k*3+1]*dirs[k*3+2]*eu;
      gradu[k][2][0] =  coef[k*3+2]*dirs[k*3+0]*eu;
      gradu[k][2][1] =  coef[k*3+2]*dirs[k*3+1]*eu;
      gradu[k][2][2] =  coef[k*3+2]*dirs[k*3+2]*eu;
   }

   int toc = 3*o*o*o;
   complex<double> (*v)[3] = new complex<double>[toc][3];
   complex<double> (*gradv)[3][3] = new complex<double>[toc][3][3];

   for(int k=0;k<toc;k+=3) {
     v[k][0] = N[k/3]; v[k][1] = 0.0; v[k][2] = 0.0;
     v[k+1][1] = N[k/3]; v[k+1][0] = 0.0; v[k+2][0] = 0.0;
     v[k+2][0] = 0.0; v[k+2][1] = 0.0; v[k+2][2] = N[k/3];
     gradv[k][0][0] = dNdx[k/3][0];
     gradv[k][0][1] = dNdx[k/3][1];
     gradv[k][0][2] = dNdx[k/3][2];
     gradv[k][1][0] = gradv[k][1][1] = gradv[k][1][2] =  0.0;
     gradv[k][2][0] = gradv[k][2][1] = gradv[k][2][2] =  0.0;

     gradv[k+1][1][0] = dNdx[k/3][0];
     gradv[k+1][1][1] = dNdx[k/3][1];
     gradv[k+1][1][2] = dNdx[k/3][2];
     gradv[k+1][0][0] = gradv[k+1][0][1] = gradv[k+1][0][2] = 0.0;
     gradv[k+1][2][0] = gradv[k+1][2][1] = gradv[k+1][2][2] = 0.0;

     gradv[k+2][2][0] = dNdx[k/3][0];
     gradv[k+2][2][1] = dNdx[k/3][1];
     gradv[k+2][2][2] = dNdx[k/3][2];
     gradv[k+2][0][0] = gradv[k+2][0][1] = gradv[k+2][0][2] = 0.0;
     gradv[k+2][1][0] = gradv[k+2][1][1] = gradv[k+2][1][2] = 0.0;
   }

   for(int l=0;l<ndir;l++) {
     for(int k=0;k<ndir;k++)
       kee[l*ndir+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradu[k], gradu[l], u[k], u[l]);
     for(int k=0;k<toc;k++) 
       kpe[l*toc+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradv[k], gradu[l], v[k], u[l]);
   }
   for(int l=0;l<toc;l++) 
     for(int k=0;k<toc;k++) 
       kpp[l*toc+k] += wdet* omega*omega*
           elasticForm(lambda, mu, omega, rho, gradv[k], gradv[l], v[k], v[l]);
   delete[] v;
   delete[] gradv;
   delete[] u;
   delete[] gradu;
 }
};


void LEDEMMatrices3d(int order, double *xy, double *materialc,
                    int ndir, complex<double> *dirs, complex<double> *coef,
                    int *nldirs, complex<double> *ldirs, complex<double> *lcoef,
                    double *xsc, double *xc,
                    complex<double> *kee, complex<double> *kel,
                    complex<double> *kpp, complex<double> *kpl,
                    complex<double> *kpe) {

 IsoParamUtils ipu(order);
 int oc = ipu.getorderc();
 int toc = 3*oc;

 int nldir =nldirs[0] + nldirs[1] + nldirs[2] +
            nldirs[3] + nldirs[4] + nldirs[5];

 for(int i=0;i<nldir*ndir;i++) kel[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) kee[i] = 0.0;
 for(int i=0;i<toc*ndir;i++) kpe[i] = 0.0;
 for(int i=0;i<toc*nldir;i++) kpl[i] = 0.0;
 for(int i=0;i<toc*toc;i++) kpp[i] = 0.0;

 int c = 0;
 for(int faceindex=1;faceindex<=6;faceindex++) if (nldirs[faceindex-1]>0) {
   LEPELFunction3d g(order,materialc,ndir,dirs,coef,
                nldirs[faceindex-1],ldirs+c*3,lcoef+c*3,
                xsc + (faceindex-1)*3, xc,
                kel+c*ndir,kpl+c*toc);
   ipu.surfInt3d(xy, faceindex, g);
   c += nldirs[faceindex-1];
 }

 LEPEFunction3d f(order,materialc,ndir,dirs, coef, xc, kee, kpe,kpp);
 ipu.volumeInt3d(xy, f);
}



void DEMLE3d::createM(complex<double>*M) {

 IsoParamUtils ipu(o);
 int os = ipu.getordersq();
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 getNodalCoord(oc,nn,xyz);

 double xref[2];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 complex<double>* cndir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) {
   complex<double> d[6];
   dir(i,d);
   cdir[3*i+0] = d[0];
   cdir[3*i+1] = d[1];
   cdir[3*i+2] = d[2];
   cndir[3*i+0] = d[3];
   cndir[3*i+1] = d[4];
   cndir[3*i+2] = d[5];
 }

 int nldir[6] = {0,0,0,0,0,0};
 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<6;i++) tnldir += nldir[i];

 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) {
   DGMLE3d_LM *l = dynamic_cast<DGMLE3d_LM*>(lm[fi]);
   if (l==0) {
     fprintf(stderr,"DGMLE3d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMLE3d_LM type specific
 double xlref[18] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<6;fi++)
// ipu.sidecenter(xyz,fi+1,xlref+fi*3);

 complex<double>* cldir = new complex<double>[3*tnldir];
 complex<double>* clndir = new complex<double>[3*tnldir];

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
     DGMLE3d_LM *l = dynamic_cast<DGMLE3d_LM*>(lm[i]);
     if (l!=0) {
       for(int j=0;j<lm[i]->nDofs();j++) {
         complex<double> d[6];
         l->ldir(j,tau1,tau2,d);
         cldir[2*cc+0] = d[0];
         cldir[2*cc+1] = d[1];
         cldir[2*cc+2] = d[2];
         clndir[3*cc+0] = d[3];
         clndir[3*cc+1] = d[4];
         clndir[3*cc+2] = d[5];
         cc++;
       }
     }
   }
 }

 double materialc[5];
 materialc[1] = getOmega();
 double E = getE();
 double nu = getNu();
 materialc[4] = getRho();
 materialc[2] = nu*E/((1.0+nu)*(1.0-2.0*nu));
 materialc[3] = E/(2.0*(1.0+nu));

 int toc = 3*oc;

 complex<double>* kpp = new complex<double>[toc*toc];
 complex<double>* kpe = new complex<double>[ndir*toc];
 complex<double>* kpl = new complex<double>[toc*tnldir];

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];


 LEDEMMatrices3d(o, xyz, materialc,
                    ndir, cdir, cndir, nldir, cldir, clndir,
                    xlref, xref, kee, kel,kpp,kpl,kpe);

 cc = 0;
 for(int i=0;i<6;i++) {
  for(int j=cc;j<cc+nldir[i];j++) {
    for(int k=0;k<ndir;k++) kel[j*ndir+k] *= sign[i];
    for(int k=0;k<toc;k++) kpl[j*toc+k] *= sign[i];
  }
  cc += nldir[i];
 }

 for(int i=0;i<toc;i++) for(int j=0;j<toc;j++)
   M[j*(toc+tnldir+ndir)+i] = kpp[j*toc+i];
 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(toc+tnldir+j)*(toc+ndir+tnldir)+(toc+tnldir+i)] = kee[j*ndir+i];
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(toc+tnldir+i)*(toc+ndir+tnldir)+toc+j] = kel[j*ndir+i];
   M[(toc+j)*(toc+ndir+tnldir)+(toc+tnldir+i)] = kel[j*ndir+i];
 }
 for(int i=0;i<toc;i++) for(int j=0;j<ndir;j++) {
   M[i*(toc+ndir+tnldir)+(toc+tnldir+j)] = kpe[j*toc+i];
   M[(toc+tnldir+j)*(toc+ndir+tnldir)+i] = kpe[j*toc+i];
 }
 for(int i=0;i<toc;i++) for(int j=0;j<tnldir;j++) {
   M[i*(toc+ndir+tnldir)+(toc+j)] = kpl[j*toc+i];
   M[(toc+j)*(toc+ndir+tnldir)+i] = kpl[j*toc+i];
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++)
   M[(toc+j)*(toc+ndir+tnldir)+toc+i] = 0.0;

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



void DEMLE3d::createRHS(complex<double>*v) { 
}
