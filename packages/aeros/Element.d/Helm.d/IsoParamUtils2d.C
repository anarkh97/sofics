#include <alloca.h>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>

#include <Element.d/Helm.d/IsoParamUtils2d.h>

using namespace std;

template<> void IsoParamUtils2d::copyOut<complex<double> >(int m,int n,complex<double> *K,
                       double *Kr, double *Ki) {
 int kk;
 for(kk=0;kk<n;kk++) {
   int l;
   for(l=0;l<m;l++) {
     Kr[kk*m+l] = real(K[kk*m+l]);
     Ki[kk*m+l] = imag(K[kk*m+l]);
   }
 }
}


template<> void IsoParamUtils2d::copyOut<double>(int m,int n, double *K, double *Kr, double *Ki)
{
 int kk;
 for(kk=0;kk<n;kk++) {
   int l;
   for(l=0;l<m;l++) {
     Kr[kk*m+l] = (K[kk*m+l]);
   }
 }
}


complex<double> IsoParamUtils2d::expdiff(complex<double> alpha,
            complex<double> eaa, complex<double> eab, double L) {

 complex<double> aL = alpha*L;

 if (abs(aL)<1e-2) {
   complex<double> s(0.0,0.0);
   s = (1.0 + aL/2.0 * (1.0 + aL/3.0 * (1.0 + aL/4.0 * (1.0 + aL/5.0*
          (1.0 + aL/6.0)))));
   return L*s*eaa;
 } else {
   return (eab-eaa)/alpha;
 }
}


double IsoParamUtils2d::sidelen(double x, double y) {
 return sqrt(x*x+y*y);
}

// Identical to 3D
int IsoParamUtils2d::getordersq() {
 return order*order;
}

// Identical to 3D
int IsoParamUtils2dTri::getordersq() {
 return (order*(order+1))/2;
}

// Identical to 3D
void IsoParamUtils2dTri::lagGalShapeFunction(int mo, int ng, double *xi,
                              double *N, int secondDerivsFlag) {
 int i,j;
 long double *x = (long double*)alloca(sizeof(long double)*order);
 long double *xx = (long double*)alloca(sizeof(long double)*order);
 long double *xix = (long double*)alloca(sizeof(long double)*order);

 for(i=0;i<order;i++) {
   x[i] = (long double)(i)/(order-1);
 }

 for(i=0;i<mo;i++)
   xx[i] = x[mo-1]-x[i];

 for(j=0;j<ng;j++) {

   if (mo==1) {
      N[j] = 1.0;
      N[j+ng] = 0.0;
      if (secondDerivsFlag)  N[j+2*ng] = 0.0;
      continue;
   }

   int xizero = -1;
   for(i=0;i<mo;i++) {
     xix[i] = xi[j]-x[i]; 
     if (fabsl(xix[i])<1e-6) xizero = i;
   }

   long double tmpx = 1.0;
   int m;
   for(m=0;m<mo-1;m++)  tmpx *= xix[m]/xx[m];

   long double tmpxx = 0.0;
   if ((xizero==-1) || (xizero==mo-1)) {
     for(m=0;m<mo-1;m++)  tmpxx += tmpx/xix[m];
   } else {
       tmpxx = 1.0;
       for(m=0;m<mo-1;m++) if (m!=xizero) tmpxx *= xix[m]; 
       for(m=0;m<mo-1;m++) tmpxx /= xx[m]; 
   } 
   N[j] = tmpx;
   N[j+ng] = tmpxx;

   if (secondDerivsFlag) {
     long double tmpxxx = 0.0;
     int q;
     for(q=0;q<mo-1;q++)  {
       int p;
       for(p=0;p<mo-1;p++) if ((p!=q)) {
         long double t = 1.0;
         long double s = 1.0;
         for(m=0;m<mo-1;m++) if ((m!=p) && (m!=q)) t *= xix[m]; 
         for(m=0;m<mo-1;m++) s *= xx[m]; 
         tmpxxx += (t/s);
       }
     }
     N[j+2*ng] = tmpxxx;
   }
 }
}


// Identical to 3D
void IsoParamUtils2d::spectralLagGalShapeFunction(int ng, double *xi, double *N,
                         int secondDerivsFlag) {
 int i,j;
 long double *x = (long double*)alloca(sizeof(long double)*order);
 long double *xx = (long double*)alloca(sizeof(long double)*order*order);
 long double *xix = (long double*)alloca(sizeof(long double)*order);

 for(i=0;i<order;i++) {
   x[i] = (long double)xi[i];
 }

 for(i=0;i<order;i++)
   for(j=0;j<order;j++)
     xx[i*order+j] = x[i]-x[j];

 for(j=0;j<ng;j++) {
   int xizero = -1;
   for(i=0;i<order;i++) {
     xix[i] = xi[j]-x[i];
     if (fabsl(xix[i])<1e-6) xizero = i;
   }

   for(i=0;i<order;i++) {
     long double tmpx = 1.0;
     int m;
     for(m=0;m<order;m++) if (m!=i) tmpx *= xix[m]/xx[i*order+m];

     long double tmpxx = 0.0;
     if ((xizero==-1) || (xizero==i)) {
       for(m=0;m<order;m++)  if (m!=i) tmpxx += tmpx/xix[m];
     } else {
       tmpxx = 1.0;
       for(m=0;m<order;m++) if ((m!=xizero) && (m!=i)) tmpxx *= xix[m];
       for(m=0;m<order;m++) if (m!=i) tmpxx /= xx[i*order+m];
     }

     N[j*order+i] = tmpx;
     N[ng*order+j*order+i] = tmpxx;

     if (secondDerivsFlag) {
       long double tmpxxx = 0.0;
       int q;
       for(q=0;q<order;q++)  if (q!=i) {
         int p;
         for(p=0;p<order;p++) if ((p!=i) && (p!=q)) {
           long double t = 1.0;
           long double s = 1.0;
           for(m=0;m<order;m++) if ((m!=i) && (m!=p) && (m!=q)) t *= xix[m];
           for(m=0;m<order;m++) if (m!=i) s *= xx[i*order+m];
           tmpxxx += (t/s);
         }
       }
       N[2*ng*order+j*order+i] = tmpxxx;
     }
   }
 }
}


// Identical to 3D
void IsoParamUtils2d::lagGalShapeFunction(int ng, double *xi, double *N, 
                         int secondDerivsFlag) {
 int i,j;
 long double *x = (long double*)alloca(sizeof(long double)*order);
 long double *xx = (long double*)alloca(sizeof(long double)*order*order);
 long double *xix = (long double*)alloca(sizeof(long double)*order);

 for(i=0;i<order;i++) {
   x[i] = -1.0+2.0*(long double)(i)/(order-1);
 }

 for(i=0;i<order;i++)
   for(j=0;j<order;j++)
     xx[i*order+j] = x[i]-x[j];

 for(j=0;j<ng;j++) {  
   int xizero = -1;
   for(i=0;i<order;i++) {
     xix[i] = xi[j]-x[i]; 
     if (fabsl(xix[i])<1e-6) xizero = i;
   }
  
   for(i=0;i<order;i++) { 
     long double tmpx = 1.0;
     int m;
     for(m=0;m<order;m++) if (m!=i) tmpx *= xix[m]/xx[i*order+m];
  
     long double tmpxx = 0.0;
     if ((xizero==-1) || (xizero==i)) {
       for(m=0;m<order;m++)  if (m!=i) tmpxx += tmpx/xix[m];
     } else {
       tmpxx = 1.0;
       for(m=0;m<order;m++) if ((m!=xizero) && (m!=i)) tmpxx *= xix[m]; 
       for(m=0;m<order;m++) if (m!=i) tmpxx /= xx[i*order+m]; 
     }
 
     N[j*order+i] = tmpx;
     N[ng*order+j*order+i] = tmpxx;

     if (secondDerivsFlag) {
       long double tmpxxx = 0.0;
       int q;
       for(q=0;q<order;q++)  if (q!=i) {
         int p;
         for(p=0;p<order;p++) if ((p!=i) && (p!=q)) {
           long double t = 1.0;
           long double s = 1.0;
           for(m=0;m<order;m++) if ((m!=i) && (m!=p) && (m!=q)) t *= xix[m]; 
           for(m=0;m<order;m++) if (m!=i) s *= xx[i*order+m]; 
           tmpxxx += (t/s);
         }
       }
       N[2*ng*order+j*order+i] = tmpxxx;
     }
   }
 }
}


void IsoParamUtils2d::jmatrix(double *N, double *xy, double (*J)[2]) {

 J[0][0] = 0.0;
 J[0][1] = 0.0;
 J[1][0] = 0.0;
 J[1][1] = 0.0;

 int ordersq = getordersq();
 int m;
 for(m=0;m<ordersq;m++) {
   J[0][0] += N[ordersq+m]*xy[m];
   J[0][1] += N[ordersq+m]*xy[m+ordersq];
   J[1][0] += N[2*ordersq+m]*xy[m];
   J[1][1] += N[2*ordersq+m]*xy[m+ordersq];
 }
}


void IsoParamUtils2dTri::jmatrix(double *N, double *xy, double (*J)[2]) {

 J[0][0] = 0.0;
 J[0][1] = 0.0;
 J[1][0] = 0.0;
 J[1][1] = 0.0;

 int ordersq = getordersq();
 int m;
 for(m=0;m<ordersq;m++) {
   J[0][0] += N[ordersq+m]*xy[m];
   J[0][1] += N[ordersq+m]*xy[m+ordersq];
   J[1][0] += N[2*ordersq+m]*xy[m];
   J[1][1] += N[2*ordersq+m]*xy[m+ordersq];
 }
}


double IsoParamUtils2d::detj(double (*J)[2]) {

 return  J[0][0]*J[1][1]-J[0][1]*J[1][0];
}


void IsoParamUtils2d::invj(double (*J)[2], double det, double (*j)[2]) {

 j[0][0] = J[1][1]/det;
 j[1][0] = -J[1][0]/det;
 j[0][1] = -J[0][1]/det;
 j[1][1] = J[0][0]/det;
}


void IsoParamUtils2d::crossj(double (*J)[2], int axis, double *cross) {

 if (axis==1) {
   cross[0] = J[1][1];
   cross[1] = -J[1][0];
 } else {
   cross[0] = J[0][1];
   cross[1] = -J[0][0];
 }
}


void IsoParamUtils2dTri::crossj(double (*J)[2], int faceindex, double *cross) {

 if (faceindex==1) {
   cross[0] = (J[1][1]-J[0][1]);
   cross[1] = -(J[1][0]-J[0][0]);
 } else if (faceindex==2) {
   cross[0] = J[1][1];
   cross[1] = -J[1][0];
 } else {
   cross[0] = J[0][1];
   cross[1] = -J[0][0];
 }
}


void IsoParamUtils2d::elementcenter(double *xy, double *cxy) {
 int ordersq = order*order;   
 int corner[4] = { 1-1, order-1, ordersq-order+1-1, ordersq-1};
 int i;
 cxy[0] = cxy[1] = 0.0;
 for(i=0;i<4;i++) {
   cxy[0] += xy[0*ordersq+corner[i]];
   cxy[1] += xy[1*ordersq+corner[i]];
 }
 cxy[0] /= 4.0;
 cxy[1] /= 4.0; 
}


void IsoParamUtils2dTri::elementcenter(double *xy, double *cxy) {
 int ordersq = getordersq();
 int corner[3] = { ordersq-1, order-1, 1-1 };
 int i;
 cxy[0] = cxy[1] = 0.0;
 for(i=0;i<3;i++) {
   cxy[0] += xy[0*ordersq+corner[i]];
   cxy[1] += xy[1*ordersq+corner[i]];
 }
 cxy[0] /= 3.0;
 cxy[1] /= 3.0; 
}


void IsoParamUtils2d::sidecenter(double *xy, int faceindex, double *scxy) {

 int ordersq = getordersq();
 int corner[5] = { 1-1, order-1, ordersq-1, ordersq-order+1-1, 1-1};
  
 scxy[0] = 0.5*(xy[0*ordersq+corner[faceindex-1]]+
                 xy[0*ordersq+corner[faceindex]]);
 scxy[1] = 0.5*(xy[1*ordersq+corner[faceindex-1]]+
                 xy[1*ordersq+corner[faceindex]]);
}


void IsoParamUtils2dTri::sidecenter(double *xy, int faceindex, double *scxy) {

 int ordersq = getordersq();
 int corner[4] = { order-1, ordersq-1, 1-1, order-1 };
  
 scxy[0] = 0.5*(xy[0*ordersq+corner[faceindex-1]]+
                 xy[0*ordersq+corner[faceindex]]);
 scxy[1] = 0.5*(xy[1*ordersq+corner[faceindex-1]]+
                 xy[1*ordersq+corner[faceindex]]);
}


void IsoParamUtils2d::faceindeces(int faceindex, int *fi) {
 int ordersq = getordersq();
 int k;
 if (faceindex==1) {
   for(k=0;k<order;k++) 
     fi[k] = k;
 }
 else if (faceindex==2) {
   for(k=0;k<order;k++)
     fi[k] = order-1+k*order;
 }
 else if (faceindex==3) {
   for(k=0;k<order;k++)
     fi[k] = ordersq-1-k;
 }
 else {
   for(k=0;k<order;k++) 
     fi[k] = order*(order-k-1);
 }
}


void IsoParamUtils2dTri::faceindeces(int faceindex, int *fi) {
 //int ordersq = getordersq();
 int k;
 if (faceindex==1) {
   fi[0] = order-1;
   for(k=1;k<order;k++) 
     fi[k] = fi[k-1] + (order-k) ;
 }
 else if (faceindex==3) {
   for(k=0;k<order;k++)
     fi[k] = k;
 }
 else {
   fi[0] = 0;
   for(k=1;k<order;k++) 
     fi[k] = fi[k-1] + (order-k+1);
 }
}


void IsoParamUtils2d::facemap(int &faceindex, int *fn, int *n, int *map) {

 int ordersq = getordersq();
 int c[5] = { 1-1, order-1, ordersq-1, ordersq-order+1-1, 1-1};
 if ( (fn[1-1]==n[c[0]] && fn[order-1]==n[c[1]]) ||
      (fn[1-1]==n[c[1]] && fn[order-1]==n[c[0]]) ) faceindex = 1;
 else 
 if ( (fn[1-1]==n[c[1]] && fn[order-1]==n[c[2]]) ||
      (fn[1-1]==n[c[2]] && fn[order-1]==n[c[1]]) ) faceindex = 2;
 else 
 if ( (fn[1-1]==n[c[2]] && fn[order-1]==n[c[3]]) ||
      (fn[1-1]==n[c[3]] && fn[order-1]==n[c[2]]) ) faceindex = 3;
 else 
 if ( (fn[1-1]==n[c[3]] && fn[order-1]==n[c[4]]) ||
      (fn[1-1]==n[c[4]] && fn[order-1]==n[c[3]]) ) faceindex = 4;
 else {
   fprintf(stderr,"Impossible! Exiting.\n"); 
   exit(-1);
 }

 int *fi = (int*)alloca(sizeof(int)*order);
 faceindeces(faceindex,fi);

 int i;
 for(i=0;i<order;i++) {
   int j;
   for(j=0;j<order;j++) {
     if (n[fi[i]]==fn[j]) map[j] = fi[i];
   } 
 }
}


void IsoParamUtils2dTri::facemap(int &faceindex, int *fn, int *n, int *map) {

 int ordersq = getordersq();
 int c[4] = { order-1, ordersq-1, 1-1, order-1 };
 if ( (fn[1-1]==n[c[0]] && fn[order-1]==n[c[1]]) ||
      (fn[1-1]==n[c[1]] && fn[order-1]==n[c[0]]) ) faceindex = 1;
 else 
 if ( (fn[1-1]==n[c[1]] && fn[order-1]==n[c[2]]) ||
      (fn[1-1]==n[c[2]] && fn[order-1]==n[c[1]]) ) faceindex = 2;
 else 
 if ( (fn[1-1]==n[c[2]] && fn[order-1]==n[c[3]]) ||
      (fn[1-1]==n[c[3]] && fn[order-1]==n[c[2]]) ) faceindex = 3;
 else {
   fprintf(stderr,"Impossible! Exiting.\n"); 
   exit(-1);
 }

 int *fi = (int*)alloca(sizeof(int)*order);
 faceindeces(faceindex,fi);

 int i;
 for(i=0;i<order;i++) {
   int j;
   for(j=0;j<order;j++) {
     if (n[fi[i]]==fn[j]) map[j] = fi[i];
   } 
 }
}


int IsoParamUtils2d::isStraight(double *xy, int faceindex) {

 int ordersq = getordersq();
 int &o = order;
 int edges[4][2] = {
        { 0,1 },
        { o-1,o },
        { o*o-1,-1 },
        { o*o-o,-o }} ;

 int fi = faceindex - 1;
 int n1 = edges[fi][0];
 int step = edges[fi][1];
 int n2 = edges[fi][0]+(o-1)*step;

 double xy1[2] = {xy[0*ordersq+n1],xy[1*ordersq+n1]};
 double xy2[2] = {xy[0*ordersq+n2],xy[1*ordersq+n2]};
 double dx[2] = {
   xy2[0]-xy1[0],
   xy2[1]-xy1[1]
 };

 int straight=1;
 double normdx = sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
 for(int j=1;j<o-1;j++) {
   int nn = n1+j*step;
   double xyn[2] = {xy[0*ordersq+nn],xy[1*ordersq+nn]};
   double dxn[2] = {
     xyn[0]-xy1[0],
     xyn[1]-xy1[1]
   };
   
   double normdxn = sqrt(dxn[0]*dxn[0]+dxn[1]*dxn[1]);
   if (fabs(dxn[0]*dx[1]-dxn[1]*dx[0])/(normdx*normdxn)>1e-4) {
     straight = 0;
     break;
   }
 }
 return straight;
}


int IsoParamUtils2dTri::isStraight(double *xy, int faceindex) {

 int ordersq = getordersq();
 int &o = order;
 int *fi = (int*)alloca(sizeof(int)*order);
 faceindeces(faceindex,fi);

 int n1 = fi[0];
 int n2 = fi[o-1];
 double xy1[2] = {xy[0*ordersq+n1],xy[1*ordersq+n1]};
 double xy2[2] = {xy[0*ordersq+n2],xy[1*ordersq+n2]};
 double dx[2] = {
   xy2[0]-xy1[0],
   xy2[1]-xy1[1]
 };

 int straight=1;
 double normdx = sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
 for(int j=1;j<o-1;j++) {
   int nn = fi[j];
   double xyn[2] = {xy[0*ordersq+nn],xy[1*ordersq+nn]};
   double dxn[2] = { xyn[0]-xy1[0], xyn[1]-xy1[1] };
   double normdxn = sqrt(dxn[0]*dxn[0]+dxn[1]*dxn[1]);
   if ((dxn[0]*dx[1]-dxn[1]*dx[0])/(normdx*normdxn)>1e-4) {
     straight = 0;
     break;
   }
 }
 return straight;
}


void IsoParamUtils2d::lineInt2d(double *xy, int faceindex, IntegFunctionL2d &f,
                                int gorder) {

 int ordersq = order*order;
 int i;

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);
 double p;
 double *sNN = (double*)alloca(sizeof(double)*2*order);

 int axis = 0;
 if (faceindex==1) {
   p = -1.0;
   axis = 2;
 } else if (faceindex==2) {
   p = 1.0;
   axis = 1;
 } else if (faceindex==3) {
   p = 1.0;
   axis = 2;
 } else {
   p = -1.0;
   axis = 1;
 }
 lagGalShapeFunction(1,&p,sNN);

 double x[6];
 double *cxy = x+2;
 elementcenter(xy,cxy);
 double *scxy = x+4;
 sidecenter(xy,faceindex,scxy);


 double *N = (double*) alloca(sizeof(double)*3*ordersq);
 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {
   int ii,jj;
   int ng1,ng2;
   double *NN1, *NN2;
   NN1 = NN2 = NN;
   if (axis==1) {
     ii = 0; jj = i;
     ng1 = 1; ng2 = ng;
     NN1 = sNN;
   } else {
     ii = i; jj = 0;
     ng1 = ng; ng2 = 1;
     NN2 = sNN;
   }

   int mx,my;
   for(mx=0;mx<order;mx++) {
     for(my=0;my<order;my++) {
        N[my*order+mx] =
          NN1[ii*order+mx]*NN2[jj*order+my];
        N[ordersq+my*order+mx] =
          NN1[ng1*order+ii*order+mx]*NN2[jj*order+my];
        N[2*ordersq+my*order+mx] =
          NN1[ii*order+mx]*NN2[ng2*order+jj*order+my];
      }
    }
         
   double w = gr.wgauss[i];

   double J[2][2];
   jmatrix(N,xy,J);

   double cross[2];
   crossj(J,axis,cross);
  
   double nsign;
   if (cross[0]*(scxy[0]-cxy[0])+
       cross[1]*(scxy[1]-cxy[1]) > 0.0) nsign = 1;
   else nsign = -1;

   x[0] = 0.0;
   x[1] = 0.0;
   int m;
   for(m=0;m<ordersq;m++) {
     x[0] += N[m]*xy[m];
     x[1] += N[m]*xy[m+ordersq];
   }
 
   f.evaluate(x,N,cross,nsign,w);
 }
}


void IsoParamUtils2d::lineLineInt2d(double *xy, IntegFunctionL2d &f, int gorder) {

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double *N = (double*) alloca(sizeof(double)*2*order);
 int ng = gr.ngauss;
 int i;
 for(i=0;i<ng;i++) {

   int m;
   for(m=0;m<order;m++) {
     N[m] = NN[i*order+m];
     N[order+m] = NN[ng*order+i*order+m];
   }
  
   double w = gr.wgauss[i];

   double J[1][2] = {{0.0,0.0}};
   for(m=0;m<order;m++) {
     J[0][0] += N[order+m]*xy[m];
     J[0][1] += N[order+m]*xy[m+order];
   }

   double cross[2] = { J[0][1], - J[0][0] };
   double nsign = 1.0;

   double x[2] = { 0.0, 0.0};
   for(m=0;m<order;m++) {
     x[0] += N[m]*xy[m];
     x[1] += N[m]*xy[m+order];
   }
 
   f.evaluate(x,N,cross,nsign,w);
 }
}


void IsoParamUtils2dTri::lineInt2d(double *xyz, int faceindex,
                                   IntegFunctionL2d &f, int gorder) {

 int ordersq = getordersq();

 double cxyz[2];
 elementcenter(xyz,cxyz);
 double scxyz[2];
 sidecenter(xyz,faceindex,scxyz);

 GaussRuleTriLine gr(gorder,GSUBDIV);
 gr.pool = (double*)alloca(sizeof(double)*4*gr.ngauss);
 gr.init(faceindex);

 int i;
 double *NN = (double*)alloca(sizeof(double)*2*3*gr.ngauss*order);
 for(i=1;i<=order;i++)
   lagGalShapeFunction(i,3*gr.ngauss,gr.xigauss,NN+(i-1)*2*3*gr.ngauss);

 double *N = (double*) alloca(sizeof(double)*3*ordersq);
 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {

   int r,s;
   int c=0;
   for(s=0;s<order;s++) {
     for(r=0;r<order-s;r++) {
       int u = order - 1 - r - s;
       int rr = r*2*3*gr.ngauss+0*gr.ngauss+i;
       int ss = s*2*3*gr.ngauss+1*gr.ngauss+i;
       int uu = u*2*3*gr.ngauss+2*gr.ngauss+i;
       int rrd = rr+3*gr.ngauss;
       int ssd = ss+3*gr.ngauss;
       int uud = uu+3*gr.ngauss;
       N[c] = NN[rr]*NN[ss]*NN[uu];
       N[1*ordersq+c] = (NN[rrd]*NN[ss]*NN[uu]-
                       NN[rr]*NN[ss]*NN[uud]);
       N[2*ordersq+c] = (NN[rr]*NN[ssd]*NN[uu]-
                       NN[rr]*NN[ss]*NN[uud]);
       c++;
     }
   }

   double w = gr.wgauss[i];

   double J[2][2];
   jmatrix(N,xyz,J);

   double cross[2];
   crossj(J,faceindex,cross);
  
   double nsign;
   if (cross[0]*(scxyz[0]-cxyz[0])+
       cross[1]*(scxyz[1]-cxyz[1]) > 0.0) nsign = 1;
   else nsign = -1;

   double x[2] = {0.0,0.0};
   int m;
   for(m=0;m<ordersq;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+ordersq];
   }

   f.evaluate(x,N,cross,nsign,w); 
 }
}


void IsoParamUtils2d::areaInt2d(double *xy, IntegFunctionA2d &f, int gorder) {

 int ordersq = getordersq();

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double x[4];
 double *cxy = x+2;
 elementcenter(xy,cxy);

 double *N = (double*) alloca(sizeof(double)*3*ordersq);
 double (*dNdx)[2] = (double(*)[2])alloca(sizeof(double)*ordersq*2);
 int ng = gr.ngauss;
 int i;
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {

     int mx,my;
     for(mx=0;mx<order;mx++) {
       for(my=0;my<order;my++) {
          N[my*order+mx] =
                NN[i*order+mx]*NN[j*order+my];
          N[ordersq+my*order+mx] =
                NN[ng*order+i*order+mx]*NN[j*order+my];
          N[2*ordersq+my*order+mx] =
                NN[i*order+mx]*NN[ng*order+j*order+my];
        }
      }
           
     double w = gr.wgauss[i]*gr.wgauss[j];

     double J[2][2];
     jmatrix(N,xy,J);
     double det = detj(J); 
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f\n",
                          det,xy[0],xy[ordersq+0]);

     double j[2][2];
     invj(J,det,j);
     int m;
     for(m=0;m<ordersq;m++) {
       dNdx[m][0] = 
        j[0][0]*N[ordersq+m]+j[0][1]*N[2*ordersq+m];
       dNdx[m][1] = 
        j[1][0]*N[ordersq+m]+j[1][1]*N[2*ordersq+m];
     }

     x[0] = 0.0;
     x[1] = 0.0;
     for(m=0;m<ordersq;m++) {
       x[0] += N[m]*xy[m];
       x[1] += N[m]*xy[m+ordersq];
     }
     f.evaluate(x,N,dNdx,w,det);
   }
 }
}


void IsoParamUtils2d::spectralAreaInt2d(double *xy, IntegFunctionA2d &f) {

 int ordersq = getordersq();

 GaussLobattoRuleLine gr(order);

 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*gr.ngauss);
 spectralLagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double x[4];
 double *cxy = x+2;
 elementcenter(xy,cxy);

 double *N = (double*) alloca(sizeof(double)*3*ordersq);
 double (*dNdx)[2] = (double(*)[2])alloca(sizeof(double)*ordersq*2);
 int ng = gr.ngauss;
 int i;
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {

     int mx,my;
     for(mx=0;mx<order;mx++) {
       for(my=0;my<order;my++) {
          N[my*order+mx] =
                NN[i*order+mx]*NN[j*order+my];
          N[ordersq+my*order+mx] =
                NN[ng*order+i*order+mx]*NN[j*order+my];
          N[2*ordersq+my*order+mx] =
                NN[i*order+mx]*NN[ng*order+j*order+my];
        }
      }
           
     double w = gr.wgauss[i]*gr.wgauss[j];

     double J[2][2];
     jmatrix(N,xy,J);
     double det = detj(J); 
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f\n",
                          det,xy[0],xy[ordersq+0]);

     double j[2][2];
     invj(J,det,j);
     int m;
     for(m=0;m<ordersq;m++) {
       dNdx[m][0] = 
        j[0][0]*N[ordersq+m]+j[0][1]*N[2*ordersq+m];
       dNdx[m][1] = 
        j[1][0]*N[ordersq+m]+j[1][1]*N[2*ordersq+m];
     }

     x[0] = 0.0;
     x[1] = 0.0;
     for(m=0;m<ordersq;m++) {
       x[0] += N[m]*xy[m];
       x[1] += N[m]*xy[m+ordersq];
     }
     f.evaluate(x,N,dNdx,w,det);
   }
 }
}


void IsoParamUtils2dTri::areaInt2d(double *xy, IntegFunctionA2d &f,
                                     int gorder) {

 int ordersq = getordersq();

 GaussRuleTriangle gr(gorder);
 gr.pool = (double*)alloca(sizeof(double)*4*gr.ngauss);
 gr.init(0);

 double *NN = (double*)alloca(sizeof(double)*2*3*gr.ngauss*order);
 int i;
 for(i=1;i<=order;i++)
   lagGalShapeFunction(i,3*gr.ngauss,gr.xigauss,NN+(i-1)*2*3*gr.ngauss);

 double *N = (double*) alloca(sizeof(double)*3*ordersq);
 double (*dNdx)[2] = (double(*)[2])alloca(sizeof(double)*ordersq*2);
 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {

   int r,s;
   int c=0;
   for(s=0;s<order;s++) {
     for(r=0;r<order-s;r++) {
       int u = order - 1 - r - s;
       int rr = r*2*3*gr.ngauss+0*gr.ngauss+i;
       int ss = s*2*3*gr.ngauss+1*gr.ngauss+i;
       int uu = u*2*3*gr.ngauss+2*gr.ngauss+i;
       int rrd = rr+3*gr.ngauss;
       int ssd = ss+3*gr.ngauss;
       int uud = uu+3*gr.ngauss;
       N[c] = NN[rr]*NN[ss]*NN[uu];
       N[1*ordersq+c] = (NN[rrd]*NN[ss]*NN[uu]-
                       NN[rr]*NN[ss]*NN[uud]);
       N[2*ordersq+c] = (NN[rr]*NN[ssd]*NN[uu]-
                       NN[rr]*NN[ss]*NN[uud]);
       c++;
     }
   }

   double J[2][2];
   jmatrix(N,xy,J);
   double det = detj(J); 
   if (det<0.0) fprintf(stderr,"Negative vol Jacobian: %f %f %f %f\n",
                        det,gr.xigauss[0*gr.ngauss+i],
                        gr.xigauss[1*gr.ngauss+i],
                        gr.xigauss[2*gr.ngauss+i]);

   double j[2][2];
   invj(J,det,j);
   int m;
   for(m=0;m<ordersq;m++) {
     dNdx[m][0] = 
      j[0][0]*N[ordersq+m]+j[0][1]*N[2*ordersq+m];
     dNdx[m][1] = 
      j[1][0]*N[ordersq+m]+j[1][1]*N[2*ordersq+m];
   }

   double x[2] = {0.0,0.0};
   for(m=0;m<ordersq;m++) {
     x[0] += N[m]*xy[m];
     x[1] += N[m]*xy[m+ordersq];
   }
   f.evaluate(x,N,dNdx,gr.wgauss[i],det);
 }
}

#ifndef _TEMPLATE_FIX_
#include <Element.d/Helm.d/IsoParamUtils2dT.C>
#endif
