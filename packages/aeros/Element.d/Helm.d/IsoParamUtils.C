#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>

#include <Element.d/Helm.d/IsoParamUtils.h>

using namespace std;

double IsoParamUtils::sidelen(double x, double y, double z) {
 return sqrt(x*x+y*y+z*z);
}


int IsoParamUtils::getorderc() {
 return order*order*order;
}


int IsoParamUtilsTetra::getorderc() {
 return (order*(order+1)*(order+2))/6;
}


int IsoParamUtilsPrism::getorderc() {
 return order*(order*(order+1))/2;
}


int IsoParamUtils::getordersq(int faceindex) {
 return order*order;
}


int IsoParamUtilsTetra::getordersq(int faceindex) {
 return (order*(order+1))/2;
}


int IsoParamUtilsPrism::getordersq(int faceindex) {
 if (faceindex==0) return 0;
 else if (faceindex<=2) return (order*(order+1))/2;
 else return order*order;
}

void IsoParamUtilsPrism::lagGalShapeFunctionTr(int mo, int ng, double *xi,
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



void IsoParamUtilsTetra::lagGalShapeFunction(int mo, int ng, double *xi,
                              double *N, int secondDerivsFlag) {
 int i,j;
 long double *x = (long double*)dbg_alloca(sizeof(long double)*order);
 long double *xx = (long double*)dbg_alloca(sizeof(long double)*order);
 long double *xix = (long double*)dbg_alloca(sizeof(long double)*order);

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


void IsoParamUtils::lagGalShapeFunction(int ng, double *xi, double *N, 
                         int secondDerivsFlag) {
 int i,j;
 long double *x = (long double*)dbg_alloca(sizeof(long double)*order);
 long double *xx = (long double*)dbg_alloca(sizeof(long double)*order*order);
 long double *xix = (long double*)dbg_alloca(sizeof(long double)*order);

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


void IsoParamUtils::spectralLagGalShapeFunction(int ng, double *xi, double *N, 
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


void IsoParamUtils::crossproduct(double *af, double *bf, double *nf) {
 nf[0] = af[1]*bf[2]-af[2]*bf[1];
 nf[1] = af[2]*bf[0]-af[0]*bf[2];
 nf[2] = af[0]*bf[1]-af[1]*bf[0];
}


void IsoParamUtils::jmatrix(double *N, double *xyz, double (*J)[3]) {

 J[0][0] = 0.0;
 J[0][1] = 0.0;
 J[0][2] = 0.0;
 J[1][0] = 0.0;
 J[1][1] = 0.0;
 J[1][2] = 0.0;
 J[2][0] = 0.0;
 J[2][1] = 0.0;
 J[2][2] = 0.0;

 int orderc = getorderc();
 for(int m=0;m<orderc;m++) {
   J[0][0] += N[orderc+m]*xyz[m];
   J[0][1] += N[orderc+m]*xyz[m+orderc];
   J[0][2] += N[orderc+m]*xyz[m+2*orderc];
   J[1][0] += N[2*orderc+m]*xyz[m];
   J[1][1] += N[2*orderc+m]*xyz[m+orderc];
   J[1][2] += N[2*orderc+m]*xyz[m+2*orderc];
   J[2][0] += N[3*orderc+m]*xyz[m];
   J[2][1] += N[3*orderc+m]*xyz[m+orderc];
   J[2][2] += N[3*orderc+m]*xyz[m+2*orderc];
 }
}

double IsoParamUtils::detj(double (*J)[3]) {

 return  J[0][0]*(J[1][1]*J[2][2]-J[2][1]*J[1][2])-
         J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+
         J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}


void IsoParamUtils::invj(double (*J)[3], double det, double (*j)[3]) {

 j[0][0] = (J[1][1]*J[2][2]-J[1][2]*J[2][1])/det;
 j[1][0] = (J[1][2]*J[2][0]-J[1][0]*J[2][2])/det;
 j[2][0] = (J[1][0]*J[2][1]-J[1][1]*J[2][0])/det;
 j[0][1] = (J[0][2]*J[2][1]-J[0][1]*J[2][2])/det;
 j[1][1] = (J[0][0]*J[2][2]-J[0][2]*J[2][0])/det;
 j[2][1] = (J[0][1]*J[2][0]-J[0][0]*J[2][1])/det;
 j[0][2] = (J[0][1]*J[1][2]-J[0][2]*J[1][1])/det;
 j[1][2] = (J[0][2]*J[1][0]-J[0][0]*J[1][2])/det;
 j[2][2] = (J[0][0]*J[1][1]-J[0][1]*J[1][0])/det;
}


void IsoParamUtils::crossj(double (*J)[3], int axis, double *cross) {

 if (axis==1) {
   cross[0] = J[1][1]*J[2][2]-J[1][2]*J[2][1]; 
   cross[1] = J[1][2]*J[2][0]-J[1][0]*J[2][2]; 
   cross[2] = J[1][0]*J[2][1]-J[1][1]*J[2][0]; 
 } else if (axis==2) {
   cross[0] = J[0][1]*J[2][2]-J[0][2]*J[2][1]; 
   cross[1] = J[0][2]*J[2][0]-J[0][0]*J[2][2]; 
   cross[2] = J[0][0]*J[2][1]-J[0][1]*J[2][0]; 
 } else {
   cross[0] = J[1][1]*J[0][2]-J[1][2]*J[0][1]; 
   cross[1] = J[1][2]*J[0][0]-J[1][0]*J[0][2]; 
   cross[2] = J[1][0]*J[0][1]-J[1][1]*J[0][0]; 
 }
}


void IsoParamUtils::crossj(double (*J)[3], int axis, double *cross,
                           double *tau1, double *tau2) {

 if (axis==1) {
   cross[0] = J[1][1]*J[2][2]-J[1][2]*J[2][1]; 
   cross[1] = J[1][2]*J[2][0]-J[1][0]*J[2][2]; 
   cross[2] = J[1][0]*J[2][1]-J[1][1]*J[2][0];
   tau1[0] = J[1][0]; tau1[1] = J[1][1]; tau1[2] = J[1][2]; 
   tau2[0] = J[2][0]; tau2[1] = J[2][1]; tau2[2] = J[2][2]; 
 } else if (axis==2) {
   cross[0] = J[0][1]*J[2][2]-J[0][2]*J[2][1]; 
   cross[1] = J[0][2]*J[2][0]-J[0][0]*J[2][2]; 
   cross[2] = J[0][0]*J[2][1]-J[0][1]*J[2][0]; 
   tau1[0] = J[0][0]; tau1[1] = J[0][1]; tau1[2] = J[0][2]; 
   tau2[0] = J[2][0]; tau2[1] = J[2][1]; tau2[2] = J[2][2]; 
 } else {
   cross[0] = J[1][1]*J[0][2]-J[1][2]*J[0][1]; 
   cross[1] = J[1][2]*J[0][0]-J[1][0]*J[0][2]; 
   cross[2] = J[1][0]*J[0][1]-J[1][1]*J[0][0]; 
   tau1[0] = J[1][0]; tau1[1] = J[1][1]; tau1[2] = J[1][2]; 
   tau2[0] = J[0][0]; tau2[1] = J[0][1]; tau2[2] = J[0][2]; 
 }
}


void IsoParamUtilsTetra::crossj(double (*J)[3], int faceindex, double *cross) {
 
 if (faceindex==1) {
   cross[0] = (J[0][1]-J[2][1])*(J[1][2]-J[2][2])-(J[0][2]-J[2][2])*(J[1][1]-J[2][1]); 
   cross[1] = (J[0][2]-J[2][2])*(J[1][0]-J[2][0])-(J[0][0]-J[2][0])*(J[1][2]-J[2][2]); 
   cross[2] = (J[0][0]-J[2][0])*(J[1][1]-J[2][1])-(J[0][1]-J[2][1])*(J[1][0]-J[2][0]); 
 } else if (faceindex==2) {
   cross[0] = -(J[1][1]*J[2][2]-J[1][2]*J[2][1]); 
   cross[1] = -(J[1][2]*J[2][0]-J[1][0]*J[2][2]); 
   cross[2] = -(J[1][0]*J[2][1]-J[1][1]*J[2][0]); 
 } else if (faceindex==3) {
   cross[0] = J[0][1]*J[2][2]-J[0][2]*J[2][1]; 
   cross[1] = J[0][2]*J[2][0]-J[0][0]*J[2][2]; 
   cross[2] = J[0][0]*J[2][1]-J[0][1]*J[2][0]; 
 } else {
   cross[0] = J[1][1]*J[0][2]-J[1][2]*J[0][1]; 
   cross[1] = J[1][2]*J[0][0]-J[1][0]*J[0][2]; 
   cross[2] = J[1][0]*J[0][1]-J[1][1]*J[0][0]; 
 }
}


void IsoParamUtilsPrism::crossj(double (*J)[3], int faceindex, double *cross) {
 
 if (faceindex==1) {
   cross[0] = J[1][1]*J[0][2]-J[1][2]*J[0][1]; 
   cross[1] = J[1][2]*J[0][0]-J[1][0]*J[0][2]; 
   cross[2] = J[1][0]*J[0][1]-J[1][1]*J[0][0]; 
 } else if (faceindex==2) {
   cross[0] = -J[1][1]*J[0][2]+J[1][2]*J[0][1]; 
   cross[1] = -J[1][2]*J[0][0]+J[1][0]*J[0][2]; 
   cross[2] = -J[1][0]*J[0][1]+J[1][1]*J[0][0]; 
 } else if (faceindex==3) {
   cross[0] = (J[0][1]-J[1][1])*(J[2][2])-(J[0][2]-J[1][2])*(J[2][1]); 
   cross[1] = (J[0][2]-J[1][2])*(J[2][0])-(J[0][0]-J[1][0])*(J[2][2]); 
   cross[2] = (J[0][0]-J[1][0])*(J[2][1])-(J[0][1]-J[1][1])*(J[2][0]); 
 } else if (faceindex==4) {
   cross[0] = -(J[1][1]*J[2][2]-J[1][2]*J[2][1]); 
   cross[1] = -(J[1][2]*J[2][0]-J[1][0]*J[2][2]); 
   cross[2] = -(J[1][0]*J[2][1]-J[1][1]*J[2][0]); 
 } else if (faceindex==5) {
   cross[0] = J[0][1]*J[2][2]-J[0][2]*J[2][1]; 
   cross[1] = J[0][2]*J[2][0]-J[0][0]*J[2][2]; 
   cross[2] = J[0][0]*J[2][1]-J[0][1]*J[2][0]; 
 }
}


void IsoParamUtils::getSurfDer(double *N, double *N2, double *xyz,
                       int axis, double (*sd)[3], double (*surfgrad)[2]) {

 double d[3][3];
 double dd[6][3];
 int m;
 for(m=0;m<3;m++)  {
   d[m][0] = 0.0;
   d[m][1] = 0.0;
   d[m][2] = 0.0;
 }
 for(m=0;m<6;m++)  {
   dd[m][0] = 0.0;
   dd[m][1] = 0.0;
   dd[m][2] = 0.0;
 }

 int orderc = getorderc();

 for(m=0;m<orderc;m++) {
   int j;
   for(j=0;j<3;j++) {
     d[j][0] += N[(j+1)*orderc+m]*xyz[m];
     d[j][1] += N[(j+1)*orderc+m]*xyz[m+orderc];
     d[j][2] += N[(j+1)*orderc+m]*xyz[m+2*orderc];
   }
   for(j=0;j<6;j++) {
     dd[j][0] += N2[j*orderc+m]*xyz[m];
     dd[j][1] += N2[j*orderc+m]*xyz[m+orderc];
     dd[j][2] += N2[j*orderc+m]*xyz[m+2*orderc];
   }
 }

 int iu,iv,iuu,iuv,ivv;
 if (axis==1) {
   iu = 1;
   iv = 2;
   iuu = 1;
   ivv = 2;
   iuv = 5;
 } else if (axis==2) {
   iu = 0;
   iv = 2;
   iuu = 0;
   ivv = 2;
   iuv = 4;
 } else {
   iu = 1;
   iv = 0;
   iuu = 1;
   ivv = 0;
   iuv = 3;
 }

 sd[0][0] = d[iu][0];
 sd[0][1] = d[iu][1];
 sd[0][2] = d[iu][2];
 sd[1][0] = d[iv][0];
 sd[1][1] = d[iv][1];
 sd[1][2] = d[iv][2];
 sd[2][0] = dd[iuu][0];
 sd[2][1] = dd[iuu][1];
 sd[2][2] = dd[iuu][2];
 sd[3][0] = dd[ivv][0];
 sd[3][1] = dd[ivv][1];
 sd[3][2] = dd[ivv][2];
 sd[4][0] = dd[iuv][0];
 sd[4][1] = dd[iuv][1];
 sd[4][2] = dd[iuv][2];

 int ll;
 for(ll=0;ll<orderc;ll++) {
   surfgrad[ll][0] = N[(iu+1)*orderc+ll];
   surfgrad[ll][1] = N[(iv+1)*orderc+ll];
 }
}


void IsoParamUtilsTetra::getSurfDer(double *N, double *N2, double *xyz,
                       int faceindex, double (*sd)[3], double (*surfgrad)[2]) {

 double d[3][3];
 double dd[6][3];
 int m;
 for(m=0;m<3;m++)  {
   d[m][0] = 0.0;
   d[m][1] = 0.0;
   d[m][2] = 0.0;
 }
 for(m=0;m<6;m++)  {
   dd[m][0] = 0.0;
   dd[m][1] = 0.0;
   dd[m][2] = 0.0;
 }

 int orderc = getorderc();

 for(m=0;m<orderc;m++) {
   int j;
   for(j=0;j<3;j++) {
     d[j][0] += N[(j+1)*orderc+m]*xyz[m];
     d[j][1] += N[(j+1)*orderc+m]*xyz[m+orderc];
     d[j][2] += N[(j+1)*orderc+m]*xyz[m+2*orderc];
   }
   for(j=0;j<6;j++) {
     dd[j][0] += N2[j*orderc+m]*xyz[m];
     dd[j][1] += N2[j*orderc+m]*xyz[m+orderc];
     dd[j][2] += N2[j*orderc+m]*xyz[m+2*orderc];
   }
 }

 int iu,iv,iuu,iuv,ivv;
 if (faceindex==1) {
   sd[0][0] = d[0][0]-d[2][0];
   sd[0][1] = d[0][1]-d[2][1];
   sd[0][2] = d[0][2]-d[2][2];
   sd[1][0] = d[1][0]-d[2][0];
   sd[1][1] = d[1][1]-d[2][1];
   sd[1][2] = d[1][2]-d[2][2];
   sd[2][0] = dd[0][0]+dd[2][0]-2.0*dd[4][0];
   sd[2][1] = dd[0][1]+dd[2][1]-2.0*dd[4][1];
   sd[2][2] = dd[0][2]+dd[2][2]-2.0*dd[4][2];
   sd[3][0] = dd[1][0]+dd[2][0]-2.0*dd[5][0];
   sd[3][1] = dd[1][1]+dd[2][1]-2.0*dd[5][1];
   sd[3][2] = dd[1][2]+dd[2][2]-2.0*dd[5][2];
   sd[4][0] = dd[3][0]+dd[2][0]-dd[4][0]-dd[5][0];
   sd[4][1] = dd[3][1]+dd[2][1]-dd[4][1]-dd[5][1];
   sd[4][2] = dd[3][2]+dd[2][2]-dd[4][2]-dd[5][2];
   int ll;
   for(ll=0;ll<orderc;ll++) {
     surfgrad[ll][0] = N[(0+1)*orderc+ll]-N[(2+1)*orderc+ll];
     surfgrad[ll][1] = N[(1+1)*orderc+ll]-N[(2+1)*orderc+ll];
   }
   return;
 } else if (faceindex==2) {
   iu = 1;
   iv = 2;
   iuu = 1;
   ivv = 2;
   iuv = 5;
 } else if (faceindex==3) {
   iu = 0;
   iv = 2;
   iuu = 0;
   ivv = 2;
   iuv = 4;
 } else {
   iu = 1;
   iv = 0;
   iuu = 1;
   ivv = 0;
   iuv = 3;
 }

 sd[0][0] = d[iu][0];
 sd[0][1] = d[iu][1];
 sd[0][2] = d[iu][2];
 sd[1][0] = d[iv][0];
 sd[1][1] = d[iv][1];
 sd[1][2] = d[iv][2];
 sd[2][0] = dd[iuu][0];
 sd[2][1] = dd[iuu][1];
 sd[2][2] = dd[iuu][2];
 sd[3][0] = dd[ivv][0];
 sd[3][1] = dd[ivv][1];
 sd[3][2] = dd[ivv][2];
 sd[4][0] = dd[iuv][0];
 sd[4][1] = dd[iuv][1];
 sd[4][2] = dd[iuv][2];
 int ll;
 for(ll=0;ll<orderc;ll++) {
   surfgrad[ll][0] = N[(iu+1)*orderc+ll];
   surfgrad[ll][1] = N[(iv+1)*orderc+ll];
 }
}


void IsoParamUtils::elementcenter(double *xyz, double *cxyz) {
 int ordersq = order*order;   
 int orderc = ordersq*order;
 int corner[8] = { 1-1, order-1, ordersq-order+1-1, ordersq-1,
                   ordersq*(order-1)+1-1, ordersq*(order-1)+order-1,
                   orderc-order+1-1, orderc-1};
 int i;
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
 for(i=0;i<8;i++) {
   cxyz[0] += xyz[0*orderc+corner[i]];
   cxyz[1] += xyz[1*orderc+corner[i]];
   cxyz[2] += xyz[2*orderc+corner[i]];
 }
 cxyz[0] /= 8.0;
 cxyz[1] /= 8.0; 
 cxyz[2] /= 8.0; 
}


void IsoParamUtilsTetra::elementcenter(double *xyz, double *cxyz) {
 int corner[4] = { 1-1 , order-1 , (order*(order+1))/2-1 ,
                   (order*(order+1)*(order+2))/6-1 };
 int orderc = getorderc();
 int i;
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
 for(i=0;i<4;i++) {
   cxyz[0] += xyz[0*orderc+corner[i]];
   cxyz[1] += xyz[1*orderc+corner[i]];
   cxyz[2] += xyz[2*orderc+corner[i]];
 }
 cxyz[0] /= 4.0;
 cxyz[1] /= 4.0; 
 cxyz[2] /= 4.0; 
}


void IsoParamUtilsPrism::elementcenter(double *xyz, double *cxyz) {
 int corner[6] = { 1-1 , order-1 , (order*(order+1))/2-1 ,
                   (order*(order+1))/2 * (order-1) + 1-1 ,
                   (order*(order+1))/2 * (order-1) + order-1 ,
                   (order*(order+1))/2 * (order-1) + (order*(order+1))/2-1 };
 int orderc = getorderc();
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
 for(int i=0;i<6;i++) {
   cxyz[0] += xyz[0*orderc+corner[i]];
   cxyz[1] += xyz[1*orderc+corner[i]];
   cxyz[2] += xyz[2*orderc+corner[i]];
 }
 cxyz[0] /= 6.0;
 cxyz[1] /= 6.0; 
 cxyz[2] /= 6.0; 
}


void IsoParamUtilsPyramid::elementcenter(double *xyz, double *cxyz) {
 int orderc = 5;
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
 for(int i=0;i<5;i++) {
   cxyz[0] += xyz[0*orderc+i];
   cxyz[1] += xyz[1*orderc+i];
   cxyz[2] += xyz[2*orderc+i];
 }
 cxyz[0] /= 5.0;
 cxyz[1] /= 5.0; 
 cxyz[2] /= 5.0; 
}


void IsoParamUtils::cornerindeces(int faceindex, int *ci) {
 int ordersq = order*order;
 int orderc = ordersq*order;
 if (faceindex==1) {
   ci[0] = 1 - 1;
   ci[1] = ordersq-order+1 - 1;
   ci[2] = ordersq - 1;
   ci[3] = order - 1;
 } else if (faceindex==2) {
   ci[0] = 1 - 1;
   ci[1] = order - 1;
   ci[2] = ordersq*(order-1)+order - 1;
   ci[3] = ordersq*(order-1)+1 - 1;
 } else if (faceindex==3) {
   ci[0] = order - 1;
   ci[1] = ordersq - 1;
   ci[2] = orderc - 1;
   ci[3] = ordersq*(order-1)+order - 1;
 } else if (faceindex==4) {
   ci[0] = ordersq-order+1 - 1;
   ci[1] = orderc-order+1 - 1;
   ci[2] = orderc - 1;
   ci[3] = ordersq - 1;
 } else if (faceindex==5) {
   ci[0] = 1 - 1;
   ci[1] = ordersq*(order-1)+1 - 1;
   ci[2] = orderc-order+1 - 1;
   ci[3] = ordersq-order+1 - 1;
 } else {
   ci[0] = ordersq*(order-1)+1 - 1;
   ci[1] = ordersq*(order-1)+order - 1;
   ci[2] = orderc - 1;
   ci[3] = orderc-order+1 - 1;
 }
}


void IsoParamUtilsTetra::cornerindeces(int faceindex, int *ci) {
// 1-1 , order-1 , (order*(order+1))/2-1 , (order*(order+1)*(order+2))/6-1
 if (faceindex==1) {
   ci[0] = order-1;
   ci[1] = (order*(order+1))/2-1;
   ci[2] = (order*(order+1)*(order+2))/6-1;
 } else if (faceindex==2) {
   ci[0] = 1 - 1;
   ci[1] = (order*(order+1)*(order+2))/6-1;
   ci[2] = (order*(order+1))/2-1;
 } else if (faceindex==3) {
   ci[0] = 1 - 1;
   ci[1] = order - 1;
   ci[2] = (order*(order+1)*(order+2))/6-1;
 } else {
   ci[0] = 1 - 1;
   ci[1] = (order*(order+1))/2-1;
   ci[2] = order - 1;
 }
}


void IsoParamUtilsPrism::cornerindeces(int faceindex, int *ci) {
 int o = order;
 int osq = (order*(order+1))/2;
 if (faceindex==1) {
   ci[0] = 1 - 1;
   ci[1] = osq - 1;
   ci[2] = o - 1;
 } else if (faceindex==2) {
   ci[0] = (o-1)*osq + 1 - 1;
   ci[1] = (o-1)*osq + o - 1;
   ci[2] = (o-1)*osq + osq - 1;
 } else if (faceindex==3) {
   ci[0] = o - 1;
   ci[1] = osq - 1;
   ci[2] = (o-1)*osq + osq - 1;
   ci[3] = (o-1)*osq + o - 1;
 } else if (faceindex==4) {
   ci[0] = osq - 1;
   ci[1] = 1 - 1;
   ci[2] = (o-1)*osq + 1 - 1;
   ci[3] = (o-1)*osq + osq - 1;
 } else {
   ci[0] = 1 - 1;
   ci[1] = o - 1;
   ci[2] = (o-1)*osq + o - 1;
   ci[3] = (o-1)*osq + 1 - 1;
 }
}


void IsoParamUtilsPyramid::cornerindeces(int faceindex, int *ci) {
 if (faceindex==1) {
   ci[0] = 0;
   ci[1] = 1;
   ci[2] = 3;
   ci[3] = 2;
 } else if (faceindex==2) {
   ci[0] = 0;
   ci[1] = 1;
   ci[2] = 4;
 } else if (faceindex==3) {
   ci[0] = 1;
   ci[1] = 3;
   ci[2] = 4;
 } else if (faceindex==4) {
   ci[0] = 3;
   ci[1] = 2;
   ci[2] = 4;
 } else {
   ci[0] = 2;
   ci[1] = 0;
   ci[2] = 4;
 }
}


void IsoParamUtils::sidecenter(double *xyz, int faceindex, double *scxyz) {
 int orderc = getorderc();
 
 int corner[4];
 cornerindeces(faceindex,corner);
  
 int i;
 scxyz[0] = scxyz[1] = scxyz[2] = 0.0;
 for(i=0;i<4;i++) {
   scxyz[0] += xyz[0*orderc+corner[i]];
   scxyz[1] += xyz[1*orderc+corner[i]];
   scxyz[2] += xyz[2*orderc+corner[i]];
 }
 scxyz[0] /= 4.0;
 scxyz[1] /= 4.0;
 scxyz[2] /= 4.0;
}


void IsoParamUtilsTetra::sidecenter(double *xyz, int faceindex, double *scxyz) {
 int orderc = getorderc();
 
 int corner[3];
 cornerindeces(faceindex,corner);
  
 int i;
 scxyz[0] = scxyz[1] = scxyz[2] = 0.0;
 for(i=0;i<3;i++) {
   scxyz[0] += xyz[0*orderc+corner[i]];
   scxyz[1] += xyz[1*orderc+corner[i]];
   scxyz[2] += xyz[2*orderc+corner[i]];
 }
 scxyz[0] /= 3.0;
 scxyz[1] /= 3.0;
 scxyz[2] /= 3.0;
}


void IsoParamUtilsPrism::sidecenter(double *xyz, int faceindex, double *scxyz) {
 int orderc = getorderc();
 
 int corner[4];
 cornerindeces(faceindex,corner);
 int nc = (faceindex<=2)?3:4;
  
 int i;
 scxyz[0] = scxyz[1] = scxyz[2] = 0.0;
 for(i=0;i<nc;i++) {
   scxyz[0] += xyz[0*orderc+corner[i]];
   scxyz[1] += xyz[1*orderc+corner[i]];
   scxyz[2] += xyz[2*orderc+corner[i]];
 }
 scxyz[0] /= nc;
 scxyz[1] /= nc;
 scxyz[2] /= nc;
}


void IsoParamUtilsPyramid::sidecenter(double *xyz, int faceindex, double *scxyz) {
 int orderc = 5;

 int corner[4];
 cornerindeces(faceindex,corner);
 int nc = (faceindex<=1)?4:3;
  
 int i;
 scxyz[0] = scxyz[1] = scxyz[2] = 0.0;
 for(i=0;i<nc;i++) {
   scxyz[0] += xyz[0*orderc+corner[i]];
   scxyz[1] += xyz[1*orderc+corner[i]];
   scxyz[2] += xyz[2*orderc+corner[i]];
 }
 scxyz[0] /= nc;
 scxyz[1] /= nc;
 scxyz[2] /= nc;
}


void IsoParamUtils::faceindeces(int faceindex, int *fi) {
 int orderc = getorderc();
 int ordersq = getordersq();
 int k,m;
 if (faceindex==1) {
   for(k=0;k<order;k++) for(m=0;m<order;m++)
     fi[k*order+m] = k*order+m;
 }
 else if (faceindex==2) {
   for(k=0;k<order;k++) for(m=0;m<order;m++)
     fi[k*order+m] = k*ordersq+m;
 }
 else if (faceindex==3) {
   for(k=0;k<order;k++) for(m=0;m<order;m++)
     fi[k*order+m] = k*ordersq+(m+1)*order-1;
 }
 else if (faceindex==4) {
   for(k=0;k<order;k++) for(m=0;m<order;m++)
     fi[k*order+m] = k*ordersq+m+ordersq-order;
 }
 else if (faceindex==5) {
   for(k=0;k<order;k++) for(m=0;m<order;m++)
     fi[k*order+m] = k*ordersq+m*order;
 }
 else {
   for(k=0;k<order;k++) for(m=0;m<order;m++)
     fi[k*order+m] = k*order+m+orderc-ordersq;
 }
}


void IsoParamUtilsTetra::faceindeces(int faceindex, int *fi) {
 //int orderc = getorderc();
 //int ordersq = getordersq();
 int c = 0;
 int fc = 0;
 int r,s,t;
 for(t=0;t<order;t++) {
   for(s=0;s<order-t;s++) {
     for(r=0;r<order-t-s;r++) {
       int u = order - 1 - r - s - t;
       if (faceindex==1 && u==0) { 
         fi[fc] = c;
         fc++;
       } else if (faceindex==2 && r==0) {
         fi[fc] = c;
         fc++;
       } else if (faceindex==3 && s==0) {
         fi[fc] = c;
         fc++;
       } else if (faceindex==4 && t==0) {
         fi[fc] = c;
         fc++;
       }
       c++;
     }
   }
 }
}


void IsoParamUtilsPrism::faceindeces(int faceindex, int *fi) {
 //int orderc = getorderc();
 //int ordersq = getordersq();
 int osq = (order*(order+1))/2;
 int c=0,fc=0;
 if (faceindex==1) 
   for(int i=0;i<osq;i++) fi[i] = i;
 else if (faceindex==2)
   for(int i=0;i<osq;i++) fi[i] = (order-1)*osq+i;
 else if (faceindex>=3) {
   for(int u=0;u<order;u++) {
     for(int s=0;s<order;s++) {
       for( int r=0;r<order-s;r++) {
         int t = order - 1 - r - s;
         if (faceindex==3 && t==0) {
           fi[fc] = c + u*osq;
           fc++;
         } else if (faceindex==4 && s==0) {
           fi[fc] = c + u*osq;
           fc++;
         } else if (faceindex==5 && r==0) {
           fi[fc] = c + u*osq;
           fc++;
         }
         c++;
       }
     }
   }
 } 
}


void IsoParamUtilsPyramid::faceindeces(int faceindex, int *fi) {
 if (faceindex==1)  {
   for(int i=0;i<4;i++) fi[i] = i;
 } else if (faceindex==2) {
   fi[0] = 0; fi[1] = 1; fi[2] = 4;
 } else if (faceindex==3) {
   fi[0] = 1; fi[1] = 3; fi[2] = 4;
 } else if (faceindex==4) {
   fi[0] = 3; fi[1] = 2; fi[2] = 4;
 } else {
   fi[0] = 2; fi[1] = 0; fi[2] = 4;
 } 
}


void IsoParamUtils::facemap(int &faceindex, int *fn, int *n, int *map) {

 int ordersq = getordersq();
 int fc[4] = { fn[1-1], fn[order-1], fn[order*order-1],
               fn[order*order-1-(order-1)] };
 int ic[4];
 int i;
 for(i=1;i<=6;i++) {
   int matches = 0;
   cornerindeces(i,ic);
   int j,k;
   for(j=0;j<4;j++) for(k=0;k<4;k++)
     if (fc[j]==n[ic[k]]) {
        matches++;
     }
   if (matches==4) { break; }
 }
 faceindex = i;

 int *fi = (int*) alloca(sizeof(int)*ordersq);
 faceindeces(faceindex,fi);

// This could be more efficient by sorting first ...
 for(i=0;i<ordersq;i++) {
   int j;
   for(j=0;j<ordersq;j++) {
     if (n[fi[i]]==fn[j]) map[j] = fi[i];
   } 
 }
}


void IsoParamUtilsTetra::facemap(int &faceindex, int *fn, int *n, int *map) {

 int ordersq = getordersq();
 int fc[3] = { fn[1-1], fn[order-1], fn[ordersq-1]};
 int ic[3];
 int i;
 for(i=1;i<=4;i++) {
   int matches = 0;
   cornerindeces(i,ic);
   int j,k;
   for(j=0;j<3;j++) for(k=0;k<3;k++)
     if (fc[j]==n[ic[k]]) {
        matches++;
     }
   if (matches==3) { break; }
 }
 faceindex = i;

 int *fi = (int*) alloca(sizeof(int)*ordersq);
 faceindeces(faceindex,fi);

// This could be more efficient by sorting first ...
 for(i=0;i<ordersq;i++) {
   int j;
   for(j=0;j<ordersq;j++) {
     if (n[fi[i]]==fn[j]) map[j] = fi[i];
   } 
 }
}


int IsoParamUtils::isFlat(double *xyz, int faceindex) {

 int orderc = getorderc();
 int corneri[4];
 cornerindeces(faceindex, corneri);

 double xf[3],af[3],bf[3];
 xf[0] = xyz[0*orderc+corneri[0]];
 xf[1] = xyz[1*orderc+corneri[0]];
 xf[2] = xyz[2*orderc+corneri[0]];
 af[0] = xyz[0*orderc+corneri[1]] - xf[0];
 af[1] = xyz[1*orderc+corneri[1]] - xf[1];
 af[2] = xyz[2*orderc+corneri[1]] - xf[2];
 bf[0] = xyz[0*orderc+corneri[3]] - xf[0];
 bf[1] = xyz[1*orderc+corneri[3]] - xf[1];
 bf[2] = xyz[2*orderc+corneri[3]] - xf[2];
// Orthogonalize [af, bf]
 double l = sidelen(af[0],af[1],af[2]);
 af[0] /= l;
 af[1] /= l;
 af[2] /= l;
 double alpha = af[0]*bf[0]+af[1]*bf[1]+af[2]*bf[2];
 bf[0] -= alpha*af[0];
 bf[1] -= alpha*af[1];
 bf[2] -= alpha*af[2];
 l = sidelen(bf[0],bf[1],bf[2]);
 bf[0] /= l;
 bf[1] /= l;
 bf[2] /= l;

 double nf[3];
 crossproduct(af,bf,nf);

 int ordersq = getordersq();
 int *fi = (int*) alloca(sizeof(int)*ordersq);
 faceindeces(faceindex, fi);

 int flat = 1;
 int i;
 for(i=0;i<ordersq;i++) {
   double len = sidelen(
                        xyz[0*orderc+fi[i]] - xf[0],
                        xyz[1*orderc+fi[i]] - xf[1],
                        xyz[2*orderc+fi[i]] - xf[2]);
   if (fabs(
        nf[0]*(xyz[0*orderc+fi[i]] - xf[0])+
        nf[1]*(xyz[1*orderc+fi[i]] - xf[1])+
        nf[2]*(xyz[2*orderc+fi[i]] - xf[2]))/len > 1e-6) {
          flat = 0;
          break;
   }
 }
 return flat;
}


int IsoParamUtils::isStraight(double *xyz, int faceindex, int edgeindex) {

 int orderc = getorderc();
 int &o = order;
 int edges[12][2] = {
        { 0,1 },
        { o-1,o },
        { o*o-1,-1 },
        { o*o-o,-o },
        { 0,o*o },
        { o-1,o*o },
        { o*o-1,o*o },
        { o*o-o,o*o },
        { o*o*(o-1)+0,1 },
        { o*o*(o-1)+o-1,o },
        { o*o*(o-1)+o*o-1,-1 },
        { o*o*(o-1)+o*o-o,-o }};
 int fe[6][4] = {
        {1,2,3,4},
        {5,9,6,1},
        {6,10,7,2},
        {3,7,11,8},
        {4,8,12,5},
        {12,11,10,9}};


 int e = fe[faceindex-1][edgeindex]-1;
 int n1 = edges[e][0];
 int step = edges[e][1];
 int n2 = edges[e][0]+(o-1)*step;

 double xyz1[3] = {xyz[0*orderc+n1],xyz[1*orderc+n1],xyz[2*orderc+n1]};
 double xyz2[3] = {xyz[0*orderc+n2],xyz[1*orderc+n2],xyz[2*orderc+n2]};
 double dx[3] = {
   xyz2[0]-xyz1[0],
   xyz2[1]-xyz1[1],
   xyz2[2]-xyz1[2]
 };

 int j;
 int straight=1;
 for(j=1;j<o-1;j++) {
   int nn = n1+j*step;
   double xyzn[3] = {xyz[0*orderc+nn],xyz[1*orderc+nn],xyz[2*orderc+nn]};
   double dxn[3] = {
     xyzn[0]-xyz1[0],
     xyzn[1]-xyz1[1],
     xyzn[2]-xyz1[2]
   };
   double cross[3];
   crossproduct(dx,dxn,cross);
   double normdx = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
   double normdxn = sqrt(dxn[0]*dxn[0]+dxn[1]*dxn[1]+dxn[2]*dxn[2]);
   double normcross = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+
                           cross[2]*cross[2]);
   if (normcross/(normdx*normdxn)>1e-6) {
     straight = 0;
     break;
   }
 }
 return straight;
}


void IsoParamUtils::adjacentfaces(int faceindex, int *baf) {
 int adjacentfaces[6][5] = { 
                            {1, 2, 3, 4, 5},
                            {1, 2, 3, 5, 6},
                            {1, 2, 3, 4, 6},
                            {1, 3, 4, 5, 6},
                            {1, 2, 4, 5, 6},
                            {2, 3, 4, 5, 6}
                           };
 
 baf[adjacentfaces[faceindex-1][0]-1] = 1;
 baf[adjacentfaces[faceindex-1][1]-1] = 1;
 baf[adjacentfaces[faceindex-1][2]-1] = 1;
 baf[adjacentfaces[faceindex-1][3]-1] = 1;
 baf[adjacentfaces[faceindex-1][4]-1] = 1;
}


void IsoParamUtils::lineInt3d(double *xyz, int faceindex, int edgeindex,
               double *nf, IntegFunctionL3d &f, int gorder) {

 //int ordersq = order*order;
 int orderc = order*order*order;

 int i;

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 int fei[6][4][4] = {
  {{-1,-1,1,3}, // face 1
   { 1,-1,2,3},
   { 1,-1,1,3},
   {-1,-1,2,3}},
  {{-1,-1,3,2}, // face 2
   {-1, 1,1,2},
   { 1,-1,3,2},
   {-1,-1,1,2}},
  {{ 1,-1,3,1}, // face 3
   { 1, 1,2,1},
   { 1, 1,3,1},
   { 1,-1,2,1}},
  {{ 1,-1,1,2}, // face 4
   { 1, 1,3,2},
   { 1, 1,1,2},
   {-1, 1,3,2}},
  {{-1,-1,2,1}, // face 5
   {-1, 1,3,1},
   {-1, 1,2,1},
   {-1,-1,3,1}},
  {{-1, 1,2,3}, // face 6
   { 1, 1,1,3},
   { 1, 1,2,3},
   {-1, 1,1,3}}};

 double p = double(fei[faceindex-1][edgeindex][0]); 
 double q = double(fei[faceindex-1][edgeindex][1]); 
 int axis = fei[faceindex-1][edgeindex][2]; 
 int faxis = fei[faceindex-1][edgeindex][3]; 

 double *sNN1 = (double*)dbg_alloca(sizeof(double)*2*order);
 double *sNN2 = (double*)dbg_alloca(sizeof(double)*2*order);
 lagGalShapeFunction(1,&p,sNN1);
 lagGalShapeFunction(1,&q,sNN2);

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {
   int ii,jj,kk;
   int ng1,ng2,ng3;
   double *NN1, *NN2, *NN3;
   NN1 = NN2 = NN3 = NN;
   if (axis==1) {
     ii = i; jj = 0; kk = 0;
     ng1 = ng; ng2 = 1; ng3 = 1;
     NN2 = sNN1;
     NN3 = sNN2;
   } else if (axis==2) {
     ii = 0; jj = i; kk = 0;
     ng1 = 1; ng2 = ng; ng3 = 1;
     NN1 = sNN1;
     NN3 = sNN2;
   } else {
     ii = 0; jj = 0; kk = i;
     ng1 = 1; ng2 = 1; ng3 = ng;
     NN1 = sNN1;
     NN2 = sNN2;
   }

   double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
   int mx,my,mz;
   for(mx=0;mx<order;mx++) {
     for(my=0;my<order;my++) {
       for(mz=0;mz<order;mz++) {
          N[mz*order*order+my*order+mx] =
            NN1[ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
          N[orderc+mz*order*order+my*order+mx] =
            NN1[ng1*order+ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
          N[2*orderc+mz*order*order+my*order+mx] =
            NN1[ii*order+mx]*NN2[ng2*order+jj*order+my]*NN3[kk*order+mz];
          N[3*orderc+mz*order*order+my*order+mx] =
            NN1[ii*order+mx]*NN2[jj*order+my]*NN3[ng3*order+kk*order+mz];
       }
     }
   }
         
   double w = gr.wgauss[i];

   double J[3][3];
   jmatrix(N,xyz,J);

   double cross[3];
   crossj(J,faxis,cross);
  
   double nsign;
   if (cross[0]*(scxyz[0]-cxyz[0])+
       cross[1]*(scxyz[1]-cxyz[1])+
       cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1.0;
   else nsign = -1.0;

   double tau[3];
   if (axis==1) {
     tau[0] = J[0][0];
     tau[1] = J[0][1];
     tau[2] = J[0][2];
   } else if (axis==2) {
     tau[0] = J[1][0];
     tau[1] = J[1][1];
     tau[2] = J[1][2];
   } else {
     tau[0] = J[2][0];
     tau[1] = J[2][1];
     tau[2] = J[2][2];
   }

   double x[3] = {0.0,0.0,0.0};
   int m;
   for(m=0;m<orderc;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+orderc];
     x[2] += N[m]*xyz[m+2*orderc];
   }

   double ns[3];
   crossproduct(cross,tau,ns);
   double tsign;
   if (ns[0]*(scxyz[0]-x[0])+
       ns[1]*(scxyz[1]-x[1])+
       ns[2]*(scxyz[2]-x[2])<0.0) tsign = 1.0;
   else tsign = -1.0;
 
   f.evaluate(x,N,cross,nsign,tau,tsign,w);
 }
}


void IsoParamUtils::surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f,
                              int gorder) {
 //int ordersq = order*order;
 int orderc = order*order*order;

 int i;

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);
 double p;
 double *sNN = (double*)dbg_alloca(sizeof(double)*2*order);
 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);

 int axis = 0;
 if (faceindex==1) {
   p = -1.0;
   axis = 3;
 } else if (faceindex==2) {
   p = -1.0;
   axis = 2;
 } else if (faceindex==3) {
   p = 1.0;
   axis = 1;
 } else if (faceindex==4) {
   p = 1.0;
   axis = 2;
 } else if (faceindex==5) {
   p = -1.0;
   axis = 1;
 } else {
   p = 1.0;
   axis = 3;
 }
 lagGalShapeFunction(1,&p,sNN);

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {
     int ii,jj,kk;
     int ng1,ng2,ng3;
     double *NN1, *NN2, *NN3;
     NN1 = NN2 = NN3 = NN;
     if (axis==1) {
       ii = 0; jj = i; kk = j;
       ng1 = 1; ng2 = ng; ng3 = ng;
       NN1 = sNN;
     } else if (axis==2) {
       ii = i; jj = 0; kk = j;
       ng1 = ng; ng2 = 1; ng3 = ng;
       NN2 = sNN;
     } else {
       ii = i; jj = j; kk = 0;
       ng1 = ng; ng2 = ng; ng3 = 1;
       NN3 = sNN;
     }

     int mx,my,mz;
     for(mx=0;mx<order;mx++) {
       for(my=0;my<order;my++) {
          for(mz=0;mz<order;mz++) {
            N[mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[orderc+mz*order*order+my*order+mx] =
              NN1[ng1*order+ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[2*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[ng2*order+jj*order+my]*NN3[kk*order+mz];
            N[3*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[ng3*order+kk*order+mz];
          }
        }
      }
           
     double w = gr.wgauss[i]*gr.wgauss[j];

     double J[3][3];
     jmatrix(N,xyz,J);

     double cross[3];
     crossj(J,axis,cross);
    
     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     double x[3] = {0.0,0.0,0.0};
     int m;
     for(m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }
   
     f.evaluate(x,N,cross,nsign,w);
   }
 }
}


void IsoParamUtils::surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f,
                              int gorder) {
 //int ordersq = order*order;
 int orderc = order*order*order;

 int i;

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);
 double p;
 double *sNN = (double*)dbg_alloca(sizeof(double)*2*order);
 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);

 int axis = 0;
 if (faceindex==1) {
   p = -1.0;
   axis = 3;
 } else if (faceindex==2) {
   p = -1.0;
   axis = 2;
 } else if (faceindex==3) {
   p = 1.0;
   axis = 1;
 } else if (faceindex==4) {
   p = 1.0;
   axis = 2;
 } else if (faceindex==5) {
   p = -1.0;
   axis = 1;
 } else {
   p = 1.0;
   axis = 3;
 }
 lagGalShapeFunction(1,&p,sNN);

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {
     int ii,jj,kk;
     int ng1,ng2,ng3;
     double *NN1, *NN2, *NN3;
     NN1 = NN2 = NN3 = NN;
     if (axis==1) {
       ii = 0; jj = i; kk = j;
       ng1 = 1; ng2 = ng; ng3 = ng;
       NN1 = sNN;
     } else if (axis==2) {
       ii = i; jj = 0; kk = j;
       ng1 = ng; ng2 = 1; ng3 = ng;
       NN2 = sNN;
     } else {
       ii = i; jj = j; kk = 0;
       ng1 = ng; ng2 = ng; ng3 = 1;
       NN3 = sNN;
     }

     int mx,my,mz;
     for(mx=0;mx<order;mx++) {
       for(my=0;my<order;my++) {
          for(mz=0;mz<order;mz++) {
            N[mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[orderc+mz*order*order+my*order+mx] =
              NN1[ng1*order+ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[2*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[ng2*order+jj*order+my]*NN3[kk*order+mz];
            N[3*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[ng3*order+kk*order+mz];
          }
        }
      }
           
     double w = gr.wgauss[i]*gr.wgauss[j];

     double J[3][3];
     jmatrix(N,xyz,J);

     double cross[3];
     double tau1[3];
     double tau2[3];
     crossj(J,axis,cross,tau1,tau2);
    
     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     double x[3] = {0.0,0.0,0.0};
     int m;
     for(m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }
   
     f.evaluate(x,N,tau1,tau2,nsign,w);
   }
 }
}


void IsoParamUtils::surfGInt3d(double *xyz, int faceindex, IntegFunctionAG3d &f,
                              int gorder) {

 //int ordersq = order*order;
 int orderc = order*order*order;

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
   axis = 3;
 } else if (faceindex==2) {
   p = -1.0;
   axis = 2;
 } else if (faceindex==3) {
   p = 1.0;
   axis = 1;
 } else if (faceindex==4) {
   p = 1.0;
   axis = 2;
 } else if (faceindex==5) {
   p = -1.0;
   axis = 1;
 } else {
   p = 1.0;
   axis = 3;
 }
 lagGalShapeFunction(1,&p,sNN);

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int ng = gr.ngauss;
 double *N = (double*) alloca(sizeof(double)*4*orderc);
 double (*dNdx)[3] = (double(*)[3])alloca(sizeof(double)*orderc*3);
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {
     int ii,jj,kk;
     int ng1,ng2,ng3;
     double *NN1, *NN2, *NN3;
     NN1 = NN2 = NN3 = NN;
     if (axis==1) {
       ii = 0; jj = i; kk = j;
       ng1 = 1; ng2 = ng; ng3 = ng;
       NN1 = sNN;
     } else if (axis==2) {
       ii = i; jj = 0; kk = j;
       ng1 = ng; ng2 = 1; ng3 = ng;
       NN2 = sNN;
     } else {
       ii = i; jj = j; kk = 0;
       ng1 = ng; ng2 = ng; ng3 = 1;
       NN3 = sNN;
     }

     int mx,my,mz;
     for(mx=0;mx<order;mx++) {
       for(my=0;my<order;my++) {
          for(mz=0;mz<order;mz++) {
            N[mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[orderc+mz*order*order+my*order+mx] =
              NN1[ng1*order+ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[2*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[ng2*order+jj*order+my]*NN3[kk*order+mz];
            N[3*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[ng3*order+kk*order+mz];
          }
        }
      }
           
     double w = gr.wgauss[i]*gr.wgauss[j];

     double J[3][3];
     jmatrix(N,xyz,J);

     double cross[3];
     crossj(J,axis,cross);
    
     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     double det = detj(J);
     double ij[3][3];
     invj(J,det,ij);
     int m;
     for(m=0;m<orderc;m++) {
       dNdx[m][0] =
        ij[0][0]*N[orderc+m]+ij[0][1]*N[2*orderc+m]+ij[0][2]*N[3*orderc+m];
       dNdx[m][1] =
        ij[1][0]*N[orderc+m]+ij[1][1]*N[2*orderc+m]+ij[1][2]*N[3*orderc+m];
       dNdx[m][2] =
        ij[2][0]*N[orderc+m]+ij[2][1]*N[2*orderc+m]+ij[2][2]*N[3*orderc+m];
     }

     double x[3] = {0.0,0.0,0.0};
     for(m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }
   
     f.evaluate(x,N,dNdx,cross,nsign,w);
   }
 }
}


void IsoParamUtils::surfSurfInt3d(double *xyz, IntegFunctionA3d &f, int gorder) {

 int ordersq = order*order;

 int i;

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double *N = (double*) alloca(sizeof(double)*3*ordersq);
 int ng = gr.ngauss;
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

     double J[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
     int m;
     for(m=0;m<ordersq;m++) {
       J[0][0] += N[ordersq+m]*xyz[m];
       J[0][1] += N[ordersq+m]*xyz[m+ordersq];
       J[0][2] += N[ordersq+m]*xyz[m+2*ordersq];
       J[1][0] += N[2*ordersq+m]*xyz[m];
       J[1][1] += N[2*ordersq+m]*xyz[m+ordersq];
       J[1][2] += N[2*ordersq+m]*xyz[m+2*ordersq];
     }

     double cross[3];
     cross[0] = J[0][1]*J[1][2]-J[0][2]*J[1][1];
     cross[1] = J[0][2]*J[1][0]-J[0][0]*J[1][2];
     cross[2] = J[0][0]*J[1][1]-J[0][1]*J[1][0];

     double x[3] = {0.0,0.0,0.0};
     for(m=0;m<ordersq;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+ordersq];
       x[2] += N[m]*xyz[m+2*ordersq];
     }
   
     double nsign = 1.0; 
     f.evaluate(x,N,cross,nsign,w);
   }
 }
}


void IsoParamUtils::spectralSurfSurfInt3d(double *xyz, IntegFunctionA3d &f) {

 int ordersq = order*order;
 GaussLobattoRuleLine gr(order);
 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*gr.ngauss);
 spectralLagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double *N = (double*) alloca(sizeof(double)*3*ordersq);
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

     double J[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
     int m;
     for(m=0;m<ordersq;m++) {
       J[0][0] += N[ordersq+m]*xyz[m];
       J[0][1] += N[ordersq+m]*xyz[m+ordersq];
       J[0][2] += N[ordersq+m]*xyz[m+2*ordersq];
       J[1][0] += N[2*ordersq+m]*xyz[m];
       J[1][1] += N[2*ordersq+m]*xyz[m+ordersq];
       J[1][2] += N[2*ordersq+m]*xyz[m+2*ordersq];
     }

     double cross[3];
     cross[0] = J[0][1]*J[1][2]-J[0][2]*J[1][1];
     cross[1] = J[0][2]*J[1][0]-J[0][0]*J[1][2];
     cross[2] = J[0][0]*J[1][1]-J[0][1]*J[1][0];

     double x[3] = {0.0,0.0,0.0};
     for(m=0;m<ordersq;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+ordersq];
       x[2] += N[m]*xyz[m+2*ordersq];
     }
   
     double nsign = 1.0; 
     f.evaluate(x,N,cross,nsign,w);
   }
 }
}


void IsoParamUtils::surfCurvInt3d(double *xyz, int faceindex, IntegFunctionAC3d &f,
                                  int gorder) {

 //int ordersq = order*order;
 int orderc = order*order*order;

 int i;

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)dbg_alloca(sizeof(double)*3*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN,1);
 double p;
 double *sNN = (double*)dbg_alloca(sizeof(double)*3*order);

 int axis = 0;
 if (faceindex==1) {
   p = -1.0;
   axis = 3;
 } else if (faceindex==2) {
   p = -1.0;
   axis = 2;
 } else if (faceindex==3) {
   p = 1.0;
   axis = 1;
 } else if (faceindex==4) {
   p = 1.0;
   axis = 2;
 } else if (faceindex==5) {
   p = -1.0;
   axis = 1;
 } else {
   p = 1.0;
   axis = 3;
 }
 lagGalShapeFunction(1,&p,sNN,1);

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int ng = gr.ngauss;
 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
 double *N2 = (double*) dbg_alloca(sizeof(double)*6*orderc);
 double (*surfgrad)[2] = (double(*)[2])dbg_alloca(sizeof(double)*2*orderc);
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {
     int ii,jj,kk;
     int ng1,ng2,ng3;
     double *NN1, *NN2, *NN3;
     NN1 = NN2 = NN3 = NN;
     if (axis==1) {
       ii = 0; jj = i; kk = j;
       ng1 = 1; ng2 = ng; ng3 = ng;
       NN1 = sNN;
     } else if (axis==2) {
       ii = i; jj = 0; kk = j;
       ng1 = ng; ng2 = 1; ng3 = ng;
       NN2 = sNN;
     } else {
       ii = i; jj = j; kk = 0;
       ng1 = ng; ng2 = ng; ng3 = 1;
       NN3 = sNN;
     }

     int mx,my,mz;
     for(mx=0;mx<order;mx++) {
       for(my=0;my<order;my++) {
         for(mz=0;mz<order;mz++) {
            N[mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[orderc+mz*order*order+my*order+mx] =
              NN1[ng1*order+ii*order+mx]*NN2[jj*order+my]*NN3[kk*order+mz];
            N[2*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[ng2*order+jj*order+my]*NN3[kk*order+mz];
            N[3*orderc+mz*order*order+my*order+mx] =
              NN1[ii*order+mx]*NN2[jj*order+my]*NN3[ng3*order+kk*order+mz];

            N2[0*orderc+mz*order*order+my*order+mx] =
                   NN1[2*ng1*order+ii*order+mx]*
                   NN2[0*ng1*order+jj*order+my]*
                   NN3[0*ng3*order+kk*order+mz];
            N2[1*orderc+mz*order*order+my*order+mx] =
                   NN1[0*ng1*order+ii*order+mx]*
                   NN2[2*ng1*order+jj*order+my]*
                   NN3[0*ng3*order+kk*order+mz];
            N2[2*orderc+mz*order*order+my*order+mx] =
                   NN1[0*ng1*order+ii*order+mx]*
                   NN2[0*ng1*order+jj*order+my]*
                   NN3[2*ng3*order+kk*order+mz];
            N2[3*orderc+mz*order*order+my*order+mx] =
                   NN1[1*ng1*order+ii*order+mx]*
                   NN2[1*ng1*order+jj*order+my]*
                   NN3[0*ng3*order+kk*order+mz];
            N2[4*orderc+mz*order*order+my*order+mx] =
                   NN1[1*ng1*order+ii*order+mx]*
                   NN2[0*ng1*order+jj*order+my]*
                   NN3[1*ng3*order+kk*order+mz];
            N2[5*orderc+mz*order*order+my*order+mx] =
                   NN1[0*ng1*order+ii*order+mx]*
                   NN2[1*ng1*order+jj*order+my]*
                   NN3[1*ng3*order+kk*order+mz];
         }
       }
     }
           
     double w = gr.wgauss[i]*gr.wgauss[j];

     double J[3][3];
     jmatrix(N,xyz,J);

     double cross[3];
     crossj(J,axis,cross);

     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     double sd[5][3];
     getSurfDer(N,N2,xyz,axis,sd,surfgrad);

     double x[3] = {0.0,0.0,0.0};
     int m;
     for(m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }
    
     f.evaluate(x,N,cross,nsign, sd, surfgrad, w);
   }
 }
}


void IsoParamUtilsTetra::surfInt3d(double *xyz, int faceindex,
                                   IntegFunctionA3d &f, int gorder ) {

 int orderc = getorderc();

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int i;

 GaussRuleTriangle gr(gorder);
 gr.pool = (double*)dbg_alloca(sizeof(double)*5*gr.ngauss);
 gr.init(faceindex);

 double *NN = (double*)dbg_alloca(sizeof(double)*2*4*gr.ngauss*order);
 for(i=1;i<=order;i++)
   lagGalShapeFunction(i,4*gr.ngauss,gr.xigauss,NN+(i-1)*2*4*gr.ngauss);

 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {

   int r,s,t;
   int c=0;
   for(t=0;t<order;t++) {
     for(s=0;s<order-t;s++) {
       for(r=0;r<order-t-s;r++) {
         int u = order - 1 - r - s - t;
         int rr = r*2*4*gr.ngauss+0*gr.ngauss+i;
         int ss = s*2*4*gr.ngauss+1*gr.ngauss+i;
         int tt = t*2*4*gr.ngauss+2*gr.ngauss+i;
         int uu = u*2*4*gr.ngauss+3*gr.ngauss+i;
         int rrd = rr+4*gr.ngauss;
         int ssd = ss+4*gr.ngauss;
         int ttd = tt+4*gr.ngauss;
         int uud = uu+4*gr.ngauss;
         double scaling = 1.0;
         N[c] = scaling*NN[rr]*NN[ss]*NN[tt]*NN[uu];
         N[1*orderc+c] = scaling*(NN[rrd]*NN[ss]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[2*orderc+c] = scaling*(NN[rr]*NN[ssd]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[3*orderc+c] = scaling*(NN[rr]*NN[ss]*NN[ttd]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         c++;
       }
     }
   }

   double w = gr.wgauss[i];

   double J[3][3];
   jmatrix(N,xyz,J);

   double cross[3];
   crossj(J,faceindex,cross);
  
   double nsign;
   if (cross[0]*(scxyz[0]-cxyz[0])+
       cross[1]*(scxyz[1]-cxyz[1])+
       cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
   else nsign = -1;

   double x[3] = {0.0,0.0,0.0};
   int m;
   for(m=0;m<orderc;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+orderc];
     x[2] += N[m]*xyz[m+2*orderc];
   }

   f.evaluate(x,N,cross,nsign,w); 
 }
}


void IsoParamUtilsTetra::surfGInt3d(double *xyz, int faceindex,
                                   IntegFunctionAG3d &f, int gorder ) {

 int orderc = getorderc();

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int i;

 GaussRuleTriangle gr(gorder);
 gr.pool = (double*)alloca(sizeof(double)*5*gr.ngauss);
 gr.init(faceindex);

 double *NN = (double*)alloca(sizeof(double)*2*4*gr.ngauss*order);
 for(i=1;i<=order;i++)
   lagGalShapeFunction(i,4*gr.ngauss,gr.xigauss,NN+(i-1)*2*4*gr.ngauss);

 double *N = (double*) alloca(sizeof(double)*4*orderc);
 int ng = gr.ngauss;
 double (*dNdx)[3] = (double(*)[3])alloca(sizeof(double)*orderc*3);
 for(i=0;i<ng;i++) {

   int r,s,t;
   int c=0;
   for(t=0;t<order;t++) {
     for(s=0;s<order-t;s++) {
       for(r=0;r<order-t-s;r++) {
         int u = order - 1 - r - s - t;
         int rr = r*2*4*gr.ngauss+0*gr.ngauss+i;
         int ss = s*2*4*gr.ngauss+1*gr.ngauss+i;
         int tt = t*2*4*gr.ngauss+2*gr.ngauss+i;
         int uu = u*2*4*gr.ngauss+3*gr.ngauss+i;
         int rrd = rr+4*gr.ngauss;
         int ssd = ss+4*gr.ngauss;
         int ttd = tt+4*gr.ngauss;
         int uud = uu+4*gr.ngauss;
         double scaling = 1.0;
         N[c] = scaling*NN[rr]*NN[ss]*NN[tt]*NN[uu];
         N[1*orderc+c] = scaling*(NN[rrd]*NN[ss]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[2*orderc+c] = scaling*(NN[rr]*NN[ssd]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[3*orderc+c] = scaling*(NN[rr]*NN[ss]*NN[ttd]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         c++;
       }
     }
   }

   double w = gr.wgauss[i];

   double J[3][3];
   jmatrix(N,xyz,J);

   double cross[3];
   crossj(J,faceindex,cross);
  
   double nsign;
   if (cross[0]*(scxyz[0]-cxyz[0])+
       cross[1]*(scxyz[1]-cxyz[1])+
       cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
   else nsign = -1;

   double det = detj(J);
   if (det<0.0) fprintf(stderr,"Negative vol Jacobian: %f %f %f %f %f\n",
                        det,gr.xigauss[0*gr.ngauss+i],
                        gr.xigauss[1*gr.ngauss+i],
                        gr.xigauss[2*gr.ngauss+i],
                        gr.xigauss[3*gr.ngauss+i]);

   double j[3][3];
   invj(J,det,j);
   int m;
   for(m=0;m<orderc;m++) {
     dNdx[m][0] =
      j[0][0]*N[orderc+m]+j[0][1]*N[2*orderc+m]+j[0][2]*N[3*orderc+m];
     dNdx[m][1] =
      j[1][0]*N[orderc+m]+j[1][1]*N[2*orderc+m]+j[1][2]*N[3*orderc+m];
     dNdx[m][2] =
      j[2][0]*N[orderc+m]+j[2][1]*N[2*orderc+m]+j[2][2]*N[3*orderc+m];
   }


   double x[3] = {0.0,0.0,0.0};
   for(m=0;m<orderc;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+orderc];
     x[2] += N[m]*xyz[m+2*orderc];
   }

   f.evaluate(x,N,dNdx,cross,nsign,w);

 }
}


void IsoParamUtilsTetra::surfSurfInt3d(double *xyz, IntegFunctionA3d &f, int gorder ) {

 int ordersq = getordersq();

 int i;

 GaussRuleTriangle gr(gorder);
 gr.pool = (double*)alloca(sizeof(double)*4*gr.ngauss);
 gr.init(0);

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

   double J[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
   int m;
   for(m=0;m<ordersq;m++) {
     J[0][0] += N[ordersq+m]*xyz[m];
     J[0][1] += N[ordersq+m]*xyz[m+ordersq];
     J[0][2] += N[ordersq+m]*xyz[m+2*ordersq];
     J[1][0] += N[2*ordersq+m]*xyz[m];
     J[1][1] += N[2*ordersq+m]*xyz[m+ordersq];
     J[1][2] += N[2*ordersq+m]*xyz[m+2*ordersq];
   }

   double cross[3];
   cross[0] = J[0][1]*J[1][2]-J[0][2]*J[1][1];
   cross[1] = J[0][2]*J[1][0]-J[0][0]*J[1][2];
   cross[2] = J[0][0]*J[1][1]-J[0][1]*J[1][0];
  
   double x[3] = {0.0,0.0,0.0};
   for(m=0;m<ordersq;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+ordersq];
     x[2] += N[m]*xyz[m+2*ordersq];
   }

   double nsign = 1.0; 
   f.evaluate(x,N,cross,nsign,w); 
 }
}


void IsoParamUtilsTetra::surfCurvInt3d(double *xyz, int faceindex,
                                       IntegFunctionAC3d &f, int gorder ) {

 int orderc = getorderc();

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 int i;

 GaussRuleTriangle gr(gorder);
 gr.pool = (double*)dbg_alloca(sizeof(double)*5*gr.ngauss);
 gr.init(faceindex);

 double *NN = (double*)dbg_alloca(sizeof(double)*3*4*gr.ngauss*order);
 for(i=1;i<=order;i++)
   lagGalShapeFunction(i,4*gr.ngauss,gr.xigauss,
     NN+(i-1)*3*4*gr.ngauss,1);

 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
 double *N2 = (double*) dbg_alloca(sizeof(double)*6*orderc);
 double (*surfgrad)[2] = (double(*)[2])dbg_alloca(sizeof(double)*2*orderc);
 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {

   int r,s,t;
   int c=0;
   for(t=0;t<order;t++) {
     for(s=0;s<order-t;s++) {
       for(r=0;r<order-t-s;r++) {
         int u = order - 1 - r - s - t;
         int rr = r*3*4*gr.ngauss+0*gr.ngauss+i;
         int ss = s*3*4*gr.ngauss+1*gr.ngauss+i;
         int tt = t*3*4*gr.ngauss+2*gr.ngauss+i;
         int uu = u*3*4*gr.ngauss+3*gr.ngauss+i;
         int rrd = rr+4*gr.ngauss;
         int ssd = ss+4*gr.ngauss;
         int ttd = tt+4*gr.ngauss;
         int uud = uu+4*gr.ngauss;
         int rrdd = rrd+4*gr.ngauss;
         int ssdd = ssd+4*gr.ngauss;
         int ttdd = ttd+4*gr.ngauss;
         int uudd = uud+4*gr.ngauss;
         N[c] = NN[rr]*NN[ss]*NN[tt]*NN[uu];
         N[1*orderc+c] = (NN[rrd]*NN[ss]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[2*orderc+c] = (NN[rr]*NN[ssd]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[3*orderc+c] = (NN[rr]*NN[ss]*NN[ttd]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N2[0*orderc+c] = (NN[rrdd]*NN[ss]*NN[tt]*NN[uu]-
                       2.0*NN[rrd]*NN[ss]*NN[tt]*NN[uud]+
                           NN[rr]*NN[ss]*NN[tt]*NN[uudd]);
         N2[1*orderc+c] = (NN[rr]*NN[ssdd]*NN[tt]*NN[uu]-
                       2.0*NN[rr]*NN[ssd]*NN[tt]*NN[uud]+
                           NN[rr]*NN[ss]*NN[tt]*NN[uudd]);
         N2[2*orderc+c] = (NN[rr]*NN[ss]*NN[ttdd]*NN[uu]-
                       2.0*NN[rr]*NN[ss]*NN[ttd]*NN[uud]+
                           NN[rr]*NN[ss]*NN[tt]*NN[uudd]);
         N2[3*orderc+c] = (NN[rrd]*NN[ssd]*NN[tt]*NN[uu]-
                           NN[rr]*NN[ssd]*NN[tt]*NN[uud]-
                           NN[rrd]*NN[ss]*NN[tt]*NN[uud]+
                           NN[rr]*NN[ss]*NN[tt]*NN[uudd]);
         N2[4*orderc+c] = (NN[rrd]*NN[ss]*NN[ttd]*NN[uu]-
                           NN[rr]*NN[ss]*NN[ttd]*NN[uud]-
                           NN[rrd]*NN[ss]*NN[tt]*NN[uud]+
                           NN[rr]*NN[ss]*NN[tt]*NN[uudd]);
         N2[5*orderc+c] = (NN[rr]*NN[ssd]*NN[ttd]*NN[uu]-
                           NN[rr]*NN[ss]*NN[ttd]*NN[uud]-
                           NN[rr]*NN[ssd]*NN[tt]*NN[uud]+
                           NN[rr]*NN[ss]*NN[tt]*NN[uudd]);
         c++;
       }
     }
   }

   double w = gr.wgauss[i];

   double J[3][3];
   jmatrix(N,xyz,J);

   double cross[3];
   crossj(J,faceindex,cross);
  
   double nsign;
   if (cross[0]*(scxyz[0]-cxyz[0])+
       cross[1]*(scxyz[1]-cxyz[1])+
       cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
   else nsign = -1;

   double x[3] = {0.0,0.0,0.0};
   int m;
   for(m=0;m<orderc;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+orderc];
     x[2] += N[m]*xyz[m+2*orderc];
   }

   double sd[5][3];
   getSurfDer(N,N2,xyz,faceindex,sd,surfgrad);
    
   f.evaluate(x,N,cross,nsign, sd, surfgrad, w);
 }
}


void IsoParamUtils::volumeInt3d(double *xyz, IntegFunctionV3d &f, int gorder ) {

 int orderc = getorderc();

 GaussRuleLine gr(gorder,GSUBDIV);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double *NN = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss*order);
 lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
 double (*dNdx)[3] = (double(*)[3])dbg_alloca(sizeof(double)*orderc*3);
 int ng = gr.ngauss;
 int i;
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {
     int k;
     for(k=0;k<ng;k++) {

       int mx,my,mz;
       for(mx=0;mx<order;mx++) {
         for(my=0;my<order;my++) {
            for(mz=0;mz<order;mz++) {
              N[mz*order*order+my*order+mx] =
                NN[i*order+mx]*NN[j*order+my]*NN[k*order+mz];
              N[orderc+mz*order*order+my*order+mx] =
                NN[ng*order+i*order+mx]*NN[j*order+my]*NN[k*order+mz];
              N[2*orderc+mz*order*order+my*order+mx] =
                NN[i*order+mx]*NN[ng*order+j*order+my]*NN[k*order+mz];
              N[3*orderc+mz*order*order+my*order+mx] =
                NN[i*order+mx]*NN[j*order+my]*NN[ng*order+k*order+mz];
            }
          }
        }
           
       double w = gr.wgauss[i]*gr.wgauss[j]*gr.wgauss[k];

       double J[3][3];
       jmatrix(N,xyz,J);
       double det = detj(J); 
       if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f %f\n",
                            det,xyz[0],xyz[orderc+0],xyz[2*orderc+0]);

       double j[3][3];
       invj(J,det,j);
       int m;
       for(m=0;m<orderc;m++) {
         dNdx[m][0] = 
          j[0][0]*N[orderc+m]+j[0][1]*N[2*orderc+m]+j[0][2]*N[3*orderc+m];
         dNdx[m][1] = 
          j[1][0]*N[orderc+m]+j[1][1]*N[2*orderc+m]+j[1][2]*N[3*orderc+m];
         dNdx[m][2] = 
          j[2][0]*N[orderc+m]+j[2][1]*N[2*orderc+m]+j[2][2]*N[3*orderc+m];
       }

       double x[3] = {0.0,0.0,0.0};
       for(m=0;m<orderc;m++) {
         x[0] += N[m]*xyz[m];
         x[1] += N[m]*xyz[m+orderc];
         x[2] += N[m]*xyz[m+2*orderc];
       }
       f.evaluate(x,N,dNdx,w,det);
     }
   }
 }
}


void IsoParamUtils::spectralVolumeInt3d(double *xyz, IntegFunctionV3d &f ) {

 int orderc = getorderc();

 GaussLobattoRuleLine gr(order);

 double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*gr.ngauss);
 spectralLagGalShapeFunction(gr.ngauss,gr.xigauss,NN);

 double *N = (double*) alloca(sizeof(double)*4*orderc);
 double (*dNdx)[3] = (double(*)[3])alloca(sizeof(double)*orderc*3);
 int ng = gr.ngauss;
 int i;
 for(i=0;i<ng;i++) {
   int j;
   for(j=0;j<ng;j++) {
     int k;
     for(k=0;k<ng;k++) {

       int mx,my,mz;
       for(mx=0;mx<order;mx++) {
         for(my=0;my<order;my++) {
            for(mz=0;mz<order;mz++) {
              N[mz*order*order+my*order+mx] =
                NN[i*order+mx]*NN[j*order+my]*NN[k*order+mz];
              N[orderc+mz*order*order+my*order+mx] =
                NN[ng*order+i*order+mx]*NN[j*order+my]*NN[k*order+mz];
              N[2*orderc+mz*order*order+my*order+mx] =
                NN[i*order+mx]*NN[ng*order+j*order+my]*NN[k*order+mz];
              N[3*orderc+mz*order*order+my*order+mx] =
                NN[i*order+mx]*NN[j*order+my]*NN[ng*order+k*order+mz];
            }
          }
        }
           
       double w = gr.wgauss[i]*gr.wgauss[j]*gr.wgauss[k];

       double J[3][3];
       jmatrix(N,xyz,J);
       double det = detj(J); 
       if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f %f\n",
                            det,xyz[0],xyz[orderc+0],xyz[2*orderc+0]);

       double j[3][3];
       invj(J,det,j);
       int m;
       for(m=0;m<orderc;m++) {
         dNdx[m][0] = 
          j[0][0]*N[orderc+m]+j[0][1]*N[2*orderc+m]+j[0][2]*N[3*orderc+m];
         dNdx[m][1] = 
          j[1][0]*N[orderc+m]+j[1][1]*N[2*orderc+m]+j[1][2]*N[3*orderc+m];
         dNdx[m][2] = 
          j[2][0]*N[orderc+m]+j[2][1]*N[2*orderc+m]+j[2][2]*N[3*orderc+m];
       }

       double x[3] = {0.0,0.0,0.0};
       for(m=0;m<orderc;m++) {
         x[0] += N[m]*xyz[m];
         x[1] += N[m]*xyz[m+orderc];
         x[2] += N[m]*xyz[m+2*orderc];
       }
       f.evaluate(x,N,dNdx,w,det);
     }
   }
 }
}


void IsoParamUtilsTetra::volumeInt3d(double *xyz, IntegFunctionV3d &f, 
                                     int gorder ) {

 int orderc = getorderc();

 GaussRuleTetra gr(gorder);
 gr.pool = (double*)dbg_alloca(sizeof(double)*5*gr.ngauss);
 gr.init();

 double *NN = (double*)dbg_alloca(sizeof(double)*2*4*gr.ngauss*order);
 int i;
 for(i=1;i<=order;i++)
   lagGalShapeFunction(i,4*gr.ngauss,gr.xigauss,NN+(i-1)*2*4*gr.ngauss);

 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
 double (*dNdx)[3] = (double(*)[3])dbg_alloca(sizeof(double)*orderc*3);
 int ng = gr.ngauss;
 for(i=0;i<ng;i++) {

   int r,s,t;
   int c=0;
   for(t=0;t<order;t++) {
     for(s=0;s<order-t;s++) {
       for(r=0;r<order-t-s;r++) {
         int u = order - 1 - r - s - t;
         int rr = r*2*4*gr.ngauss+0*gr.ngauss+i;
         int ss = s*2*4*gr.ngauss+1*gr.ngauss+i;
         int tt = t*2*4*gr.ngauss+2*gr.ngauss+i;
         int uu = u*2*4*gr.ngauss+3*gr.ngauss+i;
         int rrd = rr+4*gr.ngauss;
         int ssd = ss+4*gr.ngauss;
         int ttd = tt+4*gr.ngauss;
         int uud = uu+4*gr.ngauss;
         double scaling = 1.0;
         N[c] = scaling*NN[rr]*NN[ss]*NN[tt]*NN[uu];
         N[1*orderc+c] = scaling*(NN[rrd]*NN[ss]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[2*orderc+c] = scaling*(NN[rr]*NN[ssd]*NN[tt]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         N[3*orderc+c] = scaling*(NN[rr]*NN[ss]*NN[ttd]*NN[uu]-
                         NN[rr]*NN[ss]*NN[tt]*NN[uud]);
         c++;
       }
     }
   }

   double J[3][3];
   jmatrix(N,xyz,J);
   double det = detj(J); 
   if (det<0.0) fprintf(stderr,"Negative vol Jacobian: %f %f %f %f %f\n",
                        det,gr.xigauss[0*gr.ngauss+i],
                        gr.xigauss[1*gr.ngauss+i],
                        gr.xigauss[2*gr.ngauss+i],
                        gr.xigauss[3*gr.ngauss+i]);

   double j[3][3];
   invj(J,det,j);
   int m;
   for(m=0;m<orderc;m++) {
     dNdx[m][0] = 
      j[0][0]*N[orderc+m]+j[0][1]*N[2*orderc+m]+j[0][2]*N[3*orderc+m];
     dNdx[m][1] = 
      j[1][0]*N[orderc+m]+j[1][1]*N[2*orderc+m]+j[1][2]*N[3*orderc+m];
     dNdx[m][2] = 
      j[2][0]*N[orderc+m]+j[2][1]*N[2*orderc+m]+j[2][2]*N[3*orderc+m];
   }

   double x[3] = {0.0,0.0,0.0};
   for(m=0;m<orderc;m++) {
     x[0] += N[m]*xyz[m];
     x[1] += N[m]*xyz[m+orderc];
     x[2] += N[m]*xyz[m+2*orderc];
   }
   f.evaluate(x,N,dNdx,gr.wgauss[i],det);
 }
}



void IsoParamUtilsPrism::surfInt3d(double *xyz, int faceindex,
                                   IntegFunctionA3d &f, int gorder ) {

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);

 if (faceindex<=2) {

 GaussRuleTriangle gr(13);
 gr.pool = (double*)alloca(sizeof(double)*4*gr.ngauss);
 gr.init(0);

 double *NN = (double*)alloca(sizeof(double)*2*3*gr.ngauss*order);
 for(int i=1;i<=order;i++)
   lagGalShapeFunctionTr(i,3*gr.ngauss,gr.xigauss,NN+(i-1)*2*3*gr.ngauss);

 double *NN2 = (double*)dbg_alloca(sizeof(double)*2*1*order);
 double p = (faceindex==1)?-1.0:1.0;
 lagGalShapeFunction(1,&p,NN2);

 int orderc= getorderc();
 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);

 for(int i=0;i<1;i++) {
   for(int j=0;j<gr.ngauss;j++) {

     int c=0;
     for(int u=0;u<order;u++) {
       for(int s=0;s<order;s++) {
         for(int r=0;r<order-s;r++) {
           int t = order - 1 - r - s;
           int rr = r*2*3*gr.ngauss+0*gr.ngauss+j;
           int ss = s*2*3*gr.ngauss+1*gr.ngauss+j;
           int tt = t*2*3*gr.ngauss+2*gr.ngauss+j;
           int rrd = rr+3*gr.ngauss;
           int ssd = ss+3*gr.ngauss;
           int ttd = tt+3*gr.ngauss;
           N[c] = NN[rr]*NN[ss]*NN[tt]*NN2[i*order+u];
           N[1*orderc+c] = (NN[rrd]*NN[ss]*NN[tt]-
                           NN[rr]*NN[ss]*NN[ttd])*NN2[i*order+u];
           N[2*orderc+c] = (NN[rr]*NN[ssd]*NN[tt]-
                           NN[rr]*NN[ss]*NN[ttd])*NN2[i*order+u];
           N[3*orderc+c] = NN[rr]*NN[ss]*NN[tt]*NN2[1*order+i*order+u];
           c++;
         }
       }
     }

     double w = gr.wgauss[j];

     double J[3][3];
     jmatrix(N,xyz,J);

    double det = detj(J); 
   if (det<0.0) fprintf(stderr,"Negative vol Jacobian %d\n",faceindex);

     double cross[3];
     crossj(J,faceindex,cross);
  
     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     double x[3] = {0.0,0.0,0.0};
     for(int m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }

     f.evaluate(x,N,cross,nsign,w); 
   }
 }

 } else {

 GaussRuleLine gr2(gorder,GSUBDIV);
 gr2.pool = (double*)dbg_alloca(sizeof(double)*2*gr2.ngauss);
 gr2.init();

 double *NN2 = (double*)dbg_alloca(sizeof(double)*2*gr2.ngauss*order);
 lagGalShapeFunction(gr2.ngauss,gr2.xigauss,NN2);

 int ng = gr2.ngauss;
 double *xig = new double[3*ng];
 double *wg = new double[ng];
 for(int i=0;i<3*ng;i++) xig[i] = 0.0;
 if (faceindex==3) {
   for(int i=0;i<ng;i++) wg[i] = gr2.wgauss[i]/(2.0);
   for(int i=0;i<ng;i++) xig[0*ng+i] = (gr2.xigauss[i]+1.0)/2.0;
   for(int i=0;i<ng;i++) xig[1*ng+i] = 1.0-xig[0*ng+i];
   for(int i=0;i<ng;i++) xig[2*ng+i] = 1.0 - xig[0*ng+i] - xig[1*ng+i];
 } else if (faceindex==4) {
   for(int i=0;i<ng;i++) wg[i] = gr2.wgauss[i]/2.0;
   for(int i=0;i<ng;i++) xig[1*ng+i] = (gr2.xigauss[i]+1.0)/2.0;
   for(int i=0;i<ng;i++) xig[2*ng+i] = 1.0 - xig[1*ng+i];
 } else {
   for(int i=0;i<ng;i++) wg[i] = gr2.wgauss[i]/2.0;
   for(int i=0;i<ng;i++) xig[0*ng+i] = (gr2.xigauss[i]+1.0)/2.0;
   for(int i=0;i<ng;i++) xig[2*ng+i] = 1.0 - xig[0*ng+i];
 }

 double *NN = (double*)alloca(sizeof(double)*2*3*ng*order);
 int i;
 for(i=1;i<=order;i++)
   lagGalShapeFunctionTr(i,3*ng,xig,NN+(i-1)*2*3*ng);

 int orderc= getorderc();
 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);

 for(int i=0;i<gr2.ngauss;i++) {
   for(int j=0;j<ng;j++) {

     int c=0;
     for(int u=0;u<order;u++) {
       for(int s=0;s<order;s++) {
         for(int r=0;r<order-s;r++) {
           int t = order - 1 - r - s;
           int rr = r*2*3*ng+0*ng+j;
           int ss = s*2*3*ng+1*ng+j;
           int tt = t*2*3*ng+2*ng+j;
           int rrd = rr+3*ng;
           int ssd = ss+3*ng;
           int ttd = tt+3*ng;
           N[c] = NN[rr]*NN[ss]*NN[tt]*NN2[i*order+u];
           N[1*orderc+c] = (NN[rrd]*NN[ss]*NN[tt]-
                           NN[rr]*NN[ss]*NN[ttd])*NN2[i*order+u];
           N[2*orderc+c] = (NN[rr]*NN[ssd]*NN[tt]-
                           NN[rr]*NN[ss]*NN[ttd])*NN2[i*order+u];
           N[3*orderc+c] = NN[rr]*NN[ss]*NN[tt]*NN2[gr2.ngauss*order+i*order+u];
           c++;
         }
       }
     }

     double w = gr2.wgauss[i]*wg[j];

     double J[3][3];
     jmatrix(N,xyz,J);
    double det = detj(J); 
   if (det<0.0) fprintf(stderr,"Negative vol Jacobian %d\n",faceindex);

     double cross[3];
     crossj(J,faceindex,cross);
  
     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     double x[3] = {0.0,0.0,0.0};
     for(int m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }

     f.evaluate(x,N,cross,nsign,w); 
   }
 }
 delete[] xig;
 delete[] wg;
 }
}

void IsoParamUtilsPrism::volumeInt3d(double *xyz, IntegFunctionV3d &f, 
                                     int gorder ) {

 GaussRuleTriangle gr(13);
 gr.pool = (double*)alloca(sizeof(double)*4*gr.ngauss);
 gr.init(0);

 double *NN = (double*)alloca(sizeof(double)*2*3*gr.ngauss*order);
 int i;
 for(i=1;i<=order;i++)
   lagGalShapeFunctionTr(i,3*gr.ngauss,gr.xigauss,NN+(i-1)*2*3*gr.ngauss);

 GaussRuleLine gr2(gorder,GSUBDIV);
 gr2.pool = (double*)dbg_alloca(sizeof(double)*2*gr2.ngauss);
 gr2.init();

 double *NN2 = (double*)dbg_alloca(sizeof(double)*2*gr2.ngauss*order);
 lagGalShapeFunction(gr2.ngauss,gr2.xigauss,NN2);

 int orderc= getorderc();
 double *N = (double*) dbg_alloca(sizeof(double)*4*orderc);
 double (*dNdx)[3] = (double(*)[3])dbg_alloca(sizeof(double)*orderc*3);

 for(int i=0;i<gr2.ngauss;i++) {
   for(int j=0;j<gr.ngauss;j++) {

     int c=0;
     for(int u=0;u<order;u++) {
       for(int s=0;s<order;s++) {
         for(int r=0;r<order-s;r++) {
           int t = order - 1 - r - s;
           int rr = r*2*3*gr.ngauss+0*gr.ngauss+i;
           int ss = s*2*3*gr.ngauss+1*gr.ngauss+i;
           int tt = t*2*3*gr.ngauss+2*gr.ngauss+i;
           int rrd = rr+3*gr.ngauss;
           int ssd = ss+3*gr.ngauss;
           int ttd = tt+3*gr.ngauss;
           N[c] = NN[rr]*NN[ss]*NN[tt]*NN2[i*order+u];
           N[1*orderc+c] = (NN[rrd]*NN[ss]*NN[tt]-
                           NN[rr]*NN[ss]*NN[ttd])*NN2[i*order+u];
           N[2*orderc+c] = (NN[rr]*NN[ssd]*NN[tt]-
                           NN[rr]*NN[ss]*NN[ttd])*NN2[i*order+u];
           N[3*orderc+c] = NN[rr]*NN[ss]*NN[tt]*NN2[gr2.ngauss*order+i*order+u];
           c++;
         }
       }
     }

     double w = gr2.wgauss[i]*gr.wgauss[j];

     double J[3][3];
     jmatrix(N,xyz,J);
     double det = detj(J);
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f %f\n",
                          det,xyz[0],xyz[orderc+0],xyz[2*orderc+0]);

     double ij[3][3];
     invj(J,det,ij);

     for(int m=0;m<orderc;m++) {
       dNdx[m][0] = 
        ij[0][0]*N[orderc+m]+ij[0][1]*N[2*orderc+m]+ij[0][2]*N[3*orderc+m];
       dNdx[m][1] = 
        ij[1][0]*N[orderc+m]+ij[1][1]*N[2*orderc+m]+ij[1][2]*N[3*orderc+m];
       dNdx[m][2] = 
        ij[2][0]*N[orderc+m]+ij[2][1]*N[2*orderc+m]+ij[2][2]*N[3*orderc+m];
     }
  
     double x[3] = {0.0,0.0,0.0};
     for(int m=0;m<orderc;m++) {
       x[0] += N[m]*xyz[m];
       x[1] += N[m]*xyz[m+orderc];
       x[2] += N[m]*xyz[m+2*orderc];
     }
     f.evaluate(x,N,dNdx,w,det);
   }
 }
}


void IsoParamUtilsPyramid::surfInt3d(double *xyz, int faceindex,
                                   IntegFunctionA3d &f, int gorder ) {

 double cxyz[3];
 elementcenter(xyz,cxyz);
 double scxyz[3];
 sidecenter(xyz,faceindex,scxyz);
 int ordersq = getordersq(faceindex);
 int orderc = getorderc();
 int *fi = new int[ordersq];
 faceindeces(faceindex,fi);

 if (faceindex>1) {
 
// Adapted from Tetra:surfSurf 
   GaussRuleTriangle gr(13);
   gr.pool = (double*)alloca(sizeof(double)*4*gr.ngauss);
   gr.init(0);
  
   double *NN = (double*)alloca(sizeof(double)*2*3*gr.ngauss*order);
   for(int i=1;i<=order;i++)
     lagGalShapeFunctionTr(i,3*gr.ngauss,gr.xigauss,NN+(i-1)*2*3*gr.ngauss);
  
   double *N = (double*) alloca(sizeof(double)*3*orderc);
   int ng = gr.ngauss;
   for(int i=0;i<ng;i++) {
  
     for(int j=0;j<3*orderc;j++) N[j] = 0.0;
     int c=0;
     for(int s=0;s<order;s++) {
       for(int r=0;r<order-s;r++) {
         int u = order - 1 - r - s;
         int rr = r*2*3*gr.ngauss+0*gr.ngauss+i;
         int ss = s*2*3*gr.ngauss+1*gr.ngauss+i;
         int uu = u*2*3*gr.ngauss+2*gr.ngauss+i;
         int rrd = rr+3*gr.ngauss;
         int ssd = ss+3*gr.ngauss;
         int uud = uu+3*gr.ngauss;
         N[fi[c]] = NN[rr]*NN[ss]*NN[uu];
         N[1*orderc+fi[c]] = (NN[rrd]*NN[ss]*NN[uu]-
                         NN[rr]*NN[ss]*NN[uud]);
         N[2*orderc+fi[c]] = (NN[rr]*NN[ssd]*NN[uu]-
                         NN[rr]*NN[ss]*NN[uud]);
         c++;
       }
     }
  
     double w = gr.wgauss[i];
  
     double J[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
     for(int m=0;m<ordersq;m++) {
       int mm = fi[m];
       J[0][0] += N[orderc+mm]*xyz[mm];
       J[0][1] += N[orderc+mm]*xyz[mm+orderc];
       J[0][2] += N[orderc+mm]*xyz[mm+2*orderc];
       J[1][0] += N[2*orderc+mm]*xyz[mm];
       J[1][1] += N[2*orderc+mm]*xyz[mm+orderc];
       J[1][2] += N[2*orderc+mm]*xyz[mm+2*orderc];
     }
  
     double cross[3];
     cross[0] = J[0][1]*J[1][2]-J[0][2]*J[1][1];
     cross[1] = J[0][2]*J[1][0]-J[0][0]*J[1][2];
     cross[2] = J[0][0]*J[1][1]-J[0][1]*J[1][0];
    
     double x[3] = {0.0,0.0,0.0};
     for(int m=0;m<ordersq;m++) {
       int mm = fi[m];
       x[0] += N[mm]*xyz[mm];
       x[1] += N[mm]*xyz[mm+orderc];
       x[2] += N[mm]*xyz[mm+2*orderc];
     }

     double nsign;
     if (cross[0]*(scxyz[0]-cxyz[0])+
         cross[1]*(scxyz[1]-cxyz[1])+
         cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
     else nsign = -1;

     f.evaluate(x,N,cross,nsign,w); 
   }

 } else {
// Adapted from Hexa:surfSurf 
  
   GaussRuleLine gr(gorder,GSUBDIV);
   gr.pool = (double*)alloca(sizeof(double)*2*gr.ngauss);
   gr.init();
  
   double *NN = (double*)alloca(sizeof(double)*2*gr.ngauss*order);
   lagGalShapeFunction(gr.ngauss,gr.xigauss,NN);
  
   double *N = (double*) alloca(sizeof(double)*3*orderc);
   int ng = gr.ngauss;
   for(int i=0;i<ng;i++) {
     for(int j=0;j<ng;j++) {

       for(int k=0;k<3*orderc;k++) N[k] = 0.0;
       for(int mx=0;mx<order;mx++) {
         for(int my=0;my<order;my++) {
              N[fi[my*order+mx]] =
                NN[i*order+mx]*NN[j*order+my];
              N[orderc+fi[my*order+mx]] =
                NN[ng*order+i*order+mx]*NN[j*order+my];
              N[2*orderc+fi[my*order+mx]] =
                NN[i*order+mx]*NN[ng*order+j*order+my];
         }
       }
             
       double w = gr.wgauss[i]*gr.wgauss[j];
  
       double J[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
       for(int m=0;m<ordersq;m++) {
         int mm = fi[m];
         J[0][0] += N[orderc+mm]*xyz[mm];
         J[0][1] += N[orderc+mm]*xyz[mm+orderc];
         J[0][2] += N[orderc+mm]*xyz[mm+2*orderc];
         J[1][0] += N[2*orderc+mm]*xyz[mm];
         J[1][1] += N[2*orderc+mm]*xyz[mm+orderc];
         J[1][2] += N[2*orderc+mm]*xyz[mm+2*orderc];
       }
  
       double cross[3];
       cross[0] = J[0][1]*J[1][2]-J[0][2]*J[1][1];
       cross[1] = J[0][2]*J[1][0]-J[0][0]*J[1][2];
       cross[2] = J[0][0]*J[1][1]-J[0][1]*J[1][0];
  
       double x[3] = {0.0,0.0,0.0};
       for(int m=0;m<ordersq;m++) {
         int mm = fi[m];
         x[0] += N[mm]*xyz[mm];
         x[1] += N[mm]*xyz[mm+orderc];
         x[2] += N[mm]*xyz[mm+2*orderc];
       }
     
       double nsign;
       if (cross[0]*(scxyz[0]-cxyz[0])+
           cross[1]*(scxyz[1]-cxyz[1])+
           cross[2]*(scxyz[2]-cxyz[2])>0.0) nsign = 1;
       else nsign = -1;

       f.evaluate(x,N,cross,nsign,w);
     }
   }

 }
 delete[] fi;
}

#ifndef _TEMPLATE_FIX_
#include <Element.d/Helm.d/IsoParamUtilsT.C>
#endif
