#include <Element.d/Helm.d/GaussRules.h>
#include <cstdio>
#include <cmath>

#include <Utils.d/dbg_alloca.h>

GaussLobattoRuleLine::GaussLobattoRuleLine(int _ngauss) {
 ngauss = _ngauss;
 xigauss = new double[ngauss];
 wgauss = new double[ngauss];

 // for the exact values, see Higher-Order Numerical Methods for Transient Wave Equations
 //  by Gary Cohen, Springer 2002 (available at Physics library in Stanford
 //  , ref: QA927.C54 2002 PHYS) pp 175-176.

 if (ngauss==2) {
   xigauss[0] = -1.0;
   xigauss[1] = 1.0;
   wgauss[0] = 1.0;
   wgauss[1] = 1.0;
 } else if (ngauss==3) {
   xigauss[0] = -1.0;
   xigauss[1] = 0.0;
   xigauss[2] = 1.0;
   wgauss[0] = 0.3333333333333333333333;
   wgauss[1] = 1.3333333333333333333333;
   wgauss[2] = 0.3333333333333333333333;
 } else if (ngauss==4) {
   xigauss[0] = -1.0;
   xigauss[1] = -0.447213596;
   xigauss[2] = 0.447213596;
   xigauss[3] = 1.0;
   wgauss[0] = 0.16666666666666;
   wgauss[1] = 0.83333333333333;
   wgauss[2] = 0.83333333333333;
   wgauss[3] = 0.16666666666666;
 } else if (ngauss==5) {
   xigauss[0] = -1.0;
   xigauss[1] = -0.65465367;
   xigauss[2] = 0.0;
   xigauss[3] = 0.65465367;
   xigauss[4] = 1.0;
   wgauss[0] = 0.1;
   wgauss[1] = 0.54444444444444;
   wgauss[2] = 0.71111111111111;
   wgauss[3] = 0.54444444444444;
   wgauss[4] = 0.1;
 } else {
   fprintf(stderr,"Element.d/Helm.d/GaussRules.C: this number of gauss points has not been implemented\n");
 }
}

GaussLobattoRuleLine::~GaussLobattoRuleLine() {
 delete[] xigauss;
 delete[] wgauss;
}

GaussRuleLine::GaussRuleLine(int o ,int _nsubdiv) {
 if (o<0) {
 } else {
   ng = o;
 }
 pool = 0;
 nsubdiv = _nsubdiv;
 ngauss = ng*nsubdiv;
}


GaussRuleTriLine::GaussRuleTriLine(int o ,int _nsubdiv):
  GaussRuleLine(o,_nsubdiv) {
}


void GaussRuleLine::init() {
 wgauss = pool;
 xigauss = pool+ngauss;
 double *_xigauss, *_wgauss;
 _xigauss = (double*)dbg_alloca(sizeof(double)*ng);
 _wgauss = (double*)dbg_alloca(sizeof(double)*ng);
 if (ng==1) {
   _xigauss[0] = 0.0;
   _wgauss[0] = 2.0;
 } else if (ng==3) {
   _xigauss[0] = -0.7745966692414834;
   _xigauss[1] = 0.0;
   _xigauss[2] = 0.7745966692414834;
   _wgauss[0] = 5.0/9.0;
   _wgauss[1] = 8.0/9.0;
   _wgauss[2] = _wgauss[0];
 } else if  (ng==4) {
   _xigauss[0] = -0.861136311594053;
   _xigauss[1] = -0.339981043584856;
   _xigauss[2] = 0.339981043584856;
   _xigauss[3] = 0.861136311594053;
   _wgauss[0] = 0.347854845137454;
   _wgauss[1] = 0.652145154862546;
   _wgauss[2] = 0.652145154862546;
   _wgauss[3] = 0.347854845137454;
 } else if (ng==5) {
   _xigauss[0] = 0.0;
   _xigauss[1] =  0.538469310105683;
   _xigauss[2] = -0.538469310105683;
   _xigauss[3] =  0.906179845938664;   
   _xigauss[4] = -0.906179845938664;   
   _wgauss[0] = 0.568888888888889;
   _wgauss[1] = 0.478628670499366;
   _wgauss[2] = 0.478628670499366;
   _wgauss[3] = 0.236926885056189;
   _wgauss[4] = 0.236926885056189;

  
 } else if (ng==7) {
   _xigauss[0] = 0.949107912342758524526190;
   _xigauss[1] = -0.949107912342758524526190;
   _xigauss[2] = 0.741531185599394439863865;
   _xigauss[3] = -0.741531185599394439863865;
   _xigauss[4] = 0.405845151377397166906607;
   _xigauss[5] = -0.405845151377397166906607;
   _xigauss[6] = 0;
   _wgauss[0] = 0.12948496616886969327061;
   _wgauss[1] = 0.12948496616886969327061;
   _wgauss[2] = 0.27970539148927666790147;
   _wgauss[3] = 0.27970539148927666790147;
   _wgauss[4] = 0.38183005050511894495037;
   _wgauss[5] = 0.38183005050511894495037;
   _wgauss[6] = 0.41795918367346938775510;
 }
 int i;
 for(i=0;i<nsubdiv;i++) {
   int iw;
   for(iw=0;iw<ng;iw++) {
     wgauss[ng*i+iw] = _wgauss[iw]/double(nsubdiv);
     xigauss[ng*i+iw] = -1.0+double(2*i+1)/double(nsubdiv)+
                        _xigauss[iw]/double(nsubdiv);
   }
 }
}


void GaussRuleTriLine::init(int faceindex) {
 wgauss = pool;
 xigauss = pool+ngauss;

 double *_xigauss, *_wgauss,*__xigauss;
 _xigauss = (double*)alloca(sizeof(double)*ng);
 __xigauss = (double*)alloca(sizeof(double)*3*ngauss);
 _wgauss = (double*)alloca(sizeof(double)*ng);
 if (ng==1) {
   _xigauss[0] = 0.0;
   _wgauss[0] = 2.0;
 } else if (ng==3) {
   _xigauss[0] = -0.7745966692414834;
   _xigauss[1] = 0.0;
   _xigauss[2] = 0.7745966692414834;
   _wgauss[0] = 5.0/9.0;
   _wgauss[1] = 8.0/9.0;
   _wgauss[2] = _wgauss[0];
 } else if  (ng==4) {
   _xigauss[0] = -0.861136311594053;
   _xigauss[1] = -0.339981043584856;
   _xigauss[2] = 0.339981043584856;
   _xigauss[3] = 0.861136311594053;
   _wgauss[0] = 0.347854845137454;
   _wgauss[1] = 0.652145154862546;
   _wgauss[2] = 0.652145154862546;
   _wgauss[3] = 0.347854845137454;
 } else if (ng==7) {
   _xigauss[0] = 0.949107912342758524526190;
   _xigauss[1] = -0.949107912342758524526190;
   _xigauss[2] = 0.741531185599394439863865;
   _xigauss[3] = -0.741531185599394439863865;
   _xigauss[4] = 0.405845151377397166906607;
   _xigauss[5] = -0.405845151377397166906607;
   _xigauss[6] = 0;
   _wgauss[0] = 0.12948496616886969327061;
   _wgauss[1] = 0.12948496616886969327061;
   _wgauss[2] = 0.27970539148927666790147;
   _wgauss[3] = 0.27970539148927666790147;
   _wgauss[4] = 0.38183005050511894495037;
   _wgauss[5] = 0.38183005050511894495037;
   _wgauss[6] = 0.41795918367346938775510;
 }
 int i;
 for(i=0;i<nsubdiv;i++) {
   int iw;
   for(iw=0;iw<ng;iw++) {
     wgauss[ng*i+iw] = _wgauss[iw]/double(nsubdiv)/2.0;
//     if (faceindex==1) wgauss[ng*i+iw] *= sqrt(2.0);
     __xigauss[ng*i+iw] = -1.0+double(2*i+1)/double(nsubdiv)+
                        _xigauss[iw]/double(nsubdiv);
   }
 }
 createFaceInteg(faceindex,__xigauss);
}


void GaussRuleTriLine::createFaceInteg(int faceindex, double *_xigauss) {
 int i;
 if (faceindex==1) {
   for(i=0;i<ngauss;i++) {
     double xig = (_xigauss[i]+1.0)/2.0;
     xigauss[i] = xig;
     xigauss[ngauss+i] = 1.0-xig;
     xigauss[2*ngauss+i] = 0.0;
   }
 } else
 if (faceindex==3) {
   for(i=0;i<ngauss;i++) {
     double xig = (_xigauss[i]+1.0)/2.0;
     xigauss[i] = xig;
     xigauss[ngauss+i] = 0.0;
     xigauss[2*ngauss+i] = 1.0-xig;
   }
 } else
 if (faceindex==2) {
   for(i=0;i<ngauss;i++) {
     double xig = (_xigauss[i]+1.0)/2.0;
     xigauss[i] = 0.0;
     xigauss[ngauss+i] = xig;
     xigauss[2*ngauss+i] = 1.0-xig;
   }
 }
}


GaussRuleTetra::GaussRuleTetra(int o)  {
 if (o<0) {
 } else {
   ngauss = o;
 }
 pool = 0;
}


void GaussRuleTetra::init() {
 wgauss = pool;
 xigauss = pool+ngauss;
 if  (ngauss==1) {
   xigauss[0] = 0.3;
   xigauss[1] = 0.3;
   xigauss[2] = 0.3;
   xigauss[3] = 0.1;
 } else if  (ngauss==5) {
   xigauss[0] = 1.0/4.0;
   xigauss[1] = 1.0/3.0;
   xigauss[2] = 1.0/6.0;
   xigauss[3] = 1.0/6.0;
   xigauss[4] = 1.0/6.0;
   xigauss[5+0] = 1.0/4.0;
   xigauss[5+1] = 1.0/6.0;
   xigauss[5+2] = 1.0/3.0;
   xigauss[5+3] = 1.0/6.0;
   xigauss[5+4] = 1.0/6.0;
   xigauss[10+0] = 1.0/4.0;
   xigauss[10+1] = 1.0/6.0;
   xigauss[10+2] = 1.0/6.0;
   xigauss[10+3] = 1.0/3.0;
   xigauss[10+4] = 1.0/6.0;
   xigauss[15+0] = 1.0/4.0;
   xigauss[15+1] = 1.0/6.0;
   xigauss[15+2] = 1.0/6.0;
   xigauss[15+3] = 1.0/6.0;
   xigauss[15+4] = 1.0/3.0;
   wgauss[0] = -4.0/30.0;
   wgauss[1] = 3.0/40.0;
   wgauss[2] = 3.0/40.0;
   wgauss[3] = 3.0/40.0;
   wgauss[4] = 3.0/40.0;
 } else if (ngauss==3*3*3) {
   genFrom1D();
 } else if (ngauss==4*4*4) {
   genFrom1D();
 } else if (ngauss==7*7*7) {
   genFrom1D();
 }
}


void GaussRuleTetra::genFrom1D() {
 int ng;
 if (ngauss==3*3*3) {
   ng = 3;
 } else if (ngauss==4*4*4) {
   ng = 4;
 } else if (ngauss==7*7*7) {
   ng = 7;
 }

 GaussRuleLine gr(ng,1);
 gr.pool = (double*)dbg_alloca(2*ng*sizeof(double));
 gr.init();

 int i,j,k;
 for(i=0;i<ng;i++) for(j=0;j<ng;j++) for(k=0;k<ng;k++) {
   xigauss[0*ng*ng*ng+i*ng*ng+j*ng+k] = (1.0+gr.xigauss[k])/2.0;
   xigauss[1*ng*ng*ng+i*ng*ng+j*ng+k] = (1.0-gr.xigauss[k])/2.0*
                                        (1.0+gr.xigauss[j])/2.0;
   xigauss[2*ng*ng*ng+i*ng*ng+j*ng+k] = (1.0-gr.xigauss[k])/2.0*
                                        (1.0-gr.xigauss[j])/2.0*
                                        (1.0+gr.xigauss[i])/2.0;
   xigauss[3*ng*ng*ng+i*ng*ng+j*ng+k] = 1.0-xigauss[0*ng*ng*ng+i*ng*ng+j*ng+k]-
                                            xigauss[1*ng*ng*ng+i*ng*ng+j*ng+k]-
                                            xigauss[2*ng*ng*ng+i*ng*ng+j*ng+k];
   wgauss[i*ng*ng+j*ng+k] = (1.0-gr.xigauss[k])*
                            (1.0-gr.xigauss[k])*
                            (1.0-gr.xigauss[j])*
                            gr.wgauss[i]*gr.wgauss[j]*gr.wgauss[k]/64.0;
 }
}


GaussRuleTriangle::GaussRuleTriangle(int o)  {
 if (o<0) {
 } else {
   ngauss = o;
 }
 pool = 0;
}


void GaussRuleTriangle::createFaceInteg(int faceindex, double *_xigauss) {
 int i;
 if (faceindex==0) {
   for(i=0;i<3*ngauss;i++) xigauss[i] = _xigauss[i];
 } else 
 if (faceindex==1) {
   for(i=0;i<3*ngauss;i++) xigauss[i] = _xigauss[i];
   for(i=0;i<1*ngauss;i++) xigauss[3*ngauss+i] = 0.0;
 } else
 if (faceindex==2) {
   for(i=0;i<3*ngauss;i++) xigauss[1*ngauss+i] = _xigauss[i];
   for(i=0;i<1*ngauss;i++) xigauss[i] = 0.0;
 } else
 if (faceindex==3) {
   for(i=0;i<1*ngauss;i++) xigauss[i] = _xigauss[i];
   for(i=0;i<1*ngauss;i++) xigauss[1*ngauss+i] = 0.0;
   for(i=0;i<2*ngauss;i++) xigauss[2*ngauss+i] = _xigauss[1*ngauss+i];
 } else
 if (faceindex==4) {
   for(i=0;i<2*ngauss;i++) xigauss[i] = _xigauss[i];
   for(i=0;i<1*ngauss;i++) xigauss[2*ngauss+i] = 0.0;
   for(i=0;i<1*ngauss;i++) xigauss[3*ngauss+i] = _xigauss[2*ngauss+i];
 }
}


void GaussRuleTriangle::init(int faceindex) {
 wgauss = pool;
 xigauss = pool+ngauss;
 if (ngauss==1) {
   double _xigauss[3];
   _xigauss[0] =  0.333333333333333;
   _xigauss[1] =  0.333333333333333;
   _xigauss[2] =  0.333333333333333;

   createFaceInteg(faceindex,_xigauss);
   wgauss[0] =  1.000000000000000/2.0;
   
 } else if (ngauss==7) {
   double _xigauss[7*3];

   _xigauss[0] = 0.333333333333333; 
   _xigauss[1] = 0.797426985353087; 
   _xigauss[2] = 0.101286507323456; 
   _xigauss[3] = 0.101286507323456; 
   _xigauss[4] = 0.470142064105115; 
   _xigauss[5] = 0.470142064105115; 
   _xigauss[6] = 0.059715871789770; 

   _xigauss[7] = 0.333333333333333; 
   _xigauss[8] = 0.101286507323456; 
   _xigauss[9] = 0.797426985353087; 
   _xigauss[10] = 0.101286507323456; 
   _xigauss[11] = 0.470142064105115; 
   _xigauss[12] = 0.059715871789770; 
   _xigauss[13] = 0.470142064105115; 

   _xigauss[14] = 0.333333333333333; 
   _xigauss[15] = 0.101286507323456; 
   _xigauss[16] = 0.101286507323456; 
   _xigauss[17] = 0.797426985353087; 
   _xigauss[18] = 0.059715871789770; 
   _xigauss[19] = 0.470142064105115; 
   _xigauss[20] = 0.470142064105115; 

   createFaceInteg(faceindex,_xigauss);

   wgauss[0] = 0.225030000300000/2.0; 
   wgauss[1] = 0.125939180544827/2.0; 
   wgauss[2] = 0.125939180544827/2.0; 
   wgauss[3] = 0.125939180544827/2.0; 
   wgauss[4] = 0.132394152788506/2.0; 
   wgauss[5] = 0.132394152788506/2.0; 
   wgauss[6] = 0.132394152788506/2.0; 

   
 } else
 if (ngauss==13) {
   double _xigauss[39];
   _xigauss[0] =  0.333333333333333;

   _xigauss[1] =  0.479308067841923;
   _xigauss[2] =  0.260345966079038;
   _xigauss[3] =  0.260345966079038;

   _xigauss[4] =  0.869739794195568;
   _xigauss[5] =  0.065130102902216;
   _xigauss[6] =  0.065130102902216;

   _xigauss[7] =  0.638444188569809;
   _xigauss[8] =  0.638444188569809;
   _xigauss[9] =  0.312865496004875;
   _xigauss[10] = 0.048690315425316;
   _xigauss[11] = 0.048690315425316; 
   _xigauss[12] = 0.312865496004875; 

   _xigauss[13+0] =  0.333333333333333;

   _xigauss[13+1] =  0.260345966079038;
   _xigauss[13+2] =  0.479308067841923;
   _xigauss[13+3] =  0.260345966079038;

   _xigauss[13+4] =  0.065130102902216;
   _xigauss[13+5] =  0.869739794195568;
   _xigauss[13+6] =  0.065130102902216;

   _xigauss[13+7] =  0.312865496004875;
   _xigauss[13+8] =  0.048690315425316;
   _xigauss[13+9] =  0.638444188569809;
   _xigauss[13+10] = 0.638444188569809;
   _xigauss[13+11] = 0.312865496004875; 
   _xigauss[13+12] = 0.048690315425316; 

   _xigauss[26+0] =  0.333333333333333;

   _xigauss[26+1] =  0.260345966079038;
   _xigauss[26+2] =  0.260345966079038;
   _xigauss[26+3] =  0.479308067841923;

   _xigauss[26+4] =  0.065130102902216;
   _xigauss[26+5] =  0.065130102902216;
   _xigauss[26+6] =  0.869739794195568;

   _xigauss[26+7] =  0.048690315425316;
   _xigauss[26+8] =  0.312865496004875;
   _xigauss[26+9] =  0.048690315425316;
   _xigauss[26+10] = 0.312865496004875;
   _xigauss[26+11] = 0.638444188569809; 
   _xigauss[26+12] = 0.638444188569809; 

   createFaceInteg(faceindex,_xigauss);

   wgauss[0] = -0.149570044467670/2.0;

   wgauss[1] =  0.175615257433204/2.0;
   wgauss[2] =  0.175615257433204/2.0;
   wgauss[3] =  0.175615257433204/2.0;

   wgauss[4] =  0.053347235608839/2.0;
   wgauss[5] =  0.053347235608839/2.0;
   wgauss[6] =  0.053347235608839/2.0;

   wgauss[7] =  0.077113760890257/2.0;
   wgauss[8] =  0.077113760890257/2.0;
   wgauss[9] =  0.077113760890257/2.0;
   wgauss[10] = 0.077113760890257/2.0;
   wgauss[11] = 0.077113760890257/2.0;
   wgauss[12] = 0.077113760890257/2.0;
 } else if (ngauss==7*7) {
   genFrom1D(faceindex);
 }
}


void GaussRuleTriangle::genFrom1D(int faceindex) {
 int ng;
 if (ngauss==4*4) {
   ng = 4;
 } else if (ngauss==7*7) {
   ng = 7;
 }

 GaussRuleLine gr(ng,1);
 gr.pool = (double*)dbg_alloca(2*ng*sizeof(double));
 gr.init();
   
 double* _xigauss = (double*) dbg_alloca(sizeof(double)*3*ng*ng);
 int j,k;
 for(j=0;j<ng;j++) for(k=0;k<ng;k++) {
   _xigauss[0*ng*ng+j*ng+k] = (1.0+gr.xigauss[k])/2.0;
   _xigauss[1*ng*ng+j*ng+k] = (1.0-gr.xigauss[k])/2.0*
                              (1.0+gr.xigauss[j])/2.0;
   _xigauss[2*ng*ng+j*ng+k] = 1.0-_xigauss[0*ng*ng+j*ng+k]- 
                              _xigauss[1*ng*ng+j*ng+k];
   wgauss[j*ng+k] = (1.0-gr.xigauss[k])*gr.wgauss[j]*gr.wgauss[k]/8.0;
 }

 createFaceInteg(faceindex,_xigauss);
}
