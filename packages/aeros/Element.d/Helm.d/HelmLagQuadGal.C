#include <cstdio>
#include <cstdlib>
#include <Element.d/Helm.d/HelmLagQuadGal.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/dbg_alloca.h>
#include "GaussRule.h"


HelmLagQuadGal::HelmLagQuadGal(int o, int* nodenums) {
 int i;
 order = int(sqrt(double(o)));
 nn = new int[order*order];
 for(i=0;i<order*order;i++) nn[i] = nodenums[i];
}


HelmLagQuadGal::HelmLagQuadGal(const HelmLagQuadGal& e) {
 order = e.order;
 nn = new int[order*order];
 int i;
 for(i=0;i<order*order;i++) nn[i] = e.nn[i];
}


Element * HelmLagQuadGal::clone() {
 return new HelmLagQuadGal(*this);
}


void HelmLagQuadGal::renum(const int *table) {
 int i;
 for(i=0;i<order*order;i++) nn[i] = table[nn[i]];
}


void HelmLagQuadGal::renum(EleRenumMap& table) {
 int i;
 for(i=0;i<order*order;i++) nn[i] = table[nn[i]];
}


double HelmLagQuadGal::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmLagQuadGal::massMatrix not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmLagQuadGal::massMatrix(const CoordSet &cs, double *K, int cmflg) const
{
 int ordersq = order*order;

 double *x=(double*)dbg_alloca(sizeof(double)*ordersq);
 double *y=(double*)dbg_alloca(sizeof(double)*ordersq);

 int i;
 for(i=0;i<ordersq;i++) {
  Node nd = cs.getNode(nn[i]);
  x[i] = nd.x; y[i] = nd.y; 
 }
 
 for(i=0;i<ordersq*ordersq;i++) K[i] = 0.0;
 FullSquareMatrix ret(ordersq,K);

 GaussRule gr(order);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 // double kappa = prop ->kappaHelm;

 for(i=0;i<gr.ngauss;i++) {
   double xi = gr.xigauss[i];
   int j;
   for(j=0;j<gr.ngauss;j++) {
     double nu = gr.xigauss[j];
     double w = gr.wgauss[i]*gr.wgauss[j];

     double *N = (double*) dbg_alloca(sizeof(double)*3*ordersq);
     shapeFunctions(xi,nu,N);

     double J[2][2] = {{0.0,0.0},{0.0,0.0}};
     int k;
     for(k=0;k<ordersq;k++) {
       J[0][0] += N[ordersq+k]*x[k];
       J[0][1] += N[2*ordersq+k]*x[k];
       J[1][0] += N[ordersq+k]*y[k];
       J[1][1] += N[2*ordersq+k]*y[k];
     }
     double det = J[0][0]*J[1][1]-J[1][0]*J[0][1];
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f\n",
           det,x[0],y[0]);

     double (*dNdx)[2] = (double(*)[2])dbg_alloca(sizeof(double)*ordersq*2);
     for(k=0;k<ordersq;k++) {
       dNdx[k][0] = (N[ordersq+k]*J[1][1]-N[2*ordersq+k]*J[1][0])/det;
       dNdx[k][1] = -(N[ordersq+k]*J[0][1]-N[2*ordersq+k]*J[0][0])/det;
     }
        
     for(k=0;k<ordersq;k++) {
       int l;
       for(l=0;l<ordersq;l++) {
         ret[k][l] += w*(N[k]*N[l]
                        )*det;
       }
     }

   }
 }

 ret /= getProperty()->rho;
 return ret;
}


FullSquareMatrix HelmLagQuadGal::stiffness(const CoordSet &cs, double *K, int flg ) const
{
 int ordersq = order*order;

 double *x=(double*)dbg_alloca(sizeof(double)*ordersq);
 double *y=(double*)dbg_alloca(sizeof(double)*ordersq);

 int i;
 for(i=0;i<ordersq;i++) {
  Node nd = cs.getNode(nn[i]);
  x[i] = nd.x; y[i] = nd.y; 
 }
 
 for(i=0;i<ordersq*ordersq;i++) K[i] = 0.0;
 FullSquareMatrix ret(ordersq,K);

 GaussRule gr(order);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 // double kappa = prop ->kappaHelm;

 for(i=0;i<gr.ngauss;i++) {
   double xi = gr.xigauss[i];
   int j;
   for(j=0;j<gr.ngauss;j++) {
     double nu = gr.xigauss[j];
     double w = gr.wgauss[i]*gr.wgauss[j];

     double *N = (double*) dbg_alloca(sizeof(double)*3*ordersq);
     shapeFunctions(xi,nu,N);

     double J[2][2] = {{0.0,0.0},{0.0,0.0}};
     int k;
     for(k=0;k<ordersq;k++) {
       J[0][0] += N[ordersq+k]*x[k];
       J[0][1] += N[2*ordersq+k]*x[k];
       J[1][0] += N[ordersq+k]*y[k];
       J[1][1] += N[2*ordersq+k]*y[k];
     }
     double det = J[0][0]*J[1][1]-J[1][0]*J[0][1];
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f\n",
           det,x[0],y[0]);

     double (*dNdx)[2] = (double(*)[2])dbg_alloca(sizeof(double)*ordersq*2);
     for(k=0;k<ordersq;k++) {
       dNdx[k][0] = (N[ordersq+k]*J[1][1]-N[2*ordersq+k]*J[1][0])/det;
       dNdx[k][1] = -(N[ordersq+k]*J[0][1]-N[2*ordersq+k]*J[0][0])/det;
     }
        
     for(k=0;k<ordersq;k++) {
       int l;
       for(l=0;l<ordersq;l++) {
         ret[k][l] += w*(
                        +dNdx[k][0]*dNdx[l][0]+dNdx[k][1]*dNdx[l][1])*det;
       }
     }

   }
 }

 ret /= getProperty()->rho;
 return ret;
}


void HelmLagQuadGal::shapeFunctions(double xi, double eta, double *N) const {
 int i,j;
 long double *xx = (long double*)dbg_alloca(sizeof(long double)*order*order);
 long double *x = (long double*)dbg_alloca(sizeof(long double)*order);
 long double *xix = (long double*)dbg_alloca(sizeof(long double)*order);
 long double *etax = (long double*)dbg_alloca(sizeof(long double)*order);

 int xizero = -1;
 int etazero = -1;
 for(i=0;i<order;i++) {
  x[i] = -1.0+2.0*(long double)(i)/(order-1);
  xix[i] = xi-x[i]; 
  etax[i] = eta-x[i]; 
  if (std::abs(xix[i])<1e-6) xizero = i;
  if (std::abs(etax[i])<1e-6) etazero = i;
 }

 for(i=0;i<order;i++)
   for(j=0;j<order;j++)
     xx[i*order+j] = x[i]-x[j];

 for(i=0;i<order;i++) { 
   long double tmpx = 1.0;
   int k;
   for(k=0;k<order;k++) if (k!=i) tmpx *= xix[k]/xx[i*order+k];

   long double tmpxx = 0.0;
   if ((xizero==-1) || (xizero==i)) {
     for(k=0;k<order;k++)  if (k!=i) tmpxx += tmpx/xix[k];
   } else {
     tmpxx = 1.0;
     for(k=0;k<order;k++) if ((k!=xizero) && (k!=i)) tmpxx *= xix[k]; 
     for(k=0;k<order;k++) if (k!=i) tmpxx /= xx[i*order+k]; 
   }
   for(j=0;j<order;j++) {
 
     long double tmpe = 1.0;
     int k;
     for(k=0;k<order;k++) if (k!=j) tmpe *= etax[k]/xx[j*order+k];
     
     long double tmpee = 0.0;
     if ((etazero==-1) || (etazero==j)) {
       for(k=0;k<order;k++)  if (k!=j) tmpee += tmpe/etax[k];
     } else {
       tmpee = 1.0;
       for(k=0;k<order;k++) if ((k!=etazero) && (k!=j)) tmpee *= etax[k]; 
       for(k=0;k<order;k++) if (k!=j) tmpee /= xx[j*order+k]; 
     }
 
     N[j*order+i] = tmpx*tmpe;
     N[order*order+j*order+i] = tmpxx*tmpe;
     N[2*order*order+j*order+i] = tmpx*tmpee;
   }
 }
// for(i=0;i<order*order*3;i++) fprintf(stderr,"ff %f\n",N[i]);
}


FullSquareMatrix HelmLagQuadGal::acousticm(CoordSet &cs,double *K) {

 int ordersq = order*order;

 double *x=(double*)dbg_alloca(sizeof(double)*ordersq);
 double *y=(double*)dbg_alloca(sizeof(double)*ordersq);

 int i;
 for(i=0;i<ordersq;i++) {
  Node nd = cs.getNode(nn[i]);
  x[i] = nd.x; y[i] = nd.y; 
 }
 
 for(i=0;i<ordersq*ordersq;i++) K[i] = 0.0;
 FullSquareMatrix ret(ordersq,K);

 GaussRule gr(order);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

 double kappa = prop ->kappaHelm;

 for(i=0;i<gr.ngauss;i++) {
   double xi = gr.xigauss[i];
   int j;
   for(j=0;j<gr.ngauss;j++) {
     double nu = gr.xigauss[j];
     double w = gr.wgauss[i]*gr.wgauss[j];

     double *N = (double*) dbg_alloca(sizeof(double)*3*ordersq);
     shapeFunctions(xi,nu,N);

     double J[2][2] = {{0.0,0.0},{0.0,0.0}};
     int k;
     for(k=0;k<ordersq;k++) {
       J[0][0] += N[ordersq+k]*x[k];
       J[0][1] += N[2*ordersq+k]*x[k];
       J[1][0] += N[ordersq+k]*y[k];
       J[1][1] += N[2*ordersq+k]*y[k];
     }
     double det = J[0][0]*J[1][1]-J[1][0]*J[0][1];
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f\n",
           det,x[0],y[0]);

     double (*dNdx)[2] = (double(*)[2])dbg_alloca(sizeof(double)*ordersq*2);
     for(k=0;k<ordersq;k++) {
       dNdx[k][0] = (N[ordersq+k]*J[1][1]-N[2*ordersq+k]*J[1][0])/det;
       dNdx[k][1] = -(N[ordersq+k]*J[0][1]-N[2*ordersq+k]*J[0][0])/det;
     }
        
     for(k=0;k<ordersq;k++) {
       int l;
       for(l=0;l<ordersq;l++) {
         ret[k][l] += w*(-kappa*kappa*N[k]*N[l]
                        +dNdx[k][0]*dNdx[l][0]+dNdx[k][1]*dNdx[l][1])*det;
       }
     }

   }
 }

 ret /= getProperty()->rho;
 return ret;
}

void HelmLagQuadGal::wErrors(CoordSet&cs,
  double *l2e, double *h1e, double *l2, double *h1,
  ComplexD *u, double kappa, double *waveDir) {

 int ordersq = order*order;

 double *x=(double*)dbg_alloca(sizeof(double)*ordersq);
 double *y=(double*)dbg_alloca(sizeof(double)*ordersq);

 int i;
 for(i=0;i<ordersq;i++) {
  Node nd = cs.getNode(nn[i]);
  x[i] = nd.x; y[i] = nd.y; 
 }

 GaussRule gr(order);
 gr.pool = (double*)dbg_alloca(sizeof(double)*2*gr.ngauss);
 gr.init();

//  double kappa = prop ->kappaHelm;

 for(i=0;i<gr.ngauss;i++) {
   double xi = gr.xigauss[i];
   int j;
   for(j=0;j<gr.ngauss;j++) {
     double nu = gr.xigauss[j];
     double w = gr.wgauss[i]*gr.wgauss[j];

     double *N = (double*) dbg_alloca(sizeof(double)*3*ordersq);
     shapeFunctions(xi,nu,N);

     double J[2][2] = {{0.0,0.0},{0.0,0.0}};
     int k;
     for(k=0;k<ordersq;k++) {
       J[0][0] += N[ordersq+k]*x[k];
       J[0][1] += N[2*ordersq+k]*x[k];
       J[1][0] += N[ordersq+k]*y[k];
       J[1][1] += N[2*ordersq+k]*y[k];
     }
     double det = J[0][0]*J[1][1]-J[1][0]*J[0][1];
     if (det<0.0) fprintf(stderr,"Negative Jacobian: %f %f %f\n",
           det,x[0],y[0]);

     double (*dNdx)[2] = (double(*)[2])dbg_alloca(sizeof(double)*ordersq*2);
     for(k=0;k<ordersq;k++) {
       dNdx[k][0] = (N[ordersq+k]*J[1][1]-N[2*ordersq+k]*J[1][0])/det;
       dNdx[k][1] = -(N[ordersq+k]*J[0][1]-N[2*ordersq+k]*J[0][0])/det;
     }
        
       double xx = 0.0;
       double yy = 0.0;
       ComplexD uu = 0.0;
       ComplexD ux = 0.0;
       ComplexD uy = 0.0;
       int l;
       for(l=0;l<ordersq;l++) {
         xx += N[l]*x[l];
         yy += N[l]*y[l];
         uu += N[l]*u[l];
         ux += dNdx[l][0]*u[l];
         uy += dNdx[l][1]*u[l];
       }
       double vth = xx*waveDir[0]+yy*waveDir[1];
       *l2e += w*det*norm(exp(ComplexD(0.0,kappa*vth))- uu);
       *h1e += w*det*(
               norm(ComplexD(0.0,kappa*waveDir[0])*exp(ComplexD(0.0,kappa*vth))-
                    ux) +
               norm(ComplexD(0.0,kappa*waveDir[1])*exp(ComplexD(0.0,kappa*vth))-
                   uy));
       *l2 += w*det*norm(exp(ComplexD(0.0,kappa*vth)));
       *h1 += w*det*(
            norm(ComplexD(0.0,kappa*waveDir[0])*exp(ComplexD(0.0,kappa*vth))) +
            norm(ComplexD(0.0,kappa*waveDir[1])*exp(ComplexD(0.0,kappa*vth))));
     }
 }
}



void HelmLagQuadGal::edgeShapeFunctions(int n1, int n2, int *ng,
                                       double **gw, double **N) {

 GaussRule gr(order);
 gr.pool = new double[2*gr.ngauss];
 gr.init();

 *ng = gr.ngauss;
 *gw = gr.wgauss;

 int nc[4] = { nn[0], nn[order-1], nn[order*order-1], nn[order*(order-1)] };
 
 int i1=-1;
 int i2=-1;
 int i;
 for(i=0;i<4;i++) {
   if (nc[i]==n1) i1 = i;
   if (nc[i]==n2) i2 = i;
 }

 int *ien = (int*)dbg_alloca(sizeof(int)*order);
 int edge;
 double si;
 for(i=0;i<4;i++) {
   if (i1==0 && i2==1) {
     edge = 1;
     si = 1;
     int j;
     for(j=0;j<order;j++) ien[j] = j;
   } else if  (i1==1 && i2==0) {
     edge = 1;
     si = -1;
     int j;
     for(j=order-1;j>=0;j--) ien[order-1-j] = j;
   } else if (i1==1 && i2==2) {
     edge = 2;
     si = 1;
     int j;
     for(j=0;j<order;j++) ien[j] = order-1+j*order;
   } else if (i1==2 && i2==1) {
     edge = 2;
     si = -1;
     int j;
     for(j=order-1;j>=0;j--) ien[order-1-j] = order-1+j*order;
   } else if (i1==3 && i2==2) {
     edge = 3;
     si = 1;
     int j;
     for(j=0;j<order;j++) ien[order-1-j] = order*order-1-j;
   } else if (i1==2 && i2==3) {
     edge = 3;
     si = -1;
     int j;
     for(j=order-1;j>=0;j--) ien[j] = order*order-1-j;
   } else if (i1==0 && i2==3) {
     edge = 4;
     si = 1;
     int j;
     for(j=0;j<order;j++) ien[j] = j*order;
   } else if (i1==3 && i2==0) {
     edge = 4;
     si = -1;
     int j;
     for(j=order-1;j>=0;j--) ien[order-1-j] = j*order;
   } else {
     fprintf(stderr,"Error in HelmLagQuadGal::edgeShapeFunctions.\n");
     exit(-1);
   }
 }


 *N = new double[2*order*gr.ngauss];

 int j;
 for(j=0;j<gr.ngauss;j++) {
   int der;
   double *NN = (double*) dbg_alloca(sizeof(double)*3*order*order);
   if (edge==1) {
     shapeFunctions(gr.xigauss[j],-1,NN);
     der = 1;
   } else if (edge==2) {
     shapeFunctions(1,gr.xigauss[j],NN);
     der = 2;
   } else if (edge==3) {
     shapeFunctions(gr.xigauss[j],1,NN);
     der = 1;
   } else {
     shapeFunctions(-1,gr.xigauss[j],NN);
     der = 2;
   }
   for(i=0;i<order;i++) {
    (*N)[j*2*order+i] = NN[ien[i]];
    (*N)[j*2*order+order+i] = si*NN[ien[i]+der*order*order];
   }
 }

}


int* HelmLagQuadGal::nodes(int *p) const {

 if(p == 0) p = new int[order*order];
 int i;
 for(i=0;i<order*order;i++) p[i] = nn[i];
 return p;
}


int* HelmLagQuadGal::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 if(p == 0) p = new int[order*order];
 for(i=0;i<order*order;i++) dsa.number(nn[i],DofSet::Helm,p+i);
 return p;
}


void HelmLagQuadGal::markDofs(DofSetArray &dsa) const {

 int i;
 for(i=0;i<order*order;i++) dsa.mark(nn[i],DofSet::Helm);
}


void HelmLagQuadGal::addFaces(PolygonSet *pset) {

 fprintf(stderr,"HelmLagQuadGal::addFaces not implemented.\n");
/*
 int i;
 for(i=0;i<order-1;i++)
   pset->addLine(this,nn[i],nn[i+1]);
 for(i=0;i<order-1;i++)
    pset->addLine(this,nn[i*order+order-1],nn[(i+1)*order+order-1]);
 for(i=0;i<order-1;i++)
   pset->addLine(this,nn[order*order-1-i],nn[order*order-1-i-1]);
 for(i=0;i<order-1;i++)
   pset->addLine(this,nn[(order-1)*order-i*order],nn[(order-1)*order-(i+1)*order]);
*/
}

double HelmLagQuadGal::weight() const
{
	return order;
}

double HelmLagQuadGal::trueWeight() const
{
	return order;
}
