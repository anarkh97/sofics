#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/HelmIsoParamTri.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include "HelmIsoParamTri.h"


HelmIsoParamTri::HelmIsoParamTri(int o, int* nodenums) {

 if (o==3) order = 2;
 else if (o==6) order = 3;
 else if (o==10) order = 4;
 else if (o==15) order = 5;
 else if (o==21) order = 6;
 else {
   fprintf(stderr,"Order too high in HelmIsoParamTri::HelmIsoParamTri.\n");
   exit(-1);
 }

 int ordersq = (order*(order+1))/2;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = nodenums[i];
}


HelmIsoParamTri::HelmIsoParamTri(const HelmIsoParamTri& e) {
 order = e.order;
 int ordersq = (order*(order+1))/2;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = e.nn[i];
}


Element * HelmIsoParamTri::clone() {
 return new HelmIsoParamTri(*this);
}


void HelmIsoParamTri::renum(const int *table) {
 int i;
 int ordersq = (order*(order+1))/2;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}

void HelmIsoParamTri::renum(EleRenumMap& table) {
 int i;
 int ordersq = (order*(order+1))/2;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}

extern bool useFull;

int* HelmIsoParamTri::nodes(int *p) const {
 int ordersq = (order*(order+1))/2;
 if(useFull) {
   if(p == 0) p = new int[ordersq];
   int i;
   for(i=0;i<ordersq;i++) p[i] = nn[i];
 } else {
   if(p == 0) p = new int[3];
   p[0] = nn[0];
   p[1] = nn[order-1];
   p[2] = nn[ordersq-1];
 }
 return p;
}


int* HelmIsoParamTri::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 int ordersq = (order*(order+1))/2;
 if(p == 0) p = new int[ordersq];
 for(i=0;i<ordersq;i++) dsa.number(nn[i],DofSet::Helm,p+i);
 return p;
}


void HelmIsoParamTri::markDofs(DofSetArray &dsa) const {

 int i;
 int ordersq = (order*(order+1))/2;
 for(i=0;i<ordersq;i++) dsa.mark(nn[i],DofSet::Helm);
}


double HelmIsoParamTri::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmIsoParamTri::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmIsoParamTri::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 GalMassFunction2d<double> f(ordersq,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 int gorder = 7*7;
 if (order<=3) gorder = 4*4;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(ordersq,K);
 return ret;
}


FullSquareMatrix HelmIsoParamTri::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 GalStiffFunction2d<double> f(ordersq,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 int gorder = 7*7;
 if (order<=3) gorder = 4*4;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(ordersq,K);
 return ret;
}


FullSquareMatrixC HelmIsoParamTri::massMatrix(const CoordSet &cs,
                                               complex<double> *K) const {

 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 ipu.zeroOut<complex<double> > (ordersq*ordersq,K);
 int gorder = 7*7;
 if (order<=3) gorder = 4*4;

 if (prop->fp.PMLtype==1) {
   double Rm[2] = { -prop->fp.Rx, -prop->fp.Ry };
   double Rp[2] = { prop->fp.Rx, prop->fp.Ry };
   double Sm[2] = { -prop->fp.Sx, -prop->fp.Sy };
   double Sp[2] = { prop->fp.Sx, prop->fp.Sy };
   CartPMLGalMassFunction2d f(ordersq,Rm,Sm,Rp,Sp,prop->fp.gamma,K);
   ipu.areaInt2d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==2 || prop->fp.PMLtype==3) {
   CylPMLGalMassFunction2d f(ordersq,prop->fp.Rx,prop->fp.Sx,prop->fp.gamma, K);
   ipu.areaInt2d(xyz, f, gorder);
 } else {
   GalMassFunction2d<complex<double> > f(ordersq,K);
   ipu.areaInt2d(xyz, f, gorder);
 }

 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;

 FullSquareMatrixC ret(ordersq,K);
 return ret;
}


FullSquareMatrixC HelmIsoParamTri::stiffness(const CoordSet &cs,
                                              complex<double> *K) const {

 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 ipu.zeroOut< complex<double> > (ordersq*ordersq,K);
 int gorder = 7*7;
 if (order<=3) gorder = 4*4;

 if (prop->fp.PMLtype==1) {
   double Rm[2] = { -prop->fp.Rx, -prop->fp.Ry };
   double Rp[2] = { prop->fp.Rx, prop->fp.Ry };
   double Sm[2] = { -prop->fp.Sx, -prop->fp.Sy };
   double Sp[2] = { prop->fp.Sx, prop->fp.Sy };
   CartPMLGalStiffFunction2d f(ordersq,Rm,Sm,Rp,Sp,prop->fp.gamma,K);
   ipu.areaInt2d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==2 || prop->fp.PMLtype==3) {
   CylPMLGalStiffFunction2d f(ordersq,prop->fp.Rx,prop->fp.Sx,prop->fp.gamma,K);
   ipu.areaInt2d(xyz, f, gorder);
 } else {
   GalStiffFunction2d<complex<double> > f(ordersq,K);
   ipu.areaInt2d(xyz, f, gorder);
 }

 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;

 FullSquareMatrixC ret(ordersq,K);
 return ret;
}


int
HelmIsoParamTri::numNodes() const {
  if(useFull)
    return (order*(order+1))/2;
  else
    return(3);   // to ignore effect of mid-size nodes in dec
}

int HelmIsoParamTri::getTopNumber() const
{
    return -1;
}

double HelmIsoParamTri::weight() const
{
	return order-1;
}

double HelmIsoParamTri::trueWeight() const
{
	return order-1;
}
