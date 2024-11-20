#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/HelmIsoParamQuad.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include "HelmIsoParamQuad.h"


HelmIsoParamQuad::HelmIsoParamQuad(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/2.0)));
 int ordersq = order*order;
 nn = new int[ordersq];
 for(i=0;i<ordersq;i++) nn[i] = nodenums[i];
}


HelmIsoParamQuad::HelmIsoParamQuad(const HelmIsoParamQuad& e) {
 order = e.order;
 int ordersq = order*order;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = e.nn[i];
}


Element * HelmIsoParamQuad::clone() {
 return new HelmIsoParamQuad(*this);
}


void HelmIsoParamQuad::renum(const int *table) {
 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}

void HelmIsoParamQuad::renum(EleRenumMap& table) {
 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}

extern bool useFull;

int* HelmIsoParamQuad::nodes(int *p) const {
 if(useFull) {
   int ordersq = order*order;
   if(p == 0) p = new int[ordersq];
   int i;
   for(i=0;i<ordersq;i++) p[i] = nn[i];
 } else {
   if(p == 0) p = new int[4];
   p[0] = nn[0];
   p[1] = nn[order-1];
   p[2] = nn[order*order-1];
   p[3] = nn[order*order-order]; 
 }
 return p;
}


int* HelmIsoParamQuad::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 int ordersq = order*order;
 if(p == 0) p = new int[ordersq];
 for(i=0;i<ordersq;i++) dsa.number(nn[i],DofSet::Helm,p+i);
 return p;
}


void HelmIsoParamQuad::markDofs(DofSetArray &dsa) const {

 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) dsa.mark(nn[i],DofSet::Helm);
}


double HelmIsoParamQuad::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmIsoParamQuad::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmIsoParamQuad::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 GalMassFunction2d<double> f(ordersq,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;
 
 FullSquareMatrix ret(ordersq,K);
 return ret;
}


FullSquareMatrix HelmIsoParamQuad::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 GalStiffFunction2d<double> f(ordersq,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(ordersq,K);
 return ret;
}


FullSquareMatrixC HelmIsoParamQuad::massMatrix(const CoordSet &cs,
                                               complex<double> *K) const {

 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 ipu.zeroOut<complex<double> > (ordersq*ordersq,K);
 int gorder = 7;
 if (order<=3) gorder = 4;

 if (prop->fp.PMLtype==1) {
   double Rm[2] = { -prop->fp.Rx, -prop->fp.Ry };
   double Rp[2] = { prop->fp.Rx, prop->fp.Ry };
   double Sm[2] = { -prop->fp.Sx, -prop->fp.Sy };
   double Sp[2] = { prop->fp.Sx, prop->fp.Sy };
   CartPMLGalMassFunction2d f(ordersq,Rm,Sm,Rp,Sp,prop->fp.gamma,K);
   ipu.areaInt2d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==2 || prop->fp.PMLtype==3) {
   CylPMLGalMassFunction2d f(ordersq,prop->fp.Rx,prop->fp.Sx,prop->fp.gamma,K);
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


FullSquareMatrixC HelmIsoParamQuad::stiffness(const CoordSet &cs,
                                              complex<double> *K) const {

 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 ipu.zeroOut< complex<double> > (ordersq*ordersq,K);
 int gorder = 7;
 if (order<=3) gorder = 4;

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
HelmIsoParamQuad::numNodes() const {
  if(useFull)
    return order*order;
  else
    return(4);   // to ignore effect of mid-size nodes in dec
}

int HelmIsoParamQuad::getTopNumber() const
{
    return -1;
}

double HelmIsoParamQuad::weight() const
{
	return order;
}

double HelmIsoParamQuad::trueWeight() const
{
	return order;
}

