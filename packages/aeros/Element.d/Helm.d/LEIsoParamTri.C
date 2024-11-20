#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/LEIsoParamTri.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include "LEIsoParamTri.h"


LEIsoParamTri::LEIsoParamTri(int o, int* nodenums) {

 if (o==3) order = 2;
 else if (o==6) order = 3;
 else if (o==10) order = 4;
 else if (o==15) order = 5;
 else if (o==21) order = 6;
 else {
   fprintf(stderr,"Order too high in LEIsoParamTri::LEIsoParamTri.\n");
   exit(-1);
 }

 int ordersq = (order*(order+1))/2;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = nodenums[i];
}


LEIsoParamTri::LEIsoParamTri(const LEIsoParamTri& e) {
 order = e.order;
 int ordersq = (order*(order+1))/2;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = e.nn[i];
}


Element * LEIsoParamTri::clone() {
 return new LEIsoParamTri(*this);
}


void LEIsoParamTri::renum(const int *table) {
 int i;
 int ordersq = (order*(order+1))/2;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}


void LEIsoParamTri::renum(EleRenumMap& table) {
 int i;
 int ordersq = (order*(order+1))/2;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}


int* LEIsoParamTri::nodes(int *p) const {

 int ordersq = (order*(order+1))/2;
 if(p == 0) p = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) p[i] = nn[i];
 return p;
}


int* LEIsoParamTri::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 int ordersq = (order*(order+1))/2;
 if(p == 0) p = new int[2*ordersq];
 for(i=0;i<ordersq;i++) dsa.number(nn[i],DofSet::Xdisp|DofSet::Ydisp,p+2*i);
 return p;
}


void LEIsoParamTri::markDofs(DofSetArray &dsa) const {

 int i;
 int ordersq = (order*(order+1))/2;
 for(i=0;i<ordersq;i++) dsa.mark(nn[i],DofSet::Xdisp|DofSet::Ydisp);
}


double LEIsoParamTri::getMass(const CoordSet &cs) const {
 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 double area(0);
 AreaFunction2d f(&area);
 int gorder = 7*7;
 if (order<=3) gorder = 13;
 ipu.areaInt2d(xyz, f, gorder);

 return area*prop->eh*prop->rho; 
}


double LEIsoParamTri::getMassThicknessSensitivity(CoordSet &cs) {
 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 double area(0);
 AreaFunction2d f(&area);
 int gorder = 7*7;
 if (order<=3) gorder = 13;
 ipu.areaInt2d(xyz, f, gorder);

 return area*prop->rho; 
}


FullSquareMatrix LEIsoParamTri::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 LEMassFunction2d f(2*ordersq,prop->rho*prop->eh,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 int gorder = 7*7;
 if (order<=3) gorder = 13;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(2*ordersq,K);

 FullSquareMatrix ret(2*ordersq,K);
 return ret;
}


FullSquareMatrix LEIsoParamTri::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils2dTri ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 LEStiffFunction2d f(2*ordersq,prop->E, prop->nu,K);
 ipu.zeroOut<double> (4*ordersq*ordersq,K);
 int gorder = 7*7;
 if (order<=3) gorder = 13;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(2*ordersq,K);

 FullSquareMatrix ret(2*ordersq,K);
 return ret;
}

extern bool useFull;

int
LEIsoParamTri::numNodes() const {
  if(useFull)
    return (order*(order+1))/2;
  else
    return(3);   // to ignore effect of mid-size nodes in dec
}

int LEIsoParamTri::getTopNumber() const
{
 return -1;
}

double LEIsoParamTri::weight() const
{
	return order-1;
}

double LEIsoParamTri::trueWeight() const
{
	return order - 1;
}
