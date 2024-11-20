#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/LEIsoParamQuad.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include "LEIsoParamQuad.h"


LEIsoParamQuad::LEIsoParamQuad(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/2.0)));
 int ordersq = order*order;
 nn = new int[ordersq];
 for(i=0;i<ordersq;i++) nn[i] = nodenums[i];
}


LEIsoParamQuad::LEIsoParamQuad(const LEIsoParamQuad& e) {
 order = e.order;
 int ordersq = order*order;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = e.nn[i];
}


Element * LEIsoParamQuad::clone() {
 return new LEIsoParamQuad(*this);
}


void LEIsoParamQuad::renum(const int *table) {
 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}


void LEIsoParamQuad::renum(EleRenumMap& table) {
 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}


int* LEIsoParamQuad::nodes(int *p) const {

 int ordersq = order*order;
 if(p == 0) p = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) p[i] = nn[i];
 return p;
}


int* LEIsoParamQuad::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 int ordersq = order*order;
 if(p == 0) p = new int[2*ordersq];
 for(i=0;i<ordersq;i++) dsa.number(nn[i],DofSet::Xdisp|DofSet::Ydisp,p+2*i);
 return p;
}


void LEIsoParamQuad::markDofs(DofSetArray &dsa) const {

 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) dsa.mark(nn[i],DofSet::Xdisp|DofSet::Ydisp);
}


double LEIsoParamQuad::getMass(const CoordSet &cs) const {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 double area(0);
 AreaFunction2d f(&area);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.areaInt2d(xyz, f, gorder);

 return area*prop->eh*prop->rho;
}

double LEIsoParamQuad::getMassThicknessSensitivity(CoordSet &cs) {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 double area(0);
 AreaFunction2d f(&area);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.areaInt2d(xyz, f, gorder);

 return area*prop->rho;
}


FullSquareMatrix LEIsoParamQuad::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 LEMassFunction2d f(2*ordersq,prop->rho*prop->eh,K);
 ipu.zeroOut<double> (4*ordersq*ordersq,K);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(2*ordersq,K);

 FullSquareMatrix ret(2*ordersq,K);
 return ret;
}


FullSquareMatrix LEIsoParamQuad::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 LEStiffFunction2d f(2*ordersq,prop->E, prop->nu,K);
 ipu.zeroOut<double> (4*ordersq*ordersq,K);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.areaInt2d(xyz, f, gorder);
 ipu.symmetrize(2*ordersq,K);

 FullSquareMatrix ret(2*ordersq,K);
 return ret;
}


extern bool useFull;

int
LEIsoParamQuad::numNodes() const {
  if(useFull)
    return order*order;
  else
    return(4);   // to ignore effect of mid-size nodes in dec
}

int LEIsoParamQuad::getTopNumber() const
{
 return -1;
}

double LEIsoParamQuad::weight() const
{
	return order;
}

double LEIsoParamQuad::trueWeight() const
{
	return order;
}

