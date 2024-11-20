#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/HelmSpectralIsoParamQuad.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>


HelmSpectralIsoParamQuad::HelmSpectralIsoParamQuad(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/2.0)));
 int ordersq = order*order;
 nn = new int[ordersq];
 for(i=0;i<ordersq;i++) nn[i] = nodenums[i];
}


HelmSpectralIsoParamQuad::HelmSpectralIsoParamQuad(const HelmSpectralIsoParamQuad& e) {
 order = e.order;
 int ordersq = order*order;
 nn = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) nn[i] = e.nn[i];
}


Element * HelmSpectralIsoParamQuad::clone() {
 return new HelmSpectralIsoParamQuad(*this);
}


void HelmSpectralIsoParamQuad::renum(const int *table) {
 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}

void HelmSpectralIsoParamQuad::renum(EleRenumMap& table) {
 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) nn[i] = table[nn[i]];
}


int* HelmSpectralIsoParamQuad::nodes(int *p) const {

 int ordersq = order*order;
 if(p == 0) p = new int[ordersq];
 int i;
 for(i=0;i<ordersq;i++) p[i] = nn[i];
 return p;
}


int* HelmSpectralIsoParamQuad::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 int ordersq = order*order;
 if(p == 0) p = new int[ordersq];
 for(i=0;i<ordersq;i++) dsa.number(nn[i],DofSet::Helm,p+i);
 return p;
}


void HelmSpectralIsoParamQuad::markDofs(DofSetArray &dsa) const {

 int i;
 int ordersq = order*order;
 for(i=0;i<ordersq;i++) dsa.mark(nn[i],DofSet::Helm);
}


double HelmSpectralIsoParamQuad::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmSpectralIsoParamQuad::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmSpectralIsoParamQuad::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 GalMassFunction2d<double> f(ordersq,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 ipu.spectralAreaInt2d(xyz, f);
 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;
 
 FullSquareMatrix ret(ordersq,K);
 return ret;
}


FullSquareMatrix HelmSpectralIsoParamQuad::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils2d ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 GalStiffFunction2d<double> f(ordersq,K);
 ipu.zeroOut<double> (ordersq*ordersq,K);
 ipu.spectralAreaInt2d(xyz, f);
 ipu.symmetrize(ordersq,K);
 int i;
 for(i=0;i<ordersq*ordersq;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(ordersq,K);
 return ret;
}


extern bool useFull;

int
HelmSpectralIsoParamQuad::numNodes() const {
  if(useFull)
    return order*order;
  else
    return(4);   // to ignore effect of mid-size nodes in dec
}

double HelmSpectralIsoParamQuad::weight() const
{
	return order;
}

double HelmSpectralIsoParamQuad::trueWeight() const
{
	return order;
}

