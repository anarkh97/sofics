#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/HelmSpectralIsoParamHexa.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>

HelmSpectralIsoParamHexa::HelmSpectralIsoParamHexa(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/3.0)));
 int orderc = order*order*order;
 nn = new int[orderc];
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


HelmSpectralIsoParamHexa::HelmSpectralIsoParamHexa(const HelmSpectralIsoParamHexa& e) {
 order = e.order;
 int orderc = order*order*order;
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * HelmSpectralIsoParamHexa::clone() {
 return new HelmSpectralIsoParamHexa(*this);
}


void HelmSpectralIsoParamHexa::renum(const int *table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

void HelmSpectralIsoParamHexa::renum(EleRenumMap& table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

extern bool useFull;

int* HelmSpectralIsoParamHexa::nodes(int *p) const {
if(useFull) {
//cerr << "order = " << order << endl;//JF
 int orderc = order*order*order;
 if(p == 0) p = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) {
//cerr << "node " << i+1 << " = " << nn[i] << endl;
   p[i] = nn[i];
 }
}
else {
 if(p == 0) p = new int[8];
 p[0] = nn[0];
 p[1] = nn[order-1];
 p[2] = nn[order*order-1];
 p[3] = nn[order*order-order];                                                                    p[4] = nn[order*order*order-order*order];
 p[5] = nn[order*order*order-order*order+order-1];
 p[6] = nn[order*order*order-1];
 p[7] = nn[order*order*order-order];
}
 return p;
}


int* HelmSpectralIsoParamHexa::dofs(DofSetArray &dsa, int *p) const  {
 int i;
 int orderc = order*order*order;
 if(p == 0) p = new int[orderc];
 for(i=0;i<orderc;i++) {
   dsa.number(nn[i],DofSet::Helm,p+i);
 }
 return p;
}


void HelmSpectralIsoParamHexa::markDofs(DofSetArray &dsa) const {

 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::Helm);
}


void HelmSpectralIsoParamHexa::addFaces(PolygonSet *pset) {

 int i;
 int corner[4];
 IsoParamUtils ipu(order);
 for(i=0;i<6;i++) {
   ipu.cornerindeces(i,corner); 
   pset->addQuad(this,nn[corner[0]],nn[corner[1]],nn[corner[2]],nn[corner[3]]);
 }
}


double HelmSpectralIsoParamHexa::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmSpectralIsoParamHexa::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmSpectralIsoParamHexa::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 GalMassFunction<double> f(orderc,K);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),K);
 ipu.spectralVolumeInt3d(xyz, f);
 ipu.symmetrize(f.ncolumns(),K);
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(orderc,K);
 return ret;
}


FullSquareMatrix HelmSpectralIsoParamHexa::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 GalStiffFunction<double> f(orderc,K);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),K);
 ipu.spectralVolumeInt3d(xyz, f);
 ipu.symmetrize(f.ncolumns(),K);
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(orderc,K);
 return ret;
}

extern bool useFull;

int
HelmSpectralIsoParamHexa::numNodes() const {
  if(useFull)
    return order*order*order;
  else
    return(8);   // to ignore effect of mid-size nodes in dec
}

double HelmSpectralIsoParamHexa::weight() const
{
	return order;
}

double HelmSpectralIsoParamHexa::trueWeight() const
{
	return order;
}

