#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Tetra.d/ThermIsoParamTetra.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/dbg_alloca.h>


ThermIsoParamTetra::ThermIsoParamTetra(int o, int* nodenums) {

 if (o==4) order = 2;
 else if (o==10) order = 3;
 else if (o==20) order = 4;
 else if (o==35) order = 5;
 else if (o==56) order = 6;
 else {
   fprintf(stderr,"Order too high in ThermIsoParamTetra::ThermIsoParamTetra.\n");
   exit(-1);
 }
 
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


ThermIsoParamTetra::ThermIsoParamTetra(const ThermIsoParamTetra& e) {
 order = e.order;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * ThermIsoParamTetra::clone() {
 return new ThermIsoParamTetra(*this);
}


void ThermIsoParamTetra::renum(const int *table) {
 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

void ThermIsoParamTetra::renum(EleRenumMap& table) {
 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}


int* ThermIsoParamTetra::nodes(int *p) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 if(p == 0) p = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) p[i] = nn[i];
 return p;
}


int* ThermIsoParamTetra::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 if(p == 0) p = new int[orderc];
 for(i=0;i<orderc;i++) dsa.number(nn[i],DofSet::Temp,p+i);
 return p;
}


void ThermIsoParamTetra::markDofs(DofSetArray &dsa) const {

 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::Temp);
}


double ThermIsoParamTetra::getMass(const CoordSet&) const {
 fprintf(stderr,"ThermIsoParamTetra::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix ThermIsoParamTetra::stiffness(const CoordSet &cs, double *K, int flg ) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 GalStiffFunction<double> f(orderc,K);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(f.ncolumns(),K);
// Conductivity coefficient
 double k = prop->k;
 for(int i=0;i<orderc*orderc;i++) K[i] *= k;

 FullSquareMatrix ret(orderc,K);
 return ret;
}


FullSquareMatrix ThermIsoParamTetra::massMatrix(const CoordSet &cs,
                                               double *K, int cmflg) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 ipu.zeroOut<double> (orderc*orderc,K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;

 GalMassFunction<double> f(orderc,K);
 ipu.volumeInt3d(xyz, f, gorder);

 ipu.symmetrize(orderc,K);
 double capacitance = prop->rho*prop->Q;
 for(int i=0;i<orderc*orderc;i++) K[i] *= capacitance;

 FullSquareMatrix ret(orderc,K);
 return ret;
}


extern bool useFull;

int ThermIsoParamTetra::numNodes() const {
  //Not tested -JF
  if(useFull)
    return (order*(order+1)*(order+2))/6;
  else
    return(4);   // to ignore effect of mid-size nodes in dec
}

int ThermIsoParamTetra::getTopNumber() const {
  return 150;//5;
}


