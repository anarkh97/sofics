#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/ThermIsoParamHexa.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/dbg_alloca.h>
#define O3 4
#define O4 5


ThermIsoParamHexa::ThermIsoParamHexa(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/3.0)));
 int orderc = order*order*order;
 nn = new int[orderc];
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


ThermIsoParamHexa::ThermIsoParamHexa(const ThermIsoParamHexa& e) {
 order = e.order;
 int orderc = order*order*order;
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * ThermIsoParamHexa::clone() {
 return new ThermIsoParamHexa(*this);
}


void ThermIsoParamHexa::renum(const int *table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

void ThermIsoParamHexa::renum(EleRenumMap& table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

extern bool useFull;

int* ThermIsoParamHexa::nodes(int *p) const {
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


int* ThermIsoParamHexa::dofs(DofSetArray &dsa, int *p) const  {
 int i;
 int orderc = order*order*order;
 if(p == 0) p = new int[orderc];
 for(i=0;i<orderc;i++) {
   dsa.number(nn[i],DofSet::Temp,p+i);
 }
 return p;
}


void ThermIsoParamHexa::markDofs(DofSetArray &dsa) const {

 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::Temp);
}


double ThermIsoParamHexa::getMass(const CoordSet&) const {
 fprintf(stderr,"ThermIsoParamHexa::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix ThermIsoParamHexa::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 GalMassFunction<double> f(orderc,K);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),K);
 int gorder = 7;
 if (order<=3) gorder = O3;
 else if (order<=4) gorder = O4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(f.ncolumns(),K);
 double capacitance = prop->rho*prop->Q;
 for(int i=0;i<orderc*orderc;i++) K[i] *= capacitance;
 int j;

 FullSquareMatrix ret(orderc,K);
 if (fl==0) {
   // build a lumped mass matrix - JFD
   if (order>3) fprintf(stderr,"Warning, HRZ lumping is applied and does not converge for order>3, use spectral elements instead\n");
   // HRZ method, See Hinton, Rock & Zienkiewicz, International journal of Earthquake
   //  engineering & structural dynamics, #4, 1976, pp 245-249: "Note on mass lumping
   //  in related process in the finite element method"
   // presently wrks only for acoustics
   //fprintf(stderr,"   Mass matrix is lumped\n");
   double MM = 0.0;//total mass of element
   double MD = 0.0;//mass of the diagonal
   for (int i = 0; i < order; ++i) {
     MD += ret[i][i];
     for (j = 0; j < order; ++j) {
       MM += ret[i][j];
       if (i != j)
         ret[i][j] =0.0;
     }
   }
   for (int i = 0; i < order; ++i)
     ret[i][i] = ret[i][i]*MM/MD;
 }
 return ret;
}


FullSquareMatrix ThermIsoParamHexa::stiffness(const CoordSet &cs, double *K, int flg ) const {
 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 GalStiffFunction<double> f(orderc,K);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),K);
 int gorder = 7;
 if (order<=3) gorder = O3;
 else if (order<=4) gorder = O4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(f.ncolumns(),K);
 // Conductivity coefficient
 double k = prop->k;
 for(int i=0;i<orderc*orderc;i++) K[i] *= k;

 FullSquareMatrix ret(orderc,K);
 return ret;
}

extern bool useFull;

int
ThermIsoParamHexa::numNodes() const {
  if(useFull)
    return order*order*order;
  else
    return(8);   // to ignore effect of mid-size nodes in dec
}

double ThermIsoParamHexa::weight() const
{
	return order;
}

double ThermIsoParamHexa::trueWeight() const
{
	return order;
}

int ThermIsoParamHexa::getElementType() const
{
	return 109;
}

