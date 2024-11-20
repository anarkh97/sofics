#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/HelmIsoParamHexa.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/dbg_alloca.h>
#define O3 4
#define O4 5


HelmIsoParamHexa::HelmIsoParamHexa(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/3.0)));
 int orderc = order*order*order;
 nn = new int[orderc];
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


HelmIsoParamHexa::HelmIsoParamHexa(const HelmIsoParamHexa& e) {
 order = e.order;
 int orderc = order*order*order;
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * HelmIsoParamHexa::clone() {
 return new HelmIsoParamHexa(*this);
}


void HelmIsoParamHexa::renum(const int *table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

void HelmIsoParamHexa::renum(EleRenumMap& table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

extern bool useFull;

int* HelmIsoParamHexa::nodes(int *p) const {
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


int* HelmIsoParamHexa::dofs(DofSetArray &dsa, int *p) const  {
 int i;
 int orderc = order*order*order;
 if(p == 0) p = new int[orderc];
 for(i=0;i<orderc;i++) {
   dsa.number(nn[i],DofSet::Helm,p+i);
 }
 return p;
}


void HelmIsoParamHexa::markDofs(DofSetArray &dsa) const {

 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::Helm);
}


void HelmIsoParamHexa::addFaces(PolygonSet *pset) {

 int i;
 int corner[4];
 IsoParamUtils ipu(order);
 for(i=0;i<6;i++) {
   ipu.cornerindeces(i,corner); 
   pset->addQuad(this,nn[corner[0]],nn[corner[1]],nn[corner[2]],nn[corner[3]]);
 }
}


double HelmIsoParamHexa::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmIsoParamHexa::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmIsoParamHexa::massMatrix(const CoordSet &cs, double *K, int fl) const {

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
 int i;
 int j;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

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
   for (i = 0; i < order; ++i) {
     MD += ret[i][i];
     for (j = 0; j < order; ++j) {
       MM += ret[i][j];
       if (i != j)
         ret[i][j] =0.0;
     }
   }
   for (i = 0; i < order; ++i)
     ret[i][i] = ret[i][i]*MM/MD;
 }
 return ret;
}


FullSquareMatrix HelmIsoParamHexa::stiffness(const CoordSet &cs, double *K, int flg ) const {
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
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(orderc,K);
 return ret;
}


FullSquareMatrixC HelmIsoParamHexa::massMatrix(const CoordSet &cs,
                                               complex<double> *K) const {

 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 ipu.zeroOut<complex<double> > (orderc*orderc,K);
 int gorder = 7;
 if (order<=3) gorder = O3;
 else if (order<=4) gorder = O4;

 if (prop->fp.PMLtype==1) {
   double Rm[3] = { -prop->fp.Rx, -prop->fp.Ry, -prop->fp.Rz };
   double Rp[3] = { prop->fp.Rx, prop->fp.Ry, prop->fp.Rz };
   double Sm[3] = { -prop->fp.Sx, -prop->fp.Sy, -prop->fp.Sz };
   double Sp[3] = { prop->fp.Sx, prop->fp.Sy, prop->fp.Sz };
   CartPMLGalMassFunction f(orderc,Rm,Sm,Rp,Sp,prop->fp.gamma,K);
   ipu.volumeInt3d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==2) {
   SphPMLGalMassFunction f(orderc,prop->fp.Rx,prop->fp.Sx,prop->fp.gamma,K);
   ipu.volumeInt3d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==3) {
   CylPMLGalMassFunction f(orderc,prop->fp.Rx,prop->fp.Sx,
                        -prop->fp.Rz,prop->fp.Rz,-prop->fp.Sz,prop->fp.Sz,
                        prop->fp.gamma,K);
   ipu.volumeInt3d(xyz, f, gorder);
 } else {
   GalMassFunction<complex<double> > f(orderc,K);
   ipu.volumeInt3d(xyz, f, gorder);
 }

 ipu.symmetrize(orderc,K);
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrixC ret(orderc,K);
 return ret;
}


FullSquareMatrixC HelmIsoParamHexa::stiffness(const CoordSet &cs,
                                              complex<double> *K) const {

 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 ipu.zeroOut< complex<double> > (orderc*orderc,K);
 int gorder = 7;
 if (order<=3) gorder = O3;
 else if (order<=4) gorder = O4;

 if (prop->fp.PMLtype==1) {
   double Rm[3] = { -prop->fp.Rx, -prop->fp.Ry, -prop->fp.Rz };
   double Rp[3] = { prop->fp.Rx, prop->fp.Ry, prop->fp.Rz };
   double Sm[3] = { -prop->fp.Sx, -prop->fp.Sy, -prop->fp.Sz };
   double Sp[3] = { prop->fp.Sx, prop->fp.Sy, prop->fp.Sz };
   CartPMLGalStiffFunction f(orderc,Rm,Sm,Rp,Sp,prop->fp.gamma,K);
   ipu.volumeInt3d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==2) {
   SphPMLGalStiffFunction f(orderc,prop->fp.Rx,prop->fp.Sx,prop->fp.gamma,K);
   ipu.volumeInt3d(xyz, f, gorder);
 } else if (prop->fp.PMLtype==3) {
   CylPMLGalStiffFunction f(orderc,prop->fp.Rx,prop->fp.Sx,
                        -prop->fp.Rz,prop->fp.Rz,-prop->fp.Sz,prop->fp.Sz,
                        prop->fp.gamma,K);
   ipu.volumeInt3d(xyz, f, gorder);
 } else {
   GalStiffFunction<complex<double> > f(orderc,K);
   ipu.volumeInt3d(xyz, f, gorder);
 }

 ipu.symmetrize(orderc,K);
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrixC ret(orderc,K);
 return ret;
}


extern bool useFull;

int
HelmIsoParamHexa::numNodes() const {
  if(useFull)
    return order*order*order;
  else
    return(8);   // to ignore effect of mid-size nodes in dec
}



int HelmIsoParamHexa::getDecFace(int iFace, int *fn) {
  IsoParamUtils ipu(order);
  int ordersq = ipu.getordersq();
  ipu.faceindeces(iFace+1, fn);
  for(int i=0;i<ordersq;i++) {
    int tmp = fn[i];
    fn[i] = nn[tmp];
  }
  return ordersq;
}

double HelmIsoParamHexa::weight() const
{
	return order;
}

double HelmIsoParamHexa::trueWeight() const
{
	return order;
}
