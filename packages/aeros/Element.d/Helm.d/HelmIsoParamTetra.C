#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/HelmIsoParamTetra.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/PML.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/dbg_alloca.h>


HelmIsoParamTetra::HelmIsoParamTetra(int o, int* nodenums) {

 if (o==4) order = 2;
 else if (o==10) order = 3;
 else if (o==20) order = 4;
 else if (o==35) order = 5;
 else if (o==56) order = 6;
 else {
   fprintf(stderr,"Order too high in HelmIsoParamTetra::HelmIsoParamTetra.\n");
   exit(-1);
 }
 
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


HelmIsoParamTetra::HelmIsoParamTetra(const HelmIsoParamTetra& e) {
 order = e.order;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * HelmIsoParamTetra::clone() {
 return new HelmIsoParamTetra(*this);
}


void HelmIsoParamTetra::renum(const int *table) {
 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

void HelmIsoParamTetra::renum(EleRenumMap& table) {
 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}


int* HelmIsoParamTetra::nodes(int *p) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 if(p == 0) p = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) p[i] = nn[i];
 return p;
}


int* HelmIsoParamTetra::dofs(DofSetArray &dsa, int *p) const  {

 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 if(p == 0) p = new int[orderc];
 for(i=0;i<orderc;i++) dsa.number(nn[i],DofSet::Helm,p+i);
 return p;
}


void HelmIsoParamTetra::markDofs(DofSetArray &dsa) const {

 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::Helm);
}


void HelmIsoParamTetra::addFaces(PolygonSet *pset) {

 int i;
 int corner[3];
 IsoParamUtilsTetra ipu(order);
 for(i=0;i<4;i++) {
   ipu.cornerindeces(i,corner); 
   pset->addTri(this,nn[corner[0]],nn[corner[1]],nn[corner[2]]);
 }
}


double HelmIsoParamTetra::getMass(const CoordSet&) const {
 fprintf(stderr,"HelmIsoParamTetra::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix HelmIsoParamTetra::massMatrix(const CoordSet &cs, double *K, int fl) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 GalMassFunction<double> f(orderc,K);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(f.ncolumns(),K);
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(orderc,K);
 return ret;
}


FullSquareMatrix HelmIsoParamTetra::stiffness(const CoordSet &cs, double *K, int flg ) const {

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
 int i;
 for(i=0;i<orderc*orderc;i++) K[i] /= prop->rho;

 FullSquareMatrix ret(orderc,K);
 return ret;
}


FullSquareMatrixC HelmIsoParamTetra::massMatrix(const CoordSet &cs,
                                               complex<double> *K) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 ipu.zeroOut<complex<double> > (orderc*orderc,K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;

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


FullSquareMatrixC HelmIsoParamTetra::stiffness(const CoordSet &cs,
                                              complex<double> *K) const {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 ipu.zeroOut< complex<double> > (orderc*orderc,K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;

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
HelmIsoParamTetra::numNodes() const {
  //Not tested -JF
  if(useFull)
    return (order*(order+1)*(order+2))/6;
  else
    return(4);   // to ignore effect of mid-size nodes in dec
}

int HelmIsoParamTetra::getDecFace(int iFace, int *fn) {
  IsoParamUtilsTetra ipu(order);
  int ordersq = ipu.getordersq();
  ipu.faceindeces(iFace+1, fn);
  for(int i=0;i<ordersq;i++) {
    int tmp = fn[i];
    fn[i] = nn[tmp];
  }
  return ordersq;
}

double HelmIsoParamTetra::weight() const
{
	return order-1;
}

double HelmIsoParamTetra::trueWeight() const
{
	return order-1;
}


