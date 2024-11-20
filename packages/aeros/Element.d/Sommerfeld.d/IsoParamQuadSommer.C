#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Utils.d/dbg_alloca.h>
#include <Element.d/Sommerfeld.d/IsoParamQuadSommer.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Helm.d/HelmIsoParamHexa.h>
#include <Element.d/Helm.d/IsoParamUtils.h>

#define O3 4
#define O4 7

IsoParamQuadSommer::IsoParamQuadSommer(int o, int *nodenums, Element *_el, int etype) {

	order = int(rint(sqrt(double(o))));

	nn = new int[order * order];

	int i;
	for (i = 0; i < order * order; i++) nn[i] = nodenums[i];
	setElementType(etype);
	el = _el;
 sFlag = false;
 soundSpeed = 0.0;
}


IsoParamQuadSommer *IsoParamQuadSommer::clone() {
	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	IsoParamQuadSommer *se = new IsoParamQuadSommer(ordersq, nn, el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
	se->soundSpeed = soundSpeed;
	return se;
}


int *IsoParamQuadSommer::dofs(DofSetArray &dsa, int *p) const  {

	if (p == 0) p = new int[order * order];
	int i;
	for (i = 0; i < order * order; i++) dsa.number(nn[i], DofSet::Helm, p + i);
	return p;
}


int *IsoParamQuadSommer::solidDofs(DofSetArray &dsa, int *p) const {

	if (p == 0) p = new int[3 * order * order];
	int i;
	for (i = 0; i < order * order; i++) dsa.number(nn[i], DofSet::Xdisp | DofSet::Ydisp | DofSet::Zdisp, p + 3 * i);
	return p;
}


int *IsoParamQuadSommer::wetDofs(DofSetArray &dsa, int *p) const {

	if (p == 0) p = new int[order * order * 4];
	int i;
	for (i = 0; i < order * order; i++)
		dsa.number(nn[i],
		           DofSet::Helm | DofSet::Xdisp | DofSet::Ydisp | DofSet::Zdisp, p + 4 * i);

	return p;
}


void IsoParamQuadSommer::flipNormal() {
	int *nds = (int *) dbg_alloca(sizeof(int) * numNodes());
	nodes(nds);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) nn[i * order + j] = nds[j * order + i];
}

//#define PRTEST
#ifdef PRTEST
void IsoParamQuadSommer::neumVector(CoordSet& cs, ComplexVector& cv,
	double omega, double dx, double dy, double dz,int pflag) {

 if (el==0) {
   fprintf(stderr,"IsoParamQuadSommer::neumVector : adjacent element not defined.\n");
   exit(-1);
 }

 IsoParamUtils ipu(order);
 int ordersq = ipu.getordersq();
 double *xyz=(double*)alloca(sizeof(double)*3*ordersq);
 cs.getCoordinates(nn,ordersq,xyz,xyz+ordersq,xyz+2*ordersq);

 double dir[3] = { dx,dy, dz};
 double E = el->getProperty()->E;
 double nu = el->getProperty()->nu;
 double rho = el->getProperty()->rho;

 LENeumInterfaceBCGalFunction f(3*ordersq,omega,E,nu,rho,dir,cv.data());
 ipu.zeroOut<complex<double> > (3*ordersq,cv.data());
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.surfSurfInt3d(xyz, f, gorder);
}
#else

void IsoParamQuadSommer::neumVector(CoordSet &cs, ComplexVector &cv,
                                    double kappa, double dx, double dy, double dz, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::neumVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	HelmIsoParamHexa *e = (HelmIsoParamHexa *) el;
	int orderc = ipu.getorderc();
	int *nds = (int *) dbg_alloca(sizeof(int) * orderc);
	e->nodes(nds);

	int ordersq = ipu.getordersq();
	int faceindex;
	int *map = (int *) dbg_alloca(sizeof(int) * ordersq);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) dbg_alloca(sizeof(double) * 3 * orderc);
	cs.getCoordinates(nds, orderc, xyz, xyz + orderc, xyz + 2 * orderc);

	double incdir[3] = {dx, dy, dz};
	complex<double> *v = (complex<double> *) dbg_alloca(sizeof(complex<double>) * orderc);

	NeumannBCGalFunction f(orderc, kappa, incdir, pflag, v);
	ipu.zeroOut<complex<double> >(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfInt3d(xyz, faceindex, f, gorder);

	int i;
	for (i = 0; i < ordersq; i++) {
		cv[i] = v[map[i]] / el->getProperty()->rho;
	}
}

#endif


void IsoParamQuadSommer::neumVectorDeriv(CoordSet &cs, ComplexVector &cv,
                                         double kappa, double dx, double dy, double dz, int dero, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::neumVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	HelmIsoParamHexa *e = (HelmIsoParamHexa *) el;
	int orderc = ipu.getorderc();
	int *nds = (int *) alloca(sizeof(int) * orderc);
	e->nodes(nds);

	int ordersq = ipu.getordersq();
	int faceindex;
	int *map = (int *) alloca(sizeof(int) * ordersq);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) alloca(sizeof(double) * 3 * orderc);
	cs.getCoordinates(nds, orderc, xyz, xyz + orderc, xyz + 2 * orderc);

	double incdir[3] = {dx, dy, dz};
	complex<double> *v = (complex<double> *) alloca(sizeof(complex<double>) * orderc);

	NeumannBCGalFunction f(orderc, kappa, incdir, pflag, v, dero);
	ipu.zeroOut<complex<double> >(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfInt3d(xyz, faceindex, f, gorder);

	int i;
	for (i = 0; i < ordersq; i++) {
		cv[i] = v[map[i]] / el->getProperty()->rho;
	}
}


void
IsoParamQuadSommer::wetInterfaceVector(CoordSet &cs, ComplexVector &cv, double kappa, double dx, double dy, double dz,
                                       int dero,
                                       int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::wetInterfaceVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double dir[3] = {dx, dy, dz};
	complex<double> *tcv =
			(complex<double> *) alloca(sizeof(complex<double>) * 4 * ordersq);

	WetInterfaceBCGalFunction f(ordersq, kappa, dir, pflag, tcv, dero);
	ipu.zeroOut<complex<double> >(4 * ordersq, tcv);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfSurfInt3d(xyz, f, gorder);

	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0 && el2 == 0) {
		int j = 0;
		for (int i = 0; i < ordersq; i++) {
			cv[j++] = tcv[4 * i + 0];
			cv[j++] = tcv[4 * i + 1];
			cv[j++] = tcv[4 * i + 2];
		}
	} else if (he != 0 && el2 == 0) {
		for (int i = 0; i < ordersq; i++) cv[i] = tcv[4 * i + 3] / el->getProperty()->rho;
	} else if (he == 0) {
		for (int i = 0; i < ordersq; i++) {
			cv[4 * i + 0] = tcv[4 * i + 0];
			cv[4 * i + 1] = tcv[4 * i + 1];
			cv[4 * i + 2] = tcv[4 * i + 2];
			cv[4 * i + 3] = tcv[4 * i + 3] / el2->getProperty()->rho;
		}
	} else {
		for (int i = 0; i < ordersq; i++) {
			cv[4 * i + 0] = tcv[4 * i + 0];
			cv[4 * i + 1] = tcv[4 * i + 1];
			cv[4 * i + 2] = tcv[4 * i + 2];
			cv[4 * i + 3] = tcv[4 * i + 3] / el->getProperty()->rho;
		}
	}
}


void IsoParamQuadSommer::wetInterfaceVector(CoordSet &cs, ComplexVector &cv,
                                            complex<double> (*diri)[3], complex<double> *coefi) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::wetInterfaceVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	int *nodes = (int *) alloca(sizeof(int) * el->numNodes());
	el->nodes(nodes);
	double z = 0.0;
	int i;
	for (i = 0; i < el->numNodes(); i++) {
		Node nd = cs.getNode(nodes[i]);
		z += nd.z;
	}

	int imat;
	if (z < 0.0) imat = 1;
	else imat = 2;
	if (coefi[0] == complex<double>(0.0, 0.0) &&
	    coefi[1] == complex<double>(0.0, 0.0))
		imat = 2;

	complex<double> *tcv = (complex<double> *) alloca(sizeof(complex<double>) * 4 * ordersq);
	ipu.zeroOut<complex<double> >(4 * ordersq, tcv);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;

	WetInterfaceHeteroBCGalFunction f(ordersq, imat, diri, coefi, tcv);
	ipu.surfSurfInt3d(xyz, f, gorder);

	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0 && el2 == 0) {
		int j = 0;
		for (i = 0; i < ordersq; i++) {
			cv[j++] = tcv[4 * i + 0];
			cv[j++] = tcv[4 * i + 1];
			cv[j++] = tcv[4 * i + 2];
		}
	} else if (he != 0 && el2 == 0) {
		for (i = 0; i < ordersq; i++) cv[i] = tcv[4 * i + 3] / el->getProperty()->rho;
	} else if (he == 0) {
		for (i = 0; i < ordersq; i++) {
			cv[4 * i + 0] = tcv[4 * i + 0];
			cv[4 * i + 1] = tcv[4 * i + 1];
			cv[4 * i + 2] = tcv[4 * i + 2];
			cv[4 * i + 3] = tcv[4 * i + 3] / el2->getProperty()->rho;
		}
	} else {
		for (i = 0; i < ordersq; i++) {
			cv[4 * i + 0] = tcv[4 * i + 0];
			cv[4 * i + 1] = tcv[4 * i + 1];
			cv[4 * i + 2] = tcv[4 * i + 2];
			cv[4 * i + 3] = tcv[4 * i + 3] / el->getProperty()->rho;
		}
	}
}


void IsoParamQuadSommer::wetInterfaceVectorDeriv(CoordSet &cs,
                                                 ComplexVector &cv, complex<double> (*diri)[3],
                                                 complex<double> *coefi, complex<double> *kappai, int dero) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::wetInterfaceVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	int *nodes = (int *) alloca(sizeof(int) * el->numNodes());
	el->nodes(nodes);
	double z = 0.0;
	int i;
	for (i = 0; i < el->numNodes(); i++) {
		Node nd = cs.getNode(nodes[i]);
		z += nd.z;
	}

	int imat;
	if (z < 0.0) imat = 1;
	else imat = 2;
	if (coefi[0] == complex<double>(0.0, 0.0) &&
	    coefi[1] == complex<double>(0.0, 0.0))
		imat = 2;

	WetInterfaceHeteroBCDerivGalFunction f(ordersq, imat, diri, coefi, kappai,
	                                       dero, cv.data());
	ipu.zeroOut<complex<double> >(4 * ordersq, cv.data());
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfSurfInt3d(xyz, f, gorder);

	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0 && el2 == 0) {
		int j = 0;
		for (i = 0; i < ordersq; i++) {
			cv[j++] = cv[4 * i + 0];
			cv[j++] = cv[4 * i + 1];
			cv[j++] = cv[4 * i + 2];
		}
	} else if (he != 0 && el2 == 0) {
		for (i = 0; i < ordersq; i++) cv[i] = cv[4 * i + 3] / el->getProperty()->rho;
	} else if (he == 0) {
		for (i = 0; i < ordersq; i++) cv[4 * i + 3] /= el2->getProperty()->rho;
	} else {
		for (i = 0; i < ordersq; i++) cv[4 * i + 3] /= el->getProperty()->rho;
	}

}


/*
FullSquareMatrix IsoParamQuadSommer::sommerMatrix(CoordSet &cs, double *d) {

 if (el==0) { 
   fprintf(stderr,"IsoParamQuadSommer::sommerMatrix: adjacent element not defined.\n");
   exit(-1);
 }
 
 IsoParamUtils ipu(order);
 HelmIsoParamHexa *e = (HelmIsoParamHexa*)el;
 int orderc = ipu.getorderc();
 int *nds =(int*)dbg_alloca(sizeof(int)*orderc);
 e->nodes(nds);

 int ordersq = ipu.getordersq();
 int faceindex;
 int *map = (int*)dbg_alloca(sizeof(int)*ordersq);
 ipu.facemap(faceindex,nn,nds,map);

 double *xyz=(double*)dbg_alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nds,orderc,xyz,xyz+orderc,xyz+2*orderc);

 double* v = (double*)dbg_alloca(sizeof(double)*orderc*orderc);
 SomGalFunction f(orderc);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),v);
 int gorder = 7;
 if (order<=3) gorder = 4;
 ipu.surfInt3d<SomGalFunction>(xyz, faceindex, f, gorder);
 ipu.symmetrize(f.ncolumns(),v);

 FullSquareMatrix sommerM(ordersq,d);

 int i,j;
 for(i=0;i<ordersq;i++) for(j=0;j<ordersq;j++) {
   sommerM[i][j] = -v[map[j]*orderc+map[i]];
 }

 return sommerM;
}*/


FullSquareMatrix IsoParamQuadSommer::sommerMatrix(CoordSet &cs, double *d) const {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *v = (double *) alloca(sizeof(double) * ordersq * ordersq);
	SomGalFunction f(ordersq, v);
	ipu.zeroOut<double>(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfSurfInt3d(xyz, f, gorder);
	ipu.symmetrize(f.ncolumns(), v);

	FullSquareMatrix sommerM(ordersq, d);
	int i, j;
	for (i = 0; i < ordersq; i++)
		for (j = 0; j < ordersq; j++) {
			sommerM[i][j] = -v[j * ordersq + i] / el->getProperty()->rho;
		}

	return sommerM;
}


GenStackFSFullMatrix<double> IsoParamQuadSommer::wetInterfaceMatrix(
		CoordSet &cs, double *d) {

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	WetInterfaceGalFunction f(ordersq, d);
	ipu.zeroOut<double>(3 * ordersq * ordersq, d);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfSurfInt3d(xyz, f, gorder);

	GenStackFSFullMatrix<double> sommerM(ordersq, 3 * ordersq, d);
	return sommerM;
}


void IsoParamQuadSommer::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {
	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *d = (double *) alloca(sizeof(double) * 3 * ordersq * ordersq);
	WetInterfaceGalFunction f(ordersq, d);
	ipu.zeroOut<double>(3 * ordersq * ordersq, d);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfSurfInt3d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < ordersq; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in IsoParamQuadSommer::wetInterfaceLMPC\n");
		return;
	}

 if (!sFlag)
	for (i = 0; i < ordersq; i++) {
		LMPCTerm lmpct1(nn[i], 0, -d[i * ordersq + j]);
		lmpc->addterm(&lmpct1);
		LMPCTerm lmpct2(nn[i], 1, -d[ordersq * ordersq + i * ordersq + j]);
		lmpc->addterm(&lmpct2);
		LMPCTerm lmpct3(nn[i], 2, -d[2 * ordersq * ordersq + i * ordersq + j]);
		lmpc->addterm(&lmpct3);
	}
 else
 for(i=0;i<ordersq;i++) {
   LMPCTerm lmpct1(nn[i],0, 0.0 );
   lmpc->addterm(&lmpct1);
   LMPCTerm lmpct2(nn[i],1, 0.0 );
   lmpc->addterm(&lmpct2);
   LMPCTerm lmpct3(nn[i],2, 0.0 );
   lmpc->addterm(&lmpct3);
 } 
}


void IsoParamQuadSommer::ffp(CoordSet &cs, int numFFP, double *dirFFP,
                             complex<double> *sol, complex<double> *ffpv,
                             bool direction) {
	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::ffp: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	HelmIsoParamHexa *e = (HelmIsoParamHexa *) el;
	int orderc = ipu.getorderc();
	int *nds = (int *) alloca(sizeof(int) * orderc);
	e->nodes(nds);

	int ordersq = ipu.getordersq();
	int faceindex;
	int *map = (int *) alloca(sizeof(int) * ordersq);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) alloca(sizeof(double) * 3 * orderc);
	cs.getCoordinates(nds, orderc, xyz, xyz + orderc, xyz + 2 * orderc);

	double kappa = el->getProperty()->kappaHelm;

	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;

	if (direction) {
		FFPGalFunction f(orderc, kappa, numFFP, dirFFP, sol, ffpv);
		ipu.surfGInt3d(xyz, faceindex, f, gorder);
	} else {
		KirchhoffGalFunction f(orderc, kappa, numFFP, dirFFP, sol, ffpv);
		ipu.surfGInt3d(xyz, faceindex, f, gorder);
	}

}

FullSquareMatrix IsoParamQuadSommer::turkelMatrix(CoordSet &cs, double *d) const {

	FullSquareMatrix sommerM(order * order, d);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;
	fprintf(stderr, "IsoParamQuadSommer::turkelMatrix not implemented\n");
	return sommerM;
}


FullSquareMatrixC IsoParamQuadSommer::sommer2Matrix(CoordSet &cs,
                                                    complex<double> *d) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::sommer2Matrix: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	HelmIsoParamHexa *e = (HelmIsoParamHexa *) el;
	int orderc = ipu.getorderc();
	int *nds = (int *) dbg_alloca(sizeof(int) * orderc);
	e->nodes(nds);

	int ordersq = ipu.getordersq();
	int faceindex;
	int *map = (int *) dbg_alloca(sizeof(int) * ordersq);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) dbg_alloca(sizeof(double) * 3 * orderc);
	cs.getCoordinates(nds, orderc, xyz, xyz + orderc, xyz + 2 * orderc);

	double kappa = el->getProperty()->kappaHelm;

	complex<double> *v = (complex<double> *) dbg_alloca(sizeof(complex<double>) *
	                                                    orderc * orderc);
	Som2GalFunction f(orderc, kappa, v);
	ipu.zeroOut<complex<double> >(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = O3;
	else if (order <= 4) gorder = O4;
	ipu.surfCurvInt3d(xyz, faceindex, f, gorder);
	ipu.symmetrize(f.ncolumns(), v);

	FullSquareMatrixC sommerM(ordersq, d);

	int i, j;
	for (i = 0; i < ordersq; i++)
		for (j = 0; j < ordersq; j++) {
			sommerM[i][j] = v[map[j] * orderc + map[i]] / el->getProperty()->rho;
		}

	return sommerM;
}


void IsoParamQuadSommer::getNormal(const CoordSet &cs, double normal[3]) const {

	int corner[4] = {1 - 1, order - 1, order * order - 1, order * order - 1 - (order - 1)};
	double x[4], y[4], z[4];

	Node nd1 = cs.getNode(nn[corner[0]]);
	Node nd2 = cs.getNode(nn[corner[1]]);
	Node nd3 = cs.getNode(nn[corner[2]]);
	Node nd4 = cs.getNode(nn[corner[3]]);

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;
	x[2] = nd3.x;
	y[2] = nd3.y;
	z[2] = nd3.z;
	x[3] = nd4.x;
	y[3] = nd4.y;
	z[3] = nd4.z;

	normal[0] = (y[2] - y[0]) * (z[3] - z[1]) - (z[2] - z[0]) * (y[3] - y[1]);
	normal[1] = (z[2] - z[0]) * (x[3] - x[1]) - (x[2] - x[0]) * (z[3] - z[1]);
	normal[2] = (x[2] - x[0]) * (y[3] - y[1]) - (y[2] - y[0]) * (x[3] - x[1]);

	double l = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

	normal[0] = normal[0] / l;
	normal[1] = normal[1] / l;
	normal[2] = normal[2] / l;
}


void
IsoParamQuadSommer::markDofs(DofSetArray &dsa) const {
	for (int i = 0; i < order * order; i++) dsa.mark(nn[i], DofSet::Helm);
}
