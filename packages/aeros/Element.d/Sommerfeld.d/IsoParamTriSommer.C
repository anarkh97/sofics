#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Utils.d/dbg_alloca.h>
#include <Element.d/Sommerfeld.d/IsoParamTriSommer.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Helm.d/HelmIsoParamTetra.h>
#include <Element.d/Helm.d/IsoParamUtils.h>


IsoParamTriSommer::IsoParamTriSommer(int o, int *nodenums, Element *_el, int etype) {

	if (o == 3) order = 2;
	else if (o == 6) order = 3;
	else if (o == 10) order = 4;
	else if (o == 15) order = 5;
	else if (o == 21) order = 6;
	else {
		fprintf(stderr, "Order too high in IsoParamTriSommer::IsoParamTriSommer\n");
		exit(-1);
	}

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	nn = new int[ordersq];

	int i;
	for (i = 0; i < ordersq; i++) nn[i] = nodenums[i];
	setElementType(etype);
	el = _el;
 sFlag = false;
	soundSpeed = 0.0;
}


IsoParamTriSommer *IsoParamTriSommer::clone() {
	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	IsoParamTriSommer *se = new IsoParamTriSommer(ordersq, nn, el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
	se->soundSpeed = soundSpeed;
	return se;
}


int *IsoParamTriSommer::dofs(DofSetArray &dsa, int *p) const  {

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	if (p == 0) p = new int[ordersq];
	int i;
	for (i = 0; i < ordersq; i++) dsa.number(nn[i], DofSet::Helm, p + i);
	return p;
}


int *IsoParamTriSommer::solidDofs(DofSetArray &dsa, int *p) const {

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	if (p == 0) p = new int[3 * ordersq];
	int i;
	for (i = 0; i < ordersq; i++)
		dsa.number(nn[i],
		           DofSet::Xdisp | DofSet::Ydisp | DofSet::Zdisp, p + 3 * i);
	return p;
}


int *IsoParamTriSommer::wetDofs(DofSetArray &dsa, int *p) const {

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	if (p == 0) p = new int[ordersq * 4];
	int i;
	for (i = 0; i < ordersq; i++)
		dsa.number(nn[i],
		           DofSet::Helm | DofSet::Xdisp | DofSet::Ydisp | DofSet::Zdisp, p + 4 * i);
	return p;
}


void IsoParamTriSommer::flipNormal() {
	int *nds = (int *) dbg_alloca(sizeof(int) * numNodes());
	nodes(nds);
	int *map = (int *) dbg_alloca(sizeof(int) * order * order);
	int c = 0;
	int r, s;
	for (s = 0; s < order; s++) {
		for (r = 0; r < order - s; r++) {
			map[r * order + s] = c;
			c++;
		}
	}
	c = 0;
	for (r = 0; r < order; r++) {
		for (s = 0; s < order - r; s++) {
			nn[c] = nds[map[r * order + s]];
			c++;
		}
	}
}


void IsoParamTriSommer::neumVector(CoordSet &cs, ComplexVector &cv,
                                   double kappa, double dx, double dy, double dz, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamTriSommer::neumVector: Adjacent element not defined.\n");
		return;
	}

	HelmIsoParamTetra *e = (HelmIsoParamTetra *) el;
	IsoParamUtilsTetra ipu(order);
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
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
	ipu.surfInt3d(xyz, faceindex, f, gorder);

	int i;
	for (i = 0; i < ordersq; i++) {
		cv[i] = v[map[i]] / el->getProperty()->rho;
	}
}

void IsoParamTriSommer::neumVectorDeriv(CoordSet &cs, ComplexVector &cv,
                                        double kappa, double dx, double dy, double dz, int dero, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamQuadSommer::neumVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtilsTetra ipu(order);
	HelmIsoParamTetra *e = (HelmIsoParamTetra *) el;
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
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
	ipu.surfInt3d(xyz, faceindex, f, gorder);

	int i;
	for (i = 0; i < ordersq; i++) {
		cv[i] = v[map[i]] / el->getProperty()->rho;
	}
}


void IsoParamTriSommer::wetInterfaceVector(CoordSet &cs, ComplexVector &cv,
                                           double kappa, double dx, double dy, double dz,
                                           int dero, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamTriSommer::wetInterfaceVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double dir[3] = {dx, dy, dz};
	complex<double> *tcv = (complex<double> *) alloca(sizeof(complex<double>) * 4 * ordersq);

	WetInterfaceBCGalFunction f(ordersq, kappa, dir, pflag, tcv, dero);
	ipu.zeroOut<complex<double> >(4 * ordersq, tcv);
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
	ipu.surfSurfInt3d(xyz, f, gorder);

	int i;
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


void IsoParamTriSommer::wetInterfaceVector(CoordSet &cs, ComplexVector &cv,
                                           complex<double> (*diri)[3], complex<double> *coefi) {

	if (el == 0) {
		fprintf(stderr, "IsoParamTriSommer::wetInterfaceVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtilsTetra ipu(order);
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

	WetInterfaceHeteroBCGalFunction f(ordersq, imat, diri, coefi, cv.data());
	ipu.zeroOut<complex<double> >(4 * ordersq, cv.data());
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
	ipu.surfSurfInt3d(xyz, f, gorder);
	for (i = 0; i < ordersq; i++) cv[4 * i + 3] /= el->getProperty()->rho;
}

//#define _USE3D_EL
#ifdef _USE3D_EL
FullSquareMatrix IsoParamTriSommer::sommerMatrix(CoordSet &cs, double *d) const {

 if (el==0) {
   fprintf(stderr,"IsoParamTriSommer::sommerMatrix: adjacent element not defined.\n");
   exit(-1);
 }

 HelmIsoParamTetra *e = (HelmIsoParamTetra*)el;
 IsoParamUtilsTetra ipu(order);
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
 SomGalFunction f(orderc,v);
 ipu.zeroOut<double> (f.nrows()*f.ncolumns(),v);
 int gorder = 7*7;
 if (order<=3) gorder = 13;
 ipu.surfInt3d(xyz, faceindex, f, gorder);
 ipu.symmetrize(f.ncolumns(),v);

 FullSquareMatrix sommerM(ordersq,d);

 int i,j;
 for(i=0;i<ordersq;i++) for(j=0;j<ordersq;j++) {
   sommerM[i][j] = -v[map[j]*orderc+map[i]] / el->getProperty()->rho;
 }

 return sommerM;
}

#else

FullSquareMatrix IsoParamTriSommer::sommerMatrix(CoordSet &cs, double *d) const {

	if (el == 0) {
		fprintf(stderr, "IsoParamTriSommer::sommerMatrix: Adjacent element not defined.\n");
		return FullSquareMatrix(0);
	}

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *v = (double *) alloca(sizeof(double) * ordersq * ordersq);
	SomGalFunction f(ordersq, v);
	ipu.zeroOut<double>(f.nrows() * f.ncolumns(), v);
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
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

#endif


GenStackFSFullMatrix<double> IsoParamTriSommer::wetInterfaceMatrix(
		CoordSet &cs, double *d) {

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	WetInterfaceGalFunction f(ordersq, d);
	ipu.zeroOut<double>(3 * ordersq * ordersq, d);
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
	ipu.surfSurfInt3d(xyz, f, gorder);

	GenStackFSFullMatrix<double> sommerM(ordersq, 3 * ordersq, d);
	return sommerM;
}


void IsoParamTriSommer::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {
	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *d = (double *) alloca(sizeof(double) * 3 * ordersq * ordersq);
	WetInterfaceGalFunction f(ordersq, d);
	ipu.zeroOut<double>(3 * ordersq * ordersq, d);
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
	ipu.surfSurfInt3d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < ordersq; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in IsoParamTriSommer::wetInterfaceLMPC\n");
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
   LMPCTerm lmpct1(nn[i],0, 0.0);
   lmpc->addterm(&lmpct1);
   LMPCTerm lmpct2(nn[i],1, 0.0);
   lmpc->addterm(&lmpct2);
   LMPCTerm lmpct3(nn[i],2, 0.0);
		lmpc->addterm(&lmpct3);
	}
}


void IsoParamTriSommer::ffp(CoordSet &cs, int numFFP, double *dirFFP,
                            complex<double> *sol, complex<double> *ffpv,
                            bool direction) {
	if (el == 0) {
		fprintf(stderr, "IsoParamTriSommer::ffp: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtilsTetra ipu(order);
	HelmIsoParamTetra *e = (HelmIsoParamTetra *) el;
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

	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;

	if (direction) {
		FFPGalFunction f(orderc, kappa, numFFP, dirFFP, sol, ffpv);
		ipu.surfGInt3d(xyz, faceindex, f, gorder);
	} else {
		KirchhoffGalFunction f(orderc, kappa, numFFP, dirFFP, sol, ffpv);
		ipu.surfGInt3d(xyz, faceindex, f, gorder);
	}

}


FullSquareMatrixC IsoParamTriSommer::sommer2Matrix(CoordSet &cs,
                                                   complex<double> *d) {

	if (el == 0) {
		fprintf(stderr, "IsoParamTriSommer::sommer2Matrix: adjacent element not defined.\n");
		exit(-1);
	}

	HelmIsoParamTetra *e = (HelmIsoParamTetra *) el;
	IsoParamUtilsTetra ipu(order);
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
	int gorder = 7 * 7;
	if (order <= 3) gorder = 13;
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


FullSquareMatrix IsoParamTriSommer::turkelMatrix(CoordSet &cs, double *d) const {

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();

	FullSquareMatrix sommerM(ordersq, d);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;
	fprintf(stderr, "IsoParamTriSommer::turkelMatrix not implemented\n");
	return sommerM;
}


void IsoParamTriSommer::getNormal(const CoordSet &cs, double normal[3]) const {

	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	int corner[3] = {1 - 1, order - 1, ordersq - 1};
	double x[3], y[3], z[3];

	Node nd1 = cs.getNode(nn[corner[0]]);
	Node nd2 = cs.getNode(nn[corner[1]]);
	Node nd3 = cs.getNode(nn[corner[2]]);

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;
	x[2] = nd3.x;
	y[2] = nd3.y;
	z[2] = nd3.z;

	normal[0] = (y[2] - y[0]) * (z[1] - z[0]) - (z[2] - z[0]) * (y[1] - y[0]);
	normal[1] = (z[2] - z[0]) * (x[1] - x[0]) - (x[2] - x[0]) * (z[1] - z[0]);
	normal[2] = (x[2] - x[0]) * (y[1] - y[0]) - (y[2] - y[0]) * (x[1] - x[0]);

	double l = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

	normal[0] = normal[0] / l;
	normal[1] = normal[1] / l;
	normal[2] = normal[2] / l;
}

void
IsoParamTriSommer::markDofs(DofSetArray &dsa) const {
	IsoParamUtilsTetra ipu(order);
	int ordersq = ipu.getordersq();
	for (int i = 0; i < ordersq; i++) dsa.mark(nn[i], DofSet::Helm);
}
