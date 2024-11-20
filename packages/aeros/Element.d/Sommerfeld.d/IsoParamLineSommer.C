#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <alloca.h>
#include <Element.d/Sommerfeld.d/IsoParamLineSommer.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Helm.d/HelmIsoParamQuad.h>
#include <Element.d/Helm.d/IsoParamUtils2d.h>


IsoParamLineSommer::IsoParamLineSommer(int o, int *nodenums, Element *_el) {

	order = o;

	nn = new int[order];

	int i;
	for (i = 0; i < order; i++) nn[i] = nodenums[i];
	el = _el;
 sFlag = false;
	soundSpeed = 0.0;
}


IsoParamLineSommer *IsoParamLineSommer::clone() {
	IsoParamLineSommer *se = new IsoParamLineSommer(order, nn, el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
	se->soundSpeed = soundSpeed;
	return se;
}


int *IsoParamLineSommer::dofs(DofSetArray &dsa, int *p) const  {

	if (p == 0) p = new int[order];
	int i;
	for (i = 0; i < order; i++) dsa.number(nn[i], DofSet::Helm, p + i);
	return p;
}


int *IsoParamLineSommer::wetDofs(DofSetArray &dsa, int *p) const {

	if (p == 0) p = new int[order * 3];
	int i;
	for (i = 0; i < order; i++)
		dsa.number(nn[i],
		           DofSet::Helm | DofSet::Xdisp | DofSet::Ydisp, p + 3 * i);
	return p;
}


void IsoParamLineSommer::neumVector(CoordSet &cs, ComplexVector &cv,
                                    double kappa, double dx, double dy, double dz, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamLineSommer::neumVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils2d ipu(order);
	HelmIsoParamQuad *e = (HelmIsoParamQuad *) el;
	int ordersq = ipu.getordersq();
	int *nds = (int *) alloca(sizeof(int) * ordersq);
	e->nodes(nds);

	int faceindex;
	int *map = (int *) alloca(sizeof(int) * order);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nds, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double incdir[3] = {dx, dy, dz};
	complex<double> *v = (complex<double> *) alloca(sizeof(complex<double>) * ordersq);

	NeumannBCGalFunction2d f(ordersq, kappa, incdir, v);
	ipu.zeroOut<complex<double> >(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = 4;
	ipu.lineInt2d(xyz, faceindex, f, gorder);

	int i;
	for (i = 0; i < order; i++) {
		cv[i] = v[map[i]] / e->getProperty()->rho;
	}
}


void IsoParamLineSommer::neumVectorDeriv(CoordSet &cs, ComplexVector &cv,
                                         double kappa, double dx, double dy, double dz, int dero, int pflag) {

	if (el == 0) {
		fprintf(stderr, "IsoParamLineSommer::neumVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils2d ipu(order);
	HelmIsoParamQuad *e = (HelmIsoParamQuad *) el;
	int ordersq = ipu.getordersq();
	int *nds = (int *) alloca(sizeof(int) * ordersq);
	e->nodes(nds);

	int faceindex;
	int *map = (int *) alloca(sizeof(int) * order);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nds, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double incdir[3] = {dx, dy, dz};
	complex<double> *v = (complex<double> *) alloca(sizeof(complex<double>) * ordersq);

	NeumannBCGalFunction2d f(ordersq, kappa, incdir, v, dero);
	ipu.zeroOut<complex<double> >(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = 4;
	ipu.lineInt2d(xyz, faceindex, f, gorder);

	int i;
	for (i = 0; i < order; i++) {
		cv[i] = v[map[i]] / e->getProperty()->rho;
	}
}


void IsoParamLineSommer::wetInterfaceVector(CoordSet &cs, ComplexVector &cv,
                                            double kappa, double dx, double dy, double dz,
                                            int dero, int pflag) {

	if (el == 0) {
		fprintf(stderr,
		        "IsoParamLineSommer::wetInterfaceVector: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils2d ipu(order);
	double *xyz = (double *) alloca(sizeof(double) * 3 * order);
	cs.getCoordinates(nn, order, xyz, xyz + order, xyz + 2 * order);

	double dir[2] = {dx, dy};
	complex<double> *tcv =
			(complex<double> *) alloca(sizeof(complex<double>) * 3 * order);

	WetInterfaceBCGalFunction2d f(order, kappa, dir, cv.data(), dero);
	ipu.zeroOut<complex<double> >(3 * order, tcv);
	int gorder = 7;
	if (order <= 3) gorder = 4;
	ipu.lineLineInt2d(xyz, f, gorder);

	int i;
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0 && el2 == 0) {
		int j = 0;
		for (i = 0; i < order; i++) {
			cv[j++] = tcv[3 * i + 0];
			cv[j++] = tcv[3 * i + 1];
		}
	} else if (he != 0 && el2 == 0) {
		for (i = 0; i < order; i++) cv[i] = tcv[3 * i + 2] / el->getProperty()->rho;
	} else if (he == 0) {
		for (i = 0; i < order; i++) {
			cv[3 * i + 0] = tcv[3 * i + 0];
			cv[3 * i + 1] = tcv[3 * i + 1];
			cv[3 * i + 2] = tcv[3 * i + 2] / el2->getProperty()->rho;
		}
	} else {
		for (i = 0; i < order; i++) {
			cv[3 * i + 0] = tcv[3 * i + 0];
			cv[3 * i + 1] = tcv[3 * i + 1];
			cv[3 * i + 2] = tcv[3 * i + 2] / el->getProperty()->rho;
		}
	}

}


FullSquareMatrix IsoParamLineSommer::sommerMatrix(CoordSet &cs, double *d) const {

	if (el == 0) {
		fprintf(stderr, "IsoParamLineSommer::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils2d ipu(order);
	HelmIsoParamQuad *e = (HelmIsoParamQuad *) el;
	int ordersq = ipu.getordersq();
	int *nds = (int *) alloca(sizeof(int) * ordersq);
	e->nodes(nds);

	int faceindex;
	int *map = (int *) alloca(sizeof(int) * order);
	ipu.facemap(faceindex, nn, nds, map);

	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nds, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *v = (double *) alloca(sizeof(double) * ordersq * ordersq);
	SomGalFunction2d f(ordersq, v);
	ipu.zeroOut<double>(f.nrows() * f.ncolumns(), v);
	int gorder = 7;
	if (order <= 3) gorder = 4;
	ipu.lineInt2d(xyz, faceindex, f, gorder);
	ipu.symmetrize(f.ncolumns(), v);

	FullSquareMatrix sommerM(order, d);

	int i, j;
	for (i = 0; i < order; i++)
		for (j = 0; j < order; j++) {
			sommerM[i][j] = -v[map[j] * ordersq + map[i]] / e->getProperty()->rho;
		}

	return sommerM;
}


GenStackFSFullMatrix<double> IsoParamLineSommer::wetInterfaceMatrix(
		CoordSet &cs, double *d) {

	IsoParamUtils2d ipu(order);
	double *xyz = (double *) alloca(sizeof(double) * 3 * order);
	cs.getCoordinates(nn, order, xyz, xyz + order, xyz + 2 * order);

	WetInterfaceGalFunction2d f(order, d);
	ipu.zeroOut<double>(3 * order * order, d);
	int gorder = 7;
	if (order <= 3) gorder = 4;
	ipu.lineLineInt2d(xyz, f, gorder);

	GenStackFSFullMatrix<double> sommerM(order, 3 * order, d);
	return sommerM;
}


void IsoParamLineSommer::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {

	IsoParamUtils2d ipu(order);
	double *xyz = (double *) alloca(sizeof(double) * 3 * order);
	cs.getCoordinates(nn, order, xyz, xyz + order, xyz + 2 * order);

	double *d = (double *) alloca(sizeof(double) * 2 * order * order);
	WetInterfaceGalFunction2d f(order, d);
	ipu.zeroOut<double>(2 * order * order, d);
	int gorder = 7;
	if (order <= 3) gorder = 4;
	ipu.lineLineInt2d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < order; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in IsoParamLineSommer::wetInterfaceLMPC\n");
		return;
	}

 if (!sFlag)
	for (i = 0; i < order; i++) {
		LMPCTerm lmpct1(nn[i], 0, -d[i * order + j]);
		lmpc->addterm(&lmpct1);
		LMPCTerm lmpct2(nn[i], 1, -d[order * order + i * order + j]);
   lmpc->addterm(&lmpct2);
 }
 else
 for(i=0;i<order;i++) {
   LMPCTerm lmpct1(nn[i],0, 0.0 );
   lmpc->addterm(&lmpct1);
   LMPCTerm lmpct2(nn[i],1, 0.0 ); 
		lmpc->addterm(&lmpct2);
	}

}


FullSquareMatrix IsoParamLineSommer::turkelMatrix(CoordSet &cs, double *d) const {

	FullSquareMatrix sommerM(order, d);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;
	fprintf(stderr, "IsoParamLineSommer::turkelMatrix not implemented\n");
	return sommerM;
}


FullSquareMatrixC IsoParamLineSommer::sommer2Matrix(CoordSet &cs,
                                                    complex<double> *d) {
	FullSquareMatrixC sommerM(order, d);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;
	fprintf(stderr, "IsoParamLineSommer::sommer2Matrix not implemented\n");
	return sommerM;
}


void IsoParamLineSommer::getNormal(const CoordSet &cs, double normal[3]) const {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[order - 1]);

	double x[2], y[2];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;

	double nx, ny, l;
	nx = -y[1] + y[0];
	ny = x[1] - x[0];
	l = sqrt(nx * nx + ny * ny);
	nx = nx / l;
	ny = ny / l;

	normal[0] = nx;
	normal[1] = ny;
	normal[2] = 0.0;
}


void
IsoParamLineSommer::markDofs(DofSetArray &dsa) const {
	for (int i = 0; i < order; i++) dsa.mark(nn[i], DofSet::Helm);
}
