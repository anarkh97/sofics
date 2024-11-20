#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <alloca.h>
#include <Element.d/Sommerfeld.d/SpectralIsoParamQuadSommer.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamHexa.h>
#include <Element.d/Helm.d/IsoParamUtils.h>

#define O3 4
#define O4 5

SpectralIsoParamQuadSommer::SpectralIsoParamQuadSommer(int o, int *nodenums, Element *_el, int etype) {

	order = int(rint(sqrt(double(o))));

	nn = new int[order * order];

	int i;
	for (i = 0; i < order * order; i++) nn[i] = nodenums[i];
	setElementType(etype);
	el = _el;
}


SpectralIsoParamQuadSommer *SpectralIsoParamQuadSommer::clone() {
	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	SpectralIsoParamQuadSommer *se = new SpectralIsoParamQuadSommer(ordersq, nn, el);
	se->el2 = el2;
	se->dom = dom;
	return se;
}


int *SpectralIsoParamQuadSommer::dofs(DofSetArray &dsa, int *p) const  {

	if (p == 0) p = new int[order * order];
	int i;
	for (i = 0; i < order * order; i++) dsa.number(nn[i], DofSet::Helm, p + i);
	return p;
}


int *SpectralIsoParamQuadSommer::solidDofs(DofSetArray &dsa, int *p) const {

	if (p == 0) p = new int[3 * order * order];
	int i;
	for (i = 0; i < order * order; i++) dsa.number(nn[i], DofSet::Xdisp | DofSet::Ydisp | DofSet::Zdisp, p + 3 * i);
	return p;
}


int *SpectralIsoParamQuadSommer::wetDofs(DofSetArray &dsa, int *p) const {

	if (p == 0) p = new int[order * order * 4];
	int i;
	for (i = 0; i < order * order; i++)
		dsa.number(nn[i],
		           DofSet::Helm | DofSet::Xdisp | DofSet::Ydisp | DofSet::Zdisp, p + 4 * i);

	return p;
}


void SpectralIsoParamQuadSommer::flipNormal() {
	int *nds = (int *) alloca(sizeof(int) * numNodes());
	nodes(nds);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) nn[i * order + j] = nds[j * order + i];
}

void SpectralIsoParamQuadSommer::neumVector(CoordSet &cs, ComplexVector &cv,
                                            double kappa, double dx, double dy, double dz) {

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	fprintf(stderr, "SpectralIsoParamQuadSommer::neumVector not implemented\n");
	int i;
	for (i = 0; i < ordersq; i++) {
		cv[i] = 0.0;
	}
}


void SpectralIsoParamQuadSommer::neumVectorDeriv(CoordSet &cs, ComplexVector &cv,
                                                 double kappa, double dx, double dy, double dz, int dero) {

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	fprintf(stderr, "SpectralIsoParamQuadSommer::neumVectorDeriv not implemented\n");
	int i;
	for (i = 0; i < ordersq; i++) {
		cv[i] = 0.0;
	}
}


FullSquareMatrix SpectralIsoParamQuadSommer::sommerMatrix(CoordSet &cs, double *d) const {

	if (el == 0) {
		fprintf(stderr, "SpectralIsoParamQuadSommer::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}

	IsoParamUtils ipu(order);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *v = (double *) alloca(sizeof(double) * ordersq * ordersq);
	SomGalFunction f(ordersq, v);
	ipu.zeroOut<double>(f.nrows() * f.ncolumns(), v);
	ipu.spectralSurfSurfInt3d(xyz, f);
	ipu.symmetrize(f.ncolumns(), v);

	FullSquareMatrix sommerM(ordersq, d);
	int i, j;
	for (i = 0; i < ordersq; i++)
		for (j = 0; j < ordersq; j++) {
			sommerM[i][j] = -v[j * ordersq + i] / el->getProperty()->rho;
		}
	return sommerM;
}


void SpectralIsoParamQuadSommer::ffp(CoordSet &cs, int numFFP, double *dirFFP,
                                     complex<double> *sol, complex<double> *ffpv) {
	fprintf(stderr, "SpectralIsoParamQuadSommer::ffp not implemented\n");
}


FullSquareMatrix SpectralIsoParamQuadSommer::turkelMatrix(CoordSet &cs, double *d) const {

	FullSquareMatrix sommerM(order * order, d);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;
	fprintf(stderr, "SpectralIsoParamQuadSommer::turkelMatrix not implemented\n");
	return sommerM;
}


FullSquareMatrixC SpectralIsoParamQuadSommer::sommer2Matrix(CoordSet &cs,
                                                            complex<double> *d) {

	FullSquareMatrixC sommerM(order * order, d);
	int i, j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;
	fprintf(stderr, "SpectralIsoParamQuadSommer::sommer2Matrix not implemented\n");
	return sommerM;
}


void SpectralIsoParamQuadSommer::getNormal(const CoordSet &cs, double normal[3]) const {

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
SpectralIsoParamQuadSommer::markDofs(DofSetArray &dsa) const {
	for (int i = 0; i < order * order; i++) dsa.mark(nn[i], DofSet::Helm);
}
