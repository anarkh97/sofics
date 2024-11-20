#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Utils.d/dbg_alloca.h>
#include <Element.d/Sommerfeld.d/LagLineSommer.h>
#include <Element.d/Helm.d/HelmElement.h>


LagLineSommer::LagLineSommer(int o, int *nodenums, Element *_el) {

	order = o;

	nn = new int[order];

	int i;
	for (i = 0; i < order; i++) nn[i] = nodenums[i];
	el = _el;
 sFlag = false;
 soundSpeed = 0.0;
}


LagLineSommer *LagLineSommer::clone() {
	LagLineSommer *se = new LagLineSommer(order, nn, el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
 se->soundSpeed = soundSpeed;
	return se;
}


int *LagLineSommer::dofs(DofSetArray &dsa, int *p) const  {

	if (p == 0) p = new int[order];
	int i;
	for (i = 0; i < order; i++) dsa.number(nn[i], DofSet::Helm, p + i);
	return p;
}


void LagLineSommer::neumVector(CoordSet &cs, ComplexVector &cv,
                               double kappa, double dx, double dy, double dz) {

	double *x = (double *) dbg_alloca(sizeof(double) * order);
	double *y = (double *) dbg_alloca(sizeof(double) * order);

	int i;
	for (i = 0; i < order; i++) {
		Node nd = cs.getNode(nn[i]);
		x[i] = nd.x;
		y[i] = nd.y;
	}

	for (i = 0; i < order; i++) cv[i] = ComplexD(0.0, 0.0);

	double normal[3];
	getNormal(cs, normal);

	int ng; // number of gauss points
	double *N; // shape functions and tangential derivative in gauss points
	double *gw; //gauss weights
	if (el == 0) {
		fprintf(stderr, "LagLineSommer::neumVector: adjacent element not defined.\n");
		exit(-1);
	}
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0) {
		fprintf(stderr, "A non-Helmholtz element found.\n");
		exit(-1);
	}
	he->edgeShapeFunctions(nn[0], nn[order - 1], &ng, &gw, &N);

	double rho = el->getProperty()->rho;

	for (i = 0; i < ng; i++) {
		double tau[2] = {0, 0};
		double xx[2] = {0, 0};
		int k;
		for (k = 0; k < order; k++) {
			tau[0] += N[i * order * 2 + k + order] * x[k];
			tau[1] += N[i * order * 2 + k + order] * y[k];
			xx[0] += N[i * order * 2 + k] * x[k];
			xx[1] += N[i * order * 2 + k] * y[k];
		}
//   double si = 1.0;
//   if (normal[0]*(-tau[1])+normal[1]*tau[0]<0.0) { 
//     si =-1.0;
//   }
		int j;
		for (j = 0; j < order; j++)
//     cv[j]+= ComplexD(0.0,si*gw[i]*kappa*(-tau[1]*dx+tau[0]*dy))*
			cv[j] += ComplexD(0.0, gw[i] * kappa *
			                       ((-tau[1] * dx + tau[0] * dy) - sqrt(tau[0] * tau[0] + tau[1] * tau[1]))) *
			         exp(ComplexD(0.0, kappa * (xx[0] * dx + xx[1] * dy))) * N[i * order * 2 + j] / rho;
	}
	delete[] gw;
	delete[] N;
}


FullSquareMatrix LagLineSommer::sommerMatrix(CoordSet &cs, double *d) const {

	double *x = (double *) dbg_alloca(sizeof(double) * order);
	double *y = (double *) dbg_alloca(sizeof(double) * order);

	int i;
	for (i = 0; i < order; i++) {
		Node nd = cs.getNode(nn[i]);
		x[i] = nd.x;
		y[i] = nd.y;
	}

	FullSquareMatrix sommerM(order, d);
	int j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;

	int ng; // number of gauss points
	double *N; // shape functions and tangential derivative in gauss points
	double *gw; //gauss weights
	if (el == 0) {
		fprintf(stderr, "LagLineSommer::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0) {
		fprintf(stderr, "A non-Helmholtz element found.\n");
		exit(-1);
	}
	he->edgeShapeFunctions(nn[0], nn[order - 1], &ng, &gw, &N);

	double rho = el->getProperty()->rho;

	for (i = 0; i < ng; i++) {
		double tau[2] = {0, 0};
		int k;
		for (k = 0; k < order; k++) {
			tau[0] += N[i * order * 2 + k + order] * x[k];
			tau[1] += N[i * order * 2 + k + order] * y[k];
		}
		int j;
		for (j = 0; j < order; j++)
			for (k = 0; k < order; k++)
				sommerM[j][k] -= gw[i] * N[i * order * 2 + j] * N[i * order * 2 + k] *
				                 sqrt(tau[0] * tau[0] + tau[1] * tau[1]) / rho;
	}
	delete[] gw;
	delete[] N;

	return sommerM;
}


FullSquareMatrix LagLineSommer::turkelMatrix(CoordSet &cs, double *d) const {

	double *x = (double *) dbg_alloca(sizeof(double) * order);
	double *y = (double *) dbg_alloca(sizeof(double) * order);

	int i;
	for (i = 0; i < order; i++) {
		Node nd = cs.getNode(nn[i]);
		x[i] = nd.x;
		y[i] = nd.y;
	}

	FullSquareMatrix sommerM(order, d);
	int j;
	for (i = 0; i < order; i++) for (j = 0; j < order; j++) sommerM[i][j] = 0.0;

	int ng; // number of gauss points
	double *N; // shape functions and tangential derivative in gauss points
	double *gw; //gauss weights
	if (el == 0) {
		fprintf(stderr, "LagLineSommer::turkelMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0) {
		fprintf(stderr, "A non-Helmholtz element found.\n");
		exit(-1);
	}
	he->edgeShapeFunctions(nn[0], nn[order - 1], &ng, &gw, &N);

	double rho = el->getProperty()->rho;

	for (i = 0; i < ng; i++) {
		double tau[2] = {0, 0};
		double xx[2] = {0, 0};
		int k;
		for (k = 0; k < order; k++) {
			tau[0] += N[i * order * 2 + k + order] * x[k];
			tau[1] += N[i * order * 2 + k + order] * y[k];
			xx[0] += N[i * order * 2 + k] * x[k];
			xx[1] += N[i * order * 2 + k] * y[k];
		}
		int j;
		for (j = 0; j < order; j++)
			for (k = 0; k < order; k++)
				sommerM[j][k] -= gw[i] * N[i * order * 2 + order + j] * N[i * order * 2 + order + k] /
				                 (sqrt(tau[0] * tau[0] + tau[1] * tau[1]) * rho);
	}
	delete[] gw;
	delete[] N;

	return sommerM;
}


void LagLineSommer::getNormal(const CoordSet &cs, double normal[3]) const {
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
