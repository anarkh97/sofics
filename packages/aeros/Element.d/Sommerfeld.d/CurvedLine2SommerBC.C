#include        <cstdio>
#include        <cstdlib>
#include        <cmath>
#include        <Element.d/Sommerfeld.d/CurvedLine2SommerBC.h>
#include        <Element.d/Helm.d/HelmElement.h>


CurvedLine2SommerBC::CurvedLine2SommerBC(int n1, int n2, int n3, Element *_el) {
	nn[0] = n1;
	nn[1] = n2;
	nn[2] = n3;
	el = _el;
 sFlag = false;
 soundSpeed = 0.0;

}


int *
CurvedLine2SommerBC::dofs(DofSetArray &dsa, int *p) const  {
	if (p == 0) p = new int[3];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p + 1);
	dsa.number(nn[2], DofSet::Helm, p + 2);

	return p;
}


CurvedLine2SommerBC *CurvedLine2SommerBC::clone() {
	CurvedLine2SommerBC *se = new CurvedLine2SommerBC(nn[0], nn[1], nn[2], el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
	return se;
}


void CurvedLine2SommerBC::neumVector(CoordSet &cs, ComplexVector &cv,
                                     double kappa, double dx, double dy, double dz, int pflag) {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;
	x[2] = nd3.x;
	y[2] = nd3.y;

	int i;
	for (i = 0; i < 3; i++) cv[i] = ComplexD(0.0, 0.0);

	double normal[3];
	getNormal(cs, normal);

	int ng; // number of gauss points
	double *N; // shape functions and tangential derivative in gauss points
	double *gw; //gauss weights
	if (el == 0) {
		fprintf(stderr, "CurvedLine2SommerBC::neumVector: adjacent element not defined.\n");
		exit(-1);
	}
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0) {
		fprintf(stderr, "A non-Helmholtz element found.\n");
		exit(-1);
	}
	he->edgeShapeFunctions(nn[0], nn[2], &ng, &gw, &N);
	double rho = el->getProperty()->rho;

	for (i = 0; i < ng; i++) {
		double tau[2] = {0, 0};
		double xx[2] = {0, 0};
		int k;
		for (k = 0; k < 3; k++) {
			tau[0] += N[i * 6 + k + 3] * x[k];
			tau[1] += N[i * 6 + k + 3] * y[k];
			xx[0] += N[i * 6 + k] * x[k];
			xx[1] += N[i * 6 + k] * y[k];
		}
		double si = 1.0;
		if (normal[0] * (-tau[1]) + normal[1] * tau[0] < 0.0) {
			si = -1.0;
		}
		int j;
		for (j = 0; j < 3; j++)
			cv[j] += ComplexD(0.0, si * gw[i] * kappa * (-tau[1] * dx + tau[0] * dy)) *
			         exp(ComplexD(0.0, kappa * (xx[0] * dx + xx[1] * dy))) * N[i * 6 + j] /
			         rho;
	}
	delete[] gw;
	delete[] N;
}


FullSquareMatrix
CurvedLine2SommerBC::sommerMatrix(CoordSet &cs, double *d) const {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;
	x[2] = nd3.x;
	y[2] = nd3.y;

	FullSquareMatrix sommerM(3, d);
	int i, j;
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) sommerM[i][j] = 0.0;

	int ng; // number of gauss points
	double *N; // shape functions and tangential derivative in gauss points
	double *gw; //gauss weights
	if (el == 0) {
		fprintf(stderr, "CurvedLine2SommerBC::neumVector: adjacent element not defined.\n");
		exit(-1);
	}
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0) {
		fprintf(stderr, "A non-Helmholtz element found.\n");
		exit(-1);
	}
	he->edgeShapeFunctions(nn[0], nn[2], &ng, &gw, &N);
	double rho = el->getProperty()->rho;

	for (i = 0; i < ng; i++) {
		double tau[2] = {0, 0};
		int k;
		for (k = 0; k < 3; k++) {
			tau[0] += N[i * 6 + k + 3] * x[k];
			tau[1] += N[i * 6 + k + 3] * y[k];
		}
		int j;
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				sommerM[j][k] -= gw[i] * N[i * 6 + j] * N[i * 6 + k] * sqrt(tau[0] * tau[0] + tau[1] * tau[1]) /
				                 rho;
	}
	delete[] gw;
	delete[] N;

	return sommerM;

}

FullSquareMatrix
CurvedLine2SommerBC::turkelMatrix(CoordSet &cs, double *d) const {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;
	x[2] = nd3.x;
	y[2] = nd3.y;

	FullSquareMatrix sommerM(3, d);
	int i, j;
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) sommerM[i][j] = 0.0;

	int ng; // number of gauss points
	double *N; // shape functions and tangential derivative in gauss points
	double *gw; //gauss weights
	if (el == 0) {
		fprintf(stderr, "CurvedLine2SommerBC::neumVector: adjacent element not defined.\n");
		exit(-1);
	}
	HelmElement *he = dynamic_cast<HelmElement *>(el);
	if (he == 0) {
		fprintf(stderr, "A non-Helmholtz element found.\n");
		exit(-1);
	}
	he->edgeShapeFunctions(nn[0], nn[2], &ng, &gw, &N);
	double rho = el->getProperty()->rho;

	for (i = 0; i < ng; i++) {
		double tau[2] = {0, 0};
		double xx[2] = {0, 0};
		int k;
		for (k = 0; k < 3; k++) {
			tau[0] += N[i * 6 + k + 3] * x[k];
			tau[1] += N[i * 6 + k + 3] * y[k];
			xx[0] += N[i * 6 + k] * x[k];
			xx[1] += N[i * 6 + k] * y[k];
		}
		int j;
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				sommerM[j][k] -= gw[i] * N[i * 6 + 3 + j] * N[i * 6 + 3 + k] /
				                 (sqrt(tau[0] * tau[0] + tau[1] * tau[1]) * rho);
	}
	delete[] gw;
	delete[] N;

	return sommerM;
}


void
CurvedLine2SommerBC::getNormal(const CoordSet &cs, double normal[3]) const {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[2]);

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
