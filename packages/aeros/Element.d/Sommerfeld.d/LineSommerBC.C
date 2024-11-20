#include        <cstdio>
#include        <cmath>
#include        <Element.d/Sommerfeld.d/LineSommerBC.h>
#include        <Element.d/Helm.d/IsoParamUtils2d.h>


LineSommerBC::LineSommerBC(int n1, int n2, Element *_el, int etype) {
	nn[0] = n1;
	nn[1] = n2;
	el = _el;
	setElementType(etype);
 sFlag = false;
 soundSpeed = 0.0;

}


LineSommerBC *LineSommerBC::clone() {
	LineSommerBC *se = new LineSommerBC(nn[0], nn[1], el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
 se->soundSpeed = soundSpeed;
	return se;
}


int *
LineSommerBC::dofs(DofSetArray &dsa, int *p) const  {
	if (p == 0) p = new int[2];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p + 1);

	return p;
}


FullSquareMatrix
LineSommerBC::sommerMatrix(CoordSet &cs, double *d) const {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);

	double x[2], y[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;

	// Get the length of the edge that sees the Sommerfeld b.c.
	double l = sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]));

	if (el == 0) {
		fprintf(stderr,
		        "LineSommerBC::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}

	double r = el->getProperty()->rho;

	double dd = (-1.0 / 3.0) * l / r;
	double ee = (-1.0 / 6.0) * l / r;

	FullSquareMatrix sommerM(2, d);

	sommerM[0][0] = dd;
	sommerM[0][1] = ee;
	sommerM[1][0] = ee;
	sommerM[1][1] = dd;

	return sommerM;
}


FullSquareMatrix
LineSommerBC::turkelMatrix(CoordSet &cs, double *d) const {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);

	double x[2], y[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;

	double l = sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]));


/*
  Note that the stiffness is given by:

              |  1.0  -1.0 |
    K = (1/h) |            | 
              | -1.0   1.0 |

  and the sign is flipped to keep the same 
  convention as the mass that also appears 
         in the integration by parts term.
*/

	FullSquareMatrix sommerK(2, d);

	if (el == 0) {
		fprintf(stderr,
		        "LineSommerBC::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;

	double dd = -1.0 / l / r;
	double ee = 1.0 / l / r;

	sommerK[0][0] = dd;
	sommerK[0][1] = ee;
	sommerK[1][0] = ee;
	sommerK[1][1] = dd;

	return sommerK;
}


void LineSommerBC::sommerVector(CoordSet &cs, ComplexVector &cv,
                                ComplexVector &cf) {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	double x[2], y[2];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;

	double nx, ny, l;
	nx = -y[1] + y[0];
	ny = x[1] - x[0];
	l = sqrt(nx * nx + ny * ny);

	cf[0] = l / 3.0 * cv[0] + l / 6.0 * cv[1];
	cf[1] = l / 6.0 * cv[0] + l / 3.0 * cv[1];
}

void LineSommerBC::btVector(CoordSet &cs, ComplexVector &cv,
                            ComplexVector &cf) {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	double x[2], y[2];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;

	double nx, ny, l;
	nx = -y[1] + y[0];
	ny = x[1] - x[0];
	l = sqrt(nx * nx + ny * ny);

	cf[0] = 1.0 / l * cv[0] - 1.0 / l * cv[1];
	cf[1] = -1.0 / l * cv[0] + 1.0 / l * cv[1];
}


void LineSommerBC::ffpDir(int ndir, ComplexD *ffp,
                          CoordSet &cs,
                          ComplexD *u, ComplexD *dudn,
                          double k, double (*dir)[3], double *idir) {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);

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

	ComplexD cc = exp(ComplexD(0.0, M_PI / 4.0)) / sqrt(8.0 * M_PI * k);

	for (int i = 0; i < ndir; i++) {
		double &dirx = dir[i][0], &diry = dir[i][1];
		ComplexD u1 = ComplexD(0.0, k * (dirx * nx + diry * ny)) * u[0] + dudn[0] / l;
		ComplexD u2 = ComplexD(0.0, k * (dirx * nx + diry * ny)) * u[1] + dudn[1] / l;

		ComplexD v1 = exp(-ComplexD(0.0, k * (dirx * x[0] + diry * y[0])));
		ComplexD v2 = exp(-ComplexD(0.0, k * (dirx * x[1] + diry * y[1])));

		ffp[i] += cc * l *
		          (u1 * (1.0 / 3.0 * v1 + 1.0 / 6.0 * v2) + u2 * (1.0 / 6.0 * v1 + 1.0 / 3.0 * v2));
	}
}


double
LineSommerBC::getSize(CoordSet &cs) {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);

	double x[2], y[2];
	x[0] = nd1.x;
	y[0] = nd1.y;
	x[1] = nd2.x;
	y[1] = nd2.y;
	double nx, ny;
	nx = -y[1] + y[0];
	ny = x[1] - x[0];
	return sqrt(nx * nx + ny * ny);
}


void
LineSommerBC::getNormal(const CoordSet &cs, double normal[3]) const {
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);

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
LineSommerBC::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {
	IsoParamUtils2d ipu(2);
	double *xyz = (double *) alloca(sizeof(double) * 6);
	cs.getCoordinates(nn, 2, xyz, xyz + 2, xyz + 4);

	double *d = (double *) alloca(sizeof(double) * 8);
	WetInterfaceGalFunction2d f(2, d);
	ipu.zeroOut<double>(8, d);
	int gorder = 4;
	ipu.lineLineInt2d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < 2; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in LineSommerBC::wetInterfaceLMPC\n");
		return;
	}

 if (!sFlag)
	for (i = 0; i < 2; i++) {
		LMPCTerm lmpct1(nn[i], 0, -d[i * 2 + j]);
		lmpc->addterm(&lmpct1);
		LMPCTerm lmpct2(nn[i], 1, -d[4 + i * 2 + j]);
   lmpc->addterm(&lmpct2);
 }
 else
 for(i=0;i<2;i++) {
   LMPCTerm lmpct1(nn[i],0, 0.0);
   lmpc->addterm(&lmpct1);
   LMPCTerm lmpct2(nn[i],1, 0.0);
		lmpc->addterm(&lmpct2);
	}
}
