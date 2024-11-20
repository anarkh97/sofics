#include        <cstdio>
#include        <cmath>
#include        <Element.d/Sommerfeld.d/Line2SommerBC.h>
#include        <Element.d/Helm.d/IsoParamUtils2d.h>


Line2SommerBC::Line2SommerBC(int n1, int n2, int n3, Element *_el, int etype) 
{
	nn[0] = n1;
	nn[1] = n2;
	nn[2] = n3;
	el = _el;
	setElementType(etype);
 sFlag = false;
 soundSpeed = 0.0;

}


Line2SommerBC *Line2SommerBC::clone() 
{
	Line2SommerBC *se = new Line2SommerBC(nn[0], nn[1], nn[2], el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
 se->soundSpeed = soundSpeed;
	return se;
}


int *
Line2SommerBC::dofs(DofSetArray &dsa, int *p) const  
{
	if (p == 0) p = new int[3];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p + 1);
	dsa.number(nn[2], DofSet::Helm, p + 2);

	return p;
}


FullSquareMatrix
Line2SommerBC::sommerMatrix(CoordSet &cs, double *d) const 
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;
	x[2] = nd3.x;
	y[2] = nd3.y;
	z[2] = nd3.z;

	// Get the length of the edge that sees the Sommerfeld b.c.
	double l = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]) +
	                (z[2] - z[0]) * (z[2] - z[0]));

	if (el == 0) {
		fprintf(stderr, "Line2SommerBC::sommerMatrix: adjacent element not defined.\n");
	}

	double r = (el == 0) ? 1.0 : el->getProperty()->rho;

	double aa = (4.0 / 15.0) * (l / 2.0) / r;
	double bb = (2.0 / 15.0) * (l / 2.0) / r;
	double cc = (-1.0 / 15.0) * (l / 2.0) / r;
	double dd = (16.0 / 15.0) * (l / 2.0) / r;
	double ee = (2.0 / 15.0) * (l / 2.0) / r;
	double ff = (4.0 / 15.0) * (l / 2.0) / r;

	FullSquareMatrix sommerM(3, d);

	// Consistent Matrix
	sommerM[0][0] = -aa;
	sommerM[0][1] = -bb;
	sommerM[0][2] = -cc;
	sommerM[1][0] = -bb;
	sommerM[1][1] = -dd;
	sommerM[1][2] = -ee;
	sommerM[2][0] = -cc;
	sommerM[2][1] = -ee;
	sommerM[2][2] = -ff;

	return sommerM;
}

FullSquareMatrix
Line2SommerBC::turkelMatrix(CoordSet &cs, double *d) const {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;
	x[2] = nd3.x;
	y[2] = nd3.y;
	z[2] = nd3.z;

	// Get the length of the edge that sees the Sommerfeld b.c.
	double l = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]) +
	                (z[2] - z[0]) * (z[2] - z[0]));


	// This is for the Bayliss-Turkel Radiation condition
	FullSquareMatrix sommerK(3, d);

	if (el == 0) {
		fprintf(stderr, "Line2SommerBC::turkelMatrix: adjacent element not defined.\n");
	}

	double r = (el == 0) ? 1.0 : el->getProperty()->rho;

	double aa = (7.0 / 6.0) * (2.0 / l) / r;
	double bb = (-4.0 / 3.0) * (2.0 / l) / r;
	double cc = (1.0 / 6.0) * (2.0 / l) / r;
	double dd = (8.0 / 3.0) * (2.0 / l) / r;
	double ee = (-4.0 / 3.0) * (2.0 / l) / r;
	double ff = (7.0 / 6.0) * (2.0 / l) / r;


	sommerK[0][0] = -aa;
	sommerK[0][1] = -bb;
	sommerK[0][2] = -cc;
	sommerK[1][0] = -bb;
	sommerK[1][1] = -dd;
	sommerK[1][2] = -ee;
	sommerK[2][0] = -cc;
	sommerK[2][1] = -ee;
	sommerK[2][2] = -ff;

	return sommerK;
}


void
Line2SommerBC::getNormal(const CoordSet &cs, double normal[3]) const {
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


void Line2SommerBC::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {
	int nn[3] = {Line2SommerBC::nn[0], Line2SommerBC::nn[2], Line2SommerBC::nn[1]};
	IsoParamUtils2d ipu(3);
	double *xyz = (double *) alloca(sizeof(double) * 9);
	cs.getCoordinates(nn, 3, xyz, xyz + 3, xyz + 6);

	double *d = (double *) alloca(sizeof(double) * 18);
	WetInterfaceGalFunction2d f(3, d);
	ipu.zeroOut<double>(18, d);
	int gorder = 4;
	ipu.lineLineInt2d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < 3; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in Line2SommerBC::wetInterfaceLMPC\n");
		return;
	}

 if (!sFlag)
	for (i = 0; i < 3; i++) {
		LMPCTerm lmpct1(nn[i], 0, -d[i * 3 + j]);
		lmpc->addterm(&lmpct1);
		LMPCTerm lmpct2(nn[i], 1, -d[9 + i * 3 + j]);
   lmpc->addterm(&lmpct2);
 }
 else
 for(i=0;i<3;i++) {
   LMPCTerm lmpct1(nn[i],0, 0.0);
   lmpc->addterm(&lmpct1);
   LMPCTerm lmpct2(nn[i],1, 0.0);
		lmpc->addterm(&lmpct2);
	}
}

