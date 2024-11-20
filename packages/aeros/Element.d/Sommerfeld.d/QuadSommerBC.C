#include        <cstdio>
#include        <cmath>
#include        <Element.d/Sommerfeld.d/QuadSommerBC.h>
#include        <Utils.d/linkfc.h>
#include        <Element.d/Helm.d/IsoParamUtils.h>


extern "C" {
void _FORTRAN(quad1dofm)(double *, double *, const int &, double *, const int &);
void    _FORTRAN(q4d1sommas)(double *, double *, const int &, double *, const int &);
};


QuadSommerBC::QuadSommerBC(int n1, int n2, int n3, int n4, Element *_el, int eType) : SommerElement(_el) {
	nn[0] = n1;
	nn[1] = n2;
	nn[2] = n3;
	nn[3] = n4;
	setElementType(eType);
	dom = 0;
 sFlag = false;
 soundSpeed = 0.0;

}


QuadSommerBC *QuadSommerBC::clone() {
	QuadSommerBC *se = new QuadSommerBC(nn[0], nn[1], nn[2], nn[3], el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
 se->soundSpeed = soundSpeed;
	return se;
}


FullSquareMatrix
QuadSommerBC::sommerMatrix(CoordSet &cs, double *d) const {
	double x[4], y[4], z[4];

	// Next two lines for the consistent mass
	// Note that the global coordinates are get
	// inside getLocalCoordinates(*) and after
	// that x, y and z will contain the local
	// coordinates

	getLocalCoordinates(cs, x, y, z);

	_FORTRAN(q4d1sommas)(x, y, 2, d, 4);

	if (el == 0) {
		fprintf(stderr,
		        "QuadSommerBC::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;

	for (int i = 0; i < 16; ++i)
		d[i] = -d[i] / r;

	FullSquareMatrix sommerM(4, d);

	return sommerM;
}


FullSquareMatrix
QuadSommerBC::turkelMatrix(CoordSet &cs, double *d) const {
	double x[4], y[4], z[4];

	// Next two lines for the consistent mass
	// Note that the global coordinates are get
	// inside getLocalCoordinates(*) and after
	// that x, y and z will contain the local
	// coordinates

	getLocalCoordinates(cs, x, y, z);

	const int numgauss = 2;
	const int numdof = 4;
	_FORTRAN(quad1dofm)(x, y, numgauss, d, numdof);

	if (el == 0) {
		fprintf(stderr,
		        "QuadSommerBC::turkelMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;

	for (int i = 0; i < 16; ++i)
		d[i] /= r;

	FullSquareMatrix sommerM(4, d);
	return sommerM;

}

void
QuadSommerBC::getLocalCoordinates(CoordSet &cs, double xx[4], double yy[4], double zz[4]) const {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);

	double x[4], y[4], z[4];
	double xi[3], eta[3], nu[3];

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

	// Make xi local axis connecting nodes 1 and 2
	double l12 = sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0])
	                  + (z[1] - z[0]) * (z[1] - z[0]));


	xi[0] = (x[1] - x[0]) / l12;
	xi[1] = (y[1] - y[0]) / l12;
	xi[2] = (z[1] - z[0]) / l12;

	// Make eta local axis connecting nodes 1 and 4
	double l14 = sqrt((x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0])
	                  + (z[3] - z[0]) * (z[3] - z[0]));

	eta[0] = (x[3] - x[0]) / l14;
	eta[1] = (y[3] - y[0]) / l14;
	eta[2] = (z[3] - z[0]) / l14;

	// Now orthogonalize eta with respect to xi
	double alpha = xi[0] * eta[0] + xi[1] * eta[1] + xi[2] * eta[2];

	eta[0] = eta[0] - alpha * xi[0];
	eta[1] = eta[1] - alpha * xi[1];
	eta[2] = eta[2] - alpha * xi[2];

	// Normalize eta
	double l = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
	eta[0] = eta[0] / l;
	eta[1] = eta[1] / l;
	eta[2] = eta[2] / l;

	// Get nu = xi (crossproduct) eta
	nu[0] = xi[1] * eta[2] - xi[2] * eta[1];
	nu[1] = xi[2] * eta[0] - xi[0] * eta[2];
	nu[2] = xi[0] * eta[1] - xi[1] * eta[0];


	//double dot1 = xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2];

	//double dot1 = xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2];
	//double dot2 = nu[0]*eta[0]+nu[1]*eta[1]+nu[2]*eta[2];
	//double dot3 = nu[0]*xi[0]+nu[1]*xi[1]+nu[2]*xi[2];
	//fprintf (stderr, "\n\n DOT1 = %e   DOT2 = %e  DOT2 = %e \n\n", dot1, dot2, dot3);

	// Now we perform a transformation in the type: xx = T x

	/*

	Assuming we have a local system defined with a
	set to unit vectors (i, j, k) and the global
	system as a set of unit vectors (i', j', k'),
	the transformation T is defined as:

             | i.i'  i.j'  i.k' |
	T =  | j.i'  j.j'  j.k' |
             | k.i'  k.j'  k.k' |


	Note that in our cartesian global system:
	   i' = (1, 0, 0)
	   j' = (0, 1, 0)
	   k' = (0, 0, 1)

	*/

	double T[3][3];
	T[0][0] = xi[0];
	T[0][1] = xi[1];
	T[0][2] = xi[2];
	T[1][0] = eta[0];
	T[1][1] = eta[1];
	T[1][2] = eta[2];
	T[2][0] = nu[0];
	T[2][1] = nu[1];
	T[2][2] = nu[2];

	x[1] -= x[0];
	y[1] -= y[0];
	z[1] -= z[0];
	x[2] -= x[0];
	y[2] -= y[0];
	z[2] -= z[0];
	x[3] -= x[0];
	y[3] -= y[0];
	z[3] -= z[0];
	x[0] = 0.0;
	y[0] = 0.0;
	z[0] = 0.0;


	xx[0] = T[0][0] * x[0] + T[0][1] * y[0] + T[0][2] * z[0];
	yy[0] = T[1][0] * x[0] + T[1][1] * y[0] + T[1][2] * z[0];
	zz[0] = T[2][0] * x[0] + T[2][1] * y[0] + T[2][2] * z[0];

	xx[1] = T[0][0] * x[1] + T[0][1] * y[1] + T[0][2] * z[1];
	yy[1] = T[1][0] * x[1] + T[1][1] * y[1] + T[1][2] * z[1];
	zz[1] = T[2][0] * x[1] + T[2][1] * y[1] + T[2][2] * z[1];

	xx[2] = T[0][0] * x[2] + T[0][1] * y[2] + T[0][2] * z[2];
	yy[2] = T[1][0] * x[2] + T[1][1] * y[2] + T[1][2] * z[2];
	zz[2] = T[2][0] * x[2] + T[2][1] * y[2] + T[2][2] * z[2];

	xx[3] = T[0][0] * x[3] + T[0][1] * y[3] + T[0][2] * z[3];
	yy[3] = T[1][0] * x[3] + T[1][1] * y[3] + T[1][2] * z[3];
	zz[3] = T[2][0] * x[3] + T[2][1] * y[3] + T[2][2] * z[3];

}


void
QuadSommerBC::getNormal(const CoordSet &cs, double normal[3]) const {
	double x[4], y[4], z[4];

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);

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

	// Define normal vector as vector product of diagonals

	double nx, ny, nz;

	nx = (y[2] - y[0]) * (z[3] - z[1]) - (z[2] - z[0]) * (y[3] - y[1]);
	ny = (z[2] - z[0]) * (x[3] - x[1]) - (x[2] - x[0]) * (z[3] - z[1]);
	nz = (x[2] - x[0]) * (y[3] - y[1]) - (y[2] - y[0]) * (x[3] - x[1]);

	double l = sqrt(nx * nx + ny * ny + nz * nz);

	normal[0] = nx / l;
	normal[1] = ny / l;
	normal[2] = nz / l;
}

void QuadSommerBC::GaussCoordinates(int Ngp, double *Pg, double *weight) const {
	if (Ngp == 2) {
		Pg[0] = -0.5773502691;
		Pg[1] = 0.5773502691;
		weight[0] = 1.0;
		weight[1] = 1.0;
	}
	if (Ngp == 3) {
		Pg[0] = -0.77459666924148337703585307995648;
		Pg[1] = 0.0;
		Pg[2] = 0.77459666924148337703585307995648;
		weight[0] = 0.5555555555555555555555;
		weight[1] = 0.8888888888888888888888;
		weight[2] = 0.5555555555555555555555;
	} else if (Ngp == 4) {
		Pg[0] = -0.86113631159405257522394648889281;
		Pg[1] = -0.33998104358485626480266575910324;
		Pg[2] = 0.33998104358485626480266575910324;
		Pg[3] = 0.86113631159405257522394648889281;
		weight[0] = 0.3478548451374538;
		weight[1] = 0.6521451548625461;
		weight[2] = 0.6521451548625461;
		weight[3] = 0.3478548451374538;
	} else if (Ngp == 5) {
		Pg[0] = -0.90617984593866399279762687829939;
		Pg[1] = -0.53846931010568309103631442070021;
		Pg[2] = 0.0;
		Pg[3] = 0.53846931010568309103631442070021;
		Pg[4] = 0.90617984593866399279762687829939;
		weight[0] = 0.2369268850561891;
		weight[1] = 0.4786286704993663;
		weight[2] = 0.5688888888888890;
		weight[3] = 0.4786286704993663;
		weight[4] = 0.2369268850561891;
	} else {
		fprintf(stderr, "QuadSommerBC.C, case not implemented");
		for (int i = 0; i < Ngp; i++) {
			Pg[i] = 0.0;
			weight[i] = 1.0;
		}
	}
}

//Functions that depend on the absorbing boundary condition order and require the domain information
#include<Driver.d/Domain.h>

int *
QuadSommerBC::dofs(DofSetArray &dsa, int *p) const  {
	double numD = -2.0;
	if (dom) numD = dom->solInfo().ATDARBFlag;
	if (numD == 1.5) {
		if (p == 0) p = new int[8];
		dsa.number(nn[0], DofSet::Helm, p);
		dsa.number(nn[0], DofSet::IntPress, p + 4);
		dsa.number(nn[1], DofSet::Helm, p + 1);
		dsa.number(nn[1], DofSet::IntPress, p + 5);
		dsa.number(nn[2], DofSet::Helm, p + 2);
		dsa.number(nn[2], DofSet::IntPress, p + 6);
		dsa.number(nn[3], DofSet::Helm, p + 3);
		dsa.number(nn[3], DofSet::IntPress, p + 7);
	} else {
		if (p == 0) p = new int[4];
		dsa.number(nn[0], DofSet::Helm, p);
		dsa.number(nn[1], DofSet::Helm, p + 1);
		dsa.number(nn[2], DofSet::Helm, p + 2);
		dsa.number(nn[3], DofSet::Helm, p + 3);
	}
	return p;
}

int QuadSommerBC::numDofs() const {
	double numD = -2.0;
	if (dom) numD = dom->solInfo().ATDARBFlag;
	if (numD == 1.5)
		return 8;
	else
		return 4;
}

void
QuadSommerBC::markDofs(DofSetArray &dsa) const {
	double numD = (dom) ? dom->solInfo().ATDARBFlag : 0.0;
	if (numD == 1.5) {
		dsa.mark(nn[0], DofSet::Helm);
		dsa.mark(nn[1], DofSet::Helm);
		dsa.mark(nn[2], DofSet::Helm);
		dsa.mark(nn[3], DofSet::Helm);
		dsa.mark(nn[0], DofSet::IntPress);
		dsa.mark(nn[1], DofSet::IntPress);
		dsa.mark(nn[2], DofSet::IntPress);
		dsa.mark(nn[3], DofSet::IntPress);
	} else {
		dsa.mark(nn[0], DofSet::Helm);
		dsa.mark(nn[1], DofSet::Helm);
		dsa.mark(nn[2], DofSet::Helm);
		dsa.mark(nn[3], DofSet::Helm);
	}
}

void QuadSommerBC::SurfaceRefinement(int nNo, double *x, double *y, double *z, double *xx, double *yy, double *zz) const {
	int ttnNo = 2 * nNo;
	int i;
	//fprintf(stderr,"  QuadSommer.C - refined interpolation for surfacic quad\n");
	// the point is for each node: 1 get the normal 2 create a point on the surface 3 make a 8 nodes interpolation
	// 1. access to the normals at all the nodes
	double *nx = new double[nNo], *ny = new double[nNo], *nz = new double[nNo];
	for (i = 0; i < nNo; i++) {
		int Snn = dom->nodeToSommerNodeMap[nn[i]];
		nx[i] = dom->curvatures_normal[Snn][0];
		ny[i] = dom->curvatures_normal[Snn][1];
		nz[i] = dom->curvatures_normal[Snn][2];
	}

	// 2. create artificial points on the surface
	double *xart = new double[nNo], *yart = new double[nNo], *zart = new double[nNo];
	double dd = 0.01;//small parameter
	xart[0] = xx[0] + xx[1] * (-1 + dd) + xx[2] * (-1 + dd) + xx[3] * (-1 + dd) * (-1 + dd);
	xart[1] = xx[0] + xx[1] * (1 - dd) + xx[2] * (-1 + dd) + xx[3] * (1 - dd) * (-1 + dd);
	xart[2] = xx[0] + xx[1] * (1 - dd) + xx[2] * (1 - dd) + xx[3] * (1 - dd) * (1 - dd);
	xart[3] = xx[0] + xx[1] * (-1 + dd) + xx[2] * (1 - dd) + xx[3] * (-1 + dd) * (1 - dd);
	yart[0] = yy[0] + yy[1] * (-1 + dd) + yy[2] * (-1 + dd) + yy[3] * (-1 + dd) * (-1 + dd);
	yart[1] = yy[0] + yy[1] * (1 - dd) + yy[2] * (-1 + dd) + yy[3] * (1 - dd) * (-1 + dd);
	yart[2] = yy[0] + yy[1] * (1 - dd) + yy[2] * (1 - dd) + yy[3] * (1 - dd) * (1 - dd);
	yart[3] = yy[0] + yy[1] * (-1 + dd) + yy[2] * (1 - dd) + yy[3] * (-1 + dd) * (1 - dd);
	zart[0] = zz[0] + zz[1] * (-1 + dd) + zz[2] * (-1 + dd) + zz[3] * (-1 + dd) * (-1 + dd);
	zart[1] = zz[0] + zz[1] * (1 - dd) + zz[2] * (-1 + dd) + zz[3] * (1 - dd) * (-1 + dd);
	zart[2] = zz[0] + zz[1] * (1 - dd) + zz[2] * (1 - dd) + zz[3] * (1 - dd) * (1 - dd);
	zart[3] = zz[0] + zz[1] * (-1 + dd) + zz[2] * (1 - dd) + zz[3] * (-1 + dd) * (1 - dd);
	//put the points in the tangenetial plan: Mart st MMart.norm=0 projection
	double alpha, beta;
	for (i = 0; i < nNo; i++) {
		alpha = (xart[i] - x[i]) * nx[i] + (yart[i] - y[i]) * ny[i] + (zart[i] - z[i]) * nz[i];
		beta = nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i];
		xart[i] -= alpha * nx[i] / beta;
		yart[i] -= alpha * ny[i] / beta;
		zart[i] -= alpha * nz[i] / beta;
	}

	// 8 nodes interpolation with x=x0+x1*u+x2*v+x3*u*v+(2-u*u-v*v)*(x4+x4*u+x6*v+x7*u*v), idem for y and z
	double *temp = new double[ttnNo * ttnNo];
	FullSquareMatrix iM(ttnNo, (double *) temp);
	double *u = new double[ttnNo], *v = new double[ttnNo];
	int *vectemp = new int[ttnNo];
	u[0] = -1.0;
	u[1] = 1.0;
	u[2] = 1.0;
	u[3] = -1.0;
	u[4] = -1.0 + dd;
	u[5] = 1.0 - dd;
	u[6] = 1.0 - dd;
	u[7] = -1.0 + dd;
	v[0] = -1.0;
	v[1] = -1.0;
	v[2] = 1.0;
	v[3] = 1.0;
	v[4] = -1.0 + dd;
	v[5] = -1.0 + dd;
	v[6] = 1.0 - dd;
	v[7] = 1.0 - dd;
	for (i = 0; i < nNo; i++) {
		xx[i] = x[i];
		xx[i + nNo] = xart[i];
		yy[i] = y[i];
		yy[i + nNo] = yart[i];
		zz[i] = z[i];
		zz[i + nNo] = zart[i];
	}
	for (i = 0; i < ttnNo; i++) {
		iM[i][0] = 1;
		iM[i][1] = u[i];
		iM[i][2] = v[i];
		iM[i][3] = u[i] * v[i];
		iM[i][4] = 2 - u[i] * u[i] - v[i] * v[i];
		iM[i][5] = (2 - u[i] * u[i] - v[i] * v[i]) * u[i];
		iM[i][6] = (2 - u[i] * u[i] - v[i] * v[i]) * v[i];
		iM[i][7] = (2 - u[i] * u[i] - v[i] * v[i]) * u[i] * v[i];
	}

	//solve iM*x = xx;
	int fl = ludcmp(iM, ttnNo, vectemp); // LU decomposition
	if (fl == 0) {
		lubksb(iM, ttnNo, vectemp, xx); // LU backsubstitution
		lubksb(iM, ttnNo, vectemp, yy);
		lubksb(iM, ttnNo, vectemp, zz);
	} else
		for (i = 0; i < ttnNo; i++) {
			xx[i] = 0.0;
			yy[i] = 0.0;
			zz[i] = 0.0;
		}
}

FullSquareMatrix
QuadSommerBC::refinedSommerMatrix(CoordSet &cs, double *d) {
	int i, j;
	int nNo = numNodes();
	int ttnNo = 2 * nNo;
	// reference: do Carmo, differential geometries of curves and surfaces p. 101
	// get the nodes coordinates
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	double *x = new double[nNo], *y = new double[nNo], *z = new double[nNo];
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

	double *xx = new double[ttnNo], *yy = new double[ttnNo], *zz = new double[ttnNo];
	// Express the surface depending on 2 parameters (u,v)
	// x = xx0 + xx1*u + xx2*v + xx3*u*v
	xx[0] = (x[0] + x[1] + x[2] + x[3]) / 4;
	xx[1] = (-x[0] + x[1] + x[2] - x[3]) / 4;
	xx[2] = (-x[0] - x[1] + x[2] + x[3]) / 4;
	xx[3] = (x[0] - x[1] + x[2] - x[3]) / 4;
	yy[0] = (y[0] + y[1] + y[2] + y[3]) / 4;
	yy[1] = (-y[0] + y[1] + y[2] - y[3]) / 4;
	yy[2] = (-y[0] - y[1] + y[2] + y[3]) / 4;
	yy[3] = (y[0] - y[1] + y[2] - y[3]) / 4;
	zz[0] = (z[0] + z[1] + z[2] + z[3]) / 4;
	zz[1] = (-z[0] + z[1] + z[2] - z[3]) / 4;
	zz[2] = (-z[0] - z[1] + z[2] + z[3]) / 4;
	zz[3] = (z[0] - z[1] + z[2] - z[3]) / 4;
	for (i = nNo; i < ttnNo; i++) {
		xx[i] = 0.0;
		yy[i] = 0.0;
		zz[i] = 0.0;
	}

	bool fast = false;//if fast==true, surfStiffMatrix ~= turkelMatrix
	if (!fast) {//refine the interpolation
		SurfaceRefinement(nNo, x, y, z, xx, yy, zz);
	}

	// determine the Gauss points in the local coordinates
	const int Ngp = 4;//Number of Gauss Points
	double Pg[Ngp];
	double weight[Ngp];
	GaussCoordinates(Ngp, Pg, weight);

	//define the nodes' caracteristic functions: no=(1+coefu[0]*u)*(1+coefv[0]*v)/4, ...
	double *coefu = new double[nNo], *coefv = new double[nNo];
	coefu[0] = -1.0;
	coefu[1] = 1.0;
	coefu[2] = 1.0;
	coefu[3] = -1.0;
	coefv[0] = -1.0;
	coefv[1] = -1.0;
	coefv[2] = 1.0;
	coefv[3] = 1.0;

	for (i = 0; i < nNo * nNo; i++) d[i] = 0.0;//initialize

	double val = 0.0;
	double jj, j0, j1, j2;
	double n1, n2;
	double *UU = new double[ttnNo];//Xu = [xx;yy;zz]*UU
	double *VV = new double[ttnNo];//Xv = [xx;yy;zz]*VV
	double Xux, Xuy, Xuz, Xvx, Xvy, Xvz;
	UU[0] = 0.0;
	UU[1] = 1.0;
	UU[2] = 0.0;
	VV[0] = 0.0;
	VV[1] = 0.0;
	VV[2] = 1.0;
	for (int igp = 0; igp < Ngp; igp++) {
		for (int jgp = 0; jgp < Ngp; jgp++) {
			//build paraVect = vector of the values of the parameter
			UU[3] = Pg[jgp];
			UU[4] = -2 * Pg[igp];
			UU[6] = -2 * Pg[igp] * Pg[jgp];
			UU[5] = -3 * Pg[igp] * Pg[igp] + 2 - Pg[jgp] * Pg[jgp];
			UU[7] = Pg[jgp] * (-3 * Pg[igp] * Pg[igp] + 2 - Pg[jgp] * Pg[jgp]);
			VV[3] = Pg[igp];
			VV[4] = -2 * Pg[jgp];
			VV[5] = -2 * Pg[igp] * Pg[jgp];
			VV[6] = -3 * Pg[jgp] * Pg[jgp] + 2 - Pg[igp] * Pg[igp];
			VV[7] = Pg[igp] * (-3 * Pg[jgp] * Pg[jgp] + 2 - Pg[igp] * Pg[igp]);
			Xux = 0.0;
			Xuy = 0.0;
			Xuz = 0.0;
			Xvx = 0.0;
			Xvy = 0.0;
			Xvz = 0.0;
			for (i = 0; i < ttnNo; i++) {
				Xux += xx[i] * UU[i];
				Xuy += yy[i] * UU[i];
				Xuz += zz[i] * UU[i];
				Xvx += xx[i] * VV[i];
				Xvy += yy[i] * VV[i];
				Xvz += zz[i] * VV[i];
			}
			j0 = Xuy * Xvz - Xuz * Xvy;
			j1 = Xuz * Xvx - Xux * Xvz;
			j2 = Xux * Xvy - Xuy * Xvx;
			jj = sqrt(j0 * j0 + j1 * j1 + j2 * j2);
			//fill the d matrix
			for (i = 0; i < nNo; i++) {
				for (j = 0; j < nNo; j++) {
					// values of the function to be integrated at the gauss point (same order)
					n1 = (1 + coefu[i] * Pg[igp]) * (1 + coefv[i] * Pg[jgp]) / 4;
					n2 = (1 + coefu[j] * Pg[igp]) * (1 + coefv[j] * Pg[jgp]) / 4;
					val = n1 * n2 * jj;
					// gauss integration: build d
					d[j + i * nNo] +=
							val * weight[igp] * weight[jgp];//order of d does not matter as the matrix is symetric
				}
			}
		}
	}
	FullSquareMatrix refiSommerM(nNo, d);

	delete[] x;
	x = NULL;
	delete[] y;
	y = NULL;
	delete[] z;
	z = NULL;
	delete[] xx;
	xx = NULL;
	delete[] yy;
	yy = NULL;
	delete[] zz;
	zz = NULL;
	delete[] coefu;
	coefu = NULL;
	delete[] coefv;
	coefv = NULL;
	delete[] UU;
	UU = NULL;
	delete[] VV;
	VV = NULL;

	return refiSommerM;
}

/*
// not used anymore - developped for ABC time domain - JFD
FullSquareMatrix
QuadSommerBC::surfStiffMatrix(CoordSet &cs, double *d)
{ 
  int i, j;
  int nNo = numNodes();
  int ttnNo = 2*nNo;
  // reference: do Carmo, differential geometries of curves and surfaces p. 101
  // get the nodes coordinates
  Node nd1 = cs.getNode(nn[0]);
  Node nd2 = cs.getNode(nn[1]);
  Node nd3 = cs.getNode(nn[2]);
  Node nd4 = cs.getNode(nn[3]);
  double *x=new double[nNo], *y=new double[nNo], *z=new double[nNo];
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

  double *xx=new double[ttnNo], *yy=new double[ttnNo], *zz=new double[ttnNo];
  // Express the surface depending on 2 parameters (u,v)
  // x = xx0 + xx1*u + xx2*v + xx3*u*v
  xx[0]=(x[0]+x[1]+x[2]+x[3])/4;
  xx[1]=(-x[0]+x[1]+x[2]-x[3])/4;
  xx[2]=(-x[0]-x[1]+x[2]+x[3])/4;
  xx[3]=(x[0]-x[1]+x[2]-x[3])/4;
  yy[0]=(y[0]+y[1]+y[2]+y[3])/4;
  yy[1]=(-y[0]+y[1]+y[2]-y[3])/4;
  yy[2]=(-y[0]-y[1]+y[2]+y[3])/4;
  yy[3]=(y[0]-y[1]+y[2]-y[3])/4;
  zz[0]=(z[0]+z[1]+z[2]+z[3])/4;
  zz[1]=(-z[0]+z[1]+z[2]-z[3])/4;
  zz[2]=(-z[0]-z[1]+z[2]+z[3])/4;
  zz[3]=(z[0]-z[1]+z[2]-z[3])/4;
  for (i = nNo ; i < ttnNo ; i++) {
    xx[i] = 0.0;
    yy[i] = 0.0;
    zz[i] = 0.0;
  }

  bool fast = false;//if fast==true, surfStiffMatrix ~= turkelMatrix
  if (!fast) {//refine the interpolation
    SurfaceRefinement(nNo, x, y, z, xx, yy, zz);
  }

  // determine the Gauss points in the local coordinates
  const int Ngp = 4;//Number of Gauss Points
  double Pg[Ngp];
  double weight[Ngp];
  GaussCoordinates(Ngp, Pg, weight);
  
  //define the nodes' caracteristic functions: no=(1+coefu[0]*u)*(1+coefv[0]*v)/4, ...
  double *coefu=new double[nNo], *coefv=new double[nNo];
  coefu[0] = -1.0 ; coefu[1] = 1.0 ; coefu[2] = 1.0 ; coefu[3] = -1.0 ;
  coefv[0] = -1.0 ; coefv[1] = -1.0 ; coefv[2] = 1.0 ; coefv[3] = 1.0 ;

  for (i = 0 ; i < nNo*nNo ; i++) d[i] = 0.0;//initialize

  double val = 0.0;
  double jj, j0, j1, j2, E, F, G;
  double dn1u, dn1v, dn2u, dn2v;
  double *UU=new double[ttnNo];//Xu = [xx;yy;zz]*UU
  double *VV=new double[ttnNo];//Xv = [xx;yy;zz]*VV
  double Xux, Xuy, Xuz, Xvx, Xvy, Xvz;
  UU[0] = 0.0; UU[1] = 1.0; UU[2] = 0.0;
  VV[0] = 0.0; VV[1] = 0.0; VV[2] = 1.0;
  for (int igp = 0 ; igp < Ngp ; igp++) {
    for (int jgp = 0 ; jgp < Ngp ; jgp++) {
      //build paraVect = vector of the values of the parameter
      UU[3] = Pg[jgp]; UU[4] = -2*Pg[igp]; UU[6] = -2*Pg[igp]*Pg[jgp];
        UU[5] = -3*Pg[igp]*Pg[igp]+2-Pg[jgp]*Pg[jgp]; UU[7] = Pg[jgp]*(-3*Pg[igp]*Pg[igp]+2-Pg[jgp]*Pg[jgp]);
      VV[3] = Pg[igp]; VV[4] = -2*Pg[jgp]; VV[5] = -2*Pg[igp]*Pg[jgp];
        VV[6] = -3*Pg[jgp]*Pg[jgp]+2-Pg[igp]*Pg[igp]; VV[7] = Pg[igp]*(-3*Pg[jgp]*Pg[jgp]+2-Pg[igp]*Pg[igp]);
      Xux = 0.0; Xuy = 0.0; Xuz = 0.0;
      Xvx = 0.0; Xvy = 0.0; Xvz = 0.0;
      for (i = 0 ; i < ttnNo ; i++) {
        Xux += xx[i]*UU[i];
        Xuy += yy[i]*UU[i];
        Xuz += zz[i]*UU[i];
        Xvx += xx[i]*VV[i];
        Xvy += yy[i]*VV[i];
        Xvz += zz[i]*VV[i];
      }
      E = Xux*Xux + Xuy*Xuy + Xuz*Xuz;
      F = Xvx*Xux + Xvy*Xuy + Xvz*Xuz;
      G = Xvx*Xvx + Xvy*Xvy + Xvz*Xvz;
      j0 = Xuy*Xvz-Xuz*Xvy;
      j1 = Xuz*Xvx-Xux*Xvz;
      j2 = Xux*Xvy-Xuy*Xvx;
      jj = sqrt(j0*j0+j1*j1+j2*j2);
      //fill the d matrix
      for (i = 0 ; i < nNo ; i++) {
        for (j = 0 ; j < nNo ; j++) {
          // values of the function to be integrated at the gauss point (same order)
          dn1u = coefu[i]*(1+coefv[i]*Pg[jgp])/4 ; dn1v = (1+coefu[i]*Pg[igp])*coefv[i]/4 ;
          dn2u = coefu[j]*(1+coefv[j]*Pg[jgp])/4 ; dn2v = (1+coefu[j]*Pg[igp])*coefv[j]/4 ;
          val = (dn1u*G*dn2u-dn1u*F*dn2v-dn1v*F*dn2u+dn1v*E*dn2v)/(E*G-F*F)*jj;
          // gauss integration: build d
          d[j+i*nNo] += val*weight[igp]*weight[jgp];//order of d does not matter as the matrix is symetric
        } 
      }
    }
  }
  FullSquareMatrix surfStiffM(nNo,d);
  
  delete [] x; x = NULL;
  delete [] y; y = NULL;
  delete [] z; z = NULL;
  delete [] xx; xx = NULL;
  delete [] yy; yy = NULL;
  delete [] zz; zz = NULL;
  delete [] coefu; coefu = NULL;
  delete [] coefv; coefv = NULL;
  delete [] UU; UU = NULL;
  delete [] VV; VV = NULL;

  return surfStiffM;
}
*/

FullSquareMatrix
QuadSommerBC::HSommerMatrix(const CoordSet &cs, double *d) const {
	int i, j, k;
	int nNo = numNodes();
	int ttnNo = 2 * nNo;
	// get the nodes coordinates
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	double *x = new double[nNo], *y = new double[nNo], *z = new double[nNo];
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
	// Express the surface depending on 2 parameters (u,v)
	// x = xx0 + xx1*u + xx2*v + xx3*u*v, etc
	double *xx = new double[ttnNo], *yy = new double[ttnNo], *zz = new double[ttnNo];
	xx[0] = (x[0] + x[1] + x[2] + x[3]) / 4;
	xx[1] = (-x[0] + x[1] + x[2] - x[3]) / 4;
	xx[2] = (-x[0] - x[1] + x[2] + x[3]) / 4;
	xx[3] = (x[0] - x[1] + x[2] - x[3]) / 4;
	yy[0] = (y[0] + y[1] + y[2] + y[3]) / 4;
	yy[1] = (-y[0] + y[1] + y[2] - y[3]) / 4;
	yy[2] = (-y[0] - y[1] + y[2] + y[3]) / 4;
	yy[3] = (y[0] - y[1] + y[2] - y[3]) / 4;
	zz[0] = (z[0] + z[1] + z[2] + z[3]) / 4;
	zz[1] = (-z[0] + z[1] + z[2] - z[3]) / 4;
	zz[2] = (-z[0] - z[1] + z[2] + z[3]) / 4;
	zz[3] = (z[0] - z[1] + z[2] - z[3]) / 4;
	for (i = nNo; i < ttnNo; i++) {
		xx[i] = 0.0;
		yy[i] = 0.0;
		zz[i] = 0.0;
	}
	bool fast = false;//if fast==true, surfStiffMatrix ~= turkelMatrix
	if (!fast) {//refine the interpolation
		SurfaceRefinement(nNo, x, y, z, xx, yy, zz);
	}

	// determine the Gauss points in the local coordinates
	const int Ngp = 4;//or 4 or ...
	double Pg[Ngp];
	double weight[Ngp];
	GaussCoordinates(Ngp, Pg, weight);

	//define the nodes' caracteristic functions: no=(1+coefu[0]*u)*(1+coefv[0]*v)/4, ...
	double *coefu = new double[nNo], *coefv = new double[nNo];
	coefu[0] = -1.0;
	coefu[1] = 1.0;
	coefu[2] = 1.0;
	coefu[3] = -1.0;
	coefv[0] = -1.0;
	coefv[1] = -1.0;
	coefv[2] = 1.0;
	coefv[3] = 1.0;

	//detemine the curvatures at each node
	double *HH = new double[nNo];
	for (i = 0; i < nNo; i++) {
		int Snn = dom->nodeToSommerNodeMap[nn[i]];
		HH[i] = dom->curvaturesH[Snn];
	}

	for (i = 0; i < nNo * nNo; i++) d[i] = 0.0;//initialize

	double val = 0.0;
	double jj, j0, j1, j2;
	double n1, n2, n3;
	double *UU = new double[ttnNo];//Xu = [xx;yy;zz]*UU
	double *VV = new double[ttnNo];//Xv = [xx;yy;zz]*VV
	double Xux, Xuy, Xuz, Xvx, Xvy, Xvz;
	UU[0] = 0.0;
	UU[1] = 1.0;
	UU[2] = 0.0;
	VV[0] = 0.0;
	VV[1] = 0.0;
	VV[2] = 1.0;
	for (int igp = 0; igp < Ngp; igp++) {
		for (int jgp = 0; jgp < Ngp; jgp++) {
			UU[3] = Pg[jgp];
			UU[4] = -2 * Pg[igp];
			UU[6] = -2 * Pg[igp] * Pg[jgp];
			UU[5] = -3 * Pg[igp] * Pg[igp] + 2 - Pg[jgp] * Pg[jgp];
			UU[7] = Pg[jgp] * (-3 * Pg[igp] * Pg[igp] + 2 - Pg[jgp] * Pg[jgp]);
			VV[3] = Pg[igp];
			VV[4] = -2 * Pg[jgp];
			VV[5] = -2 * Pg[igp] * Pg[jgp];
			VV[6] = -3 * Pg[jgp] * Pg[jgp] + 2 - Pg[igp] * Pg[igp];
			VV[7] = Pg[igp] * (-3 * Pg[jgp] * Pg[jgp] + 2 - Pg[igp] * Pg[igp]);
			Xux = 0.0;
			Xuy = 0.0;
			Xuz = 0.0;
			Xvx = 0.0;
			Xvy = 0.0;
			Xvz = 0.0;
			for (i = 0; i < ttnNo; i++) {
				Xux += xx[i] * UU[i];
				Xuy += yy[i] * UU[i];
				Xuz += zz[i] * UU[i];
				Xvx += xx[i] * VV[i];
				Xvy += yy[i] * VV[i];
				Xvz += zz[i] * VV[i];
			}
			j0 = Xuy * Xvz - Xuz * Xvy;
			j1 = Xuz * Xvx - Xux * Xvz;
			j2 = Xux * Xvy - Xuy * Xvx;
			jj = sqrt(j0 * j0 + j1 * j1 + j2 * j2);
			//fill the d matrix
			for (i = 0; i < nNo; i++) {
				for (j = 0; j < nNo; j++) {
					for (k = 0; k < nNo; k++) {
						// values of the function to be integrated at the gauss point (same order)
						n1 = (1 + coefu[i] * Pg[igp]) * (1 + coefv[i] * Pg[jgp]) / 4;
						n2 = (1 + coefu[j] * Pg[igp]) * (1 + coefv[j] * Pg[jgp]) / 4;
						n3 = (1 + coefu[k] * Pg[igp]) * (1 + coefv[k] * Pg[jgp]) / 4;
						val = HH[k] * n1 * n2 * n3 * jj;
						// gauss integration: build d
						d[j + i * nNo] +=
								val * weight[igp] * weight[jgp];//order of d does not matter as the matrix is symetric
					}
				}
			}
		}
	}
	FullSquareMatrix HSommerM(nNo, d);

	delete[] x;
	x = NULL;
	delete[] y;
	y = NULL;
	delete[] z;
	z = NULL;
	delete[] xx;
	xx = NULL;
	delete[] yy;
	yy = NULL;
	delete[] zz;
	zz = NULL;
	delete[] coefu;
	coefu = NULL;
	delete[] coefv;
	coefv = NULL;
	delete[] HH;
	HH = NULL;
	delete[] UU;
	UU = NULL;
	delete[] VV;
	VV = NULL;

	return HSommerM;
}

/*
// this function is no longer used - has been developped for ABC time domain - JFD
FullSquareMatrix
QuadSommerBC::HKSommerMatrix(CoordSet &cs, double *d)
{ 
  int i, j, k;
  int nNo = numNodes();
  int ttnNo = 2*nNo;
  // get the nodes coordinates
  Node nd1 = cs.getNode(nn[0]);
  Node nd2 = cs.getNode(nn[1]);
  Node nd3 = cs.getNode(nn[2]);
  Node nd4 = cs.getNode(nn[3]);
  double *x=new double[nNo], *y=new double[nNo], *z=new double[nNo];
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  // Express the surface depending on 2 parameters (u,v)
  // x = xx0 + xx1*u + xx2*v + xx3*u*v, etc
  double *xx=new double[ttnNo], *yy=new double[ttnNo], *zz=new double[ttnNo];
  xx[0]=(x[0]+x[1]+x[2]+x[3])/4;
  xx[1]=(-x[0]+x[1]+x[2]-x[3])/4;
  xx[2]=(-x[0]-x[1]+x[2]+x[3])/4;
  xx[3]=(x[0]-x[1]+x[2]-x[3])/4;
  yy[0]=(y[0]+y[1]+y[2]+y[3])/4;
  yy[1]=(-y[0]+y[1]+y[2]-y[3])/4;
  yy[2]=(-y[0]-y[1]+y[2]+y[3])/4;
  yy[3]=(y[0]-y[1]+y[2]-y[3])/4;
  zz[0]=(z[0]+z[1]+z[2]+z[3])/4;
  zz[1]=(-z[0]+z[1]+z[2]-z[3])/4;
  zz[2]=(-z[0]-z[1]+z[2]+z[3])/4;
  zz[3]=(z[0]-z[1]+z[2]-z[3])/4;
  for (i = nNo ; i < ttnNo ; i++) {
    xx[i] = 0.0;
    yy[i] = 0.0;
    zz[i] = 0.0;
  }
  bool fast = false;//if fast==true, surfStiffMatrix ~= turkelMatrix
  if (!fast) {//refine the interpolation
    SurfaceRefinement(nNo, x, y, z, xx, yy, zz);
  }
  // determine the Gauss points in the local coordinates
  const int Ngp = 4;//or 4 or ...
  double Pg[Ngp];
  double weight[Ngp];
  GaussCoordinates(Ngp, Pg, weight);
  
  //define the nodes' caracteristic functions: no=(1+coefu[0]*u)*(1+coefv[0]*v)/4, ...
  double *coefu=new double[nNo], *coefv=new double[nNo];
  coefu[0] = -1.0 ; coefu[1] = 1.0 ; coefu[2] = 1.0 ; coefu[3] = -1.0 ;
  coefv[0] = -1.0 ; coefv[1] = -1.0 ; coefv[2] = 1.0 ; coefv[3] = 1.0 ;

  //detemine the curvatures at each node
  double *HH=new double[nNo], *KK=new double[nNo];
  for (i = 0 ; i < nNo ; i++) {
    int Snn = dom->nodeToSommerNodeMap[nn[i]];
    HH[i]=dom->curvaturesH[Snn];
    KK[i]=dom->curvaturesK[Snn];
  }

  for (i = 0 ; i < nNo*nNo ; i++) d[i] = 0.0;//initialize

  double val = 0.0;
  double jj, j0, j1, j2;
  double n1, n2, n3;
  double *UU=new double[ttnNo];//Xu = [xx;yy;zz]*UU
  double *VV=new double[ttnNo];//Xv = [xx;yy;zz]*VV
  double Xux, Xuy, Xuz, Xvx, Xvy, Xvz;
  UU[0] = 0.0; UU[1] = 1.0; UU[2] = 0.0;
  VV[0] = 0.0; VV[1] = 0.0; VV[2] = 1.0;
  for (int igp = 0 ; igp < Ngp ; igp++) {
    for (int jgp = 0 ; jgp < Ngp ; jgp++) {
      UU[3] = Pg[jgp]; UU[4] = -2*Pg[igp]; UU[6] = -2*Pg[igp]*Pg[jgp];
        UU[5] = -3*Pg[igp]*Pg[igp]+2-Pg[jgp]*Pg[jgp]; UU[7] = Pg[jgp]*(-3*Pg[igp]*Pg[igp]+2-Pg[jgp]*Pg[jgp]);
      VV[3] = Pg[igp]; VV[4] = -2*Pg[jgp]; VV[5] = -2*Pg[igp]*Pg[jgp];
        VV[6] = -3*Pg[jgp]*Pg[jgp]+2-Pg[igp]*Pg[igp]; VV[7] = Pg[igp]*(-3*Pg[jgp]*Pg[jgp]+2-Pg[igp]*Pg[igp]);
      Xux = 0.0; Xuy = 0.0; Xuz = 0.0; Xvx = 0.0; Xvy = 0.0; Xvz = 0.0;
      for (i = 0 ; i < ttnNo ; i++) {
        Xux += xx[i]*UU[i];
        Xuy += yy[i]*UU[i];
        Xuz += zz[i]*UU[i];
        Xvx += xx[i]*VV[i];
        Xvy += yy[i]*VV[i];
        Xvz += zz[i]*VV[i];
      }
      j0 = Xuy*Xvz-Xuz*Xvy;
      j1 = Xuz*Xvx-Xux*Xvz;
      j2 = Xux*Xvy-Xuy*Xvx;
      jj = sqrt(j0*j0+j1*j1+j2*j2);
      //fill the d matrix
      for (i = 0 ; i < nNo ; i++) {
        for (j = 0 ; j < nNo ; j++) {
          for (k = 0 ; k < nNo ; k++) {
            // values of the function to be integrated at the gauss point (same order)
            n1 = (1+coefu[i]*Pg[igp])*(1+coefv[i]*Pg[jgp])/4 ;
            n2 = (1+coefu[j]*Pg[igp])*(1+coefv[j]*Pg[jgp])/4 ;
            n3 = (1+coefu[k]*Pg[igp])*(1+coefv[k]*Pg[jgp])/4 ;
            val = (HH[k]*HH[k]-KK[k])*n1*n2*n3*jj;
            // gauss integration: build d
            d[j+i*nNo] += val*weight[igp]*weight[jgp];
          }
        } 
      }
    }
  }
  FullSquareMatrix HKSommerM(nNo,d);
  
  delete [] x; x = NULL;
  delete [] y; y = NULL;
  delete [] z; z = NULL;
  delete [] xx; xx = NULL;
  delete [] yy; yy = NULL;
  delete [] zz; zz = NULL;
  delete [] coefu; coefu = NULL;
  delete [] coefv; coefv = NULL;
  delete [] HH; HH = NULL;
  delete [] KK; KK = NULL;
  delete [] UU; UU = NULL;
  delete [] VV; VV = NULL;

  return HKSommerM;
}
*/

void QuadSommerBC::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {
	int nn[4] = {QuadSommerBC::nn[0], QuadSommerBC::nn[1], QuadSommerBC::nn[3], QuadSommerBC::nn[2]};
	IsoParamUtils ipu(2);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *d = (double *) alloca(sizeof(double) * 3 * ordersq * ordersq);
	WetInterfaceGalFunction f(ordersq, d);
	ipu.zeroOut<double>(3 * ordersq * ordersq, d);
	int gorder = 4;
	ipu.surfSurfInt3d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < ordersq; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in QuadSommerBC::wetInterfaceLMPC\n");
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

