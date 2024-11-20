#include        <cstdio>
#include        <cmath>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Helm.d/IsoParamUtils.h>

TriangleSommerBC::TriangleSommerBC(int n1, int n2, int n3, Element *_el, int eType) : SommerElement(_el) {
	nn[0] = n1;
	nn[1] = n2;
	nn[2] = n3;
	setElementType(eType);
	dom = 0;
 sFlag = false;
 soundSpeed = 0.0;

}


TriangleSommerBC *TriangleSommerBC::clone() {
	TriangleSommerBC *se = new TriangleSommerBC(nn[0], nn[1], nn[2], el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
 se->soundSpeed = soundSpeed;
	return se;
}


FullSquareMatrix
TriangleSommerBC::sommerMatrix(CoordSet &cs, double *d) const 
{
	// This is for the lumped mass
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

	// Define vectors:
	//                u: vector connecting nd1 and nd2
	//                v: vector connecting nd1 and nd3
	//                w: cross product between u and v

	double u1, u2, u3, v1, v2, v3, w1, w2, w3;

	u1 = x[1] - x[0];
	u2 = y[1] - y[0];
	u3 = z[1] - z[0];

	v1 = x[2] - x[0];
	v2 = y[2] - y[0];
	v3 = z[2] - z[0];

	w1 = u2 * v3 - u3 * v2;
	w2 = u3 * v1 - u1 * v3;
	w3 = u1 * v2 - u2 * v1;

	double area = 0.5 * sqrt(w1 * w1 + w2 * w2 + w3 * w3);

	FullSquareMatrix sommerM(3, d);

	if (el == 0) {
		fprintf(stderr,
		        "TriangleSommerBC::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;

	// Consistent mass
	double dd = -area / 6.0 / r;
	double ee = -area / 12.0 / r;

	// Lumped mass
	//double dd = area/3.0;
	//double ee = 0.0;

	sommerM[0][0] = dd;
	sommerM[0][1] = ee;
	sommerM[0][2] = ee;
	sommerM[1][0] = ee;
	sommerM[1][1] = dd;
	sommerM[1][2] = ee;
	sommerM[2][0] = ee;
	sommerM[2][1] = ee;
	sommerM[2][2] = dd;


	return sommerM;
}


FullSquareMatrix
TriangleSommerBC::turkelMatrix(CoordSet &cs, double *d) const {
	double x[3], y[3], z[3];

	// Transformation for the 3D "space" plane
	getLocalCoordinates(cs, x, y, z);

	double area = 0.5 * ((x[1] * y[2] - x[2] * y[1]) +
	                     (x[2] * y[0] - x[0] * y[2]) +
	                     (x[0] * y[1] - x[1] * y[0]));

	double x21 = x[1] - x[0];
	double x32 = x[2] - x[1];
	double x13 = x[0] - x[2];

	double y12 = y[0] - y[1];
	double y23 = y[1] - y[2];
	double y31 = y[2] - y[0];


	FullSquareMatrix K(3, d);

	K.zero();

	double ke = 0.25 / area;
	//double ke = 0.50/(x12*y13-y12*x13);

	// Mine
	K[0][0] = ke * (x32 * x32 + y23 * y23);
	K[0][1] = ke * (x13 * x32 + y23 * y31);
	K[0][2] = ke * (x21 * x32 + y12 * y23);

	K[1][0] = K[0][1];
	K[1][1] = ke * (x13 * x13 + y31 * y31);
	K[1][2] = ke * (x13 * x21 + y12 * y31);

	K[2][0] = K[0][2];
	K[2][1] = K[1][2];
	K[2][2] = ke * (x21 * x21 + y12 * y12);

	if (el == 0) {
		fprintf(stderr,
		        "TriangleSommerBC::turkelMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;
	int i;
	for (i = 0; i < 9; i++) d[i] /= r;

	return K;

}


void
TriangleSommerBC::getLocalCoordinates(CoordSet &cs, double xx[3], double yy[3], double zz[3]) const {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3];
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

	// Make xi local axis connecting nodes 1 and 2
	double l12 = sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0])
	                  + (z[1] - z[0]) * (z[1] - z[0]));


	xi[0] = (x[1] - x[0]) / l12;
	xi[1] = (y[1] - y[0]) / l12;
	xi[2] = (z[1] - z[0]) / l12;

	// Make eta local axis connecting nodes 1 and 3
	double l13 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0])
	                  + (z[2] - z[0]) * (z[2] - z[0]));

	eta[0] = (x[2] - x[0]) / l13;
	eta[1] = (y[2] - y[0]) / l13;
	eta[2] = (z[2] - z[0]) / l13;

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

}


double
TriangleSommerBC::getSize(CoordSet &cs) {
	return getArea(cs, nn);
}


void
TriangleSommerBC::getNormal(const CoordSet &cs, double normal[3]) const {

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

	double u1, u2, u3, v1, v2, v3, w1, w2, w3;
	u1 = x[1] - x[0];
	u2 = y[1] - y[0];
	u3 = z[1] - z[0];

	v1 = x[2] - x[0];
	v2 = y[2] - y[0];
	v3 = z[2] - z[0];

	w1 = u2 * v3 - u3 * v2;
	w2 = u3 * v1 - u1 * v3;
	w3 = u1 * v2 - u2 * v1;

	double l = sqrt(w1 * w1 + w2 * w2 + w3 * w3);

	if (l == 0) l = 1.0;

	normal[0] = w1 / l;
	normal[1] = w2 / l;
	normal[2] = w3 / l;
}


void TriangleSommerBC::get_basis(int r, int s, int t, double (*a)[3], double *e1, double *e2) {
	int i;
	for (i = 0; i < 3; i++) {
		e1[i] = a[s][i] - a[r][i];
		e2[i] = a[r][i] - a[t][i];
	}
}

double TriangleSommerBC::getArea(CoordSet &cs, int *nn) {

	// This is for the lumped mass
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



	// Define vectors:
	//                u: vector connecting nd1 and nd2
	//                v: vector connecting nd1 and nd3
	//                w: cross product between u and v

	double u1, u2, u3, v1, v2, v3, w1, w2, w3;

	u1 = x[1] - x[0];
	u2 = y[1] - y[0];
	u3 = z[1] - z[0];

	v1 = x[2] - x[0];
	v2 = y[2] - y[0];
	v3 = z[2] - z[0];

	w1 = u2 * v3 - u3 * v2;
	w2 = u3 * v1 - u1 * v3;
	w3 = u1 * v2 - u2 * v1;

	double area = 0.5 * sqrt(w1 * w1 + w2 * w2 + w3 * w3);

	return area;

}

void TriangleSommerBC::BT2(CoordSet &cs, double *e, double *f, double *g,
                           double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d) {

	int i, j, l, m;

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double a[3][3];
	a[0][0] = nd1.x;
	a[0][1] = nd1.y;
	a[0][2] = nd1.z;
	a[1][0] = nd2.x;
	a[1][1] = nd2.y;
	a[1][2] = nd2.z;
	a[2][0] = nd3.x;
	a[2][1] = nd3.y;
	a[2][2] = nd3.z;

	double Pet[3][3][2][2];
	double Pte[3][3][2][2];
	double b[3][2];

	for (m = 0; m < 3; m++)
		for (l = 0; l < 3; l++) {
			double e1[3], e2[3];
			if (l == 0) {
				int r = 0, s = 1, t = 2;
				get_basis(r, s, t, a, e1, e2);
			}
			if (l == 1) {
				int r = 1, s = 2, t = 0;
				get_basis(r, s, t, a, e1, e2);
			}
			if (l == 2) {
				int r = 2, s = 0, t = 1;
				get_basis(r, s, t, a, e1, e2);
			}

			double l1, l2;
			l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
			l2 = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
			for (i = 0; i < 3; i++) {
				e1[i] /= l1;
				e2[i] /= l2;
			}

			double dd = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
			b[l][0] = -(l2 + l1 * dd) / l1 / l2 / (1.0 - dd * dd);
			b[l][1] = (l1 + l2 * dd) / l1 / l2 / (1.0 - dd * dd);

			Pet[l][m][0][0] = tau1[m][0] * e1[0] + tau1[m][1] * e1[1] + tau1[m][2] * e1[2];
			Pet[l][m][1][0] = tau2[m][0] * e1[0] + tau2[m][1] * e1[1] + tau2[m][2] * e1[2];
			Pet[l][m][0][1] = tau1[m][0] * e2[0] + tau1[m][1] * e2[1] + tau1[m][2] * e2[2];
			Pet[l][m][1][1] = tau2[m][0] * e2[0] + tau2[m][1] * e2[1] + tau2[m][2] * e2[2];

			Pte[m][l][0][0] = Pet[l][m][0][0];
			Pte[m][l][1][0] = Pet[l][m][0][1];
			Pte[m][l][0][1] = Pet[l][m][1][0];
			Pte[m][l][1][1] = Pet[l][m][1][1];

		}

	ComplexD IRinv[3][2][2];
	for (m = 0; m < 3; m++) {
		ComplexD det = ComplexD(1.0, e[m] / k) * ComplexD(1.0, g[m] / k) -
		               ComplexD(0.0, f[m] / k) * ComplexD(0.0, f[m] / k);
		IRinv[m][0][0] = ComplexD(1.0, g[m] / k) / det;
		IRinv[m][1][1] = ComplexD(1.0, e[m] / k) / det;
		IRinv[m][0][1] = -ComplexD(0.0, f[m] / k) / det;
		IRinv[m][1][0] = -ComplexD(0.0, f[m] / k) / det;
	}

	ComplexD R[3][3];
	for (m = 0; m < 3; m++)
		for (j = 0; j < 3; j++) {
			R[m][j] = 0.0;
			for (l = 0; l < 3; l++) {
				// R^mj += b^j*Pte^lj*IRinv^l*RPet^ml*b^m
				double tmp1[2];
				ComplexD tmp2[2], tmp3[2];
				tmp1[0] = Pet[m][l][0][0] * b[m][0] + Pet[m][l][0][1] * b[m][1];
				tmp1[1] = Pet[m][l][1][0] * b[m][0] + Pet[m][l][1][1] * b[m][1];
				tmp2[0] = IRinv[l][0][0] * tmp1[0] + IRinv[l][0][1] * tmp1[1];
				tmp2[1] = IRinv[l][1][0] * tmp1[0] + IRinv[l][1][1] * tmp1[1];
				tmp3[0] = Pte[l][j][0][0] * tmp2[0] + Pte[l][j][0][1] * tmp2[1];
				tmp3[1] = Pte[l][j][1][0] * tmp2[0] + Pte[l][j][1][1] * tmp2[1];
				R[m][j] += b[j][0] * tmp3[0] + b[j][1] * tmp3[1];
			}
		}

	double area = getArea(cs, nn);
	for (m = 0; m < 3; m++)
		for (j = 0; j < 3; j++) {
			d[m * 3 + j] = area / 3.0 * R[m][j];
		}
/*
 fprintf(stderr, "BT2 term:\n");
 int jj=0;
 for(m=0;m<3;m++) {
   for(j=0;j<3;j++) {
     fprintf(stderr, "%f ", imag(d[jj]));
     jj++;
   }
   fprintf (stderr, "\n");
 }
 jj=0;
 for(m=0;m<3;m++) {
   for(j=0;j<3;j++) {
     fprintf(stderr, "%f ", real(d[jj]));
     jj++;
   }
   fprintf (stderr, "\n");
 }
*/
	if (el == 0) {
		fprintf(stderr,
		        "TriangleSommerBC::BC2: adjacent element not defined.\n");
		exit(-1);
	}
	double rho = el->getProperty()->rho;
	for (i = 0; i < 9; i++) d[i] /= rho;
}

void TriangleSommerBC::BT2n(CoordSet &cs, double *e, double *f, double *g,
                            double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d, int n) {

	int i, j, l, m;

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double a[3][3];
	a[0][0] = nd1.x;
	a[0][1] = nd1.y;
	a[0][2] = nd1.z;
	a[1][0] = nd2.x;
	a[1][1] = nd2.y;
	a[1][2] = nd2.z;
	a[2][0] = nd3.x;
	a[2][1] = nd3.y;
	a[2][2] = nd3.z;

	double Pet[3][3][2][2];
	double Pte[3][3][2][2];
	double b[3][2];

	for (m = 0; m < 3; m++)
		for (l = 0; l < 3; l++) {
			double e1[3], e2[3];
			if (l == 0) {
				int r = 0, s = 1, t = 2;
				get_basis(r, s, t, a, e1, e2);
			}
			if (l == 1) {
				int r = 1, s = 2, t = 0;
				get_basis(r, s, t, a, e1, e2);
			}
			if (l == 2) {
				int r = 2, s = 0, t = 1;
				get_basis(r, s, t, a, e1, e2);
			}

			double l1, l2;
			l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
			l2 = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
			for (i = 0; i < 3; i++) {
				e1[i] /= l1;
				e2[i] /= l2;
			}

			double dd = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
			b[l][0] = -(l2 + l1 * dd) / l1 / l2 / (1.0 - dd * dd);
			b[l][1] = (l1 + l2 * dd) / l1 / l2 / (1.0 - dd * dd);

			Pet[l][m][0][0] = tau1[m][0] * e1[0] + tau1[m][1] * e1[1] + tau1[m][2] * e1[2];
			Pet[l][m][1][0] = tau2[m][0] * e1[0] + tau2[m][1] * e1[1] + tau2[m][2] * e1[2];
			Pet[l][m][0][1] = tau1[m][0] * e2[0] + tau1[m][1] * e2[1] + tau1[m][2] * e2[2];
			Pet[l][m][1][1] = tau2[m][0] * e2[0] + tau2[m][1] * e2[1] + tau2[m][2] * e2[2];

			Pte[m][l][0][0] = Pet[l][m][0][0];
			Pte[m][l][1][0] = Pet[l][m][0][1];
			Pte[m][l][0][1] = Pet[l][m][1][0];
			Pte[m][l][1][1] = Pet[l][m][1][1];

		}

	ComplexD IRinv[3][2][2];
	for (m = 0; m < 3; m++) {
		ComplexD det = ComplexD(1.0, e[m] / k) * ComplexD(1.0, g[m] / k) -
		               ComplexD(0.0, f[m] / k) * ComplexD(0.0, f[m] / k);
		IRinv[m][0][0] = ComplexD(1.0, g[m] / k) / det;
		IRinv[m][1][1] = ComplexD(1.0, e[m] / k) / det;
		IRinv[m][0][1] = -ComplexD(0.0, f[m] / k) / det;
		IRinv[m][1][0] = -ComplexD(0.0, f[m] / k) / det;
	}

	ComplexD IRinv_n[3][2][2];
	for (m = 0; m < 3; m++) for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) IRinv_n[m][i][j] = IRinv[m][i][j];
	for (m = 0; m < 3; m++)
		for (k = 1; k <= n; ++k) {
			ComplexD c00 = IRinv_n[m][0][0];
			ComplexD c01 = IRinv_n[m][0][1];
			ComplexD c10 = IRinv_n[m][1][0];
			ComplexD c11 = IRinv_n[m][1][1];
			IRinv_n[m][0][0] = c00 * IRinv[m][0][0] + c01 * IRinv[m][1][0];
			IRinv_n[m][0][1] = c00 * IRinv[m][0][1] + c01 * IRinv[m][1][1];
			IRinv_n[m][1][0] = c10 * IRinv[m][0][0] + c11 * IRinv[m][1][0];
			IRinv_n[m][1][1] = c10 * IRinv[m][0][1] + c11 * IRinv[m][1][1];
		}

	ComplexD R[3][3];
	for (m = 0; m < 3; m++)
		for (j = 0; j < 3; j++) {
			R[m][j] = 0.0;
			for (l = 0; l < 3; l++) {
				// R^mj += b^j*Pte^lj*IRinv_n^l*Pet^ml*b^m
				double tmp1[2];
				ComplexD tmp2[2], tmp3[2];
				tmp1[0] = Pet[m][l][0][0] * b[m][0] + Pet[m][l][0][1] * b[m][1];
				tmp1[1] = Pet[m][l][1][0] * b[m][0] + Pet[m][l][1][1] * b[m][1];
				tmp2[0] = IRinv_n[l][0][0] * tmp1[0] + IRinv_n[l][0][1] * tmp1[1];
				tmp2[1] = IRinv_n[l][1][0] * tmp1[0] + IRinv_n[l][1][1] * tmp1[1];
				tmp3[0] = Pte[l][j][0][0] * tmp2[0] + Pte[l][j][0][1] * tmp2[1];
				tmp3[1] = Pte[l][j][1][0] * tmp2[0] + Pte[l][j][1][1] * tmp2[1];
				R[m][j] += b[j][0] * tmp3[0] + b[j][1] * tmp3[1];
			}
		}

	double area = getArea(cs, nn);
	for (m = 0; m < 3; m++)
		for (j = 0; j < 3; j++) {
			d[m * 3 + j] = area / 3.0 * R[m][j];
		}
/*
 fprintf(stderr, "BT2 term:\n");
 int jj=0;
 for(m=0;m<3;m++) {
   for(j=0;j<3;j++) {
     fprintf(stderr, "%f ", imag(d[jj]));
     jj++;
   }
   fprintf (stderr, "\n");
 }
 jj=0;
 for(m=0;m<3;m++) {
   for(j=0;j<3;j++) {
     fprintf(stderr, "%f ", real(d[jj]));
     jj++;
   }
   fprintf (stderr, "\n");
 }
*/
	if (el == 0) {
		fprintf(stderr,
		        "TriangleSommerBC::BT2n: adjacent element not defined.\n");
		exit(-1);
	}
	double rho = el->getProperty()->rho;
	for (i = 0; i < 9; i++) d[i] /= rho;
}

// Better BT2 for sphere
void TriangleSommerBC::sphereBT2(CoordSet &cs, double r, double k, ComplexD *d) {
	double x[3], y[3], z[3];
	getLocalCoordinates(cs, x, y, z);
	double dd[9];
	// _FORTRAN(trig6stif1)(x, y, 3, dd, 6);

	turkelMatrix(cs, dd);

	int j, m;
	for (m = 0; m < 3; m++)
		for (j = 0; j < 3; j++) {
			d[m * 3 + j] = 1.0 / ComplexD(1.0, 1.0 / r / k) * dd[m * 3 + j];
		}
	if (el == 0) {
		fprintf(stderr,
		        "TriangleSommerBC::sphereBT2: adjacent element not defined.\n");
		exit(-1);
	}
	double rho = el->getProperty()->rho;
	int i;
	for (i = 0; i < 9; i++) d[i] /= rho;
}

void TriangleSommerBC::GaussCoordinates(int Ngp, double *uPg, double *vPg, double *weight) const {
	if (Ngp == 4) {
		uPg[0] = 0.333333333333333;
		vPg[0] = 0.33333333333333;
		uPg[1] = 0.6;
		vPg[1] = 0.2;
		uPg[2] = 0.2;
		vPg[2] = 0.6;
		uPg[3] = 0.2;
		vPg[3] = 0.2;
		weight[0] = -0.28125;
		weight[1] = weight[2] = weight[3] = 0.260416666666666;
	} else if (Ngp == 7) {
		uPg[0] = 0.333333333333333;
		vPg[0] = 0.333333333333333;
		uPg[1] = 0.059715871789770;
		vPg[1] = 0.470142064105115;
		uPg[2] = 0.470142064105115;
		vPg[2] = 0.059715871789770;
		uPg[3] = 0.470142064105115;
		vPg[3] = 0.470142064105115;
		uPg[4] = 0.797426985353087;
		vPg[4] = 0.101286507323456;
		uPg[5] = 0.101286507323456;
		vPg[5] = 0.797426985353087;
		uPg[6] = 0.101286507323456;
		vPg[6] = 0.101286507323456;
		weight[0] = 0.1125;
		weight[1] = weight[2] = weight[3] = 0.066197076394253;
		weight[4] = weight[5] = weight[6] = 0.0629695902724135;
	} else if (Ngp == 13) {
		uPg[0] = 0.333333333333333;
		vPg[0] = 0.333333333333333;
		uPg[1] = 0.479308067841920;
		vPg[1] = 0.260345966079040;
		uPg[2] = 0.260345966079040;
		vPg[2] = 0.479308067841920;
		uPg[3] = 0.260345966079040;
		vPg[3] = 0.260345966079040;
		uPg[4] = 0.869739794195568;
		vPg[4] = 0.065130102902216;
		uPg[5] = 0.065130102902216;
		vPg[5] = 0.869739794195568;
		uPg[6] = 0.065130102902216;
		vPg[6] = 0.065130102902216;
		uPg[7] = 0.048690315425316;
		vPg[7] = 0.312865496004874;
		uPg[8] = 0.312865496004874;
		vPg[8] = 0.048690315425316;
		uPg[9] = 0.048690315425316;
		vPg[9] = 0.638444188569810;
		uPg[10] = 0.638444188569810;
		vPg[10] = 0.048690315425316;
		uPg[11] = 0.312865496004874;
		vPg[11] = 0.638444188569810;
		uPg[12] = 0.638444188569810;
		vPg[12] = 0.312865496004874;
		weight[0] = -0.07478502223384099;
		weight[1] = weight[2] = weight[3] = 0.087807628716604;
		weight[4] = weight[5] = weight[6] = 0.026673617804419;
		weight[7] = weight[8] = weight[9] = weight[10] = weight[11] = weight[12] = 0.0385568804451285;
	} else if (Ngp == 25) {
		double c11 = 0.028844733232685;
		double c12 = 0.485577633383657;
		double c21 = 0.781036849029926;
		double c22 = 0.109481575485037;
		double c31 = 0.141909219414880;
		double c32 = 0.307939838764121;
		double c33 = 0.550352941820999;
		double c41 = 0.025003534762686;
		double c42 = 0.246672560639903;
		double c43 = 0.728323904597411;
		double c51 = 0.009540815400299;
		double c52 = 0.066803251012200;
		double c53 = 0.923655933587500;
		uPg[0] = 0.333333333333333;
		vPg[0] = 0.333333333333333;
		uPg[1] = c11;
		vPg[1] = c12;
		uPg[2] = c12;
		vPg[2] = c11;
		uPg[3] = c12;
		vPg[3] = c12;
		uPg[4] = c21;
		vPg[4] = c22;
		uPg[5] = c22;
		vPg[5] = c21;
		uPg[6] = c22;
		vPg[6] = c22;
		uPg[7] = c31;
		vPg[7] = c32;
		uPg[8] = c32;
		vPg[8] = c31;
		uPg[9] = c31;
		vPg[9] = c33;
		uPg[10] = c33;
		vPg[10] = c31;
		uPg[11] = c32;
		vPg[11] = c33;
		uPg[12] = c33;
		vPg[12] = c32;
		uPg[13] = c41;
		vPg[13] = c42;
		uPg[14] = c42;
		vPg[14] = c41;
		uPg[15] = c41;
		vPg[15] = c43;
		uPg[16] = c43;
		vPg[16] = c41;
		uPg[17] = c42;
		vPg[17] = c43;
		uPg[18] = c43;
		vPg[18] = c42;
		uPg[19] = c51;
		vPg[19] = c52;
		uPg[20] = c52;
		vPg[20] = c51;
		uPg[21] = c51;
		vPg[21] = c53;
		uPg[22] = c53;
		vPg[22] = c51;
		uPg[23] = c52;
		vPg[23] = c53;
		uPg[24] = c53;
		vPg[24] = c52;
		weight[0] = 0.045408995191377;
		weight[1] = weight[2] = weight[3] = 0.0183629788782335;
		weight[4] = weight[5] = weight[6] = 0.022660529717764;
		weight[7] = weight[8] = weight[9] = weight[10] = weight[11] = weight[12] = 0.03637895842271;
		weight[13] = weight[14] = weight[15] = weight[16] = weight[17] = weight[18] = 0.0141636212655285;
		weight[19] = weight[20] = weight[21] = weight[22] = weight[23] = weight[24] = 0.0047108334818665;
	} else {
		fprintf(stderr, "  TriangleSommerBC.C number of gauss points non implemented\n");
		for (int i = 0; i < Ngp; i++) {
			uPg[i] = 0.0;
			vPg[i] = 0.0;
			weight[i] = 1.0;
		}
	}
}

//function that depend on the order of the absorbing boundary condition
#include<Driver.d/Domain.h>

int *
TriangleSommerBC::dofs(DofSetArray &dsa, int *p) const  {
	double numD = -2.0;
	if (dom) numD = dom->solInfo().ATDARBFlag;
	if (numD == 1.5) {
		if (p == 0) p = new int[8];
		dsa.number(nn[0], DofSet::Helm, p);
		dsa.number(nn[0], DofSet::IntPress, p + 3);
		dsa.number(nn[1], DofSet::Helm, p + 1);
		dsa.number(nn[1], DofSet::IntPress, p + 4);
		dsa.number(nn[2], DofSet::Helm, p + 2);
		dsa.number(nn[2], DofSet::IntPress, p + 5);
	} else {
		if (p == 0) p = new int[3];
		dsa.number(nn[0], DofSet::Helm, p);
		dsa.number(nn[1], DofSet::Helm, p + 1);
		dsa.number(nn[2], DofSet::Helm, p + 2);
	}
	return p;
}

int TriangleSommerBC::numDofs() const {
	double numD = -2.0;
	if (dom) numD = dom->solInfo().ATDARBFlag;
	if (numD == 1.5)
		return 6;
	else
		return 3;
}

void TriangleSommerBC::markDofs(DofSetArray &dsa) const {
	double numD = (dom) ? dom->solInfo().ATDARBFlag : 0.0;
	if (numD == 1.5) {
		dsa.mark(nn[0], DofSet::Helm);
		dsa.mark(nn[1], DofSet::Helm);
		dsa.mark(nn[2], DofSet::Helm);
		dsa.mark(nn[0], DofSet::IntPress);
		dsa.mark(nn[1], DofSet::IntPress);
		dsa.mark(nn[2], DofSet::IntPress);
	} else {
		dsa.mark(nn[0], DofSet::Helm);
		dsa.mark(nn[1], DofSet::Helm);
		dsa.mark(nn[2], DofSet::Helm);
	}
}

void TriangleSommerBC::SurfaceRefinement(int nNo, double *x, double *y, double *z, double *xx, double *yy, double *zz) const {
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
	xart[0] = xx[0] + xx[1] * dd + xx[2] * dd;
	xart[1] = xx[0] + xx[1] * (1 - dd) + xx[2] * (dd / 2);
	xart[2] = xx[0] + xx[1] * (dd / 2) + xx[2] * (1 - dd);
	yart[0] = yy[0] + yy[1] * dd + yy[2] * dd;
	yart[1] = yy[0] + yy[1] * (1 - dd) + yy[2] * (dd / 2);
	yart[2] = yy[0] + yy[1] * (dd / 2) + yy[2] * (1 - dd);
	zart[0] = zz[0] + zz[1] * dd + zz[2] * dd;
	zart[1] = zz[0] + zz[1] * (1 - dd) + zz[2] * (dd / 2);
	zart[2] = zz[0] + zz[1] * (dd / 2) + zz[2] * (1 - dd);
	//put the points in the tangenetial plan: Mart st MMart.norm=0 projection
	double alpha, beta;
	for (i = 0; i < nNo; i++) {
		alpha = (xart[i] - x[i]) * nx[i] + (yart[i] - y[i]) * ny[i] + (zart[i] - z[i]) * nz[i];
		beta = nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i];
		xart[i] -= alpha * nx[i] / beta;
		yart[i] -= alpha * ny[i] / beta;
		zart[i] -= alpha * nz[i] / beta;
	}

	// 6 nodes interpolation with x=x0+x1*u+x2*v+(u*(1-u)+v*(1-v))*(x3+x4*u+x5*v), idem for y and z
	double *temp = new double[ttnNo * ttnNo];
	FullSquareMatrix iM(ttnNo, (double *) temp);
	double *u = new double[ttnNo], *v = new double[ttnNo];
	int *vectemp = new int[ttnNo];
	u[0] = 0.0;
	u[1] = 1.0;
	u[2] = 0.0;
	u[3] = dd;
	u[4] = 1.0 - dd;
	u[5] = dd / 2;
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 1.0;
	v[3] = dd;
	v[4] = dd / 2;
	v[5] = 1.0 - dd;
	for (i = 0; i < nNo; i++) {
		xx[i] = x[i];
		xx[i + nNo] = xart[i];
		yy[i] = y[i];
		yy[i + nNo] = yart[i];
		zz[i] = z[i];
		zz[i + nNo] = zart[i];
	}
	double co = 0;
	for (i = 0; i < ttnNo; i++) {
		iM[i][0] = 1.0;
		iM[i][1] = u[i];
		iM[i][2] = v[i];
		co = u[i] * (1.0 - u[i]) + v[i] * (1.0 - v[i]);
		iM[i][3] = co;
		iM[i][4] = co * u[i];
		iM[i][5] = co * v[i];
	}
	//solve iM*x = xx;
	int fl = ludcmp(iM, ttnNo, vectemp); // LU decomposition
	if (fl == 0) {
		lubksb(iM, ttnNo, vectemp, xx); // LU backsubstitution
		lubksb(iM, ttnNo, vectemp, yy);
		lubksb(iM, ttnNo, vectemp, zz);
	} else
		for (i = 0; i < ttnNo; i++) {
			fprintf(stderr, "  TriangleSommer.C, singular matrix");
			xx[i] = 0.0;
			yy[i] = 0.0;
			zz[i] = 0.0;
		}
	delete[] nx;
	nx = NULL;
	delete[] ny;
	ny = NULL;
	delete[] nz;
	nz = NULL;
	delete[] xart;
	xart = NULL;
	delete[] yart;
	yart = NULL;
	delete[] zart;
	zart = NULL;
	delete[] temp;
	temp = NULL;
	delete[] u;
	u = NULL;
	delete[] v;
	v = NULL;
	delete[] vectemp;
	vectemp = NULL;
}

FullSquareMatrix
TriangleSommerBC::refinedSommerMatrix(CoordSet &cs, double *d) {
	int i, j;
	int nNo = numNodes();
	int ttnNo = 2 * nNo;
	// reference: do Carmo, differential geometries of curves and surfaces p. 101
	// get the nodes coordinates
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
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

	double *xx = new double[ttnNo], *yy = new double[ttnNo], *zz = new double[ttnNo];
	// Express the surface depending on 2 parameters (u,v)
	// x = xx0 + xx1*u + xx2*v
	xx[0] = x[0];
	xx[1] = (-x[0] + x[1]);
	xx[2] = (-x[0] + x[2]);
	yy[0] = y[0];
	yy[1] = (-y[0] + y[1]);
	yy[2] = (-y[0] + y[2]);
	zz[0] = z[0];
	zz[1] = (-z[0] + z[1]);
	zz[2] = (-z[0] + z[2]);
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
	const int Ngp = 13;//Number of Gauss Points, choice: 4, 7, 13, 25
	double uPg[Ngp], vPg[Ngp], weight[Ngp];
	GaussCoordinates(Ngp, uPg, vPg, weight);

	//define the nodes' caracteristic functions: no=coef0 + coefu0*u + coefv0*v, ...
	double *coef = new double[nNo], *coefu = new double[nNo], *coefv = new double[nNo];
	coef[0] = 1.0;
	coef[1] = 0.0;
	coef[2] = 0.0;
	coefu[0] = -1.0;
	coefu[1] = 1.0;
	coefu[2] = 0.0;
	coefv[0] = -1.0;
	coefv[1] = 0.0;
	coefv[2] = 1.0;

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
		//build paraVect = vector of the values of the parameter
		UU[3] = -2 * uPg[igp];
		UU[5] = -2 * uPg[igp] * vPg[igp];
		UU[4] = -3 * uPg[igp] * uPg[igp] + 2 * uPg[igp] + vPg[igp] * (1 - vPg[igp]);
		VV[3] = -2 * vPg[igp];
		VV[4] = -2 * vPg[igp] * uPg[igp];
		VV[5] = uPg[igp] * (1 - uPg[igp]) + 2 * vPg[igp] - 3 * vPg[igp] * vPg[igp];
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
				n1 = coef[i] + coefu[i] * uPg[igp] + coefv[i] * vPg[igp];
				n2 = coef[j] + coefu[j] * uPg[igp] + coefv[j] * vPg[igp];
				val = n1 * n2 * jj;
				// gauss integration: build d
				d[j + i * nNo] += val * weight[igp];//order of d does not matter as the matrix is symetric
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
	delete[] coef;
	coef = NULL;
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
TriangleSommerBC::surfStiffMatrix(CoordSet &cs, double *d)
{ 
  int i, j;
  int nNo = numNodes();
  int ttnNo = 2*nNo;
  // reference: do Carmo, differential geometries of curves and surfaces p. 101
  // get the nodes coordinates
  Node nd1 = cs.getNode(nn[0]);
  Node nd2 = cs.getNode(nn[1]);
  Node nd3 = cs.getNode(nn[2]);
  double *x=new double[nNo], *y=new double[nNo], *z=new double[nNo];
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double *xx=new double[ttnNo], *yy=new double[ttnNo], *zz=new double[ttnNo];
  // Express the surface depending on 2 parameters (u,v)
  // x = xx0 + xx1*u + xx2*v 
  xx[0]=x[0];
  xx[1]=(-x[0]+x[1]);
  xx[2]=(-x[0]+x[2]);
  yy[0]=y[0];
  yy[1]=(-y[0]+y[1]);
  yy[2]=(-y[0]+y[2]);
  zz[0]=z[0];
  zz[1]=(-z[0]+z[1]);
  zz[2]=(-z[0]+z[2]);
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
  const int Ngp = 13;//Number of Gauss Points, choice: 4, 7, 13, 25
  double uPg[Ngp], vPg[Ngp], weight[Ngp];
  GaussCoordinates(Ngp,uPg,vPg,weight);
  
  //define the nodes' caracteristic functions: no=coef0 + coefu0*u + coefv0*v, ...
  double *coef=new double[nNo],*coefu=new double[nNo], *coefv=new double[nNo];
  coef[0] = 1.0 ; coef[1] = 0.0 ; coef[2] = 0.0 ;
  coefu[0] = -1.0 ; coefu[1] = 1.0 ; coefu[2] = 0.0 ;
  coefv[0] = -1.0 ; coefv[1] = 0.0 ; coefv[2] = 1.0 ;

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
    //build paraVect = vector of the values of the parameter
    UU[3] = -2*uPg[igp]; UU[5] = -2*uPg[igp]*vPg[igp];
      UU[4] = -3*uPg[igp]*uPg[igp]+2*uPg[igp]+vPg[igp]*(1-vPg[igp]);
    VV[3] = -2*vPg[igp]; VV[4] = -2*vPg[igp]*uPg[igp];
     VV[5] = uPg[igp]*(1-uPg[igp])+2*vPg[igp]-3*vPg[igp]*vPg[igp];
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
        dn1u = coefu[i] ; dn1v = coefv[i] ;
        dn2u = coefu[j] ; dn2v = coefv[j] ;
        val = (dn1u*G*dn2u-dn1u*F*dn2v-dn1v*F*dn2u+dn1v*E*dn2v)/(E*G-F*F)*jj;
        // gauss integration: build d
        d[j+i*nNo] += val*weight[igp];//order of d does not matter as the matrix is symetric
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
  delete [] coef; coef = NULL;
  delete [] coefu; coefu = NULL;
  delete [] coefv; coefv = NULL;
  delete [] UU; UU = NULL;
  delete [] VV; VV = NULL;

  return surfStiffM;
}
*/

FullSquareMatrix
TriangleSommerBC::HSommerMatrix(const CoordSet &cs, double *d) const {
	int i, j, k;
	int nNo = numNodes();
	int ttnNo = 2 * nNo;
	// get the nodes coordinates
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
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
	// Express the surface depending on 2 parameters (u,v)
	// x = xx0 + xx1*u + xx2*v
	double *xx = new double[ttnNo], *yy = new double[ttnNo], *zz = new double[ttnNo];
	xx[0] = x[0];
	xx[1] = (-x[0] + x[1]);
	xx[2] = (-x[0] + x[2]);
	yy[0] = y[0];
	yy[1] = (-y[0] + y[1]);
	yy[2] = (-y[0] + y[2]);
	zz[0] = z[0];
	zz[1] = (-z[0] + z[1]);
	zz[2] = (-z[0] + z[2]);
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
	const int Ngp = 13;//Number of Gauss Points, choice: 4, 7, 13, 25
	double uPg[Ngp], vPg[Ngp], weight[Ngp];
	GaussCoordinates(Ngp, uPg, vPg, weight);

	//define the nodes' caracteristic functions: no=(1+coefu[0]*u)*(1+coefv[0]*v)/4, ...
	double *coef = new double[nNo], *coefu = new double[nNo], *coefv = new double[nNo];
	coef[0] = 1.0;
	coef[1] = 0.0;
	coef[2] = 0.0;
	coefu[0] = -1.0;
	coefu[1] = 1.0;
	coefu[2] = 0.0;
	coefv[0] = -1.0;
	coefv[1] = 0.0;
	coefv[2] = 1.0;

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
		//build paraVect = vector of the values of the parameter
		UU[3] = -2 * uPg[igp];
		UU[5] = -2 * uPg[igp] * vPg[igp];
		UU[4] = -3 * uPg[igp] * uPg[igp] + 2 * uPg[igp] + vPg[igp] * (1 - vPg[igp]);
		VV[3] = -2 * vPg[igp];
		VV[4] = -2 * vPg[igp] * uPg[igp];
		VV[5] = uPg[igp] * (1 - uPg[igp]) + 2 * vPg[igp] - 3 * vPg[igp] * vPg[igp];
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
					n1 = coef[i] + coefu[i] * uPg[igp] + coefv[i] * vPg[igp];
					n2 = coef[j] + coefu[j] * uPg[igp] + coefv[j] * vPg[igp];
					n3 = coef[k] + coefu[k] * uPg[igp] + coefv[k] * vPg[igp];
					val = HH[k] * n1 * n2 * n3 * jj;
					// gauss integration: build d
					d[j + i * nNo] += val * weight[igp];
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
	delete[] coef;
	coef = NULL;
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
// not used anymore - developped for ABC time domain - JFD 
FullSquareMatrix
TriangleSommerBC::HKSommerMatrix(CoordSet &cs, double *d)
{
  int i, j, k;
  int nNo = numNodes();
  int ttnNo = 2*nNo;
  // get the nodes coordinates
  Node nd1 = cs.getNode(nn[0]);
  Node nd2 = cs.getNode(nn[1]);
  Node nd3 = cs.getNode(nn[2]);
  double *x=new double[nNo], *y=new double[nNo], *z=new double[nNo];
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  // Express the surface depending on 2 parameters (u,v)
  // x = xx0 + xx1*u + xx2*v 
  double *xx=new double[ttnNo], *yy=new double[ttnNo], *zz=new double[ttnNo];
  xx[0]=x[0];
  xx[1]=(-x[0]+x[1]);
  xx[2]=(-x[0]+x[2]);
  yy[0]=y[0];
  yy[1]=(-y[0]+y[1]);
  yy[2]=(-y[0]+y[2]);
  zz[0]=z[0];
  zz[1]=(-z[0]+z[1]);
  zz[2]=(-z[0]+z[2]);
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
  const int Ngp = 13;//Number of Gauss Points, choice: 4, 7, 13, 25
  double uPg[Ngp], vPg[Ngp], weight[Ngp];
  GaussCoordinates(Ngp,uPg,vPg,weight);

  //define the nodes' caracteristic functions: no=(1+coefu[0]*u)*(1+coefv[0]*v)/4, ...
  double *coef=new double[nNo],*coefu=new double[nNo], *coefv=new double[nNo];
  coef[0] = 1.0 ; coef[1] = 0.0 ; coef[2] = 0.0 ;
  coefu[0] = -1.0 ; coefu[1] = 1.0 ; coefu[2] = 0.0 ;
  coefv[0] = -1.0 ; coefv[1] = 0.0 ; coefv[2] = 1.0 ;

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
    //build paraVect = vector of the values of the parameter
    UU[3] = -2*uPg[igp]; UU[5] = -2*uPg[igp]*vPg[igp];
      UU[4] = -3*uPg[igp]*uPg[igp]+2*uPg[igp]+vPg[igp]*(1-vPg[igp]);
    VV[3] = -2*vPg[igp]; VV[4] = -2*vPg[igp]*uPg[igp];
      VV[5] = uPg[igp]*(1-uPg[igp])+2*vPg[igp]-3*vPg[igp]*vPg[igp];
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
    j0 = Xuy*Xvz-Xuz*Xvy;
    j1 = Xuz*Xvx-Xux*Xvz;
    j2 = Xux*Xvy-Xuy*Xvx;
    jj = sqrt(j0*j0+j1*j1+j2*j2);
    //fill the d matrix
    for (i = 0 ; i < nNo ; i++) {
      for (j = 0 ; j < nNo ; j++) {
        for (k = 0 ; k < nNo ; k++) {
          // values of the function to be integrated at the gauss point (same order)
          n1 = coef[i]+coefu[i]*uPg[igp]+coefv[i]*vPg[igp];
          n2 = coef[j]+coefu[j]*uPg[igp]+coefv[j]*vPg[igp];
          n3 = coef[k]+coefu[k]*uPg[igp]+coefv[k]*vPg[igp];
          val = (HH[k]*HH[k]-KK[k])*n1*n2*n3*jj;
          // gauss integration: build d
          d[j+i*nNo] += val*weight[igp];
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
  delete [] coef; coef = NULL;
  delete [] coefu; coefu = NULL;
  delete [] coefv; coefv = NULL;
  delete [] HH; HH = NULL;
  delete [] KK; KK = NULL;
  delete [] UU; UU = NULL;
  delete [] VV; VV = NULL;

  return HKSommerM;
}
*/

void TriangleSommerBC::wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) {
	IsoParamUtilsTetra ipu(2);
	int ordersq = ipu.getordersq();
	double *xyz = (double *) alloca(sizeof(double) * 3 * ordersq);
	cs.getCoordinates(nn, ordersq, xyz, xyz + ordersq, xyz + 2 * ordersq);

	double *d = (double *) alloca(sizeof(double) * 3 * ordersq * ordersq);
	WetInterfaceGalFunction f(ordersq, d);
	ipu.zeroOut<double>(3 * ordersq * ordersq, d);
	int gorder = 13;
	ipu.surfSurfInt3d(xyz, f, gorder);

	int i, j = -1;
	for (i = 0; i < ordersq; i++) {
		if (nn[i] == nd) {
			j = i;
			break;
		}
	}
	if (j == -1) {
		fprintf(stderr, "Error in TriangleSommerBC::wetInterfaceLMPC\n");
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

