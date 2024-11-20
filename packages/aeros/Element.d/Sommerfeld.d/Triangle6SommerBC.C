#include        <cstdio>
#include        <cmath>
#include        <Utils.d/dbg_alloca.h>

#include        <Element.d/Sommerfeld.d/Triangle6SommerBC.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/linkfc.h>


extern "C" {
void    _FORTRAN(trig6stif1)(double *, double *, const int &, double *, const int &);
};
extern "C" {
void    _FORTRAN(trig6mass1)(double *, double *, const int &, double *, const int &);
};


// Sommerfeld b.c. contribution for 2-D elements
// Note that the boundary condition will be always 
// on an element edge. 

Triangle6SommerBC::Triangle6SommerBC(int n1, int n2, int n3, int n4, int n5, int n6, Element *_el, int etype)
		: SommerElement(_el) {
	nn[0] = n1;
	nn[1] = n2;
	nn[2] = n3;
	nn[3] = n4;
	nn[4] = n5;
	nn[5] = n6;
	setElementType(etype);
 sFlag = false;
 soundSpeed = 0.0;

}


Triangle6SommerBC *Triangle6SommerBC::clone() {
	Triangle6SommerBC *se = new Triangle6SommerBC(nn[0], nn[1], nn[2],
	                                              nn[3], nn[4], nn[5], el);
	se->el2 = el2;
	se->dom = dom;
 se->sFlag = sFlag;
	return se;
}


void Triangle6SommerBC::flipNormal() {
	int *nds = (int *) dbg_alloca(sizeof(int) * numNodes());
	nodes(nds);
	nn[0] = nds[0];
	nn[1] = nds[2];
	nn[2] = nds[1];
	nn[3] = nds[5];
	nn[4] = nds[4];
	nn[5] = nds[3];
}

void Triangle6SommerBC::computedxdxi(double x[6], double y[6], Matrix22 *dxdxi, double *det) {

	int i, j, k, nint = 12;

	double coord[6][2];
	for (i = 0; i < 6; i++) {
		coord[i][0] = x[i];
		coord[i][1] = y[i];
	}

	for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++) {
			for (i = 0; i < nint; i++) {
				dxdxi[i][j][k] = 0.0;
				int m;
				for (m = 0; m < 6; m++)
					dxdxi[i][j][k] += tri6_derivatives[i][m][k] * coord[m][j];
			}
		}

	for (i = 0; i < nint; i++) {
		Matrix22 &x = dxdxi[i];
		det[i] = x[0][0] * x[1][1] - x[0][1] * x[1][0];
	}

}

void
Triangle6SommerBC::sommerMatrixEllipsoid(CoordSet &cs, double kappa, double H[3], double K[3], ComplexD *d) {
/*
     double a = 1.5; double c = 5.5;
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
     double xx = (nd1.x+nd2.x+nd3.x)/3.0;
     double phi = acos(xx/5.5);
double bl =
(a*a*cos(phi)*cos(phi) + c*c*sin(phi)*sin(phi))*
    ((-3.0*c*
(-2.0*a*a*cos(phi)*sin(phi) +
           2.0*c*c*cos(phi)*sin(phi))*(-2.0*a*a*cos(phi)*sin(phi) +
           2.0*c*c*cos(phi)*sin(phi))
)/
       (2.*a*pow(a*a*cos(phi)*cos(phi) +
           c*c*sin(phi)*sin(phi),2.5)) +
      (c*(-2.0*a*a*cos(phi)*cos(phi) + 2.0*c*c*cos(phi)*cos(phi) +
           2.0*a*a*sin(phi)*sin(phi) - 2.0*c*c*sin(phi)*sin(phi)))/
       (2.*a*pow(a*a*cos(phi)*cos(phi) +
           c*c*sin(phi)*sin(phi),1.5)) +
      (15*c*
(-2*a*a*cos(phi)*sin(phi) +
           2.0*c*c*cos(phi)*sin(phi))*(-2.0*a*a*cos(phi)*sin(phi) +
           2.0*c*c*cos(phi)*sin(phi))
*
         (a*a + a*a*cos(phi)*cos(phi) +
           c*c*sin(phi)*sin(phi)))/
       (8.*a*pow(a*a*cos(phi)*cos(phi) +
           c*c*sin(phi)*sin(phi),3.5)) -
      (3.0*c*(-2.0*a*a*cos(phi)*cos(phi) +
           2.0*c*c*cos(phi)*cos(phi) + 2.0*a*a*sin(phi)*sin(phi) -
           2.0*c*c*sin(phi)*sin(phi))*
         (a*a + a*a*cos(phi)*cos(phi) +
           c*c*sin(phi)*sin(phi)))/
       (4.*a*pow(a*a*cos(phi)*cos(phi) +
           c*c*sin(phi)*sin(phi),2.5)));
*/
//    fprintf(stderr,"%f %f %f\n",bl/4.0/kappa/kappa,cos(phi),sin(phi));

	double x[6], y[6], z[6];
	getLocalCoordinates(cs, x, y, z);

	double dxdxi[12][2][2], det[12];
	computedxdxi(x, y, dxdxi, det);

	FullSquareMatrixC sm(6, d);
	int i, j, k, nint = 12;
	for (j = 0; j < 6; j++)
		for (k = j; k < 6; k++) {
			sm[j][k] = 0.0;
			for (i = 0; i < nint; i++) {
				double tmp1 = 0.5 * tri6_weights[i] * tri6_values[i][j] * tri6_values[i][k] * det[i];
				int m;
				double HH = 0.0;
				double KK = 0.0;
				for (m = 0; m < 3; m++) {
					HH += H[m] * tri3_values[i][m];
					KK += K[m] * tri3_values[i][m];
				}

				ComplexD tmp2 = ComplexD(HH, 0.0) + ComplexD(0.0, 1.0 / 2.0 / kappa) /
				                                    ComplexD(1.0, 2 * HH / kappa) *
				                                    (KK - HH * HH)/*- ComplexD(bl/4.0/kappa/kappa,0.0)*/;
				sm[j][k] += tmp1 * tmp2;
			}
		}

	for (j = 0; j < 6; j++) for (k = 0; k < j; k++) sm[j][k] = sm[k][j];

	if (el == 0) {
		fprintf(stderr,
		        "Triangle6SommerBC::sommerMatrixEllipsoid: adjacent element not defined.\n");
		exit(-1);
	}
	double rho = el->getProperty()->rho;
	for (j = 0; j < 6; j++) for (k = 0; k < 6; k++) sm[j][k] /= rho;


}


void
Triangle6SommerBC::getLocalCoordinatesNew(CoordSet &cs, double T[3][3], double xi1[3], double xi2[3], double xi3[3]) {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	Node nd5 = cs.getNode(nn[4]);
	Node nd6 = cs.getNode(nn[5]);

	double x[6], y[6], z[6];
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
	x[4] = nd5.x;
	y[4] = nd5.y;
	z[4] = nd5.z;
	x[5] = nd6.x;
	y[5] = nd6.y;
	z[5] = nd6.z;

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

	T[0][0] = xi[0];
	T[0][1] = xi[1];
	T[0][2] = xi[2];
	T[1][0] = eta[0];
	T[1][1] = eta[1];
	T[1][2] = eta[2];
	T[2][0] = nu[0];
	T[2][1] = nu[1];
	T[2][2] = nu[2];

	xi1[0] = xi[0];
	xi1[1] = xi[1];
	xi1[2] = xi[2];
	xi2[0] = eta[0];
	xi2[1] = eta[1];
	xi2[2] = eta[2];
	xi3[0] = nu[0];
	xi3[1] = nu[1];
	xi3[2] = nu[2];
}


FullSquareMatrix
Triangle6SommerBC::sommerMatrix(CoordSet &cs, double *d) const {
	double x[6], y[6], z[6];

	// Next two lines for the consistent mass
	// Note that the global coordinates are get 
	// inside getLocalCoordinates(*) and after 
	// that x, y and z will contain the local 
	// coordinates

	getLocalCoordinates(cs, x, y, z);

	_FORTRAN(trig6mass1)(x, y, 7, d, 6);

	if (el == 0) {
		fprintf(stderr,
		        "Triangle6SommerBC::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;

	int i;
	for (i = 0; i < 36; i++) d[i] = -d[i] / r;

	FullSquareMatrix sommerM(6, d);

	return sommerM;
}

FullSquareMatrix
Triangle6SommerBC::turkelMatrix(CoordSet &cs, double *d) const {
	double x[6], y[6], z[6];

	// Next two lines for the consistent mass
	// Note that the global coordinates are get
	// inside getLocalCoordinates(*) and after
	// that x, y and z will contain the local
	// coordinates

	getLocalCoordinates(cs, x, y, z);

	_FORTRAN(trig6mass1)(x, y, 3, d, 6);

	if (el == 0) {
		fprintf(stderr,
		        "Triangle6SommerBC::sommerMatrix: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;
	int i;
	for (i = 0; i < 36; i++) d[i] /= r;

	FullSquareMatrix sommerM(6, d);

	return sommerM;

}


void
Triangle6SommerBC::getLocalCoordinates(CoordSet &cs, double xx[6], double yy[6], double zz[6]) const {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	Node nd5 = cs.getNode(nn[4]);
	Node nd6 = cs.getNode(nn[5]);

	double x[6], y[6], z[6];
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
	x[4] = nd5.x;
	y[4] = nd5.y;
	z[4] = nd5.z;
	x[5] = nd6.x;
	y[5] = nd6.y;
	z[5] = nd6.z;

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
	x[3] -= x[0];
	y[3] -= y[0];
	z[3] -= z[0];
	x[4] -= x[0];
	y[4] -= y[0];
	z[4] -= z[0];
	x[5] -= x[0];
	y[5] -= y[0];
	z[5] -= z[0];
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

	xx[4] = T[0][0] * x[4] + T[0][1] * y[4] + T[0][2] * z[4];
	yy[4] = T[1][0] * x[4] + T[1][1] * y[4] + T[1][2] * z[4];
	zz[4] = T[2][0] * x[4] + T[2][1] * y[4] + T[2][2] * z[4];

	xx[5] = T[0][0] * x[5] + T[0][1] * y[5] + T[0][2] * z[5];
	yy[5] = T[1][0] * x[5] + T[1][1] * y[5] + T[1][2] * z[5];
	zz[5] = T[2][0] * x[5] + T[2][1] * y[5] + T[2][2] * z[5];

}


void
Triangle6SommerBC::getNormal(const CoordSet &cs, double normal[3]) const {

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


void Triangle6SommerBC::get_basis(int r, int s, int t, double (*a)[3], double *e1, double *e2) {
	int i;
	for (i = 0; i < 3; i++) {
		e1[i] = a[s][i] - a[r][i];
		e2[i] = a[r][i] - a[t][i];
	}
}

double Triangle6SommerBC::getArea(CoordSet &cs, int *nn) {

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

void Triangle6SommerBC::BT2(CoordSet &cs, double *e, double *f, double *g,
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
	for (m = 0; m < 36; m++) d[m] = ComplexD(0.0, 0.0);
	for (m = 0; m < 3; m++)
		for (j = 0; j < 3; j++) {
			d[m * 6 + j] = area / 3.0 * R[m][j];
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
		        "Triangle6SommerBC::BT2: adjacent element not defined.\n");
		exit(-1);
	}
	double r = el->getProperty()->rho;
	for (i = 0; i < 36; i++) d[i] /= r;
}

// Better BT2 for sphere
void Triangle6SommerBC::sphereBT2(CoordSet &cs, double r, double k, ComplexD *d) {
	double x[6], y[6], z[6];
	getLocalCoordinates(cs, x, y, z);
	double dd[36];
	_FORTRAN(trig6stif1)(x, y, 3, dd, 6);

	int j, m;
	for (m = 0; m < 6; m++)
		for (j = 0; j < 6; j++) {
			d[m * 6 + j] = 1.0 / ComplexD(1.0, 1.0 / r / k) * dd[m * 6 + j];
		}

	if (el == 0) {
		fprintf(stderr,
		        "Triangle6SommerBC::sphereBT2: adjacent element not defined.\n");
		exit(-1);
	}
	double rho = el->getProperty()->rho;
	int i;
	for (i = 0; i < 36; i++) d[i] /= rho;
}

// Better BT2 for ellipsoid
void Triangle6SommerBC::ellipsoidBT2(CoordSet &cs, double a, double b, double k, ComplexD *d) {
	double e[12];
	double f[12];
	double g[12];
	double tau1[12][3], tau2[12][3];

	int i, j, l, m;
// Basis for gradients
	double e1[3], e2[3], e3[3];
	double T[3][3];
	getLocalCoordinatesNew(cs, T, e1, e2, e3);
	double xl[6], yl[6], zl[6];
	getLocalCoordinates(cs, xl, yl, zl);

// Compute e, f, g, tau1, tau2 at each integration point
	for (i = 0; i < 12; i++) {

		double xli = 0.0, yli = 0.0;
		for (j = 0; j < 6; j++) {
			xli += tri6_values[i][j] * xl[j];
			yli += tri6_values[i][j] * yl[j];
		}
		double x = 0.0, y = 0.0, z = 0.0;
		Node nd1 = cs.getNode(nn[0]);

		x = nd1.x + xli * e1[0] + yli * e2[0];
		y = nd1.y + xli * e1[1] + yli * e2[1];
		z = nd1.z + xli * e1[2] + yli * e2[2];

//   fprintf(stderr,"%f\n",x*x/b/b+y*y/a/a+z*z/a/a);
//   fprintf(stderr,"%f\n",sqrt((x-nd1.x)*(x-nd1.x)+(y-nd1.y)*(y-nd1.y)+(z-nd1.z)*(z-nd1.z)));

		double sinasq = x * x / b / b;

		double k1 = b / a / sqrt(a * a * sinasq + b * b - b * b * sinasq);
		double k2 = a * b / pow(a * a * sinasq + b * b - b * b * sinasq, 1.5);
		e[i] = k1;
		f[i] = 0.0;
		g[i] = k2;
		if (fabs(fabs(x) - b) > 1e-4) {
			tau1[i][0] = 0.0;
			tau1[i][1] = -z;
			tau1[i][2] = y;
			double l = sqrt(y * y + z * z);
			tau1[i][0] = tau1[i][0] / l;
			tau1[i][1] = tau1[i][1] / l;
			tau1[i][2] = tau1[i][2] / l;
			tau2[i][0] = sqrt(b * b - x * x);
			tau2[i][1] = -x * y / sqrt(b * b - x * x);
			tau2[i][2] = -x * z / sqrt(b * b - x * x);
			l = sqrt(b * b - x * x + x * y * x * y / (b * b - x * x) + x * z * x * z / (b * b - x * x));
			tau2[i][0] = tau2[i][0] / l;
			tau2[i][1] = tau2[i][1] / l;
			tau2[i][2] = tau2[i][2] / l;
		} else {
			tau1[i][0] = 0.0;
			tau1[i][1] = 1.0;
			tau1[i][2] = 0.0;
			tau2[i][0] = 0.0;
			tau2[i][1] = 0.0;
			tau2[i][2] = 1.0;
		}
		/*  double w[3];
			   w[0] = tau1[i][1]*tau2[i][2]-tau2[i][1]*tau1[i][2];
			   w[1] = tau1[i][2]*tau2[i][0]-tau2[i][2]*tau1[i][0];
			   w[2] = tau1[i][0]*tau2[i][1]-tau2[i][0]*tau1[i][1];

		 fprintf(stderr,"!!!!  %f\n",w[0]*e3[0]+w[1]*e3[1]+w[2]*e3[2]); */
	}

// Compute gradients with respect to x,y
	double grad[6][12][2];
	double dxdxi[12][2][2], det[12];
	computedxdxi(xl, yl, dxdxi, det);
	for (i = 0; i < 12; i++) {
		Matrix22 &x = dxdxi[i];
		double inv[2][2];
		inv[0][0] = x[1][1];
		inv[0][1] = -x[0][1];
		inv[1][0] = -x[1][0];
		inv[1][1] = x[0][0];
// fprintf(stderr,"## %f %f %f %f\n",inv[0][0],inv[0][1],inv[1][0],inv[1][1]);

		int l, m, n;
//   fprintf(stderr,"+# %f\n",det[i]);
		for (l = 0; l < 6; l++)
			for (m = 0; m < 2; m++) {
				grad[l][i][m] = 0.0;
				for (n = 0; n < 2; n++) grad[l][i][m] += inv[n][m] * tri6_derivatives[i][l][n] / det[i];
//   fprintf(stderr,"@@ %f\n",grad[l][i][m]);
//   fprintf(stderr,"!! %f\n",tri6_derivatives[i][i][m]);
			}

	}


	double Pet[12][2][2];
	double Pte[12][2][2];
	for (m = 0; m < 12; m++) {
		Pet[m][0][0] = tau1[m][0] * e1[0] + tau1[m][1] * e1[1] + tau1[m][2] * e1[2];
		Pet[m][1][0] = tau2[m][0] * e1[0] + tau2[m][1] * e1[1] + tau2[m][2] * e1[2];
		Pet[m][0][1] = tau1[m][0] * e2[0] + tau1[m][1] * e2[1] + tau1[m][2] * e2[2];
		Pet[m][1][1] = tau2[m][0] * e2[0] + tau2[m][1] * e2[1] + tau2[m][2] * e2[2];

		Pte[m][0][0] = Pet[m][0][0];
		Pte[m][1][0] = Pet[m][0][1];
		Pte[m][0][1] = Pet[m][1][0];
		Pte[m][1][1] = Pet[m][1][1];
// fprintf(stderr,"## %f %f %f %f\n",Pet[m][0][0],Pet[m][0][1],Pet[m][1][0],Pet[m][1][1]);
	}


	ComplexD IRinv[12][2][2];
	for (m = 0; m < 12; m++) {
		ComplexD det = ComplexD(1.0, e[m] / k) * ComplexD(1.0, g[m] / k) -
		               ComplexD(0.0, f[m] / k) * ComplexD(0.0, f[m] / k);
		IRinv[m][0][0] = ComplexD(1.0, g[m] / k) / det;
		IRinv[m][1][1] = ComplexD(1.0, e[m] / k) / det;
		IRinv[m][0][1] = -ComplexD(0.0, f[m] / k) / det;
		IRinv[m][1][0] = -ComplexD(0.0, f[m] / k) / det;
// fprintf(stderr,"## %f %f %f %f\n",real(IRinv[m][0][0]),real(IRinv[m][0][1]),real(IRinv[m][1][0]),real(IRinv[m][1][1]));
// fprintf(stderr,"## %f %f %f %f\n",imag(IRinv[m][0][0]),imag(IRinv[m][0][1]),imag(IRinv[m][1][0]),imag(IRinv[m][1][1]));
	}

	ComplexD R[6][6];
	for (m = 0; m < 6; m++)
		for (j = 0; j < 6; j++) {
			R[m][j] = ComplexD(0.0, 0.0);
			for (l = 0; l < 12; l++) {
				// R^mj += b^jl*Pte^l*IRinv^l*RPet^l*b^ml
				double tmp1[2];
				ComplexD tmp2[2], tmp3[2];
				tmp1[0] = Pet[l][0][0] * grad[m][l][0] + Pet[l][0][1] * grad[m][l][1];
				tmp1[1] = Pet[l][1][0] * grad[m][l][0] + Pet[l][1][1] * grad[m][l][1];
				tmp2[0] = IRinv[l][0][0] * tmp1[0] + IRinv[l][0][1] * tmp1[1];
				tmp2[1] = IRinv[l][1][0] * tmp1[0] + IRinv[l][1][1] * tmp1[1];
				tmp3[0] = Pte[l][0][0] * tmp2[0] + Pte[l][0][1] * tmp2[1];
				tmp3[1] = Pte[l][1][0] * tmp2[0] + Pte[l][1][1] * tmp2[1];
				R[m][j] += 0.5 * tri6_weights[l] * (grad[j][l][0] * tmp3[0] + grad[j][l][1] * tmp3[1]) * det[l];

			}
		}

	for (m = 0; m < 36; m++) d[m] = R[m / 6][m % 6];
// for(m=0;m<36;m++)  fprintf(stderr,"%f %f\n",real(d[m]),imag(d[m]));
	if (el == 0) {
		fprintf(stderr,
		        "Triangle6SommerBC::ellipsoidBT2: adjacent element not defined.\n");
		exit(-1);
	}
	double rho = el->getProperty()->rho;
	for (i = 0; i < 36; i++) d[i] /= rho;
}
/*
int* Triangle6SommerBC::nodes(int *p) const
{
        if(p == 0) p = new int[6];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        p[3] = nn[3];
        p[4] = nn[4];
        p[5] = nn[5];
        return p;
}*/

/*
void Triangle6SommerBC::renum(int *table)
{
  int i;
  int *nn = getNodes();
                                                                                                                       
  for(i=0; i < numNodes(); ++i)
     nn[i] = table[nn[i]];
}*/

//functions that depend on the order of the absorbing boundary condition

int *
Triangle6SommerBC::dofs(DofSetArray &dsa, int *p) const  {
	if (p == 0) p = new int[6];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p + 1);
	dsa.number(nn[2], DofSet::Helm, p + 2);
	dsa.number(nn[3], DofSet::Helm, p + 3);
	dsa.number(nn[4], DofSet::Helm, p + 4);
	dsa.number(nn[5], DofSet::Helm, p + 5);

	return p;
}

int Triangle6SommerBC::numDofs() const {
	return 6;
}

void Triangle6SommerBC::markDofs(DofSetArray &dsa) const {
	dsa.mark(nn[0], DofSet::Helm);
	dsa.mark(nn[1], DofSet::Helm);
	dsa.mark(nn[2], DofSet::Helm);
	dsa.mark(nn[3], DofSet::Helm);
	dsa.mark(nn[4], DofSet::Helm);
	dsa.mark(nn[5], DofSet::Helm);
}
/*
FullSquareMatrix Triangle6SommerBC::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
fprintf(stderr,"massMatrix of a Triangle6SommerBC asked\n");
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
                                                                                                                       
        double x[6], y[6], z[6];
                                                                                                                       
        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
        x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;                                                                                                                
        const int numgauss = 2;

        for (int i =0; i <36 ;i++) d[i] = 0.0;
        FullSquareMatrix ret(6,d);

        return ret;
}

FullSquareMatrix Triangle6SommerBC::stiffness(const CoordSet &cs, double *Ks, int flg) const
{
fprintf(stderr,"stifness of a Traingle6SommerBC asked\n");
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);

        double x[6], y[6], z[6];
                                                                                                                       
        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
        x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;

        const int numgauss = 2;
        const int numdof   = 8;
        const int outerr   = 6;

        for (int i =0; i <36 ;i++) Ks[i] = 0.0;        
        FullSquareMatrix ret(6,Ks);

        return ret;
}*/

