#ifdef USE_EIGEN3
#ifndef _BEAMELEMENTTEMPLATE_CPP_
#define _BEAMELEMENTTEMPLATE_CPP_

#include <cmath>
#include <iostream>
#include <Element.d/Beam.d/BeamElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
doublereal
BeamElementTemplate<doublereal>
::getMassInTemplate(doublereal *x, doublereal *y, doublereal *z, doublereal rho, doublereal A) {
	using std::sqrt;

	doublereal dx = x[1] - x[0];
	doublereal dy = y[1] - y[0];
	doublereal dz = z[1] - z[0];

	doublereal length = sqrt(dx * dx + dy * dy + dz * dz);

	return length * rho * A;
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::gForce(doublereal *x, doublereal *y, doublereal *z, doublereal *gravityAcceleration,
         doublereal *_eframe, doublereal rho, doublereal A, doublereal *gravityForce, int gravflg) {
	using std::sqrt;
	doublereal massPerNode = 0.5 * getMassInTemplate(x, y, z, rho, A);

	Eigen::Map<Eigen::Matrix<doublereal, 3, 3> > eframe(_eframe);

	doublereal dx = x[1] - x[0];
	doublereal dy = y[1] - y[0];
	doublereal dz = z[1] - z[0];

	doublereal length = sqrt(dx * dx + dy * dy + dz * dz);
	doublereal t0n[3][3];

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			t0n[i][j] = eframe(j, i);
		}
	}

	int i;
	doublereal localg[3];

	for (i = 0; i < 3; ++i)
		localg[i] = 0.0;

	for (i = 0; i < 3; ++i) {
		localg[0] += t0n[0][i] * gravityAcceleration[i];
		localg[1] += t0n[1][i] * gravityAcceleration[i];
		localg[2] += t0n[2][i] * gravityAcceleration[i];
	}
	doublereal localf[3], localm[3];
	doublereal globalf[3], globalm[3];
	localf[0] = massPerNode * localg[0];
	localf[1] = massPerNode * localg[1];
	localf[2] = massPerNode * localg[2];
	if (gravflg == 2) { // consistent
		localm[0] = 0.0;
		localm[1] = -massPerNode * localg[2] * length / 6.0;
		localm[2] = massPerNode * localg[1] * length / 6.0;
	} else if (gravflg == 1) { // lumped with fixed-end moments
		localm[0] = 0.0;
		localm[1] = -massPerNode * localg[2] * length / 8.0;
		localm[2] = massPerNode * localg[1] * length / 8.0;
	} else {
		localm[0] = localm[1] = localm[2] = 0.0; // lumped without fixed-end moments
	}

	for (i = 0; i < 3; ++i) {
		globalf[i] = (t0n[0][i] * localf[0]) + (t0n[1][i] * localf[1]) + (t0n[2][i] * localf[2]);
		globalm[i] = (t0n[1][i] * localm[1]) + (t0n[2][i] * localm[2]);
	}

	gravityForce[0] = globalf[0];
	gravityForce[1] = globalf[1];
	gravityForce[2] = globalf[2];
	gravityForce[3] = globalm[0];
	gravityForce[4] = globalm[1];
	gravityForce[5] = globalm[2];
	gravityForce[6] = globalf[0];
	gravityForce[7] = globalf[1];
	gravityForce[8] = globalf[2];
	gravityForce[9] = -globalm[0];
	gravityForce[10] = -globalm[1];
	gravityForce[11] = -globalm[2];

}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::modmstif7(doublereal *_estif, doublereal A, doublereal E,
            doublereal *_eframe, doublereal Ix, doublereal Iy, doublereal Iz, doublereal alphay,
            doublereal alphaz, doublereal C1, doublereal nu, doublereal *_x, doublereal *_y, doublereal *_z,
            int flag) {
	using std::sqrt;
	using std::exp;

	Eigen::Map<Eigen::Matrix<doublereal, 9, 1> > eframe(_eframe);
	Eigen::Map<Eigen::Matrix<doublereal, 12, 12> > estif(_estif);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > x(_x), y(_y), z(_z);

//.....LOCAL VARIABLES

	int i, j, k, l;
	Eigen::Matrix<doublereal, 12, 12> Le;
	Eigen::Matrix<doublereal, 12, 12> locke;
	Eigen::Matrix<doublereal, 9, 1> T;
	doublereal zero, half, JJ, dx, dy, dz, length, G, one, two;
	doublereal c11, c22, c33, c12, c13, c23, eps;
	doublereal b, gammay, gammaz, twelve, six, xKe;
	doublereal four, bendy, bendz, bendCy, bendCz;
	doublereal s11, s22, s33, s26, s35;
	doublereal s44, s55, s66, s17, s28, s39, s59, s68;
	doublereal s77, s88, s99, s212, s311, s410, s511;
	doublereal s612, s812, s911, s1010, s1111, s1212;
	bool ortho, specialcasey, specialcasez;

//     ----
//     DATA
//     ----

	zero = 0.0;
	half = 0.5;
	one = 1.0;
	two = 2.0;
	four = 4.0;
	six = 6.0;
	twelve = 12.0;

//.....ZERO-MACHINE FOR FRAME ORTHONORMALITY CHECK

	eps = 1.0e-10;

//     -----
//     LOGIC
//     -----

//.....CLEAR THE OUTPUT STIFFNESS MATRIX
	estif.setZero();

//.....CLEAR THE LOCAL ARRAYS
	T.setZero();
	locke.setZero();
	Le.setZero();

//.....COMPUTE THE LENGTH OF THE BEAM ELEMENT

	dx = x[1] - x[0];
	dy = y[1] - y[0];
	dz = z[1] - z[0];
	length = sqrt(dx * dx + dy * dy + dz * dz);

	if (length == zero) error100();

//.....CHECK IF MATERIAL PROPERTIES ARE POSITIVE OR ZERO

	if (E <= zero) error200();
	if (A < zero) error200();
	if (nu < -1) error200();
	if (Ix < zero) error200();
	if (Iy < zero) error200();
	if (Iz < zero) error200();
	if (alphay < zero) error200();
	if (alphaz < zero) error200();
	if (C1 < zero) error200();

//.....EXTRACT THE ROTATION MATRIX FROM [EFRAME]
//.....       [ x_X y_X z_X ]
//..... [T] = [ x_Y y_Y z_Y ]
//.....       [ x_Z y_Z z_Z ]

	T = eframe;

//.....CHECK ORTHOGONALITY OF ROTATION MATRIX [T]

	c11 = T[0] * T[0] + T[1] * T[1] + T[2] * T[2];
	c22 = T[3] * T[3] + T[4] * T[4] + T[5] * T[5];
	c33 = T[6] * T[6] + T[7] * T[7] + T[8] * T[8];
	c12 = T[0] * T[3] + T[1] * T[4] + T[2] * T[5];
	c13 = T[0] * T[6] + T[1] * T[7] + T[2] * T[8];
	c23 = T[3] * T[6] + T[4] * T[7] + T[5] * T[8];

//.....ASSEMBLE THE TRANFORMATION MATRIX
//.....
//.....        [ [T]  0   0   0  ]
//.....        [  0  [T]  0   0  ]
//..... [Le] = [  0   0  [T]  0  ]
//.....        [  0   0   0  [T] ]

	for (k = 0; k < 4; ++k)
		for (j = 0; j < 3; ++j)
			for (i = 0; i < 3; ++i)
				Le(3 * k + i, 3 * k + j) = T[3 * j + i];

//.....INITIALIZE THE VARIABLE [JJ] FOR AN ARBITRARY CROSS SECTION

	JJ = Ix;

//.....INITIALIZE THE TRANSVERSE SHEAR MODULUS FOR ISOTROPIC MATERIAL

	G = E / (two * (one + nu));

//.....INITIALIZE THE ENTRIES OF THE LOCAL STIFFNESS MATRIX

	s11 = zero;
	s22 = zero;
	s33 = zero;
	s26 = zero;
	s35 = zero;
	s44 = zero;
	s55 = zero;
	s66 = zero;
	s17 = zero;
	s28 = zero;
	s39 = zero;
	s59 = zero;
	s68 = zero;
	s77 = zero;
	s88 = zero;
	s99 = zero;
	s212 = zero;
	s311 = zero;
	s410 = zero;
	s511 = zero;
	s612 = zero;
	s812 = zero;
	s911 = zero;
	s1010 = zero;
	s1111 = zero;
	s1212 = zero;

//.....INITIALIZE LOGICALS FOR PARTICULAR CASES

	specialcasey = false;
	specialcasez = false;

//.....INITIALIZE IF [ALPHA_Y] IS ZERO

	if (alphay == zero) {
		s33 = twelve * E * Iy / (length * length * length);
		s35 = -six * E * Iy / (length * length);
		s55 = four * E * Iy / length;
		s39 = -twelve * E * Iy / (length * length * length);
		s59 = six * E * Iy / (length * length);
		s99 = twelve * E * Iy / (length * length * length);
		s311 = -six * E * Iy / (length * length);
		s511 = two * E * Iy / length;
		s911 = six * E * Iy / (length * length);
		s1111 = four * E * Iy / length;
		specialcasey = true;
	}

//.....INITIALIZE IF [ALPHA_Z] IS ZERO

	if (alphaz == zero) {
		s22 = twelve * E * Iz / (length * length * length);
		s26 = six * E * Iz / (length * length);
		s66 = four * E * Iz / length;
		s28 = -twelve * E * Iz / (length * length * length);
		s68 = -six * E * Iz / (length * length);
		s88 = twelve * E * Iz / (length * length * length);
		s212 = six * E * Iz / (length * length);
		s612 = two * E * Iz / length;
		s812 = -six * E * Iz / (length * length);
		s1212 = four * E * Iz / length;
		specialcasez = true;
	}

//.....INITIALIZE IF [ALPHA_Y] IS NONZERO BUT [A] IS ZERO

	if ((alphay != zero) && (A == zero)) {
		s55 = E * Iy / length;
		s511 = -E * Iy / length;
		s1111 = E * Iy / length;
		specialcasey = true;
	}

//.....INITIALIZE IF [ALPHA_Z] IS NONZERO BUT [A] IS ZERO

	if ((alphaz != zero) && (A == zero)) {
		s66 = E * Iz / length;
		s612 = -E * Iz / length;
		s1212 = E * Iz / length;
		specialcasez = true;
	}

//.....INITIALIZE IF [I_Y] IS ZERO

	if (Iy == zero)
		specialcasey = true;

//.....INITIALIZE IF [I_Z] IS ZERO

	if (Iz == zero)
		specialcasez = true;

//.....INITIALIZE THE STIFFNESS ENTRIES IN ALL OTHER CASES

	s11 = E * A / length;
	s17 = -E * A / length;
	s77 = E * A / length;

	if (JJ == zero) {
		s44 = zero;
		s410 = zero;
		s1010 = zero;
	} else {
		if (C1 == zero) {
			s44 = G * JJ / length;
			s410 = -G * JJ / length;
			s1010 = G * JJ / length;
		} else {
			b = half * length * sqrt(G * JJ / C1);
			doublereal tanhb = (1. - exp(-2 * b)) / (1. + exp(-2 * b));
			s44 = (G * JJ) / (length * (one - (tanhb / b)));
			s410 = -(G * JJ) / (length * (one - (tanhb / b)));
			s1010 = (G * JJ) / (length * (one - (tanhb / b)));
		}
	}

	if (!specialcasey) {
		gammay = (G * A * length * length) / (twelve * E * Iy);
		bendy = (alphay + four * gammay) / (twelve * gammay * (alphay + gammay));
		bendCy = (-alphay + two * gammay) / (twelve * gammay * (alphay + gammay));
		s33 = G * A / (length * (alphay + gammay));
		s35 = -G * A / (two * (alphay + gammay));
		s55 = bendy * G * A * length;
		s39 = -G * A / (length * (alphay + gammay));
		s59 = G * A / (two * (alphay + gammay));
		s99 = G * A / (length * (alphay + gammay));
		s311 = -G * A / (two * (alphay + gammay));
		s511 = bendCy * G * A * length;
		s911 = G * A / (two * (alphay + gammay));
		s1111 = bendy * G * A * length;
	}

	if (!specialcasez) {
		gammaz = (G * A * length * length) / (twelve * E * Iz);
		bendz = (alphaz + four * gammaz) / (twelve * gammaz * (alphaz + gammaz));
		bendCz = (-alphaz + two * gammaz) / (twelve * gammaz * (alphaz + gammaz));
		s22 = G * A / (length * (alphaz + gammaz));
		s26 = G * A / (two * (alphaz + gammaz));
		s66 = bendz * G * A * length;
		s28 = -G * A / (length * (alphaz + gammaz));
		s68 = -G * A / (two * (alphaz + gammaz));
		s88 = G * A / (length * (alphaz + gammaz));
		s212 = G * A / (two * (alphaz + gammaz));
		s612 = bendCz * G * A * length;
		s812 = -G * A / (two * (alphaz + gammaz));
		s1212 = bendz * G * A * length;
	}

//.....INITIALIZE THE LOCAL STIFFNESS MATRIX

	locke(0, 0) = s11;
	locke(1, 1) = s22;
	locke(2, 2) = s33;
	locke(1, 5) = s26;
	locke(2, 4) = s35;
	locke(3, 3) = s44;
	locke(4, 4) = s55;
	locke(5, 5) = s66;
	locke(0, 6) = s17;
	locke(1, 7) = s28;
	locke(2, 8) = s39;
	locke(4, 8) = s59;
	locke(5, 7) = s68;
	locke(6, 6) = s77;
	locke(7, 7) = s88;
	locke(8, 8) = s99;
	locke(1, 11) = s212;
	locke(2, 10) = s311;
	locke(3, 9) = s410;
	locke(4, 10) = s511;
	locke(5, 11) = s612;
	locke(7, 11) = s812;
	locke(8, 10) = s911;
	locke(9, 9) = s1010;
	locke(10, 10) = s1111;
	locke(11, 11) = s1212;

	for (j = 0; j < 11; ++j)
		for (i = (j + 1); i < 12; ++i)
			locke(i, j) = locke(j, i);

//.....ASSEMBLE THE OUTPUT ELEMENTAL STIFFNESS MATRIX
//.....
//..... [ESTIF] = [Le] * [LOCKE] * [Le]^T
//.....

// if flag equals 1, return the transformed element stiffness matrix
// else, return the local element stiffness matrix

	if (flag == 1) {
		for (l = 0; l < 12; ++l) {
			for (k = 0; k < 12; ++k) {
				xKe = locke(k, l);
				for (j = 0; j < 12; ++j)
					for (i = 0; i < 12; ++i)
						estif(i, j) = estif(i, j) + Le(i, k) * xKe * Le(j, l);
			}
		}
	} else {
		for (i = 0; i < 12; ++i)
			for (j = 0; j < 12; ++j)
				estif(i, j) = locke(i, j);
	}
}


template<typename doublereal>
void
BeamElementTemplate<doublereal>
::error100() {
	using std::cerr;
	using std::endl;

	cerr << "*** FATAL ERROR in Routine MSTIF7 ***" << endl;
	cerr << "*** The Timoschenko Beam Element  ***" << endl;
	cerr << "*** Has Zero Length!              ***" << endl;
	cerr << "*** ... All Treatments Terminated ***" << endl;
	exit(-1);
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::error200() {
	using std::cerr;
	using std::endl;

	cerr << "*** FATAL ERROR in Routine MSTIF7      ***" << endl;
	cerr << "*** A Material or Geometrical Property ***" << endl;
	cerr << "*** is Not Positive: Check Input Data  ***" << endl;
	cerr << "*** ... All Treatments Terminated      ***" << endl;
	exit(-1);
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::error300() {
	using std::cerr;
	using std::endl;

	cerr << "*** FATAL ERROR in Routine MSTIF7      ***" << endl;
	cerr << "*** The Rotation Matrix Computed For a ***" << endl;
	cerr << "*** Timoshenko Beam is Not Orthogonal! ***" << endl;
	cerr << "*** ... All Treatments Terminated      ***" << endl;
	exit(-1);
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::error400() {
	using std::cerr;
	using std::endl;

	cerr << "*** FATAL ERROR in Routine MSTIF7    ***" << endl;
	cerr << "*** There is No Rotation Matrix      ***" << endl;
	cerr << "*** Available For a Timoshenko Beam  ***" << endl;
	cerr << "*** With Non-Symmetric Cross-section ***" << endl;
	cerr << "*** ... All Treatments Terminated    ***" << endl;
	exit(-1);
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::sands7(int elm, doublereal A, doublereal E,
         doublereal *_eframe, doublereal Ix, doublereal Iy, doublereal Iz,
         doublereal alphay, doublereal alphaz, doublereal C1, doublereal nu,
         doublereal *X, doublereal *Y, doublereal *Z, doublereal *ug,
         doublereal *_stress, int numel, int maxgus, int maxstr,
         int msize, doublereal alpha, doublereal tref, doublereal *_temp) {
	using std::sqrt;
	using std::exp;

	/* Initialized data */

	doublereal zero = 0.;
	doublereal half = .5;
	doublereal one = 1.;
	doublereal two = 2.;
	doublereal four = 4.;
	doublereal six = 6.;
	doublereal twelve = 12.;
	doublereal eps = 1e-10;

	/* Local variables */
	doublereal b, G;
	doublereal c11, c12, c22, c13, c33, c23,
			JJ, s11, s22, s33, s26, dx, dy, dz,
			s35, s44, s55, s66, s17, s28, s39, s59, s68, s77, s88, s99,
			s212, s311, s410, s511, s612, s812, s911, s1010, s1111, s1212;
	bool specialcasey, specialcasez;
	doublereal bendy, bendz;
	doublereal bendCy, bendCz, gammay, gammaz, length;
	doublereal tstress;

	Eigen::Matrix<doublereal, 12, 12> estif, locke, Le;
	Eigen::Matrix<doublereal, 12, 1> Fe, Ue;
	Eigen::Matrix<doublereal, 9, 1> T;
	Eigen::Matrix<doublereal, 2, 1> tl;
	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > stress(_stress, maxstr, maxgus);
	Eigen::Map<Eigen::Matrix<doublereal, 12, 1> > Ug(ug);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > temp(_temp);
	Eigen::Map<Eigen::Matrix<doublereal, 9, 1> > eframe(_eframe);

/* =====================================================================C */
/*                                                                     C */
/*     This Routine Retreives the Stresses for the 3D Timoshenko Beam  C */
/*     Element.                                                        C */
/*                                                                     C */
/*     Francois M. Hemez - July 12th 1994 - Version 1.0                C */
/*     transformed to C++ by Youngsoo Choi - Jan.2014                  C */
/*                                                                     C */
/* =====================================================================C */

/* .....CLEAR THE LOCAL STIFFNESS MATRIX */
	estif.setZero();

/* .....CLEAR THE LOCAL ARRAYS */
	T.setZero();
	locke.setZero();
	Le.setZero();
	Ue.setZero();
	Fe.setZero();


/* .....COMPUTE THE LENGTH OF THE BEAM ELEMENT */

	dx = X[1] - X[0];
	dy = Y[1] - Y[0];
	dz = Z[1] - Z[0];
	length = sqrt(dx * dx + dy * dy + dz * dz);

	if (length == zero) zeroLengthError();

/* .....CHECK IF MATERIAL PROPERTIES ARE POSITIVE OR ZERO */

	if (E <= zero || A < zero || nu < zero || Ix < zero || Iy < zero || Iz < zero || alphay < zero || alphaz < zero ||
	    C1 < zero) {
		wrongGeometricPropertyError();
	}

/* .....EXTRACT THE ROTATION MATRIX FROM [EFRAME] */
/* .....       [ x_X y_X z_X ] */
/* ..... [T] = [ x_Y y_Y z_Y ] */
/* ..... */

	for (int i = 0; i < 9; ++i) {
		T[i] = eframe[9 * (elm - 1) + i];
	}

/* .....CALCULATE THE ROTATION MATRIX IF THE INPUT ONE IS NILL */
/* .....ARCHTUNG!!! ASSUMES THAT THE BEAM ELEMENT IS SYMMETRIC */
/* .....(THAT IS: [Iy] = [Iz]) WITH CIRCULAR CROSS-SECTION SO */
/* .....THAT THE LOCAL [y] AND [z] AXES MAY BE ORIENTED */
/* .....ARBITRARILY IN THE PLANE NORMAL TO THE BEAM'S AXIS [x] */

/* .....CHECK ORTHOGONALITY OF ROTATION MATRIX [T] */

	c11 = T[0] * T[0] + T[1] * T[1] + T[2] * T[2];
	c22 = T[3] * T[3] + T[4] * T[4] + T[5] * T[5];
	c33 = T[6] * T[6] + T[7] * T[7] + T[8] * T[8];
	c12 = T[0] * T[3] + T[1] * T[4] + T[2] * T[5];
	c13 = T[0] * T[6] + T[1] * T[7] + T[2] * T[8];
	c23 = T[3] * T[6] + T[4] * T[7] + T[5] * T[8];

/* .....END OF TREATMENT FOR ROTATION MATRIX OF A SYMMETRIC BEAM */

/* .....ERROR-MESSAGE IF THE FRAME IS NOT AVAILABLE FOR A */
/* .....NON-SYMMETRIC TIMOSHENKO BEAM ELEMENT */

/* .....ASSEMBLE THE TRANFORMATION MATRIX */
/* ..... */
/* .....        [ [T]  0   0   0  ] */
/* .....        [  0  [T]  0   0  ] */
/* ..... [Le] = [  0   0  [T]  0  ] */
/* .....        [  0   0   0  [T] ] */

	for (int k = 0; k < 4; ++k) {
		for (int j = 0; j < 3; ++j) {
			for (int i = 0; i < 3; ++i) {
				Le(3 * k + i, 3 * k + j) = T[3 * j + i];
			}
		}
	}

/* .....INITIALIZE THE VARIABLE [JJ] FOR AN ARBITRARY CROSS SECTION */

	JJ = Ix;

/* .....INITIALIZE THE TRANSVERSE SHEAR MODULUS FOR ISOTROPIC MATERIAL */

	G = E / (two * (one + nu));

/* .....INITIALIZE THE ENTRIES OF THE LOCAL STIFFNESS MATRIX */

	s11 = zero;
	s17 = zero;
	s22 = zero;
	s26 = zero;
	s28 = zero;
	s212 = zero;
	s33 = zero;
	s35 = zero;
	s39 = zero;
	s311 = zero;
	s44 = zero;
	s410 = zero;
	s55 = zero;
	s59 = zero;
	s511 = zero;
	s66 = zero;
	s68 = zero;
	s612 = zero;
	s77 = zero;
	s88 = zero;
	s812 = zero;
	s99 = zero;
	s911 = zero;
	s1010 = zero;
	s1111 = zero;
	s1212 = zero;

/* .....INITIALIZE LOGICALS FOR PARTICULAR CASES */

	specialcasey = false;
	specialcasez = false;

/* .....INITIALIZE IF [ALPHA_Y] IS ZERO */

	if (alphay == zero) {
		s33 = twelve * E * Iy / (length * length * length);
		s35 = -six * E * Iy / (length * length);
		s55 = four * E * Iy / length;
		s39 = -twelve * E * Iy / (length * length * length);
		s59 = six * E * Iy / (length * length);
		s99 = twelve * E * Iy / (length * length * length);
		s311 = -six * E * Iy / (length * length);
		s511 = two * E * Iy / length;
		s911 = six * E * Iy / (length * length);
		s1111 = four * E * Iy / length;
		specialcasey = true;
	}

/* .....INITIALIZE IF [ALPHA_Z] IS ZERO */

	if (alphaz == zero) {
		s22 = twelve * E * Iz / (length * length * length);
		s26 = six * E * Iz / (length * length);
		s66 = four * E * Iz / length;
		s28 = -twelve * E * Iz / (length * length * length);
		s68 = -six * E * Iz / (length * length);
		s88 = twelve * E * Iz / (length * length * length);
		s212 = six * E * Iz / (length * length);
		s612 = two * E * Iz / length;
		s812 = -six * E * Iz / (length * length);
		s1212 = four * E * Iz / length;
		specialcasez = true;
	}

/* .....INITIALIZE IF [ALPHA_Y] IS NONZERO BUT [A] IS ZERO */

	if (alphay != zero && A == zero) {
		s55 = E * Iy / length;
		s511 = -(E) * Iy / length;
		s1111 = E * Iy / length;
		specialcasey = true;
	}

/* .....INITIALIZE IF [ALPHA_Z] IS NONZERO BUT [A] IS ZERO */

	if (alphaz != zero && A == zero) {
		s66 = E * Iz / length;
		s612 = -(E) * Iz / length;
		s1212 = E * Iz / length;
		specialcasez = true;
	}

/* .....INITIALIZE IF [I_Y] IS ZERO */

	if (Iy == zero) {
		specialcasey = true;
	}

/* .....INITIALIZE IF [I_Z] IS ZERO */

	if (Iz == zero) {
		specialcasez = true;
	}

/* .....INITIALIZE THE STIFFNESS ENTRIES IN ALL OTHER CASES */

	s11 = E * A / length;
	s17 = -(E) * A / length;
	s77 = E * A / length;

	if (JJ == zero) {
		s44 = zero;
		s410 = zero;
		s1010 = zero;
	} else {
		if (C1 == zero) {
			s44 = G * JJ / length;
			s410 = -G * JJ / length;
			s1010 = G * JJ / length;
		} else {
			b = half * length * sqrt(G * JJ / C1);
			doublereal tanhb = (1.0 - exp(-2.0 * b)) / (1.0 + exp(-2.0 * b));
			s44 = (G * JJ) / (length * (one - tanhb / b));
			s410 = -(G * JJ) / (length * (one - tanhb / b));
			s1010 = (G * JJ) / (length * (one - tanhb / b));
		}
	}

	if (!specialcasey) {
		gammay = G * A * length * length / (twelve * E * Iy);
		bendy = (alphay + four * gammay) / (twelve * gammay * (alphay + gammay));
		bendCy = (-(alphay) + two * gammay) / (twelve * gammay * (alphay + gammay));
		s33 = G * A / (length * (alphay + gammay));
		s35 = -G * A / (two * (alphay + gammay));
		s55 = bendy * G * A * length;
		s39 = -G * A / (length * (alphay + gammay));
		s59 = G * A / (two * (alphay + gammay));
		s99 = G * A / (length * (alphay + gammay));
		s311 = -G * A / (two * (alphay + gammay));
		s511 = bendCy * G * A * length;
		s911 = G * A / (two * (alphay + gammay));
		s1111 = bendy * G * A * length;
	}

	if (!specialcasez) {
		gammaz = G * A * length * length / (twelve * E * Iz);
		bendz = (alphaz + four * gammaz) / (twelve * gammaz * (alphaz + gammaz));
		bendCz = (-(alphaz) + two * gammaz) / (twelve * gammaz * (alphaz + gammaz));
		s22 = G * A / (length * (alphaz + gammaz));
		s26 = G * A / (two * (alphaz + gammaz));
		s66 = bendz * G * A * length;
		s28 = -G * A / (length * (alphaz + gammaz));
		s68 = -G * A / (two * (alphaz + gammaz));
		s88 = G * A / (length * (alphaz + gammaz));
		s212 = G * A / (two * (alphaz + gammaz));
		s612 = bendCz * G * A * length;
		s812 = -G * A / (two * (alphaz + gammaz));
		s1212 = bendz * G * A * length;
	}

/* .....INITIALIZE THE LOCAL STIFFNESS MATRIX */

	locke(0, 0) = s11;
	locke(1, 1) = s22;
	locke(2, 2) = s33;
	locke(1, 5) = s26;
	locke(2, 4) = s35;
	locke(3, 3) = s44;
	locke(4, 4) = s55;
	locke(5, 5) = s66;
	locke(0, 6) = s17;
	locke(1, 7) = s28;
	locke(2, 8) = s39;
	locke(4, 8) = s59;
	locke(5, 7) = s68;
	locke(6, 6) = s77;
	locke(7, 7) = s88;
	locke(8, 8) = s99;
	locke(1, 11) = s212;
	locke(2, 10) = s311;
	locke(3, 9) = s410;
	locke(4, 10) = s511;
	locke(5, 11) = s612;
	locke(7, 11) = s812;
	locke(8, 10) = s911;
	locke(9, 9) = s1010;
	locke(10, 10) = s1111;
	locke(11, 11) = s1212;

	for (int j = 0; j < 11; ++j) {
		for (int i = j; i < 12; ++i) {
			locke(i, j) = locke(j, i);
		}
	}

/* .....ASSEMBLE THE ELEMENTAL STIFFNESS MATRIX */
/* ..... */
/* ..... [ESTIF] = [Le] * [LOCKE] * [Le]^T */
/* ..... */


/* .....ROTATE THE GLOBAL DISPLACEMENT VECTOR BACK TO LOCAL FRAME */
	Ue = Le.transpose() * Ug;

/* .....COMPUTE THE INTERNAL FORCE RESULTANTS */
	Fe = locke * Ue;

/* .... COMPUTE THERMAL STRESS */
/* .... BE CAREFUL!!!! put tref to 0 in INPUT file if you don't want */
/* .... this term to be included */
	if (_temp != 0) {
		for (int i = 0; i < 2; ++i) tl[i] = temp[i] - tref;
		tstress = E * alpha * A * (tl[0] + tl[1]) / 2.0;
	} else tstress = 0;

/* .....WRITE THE RESULTANTS INTO THE STRESS ARRAY WHERE: */
/* .....FORCE_X  = STRESSXX */
/* .....FORCE_Y  = STRESSYY */
/* .....FORCE_Z  = STRESSZZ */
/* .....MOMENT_X = STRESSXY */
/* .....MOMENT_Y = STRESSXZ */
/* .....MOMENT_Z = STRESSYZ */

	stress(0, 0) = Fe[0] + tstress;
	stress(1, 0) = Fe[1];
	stress(2, 0) = Fe[2];
	stress(3, 0) = Fe[3];
	stress(4, 0) = Fe[4];
	stress(5, 0) = Fe[5];
	stress(0, 1) = Fe[6] - tstress;
	stress(1, 1) = Fe[7];
	stress(2, 1) = Fe[8];
	stress(3, 1) = Fe[9];
	stress(4, 1) = Fe[10];
	stress(5, 1) = Fe[11];

/* .....WRITE THE VON MISES RESULTANT INTO THE STRESS ARRAY */
	doublereal stress_0dxy = stress(0, 0) - stress(1, 0);
	doublereal stress_0dyz = stress(1, 0) - stress(2, 0);
	doublereal stress_0dxz = stress(0, 0) - stress(2, 0);
	doublereal stress_1dxy = stress(0, 1) - stress(1, 1);
	doublereal stress_1dyz = stress(1, 1) - stress(2, 1);
	doublereal stress_1dxz = stress(0, 1) - stress(2, 1);

	stress(6, 0) = sqrt(0.5 * (stress_0dxy * stress_0dxy + stress_0dyz * stress_0dyz + stress_0dxz * stress_0dxz) +
	                    3 * (stress(3, 0) * stress(3, 0) + stress(4, 0) * stress(4, 0) + stress(5, 0) * stress(5, 0)));
	stress(6, 1) = sqrt(0.5 * (stress_1dxy * stress_1dxy + stress_1dyz * stress_1dyz + stress_1dxz * stress_1dxz) +
	                    3 * (stress(3, 1) * stress(3, 1) + stress(4, 1) * stress(4, 1) + stress(5, 1) * stress(5, 1)));
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::transform(doublereal *_l, doublereal *_g, doublereal *_str) {
	using std::sqrt;

	Eigen::Map<Eigen::Matrix<doublereal, 9, 1> > l(_l);
	Eigen::Map<Eigen::Matrix<doublereal, 9, 1> > g(_g);

	Eigen::Matrix<doublereal, 3, 1> xl;
	xl << l.template block<3, 1>(0, 0);
	Eigen::Matrix<doublereal, 3, 1> yl;
	yl << l.template block<3, 1>(3, 0);
	Eigen::Matrix<doublereal, 3, 1> zl;
	zl << l.template block<3, 1>(6, 0);
	Eigen::Matrix<doublereal, 3, 1> xg;
	xg << g.template block<3, 1>(0, 0);
	Eigen::Matrix<doublereal, 3, 1> yg;
	yg << g.template block<3, 1>(3, 0);
	Eigen::Matrix<doublereal, 3, 1> zg;
	zg << g.template block<3, 1>(6, 0);
	Eigen::Map<Eigen::Matrix<doublereal, 6, 1> > str(_str);

	doublereal l1, l2, l3;
	doublereal m1, m2, m3;
	doublereal n1, n2, n3;
	Eigen::Matrix<doublereal, 6, 6> t;

	// Compute direction cosines
	l1 = xg.transpose() * xl;
	l2 = yg.transpose() * xl;
	l3 = zg.transpose() * xl;

	m1 = xg.transpose() * yl;
	m2 = yg.transpose() * yl;
	m3 = zg.transpose() * yl;

	n1 = xg.transpose() * zl;
	n2 = yg.transpose() * zl;
	n3 = zg.transpose() * zl;

	// Construct the 6x6 transformation matrix
	t(0, 0) = l1 * l1;
	t(0, 1) = m1 * m1;
	t(0, 2) = n1 * n1;
	t(0, 3) = 2.0 * l1 * m1;
	t(0, 4) = 2.0 * m1 * n1;
	t(0, 5) = 2.0 * n1 * l1;

	t(1, 0) = l2 * l2;
	t(1, 1) = m2 * m2;
	t(1, 2) = n2 * n2;
	t(1, 3) = 2.0 * l2 * m2;
	t(1, 4) = 2.0 * m2 * n2;
	t(1, 5) = 2.0 * n2 * l2;

	t(2, 0) = l3 * l3;
	t(2, 1) = m3 * m3;
	t(2, 2) = n3 * n3;
	t(2, 3) = 2.0 * l3 * m3;
	t(2, 4) = 2.0 * m3 * n3;
	t(2, 5) = 2.0 * n3 * l3;

	t(3, 0) = l1 * l2;
	t(3, 1) = m1 * m2;
	t(3, 2) = n1 * n2;
	t(3, 3) = l1 * m2 + l2 * m1;
	t(3, 4) = m1 * n2 + m2 * n1;
	t(3, 5) = n1 * l2 + n2 * l1;

	t(4, 0) = l2 * l3;
	t(4, 1) = m2 * m3;
	t(4, 2) = n2 * n3;
	t(4, 3) = l2 * m3 + l3 * m2;
	t(4, 4) = m2 * n3 + m3 * n2;
	t(4, 5) = n2 * l3 + n3 * l2;

	t(5, 0) = l3 * l1;
	t(5, 1) = m3 * m1;
	t(5, 2) = n3 * n1;
	t(5, 3) = l3 * m1 + l1 * m3;
	t(5, 4) = m3 * n1 + m1 * n3;
	t(5, 5) = n3 * l1 + n1 * l3;

	// Perform the multiplication
	str = t * str;

}

/* .....ERROR-MESSAGE IF THE BEAM ELEMENT HAS ZERO LENGTH */
template<typename doublereal>
void
BeamElementTemplate<doublereal>
::zeroLengthError() {
	std::cerr << "*** FATAL ERROR in Routine SANDS7 ***" << std::endl;
	std::cerr << "*** The Timoschenko Beam Element  ***" << std::endl;
	std::cerr << "*** Has Zero Length!              ***" << std::endl;
	std::cerr << "*** ... All Treatments Terminated ***" << std::endl;
	exit(-1);
}

/* .....ERROR-MESSAGE IF A MATERIAL/GEOMETRICAL PROPERTY IS WRONG */
template<typename doublereal>
void
BeamElementTemplate<doublereal>
::wrongGeometricPropertyError() {
	std::cerr << "*** FATAL ERROR in Routine SANDS7      ***" << std::endl;
	std::cerr << "*** A Material or Geometrical Property ***" << std::endl;
	std::cerr << "*** is Not Positive: Check Input Data  ***" << std::endl;
	std::cerr << "*** ... All Treatments Terminated      ***" << std::endl;
	exit(-1);
}

/* .....ERROR-MESSAGE IF THE ROTATION MATRIX IS NOT ORTHOGONAL */
template<typename doublereal>
void
BeamElementTemplate<doublereal>
::rotationMatrixError() {
	std::cerr << "*** FATAL ERROR in Routine SANDS7      ***" << std::endl;
	std::cerr << "*** The Rotation Matrix Computed For a ***" << std::endl;
	std::cerr << "*** Timoshenko Beam is Not Orthogonal! ***" << std::endl;
	std::cerr << "*** ... All Treatments Terminated      ***" << std::endl;
	exit(-1);
}

/* .....ERROR-MESSAGE IF THE ROTATION MATRIX IS NOT AVAILABLE */
/* .....FOR A BEAM WITH NON-SYMMETRIC CROSS-SECTIONAL AREA */
template<typename doublereal>
void
BeamElementTemplate<doublereal>
::missingRotationMatrix() {
	std::cerr << "*** FATAL ERROR in Routine SANDS7    ***" << std::endl;
	std::cerr << "*** There is No Rotation Matrix      ***" << std::endl;
	std::cerr << "*** Available For a Timoshenko Beam  ***" << std::endl;
	std::cerr << "*** With Non-Symmetric Cross-section ***" << std::endl;
	std::cerr << "*** ... All Treatments Terminated    ***" << std::endl;
	exit(-1);
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::vms7WRTdisp(int elm, doublereal A, doublereal E,
              doublereal *_eframe, doublereal Ix, doublereal Iy, doublereal Iz,
              doublereal alphay, doublereal alphaz, doublereal C1, doublereal nu,
              doublereal *X, doublereal *Y, doublereal *Z, doublereal *ug,
              doublereal *_vmsWRTdisp,
              doublereal alpha, doublereal tref, doublereal *_temp) {
	using std::sqrt;
	using std::exp;

	/* Initialized data */

	const doublereal zero = 0.;
	const doublereal half = .5;
	const doublereal one = 1.;
	const doublereal two = 2.;
	const doublereal four = 4.;
	const doublereal six = 6.;
	const doublereal twelve = 12.;
	const doublereal eps = 1e-10;

	/* Local variables */
	doublereal b, G;
	doublereal c11, c12, c22, c13, c33, c23,
			JJ, s11, s22, s33, s26, dx, dy, dz,
			s35, s44, s55, s66, s17, s28, s39, s59, s68, s77, s88, s99,
			s212, s311, s410, s511, s612, s812, s911, s1010, s1111, s1212;
	bool specialcasey, specialcasez;
	doublereal bendy, bendz;
	doublereal bendCy, bendCz, gammay, gammaz, length;
	doublereal tstress;

	Eigen::Matrix<doublereal, 12, 12> estif, locke, Le;
	Eigen::Matrix<doublereal, 12, 1> Fe, Ue;
	Eigen::Matrix<doublereal, 9, 1> T;
	Eigen::Matrix<doublereal, 2, 1> tl, vms;
	Eigen::Matrix<doublereal, 7, 2> stress;
	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > dvmsWRTdisp(_vmsWRTdisp, 2, 12);
	Eigen::Map<Eigen::Matrix<doublereal, 12, 1> > Ug(ug);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > temp(_temp);
	Eigen::Map<Eigen::Matrix<doublereal, 9, 1> > eframe(_eframe);

/* =====================================================================C */
/*                                                                     C */
/*     This Routine Retreives the sensitivity of Stresses wrt          C */
/*     displacement for the 3D Timoshenko Beam Element.                C */
/*                                                                     C */
/*     Youngsoo Choi - January 2014                                    C */
/*                                                                     C */
/* =====================================================================C */

/* .....CLEAR THE LOCAL STIFFNESS MATRIX */
	estif.setZero();

/* .....CLEAR THE LOCAL ARRAYS */
	vms.setZero();
	T.setZero();
	locke.setZero();
	Le.setZero();
	Ue.setZero();
	Fe.setZero();


/* .....COMPUTE THE LENGTH OF THE BEAM ELEMENT */

	dx = X[1] - X[0];
	dy = Y[1] - Y[0];
	dz = Z[1] - Z[0];
	length = sqrt(dx * dx + dy * dy + dz * dz);

	if (length == zero) zeroLengthError();

/* .....CHECK IF MATERIAL PROPERTIES ARE POSITIVE OR ZERO */

	if (E <= zero || A < zero || nu < zero || Ix < zero || Iy < zero || Iz < zero || alphay < zero || alphaz < zero ||
	    C1 < zero) {
		wrongGeometricPropertyError();
	}

/* .....EXTRACT THE ROTATION MATRIX FROM [EFRAME] */
/* .....       [ x_X y_X z_X ] */
/* ..... [T] = [ x_Y y_Y z_Y ] */
/* ..... */

	for (int i = 0; i < 9; ++i) {
		T[i] = eframe[9 * (elm - 1) + i];
	}

/* .....CALCULATE THE ROTATION MATRIX IF THE INPUT ONE IS NILL */
/* .....ARCHTUNG!!! ASSUMES THAT THE BEAM ELEMENT IS SYMMETRIC */
/* .....(THAT IS: [Iy] = [Iz]) WITH CIRCULAR CROSS-SECTION SO */
/* .....THAT THE LOCAL [y] AND [z] AXES MAY BE ORIENTED */
/* .....ARBITRARILY IN THE PLANE NORMAL TO THE BEAM'S AXIS [x] */

/* .....CHECK ORTHOGONALITY OF ROTATION MATRIX [T] */

	c11 = T[0] * T[0] + T[1] * T[1] + T[2] * T[2];
	c22 = T[3] * T[3] + T[4] * T[4] + T[5] * T[5];
	c33 = T[6] * T[6] + T[7] * T[7] + T[8] * T[8];
	c12 = T[0] * T[3] + T[1] * T[4] + T[2] * T[5];
	c13 = T[0] * T[6] + T[1] * T[7] + T[2] * T[8];
	c23 = T[3] * T[6] + T[4] * T[7] + T[5] * T[8];

/* .....END OF TREATMENT FOR ROTATION MATRIX OF A SYMMETRIC BEAM */

/* .....ERROR-MESSAGE IF THE FRAME IS NOT AVAILABLE FOR A */
/* .....NON-SYMMETRIC TIMOSHENKO BEAM ELEMENT */

/* .....ASSEMBLE THE TRANFORMATION MATRIX */
/* ..... */
/* .....        [ [T]  0   0   0  ] */
/* .....        [  0  [T]  0   0  ] */
/* ..... [Le] = [  0   0  [T]  0  ] */
/* .....        [  0   0   0  [T] ] */

	for (int k = 0; k < 4; ++k) {
		for (int j = 0; j < 3; ++j) {
			for (int i = 0; i < 3; ++i) {
				Le(3 * k + i, 3 * k + j) = T[3 * j + i];
			}
		}
	}

/* .....INITIALIZE THE VARIABLE [JJ] FOR AN ARBITRARY CROSS SECTION */

	JJ = Ix;

/* .....INITIALIZE THE TRANSVERSE SHEAR MODULUS FOR ISOTROPIC MATERIAL */

	G = E / (two * (one + nu));

/* .....INITIALIZE THE ENTRIES OF THE LOCAL STIFFNESS MATRIX */

	s11 = zero;
	s17 = zero;
	s22 = zero;
	s26 = zero;
	s28 = zero;
	s212 = zero;
	s33 = zero;
	s35 = zero;
	s39 = zero;
	s311 = zero;
	s44 = zero;
	s410 = zero;
	s55 = zero;
	s59 = zero;
	s511 = zero;
	s66 = zero;
	s68 = zero;
	s612 = zero;
	s77 = zero;
	s88 = zero;
	s812 = zero;
	s99 = zero;
	s911 = zero;
	s1010 = zero;
	s1111 = zero;
	s1212 = zero;

/* .....INITIALIZE LOGICALS FOR PARTICULAR CASES */

	specialcasey = false;
	specialcasez = false;

/* .....INITIALIZE IF [ALPHA_Y] IS ZERO */

	if (alphay == zero) {
		s33 = twelve * E * Iy / (length * length * length);
		s35 = -six * E * Iy / (length * length);
		s55 = four * E * Iy / length;
		s39 = -twelve * E * Iy / (length * length * length);
		s59 = six * E * Iy / (length * length);
		s99 = twelve * E * Iy / (length * length * length);
		s311 = -six * E * Iy / (length * length);
		s511 = two * E * Iy / length;
		s911 = six * E * Iy / (length * length);
		s1111 = four * E * Iy / length;
		specialcasey = true;
	}

/* .....INITIALIZE IF [ALPHA_Z] IS ZERO */

	if (alphaz == zero) {
		s22 = twelve * E * Iz / (length * length * length);
		s26 = six * E * Iz / (length * length);
		s66 = four * E * Iz / length;
		s28 = -twelve * E * Iz / (length * length * length);
		s68 = -six * E * Iz / (length * length);
		s88 = twelve * E * Iz / (length * length * length);
		s212 = six * E * Iz / (length * length);
		s612 = two * E * Iz / length;
		s812 = -six * E * Iz / (length * length);
		s1212 = four * E * Iz / length;
		specialcasez = true;
	}

/* .....INITIALIZE IF [ALPHA_Y] IS NONZERO BUT [A] IS ZERO */

	if (alphay != zero && A == zero) {
		s55 = E * Iy / length;
		s511 = -(E) * Iy / length;
		s1111 = E * Iy / length;
		specialcasey = true;
	}

/* .....INITIALIZE IF [ALPHA_Z] IS NONZERO BUT [A] IS ZERO */

	if (alphaz != zero && A == zero) {
		s66 = E * Iz / length;
		s612 = -(E) * Iz / length;
		s1212 = E * Iz / length;
		specialcasez = true;
	}

/* .....INITIALIZE IF [I_Y] IS ZERO */

	if (Iy == zero) {
		specialcasey = true;
	}

/* .....INITIALIZE IF [I_Z] IS ZERO */

	if (Iz == zero) {
		specialcasez = true;
	}

/* .....INITIALIZE THE STIFFNESS ENTRIES IN ALL OTHER CASES */

	s11 = E * A / length;
	s17 = -(E) * A / length;
	s77 = E * A / length;

	if (JJ == zero) {
		s44 = zero;
		s410 = zero;
		s1010 = zero;
	} else {
		if (C1 == zero) {
			s44 = G * JJ / length;
			s410 = -G * JJ / length;
			s1010 = G * JJ / length;
		} else {
			b = half * length * sqrt(G * JJ / C1);
			doublereal tanhb = (1.0 - exp(-2.0 * b)) / (1.0 + exp(-2.0 * b));
			s44 = (G * JJ) / (length * (one - tanhb / b));
			s410 = -(G * JJ) / (length * (one - tanhb / b));
			s1010 = (G * JJ) / (length * (one - tanhb / b));
		}
	}

	if (!specialcasey) {
		gammay = G * A * length * length / (twelve * E * Iy);
		bendy = (alphay + four * gammay) / (twelve * gammay * (alphay + gammay));
		bendCy = (-(alphay) + two * gammay) / (twelve * gammay * (alphay + gammay));
		s33 = G * A / (length * (alphay + gammay));
		s35 = -G * A / (two * (alphay + gammay));
		s55 = bendy * G * A * length;
		s39 = -G * A / (length * (alphay + gammay));
		s59 = G * A / (two * (alphay + gammay));
		s99 = G * A / (length * (alphay + gammay));
		s311 = -G * A / (two * (alphay + gammay));
		s511 = bendCy * G * A * length;
		s911 = G * A / (two * (alphay + gammay));
		s1111 = bendy * G * A * length;
	}

	if (!specialcasez) {
		gammaz = G * A * length * length / (twelve * E * Iz);
		bendz = (alphaz + four * gammaz) / (twelve * gammaz * (alphaz + gammaz));
		bendCz = (-(alphaz) + two * gammaz) / (twelve * gammaz * (alphaz + gammaz));
		s22 = G * A / (length * (alphaz + gammaz));
		s26 = G * A / (two * (alphaz + gammaz));
		s66 = bendz * G * A * length;
		s28 = -G * A / (length * (alphaz + gammaz));
		s68 = -G * A / (two * (alphaz + gammaz));
		s88 = G * A / (length * (alphaz + gammaz));
		s212 = G * A / (two * (alphaz + gammaz));
		s612 = bendCz * G * A * length;
		s812 = -G * A / (two * (alphaz + gammaz));
		s1212 = bendz * G * A * length;
	}

/* .....INITIALIZE THE LOCAL STIFFNESS MATRIX */

	locke(0, 0) = s11;
	locke(1, 1) = s22;
	locke(2, 2) = s33;
	locke(1, 5) = s26;
	locke(2, 4) = s35;
	locke(3, 3) = s44;
	locke(4, 4) = s55;
	locke(5, 5) = s66;
	locke(0, 6) = s17;
	locke(1, 7) = s28;
	locke(2, 8) = s39;
	locke(4, 8) = s59;
	locke(5, 7) = s68;
	locke(6, 6) = s77;
	locke(7, 7) = s88;
	locke(8, 8) = s99;
	locke(1, 11) = s212;
	locke(2, 10) = s311;
	locke(3, 9) = s410;
	locke(4, 10) = s511;
	locke(5, 11) = s612;
	locke(7, 11) = s812;
	locke(8, 10) = s911;
	locke(9, 9) = s1010;
	locke(10, 10) = s1111;
	locke(11, 11) = s1212;

	for (int j = 0; j < 11; ++j) {
		for (int i = j; i < 12; ++i) {
			locke(i, j) = locke(j, i);
		}
	}

/* .....ASSEMBLE THE ELEMENTAL STIFFNESS MATRIX */
/* ..... */
/* ..... [ESTIF] = [Le] * [LOCKE] * [Le]^T */
/* ..... */

/* .....ROTATE THE GLOBAL DISPLACEMENT VECTOR BACK TO LOCAL FRAME */
	Ue = Le.transpose() * Ug;

/* .....COMPUTE THE INTERNAL FORCE RESULTANTS */
	Fe = locke * Ue;


/* .... COMPUTE THERMAL STRESS */
/* .... BE CAREFUL!!!! put tref to 0 in INPUT file if you don't want */
/* .... this term to be included */
	if (_temp != 0) {
		for (int i = 0; i < 2; ++i) tl[i] = temp[i] - tref;
		tstress = E * alpha * A * (tl[0] + tl[1]) / 2.0;
	} else tstress = 0;

/* .....WRITE THE RESULTANTS INTO THE STRESS ARRAY WHERE: */
/* .....FORCE_X  = STRESSXX */
/* .....FORCE_Y  = STRESSYY */
/* .....FORCE_Z  = STRESSZZ */
/* .....MOMENT_X = STRESSXY */
/* .....MOMENT_Y = STRESSXZ */
/* .....MOMENT_Z = STRESSYZ */

	stress(0, 0) = Fe[0] + tstress;
	stress(1, 0) = Fe[1];
	stress(2, 0) = Fe[2];
	stress(3, 0) = Fe[3];
	stress(4, 0) = Fe[4];
	stress(5, 0) = Fe[5];
	stress(0, 1) = Fe[6] - tstress;
	stress(1, 1) = Fe[7];
	stress(2, 1) = Fe[8];
	stress(3, 1) = Fe[9];
	stress(4, 1) = Fe[10];
	stress(5, 1) = Fe[11];

/* .....WRITE THE VON MISES RESULTANT INTO THE STRESS ARRAY */
	doublereal stress_0dxy = stress(0, 0) - stress(1, 0);
	doublereal stress_0dyz = stress(1, 0) - stress(2, 0);
	doublereal stress_0dxz = stress(0, 0) - stress(2, 0);
	doublereal stress_1dxy = stress(0, 1) - stress(1, 1);
	doublereal stress_1dyz = stress(1, 1) - stress(2, 1);
	doublereal stress_1dxz = stress(0, 1) - stress(2, 1);

	stress(6, 0) = sqrt(0.5 * (stress_0dxy * stress_0dxy + stress_0dyz * stress_0dyz + stress_0dxz * stress_0dxz) +
	                    3 * (stress(3, 0) * stress(3, 0) + stress(4, 0) * stress(4, 0) + stress(5, 0) * stress(5, 0))) +
	               1.0e-10;
	stress(6, 1) = sqrt(0.5 * (stress_1dxy * stress_1dxy + stress_1dyz * stress_1dyz + stress_1dxz * stress_1dxz) +
	                    3 * (stress(3, 1) * stress(3, 1) + stress(4, 1) * stress(4, 1) + stress(5, 1) * stress(5, 1))) +
	               1.0e-10;

	Eigen::Matrix<doublereal, 2, 12> dvmsdStress;
	dvmsdStress.setZero();
	dvmsdStress(0, 0) = (2. * stress(0, 0) - stress(1, 0) - stress(2, 0)) / (2. * stress(6, 0));
	dvmsdStress(0, 1) = (2. * stress(1, 0) - stress(0, 0) - stress(2, 0)) / (2. * stress(6, 0));
	dvmsdStress(0, 2) = (2. * stress(2, 0) - stress(1, 0) - stress(0, 0)) / (2. * stress(6, 0));
	dvmsdStress(0, 3) = (3. * stress(3, 0)) / stress(6, 0);
	dvmsdStress(0, 4) = (3. * stress(4, 0)) / stress(6, 0);
	dvmsdStress(0, 5) = (3. * stress(5, 0)) / stress(6, 0);

	dvmsdStress(1, 6) = (2. * stress(0, 1) - stress(1, 1) - stress(2, 1)) / (2. * stress(6, 1));
	dvmsdStress(1, 7) = (2. * stress(1, 1) - stress(0, 1) - stress(2, 1)) / (2. * stress(6, 1));
	dvmsdStress(1, 8) = (2. * stress(2, 1) - stress(1, 1) - stress(0, 1)) / (2. * stress(6, 1));
	dvmsdStress(1, 9) = (3. * stress(3, 1)) / stress(6, 1);
	dvmsdStress(1, 10) = (3. * stress(4, 1)) / stress(6, 1);
	dvmsdStress(1, 11) = (3. * stress(5, 1)) / stress(6, 1);

	dvmsWRTdisp = dvmsdStress * locke * Le.transpose();

}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::frame6(doublereal *frame,
         doublereal *u,
         doublereal *v,
         doublereal *w) {
	u[0] = frame[0];
	u[1] = frame[1];
	u[2] = frame[2];
	v[0] = frame[3];
	v[1] = frame[4];
	v[2] = frame[5];
	w[0] = frame[6];
	w[1] = frame[7];
	w[2] = frame[8];
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::sands6(doublereal area, doublereal e, int elm,
         doublereal *_stress, int maxsze, int maxgus, int maxstr,
         doublereal *_eframe,
         doublereal ix, doublereal iy, doublereal iz,
         doublereal nu, doublereal *_x, doublereal *_y, doublereal *_z,
         doublereal *_ug, doublereal alpha, doublereal tref, doublereal *_temp) {
/*********************************************************************
*	THIS SUBROUTINE WILL COMPUTE THE INTERNAL FORCES FOR THE    *
* EULER-BERNOULLI BEAM ELEMENT AND STORE THEM IN THE STRESS ARRAY   *
*                                                                   *
*********************************************************************
*
*	AUTHOR  :  P.R. STERN
*	DATE    :  SEPTEMBER 1991
*	VERSION :  3.0 Multielement Revison
*
* Tranformed to C++ by Youngsoo Choi, January 2014
*
*********************************************************************
*
*		DEFINE THE GLOBAL VARIABLES
*
*	  AREA = CROSS-SECTIONAL AREA OF THE BEAM
*	     E = YOUNGS MODULUS FORTHE BEAM
*        ALPHA = DILATATION COEFFICIENT
*	   ELM = CURRENT ELEMENT NUMBER
*	STRESS = INTERNAL FORCE ARRAY
*       MAXSZE = LEADING DIMENSION OF STRESS AND STRAIN ARRAYS
*	MAXGUS = SECOND DIMENSION OF STRESS
*	MAXSTR = THIRD DIMENSION OF STRESS
*	EFRAME = ELEMENT REFERENCE FRAMES
*     IX,IY,IZ = BEAM MOMENTS OF INERTIA
*	    NU = POISSON'S RATIO
*        X,Y,Z = COORDINATES FOR THE BEAM
*           UG = GLOBAL DISPLACEMENT VECTOR FOR ELEMENT #ELM
*         TEMP = NODAL TEMPERATURE
*         TREF = REFERENCE TEMPERATURE
*      TSTRESS = THERMAL STRESS
*
*********************************************************************
*
*		CALLED BY : MDERIV
*
*********************************************************************/

	using std::sqrt;
	doublereal eiy, eiz, length2, length3, tstress;

//.... LOCAL VARIABLES

	int ii, jj, kk, ic;

	doublereal pi, G, J, dx, dy, dz, length;
	Eigen::Matrix<doublereal, 10, 1> vec;
	Eigen::Matrix<doublereal, 12, 12> ke, tran;
	Eigen::Matrix<doublereal, 12, 1> tug, res;
	Eigen::Matrix<doublereal, 3, 1> u, v, w;
	Eigen::Matrix<doublereal, 3, 3> t33;
	Eigen::Matrix<doublereal, 2, 1> tl;

	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > stress(_stress, maxstr, maxgus);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > x(_x, 2, 1);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > y(_y, 2, 1);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > z(_z, 2, 1);
	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > eframe(_eframe, 9, 1);
	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > ug(_ug, 12, 1);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > temp(_temp, 2, 1);



//.... INITIALIZE LOCAL CONSTANTS

	pi = 4.0 * atan(1.0);
	G = e / (2.0 * (1.0 + nu));
	J = ix;

//.... COMPUTE THE LENGTH OF THE BEAM

	dx = x[1] - x[0];
	dy = y[1] - y[0];
	dz = z[1] - z[0];
	length = sqrt(dx * dx + dy * dy + dz * dz);

//.... COMPUTE THE CONSTITUTIVE RELATIONS FOR THE BEAM

	eiy = e * iy;
	eiz = e * iz;
	length2 = length * length;
	length3 = length2 * length;
	vec[0] = e * area / length;
	vec[1] = 12.0 * eiy / length3;
	vec[2] = 12.0 * eiz / length3;
	vec[3] = 6.0 * eiy / length2;
	vec[4] = 6.0 * eiz / length2;
	vec[5] = G * J / length;
	vec[6] = 4.0 * eiy / length;
	vec[7] = 4.0 * eiz / length;
	vec[8] = 2.0 * eiy / length;
	vec[9] = 2.0 * eiz / length;

//.... INITIALIZE MATRICES AND VECTORS

	tug.setZero();
	res.setZero();
	ke.setZero();
	tran.setZero();

//.... BUILD THE ELEMENT STIFFNESS MATRIX

	ke(0, 0) = vec[0];
	ke(1, 1) = vec[2];
	ke(2, 2) = vec[1];
	ke(3, 3) = vec[5];
	ke(4, 4) = vec[6];
	ke(5, 5) = vec[7];
	ke(6, 6) = vec[0];
	ke(7, 7) = vec[2];
	ke(8, 8) = vec[1];
	ke(9, 9) = vec[5];
	ke(10, 10) = vec[6];
	ke(11, 11) = vec[7];

	ke(4, 2) = -vec[3];
	ke(2, 4) = ke(4, 2);
	ke(5, 1) = vec[4];
	ke(1, 5) = ke(5, 1);
	ke(6, 0) = -vec[0];
	ke(0, 6) = ke(6, 0);
	ke(7, 1) = -vec[2];
	ke(1, 7) = ke(7, 1);
	ke(7, 5) = -vec[4];
	ke(5, 7) = ke(7, 5);
	ke(8, 2) = -vec[1];
	ke(2, 8) = ke(8, 2);
	ke(8, 4) = vec[3];
	ke(4, 8) = ke(8, 4);
	ke(9, 3) = -vec[5];
	ke(3, 9) = ke(9, 3);
	ke(10, 2) = -vec[3];
	ke(2, 10) = ke(10, 2);
	ke(10, 4) = vec[8];
	ke(4, 10) = ke(10, 4);
	ke(10, 8) = vec[3];
	ke(8, 10) = ke(10, 8);
	ke(11, 1) = vec[4];
	ke(1, 11) = ke(11, 1);
	ke(11, 5) = vec[9];
	ke(5, 11) = ke(11, 5);
	ke(11, 7) = -vec[4];
	ke(7, 11) = ke(11, 7);

//        for(ii = 0; ii < 12; ++ii)
//          for(jj = 0; jj < 12; ++jj)
//            ke(jj,ii) = ke(ii,jj);

//.... COMPUTE THE TRANFORMATION MATRIX

	frame6(eframe.data(), u.data(), v.data(), w.data());

	t33 << u.transpose(), v.transpose(), w.transpose();

	for (kk = 0; kk < 4; ++kk) {
		ic = 3 * kk;
		for (ii = 0; ii < 3; ++ii)
			for (jj = 0; jj < 3; ++jj)
				tran(ic + ii, ic + jj) = t33(ii, jj);
	}

//.... ROTATE THE GLOBAL DISPLACEMENT VECTOR BACK TO LOCAL FRAME

	tug = tran * ug;

//.... COMPUTE THE INTERNAL FORCE RESULTANTS

	res = ke * tug;

//.... COMPUTE THERMAL STRESS
//.... BE CAREFUL!!!!! put tref to 0 in INPUT file if you don't want
//.... this term to be included

	if (_temp) {
		for (int i = 0; i < 2; ++i)
			tl[i] = temp[i] - tref;
	} else tl.setZero();

	tstress = e * alpha * area * (tl[0] + tl[1]) / 2.;

//.... WRITE THE RESULTANTS INTO THE STRESS ARRAY WHERE

//        FORX = STRESSXX      FORY=STRESSYY     FORZ=STRESSZZ
//        MOMX = STRESSXY      MOMY=STRESSXZ     MOMZ=STRESSYZ

	stress(0, 0) = res[0] + tstress;
	stress(1, 0) = res[1];
	stress(2, 0) = res[2];
	stress(3, 0) = res[3];
	stress(4, 0) = res[4];
	stress(5, 0) = res[5];
	stress(0, 1) = res[6] - tstress;
	stress(1, 1) = res[7];
	stress(2, 1) = res[8];
	stress(3, 1) = res[9];
	stress(4, 1) = res[10];
	stress(5, 1) = res[11];

//..... WRITE THE VON MISES RESULTANT INTO THE STRESS ARRAY
	doublereal stress_0dxy = stress(0, 0) - stress(1, 0);
	doublereal stress_0dyz = stress(1, 0) - stress(2, 0);
	doublereal stress_0dxz = stress(0, 0) - stress(2, 0);
	doublereal stress_1dxy = stress(0, 1) - stress(1, 1);
	doublereal stress_1dyz = stress(1, 1) - stress(2, 1);
	doublereal stress_1dxz = stress(0, 1) - stress(2, 1);

	stress(6, 0) = sqrt(0.5 * (stress_0dxy * stress_0dxy + stress_0dyz * stress_0dyz + stress_0dxz * stress_0dxz) +
	                    3 * (stress(3, 0) * stress(3, 0) + stress(4, 0) * stress(4, 0) + stress(5, 0) * stress(5, 0)));
	stress(6, 1) = sqrt(0.5 * (stress_1dxy * stress_1dxy + stress_1dyz * stress_1dyz + stress_1dxz * stress_1dxz) +
	                    3 * (stress(3, 1) * stress(3, 1) + stress(4, 1) * stress(4, 1) + stress(5, 1) * stress(5, 1)));
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::vms6WRTdisp(doublereal area, doublereal e, int elm,
              doublereal *_vmsWRTdisp, int maxsze, int maxgus, int maxstr,
              doublereal *_eframe,
              doublereal ix, doublereal iy, doublereal iz,
              doublereal nu, doublereal *_x, doublereal *_y, doublereal *_z,
              doublereal *_ug, doublereal alpha, doublereal tref, doublereal *_temp) {
/*********************************************************************
*	THIS SUBROUTINE WILL COMPUTE THE SENSITIVITY OF VON MISES STRESS  *
*	WITH RESPECT TO DISPLACEMENT FOR EULER-BERNOULLI BEAM ELEMENT     *
*                                                                   *
*********************************************************************
*
*	AUTHOR  :  YOUNGSOO CHOI
*	DATE    :  JANUARY 2014
*
*********************************************************************
*
*		DEFINE THE GLOBAL VARIABLES
*
*	  AREA = CROSS-SECTIONAL AREA OF THE BEAM
*	     E = YOUNGS MODULUS FORTHE BEAM
*  ALPHA = DILATATION COEFFICIENT
*	   ELM = CURRENT ELEMENT NUMBER
* MAXSZE = LEADING DIMENSION OF STRESS AND STRAIN ARRAYS
*	MAXGUS = SECOND DIMENSION OF STRESS
*	MAXSTR = THIRD DIMENSION OF STRESS
*	EFRAME = ELEMENT REFERENCE FRAMES
* IX,IY,IZ = BEAM MOMENTS OF INERTIA
*	    NU = POISSON'S RATIO
*  X,Y,Z = COORDINATES FOR THE BEAM
*     UG = GLOBAL DISPLACEMENT VECTOR FOR ELEMENT #ELM
*   TEMP = NODAL TEMPERATURE
*   TREF = REFERENCE TEMPERATURE
*TSTRESS = THERMAL STRESS
*
********************************************************************/

	using std::sqrt;
	doublereal eiy, eiz, length2, length3, tstress;

//.... LOCAL VARIABLES

	int ii, jj, kk, ic;

	doublereal pi, G, J, dx, dy, dz, length;
	Eigen::Matrix<doublereal, 10, 1> vec;
	Eigen::Matrix<doublereal, 12, 12> ke, tran;
	Eigen::Matrix<doublereal, 12, 1> tug, res;
	Eigen::Matrix<doublereal, 3, 1> u, v, w;
	Eigen::Matrix<doublereal, 3, 3> t33;
	Eigen::Matrix<doublereal, 2, 1> tl;
	Eigen::Matrix<doublereal, 7, 2> stress;

	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > dvmsWRTdisp(_vmsWRTdisp, 2, 12);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > x(_x, 2, 1);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > y(_y, 2, 1);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > z(_z, 2, 1);
	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > eframe(_eframe, 9, 1);
	Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > ug(_ug, 12, 1);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > temp(_temp, 2, 1);



//.... INITIALIZE LOCAL CONSTANTS

	pi = 4.0 * atan(1.0);
	G = e / (2.0 * (1.0 + nu));
	J = ix;

//.... COMPUTE THE LENGTH OF THE BEAM

	dx = x[1] - x[0];
	dy = y[1] - y[0];
	dz = z[1] - z[0];
	length = sqrt(dx * dx + dy * dy + dz * dz);

//.... COMPUTE THE CONSTITUTIVE RELATIONS FOR THE BEAM

	eiy = e * iy;
	eiz = e * iz;
	length2 = length * length;
	length3 = length2 * length;
	vec[0] = e * area / length;
	vec[1] = 12.0 * eiy / length3;
	vec[2] = 12.0 * eiz / length3;
	vec[3] = 6.0 * eiy / length2;
	vec[4] = 6.0 * eiz / length2;
	vec[5] = G * J / length;
	vec[6] = 4.0 * eiy / length;
	vec[7] = 4.0 * eiz / length;
	vec[8] = 2.0 * eiy / length;
	vec[9] = 2.0 * eiz / length;

//.... INITIALIZE MATRICES AND VECTORS

	tug.setZero();
	res.setZero();
	ke.setZero();
	tran.setZero();

//.... BUILD THE ELEMENT STIFFNESS MATRIX

	ke(0, 0) = vec[0];
	ke(1, 1) = vec[2];
	ke(2, 2) = vec[1];
	ke(3, 3) = vec[5];
	ke(4, 4) = vec[6];
	ke(5, 5) = vec[7];
	ke(6, 6) = vec[0];
	ke(7, 7) = vec[2];
	ke(8, 8) = vec[1];
	ke(9, 9) = vec[5];
	ke(10, 10) = vec[6];
	ke(11, 11) = vec[7];

	ke(4, 2) = -vec[3];
	ke(2, 4) = ke(4, 2);
	ke(5, 1) = vec[4];
	ke(1, 5) = ke(5, 1);
	ke(6, 0) = -vec[0];
	ke(0, 6) = ke(6, 0);
	ke(7, 1) = -vec[2];
	ke(1, 7) = ke(7, 1);
	ke(7, 5) = -vec[4];
	ke(5, 7) = ke(7, 5);
	ke(8, 2) = -vec[1];
	ke(2, 8) = ke(8, 2);
	ke(8, 4) = vec[3];
	ke(4, 8) = ke(8, 4);
	ke(9, 3) = -vec[5];
	ke(3, 9) = ke(9, 3);
	ke(10, 2) = -vec[3];
	ke(2, 10) = ke(10, 2);
	ke(10, 4) = vec[8];
	ke(4, 10) = ke(10, 4);
	ke(10, 8) = vec[3];
	ke(8, 10) = ke(10, 8);
	ke(11, 1) = vec[4];
	ke(1, 11) = ke(11, 1);
	ke(11, 5) = vec[9];
	ke(5, 11) = ke(11, 5);
	ke(11, 7) = -vec[4];
	ke(7, 11) = ke(11, 7);

//        for(ii = 0; ii < 12; ++ii)
//          for(jj = 0; jj < 12; ++jj)
//            ke(jj,ii) = ke(ii,jj);

//.... COMPUTE THE TRANFORMATION MATRIX

	frame6(eframe.data(), u.data(), v.data(), w.data());

	t33 << u.transpose(), v.transpose(), w.transpose();

	for (kk = 0; kk < 4; ++kk) {
		ic = 3 * kk;
		for (ii = 0; ii < 3; ++ii)
			for (jj = 0; jj < 3; ++jj)
				tran(ic + ii, ic + jj) = t33(ii, jj);
	}

//.... ROTATE THE GLOBAL DISPLACEMENT VECTOR BACK TO LOCAL FRAME

	tug = tran * ug;

//.... COMPUTE THE INTERNAL FORCE RESULTANTS

	res = ke * tug;

//.... COMPUTE THERMAL STRESS
//.... BE CAREFUL!!!!! put tref to 0 in INPUT file if you don't want
//.... this term to be included
	if (_temp) {
		for (int i = 0; i < 2; ++i)
			tl[i] = temp[i] - tref;
	} else tl.setZero();

	tstress = e * alpha * area * (tl[0] + tl[1]) / 2.;

//.... WRITE THE RESULTANTS INTO THE STRESS ARRAY WHERE

//        FORX = STRESSXX      FORY=STRESSYY     FORZ=STRESSZZ
//        MOMX = STRESSXY      MOMY=STRESSXZ     MOMZ=STRESSYZ

	stress(0, 0) = res[0] + tstress;
	stress(1, 0) = res[1];
	stress(2, 0) = res[2];
	stress(3, 0) = res[3];
	stress(4, 0) = res[4];
	stress(5, 0) = res[5];
	stress(0, 1) = res[6] - tstress;
	stress(1, 1) = res[7];
	stress(2, 1) = res[8];
	stress(3, 1) = res[9];
	stress(4, 1) = res[10];
	stress(5, 1) = res[11];

//..... WRITE THE VON MISES RESULTANT INTO THE STRESS ARRAY
	doublereal stress_0dxy = stress(0, 0) - stress(1, 0);
	doublereal stress_0dyz = stress(1, 0) - stress(2, 0);
	doublereal stress_0dxz = stress(0, 0) - stress(2, 0);
	doublereal stress_1dxy = stress(0, 1) - stress(1, 1);
	doublereal stress_1dyz = stress(1, 1) - stress(2, 1);
	doublereal stress_1dxz = stress(0, 1) - stress(2, 1);

	stress(6, 0) = sqrt(0.5 * (stress_0dxy * stress_0dxy + stress_0dyz * stress_0dyz + stress_0dxz * stress_0dxz) +
	                    3 * (stress(3, 0) * stress(3, 0) + stress(4, 0) * stress(4, 0) + stress(5, 0) * stress(5, 0)));
	stress(6, 1) = sqrt(0.5 * (stress_1dxy * stress_1dxy + stress_1dyz * stress_1dyz + stress_1dxz * stress_1dxz) +
	                    3 * (stress(3, 1) * stress(3, 1) + stress(4, 1) * stress(4, 1) + stress(5, 1) * stress(5, 1)));

	Eigen::Matrix<doublereal, 2, 12> dvmsdStress;
	dvmsdStress.setZero();
	dvmsdStress(0, 0) = (2. * stress(0, 0) - stress(1, 0) - stress(2, 0)) / (2. * stress(6, 0));
	dvmsdStress(0, 1) = (2. * stress(1, 0) - stress(0, 0) - stress(2, 0)) / (2. * stress(6, 0));
	dvmsdStress(0, 2) = (2. * stress(2, 0) - stress(1, 0) - stress(0, 0)) / (2. * stress(6, 0));
	dvmsdStress(0, 3) = (3. * stress(3, 0)) / stress(6, 0);
	dvmsdStress(0, 4) = (3. * stress(4, 0)) / stress(6, 0);
	dvmsdStress(0, 5) = (3. * stress(5, 0)) / stress(6, 0);

	dvmsdStress(1, 6) = (2. * stress(0, 1) - stress(1, 1) - stress(2, 1)) / (2. * stress(6, 1));
	dvmsdStress(1, 7) = (2. * stress(1, 1) - stress(0, 1) - stress(2, 1)) / (2. * stress(6, 1));
	dvmsdStress(1, 8) = (2. * stress(2, 1) - stress(1, 1) - stress(0, 1)) / (2. * stress(6, 1));
	dvmsdStress(1, 9) = (3. * stress(3, 1)) / stress(6, 1);
	dvmsdStress(1, 10) = (3. * stress(4, 1)) / stress(6, 1);
	dvmsdStress(1, 11) = (3. * stress(5, 1)) / stress(6, 1);

	dvmsWRTdisp = dvmsdStress * ke * tran;
}

template<typename doublereal>
void
BeamElementTemplate<doublereal>
::buildFrameInTemplate(doublereal *_x, doublereal *_y, doublereal *_z, doublereal *_eframe) {
	Eigen::Map<Eigen::Matrix<doublereal, 9, 1> > eframe(_eframe);
	Eigen::Map<Eigen::Matrix<doublereal, 2, 1> > x(_x), y(_y), z(_z);

	Eigen::Matrix<doublereal, 3, 1> frame1, frame2, frame3;
	frame2[0] = eframe[3];
	frame2[1] = eframe[4];
	frame2[2] = eframe[5];
	frame3[0] = eframe[6];
	frame3[1] = eframe[7];
	frame3[2] = eframe[8];

	frame1[0] = x[1] - x[0];
	frame1[1] = y[1] - y[0];
	frame1[2] = z[1] - z[0];
	frame1.normalize();
	frame3 = frame1.cross(frame2);
	frame3.normalize();
	frame2 = frame3.cross(frame1);
	frame2.normalize();

	eframe[0] = frame1[0];
	eframe[1] = frame1[1];
	eframe[2] = frame1[2];
	eframe[3] = frame2[0];
	eframe[4] = frame2[1];
	eframe[5] = frame2[2];
	eframe[6] = frame3[0];
	eframe[7] = frame3[1];
	eframe[8] = frame3[2];

}

#endif
#endif
