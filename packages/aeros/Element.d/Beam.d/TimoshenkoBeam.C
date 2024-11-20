#include    <Element.d/Beam.d/TimoshenkoBeam.h>
#include    <Element.d/Beam.d/TimoshenkoBeamStressWRTDisplacementSensitivity.h>
#include    <Element.d/Beam.d/TimoshenkoBeamStressWRTNodalCoordinateSensitivity.h>
#include    <Element.d/Beam.d/TimoshenkoBeamGravityForceWRTNodalCoordinateSensitivity.h>
#include    <Element.d/Beam.d/TimoshenkoBeamStiffnessWRTNodalCoordinateSensitivity.h>
#include    <Element.d/Function.d/SpaceDerivatives.h>
#include    <Math.d/FullSquareMatrix.h>
#include    <Math.d/matrix.h>
#include    <Corotational.d/BeamCorotator.h>
#include    <Corotational.d/GeomState.h>
#include    <Corotational.d/utilities.h>
#include    <Utils.d/dofset.h>
#include    <Utils.d/linkfc.h>

#include    <cstdio>
#include    <cstddef>
#include    <cmath>
#include    <cstring>

extern int verboseFlag;
extern "C" {
void _FORTRAN(modmstif7)(double *, double &, double &, double *,
                         double &, double &, double &, double &, double &,
                         double &, double &, double *, double *, double *,
                         int &);

void _FORTRAN(mass7)(int &, double *, double &, double &,
                     double *, double *, double *, double *, double *,
                     const int &, double &, const int &);

void _FORTRAN(sands7)(const int &, double &, double &, double *, double &,
                      double &, double &, double &, double &, double &,
                      double &, double *, double *, double *, double *,
                      double *, const int &, const int &, const int &, const int &,
                      double &, double &, double *);

void _FORTRAN(transform)(double *, double *, double *, double *, double *, double *, double *);
}

TimoshenkoBeam::TimoshenkoBeam(int *nodenums) {
	elemframe = 0;
	myElemFrame = false;
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	iniOr = 0;
	pbc = 0;
}

TimoshenkoBeam::~TimoshenkoBeam() {
	if (myElemFrame) delete[] elemframe;
	delete[] iniOr;
}

TimoshenkoBeam *
TimoshenkoBeam::clone() {
	TimoshenkoBeam *result = new TimoshenkoBeam(nn);

	if (elemframe) {
		std::cerr << "in TimoshenkoBeam::clone\n";
		result->elemframe = new EFrame[1];
		result->myElemFrame = true;
		std::memcpy(result->elemframe, elemframe, sizeof(elemframe[0]));
	}

	if (iniOr) {
		result->iniOr = new double[3];
		std::copy(iniOr, iniOr + 3, result->iniOr);
	}

	return result;
}

void
TimoshenkoBeam::renum(const int *table) {
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
TimoshenkoBeam::renum(EleRenumMap &table) {
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
TimoshenkoBeam::buildFrame(CoordSet &cs) {
	// store initial orientation of beam
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	iniOr = new double[3];

	iniOr[0] = nd2.x - nd1.x;
	iniOr[1] = nd2.y - nd1.y;
	iniOr[2] = nd2.z - nd1.z;

	double len = iniOr[0] * iniOr[0] + iniOr[1] * iniOr[1] + iniOr[2] * iniOr[2];

	len = std::sqrt(len);

	iniOr[0] = iniOr[0] / len;
	iniOr[1] = iniOr[1] / len;
	iniOr[2] = iniOr[2] / len;

	oeframe[0][0] = 1;
	oeframe[0][1] = 0;
	oeframe[0][2] = 0;
	oeframe[1][0] = 0;
	oeframe[1][1] = 1;
	oeframe[1][2] = 0;
	oeframe[2][0] = 0;
	oeframe[2][1] = 0;
	oeframe[2][2] = 1;

	if (nn[2] < 0) {  // PJSA 4-8-05 copied from EulerBeam
		// only 2nd axis from EFRAMES is actually used (to define local xy plane)
		if (elemframe != 0) {
			oeframe[0][0] = (*elemframe)[0][0];
			oeframe[0][1] = (*elemframe)[0][1];
			oeframe[0][2] = (*elemframe)[0][2];
			oeframe[1][0] = (*elemframe)[1][0];
			oeframe[1][1] = (*elemframe)[1][1];
			oeframe[1][2] = (*elemframe)[1][2];
			oeframe[2][0] = (*elemframe)[2][0];
			oeframe[2][1] = (*elemframe)[2][1];
			oeframe[2][2] = (*elemframe)[2][2];
			EFrame &theFrame = *elemframe;
			theFrame[0][0] = nd2.x - nd1.x;
			theFrame[0][1] = nd2.y - nd1.y;
			theFrame[0][2] = nd2.z - nd1.z;
			normalize(theFrame[0]);
			crossprod(theFrame[0], theFrame[1], theFrame[2]);
			normalize(theFrame[2]);
			crossprod(theFrame[2], theFrame[0], theFrame[1]);
		}
		return;
	} else {
		auto &nd3 = cs.getNode(nn[2]);

		delete[] elemframe;
		elemframe = new EFrame[1];
		myElemFrame = true;
		EFrame &theFrame = *elemframe;
		theFrame[0][0] = nd2.x - nd1.x;
		theFrame[0][1] = nd2.y - nd1.y;
		theFrame[0][2] = nd2.z - nd1.z;
		normalize(theFrame[0]);

		double xz[3] = {nd3.x - nd2.x, nd3.y - nd2.y, nd3.z - nd2.z};
		crossprod(xz, theFrame[0], theFrame[1]);
		normalize(theFrame[1]);
		crossprod(theFrame[0], theFrame[1], theFrame[2]);
		return;
	}
}


void
TimoshenkoBeam::getIntrnForce(Vector &elForce, CoordSet &cs,
                              double *elDisp, int forceIndex, double *ndTemps) {
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	int maxsze = 1;
	int maxstr = 7;
	int maxgus = 2;
	int elm = 1;
	int numel = 1;

	double elStress[2][7];

	_FORTRAN(sands7)(elm, prop->A, prop->E, (double *) *elemframe,
	                 prop->Ixx, prop->Iyy, prop->Izz, prop->alphaY,
	                 prop->alphaZ, prop->c, prop->nu, x, y, z, elDisp,
	                 (double *) elStress, numel, maxgus, maxstr, maxsze,
	                 prop->W, prop->Ta, ndTemps);

// forceIndex    = 0 =  FORCE_X
//               = 1 =  FORCE_Y
//               = 2 =  FORCE_Z
//               = 3 =  MOMENT_X
//               = 4 =  MOMENT_Y
//               = 5 =  MOMENT_Z

	elForce[0] = elStress[0][forceIndex];
	elForce[1] = elStress[1][forceIndex];
}

double
TimoshenkoBeam::getMass(const CoordSet &cs) const {
	// Check for phantom element, which has no mass
	if (prop == NULL)
		return 0;

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = std::sqrt(dx * dx + dy * dy + dz * dz);

	double mass = length * (prop->rho) * (prop->A);

	return mass;
}

void
TimoshenkoBeam::getMassNodalCoordinateSensitivity(CoordSet &cs, Vector &dMassdx) {
	using std::sqrt;

	if (dMassdx.size() != 6) {
		std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
		exit(-1);
	}

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt(dx * dx + dy * dy + dz * dz);
	double prefix = prop->A * prop->rho / length;
	dMassdx[0] = -prefix * dx;
	dMassdx[1] = -prefix * dy;
	dMassdx[2] = -prefix * dz;
	dMassdx[3] = prefix * dx;
	dMassdx[4] = prefix * dy;
	dMassdx[5] = prefix * dz;
}

void
TimoshenkoBeam::getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet &cs, double *gravityAcceleration) {
	if (dwdx.size() != 6) {
		std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
		exit(-1);
	}
	Vector dMassdx(6);
	getMassNodalCoordinateSensitivity(cs, dMassdx);
	double gravAccNorm = sqrt(gravityAcceleration[0] * gravityAcceleration[0] +
	                          gravityAcceleration[1] * gravityAcceleration[1] +
	                          gravityAcceleration[2] * gravityAcceleration[2]);
	dwdx = gravAccNorm * dMassdx;
}

double
TimoshenkoBeam::weight(CoordSet &cs, double *gravityAcceleration) {
	double _mass = getMass(cs);
	double gravAccNorm = sqrt(gravityAcceleration[0] * gravityAcceleration[0] +
	                          gravityAcceleration[1] * gravityAcceleration[1] +
	                          gravityAcceleration[2] * gravityAcceleration[2]);
	return _mass * gravAccNorm;
}

void
TimoshenkoBeam::getGravityForceNodalCoordinateSensitivity(CoordSet &cs, double *gravityAcceleration,
                                                          GenFullM<double> &dGfdx, int gravflg, GeomState *geomState) {
#ifdef USE_EIGEN3
	double massPerNode = 0.5 * getMass(cs);
	if (prop == NULL) {
		dGfdx.zero();
		return;
	}

	double x[2] = {cs[nn[0]]->x, cs[nn[1]]->x};
	double y[2] = {cs[nn[0]]->y, cs[nn[1]]->y};
	double z[2] = {cs[nn[0]]->z, cs[nn[1]]->z};
	Eigen::Array<double, 14, 1> dconst;
	dconst.segment<3>(0) = Eigen::Map<Eigen::Matrix<double, 3, 1> >(gravityAcceleration).segment(0, 3);
	dconst[3] = prop->rho;
	dconst[4] = prop->A;
	dconst[5] = (oeframe)[0][0];
	dconst[6] = (oeframe)[0][1];
	dconst[7] = (oeframe)[0][2];
	dconst[8] = (oeframe)[1][0];
	dconst[9] = (oeframe)[1][1];
	dconst[10] = (oeframe)[1][2];
	dconst[11] = (oeframe)[2][0];
	dconst[12] = (oeframe)[2][1];
	dconst[13] = (oeframe)[2][2];
	Eigen::Array<int, 1, 1> iconst;
	iconst[0] = gravflg;

	Eigen::Matrix<double, 6, 1> q;
	q << x[0], y[0], z[0], x[1], y[1], z[1];
	Eigen::Matrix<double, 12, 6> dGravityForcedx;

	Simo::Jacobian<double, TimoshenkoBeamGravityForceWRTNodalCoordinateSensitivity> dGdx(dconst, iconst);
	dGravityForcedx = dGdx(q, 0);
	dGfdx.copy(dGravityForcedx.data());
#endif
}

void
TimoshenkoBeam::getGravityForce(CoordSet &cs, double *gravityAcceleration,
                                Vector &gravityForce, int gravflg, GeomState *geomState) {
	double massPerNode = 0.5 * getMass(cs);

	double t0n[3][3] = {{0.0, 0.0, 0.0},
	                    {0.0, 0.0, 0.0},
	                    {0.0, 0.0, 0.0}};
	double length;

	if (geomState) {

		updTransMatrix(cs, geomState, t0n, length);

	} else {

		auto &nd1 = cs.getNode(nn[0]);
		auto &nd2 = cs.getNode(nn[1]);

		double x[2], y[2], z[2];

		x[0] = nd1.x;
		y[0] = nd1.y;
		z[0] = nd1.z;
		x[1] = nd2.x;
		y[1] = nd2.y;
		z[1] = nd2.z;

		double dx = x[1] - x[0];
		double dy = y[1] - y[0];
		double dz = z[1] - z[0];

		length = std::sqrt(dx * dx + dy * dy + dz * dz);

		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				t0n[i][j] = (*elemframe)[i][j];
			}
		}
	}

	// Consistent

	int i;
	double localg[3];

	for (i = 0; i < 3; ++i)
		localg[i] = 0.0;

	for (i = 0; i < 3; ++i) {
		localg[0] += t0n[0][i] * gravityAcceleration[i];
		localg[1] += t0n[1][i] * gravityAcceleration[i];
		localg[2] += t0n[2][i] * gravityAcceleration[i];
	}
	double localf[3], localm[3];
	double globalf[3], globalm[3];
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

FullSquareMatrix
TimoshenkoBeam::massMatrix(const CoordSet &cs, double *mel, int cmflg) const {
	// Check for phantom element, which has no mass
	if (prop == NULL) {
		FullSquareMatrix ret(12, mel);
		ret.zero();
		return ret;
	}

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;

	int elm = 0, grvflg = 0, masflg = 0;

	_FORTRAN(mass7)(elm, (double *) mel, prop->A, prop->rho, x, y, z,
	                gravityAcceleration, grvfor, grvflg, totmas, masflg);

	FullSquareMatrix ret(12, mel);

	return ret;
}

void
TimoshenkoBeam::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs) {
#ifdef USE_EIGEN3
	for (int i = 0; i < 6; ++i) {
		if (dStiffdx[i].dim() != 12) {
			std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
			exit(-1);
		} else dStiffdx[i].zero();
	}

	// Check for phantom element, which has no stiffness
	if (prop == NULL) return;

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	Eigen::Array<double, 18, 1> dconst;
	Eigen::Array<int, 0, 1> iconst;
	Eigen::Array<double, 6, 1> q;

	dconst[0] = prop->E;
	dconst[1] = prop->A;
	dconst[2] = prop->Ixx;
	dconst[3] = prop->Iyy;
	dconst[4] = prop->Izz;
	dconst[5] = prop->alphaY;
	dconst[6] = prop->alphaZ;
	dconst[7] = prop->C1;
	dconst[8] = prop->nu;
	dconst[9] = oeframe[0][0];
	dconst[10] = (oeframe)[0][1];
	dconst[11] = (oeframe)[0][2];
	dconst[12] = (oeframe)[1][0];
	dconst[13] = (oeframe)[1][1];
	dconst[14] = (oeframe)[1][2];
	dconst[15] = (oeframe)[2][0];
	dconst[16] = (oeframe)[2][1];
	dconst[17] = (oeframe)[2][2];

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	q << nd1.x, nd1.y, nd1.z, nd2.x, nd2.y, nd2.z;

	Eigen::Array<Eigen::Matrix<double, 12, 12>, 1, 6> dStiffnessdx;

	Simo::FirstPartialSpaceDerivatives<double, TimoshenkoBeamStiffnessWRTNodalCoordinateSensitivity> dKdx(dconst,
	                                                                                                      iconst);
	dStiffnessdx = dKdx(q, 0);
	for (int i = 0; i < 6; ++i) dStiffdx[i].copy(dStiffnessdx[i].data());
#endif
}

FullSquareMatrix
TimoshenkoBeam::stiffness(const CoordSet &cs, double *d, int flg) const {
	// Check for phantom element, which has no stiffness
	if (prop == NULL) {
		FullSquareMatrix ret(12, d);
		ret.zero();
		return ret;
	}

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

// 	EFrame is stored in following order. 
//			X_x, X_y, X_z, 
//			Y_x, Y_y, Y_z,
//      Z_x, Z_y, Z_z

	if (prop->A == 0.0)
		std::fprintf(stderr, "ERROR: Timoshenko beam has zero area. nodes %d %d\n",
		             nn[0] + 1, nn[1] + 1);
	if (std::sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) + (z[0] - z[1]) * (z[0] - z[1])) == 0.0)
		std::fprintf(stderr, "ERROR: Timoshenko beam has zero length. nodes %d %d\n",
		             nn[0] + 1, nn[1] + 1);

	// Check for the frame
	if (elemframe == 0) {
		std::fprintf(stderr, " ****************************************************\n");
		fprintf(stderr, " *** ERROR: Timoshenko beam lacks a frame"
				" (Nodes %d %d)\n", nn[0] + 1, nn[1] + 1);
		std::fprintf(stderr, " ****************************************************\n");
		exit(-1);
	}
	_FORTRAN(modmstif7)((double *) d, prop->A, prop->E,
	                    (double *) *elemframe,
	                    prop->Ixx, prop->Iyy, prop->Izz,
	                    prop->alphaY, prop->alphaZ, prop->C1,
	                    prop->nu, x, y, z, flg);

	FullSquareMatrix ret(12, d);

	return ret;
}

int
TimoshenkoBeam::numNodes() const {
	return 2;
}

int *
TimoshenkoBeam::nodes(int *p) const {
	if (p == 0) p = new int[2];
	p[0] = nn[0];
	p[1] = nn[1];
	return p;
}

int
TimoshenkoBeam::numDofs() const {
	return 12;
}

int *
TimoshenkoBeam::dofs(DofSetArray &dsa, int *p) const  {

	if (p == 0) p = new int[12];

	dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p);
	dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p + 6);

	return p;
}

void
TimoshenkoBeam::markDofs(DofSetArray &dsa) const {
	dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
}

Corotator *
TimoshenkoBeam::getCorotator(CoordSet &cs, double *kel, int, int fitAlgBeam) {
	int flag = 0;
	FullSquareMatrix myStiff = stiffness(cs, kel, flag);
	return new BeamCorotator(nn[0], nn[1], (*elemframe)[2], myStiff, fitAlgBeam);
}

int
TimoshenkoBeam::getTopNumber() const {
	return 107;//1;
}

void
TimoshenkoBeam::computePressureForce(CoordSet &cs, Vector &elPressureForce,
                                     GeomState *geomState, int cflg, double) {
	double normal[3], normal2[3];
	double px = 0.0;
	double py = 0.0;
	double pz = 0.0;
	double length;

	if (geomState) {
		double t0n[3][3] = {{0.0, 0.0, 0.0},
		                    {0.0, 0.0, 0.0},
		                    {0.0, 0.0, 0.0}};
		updTransMatrix(cs, geomState, t0n, length);
		normal[0] = t0n[1][0];
		normal[1] = t0n[1][1];
		normal[2] = t0n[1][2];
		normal2[0] = t0n[2][0];
		normal2[1] = t0n[2][1];
		normal2[2] = t0n[2][2];
	} else {
		// Obtain normal to beam (second vector in element frame)
		normal[0] = (*elemframe)[1][0];
		normal[1] = (*elemframe)[1][1];
		normal[2] = (*elemframe)[1][2];
		normal2[0] = (*elemframe)[2][0];
		normal2[1] = (*elemframe)[2][1];
		normal2[2] = (*elemframe)[2][2];

		// Get length of Beam
		auto &nd1 = cs.getNode(nn[0]);
		auto &nd2 = cs.getNode(nn[1]);

		double x[2], y[2], z[2];

		x[0] = nd1.x;
		y[0] = nd1.y;
		z[0] = nd1.z;
		x[1] = nd2.x;
		y[1] = nd2.y;
		z[1] = nd2.z;

		double dx = x[1] - x[0];
		double dy = y[1] - y[0];
		double dz = z[1] - z[0];

		length = std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	double pressureForce = 0.5 * pbc->val * length;
	px = pressureForce * normal[0];
	py = pressureForce * normal[1];
	pz = pressureForce * normal[2];

	// Consistent
	double localMz = pressureForce * length / 6.0;
	double mx = localMz * normal2[0];
	double my = localMz * normal2[1];
	double mz = localMz * normal2[2];
	elPressureForce[0] = px;
	elPressureForce[1] = py;
	elPressureForce[2] = pz;
	elPressureForce[3] = mx;
	elPressureForce[4] = my;
	elPressureForce[5] = mz;
	elPressureForce[6] = px;
	elPressureForce[7] = py;
	elPressureForce[8] = pz;
	elPressureForce[9] = -mx;
	elPressureForce[10] = -my;
	elPressureForce[11] = -mz;
}

void
TimoshenkoBeam::getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &elementThermalForce,
                                int glflag, GeomState *geomState) {
	if (prop == NULL) {
		elementThermalForce.zero();
		return;
	}

	double localThF[12];
	int i, j;

	double Tref = prop->Ta;
	double coeff = prop->E * prop->W * prop->A;
	double length;

	// Local Thermal Forces : There are only axial forces (see notes)
	// indices 0-5 -->node1, 6-11 -->node2

	double deltaT1 = ndTemps[0] - Tref;
	double deltaT2 = ndTemps[1] - Tref;

	for (i = 0; i < 12; i++) localThF[i] = 0.;

	localThF[0] = -coeff * (0.5 * deltaT1 + 0.5 * deltaT2);
	localThF[6] = coeff * (0.5 * deltaT1 + 0.5 * deltaT2);

	// Compute Global Thermal Forces: From local to global -->
	//                                transpose(Transform. Matrix)

	double t0n[3][3] = {{0., 0., 0.,},
	                    {0., 0., 0.},
	                    {0., 0., 0.}};

	if (geomState) {

		updTransMatrix(cs, geomState, t0n, length);

	} else {
		for (i = 0; i < 3; ++i) {
			for (j = 0; j < 3; ++j) {
				t0n[i][j] = (*elemframe)[i][j];
			}
		}
	}

	for (i = 0; i < 3; ++i) {
		elementThermalForce[i] = 0.;
		for (j = 0; j < 3; ++j) {
			elementThermalForce[i] += t0n[j][i] * localThF[j];
		}
	}
	for (i = 0; i < 3; ++i) {
		elementThermalForce[i + 3] = 0.;
		for (j = 0; j < 3; ++j) {
			elementThermalForce[i + 3] += t0n[j][i] * localThF[j + 3];
		}
	}

	for (i = 0; i < 3; ++i) {
		elementThermalForce[i + 6] = 0.;
		for (j = 0; j < 3; ++j) {
			elementThermalForce[i + 6] += t0n[j][i] * localThF[j + 6];
		}
	}
	for (i = 0; i < 3; ++i) {
		elementThermalForce[i + 9] = 0.;
		for (j = 0; j < 3; ++j) {
			elementThermalForce[i + 9] += t0n[j][i] * localThF[j + 9];
		}
	}
}

void
TimoshenkoBeam::getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                            Vector &elDisp, int strInd, int surface,
                            double *ndTemps, double ylayer, double zlayer, int avgnum) {
	weight = 1.0;
	if (strInd == -1) return;

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	// Find the internal forces and moments

	int maxsze = 1;
	int maxstr = 7;
	int maxgus = 2;
	int elm = 1;
	int numel = 1;

	double elStress[2][7];
	double elForce[3][2] = {{0.0, 0.0},
	                        {0.0, 0.0},
	                        {0.0, 0.0}};

	_FORTRAN(sands7)(elm, prop->A, prop->E, (double *) *elemframe,
	                 prop->Ixx, prop->Iyy, prop->Izz, prop->alphaY,
	                 prop->alphaZ, prop->c, prop->nu, x, y, z, elDisp.data(),
	                 (double *) elStress, numel, maxgus, maxstr, maxsze,
	                 prop->W, prop->Ta, ndTemps);

	// elForce[0] -> Axial Force (x-direction)
	// elForce[1] -> Moment around the y-axis (My)
	// elForce[2] -> Moment around the z-axis (Mz)

	elForce[0][0] = -elStress[0][0];
	elForce[1][0] = -elStress[0][4];
	elForce[2][0] = -elStress[0][5];
	elForce[0][1] = elStress[1][0];
	elForce[1][1] = elStress[1][4];
	elForce[2][1] = elStress[1][5];


	// Y is the position in the cross section where the stress is to be
	// calculated and is defined as Y = -ylayer * ymin if ylayer is < 0.0 (ymin
	// is negative in the input file) and Y = ylayer * ymax if ylayer is > 0.0
	// (ymax is positive), where ylayer is the percentage of the position ymin
	// or ymax.
	// Same thing for Z

	double Y = 0.0;
	double Z = 0.0;

	if (ylayer < 0.0) {
		Y = -ylayer * prop->ymin;
	} else if (ylayer > 0.0) {
		Y = ylayer * prop->ymax;
	}

	if (zlayer < 0.0) {
		Z = -zlayer * prop->zmin;
	} else if (zlayer > 0.0) {
		Z = zlayer * prop->zmax;
	}


	switch (avgnum) {

		case 0: // elemental
		{
			if (strInd == 0) {

				// Axial Stress
				double IY = prop->Iyy;
				double IZ = prop->Izz;
				double cA = prop->A;

				stress[0] = elForce[0][0] / cA - elForce[2][0] * Y / IZ + elForce[1][0] * Z / IY;
				stress[1] = elForce[0][1] / cA - elForce[2][1] * Y / IZ + elForce[1][1] * Z / IY;

			} else if (strInd == 6) {

				// von Mises stress resultant
				stress[0] = elStress[0][6];
				stress[1] = elStress[1][6];

			} else if (strInd == 7) {

				// Axial Strain

				double EA = prop->E * prop->A;
				double EIZ = prop->E * prop->Izz;
				double EIY = prop->E * prop->Iyy;

				double Tref = prop->Ta;
				double alpha = prop->W;

				double dT1 = ndTemps[0] - Tref;
				double dT2 = ndTemps[1] - Tref;

				double localThS;
				localThS = alpha * (0.5 * dT1 + 0.5 * dT2);

				stress[0] = elForce[0][0] / EA - elForce[2][0] * Y / EIZ + elForce[1][0] * Z / EIY + localThS;
				stress[1] = elForce[0][1] / EA - elForce[2][1] * Y / EIZ + elForce[1][1] * Z / EIY + localThS;

			} else {
				stress[0] = 0.0;
				stress[1] = 0.0;
			}
			break;
		}

		case 1:
		case 3: // nodalfull or nodalpartialgroup
		{
			if (strInd == 6) {

				// von Mises stress resultant
				stress[0] = elStress[0][6];
				stress[1] = elStress[1][6];

			} else if (strInd < 6) {
				// Axial Stress

				double IY = prop->Iyy;
				double IZ = prop->Izz;
				double cA = prop->A;

				double xg[3][3] = {{1.0, 0.0, 0.0},
				                   {0.0, 1.0, 0.0},
				                   {0.0, 0.0, 1.0}};
				double tmpStr1[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double tmpStr2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				tmpStr1[0] = elForce[0][0] / cA - elForce[2][0] * Y / IZ + elForce[1][0] * Z / IY;
				tmpStr2[0] = elForce[0][1] / cA - elForce[2][1] * Y / IZ + elForce[1][0] * Z / IY;

         _FORTRAN(transform)((double*)(*elemframe)[0], (double*)(*elemframe)[1], (double*)(*elemframe)[2],
				                    xg[0], xg[1], xg[2], tmpStr1);

         _FORTRAN(transform)((double*)(*elemframe)[0], (double*)(*elemframe)[1], (double*)(*elemframe)[2],
				                    xg[0], xg[1], xg[2], tmpStr2);

				stress[0] = tmpStr1[strInd];
				stress[1] = tmpStr2[strInd];

			} else if (strInd > 6 && strInd < 13) {

				// Axial Strain

				double EA = prop->E * prop->A;
				double EIZ = prop->E * prop->Izz;
				double EIY = prop->E * prop->Iyy;

				double Tref = prop->Ta;
				double alpha = prop->W;

				double dT1 = ndTemps[0] - Tref;
				double dT2 = ndTemps[1] - Tref;

				double localThS;
				localThS = alpha * (0.5 * dT1 + 0.5 * dT2);

				double xg[3][3] = {{1.0, 0.0, 0.0},
				                   {0.0, 1.0, 0.0},
				                   {0.0, 0.0, 1.0}};
				double tmpStr1[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double tmpStr2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				tmpStr1[0] = elForce[0][0] / EA - elForce[2][0] * Y / EIZ + elForce[1][0] * Z / EIY + localThS;
				tmpStr2[0] = elForce[0][1] / EA - elForce[2][1] * Y / EIZ + elForce[1][1] * Z / EIY + localThS;

         _FORTRAN(transform)((double*)(*elemframe)[0], (double*)(*elemframe)[1], (double*)(*elemframe)[2],
				                    xg[0], xg[1], xg[2], tmpStr1);

         _FORTRAN(transform)((double*)(*elemframe)[0], (double*)(*elemframe)[1], (double*)(*elemframe)[2],
				                    xg[0], xg[1], xg[2], tmpStr2);

				stress[0] = tmpStr1[strInd - 7];
				stress[1] = tmpStr2[strInd - 7];

			} else {
				stress[0] = 0.0;
				stress[1] = 0.0;
			}
			break;
		}

		case 2: // nodalpartial
		{
			weight = 0.0;
			stress[0] = 0.0;
			stress[1] = 0.0;
			break;
		}

		default:
			std::cerr << "avgnum = " << avgnum << " is not a valid number\n";
	}
}

#ifdef USE_EIGEN3

void
TimoshenkoBeam::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight,
                                                   GenFullM<double> *dDispDisp,
                                                   CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                   double *ndTemps, int avgnum, double ylayer, double zlayer) {
	if (strInd != 6) {
		std::cerr << " ... Error: strInd must be 6 in TimoshenkoBeam::getVonMisesDisplacementSensitivity\n";
		exit(-1);
	}
	if (dStdDisp.numRow() != 12 || dStdDisp.numCol() != 2) {
		std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
		exit(-1);
	}
	weight = 1;
	// scalar parameters
	Eigen::Array<double, 28, 1> dconst;
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	dconst[0] = nd1.x;
	dconst[1] = nd2.x; // x coordinates
	dconst[2] = nd1.y;
	dconst[3] = nd2.y; // y coordinates
	dconst[4] = nd1.z;
	dconst[5] = nd2.z; // z coordinates
	dconst[6] = prop->A;
	dconst[7] = prop->E;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			dconst[8 + 3 * i + j] = (*elemframe)[i][j];
		}
	}
	dconst[17] = prop->Ixx;
	dconst[18] = prop->Iyy;
	dconst[19] = prop->Izz;
	dconst[20] = prop->alphaY;
	dconst[21] = prop->alphaZ;
	dconst[22] = prop->c;
	dconst[23] = prop->nu;
	dconst[24] = prop->W;
	dconst[25] = prop->Ta;
	if (ndTemps) {
		dconst[26] = ndTemps[0];
		dconst[27] = ndTemps[1];
	} else {
		dconst[26] = 0.0;
		dconst[27] = 0.0;
	}

	// integer parameters
	Eigen::Array<int, 1, 1> iconst;
	iconst[0] = avgnum;
	// inputs
	Eigen::Matrix<double, 12, 1> q = Eigen::Map<Eigen::Matrix<double, 12, 1> >(elDisp.data()).segment(0,
	                                                                                                  12); // displacements

	//Jacobian evaluation
	Eigen::Matrix<double, 2, 12> dStressdDisp;
	dStressdDisp.setZero();
	Eigen::Matrix<double, 7, 3> stress;
	if (avgnum == 1 || avgnum == 0 || avgnum == 3) { // ELEMENTAL or NODALFULL or NODALPARTIALGROUP

		dStressdDisp.setZero();
		Eigen::Matrix<double, 9, 1> eframe = Eigen::Map<Eigen::Matrix<double, 28, 1> >(dconst.data()).segment(8,
		                                                                                                      9); // extract eframe
		vms7WRTdisp(1, prop->A, prop->E, eframe.data(), prop->Ixx, prop->Iyy, prop->Izz, prop->alphaY, prop->alphaZ,
		            prop->c,
		            prop->nu, x, y, z, q.data(), dStressdDisp.data(), prop->W, prop->Ta, ndTemps);
	}
	dStdDisp.copy(dStressdDisp.data());

	if (dDispDisp) dStdDisp ^= (*dDispDisp);
}

void
TimoshenkoBeam::getVonMisesDisplacementSensitivity(GenFullM<DComplex> &dStdDisp, ComplexVector &weight,
                                                   CoordSet &cs, ComplexVector &elDisp, int strInd, int surface,
                                                   double *, int avgnum, double ylayer, double zlayer) {
	std::cerr << "TimoshenkoBeam::getVonMisesDisplacementSensitivity is not implemented\n";
	exit(-1);
}

void
TimoshenkoBeam::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                                      Vector &elDisp, int strInd, int surface,
                                                      double *ndTemps, int avgnum, double ylayer, double zlayer) {
	if (strInd != 6) {
		std::cerr << " ... Error: strInd must be 6 in TimoshenkoBeam::getVonMisesNodalCoordinateSensitivity\n";
		exit(-1);
	}
	if (dStdx.numRow() != 6 || dStdx.numCol() != 2) {
		std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
		exit(-1);
	}
	weight = 1;
	// scalar parameters
	Eigen::Array<double, 34, 1> dconst;
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);

	double x[2], y[2], z[2];

	x[0] = nd1.x;
	y[0] = nd1.y;
	z[0] = nd1.z;
	x[1] = nd2.x;
	y[1] = nd2.y;
	z[1] = nd2.z;

	dconst.segment<12>(0) = Eigen::Map<Eigen::Matrix<double, 12, 1> >(elDisp.data()).segment(0, 12); // displacements
	dconst[12] = prop->A;
	dconst[13] = prop->E;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			dconst[14 + 3 * i + j] = oeframe[i][j];
		}
	}
	dconst[23] = prop->Ixx;
	dconst[24] = prop->Iyy;
	dconst[25] = prop->Izz;
	dconst[26] = prop->alphaY;
	dconst[27] = prop->alphaZ;
	dconst[28] = prop->c;
	dconst[29] = prop->nu;
	dconst[30] = prop->W;
	dconst[31] = prop->Ta;
	if (ndTemps) {
		dconst[32] = ndTemps[0];
		dconst[33] = ndTemps[1];
	} else {
		dconst[32] = 0.0;
		dconst[33] = 0.0;
	}

	// integer parameters
	Eigen::Array<int, 1, 1> iconst;
	iconst[0] = avgnum;

	// input
	Eigen::Matrix<double, 6, 1> q;
	q << nd1.x, nd1.y, nd1.z, nd2.x, nd2.y, nd2.z;

	//Jacobian evaluation
	Eigen::Matrix<double, 2, 6> dStressdx;
	Eigen::Matrix<double, 7, 3> stress;

	if (avgnum == 1 || avgnum == 0 || avgnum == 3) { // ELEMENTAL or NODALFULL or NODALPARTIALGROUP

		Simo::Jacobian<double, TimoshenkoBeamStressWRTNodalCoordinateSensitivity> dSdx(dconst, iconst);
		if (elDisp.norm() == 0) dStressdx.setZero();
		else dStressdx = dSdx(q, 0);
		dStdx.copy(dStressdx.data());
	} else dStdx.zero(); // NODALPARTIAL or GAUSS or any others
}

#endif

void
TimoshenkoBeam::updTransMatrix(CoordSet &cs, GeomState *geomState, double t0n[3][3], double &length) {
// Returns t0n[3][3] and length

	double xn[2][3];

	double zVecL[2][3];
	double (*rot[2])[3][3];

	// Get Nodes current coordinates
	NodeState &ns1 = (*geomState)[nn[0]];
	NodeState &ns2 = (*geomState)[nn[1]];

	xn[0][0] = ns1.x; // x coordinate of node state 1
	xn[0][1] = ns1.y; // y coordinate of node state 1
	xn[0][2] = ns1.z; // z coordinate of node state 1

	xn[1][0] = ns2.x; // x coordinate of node state 2
	xn[1][1] = ns2.y; // y coordinate of node state 2
	xn[1][2] = ns2.z; // z coordinate of node state 2

	double dx = xn[1][0] - xn[0][0];
	double dy = xn[1][1] - xn[0][1];
	double dz = xn[1][2] - xn[0][2];

	length = std::sqrt(dx * dx + dy * dy + dz * dz);

	rot[0] = &(ns1.R); // rotation tensor of node state 1
	rot[1] = &(ns2.R); // rotation tensor of node state 2

// Compute nodal rotated Z-axis in global coordinate system

	int i, nod;
	for (nod = 0; nod < 2; ++nod) {
		for (i = 0; i < 3; ++i) {
			zVecL[nod][i] = (*rot[nod])[i][0] * (*elemframe)[2][0]
			                + (*rot[nod])[i][1] * (*elemframe)[2][1]
			                + (*rot[nod])[i][2] * (*elemframe)[2][2];
		}
	}

/* Fitalg 1: Z-axis from node 1 */
	// We are setting fit Alg. to 2 to average z vectors.
	int fitAlg = 2;

	if (fitAlg == 1) {
		t0n[2][0] = zVecL[0][0];
		t0n[2][1] = zVecL[0][1];
		t0n[2][2] = zVecL[0][2];
	}

/* Fitalg .ne. 1: Z-axis as sum of nodal z-axis */
	else {
		t0n[2][0] = zVecL[0][0] + zVecL[1][0];
		t0n[2][1] = zVecL[0][1] + zVecL[1][1];
		t0n[2][2] = zVecL[0][2] + zVecL[1][2];
	}


	t0n[0][0] = xn[1][0] - xn[0][0];
	t0n[0][1] = xn[1][1] - xn[0][1];
	t0n[0][2] = xn[1][2] - xn[0][2];

/* X-axis along element in Cn */
	normalize(t0n[0]);

/* Y-axis as cross product between z and x */
	crossprod(t0n[2], t0n[0], t0n[1]);
	normalize(t0n[1]);

/* Z-axis as cross product between x and y */
	crossprod(t0n[0], t0n[1], t0n[2]);
	normalize(t0n[2]);
}
