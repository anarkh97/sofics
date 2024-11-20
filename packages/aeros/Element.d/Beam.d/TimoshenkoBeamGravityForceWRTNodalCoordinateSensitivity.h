#ifdef USE_EIGEN3
#ifndef _TIMOSHENKOBEAMGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_
#define _TIMOSHENKOBEAMGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of gravity force w.r.t the nodal displacements

template<typename Scalar>
class TimoshenkoBeamGravityForceWRTNodalCoordinateSensitivity
		: public VectorValuedFunction<6, 12, Scalar, 14, 1, double> {
public:
	BeamElementTemplate<Scalar> ele;
	Scalar rho, A; // material properties
	Eigen::Array<Scalar, 3, 1> gravityAcceleration;
	Eigen::Array<Scalar, 9, 1> elemframe;
	int gravflg;

public:
	TimoshenkoBeamGravityForceWRTNodalCoordinateSensitivity(const Eigen::Array<double, 14, 1> &sconst,
	                                                        const Eigen::Array<int, 1, 1> &iconst) {
		for (int i = 0; i < 3; ++i) gravityAcceleration[i] = sconst[i];
		rho = sconst[3];
		A = sconst[4];
		elemframe[0] = sconst[5];
		elemframe[1] = sconst[6];
		elemframe[2] = sconst[7];
		elemframe[3] = sconst[8];
		elemframe[4] = sconst[9];
		elemframe[5] = sconst[10];
		elemframe[6] = sconst[11];
		elemframe[7] = sconst[12];
		elemframe[8] = sconst[13];

		gravflg = iconst[0];
	}

	Eigen::Matrix<Scalar, 12, 1> operator()(const Eigen::Matrix<Scalar, 6, 1> &q, Scalar) {
		// inputs:
		// q = Global Displacements at the Nodal Joints

		Eigen::Matrix<Scalar, 2, 1> globalx, globaly, globalz;
		globalx << q[0], q[3];
		globaly << q[1], q[4];
		globalz << q[2], q[5];

		Eigen::Matrix<Scalar, 12, 1> gravityForce;
		ele.buildFrameInTemplate(globalx.data(), globaly.data(), globalz.data(), elemframe.data());
		ele.gForce(globalx.data(), globaly.data(), globalz.data(),
		           gravityAcceleration.data(), elemframe.data(), rho, A, gravityForce.data(), gravflg);

		return gravityForce;
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
