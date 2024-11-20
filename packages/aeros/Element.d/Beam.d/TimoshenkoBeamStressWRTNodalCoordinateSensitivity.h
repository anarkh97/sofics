#ifdef USE_EIGEN3
#ifndef _TIMOSHENKOBEAMSTRESSWRTNODALCOORDINATESENSITIVITY_H_
#define _TIMOSHENKOBEAMSTRESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class TimoshenkoBeamStressWRTNodalCoordinateSensitivity : public VectorValuedFunction<6, 2, Scalar, 34, 1, double> {
public:
	BeamElementTemplate<Scalar> ele;
	Eigen::Array<Scalar, 12, 1> globalu; // displacement
	Scalar A, E, Ixx, Iyy, Izz, alphaY, alphaZ, c, nu, W, Ta; // material properties
	Eigen::Array<Scalar, 9, 1> elemframe;
	Eigen::Array<Scalar, 2, 1> ndTemps;

public:
	TimoshenkoBeamStressWRTNodalCoordinateSensitivity(const Eigen::Array<double, 34, 1> &sconst,
	                                                  const Eigen::Array<int, 1, 1> &iconst) {
		globalu = sconst.segment<12>(0).cast<Scalar>();
		A = sconst[12];
		E = sconst[13];
		for (int i = 0; i < 9; ++i) {
			elemframe[i] = sconst[14 + i];
		}
		Ixx = sconst[23];
		Iyy = sconst[24];
		Izz = sconst[25];
		alphaY = sconst[26];
		alphaZ = sconst[27];
		c = sconst[28];
		nu = sconst[29];
		W = sconst[30];
		Ta = sconst[31];
		ndTemps[0] = sconst[32];
		ndTemps[1] = sconst[33];
	}

	Eigen::Matrix<Scalar, 2, 1> operator()(const Eigen::Matrix<Scalar, 6, 1> &q, Scalar) {
		// inputs:
		// q = Global Displacements at the Nodal Joints

		Eigen::Matrix<Scalar, 2, 1> globalx, globaly, globalz;
		globalx << q[0], q[3];
		globaly << q[1], q[4];
		globalz << q[2], q[5];

		Eigen::Array<Scalar, 7, 2> elStress;
		int numel = 1;
		int elm = 1;
		int maxgus = 2;
		int maxstr = 7;
		int maxsze = 1;
		ele.buildFrameInTemplate(globalx.data(), globaly.data(), globalz.data(), elemframe.data());
		ele.sands7(elm, A, E, elemframe.data(),
		           Ixx, Iyy, Izz, alphaY,
		           alphaZ, c, nu, globalx.data(), globaly.data(), globalz.data(), globalu.data(),
		           elStress.data(), numel, maxgus, maxstr, maxsze,
		           W, Ta, ndTemps.data());

		// return value:

		Eigen::Matrix<Scalar, 2, 1> v;
		v[0] = elStress(6, 0);
		v[1] = elStress(6, 1);

		return v;
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
