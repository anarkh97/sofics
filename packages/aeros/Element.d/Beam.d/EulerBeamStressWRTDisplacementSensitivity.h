#ifdef USE_EIGEN3
#ifndef _EULERBEAMSTRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _EULERBEAMSTRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class EulerBeamStressWRTDisplacementSensitivity : public VectorValuedFunction<12, 2, Scalar, 25, 1, double> {
public:
	BeamElementTemplate<Scalar> ele;
	Eigen::Array<Scalar, 2, 1> globalx, globaly, globalz; // nodal coordinates
	Scalar A, E, Ixx, Iyy, Izz, c, nu, W, Ta; // material properties
	Eigen::Array<Scalar, 9, 1> elemframe;
	Eigen::Array<Scalar, 2, 1> ndTemps;

public:
	EulerBeamStressWRTDisplacementSensitivity(const Eigen::Array<double, 25, 1> &sconst,
	                                          const Eigen::Array<int, 1, 1> &iconst) {
		globalx = sconst.segment<2>(0).cast<Scalar>();
		globaly = sconst.segment<2>(2).cast<Scalar>();
		globalz = sconst.segment<2>(4).cast<Scalar>();
		A = sconst[6];
		E = sconst[7];
		for (int i = 0; i < 9; ++i) {
			elemframe[i] = sconst[8 + i];
		}
		Ixx = sconst[17];
		Iyy = sconst[18];
		Izz = sconst[19];
		nu = sconst[20];
		W = sconst[21];
		Ta = sconst[22];
		ndTemps[0] = sconst[23];
		ndTemps[1] = sconst[24];
	}

	Eigen::Matrix<Scalar, 2, 1> operator()(const Eigen::Matrix<Scalar, 12, 1> &q, Scalar) {
		// inputs:
		// q = Global Displacements at the Nodal Joints

		Eigen::Matrix<Scalar, 12, 1> globalu = q;

		Eigen::Array<Scalar, 7, 2> elStress;
		int numel = 1;
		int elm = 1;
		int maxgus = 2;
		int maxstr = 7;
		int maxsze = 1;

		ele.sands6(A, E, elm,
		           elStress.data(),
		           maxsze, maxgus, maxstr,
		           elemframe.data(),
		           Ixx, Iyy, Izz, nu,
		           globalx.data(), globaly.data(), globalz.data(), globalu.data(),
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
