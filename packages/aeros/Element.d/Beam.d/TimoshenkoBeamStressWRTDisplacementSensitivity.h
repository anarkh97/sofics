#ifdef USE_EIGEN3
#ifndef _TIMOSHENKOBEAMSTRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _TIMOSHENKOBEAMSTRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class TimoshenkoBeamStressWRTDisplacementSensitivity : public VectorValuedFunction<12, 2, Scalar, 28, 1, double> {
public:
	BeamElementTemplate<Scalar> ele;
	Eigen::Array<Scalar, 2, 1> globalx, globaly, globalz; // nodal coordinates
	Scalar A, E, Ixx, Iyy, Izz, alphaY, alphaZ, c, nu, W, Ta; // material properties
	Eigen::Array<Scalar, 9, 1> elemframe;
	Eigen::Array<Scalar, 2, 1> ndTemps;

public:
	TimoshenkoBeamStressWRTDisplacementSensitivity(const Eigen::Array<double, 28, 1> &sconst,
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
		alphaY = sconst[20];
		alphaZ = sconst[21];
		c = sconst[22];
		nu = sconst[23];
		W = sconst[24];
		Ta = sconst[25];
		ndTemps[0] = sconst[26];
		ndTemps[1] = sconst[27];

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
