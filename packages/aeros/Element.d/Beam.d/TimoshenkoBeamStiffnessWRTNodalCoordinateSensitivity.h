#ifdef USE_EIGEN3
#ifndef _TIMOSHENKOBEAMSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_
#define _TIMOSHENKOBEAMSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the stiffness w.r.t the nodal coordinates

template<typename Scalar>
class TimoshenkoBeamStiffnessWRTNodalCoordinateSensitivity
		: public MatrixValuedFunction<6, 12, 12, Scalar, 18, 0, double> {
public:
	BeamElementTemplate<Scalar> ele;
	Eigen::Array<Scalar, 9, 1> elemframe;
	Scalar A, E, Ixx, Iyy, Izz, alphaY, alphaZ, c1, nu; // material properties

public:
	TimoshenkoBeamStiffnessWRTNodalCoordinateSensitivity(const Eigen::Array<double, 18, 1> &sconst,
	                                                     const Eigen::Array<int, 0, 1> &iconst) {
		E = sconst[0];
		A = sconst[1];
		Ixx = sconst[2];
		Iyy = sconst[3];
		Izz = sconst[4];
		alphaY = sconst[5];
		alphaZ = sconst[6];
		c1 = sconst[7];
		nu = sconst[8];
		elemframe[0] = sconst[9];
		elemframe[1] = sconst[10];
		elemframe[2] = sconst[11];
		elemframe[3] = sconst[12];
		elemframe[4] = sconst[13];
		elemframe[5] = sconst[14];
		elemframe[6] = sconst[15];
		elemframe[7] = sconst[16];
		elemframe[8] = sconst[17];
	}

	Eigen::Matrix<Scalar, 12, 12> operator()(const Eigen::Matrix<Scalar, 6, 1> &q, Scalar) {
		// inputs:
		// q = Nodal Coordinates

		Eigen::Matrix<Scalar, 2, 1> globalx, globaly, globalz;
		globalx << q[0], q[3];
		globaly << q[1], q[4];
		globalz << q[2], q[5];
/*
      std::cerr << "Area = " << A << std::endl;
      std::cerr << "E = " << E << std::endl;
      std::cerr << "elemframe = " << elemframe[0] << " " << elemframe[1] << " " << elemframe[2] << std::endl;
      std::cerr << "            " << elemframe[3] << " " << elemframe[4] << " " << elemframe[5] << std::endl;
      std::cerr << "            " << elemframe[6] << " " << elemframe[7] << " " << elemframe[8] << std::endl;
      std::cerr << "Ixx = " << Ixx << std::endl;
      std::cerr << "Iyy = " << Iyy << std::endl;
      std::cerr << "Izz = " << Izz << std::endl;
      std::cerr << "alphaY = " << alphaY << std::endl;
      std::cerr << "alphaZ = " << alphaZ << std::endl;
      std::cerr << "C1 = " << c1 << std::endl;
      std::cerr << "nu = " << nu << std::endl;
      std::cerr << "x = " << globalx[0] << " " << globalx[1]  << std::endl;
      std::cerr << "y = " << globaly[0] << " " << globaly[1]  << std::endl;
      std::cerr << "z = " << globalz[0] << " " << globalz[1]  << std::endl;
*/
		Eigen::Matrix<Scalar, 12, 12> estiff;
		ele.buildFrameInTemplate(globalx.data(), globaly.data(), globalz.data(), elemframe.data());
		ele.modmstif7(estiff.data(), A, E, elemframe.data(),
		              Ixx, Iyy, Izz, alphaY, alphaZ, c1,
		              nu, globalx.data(), globaly.data(), globalz.data(), 1);

		// return value:

		return estiff;
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
