//
// Created by Michel Lesoinne on 1/11/18.
//

#ifndef FEM_OPERATORS_H
#define FEM_OPERATORS_H

#include <memory>
#include <Eigen/Dense>

class ScatterMatrix {

};

/** Feti operators */
template <typename Scalar>
class Operators {
	/// \brief The matrix of 1 and 0 for regular Lagrange multipliers.
	ScatterMatrix B;
	/// \brief The matrix for wet interface
	ScatterMatrix Bwet;

	///
	ScatterMatrix Bctilde;
//	std::unique_ptr<GenCuCSparse<Scalar>> Krw;

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Ave;
};


#endif //FEM_OPERATORS_H
