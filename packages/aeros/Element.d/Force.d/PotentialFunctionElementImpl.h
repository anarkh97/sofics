#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <iostream>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>
#include <unsupported/Eigen/NumericalDiff>
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>

template<template <typename S> class ScalarValuedFunctionTemplate>
PotentialFunctionElement<ScalarValuedFunctionTemplate>
::PotentialFunctionElement(int _nNodes, DofSet nodalDofs, int* _nn)
		: BoundaryElement(_nNodes, nodalDofs, _nn)
{
}

template<template <typename S> class ScalarValuedFunctionTemplate>
PotentialFunctionElement<ScalarValuedFunctionTemplate>
::PotentialFunctionElement(int _nNodes, DofSet *nodalDofs, int* _nn)
		: BoundaryElement(_nNodes, nodalDofs, _nn)
{
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
PotentialFunctionElement<ScalarValuedFunctionTemplate>
::getInputs(Eigen::Matrix<double,ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            const CoordSet& c0, const GeomState *curState, const GeomState *refState) const
{
	// prepare the constraint function inputs
	if(curState == NULL) { // in this case the function will be evaluated in the undeformed configuration
		q.setZero();
	}
	else {
		for (int k = 0; k < inputs.size(); k++) {
			int i = inputs[k];
			switch(terms[i].dofnum) {
				case 0 :
					q[k] = (*curState)[terms[i].nnum].x - c0[terms[i].nnum]->x;
					break;
				case 1 :
					q[k] = (*curState)[terms[i].nnum].y - c0[terms[i].nnum]->y;
					break;
				case 2 :
					q[k] = (*curState)[terms[i].nnum].z - c0[terms[i].nnum]->z;
					break;
				case 3 : case 4 : case 5 : {
					q[k] = 0;
				} break;
			}
		}
	}
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
PotentialFunctionElement<ScalarValuedFunctionTemplate>::computePressureForce(CoordSet& c0, Vector& F,
                                                                             GeomState *c1, int, double t)
{
	// instantiate the function gradient object
	Eigen::Array<double, ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, nullptr, t);
	Simo::Jacobian<double,ScalarValuedFunctionTemplate> f(sconst,iconst);

	// prepare the function inputs
	const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, c1, NULL);

	// evaluate the function and store values terms
	Eigen::Map<Eigen::Matrix<double,N,1> > fval(F.data());
	fval = f(q,t);
}

template<template <typename S> class ScalarValuedFunctionTemplate>
FullSquareMatrix
PotentialFunctionElement<ScalarValuedFunctionTemplate>::stiffness(const CoordSet& c0, double* karray, int) const
{
	// evaluate the hessian (second partial derivatives w.r.t. the spatial variables)
	FullSquareMatrix K(numDofs(), karray);
	GeomState c1(c0);
	getHessian(NULL, c1, c0, K, 0);

	return K;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
PotentialFunctionElement<ScalarValuedFunctionTemplate>::getInternalForce(GeomState *refState, GeomState& c1, CoordSet& c0, FullSquareMatrix&,
                                                                         double* F, double, double t)
{
	// instantiate the function gradient object
	Eigen::Array<double, ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, &c1, t);
	Simo::Jacobian<double,ScalarValuedFunctionTemplate> f(sconst,iconst);

	// prepare the function inputs
	const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, &c1, refState);

	// evaluate the function and store values terms
	Eigen::Map<Eigen::Matrix<double,N,1> > fval(F);
	fval = f(q,t);
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
PotentialFunctionElement<ScalarValuedFunctionTemplate>::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan,
                                                                         double* F, double dt, double t)

{
	getStiffAndForce((GeomState*)NULL, c1, c0, Ktan, F, dt, t);
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
PotentialFunctionElement<ScalarValuedFunctionTemplate>::getStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan,
                                                                         double* F, double, double t)

{
	// instantiate the function object
	Eigen::Array<double, ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, &c1, t);
	Simo::Jacobian<double,ScalarValuedFunctionTemplate> f(sconst,iconst);

	// prepare the function inputs
	const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, &c1, refState);

	// evaluate the function (debug only)
	//ScalarValuedFunctionTemplate<double> V(sconst,iconst);
	//std::cerr << "V = " << V(q,t) << std::endl;

	// evaluate the gradient of the function and store values terms
	Eigen::Map<Eigen::Matrix<double,N,1> > fval(F);
	fval = f(q,t);

	// evaluate the hessian (second partial derivatives w.r.t. the spatial variables)
	getHessian(refState, c1, c0, Ktan, t);
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
PotentialFunctionElement<ScalarValuedFunctionTemplate>::getHessian(const GeomState *refState, const GeomState &c1,
                                                                   const CoordSet &c0,
                                                                   FullSquareMatrix &B, double t) const
{
	// instantiate the function hessian object
	Eigen::Array<double, ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, &c1, t);
	Simo::Hessian<double,ScalarValuedFunctionTemplate> d2fdq2(sconst,iconst);

	// prepare the function inputs
	const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, &c1, refState);

	// evaluate the hessian
	Eigen::Map<Eigen::Matrix<double,N,N,Eigen::RowMajor> > H(const_cast<double*>(B.data()));
	H = d2fdq2(q, t);
	//std::cerr << "here are the eigenvalues of H: " << H.eigenvalues().transpose() << std::endl;
	//Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,N,N> > es(H);
	//std::cerr << "here are the eigenvalues of H: " << es.eigenvalues().transpose() << std::endl;

#ifdef TOTO
	Eigen::Matrix<double,N,N> H;
  int flag = prop->constraint_hess;
  switch (flag) {
   
    default : case 0 :
      H.setZero();
      break;
    case 1: {
      Simo::Hessian<double,ScalarValuedFunctionTemplate> d2fdq2(sconst,iconst);
      H = d2fdq2(q, t);
    } break;
    case 2: {
      // instantiate the constraint hessian object
      SacadoReverseJacobian<ConstraintJacobian<double,ScalarValuedFunctionTemplate>,true> d2fdq2(dfdq);
      // evaluate the constraint hessian
      d2fdq2(q, H);
    } break;
    case 3: {
      // instantiate the constraint hessian object
      SacadoReverseJacobian<ConstraintJacobian<double,ScalarValuedFunctionTemplate>,false> d2fdq2(dfdq);
      // evaluate the constraint hessian
      d2fdq2(q, H);
    } break;
    case 4 : {
      // instantiate the forward difference object
      Eigen::NumericalDiff<ConstraintJacobian<double,ScalarValuedFunctionTemplate>,Eigen::Forward> fd(dfdq, prop->constraint_hess_eps);
      // evaluate the forward difference approximation to the hessian
      fd.df(q, H);
    } break;
    case 5 : {
      // instantiate the central difference object
      Eigen::NumericalDiff<ConstraintJacobian<double,ScalarValuedFunctionTemplate>,Eigen::Central> cd(dfdq, prop->constraint_hess_eps);
      // evaluate the central difference approximation to the hessian
      cd.df(q, H);
    } break;
  }
#endif
}

#endif
