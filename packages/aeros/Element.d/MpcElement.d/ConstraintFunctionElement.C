#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/TimeDerivatives.h>
#include <iostream>
#include <Element.d/Function.d/SacadoReverseJacobian.h>
#include <Math.d/FullSquareMatrix.h>
#include <unsupported/Eigen/NumericalDiff>

template<template <typename S> class ConstraintFunctionTemplate>
ConstraintFunctionElement<ConstraintFunctionTemplate>
::ConstraintFunctionElement(int _nNodes, DofSet nodalDofs, int* _nn, int _type)
		: MpcElement(_nNodes, nodalDofs, _nn)
{
	type = _type;
}

template<template <typename S> class ConstraintFunctionTemplate>
ConstraintFunctionElement<ConstraintFunctionTemplate>
::ConstraintFunctionElement(int _nNodes, DofSet *nodalDofs, int* _nn, int _type)
		: MpcElement(_nNodes, nodalDofs, _nn)
{
	type = _type;
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>
::getInputs(Eigen::Matrix<double,ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            const CoordSet& c0, const GeomState *curState, const GeomState *refState) const
{
	// prepare the constraint function inputs
	if(curState == NULL) {
		// in this case the function will be evaluated in the undeformed configuration
		q.setZero();
	}
	else {
		// eulerian description of rotations
		int k = 0;
		for (int i = 0; i < nterms; i++) {
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
					q[k] = 0; // spin
				} break;
			}
			k++;
		}
	}
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::buildFrame(CoordSet& _c0)
{
	c0 = &_c0;
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::setProp(StructProp *p, bool _myProp)
{
	Element::setProp(p, _myProp);

	// instantiate the constraint function object
	Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(*c0, sconst, iconst);
	ConstraintFunctionTemplate<double> f(sconst,iconst);

	// prepare the constraint function inputs
	const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, *c0, NULL, NULL);
	double t = 0;

	// evaluate the constraint function and store -ve value in LMPCons::rhs
	original_rhs.r_value = rhs.r_value = -f(q,t);

	// instantiate the constraint jacobian object
	Simo::Jacobian<double,ConstraintFunctionTemplate> dfdq(sconst,iconst);

	// evaluate the constraint jacobian (partial derivatives w.r.t. the spatial variables)
	// and store coefficients in LMPCons::terms array
	Eigen::Matrix<double,1,N> J;
	J = dfdq(q, t);
	for(int i = 0; i < nterms; ++i) terms[i].coef.r_value = J[i];
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::update(GeomState* refState, GeomState& c1, CoordSet& c0, double t)
{
	// instantiate the constraint function object
	Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, &c1);
	ConstraintFunctionTemplate<double> f(sconst,iconst);

	// prepare the constraint function inputs
	const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, &c1, refState);

	// evaluate the constraint function and store -ve value in LMPCons::rhs
	rhs.r_value = -f(q,t);

	// evaluate the constraint jacobian (partial derivatives w.r.t. the spatial variables)
	// and store coefficients in LMPCons::terms array
	Eigen::Matrix<double,1,N> J;

	Simo::Jacobian<double,ConstraintFunctionTemplate> dfdq(sconst,iconst);
	J = dfdq(q, t);

	for(int i = 0; i < nterms; ++i) terms[i].coef.r_value = J[i];
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::getHessian(const GeomState *refState, const GeomState& c1,
                                                                  const CoordSet& c0,
                                                                  FullSquareMatrix& B, double t) const
{
	// instantiate the constraint function object
	Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, &c1);
	ConstraintFunctionTemplate<double> f(sconst,iconst);

	// prepare the function inputs
	const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, &c1, refState);

	Eigen::Matrix<double,N,N> H;
	switch (prop->constraint_hess) {

		default : case 0 :
			H.setZero();
			break;
		case 1: {
			// evaluate the constraint hessian using the default implementation
			Simo::Hessian<double,ConstraintFunctionTemplate> d2fdq2(sconst,iconst);
			H = d2fdq2(q, t);
		} break;
#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
		#ifdef USE_SACADO
    /*case 2: {
      // evaluate the constraint hessian by forward automatic differentation of the jacobian
      Simo::SpatialView<double,ConstraintJacobian> dfdq(sconst,iconst,t);
      SacadoForwardJacobian<Simo::SpatialView<double,ConstraintJacobian>> d2fdq2(dfdq);
      d2fdq2(q, H);
    } break;*/
    case 3: {
      // evaluate the constraint hessian by reverse automatic differentation of the jacobian
      Simo::SpatialView<double,ConstraintJacobian> dfdq(sconst,iconst,t);
      SacadoReverseJacobian<Simo::SpatialView<double,ConstraintJacobian>> d2fdq2(dfdq);
      d2fdq2(q, H);
    } break;
#endif
    case 4 : {
      // evaluate the forward difference approximation of the constraint hessian 
      Simo::SpatialView<double,ConstraintJacobian> dfdq(sconst,iconst,t);
      Eigen::NumericalDiff<Simo::SpatialView<double,ConstraintJacobian>,Eigen::Forward> fd(dfdq, prop->constraint_hess_eps);
      fd.df(q, H);
    } break;
    case 5 : {
      // evaluate the central difference approximation of the constraint hessian
      Simo::SpatialView<double,ConstraintJacobian> dfdq(sconst,iconst,t);
      Eigen::NumericalDiff<Simo::SpatialView<double,ConstraintJacobian>,Eigen::Central> cd(dfdq, prop->constraint_hess_eps);
      cd.df(q, H);
    } break;
#endif
	}

	for(int i = 0; i < nterms; ++i)
		for(int j = 0; j < nterms; ++j)
			B[i][j] = H(i,j);
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::computePressureForce(CoordSet& c0, Vector& elPressureForce,
                                                                            GeomState *c1, int cflg, double t)
{
	// instantiate the constraint function object
	Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst);
	ConstraintFunctionTemplate<double> f(sconst,iconst);

	// prepare the constraint function inputs
	const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, NULL, NULL);

	// evaluate the function
	rhs.r_value = -f(q,t);

	MpcElement::computePressureForce(c0, elPressureForce, c1, cflg, t);
}

template<template <typename S> class ConstraintFunctionTemplate>
double
ConstraintFunctionElement<ConstraintFunctionTemplate>::getVelocityConstraintRhs(GeomState *refState, GeomState& c1,
                                                                                CoordSet& c0, double t)
{
	// instantiate the constraint function object
	Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
	Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
	getConstants(c0, sconst, iconst, &c1);
	ConstraintFunctionTemplate<double> f(sconst,iconst);

	// prepare the constraint function inputs
	const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
	Eigen::Matrix<double,N,1> q;
	getInputs(q, c0, &c1, refState);

	// evaluate the first partial time derivative of the constraint function
	Simo::FirstPartialTimeDerivative<double,ConstraintFunctionTemplate> dfdt(sconst,iconst);
	return -dfdt(q,t);
}

template<template <typename S> class ConstraintFunctionTemplate>
double
ConstraintFunctionElement<ConstraintFunctionTemplate>::getAccelerationConstraintRhs(GeomState *refState, GeomState& c1,
                                                                                    CoordSet& c0, double t)
{
#if defined(USE_SACADO) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
	// instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, &c1);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1, refState);

  // evaluate the second partial time derivative of the constraint function object
  Simo::SecondPartialTimeDerivative<double,ConstraintFunctionTemplate> d2fdt2(sconst,iconst);
  double y = d2fdt2(q,t);

  // evaluate the first partial time derivative of the constraint jacobian object
  Simo::FirstPartialTimeDerivative<double,ConstraintJacobian,1> dJdt(sconst,iconst);
  Eigen::Matrix<double,1,N> z = dJdt(q,t);

  // evaluate the first partial space derivatives of the constraint jacobian velocity product
  Eigen::Matrix<double,N,1> v;
  for(int i = 0; i < nterms; ++i) {
    if(terms[i].dofnum == 3 || terms[i].dofnum == 4 || terms[i].dofnum == 5) {
      // compute spatial angular velocity
      Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&c1[terms[i].nnum].R[0][0]);
      Eigen::Map<Eigen::Vector3d> Omega(&c1[terms[i].nnum].v[3]);
      Eigen::Vector3d omega = R*Omega;
      v[i] = omega[terms[i].dofnum-3];
    }
    else {
      v[i] = c1[terms[i].nnum].v[terms[i].dofnum];
    }
  }

  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants+N, 1> sconst2;
  sconst2 << sconst, v.array();
  Simo::Jacobian<double,ConstraintJacobianVectorProduct,1> DJV(sconst2,iconst);
  Eigen::Matrix<double,1,N> w = DJV(q,t);

  return -w.dot(v) - 2*z.dot(v) - y;
#else
	return 0;
#endif
}
#endif
