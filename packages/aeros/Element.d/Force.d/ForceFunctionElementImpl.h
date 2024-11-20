#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>

template<template <typename S> class VectorValuedFunctionTemplate>
ForceFunctionElement<VectorValuedFunctionTemplate>
::ForceFunctionElement(int _nNodes, DofSet nodalDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalDofs, _nn)
{
}

template<template <typename S> class VectorValuedFunctionTemplate>
ForceFunctionElement<VectorValuedFunctionTemplate>
::ForceFunctionElement(int _nNodes, DofSet *nodalDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalDofs, _nn)
{
}

template<template <typename S> class VectorValuedFunctionTemplate>
ForceFunctionElement<VectorValuedFunctionTemplate>
::ForceFunctionElement(int _nNodes, DofSet *nodalInputDofs, DofSet *nodalOutputDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalInputDofs, nodalOutputDofs, _nn)
{
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>
::getInputs(Eigen::Matrix<double,VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            const CoordSet& c0, const GeomState *curState, const GeomState *refState) const
{
  // prepare the constraint function inputs
  if(curState == NULL) {
    // in this case the function will be evaluated in the undeformed configuration
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

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>::computePressureForce(CoordSet& c0, Vector& F,
                                                                         GeomState *c1, int, double t)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, nullptr, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1, NULL);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);

  // fill F
  for(int i=0; i<numDofs(); ++i) F[i] = 0;

  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    F[k] = fval[i];
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
FullSquareMatrix
ForceFunctionElement<VectorValuedFunctionTemplate>::stiffness(const CoordSet& c0, double* karray, int) const
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q; // = Eigen::Matrix<double,N,1>::Zero();
  getInputs(q, c0, NULL, NULL);
  double t = 0;

  // evaluate the jacobian (partial derivatives w.r.t. the spatial variables)
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  FullM J(M,N);
  GeomState c1(c0);
  getJacobian(NULL, c1, c0, J, t);

  // fill K
  FullSquareMatrix K(numDofs(), karray);
  K.zero();
  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    for(int j=0; j<N; ++j) {
      int l = inputs[j];
      K[k][l] = J[i][j];
    }
  }
  return K;
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>::getInternalForce(GeomState *refState, GeomState& c1, CoordSet& c0,
                                                                     FullSquareMatrix&, double* F, double, double t)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, &c1, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q; // = Eigen::Matrix<double,N,1>::Zero();
  getInputs(q, c0, &c1, refState);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);

  // fill F
  for(int i=0; i<numDofs(); ++i) F[i] = 0;

  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    F[k] = fval[i];
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan,
                                                                     double* F, double dt, double t)

{
  getStiffAndForce((GeomState*)NULL, c1, c0, Ktan, F, dt, t);
}

template<template <typename S> class VectorValuedFunctionTemplate>
void 
ForceFunctionElement<VectorValuedFunctionTemplate>::getStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0,
                                                                     FullSquareMatrix& Ktan, double* F, double, double t)

{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, &c1, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1, refState);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);

  // evaluate the jacobian (partial derivatives w.r.t. the spatial variables)
  FullM J(M,N);
  getJacobian(refState, c1, c0, J, t);

  // fill Ktan and F
  Ktan.zero();
  for(int i=0; i<numDofs(); ++i) F[i] = 0;

  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    F[k] = fval[i];
    for(int j=0; j<N; ++j) {
      int l = inputs[j];
      Ktan[k][l] = J[i][j];
    }
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void 
ForceFunctionElement<VectorValuedFunctionTemplate>::getJacobian(const GeomState *refState, const GeomState &c1, const CoordSet& c0, FullM& B, double t) const
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, &c1, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1, refState);

  // instantiate the jacobian object
  Simo::Jacobian<double,VectorValuedFunctionTemplate> dfdq(sconst,iconst);

  // evaluate the jacobian
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Map<Eigen::Matrix<double,M,N,Eigen::RowMajor> > J(B.data());
  J = dfdq(q, t);
}
#endif
