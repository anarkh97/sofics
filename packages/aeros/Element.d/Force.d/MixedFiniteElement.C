#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <iostream>
#include <limits>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/MooneyRivlinMat.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Math.d/FullSquareMatrix.h>
#include <unsupported/Eigen/NumericalDiff>

template<template <typename S> class ScalarValuedFunctionTemplate>
MixedFiniteElement<ScalarValuedFunctionTemplate>
::MixedFiniteElement(int _nNodes, DofSet nodalDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalDofs, _nn)
{
  nIV = ScalarValuedFunctionTemplate<double>::NumberOfNodes2+ScalarValuedFunctionTemplate<double>::NumberOfNodes3;
  nLM = ScalarValuedFunctionTemplate<double>::NumberOfNodes4;
  first_time = true;
  materialType = 1;
  mat = NULL;
  epsilon = 1e7;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
MixedFiniteElement<ScalarValuedFunctionTemplate>
::MixedFiniteElement(int _nNodes, DofSet *nodalDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalDofs, _nn)
{
  nIV = ScalarValuedFunctionTemplate<double>::NumberOfNodes2+ScalarValuedFunctionTemplate<double>::NumberOfNodes3;
  nLM = ScalarValuedFunctionTemplate<double>::NumberOfNodes4;
  first_time = true;
  materialType = 1;
  mat = NULL;
  epsilon = 1e7;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
MixedFiniteElement<ScalarValuedFunctionTemplate>
::setMaterial(NLMaterial *_mat)
{
  if(StVenantKirchhoffMat *mat = dynamic_cast<StVenantKirchhoffMat*>(_mat)) {
    materialType = 1;
  }
  else if(HenckyMat *mat = dynamic_cast<HenckyMat*>(_mat)) {
    materialType = 2;
  }
  else if(ElaLinIsoMat *mat = dynamic_cast<ElaLinIsoMat*>(_mat)) {
    materialType = 0;
  }
  else if(NeoHookeanMat *mat = dynamic_cast<NeoHookeanMat*>(_mat)) {
    materialType = 3;
  }
  else if(MooneyRivlinMat *mat = dynamic_cast<MooneyRivlinMat*>(_mat)) {
    materialType = 4;
  }
  else if(OgdenMat *mat = dynamic_cast<OgdenMat*>(_mat)) {
    materialType = 5;
  }
  else {
    std::cerr << " *** ERROR: Unsupported MATLAW for mixed finite element.\n";
    materialType = -1;
  }

  mat = _mat;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
MixedFiniteElement<ScalarValuedFunctionTemplate>
::getConstants(const CoordSet &cs, Eigen::Array<typename ScalarValuedFunctionTemplate<double>::ScalarConstantType,
               ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
               Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
               const GeomState* gs) const
{
  // nodal coordinates in reference configuration
  const int nNode = ScalarValuedFunctionTemplate<double>::NumberOfNodes;
  const int nDime = ScalarValuedFunctionTemplate<double>::NumberOfDimensions;
  for(int i=0; i<nNode; ++i) {
    sconst[nDime*i  ] = cs[nn[i]]->x;
    sconst[nDime*i+1] = cs[nn[i]]->y;
    if(nDime == 2) continue;
    sconst[nDime*i+2] = cs[nn[i]]->z;
  }

  // material constants
  const int nMatc = ScalarValuedFunctionTemplate<double>::MaximumNumberOfMaterialConstants;
  if(mat) {
    std::vector<double> v(0);
    mat->getMaterialConstants(v);
    for(int i=0; i<nMatc; ++i) sconst[nDime*nNode+i] = (i < v.size()) ? v[i] : 0;
  }
  else {
    const double &E = prop->E;
    const double &nu = prop->nu;
    sconst[nDime*nNode  ] = E*nu/((1+nu)*(1-2*nu)); // Lamé's first parameter
    sconst[nDime*nNode+1] = E/(2*(1+nu));           // Lamé's second parameter (shear modulus)
    for(int i=2; i<nMatc; ++i) sconst[nDime*nNode+i] = 0;
  }

  // reference state for dilatation
  const int nNode2 = ScalarValuedFunctionTemplate<double>::NumberOfNodes2;
  sconst.template segment<nNode2>(nDime*nNode+nMatc).setOnes();

  // reference state for pressure
  const int nNode3 = ScalarValuedFunctionTemplate<double>::NumberOfNodes3;
  sconst.template segment<nNode3>(nDime*nNode+nMatc+nNode2).setZero();

  // Lagrange multipliers and penalty parameter
  if(nLM > 0) {
    double *state = (gs) ? gs->getElemState(this->getGlNum()) : 0;
    for(int i=0; i<nLM; ++i) sconst[nDime*nNode+nMatc+nNode2+nNode3+i] = (state) ? state[nIV+i] : 0;
    sconst[nDime*nNode+nMatc+nNode2+nNode3+nLM] = epsilon;
  }
  
  iconst[0] = getQuadratureOrder();
  iconst[1] = (gs) ? materialType : 0;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
MixedFiniteElement<ScalarValuedFunctionTemplate>
::getInputs(Eigen::Matrix<double,ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            const CoordSet& c0, const GeomState *curState) const
{
  // prepare the function inputs, which is the increment in the nodal values of the three fields
  // w.r.t. the reference state.
  // for the displacement field the reference state is the nodal coordinates (stored in c0).
  // for the dilitation field the reference state is 1, and for the pressure field the reference state is 0.

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
    double *state = curState->getElemState(this->getGlNum());
    for(int i=0; i<nIV; ++i) q[inputs.size()+i] = (state) ? state[i] : 0;
  }
}

template<template <typename S> class ScalarValuedFunctionTemplate>
FullSquareMatrix
MixedFiniteElement<ScalarValuedFunctionTemplate>
::stiffness(const CoordSet& c0, double* karray, int) const
{
  // evaluate the Hessian (second partial derivatives w.r.t. the spatial variables)
  FullSquareMatrix K(numDofs(), karray);

  const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,N> H;
  getHessian(NULL, NULL, c0, H, 0);

  const int N1 = N-nIV;
  const int N2 = nIV;
  // [ A   B ] [ x ] = [ f ]
  // [ B^T C ] [ y ]   [ g ]
  //  y = C^{-1} ( g - B^Tx ) --> (A - BC^{-1}B^T)x = f - C^{-1}g
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > S(karray,N1,N1);
  S = H.topLeftCorner(N1,N1);
  if(N2 > 0) {
    S -= H.topRightCorner(N1,N2)*H.bottomRightCorner(N2,N2).inverse()*H.bottomLeftCorner(N2,N1);
  }

  return K;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
MixedFiniteElement<ScalarValuedFunctionTemplate>
::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* F, double dt, double t)
{
  getStiffAndForce((GeomState*)NULL, c1, c0, Ktan, F, dt, t);
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void 
MixedFiniteElement<ScalarValuedFunctionTemplate>
::getStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* F, double, double t)
{
  // instantiate the function object
  Eigen::Array<double, ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, &c1);
  Simo::Jacobian<double,ScalarValuedFunctionTemplate> foo(sconst,iconst);

  // prepare the function inputs
  const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1);

  const int N1 = N-nIV;
  const int N2 = nIV;

  // update the internal variables (note: these should probably be states associated with internal nodes rather than
  // element states to distinguish them with history variables used for elasto-plasticity for instance).
  if(!first_time && N2 > 0) {
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > y(c1.getElemState(this->getGlNum()),N2);
    y -= Cinvg+CinvBt*(q.head(N1)-q_copy.head(N1));
    q.tail(N2) = y;
  }
  q_copy = q;

  // evaluate the gradient of the function
  Eigen::Matrix<double,N,1> G;
  G = foo(q,t);

  // evaluate the hessian of the function
  Eigen::Matrix<double,N,N> H;
  getHessian(refState, &c1, c0, H, t);

  // [ A   B ] [ x ] = [ -f ]
  // [ B^T C ] [ y ]   [ -g ]
  //  y = C^{-1} ( -g - B^Tx ) --> (A - BC^{-1}B^T)x = -f + BC^{-1}g
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > S(Ktan.data(),N1,N1);
  S = H.topLeftCorner(N1,N1);

  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > f(F,N1);
  f = G.head(N1);

  if(N2 > 0) {
#if EIGEN_VERSION_AT_LEAST(3,2,7)
    Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::NoQRPreconditioner> dec(H.bottomRightCorner(N2,N2), Eigen::ComputeThinU | Eigen::ComputeThinV);
    dec.setThreshold(10*std::numeric_limits<double>::epsilon()/dec.singularValues()[0]);
#else
    // XXX FullPivLU doesn't work since Eigen 3.2.7
    Eigen::FullPivLU<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dec(H.bottomRightCorner(N2,N2));
#endif

    CinvBt = dec.solve(H.bottomLeftCorner(N2,N1));
    Cinvg  = dec.solve(G.tail(N2));

    S -= H.topRightCorner(N1,N2)*CinvBt;
    f -= H.topRightCorner(N1,N2)*Cinvg;

    first_time = false;
  }
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void 
MixedFiniteElement<ScalarValuedFunctionTemplate>
::getHessian(GeomState *refState, GeomState *c1, CoordSet& c0, Eigen::Matrix<double,ScalarValuedFunctionTemplate<double>::
             NumberOfGeneralizedCoordinates, ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates>& H, double t) 
{
  // instantiate the function hessian object
  Eigen::Array<double, ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, c1);
  Simo::Hessian<double,ScalarValuedFunctionTemplate> d2fdq2(sconst,iconst);

  // prepare the function inputs
  const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1);

  // evaluate the hessian 
  H = d2fdq2(q, t);
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
MixedFiniteElement<ScalarValuedFunctionTemplate>
::initMultipliers(GeomState& c1)
{
  if(nLM > 0) {
    double *state = c1.getElemState(this->getGlNum());
    for(int i=0; i<nLM; ++i) state[nIV+i] = 0;
  }
}

template<template <typename S> class ScalarValuedFunctionTemplate>
double
MixedFiniteElement<ScalarValuedFunctionTemplate>
::getError(GeomState& c1)
{
  // update the internal variables
  // it is better to do it here rather than in updateStates because for augmented Lagrangian, updateStates is
  // only called at the end of the penalty iteration loop, while getError is called at each iteration.
  const int N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, *c1.getCoordSet(), &c1);

  const int N1 = N-nIV;
  const int N2 = nIV;
  if(!first_time && N2 > 0) {
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > y(c1.getElemState(this->getGlNum()),N2);
    y -= Cinvg+CinvBt*(q.head(N1)-q_copy.head(N1));
    first_time = true;
  }

  if(nLM > 0) {
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > h(c1.getElemState(this->getGlNum()),ScalarValuedFunctionTemplate<double>::NumberOfNodes2);
    return h.template lpNorm<1>();
  }
  else return 0;
}

template<template <typename S> class ScalarValuedFunctionTemplate>
void
MixedFiniteElement<ScalarValuedFunctionTemplate>
::updateMultipliers(GeomState& c1)
{
  if(nLM > 0) {
    // XXX assuming that the same shape functions are used for Theta and lambda
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > h(c1.getElemState(this->getGlNum()),ScalarValuedFunctionTemplate<double>::NumberOfNodes2);
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > lambda(c1.getElemState(this->getGlNum())+nIV,nLM);
    lambda += epsilon*h;
  }
}

#endif
