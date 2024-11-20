#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Utils.d/dofset.h>

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes>
void
AutoShapeFunction<ShapeFunctionTemplate,NumberOfNodes>::getLocalDerivatives(Tensor *_localDerivatives, double _xi[3])
{
#ifdef USE_EIGEN3
  Eigen::Matrix<double,3,1> xi;
  Eigen::Matrix<double,NumberOfNodes,1> N;      // shape function values
  Eigen::Matrix<double,NumberOfNodes,3> dNdxi;  // derivative of shape functions w.r.t. xi
  Simo::Jacobian<double, ShapeFunctionTemplate> dSF(Eigen::Array<double,0,1>::Zero(), Eigen::Array<int,0,1>::Zero());

  xi << _xi[0], _xi[1], _xi[2];

  dNdxi = dSF(xi, 0.);

  Tensor_d1s2_sparse *localDerivatives  = static_cast<Tensor_d1s2_sparse *>(_localDerivatives);
  for(int k = 0; k < NumberOfNodes; ++k)
    for(int i = 0; i < 3; ++i) {
      (*localDerivatives)[k][3*i  ] = dNdxi(k,0);
      (*localDerivatives)[k][3*i+1] = dNdxi(k,1);
      (*localDerivatives)[k][3*i+2] = dNdxi(k,2);
    }
#else
  std::cerr << " *** ERROR: Implementation of AutoShapeFunction::getLocalDerivatives in file" << std::endl
            << "            Element.d/NonLinearity.d/SolidElementTemplate.C requires AERO-S configured with Eigen library. Exiting...\n";
  exit(-1);
#endif
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes>
double
AutoShapeFunction<ShapeFunctionTemplate,NumberOfNodes>::interpolateScalar(double *_q, double _xi[3])
{
#ifdef USE_EIGEN3
  Eigen::Matrix<double,3,1> xi;
  ShapeFunctionTemplate<double> N(Eigen::Array<double,0,1>::Zero(), Eigen::Array<int,0,1>::Zero());
  Eigen::Map<Eigen::Matrix<double,NumberOfNodes,1> > q(_q);

  xi << _xi[0], _xi[1], _xi[2];
  return N(xi, 0.).dot(q);
#else
  std::cerr << " *** ERROR: Implementation of AutoShapeFunction::interpolate in file" << std::endl
            << "            Element.d/NonLinearity.d/SolidElementTemplate.C requires AERO-S configured with Eigen library. Exiting...\n";
  exit(-1);
  return 0;
#endif
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::SolidElementTemplate(int *nd) 
 : material(0)
{
  for(int i = 0; i < NumberOfNodes; ++i)
    n[i] = nd[i];
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::setMaterial(NLMaterial *m)
{
  material = m;
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
int
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::getNumGaussPoints() const
{
  return NumIntgPts;
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
const AutoShapeFunction<ShapeFunctionTemplate,NumberOfNodes>
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::shapeFunction;

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
ShapeFunction *
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::getShapeFunction() const
{
  return const_cast<AutoShapeFunction<ShapeFunctionTemplate,NumberOfNodes>* >(&shapeFunction);
}

extern LinearStrain linearStrain;
extern GreenLagrangeStrain greenLagrangeStrain;
extern DeformationGradient deformationGradient;

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
StrainEvaluator *
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::getStrainEvaluator() const
{
  return material->getStrainEvaluator();
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
NLMaterial *
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::getMaterial() const
{
  return material;
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::getNodeRefCoords(double (*_nodeRefCoords)[3])
{
  for(int i = 0; i < NumberOfNodes; ++i)
    for(int j = 0; j < 3; ++j)
      _nodeRefCoords[i][j] = nodeRefCoords[i][j];
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::getLocalNodalCoords(int i, double *coords)
{
  for(int j = 0; j < 3; ++j)
    coords[j] = nodeRefCoords[i][j];
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
int
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::numNodes() const
{
  return NumberOfNodes;
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
int 
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::numDofs() const
{
  return NumberOfNodes*3;
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::renum(const int *table)
{
  for(int i = 0; i < NumberOfNodes; ++i)
    n[i] = table[n[i]];
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::renum(EleRenumMap& table)
{
  for(int i = 0; i < NumberOfNodes; ++i)
    n[i] = table[n[i]];
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::markDofs(DofSetArray &dsa) const
{
  dsa.mark(n, NumberOfNodes, DofSet::XYZdisp);
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
int *
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::dofs(DofSetArray &dsa, int *p) const
{
  if(p == 0) p = new int[3*NumberOfNodes];
  for(int i = 0; i < NumberOfNodes; ++i)
    dsa.number(n[i], DofSet::XYZdisp, p+3*i);
  return p;
}

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
int*
SolidElementTemplate<ShapeFunctionTemplate,NumberOfNodes,NumIntgPts>::nodes(int *p) const
{
  if(p == 0) p = new int[NumberOfNodes];
  for(int i = 0; i < NumberOfNodes; ++i)
    p[i] = n[i];
  return p;
}
