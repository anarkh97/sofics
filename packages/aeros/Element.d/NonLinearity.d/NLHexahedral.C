#ifdef USE_EIGEN3
#include <Element.d/NonLinearity.d/NLHexahedral.h>

template<>
const double SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>
::nodeRefCoords[8][3] = 
  {{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0},{ 1.0,1.0,-1.0},{-1.0,1.0,-1.0},
   {-1.0,-1.0, 1.0},{1.0,-1.0, 1.0},{ 1.0,1.0, 1.0},{-1.0,1.0, 1.0}};

template<>
const double SolidElementTemplate<Hex20LagrangePolynomialShapeFunction,20,27>
::nodeRefCoords[20][3] = 
  {{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0},{ 1.0,1.0,-1.0},{-1.0,1.0,-1.0},
   {-1.0,-1.0, 1.0},{1.0,-1.0, 1.0},{ 1.0,1.0, 1.0},{-1.0,1.0, 1.0},
   { 0.0,-1.0,-1.0},{1.0, 0.0,-1.0},{ 0.0,1.0,-1.0},{-1.0,0.0,-1.0},
   { 0.0,-1.0, 1.0},{1.0, 0.0, 1.0},{ 0.0,1.0, 1.0},{-1.0,0.0, 1.0},
   { 0.0,-1.0, 0.0},{1.0, 0.0, 0.0},{ 0.0,1.0, 0.0},{-1.0,0.0, 0.0}};

template<>
const double SolidElementTemplate<Hex32LagrangePolynomialShapeFunction,32,64>
::nodeRefCoords[32][3] = 
  {{-1.   ,-1.,-1.},{ 1.   ,-1.,-1.},{ 1., 1.   ,-1.},{-1., 1.   ,-1.},
   {-1.   ,-1., 1.},{ 1.   ,-1., 1.},{ 1., 1.   , 1.},{-1., 1.   , 1.},
   {-1./3.,-1.,-1.},{ 1./3.,-1.,-1.},{ 1.,-1./3.,-1.},{ 1., 1./3.,-1.},
   { 1./3., 1.,-1.},{-1./3., 1.,-1.},{-1., 1./3.,-1.},{-1.,-1./3.,-1.},
   {-1./3.,-1., 1.},{ 1./3.,-1., 1.},{ 1.,-1./3., 1.},{ 1., 1./3., 1.},
   { 1./3., 1., 1.},{-1./3., 1., 1.},{-1., 1./3., 1.},{-1.,-1./3., 1.},
   {-1.,-1.,-1./3.},{-1.,-1., 1./3.},{ 1.,-1.,-1./3.},{ 1.,-1., 1./3.},
   { 1., 1.,-1./3.},{ 1., 1., 1./3.},{-1., 1.,-1./3.},{-1., 1., 1./3.}};

extern LinearStrain linearStrain;
extern GreenLagrangeStrain greenLagrangeStrain;

StrainEvaluator *
NLHexahedral::getStrainEvaluator() const
{
  // note: for this element the default strain evaluator (returned by NLMaterial::getStrainEvaluator) will not be used
  if(linearKinematics)
    return &linearStrain;
  else
    return &greenLagrangeStrain;
}

template<>
void
SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>
::getGaussPointAndWeight(int n, double *point, double &weight) const
{
  int i, j, k;
  i = n%2;
  n /= 2;
  j = n%2;
  n /= 2;
  k = n;

  static double xi[2] = {
    -1/std::sqrt(3),
     1/std::sqrt(3)
  };

  static double wi[2] = {
    1.,
    1.
  };
  weight = wi[i]*wi[j]*wi[k];
  point[0] = xi[i];
  point[1] = xi[j];
  point[2] = xi[k];
}

template<>
void
SolidElementTemplate<Hex20LagrangePolynomialShapeFunction,20,27>
::getGaussPointAndWeight(int n, double *point, double &weight) const
{
  int i, j, k;
  i = n%3;
  n /= 3;
  j = n%3;
  n /= 3;
  k = n;

  static double xi[3] = {
    -std::sqrt(15)/5,
     0.,
     std::sqrt(15)/5
  };

  static double wi[3] = {
    5.0/9.0,
    8.0/9.0,
    5.0/9.0
  };

  weight = wi[i]*wi[j]*wi[k];
  point[0] = xi[i];
  point[1] = xi[j];
  point[2] = xi[k];
}

template<>
void
SolidElementTemplate<Hex32LagrangePolynomialShapeFunction,32,64>
::getGaussPointAndWeight(int n, double *point, double &weight) const
{
  int i, j, k;
  i = n%4;
  n /= 4;
  j = n%4;
  n /= 4;
  k = n;

  static double xi[4] = {
    -0.861136311594053,
    -0.339981043584856,
     0.339981043584856,
     0.861136311594053
  };

  static double wi[4] = {
    0.347854845137454,
    0.652145154862546,
    0.652145154862546,
    0.347854845137454
  };

  weight = wi[i]*wi[j]*wi[k];
  point[0] = xi[i];
  point[1] = xi[j];
  point[2] = xi[k];
}
#endif
