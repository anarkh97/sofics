#ifndef _TET4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TET4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

template<typename Scalar>
class Tet4LagrangePolynomialShapeFunction : public VectorValuedFunction<3,4,Scalar,0,0,double>
{
  public:
    Tet4LagrangePolynomialShapeFunction() {}
    Tet4LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for four-node tet element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,4,1> N;

      N[0] = 1-xi-eta-zeta;
      N[1] = xi;
      N[2] = eta;
      N[3] = zeta;

      return N;
    }
};

namespace Simo {

template<>
Eigen::Matrix<double,4,3>
Jacobian<double,Tet4LagrangePolynomialShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double);

}

#endif
