#ifndef _HEX8LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _HEX8LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

template<typename Scalar>
class Hex8LagrangePolynomialShapeFunction : public VectorValuedFunction<3,8,Scalar,0,0,double>
{
  public:
    Hex8LagrangePolynomialShapeFunction() {}
    Hex8LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,8,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for eight-node hex element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,8,1> N;

      N[0] = 1/8.*(1-xi)*(1-eta)*(1-zeta);
      N[1] = 1/8.*(1+xi)*(1-eta)*(1-zeta);
      N[2] = 1/8.*(1+xi)*(1+eta)*(1-zeta);
      N[3] = 1/8.*(1-xi)*(1+eta)*(1-zeta);
      N[4] = 1/8.*(1-xi)*(1-eta)*(1+zeta);
      N[5] = 1/8.*(1+xi)*(1-eta)*(1+zeta);
      N[6] = 1/8.*(1+xi)*(1+eta)*(1+zeta);
      N[7] = 1/8.*(1-xi)*(1+eta)*(1+zeta);

      return N;
    }
};

namespace Simo {

template<>
Eigen::Matrix<double,8,3>
Jacobian<double,Hex8LagrangePolynomialShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double);

}

#endif
