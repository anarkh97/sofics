#ifndef _HEX20LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _HEX20LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Hex20LagrangePolynomialShapeFunction : public VectorValuedFunction<3,20,Scalar,0,0,double>
{
  public:
    Hex20LagrangePolynomialShapeFunction() {}
    Hex20LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,20,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for twenty-node hex element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,20,1> N;

      N[0]  = 1/8.*(1-xi)*(1-eta)*(1-zeta)*(-xi-eta-zeta-2);
      N[1]  = 1/8.*(1+xi)*(1-eta)*(1-zeta)*( xi-eta-zeta-2);
      N[2]  = 1/8.*(1+xi)*(1+eta)*(1-zeta)*( xi+eta-zeta-2);
      N[3]  = 1/8.*(1-xi)*(1+eta)*(1-zeta)*(-xi+eta-zeta-2);
      N[4]  = 1/8.*(1-xi)*(1-eta)*(1+zeta)*(-xi-eta+zeta-2);
      N[5]  = 1/8.*(1+xi)*(1-eta)*(1+zeta)*( xi-eta+zeta-2);
      N[6]  = 1/8.*(1+xi)*(1+eta)*(1+zeta)*( xi+eta+zeta-2);
      N[7]  = 1/8.*(1-xi)*(1+eta)*(1+zeta)*(-xi+eta+zeta-2);
      N[8]  = 1/4.*(1-xi*xi)*(1-eta)*(1-zeta);
      N[9]  = 1/4.*(1+xi)*(1-eta*eta)*(1-zeta);
      N[10] = 1/4.*(1-xi*xi)*(1+eta)*(1-zeta);
      N[11] = 1/4.*(1-xi)*(1-eta*eta)*(1-zeta);
      N[12] = 1/4.*(1-xi*xi)*(1-eta)*(1+zeta);
      N[13] = 1/4.*(1+xi)*(1-eta*eta)*(1+zeta);
      N[14] = 1/4.*(1-xi*xi)*(1+eta)*(1+zeta);
      N[15] = 1/4.*(1-xi)*(1-eta*eta)*(1+zeta);
      N[16] = 1/4.*(1-xi)*(1-eta)*(1-zeta*zeta);
      N[17] = 1/4.*(1+xi)*(1-eta)*(1-zeta*zeta);
      N[18] = 1/4.*(1+xi)*(1+eta)*(1-zeta*zeta);
      N[19] = 1/4.*(1-xi)*(1+eta)*(1-zeta*zeta);

      return N;
    }
};

#endif
