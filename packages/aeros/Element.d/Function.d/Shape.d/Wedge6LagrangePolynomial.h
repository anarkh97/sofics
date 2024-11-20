#ifndef _WEDGE6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _WEDGE6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Wedge6LagrangePolynomialShapeFunction : public VectorValuedFunction<3,6,Scalar,0,0,double>
{
  public:
    Wedge6LagrangePolynomialShapeFunction() {}
    Wedge6LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,6,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for six-node wedge element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,6,1> N;

      N[0] = 1/2.*(1-xi-eta)*(1-zeta);
      N[1] = 1/2.*xi*(1-zeta);
      N[2] = 1/2.*eta*(1-zeta);
      N[3] = 1/2.*(1-xi-eta)*(1+zeta);
      N[4] = 1/2.*xi*(1+zeta);
      N[5] = 1/2.*eta*(1+zeta);

      return N;
    }
};

#endif
