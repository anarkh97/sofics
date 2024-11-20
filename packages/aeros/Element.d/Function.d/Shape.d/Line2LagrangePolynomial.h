#ifndef _LINE2LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _LINE2LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Line2LagrangePolynomialShapeFunction : public VectorValuedFunction<1,2,Scalar,0,0,double>
{
  public:
    Line2LagrangePolynomialShapeFunction( {}
    Line2LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,2,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // local coordinate of point at which function is to be evaluated: q = [ξ]
      const Scalar &xi = q[0];

      // outputs:
      // shape functions for two-node line element: N(ξ)
      Eigen::Matrix<Scalar,2,1> N;

      //                   shape functions
      //       η           ---------------
      //       ↑           N₀ = 1/2(1-ξ)
      //    0-----1 → ξ    N₁ = 1/2(1+ξ)

      N[0] = 1/2.*(1-xi);
      N[1] = 1/2.*(1+xi);

      return N;
    }
};

#endif
