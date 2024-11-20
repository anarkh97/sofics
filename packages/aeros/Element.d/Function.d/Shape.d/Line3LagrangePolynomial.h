#ifndef _LINE3LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _LINE3LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Line3LagrangePolynomialShapeFunction : public VectorValuedFunction<1,3,Scalar,0,0,double>
{
  public:
    Line3LagrangePolynomialShapeFunction() {}
    Line3LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // local coordinate of point at which function is to be evaluated: q = [ξ]
      const Scalar &xi = q[0];

      // outputs:
      // shape functions for three-node line element: N(ξ)
      Eigen::Matrix<Scalar,3,1> N;

      //                   shape functions
      //       η           ---------------
      //       ↑           N₀ = 1/2ξ(ξ-1)
      //   0---1---2 → ξ   N₁ = 1-ξξ
      //                   N₂ = 1/2ξ(ξ+1)           

      N[0] = 1/2.*xi*(xi-1);
      N[1] = 1-xi*xi;
      N[2] = 1/2.*xi*(xi+1);

      return N;
    }
};

#endif
