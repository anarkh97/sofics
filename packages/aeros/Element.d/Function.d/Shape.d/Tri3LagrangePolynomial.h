#ifndef _TRI3LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TRI3LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Tri3LagrangePolynomialShapeFunction : public VectorValuedFunction<2,3,Scalar,0,0,double>
{
  public:
    Tri3LagrangePolynomialShapeFunction() {}
    Tri3LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for three-node tri element: N(ξ,η)
      Eigen::Matrix<Scalar,3,1> N;

      //    η              shape functions
      //    ↑              ---------------
      //    1              N₀ = ξ
      //    ¦\             N₁ = η
      //    ¦ \            N₂ = 1-ξ-η
      //    ¦  \
      //    2---0 → ξ

      N[0] = xi;
      N[1] = eta;
      N[2] = 1-xi-eta;

      return N;
    }
};

#endif
