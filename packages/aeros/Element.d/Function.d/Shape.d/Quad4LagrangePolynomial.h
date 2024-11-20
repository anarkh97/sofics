#ifndef _QUAD4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Quad4LagrangePolynomialShapeFunction : public VectorValuedFunction<2,4,Scalar,0,0,double>
{
  public:
    Quad4LagrangePolynomialShapeFunction() {}
    Quad4LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for four-node quad element: N(ξ,η)
      Eigen::Matrix<Scalar,4,1> N;

      //        η           shape functions
      //        ↑           ---------------
      //    3-------2       N₀ = 1/4(1-ξ)(1-η)
      //    ¦       ¦       N₁ = 1/4(1+ξ)(1-η)
      //    ¦       ¦ → ξ   N₂ = 1/4(1+ξ)(1+η)
      //    ¦       ¦       N₃ = 1/4(1-ξ)(1+η)
      //    0-------1

      N[0] = 1/4.*(1-xi)*(1-eta);
      N[1] = 1/4.*(1+xi)*(1-eta);
      N[2] = 1/4.*(1+xi)*(1+eta);
      N[3] = 1/4.*(1-xi)*(1+eta);

      return N;
    }
};

#endif
