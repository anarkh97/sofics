#ifndef _QUAD9LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD9LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Quad9LagrangePolynomialShapeFunction : public VectorValuedFunction<2,9,Scalar,0,0,double>
{
  public:
    Quad9LagrangePolynomialShapeFunction() {}
    Quad9LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,9,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for nine-node quad element: N(ξ,η)
      Eigen::Matrix<Scalar,9,1> N;

      //          η             shape functions
      //          ↑             ---------------
      //    3-----6-----2       N₀ = 1/4(ξ-1)(η-1)ξη
      //    ¦           ¦       N₁ = 1/4(ξ+1)(η-1)ξη
      //    ¦           ¦       N₂ = 1/4(ξ+1)(η+1)ξη
      //    7     8     5 → ξ   N₃ = 1/4(ξ-1)(η+1)ξη
      //    ¦           ¦       N₄ = 1/2(1-ξξ)η(η-1)
      //    ¦           ¦       N₅ = 1/2(1-ηη)ξ(ξ+1)
      //    0-----4-----1       N₆ = 1/2(1-ξξ)η(η+1)
      //                        N₇ = 1/2(1-ηη)ξ(ξ-1)
      //                        N₈ = (1-ξξ)(1-ηη)

      N[0] = 1/4.*(xi-1.0)*(eta-1.0)*xi*eta;
      N[1] = 1/4.*(xi+1.0)*(eta-1.0)*xi*eta;
      N[2] = 1/4.*(xi+1.0)*(eta+1.0)*xi*eta;
      N[3] = 1/4.*(xi-1.0)*(eta+1.0)*xi*eta;
      N[4] = 1/2.*(1.0-xi*xi)*eta*(eta-1.0);
      N[5] = 1/2.*(1.0-eta*eta)*xi*(xi+1.0);
      N[6] = 1/2.*(1.0-xi*xi)*eta*(eta+1.0);
      N[7] = 1/2.*(1.0-eta*eta)*xi*(xi-1.0);
      N[8] = (1.0-xi*xi)*(1.0-eta*eta);

      return N;
    }
};

#endif
