#ifndef _QUAD8LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD8LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Quad8LagrangePolynomialShapeFunction : public VectorValuedFunction<2,8,Scalar,0,0,double>
{
  public:
    Quad8LagrangePolynomialShapeFunction() {}
    Quad8LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,8,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for eight-node quad element: N(ξ,η)
      Eigen::Matrix<Scalar,8,1> N;

      //          η             shape functions
      //          ↑             ---------------
      //    3-----6-----2       N₀ = 1/4(1-ξ)(1-η)(-1-ξ-η)
      //    ¦           ¦       N₁ = 1/4(1+ξ)(1-η)(-1+ξ-η) 
      //    ¦           ¦       N₂ = 1/4(1+ξ)(1+η)(-1+ξ+η) 
      //    7           5 → ξ   N₃ = 1/4(1-ξ)(1+η)(-1-ξ+η) 
      //    ¦           ¦       N₄ = 1/2(1-ξξ)(1-η) 
      //    ¦           ¦       N₅ = 1/2(1+ξ)(1-ηη) 
      //    0-----4-----1       N₆ = 1/2(1-ξξ)(1+η) 
      //                        N₇ = 1/2(1-ξ)(1-ηη)

      N[0] = 1/4.*(1-xi)*(1-eta)*(-1-xi-eta);
      N[1] = 1/4.*(1+xi)*(1-eta)*(-1+xi-eta);
      N[2] = 1/4.*(1+xi)*(1+eta)*(-1+xi+eta);
      N[3] = 1/4.*(1-xi)*(1+eta)*(-1-xi+eta);
      N[4] = 1/2.*(1-xi*xi)*(1-eta);
      N[5] = 1/2.*(1+xi)*(1-eta*eta);
      N[6] = 1/2.*(1-xi*xi)*(1+eta);
      N[7] = 1/2.*(1-xi)*(1-eta*eta);

      return N;
    }
};

#endif
