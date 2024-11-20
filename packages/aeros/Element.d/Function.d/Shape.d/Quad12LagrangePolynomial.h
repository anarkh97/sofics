#ifndef _QUAD12LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD12LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Quad12LagrangePolynomialShapeFunction : public VectorValuedFunction<2,12,Scalar,0,0,double>
{
  public:
    Quad12LagrangePolynomialShapeFunction() {}
    Quad12LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,12,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for twelve-node quad element: N(ξ,η)
      Eigen::Matrix<Scalar,12,1> N;

      //          η             shape functions     
      //          ↑             ---------------
      //    3---9---8---2       N₀ = 1/32(1-ξ)(1-η)(-10+9(ξξ+ηη))
      //    ¦           ¦       N₁ = 1/32(1+ξ)(1-η)(-10+9(ξξ+ηη))
      //    10          7       N₂ = 1/32(1+ξ)(1+η)(-10+9(ξξ+ηη))
      //    ¦           ¦ → ξ   N₃ = 1/32(1-ξ)(1+η)(-10+9(ξξ+ηη))
      //    11          6       N₄ = 9/32(1-ξξ)(1-η)(1-3ξ)
      //    ¦           ¦       N₅ = 9/32(1-ξξ)(1-η)(1+3ξ)
      //    0---4---5---1       N₆ = 9/32(1-ηη)(1+ξ)(1-3η)
      //                        N₇ = 9/32(1-ηη)(1+ξ)(1+3η)
      //                        N₈ = 9/32(1-ξξ)(1+η)(1+3ξ)
      //                        N₉ = 9/32(1-ξξ)(1+η)(1-3ξ)
      //                        N₁₀= 9/32(1-ηη)(1-ξ)(1+3η)
      //                        N₁₁= 9/32(1-ηη)(1-ξ)(1-3η)

      N[0] = 1/32.*(1-xi)*(1-eta)*(-10+9*(xi*xi+eta*eta));
      N[1] = 1/32.*(1+xi)*(1-eta)*(-10+9*(xi*xi+eta*eta));
      N[2] = 1/32.*(1+xi)*(1+eta)*(-10+9*(xi*xi+eta*eta));
      N[3] = 1/32.*(1-xi)*(1+eta)*(-10+9*(xi*xi+eta*eta));
      N[4] = 9/32.*(1-xi*xi)*(1-eta)*(1-3*xi);
      N[5] = 9/32.*(1-xi*xi)*(1-eta)*(1+3*xi);
      N[6] = 9/32.*(1-eta*eta)*(1+xi)*(1-3*eta);
      N[7] = 9/32.*(1-eta*eta)*(1+xi)*(1+3*eta);
      N[8] = 9/32.*(1-xi*xi)*(1+eta)*(1+3*xi);
      N[9] = 9/32.*(1-xi*xi)*(1+eta)*(1-3*xi);
      N[10] = 9/32.*(1-eta*eta)*(1-xi)*(1+3*eta);
      N[11] = 9/32.*(1-eta*eta)*(1-xi)*(1-3*eta);

      return N;
    }
};

#endif
