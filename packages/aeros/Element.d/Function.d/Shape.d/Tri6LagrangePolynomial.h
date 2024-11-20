#ifndef _TRI6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TRI6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Tri6LagrangePolynomialShapeFunction : public VectorValuedFunction<2,6,Scalar,0,0,double>
{
  public:
    Tri6LagrangePolynomialShapeFunction() {}
    Tri6LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,6,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for six-node tri element: N(ξ,η)
      Eigen::Matrix<Scalar,6,1> N;

      //    η              shape functions
      //    ↑              ---------------
      //    1              N₀ = ξ(2ξ-1)
      //    ¦\             N₁ = η(2η-1)
      //    ¦ \            N₂ = ζ(2ζ-1), ζ = 1-ξ-η
      //    4  3           N₃ = 4ξη
      //    ¦   \          N₄ = 4ηζ
      //    ¦    \         N₅ = 4ζξ
      //    2--5--0 → ξ

      Scalar zeta = 1-xi-eta;
      N[0] = xi*(2*xi-1);
      N[1] = eta*(2*eta-1);
      N[2] = zeta*(2*zeta-1);
      N[3] = 4*xi*eta;
      N[4] = 4*eta*zeta;
      N[5] = 4*zeta*xi;

      return N;
    }
};

#endif
