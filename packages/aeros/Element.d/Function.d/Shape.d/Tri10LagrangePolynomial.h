#ifndef _TRI10LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TRI10LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Tri10LagrangePolynomialShapeFunction : public VectorValuedFunction<2,10,Scalar,0,0,double>
{
  public:
    Tri10LagrangePolynomialShapeFunction() {}
    Tri10LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)

    Eigen::Matrix<Scalar,10,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η]
      const Scalar &xi = q[0], &eta = q[1];

      // outputs:
      // shape functions for ten-node tri element: N(ξ,η)
      Eigen::Matrix<Scalar,10,1> N;

      //    η               shape functions
      //    ↑               ---------------
      //    1               N₀ = 9/2ξ(ξ(ξ-1)+2/9)
      //    ¦\              N₁ = 9/2η(η(η-1)+2/9)
      //    ¦ \             N₂ = 9/2ζ(ζ(ζ-1)+2/9), ζ = 1-ξ-η
      //    5  4            N₃ = 9/2ηξ(3ξ-1)
      //    ¦   \           N₄ = 9/2ηξ(3η-1)
      //    ¦    \          N₅ = 9/2ζη(3η-1)
      //    6     3         N₆ = 9/2ζη(3ζ-1)
      //    ¦      \        N₇ = 9/2ζξ(3ζ-1)
      //    ¦       \       N₈ = 9/2ζξ(3ξ-1)
      //    2--7--8--0 → ξ  N₉ = 27ξηζ
     
      Scalar zeta = 1-xi-eta;
      N[0] = 9/2.*xi*(xi*(xi-1)+2/9.);
      N[1] = 9/2.*eta*(eta*(eta-1)+2/9.);
      N[2] = 9/2.*zeta*(zeta*(zeta-1)+2/9.);
      N[3] = 9/2.*eta*xi*(3*xi-1);
      N[4] = 9/2.*eta*xi*(3*eta-1);
      N[5] = 9/2.*zeta*eta*(3*eta-1);
      N[6] = 9/2.*zeta*eta*(3*zeta-1);
      N[7] = 9/2.*zeta*xi*(3*zeta-1);
      N[8] = 9/2.*zeta*xi*(3*xi-1);
      N[9] = 27*xi*eta*zeta;

      return N;
    }
};

#endif
