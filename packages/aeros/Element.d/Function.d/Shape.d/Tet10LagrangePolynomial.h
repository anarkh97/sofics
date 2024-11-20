#ifndef _TET10LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TET10LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Tet10LagrangePolynomialShapeFunction : public VectorValuedFunction<3,10,Scalar,0,0,double>
{
  public:
    Tet10LagrangePolynomialShapeFunction() {}
    Tet10LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,10,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for ten-node tet element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,10,1> N;

      N[0] = (1-xi-eta-zeta)*(2*(1-xi-eta-zeta)-1);
      N[1] = xi*(2*xi-1);
      N[2] = eta*(2*eta-1);
      N[3] = zeta*(2*zeta-1);
      N[4] = 4*xi*(1-xi-eta-zeta);
      N[5] = 4*xi*eta;
      N[6] = 4*eta*(1-xi-eta-zeta);
      N[7] = 4*zeta*(1-xi-eta-zeta);
      N[8] = 4*zeta*xi;
      N[9] = 4*zeta*eta;

      return N;
    }
};

#endif
