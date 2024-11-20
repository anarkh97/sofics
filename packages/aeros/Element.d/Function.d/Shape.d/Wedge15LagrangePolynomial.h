#ifndef _WEDGE15LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _WEDGE15LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Wedge15LagrangePolynomialShapeFunction : public VectorValuedFunction<3,15,Scalar,0,0,double>
{
  public:
    Wedge15LagrangePolynomialShapeFunction() {}
    Wedge15LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,15,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for fifteen-node wedge element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,15,1> N;

      N[0] = -0.5*(1-xi-eta)*(1.0-zeta)*(2.0*xi+2.0*eta+zeta);
      N[1] = 0.5*xi*(1.0-zeta)*(2.0*xi-2.0-zeta);
      N[2] = 0.5*eta*(1.0-zeta)*(2.0*eta-2.0-zeta);
      N[3] = -0.5*(1-xi-eta)*(1.0+zeta)*(2.0*xi+2.0*eta-zeta);
      N[4] = 0.5*xi*(1.0+zeta)*(2.0*xi-2.0+zeta);
      N[5] = 0.5*eta*(1.0+zeta)*(2.0*eta-2.0+zeta);
      N[6] = 2.0*xi*(1-xi-eta)*(1.0-zeta);
      N[7] = 2.0*xi*eta*(1.0-zeta);
      N[8] = 2.0*eta*(1-xi-eta)*(1.0-zeta);
      N[9] = 2.0*xi*(1-xi-eta)*(1.0+zeta);
      N[10] = 2.0*xi*eta*(1.0+zeta);
      N[11] = 2.0*eta*(1-xi-eta)*(1.0+zeta);
      N[12] = (1-xi-eta)*(1.0-zeta*zeta);
      N[13] = xi*(1.0-zeta*zeta);
      N[14] = eta*(1.0-zeta*zeta);

      return N;
    }
};

#endif
