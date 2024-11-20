#ifndef _WEDGE26LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _WEDGE26LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Wedge26LagrangePolynomialShapeFunction : public VectorValuedFunction<3,26,Scalar,0,0,double>
{
  public:
    Wedge26LagrangePolynomialShapeFunction() {}
    Wedge26LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,26,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for 26-node wedge element: N(ξ,η,ζ)
      // see https://frg.bitbucket.org/aero-s/index.html#TPF17 for node numbering convention
      Eigen::Matrix<Scalar,26,1> N;

      N[0] = (2-3*xi-3*eta)*(1-3*xi-3*eta)*(1-xi-eta)*(1-zeta)/4.-9/8.*(1+zeta)*(zeta-1/3.)*(zeta-1)*(1-xi-eta)+9/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*(1-xi-eta);
      N[1] = (3*xi-1)*(3*xi-2)*xi*(1-zeta)/4.-9/8.*(1+zeta)*(zeta-1/3.)*(zeta-1)*xi+9/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*xi;
      N[2] = (3*eta-1)*(3*eta-2)*eta*(1-zeta)/4.-9/8.*(1+zeta)*(zeta-1/3.)*(zeta-1)*eta+9/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*eta;
      N[3] = (2-3*xi-3*eta)*(1-3*xi-3*eta)*(1-xi-eta)*(1+zeta)/4.-9/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*(1-xi-eta)+9/8.*(1+zeta)*(zeta+1/3.)*(zeta-1)*(1-xi-eta);
      N[4] = (3*xi-1)*(3*xi-2)*xi*(1+zeta)/4.-9/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*xi+9/8.*(1+zeta)*(zeta+1/3.)*(zeta-1)*xi;
      N[5] = (3*eta-1)*(3*eta-2)*eta*(1+zeta)/4.-9/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*eta+9/8.*(1+zeta)*(zeta+1/3.)*(zeta-1)*eta;
      N[6] = 9/4.*xi*(1-xi-eta)*(2-3*xi-3*eta)*(1-zeta);
      N[7] = 9/4.*xi*(1-xi-eta)*(3*xi-1)*(1-zeta);
      N[8] = 9/4.*xi*eta*(3*xi-1)*(1-zeta);
      N[9] = 9/4.*xi*eta*(3*eta-1)*(1-zeta);
      N[10] = 9/4.*(1-xi-eta)*eta*(3*eta-1)*(1-zeta);
      N[11] = 9/4.*(1-xi-eta)*eta*(2-3*xi-3*eta)*(1-zeta);
      N[12] = 9/4.*xi*(1-xi-eta)*(2-3*xi-3*eta)*(1+zeta);
      N[13] = 9/4.*xi*(1-xi-eta)*(3*xi-1)*(1+zeta);
      N[14] = 9/4.*xi*eta*(3*xi-1)*(1+zeta);
      N[15] = 9/4.*xi*eta*(3*eta-1)*(1+zeta);
      N[16] = 9/4.*(1-xi-eta)*eta*(3*eta-1)*(1+zeta);
      N[17] = 9/4.*(1-xi-eta)*eta*(2-3*xi-3*eta)*(1+zeta);
      N[18] = 27/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*(1-xi-eta);
      N[19] = -27/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*(1-xi-eta);
      N[20] = 27/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*xi;
      N[21] = -27/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*xi;
      N[22] = 27/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*eta;
      N[23] = -27/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*eta;
      N[24] = 27/2.*xi*eta*(1-xi-eta)*(1-zeta);
      N[25] = 27/2.*xi*eta*(1-xi-eta)*(1+zeta);

      return N;
    }
};

#endif
