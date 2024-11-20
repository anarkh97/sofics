#ifndef _HEX32LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _HEX32LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>

template<typename Scalar>
class Hex32LagrangePolynomialShapeFunction : public VectorValuedFunction<3,32,Scalar,0,0,double>
{
  public:
    Hex32LagrangePolynomialShapeFunction() {}
    Hex32LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,32,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for 32-node hex element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,32,1> N;

      N[0 ] = 1/64.*(1-xi)*(1-eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[1 ] = 1/64.*(1+xi)*(1-eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[2 ] = 1/64.*(1+xi)*(1+eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[3 ] = 1/64.*(1-xi)*(1+eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[4 ] = 1/64.*(1-xi)*(1-eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[5 ] = 1/64.*(1+xi)*(1-eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[6 ] = 1/64.*(1+xi)*(1+eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[7 ] = 1/64.*(1-xi)*(1+eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      N[8 ] = 9/64.*(1-xi*xi)*(1-3*xi)*(1-eta)*(1-zeta);
      N[9 ] = 9/64.*(1-xi*xi)*(1+3*xi)*(1-eta)*(1-zeta);
      N[10] = 9/64.*(1-eta*eta)*(1-3*eta)*(1+xi)*(1-zeta);
      N[11] = 9/64.*(1-eta*eta)*(1+3*eta)*(1+xi)*(1-zeta);
      N[12] = 9/64.*(1-xi*xi)*(1+3*xi)*(1+eta)*(1-zeta);
      N[13] = 9/64.*(1-xi*xi)*(1-3*xi)*(1+eta)*(1-zeta);
      N[14] = 9/64.*(1-eta*eta)*(1+3*eta)*(1-xi)*(1-zeta);
      N[15] = 9/64.*(1-eta*eta)*(1-3*eta)*(1-xi)*(1-zeta);
      N[16] = 9/64.*(1-xi*xi)*(1-3*xi)*(1-eta)*(1+zeta);
      N[17] = 9/64.*(1-xi*xi)*(1+3*xi)*(1-eta)*(1+zeta);
      N[18] = 9/64.*(1-eta*eta)*(1-3*eta)*(1+xi)*(1+zeta);
      N[19] = 9/64.*(1-eta*eta)*(1+3*eta)*(1+xi)*(1+zeta);
      N[20] = 9/64.*(1-xi*xi)*(1+3*xi)*(1+eta)*(1+zeta);
      N[21] = 9/64.*(1-xi*xi)*(1-3*xi)*(1+eta)*(1+zeta);
      N[22] = 9/64.*(1-eta*eta)*(1+3*eta)*(1-xi)*(1+zeta);
      N[23] = 9/64.*(1-eta*eta)*(1-3*eta)*(1-xi)*(1+zeta);
      N[24] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1-xi)*(1-eta);
      N[25] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1-xi)*(1-eta);
      N[26] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1+xi)*(1-eta);
      N[27] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1+xi)*(1-eta);
      N[28] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1+xi)*(1+eta);
      N[29] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1+xi)*(1+eta);
      N[30] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1-xi)*(1+eta);
      N[31] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1-xi)*(1+eta);
    
      return N;
    }
};

#endif
