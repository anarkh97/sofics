#ifdef USE_EIGEN3
#include <Element.d/Function.d/Shape.d/Hex8LagrangePolynomial.h>

namespace Simo {

template<>
Eigen::Matrix<double,8,3>
Jacobian<double,Hex8LagrangePolynomialShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double)
{
  const double &xi = q[0], &eta = q[1], &zeta = q[2];

  const double d1 = 0.5*(1.0+xi);
  const double d2 = 0.5*(1.0+eta);
  const double d3 = 0.5*(1.0+zeta);
  const double d4 = 1.0-d1;
  const double d5 = 1.0-d2;
  const double d6 = 1.0-d3;

  const double d1h = 0.5*d1;
  const double d2h = 0.5*d2;
  const double d3h = 0.5*d3;
  const double d4h = 0.5*d4;
  const double d5h = 0.5*d5;
  const double d6h = 0.5*d6;

  Eigen::Matrix<double,8,3> dShape;

  dShape(1,0) =  d5h*d6;
  dShape(0,0) = -dShape(1,0);
  dShape(2,0) =  d2h*d6;
  dShape(3,0) = -dShape(2,0);
  dShape(5,0) =  d5h*d3;
  dShape(4,0) = -dShape(5,0);
  dShape(6,0) =  d2h*d3;
  dShape(7,0) = -dShape(6,0);

  dShape(2,1) =  d1h*d6;
  dShape(3,1) =  d4h*d6;
  dShape(0,1) = -dShape(3,1);
  dShape(1,1) = -dShape(2,1);
  dShape(6,1) =  d1h*d3;
  dShape(7,1) =  d4h*d3;
  dShape(4,1) = -dShape(7,1);
  dShape(5,1) = -dShape(6,1);

  dShape(4,2) =  d4h*d5;
  dShape(5,2) =  d1h*d5;
  dShape(6,2) =  d1h*d2;
  dShape(7,2) =  d4h*d2;
  dShape(0,2) = -dShape(4,2);
  dShape(1,2) = -dShape(5,2);
  dShape(2,2) = -dShape(6,2);
  dShape(3,2) = -dShape(7,2);

  return dShape;
}

}

#endif
