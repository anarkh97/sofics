#ifndef _LINEARSHAPEFUNCTION_H_
#define _LINEARSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

template<typename Scalar>
class LinearShapeFunction : public VectorValuedFunction<3,4,Scalar,0,0,double>
{
  public:
    LinearShapeFunction() {}
    LinearShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // local coordinates of point at which function is to be evaluated: q = [ξ,η,ζ]
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];

      // outputs:
      // shape functions for interpolating pressure field in a mixed Q2-P1 eight-node hex element: N(ξ,η,ζ)
      Eigen::Matrix<Scalar,4,1> N;

      N[0] = 1;
      N[1] = xi;
      N[2] = eta;
      N[3] = zeta;

      return N;
    }
};

namespace Simo {

template<>
Eigen::Matrix<double,4,3>
Jacobian<double,LinearShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double);

} // namespace Simo

#endif
