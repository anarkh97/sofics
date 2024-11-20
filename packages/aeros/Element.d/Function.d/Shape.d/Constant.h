#ifndef _CONSTANTSHAPEFUNCTION_H_
#define _CONSTANTSHAPEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

template<typename Scalar>
class ConstantShapeFunction : public VectorValuedFunction<3,1,Scalar,0,0,double>
{
  public:
    ConstantShapeFunction() {}
    ConstantShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&) {}

    Eigen::Matrix<Scalar,1,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      return Eigen::Matrix<Scalar,1,1>::Ones();
    }
};

namespace Simo {

template<>
Eigen::Matrix<double,1,3>
Jacobian<double,ConstantShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double);

} // namespace Simo

#endif
