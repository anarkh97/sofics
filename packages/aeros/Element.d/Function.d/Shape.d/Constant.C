#ifdef USE_EIGEN3
#include <Element.d/Function.d/Shape.d/Constant.h>

namespace Simo {

template<>
Eigen::Matrix<double,1,3>
Jacobian<double,ConstantShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>&, double)
{
  return Eigen::Matrix<double,1,3>::Zero();
}

}

#endif

