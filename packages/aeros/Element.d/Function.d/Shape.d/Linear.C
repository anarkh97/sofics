#ifdef USE_EIGEN3
#include <Element.d/Function.d/Shape.d/Linear.h>

namespace Simo {

template<>
Eigen::Matrix<double,4,3>
Jacobian<double,LinearShapeFunction>
::operator() (const Eigen::Matrix<double,3,1>&, double)
{
  Eigen::Matrix<double,4,3> dShape;

  dShape <<  0.0,  0.0,  0.0,
             1.0,  0.0,  0.0,
             0.0,  1.0,  0.0,
             0.0,  0.0,  1.0;

  return dShape;
}

}

#endif

