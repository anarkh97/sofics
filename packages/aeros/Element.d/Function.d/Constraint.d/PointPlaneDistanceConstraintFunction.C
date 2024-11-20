#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/PointPlaneDistanceConstraintFunction.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for point-
// plane distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,3>
Jacobian<double,PointPlaneDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix<double,3,1> J;

  Eigen::Vector3d x0 = sconst.segment<3>(0);                                       
  Eigen::Vector3d x1 = sconst.segment<3>(3);
  Eigen::Vector3d x2 = sconst.segment<3>(6);
  Eigen::Vector3d x3 = sconst.segment<3>(9);

  // unit normal to plane
  Eigen::Vector3d nhat = (x2-x1).cross(x3-x1).normalized();

  if(iconst[0]) 
    J = -nhat;
  else
    J = nhat;

  return J.transpose();
}

// specializing the member function template of constraint hessian operator for point-
// plane distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,3,3>
Hessian<double,PointPlaneDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix<double,3,3> H;
  H.setZero();

  return H;
}

} // namespace Simo

#endif
