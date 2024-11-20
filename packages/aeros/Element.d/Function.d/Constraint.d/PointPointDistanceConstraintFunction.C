#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/PointPointDistanceConstraintFunction.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for point-
// point distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,3>
Jacobian<double,PointPointDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix<double,3,1> J;

  Eigen::Vector3d x0 = sconst.segment<3>(0);                                       
  Eigen::Vector3d x1 = sconst.segment<3>(3);

  // node's current coordinates
  x0 += q;

  Eigen::Vector3d dhat = (x0-x1).normalized();

  if(iconst[0])
    J = -dhat;
  else 
    J = dhat;

  return J.transpose();
}

// specializing the member function template of constraint hessian operator for point-
// point distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,3,3>
Hessian<double,PointPointDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix<double,3,3> H;

  Eigen::Vector3d x0 = sconst.segment<3>(0); 
  Eigen::Vector3d x1 = sconst.segment<3>(3);

  // node's current coordinates
  x0 += q;

  Eigen::Vector3d d = x0-x1;
  double l2 = d.squaredNorm();
  using std::sqrt;
  double l3 = sqrt(l2)*l2;

  if(iconst[0]) {
    H(0,0) = -(l2-d[0]*d[0])/l3;
    H(1,1) = -(l2-d[1]*d[1])/l3;
    H(2,2) = -(l2-d[2]*d[2])/l3;

    H(0,1) = H(1,0) = d[0]*d[1]/l3;
    H(0,2) = H(2,0) = d[0]*d[2]/l3;
    H(1,2) = H(2,1) = d[1]*d[2]/l3;
  }
  else {
    H(0,0) = (l2-d[0]*d[0])/l3;
    H(1,1) = (l2-d[1]*d[1])/l3;
    H(2,2) = (l2-d[2]*d[2])/l3;

    H(0,1) = H(1,0) = -d[0]*d[1]/l3;
    H(0,2) = H(2,0) = -d[0]*d[2]/l3;
    H(1,2) = H(2,1) = -d[1]*d[2]/l3;
  }

  return H;
}

} // namespace Simo

#endif
