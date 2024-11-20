#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/PointVariPointDistanceConstraintFunction.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for point to
// moving point distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,6>
Jacobian<double,PointVariPointDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t)
{
  Eigen::Matrix<double,6,1> J;

  Eigen::Vector3d x0 = sconst.segment<3>(0).matrix() + q.segment<3>(0);
  Eigen::Vector3d x1 = sconst.segment<3>(3).matrix() + q.segment<3>(3);
  Eigen::Vector3d a = x0-x1;

  double l = a.norm();
  
  if(iconst[0]) {
    for(int i = 0; i < 3; ++i) {
      J[0+i] = -a[i]/l;
      J[3+i] = a[i]/l;
    }
  }
  else {
    for(int i = 0; i < 3; ++i) {
      J[0+i] = a[i]/l;
      J[3+i] = -a[i]/l;
    }
  }

  return J.transpose();
}

// specializing the member function template of constraint hessian operator for point to
// moving point distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,6,6>
Hessian<double,PointVariPointDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t)
{
  Eigen::Matrix<double,6,6> H;

  using std::sqrt;

  Eigen::Vector3d x0 = sconst.segment<3>(0).matrix() + q.segment<3>(0);
  Eigen::Vector3d x1 = sconst.segment<3>(3).matrix() + q.segment<3>(3);
  Eigen::Vector3d a = x0-x1;

  double l2 = a.squaredNorm();
  double l3 = sqrt(l2)*l2;

  if(iconst[0]) {
    H(0,0) = H(3,3) = -(l2-a[0]*a[0])/l3;
    H(1,1) = H(4,4) = -(l2-a[1]*a[1])/l3;
    H(2,2) = H(5,5) = -(l2-a[2]*a[2])/l3;

    H(0,1) = H(1,0) = H(3,4) = H(4,3) = a[0]*a[1]/l3;
    H(0,2) = H(2,0) = H(3,5) = H(5,3) = a[0]*a[2]/l3;
    H(1,2) = H(2,1) = H(4,5) = H(5,4) = a[1]*a[2]/l3;
  }
  else {
    H(0,0) = H(3,3) = (l2-a[0]*a[0])/l3;
    H(1,1) = H(4,4) = (l2-a[1]*a[1])/l3;
    H(2,2) = H(5,5) = (l2-a[2]*a[2])/l3;

    H(0,1) = H(1,0) = H(3,4) = H(4,3) = -a[0]*a[1]/l3;
    H(0,2) = H(2,0) = H(3,5) = H(5,3) = -a[0]*a[2]/l3;
    H(1,2) = H(2,1) = H(4,5) = H(5,4) = -a[1]*a[2]/l3;
  }

  H(0,3) = H(3,0) = -H(0,0);
  H(1,4) = H(4,1) = -H(1,1);
  H(2,5) = H(5,2) = -H(2,2);

  H(0,4) = H(4,0) = H(1,3) = H(3,1) = -H(0,1);
  H(0,5) = H(5,0) = H(2,3) = H(3,2) = -H(0,2);
  H(1,5) = H(5,1) = H(2,4) = H(4,2) = -H(1,2);

  return H;
}

} // namespace Simo

#endif
