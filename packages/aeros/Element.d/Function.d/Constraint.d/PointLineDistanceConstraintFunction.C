#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/PointLineDistanceConstraintFunction.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for point-
// line distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,3>
Jacobian<double,PointLineDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix<double,3,1> J;

  Eigen::Vector3d x0 = sconst.segment<3>(0);                                       
  Eigen::Vector3d x1 = sconst.segment<3>(3);
  Eigen::Vector3d x2 = sconst.segment<3>(6);

  // node's current coordinates
  x0 += q;

  Eigen::Vector3d a = x2-x1;
  Eigen::Vector3d b = x1-x0;
  Eigen::Vector3d c = a.cross(b);

  if(iconst[0])
    J = -a.cross(c)/(a.norm()*c.norm());
  else 
    J = a.cross(c)/(a.norm()*c.norm());

  return J.transpose();
}

// specializing the member function template of constraint hessian operator for point-
// line distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,3,3>
Hessian<double,PointLineDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix<double,3,3> H;

  Eigen::Vector3d x0 = sconst.segment<3>(0); 
  Eigen::Vector3d x1 = sconst.segment<3>(3);
  Eigen::Vector3d x2 = sconst.segment<3>(6);

  // node's current coordinates
  x0 += q;

  Eigen::Vector3d a = x2-x1;
  Eigen::Vector3d b = x1-x0;
  Eigen::Vector3d c = a.cross(b);
  double numerator = c.norm();
  double denominator = a.norm();

  // d = sqrt( (a[1]*b[2] - a[2]*b[1])^2 + (a[2]*b[0] - a[0]*b[2])^2 + (a[0]*b[1] - a[1]*b[0])^2 ) / denominator
  //   = sqrt( (a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))^2 +
  //           (a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))^2 +
  //           (a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0]))^2 ) / denominator

  // 1/2 of the first derivatives of numerator^2:
  // a[1]*c[2]-a[2]*c[1] = a[1]*(a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0])) 
  //                      -a[2]*(a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))    
  // a[2]*c[0]-a[0]*c[2] = a[2]*(a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))
  //                      -a[0]*(a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0]))
  // a[0]*c[1]-a[1]*c[0] = a[0]*(a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))
  //                      -a[1]*(a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))
  double dx = a[1]*c[2]-a[2]*c[1];
  double dy = a[2]*c[0]-a[0]*c[2];
  double dz = a[0]*c[1]-a[1]*c[0];

  // 1/2 of the second derivatives of numerator^2
  double dxx = a[1]*a[1]+a[2]*a[2];
  double dyy = a[0]*a[0]+a[2]*a[2];
  double dzz = a[0]*a[0]+a[1]*a[1];
  double dxy = -a[1]*a[0];
  double dxz = -a[2]*a[0];
  double dyz = -a[2]*a[1];

  // second derivatives
  double n2 = numerator*numerator;
  double dn3 = denominator*numerator*numerator*numerator;
  if(iconst[0]) {
    H(0,0) = -(n2*dxx - dx*dx)/dn3;
    H(1,1) = -(n2*dyy - dy*dy)/dn3;
    H(2,2) = -(n2*dzz - dz*dz)/dn3;

    H(0,1) = H(1,0) = -(n2*dxy - dx*dy)/dn3;
    H(0,2) = H(2,0) = -(n2*dxz - dx*dz)/dn3;
    H(1,2) = H(2,1) = -(n2*dyz - dy*dz)/dn3;
  }
  else {
    H(0,0) = (n2*dxx - dx*dx)/dn3;
    H(1,1) = (n2*dyy - dy*dy)/dn3;
    H(2,2) = (n2*dzz - dz*dz)/dn3;

    H(0,1) = H(1,0) = (n2*dxy - dx*dy)/dn3;
    H(0,2) = H(2,0) = (n2*dxz - dx*dz)/dn3;
    H(1,2) = H(2,1) = (n2*dyz - dy*dz)/dn3;
  }

  return H;
}

} // namespace Simo

#endif
