#ifdef USE_EIGEN3

#include <Element.d/Function.d/Constraint.d/PointVariPlaneDistanceConstraintFunction.h>
#include <Element.d/Function.d/Constraint.d/pt2plane_d1.cc> // For the Jacobian
#include <Element.d/Function.d/Constraint.d/pt2plane_d2.cc> // For the Hessian

namespace Simo {

  // specializing the member function template of constraint jacobian operator for varying 
  // point-plane distance constraint function with double precision scalar
  template<>
  Eigen::Matrix<double,1,12>
  Jacobian<double,PointVariPlaneDistanceConstraintFunction>
  ::operator() (const Eigen::Matrix<double,12,1>& q, double t)
  {
    Eigen::Matrix<double,1,12> J;

    Eigen::Matrix<double,1,12> x;
    x.segment<3>(0) = sconst.segment<3>(0); x.segment<3>(0) += q.segment<3>(0);
    x.segment<3>(3) = sconst.segment<3>(3); x.segment<3>(3) += q.segment<3>(3);
    x.segment<3>(6) = sconst.segment<3>(6); x.segment<3>(6) += q.segment<3>(6);
    x.segment<3>(9) = sconst.segment<3>(9); x.segment<3>(9) += q.segment<3>(9);

    // Maple generated Jacobian function
    pt2plane_d1(x.data(), J.data());

    if(iconst[0])
      J = -J;
    else
      J = J;

    return J;
  }

  // specializing the member function template of constraint hessian operator for varying
  // point-plane distance constraint function with double precision scalar
  template<>
  Eigen::Matrix<double,12,12>
  Hessian<double,PointVariPlaneDistanceConstraintFunction>
  ::operator() (const Eigen::Matrix<double,12,1>& q, double t)
  {
    Eigen::Matrix<double,12,12> H;

    Eigen::Matrix<double,1,12> x;
    x.segment<3>(0) = sconst.segment<3>(0); x.segment<3>(0) += q.segment<3>(0);
    x.segment<3>(3) = sconst.segment<3>(3); x.segment<3>(3) += q.segment<3>(3);
    x.segment<3>(6) = sconst.segment<3>(6); x.segment<3>(6) += q.segment<3>(6);
    x.segment<3>(9) = sconst.segment<3>(9); x.segment<3>(9) += q.segment<3>(9);

    // Maple generated Hessian function
    pt2plane_d2(x.data(), H.data());

    if(iconst[0])
      H = -H;
    else
      H = H;

    return H.transpose();
    }

} // namespace Simo

#endif
