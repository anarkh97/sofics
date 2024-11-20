#ifndef _POINTPOINTDISTANCECONSTRAINTFUNCTION_H_
#define _POINTPOINTDISTANCECONSTRAINTFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <cmath>

namespace Simo {

template<typename Scalar>
class PointPointDistanceConstraintFunction : public ScalarValuedFunction<3,Scalar,11,1,double>
{
    // constrains the distance (d) between a point (x0) and a fixed point (x1) according to
    // d - (A*sin(omega*t+phi) + (B-C*t)*d0) = 0, <= 0 or >= 0
    
    Eigen::Matrix<double,3,1> x0;
    Eigen::Matrix<double,3,1> x1;
    double A, omega, phase, B, C, d0;
    bool negate;

  public:
    PointPointDistanceConstraintFunction(const Eigen::Array<double,11,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      x0 = sconst.segment<3>(0);
      x1 = sconst.segment<3>(3);
      A = sconst[6];
      omega = sconst[7];
      phase = sconst[8];
      B = sconst[9];
      C = sconst[10];
      negate = bool(iconst[0]);

      d0 = (x0-x1).norm();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // q(0) = x translation of point 1
      // q(1) = y translation of point 1
      // q(2) = z translation of point 1
      Eigen::Matrix<Scalar,3,1> x0 = PointPointDistanceConstraintFunction::x0.template cast<Scalar>() + q;
      Eigen::Matrix<Scalar,3,1> x1 = PointPointDistanceConstraintFunction::x1.template cast<Scalar>();

      Scalar d = (x0-x1).norm();
      using std::sin;
      Scalar f = d -(A*sin(omega*t+phase) + (B-C*t)*d0);
      if(negate) return -f; else return f;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,1,3>
Jacobian<double,PointPointDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t);

template<>
Eigen::Matrix<double,3,3>
Hessian<double,PointPointDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t);

} // namespace Simo

#endif
