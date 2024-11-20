#ifndef _POINTLINEDISTANCECONSTRAINTFUNCTION_H_
#define _POINTLINEDISTANCECONSTRAINTFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <cmath>

namespace Simo {

template<typename Scalar>
class PointLineDistanceConstraintFunction : public ScalarValuedFunction<3,Scalar,14,1,double>
{
    // constrains the distance (d) between a point (x0) and a fixed line (defined by two points x1 and x2) according to
    // d - (A*sin(omega*t+phi) + (B-C*t)*d0) = 0, <= 0 or >= 0
    // see: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    
    Eigen::Matrix<double,3,1> x0;
    Eigen::Matrix<double,3,1> x1;
    Eigen::Matrix<double,3,1> x2;
    double A, omega, phase, B, C, d0;
    bool negate;

  public:
    PointLineDistanceConstraintFunction(const Eigen::Array<double,14,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      x0 = sconst.segment<3>(0);
      x1 = sconst.segment<3>(3);
      x2 = sconst.segment<3>(6);
      A = sconst[9];
      omega = sconst[10];
      phase = sconst[11];
      B = sconst[12];
      C = sconst[13];
      negate = bool(iconst[0]);
      d0 = (x0-x1).cross(x0-x2).norm()/(x2-x1).norm();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // q[0] = x translation of point 1
      // q[1] = y translation of point 1
      // q[2] = z translation of point 1
      Eigen::Matrix<Scalar,3,1> x0 = PointLineDistanceConstraintFunction::x0.template cast<Scalar>() + q;
      Eigen::Matrix<Scalar,3,1> x1 = PointLineDistanceConstraintFunction::x1.template cast<Scalar>();
      Eigen::Matrix<Scalar,3,1> x2 = PointLineDistanceConstraintFunction::x2.template cast<Scalar>();

      Scalar d = (x0-x1).cross(x0-x2).norm()/(x2-x1).norm();
      using std::sin;
      Scalar f = d -(A*sin(omega*t+phase) + (B-C*t)*d0);
      if(negate) return -f; else return f;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,1,3>
Jacobian<double,PointLineDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t);

template<>
Eigen::Matrix<double,3,3>
Hessian<double,PointLineDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t);

} // namespace Simo

#endif
