#ifndef _POINTVARIPLANESEGMENTDISTANCECONSTRAINTFUNCTION2_H_
#define _POINTVARIPLANESEGMENTDISTANCECONSTRAINTFUNCTION2_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/utilities.hpp>
#include <cmath>

namespace Simo {

template<typename Scalar>
class PointVariPlaneSegmentDistanceConstraintFunction2 : public ScalarValuedFunction<12,Scalar,17,1,double>
{
    // constrains the distance (d) between a point and a variable plane segment (defined by three points x1, x2 and x3) according to
    // "3D Distance from a Point to a Triangle" Mark W. Jones, Technical Report February 1995

    Eigen::Matrix<double,3,1> x0;
    Eigen::Matrix<double,3,1> x1;
    Eigen::Matrix<double,3,1> x2;
    Eigen::Matrix<double,3,1> x3;
    double A, omega, phase, B, C, dist;
    double thickness;
    bool negate;

  private: 
    template<typename T>
    T distanceToSegment(Eigen::Matrix<T,3,1> x0, Eigen::Matrix<T,3,1> x1, Eigen::Matrix<T,3,1> x2, Eigen::Matrix<T,3,1> x3) {

      Eigen::Matrix<T,3,1> n = (x2-x1).cross(x3-x1).normalized();

      T d0;
      T d0_p = n.dot(x0-x1); // initial distance
      T d0_21 = n.cross(x2-x1).normalized().dot(x0-x1);
      T d0_13 = n.cross(x1-x3).normalized().dot(x0-x3);
      T d0_32 = n.cross(x3-x2).normalized().dot(x0-x2);

      if(d0_21 > T(0.0) && d0_13 < T(0.0) && d0_32 < T(0.0)) { // closest to node 3
        d0 = sqrt((x0-x3).dot(x0-x3)) - thickness;
      }
      else if(d0_13 > T(0.0) && d0_21 < T(0.0) && d0_32 < T(0.0)) { // closest to node 2
        d0 = sqrt((x0-x2).dot(x0-x2)) - thickness;
      }
      else if(d0_32 > T(0.0) && d0_21 < T(0.0) && d0_13 < T(0.0)) { // closest to node 1
        d0 = sqrt((x0-x1).dot(x0-x1)) - thickness;
      }
      else { // closest to edge or face
        if(d0_21 > T(0.0))
          d0_21 = T(0.0);
        if(d0_13 > T(0.0))
          d0_13 = T(0.0);
        if(d0_32 > T(0.0))
          d0_32 = T(0.0);

        d0 = sqrt(d0_p * d0_p + d0_21 * d0_21 + d0_13 * d0_13 + d0_32 * d0_32) - thickness;
      }

      return d0;
    }

  public:
    PointVariPlaneSegmentDistanceConstraintFunction2(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
        x0 = sconst.template segment<3>(0);
        x1 = sconst.template segment<3>(3);
        x2 = sconst.template segment<3>(6);
        x3 = sconst.template segment<3>(9);
        A = sconst[12];
        omega = sconst[13];
        phase = sconst[14];
        B = sconst[15];
        C = sconst[16];
        negate = bool(iconst[0]);

        double dummy = (x1-x3).norm();
        thickness = (x1-x2).norm();
        if (dummy > thickness)
            thickness = dummy;

        dummy = (x2-x3).norm();
        if(dummy > thickness)
            thickness = dummy;

        thickness *= 0.01;

      dist = distanceToSegment(x0,x1,x2,x3);
    }

    Scalar operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar t)
    {
      // q(0)  = x translation of point 1
      // q(1)  = y translation of point 1
      // q(2)  = z translation of point 1
      // q(3)  = x translation of point 2
      // q(4)  = y translation of point 2
      // q(5)  = z translation of point 2
      // q(6)  = x translation of point 3
      // q(7)  = y translation of point 3
      // q(8)  = z translation of point 3
      // q(9)  = x translation of point 4
      // q(10) = y translation of point 4
      // q(11) = z translation of point 4
      Eigen::Matrix<Scalar,3,1> x0 = PointVariPlaneSegmentDistanceConstraintFunction2::x0.template cast<Scalar>() + q.template segment<3>(0);
      Eigen::Matrix<Scalar,3,1> x1 = PointVariPlaneSegmentDistanceConstraintFunction2::x1.template cast<Scalar>() + q.template segment<3>(3);
      Eigen::Matrix<Scalar,3,1> x2 = PointVariPlaneSegmentDistanceConstraintFunction2::x2.template cast<Scalar>() + q.template segment<3>(6);
      Eigen::Matrix<Scalar,3,1> x3 = PointVariPlaneSegmentDistanceConstraintFunction2::x3.template cast<Scalar>() + q.template segment<3>(9);

      Scalar d = distanceToSegment(x0,x1,x2,x3);

      // note: the sign of d is positive if x0 is on the same side of the plane as the normal
      //       and negative if it is on the opposite side
      using std::sin;
      Scalar f = d -(A*sin(omega*t+phase) + (B-C*t)*dist);
      if(negate) return -f; else return f;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
