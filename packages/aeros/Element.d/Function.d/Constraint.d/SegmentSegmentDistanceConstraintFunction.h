#ifndef _SEGMENTSEGMENTDISTANCECONSTRAINTFUNCTION_H_
#define _SEGMENTSEGMENTDISTANCECONSTRAINTFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <cmath>

namespace Simo {

template<typename Scalar>
class SegmentSegmentDistanceConstraintFunction : public ScalarValuedFunction<6,Scalar,17,1,double>
{
    // constrains the distance (d) between two line segments (defined by points p0,p1 and q0,q1 respectively)
    // according to d - (A*sin(omega*t+phi) + (B-C*t)*d0) = 0, <= 0 or >= 0
    // see: http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment
    // see: http://mathworld.wolfram.com/Line-LineDistance.html 
    Eigen::Matrix<double,3,1> p0;
    Eigen::Matrix<double,3,1> p1;
    Eigen::Matrix<double,3,1> q0;
    Eigen::Matrix<double,3,1> q1;
    double A, omega, phase, B, C, d0;
    bool negate;

  public:
    SegmentSegmentDistanceConstraintFunction(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      p0 = sconst.segment<3>(0);
      p1 = sconst.segment<3>(3);
      q0 = sconst.segment<3>(6);
      q1 = sconst.segment<3>(9);
      A = sconst[12];
      omega = sconst[13];
      phase = sconst[14];
      B = sconst[15];
      C = sconst[16];
      negate = bool(iconst[0]);

      Eigen::Matrix<double,3,1> u = p1 - p0;
      Eigen::Matrix<double,3,1> v = q1 - q0;
      Eigen::Matrix<double,3,1> w = p0 - q0;
      double a = u.dot(u);
      double b = u.dot(v);
      double c = v.dot(v);
      double d = u.dot(w);
      double e = v.dot(w);
      double D = a*c - b*b;
      double sc, sN, sD = D;
      double tc, tN, tD = D;

      if(D < 1.0e-6) {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
      }
      else {
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if(sN < 0.0) {
          sN = 0.0;
          tN = e;
          tD = c;
        }
        else if(sN > sD) {
          sN = sD;
          tN = e + b;
          tD = c;
        }
      }
      if(tN < 0.0) {
        tN = 0.0;
        if(-d < 0.0) {
          sN = 0.0;
        }
        else if(-d > a) {
          sN = sD;
        }
        else {
          sN = -d;
          sD = a;
        }
      }
      else if(tN > tD) {
        tN = tD;
        if((-d + b) < 0.0) {
          sN = 0;
        }
        else if((-d + b) > a) {
          sN = sD;
        }
        else {
          sN = (-d + b);
          sD = a;
        }
      }
      using std::abs;
      sc = ((abs(sN) < 1.0e-6) ? 0.0 : (sN/sD));
      tc = ((abs(tN) < 1.0e-6) ? 0.0 : (tN/tD));

      d0  = (w + (sc*u) - (tc*v)).norm();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar t)
    {
      // q[0] = x translation of point 1
      // q[1] = y translation of point 1
      // q[2] = z translation of point 1
      // q[3] = x translation of point 2
      // q[4] = y translation of point 2
      // q[5] = z translation of point 2
      Eigen::Matrix<Scalar,3,1> p0 = SegmentSegmentDistanceConstraintFunction::p0.template cast<Scalar>() + q.template segment<3>(0);
      Eigen::Matrix<Scalar,3,1> p1 = SegmentSegmentDistanceConstraintFunction::p1.template cast<Scalar>() + q.template segment<3>(3);
      Eigen::Matrix<Scalar,3,1> q0 = SegmentSegmentDistanceConstraintFunction::q0.template cast<Scalar>();
      Eigen::Matrix<Scalar,3,1> q1 = SegmentSegmentDistanceConstraintFunction::q1.template cast<Scalar>();

      Eigen::Matrix<Scalar,3,1> u = p1 - p0;
      Eigen::Matrix<Scalar,3,1> v = q1 - q0;
      Eigen::Matrix<Scalar,3,1> w = p0 - q0;
      Scalar a = u.dot(u);
      Scalar b = u.dot(v);
      Scalar c = v.dot(v);
      Scalar d = u.dot(w);
      Scalar e = v.dot(w);
      Scalar D = a*c - b*b;
      Scalar sc, sN, sD = D;
      Scalar tc, tN, tD = D;
        
      if(D < 1.0e-6) {
        sN = (Scalar) 0.0;
        sD = (Scalar) 1.0;
        tN = e;
        tD = c;
      }
      else {
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if(sN < 0.0) {
          sN = (Scalar) 0.0;
          tN = e;
          tD = c;
        }
        else if(sN > sD) {
          sN = sD;
          tN = e + b;
          tD = c;
        }
      }
      if(tN < 0.0) {
        tN = (Scalar) 0.0;
        if(-d < 0.0) {
          sN = (Scalar) 0.0;
        }
        else if(-d > a) {
          sN = sD;
        }
        else {
          sN = -d;
          sD = a;
        }
      }
      else if(tN > tD) {
        tN = tD;
        if((-d + b) < 0.0) {
          sN = (Scalar) 0;
        }
        else if((-d + b) > a) {
          sN = sD;
        }
        else {
          sN = (-d + b);
          sD = a;
        }
      }
      using std::abs;
      if((abs(sN) < 1.0e-6)) {
        sc = (Scalar) 0.0;
      }
      else {
        sc = sN/sD;
      }
      if(abs(tN) < 1.0e-6) {
        tc = (Scalar) 0.0;
      }
      else {
        tc = tN/tD;
      }
      //sc = ((abs(sN) < 1.0e-6) ? (Scalar) 0.0 : (sN/sD));
      //tc = ((abs(tN) < 1.0e-6) ? (Scalar) 0.0 : (tN/tD));

      Scalar dist  = (w + (sc*u) - (tc*v)).norm();
         
      using std::sin;
      Scalar f = dist -(A*sin(omega*t+phase) + (B-C*t)*d0);
      if(negate) return -f; else return f;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
