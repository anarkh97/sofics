#ifndef _POINTVARIPLANESEGMENTDISTANCECONSTRAINTFUNCTION_H_
#define _POINTVARIPLANESEGMENTDISTANCECONSTRAINTFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/utilities.hpp>
#include <cmath>

namespace Simo {

template<typename Scalar>
class PointVariPlaneSegmentDistanceConstraintFunction : public ScalarValuedFunction<12,Scalar,17,1,double>
{
    // constrains the distance (d) between a point and a variable plane segment (defined by three points x1, x2 and x3) according to
    // "3D Distance from a Point to a Triangle" Mark W. Jones, Technical Report February 1995

    Eigen::Matrix<double,3,1> z0;
    Eigen::Matrix<double,3,1> z1;
    Eigen::Matrix<double,3,1> z2;
    Eigen::Matrix<double,3,1> z3;
    double A, omega, phase, B, C, dist;
    double thickness;
    bool negate;

  private:
    template<typename T>
    T distanceToSegment(Eigen::Matrix<T,3,1> x0, Eigen::Matrix<T,3,1> x1, Eigen::Matrix<T,3,1> x2, Eigen::Matrix<T,3,1> x3){

      T d0;

      // inward facing normals of tet formed by 4 points
      Eigen::Matrix<T,3,1> nplane = (x2-x1).cross(x3-x1).normalized();
      Eigen::Matrix<T,3,1> n1     = (x1-x2).cross(x0-x2).normalized();
      Eigen::Matrix<T,3,1> n2     = (x2-x3).cross(x0-x3).normalized();
      Eigen::Matrix<T,3,1> n3     = (x3-x1).cross(x0-x1).normalized();

      // angles between plane and sides of tet formed by 4 points
      T p1 = n1.dot(nplane);
      T p2 = n2.dot(nplane);
      T p3 = n3.dot(nplane);

      // this formulation assumes that the plane segment is an acute triangle
      
      if( p1 <= T(0.) && p2 <= T(0.) && p3 <= T(0.) ) { // if all angles are less than or equal to 90 degrees, then point projects to plane segment
	 // this is the only case for which contact can occure     

         d0 = (x2-x1).cross(x3-x1).normalized().dot(x0-x1) - thickness;
  
      } else {
        d0 = std::min((x0-x2).norm() - thickness,std::min((x0-x3).norm() - thickness,(x0-x1).norm() - thickness)); 
      } 

        /*else if (p1 > T(0.) && p2 > T(0.) && p3 <= T(0.)) { // if two angles are greater than 90, point is closest to that vertex
 
         d0 = (x0-x2).norm() - thickness;

      } else if (p1 <= T(0.) && p2 >= T(0.) && p3 >= T(0.)) { 

        d0 = (x0-x3).norm() - thickness;

      } else if (p1 >= T(0.) && p2 <= T(0.) && p3 >= T(0.)) {
  
        d0 = (x0-x1).norm() - thickness;

      } else if (p1 >= T(0.) && p2 <= T(0.) && p3 <= T(0.)) { // if one angle is greater than 90 degrees, then point is closest to that line segment or its end points
 
        if(p2 <= p3) { // closer to point opposite p2
          
          d0 = (x0-x1).dot(x2-x1);
          if(d0 <= T(0.)){
            d0 = (x0-x1).norm() - thickness;
          } else {
            d0 = sqrt((x0-x1).squaredNorm()-d0*d0) - thickness;
          }

        } else { // close to point opposite p3

          d0 = (x0-x2).dot(x1-x2);
          if(d0 <= T(0.)){
            d0 = (x0-x2).norm()  - thickness;
          } else {
            d0 = sqrt((x0-x2).squaredNorm()-d0*d0) - thickness;
          }

        }

      } else if (p1 <= T(0.) && p2 > T(0.) && p3 <= T(0.)) {

        if(p1 <= p3) { // closer to point opposite p1

          d0 = (x0-x3).dot(x2-x3);
          if(d0 <= T(0.)){
            d0 = (x0-x3).norm() - thickness;
          } else {
            d0 = sqrt((x0-x3).squaredNorm()-d0*d0) - thickness;
          }

        } else { // close to point opposite p3

          d0 = (x0-x2).dot(x3-x2);
          if(d0 <= T(0.)){
            d0 = (x0-x2).norm() - thickness;
          } else {
            d0 = sqrt((x0-x2).squaredNorm()-d0*d0) - thickness;
          }

        }

      } else if (p1 <= T(0.) && p2 <= T(0.) && p3 > T(0.)) {

        if(p2 <= p1) { // closer to point opposite p2

          d0 = (x0-x1).dot(x3-x1);
          if(d0 <= T(0.)){
            d0 = (x0-x1).norm() - thickness;
          } else {
            d0 = sqrt((x0-x1).squaredNorm()-d0*d0) - thickness;
          }

        } else { // close to point opposite p1

          d0 = (x0-x3).dot(x1-x3);
          if(d0 <= T(0.)){
            d0 = (x0-x3).norm() - thickness;
          } else {
            d0 = sqrt((x0-x3).squaredNorm()-d0*d0) - thickness;
          }

        }

      } else {
        
        d0 = (x0-x3).norm() - thickness;

      }*/

      return d0;

    }

public:
    PointVariPlaneSegmentDistanceConstraintFunction(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
        z0 = sconst.template segment<3>(0);
        z1 = sconst.template segment<3>(3);
        z2 = sconst.template segment<3>(6);
        z3 = sconst.template segment<3>(9);
        A = sconst[12];
        omega = sconst[13];
        phase = sconst[14];
        B = sconst[15];
        C = sconst[16];
        negate = bool(iconst[0]);

        double dummy = (z1-z3).norm();
        thickness = (z1-z2).norm();
        if (dummy > thickness)
            thickness = dummy;

        dummy = (z2-z3).norm();
        if(dummy > thickness)
            thickness = dummy;

        thickness *= 0.01;

        dist = distanceToSegment(z0,z1,z2,z3);
    }

    Scalar operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar t)
    {
      // q(0) = x translation of point 1
      // q(1) = y translation of point 1
      // q(2) = z translation of point 1
      // q(3) = x translation of point 2
      // q(4) = y translation of point 2
      // q(5) = z translation of point 2
      // q(6) = x translation of point 3
      // q(7) = y translation of point 3
      // q(8) = z translation of point 3
      // q(9) = x translation of point 4
      // q(10) = y translation of point 4
      // q(11) = z translation of point 4
      Eigen::Matrix<Scalar,3,1> z0 = PointVariPlaneSegmentDistanceConstraintFunction::z0.template cast<Scalar>() + q.template segment<3>(0);
      Eigen::Matrix<Scalar,3,1> z1 = PointVariPlaneSegmentDistanceConstraintFunction::z1.template cast<Scalar>() + q.template segment<3>(3);
      Eigen::Matrix<Scalar,3,1> z2 = PointVariPlaneSegmentDistanceConstraintFunction::z2.template cast<Scalar>() + q.template segment<3>(6);
      Eigen::Matrix<Scalar,3,1> z3 = PointVariPlaneSegmentDistanceConstraintFunction::z3.template cast<Scalar>() + q.template segment<3>(9);

      Scalar d = distanceToSegment(z0,z1,z2,z3);

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
