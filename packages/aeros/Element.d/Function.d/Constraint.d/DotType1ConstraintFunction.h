#ifndef _DOTTYPE1CONSTRAINTFUNCTION_H_
#define _DOTTYPE1CONSTRAINTFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

namespace Simo {

template<typename Scalar>
class DotType1ConstraintFunction : public ScalarValuedFunction<6,Scalar,7,0,double>
{
   // a0 and b0 are vectors which can be interpreted as one selected axis of each of the body-attached frames attached 
   // to two nodes in some specified configuration.
   // note that the specified configuration can be either the undeformed configuration (for total lagrangian),
   // or the reference configuration at t_n (for updated lagrangian), or the current configuration at t_{n+1-alphaf}^k
   // (where k is the newton iteration) for eulerian.
   Eigen::Matrix<double,3,1> a0, b0hat;
   double d0;

  public:
    DotType1ConstraintFunction(const Eigen::Array<double,7,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      Eigen::Matrix<double,3,1> b0;
      a0 << sconst(0), sconst(1), sconst(2);
      b0 << sconst(3), sconst(4), sconst(5);
      d0 = sconst(6);
      b0hat = b0.normalized();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar)
    {
      // inputs:
      // q[0] = 1st component of axis-angle rotation vector of node 1
      // q[1] = 2nd component of axis-angle rotation vector of node 1
      // q[2] = 3rd component of axis-angle rotation vector of node 1
      // q[3] = 1st component of axis-angle rotation vector of node 2
      // q[4] = 2nd component of axis-angle rotation vector of node 2
      // q[5] = 3rd component of axis-angle rotation vector of node 2
      // see: http://en.wikipedia.org/wiki/Axis-angle_representation
      // note that the rotation vector can either be the total increment from the undeformed configuration (for total
      // total lagrangian), or the increment from the reference configuration at t_n (for updated lagrangian), or the 
      // spin (i.e. zero/infintesimal increment) from the current configuration at t_{n+1-alphaf}^k (where k is the
      // newton iteration) for eulerian.

      // return value:
      // dot product of rotated axes (R(theta0)*a0).dot(R(thetb0)*b0)

      Eigen::Quaternion<Scalar> z1;
      z1.setFromOneVector(q.template segment<3>(0));

      Eigen::Quaternion<Scalar> z2;
      z2.setFromOneVector(q.template segment<3>(3));

      return -d0 + (z1.toRotationMatrix()*a0.template cast<Scalar>()).dot
                   (z2.toRotationMatrix()*b0hat.template cast<Scalar>());
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,1,6>
Jacobian<double,DotType1ConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t);

template<>
Eigen::Matrix<double,6,6>
Hessian<double,DotType1ConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t);

} // namespace Simo

#endif
