#ifndef _DOTTYPE3CONSTRAINTFUNCTION_H_
#define _DOTTYPE3CONSTRAINTFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

namespace Simo {

template<typename Scalar>
class DotType3ConstraintFunction : public ScalarValuedFunction<12,Scalar,10,0,double>
{
   Eigen::Matrix<double,3,1> a0,b0hat,c0hat;
   double d0;

  public:
    DotType3ConstraintFunction(const Eigen::Array<double,10,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      Eigen::Matrix<double,3,1> b0,c0;
      a0 << sconst(0), sconst(1), sconst(2);
      b0 << sconst(3), sconst(4), sconst(5);
      c0 << sconst(6), sconst(7), sconst(8);
      b0hat = b0.normalized();
      c0hat = c0.normalized();
      d0 = sconst(9);
    }

    Scalar operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar t)
    {
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = 1st component of axis-angle rotation vector of node 1
      // q[4] = 2nd component of axis-angle rotation vector of node 1
      // q[5] = 3rd component of axis-angle rotation vector of node 1
      // q[6] = x translation of node 2
      // q[7] = y translation of node 2
      // q[8] = z translation of node 2
      // q[9] = 1st component of axis-angle rotation vector of node 2
      // q[10] = 2nd component of axis-angle rotation vector of node 2
      // q[11] = 3rd component of axis-angle rotation vector of node 2

      // return value:
      // sum of scalar projection of a onto b, and the scalar projection of a onto c
      // where a = a0+u (i.e. translated), b = R_1*b0 (i.e. rotated) and c = R_2*c0
      // see: http://en.wikipedia.org/wiki/Scalar_projection

      Eigen::Matrix<Scalar,3,1> u1 = q.template segment<3>(0);
      Eigen::Quaternion<Scalar> z1;
      z1.setFromOneVector(q.template segment<3>(3));

      Eigen::Matrix<Scalar,3,1> u2 = q.template segment<3>(6);
      Eigen::Quaternion<Scalar> z2;
      z2.setFromOneVector(q.template segment<3>(9));

      // "unbiased" alternative to dot constraint type 2
      return -d0 + (a0.template cast<Scalar>() + u2 - u1).dot
             (z1.toRotationMatrix()*b0hat.template cast<Scalar>()+z2.toRotationMatrix()*c0hat.template cast<Scalar>());
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,1,12>
Jacobian<double,DotType3ConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, double t);

template<>
Eigen::Matrix<double,12,12>
Hessian<double,DotType3ConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, double t);

} // namespace Simo

#endif
