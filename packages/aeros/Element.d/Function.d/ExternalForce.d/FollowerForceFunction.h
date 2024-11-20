#ifndef _FOLLOWERFORCEFUNCTION_H_
#define _FOLLOWERFORCEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>

namespace Simo {

template<typename Scalar>
class FollowerForceFunction : public VectorValuedFunction<3,3,Scalar,12,0,double>
{
    Eigen::Matrix<double,3,1> f0;   // applied force in the initial configuration
    Eigen::Matrix<double,3,3> Rref; // nodal rotation in the reference configuration

  public:
    FollowerForceFunction(const Eigen::Array<double,12,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      f0  << sconst[0], sconst[1], sconst[2];
      Rref << sconst[3], sconst[4], sconst[5],
              sconst[6], sconst[7], sconst[8],
              sconst[9], sconst[10], sconst[11];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // q[0] = x component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[1] = y component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[2] = z component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1

      Eigen::Matrix<Scalar,3,3> incR;
      vec_to_mat<Scalar>(q, incR);
      return -incR*Rref*f0;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,FollowerForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double);

} // namespace Simo

#endif
