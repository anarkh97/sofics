#ifdef USE_EIGEN3
#include <Element.d/Function.d/ExternalForce.d/FollowerMomentForceFunction.h>

namespace Simo {

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,FollowerMomentForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double)
{
  Eigen::Matrix<double,3,3> J;

  Eigen::Matrix<double,3,1> m0;   // applied moment in the initial configuration
  Eigen::Matrix<double,3,3> Rref; // nodal rotation in the reference configuration

  m0   << sconst[0], sconst[1], sconst[2];
  Rref << sconst[3], sconst[4], sconst[5],
          sconst[6], sconst[7], sconst[8],
          sconst[9], sconst[10], sconst[11];

  Eigen::Matrix<double,3,1> m = Rref*m0;

  if((q.array() == 0).all()) {
    // the jacobian is 0.5*skew(m)
    J <<     0, -m[2],  m[1],
          m[2],     0, -m[0],
         -m[1],  m[0],     0;
    J *= 0.5;
  }
  else {
    Eigen::Matrix<double,3,3> C2;
    directional_deriv2(q, m, C2);
    J = -C2;
  }

  return J;
}

} // namespace Simo

#endif
