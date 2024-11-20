#ifdef USE_EIGEN3
#include <Element.d/Function.d/ExternalForce.d/FollowerForceFunction.h>

namespace Simo {

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,FollowerForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double)
{
  Eigen::Matrix<double,3,3> J;

  Eigen::Matrix<double,3,1> f0;   // applied force in the initial configuration
  Eigen::Matrix<double,3,3> Rref; // nodal rotation in the reference configuration

  f0   << sconst[0], sconst[1], sconst[2];
  Rref << sconst[3], sconst[4], sconst[5],
          sconst[6], sconst[7], sconst[8],
          sconst[9], sconst[10], sconst[11];

  Eigen::Matrix<double,3,1> f = Rref*f0;

  if((q.array() == 0).all()) {
    // the jacobian is skew(f)
    J <<     0, -f[2],  f[1],
          f[2],     0, -f[0],
         -f[1],  f[0],     0;
  }
  else {
    std::cerr << "Jacobian<double,FollowerForceFunction>::operator() not implemented for this case\n";
    J.setZero();
  }

  return J;
}

} // namespace Simo

#endif
