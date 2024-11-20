#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/InertialType2ForceFunction.h>

namespace Simo {

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,InertialType2ForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix3d T, Tdot, C1a, C1b, C2, C4, C5;
  Eigen::Vector3d inc_displacement, A, V, Omega, P, F;

  Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > J(sconst.data()+0);
  Eigen::Map<const Eigen::Vector3d> A_n(sconst.data()+9), V_n(sconst.data()+12), Psi_n(sconst.data()+15);
  const double &beta   = sconst[18];
  const double &gamma  = sconst[19];
  const double &alphaf = sconst[20];
  const double &alpham = sconst[21];
  const double &dt     = sconst[22];

  double s1 = gamma/(dt*beta);
  double s2 = (1-alpham)/(dt*dt*beta*(1-alphaf));

  inc_displacement = q - Psi_n;
  V = s1*inc_displacement + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
  A = s2*inc_displacement + -(1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;

  tangential_transf(q, T);
  tangential_transf_dot(q, V, Tdot);
  Omega = T*V;
  P = J*Omega;
  F = J*(T*A + Tdot*V) + Omega.cross(P);
  directional_deriv1(q, A, C1a);
  directional_deriv1(q, V, C1b);
  directional_deriv2(q, F, C2);
  directional_deriv4(q, V, C4);
  directional_deriv5(q, V, C5);

  return C2 + T.transpose()*(J*(C1a+s2*T + C5+s1*C4) + (skew(Omega)*J - skew(P))*(C1b+s1*T));
}

} // namespace Simo
#endif
