#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/InertialType1ForceFunction.h>
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <iostream>

namespace Simo {

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,InertialType1ForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Vector3d A, V, Psi, F, P;

  Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > J(sconst.data()+0), R_n(sconst.data()+15), Rref(sconst.data()+24);
  Eigen::Map<const Eigen::Vector3d> A_n(sconst.data()+9), V_n(sconst.data()+12);
  const double &beta   = sconst[33];
  const double &gamma  = sconst[34];
  const double &alphaf = sconst[35];
  const double &alpham = sconst[36];
  const double &dt     = sconst[37];

  double s1 = gamma/(dt*beta);
  double s2 = (1-alpham)/(dt*dt*beta*(1-alphaf));
  if( (q.array() == 0).all() ) {
    mat_to_vec<double>(R_n.transpose()*Rref, Psi);
  }
  else {
    Eigen::Matrix3d R;
    vec_to_mat(q, R);
    mat_to_vec<double>(R_n.transpose()*R*Rref, Psi);
  }
  V = s1*Psi + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
  A = s2*Psi - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
  F = Rref*(J*A + V.cross(J*V));
  P = J*V;

  Eigen::Array<double,18,1> sconst2 = sconst.segment<18>(15);
  Jacobian<double,IncrementalRotationVector> dPsidq(sconst2, Eigen::Array<int,0,1>::Zero());

  if( (q.array() == 0).all() ) {
    return -0.5*skew(F) + Rref*(s2*J + s1*(skew(V)*J - skew(P)))*dPsidq(q, t);
  }
  else {
    Eigen::Matrix3d T, C2;
    tangential_transf(q, T);
    directional_deriv2(q, F, C2);
    return C2 + T.transpose()*Rref*(s2*J + s1*(skew(V)*J - skew(P)))*dPsidq(q, t);
  }
}

} // namespace Simo
#endif
