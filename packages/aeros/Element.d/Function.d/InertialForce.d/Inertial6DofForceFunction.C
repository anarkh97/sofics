#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/Inertial6DofForceFunction.h>
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <iostream>

namespace Simo {

template<>
Eigen::Matrix<double,6,6>
Jacobian<double,Inertial6DofForceFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t)
{
  Eigen::Vector3d F, G, P, Q;
  Eigen::Matrix3d J, J0, C, E;
  Eigen::Matrix<double,6,1> a, v, inc_displacement;

  const double &m = sconst[0];
  Eigen::Map<const Eigen::Matrix<double,6,1> > j(sconst.data()+1), a_n(sconst.data()+10), v_n(sconst.data()+16);
  Eigen::Map<const Eigen::Vector3d> c(sconst.data()+7), u_n(sconst.data()+22);
  Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > R_n(sconst.data()+25), Rref(sconst.data()+34);
  const double &beta   = sconst[43];
  const double &gamma  = sconst[44];
  const double &alphaf = sconst[45];
  const double &alpham = sconst[46];
  const double &dt     = sconst[47];

  J << j[0],  j[5],  j[4],
       j[5],  j[1],  j[3],
       j[4],  j[3],  j[2];
  C <<  0.0, -c[2],  c[1],
       c[2],   0.0, -c[0],
      -c[1],  c[0],   0.0;
  J0 = J - m*C*C;

  double s1 = gamma/(dt*beta);
  double s2 = (1-alpham)/(dt*dt*beta*(1-alphaf));

  inc_displacement.head<3>() = q.head<3>() - u_n;
  if( (q.tail<3>().array() == 0).all() ) {
    mat_to_vec<double>(R_n.transpose()*Rref, inc_displacement.tail<3>());
  }
  else {
    Eigen::Matrix3d dR;
    vec_to_mat<double>(q.tail<3>(), dR);
    mat_to_vec<double>(R_n.transpose()*dR*Rref, inc_displacement.tail<3>());
  }

  v = s1*inc_displacement + (1-(1-alphaf)*gamma/beta)*v_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*a_n;
  a = s2*inc_displacement - (1-alpham)/(dt*beta)*v_n + ((alpham-1)/(2*beta)+1)*a_n;

  Eigen::VectorBlock<Eigen::Matrix<double,6,1>,3> udot  = v.head<3>(), Omega = v.tail<3>(), 
                                                  uddot = a.head<3>(), Alpha = a.tail<3>(); 
  E = Rref*C*Rref.transpose();
  F = m*E*uddot + Rref*(J0*Alpha + Omega.cross(J0*Omega));
  G = Rref*(C*Alpha + Omega.cross(C*Omega));
  P = J0*Omega;
  Q = C*Omega;

  Eigen::Array<double,18,1> sconst2 = sconst.segment<18>(25);
  Jacobian<double,IncrementalRotationVector> dPsidq(sconst2, Eigen::Array<int,0,1>::Zero());
  Eigen::Matrix3d D = dPsidq(q.tail<3>(), t);

  Eigen::Matrix<double,6,6> K;
  if( (q.tail<3>().array() == 0).all() ) {
    K.topLeftCorner<3,3>()     = s2*m*Eigen::Matrix3d::Identity();
    K.topRightCorner<3,3>()    = -m*(-skew(G) + Rref*(s2*C + s1*(skew(Omega)*C - skew(Q)))*D);
    K.bottomLeftCorner<3,3>()  = s2*m*E;
    K.bottomRightCorner<3,3>() = -0.5*skew(F) + m*E*skew(uddot) + Rref*(s2*J0 + s1*(skew(Omega)*J0 - skew(P)))*D;
  }
  else {
    std::cerr << "Error: Jacobian<double,Inertial6DofForceFunction>::operator() is not implemented for Δψ != 0\n";
  }

  return K;
}

} // namespace Simo
#endif
