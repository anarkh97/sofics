#ifndef _INERTIALTYPE2FORCEFUNCTION_H_
#define _INERTIALTYPE2FORCEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>

// Compute the inertial moment for a rigid body with 3 rotational degrees of freedom,
// according to the classical Euler's equation for the rigid body and the generalized-alpha 
// method applied to the total rotation vector and it's first and second
// time-derivatives.

namespace Simo {

template<typename Scalar>
class InertialType2ForceFunction : public VectorValuedFunction<3,3,Scalar,23,0,double>
{
  public:
    Eigen::Matrix<double,3,3> J;            // inertia tensor
    Eigen::Matrix<double,3,1> A_n, V_n;     // first and second time-derivatives of the total rotation vector at start of current time-step
    Eigen::Matrix<double,3,1> Psi_n;        // total rotation vector at start of current time-step
    double beta, gamma, alphaf, alpham, dt; // time integration scheme parameters

  public:
    InertialType2ForceFunction(const Eigen::Array<double,23,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      J = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
      A_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+9);
      V_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+12);
      Psi_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+15);
      beta   = sconst[18];
      gamma  = sconst[19];
      alphaf = sconst[20];
      alpham = sconst[21];
      dt     = sconst[22];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // inputs:
      // q[0] = x component of total axis-angle rotation vector
      // q[1] = y component of total axis-angle rotation vector
      // q[2] = z component of total axis-angle rotation vector

      // compute the increment in the total rotation vector
      Eigen::Matrix<Scalar,3,1> inc_displacement = q - Psi_n.template cast<Scalar>();

      // compute the updated first and second time-derivatives of the total rotation vector
      Eigen::Matrix<Scalar,3,1> V, A;
      V = gamma/(dt*beta)*inc_displacement
          + ((1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n).template cast<Scalar>();
      A = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement
          + (-(1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n).template cast<Scalar>();

      // compute the tangential transformation matrix T(q) and it's time-derivative Tdot(q, V)
      Eigen::Matrix<Scalar,3,3> T, Tdot;
      tangential_transf(q, T);
      tangential_transf_dot(q, V, Tdot);

      // compute the convected angular acceleration
      Eigen::Matrix<Scalar,3,1> Omega = T*V;

      // transformed angular momentum balance equation
      return T.transpose()*(J*(T*A + Tdot*V) + Omega.cross(J*Omega));
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,InertialType2ForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t);

} // namespace Simo

#endif
