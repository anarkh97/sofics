#ifndef _INERTIAL6DOFFORCEFUNCTION_H_
#define _INERTIAL6DOFFORCEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>

namespace Simo {

// Compute the inertial force and moment for a rigid body with 6 degrees of freedom and offset,
// according to the classical Newton-Euler equation for the rigid body and the extended Newmark 
// method (based on ALGO_1 from first reference below)
// Ref: UNCONDITIONALLY STABLE ALGORITHMS FOR RIGID BODY DYNAMICS THAT EXACTLY PRESERVE ENERGY
//      AND MOMENTUM, J. C. SIMO AND K. K. WONG, INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN
//      ENGINEERING, VOL. 31, 19-52 (1991)
// Ref: "Rigid Body Dynamics of Mechanisms 1 Theoretical Basis" by Hubert Hahn, Section 4.2.4.1

template<typename Scalar>
class Inertial6DofForceFunction : public VectorValuedFunction<6,6,Scalar,48,0,double>
{
  public:
    double m;                           // mass
    Eigen::Matrix<double,3,1> c;        // offset
    Eigen::Matrix<double,3,3> J0;       // inertia tensor, with offset: J - m*C*C
    Eigen::Matrix<double,6,1> a_n, v_n; // velocity and acceleration vectors
                                        // note: according to current convention, the 3rd, 4th and 5th components of
                                        // a_n and v_n store the "convected" angular velocities and accelerations
    Eigen::Matrix<double,3,1> u_n;
    Eigen::Matrix<double,3,3> R_n;
    Eigen::Matrix<double,3,3> Rref;
    double beta, gamma, alphaf, alpham, dt; // time integration scheme parameters

  public:
    Inertial6DofForceFunction(const Eigen::Array<double,48,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      m = sconst[0];
      Eigen::Matrix<double,6,1> j;        // voigt form of inertia tensor
      j << sconst[1], sconst[2], sconst[3], sconst[4], sconst[5], sconst[6]; // [Jxx, Jyy, Jzz, Jyz, Jxz, Jxy]
      c << sconst[7], sconst[8], sconst[9];

      Eigen::Matrix<double,3,3> J;        // symmetric inertia tensor, no offset
      J <<  j[0],  j[5],  j[4],
            j[5],  j[1],  j[3],
            j[4],  j[3],  j[2];

      Eigen::Matrix<double,3,3> C;        // offset, skew-symmetric matrix
      C <<   0.0, -c[2],  c[1],
            c[2],   0.0, -c[0],
           -c[1],  c[0],   0.0;

      J0 = J - m*C*C;

      a_n = Eigen::Map<Eigen::Matrix<double,6,1> >(const_cast<double*>(sconst.data())+10);
      v_n = Eigen::Map<Eigen::Matrix<double,6,1> >(const_cast<double*>(sconst.data())+16);
      u_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+22);
      R_n = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+25);
      Rref = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+34);
      beta   = sconst[43];
      gamma  = sconst[44];
      alphaf = sconst[45];
      alpham = sconst[46];
      dt     = sconst[47];
    }

    Eigen::Matrix<Scalar,6,1> operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar t)
    {
      // inputs:
      // q[0] = x translation
      // q[1] = y translation
      // q[2] = z translation
      // q[3] = x component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration)
      // q[4] = y component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration)
      // q[5] = z component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration)

      Eigen::Matrix<Scalar,6,1> inc_displacement, v, a, f;
      Eigen::Matrix<Scalar,3,3> dR, incR, T;

      // define a few 3-dimensional sub-vectors, for convenience
      Eigen::VectorBlock<Eigen::Matrix<Scalar,6,1>,3> udot  = v.template head<3>(), // 1st time-derivative of displacement
                                                      Omega = v.template tail<3>(), // convected angular velocity
                                                      uddot = a.template head<3>(), // 2nd time-derivative of displacement
                                                      Alpha = a.template tail<3>(); // convected angular acceleration

      // set translation part of incremental generalized displacement
      inc_displacement.template head<3>() = q.template head<3>() - u_n.template cast<Scalar>();

      // set the rotation part of incremental generalized displacement
      vec_to_mat<Scalar>(q.template tail<3>(), dR); // dR is the spatial (left) incremental rotation matrix w.r.t Rref,
      Eigen::Matrix<Scalar,3,3> R = dR*Rref;        // R is the updated total rotation 
      incR = R_n.transpose()*R;                     // incR is the material (right) incremental rotation matrix w.r.t R_n
                                                    // i.e. R = R_n*incR --> incR = R_n.transpose()*R,
      mat_to_vec<Scalar>(incR, inc_displacement.template tail<3>());

      // compute the updated generalized velocity and acceleration
      v = gamma/(dt*beta)*inc_displacement 
          + ((1-(1-alphaf)*gamma/beta)*v_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*a_n).template cast<Scalar>();

      a = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement 
          + (-(1-alpham)/(dt*beta)*v_n + ((alpham-1)/(2*beta)+1)*a_n).template cast<Scalar>();

      // compute the tangential transformation matrix T(q)
      tangential_transf<Scalar>(q.template tail<3>(), T);

      // compute forces
      // ref: "Rigid Body Dynamics of Mechanisms 1 Theoretical Basis" by Hubert Hahn, Section 4.2.4.1
      f.template head<3>() = m*(uddot - R*(c.cross(Alpha) - Omega.cross(Omega.cross(c))));

      // compute moments
      // note: Even though T(0) = I, we still multiply by T.transpose() so that the Jacobian will be correctly evaluated
      //       when this function is automatically or numerically differentiated
      f.template tail<3>() = T.transpose()*Rref*(-m*(R.transpose()*uddot).cross(c) + J0*Alpha + Omega.cross(J0*Omega));

      return f;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,6,6>
Jacobian<double,Inertial6DofForceFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t);

} // namespace Simo

#endif
