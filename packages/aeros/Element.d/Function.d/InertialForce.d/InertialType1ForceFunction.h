#ifndef _INERTIALTYPE1FORCEFUNCTION_H_
#define _INERTIALTYPE1FORCEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>

// Compute the inertial moment for a rigid body with 3 rotational degrees of freedom,
// according to the classical Euler's equation for the rigid body and the extended Newmark 
// method (based on ALGO_1 from reference below)
// Ref: UNCONDITIONALLY STABLE ALGORITHMS FOR RIGID BODY DYNAMICS THAT EXACTLY PRESERVE ENERGY
//      AND MOMENTUM, J. C. SIMO AND K. K. WONG, INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN
//      ENGINEERING, VOL. 31, 19-52 (1991)

namespace Simo {

template<typename Scalar>
class InertialType1ForceFunction : public VectorValuedFunction<3,3,Scalar,38,0,double>
{
  public:
    Eigen::Matrix<double,3,3> J;            // inertia tensor
    Eigen::Matrix<double,3,1> A_n, V_n;     // angular velocity and acceleration vectors
                                            // note: according to current convention, A_n and V_n store the
                                            // "convected" angular velocities and accelerations
    Eigen::Matrix<double,3,3> R_n;
    Eigen::Matrix<double,3,3> Rref;
    double beta, gamma, alphaf, alpham, dt; // time integration scheme parameters

  public:
    InertialType1ForceFunction(const Eigen::Array<double,38,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      J = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
      A_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+9);
      V_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+12);
      R_n = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+15);
      Rref = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+24);
      beta   = sconst[33];
      gamma  = sconst[34];
      alphaf = sconst[35];
      alpham = sconst[36];
      dt     = sconst[37];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // inputs:
      // q[0] = x component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration)
      // q[1] = y component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration)
      // q[2] = z component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration)

      // compute incRref, the spatial (left) incremental rotation matrix w.r.t Rref
      Eigen::Matrix<Scalar,3,3> incRref;
      vec_to_mat<Scalar>(q, incRref); 

      // compute incR_n, the material (right) incremental rotation matrix w.r.t R_n, and its corresponding rotation vector
      Eigen::Matrix<Scalar,3,3> incR_n = R_n.transpose()*incRref*Rref;
      Eigen::Matrix<Scalar,3,1> inc_displacement;
      mat_to_vec<Scalar>(incR_n, inc_displacement);

      // compute the updated convected angular velocity and acceleration
      Eigen::Matrix<Scalar,3,1> V, A;
      V = gamma/(dt*beta)*inc_displacement
          + ((1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n).template cast<Scalar>();
      A = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement
          + (-(1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n).template cast<Scalar>();

      // compute the tangential transformation matrix T(q)
      Eigen::Matrix<Scalar,3,3> T;
      tangential_transf(q, T);

      // convected description of the angular momentum balance equation (see Simo & Wong eq. 29)
      // premultiplied by transformation to fixed reference frame
      // note: Even though T(0) = I, we still multiply by T.transpose() so that the Jacobian will be correctly evaluated
      //       when this function is automatically or numerically differentiated
      return T.transpose()*Rref*(J*A + V.cross(J*V));
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,InertialType1ForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t);

} // namespace Simo

#endif
