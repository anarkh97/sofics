#ifndef _PSEUDOTANGENTIALMOMENTWORKFUNCTION_H_
#define _PSEUDOTANGENTIALMOMENTWORKFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>
#include <Eigen/Geometry>

namespace Simo {

template<typename Scalar>
class PseudoTangentialMomentWorkFunction : public ScalarValuedFunction<3,Scalar,16,0,double>
{
    double M0;
    Eigen::Matrix<double,3,1> l,k;
    Eigen::Matrix<double,3,1> m0; // applied moment in the initial configuration
    Eigen::Matrix<double,3,3> Rref; // nodal rotation in the reference configuration

  public:
    PseudoTangentialMomentWorkFunction(const Eigen::Array<double,16,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      M0 = sconst[0];
      l << sconst[1], sconst[2], sconst[3];
      k << sconst[4], sconst[5], sconst[6];
      Rref << sconst[7], sconst[8], sconst[9],
              sconst[10], sconst[11], sconst[12],
              sconst[13], sconst[14], sconst[15];

      m0 = M0*l.cross(k); // ref: Ziegler (1968)
                          //      Ritto-Corea and Camotin, 2003, sec. 4.9
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // inputs:
      // q[0] = x component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[1] = y component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[2] = z component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1

      Eigen::Quaternion<Scalar> z;
      z.setFromOneVector(q);
      Eigen::Matrix<Scalar,3,3> R = z.toRotationMatrix()*Rref.template cast<Scalar>();
      Eigen::Matrix<Scalar,3,1> Psi; mat_to_vec<Scalar>(R,Psi);
      Eigen::Matrix<Scalar,3,1> p; // a rotation measure
      Scalar psi2 = Psi.squaredNorm();
      if(psi2 < 5e-6) {
         p = (1 - psi2/6 + psi2*psi2/120)*Psi
            +(0.5 - psi2/24 + psi2*psi2/720)*Psi.dot(k.template cast<Scalar>())*k.template cast<Scalar>().cross(Psi);
      }
      else {
        using std::sqrt;
        using std::sin;
        using std::cos;
        Scalar psi = sqrt(psi2);
        p = (sin(psi)/psi)*Psi + (1-cos(psi))/psi2*Psi.dot(k.template cast<Scalar>())*k.template cast<Scalar>().cross(Psi);
      }

      return -p.dot(m0.template cast<Scalar>()); 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
