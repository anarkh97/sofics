#ifndef _LEFTCOMPOSITIONRULE_H_
#define _LEFTCOMPOSITIONRULE_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>

namespace Simo {

template<typename Scalar>
class LeftCompositionRule : public VectorValuedFunction<3,3,Scalar,3,0,double>
{
   Eigen::Matrix<double,3,1> Psi1;
   Eigen::Matrix<double,3,3> R1;

  public:
    LeftCompositionRule(const Eigen::Array<double,3,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      Psi1 << sconst[0], sconst[1], sconst[2];
      vec_to_mat<double>(Psi1, R1);
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // Euler rotation vector

      // outputs:
      Eigen::Matrix<Scalar,3,3> R;
      vec_to_mat<Scalar>(q,R);

      Eigen::Matrix<Scalar,3,1> Psi;
      mat_to_vec<Scalar>(R*R1, Psi);

      if(Psi1.norm() > M_PI) {
        Psi = complement_rot_vec<Scalar>(unscale_rotvec<Scalar>(Psi, complement_rot_vec<double>(Psi1).template cast<Scalar>()));
      }
      else {
        Psi = unscale_rotvec<Scalar>(Psi, Psi1.template cast<Scalar>());
      }

      return Psi;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
