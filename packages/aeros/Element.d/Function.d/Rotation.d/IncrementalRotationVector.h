#ifndef _INCREMENTALROTATIONVECTOR_H_
#define _INCREMENTALROTATIONVECTOR_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/utilities.hpp>

namespace Simo {

template<typename Scalar>
class IncrementalRotationVector : public VectorValuedFunction<3,3,Scalar,18,0,double>
{
   Eigen::Matrix3d R1, R2;

  public:
    IncrementalRotationVector(const Eigen::Array<double,18,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      R1 = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
      R2 = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+9);
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar)
    {
      // inputs:
      // Euler rotation vector

      // outputs:
      Eigen::Matrix<Scalar,3,3> R;
      vec_to_mat<Scalar>(q,R);

      Eigen::Matrix<Scalar,3,1> v;
      mat_to_vec<Scalar>(R1.transpose()*R*R2, v);
      return v;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,IncrementalRotationVector>
::operator() (const Eigen::Matrix<double,3,1>& q, double);

} // namespace Simo

#endif
