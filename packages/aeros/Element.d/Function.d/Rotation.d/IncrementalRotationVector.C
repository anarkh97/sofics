#ifdef USE_EIGEN3
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <iostream>

namespace Simo {

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,IncrementalRotationVector>
::operator() (const Eigen::Matrix<double,3,1>& q, double)
{
  if( (q.array() == 0).all() ) {
    Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > R1(sconst.data()+0), R2(sconst.data()+9);

    Eigen::Vector3d psi;
    mat_to_vec<double>(R1.transpose()*R2, psi);

    Eigen::Matrix3d Q;
    pseudorot_var(psi, Q); // note: Q = T(ψ)^{-T}

    return Q*R1.transpose();
  }
  else {
    std::cerr << "Error: Jacobian<double,IncrementalRotationVector>::operator() is not implemented for Δψ != 0\n";
    return Eigen::Matrix3d::Zero();
  }
}

} // namespace Simo

#endif

