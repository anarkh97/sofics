#ifndef _AUTODIFF_SCALAR_PLUGIN_H_
#define _AUTODIFF_SCALAR_PLUGIN_H_

#include <unsupported/Eigen/AutoDiff>

namespace Eigen {

template<typename DerType>
inline const AutoDiffScalar<Matrix<typename internal::traits<DerType>::Scalar,Dynamic,1> >
atan(const AutoDiffScalar<DerType>& a)
{
  using std::atan;
  typedef typename internal::traits<DerType>::Scalar Scalar;
  typedef AutoDiffScalar<Matrix<Scalar,Dynamic,1> > PlainADS;
  PlainADS ret;
  ret.value() = atan(a.value());
  ret.derivatives() = a.derivatives()/(Scalar(1)+a.value()*a.value());
  return ret;
}

}

#endif
