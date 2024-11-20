#ifndef _STRAINENERGYDENSITYFUNCTION_H_
#define _STRAINENERGYDENSITYFUNCTION_H_

#include <Eigen/Core>

namespace Simo {

template <typename Scalar>
class StrainEnergyDensityFunction
{
  public:
    virtual Scalar operator() (const Eigen::Matrix<Scalar,3,3>& F) const = 0;
};

} // namespace Simo

#endif
