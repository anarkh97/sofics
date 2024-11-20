#ifndef _NEOHOOKEANSTRAINENERGYDENSITYFUNCTION_H_
#define _NEOHOOKEANSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

namespace Simo {

template<typename Scalar>
class NeoHookeanStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double mu; // 2nd Lamé constant (shear modulus)

  public:
    NeoHookeanStrainEnergyDensityFunction(double _mu) : mu(_mu) {}
    NeoHookeanStrainEnergyDensityFunction(const Eigen::Array<double,1,1> &sconst, const Eigen::Array<int,0,1>&) {
      mu = sconst[0];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const {
      // inputs: deformation gradient
      // output: strain energy density

      return mu/2*((F*F.transpose()).trace() - 3);
    }

};

} // namespace Simo

#endif
