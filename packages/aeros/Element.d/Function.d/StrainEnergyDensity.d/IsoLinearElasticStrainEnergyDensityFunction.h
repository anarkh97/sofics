#ifndef _ISOLINEARELASTICSTRAINENERGYDENSITYFUNCTION_H_
#define _ISOLINEARELASTICSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

namespace Simo {

template<typename Scalar>
class IsoLinearElasticStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double lambda, mu; // Lamé constants

  public:
    IsoLinearElasticStrainEnergyDensityFunction(double _lambda, double _mu) : lambda(_lambda), mu(_mu) {}
    IsoLinearElasticStrainEnergyDensityFunction(const Eigen::Array<double,2,1> &sconst, const Eigen::Array<int,0,1>&) {
      lambda = sconst[0];
      mu     = sconst[1];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const {
      // inputs: deformation gradient
      // output: strain energy density

      //reference implementation:
      //Eigen::Matrix<Scalar,3,3> epsilon = 0.5*(F.transpose() + F) - Eigen::Matrix<Scalar,3,3>::Identity();
      //return lambda/2*epsilon.trace()*epsilon.trace() + mu*(epsilon*epsilon).trace();

      //optimized implementation:
      Scalar E00 = F(0,0)-1,
             E01 = 0.5*(F(0,1)+F(1,0)),
             E02 = 0.5*(F(0,2)+F(2,0)),
             E11 = F(1,1)-1,
             E12 = 0.5*(F(1,2)+F(2,1)),
             E22 = F(2,2)-1;
      Scalar Etrace = E00 + E11 + E22;
      using Eigen::numext::abs2;
      Scalar EEtrace = abs2(E00) + abs2(E11) + abs2(E22) + 2*(abs2(E01) + abs2(E02) + abs2(E12));
      return lambda/2*abs2(Etrace) + mu*EEtrace;
    }

};

} // namespace Simo

#endif
