#ifndef _STVENANTKIRCHHOFFSTRAINENERGYDENSITYFUNCTION_H_
#define _STVENANTKIRCHHOFFSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

// reference: Holzapfel, "Nonlinear Solid Mechanics" (2000), p 250 equation 6.151
// The Saint-Venant Kirchhoff model is a classical nonlinear model for compressible
// hyperelastic materials often used for metals.

namespace Simo {

template<typename Scalar>
class SaintVenantKirchhoffStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double lambda, mu; // Lamé constants

  public:
    SaintVenantKirchhoffStrainEnergyDensityFunction(double _lambda, double _mu) : lambda(_lambda), mu(_mu) {}
    SaintVenantKirchhoffStrainEnergyDensityFunction(const Eigen::Array<double,2,1> &sconst, const Eigen::Array<int,0,1>&) {
      lambda = sconst[0];
      mu     = sconst[1];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const {
      // inputs: deformation gradient
      // output: strain energy density

      //reference implementation:
      //Eigen::Matrix<Scalar,3,3> E = 0.5*(F.transpose()*F - Eigen::Matrix<Scalar,3,3>::Identity());
      //return lambda/2*E.trace()*E.trace() + mu*(E*E).trace();

      //optimized implementation:
      using Eigen::numext::abs2;
      Scalar E00 = 0.5*(F.col(0).squaredNorm()-1),
             E01 = 0.5*F.col(1).dot(F.col(0)),
             E02 = 0.5*F.col(2).dot(F.col(0)),
             E11 = 0.5*(F.col(1).squaredNorm()-1),
             E12 = 0.5*F.col(2).dot(F.col(1)),
             E22 = 0.5*(F.col(2).squaredNorm()-1);
      Scalar Etrace = E00 + E11 + E22;
      Scalar EEtrace = abs2(E00) + abs2(E11) + abs2(E22) + 2*(abs2(E01) +abs2(E02) + abs2(E12));
      return lambda/2*abs2(Etrace) + mu*EEtrace;
    }

};

} // namespace Simo

#endif
