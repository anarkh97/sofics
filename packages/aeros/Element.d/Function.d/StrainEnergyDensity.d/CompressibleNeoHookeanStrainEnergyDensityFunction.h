#ifndef _COMPRESSIBLENEOHOOKEANSTRAINENERGYDENSITYFUNCTION_H_
#define _COMPRESSIBLENEOHOOKEANSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

// references
// [1] : Bonet & Wood, "Nonlinear continuum mechanics for finite element analysis", 2nd edition (2008), p 162 equation 6.27
// [2] (uncoupled variant) : J.C. Simo, R.L. Taylor and K.S. Pister, "Variational and projection methods for the volume constraint
//     in finite deformation plasticity", Computer Methods in Applied Mechanics and Engineering, 51 (1985) 177-208, Equation 3.12a

namespace Simo {

template<typename Scalar, bool Uncoupled>
class CompressibleNeoHookeanStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double lambda, mu; // Lamé constants

  public:
    CompressibleNeoHookeanStrainEnergyDensityFunction(double _lambda, double _mu) : lambda(_lambda), mu(_mu) {}
    CompressibleNeoHookeanStrainEnergyDensityFunction(const Eigen::Array<double,2,1> &sconst, const Eigen::Array<int,0,1>&) {
      lambda = sconst[0];
      mu     = sconst[1];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const {
      // inputs: deformation gradient
      // output: strain energy density

      using std::log;
      Scalar J = F.determinant();
      Scalar lnJ = log(J);
      if(Uncoupled) {
        using std::pow;
        double K = lambda + 2*mu/3; // bulk modulus
        return mu/2*(pow(J,-2/3.)*(F*F.transpose()).trace() - 3) + K/2*(J-1)*(J-1);
      }
      else {
        return mu/2*((F*F.transpose()).trace() - 3) - mu*lnJ + lambda/2*lnJ*lnJ;
      }
    }

};

} // namespace Simo

#endif
