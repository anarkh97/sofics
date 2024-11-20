#ifndef _COMPRESSIBLEMOONEYRIVLINSTRAINENERGYDENSITYFUNCTION_H_
#define _COMPRESSIBLEMOONEYRIVLINSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

// coupled form of the compressible Mooney-Rivlin model
// reference: Holzapfel, "Nonlinear Solid Mechanics" (2000), p 247, equation 6.147

namespace Simo {

template<typename Scalar, bool Uncoupled>
class CompressibleMooneyRivlinStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double c1, c2, c, d;

  public:
    CompressibleMooneyRivlinStrainEnergyDensityFunction(double _c1, double _c2, double _c) : c1(_c1), c2(_c2), c(_c)
    {
      d  = 2*(c1+2*c2);
    }
    CompressibleMooneyRivlinStrainEnergyDensityFunction(const Eigen::Array<double,3,1> &sconst, const Eigen::Array<int,0,1>&) 
    {
      c1 = sconst[0];
      c2 = sconst[1];
      c  = sconst[2];
      d  = 2*(c1+2*c2);
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const 
    {
      // inputs: deformation gradient
      // output: strain energy density

      using std::log;
      Eigen::Matrix<Scalar,3,3> C = F.transpose()*F;
      Scalar I1 = C.trace();
      Scalar I2 = 0.5*(I1*I1-(C*C).trace());
      Scalar J = F.determinant();
      if(Uncoupled) {
        using std::pow;
        return c1*(pow(J,-2/3.)*I1 - 3) + c2*(pow(J,-4/3.)*I2 - 3) + c*(J-1)*(J-1);
      }
      else {
        // note: this is the variant from which the constitutive relation in Material.d/MooneyRivlin.cpp is derived.
        return c1*(I1-3) + c2*(I2-3) + c*(J-1)*(J-1) - d*log(J);
      }
    }
};

} // namespace Simo

#endif
