#ifndef _MOONEYRIVLINSTRAINENERGYDENSITYFUNCTION_H_
#define _MOONEYRIVLINSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

namespace Simo {

template<typename Scalar>
class MooneyRivlinStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double c1, c2;

  public:
    MooneyRivlinStrainEnergyDensityFunction(double _c1, double _c2) : c1(_c1), c2(_c2) {}
    MooneyRivlinStrainEnergyDensityFunction(const Eigen::Array<double,2,1> &sconst, const Eigen::Array<int,0,1>&) 
    {
      c1 = sconst[0];
      c2 = sconst[1];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const 
    {
      // inputs: deformation gradient
      // output: strain energy density

      Eigen::Matrix<Scalar,3,3> C = F.transpose()*F;
      Scalar I1 = C.trace();
      Scalar I2 = 0.5*(I1*I1-(C*C).trace());
      return c1*(I1-3) + c2*(I2-3);
    }

};

} // namespace Simo

#endif
