#ifndef _INCOMPRESSIBLEHEXAQ1P0_H_
#define _INCOMPRESSIBLEHEXAQ1P0_H_

#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Function.d/StrainEnergy.d/FourFieldStrainEnergyFunction.h>
#include <Element.d/Function.d/Shape.d/Hex8LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Constant.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Force.d/MixedFiniteElement.h>

template <typename S>
using HexaQ1P0FourFieldStrainEnergyFunction 
      = Simo::FourFieldStrainEnergyFunction<S, Hex8LagrangePolynomialShapeFunction,
                                            ConstantShapeFunction, ConstantShapeFunction, ConstantShapeFunction,
                                            RepeatedQuadratureRule<double,GaussLegendre,3,Eigen::Vector3d> >;

class IncompressibleHexaQ1P0 : public MixedFiniteElement<HexaQ1P0FourFieldStrainEnergyFunction>
{
  public:
    static const DofSet NODALDOFS[8];
    IncompressibleHexaQ1P0(int* _nn);

    int getTopNumber() const override;
    PrioInfo examine(int sub, MultiFront *mf) override;
    int getQuadratureOrder();
};
#endif
#endif
