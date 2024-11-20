#ifndef _HEXAQ1P0_H_
#define _HEXAQ1P0_H_

#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Function.d/StrainEnergy.d/ThreeFieldStrainEnergyFunction.h>
#include <Element.d/Function.d/Shape.d/Hex8LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Constant.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Force.d/MixedFiniteElement.h>

template <typename S>
using HexaQ1P0ThreeFieldStrainEnergyFunction 
      = Simo::ThreeFieldStrainEnergyFunction<S, Hex8LagrangePolynomialShapeFunction,
                                             ConstantShapeFunction, ConstantShapeFunction,
                                             RepeatedQuadratureRule<double,GaussLegendre,3,Eigen::Vector3d> >;

class HexaQ1P0 : public MixedFiniteElement<HexaQ1P0ThreeFieldStrainEnergyFunction>
{
  public:
    static const DofSet NODALDOFS[8];
    HexaQ1P0(int* _nn);

    int getTopNumber() const override;
    PrioInfo examine(int sub, MultiFront *mf) override;
    int getQuadratureOrder(); 
};
#endif
#endif
