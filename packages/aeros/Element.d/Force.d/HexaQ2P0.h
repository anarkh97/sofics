#ifndef _HEXAQ2P0_H_
#define _HEXAQ2P0_H_

#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Function.d/StrainEnergy.d/ThreeFieldStrainEnergyFunction.h>
#include <Element.d/Function.d/Shape.d/Constant.h>
#include <Element.d/Function.d/Shape.d/Hex20LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Force.d/MixedFiniteElement.h>

template <typename S>
using HexaQ2P0ThreeFieldStrainEnergyFunction 
      = Simo::ThreeFieldStrainEnergyFunction<S, Hex20LagrangePolynomialShapeFunction,
                                             ConstantShapeFunction, ConstantShapeFunction,
                                             RepeatedQuadratureRule<double,GaussLegendre,3,Eigen::Vector3d> >;

class HexaQ2P0 : public MixedFiniteElement<HexaQ2P0ThreeFieldStrainEnergyFunction>
{
  public:
    static const DofSet NODALDOFS[20];
    HexaQ2P0(int* _nn);

    int getTopNumber() const override;
    PrioInfo examine(int sub, MultiFront *mf) override;
    int getQuadratureOrder(); 
};
#endif
#endif
