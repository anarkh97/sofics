#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Force.d/HexaQ2P0.h>

const DofSet HexaQ2P0::NODALDOFS[20] = { DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp };

PrioInfo examineHex20(int sub, MultiFront *mf, int *nn);

HexaQ2P0::HexaQ2P0(int* _nn)
 : MixedFiniteElement<HexaQ2P0ThreeFieldStrainEnergyFunction>(20, const_cast<DofSet*>(NODALDOFS), _nn)
{}

int
HexaQ2P0::getTopNumber() const
{
  return 172;
}

PrioInfo
HexaQ2P0::examine(int sub, MultiFront *mf)
{
  return examineHex20(sub, mf, nn);
}

int
HexaQ2P0::getQuadratureOrder()
{
  return 3;
}

#endif
