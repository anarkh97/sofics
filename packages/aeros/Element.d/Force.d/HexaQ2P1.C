#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Force.d/HexaQ2P1.h>

const DofSet HexaQ2P1::NODALDOFS[20] = { DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp };

PrioInfo examineHex20(int sub, MultiFront *mf, int *nn);

HexaQ2P1::HexaQ2P1(int* _nn)
 : MixedFiniteElement<HexaQ2P1ThreeFieldStrainEnergyFunction>(20, const_cast<DofSet*>(NODALDOFS), _nn)
{}

int
HexaQ2P1::getTopNumber() const
{
  return 172;
}

PrioInfo
HexaQ2P1::examine(int sub, MultiFront *mf)
{
  return examineHex20(sub, mf, nn);
}

int
HexaQ2P1::getQuadratureOrder()
{
  return 3;
}

#endif
