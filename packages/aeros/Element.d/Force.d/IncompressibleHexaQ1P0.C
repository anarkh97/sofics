#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Force.d/IncompressibleHexaQ1P0.h>

const DofSet IncompressibleHexaQ1P0::NODALDOFS[8] = { DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                                      DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp };

PrioInfo examineHex8(int sub, MultiFront *mf, int *nn);

IncompressibleHexaQ1P0::IncompressibleHexaQ1P0(int* _nn)
 : MixedFiniteElement<HexaQ1P0FourFieldStrainEnergyFunction>(8, const_cast<DofSet*>(NODALDOFS), _nn)
{}

int
IncompressibleHexaQ1P0::getTopNumber() const
{
  return 117;
}

PrioInfo
IncompressibleHexaQ1P0::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

int 
IncompressibleHexaQ1P0::getQuadratureOrder()
{
  return 2;
}

#endif
