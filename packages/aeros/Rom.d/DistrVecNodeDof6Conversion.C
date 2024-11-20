#include "DistrVecNodeDof6Conversion.h"

namespace Rom {

template <int DOFS_PER_NODE>
DistrVecNodeDofConversion<DOFS_PER_NODE>::~DistrVecNodeDofConversion() {
  typedef typename ConversionContainer::const_iterator ConversionIt;
  const ConversionIt it_end = subConversions_.end();
  for (ConversionIt it = subConversions_.begin(); it != it_end; ++it) {
    delete *it;
  }

  typedef typename RestrictedConversionContainer::const_iterator RestrictedConversionIt;
  const RestrictedConversionIt jt_end = subRestrictedConversions_.end();
  for (RestrictedConversionIt jt = subRestrictedConversions_.begin(); jt != jt_end; ++jt) {
    delete *jt;
  } 
}

template
DistrVecNodeDofConversion<6>::~DistrVecNodeDofConversion();

template
DistrVecNodeDofConversion<1>::~DistrVecNodeDofConversion();

} // end namespace Rom
