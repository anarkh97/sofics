#ifndef FWK_INTERVAL_H
#define FWK_INTERVAL_H

#include "Ordinal.h"

namespace Fwk {

template<typename UnitType, typename RepType, typename DiffType>
class Interval : public Ordinal<UnitType, RepType> {
public:
  explicit Interval(RepType v = RepType()) : Ordinal<UnitType, RepType>(v) {}

  DiffType operator-(const Interval<UnitType, RepType, DiffType> & v) const { return DiffType(this->value() - v.value()); }
};

}

#endif /* FWK_INTERVAL_H */
