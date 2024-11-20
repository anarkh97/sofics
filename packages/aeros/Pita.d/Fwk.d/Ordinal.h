#ifndef FWK_ORDINAL_H
#define FWK_ORDINAL_H

#include "Nominal.h"

namespace Fwk {

template<class UnitType, class RepType>
class Ordinal : public Nominal<UnitType, RepType> {
public:
  explicit Ordinal(RepType v = RepType()) : Nominal<UnitType, RepType>(v) {}
	
	bool operator<(const Ordinal<UnitType, RepType>& v) const
	{ return Nominal<UnitType, RepType>::value_ < v.value_; }
	
	bool operator<=(const Ordinal<UnitType, RepType>& v) const
	{ return Nominal<UnitType, RepType>::value_ <= v.value_; }
	
	bool operator>(const Ordinal<UnitType, RepType>& v) const
	{ return Nominal<UnitType, RepType>::value_ > v.value_; }

	bool operator>=(const Ordinal<UnitType, RepType>& v) const
	{ return Nominal<UnitType, RepType>::value_ >= v.value_; }
};

}

#endif /* FWK_ORDINAL_H */
