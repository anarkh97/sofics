#ifndef FWK_NOMINAL_H
#define FWK_NOMINAL_H

namespace Fwk {

template<class UnitType, class RepType>
class Nominal {
public:
  explicit Nominal(RepType v = RepType()) : value_(v) {}
	
	bool operator==(const Nominal<UnitType, RepType>& v) const
	{ return value_ == v.value_; }
	
	bool operator!=(const Nominal<UnitType, RepType>& v) const
	{ return value_ != v.value_; }
	
	Nominal<UnitType, RepType>& operator=(const Nominal<UnitType, RepType>& v)
	{ value_ = v.value_; return *this; }
	
	RepType value() const
	{ return value_; }
	
protected:
	RepType value_;
};

template <class UnitType, class RepType>
inline std::ostream & operator<< (std::ostream & os, const Nominal<UnitType, RepType> & n) {
    return os << n.value();
}


} // end namespace Fwk

#endif
