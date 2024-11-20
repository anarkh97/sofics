#ifndef FWK_NUMERIC_H
#define FWK_NUMERIC_H

#include "Types.h"

namespace Fwk {

template<typename UnitType, typename RepType>
class Numeric {
public:
  explicit Numeric(RepType v = RepType()) : value_(v) {}
	
  RepType value() const
	{ return value_; }
	
  const Numeric<UnitType, RepType>& operator=(const Numeric<UnitType, RepType>& v)
	{ value_ = v.value_; return *this; }
	
	bool operator==(const Numeric<UnitType, RepType>& v) const
	{ return value_ == v.value_; }
	
	bool operator!=(const Numeric<UnitType, RepType>& v) const
	{ return value_ != v.value_; }
	
	bool operator<(const Numeric<UnitType, RepType>& v) const
	{ return Numeric<UnitType, RepType>::value_ < v.value_; }
	
	bool operator<=(const Numeric<UnitType, RepType>& v) const
	{ return Numeric<UnitType, RepType>::value_ <= v.value_; }
	
	bool operator>(const Numeric<UnitType, RepType>& v) const
	{ return Numeric<UnitType, RepType>::value_ > v.value_; }

	bool operator>=(const Numeric<UnitType, RepType>& v) const
	{ return Numeric<UnitType, RepType>::value_ >= v.value_; }
 
  const Numeric<UnitType, RepType> & operator+=(const Numeric<UnitType, RepType> & v)
  { value_ += v.value_; return *this; }

  Numeric<UnitType, RepType> operator+() const
  { return Numeric<UnitType, RepType>(+value_); }
  
  Numeric<UnitType, RepType> operator-() const
  { return Numeric<UnitType, RepType>(-value_); }
  
  const Numeric<UnitType, RepType> & operator++()
  { ++value_; return *this; }
  
  Numeric<UnitType, RepType> operator++(int)
  { Numeric<UnitType, RepType> tmp(*this); ++value_; return tmp; }
  
  const Numeric<UnitType, RepType> & operator-=(const Numeric<UnitType, RepType> & v)
  { value_ -= v.value_; return *this; }
  
  const Numeric<UnitType, RepType> & operator--()
  { --value_; return *this; }
  
  Numeric<UnitType, RepType> operator--(int)
  { Numeric<UnitType, RepType> tmp(*this); --value_; return tmp; }

  template <typename FactorType>
  const Numeric<UnitType, RepType> & operator*=(const FactorType & c)
  { value_ *= c; return *this; }

  template<typename FactorType>
  Numeric<UnitType, RepType> operator*(const FactorType & c) const
  { return Numeric<UnitType, RepType>(value_ * c); }
  
  template <typename FactorType>
  const Numeric<UnitType, RepType> & operator/=(const FactorType & c)
  { value_ /= c; return *this; }

  template <typename FactorType>
  const Numeric<UnitType, RepType> & operator%=(const FactorType & c)
  { value_ %= c; return *this; }

protected:
  RepType value_;
};

template<typename UnitType, typename RepType>
inline
Numeric<UnitType, RepType> operator+(const Numeric<UnitType, RepType> & v1, const Numeric<UnitType, RepType> & v2) {
  return Numeric<UnitType, RepType>(v1.value() + v2.value());
}

template<typename UnitType, typename RepType>
inline
Numeric<UnitType, RepType> operator-(const Numeric<UnitType, RepType> & v1, const Numeric<UnitType, RepType> & v2) {
  return Numeric<UnitType, RepType>(v1.value() - v2.value());
}

template<typename UnitType, typename RepType, typename FactorType>
inline
Numeric<UnitType, RepType> operator*(const FactorType & c, const Numeric<UnitType, RepType> & v) {
  return Numeric<UnitType, RepType>(v.value() * c);
}

template<typename UnitType, typename RepType, typename FactorType>
inline
Numeric<UnitType, RepType> operator/(const Numeric<UnitType, RepType> & v, const FactorType & c) {
  return Numeric<UnitType, RepType>(v.value() / c);
}

template<typename UnitType, typename RepType, typename FactorType>
inline
Numeric<UnitType, RepType> operator%(const Numeric<UnitType, RepType> & v, const FactorType & c) {
  return Numeric<UnitType, RepType>(v.value() % c);
}

template<typename UnitType, typename RepType>
inline
Fwk::OStream & operator<<(Fwk::OStream & out, const Numeric<UnitType, RepType> & v) {
  return out << v.value();
}

} // end namespace Fwk

#endif /* FWK_NUMERIC_H */
