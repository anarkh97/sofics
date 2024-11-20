#ifndef PITA_DYNAMSTATE_H
#define PITA_DYNAMSTATE_H

#include "Fwk.h"
#include <Math.d/Vector.h>

namespace Pita {

class DynamState {
public:
  typedef double Scalar;
  typedef GenVector<Scalar> VectorType;
    
  explicit DynamState(size_t vectorSize = 0);
  DynamState(size_t vectorSize, Scalar initialValue);
  DynamState(size_t vectorSize, const Scalar * data);
  DynamState(const GenVector<double> & disp, const GenVector<double> & vel);
  
  DynamState(const DynamState & other) :
    desc_(other.desc_)
  {}

  DynamState & operator=(const DynamState & other) {
    desc_ = other.desc_;
    return *this;
  }

  size_t vectorSize() const;

  const GenVector<Scalar> & displacement() const;
  const GenVector<Scalar> & velocity() const;
  
  GenVector<Scalar> & displacement();
  GenVector<Scalar> & velocity();

  DynamState & operator+=(const DynamState & ds); 
  DynamState & operator-=(const DynamState & ds);
  DynamState & operator*=(double coef);
  DynamState & operator/=(double coef);

  void linAdd(double coef, const DynamState & ds);

protected:
  class Desc : public PtrInterface<Desc> {
  public:
    EXPORT_PTRINTERFACE_TYPES(Desc);

    typedef double Scalar;
    typedef GenVector<Scalar> VectorType;

    explicit Desc(size_t vectorSize);
    Desc(size_t vectorSize, Scalar initvalue);
    Desc(size_t vectorSize, const Scalar * data);
    Desc(const GenVector<double> & disp, const GenVector<double> & vel);

    const VectorType & displacement() const { return disp_; }
    const VectorType & velocity() const { return vel_; }

    VectorType & displacement() { return disp_; }
    VectorType & velocity() { return vel_; }

    Desc & operator+=(const Desc & dsd);
    Desc & operator-=(const Desc & dsd);
    Desc & operator*=(double coef);
    Desc & operator/=(double coef);

    void linAdd(double coef, const Desc & dsd);

  private:
    VectorType disp_;
    VectorType vel_;

    DISALLOW_COPY_AND_ASSIGN(Desc);
};
  
  void unshareDesc();
  friend void unshare(DynamState & state);

private:
  DynamState::Desc::Ptr desc_;
};

inline
void
unshare(DynamState & state) {
  state.unshareDesc();
}

const DynamState operator+(const DynamState & op1, const DynamState & op2);
const DynamState operator-(const DynamState & op1, const DynamState & op2);
double operator*(const DynamState & op1, const DynamState & op2);

template <typename StateType>
void
bufferStateCopy(const StateType & s, typename StateType::Scalar * buffer) {
  typename StateType::Scalar * data = const_cast<typename StateType::VectorType &>(s.displacement()).getData();
  std::copy(data, data + s.vectorSize(), buffer);
  data = const_cast<typename StateType::VectorType &>(s.velocity()).getData();
  std::copy(data, data + s.vectorSize(), buffer + s.vectorSize());
}


} // end namespace Pita

#endif
