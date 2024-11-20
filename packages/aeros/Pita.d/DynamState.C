#include "DynamState.h"

#include <functional>
#include <algorithm>

namespace Pita {

DynamState::Desc::Desc(size_t vectorSize, Scalar initValue) :
  disp_(vectorSize, initValue),
  vel_(vectorSize, initValue)
{}

DynamState::Desc::Desc(size_t vectorSize) :
  disp_(vectorSize),
  vel_(vectorSize)
{}

DynamState::Desc::Desc(size_t vectorSize, const DynamState::Scalar * data) :
  disp_(vectorSize, const_cast<Scalar *>(data), true),
  vel_(vectorSize, const_cast<Scalar *>(data + vectorSize), true)
{}

DynamState::Desc::Desc(const GenVector<double> & disp, const GenVector<double> & vel) :
  disp_(disp),
  vel_(vel)
{}

DynamState::Desc &
DynamState::Desc::operator+=(const DynamState::Desc & dsd) {
  disp_ += dsd.disp_;
  vel_ += dsd.vel_;
  return *this;
}

DynamState::Desc &
DynamState::Desc::operator-=(const DynamState::Desc & dsd) {
  disp_ -= dsd.disp_;
  vel_ -= dsd.vel_;
  return *this;
}

DynamState::Desc &
DynamState::Desc::operator*=(double coef) {
  std::transform(disp_.data(), disp_.data() + disp_.size(), disp_.data(), std::bind2nd(std::multiplies<double>(), coef));
  std::transform(vel_.data(), vel_.data() + vel_.size(), vel_.data(), std::bind2nd(std::multiplies<double>(), coef));
  return *this;
}

DynamState::Desc &
DynamState::Desc::operator/=(double coef) {
  std::transform(disp_.data(), disp_.data() + disp_.size(), disp_.data(), std::bind2nd(std::divides<double>(), coef));
  std::transform(vel_.data(), vel_.data() + vel_.size(), vel_.data(), std::bind2nd(std::divides<double>(), coef));
  return *this;
}

void
DynamState::Desc::linAdd(double coef, const DynamState::Desc & dsd) {
  disp_.linAdd(coef, dsd.disp_);
  vel_.linAdd(coef, dsd.vel_);
}

//--------------------
  
DynamState::DynamState(size_t vectorSize) :
  desc_(new DynamState::Desc(vectorSize))
{}

DynamState::DynamState(size_t vectorSize, Scalar initialValue) :
  desc_(new DynamState::Desc(vectorSize, initialValue))
{}

DynamState::DynamState(size_t vectorSize, const DynamState::Scalar * data) :
  desc_(new DynamState::Desc(vectorSize, data))
{}

DynamState::DynamState(const GenVector<double> & disp, const GenVector<double> & vel) :
  desc_(new DynamState::Desc(disp, vel))
{}

size_t
DynamState::vectorSize() const {
  return displacement().size(); 
}

const DynamState::VectorType &
DynamState::displacement() const {
  return desc_->displacement();
}

const DynamState::VectorType &
DynamState::velocity() const {
  return desc_->velocity();
}

void
DynamState::unshareDesc() {
  if (desc_->references() > 1ul)
    desc_ = new DynamState::Desc(desc_->displacement(), desc_->velocity());
}

DynamState::VectorType &
DynamState::displacement() {
  unshareDesc();
  return desc_->displacement();
}

DynamState::VectorType &
DynamState::velocity() {
  unshareDesc();
  return desc_->velocity();
}

DynamState &
DynamState::operator+=(const DynamState & ds) {
  unshareDesc();
  *desc_ += *(ds.desc_);
  return *this;
}

DynamState &
DynamState::operator-=(const DynamState & ds) {
  unshareDesc();
  *desc_ -= *(ds.desc_);
  return *this;
}

DynamState &
DynamState::operator*=(double coef) {
  unshareDesc();
  *desc_ *= coef;
  return *this;
}

DynamState &
DynamState::operator/=(double coef) {
  unshareDesc();
  *desc_ /= coef;
  return *this;
}

void
DynamState::linAdd(double coef, const DynamState & ds) {
  unshareDesc();
  desc_->linAdd(coef, *ds.desc_);
}

const DynamState
operator+(const DynamState & op1, const DynamState & op2) {
  DynamState temp(op1);
  temp += op2;
  return temp;
}

const DynamState
operator-(const DynamState & op1, const DynamState & op2) {
  DynamState temp(op1);
  temp -= op2;
  return temp;
}

double
operator*(const DynamState & op1, const DynamState & op2) {
  return (op1.displacement() * op2.displacement()) + (op1.velocity() * op2.velocity());
}

} // end namespace Pita
