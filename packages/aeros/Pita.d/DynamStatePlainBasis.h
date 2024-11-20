#ifndef PITA_DYNAMSTATEPLAINBASIS_H
#define PITA_DYNAMSTATEPLAINBASIS_H

#include "DynamStateBasis.h"
#include <deque>

namespace Pita {

class DynamStatePlainBasis : public DynamStateBasis {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStatePlainBasis);

  virtual size_t stateCount() const { return state_.size(); }
  using DynamStateBasis::state;
  DynamState state(size_t index) const { return state_[index]; } // Unchecked access
  DynamState & internalState(size_t index); // Unchecked access and unsafe (returns a reference to the internals)

  void stateIs(size_t index, const DynamState & newState) { state_[index] = newState; } // Unchecked access

  virtual void lastStateIs(const DynamState & ds);
  virtual void lastStateBasisIs(const DynamStateBasis * dsb);
  virtual void lastStateBasisIs(const DynamStatePlainBasis * dsb);
  
  virtual void firstStateIs(const DynamState & ds);
  virtual void firstStateBasisIs(const DynamStateBasis * dsb);
  virtual void firstStateBasisIs(const DynamStatePlainBasis * dsb);

  void stateBasisDel();
  
  static DynamStatePlainBasis::Ptr New(size_t vectorSize) { return new DynamStatePlainBasis(vectorSize); }
  
protected:
  explicit DynamStatePlainBasis(size_t vectorSize) : DynamStateBasis(vectorSize) {}

  void prependState(const DynamState & ds);
  void prependStateBasis(const DynamStatePlainBasis * dsb);
  
  void appendState(const DynamState & ds);
  void appendStateBasis(const DynamStatePlainBasis * dsb);
  
private:
  std::deque<DynamState> state_;
};

inline
void
DynamStatePlainBasis::stateBasisDel() {
  this->state_.clear();
}

inline
void
DynamStatePlainBasis::prependState(const DynamState & ds) {
  this->state_.push_front(ds);
}

inline 
void
DynamStatePlainBasis::prependStateBasis(const DynamStatePlainBasis * dsb) {
  this->state_.insert(this->state_.begin(), dsb->state_.begin(), dsb->state_.end());
}

inline
void
DynamStatePlainBasis::appendState(const DynamState & ds) {
  this->state_.push_back(ds);
}

inline 
void
DynamStatePlainBasis::appendStateBasis(const DynamStatePlainBasis * dsb) {
  this->state_.insert(this->state_.end(), dsb->state_.begin(), dsb->state_.end());
}

} // end namespace Pita

#endif /* PITA_DYNAMSTATEPLAINBASIS_H */
