#ifndef PITA_DYNAMSTATEBASIS_H
#define PITA_DYNAMSTATEBASIS_H

#include "Fwk.h"
#include "DynamState.h"

namespace Pita {

class DynamStateBasis : public Fwk::PtrInterface<DynamStateBasis> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateBasis);
  class IteratorConst;

  size_t vectorSize() const { return vectorSize_; }
  virtual size_t stateCount() const = 0;
  
  virtual DynamState state(size_t index) const = 0; // Unsafe
  IteratorConst state() const;
  
protected:
  explicit DynamStateBasis(size_t vectorSize) : vectorSize_(vectorSize) {}
  
  DynamStateBasis(const DynamStateBasis &); // No implementation
  DynamStateBasis & operator=(const DynamStateBasis &); // No implementation
  
private:
  size_t vectorSize_;
};


class DynamStateBasis::IteratorConst {
public:
  const DynamState & operator*() const { return current_; }
  const DynamState * operator->() const { return &current_; }
  
  operator bool() const { return current_.vectorSize(); }

  IteratorConst & operator++() { ++currentRank_; updateCurrent(); return *this; }
  IteratorConst operator++(int) {
    IteratorConst temp(*this);
    ++(*this);
    return temp;
  }
  
  // IteratorConst(const IteratorConst &); // default
  // IteratorConst & operator=(const IteratorConst &); // default

protected:
  explicit IteratorConst(const DynamStateBasis * parent) :
    parent_(parent),
    currentRank_(0),
    current_()
  {
    updateCurrent();
  }

  friend class DynamStateBasis;

private:
  void updateCurrent() {
    current_ = (currentRank_ < parent_->stateCount()) ? parent_->state(currentRank_) : DynamState();
  }

  DynamStateBasis::PtrConst parent_;
  size_t currentRank_;
  DynamState current_;
};

inline
DynamStateBasis::IteratorConst
DynamStateBasis::state() const {
  return IteratorConst(this);
}

} /* end namespace Pita */

#endif /* PITA_DYNAMSTATEBASIS_H */
