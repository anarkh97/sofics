#ifndef PITA_DYNAMSTATEBASISWRAPPER_H
#define PITA_DYNAMSTATEBASISWRAPPER_H

#include "Fwk.h"
#include "DynamStateBasis.h"

namespace Pita {

class DynamStateBasisWrapper : public DynamStateBasis {
public:
  typedef Fwk::Ptr<DynamStateBasisWrapper> Ptr;
  typedef Fwk::Ptr<const DynamStateBasisWrapper> PtrConst;

  virtual size_t stateCount() const { return stateCount_; }
  using DynamStateBasis::state;
  virtual DynamState state(size_t index) const { // Unsafe
    return DynamState(vectorSize(), data_ + 2 * index * vectorSize());
  }

  static DynamStateBasisWrapper::Ptr New(size_t vectorSize, size_t stateCount, double * data) { 
    return new DynamStateBasisWrapper(vectorSize, stateCount, data); 
  }

protected:
  DynamStateBasisWrapper(size_t vectorSize, size_t stateCount, double * data);
  
private:
  size_t stateCount_;
  double * data_;
};

}

#endif /* PITA_DYNAMSTATEBASISWRAPPER_H */
