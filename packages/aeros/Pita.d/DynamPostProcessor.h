#ifndef PITA_DYNAMPOSTPROCESSOR_H
#define PITA_DYNAMPOSTPROCESSOR_H

#include "Fwk.h"
#include "DynamTimeIntegrator.h"

namespace Pita {

class DynamPostProcessor : public Fwk::PtrInterface<DynamPostProcessor> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamPostProcessor);

  enum Status {
    CLOSED = 0,
    OPEN = 1
  };

  Status status() const { return status_; }
  SliceRank sliceRank() const { return sliceRank_; } // HACK
  
  virtual void statusIs(Status s) = 0;
  virtual void sliceRankIs(SliceRank r) = 0; // HACK

  class IntegratorReactor : public DynamTimeIntegrator::NotifieeConst {
  public:
    typedef Fwk::Ptr<IntegratorReactor> Ptr;
    typedef Fwk::Ptr<const IntegratorReactor> PtrConst;
   
    DynamPostProcessor::Ptr parent() const { return getParent(); } 
    
    virtual void onCurrentSlice() { getParent()->sliceRankIs(notifier()->currentSlice()); } // HACK

  protected:
    explicit IntegratorReactor(const DynamTimeIntegrator * notifier) :
      DynamTimeIntegrator::NotifieeConst(notifier) {}

    virtual DynamPostProcessor * getParent() const = 0;
  };

  IntegratorReactor::Ptr integratorReactorNew(const DynamTimeIntegrator * notifier) {
    return getNewIntegratorReactor(notifier);
  }
 
protected:
  explicit DynamPostProcessor(SliceRank slice = SliceRank()) : 
    status_(CLOSED), 
    sliceRank_(slice) // HACK
  {}
  
  void setStatus(Status s) { status_ = s; }
  void setSliceRank(SliceRank r) { sliceRank_ = r; } // HACK
  
  virtual IntegratorReactor * getNewIntegratorReactor(const DynamTimeIntegrator * notifier) = 0;

private:
  Status status_;
  SliceRank sliceRank_; // HACK
};
  
} // end namespace Pita

#endif /* PITA_DYNAMPOSTPROCESSOR_H */
