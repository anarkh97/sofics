#ifndef PITA_LOCALCORRECTIONTIMESLICE_H
#define PITA_LOCALCORRECTIONTIMESLICE_H

#include "CorrectionPropagator.h"
#include "Seed.h"

#include "DynamPropagator.h"

namespace Pita {

class FullCorrectionPropagatorImpl : public CorrectionPropagator<DynamState> {
public:
  EXPORT_PTRINTERFACE_TYPES(FullCorrectionPropagatorImpl);
  class Manager;

  virtual void iterationIs(IterationRank i); // overriden

protected:
  friend class Manager;
 
  FullCorrectionPropagatorImpl(const String & name, DynamPropagator * prop);
  
private:
  DynamPropagator::Ptr propagator_;  
};

class FullCorrectionPropagatorImpl::Manager : public CorrectionPropagator<DynamState>::Manager, private Fwk::GenManagerImpl<FullCorrectionPropagatorImpl, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  // Instance management
  FullCorrectionPropagatorImpl * instance(const String & key) const { return Impl::instance(key); }
  InstanceCount instanceCount() const { return Impl::instanceCount(); }
  FullCorrectionPropagatorImpl * instanceNew(const String & key) { return Impl::instanceNew(key); }
  void instanceDel(const String & key) { Impl::instanceDel(key); }
  
  // Propagator is shared by all instances
  DynamPropagator * sharedPropagator() const { return sharedPropagator_.ptr(); }
  void sharedPropagatorIs(DynamPropagator * sp) { sharedPropagator_ = sp; }
  
  static Ptr New(DynamPropagator * sharedPropagator) {
    return new Manager(sharedPropagator);
  }
  
protected:
  explicit Manager(DynamPropagator * sp);
  
  virtual FullCorrectionPropagatorImpl * createNewInstance(const String & key); // overriden
  
private:
  typedef Fwk::GenManagerImpl<FullCorrectionPropagatorImpl, String> Impl;
  
  DynamPropagator::Ptr sharedPropagator_;

  DISALLOW_COPY_AND_ASSIGN(Manager);
};

} /* end namespace Pita */

#endif /* PITA_LOCALCORRECTIONTIMESLICE_H */
