#ifndef PITA_CONCURRENTBASIS_MANAGER_H
#define PITA_CONCURRENTBASIS_MANAGER_H

#include "Fwk.h"

#include "DynamStatePlainBasis.h"

#include "IntegratorPropagator.h"
#include "NlDynamTimeIntegrator.h"

#include "LinearizedPropagator.h"

#include <map>

namespace Pita {

namespace ConcurrentBasis {

typedef GenIntegratorPropagator<NlDynamTimeIntegrator> Propagator;

class PropagatorReactor;

class Manager : public Fwk::PtrInterface<Manager> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  const Propagator * referencePropagator(const DynamStatePlainBasis * concurrentBasis) const;
  void referencePropagatorIs(DynamStatePlainBasis * concurrentBasis, const Propagator * cp);

  static Ptr New() { return new Manager(); }

protected:
  Manager();

protected:
  friend class PropagatorReactor;

private:
  typedef std::map<const DynamStatePlainBasis *, Fwk::Ptr<PropagatorReactor> > ReactorMap;
  ReactorMap reactor_; 

  struct LinPropFactory : public Fwk::InstanceFactory<LinearizedPropagator, const NlDynamTimeIntegrator *> {
    LinearizedPropagator::Ptr operator()(const NlDynamTimeIntegrator * key) const {
      return LinearizedPropagator::New(key);
    }
  };

  typedef Fwk::FactoryManagerImpl<LinPropFactory> LinPropMgr;
  LinPropMgr linearizedPropagator_;

  DISALLOW_COPY_AND_ASSIGN(Manager);
};

class IntegratorReactor : public DynamTimeIntegrator::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(IntegratorReactor);

  virtual void onCurrentCondition(); // overriden

  IntegratorReactor(const NlDynamTimeIntegrator * notifier, DynamStatePlainBasis * concurrentBasis, LinearizedPropagator * stepPropagator);

private:
  DynamStatePlainBasis::Ptr basis_;
  LinearizedPropagator::Ptr stepPropagator_;
};

class PropagatorReactor : public DynamPropagator::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagatorReactor);

  virtual void onInitialState(); // overriden
  virtual void onFinalState(); // overriden

  PropagatorReactor(const Propagator * notifier, DynamStatePlainBasis * concurrentBasis, Manager * parent);

private:
  Manager * parent_;
  DynamStatePlainBasis::Ptr basis_;
  Fwk::Ptr<IntegratorReactor> stepReactor_;
};

} /* end namespace ConcurrentBasis */

} /* end namespace Pita */

#endif /* PITA_CONCURRENTBASIS_MANAGER_H */
