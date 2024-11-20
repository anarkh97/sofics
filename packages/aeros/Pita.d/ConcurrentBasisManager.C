#include "ConcurrentBasisManager.h"

namespace Pita {

namespace ConcurrentBasis {

Manager::Manager() :
  linearizedPropagator_(LinPropFactory())
{}

const Propagator *
Manager::referencePropagator(const DynamStatePlainBasis * concurrentBasis) const {
   ReactorMap::const_iterator it = reactor_.find(concurrentBasis);
   return (it != reactor_.end()) ? static_cast<const Propagator *>(it->second->notifier()) : NULL;
}

void
Manager::referencePropagatorIs(DynamStatePlainBasis * concurrentBasis, const Propagator * rp) {
  if (!concurrentBasis) return;
  if (!rp) {
    reactor_.erase(concurrentBasis);
  } else {
    reactor_[concurrentBasis] = new PropagatorReactor(rp, concurrentBasis, this);
  }
}

PropagatorReactor::PropagatorReactor(const Propagator * notifier, DynamStatePlainBasis * basis, Manager * parent) :
  DynamPropagator::Notifiee(notifier),
  parent_(parent),
  basis_(basis),
  stepReactor_(NULL)
{}

void
PropagatorReactor::onInitialState() {
  const Propagator * actualNotifier = static_cast<const Propagator *>(notifier()); // Safe cast
  NlDynamTimeIntegrator::PtrConst integrator = actualNotifier->integrator();
  
  LinearizedPropagator::Ptr stepPropagator = parent_->linearizedPropagator_.instance(integrator.ptr());
  if (!stepPropagator) {
    stepPropagator = parent_->linearizedPropagator_.instanceNew(integrator.ptr());
  }
   
  stepReactor_ = new IntegratorReactor(integrator.ptr(), basis_.ptr(), stepPropagator.ptr()); 
}

void
PropagatorReactor::onFinalState() {
  stepReactor_ = NULL;
}

IntegratorReactor::IntegratorReactor(const NlDynamTimeIntegrator * notifier,
                                     DynamStatePlainBasis * basis,
                                     LinearizedPropagator * stepPropagator) :
  DynamTimeIntegrator::NotifieeConst(notifier),
  basis_(basis),
  stepPropagator_(stepPropagator)
{}

void
IntegratorReactor::onCurrentCondition() {
  size_t stateCount = basis_->stateCount();
  for (size_t i = 0; i < stateCount; ++i) {
    // Ignore return value, use only side-effect only
    stepPropagator_->finalState(basis_->internalState(i));
  }
}

} /* end namespace ConcurrentBasis */

} /* end namespace Pita */
