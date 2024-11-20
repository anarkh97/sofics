#include "AffineBasisCollector.h"

namespace Pita { namespace Std {

// Notification (Auxiliary class)

class AffineBasisCollector::PropagationReactor : public DynamPropagator::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationReactor);

  virtual void onFinalState(); // overriden

  PropagationReactor(SliceRank sliceId, const AffineDynamPropagator * notifier, AffineBasisCollector * parent) :
    DynamPropagator::Notifiee(notifier), // Notifier can be downcasted back safely
    sliceId_(sliceId),
    parent_(parent)
  {}

private:
  SliceRank sliceId_;
  AffineBasisCollector * parent_;
}; 

void
AffineBasisCollector::PropagationReactor::onFinalState() {
  const AffineDynamPropagator * downcasted = static_cast<const AffineDynamPropagator *>(notifier());

  // Interested only in linear propagation
  if (downcasted->constantTerm() == AffineDynamPropagator::HOMOGENEOUS) {
    parent_->initialState_.push(CollectedState(sliceId_, downcasted->initialState()));
    parent_->finalState_.push(CollectedState(sliceId_, downcasted->finalState()));
  }
}


// Constructor

AffineBasisCollector::AffineBasisCollector() {}

// Destructor

AffineBasisCollector::~AffineBasisCollector() {}

// Data collection

size_t
AffineBasisCollector::sourceCount() const {
  return reactor_.size(); 
}

const AffineDynamPropagator *
AffineBasisCollector::source(SliceRank sliceId) const {
  ReactorMap::const_iterator it = reactor_.find(sliceId);
  return (it != reactor_.end()) ? static_cast<const AffineDynamPropagator *>(it->second->notifier()) : NULL;
}

void
AffineBasisCollector::sourceIs(SliceRank sliceId, const AffineDynamPropagator * source) {
  if (source == NULL) {
    // Remove current
    reactor_.erase(sliceId);
    return;
  }

  // Find insertion point 
  ReactorMap::iterator it = reactor_.lower_bound(sliceId);
  if (it != reactor_.end() && it->first == sliceId) {
    // Replace previous
    it->second->notifierIs(source);
  } else {
    // Insert new 
    PropagationReactor::Ptr newReactor = new PropagationReactor(sliceId, source, this);
    reactor_.insert(it, std::make_pair(sliceId, newReactor));
  }
}


// Data access

AffineBasisCollector::CollectedState
AffineBasisCollector::firstInitialState() const {
  return !initialState_.empty() ? initialState_.front() : CollectedState();
}

void
AffineBasisCollector::firstInitialStateDel() {
  initialState_.pop();
}

AffineBasisCollector::CollectedState
AffineBasisCollector::firstFinalState() const {
  return !finalState_.empty() ? finalState_.front() : CollectedState();
}

void
AffineBasisCollector::firstFinalStateDel() {
  finalState_.pop();
}

} /* end namespace Std */ } /* end namespace Pita */
