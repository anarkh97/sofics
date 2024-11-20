#include "AffineBasisCollector.h"

namespace Pita { namespace Hts {

AffineBasisCollector::AffineBasisCollector() {}

const AffineDynamPropagator *
AffineBasisCollector::source(const HalfSliceId & sliceId) const {
  PropagationReactorContainer::const_iterator it = propagationReactor_.find(sliceId);
  return (it != propagationReactor_.end()) ? static_cast<const AffineDynamPropagator *>(it->second->notifier()) : NULL;
}

size_t
AffineBasisCollector::sourceCount() const {
  return propagationReactor_.size();
}

void
AffineBasisCollector::sourceIs(const HalfSliceId & sliceId, const AffineDynamPropagator * source) {
  if (source == NULL) {
    // Remove current
    propagationReactor_.erase(sliceId);
    return;
  }

  PropagationReactor::Ptr newReactor = propagationReactorNew(source, sliceId);
  
  // Find insertion point
  PropagationReactorContainer::iterator it = propagationReactor_.lower_bound(sliceId);
  if (it != propagationReactor_.end() && it->first == sliceId) {
    it->second = newReactor; // Replace previous
  } else {
    propagationReactor_.insert(it, std::make_pair(sliceId, newReactor)); // Insert new
  }
}

AffineBasisCollector::CollectedState
AffineBasisCollector::firstForwardFinalState() const {
  return !forwardFinalState_.empty() ? forwardFinalState_.top() : std::make_pair(HalfSliceRank(-1), DynamState());
}

void
AffineBasisCollector::firstForwardFinalStateDel() {
  forwardFinalState_.pop();
}

AffineBasisCollector::CollectedState
AffineBasisCollector::firstBackwardFinalState() const {
  return !backwardFinalState_.empty() ? backwardFinalState_.top() : std::make_pair(HalfSliceRank(-1), DynamState());
}

void
AffineBasisCollector::firstBackwardFinalStateDel() {
  backwardFinalState_.pop();
}

void
AffineBasisCollector::finalStateIs(const HalfSliceId & sliceId, const DynamState & state) {
  Direction dir = sliceId.direction();
  StateContainer & stateContainer = (dir == FORWARD) ? forwardFinalState_ : backwardFinalState_;
  stateContainer.push(std::make_pair(sliceId.rank(), state));
}

AffineBasisCollector::PropagationReactor *
AffineBasisCollector::propagationReactorNew(const AffineDynamPropagator * notifier,
                                            const HalfSliceId & id) {
  return new PropagationReactor(notifier, id, this);
}

AffineBasisCollector::PropagationReactor::PropagationReactor(const AffineDynamPropagator * notifier,
                                                             const HalfSliceId & id,
                                                             AffineBasisCollector * parent) :
  DynamPropagator::Notifiee(notifier),
  sliceId_(id),
  parent_(parent)
{}

void
AffineBasisCollector::PropagationReactor::onFinalState() {
  const AffineDynamPropagator * downcasted = static_cast<const AffineDynamPropagator *>(notifier());

  // Interested only in linear propagation
  if (downcasted->constantTerm() == AffineDynamPropagator::HOMOGENEOUS) {
    parent_->finalStateIs(sliceId_, notifier()->finalState());
  }
}

} /* end notifier Hts */ } /* end namespace Pita */
