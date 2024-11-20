#include "PostProcessingManager.h"

// HACK ...
#include "AffineIntegratorPropagator.h"
#include "IntegratorPropagator.h"

namespace Pita {

namespace PostProcessing {

// PropagatorReactor::Manager implementation

PropagatorReactor::Manager::Manager(IntegratorReactor::Builder * reactorBuilder) :
  reactorBuilder_(reactorBuilder),
  propagatorReactor_()
{}

PostProcessor::FileSetId
PropagatorReactor::Manager::outputFileSet(const DynamPropagator * op) const {
  PropagatorReactorContainer::const_iterator it = propagatorReactor_.find(op);
  return it != propagatorReactor_.end() ? it->second->fileSet() : PostProcessor::FileSetId();
}

void
PropagatorReactor::Manager::outputFileSetIs(const DynamPropagator * op, PostProcessor::FileSetId fs) {
  PropagatorReactorContainer::iterator it = propagatorReactor_.lower_bound(op);
  if (it != propagatorReactor_.end() && it->first == op) {
    // The propagator is already associated with a FileSet
    if (fs != PostProcessor::FileSetId()) {
      // Update FileSet association
      it->second->fileSetIs(fs);
    } else {
      // Default value corresponds to no fileset, stop observing propagator 
      propagatorReactor_.erase(it);
    }
  } else {
    if (fs != PostProcessor::FileSetId()) {
      // Start observing propagator with associated FileSet
      PropagatorReactor::Ptr reactor = new PropagatorReactor(op, this, fs);
      propagatorReactor_.insert(it, std::make_pair(op, reactor));
    }
  }
}

// PropagatorReactor implementation

PropagatorReactor::PropagatorReactor(const DynamPropagator * notifier,
                                     PropagatorReactor::Manager * parent,
                                     PostProcessor::FileSetId fileSet) :
  DynamPropagator::Notifiee(NULL),
  notifier_(notifier),
  parent_(parent),
  fileSet_(fileSet),
  integratorReactor_(NULL)
{
  this->notifierIs(notifier);
}

void
PropagatorReactor::onInitialState() {
  if (const AffineIntegratorPropagator * downcasted = dynamic_cast<const AffineIntegratorPropagator *>(notifier_)) {
    if (const DynamTimeIntegrator * integrator = downcasted->integrator()) {
      integratorReactor_ = parent_->reactorBuilder()->reactorNew(integrator, fileSet_);
    }
  }
  
  if (const IntegratorPropagatorRoot * downcasted = dynamic_cast<const IntegratorPropagatorRoot *>(notifier_)) {
    if (const DynamTimeIntegrator * integrator = downcasted->integrator()) {
      integratorReactor_ = parent_->reactorBuilder()->reactorNew(integrator, fileSet_);
    }
  }
}

void
PropagatorReactor::onFinalState() {
  integratorReactor_ = NULL;
  parent_->reactorBuilder()->postProcessor()->fileStatusIs(this->fileSet(), PostProcessor::CLOSED);
}

} // end namespace PostProcessing

} // end namespace Pita
