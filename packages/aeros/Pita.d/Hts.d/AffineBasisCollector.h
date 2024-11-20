#ifndef PITA_HTS_AFFINEBASISCOLLECTOR_H
#define PITA_HTS_AFFINEBASISCOLLECTOR_H

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceId.h"

#include "../DynamState.h"
#include "../AffineDynamPropagator.h"

#include <stack>

namespace Pita { namespace Hts {

class AffineBasisCollector : public Fwk::PtrInterface<AffineBasisCollector> {
public: 
  EXPORT_PTRINTERFACE_TYPES(AffineBasisCollector);
  
  // Data collection  
  size_t sourceCount() const;
  const AffineDynamPropagator * source(const HalfSliceId & sliceId) const;
  void sourceIs(const HalfSliceId & sliceId, const AffineDynamPropagator * source);

  // Data access
  typedef std::pair<HalfSliceRank, DynamState> CollectedState;

  CollectedState firstForwardFinalState() const;
  void firstForwardFinalStateDel();
  CollectedState firstBackwardFinalState() const;
  void firstBackwardFinalStateDel();
  
  void finalStateIs(const HalfSliceId & sliceId, const DynamState & state); 

  static Ptr New() {
    return new AffineBasisCollector();
  }

protected:
  AffineBasisCollector();
 
  class PropagationReactor; 
  typedef std::map<HalfSliceId, Fwk::Ptr<PropagationReactor> > PropagationReactorContainer;
  
  const PropagationReactorContainer & propagationReactor() const {
    return propagationReactor_;
  }

  virtual PropagationReactor * propagationReactorNew(const AffineDynamPropagator * notifier,
                                                     const HalfSliceId & id);

private:
  PropagationReactorContainer propagationReactor_;
  
  typedef std::stack<CollectedState> StateContainer;
  StateContainer forwardFinalState_;
  StateContainer backwardFinalState_;
};

class AffineBasisCollector::PropagationReactor : public DynamPropagator::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationReactor);
  
  PropagationReactor(const AffineDynamPropagator * notifier, const HalfSliceId & id, AffineBasisCollector * parent);

  virtual void onFinalState(); // Overriden

  const HalfSliceId & sliceId() const { return sliceId_; }
  AffineBasisCollector * parent() const { return parent_; }

private:
  HalfSliceId sliceId_;
  AffineBasisCollector * parent_;
};

} /* end notifier Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_AFFINEBASISCOLLECTOR_H */
