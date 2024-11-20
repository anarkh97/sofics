#ifndef PITA_STD_AFFINEBASISCOLLECTOR_H
#define PITA_STD_AFFINEBASISCOLLECTOR_H

#include "Fwk.h"
#include "Types.h"

#include "../DynamState.h"
#include "../AffineDynamPropagator.h"

#include <queue>

namespace Pita { namespace Std {

class AffineBasisCollector : public Fwk::PtrInterface<AffineBasisCollector> {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineBasisCollector);

  // Data collection  
  size_t sourceCount() const;
  const AffineDynamPropagator * source(SliceRank sliceId) const;
  void sourceIs(SliceRank sliceId, const AffineDynamPropagator * source);

  // Data access
  struct CollectedState {
    SliceRank sliceId;
    DynamState state;
    
    CollectedState(SliceRank sl = SliceRank(-1), DynamState st = DynamState()) :
      sliceId(sl), state(st)
    {}
  };

  CollectedState firstInitialState() const;
  void firstInitialStateDel();
  CollectedState firstFinalState() const;
  void firstFinalStateDel();

  static Ptr New() {
    return new AffineBasisCollector();
  }

protected:
  AffineBasisCollector();
  ~AffineBasisCollector();

private:
  class PropagationReactor;
  friend class PropagationReactor;
  typedef std::map<SliceRank, Fwk::Ptr<PropagationReactor> > ReactorMap;
  ReactorMap reactor_;

  typedef std::queue<CollectedState> StateQueue;
  StateQueue initialState_;
  StateQueue finalState_;

  DISALLOW_COPY_AND_ASSIGN(AffineBasisCollector);
};


} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_AFFINEBASISCOLLECTOR_H */
