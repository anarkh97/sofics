#include "NlDriverImpl.h"
#include "PitaNonLinDynam.h"

namespace Pita {

DynamState
NlDriverImpl::initialState() const {
  DynamState result(probDesc()->solVecInfo());

  GenVector<double> dummy_acc(probDesc()->solVecInfo());
  GenVector<double> dummy_vp(probDesc()->solVecInfo());
  probDesc()->getInitState(result.displacement(), result.velocity(),
                           dummy_acc, dummy_vp);

  return result;
}

} // end namespace Pita
