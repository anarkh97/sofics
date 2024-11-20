#include "UpdatedSeedAssemblerImpl.h"

#include <memory>
#include <cassert>

namespace Pita {

/* UpdatedSeedAssemblerImpl */

UpdatedSeedAssemblerImpl::UpdatedSeedAssemblerImpl(const String & name, const Manager * manager) :
  UpdatedSeedAssembler(name),
  manager_(manager)
{}

void
UpdatedSeedAssemblerImpl::iterationIs(IterationRank ir) {
  assert(correctionComponents()->iteration() == propagatedSeed()->iteration() || correctionComponents()->status() == Seed::INACTIVE);
 
  if (correctionComponents()->status() != Seed::INACTIVE) {
    const Vector & components = correctionComponents()->state();

    DynamState result(propagatedSeed()->state().vectorSize(), 0.0);
    int rbs = static_cast<int>(correctionBasis()->stateCount());
    for (int i = 0; i < rbs; ++i) {
      if (components[i] != 0.0) {
        result.linAdd(components[i], correctionBasis()->state(i));
      }
    }

    correction()->stateIs(result);

    correction()->statusIs(correctionComponents()->status());
    correction()->iterationIs(correctionComponents()->iteration());
  }

  updateSeed();

  setIteration(ir);
}

/* Manager */

UpdatedSeedAssemblerImpl::Manager::Manager(const DynamStateBasis * dcb) :
  defaultCorrectionBasis_(dcb)
{}

UpdatedSeedAssemblerImpl *
UpdatedSeedAssemblerImpl::Manager::createNewInstance(const String & key) {
  String taskName = String("Update Seed ") + key;
  std::unique_ptr<UpdatedSeedAssemblerImpl> newInstance(new UpdatedSeedAssemblerImpl(taskName, this));
  return newInstance.release();
}

} /* end namespace Pita */
