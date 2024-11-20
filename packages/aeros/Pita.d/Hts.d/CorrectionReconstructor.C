#include "CorrectionReconstructor.h"

namespace Pita { namespace Hts {

void
CorrectionReconstructor::iterationIs(IterationRank ir) {
  if (correctionComponents()->status() != Seed::INACTIVE) {
    correction()->stateIs(reconstructor_->result(correctionComponents()->state()));
  }
  correction()->statusIs(correctionComponents()->status());
  correction()->iterationIs(correctionComponents()->iteration());
}

CorrectionReconstructor::Manager::Manager(OperatorManager * reconstructorManager) :
  impl_(Factory(reconstructorManager))
{}

CorrectionReconstructor *
CorrectionReconstructor::Manager::Factory::operator()(HalfSliceRank key) {
  DynamStateReconstructor * reconstructor = manager_->instanceNew(key);
  String taskName = String("CorrectionReconstructor ") + toString(key);
  return new CorrectionReconstructor(taskName, reconstructor);
}

} /* end namespace Hts */ } /* end namespace Pita */
