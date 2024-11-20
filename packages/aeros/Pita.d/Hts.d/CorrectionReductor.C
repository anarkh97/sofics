#include "CorrectionReductor.h"

namespace Pita { namespace Hts {

CorrectionReductor::CorrectionReductor(const String & name, const DynamStateReductor * reductor) :
  CorrectionPropagator<DynamState, Vector>(name),
  reductor_(reductor)
{}

void
CorrectionReductor::iterationIs(IterationRank ir) {
  assert(jump()->iteration().next() == ir);
  assert(correction()->iteration() == jump()->iteration() || correction()->status() == Seed::INACTIVE);
  assert(jump()->status() == Seed::CONVERGED || correction()->status() != Seed::INACTIVE);

  DynamState update = jump()->state();
  if (correction()->status() != Seed::INACTIVE) {
    update += correction()->state();
  }

  nextCorrection()->stateIs(reductor()->result(update));
  nextCorrection()->statusIs(correction()->status() == Seed::ACTIVE ? Seed::ACTIVE : jump()->status());
  nextCorrection()->iterationIs(jump()->iteration());

  setIteration(ir);
}

CorrectionReductor::Manager::Manager(OperatorManager * reductorManager) :
  impl_(Factory(reductorManager))
{}

CorrectionReductor *
CorrectionReductor::Manager::Factory::operator()(HalfSliceRank key) {
  DynamStateReductor * reductor = manager_->instanceNew(key);
  String taskName = String("CorrectionReductor ") + toString(key);
  return new CorrectionReductor(taskName, reductor);
}

} /* end namespace Hts */ } /* end namespace Pita */
