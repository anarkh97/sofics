#include "RemoteStateTask.h"

namespace Pita {

RemoteStateTask::RemoteStateTask(
    const String & name,
     RemoteState::Activity * remoteActivity) :
  NamedTask(name),
  remoteActivity_(remoteActivity)
{}

void
RemoteStateTask::iterationIs(IterationRank i) {
  remoteActivity_->statusIs(RemoteState::EXECUTING);
  setIteration(i);
}

} // end namespace Pita
