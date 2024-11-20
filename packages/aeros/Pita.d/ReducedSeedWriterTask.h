#ifndef PITA_REDUCEDSEEDWRITERTASK_H
#define PITA_REDUCEDSEEDWRITERTASK_H

#include "RemoteStateTask.h"
#include "RemoteStateMpiImpl.h"

namespace Pita {

template <typename Reconstructor>
class ReducedSeedWriterTask : public RemoteStateTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedSeedWriterTask);

  virtual void iterationIs(IterationRank i); // overriden

  static Ptr New(const String & name, RemoteState::MpiReducedSeedWriter * remoteActivity, const Reconstructor * reconstructor) {
    return new ReducedSeedWriterTask(name, remoteActivity, reconstructor);
  }

protected:
  ReducedSeedWriterTask(const String &, RemoteState::Activity *, const Reconstructor * reconstructor);

private:
  Fwk::Ptr<const Reconstructor> reconstructor_;
};

template <typename Reconstructor>
ReducedSeedWriterTask<Reconstructor>::ReducedSeedWriterTask(const String & name, RemoteState::Activity * activity, const Reconstructor * reconstructor) :
  RemoteStateTask(name, activity),
  reconstructor_(reconstructor)
{}

template <typename Reconstructor>
void
ReducedSeedWriterTask<Reconstructor>::iterationIs(IterationRank i) {
  RemoteState::MpiReducedSeedWriter * writer = static_cast<RemoteState::MpiReducedSeedWriter *>(remoteActivity()); // Safe cast
  writer->reducedStateSizeIs(reconstructor_->reducedBasisSize());

  RemoteStateTask::iterationIs(i);
}

} /* end namespace Pita */

#endif /* PITA_REMOTESTATETASK_H */
