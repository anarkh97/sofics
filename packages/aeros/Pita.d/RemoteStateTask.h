#ifndef PITA_REMOTESTATETASK_H
#define PITA_REMOTESTATETASK_H

#include "NamedTask.h"
#include "RemoteState.h"

namespace Pita {

class RemoteStateTask : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteStateTask);

  virtual void iterationIs(IterationRank i);

  const RemoteState::Activity * remoteActivity() const { return remoteActivity_.ptr(); }

  static Ptr New(const String & name, RemoteState::Activity * remoteActivity) {
    return new RemoteStateTask(name, remoteActivity);
  }

protected:
  RemoteStateTask(const String &, RemoteState::Activity *);

  RemoteState::Activity * remoteActivity() { return remoteActivity_.ptr(); }

private:
  RemoteState::Activity::Ptr remoteActivity_;
};

} /* end namespace Pita */

#endif /* PITA_REMOTESTATETASK_H */
