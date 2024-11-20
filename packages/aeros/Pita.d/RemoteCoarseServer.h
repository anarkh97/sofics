#ifndef PITA_REMOTECOARSESERVER_H
#define PITA_REMOTECOARSESERVER_H

#include "Fwk.h"

namespace Pita {

class RemoteCoarseServer : public Fwk::PtrInterface<RemoteCoarseServer> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteCoarseServer);

  enum Status {
    READY = 0,
    BUSY
  };

  Status status() const { return status_; }
  virtual void statusIs(Status s) = 0;

protected:
  explicit RemoteCoarseServer() :
    status_(READY)
  {}

  void setStatus(Status s) { status_ = s; }

private:
  Status status_;
};

} /* end namespace Pita */

#endif /* PITA_REMOTECOARSESERVER_H */
