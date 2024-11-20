#include "RemoteCoarseCorrectionServer.h"

namespace Pita { namespace Std {

RemoteCoarseCorrectionServer::RemoteCoarseCorrectionServer(
    RemoteDynamPropagatorServer * server,
    const SliceMapping * clientMapping) :
  server_(server),
  clientMapping_(clientMapping)
{}

void
RemoteCoarseCorrectionServer::statusIs(Status s) {
  if (status() == s)
    return;

  if (s == BUSY) {
    setStatus(BUSY);

    SliceRank firstInactive = clientMapping()->firstInactiveSlice();
    for (SliceRank r = clientMapping()->firstActiveSlice(); r < firstInactive; r = r.next()) {
      CpuRank clientCpu = clientMapping()->hostCpu(r);
      server_->initialStateNew(clientCpu);
    }
  }

  setStatus(READY);
}


} /* end namespace Std */ } /* end namespace Pita */
