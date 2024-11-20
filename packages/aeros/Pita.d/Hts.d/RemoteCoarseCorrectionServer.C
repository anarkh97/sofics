#include "RemoteCoarseCorrectionServer.h"

namespace Pita { namespace Hts {

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

    // Only serve the even-numbered (starting from first active) active half-slices corresponding to a complete full time-slice
    // That is, [firstActive, firstActive + 2, ..., firstActive + 2*n] where firstActive + 2*n + 1 < firstInactive and firstActive + 2*n + 3 >= firstInactive
    HalfSliceRank firstInactive = clientMapping()->firstInactiveSlice();
    for (HalfSliceRank r = clientMapping()->firstActiveSlice(); r + HalfSliceCount(1) < firstInactive; r = r + HalfSliceCount(2)) {
      CpuRank clientCpu = clientMapping()->hostCpu(r);
      server_->initialStateNew(clientCpu);
    }
  }

  setStatus(READY);
}


} /* end namespace Hts */ } /* end namespace Pita */
