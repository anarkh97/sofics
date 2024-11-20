#include "RemoteSeedInitializerServer.h"

#include <Comm.d/Communicator.h>

namespace Pita { namespace Hts {

RemoteSeedInitializerServer::RemoteSeedInitializerServer(Communicator * cc, SeedInitializer * si, const SliceMapping * m) :
  baseInitializer_(si),
  clientCommunicator_(cc),
  sBuffer_(),
  clientMapping_(m)
{}

void
RemoteSeedInitializerServer::statusIs(RemoteCoarseServer::Status s) {
  if (status() == s)
    return;

  if (s == BUSY) {
    setStatus(BUSY);

    int seedCount = (clientMapping()->activeSlices().value() / 2) + 1;
    size_t stateSize = 2 * baseInitializer_->vectorSize();
    sBuffer_.sizeIs(stateSize * seedCount);
 
    for (int s = 0; s < seedCount; ++s) {
      DynamState seed = baseInitializer_->initialSeed(SliceRank(s));
      double * buffer = sBuffer_.array() + s * stateSize;
      bufferStateCopy(seed, buffer);
     
      CpuRank backwardTarget(-1);
      if (s != 0) {
        backwardTarget = clientMapping()->hostCpu(HalfSliceRank(s * 2 - 1));
        clientCommunicator()->sendTo(backwardTarget.value(), s, buffer, stateSize);
      }

      CpuRank forwardTarget(-1);
      if (s != seedCount - 1) {
        forwardTarget = clientMapping()->hostCpu(HalfSliceRank(s * 2));
        if (forwardTarget != backwardTarget) {
          clientCommunicator()->sendTo(forwardTarget.value(), s, buffer, stateSize);
        }
      }
    }

    clientCommunicator()->waitForAllReq();
  }
  
  setStatus(READY);
}

} /* end namespace Hts */ } /* end namespace Hts */
