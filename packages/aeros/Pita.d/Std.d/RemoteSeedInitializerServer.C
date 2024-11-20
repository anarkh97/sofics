#include "RemoteSeedInitializerServer.h"

#include <Comm.d/Communicator.h>

namespace Pita { namespace Std {

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

    int seedCount = clientMapping()->activeSlices().value();
    size_t stateSize = 2 * baseInitializer_->vectorSize();
    sBuffer_.sizeIs(stateSize * seedCount);
 
    for (int s = 0; s < seedCount; ++s) {
      DynamState seed = baseInitializer_->initialSeed(SliceRank(s));
      double * buffer = sBuffer_.array() + s * stateSize;
      bufferStateCopy(seed, buffer);
    
      CpuRank client = clientMapping()->hostCpu(SliceRank(s));
      clientCommunicator()->sendTo(client.value(), s, buffer, stateSize);
    }

    clientCommunicator()->waitForAllReq();
  }
  
  setStatus(READY);
}

} /* end namespace Std */ } /* end namespace Std */
