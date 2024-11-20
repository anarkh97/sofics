#ifndef PITA_REMOTEDYNAMPROPAGATORSERVER_H
#define PITA_REMOTEDYNAMPROPAGATORSERVER_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "DynamPropagator.h"
#include <Comm.d/Communicator.h>

#include "SimpleBuffer.h"

namespace Pita {

class RemoteDynamPropagatorServer : public Fwk::PtrInterface<RemoteDynamPropagatorServer> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteDynamPropagatorServer);

  DynamPropagator * propagator() const { return propagator_.ptr(); }
  Communicator * clientCommunicator() const { return clientCommunicator_; }

  void initialStateNew(CpuRank clientCpu);

  RemoteDynamPropagatorServer(DynamPropagator * propagator,
      Communicator * clientCommunicator); 

private:
  DynamPropagator::Ptr propagator_;
  Communicator * clientCommunicator_;

  SimpleBuffer<double> sBuffer_;

  DISALLOW_COPY_AND_ASSIGN(RemoteDynamPropagatorServer);
};

} // end namespace Pita

#endif /* PITA_REMOTEDYNAMPROPAGATORSERVER_H */
