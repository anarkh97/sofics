#include "RemoteDynamPropagatorServer.h"

namespace Pita {

RemoteDynamPropagatorServer::RemoteDynamPropagatorServer(
    DynamPropagator * propagator, Communicator * clientCommunicator) :
  propagator_(propagator),
  clientCommunicator_(clientCommunicator),
  sBuffer_()
{}

void
RemoteDynamPropagatorServer::initialStateNew(CpuRank clientCpu) {
  size_t stateSize = propagator()->vectorSize() * 2;
  sBuffer_.sizeIs(stateSize);

  clientCommunicator_->recFrom(clientCpu.value(), sBuffer_.array(), stateSize);
  DynamState receivedState = DynamState(propagator()->vectorSize(), sBuffer_.array());

  propagator()->initialStateIs(receivedState);
  bufferStateCopy(propagator()->finalState(), sBuffer_.array());

  clientCommunicator_->sendTo(clientCpu.value(), clientCpu.value(), sBuffer_.array(), stateSize);
  clientCommunicator_->waitForAllReq();
}

} // end namespace Pita
