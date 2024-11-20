#include "RemoteDynamPropagatorProxy.h"

namespace Pita {

RemoteDynamPropagatorProxy::RemoteDynamPropagatorProxy(
    size_t vectorSize, Communicator * serverCommunicator, CpuRank serverCpu) :
  DynamPropagator(vectorSize),
  serverCommunicator_(serverCommunicator),
  serverCpu_(serverCpu)
{}

void
RemoteDynamPropagatorProxy::initialStateIs(const DynamState & is) {
  setInitialState(is);
  initialStateNotify();

  size_t stateSize = 2 * this->vectorSize(); 
  sBuffer_.sizeIs(stateSize);
  bufferStateCopy(is, sBuffer_.array());

  serverCommunicator_->sendTo(serverCpu_.value(), serverCommunicator_->myID(), sBuffer_.array(), stateSize);
  serverCommunicator_->waitForAllReq();
  serverCommunicator_->recFrom(serverCommunicator_->myID(), sBuffer_.array(), stateSize);

  DynamState receivedState = DynamState(vectorSize(), sBuffer_.array());
  setFinalState(receivedState);
  finalStateNotify();
}

} // end namespace Pita
