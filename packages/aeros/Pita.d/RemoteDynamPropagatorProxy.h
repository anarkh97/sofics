#ifndef PITA_REMOTEPROPAGATORPPROXY_H
#define PITA_REMOTEPROPAGATORPPROXY_H

#include "DynamPropagator.h"
#include "SimpleBuffer.h"
#include <Comm.d/Communicator.h>

namespace Pita {

class RemoteDynamPropagatorProxy : public DynamPropagator {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteDynamPropagatorProxy);

  // Added
  Communicator * serverCommunicator() const { return serverCommunicator_; }
  CpuRank serverCpu() const { return serverCpu_; }
  
  // Overriden
  virtual void initialStateIs(const DynamState & is);

  static Ptr New(size_t vectorSize, Communicator * serverCommunicator, CpuRank serverCpu) {
    return new RemoteDynamPropagatorProxy(vectorSize, serverCommunicator, serverCpu);
  }

protected:
  RemoteDynamPropagatorProxy(size_t vectorSize, Communicator * serverCommunicator, CpuRank serverCpu);
  
private:
  Communicator * serverCommunicator_;
  CpuRank serverCpu_; 

  SimpleBuffer<double> sBuffer_;
};

} // end namespace Pita

#endif /* PITA_REMOTEPROPAGATORPPROXY_H */
