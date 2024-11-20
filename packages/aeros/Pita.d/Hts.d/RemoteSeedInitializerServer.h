#ifndef PITA_HTS_REMOTESEEDINITIALIZERSERVER_H
#define PITA_HTS_REMOTESEEDINITIALIZERSERVER_H

#include "../RemoteCoarseServer.h"
#include "SliceMapping.h"

#include "../SeedInitializer.h"
#include "../SimpleBuffer.h"

class Communicator;

namespace Pita { namespace Hts {

class RemoteSeedInitializerServer : public RemoteCoarseServer {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSeedInitializerServer);

  virtual void statusIs(Status s); // overriden
 
  SeedInitializer * baseInitializer() const { return baseInitializer_.ptr(); }
  const SliceMapping * clientMapping() const { return clientMapping_.ptr(); }
  Communicator * clientCommunicator() const { return clientCommunicator_; }

  static Ptr New(Communicator * clientCommunicator, SeedInitializer * baseInitializer, const SliceMapping * mapping) {
    return new RemoteSeedInitializerServer(clientCommunicator, baseInitializer, mapping);
  }

protected:
  RemoteSeedInitializerServer(Communicator * cc, SeedInitializer * si, const SliceMapping * m);

private:
  SeedInitializer::Ptr baseInitializer_;
  Communicator * clientCommunicator_;
  SimpleBuffer<double> sBuffer_;
  
  SliceMapping::PtrConst clientMapping_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_REMOTESEEDINITIALIZERSERVER_H */
