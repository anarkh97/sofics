#ifndef PITA_STD_REMOTECOARSECORRECTIONSERVER_H
#define PITA_STD_REMOTECOARSECORRECTIONSERVER_H

#include "Fwk.h"
#include "Types.h"

#include "../RemoteCoarseServer.h" 
#include "SliceMapping.h"

#include "../RemoteDynamPropagatorServer.h"

namespace Pita { namespace Std {

class RemoteCoarseCorrectionServer : public RemoteCoarseServer {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteCoarseCorrectionServer);

  virtual void statusIs(Status s); // overriden
  
  RemoteDynamPropagatorServer * server() const { return server_.ptr(); }
  const SliceMapping * clientMapping() const { return clientMapping_.ptr(); }

  // Caution: mapping must be kept up-to-date 
  static Ptr New(RemoteDynamPropagatorServer * server, const SliceMapping * mapping) {
    return new RemoteCoarseCorrectionServer(server, mapping);
  }

protected:
  RemoteCoarseCorrectionServer(RemoteDynamPropagatorServer * server, const SliceMapping * clientMapping);

private:
  RemoteDynamPropagatorServer::Ptr server_;
  SliceMapping::PtrConst clientMapping_;

  DISALLOW_COPY_AND_ASSIGN(RemoteCoarseCorrectionServer);
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_REMOTECOARSECORRECTIONSERVER_H */
