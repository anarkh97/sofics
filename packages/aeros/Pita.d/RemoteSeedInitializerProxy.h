#ifndef PITA_REMOTESEEDINITIALIZERPROXY_H
#define PITA_REMOTESEEDINITIALIZERPROXY_H

#include "SeedInitializer.h"

#include "SimpleBuffer.h"

#include <Comm.d/Communicator.h>

namespace Pita {

class RemoteSeedInitializerProxy : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSeedInitializerProxy);

  virtual DynamState initialSeed(SliceRank rank) const; // Overriden

  const Communicator * communicator() const { return communicator_; }

  static Ptr New(Communicator * communicator, size_t vectorSize) {
    return new RemoteSeedInitializerProxy(communicator, vectorSize);
  }

protected:
  RemoteSeedInitializerProxy(Communicator * c, size_t vs);

private:
  Communicator * communicator_;
  mutable SimpleBuffer<double> sBuffer_;

  typedef std::map<SliceRank, DynamState> StateMap;
  mutable StateMap state_;
};

} // end namespace Pita

#endif /* PITA_REMOTESEEDINITIALIZERPROXY_H */
