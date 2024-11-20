#ifndef PITA_COMMSPLITTER_H
#define PITA_COMMSPLITTER_H

#include "Fwk.h"
#include "Types.h"

#include <Comm.d/Communicator.h>

namespace Pita {

class CommSplitter : public Fwk::PtrInterface<CommSplitter> {
public:
  EXPORT_PTRINTERFACE_TYPES(CommSplitter);

  enum CpuGroup {
    STANDARD    = 0,
    SPECIALIZED = 1
  };

  CpuGroup localGroup() const { return localGroup_; }

  Communicator * originComm() const { return originComm_; } // Original intracommunicator, not owned by CommSplitter
  Communicator * splitComm()  const { return splitComm_; }  // Local intracommunicator, owned by CommSplitter
  Communicator * interComm()  const { return interComm_; }  // Intercommunicator, owned by CommSplitter

  // Creates 2 intracommunicators from originComm linked by an intercommunicator
  // The communicator originComm must remain valid while CommSplitter exists
  static Ptr New(Communicator * originComm, CpuCount specializedCpus) {
    return new CommSplitter(originComm, specializedCpus);
  }

protected:
  explicit CommSplitter(Communicator * originComm, CpuCount specializedCpus);

  ~CommSplitter();

private:
  CpuGroup localGroup_;
  Communicator * originComm_;
  Communicator * splitComm_;
  Communicator * interComm_;
};

} /* end namespace Pita */

#endif /* PITA_COMMSPLITTER_H */
