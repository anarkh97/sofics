#include "CommSplitter.h"

namespace Pita {

CommSplitter::CommSplitter(Communicator * originComm, CpuCount specializedCpus) :
  localGroup_(STANDARD),
  originComm_(originComm),
  splitComm_(NULL),
  interComm_(NULL)
{
  // Determine local context
  CpuRank firstStandardCpu(0); 
  CpuRank firstSpecializedCpu = CpuRank(originComm_->numCPUs()) - specializedCpus;
  if (CpuRank(originComm_->myID()) >= firstSpecializedCpu) {
    localGroup_ = SPECIALIZED;
  }

  // Private common communicator
  MPI_Comm peerMComm;
  MPI_Comm_dup(*originComm->getCommunicator(), &peerMComm);
  int peerId;
  MPI_Comm_rank(peerMComm, &peerId);

  // Create non-overlapping communicators
  MPI_Comm splitMComm;
  MPI_Comm_split(peerMComm, localGroup_, peerId, &splitMComm);

  // Create intercommunicator
  MPI_Comm interMComm;
  CpuRank remoteLeader = (localGroup_ == SPECIALIZED) ? firstStandardCpu : firstSpecializedCpu;
  const int localLeader = 0;
  const int msgTag = 0;
  MPI_Intercomm_create(splitMComm, localLeader, peerMComm, remoteLeader.value(), msgTag, &interMComm);

  // Release private communicator
  MPI_Comm_free(&peerMComm); 

  // Publish results
  splitComm_ = new Communicator(splitMComm, stderr);
  interComm_ = new Communicator(interMComm, stderr);
}

CommSplitter::~CommSplitter() {
  // Mark MPI Communicator for deletion
  MPI_Comm_free(interComm_->getCommunicator());
  MPI_Comm_free(splitComm_->getCommunicator());
  
  delete interComm_;
  delete splitComm_;
}

} /* end namespace Pita */
