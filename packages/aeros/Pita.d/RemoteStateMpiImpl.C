#include "RemoteStateMpiImpl.h"

#include "SimpleBuffer.h"
#include <Comm.d/Communicator.h>

namespace Pita { namespace RemoteState {

/* MpiSeedReader implementation */

MpiSeedReader::MpiSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer) :
  Reader<DynamState>(targetCpu),
  communicator_(comm),
  stateBuffer_(buffer)
{}

CpuRank
MpiSeedReader::localCpu() const {
  return CpuRank(communicator_->myID());
}

void
MpiSeedReader::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiSeedReader::statusIs -- Not ready");

  setStatus(s);

  if (status() != EXECUTING)
    return;

  if (localCpu() != targetCpu() && targetCpu() != nullCpu) {
    DynamState state = origin()->state();

    size_t dim = 2 * state.vectorSize() + 2;
    bufferStateCopy(state, stateBuffer().array());

    stateBuffer()[dim - 2] = static_cast<double>(origin()->status()); 
    stateBuffer()[dim - 1] = static_cast<double>(origin()->iteration().value());

    int messageTag = localCpu().value(); // TODO determine meaningful messageId
    communicator()->sendTo(targetCpu().value(), messageTag, stateBuffer().array(), dim);
    communicator()->waitForAllReq(); // TODO delay waitForAllReq
  }

  setStatus(READY);
}

/* MpiSeedWriter implementation */

MpiSeedWriter::MpiSeedWriter(CpuRank originCpu, Communicator * comm, size_t vSize, SimpleBuffer<double> * buffer) :
  Writer<DynamState>(originCpu),
  communicator_(comm),
  vectorSize_(vSize),
  stateBuffer_(buffer)
{}

CpuRank
MpiSeedWriter::localCpu() const {
  return CpuRank(communicator_->myID());
}

void
MpiSeedWriter::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiSeedWriter::statusIs -- Not ready");

  setStatus(s);
  
  if (status() != EXECUTING)
    return;

  if (originCpu() != localCpu() && originCpu() != nullCpu) {
    size_t dim = 2 * vectorSize() + 2;
    int messageTag = originCpu().value(); // TODO messageTag
    communicator()->recFrom(messageTag, stateBuffer().array(), dim); // Blocking MPI_RECV
    target()->statusIs(Seed::Status(static_cast<int>(stateBuffer()[dim - 2])));
    target()->stateIs(DynamState(vectorSize(), stateBuffer().array()));
    target()->iterationIs(IterationRank(static_cast<int>(stateBuffer()[dim - 1])));
  }

  setStatus(READY);
}

/* MpiReducedSeedReader implementation */

MpiReducedSeedReader::MpiReducedSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer) :
  Reader<Vector>(targetCpu),
  communicator_(comm),
  stateBuffer_(buffer)
{}

CpuRank
MpiReducedSeedReader::localCpu() const {
  return CpuRank(communicator_->myID());
}

void
MpiReducedSeedReader::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiReducedSeedReader::statusIs -- Not ready");

  setStatus(s);

  if (status() != EXECUTING)
    return;

  if (targetCpu() != localCpu() && targetCpu() != nullCpu) {
    size_t dim = origin()->state().size() + 2;
    const double * stateBegin = origin()->state().data();
    std::copy(stateBegin, stateBegin + origin()->state().size(), stateBuffer().array());

    stateBuffer()[dim - 2] = static_cast<double>(origin()->status()); 
    stateBuffer()[dim - 1] = static_cast<double>(origin()->iteration().value());

    int messageTag = localCpu().value(); // TODO more meaningful messageTag
    communicator()->sendTo(targetCpu().value(), messageTag, stateBuffer().array(), dim);
    communicator()->waitForAllReq(); // TODO delayed waitForAllReq
  }

  setStatus(READY);
}

/* MpiReducedSeedWriter implementation */

MpiReducedSeedWriter::MpiReducedSeedWriter(CpuRank originCpu, Communicator * comm, size_t rSize, SimpleBuffer<double> * buffer) :
  Writer<Vector>(originCpu),
  communicator_(comm),
  reducedStateSize_(rSize),
  stateBuffer_(buffer)
{}

CpuRank
MpiReducedSeedWriter::localCpu() const {
  return CpuRank(communicator_->myID());
}

void
MpiReducedSeedWriter::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiReducedSeedWriter::statusIs -- Not ready");

  setStatus(s);
  
  if (status() != EXECUTING)
    return;

  if (originCpu() != localCpu() && originCpu() != nullCpu) {
    size_t dim = reducedStateSize() + 2;
    int messageTag = originCpu().value(); // TODO more meaningful messageTag
    communicator()->recFrom(messageTag, stateBuffer().array(), dim); // Blocking MPI_RECV
    target()->statusIs(Seed::Status(static_cast<int>(stateBuffer()[dim - 2])));
    target()->stateIs(Vector(reducedStateSize(), stateBuffer().array()));
    target()->iterationIs(IterationRank(static_cast<int>(stateBuffer()[dim - 1])));
  }

  setStatus(READY);
}

/* MpiManager implementation */

MpiManager::MpiManager(Communicator * comm, size_t vectorSize) :
  communicator_(comm),
  vectorSize_(vectorSize),
  reducedStateSize_(0),
  sharedBuffer_(computeBufferSize(vectorSize)),
  seedReaderMgr_(SeedReaderFactory(this)),
  reducedSeedReaderMgr_(ReducedSeedReaderFactory(this)),
  seedWriterMgr_(SeedWriterFactory(this)),
  reducedSeedWriterMgr_(ReducedSeedWriterFactory(this))
{}

CpuRank
MpiManager::localCpu() const {
  return CpuRank(communicator_->myID());
}

MpiSeedReader *
MpiManager::reader(const Seed * origin, CpuRank targetCpu) const {
  return seedReaderMgr_.instance(ActivityId<const Seed>(origin, targetCpu));
}

MpiSeedWriter *
MpiManager::writer(Seed * target, CpuRank originCpu) const {
  return seedWriterMgr_.instance(ActivityId<Seed>(target, originCpu));
}

MpiReducedSeedReader *
MpiManager::reader(const ReducedSeed * origin, CpuRank targetCpu) const {
  return reducedSeedReaderMgr_.instance(ActivityId<const ReducedSeed>(origin, targetCpu));
}

MpiReducedSeedWriter *
MpiManager::writer(ReducedSeed * target, CpuRank originCpu) const {
  return reducedSeedWriterMgr_.instance(ActivityId<ReducedSeed>(target, originCpu));
}

MpiSeedReader *
MpiManager::readerNew(const Seed * origin, CpuRank targetCpu) {
  return seedReaderMgr_.instanceNew(ActivityId<const Seed>(origin, targetCpu));
}

MpiSeedWriter *
MpiManager::writerNew(Seed * target, CpuRank originCpu) {
  return seedWriterMgr_.instanceNew(ActivityId<Seed>(target, originCpu));
}

MpiReducedSeedReader *
MpiManager::readerNew(const ReducedSeed * origin, CpuRank targetCpu) {
  return reducedSeedReaderMgr_.instanceNew(ActivityId<const ReducedSeed>(origin, targetCpu));
}

MpiReducedSeedWriter *
MpiManager::writerNew(ReducedSeed * target, CpuRank originCpu) {
  return reducedSeedWriterMgr_.instanceNew(ActivityId<ReducedSeed>(target, originCpu));
}

void
MpiManager::readerDel(const Seed * origin, CpuRank targetCpu) {
  seedReaderMgr_.instanceDel(ActivityId<const Seed>(origin, targetCpu));
}

void
MpiManager::writerDel(Seed * target, CpuRank originCpu) {
  seedWriterMgr_.instanceDel(ActivityId<Seed>(target, originCpu));
}

void
MpiManager::readerDel(const ReducedSeed * origin, CpuRank targetCpu) {
  reducedSeedReaderMgr_.instanceDel(ActivityId<const ReducedSeed>(origin, targetCpu));
}

void
MpiManager::writerDel(ReducedSeed * target, CpuRank originCpu) {
  reducedSeedWriterMgr_.instanceDel(ActivityId<ReducedSeed>(target, originCpu));
}

void
MpiManager::vectorSizeIs(size_t vSize) {
  sharedBuffer_.sizeIs(computeBufferSize(vSize));

  SeedWriterMgr::Iterator it_end = seedWriterMgr_.instanceEnd();
  for (SeedWriterMgr::Iterator it = seedWriterMgr_.instanceBegin(); it != it_end; ++it) {
    it->second->vectorSizeIs(vSize);
  }
  
  vectorSize_ = vSize;
}

void
MpiManager::reducedStateSizeIs(size_t rSize) {
  
  ReducedSeedWriterMgr::Iterator it_end = reducedSeedWriterMgr_.instanceEnd();
  for (ReducedSeedWriterMgr::Iterator it = reducedSeedWriterMgr_.instanceBegin(); it != it_end; ++it) {
    it->second->reducedStateSizeIs(rSize);
  }
  
  reducedStateSize_ = rSize;
}

} /* end namespace RemoteState */ } /* end namespace Pita */
