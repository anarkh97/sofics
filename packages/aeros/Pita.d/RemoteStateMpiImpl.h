#ifndef PITA_REMOTESTATEMPIIMPL_H
#define PITA_REMOTESTATEMPIIMPL_H

#include "RemoteState.h"
#include "Seed.h"

#include "SimpleBuffer.h"
#include <map>

class Communicator;

namespace Pita { namespace RemoteState {

class MpiManager;

/* Reader/Writer for Seed */

class MpiSeedReader : public Reader<DynamState> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiSeedReader);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }
  CpuRank localCpu() const;

  MpiSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer);

protected:
  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  SimpleBuffer<double> * stateBuffer_;
};

class MpiSeedWriter : public Writer<DynamState> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiSeedWriter);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }
  CpuRank localCpu() const;

  size_t vectorSize() const { return vectorSize_; }
  void vectorSizeIs(size_t vSize) { vectorSize_ = vSize; }

  MpiSeedWriter(CpuRank originCpu, Communicator * comm, size_t vSize, SimpleBuffer<double> * buffer);

protected:
  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  size_t vectorSize_;
  SimpleBuffer<double> * stateBuffer_;
};

/* Reader/Writer for ReducedSeed */

class MpiReducedSeedReader : public Reader<Vector> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiReducedSeedReader);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }
  CpuRank localCpu() const;

  MpiReducedSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer);
  
protected:
  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  SimpleBuffer<double> * stateBuffer_;
};

class MpiReducedSeedWriter : public Writer<Vector> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiReducedSeedWriter);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }
  CpuRank localCpu() const;

  size_t reducedStateSize() const { return reducedStateSize_; }
  void reducedStateSizeIs(size_t rSize) { reducedStateSize_ = rSize; }

  MpiReducedSeedWriter(CpuRank originCpu, Communicator * comm, size_t rSize, SimpleBuffer<double> * buffer);

protected:
  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  size_t reducedStateSize_;
  SimpleBuffer<double> * stateBuffer_;
};

/* MpiManager */

class MpiManager : public Manager {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiManager);

  /* Overriden */
  virtual CpuRank localCpu() const;
  
  virtual MpiSeedReader * reader(const Seed * origin, CpuRank targetCpu) const;
  virtual MpiSeedWriter * writer(Seed * target, CpuRank originCpu) const;

  virtual MpiSeedReader * readerNew(const Seed * origin, CpuRank targetCpu);
  virtual MpiSeedWriter * writerNew(Seed * target, CpuRank originCpu);
  
  virtual void readerDel(const Seed * origin, CpuRank targetCpu);
  virtual void writerDel(Seed * target, CpuRank originCpu);
  
  virtual MpiReducedSeedReader * reader(const ReducedSeed * origin, CpuRank targetCpu) const;
  virtual MpiReducedSeedWriter * writer(ReducedSeed * target, CpuRank originCpu) const;
  
  virtual MpiReducedSeedReader * readerNew(const ReducedSeed * origin, CpuRank targetCpu);
  virtual MpiReducedSeedWriter * writerNew(ReducedSeed * target, CpuRank originCpu);

  virtual void readerDel(const ReducedSeed * origin, CpuRank targetCpu);
  virtual void writerDel(ReducedSeed * target, CpuRank originCpu);
  
  /* Added */
  Communicator * communicator() const { return communicator_; }

  size_t vectorSize() const { return vectorSize_; }
  size_t reducedStateSize() const { return reducedStateSize_; } // Assumed to be smaller than 2 * vectorSize

  void vectorSizeIs(size_t vSize);
  void reducedStateSizeIs(size_t rSize);

  SimpleBuffer<double> * sharedBuffer() { return &sharedBuffer_; }
  
  static Ptr New(Communicator * comm, size_t vectorSize) {
    return new MpiManager(comm, vectorSize);
  }

protected:
  MpiManager(Communicator * comm, size_t vectorSize);

  static int computeBufferSize(size_t vectorSize) { return 2 * vectorSize + 2; }

private:
  Communicator * communicator_;
  
  size_t vectorSize_, reducedStateSize_;
  SimpleBuffer<double> sharedBuffer_;


  /* Manager implementations (4 x) */
  template <typename T>
  class ActivityId {
  public:
    ActivityId(T * state, CpuRank peerCpu) :
      state_(state), peerCpu_(peerCpu)
    {}

    T * state() const { return state_; }
    CpuRank peerCpu() const { return peerCpu_; }

    bool operator<(const ActivityId<T> & other) const {
      return state() != other.state() ? state() < other.state() : peerCpu() < other.peerCpu();
    }

    bool operator==(const ActivityId<T> & other) const {
      return state() == other.state() && peerCpu() == other.peerCpu();
    }

  private:
    T * state_;
    CpuRank peerCpu_;
  };


  class SeedReaderFactory : public Fwk::InstanceFactory<MpiSeedReader, ActivityId<const Seed> > {
  public:
    explicit SeedReaderFactory(MpiManager * parent) : parent_(parent) {}

    MpiSeedReader * operator()(const ActivityId<const Seed> & key) const {
      MpiSeedReader * result = new MpiSeedReader(key.peerCpu(), parent_->communicator(), parent_->sharedBuffer());
      result->originIs(key.state());
      return result;
    }

  private:
    MpiManager * parent_;
  };

  typedef Fwk::FactoryManagerImpl<SeedReaderFactory> SeedReaderMgr;
  SeedReaderMgr seedReaderMgr_;


  class ReducedSeedReaderFactory : public Fwk::InstanceFactory<MpiReducedSeedReader, ActivityId<const ReducedSeed> > {
  public:
    explicit ReducedSeedReaderFactory(MpiManager * parent) : parent_(parent) {}

    MpiReducedSeedReader * operator()(const ActivityId<const ReducedSeed> & key) const {
      MpiReducedSeedReader * result = new MpiReducedSeedReader(key.peerCpu(), parent_->communicator(), parent_->sharedBuffer());
      result->originIs(key.state());
      return result;
    }

  private:
    MpiManager * parent_;
  };

  typedef Fwk::FactoryManagerImpl<ReducedSeedReaderFactory> ReducedSeedReaderMgr;
  ReducedSeedReaderMgr reducedSeedReaderMgr_;
  
  
  class SeedWriterFactory : public Fwk::InstanceFactory<MpiSeedWriter, ActivityId<Seed> > {
  public:
    explicit SeedWriterFactory(MpiManager * parent) : parent_(parent) {}

    MpiSeedWriter * operator()(const ActivityId<Seed> & key) const {
      MpiSeedWriter * result = new MpiSeedWriter(key.peerCpu(), parent_->communicator(), parent_->vectorSize(), parent_->sharedBuffer());
      result->targetIs(key.state());
      return result;
    }

  private:
    MpiManager * parent_;
  };

  typedef Fwk::FactoryManagerImpl<SeedWriterFactory> SeedWriterMgr;
  SeedWriterMgr seedWriterMgr_;


  class ReducedSeedWriterFactory : public Fwk::InstanceFactory<MpiReducedSeedWriter, ActivityId<ReducedSeed> > {
  public:
    explicit ReducedSeedWriterFactory(MpiManager * parent) : parent_(parent) {}

    MpiReducedSeedWriter * operator()(const ActivityId<ReducedSeed> & key) const {
      MpiReducedSeedWriter * result = new MpiReducedSeedWriter(key.peerCpu(), parent_->communicator(), parent_->reducedStateSize(), parent_->sharedBuffer());
      result->targetIs(key.state());
      return result;
    }

  private:
    MpiManager * parent_;
  };

  typedef Fwk::FactoryManagerImpl<ReducedSeedWriterFactory> ReducedSeedWriterMgr;
  ReducedSeedWriterMgr reducedSeedWriterMgr_;
};

} /* end namespace RemoteState */ } /* end namespace Pita */

#endif /* PITA_REMOTESTATEMPIIMPL_H */
