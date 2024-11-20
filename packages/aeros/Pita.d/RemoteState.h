#ifndef PITA_REMOTESTATE_H
#define PITA_REMOTESTATE_H

#include "Fwk.h"
#include "Types.h"
#include "SharedState.h"

#include "DynamState.h"

namespace Pita { namespace RemoteState {

const CpuRank nullCpu(-1);

enum Status {
  INACTIVE = 0,
  READY,
  EXECUTING
};

class Activity : public Fwk::PtrInterface<Activity> {
public:
  EXPORT_PTRINTERFACE_TYPES(Activity);

  Status status() const { return status_; }
  virtual void statusIs(Status s) = 0;

protected:
  Activity() : status_(INACTIVE) {}

  void setStatus(Status s) { status_ = s; }

  virtual bool isReady() const = 0;

private:
  Status status_;

  DISALLOW_COPY_AND_ASSIGN(Activity);
};

class BaseReader : public Activity {
public:
  EXPORT_PTRINTERFACE_TYPES(BaseReader);

  CpuRank targetCpu() const { return targetCpu_; }

protected:
  explicit BaseReader(CpuRank targetCpu) :
    targetCpu_(targetCpu)
  {}

  void setTargetCpu(CpuRank tc) { targetCpu_ = tc; }

private:
  CpuRank targetCpu_;
};

template <typename S>
class Reader : public BaseReader {
public:
  EXPORT_PTRINTERFACE_TYPES(Reader);

  const SharedState<S> * origin() const { return origin_.ptr(); }
  virtual void originIs(const SharedState<S> * o);

protected:
  explicit Reader(CpuRank targetCpu) :
    BaseReader(targetCpu)
  {}

  virtual bool isReady() const { return origin_; }
  
  void setOrigin(const SharedState<S> * o) { origin_ = o; }

private:
  Fwk::Ptr<const SharedState<S> > origin_;
};

class BaseWriter : public Activity {
public:
  EXPORT_PTRINTERFACE_TYPES(BaseWriter);

  CpuRank originCpu() const { return originCpu_; }

protected:
  explicit BaseWriter(CpuRank originCpu) :
    originCpu_(originCpu)
  {}

  void setOriginCpu(CpuRank oc) { originCpu_ = oc; }

private:
  CpuRank originCpu_;
};

template <typename S>
class Writer : public BaseWriter {
public:
  EXPORT_PTRINTERFACE_TYPES(Writer);

  SharedState<S> * target() const { return target_.ptr(); }
  virtual void targetIs(SharedState<S> * t);

protected:
  explicit Writer(CpuRank originCpu) :
    BaseWriter(originCpu)
  {}

  virtual bool isReady() const { return target_; }

  void setTarget(SharedState<S> * t) { target_ = t; }

private:
  Fwk::Ptr<SharedState<S> > target_;
};


class Manager : public Fwk::PtrInterface<Manager> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  virtual CpuRank localCpu() const = 0;

  // Seed
  virtual Reader<DynamState> * reader(const SharedState<DynamState> * origin, CpuRank targetCpu) const = 0;
  virtual Writer<DynamState> * writer(SharedState<DynamState> * target, CpuRank originCpu) const = 0;

  virtual Reader<DynamState> * readerNew(const SharedState<DynamState> * origin, CpuRank targetCpu) = 0;
  virtual Writer<DynamState> * writerNew(SharedState<DynamState> * target, CpuRank originCpu) = 0;
  
  virtual void readerDel(const SharedState<DynamState> * origin, CpuRank targetCpu) = 0;
  virtual void writerDel(SharedState<DynamState> * target, CpuRank originCpu) = 0;

  // ReducedSeed
  virtual Reader<Vector> * reader(const SharedState<Vector> * origin, CpuRank targetCpu) const = 0;
  virtual Writer<Vector> * writer(SharedState<Vector> * target, CpuRank originCpu) const = 0;
  
  virtual Reader<Vector> * readerNew(const SharedState<Vector> * origin, CpuRank targetCpu) = 0;
  virtual Writer<Vector> * writerNew(SharedState<Vector> * target, CpuRank originCpu) = 0;

  virtual void readerDel(const SharedState<Vector> * origin, CpuRank targetCpu) = 0;
  virtual void writerDel(SharedState<Vector> * target, CpuRank originCpu) = 0;

protected:
  Manager() {}

private:
  DISALLOW_COPY_AND_ASSIGN(Manager);
};


template <typename S>
void
Reader<S>::originIs(const SharedState<S> * o) {
  setOrigin(o);
  setStatus(isReady() ? READY : INACTIVE);
}

template <typename S>
void
Writer<S>::targetIs(SharedState<S> * t) {
  setTarget(t);
  setStatus(isReady() ? READY : INACTIVE);
}

typedef Reader<DynamState> SeedReader;
typedef Reader<Vector> ReducedSeedReader;
typedef Writer<DynamState> SeedWriter;
typedef Writer<Vector> ReducedSeedWriter;

} /* end namespace RemoteState */ } /* end namespace Pita */

#endif /* PITA_REMOTESTATE_H */
