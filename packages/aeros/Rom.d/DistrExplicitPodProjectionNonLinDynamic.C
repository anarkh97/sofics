#include "DistrExplicitPodProjectionNonLinDynamic.h"

#include "DistrGalerkinProjectionSolver.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrVecBasis.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Feti.d/DistrVector.h>

#include <Driver.d/GeoSource.h>

#include <cstddef>

extern Communicator *structCom;
extern GeoSource *geoSource;

namespace Rom {

class DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler {
public:
  virtual void currentTimeIs(double t) = 0;
  virtual void stateSnapshotAdd(const DistrVector &state) = 0;
  virtual void accelerationSnapshotAdd(const DistrVector &accel) = 0;
  virtual void velocitySnapshotAdd(const DistrVector &veloc) = 0;
  virtual void forceSnapshotAdd(const DistrVector &) = 0;
 
  SnapshotHandler() {} 
  virtual ~SnapshotHandler();

private:
  // Disallow copy and assignment
  SnapshotHandler(const SnapshotHandler &);
  SnapshotHandler &operator=(const SnapshotHandler &);
};

DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler::~SnapshotHandler() {
  // Nothing to do
}


// Dummy class, used for namespace access
class DistrExplicitPodProjectionNonLinDynamicDetail : public DistrExplicitPodProjectionNonLinDynamic {
public:
  class NoOpSnapshotHandler;
  class RecordingSnapshotHandler;

private:
  // Dummy constructor
  DistrExplicitPodProjectionNonLinDynamicDetail();
};

// NoOp implementation
class DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler : public DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler {
public:
  void currentTimeIs(double t);
  void stateSnapshotAdd(const DistrVector &s);
  void accelerationSnapshotAdd(const DistrVector &accel);
  void velocitySnapshotAdd(const DistrVector &veloc);
  void forceSnapshotAdd(const DistrVector &); //overriden
};

// Recording implementation
class DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler : public DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler {
public:
  void currentTimeIs(double t);
  void stateSnapshotAdd(const DistrVector &s);
  void accelerationSnapshotAdd(const DistrVector &accel);
  void velocitySnapshotAdd(const DistrVector &veloc);
  void forceSnapshotAdd(const DistrVector &s); // overriden
  explicit RecordingSnapshotHandler(DecDomain *decDom);

private:
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  
  DistrExplicitPodProjectionNonLinDynamic *parent_;  
  DecDomain *decDomain_;

  DistrVecNodeDof6Conversion converter_;
  DistrMasterMapping masterMapping_;
  DistrVector assembledSnapshot_;
  DistrNodeDof6Buffer buffer_;
  DistrBasisOutputFile * stateOutputFile_;
  DistrBasisOutputFile * accelOutputFile_;
  DistrBasisOutputFile * velocOutputFile_;
  DistrBasisOutputFile * forceOutputFile_;

  int stateSkip_;
  int accelSkip_;
  int velocSkip_;
  int forceSkip_;  

  double currentTime_;
};


// Main class implementation

DistrExplicitPodProjectionNonLinDynamic::DistrExplicitPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitPodProjectionNonLinDynamicBase(domain),
  snapshotHandler_(nullptr)
{}

DistrExplicitPodProjectionNonLinDynamic::~DistrExplicitPodProjectionNonLinDynamic() {
  // Nothing to do
}

void
DistrExplicitPodProjectionNonLinDynamic::preProcess() {
  DistrExplicitPodProjectionNonLinDynamicBase::preProcess();

  collectState = false;
  collectAccel = false;
  collectVeloc = false;
  collectForce = false;

  if(domain->solInfo().statevectPodRom)
    collectState = true;
  if(domain->solInfo().accelvectPodRom)
    collectAccel = true;
  if(domain->solInfo().velocvectPodRom)
    collectVeloc = true;
  if(domain->solInfo().forcevectPodRom)
    collectForce = true;

  if (domain->solInfo().snapshotsPodRom) {
    snapshotHandler_.reset(new DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler(this->decDomain));
  } else {
    snapshotHandler_.reset(new DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler);
  }
}

void
DistrExplicitPodProjectionNonLinDynamic::currentTimeIs(double t) {
  snapshotHandler_->currentTimeIs(t);
}

void
DistrExplicitPodProjectionNonLinDynamic::stateSnapshotAdd(const DistrVector &f) {
  snapshotHandler_->stateSnapshotAdd(f);
}

void
DistrExplicitPodProjectionNonLinDynamic::accelerationSnapshotAdd(const DistrVector &accel) {
  snapshotHandler_->accelerationSnapshotAdd(accel);
}

void
DistrExplicitPodProjectionNonLinDynamic::velocitySnapshotAdd(const DistrVector &veloc) {
  snapshotHandler_->velocitySnapshotAdd(veloc);
}

void
DistrExplicitPodProjectionNonLinDynamic::forceSnapshotAdd(const DistrVector &f) {
  snapshotHandler_->forceSnapshotAdd(f);
}

// Implementation of auxiliary classes

void
DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler::currentTimeIs(double time) {
  // Nothing to do
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler::stateSnapshotAdd(const DistrVector &state) {
  // Nothing to do
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler::accelerationSnapshotAdd(const DistrVector &accel) {
  // Nothing to do
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler::velocitySnapshotAdd(const DistrVector &veloc) {
  // Nothing to do
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::NoOpSnapshotHandler::forceSnapshotAdd(const DistrVector &) {
  // Nothing to do
}

DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler::RecordingSnapshotHandler(DecDomain *decDom) :
  decDomain_(decDom),
  converter_(decDom->getAllSubDomains(), decDom->getAllSubDomains() + decDom->getNumSub()),
  masterMapping_(SubDomIt(decDom->getAllSubDomains()), SubDomIt(decDom->getAllSubDomains() + decDom->getNumSub())),
  buffer_(masterMapping_.localNodeBegin(), masterMapping_.localNodeEnd()),
  assembledSnapshot_(decDom->solVecInfo()),
  stateSkip_(0),
  accelSkip_(0),
  velocSkip_(0),
  forceSkip_(0)
{
  if(decDom->getDomain()->solInfo().statevectPodRom){
    stateOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::STATE, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
              buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(decDom->getDomain()->solInfo().accelvectPodRom){
    accelOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::ACCELERATION, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
              buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(decDom->getDomain()->solInfo().velocvectPodRom){
    velocOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::VELOCITY, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
              buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(decDom->getDomain()->solInfo().forcevectPodRom){
    forceOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::FORCE, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
              buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler::currentTimeIs(double time) {
  currentTime_ = time;
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler::stateSnapshotAdd(const DistrVector &state) {
  ++stateSkip_;
  if (stateSkip_ >= decDomain_->getDomain()->solInfo().skipState) {
    converter_.paddedNodeDof6(state, buffer_);
    stateOutputFile_->stateAdd(buffer_, currentTime_);
    stateSkip_ = 0;
  }
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler::accelerationSnapshotAdd(const DistrVector &accel) {
  ++accelSkip_;
  if (accelSkip_ >= decDomain_->getDomain()->solInfo().skipAccel) {
    converter_.paddedNodeDof6(accel, buffer_);
    accelOutputFile_->stateAdd(buffer_, currentTime_);
    accelSkip_ = 0;
  }
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler::velocitySnapshotAdd(const DistrVector &veloc) {
  ++velocSkip_;
  if (velocSkip_ >= decDomain_->getDomain()->solInfo().skipVeloc) {
    converter_.paddedNodeDof6(veloc, buffer_);
    velocOutputFile_->stateAdd(buffer_, currentTime_);
    velocSkip_ = 0;
  }
}

void
DistrExplicitPodProjectionNonLinDynamicDetail::RecordingSnapshotHandler::forceSnapshotAdd(const DistrVector &f) {
  ++forceSkip_;
  if (forceSkip_ >= decDomain_->getDomain()->solInfo().skipForce) {
    assembledSnapshot_ = f;
    decDomain_->getSolVecAssembler()->assemble(assembledSnapshot_);
    converter_.paddedNodeDof6(assembledSnapshot_, buffer_);
    forceOutputFile_->stateAdd(buffer_, currentTime_);
    forceSkip_ = 0;
  }
}

} // end namespace Rom
