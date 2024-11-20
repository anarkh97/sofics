#include "DistrExplicitSnapshotNonLinDynamic.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"
#include "BasisOutputFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrDomainUtils.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include <Driver.d/GeoSource.h>

#include <algorithm>
#include <cstddef>

extern Communicator *structCom;
extern GeoSource *geoSource;

namespace Rom {

class DistrExplicitSnapshotNonLinDynamic::SnapshotHandler {
public:
  void currentTimeIs(double t);
  void snapshotAdd(const DistrVector &s);
  void accelerationSnapshotAdd(const DistrVector &accel);
  void velocitySnapshotAdd(const DistrVector &veloc);
  void forceSnapshotAdd(const DistrVector &f);
  explicit SnapshotHandler(DistrExplicitSnapshotNonLinDynamic *parent);
  


private:
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  
  DistrExplicitSnapshotNonLinDynamic *parent_;

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

DistrExplicitSnapshotNonLinDynamic::DistrExplicitSnapshotNonLinDynamic(Domain *domain) :
  MultiDomainDynam(domain),
  domain_(domain),
  snapshotHandler_(NULL)
{}

DistrExplicitSnapshotNonLinDynamic::~DistrExplicitSnapshotNonLinDynamic() {
  delete snapshotHandler_;
}

void
DistrExplicitSnapshotNonLinDynamic::preProcess() {
  MultiDomainDynam::preProcess();
  snapshotHandler_ = new SnapshotHandler(this);

  collectState = false;
  collectAccel = false;
  collectVeloc = false;
  collectForce = false;

  if(domain_->solInfo().statevectPodRom)
    collectState = true;
  if(domain_->solInfo().accelvectPodRom)
    collectAccel = true;
  if(domain_->solInfo().velocvectPodRom)
    collectVeloc = true;
  if(domain_->solInfo().forcevectPodRom)
    collectForce = true;

}

void
DistrExplicitSnapshotNonLinDynamic::currentTimeIs(double t) {
  snapshotHandler_->currentTimeIs(t);
}

void
DistrExplicitSnapshotNonLinDynamic::stateSnapshotAdd(const DistrVector &state) {
  snapshotHandler_->snapshotAdd(state);
}

void
DistrExplicitSnapshotNonLinDynamic::accelerationSnapshotAdd(const DistrVector &accel) {
  if(geomState->getHaveRot()) {
    DistrVector a(accel);
    geomState->transform(a, 6, true); // transform convected angular acceleration to second time derivative of total rotation vector
    snapshotHandler_->accelerationSnapshotAdd(a);
  }
  else {
    snapshotHandler_->accelerationSnapshotAdd(accel);
  }
}

void
DistrExplicitSnapshotNonLinDynamic::velocitySnapshotAdd(const DistrVector &veloc) {
  if(geomState->getHaveRot()) {
    DistrVector v(veloc);
    geomState->transform(v, 2, true); // transform convected angular velocity to time derivative of total rotation vector
    snapshotHandler_->velocitySnapshotAdd(v);
  }
  else {
    snapshotHandler_->velocitySnapshotAdd(veloc);
  }
}

void
DistrExplicitSnapshotNonLinDynamic::forceSnapshotAdd(const DistrVector &f) {
  snapshotHandler_->forceSnapshotAdd(f);
}

DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::SnapshotHandler(DistrExplicitSnapshotNonLinDynamic *parent) :
  parent_(parent),
  converter_(parent->decDomain->getAllSubDomains(), parent->decDomain->getAllSubDomains() + parent->decDomain->getNumSub()),
  masterMapping_(SubDomIt(parent->decDomain->getAllSubDomains()), SubDomIt(parent->decDomain->getAllSubDomains() + parent->decDomain->getNumSub())),
  buffer_(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd()),
  assembledSnapshot_(parent->decDomain->solVecInfo()),
  stateSkip_(0),
  accelSkip_(0),
  velocSkip_(0),
  forceSkip_(0),
  currentTime_(0.0)
{
  if(parent_->domain->solInfo().statevectPodRom){
    stateOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::STATE, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
             buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(parent_->domain->solInfo().accelvectPodRom){
    accelOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::ACCELERATION, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
             buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(parent_->domain->solInfo().velocvectPodRom){
    velocOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::VELOCITY, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
             buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(parent_->domain->solInfo().forcevectPodRom){
    forceOutputFile_ = new DistrBasisOutputFile(BasisFileId(FileNameInfo(), BasisId::FORCE, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
             buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom, (geoSource->getCheckFileInfo()->lastRestartFile != 0));}

}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::currentTimeIs(double time) {
  currentTime_ = time;
}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::snapshotAdd(const DistrVector &state) {
  ++stateSkip_;
  if (stateSkip_ >= parent_->domain->solInfo().skipState) {
    converter_.paddedNodeDof6(state, buffer_);
    stateOutputFile_->stateAdd(buffer_, currentTime_);
    stateSkip_ = 0;
  }
}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::accelerationSnapshotAdd(const DistrVector &accel) {
  ++accelSkip_;
  if (accelSkip_ >= parent_->domain->solInfo().skipAccel) {
    converter_.paddedNodeDof6(accel, buffer_);
    accelOutputFile_->stateAdd(buffer_, currentTime_);
    accelSkip_ = 0;
  }
}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::velocitySnapshotAdd(const DistrVector &veloc) {
  ++velocSkip_;
  if (velocSkip_ >= parent_->domain->solInfo().skipVeloc) {
    converter_.paddedNodeDof6(veloc, buffer_);
    velocOutputFile_->stateAdd(buffer_, currentTime_);
    velocSkip_ = 0;
  }
}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::forceSnapshotAdd(const DistrVector &f) {
  ++forceSkip_;
  if (forceSkip_ >= parent_->domain->solInfo().skipForce) {
    assembledSnapshot_ = f;
    parent_->decDomain->getSolVecAssembler()->assemble(assembledSnapshot_);
    converter_.paddedNodeDof6(assembledSnapshot_, buffer_);
    forceOutputFile_->stateAdd(buffer_, currentTime_);
    forceSkip_ = 0;
  }
}

} // end namespace Rom
