#include "DistrSnapshotNonLinDynamic.h"

#include "DistrBasisFile.h"
#include "FileNameInfo.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/Domain.h>
#include <Driver.d/DecDomain.h>
#include <Driver.d/GeoSource.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

#include <cstddef>
#include <memory>

extern Communicator *structCom;
extern GeoSource *geoSource;

namespace Rom {

// Dummy class holding the implementation of DistrSnapshotNonLinDynamic
struct DistrSnapshotNonLinDynamicDetail : private DistrSnapshotNonLinDynamic {
  class RawImpl : public Impl {
  public:
    // Overriden
    virtual void lastMidTimeIs(double t);
    virtual void lastDeltaIs(double delta);
    virtual void stateSnapshotAdd(const DistrGeomState &);
    virtual void velocSnapshotAdd(const DistrVector &);
    virtual void accelSnapshotAdd(const DistrVector &);
    virtual void dsvarSnapshotAdd(const DistrGeomState &);
    virtual void muvarSnapshotAdd(const DistrGeomState &);
    virtual void postProcess();

    int nodeCount() const { return geoSource->getNumGlobNodes(); }

    explicit RawImpl(DecDomain *);

  private:
    typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
    
    DecDomain * decDomain_;
    DistrVecNodeDof6Conversion converter_;
    DistrMasterMapping masterMapping_;
    DistrMpcMasterMapping *mpcMasterMapping_;
    DistrMasterMapping *dualMasterMapping_;

    DistrNodeDof6Buffer snapBuffer_;
    DistrNodeDof1Buffer *dualSnapBuffer_;
    DistrNodeDof1Buffer *muSnapBuffer_;
    
    FileNameInfo fileInfo_;
    DistrBasisOutputFile *stateSnapFile_;
    DistrBasisOutputFile *velocSnapFile_;
    DistrBasisOutputFile *accelSnapFile_;
    DistrBasisOutputFile *dsvarSnapFile_;
    DistrBasisOutputFile *muvarSnapFile_;
    
    double timeStamp_;
    int stateSkip_;
    int velocSkip_;
    int accelSkip_;
    int dsvarSkip_;
    int muvarSkip_;
    int mpcOffset_;
  };

private:
  // Dummy constructor to avoid compilation failures
  DistrSnapshotNonLinDynamicDetail() :
    DistrSnapshotNonLinDynamic(NULL)
  {}
};

DistrSnapshotNonLinDynamicDetail::RawImpl::RawImpl(DecDomain *decDomain) :
  decDomain_(decDomain),
  converter_(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub()),
  masterMapping_(SubDomIt(decDomain->getAllSubDomains()), SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub())),
  mpcMasterMapping_(NULL),
  dualMasterMapping_(NULL),
  snapBuffer_(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd()),
  dualSnapBuffer_(NULL),
  fileInfo_(),
  stateSnapFile_(NULL),
  velocSnapFile_(NULL),
  accelSnapFile_(NULL),
  dsvarSnapFile_(NULL),
  muvarSnapFile_(NULL),
  timeStamp_(decDomain->getDomain()->solInfo().initialTime),
  stateSkip_(0),
  velocSkip_(0),
  accelSkip_(0),
  dsvarSkip_(0),
  muvarSkip_(0)
{
  if(decDomain->getDomain()->solInfo().statevectPodRom){

    stateSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::STATE, BasisId::SNAPSHOTS), nodeCount(),
                     snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom,
                     (geoSource->getCheckFileInfo()->lastRestartFile != 0));

  }
  if(decDomain->getDomain()->solInfo().velocvectPodRom){ 

    velocSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::VELOCITY, BasisId::SNAPSHOTS), nodeCount(),
                     snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom,
                     (geoSource->getCheckFileInfo()->lastRestartFile != 0));

  }
  if(decDomain->getDomain()->solInfo().accelvectPodRom){ 

    accelSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::ACCELERATION, BasisId::SNAPSHOTS), nodeCount(),
                     snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom,
                     (geoSource->getCheckFileInfo()->lastRestartFile != 0));

  }
  if(decDomain->getDomain()->solInfo().dsvPodRom){
    if(isFeti(decDomain->getDomain()->solInfo().solvercntl->type)) {

      mpcMasterMapping_ = new DistrMpcMasterMapping(SubDomIt(decDomain->getAllSubDomains()), SubDomIt(decDomain->getAllSubDomains()
                                                    + decDomain->getNumSub()));
      dualSnapBuffer_   = new DistrNodeDof1Buffer(mpcMasterMapping_->masterNodeBegin(), mpcMasterMapping_->masterNodeEnd());

    } else {

      dualMasterMapping_ = new DistrMasterMapping(SubDomIt(decDomain->getAllSubDomains()), SubDomIt(decDomain->getAllSubDomains()
                                                  + decDomain->getNumSub()), true);
      mpcOffset_         = decDomain->getDomain()->numNode()-decDomain->getDomain()->getNumCTC();
      dualSnapBuffer_    = new DistrNodeDof1Buffer(dualMasterMapping_->masterNodeBegin(), dualMasterMapping_->masterNodeEnd(), mpcOffset_);

    }

    dsvarSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::DUALSTATE, BasisId::SNAPSHOTS), decDomain->getDomain()->getNumCTC(),
                     dualSnapBuffer_->globalNodeIndexBegin(), dualSnapBuffer_->globalNodeIndexEnd(), structCom, 
                     (geoSource->getCheckFileInfo()->lastRestartFile != 0), 1);
  }
  if(decDomain->getDomain()->solInfo().muvPodRom){

    muSnapBuffer_  = new DistrNodeDof1Buffer(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd());

    muvarSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::MUSTATE, BasisId::SNAPSHOTS), nodeCount(),
                     muSnapBuffer_->globalNodeIndexBegin(), muSnapBuffer_->globalNodeIndexEnd(), structCom,
                     (geoSource->getCheckFileInfo()->lastRestartFile != 0), 1);
  }
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::postProcess() {
  // Nothing to do
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::lastMidTimeIs(double t) {
  // t = t_n + dt*(1 - alphaf)
  timeStamp_ = t;
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::lastDeltaIs(double delta) {
  // delta = dt/2
  // timeStamp_ = t_{n+1} = t_n + dt = t_n + dt*(1 - alphaf) + dt*alphaf 
  //            = timeStamp_ + 2*delta*alphaf
  timeStamp_ += 2*delta*::domain->solInfo().newmarkAlphaF;
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::stateSnapshotAdd(const DistrGeomState &snap) {
  ++stateSkip_;
  if(stateSnapFile_ && (stateSkip_ >= decDomain_->getDomain()->solInfo().skipState)) {
    const int subDomCount = snap.getNumSub();
    DistrMasterMapping::SubMasterMappingIt mappingIt = masterMapping_.begin();
    for (int iSub = 0; iSub < subDomCount; ++iSub) {
      const GeomState &subSnap = *snap[iSub];
      const CoordSet &refCoords = decDomain_->getSubDomain(iSub)->getNodes();
      const MasterMapping &mapping = *mappingIt++;

      typedef MasterMapping::IndexPairIterator IndexPairIt;
      const IndexPairIt nodeItEnd = mapping.end();
      for (IndexPairIt nodeIt = mapping.begin(); nodeIt != nodeItEnd; ++nodeIt) {
        // Indexing
        const int iLocalNode = nodeIt->local;
        if(!(refCoords[iLocalNode] && subSnap.getNodeFlag(iLocalNode) > 0)) continue;
        const int iGlobalNode = nodeIt->global;

        // Translational dofs
        snapBuffer_[iGlobalNode][0] = subSnap[iLocalNode].x - refCoords[iLocalNode]->x;
        snapBuffer_[iGlobalNode][1] = subSnap[iLocalNode].y - refCoords[iLocalNode]->y;
        snapBuffer_[iGlobalNode][2] = subSnap[iLocalNode].z - refCoords[iLocalNode]->z;

        // Rotational dofs
        // old method: collect the rescaled rotation vector
        //mat_to_vec(const_cast<double (*)[3]>(subSnap[iLocalNode].R), &snapBuffer_[iGlobalNode][3]);

        // new method: collect the unscaled rotation vector which has already be computed and stored in NodeState::theta
        snapBuffer_[iGlobalNode][3] = subSnap[iLocalNode].theta[0];
        snapBuffer_[iGlobalNode][4] = subSnap[iLocalNode].theta[1];
        snapBuffer_[iGlobalNode][5] = subSnap[iLocalNode].theta[2];
      }
    }

    stateSnapFile_->stateAdd(snapBuffer_, timeStamp_);
    stateSkip_ = 0;
  }
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::velocSnapshotAdd(const DistrVector &veloc) {
  ++velocSkip_;
  if (velocSnapFile_ && (velocSkip_ >= decDomain_->getDomain()->solInfo().skipVeloc)) {
    converter_.paddedNodeDof6(veloc, snapBuffer_);
    velocSnapFile_->stateAdd(snapBuffer_, timeStamp_);
    velocSkip_ = 0;
  }
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::accelSnapshotAdd(const DistrVector &accel) {
  ++accelSkip_;
  if (accelSnapFile_ && (accelSkip_ >= decDomain_->getDomain()->solInfo().skipAccel)) {
    converter_.paddedNodeDof6(accel, snapBuffer_);
    accelSnapFile_->stateAdd(snapBuffer_, timeStamp_);
    accelSkip_ = 0;
  }
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::muvarSnapshotAdd(const DistrGeomState &snap) {
  ++muvarSkip_;
  if(muvarSnapFile_ && (muvarSkip_ >= decDomain_->getDomain()->solInfo().skipMuStateVar)) {
    const int subDomCount = snap.getNumSub();
    muSnapBuffer_->zero(); // zero out snapshot buffer
    for (int iSub = 0; iSub < subDomCount; ++iSub) { // loop over each subdomain in this mpi process
      SubDomain *sd = decDomain_->getSubDomain(iSub);
      // loop through each multiplier in this subdomain
      for(std::map<std::pair<int,int>, double>::iterator it = snap.mu[iSub].begin(); it != snap.mu[iSub].end(); it++) { 
        std::pair<int,int> id = it->first; // get lmpc id and slave node pair
        double lagrangeMultiplier = it->second; // get lagrangemultiplier
        int slaveNodeGId = id.second;
        //fprintf(stderr,"slave node: %d, mu: %3.2e\n",slaveNodeGId, lagrangeMultiplier);
        const int pnId = sd->globalToLocal(slaveNodeGId); // see if its in this subdomain
        if(pnId > 0 && lagrangeMultiplier < -1e-16){
          if((*muSnapBuffer_)[slaveNodeGId] != NULL) // see if this process owns that node
            (*muSnapBuffer_)[slaveNodeGId][0] = lagrangeMultiplier; // write to buffer
        }
      }
    }
  }
  muvarSnapFile_->stateAdd(*muSnapBuffer_, timeStamp_);
  muvarSkip_ = 0;
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::dsvarSnapshotAdd(const DistrGeomState &snap) {
  ++dsvarSkip_;
  if(dsvarSnapFile_ && (dsvarSkip_ >= decDomain_->getDomain()->solInfo().skipDualStateVar)) {
    const int subDomCount = snap.getNumSub();
    if(isFeti(decDomain_->getDomain()->solInfo().solvercntl->type)) { // FETI-DP with "multipliers" constraint method
    DistrMasterMapping::SubMasterMappingIt mappingIt = mpcMasterMapping_->begin();
    for (int iSub = 0; iSub < subDomCount; ++iSub) {
      const GeomState     &subSnap = *snap[iSub];
      const MasterMapping &mapping = *mappingIt++;

      std::vector<double> lambda;
      decDomain_->getSubDomain(iSub)->getLocalMultipliers(lambda);

      typedef MasterMapping::IndexPairIterator IndexPairIt;
      const IndexPairIt nodeItEnd = mapping.end();
      for (IndexPairIt nodeIt = mapping.begin(); nodeIt != nodeItEnd; ++nodeIt) {
        // Indexing
        const int iLocalNode = nodeIt->local;
        const int iGlobalNode = nodeIt->global;
        
        (*dualSnapBuffer_)[iGlobalNode][0] = lambda[iLocalNode];
      }
    }}
    else { // Multi-domain MUMPS with "augmented" constraint method
    DistrMasterMapping::SubMasterMappingIt mappingIt = dualMasterMapping_->begin();
    for (int iSub = 0; iSub < subDomCount; ++iSub) {
      const GeomState     &subSnap = *snap[iSub];
      const MasterMapping &mapping = *mappingIt++;

      typedef MasterMapping::IndexPairIterator IndexPairIt;
      const IndexPairIt nodeItEnd = mapping.end();
      for (IndexPairIt nodeIt = mapping.begin(); nodeIt != nodeItEnd; ++nodeIt) {
        // Indexing
        const int iLocalNode = nodeIt->local;
        if(subSnap.getNodeFlag(iLocalNode) > 0) continue;
        const int iGlobalNode = nodeIt->global - mpcOffset_;

        (*dualSnapBuffer_)[iGlobalNode][0] = -subSnap[iLocalNode].x;
      }
    }}

    dsvarSnapFile_->stateAdd(*dualSnapBuffer_, timeStamp_);
    dsvarSkip_ = 0;
  }
}

DistrSnapshotNonLinDynamic::DistrSnapshotNonLinDynamic(Domain *domain) :
  MDNLDynamic(domain),
  impl_(nullptr)
{}

void
DistrSnapshotNonLinDynamic::preProcess() {
  MDNLDynamic::preProcess();
  if(!impl_.get()) 
    impl_.reset(new DistrSnapshotNonLinDynamicDetail::RawImpl(this->getDecDomain()));
}

void
DistrSnapshotNonLinDynamic::saveVelocSnapshot(DistrGeomState &state, const DistrVector &veloc) {
  if(domain->solInfo().velocvectPodRom) {
    if(state.getHaveRot()) {
      DistrVector v(veloc);
      state.transform(v, 2, true); // transform convected angular velocity to time derivative of total rotation vector
      impl_->velocSnapshotAdd(v);
    }
    else {
      impl_->velocSnapshotAdd(veloc);
    }
  }

}

void
DistrSnapshotNonLinDynamic::saveAccelSnapshot(DistrGeomState &state, const DistrVector &accel) {
  if(domain->solInfo().accelvectPodRom) {
    if(state.getHaveRot()) {
      DistrVector a(accel);
      state.transform(a, 6, true); // transform convected angular acceleration to second time derivative of total rotation vector
      impl_->accelSnapshotAdd(a);
    }
    else {
      impl_->accelSnapshotAdd(accel);
    }
  }
}

} /* end namespace Rom */
