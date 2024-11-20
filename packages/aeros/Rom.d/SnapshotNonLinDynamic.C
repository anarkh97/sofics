#include "SnapshotNonLinDynamic.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisBinaryFile.h"
#include "NodeDof6Buffer.h"
#include "FileNameInfo.h"
#include "BasisFileStream.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/utilities.hpp>

#include <deque>

#include <cstddef>
#include <algorithm>
#include <memory>

extern GeoSource * geoSource;

namespace Rom {

// Dummy class holding the implementation of SnapshotNonLinDynamic
struct SnapshotNonLinDynamicDetail : private SnapshotNonLinDynamic {
  class sttSnapImpl : public Impl {
  public:
     void lastMidTimeIs(double t);
     void lastDeltaIs(double delta);
     void stateSnapshotAdd(const GeomState &);
     void internalStateSnapshotAdd(const GeomState &);
     void velocSnapshotAdd(const Vector &);
     void accelSnapshotAdd(const Vector &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();

    int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
    int vectorSize() const { return converter_.vectorSize(); }

    explicit sttSnapImpl(Domain *, BasisId::Level level = BasisId::SNAPSHOTS);

  protected:
    template <typename VecType>
    void fillSnapBuffer(const VecType &origin);
  
    const NodeDof6Buffer &snapBuffer() const { return snapBuffer_; }
    const VecNodeDof6Conversion &converter() const { return converter_; }

    int maxSizePodRom() const { return domain_->solInfo().maxSizePodRom; }

  private:
    Domain * domain_;
    int stateSkip_;    
    VecNodeDof6Conversion converter_;
    NodeDof6Buffer snapBuffer_;

    double timeStamp_;

  protected:
    FileNameInfo fileInfo_;
    BasisBinaryOutputFile stateSnapFile_;
  };

  // Implementation with velocity snapshots
  class velSnapImpl : public Impl {
  public:
    explicit velSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);
  
    // Overriden functions
     void lastMidTimeIs(double t);
     void lastDeltaIs(double delta);
     void stateSnapshotAdd(const GeomState &state);
     void internalStateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &);
     void accelSnapshotAdd(const Vector &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();
  
  private:
    Domain * domain_;
    int velSkip_;
    SnapshotNonLinDynamic *parent_;
    VecNodeDof6Conversion vecNodeDof6Conversion_;
    FileNameInfo fileInfo_;
    BasisOutputStream<6> velocitySnapFile_;
  };
  
  // Implementation with acceleration snapshots
  class accSnapImpl : public Impl {
  public:
    explicit accSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);
  
    // Overriden functions
     void lastMidTimeIs(double t);
     void lastDeltaIs(double delta);
     void stateSnapshotAdd(const GeomState &state);
     void internalStateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &);
     void accelSnapshotAdd(const Vector &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();
  
  private:
    Domain * domain_;
    int accSkip_;
    SnapshotNonLinDynamic *parent_;
    VecNodeDof6Conversion vecNodeDof6Conversion_;
    FileNameInfo fileInfo_;
    BasisOutputStream<6> accelerationSnapFile_;
  };
  
  // Implementation with residual snapshots
  class resSnapImpl : public Impl {
  public:
    explicit resSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);
  
    // Overriden functions
     void lastMidTimeIs(double t);
     void lastDeltaIs(double delta);
     void stateSnapshotAdd(const GeomState &state);
     void internalStateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &);
     void accelSnapshotAdd(const Vector &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();
  
  private:
    Domain * domain_;
    int resSkip_;
    SnapshotNonLinDynamic *parent_;
    VecNodeDof6Conversion vecNodeDof6Conversion_;
    FileNameInfo fileInfo_;
    BasisOutputStream<6> residualSnapFile_;
  };
  
  //Implementation with jacobian snapshots
  class jacSnapImpl : public Impl {
  public:
    explicit jacSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);
  
    // Overriden functions
     void lastMidTimeIs(double t);
     void lastDeltaIs(double delta);
     void stateSnapshotAdd(const GeomState &state);
     void internalStateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &);
     void accelSnapshotAdd(const Vector &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();
  
  private:
    Domain * domain_;
    int jacSkip_;
    SnapshotNonLinDynamic *parent_;
    VecNodeDof6Conversion vecNodeDof6Conversion_;
    FileNameInfo fileInfo_;
    BasisOutputStream<6> jacobianSnapFile_;
  };
  
  //Implementation with internal state snapshots
  class isvSnapImpl : public Impl {
  public:
    explicit isvSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);
  
    // Overriden functions
     void lastMidTimeIs(double t);
     void lastDeltaIs(double delta);
     void stateSnapshotAdd(const GeomState &state);
     void internalStateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &);
     void accelSnapshotAdd(const Vector &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();
  
  private:
    Domain * domain_;
    int isvSkip_;
    SnapshotNonLinDynamic *parent_;
    //VecNodeDof6Conversion vecNodeDof6Conversion_;
    //FileNameInfo fileInfo_;
    //BasisOutputStream<6> internalStateSnapFile_;
    int internalStateSnapFile_; 
  };

  private:
  // Dummy constructor to avoid compilation failures
  SnapshotNonLinDynamicDetail(Domain *d) :
    SnapshotNonLinDynamic(d)
  {}
};


SnapshotNonLinDynamicDetail::sttSnapImpl::sttSnapImpl(Domain * domain, BasisId::Level level) :
  domain_(domain),
  stateSkip_(0),
  converter_(*domain->getCDSA()),
  snapBuffer_(dofSetNodeCount()),
  fileInfo_(),
  stateSnapFile_(BasisFileId(fileInfo_, BasisId::STATE, level), dofSetNodeCount(),
                 (geoSource->getCheckFileInfo()->lastRestartFile != 0)),
  timeStamp_(domain->solInfo().initialTime)
{}

SnapshotNonLinDynamicDetail::velSnapImpl::velSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain) :
  domain_(domain),
  velSkip_(0),
  parent_(parent),
  vecNodeDof6Conversion_(*domain->getCDSA()),
  fileInfo_(),
  velocitySnapFile_(BasisFileId(fileInfo_, BasisId::VELOCITY, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                    (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

SnapshotNonLinDynamicDetail::accSnapImpl::accSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain) :
  domain_(domain),
  accSkip_(0),
  parent_(parent),
  vecNodeDof6Conversion_(*domain->getCDSA()),
  fileInfo_(),
  accelerationSnapFile_(BasisFileId(fileInfo_, BasisId::ACCELERATION, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                        (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

SnapshotNonLinDynamicDetail::resSnapImpl::resSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain) :
  domain_(domain),
  resSkip_(0),
  parent_(parent),
  vecNodeDof6Conversion_(*domain->getCDSA()),
  fileInfo_(),
  residualSnapFile_(BasisFileId(fileInfo_, BasisId::RESIDUAL, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                    (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

SnapshotNonLinDynamicDetail::jacSnapImpl::jacSnapImpl(SnapshotNonLinDynamic *parent, Domain * domain) :
  domain_(domain),
  jacSkip_(0),
  parent_(parent),
  vecNodeDof6Conversion_(*domain->getCDSA()),
  fileInfo_(),
  jacobianSnapFile_(BasisFileId(fileInfo_, BasisId::JACOBIAN, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                    (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

SnapshotNonLinDynamicDetail::isvSnapImpl::isvSnapImpl(SnapshotNonLinDynamic *parent, Domain * domain) :
  domain_(domain),
  isvSkip_(0),
  parent_(parent)
  //vecNodeDof6Conversion_(*domain->getCDSA()),
  //fileInfo_(),
  //internalStateSnapFile_(BasisFileId(fileInfo_, BasisId::INTERNALSTATE, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{
  internalStateSnapFile_ = open(domain->solInfo().isvPodRomFile, O_WRONLY | O_CREAT, 0666);
}


void
SnapshotNonLinDynamicDetail::sttSnapImpl::postProcess()
{
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::postProcess()
{
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::postProcess()
{
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::postProcess()
{
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::postProcess()
{
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::postProcess()
{
  // Nothing to do
}

template <typename VecType>
inline
void
SnapshotNonLinDynamicDetail::sttSnapImpl::fillSnapBuffer(const VecType &snap)
{
  converter_.paddedNodeDof6(snap, snapBuffer_);
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::lastMidTimeIs(double t)
{
  // t = t_n + dt*(1 - alphaf)
  timeStamp_ = t;
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::lastDeltaIs(double delta)
{
  // delta = dt/2
  // timeStamp_ = t_{n+1} = t_n + dt = t_n + dt*(1 - alphaf) + dt*alphaf 
  //            = timeStamp_ + 2*delta*alphaf
  timeStamp_ += 2*delta*domain_->solInfo().newmarkAlphaF;
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::lastMidTimeIs(double t)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::lastDeltaIs(double delta)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::lastMidTimeIs(double t)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::lastDeltaIs(double delta)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::lastMidTimeIs(double t)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::lastDeltaIs(double delta)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::lastMidTimeIs(double t)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::lastDeltaIs(double delta)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::lastMidTimeIs(double t)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::lastDeltaIs(double delta)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::stateSnapshotAdd(const GeomState &snap)
{
  ++stateSkip_;
  if(stateSkip_ >= domain_->solInfo().skipState) {
  const CoordSet &refCoords = domain_->getNodes();

  for (int iNode = 0, iNodeEnd = dofSetNodeCount(); iNode != iNodeEnd; ++iNode) {
    double *nodeBuffer = snapBuffer_[iNode];

    const Node *refNode = refCoords[iNode];
    if (refNode) {
      const NodeState &snapNode = snap[iNode];

      // Translational dofs
      nodeBuffer[0] = snapNode.x - refNode->x;
      nodeBuffer[1] = snapNode.y - refNode->y;
      nodeBuffer[2] = snapNode.z - refNode->z;

      // Rotational dofs
      // old method: collect the rescaled rotation vector
      //mat_to_vec(const_cast<double (*)[3]>(snapNode.R), &nodeBuffer[3]);

      // new method: collect the unscaled rotation vector which has already be computed and stored in NodeState::theta
      nodeBuffer[3] = snapNode.theta[0];
      nodeBuffer[4] = snapNode.theta[1];
      nodeBuffer[5] = snapNode.theta[2];

      if(NFrameData *cd = refCoords.dofFrame(iNode)) {
        cd->transformVector6(nodeBuffer);
      }

    } else {
      // Node does not really exist, corresponds to a gap in node numbering
      std::fill_n(nodeBuffer, 6, 0.0);
    }
  }

  stateSnapFile_.stateAdd(snapBuffer_, timeStamp_);
  stateSkip_ = 0;
 }
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::accSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::resSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::internalStateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::velSnapImpl::internalStateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::accSnapImpl::internalStateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::resSnapImpl::internalStateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::internalStateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::internalStateSnapshotAdd(const GeomState &snap)
{
  ++isvSkip_;
  if(isvSkip_ >= domain_->solInfo().skipInternalStateVar) {
    int numElemStates = snap.getTotalNumElemStates();
    Vector elemStates(numElemStates);
    snap.getElemStates(elemStates.data());
    int writeSize = write(internalStateSnapFile_, elemStates.data(), numElemStates*sizeof(double));
    isvSkip_ = 0;
  }
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::velocSnapshotAdd(const Vector &veloc)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::velocSnapshotAdd(const Vector &veloc)
{
  ++velSkip_;
  if(velSkip_ >= domain_->solInfo().skipVeloc) {
    velocitySnapFile_ << veloc;
    velSkip_ = 0;
  }
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::velocSnapshotAdd(const Vector &veloc)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::velocSnapshotAdd(const Vector &veloc)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::velocSnapshotAdd(const Vector &veloc)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::velocSnapshotAdd(const Vector &veloc)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::accelSnapshotAdd(const Vector &accel)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::accelSnapshotAdd(const Vector &accel)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::accelSnapshotAdd(const Vector &accel)
{
  ++accSkip_;
  if(accSkip_ >= domain_->solInfo().skipAccel) {
    accelerationSnapFile_ << accel;
    accSkip_ = 0;
  }
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::accelSnapshotAdd(const Vector &accel)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::accelSnapshotAdd(const Vector &accel)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::accelSnapshotAdd(const Vector &accel)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::handleResidualSnapshot(const Vector &res)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::handleResidualSnapshot(const Vector &res)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::handleResidualSnapshot(const Vector &res)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::handleResidualSnapshot(const Vector &res)
{
  ++resSkip_;
  if(resSkip_ >= domain_->solInfo().skipResidual) {
    residualSnapFile_ << res;
    resSkip_ = 0;
  }
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::handleResidualSnapshot(const Vector &res)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::handleResidualSnapshot(const Vector &res)
{
  //empty
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::handleJacobianSnapshot()
{
  //empty
}

void
SnapshotNonLinDynamicDetail::velSnapImpl::handleJacobianSnapshot()
{
  //empty
}

void
SnapshotNonLinDynamicDetail::accSnapImpl::handleJacobianSnapshot()
{
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::handleJacobianSnapshot()
{
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::handleJacobianSnapshot()
{
  ++jacSkip_;
  if(jacSkip_ >= domain_->solInfo().skipJacobian) {
    Vector snap(parent_->solVecInfo());
    //expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
    jacobianSnapFile_ << snap;
    jacSkip_ = 0;
  }
}

void
SnapshotNonLinDynamicDetail::isvSnapImpl::handleJacobianSnapshot()
{
  //empty
}

SnapshotNonLinDynamic::SnapshotNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain),
  stateImpl_(nullptr),
  internalStateImpl_(nullptr),
  velocImpl_(nullptr),
  accelImpl_(nullptr),
  resImpl_(nullptr),
  jacImpl_(nullptr)
{}

void
SnapshotNonLinDynamic::preProcess()
{
  NonLinDynamic::preProcess();
  if(domain->solInfo().statevectPodRom)
    stateImpl_.reset(new SnapshotNonLinDynamicDetail::sttSnapImpl(this->domain));
  if(domain->solInfo().isvPodRom)
    internalStateImpl_.reset(new SnapshotNonLinDynamicDetail::isvSnapImpl(this,this->domain));
  if(domain->solInfo().velocvectPodRom)
    velocImpl_.reset(new SnapshotNonLinDynamicDetail::velSnapImpl(this,this->domain));
  if(domain->solInfo().accelvectPodRom)
    accelImpl_.reset(new SnapshotNonLinDynamicDetail::accSnapImpl(this,this->domain));
  if(domain->solInfo().residvectPodRom)
    resImpl_.reset(new SnapshotNonLinDynamicDetail::resSnapImpl(this,this->domain));
  if(domain->solInfo().jacobvectPodRom)
    jacImpl_.reset(new SnapshotNonLinDynamicDetail::jacSnapImpl(this,this->domain));
}

void
SnapshotNonLinDynamic::postProcess()
{
  if(domain->solInfo().statevectPodRom)
    stateImpl_->postProcess();
  if(domain->solInfo().isvPodRom)
    internalStateImpl_->postProcess();
  if(domain->solInfo().velocvectPodRom)
    velocImpl_->postProcess();
  if(domain->solInfo().accelvectPodRom)
    accelImpl_->postProcess();
  if(domain->solInfo().residvectPodRom)
    resImpl_->postProcess();
  if(domain->solInfo().jacobvectPodRom)
    jacImpl_->postProcess();
}

void
SnapshotNonLinDynamic::saveMidTime(double t)
{
  if(domain->solInfo().statevectPodRom)
    stateImpl_->lastMidTimeIs(t);
}

void
SnapshotNonLinDynamic::saveDelta(double delta)
{
  if(domain->solInfo().statevectPodRom)
    stateImpl_->lastDeltaIs(delta);
}

void
SnapshotNonLinDynamic::saveStateSnapshot(const GeomState &state)
{
  if(domain->solInfo().statevectPodRom)
    stateImpl_->stateSnapshotAdd(state);
}

void
SnapshotNonLinDynamic::saveInternalStateSnapshot(const GeomState &state)
{
  if(domain->solInfo().isvPodRom)
    internalStateImpl_->internalStateSnapshotAdd(state);
}

void
SnapshotNonLinDynamic::saveVelocitySnapshot(const GeomState &state, const Vector &veloc)
{
  if(domain->solInfo().velocvectPodRom) {
    if(state.getHaveRot()) {
      Vector v(veloc);
      state.transform(v, 2, true); // transform convected angular velocity to time derivative of total rotation vector
      velocImpl_->velocSnapshotAdd(v);
    }
    else {
      velocImpl_->velocSnapshotAdd(veloc);
    }
  }
}

void
SnapshotNonLinDynamic::saveAccelerationSnapshot(const GeomState &state, const Vector &accel)
{
  if(domain->solInfo().accelvectPodRom) {
    if(state.getHaveRot()) {
      Vector a(accel);
      state.transform(a, 6, true); // transform convected angular acceleration to second time derivative of total rotation vector
      accelImpl_->accelSnapshotAdd(a);
    }
    else {
      accelImpl_->accelSnapshotAdd(accel);
    }
  }
}

void
SnapshotNonLinDynamic::handleResidualSnapshot(const Vector &snap)
{
  if(domain->solInfo().residvectPodRom)
    resImpl_->handleResidualSnapshot(snap);
}

int
SnapshotNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time)
{
  if(domain->solInfo().jacobvectPodRom)
    jacImpl_->handleJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time);
}

} /* end namespace Rom */
