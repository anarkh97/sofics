#include "PodProjectionNonLinDynamic.h"

#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "NodeDof6Buffer.h"
#include "VecNodeDof6Conversion.h"

#include "VecBasis.h"
#include "VecBasisFile.h"

#include "BasisOps.h"
#include "VecBasisOps.h"
#include "PodProjectionSolver.h"

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/SysState.h>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Utils.d/DistHelper.h>
#include <Corotational.d/TemperatureState.h>

#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <numeric>

#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Eigen/Sparse>
#endif

extern GeoSource * geoSource;
extern int verboseFlag;

namespace Rom {

// Implementation classes

// Common implementation
class PodProjectionNonLinDynamic::Impl {
public:

  virtual void lastMidTimeIs(double) = 0;
  virtual void lastDeltaIs(double) = 0;
  virtual void stateSnapshotAdd(const GeomState &) = 0;
  virtual void velocSnapshotAdd(const Vector &) = 0;
  virtual void accelSnapshotAdd(const Vector &) = 0;
  virtual void handleResidualSnapshot(const Vector &res) = 0;
  virtual void handleJacobianSnapshot() = 0;

  virtual ~Impl();

protected:
  explicit Impl(PodProjectionNonLinDynamic *parent);

  const ConstrainedDSA &getCDSA() const { return *parent_->domain->getCDSA(); }
  int solVecInfo() const { return parent_->solVecInfo(); }
  const SolverInfo &solInfo() const { return parent_->domain->solInfo(); }
  int getNumLMPC() const { return parent_->domain->getNumLMPC(); }
  LMPCons** getLMPCs() const { return parent_->domain->getLMPC()->data(); }
  PodProjectionSolver *getSolver() { return parent_->getSolver(); }
  void getKtimesU(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                     FullSquareMatrix *kelArray=0) { return parent_->domain->getKtimesU(dsp, bcx, ext_f, eta, kelArray); }
  void createKelArray(FullSquareMatrix *& kel) { return parent_->domain->createKelArray(kel); }
  const PodProjectionSolver *getSolver() const { return parent_->getSolver(); }
    
  PodProjectionNonLinDynamic *parent_;

private:
  // Disallow copy & assignment
  Impl(const Impl &);
  Impl &operator=(const Impl &);
};

PodProjectionNonLinDynamic::Impl::Impl(PodProjectionNonLinDynamic *parent) :
  parent_(parent)
{}

PodProjectionNonLinDynamic::Impl::~Impl() {
  // Nothing to do
}

// Dummy class, used for namespace access
class PodProjectionNonLinDynamicDetail : private PodProjectionNonLinDynamic {
public:
  class BasicImpl;
  class sttSnapImpl;
  class velSnapImpl;
  class accSnapImpl;
  class resSnapImpl;
  class jacSnapImpl;

private:
  // Dummy constructor
  PodProjectionNonLinDynamicDetail();
};

// Basic implementation
class PodProjectionNonLinDynamicDetail::BasicImpl : public PodProjectionNonLinDynamic::Impl {
public:
  explicit BasicImpl(PodProjectionNonLinDynamic *parent);
  
  // Overriden functions
  virtual void lastMidTimeIs(double t);
  virtual void lastDeltaIs(double dt);
  virtual void stateSnapshotAdd(const GeomState &);
  virtual void velocSnapshotAdd(const Vector &);
  virtual void accelSnapshotAdd(const Vector &);
  virtual void handleResidualSnapshot(const Vector &res);
  virtual void handleJacobianSnapshot();

  ~BasicImpl();


protected: 
  VecNodeDof6Conversion vecNodeDof6Conversion_;
  FileNameInfo fileInfo_, fileInfo2_;

  VecNodeDof6Conversion &getDof6Conv() { return vecNodeDof6Conversion_; }

  VecBasis projectionBasis_, dualProjectionBasis_;
};

PodProjectionNonLinDynamicDetail::BasicImpl::~BasicImpl() {}

PodProjectionNonLinDynamicDetail::BasicImpl::BasicImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamic::Impl(parent),
  vecNodeDof6Conversion_(getCDSA()),
  fileInfo_()
{

  // this loop checks how many vectors are in each file for memory allocation
  // if pod size is given as 0, then use all vectors in basis file
  std::vector<int> locBasisVec; 
  for(int j=0; j<solInfo().readInROBorModes.size(); ++j){

     std::string fileName = BasisFileId(fileInfo_, BasisId::STATE, BasisId::POD, j);
     if(solInfo().useMassNormalizedBasis) {
       fileName.append(".massorthonormalized");
     }
     BasisInputStream<6> projectionBasisInput(fileName, vecNodeDof6Conversion_);
     if(solInfo().localBasisSize[j] <=0 || solInfo().localBasisSize[j] > projectionBasisInput.size())
       locBasisVec.push_back(projectionBasisInput.size());
     else
       locBasisVec.push_back(solInfo().localBasisSize[j]);
  
     const_cast<int&>(solInfo().localBasisSize[j]) = locBasisVec[j];
     if(solInfo().readInROBorModes.size() > 1)
       filePrint(stderr, " ... Local Basis %d size %-3d         ...\n", j, locBasisVec[j]);
  }

  // this loop checks how many dual vectors are in each file for memory allocation
  std::vector<int> locDualBasisVec;
  if(!solInfo().readInDualROB.empty()){
    VecNodeDof1Conversion vecNodeDof1Conversion(geoSource->getNumConstraintElementsIeq());
    for(int j = 0; j < solInfo().readInDualROB.size(); ++j){
      std::string fileName = BasisFileId(fileInfo_, BasisId::DUALSTATE, BasisId::POD,j);
      BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);

      if(solInfo().localDualBasisSize[j] <= 0)
        locDualBasisVec.push_back(dualProjectionBasisInput.size());
      else
        locDualBasisVec.push_back(solInfo().localDualBasisSize[j]);

      const_cast<int&>(solInfo().localDualBasisSize[j]) = locDualBasisVec[j];
      if(solInfo().readInDualROB.size() > 1)
        filePrint(stderr, " ... Local Dual Basis %d size %-3d    ...\n", j, locDualBasisVec[j]);
    }
  }

  // this loop reads in vectors and stores them in a single data structure
  for(int j=0; j<solInfo().readInROBorModes.size(); ++j) {
    // Load projection basis
    std::string fileName = BasisFileId(fileInfo_, BasisId::STATE, BasisId::POD, j);
    if(solInfo().useMassNormalizedBasis) {
      if(j==0) filePrint(stderr, " ... Using Mass-normalized Basis    ...\n");
      fileName.append(".massorthonormalized");
    }
    filePrint(stderr," ... Reading %s ...\n",fileName.c_str());
    BasisInputStream<6> projectionBasisInput(fileName, vecNodeDof6Conversion_);


    const int projectionSubspaceSize = locBasisVec[j] ?
                                       std::min(locBasisVec[j], projectionBasisInput.size()) :
                                       projectionBasisInput.size();
  
    readVectors(projectionBasisInput, projectionBasis_, 
                std::accumulate(locBasisVec.begin(), locBasisVec.end(), 0),
                projectionSubspaceSize,
                std::accumulate(locBasisVec.begin(), locBasisVec.begin()+j, 0));
  }
 
 
  filePrint(stderr, " ... Proj. Subspace Dimension = %-3d ...\n", projectionBasis_.vectorCount());
  if(solInfo().readInROBorModes.size() > 1) {
    filePrint(stderr, " ... Number of Local Bases = %-3d    ...\n", solInfo().readInROBorModes.size());
    parent->readLocalBasesCent(vecNodeDof6Conversion_);
    parent->readLocalBasesAuxi();
  }

  if(!solInfo().readInDualROB.empty()) { // this loop reads in dual vectors and stores them in a single data structure if not using modal LMPCs
    VecNodeDof1Conversion vecNodeDof1Conversion(geoSource->getNumConstraintElementsIeq());
    for(int j = 0 ; j < solInfo().readInDualROB.size(); j++) {
      // Load dual projection basis    
      std::string fileName = BasisFileId(fileInfo_, BasisId::DUALSTATE, BasisId::POD,j);
      filePrint(stderr," ... Reading %s ...\n",fileName.c_str());
      BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);
      const int dualProjectionSubspaceSize = locDualBasisVec[j] ?
                                             std::min(locDualBasisVec[j], dualProjectionBasisInput.size()) :
                                             dualProjectionBasisInput.size();
  
      readVectors(dualProjectionBasisInput, dualProjectionBasis_,
                  std::accumulate(locDualBasisVec.begin(), locDualBasisVec.end(), 0),
                  dualProjectionSubspaceSize,
                  std::accumulate(locDualBasisVec.begin(), locDualBasisVec.begin()+j, 0));

    }
    filePrint(stderr, " ... Total Dual Proj. Subspace Dim. = %-3d ...\n", dualProjectionBasis_.vectorCount());
  }

  // Setup solver
  PodProjectionSolver *solver = getSolver();
  solver->projectionBasisIs(projectionBasis_);

  if(geoSource->getNumConstraintElementsIeq() > 0 && !solInfo().modalLMPC) {
    solver->dualProjectionBasisIs(dualProjectionBasis_);
    double dt = solInfo().getTimeStep(), beta = solInfo().newmarkBeta;
    double Kcoef = dt*dt*beta;
    Elemset &eset = geoSource->getPackedEsetConstraintElementIeq();
    ResizeArray<LMPCons *> lmpc(0);
    int numLMPC = 0;
    for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) {
      Element *ele = eset[i];
      int n = ele->getNumMPCs();
      LMPCons **l = ele->getMPCs();
      for(int j = 0; j < n; ++j) {
        lmpc[numLMPC++] = l[j];
      }
      delete [] l;
    }
    solver->addLMPCs(numLMPC, lmpc.data(), Kcoef);
    for(int i=0; i<numLMPC; ++i) delete lmpc[i];
  }

  if(solInfo().modalLMPC) {
    for(int b = 0; b < solInfo().localDualBasisSize.size(); b++)
      std::cout << " ... Dual Basis " << b << " is size " << solInfo().localDualBasisSize[b] << " ... " << std::endl;
    double dt = solInfo().getTimeStep(), beta = solInfo().newmarkBeta;
    double Kcoef = dt*dt*beta;
    int numCols = std::accumulate(solInfo().localDualBasisSize.begin(),solInfo().localDualBasisSize.end(),0);
    std::cout << " ... Total Dual Basis Size is " << numCols << std::endl;
    solver->addModalLMPCs(Kcoef,numCols,geoSource->ROMLMPCVecBegin(),geoSource->ROMLMPCVecEnd());
  }

  solver->factor(); // Delayed factorization
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::stateSnapshotAdd(const GeomState &snap) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::velocSnapshotAdd(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::accelSnapshotAdd(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::handleResidualSnapshot(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::handleJacobianSnapshot() {
  // Nothing to do
}

class PodProjectionNonLinDynamicDetail::sttSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
  public:
     void lastMidTimeIs(double t);
     void lastDeltaIs(double dt);
     void stateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &res);
     void accelSnapshotAdd(const Vector &res);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot(); 

    int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
    int vectorSize() const { return converter_.vectorSize(); }

    explicit sttSnapImpl(Domain * domain, PodProjectionNonLinDynamic * parent);


  protected:

    Domain * domain_;
   
    template <typename VecType>
    void fillSnapBuffer(const VecType &origin);

    const NodeDof6Buffer &snapBuffer() const { return snapBuffer_; }
    const VecNodeDof6Conversion &converter() const { return converter_; }

  private:

    VecNodeDof6Conversion converter_;
    NodeDof6Buffer snapBuffer_;

    double timeStamp_;

  protected:
    BasisBinaryOutputFile stateSnapFile_;
};

// Implementation with velocity snapshots
class PodProjectionNonLinDynamicDetail::velSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit velSnapImpl(PodProjectionNonLinDynamic *parent);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream<6> velocitySnapFile_;
};

// Implementation with acceleration snapshots
class PodProjectionNonLinDynamicDetail::accSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit accSnapImpl(PodProjectionNonLinDynamic *parent);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream<6> accelerationSnapFile_;
};

// Implementation with residual snapshots
class PodProjectionNonLinDynamicDetail::resSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit resSnapImpl(PodProjectionNonLinDynamic *parent);
  
  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream<6> residualSnapFile_;
};

//Implementation with jacobian snapshots
class PodProjectionNonLinDynamicDetail::jacSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit jacSnapImpl(PodProjectionNonLinDynamic *parent);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream<6> jacobianSnapFile_;
};

PodProjectionNonLinDynamicDetail::sttSnapImpl::sttSnapImpl(Domain * domain, PodProjectionNonLinDynamic * parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  domain_(domain),
  converter_(*domain->getCDSA()),
  snapBuffer_(dofSetNodeCount()),
  stateSnapFile_(BasisFileId(fileInfo_, BasisId::STATE, BasisId::SNAPSHOTS), dofSetNodeCount(),
                 (geoSource->getCheckFileInfo()->lastRestartFile != 0)),
  timeStamp_(domain->solInfo().initialTime)
{}

template <typename VecType>
inline
void
PodProjectionNonLinDynamicDetail::sttSnapImpl::fillSnapBuffer(const VecType &snap) {
  converter_.paddedNodeDof6(snap, snapBuffer_);
}

PodProjectionNonLinDynamicDetail::velSnapImpl::velSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  velocitySnapFile_(BasisFileId(fileInfo_, BasisId::VELOCITY, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                    (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

PodProjectionNonLinDynamicDetail::accSnapImpl::accSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  accelerationSnapFile_(BasisFileId(fileInfo_, BasisId::ACCELERATION, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                        (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

PodProjectionNonLinDynamicDetail::resSnapImpl::resSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  residualSnapFile_(BasisFileId(fileInfo_, BasisId::RESIDUAL, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                    (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

PodProjectionNonLinDynamicDetail::jacSnapImpl::jacSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  jacobianSnapFile_(BasisFileId(fileInfo_, BasisId::JACOBIAN, BasisId::SNAPSHOTS), vecNodeDof6Conversion_,
                    (geoSource->getCheckFileInfo()->lastRestartFile != 0))
{}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::lastMidTimeIs(double t) {
  timeStamp_ = t;
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::lastDeltaIs(double dt) {
  timeStamp_ += dt;
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::stateSnapshotAdd(const GeomState &snap) {
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
    } else {
      // Node does not really exist, corresponds to a gap in node numbering
      std::fill_n(nodeBuffer, 6, 0.0);
    }
  }

  stateSnapFile_.stateAdd(snapBuffer_, timeStamp_);
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::stateSnapshotAdd(const GeomState &snap) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::stateSnapshotAdd(const GeomState &snap) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::stateSnapshotAdd(const GeomState &snap) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::stateSnapshotAdd(const GeomState &snap) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  velocitySnapFile_ << veloc;
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::accelSnapshotAdd(const Vector &accel) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::accelSnapshotAdd(const Vector &accel) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::accelSnapshotAdd(const Vector &accel) {
  accelerationSnapFile_ << accel;
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::accelSnapshotAdd(const Vector &accel) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::accelSnapshotAdd(const Vector &accel) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::handleResidualSnapshot(const Vector &res) {
  residualSnapFile_ << res;
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::handleJacobianSnapshot() {
  Vector snap(solVecInfo()); // XXX this is wrong
  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
  jacobianSnapFile_ << snap;
}


// Main class members

PodProjectionNonLinDynamic::PodProjectionNonLinDynamic(Domain *d) :
  NonLinDynamic(d),
  impl_(nullptr),
  sttImpl_(nullptr),
  velImpl_(nullptr),
  accImpl_(nullptr),
  resImpl_(nullptr),
  jacImpl_(nullptr),
  podPostPro(nullptr),
  d0_Big(nullptr),
  v0_Big(nullptr),
  geomState_Big(nullptr),
  refState_Big(nullptr),
  solver_(nullptr),
  localBasisId(0)
{}

PodProjectionNonLinDynamic::~PodProjectionNonLinDynamic() {
  if(podPostPro) delete podPostPro;
  if(d0_Big) delete d0_Big;
  if(v0_Big) delete v0_Big;
  if(geomState_Big) delete geomState_Big;
  if(refState_Big) delete refState_Big;
}

void
PodProjectionNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  
  if (!(solver_ = dynamic_cast<PodProjectionSolver *>(NonLinDynamic::getSolver()))) {
    throw std::runtime_error("Solver must be a PodProjectionSolver");
  }

  if (domain->solInfo().snapshotsPodRom) {
   if(domain->solInfo().statevectPodRom)
    sttImpl_.reset(new PodProjectionNonLinDynamicDetail::sttSnapImpl(this->domain,this)); 
   if(domain->solInfo().velocvectPodRom)
    velImpl_.reset(new PodProjectionNonLinDynamicDetail::velSnapImpl(this));
   if(domain->solInfo().accelvectPodRom)
    accImpl_.reset(new PodProjectionNonLinDynamicDetail::accSnapImpl(this));
   if(domain->solInfo().residvectPodRom) 
    resImpl_.reset(new PodProjectionNonLinDynamicDetail::resSnapImpl(this));
   if(domain->solInfo().jacobvectPodRom)
    jacImpl_.reset(new PodProjectionNonLinDynamicDetail::jacSnapImpl(this));
  } else {
    impl_.reset(new PodProjectionNonLinDynamicDetail::BasicImpl(this));
  }

  podPostPro = new SDDynamPodPostProcessor(domain, bcx, vcx, acx, times);
  GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
  if(domain->solInfo().readInROBorModes.size() == 1)
    podPostPro->printPODSize(projectionBasis.numVectors());
  else
    podPostPro->printPODSize(std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0));
  podPostPro->makeSensorBasis(&projectionBasis);

  if(domain->solInfo().getNLInfo().linearelastic) {
    // there isn't any good reason to use ECSW with "linearelastic"; this is just for testing...
    if(!domain->solInfo().elemLumpPodRom) calculateReducedStiffness(*K, projectionBasis, K_reduced);
    delete K; K = 0;

    // calculate the reduced time-dependent forces (default loadset only) TODO add support for multiple loadsets
    Vector forceFull(NonLinDynamic::solVecInfo(), 0.0);
    domain->computeUnamplifiedExtForce(forceFull, 0);
    if(forceFull.norm() != 0) {
      Vector forceRed(projectionBasis.numVectors());
      reduce(projectionBasis, forceFull, forceRed);
      BCond *nbcModal = new BCond[projectionBasis.numVectors()];
      for(int i=0; i<projectionBasis.numVectors(); ++i) {
        nbcModal[i].setData(i, -1, forceRed[i], BCond::Undefined, 0);
      }
      int numNeumanModal = domain->nNeumannModal();
      domain->setNeumanModal(projectionBasis.numVectors(), nbcModal);
      if(numNeumanModal) delete [] nbcModal;
    }
  }

  if(!domain->solInfo().useMassNormalizedBasis) {
#ifdef USE_EIGEN3
    if(domain->solInfo().modalDIMASS) {
      filePrint(stderr, " ... Reading Reduced Mass Matrix    ...\n");
      std::ifstream matrixin(domain->solInfo().reducedMassFile);
      int n = solver_->projectionBasis().vectorCount();
      VtMV.resize(n,n);
      for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
          matrixin >>VtMV(i,j);
      matrixin.close();
    }
#endif
  }
}

const PodProjectionSolver *
PodProjectionNonLinDynamic::getSolver() const {
  return solver_;
}

PodProjectionSolver *
PodProjectionNonLinDynamic::getSolver() {
  return solver_;
}

int
PodProjectionNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time) {
  if(domain->solInfo().jacobvectPodRom)
  jacImpl_->handleJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time); 
}

double
PodProjectionNonLinDynamic::getResidualNorm(const Vector &residual, ModalGeomState &, double) {
#ifdef USE_EIGEN3
  if(!domain->solInfo().readInDualROB.empty() || domain->solInfo().modalLMPC) {
    const Eigen::VectorXd &fc = solver_->lastReducedConstraintForce();
    return (Eigen::Map<const Eigen::VectorXd>(residual.data(), residual.size()) - fc).norm();
  } else
#endif
  return solver->getResidualNorm(residual);
}

void
PodProjectionNonLinDynamic::handleResidualSnapshot(const Vector &snap) {
  if(domain->solInfo().residvectPodRom)
  resImpl_->handleResidualSnapshot(snap);
}

bool
PodProjectionNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}

void
PodProjectionNonLinDynamic::saveMidTime(double t) {
  if(domain->solInfo().statevectPodRom)
  sttImpl_->lastMidTimeIs(t);
}

void
PodProjectionNonLinDynamic::saveDelta(double dt) {
  if(domain->solInfo().statevectPodRom)
  sttImpl_->lastDeltaIs(dt);
}

void
PodProjectionNonLinDynamic::saveStateSnapshot(const GeomState &state) {
  if(domain->solInfo().statevectPodRom)
  sttImpl_->stateSnapshotAdd(state);
}

void
PodProjectionNonLinDynamic::saveVelocitySnapshot(const Vector &veloc) {
  if(domain->solInfo().velocvectPodRom)
  velImpl_->velocSnapshotAdd(veloc);
}

void
PodProjectionNonLinDynamic::saveAccelerationSnapshot(const Vector &accel) {
  if(domain->solInfo().accelvectPodRom)
  accImpl_->accelSnapshotAdd(accel);
}

int
PodProjectionNonLinDynamic::solVecInfo() const
{
  return getSolver()->basisSize();
}

int
PodProjectionNonLinDynamic::getInitState(Vector &d, Vector &v, Vector &a, Vector &v_p)
{
  // need a copy of d0_Big, v0_Big to use later
  v0_Big = new Vector(NonLinDynamic::solVecInfo(), 0.0);
  d0_Big = new Vector(NonLinDynamic::solVecInfo(), 0.0);

  // d, v, a and v_p are on entry are already initialized to zero
  int numIDisModal = domain->numInitDispModal();
  if(numIDisModal) {
    filePrint(stderr, " ... Using Modal IDISPLACEMENTS     ...\n");
    BCond* iDisModal = domain->getInitDispModal();
    for(int i = 0; i < numIDisModal; ++i) {
      if(iDisModal[i].nnum < d.size())
        d[iDisModal[i].nnum] = iDisModal[i].val;
    }
    const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    projectionBasis.expand(d, *d0_Big);
    initLocalBasis(d);
  }

  int numIVelModal = domain->numInitVelocityModal();
  if(numIVelModal) {
    filePrint(stderr, " ... Using Modal IVELOCITIES        ...\n");
    BCond* iVelModal = domain->getInitVelocityModal();
    const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    if(numIDisModal == 0) initLocalBasis(d);
    // note: It is currently assumed that the modal ivel is defined in the local basis corresponding to the initial state.
    //       Any modal ivel defined in other local bases will simply be ignored.
    for(int i = 0; i < numIVelModal; ++i) {
      if((iVelModal[i].nnum >= projectionBasis.startCol()) && 
         (iVelModal[i].nnum < projectionBasis.startCol()+projectionBasis.numVectors())) {
        v[iVelModal[i].nnum] = iVelModal[i].val;
      }
    }
    projectionBasis.expand(v, *v0_Big);
  }

  // currently, if modal initial conditions are defined then any non-modal initial conditions are ignored
  if(numIDisModal == 0 && numIVelModal == 0) {
    Vector a_Big(NonLinDynamic::solVecInfo(), 0.0),
           v_p_Big(NonLinDynamic::solVecInfo(), 0.0);
  
    NonLinDynamic::getInitState(*d0_Big, *v0_Big, a_Big, v_p_Big);

    if(d0_Big->norm() != 0) reduceDisp(*d0_Big, d);
    initLocalBasis(d);
    // XXX the initial angular velocities are convected, by convention, so they should be transformed
    //     to the first time-derivative of the total rotation vector before reducing
    if(v0_Big->norm() != 0) reduceDisp(*v0_Big, v);
  }

  return domain->solInfo().aeroFlag;
}

void
PodProjectionNonLinDynamic::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n, Vector &v_p, ModalGeomState &geomState)
{
  if(geoSource->getCheckFileInfo()->lastRestartFile) {
    Vector d_n_Big(NonLinDynamic::solVecInfo()),
           v_n_Big(NonLinDynamic::solVecInfo()),
           a_n_Big(NonLinDynamic::solVecInfo()),
           v_p_Big(NonLinDynamic::solVecInfo()),
           q_Big(NonLinDynamic::solVecInfo());

    NonLinDynamic::readRestartFile(d_n_Big, v_n_Big, a_n_Big, v_p_Big, *geomState_Big);
    geomState_Big->setVelocityAndAcceleration(v_n_Big, a_n_Big);

    reduceDisp(d_n_Big, d_n);
    reduceDisp(v_n_Big, v_n);
    reduceDisp(a_n_Big, a_n);
    reduceDisp(v_p_Big, v_p);
    geomState_Big->get_tot_displacement(q_Big);
    reduceDisp(q_Big, geomState.q);
  }
  else if(domain->solInfo().initialTime == 0.0) {
    if(domain->numInitDispModal() || domain->numInitVelocityModal()) {
      geomState.q = d_n;
      geomState_Big->explicitUpdate(domain->getNodes(), *d0_Big);
      geomState_Big->setVelocity(*v0_Big);
    }
  }
}

void
PodProjectionNonLinDynamic::updatePrescribedDisplacement(ModalGeomState *geomState)
{
  if(domain->solInfo().initialTime == 0.0) {
    if(domain->numInitDispModal() == 0 && domain->numInitVelocityModal() == 0) {
      Vector q_Big(NonLinDynamic::solVecInfo());
      NonLinDynamic::updatePrescribedDisplacement(geomState_Big);
      geomState_Big->get_tot_displacement(q_Big);
      reduceDisp(q_Big, geomState->q);
      geomState_Big->setVelocity(*v0_Big);
    }
  }
}

void
PodProjectionNonLinDynamic::getConstForce(Vector &constantForce)
{
  int numNeumanModal = domain->nNeumannModal();
  if(numNeumanModal) {
    filePrint(stderr, " ... Using Reduced Constant Force   ...\n");
    BCond* nbcModal = domain->getNBCModal();
    constantForce.zero();
    for(int i = 0; i < numNeumanModal; ++i) {
      if(nbcModal[i].nnum < constantForce.size()) {
        if(!domain->getMFTT(nbcModal[i].loadsetid)) {
          double loadFactor = (nbcModal[i].loadsetid == -1) ? 1.0 : domain->getLoadFactor(nbcModal[i].loadsetid); // -1 is gravity
          constantForce[nbcModal[i].nnum] += loadFactor*nbcModal[i].val;
        }
      }
    }
  }
  else { // XXX currently, if modal forces are defined then any non-modal forces are ignored
    Vector constantForce_Big(NonLinDynamic::solVecInfo());

    NonLinDynamic::getConstForce(constantForce_Big);

    const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    projectionBasis.reduceAll(constantForce_Big, constantForce);
  }
}

void
PodProjectionNonLinDynamic::getExternalForce(Vector &rhs, Vector &constantForce, int tIndex, double t,
                                             ModalGeomState *geomState, Vector &elementInternalForce,
                                             Vector &aeroForce, double localDelta)
{
  int numNeumanModal = domain->nNeumannModal();
  const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();

  if(numNeumanModal) {
    BCond* nbcModal = domain->getNBCModal();
    rhs = 0.;
    projectionBasis.addLocalPart(constantForce, rhs);
    if(domain->getNumMFTT() > 0) {
      for(int i = 0; i < numNeumanModal; ++i) {
        if((nbcModal[i].nnum >= projectionBasis.startCol()) &&
           (nbcModal[i].nnum < projectionBasis.startCol()+projectionBasis.numVectors())) {
          if(MFTTData *mftt = domain->getMFTT(nbcModal[i].loadsetid)) rhs[nbcModal[i].nnum] += mftt->getVal(t)*nbcModal[i].val;
        }
      }
    }
  }
  else { // XXX currently, if modal forces are defined then any non-modal forces are ignored
    Vector rhs_Big(NonLinDynamic::solVecInfo()),
           constantForce_Big(NonLinDynamic::solVecInfo(), 0.0),
           aeroForce_Big(NonLinDynamic::solVecInfo());

    NonLinDynamic::getExternalForce(rhs_Big, constantForce_Big, tIndex, t, geomState_Big, elementInternalForce, 
                                    aeroForce_Big, localDelta);

    if(rhs_Big.norm() != 0) {
      projectionBasis.reduce(rhs_Big, rhs, false);
      projectionBasis.addLocalPart(constantForce, rhs);
    }
    else {
      rhs = 0.;
      projectionBasis.addLocalPart(constantForce, rhs);
    }
  }
}

void
PodProjectionNonLinDynamic::getIncDisplacement(ModalGeomState *geomState, Vector &du, ModalGeomState *refState,
                                               bool zeroRot)
{
  geomState->get_inc_displacement(du, *refState, zeroRot);
}

double
PodProjectionNonLinDynamic::formRHScorrector(Vector &inc_displacement, Vector &velocity, Vector &acceleration,
                                             Vector &residual, Vector &rhs, ModalGeomState *geomState, double localDelta)
{
  double beta, gamma, alphaf, alpham, dt = 2*localDelta;
  getNewmarkParameters(beta, gamma, alphaf, alpham);

  if(domain->solInfo().useMassNormalizedBasis && C == 0) {
    if(domain->solInfo().order == 1)
      rhs = -1.0*inc_displacement;
    else {
      if(domain->solInfo().quasistatic) rhs = 0.;
      else rhs = -(1-alpham)/(1-alphaf)*inc_displacement + dt*(1-alpham)*velocity + dt*dt*((1-alpham)/2-beta)*acceleration;
    }
  }
#ifdef USE_EIGEN3
  else if(!domain->solInfo().useMassNormalizedBasis && domain->solInfo().modalDIMASS && C == 0) {
    if(domain->solInfo().order == 1)
      rhs = -1.0*inc_displacement;
    else
      rhs = -(1-alpham)/(1-alphaf)*inc_displacement + dt*(1-alpham)*velocity + dt*dt*((1-alpham)/2-beta)*acceleration;
    Eigen::Map<Eigen::VectorXd> rhsMap(rhs.data(),rhs.size());
    rhsMap = VtMV*rhsMap;
  }
#endif
  else {
    // this can be improved by pre-computing V^T*M*V and V^T*C*V
    Vector inc_displacement_Big(NonLinDynamic::solVecInfo()),
           velocity_Big(NonLinDynamic::solVecInfo()),
           acceleration_Big(NonLinDynamic::solVecInfo()),
           residual_Big(NonLinDynamic::solVecInfo(), 0.0),
           rhs_Big(NonLinDynamic::solVecInfo());

    const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    projectionBasis.fullExpand(inc_displacement, inc_displacement_Big);
    projectionBasis.fullExpand(velocity, velocity_Big);
    projectionBasis.fullExpand(acceleration, acceleration_Big);

    NonLinDynamic::formRHScorrector(inc_displacement_Big, velocity_Big, acceleration_Big,
                                    residual_Big, rhs_Big, geomState_Big, localDelta);

    projectionBasis.reduce(rhs_Big, rhs);
  }
  if(domain->solInfo().order == 1)
    rhs += localDelta*residual;
  else
    rhs += (dt*dt*beta)*residual;
  return rhs.norm();
}

void
PodProjectionNonLinDynamic::formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs,
                                             ModalGeomState &geomState, double midtime, double localDelta)
{
  std::cerr << "PodProjectionNonLinDynamic::formRHSpredictor is not implemented\n";
}

void
PodProjectionNonLinDynamic::formRHSinitializer(Vector &fext, Vector &velocity, Vector &elementInternalForce,
                                               ModalGeomState &geomState, Vector &rhs, ModalGeomState *refState)
{
  // XXX For the case of initial values for the rotation, the transformed mass matrix
  //     should be used to solve for the initial accelerations.
  Vector fext_Big(NonLinDynamic::solVecInfo(), 0.0),
         velocity_Big(NonLinDynamic::solVecInfo()),
         rhs_Big(NonLinDynamic::solVecInfo());

  const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
  projectionBasis.expand(velocity, velocity_Big);
  NonLinDynamic::formRHSinitializer(fext_Big, velocity_Big, elementInternalForce, *geomState_Big, rhs_Big, refState_Big);

  projectionBasis.reduce(rhs_Big, rhs);
  rhs += fext;
}

ModalGeomState*
PodProjectionNonLinDynamic::createGeomState()
{
  if(domain->solInfo().soltyp == 2)
    geomState_Big = new TemperatureState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes());
  else
    geomState_Big = new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes(), &domain->getElementSet(),
                                  domain->getNodalTemperatures());

  return new ModalGeomState(solVecInfo());
}

ModalGeomState*
PodProjectionNonLinDynamic::copyGeomState(ModalGeomState *geomState)
{
  if(domain->solInfo().soltyp == 2)
    refState_Big = new TemperatureState(* (TemperatureState *) geomState_Big);
  else
    refState_Big = new GeomState(*geomState_Big);
  return new ModalGeomState(*geomState);
}

void
PodProjectionNonLinDynamic::updateStates(ModalGeomState *refState, ModalGeomState &geomState, double time)
{
  if((!domain->solInfo().getNLInfo().linearelastic && (geomState_Big->getHaveRot() || geomState_Big->getTotalNumElemStates() > 0))
     || domain->solInfo().readInROBorModes.size() > 1) {
    // updateStates is called after midpoint update (i.e. once per timestep)
    // so it is a convenient place to update and copy geomState_Big, if necessary
    const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    if(domain->solInfo().readInROBorModes.size() == 1) {
      // note: for local bases method, geomState_Big has already been updated in setLocalBasis
      Vector q_Big(NonLinDynamic::solVecInfo());
      projectionBasis.expand(geomState.q, q_Big);
      geomState_Big->explicitUpdate(domain->getNodes(), q_Big);
    }

    if(geomState_Big->getHaveRot()) {
      Vector vel_Big(NonLinDynamic::solVecInfo()), acc_Big(NonLinDynamic::solVecInfo());
      projectionBasis.expand(geomState.vel, vel_Big);
      geomState_Big->setVelocity(vel_Big);
      projectionBasis.expand(geomState.acc, acc_Big);
      geomState_Big->setAcceleration(acc_Big);
    }

    if(geomState_Big->getTotalNumElemStates() > 0) 
      domain->updateStates(refState_Big, *geomState_Big, allCorot, time);

    *refState_Big = *geomState_Big;
  }
}

int
PodProjectionNonLinDynamic::selectLocalBasis(Vector &q)
{// this function determines which basis to use
#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && HAS_CXX11_LAMBDA
  if(d.size() > 0 && w.size() > 0) {
    // modelIII: fast implementation using pre-computed auxiliary quantities
    const int Nv = domain->solInfo().readInROBorModes.size();
    const int k  = q.size();
    Eigen::Map<Eigen::VectorXd> qi(q.data(),k);
    std::vector<int> s(Nv); for(int m=0; m<Nv; ++m) s[m] = m;
    std::sort(s.begin(), s.end(), [&](int m, int p) {
      return (p>m && (w(m,p).dot(qi) + d(m,p)) < 0) || (m>p && (w(p,m).dot(qi) + d(p,m)) > 0);
    });
    if(verboseFlag) std::cerr << " ... Selecting Local Basis # " << s[0] << "     ...\n";
    return s[0];
  } else if(uc.size() > 0) {
    // modelII: slow implementation of using cluster centroids.
    Vector q_Big(NonLinDynamic::solVecInfo());
    GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    projectionBasis.fullExpand(q, q_Big);

    Eigen::Map<Eigen::VectorXd> u(q_Big.data(),NonLinDynamic::solVecInfo());
    Eigen::MatrixXd::Index j;
    double minNorm = (uc.colwise()-u).colwise().norm().minCoeff(&j); 
    if(verboseFlag) std::cerr << " ... Selecting Local Basis # " << j << "     ...\n";
    return int(j);
  }
  else {
    std::cerr << " *** Error: cluster centroids or auxiliary quantities required to select local basis.\n";
    exit(-1);
  }
#else
  std::cerr << " *** Error: Local bases requires Aero-S to be built with Eigen library and C++11 support.\n"; 
  exit(-1);
#endif
}

void
PodProjectionNonLinDynamic::initLocalBasis(Vector &q0)
{
#ifdef USE_EIGEN3
  // Local bases
  if(domain->solInfo().readInROBorModes.size() > 1) {
    int oldLBI = localBasisId;
    localBasisId = selectLocalBasis(q0);
    if(oldLBI != localBasisId) 
      std::cerr << " selecting local basis " << localBasisId << "               " << std::endl;
    GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
    int blockCols = domain->solInfo().localBasisSize[localBasisId];
    int startCol = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.begin()+localBasisId, 0);
    getSolver()->setLocalBasis(startCol, blockCols);
    projectionBasis.localBasisIs(startCol, blockCols);
    setLocalReducedMesh(localBasisId);
    // if multiple dual basis are given, select corresponding dual basis, obviously local dual basis can only be done with local primal basis
    if(domain->solInfo().localDualBasisSize.size() > 1) {
      blockCols = domain->solInfo().localDualBasisSize[localBasisId];
      startCol = std::accumulate(domain->solInfo().localDualBasisSize.begin(), domain->solInfo().localDualBasisSize.begin()+localBasisId, 0);
      getSolver()->setLocalDualBasis(startCol,blockCols);
      if(!domain->solInfo().modalLMPC){ //if using modal LMPCs, skip this since W is never explicitly read in
        GenVecBasis<double> &dualProjectionBasis = solver_->dualProjectionBasis();
        dualProjectionBasis.localBasisIs(startCol, blockCols);
      }
    }
  }
#endif
}

void
PodProjectionNonLinDynamic::setLocalBasis(ModalGeomState *refState, ModalGeomState *geomState, Vector &q_n, Vector &vel, Vector &acc)
{
#ifdef USE_EIGEN3
  // Local bases: this block set new local basis and orthogonaly projects onto new basis
  if(domain->solInfo().readInROBorModes.size() > 1) {

    int j = selectLocalBasis(geomState->q);

    Vector dq_Big(NonLinDynamic::solVecInfo());
    GenVecBasis<double> &projectionBasis = solver_->projectionBasis(); 

    Vector dq = geomState->q - q_n;
    projectionBasis.expand(dq, dq_Big, false); // XXX not using compressed basis
    geomState_Big->update(*refState_Big, dq_Big, 2); 

    if(j != localBasisId) { // if a new local basis has been selected, update the things
      std::cerr << "\r selecting local basis " << j << "               " << std::endl;
      Vector vel_Big(NonLinDynamic::solVecInfo()),
             acc_Big(NonLinDynamic::solVecInfo());

      if(VtV.size() == 0) {
        projectionBasis.expand(vel, vel_Big, false);
        projectionBasis.expand(acc, acc_Big, false);
      }

      int blockCols = domain->solInfo().localBasisSize[j];
      int startCol = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.begin()+j, 0);
      getSolver()->setLocalBasis(startCol, blockCols);
      projectionBasis.localBasisIs(startCol, blockCols);
      setLocalReducedMesh(j);

      if(domain->solInfo().localDualBasisSize.size() > 1) { // if multiple dual bases are given, execute this block
        blockCols = domain->solInfo().localDualBasisSize[j];
        startCol = std::accumulate(domain->solInfo().localDualBasisSize.begin(), domain->solInfo().localDualBasisSize.begin()+j, 0); 
        getSolver()->setLocalDualBasis(startCol, blockCols);
        if(!domain->solInfo().modalLMPC){// we don't have W if using modal LMPCs
          GenVecBasis<double> &dualProjectionBasis = solver_->dualProjectionBasis();
          dualProjectionBasis.localBasisIs(startCol, blockCols);
        }
      }

      if(VtV.size() == 0) {
        vel = 0.; reduceDisp(vel_Big, vel);
        acc = 0.; reduceDisp(acc_Big, acc);
      } else {
        // using precomputed Vi^T*Vj (or Vi^T*M*Vj for M-orthogonal local bases)
        projectLocalBases(localBasisId, j, vel);
        projectLocalBases(localBasisId, j, acc);
      }

      refState->setVelocityAndAcceleration(vel, acc);
      geomState->setVelocityAndAcceleration(vel, acc);

      localBasisId = j;
    }
  }
#endif
}

void
PodProjectionNonLinDynamic::readLocalBasesCent(const VecNodeDof6Conversion &vecNodeDof6Conversion)
{
// this block reads in the centriods of the local bases
#ifdef USE_EIGEN3
  uc.resize(NonLinDynamic::solVecInfo(), domain->solInfo().readInLocalBasesCent.size());
  for(int i=0; i<domain->solInfo().readInLocalBasesCent.size(); ++i) {
    BasisInputStream<6> in(domain->solInfo().readInLocalBasesCent[i], vecNodeDof6Conversion);
    in >> uc.col(i).data();
  }
#endif
}

void
PodProjectionNonLinDynamic::readLocalBasesAuxi()
{ // this block is only executed when the precomputed auxilary matrices are given in the input file
  if(domain->solInfo().readInLocalBasesAuxi.size() > 0) {
#ifdef USE_EIGEN3
    const int Nv = domain->solInfo().readInROBorModes.size();
    const int k = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0);
    VtV.resize(Nv,Nv); //grammian
    d.resize(Nv,Nv);
    w.resize(Nv,Nv);
    for(int i=0; i<Nv; ++i) {
      int ki = domain->solInfo().localBasisSize[i];
      for(int j=i+1; j<Nv; ++j) {
        int kj = domain->solInfo().localBasisSize[j];
        std::string fileName = domain->solInfo().readInLocalBasesAuxi[std::make_pair(i,j)];
        std::ifstream file(fileName.c_str());
        VtV(i,j).resize(ki,kj);
        for(int irow = 0; irow < ki; ++irow) {
          for(int jcol = 0; jcol < kj; ++jcol) {
            file >> VtV(i,j)(irow,jcol);
          }
        }
        file >> d(i,j);
        w(i,j).resize(k);
        for(int irow=0; irow<k; ++irow) {
          file >> w(i,j)[irow];
        }
        file.close();
      }
    }
#endif
  }
}

void
PodProjectionNonLinDynamic::projectLocalBases(int i, int j, Vector &q)
{
#ifdef USE_EIGEN3
  int ki = domain->solInfo().localBasisSize[i];
  int pi = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.begin()+i, 0);
  Eigen::Map<Eigen::VectorXd> qi(q.data()+pi,ki);

  int kj = domain->solInfo().localBasisSize[j];
  int pj = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.begin()+j, 0);
  Eigen::Map<Eigen::VectorXd> qj(q.data()+pj,kj);

  if(j < i) qj = VtV(j,i)*qi;
  else qj = VtV(i,j).transpose()*qi;

  qi.setZero();
#endif
}

double
PodProjectionNonLinDynamic::getStiffAndForce(ModalGeomState &geomState, Vector &residual,
                                             Vector &elementInternalForce, double t, ModalGeomState *refState, bool forceOnly)
{
  Vector r(solVecInfo(),0.);
  if(domain->solInfo().getNLInfo().linearelastic == 1 && domain->getFollowedElemList().empty() && !domain->solInfo().elemLumpPodRom) {
    K_reduced.mult(geomState.q, r);
    residual -= r;
  }
  else {
    Vector q_Big(NonLinDynamic::solVecInfo()),
           residual_Big(NonLinDynamic::solVecInfo(), 0.0);
    const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();

    if(domain->solInfo().readInROBorModes.size() > 1) {
      Vector dq = geomState.q - refState->q;
      projectionBasis.expand(dq, q_Big);
      geomState_Big->update(*refState_Big, q_Big, 2);
    }
    else {
      projectionBasis.expand(geomState.q, q_Big);
      geomState_Big->explicitUpdate(domain->getNodes(), q_Big);
    }

    NonLinDynamic::getStiffAndForce(*geomState_Big, residual_Big, elementInternalForce, t, refState_Big, forceOnly);

    projectionBasis.reduce(residual_Big, r);
    residual += r;
  }

#ifdef USE_EIGEN3
  if(geoSource->getNumConstraintElementsIeq() || domain->solInfo().modalLMPC) {
    solver_->activateContact();
    if(geoSource->getLmpcFlag() || domain->solInfo().modalLMPC) { // linear
      solver_->updateLMPCs(geomState.q);
    }
    else { // nonlinear, or mixed
      Elemset &eset = geoSource->getPackedEsetConstraintElementIeq();
      ResizeArray<LMPCons *> lmpc(0);
      int numLMPC = 0;
      for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) {
        Element *ele = eset[i];
        static_cast<MpcElement*>(ele)->update(refState_Big, *geomState_Big, domain->getNodes(), t);
        int n = ele->getNumMPCs();
        LMPCons **l = ele->getMPCs();
        for(int j = 0; j < n; ++j) {
          lmpc[numLMPC++] = l[j];
        }
        delete [] l;
      }
      double dt = domain->solInfo().getTimeStep(), beta = domain->solInfo().newmarkBeta;
      double Kcoef = dt*dt*beta;
      solver_->addLMPCs(numLMPC, lmpc.data(), Kcoef);
      for(int i=0; i<numLMPC; ++i) delete lmpc[i];
    }
  }
#endif
  return residual.norm();
}

void
PodProjectionNonLinDynamic::reBuild(ModalGeomState &geomState, int iteration, double localDelta, double t)
{
  NonLinDynamic::reBuild(*geomState_Big, iteration, localDelta, t);
}

void
PodProjectionNonLinDynamic::dynamCommToFluid(ModalGeomState *geomState, ModalGeomState *bkGeomState, Vector &velocity,
                                             Vector &bkVelocity, Vector &vp, Vector &bkVp, int step, int parity,
                                             int aeroAlg, double time)
{
  // XXX this could be more implemented more efficiently (only called for AERO)
  Vector velocity_Big(NonLinDynamic::solVecInfo()),
         bkVelocity_Big(NonLinDynamic::solVecInfo()),
         vp_Big(NonLinDynamic::solVecInfo()),
         bkVp_Big(NonLinDynamic::solVecInfo());
  GeomState *bkGeomState_Big = nullptr;
  const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
  projectionBasis.expand(velocity, velocity_Big);
  projectionBasis.expand(bkVelocity, bkVelocity_Big);
  projectionBasis.expand(vp, vp_Big);
  projectionBasis.expand(bkVp, bkVp_Big);

  NonLinDynamic::dynamCommToFluid(geomState_Big, bkGeomState_Big, velocity_Big, bkVelocity_Big, vp_Big, bkVp_Big, step, parity, aeroAlg, time);
}

void
PodProjectionNonLinDynamic::dynamOutput(ModalGeomState *geomState, Vector &velocity, Vector &vp, double time,
                                        int step, Vector &force, Vector &aeroF, Vector &acceleration,
                                        ModalGeomState *refState) const
{
  DynamMat dynOps;
  SysState<Vector> systemState(geomState->q, velocity, acceleration, vp);
  podPostPro->dynamOutput(step+1, time, dynOps, force, &aeroF, systemState);
}

void
PodProjectionNonLinDynamic::getConstraintMultipliers(ModalGeomState &geomState)
{
  NonLinDynamic::getConstraintMultipliers(*geomState_Big);
}

void
PodProjectionNonLinDynamic::initializeParameters(int step, ModalGeomState *geomState)
{
  NonLinDynamic::initializeParameters(step, geomState_Big);
}

void
PodProjectionNonLinDynamic::updateParameters(ModalGeomState *geomState)
{
  NonLinDynamic::updateParameters(geomState_Big);
}

bool
PodProjectionNonLinDynamic::checkConstraintViolation(double &err, ModalGeomState *geomState)
{
  return NonLinDynamic::checkConstraintViolation(err, geomState_Big);
}

void
PodProjectionNonLinDynamic::expandForce(Vector &fr, Vector &f) const
{
  const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
  if(domain->solInfo().useMassNormalizedBasis) {
    // note: this case should not be relied upon for modelIII with reduced mesh, because the full mass matrix is not available
    Vector Vfr(NonLinDynamic::solVecInfo());
    projectionBasis.fullExpand(fr, Vfr);
    M->mult(Vfr, f);
  }
  else {
    projectionBasis.fullExpand(fr, f);
  }
}

void
PodProjectionNonLinDynamic::reduceDisp(Vector &d, Vector &dr) const
{
  const GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
  if(domain->solInfo().useMassNormalizedBasis) {
    Vector Md(NonLinDynamic::solVecInfo());
    M->mult(d, Md);
    projectionBasis.reduce(Md, dr, false);
  }
  else {
    projectionBasis.reduce(d, dr, false);
  }
}

} /* end namespace Rom */
