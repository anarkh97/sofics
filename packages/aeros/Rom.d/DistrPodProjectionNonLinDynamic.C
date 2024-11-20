#ifdef USE_EIGEN3
#include "DistrPodProjectionNonLinDynamic.h"

#include "DistrBasisFile.h"
#include "FileNameInfo.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "VecBasis.h"
#include "DistrVecBasis.h"

#include <Driver.d/Domain.h>
#include <Driver.d/DecDomain.h>
#include <Driver.d/GeoSource.h>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Threads.d/PHelper.h>
#include <Threads.d/Paral.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/SysState.h>

#include <cstddef>
#include <memory>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

//Constructor
DistrPodProjectionNonLinDynamic::DistrPodProjectionNonLinDynamic(Domain *domain) :
  MDNLDynamic(domain),
  podPostPro(NULL),
  d0_Big(NULL),
  v0_Big(NULL),
  geomState_Big(NULL),
  refState_Big(NULL),
  solver_(NULL),  
  localBasisId(0),
  resetFromClean(false)
{}

//Destructor
DistrPodProjectionNonLinDynamic::~DistrPodProjectionNonLinDynamic(){
  if(podPostPro)    delete podPostPro;
  if(d0_Big)        delete d0_Big;
  if(v0_Big)        delete v0_Big;
  if(geomState_Big) delete geomState_Big;
  if(refState_Big)  delete refState_Big;
}

// distributed reduced vector info
DistrInfo&
DistrPodProjectionNonLinDynamic::solVecInfo()
{
  return reducedInfo;
}

void
DistrPodProjectionNonLinDynamic::preProcess() {
  MDNLDynamic::preProcess();

  if (!(solver_ = dynamic_cast<GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *>(MDNLDynamic::getSolver()))) {
    throw std::runtime_error("Solver must be EiGalerkinProjectionSolver");
  }


  if(projectionBasis_.data() != NULL) {
    resetFromClean = true;
  }
  // read in primal and dual bases
  {
    std::vector<int> locBasisVec;
    for(int rob=0; rob<domain->solInfo().readInROBorModes.size(); ++rob) {
      FileNameInfo fileInfo;
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD,rob);
      fileName.append(".massorthonormalized");
      DistrBasisInputFile BasisFile(fileName);
      int locSize = domain->solInfo().localBasisSize[rob] ?
                    std::min(domain->solInfo().localBasisSize[rob], BasisFile.stateCount()) :
                    BasisFile.stateCount();
      locBasisVec.push_back(locSize);
      if(verboseFlag && domain->solInfo().readInROBorModes.size()>1 && !resetFromClean) 
        filePrint(stderr, " ... Local Basis %d size %-3d ...\n",rob,locBasisVec[rob]);
    }

    //initialize helper objects for reading in distributed basis vectors
    const int projectionSubspaceSize = std::accumulate(locBasisVec.begin(),locBasisVec.end(),0);
    if(!resetFromClean) filePrint(stderr, " ... Proj. Subspace Dimension = %-3d ...\n", projectionSubspaceSize);
    projectionBasis_.dimensionIs(projectionSubspaceSize, MDNLDynamic::solVecInfo());

    DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

    typedef PtrPtrIterAdapter<SubDomain> SubDomIt;

    DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                     SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
    DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  

    // this loop reads in vectors and stores them in a single Distributed Basis structure
    DistrVecBasis::iterator it = projectionBasis_.begin();
    for(int rob=0; rob<domain->solInfo().readInROBorModes.size(); ++rob){
      FileNameInfo fileInfo;
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD,rob);
      fileName.append(".massorthonormalized");
      DistrBasisInputFile BasisFile(fileName); //read in mass-normalized basis

      if(verboseFlag) filePrint(stderr, " ... Reading basis from file %s ...\n", fileName.c_str());
      for (int currentVec = 0; currentVec<locBasisVec[rob]; ++currentVec) {
        assert(BasisFile.validCurrentState());

        BasisFile.currentStateBuffer(buffer);
        converter.vector(buffer, *it); 
        it++;  
        BasisFile.currentStateIndexInc();
      }
    }
   
    solver_->projectionBasisIs(projectionBasis_);

    if(!domain->solInfo().readInDualROB.empty()) { 
      // this loop reads in dual vectors and stores them in a single data structure if not using modal LMPCs 
      // check size of each local basis and allocate memory
      std::vector<int> dualLocBasisVec;
      for(int rob=0; rob<domain->solInfo().readInROBorModes.size(); ++rob) {
        FileNameInfo fileInfo;
        std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD,rob);
        DistrBasisInputFileTemplate<1> BasisFile(fileName);
        int locSize = domain->solInfo().localDualBasisSize[rob] ?
                      std::min(domain->solInfo().localDualBasisSize[rob], BasisFile.stateCount()) :
                      BasisFile.stateCount();
        dualLocBasisVec.push_back(locSize);
        if(verboseFlag && domain->solInfo().readInDualROB.size()>1 && !resetFromClean) 
          filePrint(stderr, " ... Local Dual Basis %d size %-3d ...\n",rob,dualLocBasisVec[rob]);
      }

      //initialize helper objects for reading in distributed basis vectors
      const int dualProjectionSubspaceSize = std::accumulate(dualLocBasisVec.begin(),dualLocBasisVec.end(),0);
      if(!resetFromClean) filePrint(stderr, " ... Dual Proj. Subspace Dimension = %-3d ...\n", dualProjectionSubspaceSize);
      dualMapVectors_.resize(dualProjectionSubspaceSize); 
      
      std::vector<std::map<int,double> >::iterator dualMapIt = dualMapVectors_.begin();
      for(int rob=0; rob<domain->solInfo().readInDualROB.size(); ++rob){
        FileNameInfo fileInfo;
        std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD,rob);
        DistrBasisInputFileTemplate<1> BasisFile(fileName); //read in dual variable basis
        if(verboseFlag) filePrint(stderr, " ... Reading dual basis from file %s ...\n", fileName.c_str());
        for (int currentVec = 0; currentVec<dualLocBasisVec[rob]; ++currentVec) {
          assert(BasisFile.validCurrentState());

          BasisFile.currentStateBuffer(*dualMapIt);
          dualMapIt++;
 
          BasisFile.currentStateIndexInc();
        }
      }
      solver_->dualProjectionBasisIs(dualMapVectors_);
    } 

    if(domain->solInfo().modalLMPC) {
      for(int b = 0; b < domain->solInfo().localDualBasisSize.size(); b++)
        filePrint(stderr," ... Dual Basis %d is size %d ... \n",b,domain->solInfo().localDualBasisSize[b]);
      double dt = domain->solInfo().getTimeStep(), beta = domain->solInfo().newmarkBeta;
      double Kcoef = dt*dt*beta;
      int numCols = std::accumulate(domain->solInfo().localDualBasisSize.begin(),domain->solInfo().localDualBasisSize.end(),0);
      filePrint(stderr," ... Total Dual Basis Size is %d ...\n",numCols);
      solver_->addModalLMPCs(Kcoef,numCols,geoSource->ROMLMPCVecBegin(),geoSource->ROMLMPCVecEnd());
    }

    // if performing local basis analysis and centroids are provided, read in centroids for basis switching
    readLocalBasesCent(converter, buffer, locBasisVec);

    // if auxillary bases are provided, read in data structures for fast switching
    // every process owns a copy of VtV, w, and d
    readLocalBasesAuxi();
  }

  //preProcessing for reduced solution vector information
  if(!resetFromClean){
    reducedInfo.domLen = new int[MDNLDynamic::solVecInfo().numDom];
    reducedInfo.numDom = MDNLDynamic::solVecInfo().numDom;
    int totLen = 0;
    for(int iSub = 0; iSub < MDNLDynamic::solVecInfo().numDom; ++iSub) {
      if(domain->solInfo().readInROBorModes.size() == 1)
        reducedInfo.domLen[iSub] = (iSub==0) ? projectionBasis_.numVec() : 0;
      else
        reducedInfo.domLen[iSub] = (iSub==0) ? std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0) : 0;
      totLen += reducedInfo.domLen[iSub];
    }

    reducedInfo.len = totLen;
    reducedInfo.setMasterFlag();
  }
 
  // initialized postprocessor for outputing reduced coordinates
  if(!resetFromClean) {
    if(!(podPostPro = new MultiDomDynPodPostProcessor(decDomain, times, geomState_Big, allCorot)))
      throw std::runtime_error("Pod Post Processor did not bind\n");

    if(domain->solInfo().readInROBorModes.size() == 1){
      podPostPro->printPODSize(projectionBasis_.numVec());
    } else {
      int podSize = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0);
      podPostPro->printPODSize(podSize);
    }
    podPostPro->makeSensorBasis(&projectionBasis_);

    if(domain->solInfo().getNLInfo().linearelastic) {
      // TODO for doing linear elastic analysis with nonlinear driver
    }
  }

  solver_->refactor();
}

void
DistrPodProjectionNonLinDynamic::readRestartFile(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                                                 DistrVector &v_p, DistrModalGeomState &geomState){
  if(geoSource->getCheckFileInfo()->lastRestartFile) {
    DistrVector d_n_Big(MDNLDynamic::solVecInfo()),
                v_n_Big(MDNLDynamic::solVecInfo()),
                a_n_Big(MDNLDynamic::solVecInfo()),
                v_p_Big(MDNLDynamic::solVecInfo()),
                  q_Big(MDNLDynamic::solVecInfo());
           
    MDNLDynamic::readRestartFile(d_n_Big, v_n_Big, a_n_Big, v_p_Big, *geomState_Big);
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
      geomState_Big->explicitUpdate(decDomain, *d0_Big);
      geomState_Big->setVelocity(*v0_Big);
    } 
  }
}

int
DistrPodProjectionNonLinDynamic::getInitState(DistrVector &d, DistrVector& v, DistrVector &a, DistrVector &v_p)
{
  v0_Big = new DistrVector(MDNLDynamic::solVecInfo()); v0_Big->zero();
  d0_Big = new DistrVector(MDNLDynamic::solVecInfo()); d0_Big->zero();

  int numIDisModal = domain->numInitDispModal();
  if(numIDisModal) {
    filePrint(stderr, " ... Using Modal IDISPLACEMENTS     ...\n");
    BCond* iDisModal = domain->getInitDispModal();
    for(int i = 0; i < numIDisModal; ++i) {
      if(iDisModal[i].nnum < d.size())
        d[iDisModal[i].nnum] = iDisModal[i].val;
    }
    GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
    projectionBasis.expand(d, *d0_Big);
    initLocalBasis(d);
  }

  int numIVelModal = domain->numInitVelocityModal();
  if(numIVelModal) {
    filePrint(stderr, " ... Using Modal IVELOCITIES        ...\n");
    BCond* iVelModal = domain->getInitVelocityModal();
    GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
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

  if(numIDisModal == 0 && numIVelModal == 0) {
    DistrVector a_Big(MDNLDynamic::solVecInfo()),
                v_p_Big(MDNLDynamic::solVecInfo());
    a_Big.zero();
    v_p_Big.zero();

    MDNLDynamic::getInitState(*d0_Big, *v0_Big, a_Big, v_p_Big);

    if(d0_Big->norm() != 0) reduceDisp(*d0_Big, d);
    initLocalBasis(d);
    // XXX the initial angular velocities are convected, by convention, so they should be transformed
    //     to the first time-derivative of the total rotation vector before reducing
    if(v0_Big->norm() != 0) reduceDisp(*v0_Big, v);
  }
  return domain->solInfo().aeroFlag;
}

void
DistrPodProjectionNonLinDynamic::updatePrescribedDisplacement(DistrModalGeomState *geomState)
{
  if(domain->solInfo().initialTime == 0.0) {
    if(domain->numInitDispModal() == 0 && domain->numInitVelocityModal() == 0) {
      DistrVector q_Big(MDNLDynamic::solVecInfo());
      MDNLDynamic::updatePrescribedDisplacement(geomState_Big);
      geomState_Big->get_tot_displacement(q_Big);
      reduceDisp(q_Big, geomState->q);
      geomState_Big->setVelocity(*v0_Big);
    }
  }
}

void
DistrPodProjectionNonLinDynamic::getConstForce(DistrVector &constantForce)
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
    DistrVector constantForce_Big(MDNLDynamic::solVecInfo());

    MDNLDynamic::getConstForce(constantForce_Big);

    const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
    projectionBasis.reduceAll(constantForce_Big, constantForce);
  }
}

void
DistrPodProjectionNonLinDynamic::getExternalForce(DistrVector &rhs, DistrVector &constantForce, 
                        int tIndex, double t, DistrModalGeomState *geomState, 
                        DistrVector &elementInternalForce, DistrVector &aeroF, double localDelta)
{
  int numNeumanModal = domain->nNeumannModal();
  const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();

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
    DistrVector           rhs_Big(MDNLDynamic::solVecInfo()),
                constantForce_Big(MDNLDynamic::solVecInfo()),
                    aeroForce_Big(MDNLDynamic::solVecInfo());
    constantForce_Big.zero();

    MDNLDynamic::getExternalForce(rhs_Big, constantForce_Big, tIndex, t, geomState_Big, elementInternalForce,
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
DistrPodProjectionNonLinDynamic::getIncDisplacement(DistrModalGeomState *geomState, DistrVector &du, DistrModalGeomState *refState, bool zeroRot)
{
  geomState->get_inc_displacement(du, *refState, zeroRot);
}

double
DistrPodProjectionNonLinDynamic::formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                                                  DistrVector& residual, DistrVector& rhs, DistrModalGeomState *geomState, double localDelta)
{
  double beta, gamma, alphaf, alpham, dt = 2*localDelta;
  getNewmarkParameters(beta, gamma, alphaf, alpham);

  if(domain->solInfo().useMassNormalizedBasis && C == 0) {
    if(domain->solInfo().order == 1)
      rhs = -1.0*inc_displacement;
    else
      rhs = -(1-alpham)/(1-alphaf)*inc_displacement + dt*(1-alpham)*velocity + dt*dt*((1-alpham)/2-beta)*acceleration;
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
    DistrVector inc_displacement_Big(MDNLDynamic::solVecInfo()),
                        velocity_Big(MDNLDynamic::solVecInfo()),
                    acceleration_Big(MDNLDynamic::solVecInfo()),
                        residual_Big(MDNLDynamic::solVecInfo()),
                             rhs_Big(MDNLDynamic::solVecInfo());
    residual_Big.zero();

    const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
    projectionBasis.fullExpand(inc_displacement, inc_displacement_Big);
    projectionBasis.fullExpand(velocity, velocity_Big);
    projectionBasis.fullExpand(acceleration, acceleration_Big);

    MDNLDynamic::formRHScorrector(inc_displacement_Big, velocity_Big, acceleration_Big,
                                  residual_Big, rhs_Big, geomState_Big, localDelta);

    projectionBasis.reduce(rhs_Big, rhs);
  }
  if(domain->solInfo().order == 1)
    rhs += localDelta*residual;
  else
    rhs += (dt*dt*beta)*residual;
  return solver_->getResidualNorm(rhs);
}

void
DistrPodProjectionNonLinDynamic::formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                                                  DistrVector& rhs, DistrModalGeomState &, double mid, double localDelta)
{
   filePrint(stderr,"PodProjectionNonLinDynamic::formRHSpredictor is not implemented");
}

void
DistrPodProjectionNonLinDynamic::formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce,
                                                    DistrModalGeomState &geomState, DistrVector &rhs, DistrModalGeomState *refState)
{
  // XXX For the case of initial values for the rotation, the transformed mass matrix
  //     should be used to solve for the initial accelerations.
  DistrVector     fext_Big(MDNLDynamic::solVecInfo()),
              velocity_Big(MDNLDynamic::solVecInfo()),
                   rhs_Big(MDNLDynamic::solVecInfo());
  fext_Big.zero(); rhs_Big.zero(); velocity_Big.zero();
  
  const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
  projectionBasis.expand(velocity, velocity_Big);
  MDNLDynamic::formRHSinitializer(fext_Big, velocity_Big, elementInternalForce, *geomState_Big, rhs_Big, refState_Big);

  projectionBasis.reduce(rhs_Big, rhs);
  rhs += fext;  
}

bool
DistrPodProjectionNonLinDynamic::factorWhenBuilding() const 
{
  // Delayed factorization
  return false;
}

void
DistrPodProjectionNonLinDynamic::updateStates(DistrModalGeomState *refState, DistrModalGeomState& geomState, double time)
{
  if(domain->solInfo().piecewise_contact) updateCS = true; // update that shit
  if((!domain->solInfo().getNLInfo().linearelastic && (geomState_Big->getHaveRot() || geomState_Big->getTotalNumElemStates() > 0))
     || domain->solInfo().readInROBorModes.size() > 1) {
    // updateStates is called after midpoint update (i.e. once per timestep)
    // so it is a convenient place to update and copy geomState_Big, if necessary
    const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
    if(domain->solInfo().readInROBorModes.size() == 1) {
      // note: for local bases method, geomState_Big has already been updated in setLocalBasis
      DistrVector q_Big(MDNLDynamic::solVecInfo());
      projectionBasis.expand(geomState.q, q_Big);
      geomState_Big->explicitUpdate(decDomain, q_Big);
    } 
    
    if(geomState_Big->getHaveRot()) {
      DistrVector vel_Big(MDNLDynamic::solVecInfo()), acc_Big(MDNLDynamic::solVecInfo());
      projectionBasis.expand(geomState.vel, vel_Big);
      geomState_Big->setVelocity(vel_Big);
      projectionBasis.expand(geomState.acc, acc_Big);
      geomState_Big->setAcceleration(acc_Big);
    }

    if(geomState_Big->getTotalNumElemStates() > 0){
      MDNLDynamic::updateStates(refState_Big, *geomState_Big, time);
    }

    *refState_Big = *geomState_Big;
  }
}

double
DistrPodProjectionNonLinDynamic::getStiffAndForce(DistrModalGeomState& geomState, DistrVector& residual, DistrVector& elementInternalForce,
                                                  double t, DistrModalGeomState *refState, bool forceOnly)
{
  DistrVector r(solVecInfo()); r.zero();
  if(domain->solInfo().getNLInfo().linearelastic == 1 && domain->getFollowedElemList().empty() && !domain->solInfo().elemLumpPodRom) {
    // TODO
    //K_reduced.mult(geomState.q, r);
    //residual -= r;
  }
  else {
    DistrVector        q_Big(MDNLDynamic::solVecInfo()),
                residual_Big(MDNLDynamic::solVecInfo());
    residual_Big.zero();
    const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
                                             
    if(domain->solInfo().readInROBorModes.size() > 1) {
      DistrVector dq = geomState.q;
      dq -= refState->q;
      projectionBasis.expand(dq, q_Big);
      geomState_Big->update(*refState_Big, q_Big, 2);
    }
    else {
      projectionBasis.expand(geomState.q, q_Big);
      geomState_Big->explicitUpdate(decDomain, q_Big);
    }

    MDNLDynamic::getStiffAndForce(*geomState_Big, residual_Big, elementInternalForce, t, refState_Big, forceOnly); 

    projectionBasis.reduce(residual_Big, r);
    residual += r;
  }

  int glNumMPCs = 0; 
  for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
        glNumMPCs += decDomain->getSubDomain(iSub)->numMPCs();
  }
  
  if(structCom)
    glNumMPCs = structCom->globalSum(glNumMPCs);

#ifdef USE_EIGEN3
  if((!domain->solInfo().readInDualROB.empty() && glNumMPCs > 0) || domain->solInfo().modalLMPC) {
    solver_->activateContact();
    if(domain->solInfo().modalLMPC) { // linear
      solver_->updateLMPCs(geomState.q);
    } else { // nonlinear, or mixed
     
      // get primal and dual bases
      const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis(); 
      int startDualCol = 0; int blockDualCols = 0; 
      solver_->getLocalDualBasisInfo(startDualCol, blockDualCols); 

      // allocate space for projection of constraints and their rhs onto dual basis
      DistrVecBasis CtW(blockDualCols, decDomain->masterSolVecInfo());
      Eigen::Matrix<double,Eigen::Dynamic,1> WtRhs(blockDualCols); WtRhs.setZero();
      // get MPCs from each subDomain
      for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
        //decDomain->getSubDomain(iSub)->dualConstraintProjection( dualBasis, CtW, WtRhs); 
        decDomain->getSubDomain(iSub)->dualConstraintProjection(dualMapVectors_, CtW, WtRhs, startDualCol, blockDualCols);
      } 
    
      // unify reduced constraint matrix and right hand side
      if(structCom)
        structCom->globalSum(WtRhs.size(), WtRhs.data());
 
      // project constraints and rhs onto primal basis
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> WtCV(blockDualCols,projectionBasis.numVec()); 
      for(int i = 0; i < blockDualCols; ++i) {
        Eigen::Matrix<double,Eigen::Dynamic,1> buffer(projectionBasis_.numVec()); buffer.setZero();
        projectionBasis_.sparseVecReduce(CtW[i],buffer.data()); // force through specialized tempalate function
        WtCV.row(i).array() = buffer.array();
      }

      // send to solver
      double dt = domain->solInfo().getTimeStep(), beta = domain->solInfo().newmarkBeta;
      double Kcoef = dt*dt*beta;
      solver_->addMPCs(WtCV, WtRhs, Kcoef); 
    }
  }
#endif
  return solver_->getResidualNorm(residual);
}

void
DistrPodProjectionNonLinDynamic::reBuild(DistrModalGeomState& geomState, int iter, double localDelta, double t)
{
  MDNLDynamic::reBuild(*geomState_Big, iter, localDelta, t);
}

void
DistrPodProjectionNonLinDynamic::dynamCommToFluid(DistrModalGeomState* geomState, DistrModalGeomState* bkGeomState,
                                                  DistrVector& velocity, DistrVector& bkVelocity,
                                                  DistrVector& vp, DistrVector& bkVp, int step, int parity,
                                                  int aeroAlg, double time)
{
  // XXX this could be more implemented more efficiently (only called for AERO)
  DistrVector   velocity_Big(MDNLDynamic::solVecInfo()),
              bkVelocity_Big(MDNLDynamic::solVecInfo()),
                      vp_Big(MDNLDynamic::solVecInfo()),
                    bkVp_Big(MDNLDynamic::solVecInfo());
  DistrGeomState *bkGeomState_Big = NULL;
  const GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
  projectionBasis.expand(velocity  , velocity_Big);
  projectionBasis.expand(bkVelocity, bkVelocity_Big);
  projectionBasis.expand(vp        , vp_Big);
  projectionBasis.expand(bkVp      , bkVp_Big);

  MDNLDynamic::dynamCommToFluid(geomState_Big, bkGeomState_Big, velocity_Big, bkVelocity_Big, vp_Big, bkVp_Big, step, parity, aeroAlg, time);
}

void
DistrPodProjectionNonLinDynamic::dynamOutput(DistrModalGeomState* geomState, DistrVector& velocity, DistrVector &vp,
                                             double time, int timestep, DistrVector& force, DistrVector &aeroF, DistrVector &acceleration,
                                             DistrModalGeomState *refState)
{
  MDDynamMat dynOps;
  SysState<DistrVector> systemState(geomState->q, velocity, acceleration, vp); // load data to be output
  podPostPro->dynamOutput(timestep+1, time, dynOps, force, &aeroF, systemState);
}

void
DistrPodProjectionNonLinDynamic::getConstraintMultipliers(DistrModalGeomState &geomState)
{
  MDNLDynamic::getConstraintMultipliers(*geomState_Big);
}

void
DistrPodProjectionNonLinDynamic::initializeParameters(int step, DistrModalGeomState *geomState)
{
  MDNLDynamic::initializeParameters(step, geomState_Big);
}

void
DistrPodProjectionNonLinDynamic::updateParameters(DistrModalGeomState *geomState)
{
  MDNLDynamic::updateParameters(geomState_Big);
}

bool
DistrPodProjectionNonLinDynamic::checkConstraintViolation(double &err, DistrModalGeomState *geomState)
{
  return MDNLDynamic::checkConstraintViolation(err, geomState_Big);
}

double
DistrPodProjectionNonLinDynamic::getResidualNorm(DistrVector &residual, DistrModalGeomState &, double)
{
#ifdef USE_EIGEN3
  if(!domain->solInfo().readInDualROB.empty() || domain->solInfo().modalLMPC) {
    const Eigen::VectorXd &fc = solver_->lastReducedConstraintForce();
    return (Eigen::Map<const Eigen::VectorXd>(residual.data(), residual.size()) - fc).norm();
  } else
#endif
  return solver_->getResidualNorm(residual);
}

void
DistrPodProjectionNonLinDynamic::expandForce(DistrVector &fr, DistrVector &f)
{
  GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
  if(domain->solInfo().useMassNormalizedBasis) {
    // note: this case should not be relied upon for modelIII with reduced mesh, 
    // because the full mass matrix is not available
    DistrVector Vfr(MDNLDynamic::solVecInfo());
    projectionBasis.fullExpand(fr, Vfr);
    M->mult(Vfr, f);
  }   
  else {
    projectionBasis.fullExpand(fr, f);
  }
}

void
DistrPodProjectionNonLinDynamic::reduceDisp(DistrVector &d, DistrVector &dr)
{
  GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
  if(domain->solInfo().useMassNormalizedBasis) {
    DistrVector Md(MDNLDynamic::solVecInfo());
    M->mult(d, Md);
    projectionBasis.reduce(Md, dr, false);
  }
  else {
    projectionBasis.reduce(d, dr, false);
  }
}

DistrModalGeomState*
DistrPodProjectionNonLinDynamic::createGeomState(){

  if(!(geomState_Big = new DistrGeomState(decDomain)))
    throw std::runtime_error("Geomstate didn't bind");

  return new DistrModalGeomState(solVecInfo());
}

DistrModalGeomState*
DistrPodProjectionNonLinDynamic::copyGeomState(DistrModalGeomState* geomState){

 refState_Big = new DistrGeomState(*geomState_Big);

 return new DistrModalGeomState(*geomState);
}

GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *
DistrPodProjectionNonLinDynamic::getSolver(){
  return solver_;
}

const GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *
DistrPodProjectionNonLinDynamic::getSolver() const {
  return solver_;
}

int 
DistrPodProjectionNonLinDynamic::selectLocalBasis(DistrVector &q)
{
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
  } else if(centroids.numVec() > 0) {
    // modelII: slow implementation of using cluster centroids.
    int    cc      = std::numeric_limits<int>::max();
    double minDist = std::numeric_limits<double>::max(); 
    GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
    for(int c=0; c< centroids.numVec(); ++c){
      DistrVector centroidVec = centroids[c]; 
      DistrVector q_Big(MDNLDynamic::solVecInfo());
      projectionBasis.expand(q,q_Big);
      centroidVec -= q_Big; 
      double distNorm = centroidVec.norm();
      if(distNorm < minDist) { // check proximity to centroid, update if closer
        minDist = distNorm;
        cc      = c; 
      } 
    }
    return cc; 
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
DistrPodProjectionNonLinDynamic::initLocalBasis(DistrVector &q0){
#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11))
  if(domain->solInfo().readInROBorModes.size() > 1) {
    int oldLBI = localBasisId;
    localBasisId = selectLocalBasis(q0);
    if(oldLBI != localBasisId) 
      filePrint(stderr," ... Selecting Local Basis # %d     ...\n",localBasisId);
    GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
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
    }
  }
#else
  filePrint(stderr,"*** Error: Local auxiliary bases requres Aero-S to be built with eigen library and c++11 support.\n");
  exit(-1);
#endif
}

void
DistrPodProjectionNonLinDynamic::setLocalBasis(DistrModalGeomState *refState, DistrModalGeomState *geomState, DistrVector &q_n, DistrVector &vel, DistrVector &acc)
{
#ifdef USE_EIGEN3
  // Local bases: this block set new local basis and orthogonaly projects onto new basis
  if(domain->solInfo().readInROBorModes.size() > 1) {

    int j = selectLocalBasis(geomState->q);

    DistrVector dq_Big(MDNLDynamic::solVecInfo());
    GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();

    DistrVector dq = geomState->q;
    dq -= q_n;
    projectionBasis.expand(dq, dq_Big, false); // XXX not using compressed basis
    geomState_Big->update(dq_Big);

    if(j != localBasisId) { // if a new local basis has been selected, update the things
      filePrint(stderr," ... Selecting Local Basis # %d \n     ... ",j);
      DistrVector vel_Big(MDNLDynamic::solVecInfo()),
                  acc_Big(MDNLDynamic::solVecInfo());
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
      }
      if(VtV.size() == 0) {
        reduceDisp(vel_Big, vel);
        reduceDisp(acc_Big, acc);
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
DistrPodProjectionNonLinDynamic::readLocalBasesCent(DistrVecNodeDof6Conversion &converter, DistrNodeDof6Buffer &buffer, std::vector<int> &locBasisVec)
{
  if(domain->solInfo().readInLocalBasesCent.size() >= 1) {
    centroids.dimensionIs(locBasisVec.size(), decDomain->masterSolVecInfo());
    DistrVecBasis::iterator it = centroids.begin();
    for(int rob=0; rob<locBasisVec.size(); ++rob) {
      std::string fileName = domain->solInfo().readInLocalBasesCent[rob].c_str();
      DistrBasisInputFile CentroidFile(fileName);
      if(verboseFlag) filePrint(stderr, " ... Reading centroid from file %s ...\n", fileName.c_str());
      assert(CentroidFile.validCurrentState());
      CentroidFile.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      it++;
    }
  }
}

void
DistrPodProjectionNonLinDynamic::readLocalBasesAuxi()
{
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
        if(verboseFlag) filePrint(stderr, " ... Reading aux data from file %s ...\n", fileName.c_str());
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
DistrPodProjectionNonLinDynamic::projectLocalBases(int i, int j, DistrVector &q)
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

} /* end namespace Rom */
#endif
