#ifdef USE_EIGEN3
#include "DistrLumpedPodProjectionNonLinDynamic.h"
#include "PodProjectionSolver.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Timers.d/StaticTimers.h>

#include <fstream>
#include <utility>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

DistrLumpedPodProjectionNonLinDynamic::DistrLumpedPodProjectionNonLinDynamic(Domain *domain) :
  DistrPodProjectionNonLinDynamic(domain),
  localReducedMeshId_(0)
{}

void
DistrLumpedPodProjectionNonLinDynamic::preProcess() {
  DistrPodProjectionNonLinDynamic::preProcess();

  //if(!resetFromClean)
    buildPackedElementWeights();
}

double
DistrLumpedPodProjectionNonLinDynamic::getStiffAndForce(DistrModalGeomState& geomState, DistrVector& residual, DistrVector& elementInternalForce,
                                                  double t, DistrModalGeomState *refState, bool forceOnly)
{
  times->buildStiffAndForce -= getTime();
  DistrVector r(solVecInfo()); r.zero();
  {
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

    // update the contact surfaces
    if(t != domain->solInfo().initialTime) {
      if(domain->GetnContactSurfacePairs()) {
        if(!domain->solInfo().piecewise_contact || updateCS) {
          updateContactSurfaces(*geomState_Big, refState_Big);
          updateCS = false;
        }
        elementInternalForce.resize(MDNLDynamic::elemVecInfo());
        residual_Big.conservativeResize(MDNLDynamic::solVecInfo());
        localTemp->resize(MDNLDynamic::solVecInfo());
      }
      // set the gap for the linear constraints
      if(fetiSolver) decDomain->setConstraintGap(geomState_Big, refState_Big, fetiSolver, t);
    }

    execParal(decDomain->getNumSub(), this, &DistrLumpedPodProjectionNonLinDynamic::subGetStiffAndForce, *geomState_Big, residual_Big, elementInternalForce, t, refState_Big, forceOnly);

    projectionBasis.sparseVecReduce(residual_Big, r);
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
    }
    else { // nonlinear, or mixed

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
  times->buildStiffAndForce += getTime();
  return solver_->getResidualNorm(residual);
}

void 
DistrLumpedPodProjectionNonLinDynamic::subGetStiffAndForce(int iSub,  DistrGeomState &geomState, DistrVector &res,
                                                           DistrVector &elemIntForce, double t, DistrGeomState *refState, bool forceOnly)
{
  std::map<int, double> &subElementWeights = packedElementWeights_[localReducedMeshId_][iSub];
  SubDomain *sd = decDomain->getSubDomain(iSub);
  StackVector residual(res.subData(iSub), res.subLen(iSub));
  // eIF = element internal force
  StackVector eIF(elemIntForce.subData(iSub), elemIntForce.subLen(iSub));
  GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
  if(forceOnly) {
    sd->getWeightedInternalForceOnly(subElementWeights, *geomState[iSub], eIF, allCorot[iSub], kelArray[iSub], residual, 1.0, t, subRefState, melArray[iSub]);
  }
  else {
    sd->getWeightedStiffAndForceOnly(subElementWeights, *geomState[iSub], eIF, allCorot[iSub], kelArray[iSub], residual, 1.0, t, subRefState, melArray[iSub]);
  }
}

void
DistrLumpedPodProjectionNonLinDynamic::updateStates(DistrModalGeomState *refState, DistrModalGeomState& geomState, double time)
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
      DistrVector vel_Big(MDNLDynamic::solVecInfo()),
                  acc_Big(MDNLDynamic::solVecInfo());
      projectionBasis.expand(geomState.vel, vel_Big);
      geomState_Big->setVelocity(vel_Big);
      projectionBasis.expand(geomState.acc, acc_Big);
      geomState_Big->setAcceleration(acc_Big);
    }

    if(geomState_Big->getTotalNumElemStates() > 0)
      execParal(decDomain->getNumSub(), this, &DistrLumpedPodProjectionNonLinDynamic::subUpdateStates, refState_Big, geomState_Big, time);

    *refState_Big = *geomState_Big;
  }
}

void
DistrLumpedPodProjectionNonLinDynamic::subUpdateStates(int iSub, DistrGeomState *refState, DistrGeomState *geomState, double time)
{
  std::map<int, double> &subWeightedElems = packedElementWeights_[localReducedMeshId_][iSub];
  SubDomain *sd = decDomain->getSubDomain(iSub);
  GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
  sd->updateWeightedElemStatesOnly(subWeightedElems, subRefState, *(*geomState)[iSub], allCorot[iSub], time);
}

void
DistrLumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  
  packedElementWeights_.resize(geoSource->elementLumpingWeightSize());    // resize packedElementWeights_ to number of local meshes provided
  packedWeightedElems_.resize(decDomain->getNumSub());
  localPackedWeightedNodes_.resize(geoSource->elementLumpingWeightSize());// resize localPackedWeightedNodes_ to number of local meshes provided
  packedWeightedNodes_.resize(decDomain->getNumSub());
  int elemCounter;
  for (int j=0; j<geoSource->elementLumpingWeightSize(); ++j) {           // loop over each reduced Mesh

    packedElementWeights_[j].resize(decDomain->getNumSub());
    localPackedWeightedNodes_[j].resize(decDomain->getNumSub());
    localReducedMeshId_ = j; 
    execParal(decDomain->getNumSub(),
              this, &DistrLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights);
  
    elemCounter = 0;
    for(int i=0; i<decDomain->getNumSub(); ++i) elemCounter += packedElementWeights_[j][i].size();
    if(structCom) elemCounter = structCom->globalSum(elemCounter);
 
    if(!resetFromClean) {
      if(geoSource->elementLumpingWeightSize() == 1)
        filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", elemCounter);
      else
        filePrint(stderr, " ... # Elems. in Reduced Mesh %d = %-4d...\n", j,elemCounter);
    } else {
      if(verboseFlag) filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", elemCounter);
    }
  }

  int numElems = 0;
  for(int ns = 0; ns < packedWeightedElems_.size(); ++ns)
    numElems += packedWeightedElems_[ns].size(); 

  if(structCom) numElems = structCom->globalSum(numElems);

  if(geoSource->elementLumpingWeightSize() > 1 && !resetFromClean)
    filePrint(stderr, " ... # of unique Elems in Mesh = %-4d ... \n",numElems);

  if(geoSource->elementLumpingWeightSize() == 1 && numElems < domain->numElements()) {
    if(domain->solInfo().useMassNormalizedBasis) { // don't compress if using Local mesh
      if(!resetFromClean) filePrint(stderr, " ...       Compressing Basis        ...\n");
      DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
      for(int i=0; i<decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();
      GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
      projectionBasis.makeSparseBasis(packedWeightedNodes_, all_cdsa);
      delete [] all_cdsa;
    }
  }
  else {
    if(!domain->solInfo().useMassNormalizedBasis && !resetFromClean) {
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
        solver_->storeReducedMass(VtMV);
      }
      else
#endif
        filePrint(stderr, " *** WARNING: \"use_mass_normalized_basis off\" is not supported for\n"
                          "     for model III when \"samplmsh.elementmesh.inc\" file is used   \n"
                          "     unless a modal DIMASS file is specified containing the reduced \n"
                          "     mass matrix.\n");
    }
  }
}

void
DistrLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights(int iSub)
{
  SubDomain *sd = decDomain->getSubDomain(iSub);
  std::map<int, double> &subElementWeights = packedElementWeights_[localReducedMeshId_][iSub];
  std::set<int>         &subWeightedElems  = packedWeightedElems_[iSub];
  std::vector<int>      &subWeightedNodes  = packedWeightedNodes_[iSub];
  std::vector<int> &subLocalWeightedNodes  = localPackedWeightedNodes_[localReducedMeshId_][iSub];

  // clear lists to account for shrinking element set
  subElementWeights.clear();
  subWeightedElems.clear();
  subWeightedNodes.clear();
  subLocalWeightedNodes.clear();

  //fprintf(stderr," num MPCs %d, num MPCs_primal %d \n", sd->numMPCs(), sd->numMPCs_primal());
  for (GeoSource::ElementWeightMap::const_iterator it     = geoSource->elementLumpingWeightBegin(localReducedMeshId_),
                                                   it_end = geoSource->elementLumpingWeightEnd(localReducedMeshId_);
                                                   it    != it_end; ++it) {

    const int elemId = it->first;

    const int packedId = sd->glToPackElem(elemId);
    if (packedId < 0) {
      continue;
    }

    const double weight = it->second;
    if (weight > 1e-16) {
      Element *ele = sd->getElementSet()[packedId]; // get weighted element data
      std::vector<int> node_buffer(ele->numNodes());
      subElementWeights.insert(subElementWeights.end(), std::make_pair(packedId, weight)); //pack element weight
      //put nodes for weighted element into dummy vector and insert into packed node vector
      ele->nodes(node_buffer.data());
      subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
      if(geoSource->elementLumpingWeightSize() > 1)
        subLocalWeightedNodes.insert(subLocalWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
      subWeightedElems.insert(packedId);
    }
  }

  if(!domain->solInfo().reduceFollower) {
    // add the nodes for all elements to which follower forces have been applied
    std::vector<int> &followedElemList = sd->getFollowedElemList();
    for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) {
      Element *ele = sd->getElementSet()[*it];
      std::vector<int> node_buffer(ele->numNodes());
      ele->nodes(node_buffer.data());
      subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
      if(geoSource->elementLumpingWeightSize() > 1)
        subLocalWeightedNodes.insert(subLocalWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
    }
  }

  int numSurf = domain->getNumSurfs();
  if(numSurf > 0) { // add nodes associated with surface topologies
    for(int surf = 0; surf < numSurf; ++surf) { // loop through each surface
      int nNodes = domain->GetSurfaceEntity(surf)->GetnNodes(); // get number of nodes in that surface
      for(int nn = 0; nn < nNodes; ++nn){ // loop throught the nodes
        int newNode = domain->GetSurfaceEntity(surf)->GetGlNodeId(nn); //get node from surface
        const int pnId = sd->globalToLocal(newNode); // see if its in this subdomain
        if(pnId > 0)
          subWeightedNodes.push_back(pnId);
      }
    }
  }

  //sort nodes in ascending order and erase redundant nodes
  std::sort(subWeightedNodes.begin(), subWeightedNodes.end());
  std::vector<int>::iterator packedNodeIt = std::unique(subWeightedNodes.begin(),subWeightedNodes.end());
  subWeightedNodes.resize(packedNodeIt-subWeightedNodes.begin());
  if(geoSource->elementLumpingWeightSize() > 1) {
    std::sort(subLocalWeightedNodes.begin(), subLocalWeightedNodes.end());
    std::vector<int>::iterator packedNodeIt = std::unique(subLocalWeightedNodes.begin(),subLocalWeightedNodes.end());
    subLocalWeightedNodes.resize(packedNodeIt-subLocalWeightedNodes.begin());
  }
}

void
DistrLumpedPodProjectionNonLinDynamic::setLocalReducedMesh(int j)
{

  // first zero out the old element stiffness arrays
  execParal(decDomain->getNumSub(), this, &DistrLumpedPodProjectionNonLinDynamic::subZeroStiffMats); 

  // then set new ID number
  localReducedMeshId_ = std::min(geoSource->elementLumpingWeightSize()-1, j);

  DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
  for(int i=0; i<decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();

  // make new sparse basis
  if(verboseFlag) filePrint(stderr," ... Compressing Local Basis       ... \n");
  GenVecBasis<double,GenDistrVector> &projectionBasis = solver_->projectionBasis();
  projectionBasis.makeSparseBasis(localPackedWeightedNodes_[localReducedMeshId_], all_cdsa); // these could be computed once, stored and then switched between
  delete [] all_cdsa;
}

void
DistrLumpedPodProjectionNonLinDynamic::subZeroStiffMats(int iSub)
{
  std::map<int, double> &subElementWeights = packedElementWeights_[localReducedMeshId_][iSub];

  // zero out current element stiffness matrices
  for(std::map<int, double>::const_iterator it = subElementWeights.begin(),
                                        it_end = subElementWeights.end(); it != it_end; ++it) {
    const int iele = it->first;
    kelArray[iSub][iele].zero();
  }
}

} /* end namespace Rom */
#endif
