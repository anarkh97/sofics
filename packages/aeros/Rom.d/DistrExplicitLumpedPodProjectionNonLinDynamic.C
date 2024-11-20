#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>

#include <Driver.d/GeoSource.h>

#include <utility>
#include <cstddef>

#include <sys/time.h>
#include <Utils.d/DistHelper.h>
#include <Corotational.d/DistrGeomState.h>


extern GeoSource *geoSource;

namespace Rom {

DistrExplicitLumpedPodProjectionNonLinDynamic::DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitPodProjectionNonLinDynamicBase(domain), K(NULL),
  localReducedMeshId_(0)
{}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess() {
  
  DistrExplicitPodProjectionNonLinDynamicBase::preProcess();

  buildPackedElementWeights();

  if(domain->solInfo().stable) {
    execParal(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subInitWeightedStiffOnly);
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n1) {

  setLocalBasis(d_n1,v_n_h); // this also sets the reduced Mesh and resets sparse basis if needed
  DistrVector temp1(solVecInfo());
  temp1 = dt_n_h*v_n_h;
  d_n->zero();
  normalizedBasis_.expand( temp1, *d_n);
  execParal(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subUpdateWeightedNodesOnly,*d_n);
  normalizedBasis_.addLocalPart(temp1,d_n1);

  if(haveRot) {
    v_n->zero();
    normalizedBasis_.expand(v_n_h, *v_n);
    execParal(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subSetVelocityWeightedNodesOnly,*v_n);
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::setLocalReducedMesh(int j) {

  // set new ID number
  localReducedMeshId_ = j; 

  // make new sparse basis
/*  DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
  for(int i=0; i<decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();
  normalizedBasis_.makeSparseBasis(localPackedWeightedNodes_[localReducedMeshId_], all_cdsa);
  delete [] all_cdsa;*/
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {

  execParal(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly,*fInt,t,tIndex);

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    GenMDDynamMat<double> ops;
    ops.K = K;
    decDomain->rebuildOps(ops, 0.0, 0.0, 0.0, kelArray);
  }
  
  if (domain->solInfo().filterFlags) {
    trProject(*fInt);
  }

  *tempVec = *fInt - *fExt;
  f.zero();
  normalizedBasis_.sparseVecReduce(*tempVec, f);
  //  the residual is computed in this step to avoid projecting into the reduced coordinates twice
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  packedElementWeights_.resize(geoSource->elementLumpingWeightSize());
  localPackedWeightedNodes_.resize(geoSource->elementLumpingWeightSize());
  packedWeightedNodes_.resize(decDomain->getNumSub());
  int elemCounter;
  for(int j=0; j<geoSource->elementLumpingWeightSize(); ++j){ // loop over each reduced mesh and load data structures
    packedElementWeights_[j].resize(decDomain->getNumSub());
    localPackedWeightedNodes_[j].resize(decDomain->getNumSub());
    localReducedMeshId_ = j;
    execParal(decDomain->getNumSub(),
              this, &DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights);

    elemCounter = 0;
    for(int i=0; i<decDomain->getNumSub(); ++i) elemCounter += packedElementWeights_[j][i].size();
    if(structCom) elemCounter = structCom->globalSum(elemCounter);

    if(geoSource->elementLumpingWeightSize() == 1)
      filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", elemCounter);
    else
      filePrint(stderr, " ... # Elems. in Reduced Mesh %d = %-4d...\n", j,elemCounter);
  }

  if(elemCounter < domain->numElements() && geoSource->elementLumpingWeightSize() == 1) {
    filePrint(stderr, " ... Compressing Basis              ...\n");
    DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
    for(int i=0; i<decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();
    normalizedBasis_.makeSparseBasis(packedWeightedNodes_, all_cdsa);
    delete [] all_cdsa;
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subUpdateWeightedNodesOnly(int iSub, DistrVector &v) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->update(vec, packedWeightedNodes_[iSub], 2);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subSetVelocityWeightedNodesOnly(int iSub, DistrVector &v) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->setVelocity(packedWeightedNodes_[iSub].size(), packedWeightedNodes_[iSub].data(), vec, 2);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subTransformWeightedNodesOnly(int iSub, DistrVector &v, int type) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->transform(vec, packedWeightedNodes_[iSub], type, true);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    sd->getWeightedStiffAndForceOnly(packedElementWeights_[localReducedMeshId_][iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub], melArray[iSub]); // residual -= internal force);
  }
  else {
    sd->getWeightedInternalForceOnly(packedElementWeights_[localReducedMeshId_][iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub], melArray[iSub]); // residual -= internal force);
  }
  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, -1.0); // f = -residual
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subInitWeightedStiffOnly(int iSub) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(sd->numUncon(), 0.0);
  Vector eIF(sd->maxNumDOF());

  sd->getWeightedStiffAndForceOnly(packedElementWeights_[localReducedMeshId_][iSub], *(*geomState)[iSub], eIF,
                                   allCorot[iSub], kelArray[iSub], residual,
                                   1.0, 0.0, (*geomState)[iSub], melArray[iSub]);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights(int iSub) {
  //need to sort out issue with subdomains for packed node vector
  SubDomain *sd = decDomain->getSubDomain(iSub);
  std::map<int, double> &subElementWeights = packedElementWeights_[localReducedMeshId_][iSub];
  std::vector<int> &subWeightedNodes       = packedWeightedNodes_[iSub];
  std::vector<int> &subLocalWeightedNodes  = localPackedWeightedNodes_[localReducedMeshId_][iSub];
  
  for (GeoSource::ElementWeightMap::const_iterator it = geoSource->elementLumpingWeightBegin(localReducedMeshId_),
                                                   it_end = geoSource->elementLumpingWeightEnd(localReducedMeshId_);
                                                   it != it_end; ++it) {
    const int elemId = it->first;

    const int packedId = sd->glToPackElem(elemId);
    if (packedId < 0) {
      continue;
    }

    const double weight = it->second;
    if (weight != 0.0) {
      Element *ele = sd->getElementSet()[packedId]; // get weighted element data
      std::vector<int> node_buffer(ele->numNodes());
      subElementWeights.insert(subElementWeights.end(), std::make_pair(packedId, weight)); //pack element weight
      //put nodes for weighted element into dummy vector and insert into packed node vector
      ele->nodes(node_buffer.data());
      int numNodes = node_buffer.size();
      subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
      if(geoSource->elementLumpingWeightSize() > 1)
        subLocalWeightedNodes.insert(subLocalWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
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
    for(int surf = 0; surf < numSurf; ++surf) {
      int nNodes = domain->GetSurfaceEntity(surf)->GetnNodes();
      for(int nn = 0; nn < nNodes; ++nn){
        int newNode = domain->GetSurfaceEntity(surf)->GetGlNodeId(nn); //get node from surface
        const int pnId = sd->globalToLocal(newNode); // see if its in this subdomain
        if(pnId > 0)
          subWeightedNodes.push_back(pnId); 
      }
    }
  }

/* XXX consider whether to also add the nodes with non-follower external forces
  for(int i = 0; i < sd->nNeumann(); ++i) {
    subWeightedNodes.push_back(sd->getNBC()[i].nnum);
  }
*/
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

MDDynamMat *
DistrExplicitLumpedPodProjectionNonLinDynamic::buildOps(double mCoef, double cCoef, double kCoef) {

  MDDynamMat *result = DistrExplicitPodProjectionNonLinDynamicBase::buildOps(mCoef, cCoef, kCoef);
  K = result->K;

  return result;
}

} // end namespace Rom
