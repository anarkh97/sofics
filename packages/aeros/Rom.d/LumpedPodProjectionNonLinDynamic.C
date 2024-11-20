#include "LumpedPodProjectionNonLinDynamic.h"
#include "PodProjectionSolver.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

#include <fstream>
#include <utility>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

LumpedPodProjectionNonLinDynamic::LumpedPodProjectionNonLinDynamic(Domain *domain) :
  PodProjectionNonLinDynamic(domain),
  localReducedMeshId_(0)
{}

void
LumpedPodProjectionNonLinDynamic::preProcess() {
  PodProjectionNonLinDynamic::preProcess();

  buildPackedElementWeights();
}

void
LumpedPodProjectionNonLinDynamic::getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                                             Corotator **allCorot, FullSquareMatrix *kelArray,
                                                             Vector &residual, double lambda, double time, GeomState *refState,
                                                             FullSquareMatrix *_melArray, bool forceOnly) {
  FullSquareMatrix *melArray = (domain->solInfo().quasistatic) ? (FullSquareMatrix*) NULL : _melArray;
  if(forceOnly) {
    domain->getWeightedInternalForceOnly(packedElementWeights_[localReducedMeshId_],
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray);
  }
  else {
    domain->getWeightedStiffAndForceOnly(packedElementWeights_[localReducedMeshId_],
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray);
  }
}

void
LumpedPodProjectionNonLinDynamic::updateStates(ModalGeomState *refState, ModalGeomState& geomState, double time)
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
      Vector vel_Big(NonLinDynamic::solVecInfo()),
             acc_Big(NonLinDynamic::solVecInfo());
      projectionBasis.expand(geomState.vel, vel_Big);
      geomState_Big->setVelocity(vel_Big);
      projectionBasis.expand(geomState.acc, acc_Big);
      geomState_Big->setAcceleration(acc_Big);
    }

    if(geomState_Big->getTotalNumElemStates() > 0)
      domain->updateWeightedElemStatesOnly(packedWeightedElems_, refState_Big, *geomState_Big, allCorot, time);

    *refState_Big = *geomState_Big;
  }
}

void
LumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  packedElementWeights_.resize(geoSource->elementLumpingWeightSize());    // resize packedElementWeights_ to number of local meshes provided
  localPackedWeightedNodes_.resize(geoSource->elementLumpingWeightSize());// resize localPackedWeightedNodes_ to number of local meshes provided
  for (int j=0; j<geoSource->elementLumpingWeightSize(); ++j) {           // loop over each reduced Mesh
    for (GeoSource::ElementWeightMap::const_iterator it = geoSource->elementLumpingWeightBegin(j),
                                                     it_end = geoSource->elementLumpingWeightEnd(j);
         it != it_end; ++it) {
      const int elemId = it->first;

      const int packedId = geoSource->glToPackElem(elemId);
      if (packedId < 0) {
        continue;
      }

      const double weight = it->second;
      if (weight != 0.0) {
        Element *ele = domain->getElementSet()[packedId]; // get weighted element data
        std::vector<int> node_buffer(ele->numNodes());
        packedElementWeights_[j].insert(packedElementWeights_[j].end(), std::make_pair(packedId, weight));
        //put nodes for weighted element into dummy vector and insert into packed node vector
        ele->nodes(node_buffer.data());
        packedWeightedNodes_.insert(packedWeightedNodes_.end(), node_buffer.begin(), node_buffer.end());
        if(geoSource->elementLumpingWeightSize() > 1) {
          localPackedWeightedNodes_[j].insert(localPackedWeightedNodes_[j].end(), node_buffer.begin(), node_buffer.end());
        }
        packedWeightedElems_.insert(packedId);
      }
    }
    filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", packedElementWeights_[j].size());
  }

  // XXX also need to add to packedWeightedNodes the nodes of any elements to which a follower force has been applied
  // if the follower forces are not reduced. See: DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights

  //sort nodes in ascending order and erase redundant nodes
  std::sort(packedWeightedNodes_.begin(), packedWeightedNodes_.end());
  std::vector<int>::iterator packedNodeIt = std::unique(packedWeightedNodes_.begin(), packedWeightedNodes_.end());
  packedWeightedNodes_.resize(packedNodeIt-packedWeightedNodes_.begin());

  if(geoSource->elementLumpingWeightSize() == 1 && packedWeightedElems_.size() < domain->numElements()) {
    if(domain->solInfo().useMassNormalizedBasis) { // don't compress yet if using Local mesh
      filePrint(stderr, " ... Compressing Basis              ...\n");
      GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
      projectionBasis.makeSparseBasis(packedWeightedNodes_, domain->getCDSA());
    }
  }
  else {
    if(geoSource->elementLumpingWeightSize() > 1) {
      GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
      for (int j=0; j<geoSource->elementLumpingWeightSize(); ++j) {    
        std::sort(localPackedWeightedNodes_[j].begin(), localPackedWeightedNodes_[j].end());
        std::vector<int>::iterator packedNodeIt = std::unique(localPackedWeightedNodes_[j].begin(), localPackedWeightedNodes_[j].end());
        localPackedWeightedNodes_[j].resize(packedNodeIt-localPackedWeightedNodes_[j].begin());
        filePrint(stderr, " ... # Nodes in Reduced Mesh = %-5d...\n", localPackedWeightedNodes_[j].size());
      }
    }
    if(!domain->solInfo().useMassNormalizedBasis && !domain->solInfo().modalDIMASS) {
      filePrint(stderr, " *** WARNING: \"use_mass_normalized_basis off\" is not supported for\n"
                        "     for model III when \"samplmsh.elementmesh.inc\" file is used   \n"
                        "     unless a modal DIMASS file is specified containing the reduced \n"
                        "     mass matrix.\n");
    }
  }
}

void
LumpedPodProjectionNonLinDynamic::setLocalReducedMesh(int j)
{
  // zero out current element stiffness matrices
  for(std::map<int, double>::const_iterator it = packedElementWeights_[localReducedMeshId_].begin(),
                                        it_end = packedElementWeights_[localReducedMeshId_].end(); it != it_end; ++it) {
    const int iele = it->first;
    kelArray[iele].zero();
  }

  // set new ID number
  localReducedMeshId_ = std::min(geoSource->elementLumpingWeightSize()-1, j);

  // make new sparse basis
  GenVecBasis<double> &projectionBasis = solver_->projectionBasis();
  projectionBasis.makeSparseBasis(localPackedWeightedNodes_[j], domain->getCDSA()); // these could be computed once, stored and then switched between
}

} /* end namespace Rom */
