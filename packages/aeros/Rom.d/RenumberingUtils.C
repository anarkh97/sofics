#include "RenumberingUtils.h"

#include "ConnectivityUtils.h"

#include <iterator>

#include <Mortar.d/MortarDriver.d/MortarHandler.h>

// Debug
#include <iostream>

namespace Rom {

void
MeshRenumbering::init(const Connectivity &elemToNode, bool verboseFlag) {
  if(verboseFlag) {
    std::cerr << "Elements = ";
    for (std::vector<int>::const_iterator it = reducedElemIds_.begin(); it != reducedElemIds_.end(); ++it) {
      std::cerr << *it + 1 << " ";
    }
    std::cerr << std::endl;
  }
  
  inverse_numbering(reducedElemIds_.begin(), reducedElemIds_.end(), std::inserter(elemRenumbering_, elemRenumbering_.end()));
 
  if(verboseFlag) {
    std::cerr << "Element renumbering = ";
    for (std::map<int, int>::const_iterator it = elemRenumbering_.begin(); it != elemRenumbering_.end(); ++it) {
      std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
    }
    std::cerr << std::endl;
  }
  
  connections(elemToNode, reducedElemIds_.begin(), reducedElemIds_.end(), std::back_inserter(reducedNodeIds_));

  if(verboseFlag) {
    std::cerr << "Nodes = ";
    for (std::vector<int>::const_iterator it = reducedNodeIds_.begin(); it != reducedNodeIds_.end(); ++it) {
      std::cerr << *it + 1 << " ";
    }
    std::cerr << std::endl;
  }

  inverse_numbering(reducedNodeIds_.begin(), reducedNodeIds_.end(), std::inserter(nodeRenumbering_, nodeRenumbering_.end()));

  if(verboseFlag) {
    std::cerr << "Node renumbering = ";
    for (std::map<int, int>::const_iterator it = nodeRenumbering_.begin(); it != nodeRenumbering_.end(); ++it) {
      std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
    }
    std::cerr << std::endl;
  }
}

void
MeshRenumbering::init(const Connectivity &elemToNode, Domain * domain, bool verboseFlag) {
 
  if(verboseFlag) {
    std::cerr << "Elements = ";
    for (std::vector<int>::const_iterator it = reducedElemIds_.begin(); it != reducedElemIds_.end(); ++it) {
      std::cerr << *it + 1 << " ";
    }
    std::cerr << std::endl;
  }

  inverse_numbering(reducedElemIds_.begin(), reducedElemIds_.end(), std::inserter(elemRenumbering_, elemRenumbering_.end()));

  if(verboseFlag) {
    std::cerr << "Element renumbering = ";
    for (std::map<int, int>::const_iterator it = elemRenumbering_.begin(); it != elemRenumbering_.end(); ++it) {
      std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
    }
    std::cerr << std::endl;
  }

  connections(elemToNode, reducedElemIds_.begin(), reducedElemIds_.end(), std::back_inserter(reducedNodeIds_));

/*
  //Extract Nodes belonging to Surface Topologies
  int numMC = domain->GetnMortarConds();
  if(numMC > 0) {
    std::set<int> activeSurfs; 
    //First: loop over Mortar Constraints to identify active surfaces
    for(int mc = 0; mc<numMC; ++mc){
      if(domain->GetMortarCond(mc)->GetInteractionType() == MortarHandler::CTC){
        activeSurfs.insert(domain->GetMortarCond(mc)->GetMasterEntityId()); 
        activeSurfs.insert(domain->GetMortarCond(mc)->GetSlaveEntityId());
      }
    }

    //Second: loop over list of surfaces to extract nodes using connectivity 
    std::set<int> surfNodeSet; 
    for(std::set<int>::iterator it = activeSurfs.begin(); it != activeSurfs.end(); ++it){
      for(int surf = 0; surf < domain->getNumSurfs(); ++surf) { // loop through surfaces and check for matching ID
        if(*it == domain->GetSurfaceEntity(surf)->GetId()) {
          int nNodes = domain->GetSurfaceEntity(surf)->GetnNodes(); 
          for(int nn = 0; nn < nNodes; ++nn) {
            surfNodeSet.insert(domain->GetSurfaceEntity(surf)->GetGlNodeId(nn));
          }
        }
      }
    }
    // copy surface topology nodes to container and removed duplicates
    std::copy(surfNodeSet.begin(),surfNodeSet.end(),std::back_inserter(reducedNodeIds_));
    std::sort(reducedNodeIds_.begin(), reducedNodeIds_.end());
    std::vector<int>::iterator it = std::unique(reducedNodeIds_.begin(), reducedNodeIds_.end());
    reducedNodeIds_.resize(std::distance(reducedNodeIds_.begin(),it));
  }
  */

  if(verboseFlag) {
    std::cerr << "Nodes = ";
    for (std::vector<int>::const_iterator it = reducedNodeIds_.begin(); it != reducedNodeIds_.end(); ++it) {
      std::cerr << *it + 1 << " ";
    }
    std::cerr << std::endl;
  }

  inverse_numbering(reducedNodeIds_.begin(), reducedNodeIds_.end(), std::inserter(nodeRenumbering_, nodeRenumbering_.end()));

  if(verboseFlag) {
    std::cerr << "Node renumbering = ";
    for (std::map<int, int>::const_iterator it = nodeRenumbering_.begin(); it != nodeRenumbering_.end(); ++it) {
      std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
    }
    std::cerr << std::endl;
  }
}

void
SampledMeshRenumbering::init(const Connectivity &nodeToNode,
                             const Connectivity &nodeToElem) {
  std::cerr << "Sample nodes = ";
  for (std::vector<int>::const_iterator it = sampleNodeIds_.begin(); it != sampleNodeIds_.end(); ++it) {
    std::cerr << *it + 1 << " ";
  }
  
  std::cerr << std::endl;
  neighborhood(nodeToNode, sampleNodeIds_.begin(), sampleNodeIds_.end(), std::back_inserter(reducedNodeIds_));

  std::cerr << "Nodes = ";
  for (std::vector<int>::const_iterator it = reducedNodeIds_.begin(); it != reducedNodeIds_.end(); ++it) {
    std::cerr << *it + 1 << " ";
  }
  std::cerr << std::endl;

  inverse_numbering(reducedNodeIds_.begin(), reducedNodeIds_.end(), std::inserter(nodeRenumbering_, nodeRenumbering_.end()));

  renumber(nodeRenumbering_, sampleNodeIds_.begin(), sampleNodeIds_.end(), std::back_inserter(reducedSampleNodeIds_));

  std::cerr << "Node renumbering = ";
  for (std::map<int, int>::const_iterator it = nodeRenumbering_.begin(); it != nodeRenumbering_.end(); ++it) {
    std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
  }
  std::cerr << std::endl;

  connections(nodeToElem, sampleNodeIds_.begin(), sampleNodeIds_.end(), std::back_inserter(reducedElemIds_));

  std::cerr << "Elements = ";
  for (std::vector<int>::const_iterator it = reducedElemIds_.begin(); it != reducedElemIds_.end(); ++it) {
    std::cerr << *it + 1 << " ";
  }
  std::cerr << std::endl;

  inverse_numbering(reducedElemIds_.begin(), reducedElemIds_.end(), std::inserter(elemRenumbering_, elemRenumbering_.end()));
  
  std::cerr << "Element renumbering = ";
  for (std::map<int, int>::const_iterator it = elemRenumbering_.begin(); it != elemRenumbering_.end(); ++it) {
    std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
  }
  std::cerr << std::endl;
}

} /* end namespace Rom */
