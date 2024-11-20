#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <iostream>

#include <Driver.d/Domain.h>
#include <Math.d/VectorSet.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/SubDomain.h>
#include <Corotational.d/DistrGeomState.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Utils.d/pstress.h>
#include <Utils.d/BinFileHandler.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/CornerMaker.h>
#include <Feti.d/Feti.h>
#include <Driver.d/Mpc.h>
#include <Solvers.d/DiagParallelSolver.h>
#include <Solvers.d/GmresSolver.h>
#include <Paral.d/MDDynam.h>
#include <Paral.d/GenMS.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Paral.d/DomainGroupTask.h>
#include <Solvers.d/MultiDomainRbm.h>
#include <Solvers.d/Mumps.h>
#include <Driver.d/SysState.h>
#include <Rom.d/BlockCyclicMap.h>
#include <Rom.d/EiGalerkinProjectionSolver.h>
#ifdef USE_MPI
#include <Comm.d/Communicator.h>
#endif
extern Connectivity *procMpcToMpc;

#include <Driver.d/SubDomainFactory.h>
#include "Feti.d/FetiSub.h"

template<class Scalar>
GenDecDomain<Scalar>::GenDecDomain(Domain *d, Communicator *structCom, bool soweredInput, bool coarseLevel)
 : mt(d->getTimers()), soweredInput(soweredInput), coarseLevel(coarseLevel)
{
  domain = d;
  initialize(); 
#ifdef USE_MPI
  communicator = new FSCommunicator(structCom);
#else
  communicator = new FSCommunicator;
#endif
  myCPU = communicator->cpuNum();
}

template<class Scalar>
void 
GenDecDomain<Scalar>::initialize()
{
  mpcToSub_primal = 0;
  primalFile = 0;
  stress     = 0;
  weight     = 0;
  globalStress = 0;
  globalWeight = 0;
  grToSub    = 0;
  wetInterfaceNodes = 0;
  numWetInterfaceNodes = 0;
  outFreqCount = 0;
  outEigCount = 0;
  glSubToLocal = 0;
  communicator = 0;
  cpuToCPU = 0; 
  numPrimalMpc = 0;
  numDualMpc = 0;
  firstOutput = true;
  internalInfo = 0;
  internalInfo2 = 0;
  masterSolVecInfo_ = 0;
  nodeInfo = 0;
  nodeVecInfo = 0;
  eleVecInfo = 0;
  bcVecInfo = 0;
  wiPat = 0;
  ba = 0;
  ba2 = 0;
} 

template<class Scalar>
GenDecDomain<Scalar>::~GenDecDomain()
{
	if(mpcToSub_primal) { delete mpcToSub_primal; mpcToSub_primal = 0; }
	for(auto sub :subDomain)
		delete sub;

	if(primalFile) { delete primalFile; primalFile = 0; }
	if(stress) { delete stress; stress = 0; }
	if(weight) { delete weight; weight = 0; }
	if(globalStress) { delete [] globalStress; globalStress = 0;}
	if(globalWeight) {delete [] globalWeight; globalWeight = 0;}

	if(wetInterfaceNodes) { delete [] wetInterfaceNodes; }
	if(communicator) { delete communicator; communicator = 0; }
	if(glSubToLocal) { delete [] glSubToLocal; glSubToLocal = 0; }
	if(internalInfo) delete internalInfo;
	if(internalInfo2) delete internalInfo2;
	if(nodeInfo) delete nodeInfo;
	if(nodeVecInfo) delete nodeVecInfo;
	if(eleVecInfo) delete eleVecInfo;
	if(bcVecInfo) delete bcVecInfo;
	if(masterSolVecInfo_) delete masterSolVecInfo_;
	if(ba) delete ba;
	if(ba2) delete ba2;
	for(std::vector<DistrInfo*>::iterator it = vecInfoStore.begin(); it != vecInfoStore.end(); ++it) delete *it;
	geoSource->deleteMatchArrays(numSub);
}

template<class Scalar>
void GenDecDomain<Scalar>::addBMPCs()
{
  // step 1. create all the candidate bmpcs 
  std::vector<LMPCons *> bmpcs; // candidate bmpcs
  int nbmpc = 0; // number of bmpcs
  int nsn; // number of subdomains sharing node
  for(int i=0; i<nodeToSub->csize(); ++i) {
    if((nsn = nodeToSub->num(i)) < 2) continue; 
    for(int dof=0; dof<DofSet::max_known_nonL_dof; ++dof) { // x,y,z translations and rotations, temperature, helmholtz and intpress
      for(int j=0; j<nsn; ++j) {
        for(int k=j+1; k<nsn; ++k) {
          LMPCons *bmpc = new LMPCons(nbmpc, 0.0, new LMPCTerm(i,dof,0.0)); // dummy mpc, set coef later in distributeMPCs
          bmpc->type = 2;
          int subA = (*nodeToSub)[i][j];
          int subB = (*nodeToSub)[i][k];
          bmpc->psub = (subA > subB) ? subA : subB;
          bmpc->nsub = (subA > subB) ? subB : subA;
          bmpcs.push_back(bmpc);
          nbmpc++;
        }
      }
    }
  }

  // step 2. check which of the candidates are active (ie both subdomains share an active "r" dof)
  int *pstatus = new int[2*nbmpc];
  int *nstatus = pstatus+nbmpc;
  for(int i=0; i<nbmpc; ++i) pstatus[i] = nstatus[i] = 0;
  paralApply(subDomain, &GenSubDomain<Scalar>::bmpcQualify, &bmpcs, pstatus, nstatus);
#ifdef DISTRIBUTED
  communicator->globalSum(2*nbmpc, pstatus);
#endif

  // step 3. now add the active bmpcs to the global mpc collection
  // (check for and eliminate redundant bmpcs)
  nbmpc = 0; // reset counter
  int total = 0;
  for(int i=0; i<nodeToSub->csize(); ++i) {
    if((nsn = nodeToSub->num(i)) < 2) continue;
    for(int dof=0; dof<DofSet::max_known_nonL_dof; ++dof) { // x,y,z translations and rotations, temperature, helmholtz and intpress
      bool got_enough = false;
      for(int j=0; j<nsn; ++j) {
        int count = 0;
        for(int k=j+1; k<nsn; ++k) {
          //if(pstatus[nbmpc] && nstatus[nbmpc] && (*nodeToSub)[j][0] >= (*nodeToSub)[k][0]) {
          if(pstatus[nbmpc] && nstatus[nbmpc] && !got_enough) {
            bmpcs[nbmpc]->lmpcnum = 0; // shouldn't be used 
            domain->addLMPC(bmpcs[nbmpc],false);
            count++; total++;
          }
          nbmpc++;
        }
        if(count > 0) { 
          got_enough = true;
        }
      }
    }
  }
}

template<class Scalar>
template<class ConnectivityType1, class ConnectivityType2>
void GenDecDomain<Scalar>::buildSharedNodeComm(const ConnectivityType1 *nodeToSub, const ConnectivityType2 *subToNode) const
{
  // Start timer
  startTimerMemory(mt.makeConnectivity, mt.memoryConnect);

  map<int, int> NESubMap;
  int numNESub = 0; 
  for(int iSub = 0; iSub < subToNode->csize(); iSub++)
    if(subToNode->num(iSub))
      NESubMap[iSub] = numNESub++;

  // PJSA 9-1-04 for coupled_dph all the wet interface nodes need to be included in the sharedNodes list in SComm
  bool coupled_dph = false;
  std::vector<int> wetInterfaceNodeMap;
  if(domain->solInfo().isCoupled) {
    coupled_dph = true;
    //wetInterfaceNodes = domain->getAllWetInterfaceNodes(numWetInterfaceNodes);
    wetInterfaceNodeMap.assign(nodeToSub->csize(), -1);
    for(int i=0; i<numWetInterfaceNodes; ++i) wetInterfaceNodeMap[wetInterfaceNodes[i]] = i;
  }

  // create each subdomain's interface lists
  int iSub, jSub, subJ, iNode;
  int totConnect;
  std::vector<int> nConnect(numNESub, 0);
  std::vector<int> flag(numNESub, -1);

  // Count connectivity
  totConnect = 0;
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) { // loop over the nodes
      typename ConnectivityType2::IndexType thisNode = (*subToNode)[iSub][iNode];
      bool isWetInterfaceNode = (coupled_dph && (wetInterfaceNodeMap[thisNode] != -1));
      for(jSub = 0; jSub < nodeToSub->num(thisNode); ++jSub) {
        // loop over the subdomains connected to this node
        subJ = (*nodeToSub)[thisNode][jSub];
        bool neighbWithSelf =
		        (subJ == iSub) && isWetInterfaceNode && (domain->solInfo().solvercntl->fetiInfo.fsi_corner == 0);
        if((subJ > iSub) || neighbWithSelf) {
          // only deal with connection to highered numbered subdomains, guarantees symmetry of lists
	  if(!subToNode->num(subJ)) {
	    fprintf(stderr, "reference to empty subdomain 1\n");
	    exit(1);
	  }
          if(flag[NESubMap[subJ]] != iSub) {
            flag[NESubMap[subJ]] = iSub;
            nConnect[NESubMap[iSub]]++;
            if(!neighbWithSelf) {
              nConnect[NESubMap[subJ]]++;
              totConnect += 2;
            }
            else totConnect += 1;
          }
        }
      }
    }
  }

  int **nodeCount = new int*[numNESub];
	std::vector<std::vector<gl_sub_idx>> connectedDomain(numNESub);
	std::vector<std::vector<lc_sub_idx>> remoteID(numNESub);

  // Allocate memory for list of connected subdomains
  // (de-allocated in ~SComm)
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    if(subToNode->num(iSub)) {
      int size = nConnect[NESubMap[iSub]];
      connectedDomain[NESubMap[iSub]].resize(size);
      remoteID[NESubMap[iSub]].resize(size);
      nodeCount[NESubMap[iSub]] = new int[size];
      flag[NESubMap[iSub]] = -1;
      nConnect[NESubMap[iSub]] = 0;
    }
  }

	std::vector<int> whichLocal(numNESub);
	std::vector<int> whichRemote(numNESub);

  for(iSub=0; iSub < subToNode->csize(); ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
      typename ConnectivityType2::IndexType nd = (*subToNode)[iSub][iNode];
      bool isWetInterfaceNode = coupled_dph && (wetInterfaceNodeMap[nd] != -1);
      for(jSub = 0; jSub < nodeToSub->num(nd); ++jSub) {
        subJ = (*nodeToSub)[nd][jSub];
        bool neighbWithSelf =
            (subJ == iSub) && isWetInterfaceNode && (domain->solInfo().solvercntl->fetiInfo.fsi_corner == 0);
        if((subJ > iSub) || neighbWithSelf) {
          if(flag[NESubMap[subJ]] != iSub) { // attribute location for this sub
            flag[NESubMap[subJ]] = iSub;
            connectedDomain[NESubMap[subJ]][nConnect[NESubMap[subJ]]] = iSub;
            connectedDomain[NESubMap[iSub]][nConnect[NESubMap[iSub]]] = subJ;
            remoteID[NESubMap[subJ]][nConnect[NESubMap[subJ]]] = nConnect[NESubMap[iSub]];
            remoteID[NESubMap[iSub]][nConnect[NESubMap[iSub]]] = nConnect[NESubMap[subJ]];
            if(!neighbWithSelf) whichLocal[NESubMap[subJ]] = nConnect[NESubMap[iSub]]++;
            else whichLocal[NESubMap[subJ]] = nConnect[NESubMap[iSub]];
            whichRemote[NESubMap[subJ]] = nConnect[NESubMap[subJ]]++;
            nodeCount[NESubMap[iSub]][whichLocal[NESubMap[subJ]]] = 1;
            nodeCount[NESubMap[subJ]][whichRemote[NESubMap[subJ]]] = 1;
          }
          else {
            nodeCount[NESubMap[iSub]][whichLocal[NESubMap[subJ]]]++;
            if(!neighbWithSelf) nodeCount[NESubMap[subJ]][whichRemote[NESubMap[subJ]]]++;
          }
        }
      }
    }
  }

  // allocate memory for interface node lists
  std::vector<std::unique_ptr<Connectivity>> interfNode(numNESub);
  for(iSub=0; iSub < subToNode->csize(); ++iSub)
    if(subToNode->num(iSub))
      interfNode[NESubMap[iSub]] = std::make_unique<Connectivity>(nConnect[NESubMap[iSub]], nodeCount[NESubMap[iSub]]);

  // fill the lists
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    if(subToNode->num(iSub)) {
      flag[NESubMap[iSub]]     = -1;
      nConnect[NESubMap[iSub]] =  0;
    }
  }

  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
      auto nd = (*subToNode)[iSub][iNode];
      bool isWetInterfaceNode = coupled_dph && (wetInterfaceNodeMap[nd] != -1);
      for(jSub = 0; jSub < nodeToSub->num(nd); ++jSub) {
        subJ = (*nodeToSub)[nd][jSub];
        bool neighbWithSelf =
		        (subJ == iSub) && isWetInterfaceNode && (domain->solInfo().solvercntl->fetiInfo.fsi_corner == 0);
        if((subJ > iSub) || neighbWithSelf) {
          if(flag[NESubMap[subJ]] != iSub) { // attribute location for this sub
            flag[NESubMap[subJ]] = iSub;
            if(!neighbWithSelf) whichLocal[NESubMap[subJ]] = nConnect[NESubMap[iSub]]++;
            else whichLocal[NESubMap[subJ]] = nConnect[NESubMap[iSub]];
            whichRemote[NESubMap[subJ]] = nConnect[NESubMap[subJ]]++;
            (*interfNode[NESubMap[iSub]])[whichLocal[NESubMap[subJ]]][0] = (glSubToLocal[iSub] >= 0) ? subDomain[glSubToLocal[iSub]]->globalToLocal(nd) : nd;
            (*interfNode[NESubMap[subJ]])[whichRemote[NESubMap[subJ]]][0] = (glSubToLocal[subJ] >= 0) ? subDomain[glSubToLocal[subJ]]->globalToLocal(nd) : nd;

            nodeCount[NESubMap[iSub]][whichLocal[NESubMap[subJ]]]=1;
            nodeCount[NESubMap[subJ]][whichRemote[NESubMap[subJ]]]=1;
          }
          else {
            int il = nodeCount[NESubMap[iSub]][whichLocal[NESubMap[subJ]]]++;
            (*interfNode[NESubMap[iSub]])[whichLocal[NESubMap[subJ]]][il] = (glSubToLocal[iSub] >= 0) ? subDomain[glSubToLocal[iSub]]->globalToLocal(nd) : nd;

            if(!neighbWithSelf) {
              int jl = nodeCount[NESubMap[subJ]][whichRemote[NESubMap[subJ]]]++;
              (*interfNode[NESubMap[subJ]])[whichRemote[NESubMap[subJ]]][jl] = (glSubToLocal[subJ] >= 0) ? subDomain[glSubToLocal[subJ]]->globalToLocal(nd) : nd;
            }
          }
        }
      }
    }
  }

  for(iSub = 0; iSub < numNESub; ++iSub)
    delete [] nodeCount[iSub];
  delete [] nodeCount;

  for(iSub = 0; iSub < numSub; ++iSub) {
	  gl_sub_idx subI = (localSubToGl.size() > 0) ? localSubToGl[iSub] : iSub;
    if(!subToNode->num(subI)) {
      fprintf(stderr, "reference to empty subdomain 2\n");
      exit(1);
    }
    SComm *sc = new SComm(nConnect[NESubMap[subI]], std::move(connectedDomain[NESubMap[subI]]),
                          std::move(remoteID[NESubMap[subI]]), std::move(interfNode[NESubMap[subI]]));
    sc->locSubNum = iSub;
    subDomain[iSub]->setSComm(sc);
  }

  stopTimerMemory(mt.makeConnectivity, mt.memoryConnect);
}

template<class Scalar>
void
GenDecDomain<Scalar>::preProcessBCsEtc()
{
  if(!soweredInput) {
    distributeBCs();
    distributeControlLawData();
    distributeDiscreteMass();
    paralApply(subDomain, &GenSubDomain<Scalar>::renumberBCsEtc); 
    paralApply(subDomain, &GenSubDomain<Scalar>::renumberControlLaw); 
  }
  else {
    if (geoSource->binaryInputControlLeft) {
      distributeControlLawData();
      paralApply(subDomain, &GenSubDomain<Scalar>::renumberControlLaw);
    }
  }

  paralApply(subDomain, &BaseSub::makeCDSA); 
  paralApply(subDomain, &Domain::makeNodeToNode_sommer);
}

template<class Scalar>
void
GenDecDomain<Scalar>::preProcessMPCs()
{
#ifdef SOWER_SURFS
	if(soweredInput) {
		//HB compute mortar LMPCs
		domain->SetMortarPairing();
		domain->SetUpSurfaces(); //domain->SetUpSurfaces(geoSource->sower_nodes);
		if(domain->solInfo().newmarkBeta != 0.0) { // not for explicit dynamics
			domain->ComputeMortarLMPC();
			domain->computeMatchingWetInterfaceLMPC();
			domain->CreateMortarToMPC();
		}
#ifdef MORTAR_DEBUG
		domain->PrintSurfaceEntities();
    domain->PrintMortarConds();
#endif
	}
#endif
	if(domain->solInfo().solvercntl->fetiInfo.bmpc ||
	   (!isFeti(domain->solInfo().solvercntl->type) &&
	    (domain->solInfo().filterFlags || domain->solInfo().rbmflg))) addBMPCs();
	if(domain->getNumLMPC() > 0) {
		if(verboseFlag) filePrint(stderr, " ... Applying the MPCs              ...\n");
		// check for mpcs involving bad nodes and constrained DOFs
		if(nodeToSub) domain->checkLMPCs(nodeToSub.get());
		// select which mpcs are to be included in the coarse problem
		domain->setPrimalLMPCs(numDualMpc, numPrimalMpc);
		// distribute mpcs
		execParal(numSub, this, &GenDecDomain<Scalar>::extractSubDomainMPCs);
		makeMpcToSub(); // new version works for distributed data
		// renumber local mpcs
		paralApplyToAll(numSub, subDomain, &GenSubDomain<Scalar>::renumberMPCs);
		// locate and store dsa, c_dsa and cc_dsa numbering for future reference
		paralApply(subDomain, &GenSubDomain<Scalar>::locateMpcDofs);
		// make mpcToMpc connectivity
		makeMpcToMpc(); // this should be in the Feti-DP constructor
	}
	else domain->solInfo().getFetiInfo().mpcflag = 0;
}

template<class Scalar>
void
GenDecDomain<Scalar>::deleteMPCs()
{
  paralApply(subDomain, &GenSubDomain<Scalar>::deleteMPCs);
  mpcToSub_dual.reset();
  mpcToMpc.reset();
  mpcToCpu.reset();
  numDualMpc = 0;
}

template<class Scalar>
void
GenDecDomain<Scalar>::reProcessMPCs()
{
  deleteMPCs();
  preProcessMPCs();
  getSharedMPCs();
  paralApply(subDomain, &BaseSub::mergeInterfaces);
  paralApply(subDomain, &GenSubDomain<Scalar>::applyMpcSplitting);
}

template<class Scalar>
void
GenDecDomain<Scalar>::extractSubDomainMPCs(int iSub)
{
  if(numDualMpc) subDomain[iSub]->extractMPCs(domain->numLMPC, domain->lmpc);
  if(numPrimalMpc) subDomain[iSub]->extractMPCs_primal(domain->numLMPC, domain->lmpc);
}

/** \brief Build connectivity data for the DecDomain.
 * \details Builds:
 *  - elemToNode
 *  - subToNode
 *  - nodeToSub
 *  - subToSub
 * @tparam Scalar
 */
template<class Scalar>
void
GenDecDomain<Scalar>::makeSubToSubEtc()
{
	if(soweredInput) {
		geoSource->setNumNodes(communicator->globalMax(geoSource->nodeToSub_sparse->csize()));
	}
	else {
		mt.memoryElemToNode -= memoryUsed();
		if(!elemToNode) elemToNode = std::make_unique<Connectivity>(domain->packedEset.asSet());
		mt.memoryElemToNode += memoryUsed();

		mt.memorySubToNode -= memoryUsed();
		subToNode = std::make_unique<Connectivity>( subToElem->transcon(*elemToNode) );
		if(!domain->GetnContactSurfacePairs()) subToNode->sortTargets();
		mt.memorySubToNode += memoryUsed();

		mt.memoryNodeToSub -= memoryUsed();
		nodeToSub = std::make_unique<Connectivity>( subToNode->reverse() );
		mt.memoryNodeToSub += memoryUsed();
		domain->setNumNodes(nodeToSub->csize());

		mt.memorySubToNode -= memoryUsed();
		subToSub = std::make_unique<Connectivity>( subToNode->transcon(*nodeToSub) );
		mt.memorySubToNode += memoryUsed();

#ifdef USE_MPI
		if(domain->numSSN() || domain->solInfo().isCoupled
		   || domain->solInfo().solvercntl->type == SolverSelection::Direct
		   || domain->solInfo().aeroFlag > -1 || geoSource->binaryOutput == 0) {
#else
		  // TODO Find out why this would be different from the MPI case.
			if(domain->numSSN() || domain->solInfo().isCoupled
			|| domain->solInfo().solvercntl->type == SolverSelection::Direct || domain->solInfo().aeroFlag > -1) {
#endif
			// sommerfeld, scatter, wet, distributed neum PJSA 6/28/2010 multidomain mumps PJSA 12/01/2010 non-binary output for mpi
			mt.memoryNodeToElem -= memoryUsed();
			if(!domain->nodeToElem) domain->nodeToElem = elemToNode->alloc_reverse();
			mt.memoryNodeToElem += memoryUsed();

			mt.memoryElemToSub -= memoryUsed();
			elemToSub = std::make_unique<Connectivity>( subToElem->reverse() );
			mt.memoryElemToSub += memoryUsed();
			if(domain->numSSN() > 0) domain->checkSommerTypeBC(domain, elemToNode.get(), domain->nodeToElem); // flip normals if necessary
		}
	}
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSubDomains() 
{
  makeSubDMaps();
  if(!soweredInput && geoSource->getGlob()) geoSource->computeClusterInfo(localSubToGl[0], subToNode.get());
  subDomain.resize(numSub);

  startTimerMemory(mt.makeSubDomains, mt.memorySubdomain);
  if(soweredInput) {
    subDomain = geoSource->template readDistributedInputFiles<Scalar>(localSubToGl);

    ControlLawInfo *claw = geoSource->getControlLaw();
    if(claw) {
#ifdef DISTRIBUTED
      claw->numActuator = communicator->globalSum(claw->numActuator);
      claw->numSensor = communicator->globalSum(claw->numSensor);
      claw->numUserForce = communicator->globalSum(claw->numUserForce);
      claw->numUserDisp = communicator->globalSum(claw->numUserDisp);
#endif
      domain->setClaw(claw);
    }
#ifdef SOWER_SURFS
    geoSource->readDistributedSurfs(localSubToGl[0]); //pass dummy sub number
#endif
  }
  else {
    execParal(numSub, this, &GenDecDomain<Scalar>::constructSubDomains);
    if(domain->solInfo().isCoupled) {
      for(int iSub = 0; iSub < numSub; ++iSub) subDomain[iSub]->setnodeToSubConnectivity(nodeToSub.get());
      addFsiElements();
    }
	  threadManager->callParal(subDomain.size(), [this](int iSub) {
	  	subDomain[iSub]->renumberElements();
	  });
  }

  paralApply(subDomain, &BaseSub::makeDSA);
  stopTimerMemory(mt.makeSubDomains, mt.memorySubdomain);
}

template<class Scalar>
void
GenDecDomain<Scalar>::distributeDiscreteMass()
{
  // Distribute masses to the subdomains
  DMassData *cmass = domain->firstDiMass;
  while(cmass != 0) {
   int node = cmass->node;
   for(int iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
     int tSub = glSubToLocal[(*nodeToSub)[node][iSub]];
     if(tSub >= 0)
        subDomain[tSub]->addDMass(node, cmass->dof, cmass->diMass);
   }
   cmass = cmass->next;
  }
}

template<class Scalar>
GenParallelSolver<Scalar> *
GenDecDomain<Scalar>::getFetiSolver(GenDomainGroupTask<Scalar> &dgt)
{
	FetiInfo *finfo = &domain->solInfo().getFetiInfo();
	int verboseFlag = domain->solInfo().solvercntl->verbose;
	if(finfo->version == FetiInfo::fetidp) {
		bool rbmFlag = ((domain->solInfo().isStatic() || domain->probType() == SolverInfo::Modal) && !geoSource->isShifted());
		bool geometricRbms = (domain->solInfo().rbmflg && !domain->solInfo().isNonLin());
		std::vector<std::unique_ptr<GenSolver<Scalar>>> dynMats;
		dynMats.reserve(numSub);
		for(int i = 0; i < numSub; ++i)
			dynMats.emplace_back(dgt.dynMats[i]);
		return new GenFetiDPSolver<Scalar>(numSub, globalNumSub,
			{subDomain.begin(), subDomain.end()},
		                                   subToSub.get(), finfo, communicator, glSubToLocal,
		                                   mpcToSub_dual.get(),
		                                   mpcToSub_primal,
		                                   mpcToMpc.get(),
		                                   mpcToCpu.get(), cpuToSub.get(), grToSub.get(),
		                                   std::move(dynMats), dgt.spMats, dgt.rbms, rbmFlag, geometricRbms, verboseFlag);
	}
	else {
		return new GenFetiSolver<Scalar>(numSub, subDomain.data(), subToSub.get(), finfo, communicator,
		                                 glSubToLocal, mpcToSub_dual.get(), cpuToSub.get(),
		                                 dgt.dynMats, dgt.spMats, dgt.rbms, verboseFlag);
	}
}

template<class Scalar>
DiagParallelSolver<Scalar> *
GenDecDomain<Scalar>::getDiagSolver(int nSub, GenSubDomain<Scalar> **sd,
                                    GenSolver<Scalar> **sol)
{
 return new DiagParallelSolver<Scalar>(nSub, sd, sol, cpuToSub.get(), communicator);
}

template<class Scalar>
void
GenDecDomain<Scalar>::getCPUMap()
{
  mt.memoryCPUMAP -= memoryUsed();
#ifdef DISTRIBUTED
  const char *mapName = geoSource->getCpuMapFile();
  FILE *f = fopen(mapName,"r");
  numCPU = geoSource->getCPUMap(f, nullptr, globalNumSub, structCom->numCPUs());
  cpuToCPU = geoSource->getCpuTOCPU();
  if(f) fclose(f);
#else
  numCPU = 1; myCPU = 0;
  geoSource->createSingleCpuToSub(globalNumSub);
#endif
  cpuToSub = geoSource->getCpuToSub();
  mt.memoryCPUMAP += memoryUsed();
}

template<class Scalar>
void
GenDecDomain<Scalar>::setCPUMap(std::unique_ptr<Connectivity> _cpuToSub)
{
  cpuToSub = std::move(_cpuToSub);
  numCPU = cpuToSub->csize();
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSubDMaps()
{
  glSubToLocal = new int[globalNumSub];
  if(cpuToSub) {
    numSub = cpuToSub->num(myCPU);
    for(int iSub = 0; iSub < globalNumSub; ++iSub)
      glSubToLocal[iSub] = -1;
    auto glSubIndices = (*cpuToSub)[myCPU];
    for(int iSub = 0; iSub < numSub; ++iSub) 
      glSubToLocal[ glSubIndices[iSub] ] = iSub;

    localSubToGl = { glSubIndices.begin(), glSubIndices.end() };
  }
  else { // shared memory, all subds on single cpu
    for(int i=0; i<globalNumSub; ++i) glSubToLocal[i] = i;
    numSub = globalNumSub;
    localSubToGl.resize(numSub);
    for(int i=0; i<numSub; ++i) localSubToGl[i] = i;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::preProcess()
{
	if(domain->solInfo().solvercntl->type == SolverSelection::Direct) { // Makes global renumbering, connectivities and dofsets
		domain->preProcessing();

		int numdof = domain->numdof();
		std::vector<int> bc(numdof);
		std::vector<double> bcx(numdof);
		domain->make_bc(bc, bcx);

		domain->make_constrainedDSA();
		domain->makeAllDOFs();
	}
	//soweredInput = geoSource->binaryInput;

	if(!subToElem) {
		if(verboseFlag) filePrint(stderr, " ... Reading Decomposition File     ...\n");
		if(soweredInput)
			geoSource->getBinaryDecomp();
		else
			subToElem = geoSource->getDecomposition();
		//subToElem->sortTargets(); // JAT 021915 // PJSA 11-16-2006
	}

	makeSubToSubEtc();

	if(subToSub)
		globalNumSub = subToSub->csize(); // JAT 021915 // JAT 021916

	if(soweredInput)
		globalNumSub = communicator->globalMax(geoSource->subToNode_sparse->csize());

	if(!cpuToSub) getCPUMap();

	if(verboseFlag) filePrint(stderr, " ... Making the Subdomains          ...\n");
	makeSubDomains();

	preProcessBCsEtc();

	preProcessFSIs();// FLuid-Structure Interaction

	if(soweredInput)
		buildSharedNodeComm(geoSource->nodeToSub_sparse, geoSource->subToNode_sparse);
	else
		buildSharedNodeComm(nodeToSub.get(), subToNode.get());

	makeCorners();// Corners for FETI-DP

	getSharedDOFs();

	preProcessMPCs();//Multi-Point Constraint

	getSharedFSIs();

	getSharedMPCs();

	paralApply(subDomain, &BaseSub::mergeInterfaces);
	paralApply(subDomain, &GenSubDomain<Scalar>::applySplitting);

	//paralApply(subDomain, &GenSubDomain<Scalar>::initSrc);
	makeInternalInfo();

	makeNodeInfo();

#ifdef DISTRIBUTED
	if( ! coarseLevel ) {
		geoSource->setNumNodalOutput();
		if (geoSource->getNumNodalOutput()) {
			for (int i = 0; i < numSub; ++i)
				geoSource->distributeOutputNodesX(subDomain[i], (nodeToSub_copy) ? nodeToSub_copy.get()
				                                                                 : nodeToSub.get()); // make sure each node always gets
			// assigned to the same subdomain.
		}
	}
#endif

	// compute the number of unconstrained dofs for timing file and screen output
	GenDistrVector<int> toto(masterSolVecInfo());
	toto = 1;
	domain->setNumDofs(toto.sqNorm()+domain->nDirichlet()+domain->nCDirichlet());

	// free up some memory
	if(domain->solInfo().solvercntl->type != SolverSelection::Direct && domain->solInfo().aeroFlag < 0) {
		elemToSub.reset();
		//if(!geoSource->elemOutput() && elemToNode) { delete elemToNode; elemToNode = 0; }
	}
}

template<class Scalar>
void
GenDecDomain<Scalar>::scaleDisp(GenDistrVector<Scalar> &u)
{
  execParal(numSub, this, &GenDecDomain<Scalar>::scaleSubDisp, u);
}

template<class Scalar>
void
GenDecDomain<Scalar>::scaleSubDisp(int iSub, GenDistrVector<Scalar> &u)
{
  subDomain[iSub]->scaleDisp(u.subData(iSub));
}


template<class Scalar>
void
GenDecDomain<Scalar>::scaleInvDisp(GenDistrVector<Scalar> &u)
{
  execParal(numSub, this, &GenDecDomain<Scalar>::scaleInvSubDisp, u);
}


template<class Scalar>
void
GenDecDomain<Scalar>::scaleInvSubDisp(int iSub, GenDistrVector<Scalar> &u)
{
  subDomain[iSub]->scaleInvDisp(u.subData(iSub));
}


template<class Scalar>
void
GenDecDomain<Scalar>::scaleDisp(GenDistrVector<Scalar> &u, double alpha)
{
  execParal(numSub, this, &GenDecDomain<Scalar>::scaleSubDisp_2, u, alpha);
}

template<class Scalar>
void
GenDecDomain<Scalar>::scaleSubDisp_2(int iSub, GenDistrVector<Scalar> &u, double alpha) const
{
  subDomain[iSub]->scaleDisp(u.subData(iSub), alpha);
}


template<class Scalar>
void
GenDecDomain<Scalar>::forceContinuity(GenDistrVector<Scalar> &u)
{

  int numNodes = geoSource->numNode();

  Scalar (*mergedDis)[11] = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
  for(int i = 0; i < numNodes; ++i)
    for(int j=0; j<11; j++) mergedDis[i][j] = 0.0;//DofSet::max_known_nonL_dof

  for(int iSub = 0; iSub < numSub; ++iSub)
    subDomain[iSub]->mergeAllDisp(mergedDis, u.subData(iSub));

  for(int iSub = 0; iSub < numSub; ++iSub)
    subDomain[iSub]->forceContinuity(u.subData(iSub),mergedDis);

}


template<class Scalar>
void
GenDecDomain<Scalar>::forceAssemble(GenDistrVector<Scalar> &u)
{
 ba2->assemble(u,0);
}


template<class Scalar>
void
GenDecDomain<Scalar>::postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f,
                                     double eigV, GenDistrVector<Scalar> *aeroF, int x, 
                                     GenMDDynamMat<Scalar> *dynOps, SysState<GenDistrVector<Scalar> > *distState, int ndflag)
{
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // check if there are any output files which need to be processed now
  if(geoSource->noOutput(x) && x != domain->solInfo().initialTimeIndex) return;

  if(verboseFlag && x == 0 && ndflag == 0 && !(domain->solInfo().isDynam() || domain->solInfo().timeIntegration == 1))
    filePrint(stderr," ... Postprocessing                 ...\n");

  Scalar *globVal = 0;  
  if(domain->outFlag && domain->nodeTable == 0) domain->makeNodeTable(domain->outFlag);
  int numNodes = (domain->outFlag) ? domain->exactNumNodes : geoSource->numNode(); 
  int i, j, iSub, inode;

  // initialize and merge displacements from subdomains into global array
  Scalar (*glMergedDis)[11] = new Scalar[numNodes][11];
  Scalar (*locMergedDis)[11] = (domain->solInfo().basicDofCoords) ? 0 : new Scalar[numNodes][11];
  for(i = 0; i < numNodes; ++i)
    for(j=0; j<11; j++) glMergedDis[i][j] = 0.0;
  for(iSub = 0; iSub < numSub; ++iSub)
    subDomain[iSub]->mergeAllDisp(glMergedDis, u.subData(iSub), locMergedDis);

  // intialize and merge aeroelastic forces from subdomains into global array
  Scalar (*mergedAeroF)[6] = 0;
  if(domain->solInfo().aeroFlag > -1 && aeroF) {
    mergedAeroF = new Scalar[numNodes][6];
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<6; ++j) mergedAeroF[i][j] = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub) 
      subDomain[iSub]->mergeForces(mergedAeroF, aeroF->subData(iSub));
  }

  // initialize and merge velocities and accelerations from subdomains into global array
  GenDistrVector<Scalar> *v_n = 0, *a_n = 0;
  Scalar (*glMergedVel)[11] = 0, (*glMergedAcc)[11] = 0;
  Scalar (*locMergedVel)[11] = 0, (*locMergedAcc)[11] = 0;
  if(distState) {
    v_n = &distState->getVeloc();
    a_n = &distState->getAccel();
    glMergedVel = new Scalar[numNodes][11];
    glMergedAcc = new Scalar[numNodes][11];
    if(!domain->solInfo().basicDofCoords) {
      locMergedVel = new Scalar[numNodes][11];
      locMergedAcc = new Scalar[numNodes][11];
    }
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<11; ++j) glMergedVel[i][j] = glMergedAcc[i][j] = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->mergeAllVeloc(glMergedVel, v_n->subData(iSub), locMergedVel);
      subDomain[iSub]->mergeAllAccel(glMergedAcc, a_n->subData(iSub), locMergedAcc);
    }
  }

  // compute current time (or frequency in the case of a helmholtz problem)
  double time;
  if(geoSource->isShifted() && domain->probType() != SolverInfo::Modal      
                            && domain->probType() != SolverInfo::Dynamic) {
    time = domain->getFrequencyOrWavenumber();
    if(domain->solInfo().doFreqSweep) x = outFreqCount++;
  } else if(domain->probType() == SolverInfo::Modal) {
    time = eigV;
    if(domain->solInfo().doEigSweep) x = outEigCount++; 
  }
  else time = eigV;
  if (domain->solInfo().loadcases.size() > 0 && !domain->solInfo().doFreqSweep) time = domain->solInfo().loadcases.front();

  // open output files
  if(x == domain->solInfo().initialTimeIndex && firstOutput) geoSource->openOutputFiles();

  Scalar Wext = 0, Waero = 0, Wela = 0, Wkin = 0, Wdmp = 0;

  for(i = 0; i < numOutInfo; i++) { 
    if(oinfo[i].ndtype != ndflag) continue;
    if(ndflag !=0 && oinfo[i].type != OutputInfo::Disp6DOF && oinfo[i].type !=  OutputInfo::Displacement) continue;
    if(oinfo[i].interval != 0 && x % oinfo[i].interval == 0) {

      Scalar (*mergedDis)[11] = (oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords) ? glMergedDis : locMergedDis;
      Scalar (*mergedVel)[11] = (oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords) ? glMergedVel : locMergedVel;
      Scalar (*mergedAcc)[11] = (oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords) ? glMergedAcc : locMergedAcc;

      switch(oinfo[i].type) {
        case OutputInfo::EigenPair:
        case OutputInfo::FreqRespModes:
        case OutputInfo::Displacement:
          getPrimalVector(i, mergedDis, numNodes, 3, time);
          break;
        case OutputInfo::Velocity:
          if(distState) getPrimalVector(i, mergedVel, numNodes, 3, time);
          break;
        case OutputInfo::Acceleration:
          if(distState) getPrimalVector(i, mergedAcc, numNodes, 3, time);
          break;
        case OutputInfo::EigenPair6:
        case OutputInfo::Disp6DOF:
          getPrimalVector(i, mergedDis, numNodes, 6, time);
          break;
        case OutputInfo::Velocity6:
          if(distState) getPrimalVector(i, mergedVel, numNodes, 6, time);
          break;
        case OutputInfo::Accel6:
          if(distState) getPrimalVector(i, mergedAcc, numNodes, 6, time);
          break;
        case OutputInfo::Temperature:
          getPrimalScalar(i, mergedDis, numNodes, 6, time);
          break;
        case OutputInfo::TemperatureFirstTimeDerivative:
          if(distState) getPrimalScalar(i, mergedVel, numNodes, 6, time);
          break;
        case OutputInfo::AcousticPressure:
        case OutputInfo::EigenPressure:
        case OutputInfo::HelmholtzModes:
        case OutputInfo::Helmholtz:
        case OutputInfo::EigenSlosh:
          getPrimalScalar(i, mergedDis, numNodes, 7, time);
          break;
        case OutputInfo::PressureFirstTimeDerivative:
          if(distState) getPrimalScalar(i, mergedVel, numNodes, 7, time);
          break;
        case OutputInfo::PressureSecondTimeDerivative:
          if(distState) getPrimalScalar(i, mergedAcc, numNodes, 7, time);
          break;
        case OutputInfo::StressXX:
          getStressStrain(u, i, SXX, time);
          break;
        case OutputInfo::StressYY:
          getStressStrain(u, i, SYY, time);
          break;
        case OutputInfo::StressZZ:
          getStressStrain(u, i, SZZ, time);
          break;
        case OutputInfo::StressXY:
          getStressStrain(u, i, SXY, time);
          break; 
        case OutputInfo::StressYZ:
          getStressStrain(u, i, SYZ, time);
          break;
        case OutputInfo::StressXZ:
          getStressStrain(u, i, SXZ, time);
          break;
        case OutputInfo::StrainXX:
          getStressStrain(u, i, EXX, time);
          break;
        case OutputInfo::StrainYY:
          getStressStrain(u, i, EYY, time);
          break;
        case OutputInfo::StrainZZ:
          getStressStrain(u, i, EZZ, time);
          break;
        case OutputInfo::StrainXY:
          getStressStrain(u, i, EXY, time);
          break;
        case OutputInfo::StrainYZ:
          getStressStrain(u, i, EYZ, time);
          break;
        case OutputInfo::StrainXZ:
          getStressStrain(u, i, EXZ, time);
          break;
        case OutputInfo::StressVM:
          getStressStrain(u, i, VON, time);
          break;
        case OutputInfo::StrainVM:
          getStressStrain(u,i,STRAINVON, time);
          break;
        case OutputInfo::ContactPressure: {
          if(!domain->tdenforceFlag())
            getStressStrain(u, i, CONPRESS, time);
          else 
            filePrint(stderr," *** WARNING: Output case %d not supported \n", i);
        } break;
        case OutputInfo::Damage:
          getStressStrain(u, i, DAMAGE, time);
          break;
        case OutputInfo::EquivalentPlasticStrain:
          getStressStrain(u, i, EQPLSTRN, time);
          break;
        case OutputInfo::StressPR1:
          getPrincipalStress(u, i, PSTRESS1, time);
          break;
        case OutputInfo::StressPR2:
          getPrincipalStress(u, i, PSTRESS2, time);
          break;
        case OutputInfo::StressPR3:
          getPrincipalStress(u, i, PSTRESS3, time);
          break;
        case OutputInfo::StrainPR1:
          getPrincipalStress(u, i, PSTRAIN1, time);
          break;
        case OutputInfo::StrainPR2:
          getPrincipalStress(u, i, PSTRAIN2, time);
          break;
        case OutputInfo::StrainPR3:
          getPrincipalStress(u, i, PSTRAIN3, time);
          break;
        case OutputInfo::StressPR1Direc:
          getPrincipalStress(u, i, PSTRESS1DIREC, time);
          break;
        case OutputInfo::StressPR2Direc:
          getPrincipalStress(u, i, PSTRESS2DIREC, time);
          break;
        case OutputInfo::StressPR3Direc:
          getPrincipalStress(u, i, PSTRESS3DIREC, time);
          break;
        case OutputInfo::StrainPR1Direc:
          getPrincipalStress(u, i, PSTRAIN1DIREC, time);
          break;
        case OutputInfo::StrainPR2Direc:
          getPrincipalStress(u, i, PSTRAIN2DIREC, time);
          break;
        case OutputInfo::StrainPR3Direc:
          getPrincipalStress(u, i, PSTRAIN3DIREC, time);
          break;
        case OutputInfo::InXForce:
          getElementForce(u, i, INX, time);
          break;
        case OutputInfo::InYForce:
          getElementForce(u, i, INY, time);
          break;
        case OutputInfo::InZForce:
          getElementForce(u, i, INZ, time);
          break;
        case OutputInfo::AXMoment:
          getElementForce(u, i, AXM, time);
          break;
        case OutputInfo::AYMoment:
          getElementForce(u, i, AYM, time);
          break;
        case OutputInfo::AZMoment:
          getElementForce(u, i, AZM, time);
          break;
        case OutputInfo::DispX:
          getPrimalScalar(i, mergedDis, numNodes, 0, time);
          break;
        case OutputInfo::DispY:
          getPrimalScalar(i, mergedDis, numNodes, 1, time);
          break;
        case OutputInfo::DispZ:
          getPrimalScalar(i, mergedDis, numNodes, 2, time);
          break;
        case OutputInfo::RotX:
          getPrimalScalar(i, mergedDis, numNodes, 3, time);
          break;
        case OutputInfo::RotY:
          getPrimalScalar(i, mergedDis, numNodes, 4, time);
          break;
        case OutputInfo::RotZ:
          getPrimalScalar(i, mergedDis, numNodes, 5, time);
          break;
        case OutputInfo::DispMod:
          if(oinfo[i].nodeNumber == -1) {
            if(!globVal) globVal = new Scalar[numNodes];
            for(inode=0; inode<numNodes; ++inode) {
              globVal[inode] = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                                 mergedDis[inode][1]*mergedDis[inode][1] +
                                                 mergedDis[inode][2]*mergedDis[inode][2]);
            }
            geoSource->outputNodeScalars(i, globVal, numNodes, time);
          }
          else {
            inode = oinfo[i].nodeNumber;
            Scalar dm = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                          mergedDis[inode][1]*mergedDis[inode][1] +
                                          mergedDis[inode][2]*mergedDis[inode][2]);
            geoSource->outputNodeScalars(i, &dm, 1, time);
          }
          break;
	case OutputInfo::RotMod:
          if(oinfo[i].nodeNumber == -1) {
            if(!globVal) globVal = new Scalar[numNodes];
            for(inode=0; inode<numNodes; ++inode) {
              globVal[inode] = ScalarTypes::sqrt(mergedDis[inode][3]*mergedDis[inode][3] +
                                                 mergedDis[inode][4]*mergedDis[inode][4] +
                                                 mergedDis[inode][5]*mergedDis[inode][5]);
            }
            geoSource->outputNodeScalars(i, globVal, numNodes, time);
          }
          else {
            inode = oinfo[i].nodeNumber;
            Scalar rm = ScalarTypes::sqrt(mergedDis[inode][3]*mergedDis[inode][3] +
                                          mergedDis[inode][4]*mergedDis[inode][4] +
                                          mergedDis[inode][5]*mergedDis[inode][5]);
            geoSource->outputNodeScalars(i, &rm, 1, time);
          }
          break;
        case OutputInfo::TotMod:
          if(oinfo[i].nodeNumber == -1) {
            if(!globVal) globVal = new Scalar[numNodes];
            for(inode=0; inode<numNodes; ++inode) {
              globVal[inode] = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                                 mergedDis[inode][1]*mergedDis[inode][1] +
                                                 mergedDis[inode][2]*mergedDis[inode][2] +
                                                 mergedDis[inode][3]*mergedDis[inode][3] +
                                                 mergedDis[inode][4]*mergedDis[inode][4] +
                                                 mergedDis[inode][5]*mergedDis[inode][5]);
            }
            geoSource->outputNodeScalars(i, globVal, numNodes, time);
          }
          else {
            inode = oinfo[i].nodeNumber;
            Scalar tm = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                          mergedDis[inode][1]*mergedDis[inode][1] +
                                          mergedDis[inode][2]*mergedDis[inode][2] +
                                          mergedDis[inode][3]*mergedDis[inode][3] +
                                          mergedDis[inode][4]*mergedDis[inode][4] +
                                          mergedDis[inode][5]*mergedDis[inode][5]);
            geoSource->outputNodeScalars(i, &tm, 1, time);
          }
          break;
/*
        case OutputInfo::Rigid: {
#ifdef DISTRIBUTED
          int glNumRBM = (fetiSolver) ? fetiSolver->numRBM() : 0;
          if (glNumRBM) {
            int rbmSize = glNumRBM*subDomain[0]->numUncon();
            Scalar *localRBMs = new Scalar[rbmSize];
            fetiSolver->getRBMs(localRBMs);
            VectorSet globalRBMs(glNumRBM,numNodes*6, 0.0);
            subDomain[0]->expandRBM(localRBMs,globalRBMs);
            int iRBM;
            for (iRBM=0; iRBM<glNumRBM; ++iRBM)
              communicator->globalSum(numNodes*6,globalRBMs[iRBM].data());

            for (iRBM=0; iRBM<glNumRBM; ++iRBM) {
              filePrint(oinfo[i].filptr,"%e\n",0.0);
              for (iNode=0; iNode<numNodes; ++iNode)
                filePrint(oinfo[i].filptr,"%f %f %f\n",
                   globalRBMs[iRBM][6*iNode+0], globalRBMs[iRBM][6*iNode+1],
                                                globalRBMs[iRBM][6*iNode+2]);
            }
          }
#endif
          }
          break;
*/
        case OutputInfo::Energies:
          getEnergies(u, f, i, time, distState, dynOps, aeroF);
          break;
        case OutputInfo::Farfield: 
          domain->nffp = oinfo[i].interval;
          buildFFP(u,oinfo[i].filptr, true);
          break;
        case OutputInfo::Kirchhoff:
          domain->nffp = oinfo[i].interval;
          buildFFP(u,oinfo[i].filptr, false);
          break;
        case OutputInfo::AeroForce: break; // this is done in DistFlExchange.C
        case OutputInfo::AeroXForce:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 0, time);
          break;
        case OutputInfo::AeroYForce:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 1, time);
          break;
        case OutputInfo::AeroZForce:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 2, time);
          break;
        case OutputInfo::AeroXMom:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 3, time);
          break;
        case OutputInfo::AeroYMom:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 4, time);
          break;
        case OutputInfo::AeroZMom:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 5, time);
          break;
/* TODO
        case OutputInfo::Reactions:
          break;
        case OutputInfo::Reactions6:
          break;
*/
        case OutputInfo::YModulus:
          this->getElementAttr(i,YOUNG, time);
          break;
        case OutputInfo::MDensity:
          this->getElementAttr(i,MDENS, time);
          break;
        case OutputInfo::Thicknes:
          this->getElementAttr(i,THICK, time);
          break;
        case OutputInfo::TDEnforcement: {
          if(domain->tdenforceFlag()) {
            double *plot_data = new double[numNodes]; 
            if(oinfo[i].tdenforc_var == 1) // CONFACE
              for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.5;
            else
              for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.0;
            for(int iMortar=0; iMortar<domain->GetnMortarConds(); iMortar++) {
              domain->GetMortarCond(iMortar)->get_plot_variable(oinfo[i].tdenforc_var,plot_data);
            }
            if(oinfo[i].nodeNumber == -1) 
              geoSource->outputNodeScalars(i, plot_data, numNodes, time);
            else 
              geoSource->outputNodeScalars(i, &plot_data[oinfo[i].nodeNumber], 1, time);
            delete [] plot_data;
          }
          else filePrint(stderr," *** WARNING: Output case %d not supported \n", i);
        } break;
        case OutputInfo::Statevector:
        case OutputInfo::Velocvector:
        case OutputInfo::Accelvector:
        case OutputInfo::InternalStateVar:
        case OutputInfo::DualStateVar:
        case OutputInfo::Forcevector:
        case OutputInfo::Constraintvector:
        case OutputInfo::Constraintviolation:
        case OutputInfo::Residual:
        case OutputInfo::Jacobian:
        case OutputInfo::RobData:
        case OutputInfo::SampleMesh:
        case OutputInfo::ModalDsp:
        case OutputInfo::ModalExF:
        case OutputInfo::ModalMass:
        case OutputInfo::ModalStiffness:
        case OutputInfo::ModalDamping:
        case OutputInfo::ModalDynamicMatrix:
        case OutputInfo::ModalMatrices:
          break;
        default:
          filePrint(stderr," *** WARNING: Output case %d not implemented \n", i);
          break;
      }
    }
    if(globVal) { delete [] globVal; globVal = 0; }
  }
  firstOutput = false; 

  if(glMergedDis) delete [] glMergedDis;
  if(locMergedDis) delete [] locMergedDis;
  if(aeroF) delete [] mergedAeroF;
  if(distState) { 
    delete [] glMergedVel;
    delete [] glMergedAcc; 
    if(locMergedVel) delete [] locMergedVel;
    if(locMergedAcc) delete [] locMergedAcc; 
  }
  if(globVal) delete [] globVal; 

}

template<class Scalar>
void
GenDecDomain<Scalar>::getPrimalVector(int fileNumber, Scalar (*xyz)[11], int numNodes,
                                      int ndof, double time)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];

  int inode;
  if (ndof == 6) {
    if (oinfo.nodeNumber == -1) {
      geoSource->outputNodeVectors6(fileNumber, xyz, numNodes, time);
    }
    else  {
      inode = (domain->outFlag == 1) ? domain->nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
      geoSource->outputNodeVectors6(fileNumber, xyz+inode, 1, time);
    }
  }
  else {
    if (oinfo.nodeNumber == -1) {
      geoSource->outputNodeVectors(fileNumber, xyz, numNodes, time);
    }
    else  {
      inode = (domain->outFlag == 1) ? domain->nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
      geoSource->outputNodeVectors(fileNumber, xyz+inode, 1, time);
    }
  }

}

template<class Scalar>
void
GenDecDomain<Scalar>::getPrimalScalar(int fileNumber, Scalar (*xyz)[11], int numNodes, 
                                      int dof, double time)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber]; 

  int inode;
  if(oinfo.nodeNumber == -1) {
    Scalar *xyz_dof = new Scalar[numNodes];
    for(inode = 0; inode < numNodes; ++inode) xyz_dof[inode] = xyz[inode][dof];
    geoSource->outputNodeScalars(fileNumber, xyz_dof, numNodes, time);
  }
  else {
    inode = (domain->outFlag == 1) ? domain->nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
    geoSource->outputNodeScalars(fileNumber, xyz[inode]+dof, 1, time);
  }
} 

template<class Scalar>
void
GenDecDomain<Scalar>::getAeroForceScalar(int fileNumber, Scalar (*mergedAeroF)[6],
                                         int numNodes, int dof, double time)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];

  int inode;
  if(oinfo.nodeNumber == -1) {
    geoSource->outputNodeScalars(fileNumber, (double*)0, 0, time); // output time
    for(inode = 0; inode < numNodes; ++inode)
      geoSource->outputNodeScalars(fileNumber, mergedAeroF[inode]+dof, 1);
  }
  else {
    inode = (domain->outFlag == 1) ? domain->nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
    geoSource->outputNodeScalars(fileNumber, mergedAeroF[inode]+dof, 1, time);
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::computeSubdElemForce(int iSub, Scalar *globForce,
                                          GenDistrVector<Scalar> *u, int fileNumber, int Findex)
{
  Scalar *locForce = new Scalar[subDomain[iSub]->countElemNodes()];
  subDomain[iSub]->computeElementForce(fileNumber, u->subData(iSub), Findex, locForce);
  subDomain[iSub]->mergeElemStress(locForce, globForce, elemToNode.get());
  delete [] locForce;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementForce(GenDistrVector<Scalar> &u, int fileNumber,
                                           int Findex, double time)
{
  int numElemNodes = elemToNode->numConnect();
  Scalar *globForce = new Scalar[numElemNodes];
  execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemForce,
            globForce, &u, fileNumber, Findex);
  geoSource->outputElemStress(fileNumber, globForce, elemToNode->csize(),
                              elemToNode->ptr(), time);
  delete [] globForce;
}

template<class Scalar>
void
GenDecDomain<Scalar>::computeSubdStress(int iSub, GenDistrVector<Scalar> *globStress, 
                                        GenDistrVector<Scalar> *globWeight, GenDistrVector<Scalar> *u, 
                                        int fileNumber, int Findex) const
{
  GenStackVector<Scalar> gstress(globStress->subData(iSub),globStress->subLen(iSub));
  GenStackVector<Scalar> gweight(globWeight->subData(iSub),globWeight->subLen(iSub));
  subDomain[iSub]->computeStressStrain(fileNumber, u->subData(iSub), 
	                               Findex, gstress.data(), gweight.data());
}

template<class Scalar>
void GenDecDomain<Scalar>::computeSubdElemStress(int iSub, Scalar *glElemStress,
                                                 GenDistrVector<Scalar> *u, int fileNumber, int Findex)  const
{
  Scalar *locStress = new Scalar[subDomain[iSub]->countElemNodes()];
  subDomain[iSub]->computeStressStrain(fileNumber, u->subData(iSub),
                                       Findex, locStress);
  subDomain[iSub]->mergeElemStress(locStress, glElemStress, elemToNode.get());
  delete [] locStress;
}

template<class Scalar>
void
GenDecDomain<Scalar>::computeSubdStress_NL(int iSub, GenDistrVector<Scalar> *globStress,
                                        GenDistrVector<Scalar> *globWeight, DistrGeomState *u, 
                                        Corotator ***allCorot, int *fileNumber, int *Findex,
                                        DistrGeomState *refState) const
{
  // Non-linear version of computeSubdStress
  GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
  subDomain[iSub]->computeStressStrain((*u)[iSub], allCorot[iSub],
                                       *fileNumber, *Findex, globStress->subData(iSub),
                                       globWeight->subData(iSub), subRefState);
}

template<class Scalar>
void GenDecDomain<Scalar>::computeSubdElemStress_NL(int iSub, Scalar *glElemStress,
                                                 DistrGeomState *u, Corotator ***allCorot, 
                                                 int fileNumber, int Findex,
                                                 DistrGeomState *refState) const
{
  // Non-linear version of computeSubdElemStress
  Scalar *locStress = new Scalar[subDomain[iSub]->countElemNodes()];
  GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
  subDomain[iSub]->computeStressStrain((*u)[iSub], allCorot[iSub], fileNumber,
                                       Findex, locStress, (Scalar *) 0, subRefState);
  subDomain[iSub]->mergeElemStress(locStress, glElemStress, elemToNode.get());
  delete [] locStress;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot, 
                                                  int fileNumber, int Findex, double time,
                                                  DistrGeomState *refState) const
{
  // Non-linear version of getElementStressStrain
  int numElemNodes = elemToNode->numConnect();
  Scalar *globStress = new Scalar[numElemNodes];
  execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress_NL,
            globStress, gs, allCorot, fileNumber, Findex, refState);
  geoSource->outputElemStress(fileNumber, globStress, elemToNode->csize(),
                              elemToNode->ptr(), time);
  delete [] globStress;
}

template<class Scalar>
void
GenDecDomain<Scalar>::getStressStrain(DistrGeomState *gs, Corotator ***allCorot,
                                      int fileNumber, int Findex, double time,
                                      DistrGeomState *refState) 
{
 // Non-linear version of getStressStrain
 OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
 if(oinfo.averageFlg == 0) {
   getElementStressStrain(gs, allCorot, fileNumber, Findex, time, refState);
   return;
 }

 // Allocate a distributed vector and initialize it to zero
 // if it hasn't already been allocated.
 if(stress == 0) stress = new GenDistrVector<Scalar>(*nodeInfo);
 if(weight == 0) weight = new GenDistrVector<Scalar>(*nodeInfo);

 stress->zero();
 weight->zero();

 // each subdomain computes its stress vector
 if(Findex != 16) {
   execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress_NL,
             stress, weight, gs, allCorot, &fileNumber, &Findex, refState);
 }

 int numNodes = (domain->outFlag) ? domain->exactNumNodes : geoSource->numNode();

 if(globalStress == 0) globalStress = new Scalar[numNodes]; 
 if(globalWeight == 0) globalWeight = new Scalar[numNodes];

 int i;
 for(i=0; i<numNodes; ++i)
   globalStress[i] = globalWeight[i] = 0.0;

 int iSub;
 for(iSub=0; iSub<numSub; ++iSub) {
   if(Findex != 16) {
     subDomain[iSub]->mergeStress(stress->subData(iSub), weight->subData(iSub),
                                  globalStress, globalWeight, numNodes);
   }
   else {
     subDomain[iSub]->computeContactPressure(globalStress, globalWeight);
   }
 }

 for(i=0; i < numNodes; ++i)  {
   if(globalWeight[i] == 0.0)
     globalStress[i] = 0.0;
   else
     globalStress[i] = globalStress[i]/globalWeight[i];
 }

 if(oinfo.nodeNumber == -1)
   geoSource->outputNodeScalars(fileNumber, globalStress, numNodes, time);
 else
   geoSource->outputNodeScalars(fileNumber, globalStress+oinfo.nodeNumber, 1, time);

 delete [] globalWeight; globalWeight=0;
 delete [] globalStress; globalStress=0;

 delete stress; stress = 0;
 delete weight; weight = 0;
}

template<class Scalar>
void GenDecDomain<Scalar>::setsizeSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1)  sizeSfemStress = geoSource->numNode();  // node-based output
  else {
   std::cerr << "avgnum = " << avgnum << " not implemented in Domain::setsizeSfemStress()" << std::endl;
   sizeSfemStress = 0;
  }
}

template<class Scalar>
void GenDecDomain<Scalar>::updateSfemStress(Scalar* str, int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  int numNodes = geoSource->numNode();
  if(avgnum == 1)  for (int i=0;i<numNodes;++i) globalStress[i] = str[i];
  else std::cerr << "avgnum = " << avgnum << " not implemented in Domain::updateSfemStress()" << std::endl;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementStressStrain(GenDistrVector<Scalar> &u, int fileNumber,
                                                  int Findex, double time, int printFlag)  
{
  // allocate arrays
  int numElemNodes = elemToNode->numConnect();
  Scalar *globStress = new Scalar[numElemNodes];
  // each subdomain computes its stress vector
  execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress,
            globStress, &u, fileNumber, Findex);
  geoSource->outputElemStress(fileNumber, globStress, elemToNode->csize(),
                              elemToNode->ptr(), time);
}

template<class Scalar>
void GenDecDomain<Scalar>::getStressStrain(GenDistrVector<Scalar> &u, int fileNumber, 
           				   int Findex, double time, int printFlag)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementStressStrain(u, fileNumber, Findex, time, printFlag);
    return;
  }

  // Allocate a distributed vector and initialize it to zero
  // if it hasn't already been allocated.
  if(stress == 0) stress = new GenDistrVector<Scalar>(*nodeInfo);
  if(weight == 0) weight = new GenDistrVector<Scalar>(*nodeInfo);

  stress->zero();
  weight->zero();
  

  if(printFlag != 2) {
    // each subdomain computes its stress vector
    if(Findex != 16)
      execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress,
                stress, weight, &u, fileNumber, Findex);
  }


  // allocate global stress and weight arrays 
  int numNodes = (domain->outFlag) ? domain->exactNumNodes : geoSource->numNode();
  if(globalStress == 0) globalStress = new Scalar[numNodes]; 
  if(globalWeight == 0) globalWeight = new Scalar[numNodes];

  if (printFlag != 2) { 
    int i;
    for (i = 0; i < numNodes; ++i)
      globalStress[i] = globalWeight[i] = 0.0;
 

    int iSub;
    for(iSub=0; iSub < numSub; ++iSub) {
      if(Findex != 16) {
        subDomain[iSub]->mergeStress(stress->subData(iSub), weight->subData(iSub),
                                     globalStress, globalWeight, numNodes);
      }
      else {
        subDomain[iSub]->computeContactPressure(globalStress, globalWeight);
      }
    }
    for(i = 0; i < numNodes; ++i)  {
      if(globalWeight[i] == 0.0)
        globalStress[i] = 0.0;
      else
        globalStress[i] = globalStress[i]/globalWeight[i];
    }

  }

  if(printFlag != 1) {
    if(oinfo.nodeNumber == -1)
      geoSource->outputNodeScalars(fileNumber, globalStress, numNodes, time);
    else
      geoSource->outputNodeScalars(fileNumber, globalStress+oinfo.nodeNumber, 1, time);
  }

  delete stress; stress = 0;
  delete weight; weight = 0;
}

template<class Scalar>
void
GenDecDomain<Scalar>::getPrincipalStress(DistrGeomState *gs, Corotator ***allCorot,
                                         int fileNumber, int strIndex, double time,
                                         DistrGeomState *refState)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(gs, allCorot, fileNumber, strIndex, time, refState);
    return;
  }

  // set stress VS. strain for element subroutines (put this in a function!!)
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
  bool direction;

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = false;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = false;
  }
  else if((strIndex==6) || (strIndex==7) || (strIndex==8)) {
    strInd = strIndex-6;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = true;
  }
  else if((strIndex==9) || (strIndex==10) || (strIndex==11)) {
    strInd = strIndex-9;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = true;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // Allocate a distributed vector for stress if it hasn't already been allocated
  if(stress == 0) stress = new GenDistrVector<Scalar>(*nodeInfo);
  if(weight == 0) weight = new GenDistrVector<Scalar>(*nodeInfo);

  // stress storage
  int numNodes = (domain->outFlag) ? domain->exactNumNodes : geoSource->numNode();
  Scalar (*globalAllStress)[6] = new Scalar [numNodes][6];

  // Compute Each Required Stress (all 6) using same routines as for 
  // individual stresses
  int str_loop;
  int Findex;
  for (str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];
    // Initialize distributed vector to zero
    stress->zero();
    weight->zero();
    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress_NL,
             stress, weight, gs, allCorot, &fileNumber, &Findex, refState);
    Scalar *globalStress = new Scalar[numNodes]; 
    Scalar *globalWeight = new Scalar[numNodes];
    for(i = 0; i < numNodes; ++i)
      globalStress[i] = globalWeight[i] = 0.0;
    int iSub;
    for(iSub=0; iSub<numSub; ++iSub)
      subDomain[iSub]->mergeStress(stress->subData(iSub),
                                   weight->subData(iSub),
                                   globalStress,globalWeight,numNodes);
    for(i = 0; i < numNodes; ++i)
      if(globalWeight[i] != 0.0)
        globalAllStress[i][str_loop] = globalStress[i]/globalWeight[i];
      else
        globalAllStress[i][str_loop] = 0.0;
    delete [] globalWeight; globalWeight=0;
    delete [] globalStress; globalStress=0;
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  // PJSA 3-24-05: modified to compute principal direction if required
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = 0; 
  Scalar (*globalPDir)[3] = 0;//DofSet::max_known_nonL_dof 
  Scalar *pdir = 0;
  if(direction) {
    globalPDir = new Scalar[numNodes][3];//DofSet::max_known_nonL_dof 
    pdir = (Scalar *) dbg_alloca(sizeof(Scalar)*3);
  }
  else globalPVec = new Scalar[numNodes];

  for(i = 0; i < numNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globalAllStress[i][j];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec, pdir);
    if(direction) for(j=0;j<3;++j) globalPDir[i][j] = pdir[strInd*3+j];
    else globalPVec[i] = pvec[strInd];
  }
  if(direction) {
    if(oinfo.nodeNumber == -1) 
       geoSource->outputNodeVectors(fileNumber, globalPDir, numNodes, time);
    else 
      geoSource->outputNodeVectors(fileNumber, globalPDir+oinfo.nodeNumber, 1, time);
  }
  else {
    if(oinfo.nodeNumber == -1)
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodes, time);
    else
      geoSource->outputNodeScalars(fileNumber, globalPVec+oinfo.nodeNumber, 1, time);
  }

  delete [] globalAllStress;
  if(globalPVec) delete [] globalPVec;
  if(globalPDir) delete [] globalPDir;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementPrincipalStress(GenDistrVector<Scalar> &u, int fileNumber,
                                                      int strIndex, double time)
{
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
                                                                                                                               
  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  int numElemNodes = elemToNode->numConnect();
  Scalar **globAllStress = new Scalar * [6];
  for(i=0; i<6; ++i) globAllStress[i] = new Scalar[numElemNodes];

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int str_loop;
  int Findex;
  for(str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];

    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress,
              globAllStress[str_loop], &u, fileNumber, Findex);
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = new Scalar[numElemNodes];
  for(i = 0; i < numElemNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globAllStress[j][i];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec);
    globalPVec[i] = pvec[strInd];
  }
  geoSource->outputElemStress(fileNumber, globalPVec, elemToNode->csize(),
                              elemToNode->ptr(), time);

  for(i=0; i<6; ++i) delete [] globAllStress[i];
  delete [] globAllStress;
  delete [] globalPVec;
}

template<class Scalar>
void
GenDecDomain<Scalar>::getElementPrincipalStress(DistrGeomState *gs, Corotator ***allCorot,
                                                int fileNumber, int strIndex, double time,
                                                DistrGeomState *refState)
{
  // Non-linear version of getElementPrincipalStress
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  int numElemNodes = elemToNode->numConnect();
  Scalar **globAllStress = new Scalar * [6];
  for(i=0; i<6; ++i) globAllStress[i] = new Scalar[numElemNodes];

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int str_loop;
  int Findex;
  for(str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];

    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress_NL,
              globAllStress[str_loop], gs, allCorot, fileNumber, Findex, refState);
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = new Scalar[numElemNodes];
  for(i = 0; i < numElemNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globAllStress[j][i];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec);
    globalPVec[i] = pvec[strInd];
  }
  geoSource->outputElemStress(fileNumber, globalPVec, elemToNode->csize(),
                              elemToNode->ptr(), time);

  for(i=0; i<6; ++i) delete [] globAllStress[i];
  delete [] globAllStress;
  delete [] globalPVec;
}


template<class Scalar>
void
GenDecDomain<Scalar>::getPrincipalStress(GenDistrVector<Scalar> &u, int fileNumber, int strIndex, 
			     	         double time)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(u, fileNumber, strIndex, time);
    return;
  }

  // set stress VS. strain for element subroutines (put this in a function!!)
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
  bool direction;

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = false;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = false;
  }
  else if((strIndex==6) || (strIndex==7) || (strIndex==8)) {
    strInd = strIndex-6;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = true;
  }
  else if((strIndex==9) || (strIndex==10) || (strIndex==11)) {
    strInd = strIndex-9;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = true;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // Allocate a distributed vector for stress if it hasn't already been allocated
  if(stress == 0) stress = new GenDistrVector<Scalar>(*nodeInfo);
  if(weight == 0) weight = new GenDistrVector<Scalar>(*nodeInfo);

  // stress storage
  int numNodes = (domain->outFlag) ? domain->exactNumNodes : geoSource->numNode();
  Scalar (*globalAllStress)[6] = new Scalar [numNodes][6];

  // Compute Each Required Stress (all 6) using same routines as for 
  // individual stresses
  int str_loop;
  int Findex;
  for (str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];
    // Initialize distributed vector to zero
    stress->zero();
    weight->zero();
    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress,
              stress, weight, &u, fileNumber, Findex);
    Scalar *globalStress = new Scalar[numNodes]; 
    Scalar *globalWeight = new Scalar[numNodes];
    for(i = 0; i < numNodes; ++i)
      globalStress[i] = globalWeight[i] = 0.0;
    int iSub;
    for(iSub=0; iSub<numSub; ++iSub)
      subDomain[iSub]->mergeStress(stress->subData(iSub),
                                   weight->subData(iSub),
                                   globalStress,globalWeight,numNodes);
    for(i = 0; i < numNodes; ++i)
      if(globalWeight[i] != 0.0)
        globalAllStress[i][str_loop] = globalStress[i]/globalWeight[i];
      else
        globalAllStress[i][str_loop] = 0.0;
    delete [] globalWeight; globalWeight=0;
    delete [] globalStress; globalStress=0;
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  // PJSA 3-24-05: modified to compute principal direction if required
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = 0; 
  Scalar (*globalPDir)[3] = 0;//DofSet::max_known_nonL_dof 
  Scalar *pdir = 0;
  if(direction) {
    globalPDir = new Scalar[numNodes][3];//DofSet::max_known_nonL_dof 
    pdir = (Scalar *) dbg_alloca(sizeof(Scalar)*3);
  }
  else globalPVec = new Scalar[numNodes];

  for(i = 0; i < numNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globalAllStress[i][j];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec, pdir);
    if(direction) for(j=0;j<3;++j) globalPDir[i][j] = pdir[strInd*3+j];
    else globalPVec[i] = pvec[strInd];
  }
  if(direction) {
    if(oinfo.nodeNumber == -1) 
       geoSource->outputNodeVectors(fileNumber, globalPDir, numNodes, time);
    else 
      geoSource->outputNodeVectors(fileNumber, globalPDir+oinfo.nodeNumber, 1, time);
  }
  else {
    if(oinfo.nodeNumber == -1)
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodes, time);
    else
      geoSource->outputNodeScalars(fileNumber, globalPVec+oinfo.nodeNumber, 1, time);
  }

  delete [] globalAllStress;
  if(globalPVec) delete [] globalPVec;
  if(globalPDir) delete [] globalPDir;
}

// -----------------------------
// Nonlinear DecDomain functions
// -----------------------------

template<class Scalar>
void
GenDecDomain<Scalar>::postProcessing(DistrGeomState *geomState, GenDistrVector<Scalar> &extF, Corotator ***allCorot, double x,
                                     SysState<GenDistrVector<Scalar> > *distState, GenDistrVector<Scalar> *aeroF, DistrGeomState *refState,
                                     GenDistrVector<Scalar> *reactions, GenMDDynamMat<Scalar> *dynOps, GenDistrVector<Scalar> *resF)
{
  // NOTE: for dynamic runs, x represents the time
  //       for static runs, x represents the load parameter, lambda
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // check if there are any output files which need to be processed now
  int step = (domain->solInfo().isDynam()) ? int(x/domain->solInfo().getTimeStep()+0.5) : int(x/domain->solInfo().getNLInfo().dlambda+0.5);
  if(geoSource->noOutput(step) && x != domain->solInfo().initialTime) return;

  if(verboseFlag && x == 0)
    filePrint(stderr," ... Postprocessing                 ...\n");

  if(domain->outFlag && domain->nodeTable == 0) domain->makeNodeTable(domain->outFlag);
  int numNodes = (domain->outFlag) ? domain->exactNumNodes : geoSource->numNode();
  Scalar (*glMergedDis)[11] = new Scalar[numNodes][11];
  Scalar (*locMergedDis)[11] = (domain->solInfo().basicDofCoords) ? 0 : new Scalar[numNodes][11];
  Scalar *globVal = 0;  // for output

  int i,j,iSub;
  for(i = 0; i < numNodes; ++i)
    for (j = 0 ; j < 11 ; j++)
    glMergedDis[i][j] = 0.0;

  int isub;
  for(isub = 0; isub < numSub; ++isub)
    subDomain[isub]->mergeDisp(glMergedDis, (*geomState)[isub], locMergedDis);

  // intialize and merge aeroelastic forces from subdomains into global array
  Scalar (*mergedAeroF)[6] = 0;
  if(domain->solInfo().aeroFlag > -1 && aeroF) {
    mergedAeroF = new Scalar[numNodes][6];
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<6; ++j) mergedAeroF[i][j] = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub)
      subDomain[iSub]->mergeForces(mergedAeroF, aeroF->subData(iSub));
  }

  // for nonlinear dynamics: initialize and merge velocities and accelerations from subdomains into global array
  GenDistrVector<Scalar> *v_n = 0, *a_n = 0;
  Scalar (*glMergedVel)[11] = 0, (*glMergedAcc)[11] = 0;
  Scalar (*locMergedVel)[11] = 0, (*locMergedAcc)[11] = 0;
  if(distState) {
    v_n = &distState->getVeloc();
    a_n = &distState->getAccel();
    glMergedVel = new Scalar[numNodes][11];
    glMergedAcc = new Scalar[numNodes][11];
    if(!domain->solInfo().basicDofCoords) {
      locMergedVel = new Scalar[numNodes][11];
      locMergedAcc = new Scalar[numNodes][11];
    }
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<11; ++j) glMergedVel[i][j] = glMergedAcc[i][j] = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->mergeAllVeloc(glMergedVel, v_n->subData(iSub), locMergedVel);
      subDomain[iSub]->mergeAllAccel(glMergedAcc, a_n->subData(iSub), locMergedAcc);
    }
  }

  // merge reaction forces from subdomains into global array
  Scalar (*mergedReactions)[11] = 0;
  if(reactions) {
    mergedReactions = new Scalar[numNodes][11];
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<11; ++j) mergedReactions[i][j] = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub)
      subDomain[iSub]->mergeReactions(mergedReactions, reactions->subData(iSub));
  }

  if(x == domain->solInfo().initialTime) {
    geoSource->openOutputFiles();
  }

  int inode;
  for(i = 0; i < numOutInfo; i++) {
   if(oinfo[i].interval != 0 && step % oinfo[i].interval == 0) {

    Scalar (*xyz)[11] = (oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords) ? glMergedDis : locMergedDis;
    Scalar (*mergedVel)[11] = (oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords) ? glMergedVel : locMergedVel;
    Scalar (*mergedAcc)[11] = (oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords) ? glMergedAcc : locMergedAcc;

    switch(oinfo[i].type) {
     case OutputInfo::FreqRespModes:
     case OutputInfo::Displacement:
       getPrimalVector(i, xyz, numNodes, 3, x);
       break;
     case OutputInfo::Velocity:
       if(distState) getPrimalVector(i, mergedVel, numNodes, 3, x);
       break;
     case OutputInfo::Acceleration:
       if(distState) getPrimalVector(i, mergedAcc, numNodes, 3, x);
       break;
     case OutputInfo::Disp6DOF:
       if(oinfo[i].rotvecouttype != OutputInfo::Euler || !oinfo[i].rescaling) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       getPrimalVector(i, xyz, numNodes, 6, x);
       break;
     case OutputInfo::Velocity6:
       if(oinfo[i].angularouttype != OutputInfo::convected) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       if(distState) getPrimalVector(i, mergedVel, numNodes, 6, x);
       break;
     case OutputInfo::Accel6:
       if(oinfo[i].angularouttype != OutputInfo::convected) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       if(distState) getPrimalVector(i, mergedAcc, numNodes, 6, x);
       break;
     case OutputInfo::Temperature:
       getPrimalScalar(i, xyz, numNodes, 0, x);
       break;
     case OutputInfo::TemperatureFirstTimeDerivative:
       if(distState) getPrimalScalar(i, mergedVel, numNodes, 6, x);
       break;
     case OutputInfo::PressureFirstTimeDerivative:
       if(distState) getPrimalScalar(i, mergedVel, numNodes, 7, x);
       break;
     case OutputInfo::PressureSecondTimeDerivative:
       if(distState) getPrimalScalar(i, mergedAcc, numNodes, 7, x);
       break;
     case OutputInfo::StressXX:
       getStressStrain(geomState, allCorot, i, SXX, x, refState);
       break;
     case OutputInfo::StressYY:
       getStressStrain(geomState, allCorot, i, SYY, x, refState);
       break;
     case OutputInfo::StressZZ:
       getStressStrain(geomState, allCorot, i, SZZ, x, refState);
       break;
     case OutputInfo::StressXY:
       getStressStrain(geomState, allCorot, i, SXY, x, refState);
       break;
     case OutputInfo::StressYZ:
       getStressStrain(geomState, allCorot, i, SYZ, x, refState);
       break;
     case OutputInfo::StressXZ:
       getStressStrain(geomState, allCorot, i, SXZ, x, refState);
       break;
     case OutputInfo::StrainXX:
       getStressStrain(geomState, allCorot, i, EXX, x, refState);
       break;
     case OutputInfo::StrainYY:
       getStressStrain(geomState, allCorot, i, EYY, x, refState);
       break;
     case OutputInfo::StrainZZ:
       getStressStrain(geomState, allCorot, i, EZZ, x, refState);
       break;
     case OutputInfo::StrainXY:
       getStressStrain(geomState, allCorot, i, EXY, x, refState);
       break;
     case OutputInfo::StrainYZ:
       getStressStrain(geomState, allCorot, i, EYZ, x, refState);
       break;
     case OutputInfo::StrainXZ:
       getStressStrain(geomState, allCorot, i, EXZ, x, refState);
       break;
     case OutputInfo::StressVM:
       getStressStrain(geomState, allCorot, i, VON, x, refState);
       break;
     case OutputInfo::StrainVM:
       getStressStrain(geomState, allCorot, i, STRAINVON, x, refState);
       break;
     case OutputInfo::ContactPressure: {
       if(!domain->tdenforceFlag()) 
         getStressStrain(geomState, allCorot, i, CONPRESS, x, refState);
       else
         filePrint(stderr," *** WARNING: Output case %d not supported \n", i);
     } break;
     case OutputInfo::Damage:
       getStressStrain(geomState, allCorot, i, DAMAGE, x, refState);
       break;
     case OutputInfo::EquivalentPlasticStrain:
       getStressStrain(geomState, allCorot, i, EQPLSTRN, x, refState);
       break;
     case OutputInfo::Energies:
       getEnergies_b(geomState, extF, allCorot, i, x, distState, dynOps, aeroF);
       break;
     case OutputInfo::DissipatedEnergy:
       getDissipatedEnergy(geomState, allCorot, i, x);
       break;
     case OutputInfo::StressPR1:
       getPrincipalStress(geomState, allCorot, i, PSTRESS1, x, refState);
       break;
     case OutputInfo::StressPR2:
       getPrincipalStress(geomState, allCorot, i, PSTRESS2, x, refState);
       break;
     case OutputInfo::StressPR3:
       getPrincipalStress(geomState, allCorot, i, PSTRESS3, x, refState);
       break;
     case OutputInfo::StrainPR1:
       getPrincipalStress(geomState, allCorot, i, PSTRAIN1, x, refState);
       break;
     case OutputInfo::StrainPR2:
       getPrincipalStress(geomState, allCorot, i, PSTRAIN2, x, refState);
       break;
     case OutputInfo::StrainPR3:
       getPrincipalStress(geomState, allCorot, i, PSTRAIN3, x, refState);
       break;
     case OutputInfo::DispX:
       getPrimalScalar(i, xyz, numNodes, 0, x);
       break;
     case OutputInfo::DispY:
       getPrimalScalar(i, xyz, numNodes, 1, x);
       break;
     case OutputInfo::DispZ:
       getPrimalScalar(i, xyz, numNodes, 2, x);
       break;
     case OutputInfo::RotX:
       if(oinfo[i].rotvecouttype != OutputInfo::Euler || !oinfo[i].rescaling) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       getPrimalScalar(i, xyz, numNodes, 3, x);
       break;
     case OutputInfo::RotY:
       if(oinfo[i].rotvecouttype != OutputInfo::Euler || !oinfo[i].rescaling) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       getPrimalScalar(i, xyz, numNodes, 4, x);
       break;
     case OutputInfo::RotZ:
       if(oinfo[i].rotvecouttype != OutputInfo::Euler || !oinfo[i].rescaling) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       getPrimalScalar(i, xyz, numNodes, 5, x);
       break;
     case OutputInfo::DispMod:
       if(!globVal) globVal = new Scalar[numNodes];
       for(inode=0; inode<numNodes; ++inode) {
         globVal[inode] = ScalarTypes::sqrt(xyz[inode][0]*xyz[inode][0] +
                                            xyz[inode][1]*xyz[inode][1] +
                                            xyz[inode][2]*xyz[inode][2]);
       }
       geoSource->outputNodeScalars(i, globVal, numNodes, x);
       break;
     case OutputInfo::RotMod:
       if(oinfo[i].rotvecouttype != OutputInfo::Euler || !oinfo[i].rescaling) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       if(!globVal) globVal = new Scalar[numNodes];
       for(inode=0; inode<numNodes; ++inode) {
         globVal[inode] = ScalarTypes::sqrt(xyz[inode][3]*xyz[inode][3] +
                                            xyz[inode][4]*xyz[inode][4] +
                                            xyz[inode][5]*xyz[inode][5]);
       }
       geoSource->outputNodeScalars(i, globVal, numNodes, x);
       break;
     case OutputInfo::TotMod:
       if(oinfo[i].rotvecouttype != OutputInfo::Euler || !oinfo[i].rescaling) {
         filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
         break;
       }
       if(!globVal) globVal = new Scalar[numNodes];
       for(inode=0; inode<numNodes; ++inode) {
         globVal[inode] = ScalarTypes::sqrt(xyz[inode][0]*xyz[inode][0] +
                                            xyz[inode][1]*xyz[inode][1] +
                                            xyz[inode][2]*xyz[inode][2] +
                                            xyz[inode][3]*xyz[inode][3] +
                                            xyz[inode][4]*xyz[inode][4] +
                                            xyz[inode][5]*xyz[inode][5]);
       }
       geoSource->outputNodeScalars(i, globVal, numNodes, x);
       break;
     case OutputInfo::AeroForce: break; // this is done in DistFlExchange.C
     case OutputInfo::AeroXForce:
       if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 0, x);
       break;
     case OutputInfo::AeroYForce:
       if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 1, x);
       break;
     case OutputInfo::AeroZForce:
       if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 2, x);
       break;
     case OutputInfo::AeroXMom:
       if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 3, x);
       break;
     case OutputInfo::AeroYMom:
       if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 4, x);
       break;
     case OutputInfo::AeroZMom:
       if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 5, x);
       break;
     case OutputInfo::Reactions:
       if(reactions) getPrimalVector(i, mergedReactions, numNodes, 3, x);
       break;
     case OutputInfo::Reactions6:
       if(reactions) getPrimalVector(i, mergedReactions, numNodes, 6, x);
       break;
     case OutputInfo::TDEnforcement: {
       if(domain->tdenforceFlag()) {
         double *plot_data = new double[numNodes];
         if(oinfo[i].tdenforc_var == 1) // CONFACE
           for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.5;
         else
           for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.0;
         for(int iMortar=0; iMortar<domain->GetnMortarConds(); iMortar++) {
           domain->GetMortarCond(iMortar)->get_plot_variable(oinfo[i].tdenforc_var,plot_data);
         }
         if(oinfo[i].nodeNumber == -1)
           geoSource->outputNodeScalars(i, plot_data, numNodes, x);
         else
           geoSource->outputNodeScalars(i, &plot_data[oinfo[i].nodeNumber], 1, x);
         delete [] plot_data;
       }
       else filePrint(stderr," *** WARNING: Output case %d not supported \n", i);
     } break;
     case OutputInfo::DeletedElements: {
       for(int iSub = 0; iSub < numSub; ++iSub) {
         std::vector<std::pair<double,int> > &deletedElements = subDomain[iSub]->getDeletedElements();
         for(std::vector<std::pair<double,int> >::iterator it = deletedElements.begin(); it != deletedElements.end(); ++it) {
           filePrint(oinfo[i].filptr, " %12.6e  %9d          Undetermined\n", it->first, it->second+1);
           fflush(oinfo[i].filptr);
         }
         deletedElements.clear();
       }
     } break;
     case OutputInfo::Statevector:
     case OutputInfo::Velocvector:
     case OutputInfo::Accelvector:
     case OutputInfo::InternalStateVar:
     case OutputInfo::DualStateVar:
     case OutputInfo::Forcevector:
     case OutputInfo::Constraintvector:
     case OutputInfo::Constraintviolation:
     case OutputInfo::Residual:
     case OutputInfo::Jacobian:
     case OutputInfo::RobData:
     case OutputInfo::SampleMesh:
        break;
     default:
       filePrint(stderr," *** WARNING: Output case %d not implemented\n", i);
       break;
   }
  }
 }
 if(globVal) delete [] globVal;
 if(aeroF) delete [] mergedAeroF;
 if(glMergedDis) delete [] glMergedDis;
 if(locMergedDis) delete [] locMergedDis;
 if(glMergedVel) delete [] glMergedVel;
 if(glMergedAcc) delete [] glMergedAcc;
 if(locMergedVel) delete [] locMergedVel;
 if(locMergedAcc) delete [] locMergedAcc;
 if(mergedReactions) delete [] mergedReactions;
}

// element vector distributed vector info
// each element has a certain length based on its element stiffness matrix
template<class Scalar>
DistrInfo*
GenDecDomain<Scalar>::elementVectorInfo()
{
  if(!eleVecInfo) {
    eleVecInfo = new DistrInfo;
    makeBasicDistrInfo(*eleVecInfo, &Domain::maxNumDOF);
  }

  return eleVecInfo;
}

// prescribed boundary condition distributed vector info
template<class Scalar>
DistrInfo*
GenDecDomain<Scalar>::pbcVectorInfo()
{
 if(!bcVecInfo) {
   bcVecInfo = new DistrInfo;
   makeBasicDistrInfo(*bcVecInfo, &Domain::nDirichlet);
 }
 return bcVecInfo;
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeInternalInfo()
{
 startTimerMemory(mt.makeInternalInfo, mt.memoryInternal);

 makeSolVecInfo();
 makeSysVecInfo();
 
 stopTimerMemory(mt.makeInternalInfo, mt.memoryInternal);
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSolVecInfo()
{
 // Create internal Distributed information, only for unconstrained dofs
 if(!internalInfo) {
   internalInfo = new DistrInfo();
   makeBasicDistrInfo(*internalInfo, &Domain::numUncon);

   if(domain->solInfo().inpc || domain->solInfo().timeIntegration == SolverInfo::Qstatic) {
     setNonTrivialMasterFlag(*internalInfo);
   } else {
     internalInfo->setMasterFlag();
   }
 }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSysVecInfo()
{
 // Create internal Distributed information for all dofs, both constrained and unconstrained
 if(!internalInfo2) {
   internalInfo2 = new DistrInfo();
   makeBasicDistrInfo(*internalInfo2, &Domain::numdof);
   internalInfo2->setMasterFlag();
 }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeNodeInfo()
{
 startTimerMemory(mt.makeInternalInfo, mt.memoryInternal);

 // Create nodal Distributed information (used for nodal stress output)
 nodeInfo = new DistrInfo;
 makeBasicDistrInfo(*nodeInfo, &Domain::numNodes);
#ifdef DISTRIBUTED
 nodeInfo->computeOffsets();
#else
 nodeInfo->setMasterFlag();
#endif
 stopTimerMemory(mt.makeInternalInfo, mt.memoryInternal);
}

template<class Scalar>
const DistrInfo&
GenDecDomain<Scalar>::masterSolVecInfo() const
{
 if(!masterSolVecInfo_) {
   GenDecDomain<Scalar> *self = const_cast<GenDecDomain<Scalar> *>(this);
   self->masterSolVecInfo_ = new DistrInfo;
   self->makeBasicDistrInfo(*masterSolVecInfo_, &Domain::numUncon);
   self->setNonTrivialMasterFlag(*masterSolVecInfo_);
 }
 return *masterSolVecInfo_;
}

template<class Scalar>
DistrInfo&
GenDecDomain<Scalar>::ndVecInfo()
{
 // Create nodal Distributed information (used for nodal stress output)
 if(!nodeVecInfo) {
   nodeVecInfo = new DistrInfo;
   makeBasicDistrInfo(*nodeVecInfo, &Domain::numNodes);
   nodeVecInfo->setMasterFlag();
 }
 return *nodeVecInfo;
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeBasicDistrInfo(DistrInfo &info, int(Domain::*countFunc)() const) {
 info.domLen = new int[numSub];
 info.numDom = numSub;
 int totLen = 0;
 for(int iSub = 0; iSub < numSub; ++iSub) {
   info.domLen[iSub] = (subDomain[iSub]->*countFunc)();
   totLen += info.domLen[iSub];
 }
 info.len = totLen;
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeBlockCyclicDistrInfo(DistrInfo &info, int globalLen, int blockSize)
{
  info.domLen = new int[numSub];
  info.numDom = numSub;
  int totLen = 0;
  Rom::BlockCyclicMap bcMap(globalLen, blockSize, numCPU, numSub);
  for(int iSub = 0; iSub < numSub; ++iSub) {
    info.domLen[iSub] = bcMap.subLen(myCPU, iSub);
    totLen += info.domLen[iSub];
  }
  info.len = totLen;
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeNonOverlappingDistrInfo(DistrInfo &info)
{
  info.domLen = new int[numSub];
  info.numDom = numSub;
  int totLen = 0;
  for(int iSub = 0; iSub < numSub; ++iSub) {
    const bool *subMasterFlag = subDomain[iSub]->getInternalMasterFlag();
    info.domLen[iSub] = 0;
    for(int i=0; i<subDomain[iSub]->numUncon(); ++i) if(subMasterFlag[i]) info.domLen[iSub]++;
    totLen += info.domLen[iSub];
  }
  info.len = totLen;
}

template<class Scalar>
void
GenDecDomain<Scalar>::setNonTrivialMasterFlag(DistrInfo &info)
{
  bool *internalMasterFlag = new bool[info.len];
  info.computeOffsets();
  for(int iSub = 0; iSub < numSub; ++iSub) {
    const bool *subMasterFlag = subDomain[iSub]->getInternalMasterFlag();
    const int subDofCount = info.domLen[iSub];
    const int subOffset = info.subOffset[iSub];
    std::copy(subMasterFlag, subMasterFlag + subDofCount, internalMasterFlag + subOffset);
  }
  info.setMasterFlag(internalMasterFlag);
}

template<class Scalar>
void
GenDecDomain<Scalar>::constructSubDomains(int iSub)
{
  subDomain[iSub] = GenSubDomainFactory<Scalar>::getFactory()->
    createSubDomain(*domain, iSub, *subToElem, *subToNode, localSubToGl[iSub]);
}

template<class Scalar>
void
GenDecDomain<Scalar>::getSharedDOFs()
{
  startTimerMemory(mt.makeInterface, mt.memoryInterface);

  FSCommPattern<int> nodeIntPat(communicator, cpuToSub.get(), myCPU, FSCommPattern<int>::CopyOnSend);
  for(int i=0; i<numSub; ++i) subDomain[i]->setNodeCommSize(&nodeIntPat);
  nodeIntPat.finalize();

  paralApplyToAll(numSub, subDomain, &FetiBaseSub::sendDOFList, &nodeIntPat);
  nodeIntPat.exchange();
  paralApply(subDomain, &FetiSub<Scalar>::gatherDOFList, &nodeIntPat);
  paralApply(subDomain, &FetiBaseSub::gatherDOFListPlus, &nodeIntPat);


  stopTimerMemory(mt.makeInterface, mt.memoryInterface);
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeCorners()
{
  if(!( isFeti(domain->solInfo().solvercntl->type)
       && domain->solInfo().solvercntl->fetiInfo.version == FetiInfo::fetidp)) return;
  if(verboseFlag) filePrint(stderr, " ... Selecting the Corners          ...\n");

  FSCommPattern<int> cpat(communicator, cpuToSub.get(), myCPU, FSCommPattern<int>::CopyOnSend);
  for(int i=0; i<numSub; ++i) subDomain[i]->setNodeCommSize(&cpat);
  cpat.finalize();
 
  std::vector<FetiSubCornerHandler *>cornerHandler(numSub); // deleted by cornerSelector
  execParal(numSub, this, &GenDecDomain<Scalar>::makeCornerHandler, cornerHandler.data());
  CornerSelector cornerSelector(globalNumSub, numSub, std::move(cornerHandler), &cpat, communicator);
  cornerSelector.makeCorners();
  grToSub = cornerSelector.yieldGrToSub();
  execParal(numSub, this, &GenDecDomain<Scalar>::setLocalCorners, cornerSelector.handlers().data());

  paralApply(subDomain, &BaseSub::makeCCDSA);
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeCornerHandler(int iSub, FetiSubCornerHandler **cornerHandler)
{
  cornerHandler[iSub] = subDomain[iSub]->getCornerHandler();
}

template<class Scalar>
void
GenDecDomain<Scalar>::setLocalCorners(int iSub, FetiSubCornerHandler **cornerHandler)
{
  subDomain[iSub]->setCorners(cornerHandler[iSub]->getCorners());
}

template<class Scalar>
void GenDecDomain<Scalar>::distributeBCs() 
{
  startTimerMemory(mt.distributeBCs,  mt.memoryDistBC);

  int i, iSub, subI;

  std::vector<int> nDirichletPerSub (numSub, 0);
  std::vector<int> nNeumannPerSub   (numSub, 0);
  std::vector<int> nIDisPerSub      (numSub, 0);
  std::vector<int> nIDis6PerSub     (numSub, 0);
  std::vector<int> nIVelPerSub      (numSub, 0);
 
  // get bc's from domain
  BCond* dbc = 0;
  BCond* nbc = 0;
  BCond* cvbc = 0;
  BCond* iDis = 0;
  BCond* iDis6 = 0;
  BCond* iVel = 0;

  int numDirichlet = domain->nDirichlet();   dbc = domain->getDBC();
  int numNeuman    = domain->nNeumann();     nbc = domain->getNBC();
  int numIDis      = domain->numInitDisp();  iDis = domain->getInitDisp();
  int numIDis6     = domain->numInitDisp6(); iDis6 = domain->getInitDisp6();
  int numIVel   = domain->numInitVelocity(); iVel = domain->getInitVelocity();

  // Count the number of boundary conditions per subdomain
  int numDispDirichlet = 0; // number of displacement dirichlet BCs
  for(i = 0; i < numDirichlet; ++i) {
    int node = dbc[i].nnum;
    if(dbc[i].dofnum < 6) numDispDirichlet++; 
    for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nDirichletPerSub[ subI ]++;
    }
  }
  domain->setNumDispDirichlet(numDispDirichlet);

  for(i = 0; i < numNeuman; ++i) {
    int node = nbc[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nNeumannPerSub[subI]++;
  }

  for(i = 0; i < numIDis; ++i) {
    int node = iDis[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nIDisPerSub[subI]++;
  }

  for(i = 0; i < numIDis6; ++i) {
    int node = iDis6[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nIDis6PerSub[subI]++;
  }

  for(i = 0; i < numIVel; ++i) {
    int node = iVel[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nIVelPerSub[subI]++;
  }

  // Set BC's for all subdomains
  BCond **subBC = new BCond *[numSub];
  for(iSub = 0; iSub < numSub; ++iSub) {
    subBC[iSub] = new BCond[nDirichletPerSub[iSub]];
    nDirichletPerSub[iSub] = 0;
  }
  for(i = 0; i < numDirichlet; ++i) {
    int node = dbc[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nDirichletPerSub[subI]] = dbc[i];
        nDirichletPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setDirichlet(nDirichletPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nNeumannPerSub[iSub]];
    nNeumannPerSub[iSub] = 0;
  }
  for(i = 0; i < numNeuman; ++i) {
    int node = nbc[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nNeumannPerSub[subI]] = nbc[i];
        nNeumannPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setNeuman(nNeumannPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nIDisPerSub[iSub]];
    nIDisPerSub[iSub] = 0;
  }
  for(i = 0; i < numIDis; ++i) {
    int node = iDis[i].nnum;
    for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nIDisPerSub[subI]] = iDis[i];
        nIDisPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setIDis(nIDisPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nIVelPerSub[iSub]];
    nIVelPerSub[iSub] = 0;
  }
  for(i = 0; i < numIVel; ++i) {
    int node = iVel[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nIVelPerSub[subI]] = iVel[i];
        nIVelPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setIVel(nIVelPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nIDis6PerSub[iSub]];
    nIDis6PerSub[iSub] = 0;
  }
  for(i = 0; i < numIDis6; ++i) {
    int node = iDis6[i].nnum;
    for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nIDis6PerSub[subI]] = iDis6[i];
        nIDis6PerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setIDis6(nIDis6PerSub[iSub], subBC[iSub]);
  }
  delete [] subBC;

  if(domain->numSSN() > 0) {

    int totEle = (domain->getEset())->last();
    std::vector<int> eleTouch(totEle);
    std::vector<int> eleCount(totEle);

    // Sommerfeld
    for (i=0;i<totEle;i++) eleTouch[i] = -1;
    int *somToSub = new int[domain->numSommer];
    for (i=0;i<domain->numSommer;i++) somToSub[i] = -1;
    for(int iSommEle =0; iSommEle < domain->numSommer; ++iSommEle) {
      int iele = domain->sommer[iSommEle]->findEle(domain->nodeToElem,
                                                   eleTouch.data(), eleCount.data(), iSommEle,
                                                   &domain->getElementSet());
      if(iele < 0) {
        fprintf (stderr, "Error in the Sommerfeld b.c.'s - aborting\n");
        exit(0);
      }
      else {
        subI = glSubToLocal[(*elemToSub)[iele][0]];
        if(subI >= 0) somToSub[iSommEle] = subI;
      }
    }

    // Scatterer
    for (i=0;i<totEle;i++) eleTouch[i] = -1;
    int *scaToSub = new int[domain->numScatter];
    for (i=0;i<domain->numScatter;i++) scaToSub[i] = -1;
    int *sBoundFlag = new int[domain->numNode()];
    for(i=0;i<domain->numNode();i++) sBoundFlag[i] = 0;
    for(int iScatter=0; iScatter < domain->numScatter; ++iScatter) {
      for(int iNode = 0;iNode<domain->scatter[iScatter]->numNodes();iNode++) {
        int ndNum = (domain->scatter[iScatter]->getNodes())[iNode];
        sBoundFlag[ndNum] = 1;
      }
      int iele = domain->scatter[iScatter]->findEle(domain->nodeToElem,
                                                    eleTouch.data(), eleCount.data(), iScatter,
                                                    &domain->getElementSet(),1);
      if(iele < 0) {
        fprintf (stderr, "Error in the scatterer b.c.'s - aborting %d \n",iScatter+1);
        exit(0);
      }
      else {
        subI = glSubToLocal[(*elemToSub)[iele][0]];
        if(subI >= 0) scaToSub[iScatter] = subI;
      }
    }

   // Wet
   for (i=0;i<totEle;i++) eleTouch[i] = -1;
   int (*wetToSub)[2] = new int[domain->numWet][2];
   for (i=0;i<domain->numWet;i++) wetToSub[i][0] = wetToSub[i][1] = -1;
   int iWetEle;
   for(iWetEle =0; iWetEle < domain->numWet; ++iWetEle) {
     int iele[2];
     domain->wet[iWetEle]->findBothEle(domain->nodeToElem, eleTouch.data(),
        eleCount.data(), iWetEle, &domain->getElementSet(),iele);
     if(iele[0]<0 || iele[1]<0) {
       fprintf (stderr, "Error in the wet b.c.'s - aborting\n");
       exit(0);
     }
     else {
       int subI = glSubToLocal[(*elemToSub)[iele[0]][0]];
       if(subI >= 0) wetToSub[iWetEle][0] = subI;
       int subI2 = glSubToLocal[(*elemToSub)[iele[1]][0]];
       if(subI2 >= 0) wetToSub[iWetEle][1] = subI2;
     }
   }

    // Implicit Neumann
    for (i=0;i<totEle;i++) eleTouch[i] = -1;
    int *neumToSub = new int[domain->numNeum];
    for (i=0;i<domain->numNeum;i++) neumToSub[i] = -1;
    for(int iNeum=0; iNeum < domain->numNeum; ++iNeum) {
      int iele = domain->neum[iNeum]->findEle(domain->nodeToElem,
                                              eleTouch.data(), eleCount.data(), iNeum);
      if(iele < 0) {
        fprintf (stderr, "Error in the Neumann b.c.'s - aborting\n");
        exit(0);
      }
      else {
        subI = glSubToLocal[(*elemToSub)[iele][0]];
        if(subI >= 0) neumToSub[iNeum] = subI;
      }
    }
 
    execParal(numSub, this, &GenDecDomain<Scalar>::distribBC, subDomain, domain,
              somToSub, scaToSub, neumToSub, wetToSub, sBoundFlag);

    delete[] sBoundFlag;
    delete[] somToSub;
    delete[] wetToSub;
    delete[] scaToSub;
    delete[] neumToSub;
  }

  // complex nodal boundary conditions 
  if((domain->numComplexDirichlet > 0) || (domain->numComplexNeuman > 0)) {
    int *nComplexDirichletPerSub = new int[numSub]; 
    int *nComplexNeumannPerSub = new int[numSub];
    ComplexBCond **subCBC = new ComplexBCond *[numSub];
    for(iSub = 0; iSub < numSub; ++iSub)
      nComplexDirichletPerSub[iSub] = nComplexNeumannPerSub[iSub] = 0;

    // count Complex Dirichlet
    for(i =0; i < domain->numComplexDirichlet; ++i) {
      int node = domain->cdbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
          nComplexDirichletPerSub[ subI ] += 1;
      }
    }

    // count Complex Neuman
    for(i = 0; i < domain->numComplexNeuman; ++i) {
      int node = domain->cnbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
          nComplexNeumannPerSub[ subI ] += 1;
      }
    }

    // distribute Complex Dirichlet
    for(iSub = 0; iSub < numSub; ++iSub) {
      subCBC[iSub] = new ComplexBCond[nComplexDirichletPerSub[iSub]];
      nComplexDirichletPerSub[iSub] = 0;
    }
    for(i =0; i < domain->numComplexDirichlet; ++i) {
      int node = domain->cdbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
          subCBC[subI][ nComplexDirichletPerSub[subI] ] = domain->cdbc[i];
          nComplexDirichletPerSub[subI] += 1;
        }
      }
    }
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->setComplexDirichlet(nComplexDirichletPerSub[iSub], subCBC[iSub]);
    }

    // Distribute Complex Neuman
    for(iSub = 0; iSub < numSub; ++iSub) {
      subCBC[iSub] = new ComplexBCond[nComplexNeumannPerSub[iSub]];
      nComplexNeumannPerSub[iSub] = 0;
    }
    for(i =0; i < domain->numComplexNeuman; ++i) {
      int node = domain->cnbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
          subCBC[subI][ nComplexNeumannPerSub[subI] ] = domain->cnbc[i];
          nComplexNeumannPerSub[ subI ] += 1;
        }
      }
    }
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->setComplexNeuman(nComplexNeumannPerSub[iSub], subCBC[iSub]);
    }
  
    delete [] nComplexDirichletPerSub;
    delete [] nComplexNeumannPerSub; 
    delete [] subCBC;
  }

  stopTimerMemory(mt.distributeBCs, mt.memoryDistBC);
}

template<class Scalar>
void
GenDecDomain<Scalar>::distribBC(int iSub, gsl::span<GenSubDomain<Scalar> *> sd, Domain *domain,
                                int *somToSub, int *scaToSub, int *neumToSub, int (*wetToSub)[2], int *sBoundFlag)
{
 int iS;
 for(iS=0;iS<domain->numSommer;iS++) {
   if (somToSub[iS] == iSub) sd[iSub]->addSommer(domain->sommer[iS]);
 }
 for(iS=0;iS<domain->numScatter;iS++)
   if (scaToSub[iS] == iSub) sd[iSub]->addScatter(domain->scatter[iS]);
 for(iS=0;iS<domain->numNeum;iS++)
   if (neumToSub[iS] == iSub) sd[iSub]->addNeum(domain->neum[iS]);
 for(iS=0;iS<domain->numWet;iS++) {
   if (wetToSub[iS][0] == iSub && wetToSub[iS][1] == iSub) {
      sd[iSub]->addWet(domain->wet[iS]);
   }
   else if (wetToSub[iS][0] == iSub && (wetToSub[iS][0] != wetToSub[iS][1])) {
      SommerElement *se = domain->wet[iS]->clone();
      se->el2 = 0;
      sd[iSub]->addWet(se);
   }
   else if (wetToSub[iS][1] == iSub && (wetToSub[iS][0] != wetToSub[iS][1])) {
      SommerElement *se = domain->wet[iS]->clone();
      se->el = se->el2;
      se->el2 = 0;
      sd[iSub]->addWet(se);
   }
 }

 int i;
 for(i=0;i<domain->numNode();i++) {
   if(sBoundFlag[i]) 
     for(int jSub =0; jSub < nodeToSub->num(i); ++jSub) {
       int subJ = glSubToLocal[(*nodeToSub)[i][jSub]];
       if(subJ < 0) continue;
       if(subJ==iSub) sd[iSub]->addSBoundNode(i);
     }
 }

 // RT: Added for SO3
 int c=0;
 sd[iSub]->subScaToSca = new int[sd[iSub]->numScatter];
 for(iS=0;iS<domain->numScatter;iS++)
   if (scaToSub[iS] == iSub) sd[iSub]->subScaToSca[c++] = iS;
}

template<class Scalar>
void GenDecDomain<Scalar>::addUserForce(GenDistrVector<Scalar> &f, Scalar *userDefineForce) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->addUserForce(f.subData(iSub), userDefineForce);
}

template<class Scalar>
void GenDecDomain<Scalar>::addCtrl(GenDistrVector<Scalar> &force, Scalar *ctrfrc) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->addCtrl(force.subData(iSub), ctrfrc);
}

template<class Scalar>
void GenDecDomain<Scalar>::extractControlData(GenDistrVector<Scalar> &disp, GenDistrVector<Scalar> &vel,
                                              GenDistrVector<Scalar> &acc, Scalar *ctrdisp,
                                              Scalar *ctrvel, Scalar *ctracc) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->extractControlData(disp.subData(iSub),
        vel.subData(iSub), acc.subData(iSub), ctrdisp, ctrvel, ctracc);
}

template<class Scalar>
void GenDecDomain<Scalar>::distributeControlLawData() 
{
  // get Global Control Law
  ControlLawInfo *claw = geoSource->getControlLaw();
  if(!claw) return;

  int i, iSub, subI;

  // allocate Control Law Data
  int *nSensorsPerSub = new int[numSub];
  int *nActuatorsPerSub = new int[numSub];
  int *nUserDispPerSub = new int[numSub];
  int *nUserForcePerSub = new int[numSub];

  // zero all control law counters
  for (iSub = 0; iSub < numSub; ++iSub)  {
    nSensorsPerSub[iSub] = 0;
    nActuatorsPerSub[iSub] = 0;
    nUserDispPerSub[iSub] = 0;
    nUserForcePerSub[iSub] = 0;
  }

  // Count the number of control law data per subdomain
  for(i = 0; i < claw->numSensor; ++i) {
    int node = claw->sensor[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nSensorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numActuator; ++i) {
    int node = claw->actuator[i].nnum;
    for (iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nActuatorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserDisp; ++i) {
    int node = claw->userDisp[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nUserDispPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserForce; ++i) {
    int node = claw->userForce[i].nnum;
    for (iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nUserForcePerSub[subI]++;
    }
  }

  // Set Control Law for all subdomains
  BCond **subSensor = new BCond *[numSub];
  BCond **subActuator = new BCond *[numSub];
  BCond **subUserDisp = new BCond *[numSub];
  BCond **subUserForce = new BCond *[numSub];

  int **locToGlSensorData = new int*[numSub];
  int **locToGlActuatorData = new int*[numSub];
  int **locToGlUserDispData = new int*[numSub];
  int **locToGlUserForceData = new int*[numSub];

  for(iSub = 0; iSub < numSub; ++iSub) {
    subSensor[iSub] = new BCond[nSensorsPerSub[iSub]];
    subActuator[iSub] = new BCond[nActuatorsPerSub[iSub]];
    subUserDisp[iSub] = new BCond[nUserDispPerSub[iSub]];
    subUserForce[iSub] = new BCond[nUserForcePerSub[iSub]];
    locToGlSensorData[iSub] = new int[ nSensorsPerSub[iSub] ];
    locToGlActuatorData[iSub] = new int[ nActuatorsPerSub[iSub] ];
    locToGlUserDispData[iSub] = new int[ nUserDispPerSub[iSub] ];
    locToGlUserForceData[iSub] = new int[ nUserForcePerSub[iSub] ];

    // reinitialize the counters
    nSensorsPerSub[iSub] = 0;
    nActuatorsPerSub[iSub] = 0;
    nUserDispPerSub[iSub] = 0;
    nUserForcePerSub[iSub] = 0;
  }

  for(i = 0; i < claw->numSensor; ++i) {
    int node = claw->sensor[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subSensor[subI][ nSensorsPerSub[subI] ] = claw->sensor[i];
      locToGlSensorData[subI][ nSensorsPerSub[subI] ] = i;
      nSensorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numActuator; ++i) {
    int node = claw->actuator[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subActuator[subI][ nActuatorsPerSub[subI] ] = claw->actuator[i];
      locToGlActuatorData[subI][ nActuatorsPerSub[subI] ] = i;
      nActuatorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserDisp; ++i) {
    int node = claw->userDisp[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subUserDisp[subI][ nUserDispPerSub[subI] ] = claw->userDisp[i];
      locToGlUserDispData[subI][ nUserDispPerSub[subI] ] = i;
      nUserDispPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserForce; ++i) {
    int node = claw->userForce[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      subI = (*nodeToSub)[node][iSub];
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subUserForce[subI][ nUserForcePerSub[subI] ] = claw->userForce[i];
      locToGlUserForceData[subI][ nUserForcePerSub[subI] ] = i;
      nUserForcePerSub[subI]++;
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub)  {
    ControlLawInfo *locCLaw = new ControlLawInfo();
    locCLaw->fileName = claw->fileName;
    locCLaw->routineName = claw->routineName;
    locCLaw->numSensor = nSensorsPerSub[iSub];
    locCLaw->sensor = subSensor[iSub];
    locCLaw->numActuator = nActuatorsPerSub[iSub];
    locCLaw->actuator = subActuator[iSub];
    locCLaw->numUserDisp = nUserDispPerSub[iSub];
    locCLaw->userDisp = subUserDisp[iSub];
    locCLaw->numUserForce = nUserForcePerSub[iSub];
    locCLaw->userForce = subUserForce[iSub];

    subDomain[iSub]->setControlData(locCLaw, locToGlSensorData[iSub],
                locToGlActuatorData[iSub], locToGlUserDispData[iSub],
                locToGlUserForceData[iSub]);
  }
  delete [] locToGlSensorData;
  delete [] locToGlActuatorData;
  delete [] locToGlUserDispData;
  delete [] locToGlUserForceData;
  delete [] subUserForce;
  delete [] subUserDisp;
  delete [] subActuator;
  delete [] subSensor;
  delete [] nUserForcePerSub;
  delete [] nUserDispPerSub;
  delete [] nActuatorsPerSub;
  delete [] nSensorsPerSub;
}

// NOTE: this creates global element connectivity in shared
//       and cluster element connectivity in distributed
template<class Scalar>
void GenDecDomain<Scalar>::createElemToNode() 
{
  mt.memoryElemToNode -= memoryUsed();
  int iSub;

  // get total number of elements
  int size = 0;
  for(iSub = 0; iSub < numSub; iSub++)
    size += subDomain[iSub]->numElements();

  // allocate connectivity pointers
  int *ptr = new int[size+1];

  // count number of targets
  int numTargetNodes = 0;
  for(iSub = 0; iSub < numSub; iSub++) {
    Elemset &elems = subDomain[iSub]->getElementSet();
    for(int iElem = 0; iElem < subDomain[iSub]->numElements(); iElem++) {
      numTargetNodes += elems[iElem]->numNodes();
      ptr[ elems[iElem]->getGlNum() ] = elems[iElem]->numNodes();
    }
  }

  // allocate targets
  int *nodeTargets = new int[numTargetNodes];

  // populate targets
  int lastPtr = 0;
  int count = 0;
  for(iSub = 0; iSub < numSub; iSub++) {
    Elemset &elems = subDomain[iSub]->getElementSet();
    for(int iElem = 0; iElem < subDomain[iSub]->numElements(); iElem++) {
      ptr[count++] = lastPtr;
      lastPtr += elems[iElem]->numNodes();  
      elems[iElem]->nodes(nodeTargets + elems[iElem]->getGlNum());
    }
  }
  ptr[size] = lastPtr;
  elemToNode = std::make_unique<Connectivity>(size, ptr, nodeTargets);
      
  mt.memoryElemToNode += memoryUsed();
}

template<class Scalar>
void
GenDecDomain<Scalar>::getSharedMPCs()
{
	if(numDualMpc) {
		Connectivity *subToMpc = mpcToSub_dual->altReverse(); // reverse without reordering
		Connectivity *subToSub_mpc = subToMpc->altTranscon(*mpcToSub_dual); // remove connection with self unless pure internal
		paralApply(subDomain, &BaseSub::makeMpcInterface, subToMpc, *mpcToSub_dual, subToSub_mpc);
		delete subToMpc;
		delete subToSub_mpc;
	}
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeMpcToMpc()
{
// This function constructs the mpcToMpc connectivity which will be used 
// to assemble the matrix CC^T for the generalized preconditioner (rixen/dual method)
 if(numDualMpc) {

  // build local mpcToMpc connectivities
  paralApplyToAll(numSub, subDomain, &GenSubDomain<Scalar>::makeLocalMpcToMpc);

  int i, j, iSub;
  int size = mpcToSub_dual->csize();
  int *flags = new int[size];
  for(i=0; i<size; ++i) flags[i] = -1;

  // identify the number of connections from MPC i
  int *np = new int[size+1];
  int cp = 0;
  for(i=0; i<size; ++i) {
    np[i] = cp;
    for(iSub=0; iSub < mpcToSub_dual->num(i); ++iSub) {
      int subi = glSubToLocal[(*mpcToSub_dual)[i][iSub]];
      if(subi == -1) continue;
      int ln = subDomain[subi]->globalToLocalMPC[i];
      Connectivity *locMpcToMpc = subDomain[subi]->localMpcToMpc;
      for(j=0; j < locMpcToMpc->num(ln); ++j) {
        int lj = (*locMpcToMpc)[ln][j];
        int gj = subDomain[subi]->localToGlobalMPC[lj];
        if(flags[gj] != i) {
          cp++;
          flags[gj] = i;
        }
      }
    }
  }
  np[size] = cp;

  // Now allocate and fill the new target
  for(i=0; i<size; ++i) flags[i] = -1;
  int *ntg = new int[cp];
  cp = 0;
  for(i=0; i<size; ++i) {
    for(iSub=0; iSub < mpcToSub_dual->num(i); ++iSub) {
      int subi = glSubToLocal[(*mpcToSub_dual)[i][iSub]];
      if(subi == -1) continue;
      int ln = subDomain[subi]->globalToLocalMPC[i];
      Connectivity *locMpcToMpc = subDomain[subi]->localMpcToMpc;
      for(j=0; j < locMpcToMpc->num(ln); ++j) {
        int lj = (*locMpcToMpc)[ln][j];
        int gj = subDomain[subi]->localToGlobalMPC[lj];
        if(flags[gj] != i) {
          ntg[cp] = gj;
          cp++;
          flags[gj] = i;
        }
      }
    }
  }
 
  delete [] flags;
  mpcToMpc = std::make_unique<Connectivity>(size, np, ntg); // for all mpcs on this cpu
  if(domain->solInfo().solvercntl->fetiInfo.bmpc) {
    mpcToMpc = std::make_unique<Connectivity>( mpcToMpc->modify() );
  }

#ifdef DISTRIBUTED
  makeGlobalMpcToMpc(mpcToMpc.get()); // merge all cpus
#endif

  FetiInfo *finfo = &domain->solInfo().getFetiInfo();
  if(finfo->mpc_precno == FetiInfo::autoSelectCCt) {
    if(mpcToMpc->isDiagonal()) finfo->mpc_precno = FetiInfo::diagCCt;
    else if((numCPU > 1) || (numDualMpc < 100)) finfo->mpc_precno = FetiInfo::globalCCt;
    else finfo->mpc_precno = FetiInfo::superBlockDiagCCt;
  }
 }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeGlobalMpcToMpc(Connectivity *_procMpcToMpc)
{
// This function constructs the global mpcToMpc connectivity which will be used
// to assemble the matrix CC^T for the generalized preconditioner (rixen/dual method)
	int i, j, k;
	int pid = myCPU;
	std::vector<int> size(numCPU);
	std::vector<int> numtarget(numCPU);
	for(i=0; i<numCPU; ++i) {
		if(i == pid) {
			size[i] = _procMpcToMpc->csize();
			numtarget[i] = _procMpcToMpc->numConnect();
		}
		else {
			size[i] = 0;
			numtarget[i] = 0;
		}
	}
#ifdef DISTRIBUTED
	communicator->globalSum(size);
	communicator->globalSum(numtarget);
#endif

	int totSize = 0;
	size_t totNumtarget = 0;
	for(i=0; i<numCPU; ++i) {
		totSize += size[i];
		totNumtarget += numtarget[i];
	}

	std::vector<int> pointer(totSize+numCPU);
	std::vector<int> target(totNumtarget);
	int startp = 0;
	size_t startt = 0;
	for(i=0; i<pid; ++i) {
		startp += size[i]+1;
		startt += numtarget[i];
	}
	for(i=0; i<totSize+numCPU; ++i) pointer[i] = 0;
	for(i=0; i<totNumtarget; ++i) target[i] = 0;
	for(i=0; i<size[pid]; ++i)
		pointer[startp + i] = _procMpcToMpc->offset(i);
	pointer[startp + size[pid]] = _procMpcToMpc->numConnect();
	for(i=0; i<numtarget[pid]; ++i)
		target[startt + i] = _procMpcToMpc->getTargetValue(i);
#ifdef DISTRIBUTED
	communicator->globalSum(pointer);
	communicator->globalSum(target);
#endif

	// now all the _procMpcToMpc connectivity data is stored in size, pointer and target
	std::vector<Connectivity> tmpMpcToMpc;
	tmpMpcToMpc.reserve(numCPU);
	for(i=0; i<numCPU; ++i) {
		startp = 0; startt = 0;
		for(j=0; j<i; ++j) {
			startp += size[j]+1;
			startt += numtarget[j];
		}
		tmpMpcToMpc.emplace_back(size[i], &pointer[startp], &target[startt], 0);
	}
	// now each processor has the _procMpcToMpc connectivities for all other processors
	Connectivity subToCpu = cpuToSub->reverse();
	mpcToCpu = std::make_unique<Connectivity>( mpcToSub_dual->transcon(subToCpu) );
	int csize = mpcToCpu->csize();
	std::vector<int> flags(csize);
	for(i=0; i<csize; ++i) flags[i] = -1;
	// identify the number of connections from MPC i
	std::vector<size_t> np(csize+1);
	size_t cp = 0;
	for(i=0; i<csize; ++i) {
		np[i] = cp;
		for(j=0; j < mpcToCpu->num(i); ++j) {
			int cpu = (*mpcToCpu)[i][j];
			for(k=0; k < tmpMpcToMpc[cpu].num(i); ++k) {
				int mpck = tmpMpcToMpc[cpu][i][k];
				if(flags[mpck] != i) {
					cp++;
					flags[mpck] = i;
				}
			}
		}
	}
	np[csize] = cp;

	// Now allocate and fill the new target
	for(i=0; i<csize; ++i) flags[i] = -1;
	std::vector<int> ntg(cp);
	cp = 0;
	for(i=0; i<csize; ++i) {
		np[i] = cp;
		for(j=0; j < mpcToCpu->num(i); ++j) {
			int cpu = (*mpcToCpu)[i][j];
			for(k=0; k < tmpMpcToMpc[cpu].num(i); ++k) {
				int mpck = tmpMpcToMpc[cpu][i][k];
				if(flags[mpck] != i) {
					ntg[cp] = mpck;
					cp++;
					flags[mpck] = i;
				}
			}
		}
	}

#ifdef USE_MUMPS
	if(domain->solInfo().solvercntl->fetiInfo.cct_cntl->subtype == FetiInfo::mumps && domain->solInfo().solvercntl->fetiInfo.cct_cntl->mumps_icntl[18] == 3) {
    procMpcToMpc = _procMpcToMpc;
  } else
#endif
	delete _procMpcToMpc;
	mpcToMpc = std::make_unique<Connectivity>(csize, std::move(np), std::move(ntg));
}

template<class Scalar>
void GenDecDomain<Scalar>::addFsiElements()
{
 if ( domain->solInfo().isCoupled && isFeti(domain->solInfo().solvercntl->type) &&
      domain->solInfo().isMatching && domain->solInfo().solvercntl->fetiInfo.fsi_corner != 0 ) {
   // JLchange: replace addSubFsiElem() such that fsi elements are added only to structure elements. 
   for (int i=0; i< domain->getNumFSI(); ++i) {
     LMPCons *thisGlFSI = domain->getFsi(i);
     int glFluidNode = thisGlFSI->lmpcnum;
     for(int j=0; j< thisGlFSI->nterms; j++) {
       int glStrucNode = thisGlFSI->terms[j].nnum;
       // Find all subdomains which have a structure element covering these two nodes 

       int weight = 0;
       for (int fSub=0; fSub< nodeToSub->num(glFluidNode); ++fSub) {
         int thisSubNum = (*nodeToSub)[glFluidNode][fSub];
         if (nodeToSub->locate(glStrucNode,thisSubNum)) {
           bool foundOne = false;
           for (int iElem=0; iElem< (domain->nodeToElem)->num(glFluidNode); ++iElem) {
             if (foundOne == true) break;
             int thisElemNum = (*(domain->nodeToElem))[glFluidNode][iElem];
             if ((domain->nodeToElem)->locate(glStrucNode,thisElemNum) &&
                 domain->isStructureElement(thisElemNum) && subToElem->locate(thisSubNum,thisElemNum)) {
               weight++;
               foundOne = true;
             }
           }
         }
       }

       for (int fSub=0; fSub< nodeToSub->num(glFluidNode); ++fSub) {
         int thisSubNum = (*nodeToSub)[glFluidNode][fSub];
         int localSubNumber = glSubToLocal[thisSubNum];
         if ((localSubNumber >= 0) && nodeToSub->locate(glStrucNode,thisSubNum)) {
           bool foundOne = false;
           for (int iElem=0; iElem< (domain->nodeToElem)->num(glFluidNode); ++iElem) {
             if (foundOne == true) break;
             int thisElemNum = (*(domain->nodeToElem))[glFluidNode][iElem];
             if ((domain->nodeToElem)->locate(glStrucNode,thisElemNum) &&
                 domain->isStructureElement(thisElemNum) && subToElem->locate(thisSubNum,thisElemNum)) {
               LMPCons *localFsi = new LMPCons(glFluidNode, 0.0);
               LMPCTerm thisLmpcTerm(thisGlFSI->terms[j], 1.0/weight);
               localFsi->addterm(&(thisLmpcTerm));
               subDomain[localSubNumber]->addSingleFsi(localFsi);
               foundOne = true;
             }
           }
         }
       }
     }
   }
 }
}

template<class Scalar>
void GenDecDomain<Scalar>::preProcessFSIs()
{
  if(!domain->solInfo().isCoupled) return;
  if(verboseFlag) filePrint(stderr, " ... Applying the Fluid-Structure Interactions");
  domain->computeCoupledScaleFactors();
  domain->makeFsiToNode();
  wetInterfaceNodes = domain->getAllWetInterfaceNodes(numWetInterfaceNodes);
  distributeWetInterfaceNodes();
  if(verboseFlag) filePrint(stderr, " ...\n");
}

template<class Scalar>
void GenDecDomain<Scalar>::distributeWetInterfaceNodes()
{
  int i, iSub, subI;
  int *nWetInterfaceNodesPerSub = new int[numSub];
  for (iSub = 0; iSub < numSub; ++iSub)  {
    nWetInterfaceNodesPerSub[iSub] = 0;
  }
  for(i = 0; i < numWetInterfaceNodes; ++i) {
    int node = wetInterfaceNodes[i];
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub){
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) 
        nWetInterfaceNodesPerSub[subI]++;
    }
  }

  int **subWetInterfaceNodes = new int *[numSub];
  for(iSub = 0; iSub < numSub; ++iSub) {
    subWetInterfaceNodes[iSub] = new int[nWetInterfaceNodesPerSub[iSub]];
    nWetInterfaceNodesPerSub[iSub] = 0;
  }
  for(i = 0; i < numWetInterfaceNodes; ++i) {
    int node = wetInterfaceNodes[i];
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) < 0) continue;
      subWetInterfaceNodes[subI][ nWetInterfaceNodesPerSub[subI] ] = node;  // global number
      nWetInterfaceNodesPerSub[subI]++;
    }
  }

  if ((domain->solInfo().isCoupled) && (domain->solInfo().getFetiInfo().fsi_corner != 0)) {
    execParal(numSub, this, &GenDecDomain<Scalar>::markSubWetInterface,
              nWetInterfaceNodesPerSub, subWetInterfaceNodes);
  }  

  execParal(numSub, this, &GenDecDomain<Scalar>::setSubWetInterface, nWetInterfaceNodesPerSub, subWetInterfaceNodes);
  delete [] nWetInterfaceNodesPerSub;
  for(iSub = 0; iSub < numSub; ++iSub) delete [] subWetInterfaceNodes[iSub];
  delete [] subWetInterfaceNodes;
}

template<class Scalar>
void GenDecDomain<Scalar>::markSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes)
{
  subDomain[iSub]->markWetInterface(nWetInterfaceNodesPerSub[iSub], subWetInterfaceNodes[iSub], soweredInput);
}

template<class Scalar>
void GenDecDomain<Scalar>::setSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes)
{
  subDomain[iSub]->setWetInterface(nWetInterfaceNodesPerSub[iSub], subWetInterfaceNodes[iSub], soweredInput);
}

template<class Scalar>
void GenDecDomain<Scalar>::getSharedFSIs()
{
  if(domain->solInfo().isCoupled) {
    Connectivity fsiToSub = domain->getFsiToNode()->transcon(*nodeToSub);
    Connectivity *subToFsi = fsiToSub.altReverse(); // reverse without reordering
    Connectivity subToSub_fsi = subToFsi->transcon(fsiToSub);
    paralApply(subDomain, &BaseSub::makeFsiInterface, subToFsi, fsiToSub, &subToSub_fsi);
    delete subToFsi;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeMpcToSub()
{
  if(numDualMpc) {
    int *pointer = new int[globalNumSub+1];
    for(int i =0; i < globalNumSub; ++i) pointer[i] = 0;
    paralApply(subDomain, &BaseSub::putNumMPC, pointer);
#ifdef DISTRIBUTED
    communicator->globalSum(globalNumSub, pointer);
#endif

    int total = 0;
    for(int i = 0; i < globalNumSub; ++i) {
      int tmp = pointer[i];
      pointer[i] = total;
      total += tmp;
    }
    pointer[globalNumSub] = total;

    int *target = new int[total];
    for(int i = 0; i < total; ++i) target[i] = 0;
    paralApply(subDomain, &BaseSub::putLocalToGlobalMPC, pointer, target);
#ifdef DISTRIBUTED
    communicator->globalSum(total, target);
#endif

    auto subToMpc = std::make_unique<Connectivity>(globalNumSub, pointer, target);
    mpcToSub_dual.reset(subToMpc->alloc_reverse());
  }
  if(numPrimalMpc) {
    int *pointer = new int[globalNumSub+1];
    for(int i =0; i < globalNumSub; ++i) pointer[i] = 0;
    paralApply(subDomain, &BaseSub::putNumMPC_primal, pointer);
#ifdef DISTRIBUTED
    communicator->globalSum(globalNumSub, pointer);
#endif

    int total = 0;
    for(int i = 0; i < globalNumSub; ++i) {
      int tmp = pointer[i];
      pointer[i] = total;
      total += tmp;
    }
    pointer[globalNumSub] = total;

    int *target = new int[total];
    for(int i = 0; i < total; ++i) target[i] = 0;
    paralApply(subDomain, &BaseSub::putLocalToGlobalMPC_primal, pointer, target);
#ifdef DISTRIBUTED
    communicator->globalSum(total, target);
#endif

    Connectivity *subToMpc = new Connectivity(globalNumSub, pointer, target);
    if(mpcToSub_primal) delete mpcToSub_primal;
    mpcToSub_primal = subToMpc->alloc_reverse();
    delete subToMpc;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::buildOps(GenMDDynamMat<Scalar> &res, double coeM, double coeC, double coeK,
                               Rbm **rbms, FullSquareMatrix **kelArray, bool make_feti,
                               FullSquareMatrix **melArray, FullSquareMatrix **celArray, bool factor)
{
 GenDomainGroupTask<Scalar> dgt(numSub, subDomain.data(), coeM, coeC, coeK, rbms, kelArray,
                                domain->solInfo().alphaDamp, domain->solInfo().betaDamp,
                                domain->numSommer, domain->solInfo().getFetiInfo().local_cntl->subtype,
                                communicator, melArray, celArray, mt);
 if(domain->solInfo().solvercntl->type == SolverSelection::Direct) {
   switch(domain->solInfo().solvercntl->subtype) {
#ifdef USE_MUMPS
     case 9 : {
       if(soweredInput) { filePrint(stderr, " *** ERROR: MUMPS solver is not supported with binary input. Exiting...\n"); exit(-1); }
       GenMumpsSolver<Scalar> *msmat;
       std::map<int,int>::iterator it = domain->solInfo().solvercntl->mumps_icntl.find(18);
       if(it != domain->solInfo().solvercntl->mumps_icntl.end() && it->second == 3) {
         Connectivity *subToCpu = cpuToSub->alloc_reverse();
         Connectivity *elemToCpu = elemToSub->transcon(subToCpu);
         Connectivity *nodeToNodeLocal = domain->getNodeToElem()->transcon(elemToNode.get(), elemToCpu->getTarget(), communicator->cpuNum());
         msmat = new GenMumpsSolver<Scalar>(nodeToNodeLocal, domain->getDSA(), domain->getCDSA(), numSub, subDomain, 
                                            *domain->solInfo().solvercntl, communicator);
         delete nodeToNodeLocal;
         delete elemToCpu;
         delete subToCpu;
       }
       else {
         msmat = new GenMumpsSolver<Scalar>(domain->getNodeToNode(), domain->getDSA(), domain->getCDSA(), numSub, subDomain,
                                            *domain->solInfo().solvercntl, communicator);
       }
       for(int i = 0; i < numSub; ++i) {
         dgt.dynMats[i] = dynamic_cast<GenSolver<Scalar>* >(msmat);
         dgt.spMats[i] = dynamic_cast<GenSparseMatrix<Scalar>* >(msmat);
       }
     } break;
#endif
#ifdef USE_EIGEN3
     case 13 : { // eisgal for implicit POD with parallel assembly and serial direct inversion
       // do nothing
     } break;
#endif
     default :
       filePrint(stderr, " *** ERROR: subtype %d not supported here in GenDecDomain::buildOps\n", domain->solInfo().solvercntl->subtype);
       exit(-1);
   }
 }

 if(verboseFlag) filePrint(stderr," ... Assemble Subdomain Matrices    ... \n");
 execParal(numSub, &dgt, &GenDomainGroupTask<Scalar>::runForWB, make_feti);

 GenAssembler<Scalar> * assembler = 0;
 if(domain->solInfo().inpc || domain->solInfo().aeroFlag > -1 || domain->solInfo().solvercntl->type == SolverSelection::Iterative) {
   assembler = getSolVecAssembler(); 
 }
// RT0212
{
   FSCommPattern<Scalar> *pat = new FSCommPattern<Scalar>(communicator, cpuToSub.get(), myCPU,
                                                          FSCommPattern<Scalar>::CopyOnSend);
   for(int i=0; i<numSub; ++i) subDomain[i]->setDofPlusCommSize(pat);
   pat->finalize();
   ba2 = new GenBasicAssembler<Scalar>(numSub, subDomain.data(), pat);
}
// RT0212


 if(domain->solInfo().inpc) res.K = new GenSubDOp<Scalar>(numSub, dgt.K, assembler);
 else res.K = new GenSubDOp<Scalar>(numSub, dgt.K);
 res.Kuc = new GenSubDOp<Scalar>(numSub, dgt.Kuc);

 if(dgt.makeC) {
   res.C = new GenSubDOp<Scalar>(numSub, dgt.C);
   res.Cuc = new GenSubDOp<Scalar>(numSub, dgt.Cuc);
   res.Ccc = new GenSubDOp<Scalar>(numSub, dgt.Ccc);
 }
 else {
   res.C   = 0; delete [] dgt.C;
   res.Cuc = 0; delete [] dgt.Cuc;
   res.Ccc = 0; delete [] dgt.Ccc;
 }
 res.M   = new GenSubDOp<Scalar>(numSub, dgt.M);
 res.Muc = new GenSubDOp<Scalar>(numSub, dgt.Muc);
 res.Mcc = new GenSubDOp<Scalar>(numSub, dgt.Mcc);

// RT
 if(dgt.makeC_deriv) {
   res.C_deriv = new GenSubDOp<Scalar>*[1];
   (res.C_deriv)[0] = new GenSubDOp<Scalar>(numSub, dgt.C_deriv,0);
   res.Cuc_deriv = new GenSubDOp<Scalar>*[1];
   res.Cuc_deriv[0] = new GenSubDOp<Scalar>(numSub, dgt.Cuc_deriv,0);
   res.num_K_deriv = dgt.num_K_deriv;
   res.K_deriv = new GenSubDOp<Scalar>*[res.num_K_deriv+1];
   for(int i=0;i<=res.num_K_deriv;i++)
     res.K_deriv[i] = new GenSubDOp<Scalar>(numSub, dgt.K_deriv,i);
   for(int i=0;i<=res.num_K_deriv;i++) for(int j=0;j<numSub;j++) {
      GenSparseMatrix<Scalar>* p = (*(res.K_deriv[i]))[j];
   }
   if (dgt.Kuc_deriv) {
      res.Kuc_deriv = new GenSubDOp<Scalar>*[res.num_K_deriv+1];
      for(int i=0;i<=res.num_K_deriv;i++)
       res.Kuc_deriv[i] = new GenSubDOp<Scalar>(numSub, dgt.Kuc_deriv,i);
   } else { 
     res.Kuc_deriv = 0;
   }
   res.num_K_arubber = dgt.num_K_arubber;
   res.K_arubber_l = new GenSubDOp<Scalar>*[res.num_K_arubber];
   res.K_arubber_m = new GenSubDOp<Scalar>*[res.num_K_arubber];
   for(int i=0;i<res.num_K_arubber;i++) {
     res.K_arubber_l[i] = new GenSubDOp<Scalar>(numSub, dgt.K_arubber_l,i);
     res.K_arubber_m[i] = new GenSubDOp<Scalar>(numSub, dgt.K_arubber_m,i);
   }
   if (dgt.Kuc_arubber_l) {
      res.Kuc_arubber_l = new GenSubDOp<Scalar>*[res.num_K_arubber];
      res.Kuc_arubber_m = new GenSubDOp<Scalar>*[res.num_K_arubber];
      for(int i=0;i<res.num_K_arubber;i++) {
       res.Kuc_arubber_l[i] = new GenSubDOp<Scalar>(numSub, dgt.Kuc_arubber_l,i);     
       res.Kuc_arubber_m[i] = new GenSubDOp<Scalar>(numSub, dgt.Kuc_arubber_m,i);     
      }
   } else { 
     res.Kuc_arubber_l = 0;
     res.Kuc_arubber_m = 0;
   }
 } else {
   res.C_deriv = 0;
   if (dgt.C_deriv) delete [] dgt.C_deriv;
   res.Cuc_deriv = 0;
   if (dgt.Cuc_deriv) delete [] dgt.Cuc_deriv;
   res.K_deriv = 0;
   if (dgt.K_deriv) delete [] dgt.K_deriv;
   res.Kuc_deriv = 0;
   if (dgt.Kuc_deriv) delete [] dgt.Kuc_deriv;
   res.K_arubber_l = 0;
   if (dgt.K_arubber_l) delete [] dgt.K_arubber_l;
   res.Kuc_arubber_l = 0;
   if (dgt.Kuc_arubber_l) delete [] dgt.Kuc_arubber_l;
   res.K_arubber_m = 0;
   if (dgt.K_arubber_m) delete [] dgt.K_arubber_m;
   res.Kuc_arubber_m = 0;
   if (dgt.Kuc_arubber_m) delete [] dgt.Kuc_arubber_m;
 }
// RT end
 switch(domain->solInfo().solvercntl->type) {
   case SolverSelection::Direct : {
     if(domain->solInfo().solvercntl->subtype == 13) { // if ROM, dont need dgt.dynMats, keep dgt.spMats for V^T*K*V, factor after pod basis has been buffered
       if(!(res.spMat  = new GenSubDOp<Scalar>(numSub, dgt.spMats)))
         throw std::runtime_error("subdomain matrices did not bind");
       res.spMat->zeroAll();
#ifdef USE_EIGEN3
       res.dynMat = 
           new Rom::GenEiSparseGalerkinProjectionSolver<Scalar,GenDistrVector,GenParallelSolver<Scalar> >
              (domain->getNodeToNode(), domain->getDSA(), domain->getCDSA(), 
               numSub, dgt.spMats, !domain->solInfo().unsym(), 1e-6, 
               domain->solInfo().numberOfRomCPUs);
#endif
       delete [] dgt.dynMats;
     } else { // else mumps
       if(dgt.dynMats[0]) {
         dgt.dynMats[0]->unify(communicator);
         res.dynMat = dynamic_cast<GenParallelSolver<Scalar>* >(dgt.dynMats[0]);
       }
       delete [] dgt.dynMats;
       delete [] dgt.spMats;
     }
     if(factor) res.dynMat->refactor();
   } break;
   case SolverSelection::Iterative : {
     switch(domain->solInfo().solvercntl->iterType) {
       case 1: {
         if(myCPU == 0) std::cerr << " ... GMRES Solver is Selected       ... \n";
         res.spMat = new GenSubDOp<Scalar>(numSub, dgt.spMats, assembler);
         if(domain->solInfo().solvercntl->precond == 1) res.prec = getDiagSolver(numSub, dgt.sd, dgt.sps);
         GmresSolver<Scalar, GenDistrVector<Scalar>, GenSubDOp<Scalar>, GenParallelSolver<Scalar>, GenParallelSolver<Scalar> > *gmresSolver
           = new GmresSolver<Scalar, GenDistrVector<Scalar>, GenSubDOp<Scalar>, GenParallelSolver<Scalar>, GenParallelSolver<Scalar> >
             (domain->solInfo().solvercntl->maxit, domain->solInfo().solvercntl->tol, res.spMat, &GenSubDOp<Scalar>::mult, res.prec,
              &GenParallelSolver<Scalar>::solve, NULL, &GenParallelSolver<Scalar>::solve, communicator); 
         if(domain->solInfo().solvercntl->maxvecsize > 0) gmresSolver->maxortho = domain->solInfo().solvercntl->maxvecsize;
         gmresSolver->verbose = verboseFlag;
         gmresSolver->printNumber = domain->solInfo().solvercntl->fetiInfo.printNumber;
         res.dynMat = gmresSolver;
       } break;
       default:
         std::cerr << " *** ERROR: iterType " << domain->solInfo().solvercntl->iterType << " not supported here in GenDecDomain::buildOps\n";
         exit(-1);
     }
   } break;
   case SolverSelection::Feti : {
     if(myCPU == 0) std::cerr << " ... FETI-DP Solver is Selected     ... \n";
     res.dynMat = getFetiSolver(dgt);
     delete [] dgt.dynMats;
     delete [] dgt.spMats;
   } break;
   case SolverSelection::BlockDiag : {
     if(myCPU == 0) std::cerr << " ... Diagonal Solver is Selected    ... \n";
     res.dynMat = getDiagSolver(numSub, dgt.sd, dgt.dynMats);
   } break;
	 case SolverSelection::FetiLib : {
		 if(myCPU == 0) std::cerr << " ... FetiLib Solver is Selected     ... \n";
		 res.dynMat = getFetiSolver(dgt);
		 delete [] dgt.dynMats;
		 delete [] dgt.spMats;
	 }
	 break;
 }
}

template<class Scalar>
void
GenDecDomain<Scalar>::rebuildOps(GenMDDynamMat<Scalar> &res, double coeM, double coeC, double coeK, 
                                 FullSquareMatrix **kelArray, FullSquareMatrix **melArray, FullSquareMatrix **celArray)
{
 if(res.dynMat) res.dynMat->reconstruct(); // do anything that needs to be done before zeroing and assembling the matrices

 execParal(numSub, this, &GenDecDomain<Scalar>::subRebuildOps, res, coeM, coeC, coeK, kelArray, melArray, celArray);

 if(domain->solInfo().solvercntl->type == SolverSelection::Direct && res.dynMat) {
   GenSolver<Scalar> *dynmat = dynamic_cast<GenSolver<Scalar>*>(res.dynMat);
   if(!verboseFlag) dynmat->setPrintNullity(false);
   dynmat->unify(communicator);
 }
 if(res.dynMat) res.dynMat->refactor(); // do anything that needs to be done after zeroing and assembling the matrices
}


template<>
void
GenDecDomain<double>::rebuildOps(GenMDDynamMat<double> &res, double coeM, double coeC, double coeK, 
                                 FullSquareMatrix **kelArray, FullSquareMatrix **melArray, FullSquareMatrix **celArray)
{
 if(domain->solInfo().galerkinPodRom && domain->solInfo().newmarkBeta != 0) {
#ifdef USE_EIGEN3

   // rebuild the dynamic stiffness matrix
   if (domain->solInfo().useMassNormalizedBasis || domain->solInfo().modalDIMASS)
     execParal(numSub, this, &GenDecDomain<double>::subRebuildOps, res, 0.0, coeC, coeK, kelArray, melArray, celArray);
   else 
     execParal(numSub, this, &GenDecDomain<double>::subRebuildOps, res, coeM, coeC, coeK, kelArray, melArray, celArray);

   // every process now has the reduced dynamic stiffness matrix, add it to the the galerkin solver
   if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().modalDIMASS)
     dynamic_cast<Rom::GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *>(res.sysSolver)->addReducedMass(coeM);

#endif
 } else {
   if(res.dynMat) res.dynMat->reconstruct(); // do anything that needs to be done before zeroing and assembling the matrices

   execParal(numSub, this, &GenDecDomain<double>::subRebuildOps, res, coeM, coeC, coeK, kelArray, melArray, celArray);

   if(domain->solInfo().solvercntl->type == SolverSelection::Direct && res.dynMat) {
     GenSolver<double> *dynmat = dynamic_cast<GenSolver<double>*>(res.dynMat);
     if(!verboseFlag) dynmat->setPrintNullity(false);
     dynmat->unify(communicator);
   }
 }

 if(res.dynMat) res.dynMat->refactor(); // do anything that needs to be done after zeroing and assembling the matrices
}

template<class Scalar>
void
GenDecDomain<Scalar>::subRebuildOps(int iSub, GenMDDynamMat<Scalar> &res, double coeM, double coeC, double coeK, 
                                    FullSquareMatrix **kelArray, FullSquareMatrix **melArray, FullSquareMatrix **celArray)
{
  AllOps<Scalar> allOps;

  if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
   allOps.K = new GenMultiSparse<Scalar>((res.K) ? (*res.K)[iSub] : 0, subDomain[iSub]->KiiSparse.get(),
                                         subDomain[iSub]->Kbb.get(), subDomain[iSub]->Kib.get());
  else
    allOps.K = (res.K) ? (*res.K)[iSub] : 0;
  if(res.C)    allOps.C   = (*res.C)[iSub];
  if(res.Cuc)  allOps.Cuc = (*res.Cuc)[iSub];
  if(res.Ccc)  allOps.Ccc = (*res.Ccc)[iSub];
  if(res.M)    allOps.M   = (*res.M)[iSub];
  if(res.Muc)  allOps.Muc = (*res.Muc)[iSub];
  if(res.Mcc)  allOps.Mcc = (*res.Mcc)[iSub];
  if(res.Kuc)  allOps.Kuc = (*res.Kuc)[iSub];
// RT: 053013 : not finished
  if(res.C_deriv) {

     allOps.C_deriv = new GenSparseMatrix<Scalar>*[1];
     allOps.C_deriv[0] = (*(res.C_deriv[0]))[iSub];
  }
  if(res.Cuc_deriv) {
     allOps.Cuc_deriv = new GenSparseMatrix<Scalar>*[1];
     allOps.Cuc_deriv[0] = (*(res.Cuc_deriv[0]))[iSub];
  }
  if(res.K_deriv) {
     allOps.K_deriv = new GenSparseMatrix<Scalar>*[res.num_K_deriv+1];
     allOps.n_Kderiv = res.num_K_deriv;
     for(int i=0;i<=res.num_K_deriv;i++)
        allOps.K_deriv[i] = (*(res.K_deriv[i]))[iSub];
  }
  if(res.Kuc_deriv) {
     allOps.Kuc_deriv = new GenSparseMatrix<Scalar>*[res.num_K_deriv+1];
     for(int i=0;i<=res.num_K_deriv;i++)
       allOps.Kuc_deriv[i] = (*(res.Kuc_deriv[i]))[iSub];
  }
  if(res.K_arubber_l) {
     allOps.K_arubber_l = new GenSparseMatrix<Scalar>*[res.num_K_arubber];
     allOps.K_arubber_m = new GenSparseMatrix<Scalar>*[res.num_K_arubber];
     allOps.num_K_arubber = res.num_K_arubber;
     for(int i=0;i<res.num_K_arubber;i++) {
        allOps.K_arubber_l[i] = (*(res.K_arubber_l[i]))[iSub];
        allOps.K_arubber_m[i] = (*(res.K_arubber_m[i]))[iSub];
     }
  }
  if(res.Kuc_arubber_l) {
     allOps.Kuc_arubber_l = new GenSparseMatrix<Scalar>*[res.num_K_arubber];
     allOps.Kuc_arubber_m = new GenSparseMatrix<Scalar>*[res.num_K_arubber];
     for(int i=0;i<res.num_K_arubber;i++) {
       allOps.Kuc_arubber_l[i] = (*(res.Kuc_arubber_l[i]))[iSub];
       allOps.Kuc_arubber_m[i] = (*(res.Kuc_arubber_m[i]))[iSub];
     }
  }

  allOps.zero();

  if(domain->solInfo().solvercntl->type == SolverSelection::Direct) {
    if(domain->solInfo().solvercntl->subtype == 13) { // if performing reduced order model, point to subdomain matrices and pass these to makeSparseOps
#ifdef USE_EIGEN3
      GenSparseMatrix<Scalar> *spmat = 
      dynamic_cast<Rom::GenEiSparseGalerkinProjectionSolver<Scalar,GenDistrVector,GenParallelSolver<Scalar> > *>(res.dynMat)->getSpMat(iSub);
      dynamic_cast<Rom::GenEiSparseGalerkinProjectionSolver<Scalar,GenDistrVector,GenParallelSolver<Scalar> > *>(res.dynMat)->zeroAll();
      subDomain[iSub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, spmat, (kelArray) ? kelArray[iSub] : 0,
                                                     (melArray) ? melArray[iSub] : 0, (celArray) ? celArray[iSub] : 0);
#endif
    } else {
      GenSparseMatrix<Scalar> *spmat = (res.dynMat) ? dynamic_cast<GenSparseMatrix<Scalar>*>(res.dynMat) : NULL;
      if(iSub == 0 && spmat) spmat->zeroAll();
#if defined(_OPENMP)
    #pragma omp barrier
#endif
      subDomain[iSub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, spmat, (kelArray) ? kelArray[iSub] : 0,
                                                     (melArray) ? melArray[iSub] : 0, (celArray) ? celArray[iSub] : 0);
    }
  }
  else if(domain->solInfo().solvercntl->type == SolverSelection::BlockDiag) {
    GenSparseMatrix<Scalar> *spmat = (res.dynMat) ? dynamic_cast<GenSparseMatrix<Scalar>*>(res.dynMat) : NULL;
    if(spmat) spmat->zeroAll();
    subDomain[iSub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, spmat, (kelArray) ? kelArray[iSub] : 0,
                                                    (melArray) ? melArray[iSub] : 0, (celArray) ? celArray[iSub] : 0);
  }
  else {
    GenMultiSparse<Scalar> *allMats;
    if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
      allMats = new GenMultiSparse<Scalar>(subDomain[iSub]->KrrSparse, subDomain[iSub]->Krc.get(), subDomain[iSub]->Kcc.get());
    else 
      allMats = new GenMultiSparse<Scalar>(subDomain[iSub]->KrrSparse, subDomain[iSub]->KiiSparse.get(), subDomain[iSub]->Kbb.get(),
                                           subDomain[iSub]->Kib.get(), subDomain[iSub]->Krc.get(), subDomain[iSub]->Kcc.get());
    allMats->zeroAll();
    subDomain[iSub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, allMats, (kelArray) ? kelArray[iSub] : 0, 
                                                    (melArray) ? melArray[iSub] : 0);
    delete allMats;
  }

  if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted) delete allOps.K;
}

template<>
void
GenDecDomain<double>::buildFFP(GenDistrVector<double> &u, FILE *fffp, bool direction);

template<>
void
GenDecDomain<DComplex>::buildFFP(GenDistrVector<DComplex> &u, FILE *fffp, bool direction);

template<>
void
GenDecDomain<DComplex>::buildLocalFFP(int iSub, GenDistrVector<DComplex> *u,
                                      DComplex **ffp, int *numSample, double (*dir)[3], bool direction);

template<>
void
GenDecDomain<double>::buildLocalFFP(int iSub, GenDistrVector<double> *u,
                                    double **ffp, int *numSample, double (*dir)[3], bool direction);

template<class Scalar>
void
GenDecDomain<Scalar>::setNewProperties(int i)
{
  paralApply(subDomain, &Domain::setNewProperties, i);
}

template<class Scalar>
void
GenDecDomain<Scalar>::assignRandMat()
{
  paralApply(subDomain, &Domain::assignRandMat);
}

template<class Scalar>
void
GenDecDomain<Scalar>::retrieveElemset()
{
  paralApply(subDomain, &Domain::retrieveElemset);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain<Scalar>::getElementAttr(int fileNumber,int iAttr, double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  int numNodes = geoSource->numNode();

  double *weights = (double *) dbg_alloca(numNodes*sizeof(double));
  double *props = (double *) dbg_alloca(numNodes*sizeof(double));
  for(int i=0; i<numNodes; ++i) weights[i] = props[i] = 0.0;

  // ... WRITE THE TIME VALUE
  //filePrint(oinfo[fileNumber].filptr,"%20.10e\n",domain->solInfo().getTimeStep());
  // ... OUTPUT precision
  //int p = oinfo[fileNumber].precision;

  for(int i=0; i<this->getNumSub(); ++i)
    { this->getSubDomain(i)->mergeElemProps(props, weights, iAttr); }

  if(avgnum == 1 ) 
    {
      for(int k=0; k<numNodes; ++k)  
	{
	  if(weights[k] > 0)
	    { props[k] /= weights[k]; }
	}      
      geoSource->outputNodeScalars(fileNumber, props, numNodes, time);
    }
  else
    {
      // not implemented
      assert(0);
    }
  return;
}

template <>
void GenDecDomain<std::complex<double>>::setConstraintGap(DistrGeomState *geomState, DistrGeomState *refState,
                                            GenFetiSolver<std::complex<double>> *fetiSolver, double t) {
	throw "MLX Calling an unimplemented GenDecDomain<std::complex<double>>::setConstraintGap(...)";
}

template<class Scalar>
void GenDecDomain<Scalar>::setConstraintGap(DistrGeomState *geomState, DistrGeomState *refState,
                                            GenFetiSolver<Scalar> *fetiSolver, double t)
{
  // note: for nonlinear statics t is pseudo time (i.e. load factor)
  if(numDualMpc) {
    GenDistrVector<Scalar> cu(fetiSolver->interfInfo());
    GenDistrVector<Scalar> u(fetiSolver->localInfo());
    geomState->get_tot_displacement(u);
    ((GenFetiDPSolver<Scalar> *)fetiSolver)->multC(u, cu); // cu = C*u
    execParal(this->numSub, this, &GenDecDomain<Scalar>::setMpcRhs, cu, t, 0);
    if(domain->GetnContactSurfacePairs() && domain->solInfo().piecewise_contact) {
      GenDistrVector<Scalar> u0(fetiSolver->localInfo());
      refState->get_tot_displacement(u0);
      u -= u0;
      ((GenFetiDPSolver<Scalar> *)fetiSolver)->multC(u, cu);
      execParal(this->numSub, this, &GenDecDomain<Scalar>::setMpcRhs, cu, t, 1);
    }
  }
}


template <>
void GenDecDomain<std::complex<double>>::setConstraintGap(DistrGeomState *geomState, DistrGeomState *refState,
                                                          FetiBaseClass<std::complex<double>> *fetiSolver, double t) {
	throw "MLX Calling an unimplemented GenDecDomain<std::complex<double>>::setConstraintGap(...)";
}

template<class Scalar>
void GenDecDomain<Scalar>::setConstraintGap(DistrGeomState *geomState, DistrGeomState *refState,
                                            FetiBaseClass<Scalar> *fetiSolver, double t)
{
	// note: for nonlinear statics t is pseudo time (i.e. load factor)
	if(numDualMpc) {
		GenDistrVector<Scalar> cu(fetiSolver->interfInfo());
		GenDistrVector<Scalar> u(fetiSolver->localInfo());
		geomState->get_tot_displacement(u);
		((GenFetiDPSolver<Scalar> *)fetiSolver)->multC(u, cu); // cu = C*u
		execParal(this->numSub, this, &GenDecDomain<Scalar>::setMpcRhs, cu, t, 0);
		if(domain->GetnContactSurfacePairs() && domain->solInfo().piecewise_contact) {
			GenDistrVector<Scalar> u0(fetiSolver->localInfo());
			refState->get_tot_displacement(u0);
			u -= u0;
			((GenFetiDPSolver<Scalar> *)fetiSolver)->multC(u, cu);
			execParal(this->numSub, this, &GenDecDomain<Scalar>::setMpcRhs, cu, t, 1);
		}
	}
}

template<>
void
GenDecDomain<std::complex<double>>::extractPosition(int iSub, DistrGeomState &geomState, GenDistrVector<std::complex<double>> &x)
{
	throw "MLX Unimplemented extractPosition.";
}

template<class Scalar>
void
GenDecDomain<Scalar>::extractPosition(int iSub, DistrGeomState &geomState, GenDistrVector<Scalar> &x)
{
  geomState[iSub]->extract(x.subData(subDomain[iSub]->localSubNum()));
}

template<class Scalar>
void
GenDecDomain<Scalar>::setMpcRhs(int iSub, GenDistrVector<Scalar> &cu, double t, int flag)
{
  subDomain[iSub]->setMpcRhs(cu.subData(subDomain[iSub]->localSubNum()), t, flag);
}

template<class Scalar>
void
GenDecDomain<Scalar>::printLMPC()
{
  for(int i=0; i<numSub; ++i) subDomain[i]->printLMPC();
}

template<class Scalar>
FSCommPattern<Scalar> *
GenDecDomain<Scalar>::getWiCommPattern()
{
  if(!wiPat) {
    wiPat = new FSCommPattern<Scalar>(communicator, cpuToSub.get(), myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                      FSCommPattern<Scalar>::NonSym);
    for(int iSub = 0; iSub < numSub; ++iSub) subDomain[iSub]->setWICommSize(wiPat);
    wiPat->finalize();
  }
  return wiPat;
}

template<class Scalar>
GenAssembler<Scalar> *
GenDecDomain<Scalar>::getSolVecAssembler() {
  if (!ba) {
    ba = solVecAssemblerNew(); 
  }
  return ba;
}

template<class Scalar>
GenBasicAssembler<Scalar> *
GenDecDomain<Scalar>::solVecAssemblerNew() {
  FSCommPattern<Scalar> *pat = new FSCommPattern<Scalar>(communicator, cpuToSub.get(), myCPU,
                                                         FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<numSub; ++i) {
    subDomain[i]->setDofPlusCommSize(pat);
  }
  pat->finalize();

  return new GenBasicAssembler<Scalar>(numSub, subDomain.data(), pat);
}

template<class Scalar>
void GenDecDomain<Scalar>::getEnergies(GenDistrVector<Scalar> &disp, GenDistrVector<Scalar> &extF, int fileNumber, double time,
                                       SysState<GenDistrVector<Scalar> > *distState, GenMDDynamMat<Scalar> *dynOps,
                                       GenDistrVector<Scalar> *aeroF)
{
  double Wext = 0.0, Waero = 0.0, Wdmp = 0.0, Wela = 0.0, Wkin = 0.0, error = 0.0;
  if(domain->solInfo().isDynam()) {
    double *subW = new double[6*numSub];
    execParal(numSub, this, &GenDecDomain<Scalar>::subGetEnergies, disp, extF, time, distState, dynOps, aeroF, subW);

    for(int i=0; i<numSub; ++i) {
      Wext  += subW[6*i  ];
      Waero += subW[6*i+1];
      Wdmp  += subW[6*i+2];
      Wela  += subW[6*i+3];
      Wkin  += subW[6*i+4];
      error += subW[6*i+5];
    }
    delete [] subW;
#ifdef DISTRIBUTED
    communicator->reduce(1, &Wext);
    communicator->reduce(1, &Waero);
    communicator->reduce(1, &Wdmp);
    communicator->reduce(1, &Wela);
    communicator->reduce(1, &Wkin);
    communicator->reduce(1, &error);
#endif
  }
  else {
    Wext = ScalarTypes::Real(extF * disp);
    Wela = 0.5 * Wext;
  }

#ifdef DISTRIBUTED
  if(myCPU == 0)
#endif
  geoSource->outputEnergies(fileNumber, time, Wext, Waero, Wela, Wkin, Wdmp, error);
}

template<class Scalar>
void GenDecDomain<Scalar>::getEnergies_b(DistrGeomState *geomState, GenDistrVector<Scalar> &extF,
                                       Corotator ***allCorot, int fileNumber, double time,
                                       SysState<GenDistrVector<Scalar> > *distState, GenMDDynamMat<Scalar> *dynOps,
                                       GenDistrVector<Scalar> *aeroF)
{
  double Wext = 0.0, Waero = 0.0, Wdmp = 0.0, Wela = 0.0, Wkin = 0.0, Wdis = 0.0, error = 0.0;
  double *subW = new double[7*numSub];
  execParal(numSub, this, &GenDecDomain<Scalar>::subGetEnergies_b, geomState, extF, allCorot, time,
              distState, dynOps, aeroF, subW);

  for(int i=0; i<numSub; ++i) {
    Wext  += subW[7*i  ];
    Waero += subW[7*i+1];
    Wdmp  += subW[7*i+2];
    Wela  += subW[7*i+3];
    Wkin  += subW[7*i+4];
    Wdis  += subW[7*i+5];
    error += subW[7*i+6];
  }
  delete [] subW;
#ifdef DISTRIBUTED
  communicator->reduce(1, &Wext);
  communicator->reduce(1, &Waero);
  communicator->reduce(1, &Wdmp);
  communicator->reduce(1, &Wela);
  communicator->reduce(1, &Wkin);
  communicator->reduce(1, &Wdis);
  communicator->reduce(1, &error);

  if(myCPU == 0)
#endif
  geoSource->outputEnergies(fileNumber, time, Wext, Waero, Wela, Wkin, Wdmp+Wdis, error);
}

template<>
inline void GenDecDomain<double>::subGetEnergies(int iSub, GenDistrVector<double> &disp, GenDistrVector<double> &extF, double time,
                                                 SysState<GenDistrVector<double> > *distState, GenMDDynamMat<double> *dynOps,
                                                 GenDistrVector<double> *aeroF, double *subW)
{
  GenStackVector<double> subDisp(disp.subData(iSub), disp.subLen(iSub));
  GenStackVector<double> subExtF(extF.subData(iSub), extF.subLen(iSub));
  GenStackVector<double> *subAeroF = (aeroF) ? new GenStackVector<double>(aeroF->subData(iSub), aeroF->subLen(iSub)) : NULL;
  GenStackVector<double> *subVel = (distState) ? new GenStackVector<double>(distState->getVeloc().subData(iSub), distState->getVeloc().subLen(iSub)) : NULL;
  GenSparseMatrix<double> *subK = (dynOps && dynOps->K) ? (*dynOps->K)[iSub] : NULL;
  GenSparseMatrix<double> *subM = (dynOps && dynOps->M) ? (*dynOps->M)[iSub] : NULL;
  GenSparseMatrix<double> *subC = (dynOps && dynOps->C) ? (*dynOps->C)[iSub] : NULL;
  subDomain[iSub]->computeEnergies(subDisp, subExtF, time, subAeroF, subVel, subK, subM, subC,
                                   subW[6*iSub+3], subW[6*iSub+4], subW[6*iSub+5]);
  subW[6*iSub+0] = subDomain[iSub]->getWext();
  subW[6*iSub+1] = subDomain[iSub]->getWaero();
  subW[6*iSub+2] = subDomain[iSub]->getWdmp();
  if(subAeroF) delete subAeroF;
  if(subVel) delete subVel;
}

template<>
inline void GenDecDomain<complex<double> >::subGetEnergies(int iSub, GenDistrVector<complex<double> > &, GenDistrVector<complex<double> > &,
                                                                double, SysState<GenDistrVector<complex<double> > > *,
                                                                GenMDDynamMat<complex<double> > *, GenDistrVector<complex<double> > *, double *subW)
{
  for(int i=0; i<6; ++i) subW[6*iSub+i] = 0;
}

template<>
inline void GenDecDomain<double>::subGetEnergies_b(int iSub, DistrGeomState *geomState, GenDistrVector<double> &extF,
                                                 Corotator ***allCorot, double time, SysState<GenDistrVector<double> > *distState,
                                                 GenMDDynamMat<double> *dynOps, GenDistrVector<double> *aeroF, double *subW)
{
  GenStackVector<double> subExtF(extF.subData(iSub), extF.subLen(iSub));
  GenStackVector<double> *subAeroF = (aeroF) ? new GenStackVector<double>(aeroF->subData(iSub), aeroF->subLen(iSub)) : NULL;
  double *subVel = (distState) ? distState->getVeloc().subData(iSub) : NULL;
  GenSparseMatrix<double> *subM = (dynOps && dynOps->M) ? (*dynOps->M)[iSub] : NULL;
  GenSparseMatrix<double> *subC = (dynOps && dynOps->C) ? (*dynOps->C)[iSub] : NULL;
  subDomain[iSub]->computeEnergies((*geomState)[iSub], subExtF, time, subAeroF, subVel, allCorot[iSub], subM, subC,
                                   subW[7*iSub+3], subW[7*iSub+4], subW[7*iSub+5], subW[7*iSub+6]);
  subW[7*iSub+0] = subDomain[iSub]->getWext();
  subW[7*iSub+1] = subDomain[iSub]->getWaero();
  subW[7*iSub+2] = subDomain[iSub]->getWdmp();
  if(subAeroF) delete subAeroF;
}

template<>
inline void GenDecDomain<complex<double> >::subGetEnergies_b(int iSub, DistrGeomState *, GenDistrVector<complex<double> > &,
                                                                Corotator ***, double, SysState<GenDistrVector<complex<double> > > *,
                                                                GenMDDynamMat<complex<double> > *, GenDistrVector<complex<double> > *, double *subW)
{
  for(int i=0; i<7; ++i) subW[7*iSub+i] = 0;
}

template<class Scalar>
void GenDecDomain<Scalar>::getDissipatedEnergy(DistrGeomState *geomState, Corotator ***allCorot, 
                                               int fileNumber, double time)
{
  int groupIndex = geoSource->getGroupNumber(fileNumber);
  if(groupIndex<0) { // output for all elements

    double *subD = new double[numSub];
    execParal(numSub, this, &GenDecDomain<Scalar>::subGetDissipatedEnergy, geomState, allCorot, subD);
    double D = 0;
    for(int i=0; i<numSub; ++i) D += subD[i];
#ifdef DISTRIBUTED
    communicator->reduce(1, &D);
    if(myCPU == 0)
#endif
    geoSource->outputEnergy(fileNumber, time, D);
    delete [] subD;
  }
  else {

    std::map<int, Group> &group = geoSource->group;
    int numAttributesInGrp = group[groupIndex].attributes.size();

    // initialize array to store data
    double** subDPerAttrib = new double*[numSub];
    for(int i=0; i<numSub; i++) subDPerAttrib[i] = new double[numAttributesInGrp];
    
    execParal(numSub, this, &GenDecDomain<Scalar>::subGetDissipEnergyPerAttributes, geomState, allCorot, groupIndex, subDPerAttrib);
    
    double D[numAttributesInGrp];
    for(int i=0; i<numAttributesInGrp; i++) {
      D[i] = 0;
      for(int j=0; j<numSub; j++)
        D[i] += subDPerAttrib[j][i];
    }
#ifdef DISTRIBUTED
    communicator->reduce(numAttributesInGrp, D);
    if(myCPU == 0)
#endif
    geoSource->outputEnergyPerAttribute(fileNumber, time, D, numAttributesInGrp);

    // free memory
    for(int i=0; i<numSub; i++) delete [] subDPerAttrib[i];
    delete [] subDPerAttrib;
    
  }

}

template<class Scalar>
void GenDecDomain<Scalar>::subGetDissipatedEnergy(int iSub, DistrGeomState *geomState, Corotator ***allCorot, double *D)
{
  D[iSub] = subDomain[iSub]->getDissipatedEnergy((*geomState)[iSub], allCorot[iSub]);
}

template<class Scalar>
void GenDecDomain<Scalar>::subGetDissipEnergyPerAttributes(int iSub, DistrGeomState *geomState, Corotator ***allCorot, int gIndex,
                                                           double **subD)
{
  subD[iSub] = subDomain[iSub]->getDissipatedEnergy((*geomState)[iSub], allCorot[iSub], gIndex);
}

template<class Scalar>
void
GenDecDomain<Scalar>::exchangeInterfaceGeomState(DistrGeomState *geomState)
{
	FSCommPattern<double> *geomStatePat = new FSCommPattern<double>(communicator, cpuToSub.get(),
	                                                                myCPU, FSCommPattern<Scalar>::CopyOnSend,
	                                                                FSCommPattern<Scalar>::NonSym);
  int len = (domain->solInfo().isDynam()) ? 28 : 16;
  for(int i=0; i<numSub; ++i) subDomain[i]->setNodeCommSize(geomStatePat, len);
  geomStatePat->finalize();

  execParal(numSub, this, &GenDecDomain<Scalar>::dispatchInterfaceGeomState, geomStatePat, geomState);
  geomStatePat->exchange();
  execParal(numSub, this, &GenDecDomain<Scalar>::collectInterfaceGeomState, geomStatePat, geomState);

  delete geomStatePat;
}

template<class Scalar>
void
GenDecDomain<Scalar>::dispatchInterfaceGeomState(int isub, FSCommPattern<double> *geomStatePat, DistrGeomState *geomState)
{
  if(domain->solInfo().isDynam())
    subDomain[isub]->dispatchInterfaceGeomStateDynam(geomStatePat, (*geomState)[isub]);
  else
    subDomain[isub]->dispatchInterfaceGeomState(geomStatePat, (*geomState)[isub]);
}

template<class Scalar>
void
GenDecDomain<Scalar>::collectInterfaceGeomState(int isub, FSCommPattern<double> *geomStatePat, DistrGeomState *geomState)
{
  if(domain->solInfo().isDynam())
    subDomain[isub]->collectInterfaceGeomStateDynam(geomStatePat, (*geomState)[isub]);
  else
    subDomain[isub]->collectInterfaceGeomState(geomStatePat, (*geomState)[isub]);
}

template<class Scalar>
void
GenDecDomain<Scalar>::clean()
{
  // this function should be used before re-calling preProcess() to prevent memory leaks
  if(ba2) delete ba2;
  for(auto sub : subDomain)
  	delete sub;
  subDomain.clear();
  elemToNode.reset();
  subToNode.reset();
  if(nodeToSub) {
    if(!nodeToSub_copy)
    	nodeToSub_copy = std::move(nodeToSub);
    nodeToSub.reset();
  }
  subToSub.reset();
  elemToSub.reset();
  delete glSubToLocal;
  glSubToLocal = nullptr;

  if(internalInfo) { vecInfoStore.push_back(internalInfo); internalInfo = 0; }
  if(internalInfo2) { vecInfoStore.push_back(internalInfo2); internalInfo2 = 0; }
  if(masterSolVecInfo_) { vecInfoStore.push_back(masterSolVecInfo_); masterSolVecInfo_ = 0; }
  if(nodeInfo) { vecInfoStore.push_back(nodeInfo); nodeInfo = 0; }
  if(nodeVecInfo) { vecInfoStore.push_back(nodeVecInfo); nodeVecInfo = 0; }
  if(eleVecInfo) { vecInfoStore.push_back(eleVecInfo); eleVecInfo = 0; }
  if(bcVecInfo) { vecInfoStore.push_back(bcVecInfo); bcVecInfo = 0; }
}

template<class Scalar>
void
GenDecDomain<Scalar>::assembleNodalInertiaTensors(FullSquareMatrix **melArray)
{
  FSCommPattern<double> *pat = new FSCommPattern<double>(communicator, cpuToSub.get(), myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<numSub; ++i) subDomain[i]->setNodeCommSize(pat, 9);
  pat->finalize();

  execParal(numSub, this, &GenDecDomain<Scalar>::dispatchInterfaceNodalInertiaTensors, pat, melArray);
  pat->exchange();
  execParal(numSub, this, &GenDecDomain<Scalar>::collectInterfaceNodalInertiaTensors, pat);

  delete pat;
}

template<class Scalar>
void
GenDecDomain<Scalar>::dispatchInterfaceNodalInertiaTensors(int isub, FSCommPattern<double> *pat, FullSquareMatrix **melArray)
{
  subDomain[isub]->assembleNodalInertiaTensors(melArray[isub]);
  subDomain[isub]->dispatchInterfaceNodalInertiaTensors(pat);
}

template<class Scalar>
void
GenDecDomain<Scalar>::collectInterfaceNodalInertiaTensors(int isub, FSCommPattern<double> *pat)
{
  subDomain[isub]->collectInterfaceNodalInertiaTensors(pat);
}

template<class Scalar>
MultiDomainRbm<Scalar> *
GenDecDomain<Scalar>::constructRbm(bool printFlag)
{
  if(printFlag) filePrint(stderr," ... Using Geometric RBM Method     ...\n");
  MultiDomainRbm<Scalar> *rbm = new MultiDomainRbm<Scalar>(this, domain->solInfo().tolsvd);
  if(printFlag)
    filePrint(stderr, " ... GRBM algorithm detected %d rigid body or zero energy modes ...\n", rbm->numRBM());
  return rbm;
}

