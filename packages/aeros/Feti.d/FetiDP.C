#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include <Driver.d/DecDomain.h>
#include <Driver.d/SubDomain.h>
#include <Feti.d/Feti.h>
#include <Feti.d/GMRESOrthoSet.h>
#include <Feti.d/GCROrthoSet.h>
#include <Feti.d/CGOrthoSet.h>
#include <Threads.d/PHelper.h>
#include <Math.d/matrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Solvers.d/Spooles.h>
#include <Utils.d/linkfc.h>
#include <Feti.d/CCtSolver.d/GlobalCCt.h>
#include <Feti.d/CCtSolver.d/BlockCCt.h>
#include <Feti.d/CCtSolver.d/SuperBlockCCt.h>
#include <Feti.d/CCtSolver.d/SubBlockCCt.h>
#include <Feti.d/FetiOp.h>
#include <Driver.d/Mpc.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/GetTime.h>
#include <Math.d/VecSet.h>
#include <Math.d/FullMatrix.h>
#include <Utils.d/print_debug.h>
#include <Solvers.d/SolverFactory.h>
#include <Element.d/MatrixElement.d/MatrixElement.h>
#include <Paral.d/MDDynam.h>
#include <Math.d/BLAS.h>

extern double t0;
extern double t1;
extern double t2;
extern double t4;
extern double t5;
extern double t6;

// number of main loop iterations + number of friction orthogonalisation loop iterations
extern int iterTotal;

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

#include <Utils.d/Memory.h>

inline double DABS(double x) { return (x>0.0) ? x : -x; }

// New constructor for both shared and distributed memory
template<class Scalar>
GenFetiDPSolver<Scalar>::GenFetiDPSolver(int _nsub, int _glNumSub, std::vector<FetiSub<Scalar> *> subdomains,
                                         const Connectivity *_subToSub, FetiInfo *_fetiInfo, FSCommunicator *_fetiCom,
                                         int *_glSubToLoc, const Connectivity *_mpcToSub,
                                         const Connectivity *_mpcToSub_primal,
                                         Connectivity *_mpcToMpc,
                                         const Connectivity *_mpcToCpu, const Connectivity *_cpuToSub,
                                         const Connectivity *_bodyToSub,
                                         std::vector<std::unique_ptr<GenSolver<Scalar>>> sysMatrices,
                                         GenSparseMatrix<Scalar> **sysSparse,
                                         Rbm **, bool _rbmFlag, bool _geometricRbms, int _verboseFlag)
	: FetiBaseClass<Scalar>(std::move(subdomains), threadManager->numThr(), _verboseFlag),
	  internalR(_nsub), internalC(_nsub), internalWI(_nsub)
{

	// Compute memory used by FETI Solver
	t6 -= getTime();
#ifdef DISTRIBUTED
	this->times.memoryFETI = -threadManager->memoryUsed();
#else
	this->times.memoryFETI -= memoryUsed();
#endif

	this->nsub       = _nsub;        // Number of subdomains
	this->glNumSub   = _glNumSub;
	this->subToSub   = _subToSub;    // subdomain to subdomain connectivity
	//std::cerr << "this->nsub = " << this->nsub << ", this->subToSub =\n"; this->subToSub->print();
	this->mpcToSub   = _mpcToSub;    // MPC to subdomain connectivity
	this->glNumMpc = (this->mpcToSub) ? this->mpcToSub->csize() : 0;
	this->mpcToSub_primal = _mpcToSub_primal;
	this->glNumMpc_primal = (this->mpcToSub_primal) ? this->mpcToSub_primal->csize() : 0;
	mpcToMpc   = _mpcToMpc;    // MPC to MPC connectivity (used for CC^t preconditioner)
	mpcToCpu   = _mpcToCpu;
	this->cpuToSub   = _cpuToSub;
	fetiInfo   = _fetiInfo;    // Feti solver information
	this->fetiCom    = _fetiCom;
	this->glSubToLoc = _glSubToLoc;

	globalFlagCtc = domain->getNumCTC();
#ifdef DISTRIBUTED
	globalFlagCtc = this->fetiCom->globalMax((int) globalFlagCtc);
#endif

	rbmFlag = _rbmFlag; // if this is true then check for singularities in Kcc^* and each Krr^{(s)}
	// and deal with the rigid body modes using projection method
	geometricRbms = _geometricRbms; // if this is true then use geometric method
	// to compute the rigid body modes of Kcc^* and Krr^{(s)}
	// when rbmFlag is true, otherwise use algebraic null space

	this->myCPU = this->fetiCom->cpuNum();
	this->numCPUs = this->fetiCom->size();
	bodyToSub = _bodyToSub;
	subToBody = bodyToSub->alloc_reverse();

	// Define FETI tolerance and maximum number of iterations
	double fetiTolerance = fetiInfo->tol;
	this->epsilon2 = fetiTolerance*fetiTolerance;
	this->maxiter  = fetiInfo->maxit;

	int iSub;
	// create this->vPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
	this->vPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
	for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setDofCommSize(this->vPat);
	this->vPat->finalize();

	// create this->sPat FSCommPattern objects, used to send/receive a single integer
	this->sPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
	for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setCommSize(this->sPat, 1);
	this->sPat->finalize();

	if(globalFlagCtc) {
		mpcPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
		for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setMpcCommSize(mpcPat);
		mpcPat->finalize();
	}

	// Classes to organize parallel execution of tasks
	this->fetiOps   = new GenFetiOp<Scalar> *[this->nsub];

	// Compute total this->interface length, total internal length
	// and total half this->interface length
	int tInterfLen    = 0;
	int tLocalLen     = 0;
	int tLocalRLen    = 0;
	this->halfSize = 0;
	for(iSub = 0; iSub < this->nsub; ++iSub) {
		this->fetiOps[iSub] = new GenFetiOp<Scalar>;
		this->interface.domLen[iSub] = this->subdomains[iSub]->interfLen();
		this->internalDI.domLen[iSub]  = this->subdomains[iSub]->getNumUncon();
		internalR.domLen[iSub] = this->subdomains[iSub]->localRLen();

		this->subdomains[iSub]->computeMasterFlag(*this->mpcToSub);
		this->fetiOps[iSub]->setHalfOffset(this->halfSize);
		this->halfSize      += this->subdomains[iSub]->halfInterfLen();
		tInterfLen    += this->interface.domLen[iSub];
		tLocalLen     += this->internalDI.domLen[iSub];
		tLocalRLen    += internalR.domLen[iSub];
	}
	this->interface.len = tInterfLen;
	this->internalDI.len  = tLocalLen;
	internalR.len = tLocalRLen;

	// compute the masterFlags
	bool *interfaceMasterFlag = new bool[tInterfLen];
	this->interface.computeOffsets();
	for(iSub = 0; iSub < this->nsub; ++iSub) {
		auto &subMasterFlag = this->subdomains[iSub]->getMasterFlag();
		int subOffset = this->interface.subOffset[iSub];
		int j;
		for(j=0; j<this->interface.domLen[iSub]; ++j)
			interfaceMasterFlag[subOffset+j] = subMasterFlag[j];
	}
	this->interface.setMasterFlag(interfaceMasterFlag);
	this->internalDI.setMasterFlag();
	internalR.setMasterFlag();
	// don't delete interfaceMasterFlag

	// Allocate space for reorthogonalization set
	this->times.memoryOSet -= memoryUsed();
	if((fetiInfo->outerloop == 0) || (fetiInfo->outerloop == 3))
		this->oSetCG = (fetiInfo->maxortho > 0) ? new GenCGOrthoSet<Scalar>(this->halfSize, fetiInfo->maxortho, this->fetiCom) : 0;
	else if(fetiInfo->outerloop == 1)
		this->oSetGMRES = new GenGMRESOrthoSet<Scalar>(this->halfSize, fetiInfo->maxortho, this->fetiCom);
	else
		this->oSetGCR = new GenGCROrthoSet<Scalar>(this->halfSize, fetiInfo->maxortho, this->fetiCom);
	this->times.memoryOSet += memoryUsed();

	if(sysMatrices.size() != 0) {
		for(iSub = 0; iSub < this->nsub; ++iSub) {
			this->subdomains[iSub]->Krr = std::move(sysMatrices[iSub]);
			this->subdomains[iSub]->KrrSparse = sysSparse[iSub];
			if(fetiInfo->printMatLab) {
				std::stringstream filename;
				filename << "localmat" << this->subdomains[iSub]->subNum();
				this->subdomains[iSub]->KrrSparse->printSparse(filename.str());
			}
			this->fetiOps[iSub]->setSysMatrix(this->subdomains[iSub]->Krr.get(), sysSparse[iSub]);
		}
	}

	if(verboseFlag) filePrint(stderr," ... Build Edge Augmentation (Q)    ... \n");
	computeLocalWaveNumbers();
	paralApplyToAll(this->nsub, this->subdomains, &FetiSub<Scalar>::makeQ);  // build augmentation matrix
	if(fetiInfo->augment == FetiInfo::Gs) {
		// exchange number of each neighbors rbms
		paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::sendNumNeighbGrbm, this->sPat);
		this->sPat->exchange();
		paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::recvNumNeighbGrbm, this->sPat);
	}

	// Compute stiffness scaling if required
	// this also perform the LMPCs stiffness scaling/splitting for the primal method
	paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::initScaling);
	if((fetiInfo->scaling == FetiInfo::kscaling) || ((fetiInfo->mpc_scaling == FetiInfo::kscaling) && (this->glNumMpc_primal > 0))
	   || (fetiInfo->augment == FetiInfo::WeightedEdges)) {
		execParal(this->nsub, this, &FetiBaseClass<Scalar>::sendScale);
		this->vPat->exchange();
		paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectScaling, this->vPat);
	}
	if(domain->solInfo().isCoupled)
		paralApplyToAll(this->subdomains, &FetiSub<Scalar>::scaleAndSplitKww);

	if(fetiInfo->augment == FetiInfo::WeightedEdges)
		paralApplyToAll(this->subdomains, &FetiSub<Scalar>::weightEdgeGs); // W*Q

	// MPCs
	mpcPrecon = false;
	if(this->glNumMpc > 0) {
		if(fetiInfo->mpc_scaling == FetiInfo::kscaling) { // MPC stiffness scaling
			auto mpcDiagPat = std::make_unique<FSCommPattern<Scalar>>(this->fetiCom, this->cpuToSub, this->myCPU,
			                                                              FSCommPattern<Scalar>::CopyOnSend);
			for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setMpcDiagCommSize(mpcDiagPat.get());
			mpcDiagPat->finalize();
			paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::sendMpcDiag, mpcDiagPat.get());
			mpcDiagPat->exchange();
			paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectMpcDiag, mpcDiagPat.get());
		}

		if(fetiInfo->c_normalize) normalizeC();

		if(fetiInfo->mpc_precno == FetiInfo::diagCCt) {
			// use W scaling for preconditioning mpcs, don't need to build & invert CC^t
			paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::sendMpcScaling, this->vPat);
			this->vPat->exchange();  // neighboring subs mpc weights
			paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectMpcScaling, this->vPat);
		}
		else if(fetiInfo->mpc_precno != FetiInfo::noMpcPrec) {
			// used generalized proconditioner for mpcs, need to build and invert CC^t
			mpcPrecon = true;
			paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::initMpcScaling);
			buildCCt();
		}
	}

	// Factor matrices: K and Kii (if dirichlet preconditioner))
	if(verboseFlag) filePrint(stderr," ... Factor Subdomain Matrices      ... \n");
	startTimerMemory(this->times.factor, this->times.memoryFactor);
	if(fetiInfo->local_cntl->subtype != FetiInfo::spooles && fetiInfo->local_cntl->subtype != FetiInfo::mumps) // spooles/mumps factor is apparently not thread-safe
		timedParal(this->times.factorMat, this->nsub, this, &GenFetiDPSolver<Scalar>::factorLocalMatrices);
	else
		for(iSub=0; iSub<this->nsub; ++iSub) factorLocalMatrices(iSub);

	stopTimerMemory(this->times.factor, this->times.memoryFactor);
	if(fetiInfo->augment == FetiInfo::Gs) {
		this->makeRbmPat();
		paralApplyToAll(this->subdomains, &FetiSub<Scalar>::precondGrbm);
		// Get all the numbers of rigid body modes and dispatch RBMs to neighbors
		paralApplyToAll(this->subdomains, &FetiSub<Scalar>::sendInterfaceGrbm, this->rbmPat);
		this->rbmPat->exchange();
		paralApplyToAll(this->subdomains, &FetiSub<Scalar>::receiveInterfaceGrbm, this->rbmPat);
	}

	// Make coarse problems (Kcc^* and GtG)
	makeKcc();
	if(ngrbms) makeGtG();  // currently G = C^T*R (ie restriction of R to mpc interface)

	// build CC^t for preconditioning mpc residual if necessary
	// done earlier if(mpcPrecon) buildCCt();

	if(domain->solInfo().isCoupled) wetInterfaceComms();

	int tLocalCLen = 0;
	for(iSub = 0; iSub < this->nsub; ++iSub) {
		internalC.domLen[iSub] = this->subdomains[iSub]->numCoarseDofs();
		tLocalCLen += internalC.domLen[iSub];
	}
	internalC.len = tLocalCLen;

	if(domain->solInfo().isCoupled) {
		int tLocalWILen = 0;
		for(iSub = 0; iSub < this->nsub; ++iSub) {
			internalWI.domLen[iSub] = this->subdomains[iSub]->numWetInterfaceDofs();
			tLocalWILen += internalWI.domLen[iSub];
		}
		internalWI.len = tLocalWILen;
		internalWI.setMasterFlag();
	}

	// Allocate Distributed Vectors necessary for FETI solve loop
	this->times.memoryDV -= memoryUsed();
	int numC = (KccSolver) ? KccSolver->neqs() : 0;
	if(KccParallelSolver) numC = coarseInfo->totLen();
	this->wksp = new GenFetiWorkSpace<Scalar>(this->interface, internalR, internalWI, ngrbms, numC, globalFlagCtc);
	this->times.memoryDV += memoryUsed();

	threadManager->callParal(this->nsub, [this](int iSub) { this->subdomains[iSub]->makeBs(); });

	paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::initMpcStatus);

	t6 += getTime();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeLocalWaveNumbers()
{
  int iSub;
  if(fetiInfo->numdir > 0) { // support for multiple fluids 
    if(fetiInfo->waveMethod == FetiInfo::uniform) {
      paralApply(this->subdomains, &FetiBaseSub::computeWaveNumbers);
    }
    else if(fetiInfo->waveMethod == FetiInfo::averageK) {
      paralApply(this->subdomains, &FetiBaseSub::computeWaveNumbers);
      // send/receive wave numbers for FETI-DPH EdgeWs augmentation
      FSCommPattern<double> *kPat = new FSCommPattern<double>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<double>::CopyOnSend);
      for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setCommSize(kPat, 4);
      kPat->finalize();
      paralApply(this->subdomains, &FetiBaseSub::sendWaveNumbers, kPat);
      kPat->exchange();
      paralApply(this->subdomains, &FetiBaseSub::collectWaveNumbers,kPat);
      delete kPat;
    }
    else {
      paralApply(this->subdomains, &FetiBaseSub::averageMatProps);
      // send/receive neighb material properties for FETI-DPH EdgeWs augmentation
      FSCommPattern<double> *matPat = new FSCommPattern<double>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<double>::CopyOnSend);
      for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setCommSize(matPat, 5);
      matPat->finalize();
      paralApply(this->subdomains, &FetiBaseSub::sendMatProps, matPat);
      matPat->exchange();
      paralApply(this->subdomains, &FetiBaseSub::collectMatProps, matPat);
      delete matPat;
    }
  }
}


template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeKcc()
{
	startTimerMemory(this->times.coarse1, this->times.memoryGtG);

	if(verboseFlag) filePrint(stderr," ... Number of Subdomains    %5d  ...\n", this->glNumSub);

	// STEP 1. count number of corner nodes and make subToCorner connectivity
	// A global subdomain to corner connectivity. Does not include primal augmentation.
	std::unique_ptr<Connectivity>subToCorner = 0;
	if(!cornerToSub) {
		subToCorner = std::make_unique<Connectivity>( makeCornerToSub() );
	}
	else {
		subToCorner = std::make_unique<Connectivity>( cornerToSub->reverse() );
	}

	// STEP 2. make subToCoarse connectivity (includes primal mpcs and augmentation) and coarseConnectivity
	Connectivity *coarseToSub = 0;
	Connectivity *subToCoarse = 0;
	Connectivity *coarseConnectivity = 0;
	if(isFeti(fetiInfo->coarse_cntl->type)) {
		if(this->glNumMpc_primal > 0) {
			filePrint(stderr, " *** ERROR: Mpc not supported with multi-level solver *** \n");
			exit(-1);
		}
		augOffset = glNumCorners;
		if((fetiInfo->augment == FetiInfo::WeightedEdges ||
		    fetiInfo->augment == FetiInfo::Edges) && !this->edgeToSub)
			makeEdgeConnectivity();
	} else {
		coarseToSub = cornerToSub;
		mpcOffset = coarseToSub->csize();
		if(this->glNumMpc_primal > 0) {
			coarseToSub = new Connectivity(coarseToSub->append(*this->mpcToSub_primal));
		}
		augOffset = coarseToSub->csize();
		switch(fetiInfo->augment) {
			case FetiInfo::Gs: {
				Connectivity *augcoarseToSub = new Connectivity(coarseToSub->append(*this->subToSub));
				if(coarseToSub != cornerToSub) delete coarseToSub;
				coarseToSub = augcoarseToSub;
			} break;
			case FetiInfo::WeightedEdges:
			case FetiInfo::Edges: {
				if(!this->edgeToSub) makeEdgeConnectivity();
				Connectivity *augcoarseToSub = new Connectivity(coarseToSub->append(*this->edgeToSub));
				if(coarseToSub != cornerToSub) delete coarseToSub;
				coarseToSub = augcoarseToSub;
				// filePrint(stderr,"coarseToSub %d kb\n",(int)(coarseToSub->memsize()/1024));
			} break;
			default:
				break;
		}
		subToCoarse = (coarseToSub != cornerToSub) ? coarseToSub->alloc_reverse() : subToCorner.get();
		// filePrint(stderr,"subToCoarse %d kb\n",(int)(subToCoarse->memsize()/1024));
		coarseConnectivity = coarseToSub->transcon(subToCoarse);
		// filePrint(stderr,"coarseConnectivity %d kb\n",(int)(coarseConnectivity->memsize()/1024));
	}

	// STEP 3. make the coarse problem equation numberer
	int renumFlag = (fetiInfo->coarse_cntl->subtype == FetiInfo::skyline || fetiInfo->coarse_cntl->subtype == FetiInfo::blocksky) ? domain->solInfo().renum : 0;
	if(cornerEqs) delete cornerEqs;
	if(isFeti(fetiInfo->coarse_cntl->type)) {
		cornerEqs = 0;
	} else {
		compStruct renumber = coarseConnectivity->renumByComponent(renumFlag);
		cornerEqs = new DofSetArray(coarseConnectivity->csize(), renumber.renum, 1);
		delete [] renumber.xcomp;
	}

#if defined(USE_MUMPS) && defined(DISTRIBUTED)
	if(fetiInfo->coarse_cntl->subtype == FetiInfo::mumps && fetiInfo->coarse_cntl->mumps_icntl[18] == 3) { // matrix is distributed, use local graph for matrix structure
   delete coarseConnectivity;
   auto subToCPU = this->cpuToSub->reverse();
   coarseConnectivity = coarseToSub->transcon(subToCoarse, subToCPU.getTarget(), this->myCPU);
 }
#endif
	if(coarseToSub && coarseToSub != cornerToSub) delete coarseToSub;
	if(subToCoarse != subToCorner.get()) delete subToCoarse;

	if ( !isFeti(fetiInfo->coarse_cntl->type) ) {
		std::vector<int> glCornerDofs(glNumCorners, 0);
		for(int iSub = 0; iSub < this->nsub; ++iSub)
			this->subdomains[iSub]->markCornerDofs(glCornerDofs);
#ifdef DISTRIBUTED
		this->fetiCom->globalMpiOp(glNumCorners, glCornerDofs.data(), MPI_BOR); // MPI_BOR is an mpi bitwise or
#endif
		for(int i = 0; i < glNumCorners; ++i) cornerEqs->mark(i, glCornerDofs[i]);

		if(this->glNumMpc_primal > 0) {
			for(int iMPC=0; iMPC<this->glNumMpc_primal; ++iMPC)
				cornerEqs->setWeight(mpcOffset+iMPC, 1);
			this->times.numMPCs = this->glNumMpc_primal;
		}
	}

	if(fetiInfo->augment == FetiInfo::Gs) {
		this->times.numRBMs = 0;
		std::vector<int> numRBMPerSub(this->glNumSub);
		for(int iSub=0; iSub<this->glNumSub; ++iSub) numRBMPerSub[iSub] = 0;
		for(int iSub=0; iSub<this->nsub; ++iSub)
			numRBMPerSub[this->subdomains[iSub]->subNum()] = this->subdomains[iSub]->numRBM();
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(numRBMPerSub);
#endif
		if(!isFeti(fetiInfo->coarse_cntl->type)) {
			for(int iSub=0; iSub<this->glNumSub; ++iSub) {
				cornerEqs->setWeight(augOffset+iSub, numRBMPerSub[iSub]);
				this->times.numRBMs += numRBMPerSub[iSub];
			}
		}
	}

	if(fetiInfo->isEdgeAugmentationOn()) {
		this->times.numEdges = 0;
		std::vector<int> edgeWeights(this->edgeToSub->csize());
		for(int i=0; i<this->edgeToSub->csize(); ++i) edgeWeights[i] = 0;
		for(int iSub=0; iSub<this->nsub; ++iSub) {
			int numNeighbor = this->subdomains[iSub]->numNeighbors();
			int myNum = this->subdomains[iSub]->subNum();
			int jEdgeN = 0;
			for(int jSub=0; jSub<numNeighbor; ++jSub) {
				int subJ = this->subdomains[iSub]->getSComm()->subNums[jSub];
				if(this->subdomains[iSub]->isEdgeNeighbor(jSub)) {
					if(myNum < subJ) {
						int nEdge = this->subdomains[iSub]->numEdgeDofs(jSub);
						edgeWeights[(*this->subToEdge)[this->subdomains[iSub]->subNum()][jEdgeN]] = nEdge;
					}
					jEdgeN++;
				}
			}
		}
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(edgeWeights);
#endif
		if(!isFeti(fetiInfo->coarse_cntl->type)) {
			for(int iEdge=0; iEdge<this->edgeToSub->csize(); ++iEdge) {
				cornerEqs->setWeight(augOffset+iEdge, edgeWeights[iEdge]);
				this->times.numEdges += edgeWeights[iEdge];
			}
		}
	}

	if(!isFeti(fetiInfo->coarse_cntl->type)) {
		cornerEqs->finish();
		this->times.numCRNs = cornerEqs->size();
	} else
		this->times.numCRNs = 0;

	if(verboseFlag) {
		filePrint(stderr, " ... Size of Interface  %7d     ...\n", this->interface.len);
		filePrint(stderr, " ... Size of Global Kcc %7d     ...\n", this->times.numCRNs);
		filePrint(stderr, " ... global_rbm_tol     = %9.3e ...\n", fetiInfo->coarse_cntl->trbm);
	}

	// STEP 4. make the rigid body modes and condensed equation numberer
	// PJSA: start new code **************************************************

	int nBodies = (bodyToSub) ? bodyToSub->csize() : 0;
	bool mbgflag = false;
	Connectivity *groupToBody = 0;
	if(this->glNumMpc_primal > 0) {

		// Find out which bodies are connected together by mpcs
		// a collection of inter-connected bodies is referred to as a group
		Connectivity subToMpc = this->mpcToSub_primal->reverse();
		Connectivity bodyToMpc = bodyToSub->transcon(subToMpc);
		Connectivity mpcToBody = bodyToMpc.reverse();
		Connectivity bodyToBodyTmp = bodyToMpc.transcon(mpcToBody);
		Connectivity bodyToBody = bodyToBodyTmp.withSelfConnection();
		compStruct renumber = bodyToBody.renumByComponent(1);  // 1 = sloan renumbering
		nGroups = renumber.numComp;
		if(nGroups < nBodies) { // at least one multi-body group exists
			mbgflag = true;
			// make groupToBody connectivity
			renumber.order = new int[nBodies];
			for(int i = 0; i < nBodies; ++i)
				renumber.order[renumber.renum[i]] = i;
			int *pointer = new int[nGroups + 1];
			pointer[0] = 0;
			for(int i = 0; i < nGroups; ++i) {
				int nbod = renumber.xcomp[i + 1] - renumber.xcomp[i];
				pointer[i + 1] = pointer[i] + nbod;
			}
			groupToBody = new Connectivity(nGroups, pointer, renumber.order);
			groupToSub = groupToBody->transcon(bodyToSub);
			subToGroup = groupToSub->alloc_reverse();
		}
		delete [] renumber.xcomp;
		delete [] renumber.renum;
	}
	if(!mbgflag) { // one body per group
		groupToSub = bodyToSub;
		subToGroup = subToBody;
		nGroups = nBodies;
	}

	// tell each subDomain what group it is in and find groups on each processor
	paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::setGroup, subToGroup);
	if(groups) delete [] groups;
	groups = new int[nGroups];  // groups represented on this processor
#ifdef DISTRIBUTED
	if(this->subdomains.size() > 0) {
		groups[0] = (*subToGroup)[this->subdomains[0]->subNum()][0];
		int n = 1;
		for(int i = 1; i < this->nsub; ++i) {
			int group = (*subToGroup)[this->subdomains[i]->subNum()][0];
			int j;
			for(j = 0; j < n; ++j) if(group == groups[j]) j = n+1;
			if(j == n) groups[n++] = group;
		}
		nGroups1 = n;  // number of groups represented on this processor
	}
	else nGroups1 = 0;
#else
	for(int i = 0; i < nGroups; ++i) groups[i] = i;
 nGroups1 = nGroups;
#endif
	if(ngrbmGr) delete [] ngrbmGr;
	ngrbmGr = new int[nGroups];
	for(int i = 0; i < nGroups; ++i) ngrbmGr[i] = 0;
	if(verboseFlag) filePrint(stderr, " ... Number of bodies = %3d         ...\n", nGroups);

	int ngrbm = 0;
	ngrbms = 0;
	if(rbmFlag && geometricRbms) {
		if((fetiInfo->corners == FetiInfo::noCorners) && (this->glNumMpc_primal > 0)) {
			// subdomain ZEMs need to be projected or eliminated
			// I think adding the dofs of primal mpcs to the "c" dofs would resolve this issue
			filePrint(stderr, " *** ERROR: mpc_type 2 is not supported with no corners *** \n");
			exit(-1);
		}
		if(verboseFlag) filePrint(stderr, " ... Computing multi-body GRBMs     ...\n");
		// calculate the centroid of each group
		std::vector<double> centroid(nGroups*3);   // pseudo centroid of each group
		std::vector<double> nNodes(nGroups);       // number of nodes in each group;
		for(int i = 0; i < nGroups; ++i) {
			nNodes[i] = 0.0;
			for(int j = 0; j < 3; ++j) centroid[3*i+j] = 0.0;
		}

		// groups could be done in parallel
		for(int i = 0; i < this->nsub; ++i) {
			auto csum = this->subdomains[i]->getNodeSet().computeSums();
			for (int j = 0; j < 3; ++j)
				centroid[this->subdomains[i]->group*3+j] = csum.first[j];
			nNodes[this->subdomains[i]->group] += csum.second;
		}

#ifdef DISTRIBUTED
		this->fetiCom->globalSum(nGroups, nNodes.data());
		this->fetiCom->globalSum(nGroups*3, centroid.data());
#endif
		for(int i = 0; i < nGroups; ++i) {
			if(nNodes[i] > 0)
				for(int j = 0; j < 3; ++j)
					centroid[3*i+j] = centroid[3*i+j]/nNodes[i];
		}
		// note: this centroid calculation assumes nodes are equally spaced, and also
		// interface nodes from a group split over more than one interface are used more than once.
		// however, the objective is only to find a point inside the group to be used as
		// a reference for the geometric rbm calculation.  it is not necessary to use the exact
		// geometric centroid.

		// make Zstar and R matrices for each subdomain
		paralApply(this->subdomains, &FetiSub<Scalar>::makeZstarAndR, centroid.data());

		Connectivity groupToMpc;
		if(this->glNumMpc_primal > 0) {
			Connectivity subToMpc = this->mpcToSub_primal->reverse();
			groupToMpc = groupToSub->transcon(subToMpc);
			paralApply(this->subdomains, &FetiBaseSub::makeLocalToGroupMPC, groupToMpc);
		}

		// assemble global Zstar matrix for each body
		std::vector<FullM> globalZstar;
		globalZstar.reserve(nGroups);
		std::vector<int> zRow(nGroups);
		std::vector<int> zRowDim(nGroups);
		std::vector<int> zColDim(nGroups);
		std::vector<int> zColOffset(nBodies);
		int zColDim1 = (this->subdomains.size() > 0) ? this->subdomains[0]->zColDim() : 0;  // (6 for 3D, 3 for 2D)
#ifdef DISTRIBUTED
		zColDim1 = this->fetiCom->globalMax(zColDim1);  // enforce it to be the same
#endif
		for(int i = 0; i < nGroups; ++i) {
			zRowDim[i] = 0;
			if(mbgflag) {
				int cOff = 0;
				for(int j = 0; j < groupToBody->num(i); ++j) {
					int body = (*groupToBody)[i][j];
					zColOffset[body] = cOff;
					cOff += zColDim1;
				}
				zColDim[i] = cOff;
			}
			else {
				zColOffset[i] = 0;
				zColDim[i] = zColDim1;
			}
		}
		for(int iSub = 0; iSub < this->nsub; ++iSub) {
			int subGroup = (*subToGroup)[this->subdomains[iSub]->subNum()][0];
			zRowDim[subGroup] += this->subdomains[iSub]->zRowDim();
		}
		if(this->glNumMpc_primal > 0) {
			for(int i = 0; i < nGroups; ++i) zRowDim[i] += groupToMpc.num(i);
		}
		if(groupToBody) delete groupToBody;

#ifdef DISTRIBUTED
		std::vector<int> zRowOffset(this->numCPUs*nGroups);
		for(int i = 0; i < this->numCPUs*nGroups; ++i) zRowOffset[i] = 0;
		for(int i = 0; i < nGroups1; ++i) {
			int iGroup = groups[i];
			for(int j = this->myCPU+1; j < this->numCPUs; ++j) zRowOffset[iGroup*this->numCPUs +j] = zRowDim[iGroup];
		}
		this->fetiCom->globalSum(nGroups, zRowDim.data());
		this->fetiCom->globalSum(zRowOffset);
		for(int i = 0; i < nGroups; ++i) zRow[i] = zRowOffset[i*this->numCPUs + this->myCPU];
#else
		for(int i = 0; i < nGroups; ++i) zRow[i] = 0;
#endif
		for(int i = 0; i < nGroups; ++i) {
			globalZstar.emplace_back(zRowDim[i], zColDim[i]);
			globalZstar[i].zero();
		}
		// could do this in parallel (by groups)
		for(int iSub = 0; iSub < this->nsub; ++iSub) {
			int subBody = (*subToBody)[this->subdomains[iSub]->subNum()][0];
			int subGroup = (*subToGroup)[this->subdomains[iSub]->subNum()][0];
			if(this->subdomains[iSub]->zRowDim() > 0)
				this->subdomains[iSub]->addSPCsToGlobalZstar(globalZstar[subGroup], zRow[subGroup], zColOffset[subBody]);
			if(this->subdomains[iSub]->numMPCs_primal() > 0) {
				int startRow = zRowDim[subGroup] - groupToMpc.num(subGroup);
				this->subdomains[iSub]->addMPCsToGlobalZstar(globalZstar[subGroup], startRow, zColOffset[subBody], zColDim1);
			}
			if(this->glNumMpc_primal > 0) {
				int subBody = (*subToBody)[this->subdomains[iSub]->subNum()][0];
				this->subdomains[iSub]->setBodyRBMoffset(zColOffset[subBody]);
			}
		}

		std::vector<int> groupProc(nGroups);
#ifdef DISTRIBUTED
		for(int i = 0; i < nGroups; ++i) {
			this->fetiCom->globalSum(zRowDim[i]*zColDim[i], globalZstar[i].data());
			groupProc[i] = -1;
		}
		for(int i = 0; i < nGroups1; ++i) groupProc[groups[i]] = this->myCPU;
		for(int i = 0; i < nGroups; ++i) groupProc[i] = this->fetiCom->globalMax(groupProc[i]);
#else
		for(int i = 0; i < nGroups; ++i) groupProc[i] = this->myCPU;
#endif

		// now do svd on globalZstar for each group to get globalQ for each group
		ngrbm = 0;  // total of all groups
		FullM  **Qtranspose;
		Qtranspose = new FullM * [nGroups];
		for(int i = 0; i < nGroups1; ++i) {
			int iGroup = groups[i];
			int ncol = zColDim[iGroup];
			int nrow = zRowDim[iGroup];
			FullM U(ncol,ncol); U.zero();
			int rank = 0;
			singularValueDecomposition(globalZstar[iGroup], U, ncol, nrow, rank, domain->solInfo().tolsvd);
			int ngrbmGrTmp = ncol - rank;
			globalZstar[iGroup].clean_up();
			if(groupProc[iGroup] == this->myCPU) {
				ngrbmGr[iGroup] = ngrbmGrTmp;
				ngrbm += ngrbmGr[iGroup];
				fprintf(stderr, " ... Number of GRBMs for body %d: %d ...\n", iGroup, ngrbmGrTmp);
			}
			Qtranspose[iGroup] = new FullM(U, ngrbmGrTmp, rank, ncol, 0);
		}
#ifdef DISTRIBUTED
		ngrbms = this->fetiCom->globalSum(ngrbm);  // total number of rigid body modes for all processes
#else
		ngrbms = ngrbm;
#endif
		if(verboseFlag) filePrint(stderr, " ... total number of GRBMs = %5d  ...\n", ngrbms);

		// make local Rstar
		paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::makeLocalRstar, Qtranspose);
		for(int i = 0; i < nGroups1; ++i) delete Qtranspose[groups[i]];
		delete [] Qtranspose;
	}
	if(rbmFlag && !geometricRbms && fetiInfo->corners == FetiInfo::noCorners) {
		paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::useKrrNullspace);
		ngrbmGr = new int[nGroups];
		ngrbm = 0;
		for(int i = 0; i < nGroups; ++i) ngrbmGr[i] = 0;
		for(int iSub = 0; iSub < this->nsub; ++iSub) {
			int subGroup = (*subToGroup)[this->subdomains[iSub]->subNum()][0];
			ngrbmGr[subGroup] = this->subdomains[iSub]->Krr->numRBM();
			ngrbm += ngrbmGr[subGroup];
		}
#ifdef DISTRIBUTED
		ngrbms = this->fetiCom->globalSum(ngrbm);  // total number of rigid body modes for all processes
#else
		ngrbms = ngrbm;
#endif
	}

	if(ngrbms) {
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(nGroups, ngrbmGr);
#endif
		paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::setNumGroupRBM, ngrbmGr);
		paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::deleteLocalRBMs);
	}

	// end new code ********************************************************

	KccSparse = 0;

//  if(cornerEqs->size() > 0) {
	if(1) {

		// assemble the coarse problem: Kcc^* -> Kcc - Krc^T Krr^-1 Krc
		if(verboseFlag) filePrint(stderr, " ... Assemble Kcc solver            ...\n");
		t5 -= getTime();
		paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::formKccStar); // create the local Kcc^*
		t5 += getTime();

		// EXPERIMENTAL CODE to use feti dp for the coarse problem
		if(fetiInfo->coarse_cntl->type == SolverSelection::Feti) {
			// build coarse solver
			makeMultiLevelDP(std::move(subToCorner));
		}
		else {

			this->times.memoryGtGsky -= memoryUsed();
			int sparse_ngrbms = (geometricRbms) ? ngrbms : 0; // TODO pass Rbm object, not just ngrbms
#ifdef DISTRIBUTED
			if(this->glNumSub == this->numCPUs && fetiInfo->type != FetiInfo::nonlinear &&
			   !domain->solInfo().doEigSweep && !domain->solInfo().doFreqSweep &&
			   (fetiInfo->coarse_cntl->type == SolverSelection::Direct && (fetiInfo->coarse_cntl->subtype == 0 || fetiInfo->coarse_cntl->subtype == 1)))
				KccSolver = GenSolverFactory<Scalar>::getFactory()->createDistSolver(coarseConnectivity, cornerEqs, *fetiInfo->coarse_cntl,
				                                                                     KccSparse, this->fetiCom);
			else
#endif
				KccSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(coarseConnectivity, cornerEqs, *fetiInfo->coarse_cntl,
				                                                                 KccSparse, sparse_ngrbms, this->fetiCom);

			this->times.memoryGtGsky += memoryUsed();


			t0 -= getTime();
			paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::makeKccDofs, cornerEqs, augOffset, this->subToEdge, mpcOffset);
			if(KccSparse)
				for (const auto &subdomain : this->subdomains)
					KccSparse->add(*subdomain->Kcc, subdomain->getCornerEqNums().data());

			t0 += getTime();

			// Factor coarse solver
			startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
#ifdef DISTRIBUTED
			if(verboseFlag) filePrint(stderr, " ... Unify Kcc                      ...\n");
			KccSolver->unify(this->fetiCom);
#endif
			if(fetiInfo->printMatLab) {
				KccSparse->printSparse("coarsemat");
			}

			if(verboseFlag) filePrint(stderr, " ... Factor Kcc solver              ...\n");
			KccSolver->setPrintNullity(fetiInfo->contactPrintFlag && this->myCPU == 0);
			KccSolver->parallelFactor();
			stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);

			if(rbmFlag && geometricRbms && this->myCPU == 0) {
				if(KccSolver->neqs() > 0 && KccSolver->numRBM() != ngrbms) {
					std::cerr << " *** WARNING: number of singularities in Kcc (" << KccSolver->numRBM() << ")" << std::endl
					          << "     does not match the number of Geometric RBMs (" << ngrbms << ")" << std::endl
					          << " *** try adjusting global_cor_rbm_tol or use TRBM method" << std::endl;
				}
			}

			if(rbmFlag && !geometricRbms && (ngrbms = KccSolver->numRBM()) > 0) {
				kccrbms = new Scalar[KccSolver->neqs()*KccSolver->numRBM()];
				KccSolver->getNullSpace(kccrbms);
				if(fetiInfo->nullSpaceFilterTol > 0.0) {
					for(int i= 0; i<KccSolver->numRBM(); ++i)
						for(int j=0; j<KccSolver->neqs(); ++j)
							if(ScalarTypes::norm(kccrbms[i*KccSolver->neqs()+j]) < fetiInfo->nullSpaceFilterTol) kccrbms[i*KccSolver->neqs()+j] = 0.0; // FILTER
				}
			}
		}

	} else
		KccSolver = 0;

	if(coarseConnectivity) delete coarseConnectivity;
	stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

/** \details The corners are given global unique numbers from 0 to glNumCorners-1.
 *
 * Depends only on non 'Scalar' related members:
 * glNumSub, subdomains, fetiCom,
 *
 * Also builds glNumCorners.
 * @tparam Scalar
 * @return
 */
template <typename Scalar>
Connectivity GenFetiDPSolver<Scalar>::makeCornerToSub() {
	std::vector<size_t> pointer(this->glNumSub + 1, 0);
	for(int iSub=0; iSub < this->nsub; ++iSub)
			pointer[this->subdomains[iSub]->subNum()] = this->subdomains[iSub]->numCorners();
#ifdef DISTRIBUTED
	this->fetiCom->globalSum(this->glNumSub, pointer.data());
#endif
	size_t total = 0;
	for(int iSub=0; iSub < this->glNumSub; ++iSub) {
			int tmp = pointer[iSub];
			pointer[iSub] = total;
			total += tmp;
		}
	pointer[this->glNumSub] = total;

	int *glCornerNodes = new int[total];
	for(int i=0; i<total; ++i) glCornerNodes[i] = 0;
	for(int iSub=0; iSub < this->nsub; ++iSub) {
			int numCorner = this->subdomains[iSub]->numCorners();
			const auto &localCornerNodes = this->subdomains[iSub]->getLocalCornerNodes();
			auto glN = this->subdomains[iSub]->getGlNodes();
			for(int iCorner=0; iCorner<numCorner; ++iCorner)
				glCornerNodes[pointer[this->subdomains[iSub]->subNum()] + iCorner] = glN[localCornerNodes[iCorner]];
		}
#ifdef DISTRIBUTED
	this->fetiCom->globalSum(total, glCornerNodes);
#endif
	glNumCorners = 0;
#ifdef TFLOP
	std::map<int, int, std::less<int> > glCornerMap;
#else
	map<int, int> glCornerMap;
#endif
	for(int iCorner=0; iCorner<total; ++iCorner)
			if(glCornerMap.find(glCornerNodes[iCorner]) == glCornerMap.end() )
				glCornerMap[ glCornerNodes[iCorner] ] = glNumCorners++;
	delete [] glCornerNodes;
	if(verboseFlag) filePrint(stderr, " ... Total Number of Corners %5d  ...\n", glNumCorners);

	std::vector<int> target(total, 0);
	for(int iSub=0; iSub < this->nsub; ++iSub) {
			int numCorner    = this->subdomains[iSub]->numCorners();
			auto &cornerNodes = this->subdomains[iSub]->getCornerNodes();
			const auto &localCornerNodes = this->subdomains[iSub]->getLocalCornerNodes();
			auto glN = this->subdomains[iSub]->getGlNodes();
			for(int iCorner=0; iCorner<numCorner; ++iCorner) {
				cornerNodes[iCorner] = glCornerMap[glN[localCornerNodes[iCorner]]];
				target[iCorner+pointer[this->subdomains[iSub]->subNum()]] = cornerNodes[iCorner];
			}
		}
#ifdef DISTRIBUTED
	this->fetiCom->globalSum(total, target.data());
#endif

	Connectivity subToCorner(this->glNumSub, std::move(pointer), std::move(target));
	this->cornerToSub = subToCorner.alloc_reverse();
	return subToCorner;
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::KrrReSolve(int iSub, GenDistrVector<Scalar> &ur)
{
 this->subdomains[iSub]->Krr->reSolve(ur.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::extractFc(int iSub, const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fc) const
{
 this->subdomains[iSub]->getFc(f.subData(iSub), fc.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::getFc(const GenDistrVector<Scalar> &f, GenVector<Scalar> &fc) const
{
 GenDistrVector<Scalar> distFc(internalC);
 distFc.zero(); 

 execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::extractFc, f, distFc);

 // Assemble fc (force at the corners)
 fc.zero();
 for(int iSub=0; iSub<this->nsub; ++iSub) {
    const auto &subdCDofs = this->subdomains[iSub]->getCornerEqNums();
    Scalar *subFc = distFc.subData(iSub);
    for(int iDof = 0; iDof < this->subdomains[iSub]->numCoarseDofs(); ++iDof)
      if(subdCDofs[iDof] > -1)
        fc[subdCDofs[iDof] ] += subFc[iDof];
 }
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::updateActiveSet(GenDistrVector<Scalar> &v, int flag, double tol)
{
  // flag = 0 dual planing
  // flag = 1 primal planing
  paralApply(this->subdomains, &FetiSub<Scalar>::saveMpcStatus2);

  bool *local_status_change = new bool[this->nsub];
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::subUpdateActiveSet, v, tol, flag, local_status_change);
  bool status_change1 = false;
  for(int i=0; i<this->nsub; ++i) if(local_status_change[i]) { status_change1 = true; break; }
#ifdef DISTRIBUTED
  status_change1 = this->fetiCom->globalMax((int) status_change1);
#endif

  if(status_change1) {
    paralApply(this->subdomains, &FetiSub<Scalar>::sendMpcStatus, mpcPat, flag);
    mpcPat->exchange();
    execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::subRecvMpcStatus, mpcPat, flag, local_status_change);
    bool status_change2 = false;
    for(int i=0; i<this->nsub; ++i) if(local_status_change[i]) { status_change2 = true; break; }
#ifdef DISTRIBUTED
    status_change2 = this->fetiCom->globalMax((int) status_change2);
#endif
    if(status_change2 && ngrbms) rebuildGtGtilda();
    if(fetiInfo->contactPrintFlag && this->myCPU == 0) std::cerr << " ";
    delete [] local_status_change;
    return status_change2;
  }
  else {
    delete [] local_status_change;
    return false;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subUpdateActiveSet(int iSub, GenDistrVector<Scalar> &lambda, double tol, int flag, bool *statusChange)
{
  this->subdomains[iSub]->updateActiveSet(lambda.subData(this->subdomains[iSub]->localSubNum()), tol, flag, statusChange[iSub]);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subRecvMpcStatus(int iSub, FSCommPattern<int> *mpcPat, int flag, bool *statusChange)
{
  this->subdomains[iSub]->recvMpcStatus(mpcPat, flag, statusChange[iSub]);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::update(Scalar nu, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p, 
                                GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fp, 
                                GenDistrVector<Scalar> &ur, GenDistrVector<Scalar> &dur, 
                                GenVector<Scalar> &uc, GenVector<Scalar> &duc) const
{
  if(globalFlagCtc) saveStep();
  GenDistrVector<Scalar> &lambda_k = this->wksp->ret_lambda_copy();
  GenDistrVector<Scalar> &r_k = this->wksp->ret_r_copy();

  int i;
  for(i = 0; i < fetiInfo->linesearch_maxit+1; ++i, ++nLinesearchIter) {

    // update lambda
    lambda.linAdd(nu,p);
    if(globalFlagCtc) project(lambda, lambda, true);

    // Update residual (r)
    if(dualStatusChange) {
      p.linC(1.0/nu, lambda, -1.0/nu, lambda_k); // reduced search direction p = (lambda-lambda_copy)/nu
      localSolveAndJump(p, dur, duc, Fp); // recompute Fp using reduced search direction
    }
    r.linAdd(nu, Fp); 

    if(i == 0 && !dualStatusChange) break; // CG step
    else { // gradient projection step
      Scalar rp = r_k*p;
      Scalar pFp = p*Fp; 
      //if(ScalarTypes::lessThan(pFp, 0)) throw std::runtime_error("FETI operator is not positive semidefinite");
      Scalar delta_f = nu*nu/2.0*pFp + nu*rp;
      if(fetiInfo->contactPrintFlag >= 2 && this->myCPU == 0)
        std::cerr << " linesearch: iteration = " << i << ", delta_f = " << delta_f << ", pFp = " << pFp << ", nu = " << nu << ", rp = " << rp << std::endl;
      if(ScalarTypes::lessThanEq(delta_f, 0)) break; // sequence is monotonic (note: check for gcr and gmres)
      else {
        if(i < fetiInfo->linesearch_maxit) { 
          restoreStep();
          nu *= fetiInfo->linesearch_tau;
        }
        else {
          throw std::runtime_error("linesearch did not converge");
        }
      }
    }
  }

  stepLengthChange = (i > 0);
  if(stepLengthChange) nLinesearch++;
  
  // Update primal (ur, uc)
  ur.linAdd(nu,dur); 
  uc += nu*duc;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::saveStep() const
{
  this->wksp->save();
  if(globalFlagCtc) paralApply(this->subdomains, &FetiSub<Scalar>::saveMpcStatus);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::restoreStep() const
{
  this->wksp->restore();
  paralApply(this->subdomains, &FetiSub<Scalar>::restoreMpcStatus);
  if(ngrbms)
	  const_cast<GenFetiDPSolver<Scalar> *>(this)->rebuildGtGtilda();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::solve(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
	if (verboseFlag) filePrint(stderr," ... Begin FETI-DP Solve            ...\n");

	switch(fetiInfo->outerloop) {
		default:
		case FetiInfo::OuterloopType::CG:
			if(verboseFlag) {
				filePrint (stderr, " ... CG solver selected             ...\n\n");
			}
			if(this->myCPU == 0 && typeid(Scalar)==typeid(DComplex) && !fetiInfo->complex_hermitian)
				std::cerr << " *** WARNING: CG is not valid for all complex symmetric matrices even if positive definite.\n"
				          << " *** If your matrix is complex hermitian include \"interf_solver CG hermitian\" under FETI in your input file.\n"
				          << " *** Otherwise, \"interf_solver GMRES\" or \"interf_solver GCR\" are safer choices.\n";
			if(verboseFlag) {
				filePrint (stderr, "----------------------------------------------\n");
				filePrint (stderr, "          Iterations loop monitoring          \n");
				filePrint (stderr, "----------------------------------------------\n");
			}
			solveCG(f, u);
			break;
		case FetiInfo::OuterloopType::GMRES:
			if(verboseFlag) {
				filePrint (stderr, " ... GMRES solver selected          ...\n");
				filePrint (stderr, "----------------------------------------------\n");
				filePrint (stderr, "          Iterations loop monitoring          \n");
				filePrint (stderr, "----------------------------------------------\n");
			}
			solveGMRES(f, u);
			break;
		case FetiInfo::OuterloopType::GCR:
			if(verboseFlag) {
				filePrint (stderr, " ... GCR solver selected            ...\n\n");
				filePrint (stderr, "----------------------------------------------\n");
				filePrint (stderr, "          Iterations loop monitoring          \n");
				filePrint (stderr, "----------------------------------------------\n");
			}
			solveGCR(f, u);
			break;
	}
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solveCG(const GenDistrVector<Scalar> &_f, GenDistrVector<Scalar> &u)
{
 t7 = -getTime(); this->times.solve -= getTime();
    GenDistrVector<Scalar> f{_f};
 // vectors
 GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();  
 GenDistrVector<Scalar> &ur      = this->wksp->ret_ur(); 
 GenDistrVector<Scalar> &dur     = this->wksp->ret_du(); 
 GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda();      // lagrange multipliers
 GenDistrVector<Scalar> &r       = this->wksp->ret_r();           // residual
 GenDistrVector<Scalar> &w       = this->wksp->ret_w();           // projected residual
 GenDistrVector<Scalar> &y       = this->wksp->ret_y();           // re-projected residual
 GenDistrVector<Scalar> &z       = this->wksp->ret_z();           // preconditioned residual
 GenDistrVector<Scalar> &p       = this->wksp->ret_p();           // search direction
 GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp(); 
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); 
 GenDistrVector<Scalar> &fw      = this->wksp->ret_fw(); 
 GenVector<Scalar> &fc           = this->wksp->ret_fc(); 
 GenVector<Scalar> &uc           = this->wksp->ret_uc();  
 GenVector<Scalar> &duc          = this->wksp->ret_duc(); 

 double ww, ww0, ff, error;
 int iter = 0;
 Scalar pFp, nu;

 nMatVecProd = nRebuildGtG = nRebuildCCt = nLinesearch = nLinesearchIter = nSubIterDual = nSubIterPrimal = nStatChDual
	 = nStatChPrimal = 0;
 if(globalFlagCtc && this->numSystems > 0)
	 const_cast<GenFetiDPSolver<Scalar> *>(this)->resetOrthoSet();

 // extract fr, fc and fw from f
 ff = extractForceVectors(f, fr, fc, fw);
 
 // Compute initial lagrange multipliers: lambda^0 = P * lambda^00 - G*(G^T*G)^-1 * e, where e = R^T*f ... also gamma = G^T*lambda+e
 // also for feti-dpc: expand active set and re-compute initial lagrange multipliers if lambda^0 is not feasible
 computeL0(lambda, f); 

 // Compute initial residual: r^0 = F*lambda^0 + d ... also uc = Kcc^-1(fc+Krr^-1*Krc^T*lambda), ur = Krr^-1*(fr-Br^T*lambda-Krc*uc)
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw); 

 // Compute initial projected residual: w^0 = P^T * r^0 
 // also for feti-dpc: contract active set if w^0 is not proportional (primal planing)
 ww0 = ww = tProject(r, w); 
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(ww0));

 //if(ww == 0.0) { 
 //  error = 0.0;
 //  this->times.iterations[this->numSystems].stagnated = 0;
 //} 
 //else { 
   // Multiple rhs prediction (note: not used for contact)
   if(this->predict(w, lambda)) {
     localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);
     ww = tProject(r, w); 
     if(verboseFlag) filePrint(stderr," ... Initial residual norm after MRHS prediction %e\n", sqrt(ww));
   }

   if(verboseFlag) filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

   for(iter = 0; iter < this->maxiter; ++iter, ++iterTotal) {
     dbg_alloca(0);

     // Precondition: z = M^-1 * w
     error = preCondition(w, z);

     // Check stopping criteria
     bool stop = checkStoppingCriteria(iter, error, ff);

     // Print errors
     if(verboseFlag && (stop || ((fetiInfo->numPrint() > 0) && (iter % fetiInfo->numPrint() == 0))))
       filePrint(stderr," %4d %23.6e %21.6e\n", iter, (ff == 0.0) ? sqrt(error) : sqrt(error/ff), sqrt(ww/ww0));

     // Exit CG iterations loop if one of the stopping criteria has been satisfied
     if(stop) break;

     // Krylov acceleration
     if(fetiInfo->nlPrecFlg) this->nlPreCondition(w, z);

     // Re-project: y = P * z
     project(z, y); 

     // Search direction
     this->orthogonalize(y, p);
   
     // Matrix vector product
     localSolveAndJump(p, dur, duc, Fp);

     // Step length
     if(!fetiInfo->complex_hermitian) // Note: CG by default uses the indefinite vector inner product i.e. x^y = x^T y
      { pFp = p^Fp; nu = -(w^p)/pFp; }      // CG in this form can be used for complex symmetric matrices but is NOT valid for
     else                                   // all such matrices, so we don't recommed using it. Try GCR or GMRES instead!!!
      { pFp = p*Fp; nu = -((w*p) / pFp); }  // CG is valid for positive definite Hermitian matrices but to set this up you need 
                                          // to use to the standard inner product x*y = y^H x for pFp, nu and beta (in Feti.d/CGOrthoset.C)
                                          // this option hasn't really been tested since our matrices are typically NOT Hermitian

     // Update solution: lambda += nu*p, r += nu*Fp, ur += nu*dur, uc += nu*duc 
     // optional linesearch to adjust step length nu if sequence is not monotonic
     // also for feti-dpc: reduce step and expand active set if lambda is not feasible (dual planing)
     update(nu, lambda, p, r, Fp, ur, dur, uc, duc);

     // Project: w = P^T * r 
     // also for feti-dpc: contract active set if w is not proportional (primal planing)
     ww = tProject(r, w); 

     // add search direction to orthoset or reset if necessary
     if(globalFlagCtc && (dualStatusChange || primalStatusChange || stepLengthChange))
	     const_cast<GenFetiDPSolver<Scalar> *>(this)->resetOrthoSet();
     else this->orthoAddCG(p, Fp, pFp);
   }

   ur += deltaU; // make solution compatible ur += deltaU
 //}

 // Assemble and store primal solution u
 mergeSolution(ur, uc, u, lambda);

 // Store number of iterations, errors, timings and memory used
 this->setAndStoreInfo(iter+1, (ff == 0.0) ? error : (error/ff), (ww0 == 0.0) ? ww : (ww/ww0));
 if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
 this->times.solve += getTime(); t7 += getTime();
 this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;
 printSummary(iter);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::solveGMRES(const GenDistrVector<Scalar> &_f, GenDistrVector<Scalar> &u)
{
	GenDistrVector<Scalar> f{_f}; // The force gets redistributed between subdomains.
	this->times.solve -= getTime();
	t7 = -getTime();
	if(this->oSetGMRES->numDir() > 0) this->oSetGMRES->reInit(); // just in case someone tries to use MRHS with GMRES

	GenDistrVector<Scalar> &ur = this->wksp->ret_ur(); ur.zero();
	GenDistrVector<Scalar> &dur = this->wksp->ret_du(); dur.zero();
	GenDistrVector<Scalar> &lambda0 = this->wksp->ret_lambda();
	GenDistrVector<Scalar> &r = this->wksp->ret_r(); r.zero();
	GenDistrVector<Scalar> &z = this->wksp->ret_z(); z.zero();
	GenDistrVector<Scalar> &Fp = this->wksp->ret_Fp(); Fp.zero();
	GenVector<Scalar> &uc = this->wksp->ret_uc();
	GenVector<Scalar> &duc = this->wksp->ret_duc();
	GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); deltaU.zero();

	GenDistrVector<Scalar> &fr = this->wksp->ret_fr();
	GenVector<Scalar> &fc = this->wksp->ret_fc();
	GenDistrVector<Scalar> &fw = this->wksp->ret_fw();

	double ff = extractForceVectors(f, fr, fc, fw);

	computeL0(lambda0, f);

	// Solve uc = Kcc^-1(  fc + Krr^-1 Krc^T lambda0)
	//          = Kcc^-1(  fc )
	// Solve ur = Krr^-1 ( fr - Br^T lambda0 - Krc uc)
	//   and  r = Br ur
	localSolveAndJump(fr, lambda0, ur, fc, uc, r, fw);

	double r0Norm2 = r.sqNorm();
	if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(r0Norm2) );
	if(r0Norm2 == 0.0) {  mergeSolution(ur, uc, u, lambda0); return; }

	// Precondition r,
	double error = preCondition(r, z);  // z = F^-1 r

// if((fetiInfo->numPrint() > 0) && (fetiInfo->numPrint() < 10))
	if(verboseFlag)
		filePrint(stderr," Iteration  Relative Primal Error  Relative Preconditioned Dual Error\n");

	GenDistrVector<Scalar> rzero(this->interface);
	GenDistrVector<Scalar> zzero(this->interface);
	GenDistrVector<Scalar> medvec(this->interface);
	GenDistrVector<Scalar> lambda(this->interface);

	rzero = r;
	zzero = z;

	this->initGMRES(z);

	bool primalresidual = fetiInfo->gmresResidual;
	int J = 0;
	const int ReStep = fetiInfo->maxortho; // for restarted GMRES

	for(int iter = 0; true; ++iter) {
		// Arnoldi iteration (Algorithm see Saad SISC)
		for (int j=0; j<ReStep; j++, J++, ++iterTotal) {

			localSolveAndJump(z, dur, duc, Fp); // Fp = F*z

			error = preCondition(Fp, medvec);   // medvec = M^-1*Fp

			// Do Arnoldi step
			double resGMRES = this->orthoAddGMRES(z, medvec);

			if((fabs(resGMRES)<=sqrt(this->epsilon2*ff)) || (J == this->maxiter-1) || primalresidual) {

				primalresidual = true; // Since now we compute the primal residual in each step

				this->GMRESSolution(lambda);

				localSolveAndJump(lambda, dur, duc, Fp);  // Fp = F*lambda

				r.linC(1.0,rzero,1.0,Fp); // r = r0 + Fp

				error = preCondition(r, medvec); // medvec = M^-1*r

				bool isConverged = ((sqrt(error) < sqrt(this->epsilon2*ff)) || (J == this->maxiter-1));

				if(verboseFlag && (isConverged || ((fetiInfo->numPrint() > 0) && (J % fetiInfo->numPrint() == 0))))
					filePrint(stderr, "%4d %23.6e %23.6e\n", iter*ReStep+j,sqrt(error/ff),fabs(resGMRES)/sqrt(ff));

				// Determine convergence
				if (isConverged) {
					// For timing and records
					this->times.iterations[this->numSystems].stagnated = 0;
					// Store number of iterations, primal error and dual error
					this->setAndStoreInfo(iter*ReStep+j, error/ff, resGMRES*resGMRES/ff);
					if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
					this->times.solve += getTime();
					this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;

					t7 += getTime();
					printSummary(iter*ReStep+j+1);

					lambda += lambda0;
					ur += dur;
					uc += duc;
					ur += deltaU; // make solution compatible u += deltaU
					mergeSolution(ur, uc, u, lambda);

					return;
				}
			}
			else
			if(verboseFlag && ((fetiInfo->numPrint() > 0) && (j % fetiInfo->numPrint() == 0)))
				filePrint(stderr, "%4d %47.6e\n", iter*ReStep+j, fabs(resGMRES)/sqrt(ff));
		}

		// restart GMRES
		if(verboseFlag) filePrint(stderr, " *** Krylov Space Full - Restarting GMRES \n");
		if(!primalresidual) {
			this->GMRESSolution(lambda);  // compute incremental solution lambda
			localSolveAndJump(lambda, dur, duc, Fp); // Fp = F*lambda
			r.linC(1.0,rzero,1.0,Fp); // r = rzero + Fp;
		}
		rzero       = r;
		uc += duc;
		ur += dur;

		error = preCondition(Fp, medvec);  // medvec = M^-1*Fp
		z.linC(1.0,zzero,1.0,medvec); // z = zzero + medvec;
		zzero = z;

		lambda0 += lambda;
		primalresidual = fetiInfo->gmresResidual;  // primalresidual might not be reached after restart

		this->oSetGMRES->reInit(); // Reinitialize Krylov space and set z of last step as initial vector
		this->initGMRES(z);
	}
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::solveGCR(const GenDistrVector<Scalar> &_f, GenDistrVector<Scalar> &u)
{
 t7 = -getTime();
 this->times.solve -= getTime();
	GenDistrVector<Scalar> f(_f);
 GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();
 GenDistrVector<Scalar> &ur      = this->wksp->ret_ur();
 GenDistrVector<Scalar> &work1   = this->wksp->ret_du();
 GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda();  // Lagrange multipliers
 GenDistrVector<Scalar> &r       = this->wksp->ret_r();       // residual
 GenDistrVector<Scalar> &z       = this->wksp->ret_z();       // preconditioned residual
 GenDistrVector<Scalar> &p       = this->wksp->ret_p();       // search direction
 GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp();
 GenDistrVector<Scalar> &Fz      = this->wksp->ret_y();
 GenDistrVector<Scalar> &fw      = this->wksp->ret_fw();
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU();
 GenVector<Scalar> &fc           = this->wksp->ret_fc();
 GenVector<Scalar> &uc           = this->wksp->ret_uc();
 GenVector<Scalar> &work2        = this->wksp->ret_duc();

 double rr, rr0, ff, error;
 int iter = 0;

 // extract fr, fc and fw from f
 ff = extractForceVectors(f, fr, fc, fw);

 // compute the initial guess for for the Lagrange Multipliers (zero if no projection)
 computeL0(lambda, f);

 // compute the initial residual
 // Solve uc = Kcc^-1(  fc + Krr^-1 Krc^T lambda)
 //       ur = Krr^-1 ( fr - Br^T lambda - Krc uc)
 //       r = Br ur
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);
 rr0 = rr = r.sqNorm();
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(rr0));

 // multiple rhs prediction
 if(this->predictGCR(r, lambda)) {
   localSolveAndJump(fr, lambda, ur, fc, uc, r, fw); 
   rr = r.sqNorm();
   if(verboseFlag) filePrint(stderr," ... Initial residual norm after MRHS prediction %e\n", sqrt(rr));
 }

 if(verboseFlag) filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

 for(iter = 0; iter < this->maxiter; ++iter, ++iterTotal) {
   dbg_alloca(0);

   // Precondition: z = F^-1 r
   error = preCondition(r, z);

   // Check stopping criteria
   bool stop = checkStoppingCriteria(iter, error, ff);

   // Print errors
   if(verboseFlag && (stop || ((fetiInfo->numPrint() > 0) && (iter % fetiInfo->numPrint() == 0))))
     filePrint(stderr," %4d %23.6e %21.6e\n", iter, (ff == 0.0) ? sqrt(error) : sqrt(error/ff), sqrt(rr/rr0));

   // Exit CG iterations loop if one of the stopping criteria has been satisfied
   if(stop) break;

   localSolveAndJump(z, work1, work2, Fz);

   this->orthogonalizeGCR(z, Fz, p, Fp);  // computes new p, Fp

   Scalar FpFp = Fp * Fp;

   Scalar rFp = r * Fp;
   Scalar nu = -(rFp/FpFp);

   // Update solution: lambda += nu*p, r += nu*Fp
   lambda.linAdd(nu, p);
   r.linAdd(nu, Fp);
   rr = r.sqNorm();

   this->orthoAddGCR(p, Fp, FpFp);
 }

 // get primal solution 
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);

 ur += deltaU;

 // merge ur and uc into u
 mergeSolution(ur, uc, u, lambda);

 // Store number of iterations, errors, timings and memory used
 this->setAndStoreInfo(iter+1, error/ff, sqrt(rr/rr0));
 if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
 this->times.solve += getTime(); t7 += getTime();
 this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;
 printSummary(iter);
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::extractForceVectors(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr,
                                             GenVector<Scalar> &fc, GenDistrVector<Scalar> &fw) const
{
	// distribute force  for sfem inpc
	std::unique_ptr<GenDistrVector<Scalar>> f_copy;
	if(domain->solInfo().inpc) {
		f_copy = std::make_unique<GenDistrVector<Scalar>>(f);
		execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::fSplit, f);
	}

	// extract fr and fw from f
	auto perSubExtract = [&, isCoupled = domain->solInfo().isCoupled](int iSub) {
		fr.subVec(iSub).setZero();
		this->subdomains[iSub]->getFr(f.subData(iSub), fr.subData(iSub));
		if(isCoupled) {
			fw.subVec(iSub).setZero();
			this->subdomains[iSub]->getFw(f.subData(iSub), fw.subData(iSub));
		}
	};
	threadManager->callParal(this->nsub, perSubExtract);

	// Assemble and split fr and fw on subdomain interface (note: f for first system is already split by topological scaling)
	if((this->numSystems == 0 && fetiInfo->scaling == FetiInfo::kscaling) || (this->numSystems > 0 && fetiInfo->rescalef)) {
		if(domain->solInfo().isCoupled) this->distributeForce(fr, fw);
		else this->distributeForce(fr);
	}
	double ffr = fr.sqNorm();
	double ffw = (domain->solInfo().isCoupled) ? fw.sqNorm() : 0.0;

	// extract fc from f
	getFc(f, fc);
	double ffc;
	if(KccParallelSolver) {
		GenStackDistVector<Scalar> DistFc(*coarseInfo, fc.data());
		ffc = DistFc.sqNorm();
	}
	else {
#ifdef DISTRIBUTED
		GenVector<Scalar> fc_copy(fc);
    this->fetiCom->globalSum(fc_copy.size(), fc_copy.data());
    ffc = fc_copy.sqNorm();
#else
		ffc = fc.sqNorm();
#endif
	}

	double ff = f.sqNorm();

	// print norms
	if((fetiInfo->numPrint() > 0) && (fetiInfo->numPrint() < 10) && verboseFlag) {
		filePrint(stderr, " ... f*f   %e             ...\n", ff);
		filePrint(stderr, " ... fr*fr %e             ...\n", ffr);
		filePrint(stderr, " ... fc*fc %e             ...\n", ffc);
		if(domain->solInfo().isCoupled)
			filePrint(stderr, " ... fw*fw %e             ...\n", ffw);
	}

	if(domain->solInfo().inpc) f = (*f_copy);

	return ff == 0.0 ? 1.0 : ff;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::printSummary(int iter) const
{
  if(verboseFlag) {
    if(globalFlagCtc) {
      filePrint(stderr," -----------------------------------------------\n");
      filePrint(stderr," number of main iter.                = %5d ...\n", iter);
      filePrint(stderr," number of line searches             = %5d ...\n", nLinesearch);
      filePrint(stderr," number of line search sub-iter.     = %5d ...\n", nLinesearchIter);
      filePrint(stderr," number of dual planings             = %5d ...\n", nStatChDual); 
      filePrint(stderr," number of dual planing sub-iter.    = %5d ...\n", nSubIterDual);
      filePrint(stderr," number of primal planings           = %5d ...\n", nStatChPrimal);
      filePrint(stderr," number of primal planing sub-iter.  = %5d ...\n", nSubIterPrimal);
      filePrint(stderr," total number of matrix-vector op.s  = %5d ...\n", nMatVecProd); 
      filePrint(stderr," number of GtG rebuilds              = %5d ...\n", nRebuildGtG);
      filePrint(stderr," -----------------------------------------------\n");
    }
    else {
      filePrint(stderr," --------------------------------------\n");
      filePrint(stderr," number of main iter.       = %5d ...\n", iter);
      filePrint(stderr," --------------------------------------\n");
    }
    filePrint(stderr," Assembly of Kcc*   %13.4f s ...\n",t0/1000.0);
    filePrint(stderr," Krc^T Krr^-1 f     %13.4f s ...\n",t1/1000.0);
    filePrint(stderr," Assembly of fc*    %13.4f s ...\n",t2/1000.0);
    filePrint(stderr," Kcc*^-1 fc*        %13.4f s ...\n",this->times.project/1000.0);
    filePrint(stderr," Local Solve        %13.4f s ...\n",t4/1000.0);
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," Kcc - Krc^t Krr Krc %12.4f s ...\n",t5/1000.0);
    filePrint(stderr," Total Making Kcc    %12.4f s ...\n",this->times.coarse1/1000.0);
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," Total Making FETI  %13.4f s ...\n",t6/1000.0);
    filePrint(stderr," Total Solve RHS#   %13.4f s ...\n",t7/1000.0);
  }
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Mv, bool errorFlag) const
{
  // this function does dual mpc (CCt) preconditioning in addition to the usual feti preconditioning
  double error;
  if(mpcPrecon) {
    if(globalFlagCtc && (dualStatusChange || primalStatusChange) && fetiInfo->rebuildcct && CCtsolver)
	    const_cast<GenFetiDPSolver<Scalar>*>(this)->rebuildCCt();
    if(&Mv != &v) Mv = v; cctSolveMpc(Mv); // Mv = CCt^-1*v
    error = FetiBaseClass<Scalar>::preCondition(Mv, Mv, errorFlag);
    cctSolveMpc(Mv); 
  }
  else error = FetiBaseClass<Scalar>::preCondition(v, Mv, errorFlag);
  return (errorFlag) ? error : std::numeric_limits<double>::max();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::mergeSolution(GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc,
                                       GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &lambda) const
{
  if(ngrbms) {
    if(!geometricRbms && fetiInfo->corners != FetiInfo::noCorners) {
      GenVector<Scalar> &alpha = this->wksp->ret_alpha();
      for(int i= 0; i<KccSolver->numRBM(); ++i) 
        for(int j=0; j<KccSolver->neqs(); ++j) 
          uc[j] += alpha[i]*kccrbms[i*KccSolver->neqs()+j];
      for(int i=0; i<this->nsub;++i)
	      this->subdomains[i]->addTrbmRalpha(kccrbms, KccSolver->numRBM(), KccSolver->neqs(), alpha.data(), ur.subData(i));
    }
  }

  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, ur, uc, u, lambda);
  if(ngrbms) { // add rigid body modes
    if(geometricRbms || fetiInfo->corners == FetiInfo::noCorners) {
      GenVector<Scalar> &alpha = this->wksp->ret_alpha(); // amplitudes of the rigid body modes
      execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::addRalpha, u, alpha);
    }
    if(fetiInfo->uproj) computeProjectedDisplacement(u);
    else filePrint(stderr, " ... Do not project the displacement ...\n");
  }
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::mergeUr(int iSub, GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, 
                                 GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &lambda) const
{
  this->subdomains[iSub]->mergeUr(ur.subData(iSub), uc.data(), u.subData(iSub), lambda.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda,
                                           GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
                                           GenVector<Scalar> &uc, GenDistrVector<Scalar> &r,
                                           GenDistrVector<Scalar> &fw) const
{
 startTimerMemory(this->times.sAndJ, this->times.memorySAndJ);
 r.zero();

 GenVector<Scalar> &FcStar(uc); 

 if(KccSolver || KccParallelSolver) { 
   // Step 1: fc^*(s) = fc^(s) - (Krc^T Krr^-1)^(s) (fr^(s) - Br^(s)T lambda) - Bc^(s)tilde^T lambda
   t1 -= getTime();
   execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::makeFc, fr, lambda);
   FcStar = fc;
   t1 += getTime();

   // Step 2: Assemble local fc^* into global FcStar
   t2 -= getTime();
   assembleFcStar(FcStar);
   t2 += getTime();

   if(KccParallelSolver) {
     GenStackDistVector<Scalar> DistFcStar(*coarseInfo, FcStar.data());
     // Step 3: Solve uc = Kcc^-1 FcStar
     this->times.project -= getTime();
     KccParallelSolver->reSolve(DistFcStar);
     this->times.project += getTime();
   }
   else {
#ifdef DISTRIBUTED
     this->fetiCom->globalSum(FcStar.size(), FcStar.data());
#endif

     // Step 3: Solve uc = Kcc^-1 FcStar
     this->times.project -= getTime();
     if(this->glNumMpc_primal > 0) execParal(this->mpcToSub_primal->csize(), this, &GenFetiDPSolver<Scalar>::addMpcRHS, FcStar.data()); 
     KccSolver->reSolve(FcStar);  // now Fcstar is uc;
     this->times.project += getTime();
   }
 }
 else FcStar.zero();

 // Step 3.5: fr^(s) = fr^(s) - Krc^(s) uc^(s)
 t4 -= getTime();
 GenDistrVector<Scalar> &fr2 = this->wksp->ret_fr2(); fr2 = fr;
 execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, fr2, FcStar);

 // Step 4:
 if(domain->solInfo().isCoupled && !fetiInfo->fsi_element) {
   timedParal(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled1,
                ur, r, fr2, lambda, FcStar);

   if (fetiInfo->fsi_corner == 0) this->wiPat->exchange();

   timedParal(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled2,
                ur, r, fr2, lambda, FcStar, fw);
 }
 else 
   timedParal(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolve,
                ur, r, fr2, lambda, FcStar); // FcStar is necessary for contact

 t4 += getTime();
 this->vPat->exchange();
 timedParal(this->times.solveAndJump, this->nsub, this, &FetiBaseClass<Scalar>::interfaceDiff, r);

 // Step 5:
 if(this->glNumMpc > 0) execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::subtractMpcRhs, r);

 nMatVecProd++;
 stopTimerMemory(this->times.sAndJ, this->times.memorySAndJ);
}

template<class Scalar> 
Scalar
GenFetiDPSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur, 
                                           GenVector<Scalar> &duc, GenDistrVector<Scalar> &Fp) const
{
 startTimerMemory(this->times.sAndJ, this->times.memorySAndJ);

 GenVector<Scalar> &FcStar(duc);
 FcStar.zero();
 if(KccSolver || KccParallelSolver) {
   // Step 1: fc^(s) = - (Krc^T Krr^-1)^(s) ( Br^(s)T * p) + Bc^(s)T * p
   t1 -= getTime();
   execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::makeFcB, p);
   t1 += getTime();

   // Step 2: Assemble fc^(s) into FcStar
   t2 -= getTime();
   assembleFcStar(FcStar);
   t2 += getTime();

   if(KccParallelSolver) {
     GenStackDistVector<Scalar> DistFcStar(*coarseInfo, FcStar.data());
     // Step 3: Solve uc = Kcc^-1 FcStar
     this->times.project -= getTime();
     KccParallelSolver->reSolve(DistFcStar);
     this->times.project += getTime();
   }
   else {
#ifdef DISTRIBUTED
     this->fetiCom->globalSum(FcStar.size(), FcStar.data());
#endif

     // Step 3: Solve uc = Kcc^-1 FcStar
     this->times.project -= getTime();
     KccSolver->reSolve(FcStar);
     this->times.project += getTime();
   }
 }

 // Step 3.5: fr2^(s) = - Krc^(s) uc^(s)
 t4 -= getTime();
 GenDistrVector<Scalar> &fr2 = this->wksp->ret_fr2(); fr2.zero();
 execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, fr2, FcStar);

 // Step 4:
 // CKT : fr=Krc (-uc) 
 // CKT : Fp = B Krr^-1 [B^T p + Krc (-uc)] - Bc (-uc)
 // B contains both usual Br u = +- u and ctc restriction B u = +- u.n
 // Bc uc < > 0 if some corners nodes are in ctc. Then Bc uc = +- uc.n
 if(domain->solInfo().isCoupled && !fetiInfo->fsi_element) {
   timedParal(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled1,
                dur, Fp, fr2, p, FcStar);

   if (fetiInfo->fsi_corner == 0) this->wiPat->exchange();
   timedParal(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled2b,
                dur, Fp, fr2, p, FcStar);
 }
 else
   timedParal(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolve,
                dur, Fp, fr2, p, FcStar);
 t4 += getTime();
 this->vPat->exchange();
 timedParal(this->times.solveAndJump, this->nsub, this, &FetiBaseClass<Scalar>::interfaceDiff, Fp);
 Scalar ret = p*Fp;
#ifdef DEBUG_FETI
 Scalar pHFp = ScalarTypes::conj(ret); // note: p*Fp = (Fp)^H p therefore p^H Fp = conj(p*Fp)
 if(this->myCPU == 0 && (fetiInfo->outerloop == FetiInfo::CG)
    && (ScalarTypes::Real(pHFp) < 0.0 || fabs(ScalarTypes::Imag(pHFp)) > 1.0e-10))
   std::cerr << " *** WARNING: x^H F x = " << pHFp << ", must be positive and real for any x when F is Hermitian and positive definite. CG may not work \n";
#endif
 nMatVecProd++;
 stopTimerMemory(this->times.sAndJ, this->times.memorySAndJ);
 return ret;
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::assembleFcStar(GenVector<Scalar> &FcStar) const
{
 int i, iSub;
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   int numCornerDofs = this->subdomains[iSub]->numCoarseDofs();
   const auto &dofs = this->subdomains[iSub]->getCornerEqNums();
   const auto&fc = this->subdomains[iSub]->getfc();  // returns sub fcstar, not condensed
   for(i = 0; i < numCornerDofs; ++i) {  // assemble global condensed fcstar
     if(dofs[i] != -1) {
       FcStar[dofs[i]] += fc[i];
     }
   }
 }
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::multKrc(int iSub, GenDistrVector<Scalar> &fr, const GenVector<Scalar> &uc) const
{
  this->subdomains[iSub]->multKrc(fr.subData(iSub), uc.data());
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeFc(int iSub, const GenDistrVector<Scalar> &fr, /*GenVector<Scalar> &fc,*/
                                const GenDistrVector<Scalar> &lambda) const
{
  this->subdomains[iSub]->multfc(fr.subVec(this->subdomains[iSub]->localSubNum()), /*fc.data(),*/
                         lambda.subVec(this->subdomains[iSub]->localSubNum()));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeFcB(int iSub, GenDistrVector<Scalar> &p) const
{
  this->subdomains[iSub]->multFcB(p.subData(this->subdomains[iSub]->localSubNum()));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeEdgeConnectivity()
{ 
  int iSub;
//  int glNumSub = this->subToSub->csize();

  paralApply(this->subdomains, &FetiBaseSub::findEdgeNeighbors);

  // First count number of edges per subdomain
  int *cx = new int[this->glNumSub+1];
  int *cxx = new int[this->glNumSub+1];
  for(iSub=0; iSub<this->glNumSub+1; ++iSub) { cx[iSub] = 0; cxx[iSub] = 0; }

  // count edges in parallel
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::countEdges, cx);

#ifdef DISTRIBUTED
  this->fetiCom->globalSum(this->glNumSub+1, cx);
#endif

  // Sum each subdomain's edge count to compute total number of edges
  // We modify 'edges' at the same time so that it has the first edge number
  // for each subdomain
  int numEdges      = 0;
  int totalNumEdges = 0;
#ifdef DISTRIBUTED
  int i;
  int *numEdgesPerSub = new int[this->glNumSub];
  for(iSub=0; iSub<this->glNumSub; ++iSub)
    numEdgesPerSub[iSub] = 0;

  for(iSub=0; iSub<this->nsub; ++iSub)  
    numEdgesPerSub[this->subdomains[iSub]->subNum()] = this->subdomains[iSub]->numEdgeNeighbors();  // don't include virtual neighbors

  this->fetiCom->globalSum(this->glNumSub, numEdgesPerSub);

  for(iSub=0; iSub<this->glNumSub; ++iSub) {
    int tmp = numEdges;
    numEdges += cx[iSub];
    cx[iSub] = tmp;

    cxx[iSub] = totalNumEdges;
    totalNumEdges += numEdgesPerSub[iSub];
  }
  // make sure I can delete this
  delete [] numEdgesPerSub;
#else
  for(iSub=0; iSub<this->glNumSub; ++iSub) {
    int tmp = numEdges;
    numEdges += cx[iSub];
    cx[iSub] = tmp;

    cxx[iSub] = totalNumEdges;
    totalNumEdges += this->subdomains[iSub]->numEdgeNeighbors(); // don't include virtual neighbors
  }
#endif
  // Find total number of edges including duplicated ones.
  // totalNumEdges = 2 * numEdges (i.e. totalNumEdges counts the duplicate
  // edges since 2 subdomains see the same edge)
  // filePrint(stderr,"nE %d tnE %d\n",numEdges, totalNumEdges);

  int *connect = new int[totalNumEdges];
  cxx[this->glNumSub] = totalNumEdges;

#ifdef DISTRIBUTED
  for(i=0; i<totalNumEdges; ++i) connect[i] = -1;
#endif

  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::numberEdges, cx, cxx, connect, this->sPat);
  delete [] cx;

 // Due to the -1, how could I do a globalSum?
#ifdef DISTRIBUTED
  for(i=0; i<totalNumEdges; ++i) connect[i] += 1;
  this->fetiCom->globalSum(totalNumEdges, connect);
  for(i=0; i<totalNumEdges; ++i) connect[i] -= 1;
#endif

  this->sPat->exchange();
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::receiveNeighbEdgeNums, cxx, connect, this->sPat);

#ifdef DISTRIBUTED
  int *connect2 = new int[totalNumEdges];
  for(i=0; i<totalNumEdges; ++i) {
    if(connect[i] == -1) 
      connect2[i] = 0;
    else
      connect2[i] = connect[i]; 
  }
  this->fetiCom->globalSum(totalNumEdges, connect2);
  for(i=0; i<totalNumEdges; ++i) {
    if(connect[i] == -1) 
     connect[i] = connect2[i];
  }
  delete [] connect2;
#endif
  this->subToEdge = new Connectivity(this->glNumSub, cxx, connect);
  // create the edge to subdomain connectivity
  this->edgeToSub = this->subToEdge->alloc_reverse();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::countEdges(int iSub, int *edges) const
{
  // get my global subdomain number
  int myNum = this->subdomains[iSub]->subNum();

  // count a subdomain's edges
  edges[myNum] = 0;
  int j;
  int numNeighbor = this->subdomains[iSub]->numNeighbors();
  for(j=0; j<numNeighbor; ++j) {
    int subI = this->subdomains[iSub]->getSComm()->subNums[j];
    if(this->subdomains[iSub]->isEdgeNeighbor(j) && (myNum < subI))
      edges[myNum] += 1;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::numberEdges(int iSub, int *eP, int *ep2, int* edges, FSCommPattern<int> *sPat)
{
  // first compute each subdomain's starting edge number
  int myNum = this->subdomains[iSub]->subNum();
  int startI = eP[myNum];
  int fP = ep2[myNum];

  // number all subdomain's edges
  int jSub;
  int numNeighbor = this->subdomains[iSub]->numNeighbors();
  for(jSub=0; jSub<numNeighbor; ++jSub) {
    int subJ = this->subdomains[iSub]->getSComm()->subNums[jSub];
    FSSubRecInfo<int> sInfo = this->sPat->getSendBuffer(myNum, subJ);
    if(this->subdomains[iSub]->isEdgeNeighbor(jSub)) {
      if(myNum < subJ) {
        edges[fP] = startI;
        startI++;
      }
      else edges[fP] = -1;
      // Send the numbered edges to all neighbors so that they can complete
      // the edge vector
      sInfo.data[0] = edges[fP];
      fP++;
    }
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::receiveNeighbEdgeNums(int iSub, int *eP, int* edges, FSCommPattern<int> *sPat)
{
  int myNum = this->subdomains[iSub]->subNum();
  // number all subdomain's edges
  int jSub;
  int numNeighbor = this->subdomains[iSub]->numNeighbors();
  int jEdgeN = 0;
  for(jSub=0; jSub<numNeighbor; ++jSub) {
    if(this->subdomains[iSub]->isEdgeNeighbor(jSub)) {
      FSSubRecInfo<int> rInfo = this->sPat->recData(this->subdomains[iSub]->getSComm()->subNums[jSub], myNum);
      int en =  rInfo.data[0];
      if(en >= 0) edges[eP[myNum]+jEdgeN] = en;
      jEdgeN++;
    }
  }
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::factorLocalMatrices(int iSub)
{
  auto &K = this->subdomains[iSub]->Krr;
  if(K) {
    K->setPrintNullity(false);
    K->factor(); 
    if(K->numRBM() != 0) {
      filePrint(stderr," ... Subdomain %3d found %3d ZEMs   ...\n",
                this->subdomains[iSub]->localSubNum()+1,K->numRBM());
    }
  }
  this->subdomains[iSub]->factorKii();
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subdomainSolve(int iSub, GenDistrVector<Scalar> &v1, 
                                        GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                        GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5) const
{
  int sn = this->subdomains[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  Scalar *localsrc  = v3.subData(sn);
  Scalar *interfsrc = v4.subData(sn);
  Scalar *uc        = v5.data(); // Necessary for contact

  int localLen = this->subdomains[iSub]->localRLen();
  int i;

  for(i = 0; i < localLen; ++i) 
     localvec[i] = localsrc[i];

  int interfaceLen = this->subdomains[iSub]->interfLen();
  for(i = 0; i < interfaceLen; ++i) 
     interfvec[i] = interfsrc[i];

  this->subdomains[iSub]->fetiBaseOp(uc, this->subdomains[iSub]->Krr.get(), localvec, interfvec);
  
  this->subdomains[iSub]->sendInterf(interfvec, this->vPat);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subdomainSolveCoupled1(int iSub, GenDistrVector<Scalar> &v1,
                                                GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                                GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5) const
{
  int sn = this->subdomains[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  Scalar *localsrc  = v3.subData(sn);
  Scalar *interfsrc = v4.subData(sn);

  int localLen = this->subdomains[iSub]->localRLen();
  int i;

  for(i = 0; i < localLen; ++i)
     localvec[i] = localsrc[i];

  int interfaceLen = this->subdomains[iSub]->interfLen();
  for(i = 0; i < interfaceLen; ++i)
     interfvec[i] = interfsrc[i];

  this->subdomains[iSub]->fetiBaseOpCoupled1(this->subdomains[iSub]->Krr.get(), localvec, interfvec, this->wiPat);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1,
                                                GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                                GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5,
                                                GenDistrVector<Scalar> &fw) const
{
  int sn = this->subdomains[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  Scalar *uc        = v5.data();
  Scalar *localfw   = fw.subData(sn);

  this->subdomains[iSub]->fetiBaseOpCoupled2(uc, localvec, interfvec, this->wiPat, localfw);
  this->subdomains[iSub]->sendInterf(interfvec, this->vPat);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subdomainSolveCoupled2b(int iSub, GenDistrVector<Scalar> &v1,
                                                GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                                GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5) const
{
  int sn = this->subdomains[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  Scalar *uc        = v5.data();
  Scalar *localfw   = 0;

  this->subdomains[iSub]->fetiBaseOpCoupled2(uc, localvec, interfvec, this->wiPat, localfw);
  this->subdomains[iSub]->sendInterf(interfvec, this->vPat);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::clean_up()
{
}

template<class Scalar> 
double
GenFetiDPSolver<Scalar>::getFNormSq(GenDistrVector<Scalar> &f)
{
	// this is used by nonlinear analysis, need to assemble force residual on subdomain interface for error norm
	GenDistrVector<Scalar> &fr = this->wksp->ret_fr();
	// extract fr from f
	auto perSubExtract = [&](int iSub) {
		fr.subVec(iSub).setZero();
		this->subdomains[iSub]->getFr(f.subData(iSub), fr.subData(iSub));
	};
	threadManager->callParal(this->nsub, perSubExtract);
	this->distributeForce(fr);
	GenVector<Scalar> &fc  = this->wksp->ret_fc();
	getFc(f, fc);
#ifdef DISTRIBUTED
	this->fetiCom->globalSum(fc.size(), fc.data());
#endif
	double mpcerr = 0.0;

	for(int i=0; i<this->nsub; ++i) mpcerr += this->subdomains[i]->getMpcError();
#ifdef DISTRIBUTED
	mpcerr = this->fetiCom->globalSum(mpcerr);
#endif

	return (fr.sqNorm() + fc.sqNorm() + mpcerr);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subtractMpcRhs(int iSub, GenDistrVector<Scalar> &dv1) const
{
  this->subdomains[iSub]->subtractMpcRhs(dv1.subData(this->subdomains[iSub]->localSubNum()));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V)
{
  int info = 0;
  int mindim = std::min(nrow,ncol);
  int maxdim = std::max(nrow,ncol);
  double max_value = A.maxAbs();
#ifdef FILERING
  for(int i=0; i<A.numCol()*A.numRow(); i++)
    if(fabs((A.data())[i])<1.E-16*max_values+1.E-20) (A.data())[i] = 0.0;
#endif
  double *w    = (double*) dbg_alloca(sizeof(double)* maxdim);
  double *e    = (double*) dbg_alloca(sizeof(double)* maxdim);
  double *work = (double*) dbg_alloca(sizeof(double)* maxdim);

  for(int i=0; i<maxdim; i++) { w[i] = 0.0; e[i] = 0.0; work[i] = 0.0; }
                                                                                                  
  if(V == 0)
    Tsvdc(A.data(), ncol, ncol, nrow, w, e, U.data(), ncol, U.data(), ncol, work, 10, info);
  else 
    Tsvdc(A.data(), ncol, ncol, nrow, w, e, U.data(), ncol, V->data(), nrow, work, 11, info);

  //  ... CHECK RETURN STATUS OF DSVDC:
  if(info <  0)
    filePrint(stderr," *** ERROR: Illegal value in argument #%d in FetiDPSolver::singularValueDecomposition()\n",info);
  if(info >  0)
    filePrint(stderr," *** WARNING: %d diagonals did not converge in FetiDPSolver::singularValueDecomposition(). Check your result.\n",info);

  //  ... DETERMINE RANK
  rank = 0;
  double tolerance = max_value*tol;
  for(int i=0; i<mindim; ++i) {
    if(fabs(w[i]) > tolerance) rank += 1;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::rebuildGtGtilda()
{
  nRebuildGtG++;
  startTimerMemory(this->times.coarse1, this->times.memoryGtG);

  if(GtGtilda == NULL) {
	  GtGtilda = GenSolverFactory<Scalar>::getFactory()->createSolver(
		  &coarseConnectGtG, eqNumsGtG, *fetiInfo->auxcoarse_cntl, GtGsparse);
    GtGtilda->setPrintNullity(fetiInfo->contactPrintFlag && this->myCPU == 0);
  } else
  GtGtilda->zeroAll();
  execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::assembleGtG);
#ifdef DISTRIBUTED
  GtGtilda->unify(this->fetiCom);
#endif
  startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  GtGtilda->parallelFactor();
  //int* pivnull_list = ((GenMumpsSolver<Scalar>*) GtGtilda)->getPivnull_list();
  stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);

  stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::assembleGtG(int iGroup)
{
  // assembles groups in parallel, subdomains with same group sequentially 
  // threadsafe implementation - avoids simultaneous writing to same memory
  // note: the distributed version will work for shared memory too, but the 
  // alternative code is a bit more efficient
  int i;
#ifdef DISTRIBUTED
  for(i = 0; i < this->nsub; ++i) {
    if(this->subdomains[i]->getGroup() == groups[iGroup])
      this->subdomains[i]->assembleGtGsolver(GtGsparse);
  }
#else
  auto grsubs = (*groupToSub)[iGroup];
  for(i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->subdomains[iSub]->assembleGtGsolver(GtGsparse);
  }
#endif
}

template<class Scalar>
void 
GenFetiDPSolver<Scalar>::buildCCt()
{
  // build CC^t for preconditioning mpc residual
  startTimerMemory(this->times.buildCCt, this->times.memoryBuildCCt);
  // find subs with mpcs
  mpcSubMap = new int[this->nsub];
  numSubsWithMpcs = 0;
  for(int i=0; i<this->nsub; ++i) {
    if(this->subdomains[i]->numMPC > 0) {
      mpcSubMap[i] = numSubsWithMpcs;
    }
    else mpcSubMap[i] = -1;
  }

  switch(fetiInfo->mpc_precno) {
    case (FetiInfo::globalCCt) :
      CCtsolver = new GlobalCCtSolver<Scalar>(mpcToMpc, mpcToCpu, numSubsWithMpcs, this->subdomains,
                                              fetiInfo, this->fetiCom);
      break;
    case (FetiInfo::blockDiagCCt) : {
      Connectivity *blockToMpc = getBlockToMpc();
      CCtsolver = new BlockCCtSolver<Scalar>(blockToMpc, mpcToMpc, this->mpcToSub, mpcToCpu, this->numSubsWithMpcs,
                                             this->subdomains,
                                             mpcSubMap, fetiInfo, this->fetiCom);
      } break;
    case (FetiInfo::subBlockDiagCCt) :
      CCtsolver = new SubBlockCCtSolver<Scalar>(mpcToMpc, this->mpcToSub, numSubsWithMpcs, this->subdomains,
                                                this->fetiCom, this->cpuToSub);
      break;
    case (FetiInfo::superBlockDiagCCt) : {
      Connectivity *blockToMpc = getBlockToMpc();
      bool super_flag = (fetiInfo->mpc_block == FetiInfo::subBlock) ? false : true;
      bool sub_flag = (fetiInfo->mpc_block == FetiInfo::mortarBlock) ? true : false;
      CCtsolver = new SuperBlockCCtSolver<Scalar>(blockToMpc, mpcToMpc, this->mpcToSub, mpcToCpu, numSubsWithMpcs,
                                                  this->subdomains,
                                                  fetiInfo, this->fetiCom, super_flag, sub_flag);
    } break;
    default :
      std::cerr << " *** ERROR: don't know mpc_precno = " << fetiInfo->mpc_precno << std::endl;
      break;
  }
  CCtsolver->assemble();
  CCtsolver->factor();
  stopTimerMemory(this->times.buildCCt, this->times.memoryBuildCCt);
}

template<class Scalar>
Connectivity *
GenFetiDPSolver<Scalar>::getBlockToMpc()
{
 // return a decomposition of global mpcs into blocks for block solver
 Connectivity *blockToMpc;
 switch(fetiInfo->mpc_block) {
   case (FetiInfo::topoBlock) :
     {
       compStruct renumber = mpcToMpc->renumByComponent(-1);  // -1 = no renumbering (optimized implementation)
       int nMpcBlocks = renumber.numComp;
       renumber.order = new int[this->glNumMpc];
       for(int i=0; i<this->glNumMpc; ++i) renumber.order[renumber.renum[i]] = i;
       int *target = new int[this->glNumMpc];
       for(int i=0; i<this->glNumMpc; ++i) target[renumber.renum[i]] = i;
       int *pointer = new int[nMpcBlocks+1];
       for(int i=0; i<=nMpcBlocks; ++i) pointer[i] = renumber.xcomp[i];
       renumber.clearMemory();
       blockToMpc = new Connectivity(nMpcBlocks, pointer, target);
     } 
     break;
   case (FetiInfo::subBlock) :
     { 
       //HB: WARNING: HARD CODED FOR SHARED MEMORY !!!
       Connectivity *subToMpc = this->mpcToSub->alloc_reverse();
       int nMpcBlocks = 0; 
       int ntarget = 0;
       for(int s=0;s<subToMpc->csize();s++){
         if(subToMpc->num(s)) { 
           nMpcBlocks++; 
           ntarget += subToMpc->num(s);
         }
       }
       int *ptr = new int[nMpcBlocks+1]; ptr[0] = 0;
       int *target = new int[ntarget];
       int iblk = 0;
       for(int s=0;s<subToMpc->csize();s++){
         if(subToMpc->num(s)) { 
           ptr[iblk+1] = ptr[iblk];
           for(int i=0;i<subToMpc->num(s);i++){
             target[ptr[iblk]+i] = (*subToMpc)[s][i];
             ptr[iblk+1]++;
           } 
           iblk++;
         }
       }
       blockToMpc = new Connectivity(nMpcBlocks,ptr,target);
       delete subToMpc;
     } 
     break;
   case (FetiInfo::mortarBlock) :
     {
       // make initial blockToMpc connectivity: make one block for each set of MPCs representing 
       // a mortar face, and put all remaining (non-mortar) MPCs together in one block 
       int numMortarMpcs = domain->GetnMortarLMPCs();
       int numOtherMpcs = this->glNumMpc - numMortarMpcs;
       if(numOtherMpcs > 0) {
         int *pointer = new int[2]; pointer[0] = 0; pointer[1] = numOtherMpcs;
         int *target = new int[numOtherMpcs];
         for(int i=0; i<numOtherMpcs; ++i) target[i] = i;
         Connectivity *otherToMpc = new Connectivity(1, pointer, target);
         if(numMortarMpcs == 0) blockToMpc = otherToMpc;
         else blockToMpc = new Connectivity{ otherToMpc->append(*domain->GetMortarToMPC()) };
       }
       else blockToMpc = new Connectivity(*domain->GetMortarToMPC());
     } 
     break;
 }

 // block overlapping if requested
 if(blockToMpc->csize() > 1) {
   int blockOverlapLevel = fetiInfo->mpcBlkOverlap;
   filePrint(stderr, " ... Block Overlap Level %d ... \n", blockOverlapLevel);
   for(int i = 0; i < blockOverlapLevel; ++i) {
     Connectivity *tmpBlockToMpc = blockToMpc->transcon(mpcToMpc);
     delete blockToMpc; blockToMpc = tmpBlockToMpc;
   }
 }

 return blockToMpc;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::rebuildCCt()
{
  CCtsolver->zeroAll();
  CCtsolver->assemble();
  CCtsolver->factor();
  nRebuildCCt++;
}

template<class Scalar>
void 
GenFetiDPSolver<Scalar>::cctSolveMpc(GenDistrVector<Scalar> &v) const
{
  // this function is used when applying the generalized preconditioner
  // to the mpc part of the residual
  // computes:  v = (CC^t)^-1 v
  startTimerMemory(this->times.precond, this->times.memoryPrecond);
  startTimerMemory(this->times.solveCCt, this->times.memorySolveCCt);
  
  CCtsolver->reSolve(v); 

  stopTimerMemory(this->times.solveCCt, this->times.memorySolveCCt);
  if(this->times.memorySolveCCt < 0) this->times.memorySolveCCt = 0;
  stopTimerMemory(this->times.precond, this->times.memoryPrecond); 
}

template<class Scalar>
GenFetiDPSolver<Scalar>::~GenFetiDPSolver() 
{
  if(ngrbmGr) { delete [] ngrbmGr; ngrbmGr = 0; }
  if(groups) { delete [] groups; groups = 0; }
  if(groupToSub != bodyToSub) { delete groupToSub; groupToSub = 0; }
  if(subToGroup && (subToGroup != subToBody)) { delete subToGroup; subToGroup = 0; }
  if(subToBody) { delete subToBody; subToBody = 0; }
  if(GtGtilda) { delete GtGtilda; GtGtilda = 0; }
  //if(glCrnGroup) { delete [] glCrnGroup; glCrnGroup = 0; }
  
  if(KccSparse) { delete KccSparse; KccSolver = 0; KccSparse = 0; }
  if(cornerToSub) { delete cornerToSub; cornerToSub = 0; }
  if(cornerEqs) { delete cornerEqs; cornerEqs = 0; }
  if(eqNumsGtG) { delete eqNumsGtG; eqNumsGtG = 0; }
  
  if(CCtsolver) { delete CCtsolver; CCtsolver = 0; }
  if(mpcSubMap) { delete [] mpcSubMap; mpcSubMap = 0; }
  if(mpcPat) { delete mpcPat; mpcPat = 0; }
  if(kccrbms) { delete [] kccrbms; kccrbms = 0; }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getLocalMpcForces(int iSub, double *mpcLambda)
{
  // mpcLambda is the local MPC forces for subdomain iSub (required for salinas)
  GenVector<Scalar> &uc = this->wksp->ret_uc();
  this->subdomains[iSub]->getLocalMpcForces(mpcLambda, cornerEqs, mpcOffset, uc);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::addMpcRHS(int iMPC, Scalar *fcstar) const
{
 int dof = cornerEqs->firstdof(mpcOffset+iMPC);
 int myNum = this->glSubToLoc[(*this->mpcToSub_primal)[iMPC][0]];
 if(myNum >= 0)
   fcstar[dof] += this->subdomains[myNum]->getMpcRhs_primal(iMPC);
}

template<class Scalar>
int
GenFetiDPSolver<Scalar>::numRBM()
{
  bool useKccSolver = (this->glNumMpc == 0 && !geometricRbms);
  if(GtGtilda && !useKccSolver) {
    return GtGtilda->numRBM();
  }
  else return (KccSolver) ? KccSolver->numRBM() : 0;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::wetInterfaceComms()
{
  // send and receive numWIdof
  FSCommPattern<int> *wiOnePat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
  for(int iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setWIoneCommSize(wiOnePat);
  wiOnePat->finalize();
  paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::sendNumWIdof, wiOnePat);
  wiOnePat->exchange();
  paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::recvNumWIdof, wiOnePat);
  delete wiOnePat;

  // send and receive glToLocalWImaps
  FSCommPattern<int> *wiMapPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend,
                                                        FSCommPattern<int>::NonSym);
  for(int iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setWImapCommSize(wiMapPat);
  wiMapPat->finalize();
  paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::sendWImap, wiMapPat);
  wiMapPat->exchange();
  paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::recvWImap, wiMapPat);
  delete wiMapPat;

  // create this->wiPat FSCommPattern object, used to send/receive wet this->interface interaction vectors
  this->wiPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                    FSCommPattern<Scalar>::NonSym);
  for(int iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setWICommSize(this->wiPat);
  this->wiPat->finalize();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::reconstruct()
{
  // 1. reset the orthoset
  this->resetOrthoSet();
 
  if(domain->solInfo().isCoupled) domain->computeCoupledScaleFactors();

  if(fetiInfo->isEdgeAugmentationOn() && fetiInfo->numdir > 0) { // augmentation depends on freq/k
    // 2. reconstruct local Kcc etc. since size of Kcc may have changed
    startTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
    paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::constructKcc);
    if(domain->solInfo().isCoupled)  paralApply(this->subdomains, &FetiSub<Scalar>::constructKcw);
    stopTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);

    // 3. rebuild augmentation Q
    if(verboseFlag) filePrint(stderr," ... Rebuild Edge Augmentation (Q)  ... \n");
    if(fetiInfo->waveMethod != FetiInfo::averageMat) computeLocalWaveNumbers();
    paralApplyToAll(this->subdomains, &FetiBaseSub::zeroEdgeDofSize);
    paralApplyToAll(this->subdomains, &FetiSub<Scalar>::makeQ);  // rebuild augmentation matrix
  }

  geometricRbms = 0;
  ngrbms = 0;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::refactor()
{
  // 4.  stiffness scaling && edge weighting -> check !!
  if(fetiInfo->scaling == FetiInfo::kscaling) {
    execParal(this->nsub, this, &FetiBaseClass<Scalar>::sendScale);
    this->vPat->exchange();
    paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectScaling, this->vPat);
  }

  if(domain->solInfo().isCoupled) {
    paralApplyToAll(this->subdomains, &FetiSub<Scalar>::reScaleAndReSplitKww);
  }

  if(fetiInfo->augment == FetiInfo::WeightedEdges)
    paralApplyToAll(this->subdomains, &FetiSub<Scalar>::weightEdgeGs); // W*Q

 // MPCs 
 mpcPrecon = false;
 if(this->glNumMpc > 0) {
   if(fetiInfo->mpc_scaling == FetiInfo::kscaling) { // MPC stiffness scaling
     FSCommPattern<Scalar> *mpcDiagPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU,
                                                                   FSCommPattern<Scalar>::CopyOnSend);
     for(int iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setMpcDiagCommSize(mpcDiagPat);
     mpcDiagPat->finalize();
     paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::sendMpcDiag, mpcDiagPat);
     mpcDiagPat->exchange();
     paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectMpcDiag, mpcDiagPat);
     delete mpcDiagPat;
   }

   if(fetiInfo->c_normalize) normalizeC();

   if(fetiInfo->mpc_precno == FetiInfo::diagCCt) {
     // use W scaling for preconditioning mpcs, don't need to build & invert CC^t
     paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::sendMpcScaling, this->vPat);
     this->vPat->exchange();  // neighboring subs mpc weights
     paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectMpcScaling, this->vPat);
   }
   else if(fetiInfo->mpc_precno != FetiInfo::noMpcPrec) {
     // used generalized proconditioner for mpcs, need to build and invert CC^t
     mpcPrecon = true;
     paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::initMpcScaling);
     deleteCCt();
     buildCCt();
   }
 }

  // 5. factor local matrices
  if(verboseFlag) filePrint(stderr," ... Factor Subdomain Matrices      ... \n");
  startTimerMemory(this->times.factor, this->times.memoryFactor);
  if(fetiInfo->local_cntl->subtype != FetiInfo::spooles && fetiInfo->local_cntl->subtype != FetiInfo::mumps)
    timedParal(this->times.factorMat, this->nsub, this, &GenFetiDPSolver<Scalar>::factorLocalMatrices);
  else // spooles and mumps are not thread safe
    for(int iSub=0; iSub<this->nsub; ++iSub) factorLocalMatrices(iSub);
  stopTimerMemory(this->times.factor, this->times.memoryFactor);

  if(verboseFlag) filePrint(stderr, " ... Reconstruct Kcc solver         ... \n");
  delete KccSparse; KccSparse = 0; KccSolver = 0;
  makeKcc();

  if(GtGtilda) { delete GtGtilda; GtGtilda = 0; } 
  if(ngrbms) makeGtG();

  delete this->wksp;
  this->times.memoryDV -= memoryUsed();
  int numC = (KccSolver) ? KccSolver->neqs() : 0;
  this->wksp = new GenFetiWorkSpace<Scalar>(this->interface, internalR, internalWI, ngrbms, numC, globalFlagCtc);
  this->times.memoryDV += memoryUsed();

  int tLocalCLen = 0;
  for(int iSub = 0; iSub < this->nsub; ++iSub) {
    internalC.domLen[iSub] = this->subdomains[iSub]->numCoarseDofs();
    tLocalCLen += internalC.domLen[iSub];
  }
  internalC.len = tLocalCLen;

  threadManager->callParal(this->nsub, [this](int iSub) { this->subdomains[iSub]->makeBs(); });

  paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::initMpcStatus);

}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::projectActiveIneq(const GenDistrVector<Scalar> &x, GenDistrVector<Scalar> &y) const
{
	if(&x != &y) y = x;
	threadManager->callParal(this->nsub, [&](int iSub) {
		this->subdomains[iSub]->projectActiveIneq(y.subData(this->subdomains[iSub]->localSubNum()));
	});
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::multG(const GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, double beta) const
{
	// y = alpha*G*x + beta*y
	if(beta == 0) y.zero(); else y *= beta;
	threadManager->callParal(this->nsub, [&](int iSub) {
		this->subdomains[iSub]->multG(x, y.subData(this->subdomains[iSub]->localSubNum()), alpha);
	});
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::trMultG(const GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, double beta) const
{
  // y = alpha*G^T*x + beta*y
  if(beta == 0) y.zero(); else y *= beta;
  GenVector<Scalar> v(ngrbms, 0.0);
  execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::subTrMultG, x, v, alpha); // v += alpha*G^T*x
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(ngrbms, v.data());
#endif
  y += v;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subTrMultG(int iGroup, const GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha) const
{
#ifdef DISTRIBUTED
  for(int i=0; i<this->nsub; ++i) {
    if(this->subdomains[i]->getGroup() == groups[iGroup])
      this->subdomains[i]->trMultG(x.subData(this->subdomains[i]->localSubNum()), y, alpha);
  }
#else
  auto grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->subdomains[iSub]->trMultG(x.subData(this->subdomains[iSub]->localSubNum()), y, alpha);
  }
#endif
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::project(GenDistrVector<Scalar> &z, GenDistrVector<Scalar> &y, int eflag) const
{
  startTimerMemory(this->times.project, this->times.memoryProject1);
  // if eflag is true this function computes y as the projection of z on to the feasible subset in which all active constraints remain active
  // if eflag is false this function computes y as the projection of z on to the tangent subspace

  paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::saveMpcStatus1);

  GenDistrVector<Scalar> x(z);
  bool status_change = false;

  int i;
  for(i = 0; i < fetiInfo->dual_plan_maxit+1; ++i) {
    // y = P_i*x
    projectActiveIneq(x, y);

    // res = -e-G^T*x
    GenVector<Scalar> res(ngrbms);
    if(ngrbms) { 
      trMultG(y, res, -1.0, 0.0);
      if(eflag) { 
        GenVector<Scalar> &e = this->wksp->ret_e();
        res -= e; 
      }
    }

    // check stopping criteria
    if(i > 0) {
      double resnorm = (eflag && ngrbms) ? res.norm() : 0;
      if(fetiInfo->contactPrintFlag && this->myCPU == 0) std::cerr << "dual planing: iteration = " << i << ", residual = " << resnorm << std::endl;
      if(/*resnorm <= fetiInfo->dual_proj_tol ||*/ !status_change) break;
      else if(i == std::max(1,fetiInfo->dual_plan_maxit)) {
        if(this->myCPU == 0) std::cerr << "warning: dual planing did not converge after " << i << " iterations. Error = " << resnorm << std::endl;
        // note: if we break the loop here then y will not be feasible wrt the equality constraints (i.e. G^T*y != e)
        // if we don't break here then y will not be feasible wrt the inequality constraints
        break;
      }
    }
  
    // res = (G^T*P_i*G)^(-1)*res
    if(GtGtilda) GtGtilda->reSolve(res);

    // x = x + G*res
    if(ngrbms) multG(res, x, 1.0, 1.0);

    // update the active set
    if(eflag && globalFlagCtc)
      status_change =
	      const_cast<GenFetiDPSolver<Scalar> *>(this)->updateActiveSet(x, 0, -fetiInfo->dual_plan_tol);
  }

  if((dualStatusChange = (i > 1))) {
    nSubIterDual += (i-1);
    nStatChDual++;
    if(fetiInfo->contactPrintFlag && this->myCPU == 0) std::cerr << std::endl;
  }

  stopTimerMemory(this->times.project, this->times.memoryProject1);
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::tProject(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w) const
{
  startTimerMemory(this->times.project, this->times.memoryProject1);
  // this function computes w as the projection of r on to the tangent cone at a feasible point
  // unless the "chopped gradient" error term is proportional, in which case this function
  // computes w as the projection of r on to the tangent subspace
  // note: this function also computes alpha which is stored in the feti workspace

  paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::saveMpcStatus1);

  GenDistrVector<Scalar> x(r);
  GenVector<Scalar> &alpha = this->wksp->ret_alpha();
  GenDistrVector<Scalar> &gc = this->wksp->ret_gc(); // chopped gradient
  GenDistrVector<Scalar> &gf = this->wksp->ret_gf(); // free gradient
  alpha.zero();
  bool status_change = false;
  bool proportional = false;

  int i;
	for(i = 0; i < fetiInfo->primal_plan_maxit+1; ++i) {
		// w = P_i*x
		projectActiveIneq(x, w);

		// res = -G^T*x
		GenVector<Scalar> res(ngrbms);
		if(ngrbms) trMultG(w, res, -1.0, 0.0);

		// check stopping criteria
		if(i > 0) {
			double resnorm = (ngrbms) ? res.norm() : 0;
			if(fetiInfo->contactPrintFlag && this->myCPU == 0) std::cerr << "primal planing: iteration " << i << ", residual = " << resnorm << std::endl;
			if(/*resnorm <= fetiInfo->primal_proj_tol ||*/ !status_change) break;
			else if(i == std::max(1,fetiInfo->primal_plan_maxit)) {
				if(this->myCPU == 0) std::cerr << "warning: primal planing did not converge after " << i << " iterations. " << std::endl;
				break;
			}
		}

		// res = (G^T*P_i*G)^(-1)*res
		if(GtGtilda) GtGtilda->reSolve(res);

		// x += G*res, alpha += res
		if(ngrbms) multG(res, x, 1.0, 1.0);
		alpha += res;

		// update active set unless error is proportional
		if(globalFlagCtc) {
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::split, x, gf, gc);
			proportional = (i == 0 && (gc.norm() <= fetiInfo->gamma*gf.norm()));
			if(!proportional)
				status_change =
					const_cast<GenFetiDPSolver<Scalar> *>(this)->updateActiveSet(x, 1, -fetiInfo->primal_plan_tol);
		}
	}

  if((primalStatusChange = (i > 1))) {
    nSubIterPrimal += (i-1);
    nStatChPrimal++;
    if(fetiInfo->contactPrintFlag && this->myCPU == 0) std::cerr << std::endl;
  }

  stopTimerMemory(this->times.project, this->times.memoryProject1);

  return w.sqNorm() + (proportional ? gc.sqNorm() : 0);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeL0(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f) const
{
  lambda.zero(); // XXXX consider starting from previous iteration or time (load) step for nonlinear dynamics (statics)
                 // currently lambda is initialized to zero in workspace constructor for feti-dp and a 
                 // new workspace is constructed every time the solver is refactored
  if(this->glNumMpc) {
    if(ngrbms > 0) makeE(f); // compute e = R^T*f
    project(lambda, lambda, true); 
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::normalizeC()
{
  GenVector<Scalar> cnorm(this->glNumMpc, 0.0);
  for(int i=0; i<this->nsub; ++i) this->subdomains[i]->normalizeCstep1(cnorm.data());
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(this->glNumMpc, cnorm.data());
#endif
  for(int i=0; i<this->glNumMpc; ++i) cnorm[i] = ScalarTypes::sqrt(cnorm[i]);
  paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::normalizeCstep2, cnorm.data());
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::trMultC(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f)
{
  // compute f = C^T*lambda
  f.zero();
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::subTrMultC, lambda, f);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subTrMultC(int iSub, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f)
{ 
  this->subdomains[iSub]->multAddCT(lambda.subData(this->subdomains[iSub]->localSubNum()), f.subData(this->subdomains[iSub]->localSubNum()));
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::multC(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu)
{
  // compute cu = C*u
  cu.zero();
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::subMultC, u, cu);

  execParal(this->nsub, this, &FetiBaseClass<Scalar>::interfSend, cu);
  this->vPat->exchange();
  execParal(this->nsub, this, &FetiBaseClass<Scalar>::interfaceDiff, cu);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subMultC(int iSub, GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu) const
{
  this->subdomains[iSub]->multC(u.subData(this->subdomains[iSub]->localSubNum()), cu.subData(this->subdomains[iSub]->localSubNum()));
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::reconstructMPCs(Connectivity *_mpcToSub, Connectivity *_mpcToMpc, Connectivity *_mpcToCpu)
{
#ifdef AS_LIBRARY
#else
 this->mpcToSub   = _mpcToSub;    // MPC to subdomain connectivity
 this->glNumMpc = (this->mpcToSub) ? this->mpcToSub->csize() : 0;
 mpcToMpc   = _mpcToMpc;    // MPC to MPC connectivity (used for CC^t preconditioner)
 mpcToCpu   = _mpcToCpu;
 globalFlagCtc = domain->getNumCTC();
#ifdef DISTRIBUTED
 globalFlagCtc = this->fetiCom->globalMax((int) globalFlagCtc);
#endif
 int iSub;

 // create this->vPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
 delete this->vPat;
 this->vPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
 for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setDofCommSize(this->vPat);
 this->vPat->finalize();
 
 // create this->sPat FSCommPattern objects, used to send/receive a single integer
 delete this->sPat;
 this->sPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
 for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setCommSize(this->sPat, 1);
 this->sPat->finalize();

 if(mpcPat) { delete mpcPat; mpcPat = 0; }
 if(globalFlagCtc) {
   mpcPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
   for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setMpcCommSize(mpcPat);
   mpcPat->finalize();
 }
 
 int tInterfLen    = 0;
 this->halfSize = 0;
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   this->interface.domLen[iSub] = this->subdomains[iSub]->interfLen();

   this->subdomains[iSub]->computeMasterFlag(*this->mpcToSub);
   this->fetiOps[iSub]->setHalfOffset(this->halfSize);
   this->halfSize      += this->subdomains[iSub]->halfInterfLen();
   tInterfLen    += this->interface.domLen[iSub];
 }
 this->interface.len = tInterfLen;

 // compute the masterFlags
 bool *interfaceMasterFlag = new bool[tInterfLen];
 this->interface.recomputeOffsets();
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   auto &subMasterFlag = this->subdomains[iSub]->getMasterFlag();
   int subOffset = this->interface.subOffset[iSub];
   int j;
   for(j=0; j<this->interface.domLen[iSub]; ++j)
     interfaceMasterFlag[subOffset+j] = subMasterFlag[j];
 }
 this->interface.setMasterFlag(interfaceMasterFlag);
 // don't delete interfaceMasterFlag

 // Allocate space for reorthogonalization set
 this->times.memoryOSet -= memoryUsed();
 if(fetiInfo->outerloop == 0) {
   delete this->oSetCG;
   this->oSetCG = (fetiInfo->maxortho > 0) ? new GenCGOrthoSet<Scalar>(this->halfSize, fetiInfo->maxortho, this->fetiCom) : 0;
 }
 else if(fetiInfo->outerloop == 1) {
   delete this->oSetGMRES;
   this->oSetGMRES = new GenGMRESOrthoSet<Scalar>(this->halfSize, fetiInfo->maxortho, this->fetiCom);
 }
 else {
   delete this->oSetGCR;
   this->oSetGCR = new GenGCROrthoSet<Scalar>(this->halfSize, fetiInfo->maxortho, this->fetiCom);
 }
 this->times.memoryOSet += memoryUsed();

 paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::rebuildKbb);

 paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::initScaling);
 if((fetiInfo->scaling == FetiInfo::kscaling) || ((fetiInfo->mpc_scaling == FetiInfo::kscaling) && (this->glNumMpc_primal > 0))
    || (fetiInfo->augment == FetiInfo::WeightedEdges)) {
   execParal(this->nsub, this, &FetiBaseClass<Scalar>::sendScale);
   this->vPat->exchange();
   paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::collectScaling, this->vPat);
 }

 deleteCCt(); mpcPrecon = false;
 if(this->glNumMpc > 0) {
   if(fetiInfo->mpc_scaling == FetiInfo::kscaling) { // MPC stiffness scaling
     FSCommPattern<Scalar> *mpcDiagPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU,
                                                                   FSCommPattern<Scalar>::CopyOnSend);
     for(iSub=0; iSub<this->nsub; ++iSub) this->subdomains[iSub]->setMpcDiagCommSize(mpcDiagPat);
     mpcDiagPat->finalize();
     paralApply(this->subdomains, &FetiSub<Scalar>::sendMpcDiag, mpcDiagPat);
     mpcDiagPat->exchange();
     paralApply(this->subdomains, &FetiSub<Scalar>::collectMpcDiag, mpcDiagPat);
     delete mpcDiagPat;
   }

   if(fetiInfo->c_normalize) normalizeC();

   if(fetiInfo->mpc_precno == FetiInfo::diagCCt) {
     // use W scaling for preconditioning mpcs, don't need to build & invert CC^t
     paralApply(this->subdomains, &FetiSub<Scalar>::sendMpcScaling, this->vPat);
     this->vPat->exchange();  // neighboring subs mpc weights
     paralApply(this->subdomains, &FetiSub<Scalar>::collectMpcScaling, this->vPat);
   }
   else if(fetiInfo->mpc_precno != FetiInfo::noMpcPrec) {
     // used generalized proconditioner for mpcs, need to build and invert CC^t
     mpcPrecon = true;
     paralApply(this->subdomains, &FetiSub<Scalar>::initMpcScaling);
   }
   if(mpcSubMap) { delete [] mpcSubMap; mpcSubMap = 0; }
 }

 paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::cleanMpcData);
#endif  // As_LIBRARY
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::checkStoppingCriteria(int iter, double error, double ff) const
{
  // 1. check if maximum number of iterations reached
  if(iter == this->maxiter) {
    this->times.iterations[this->numSystems].stagnated = 0;
    return true;
  }
 
  // 2. check for convergence
  if(sqrt(error) <= std::max(fetiInfo->tol*sqrt(ff), fetiInfo->absolute_tol)) {
    this->times.iterations[this->numSystems].stagnated = 0;
    return true;
  }

  // 3. check for stagnation
  if(iter > 0 && (std::fabs(sqrt(error)-sqrt(lastError)) < std::max(fetiInfo->stagnation_tol*sqrt(lastError), fetiInfo->absolute_stagnation_tol))) {
     this->times.iterations[this->numSystems].stagnated = 1;
     return true;
  }

  lastError = error;
  return false;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::zeroG()
{
  paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::zeroG);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::deleteG()
{
  paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::deleteG);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::addRstar_gT(int iGroup, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta) const
{
#ifdef DISTRIBUTED
	for(int i=0; i<this->nsub; ++i) {
    if(this->subdomains[i]->getGroup() == groups[iGroup])
      this->subdomains[i]->addRstar_gT(u.subData(this->subdomains[i]->localSubNum()), beta);
  }
#else
	auto grsubs = (*groupToSub)[iGroup];
	for(int i = 0; i < groupToSub->num(iGroup); ++i) {
		int iSub = grsubs[i];
		this->subdomains[iSub]->addRstar_gT(u.subData(this->subdomains[iSub]->localSubNum()), beta);
	}
#endif
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subtractRstar_g(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta) const
{
	this->subdomains[iSub]->subtractRstar_g(u.subData(this->subdomains[iSub]->localSubNum()), beta);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeProjectedDisplacement(GenDistrVector<Scalar> &u) const
{
	int numGtGsing = (GtGtilda) ? GtGtilda->numRBM() : 0;
	if(numGtGsing > 0) {
		// get null space of GtGtilda
		Scalar *zem = new Scalar[numGtGsing*ngrbms];
		GtGtilda->getNullSpace(zem);
		GenFullM<Scalar> X(zem, ngrbms, numGtGsing, 1);
		// build global RBMs
		if(geometricRbms || fetiInfo->corners == FetiInfo::noCorners)
			paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::buildGlobalRBMs, X, cornerToSub);
		else {
			// TODO: FetiSub<Scalar>::Rstar first needs to be filled from the nullspace of Kcc^*
			std::cerr << " *** WARNING: FetiDPSolver::computeProjectedDisplacement requires GRBM or \"corners none\"\n";
			return;
		}

		GenVector<Scalar> beta(numGtGsing, 0.0);
		// 1. compute beta = Rstar_g^t * u
		execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::addRstar_gT, u, beta);
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(numGtGsing, beta.data());
#endif
		// build RtR
		GenFullM<Scalar> RtR(numGtGsing, numGtGsing);
		RtR.zero();
		for(int i=0; i<this->nsub; ++i) this->subdomains[i]->assembleRtR(RtR); // can't be done in parallel
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(numGtGsing*numGtGsing, RtR.getData());
#endif
		RtR.factor();

		// 3. compute beta = (Rstar_g^t Rstar_g)^(-1) * beta
		RtR.reSolve(beta.getData());
		// 4. compute u = u - Rstar_g * beta
		execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::subtractRstar_g, u, beta);
	}
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::makeGtG()
{
	if(verboseFlag) filePrint(stderr, " ... Build G matrix and GtG solver  ...\n");
	startTimerMemory(this->times.coarse1, this->times.memoryGtG);

	// 0. delete previous G if it exists
	deleteG();

	// 1. make local G = Bbar * R
	if(geometricRbms || fetiInfo->corners == FetiInfo::noCorners)
		paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::makeG);
	else {
		paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::makeTrbmG, kccrbms, KccSolver->numRBM(), KccSolver->neqs());
		paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::getNumGroupRBM, ngrbmGr);
#ifdef DISTRIBUTED
		this->fetiCom->globalMax(nGroups, ngrbmGr);
#endif
	}

	// 2. exchange G's between neighboring subdomains
	FSCommPattern<int> *sPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
	for(int i = 0; i < this->nsub; ++i) this->subdomains[i]->setMpcNeighbCommSize(sPat, 2);
	sPat->finalize();
	paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::sendNeighbGrbmInfo, sPat);
	sPat->exchange();  // neighboring subs number of group RBMs and offset
	paralApply(this->nsub, this->subdomains.data(), &FetiBaseSub::receiveNeighbGrbmInfo, sPat);
	delete sPat;
	// create bodyRbmPat FSCommPattern object, used to send/receive G matrix
	FSCommPattern<Scalar> *gPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU,
	                                                        FSCommPattern<Scalar>::CopyOnSend, FSCommPattern<Scalar>::NonSym);
	for(int i = 0; i < this->nsub; ++i) this->subdomains[i]->setGCommSize(gPat);
	gPat->finalize();
	paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::sendG, gPat);
	gPat->exchange();
	paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::receiveG, gPat);
	delete gPat;

	// 3. build coarse connectivity and equation numberer
	if(this->mpcToSub) {
		Connectivity mpcToBody = this->mpcToSub->transcon(*subToGroup); // PJSA 6-15-06
		Connectivity bodyToMpc = mpcToBody.reverse();
		Connectivity bodyToBody_mpc = bodyToMpc.transcon(mpcToBody);
		coarseConnectGtG = bodyToBody_mpc.modify();
	}
	else {
		Connectivity *bodyToSub = subToBody->alloc_reverse();
		coarseConnectGtG = bodyToSub->transcon(*subToBody);
		delete bodyToSub;
	}
	eqNumsGtG = new SimpleNumberer(nGroups);
	for(int i = 0; i < nGroups; ++i) eqNumsGtG->setWeight(i, ngrbmGr[i]);
	eqNumsGtG->makeOffset();

	// 4. create, assemble and factorize GtG
	rebuildGtGtilda();

	// 5. check for singularities in GtGstar (representing global RBMs)
	//if(GtG->numRBM() > 0)
	//  filePrint(stderr, " ... GtG has %d singularities for tol %e ...\n", GtG->numRBM(), fetiInfo->grbm_tol);

	stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::makeE(GenDistrVector<Scalar> &f) const
{
	// Compute e = R^t * f
	GenVector<Scalar> &e = this->wksp->ret_e();
	e.zero();

	if(geometricRbms || fetiInfo->corners == FetiInfo::noCorners) {
		execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::assembleE, e, f);
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(ngrbms, e.data());
#endif
	}
	else {
		GenVector<Scalar> &fc = this->wksp->ret_fc();
		GenDistrVector<Scalar> &fr = this->wksp->ret_fr();
		for(int i = 0; i < this->nsub; ++i)
			this->subdomains[i]->assembleTrbmE(kccrbms, KccSolver->numRBM(), KccSolver->neqs(), e.data(), fr.subData(this->subdomains[i]->localSubNum()));
		for(int i = 0; i < KccSolver->numRBM(); ++i)
			for(int j = 0; j < KccSolver->neqs(); ++j)
				e[i] += fc[j]*kccrbms[i*KccSolver->neqs()+j];
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(ngrbms, e.data());
#endif
	}
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::assembleE(int iGroup, GenVector<Scalar> &e, GenDistrVector<Scalar> &f) const
{
#ifdef DISTRIBUTED
	for(int i = 0; i < this->nsub; ++i) {
    if(this->subdomains[i]->getGroup() == groups[iGroup])
      this->subdomains[i]->assembleE(e, f.subData(this->subdomains[i]->localSubNum()));
  }
#else
	auto grsubs = (*groupToSub)[iGroup];
	for(int i = 0; i < groupToSub->num(iGroup); ++i) {
		int iSub = grsubs[i];
		this->subdomains[iSub]->assembleE(e, f.subData(this->subdomains[iSub]->localSubNum()));
	}
#endif
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::addRalpha(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &alpha) const
{
	this->subdomains[iSub]->addRalpha(u.subData(this->subdomains[iSub]->localSubNum()), alpha);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getRBMs(Scalar *globRBM)
{
	bool useKccSolver = (this->glNumMpc == 0 && !geometricRbms);
	if(GtGtilda && !useKccSolver) {
		int nRBM = numRBM();
		int iRBM;
		for(iRBM = 0; iRBM < nRBM; ++iRBM) {
			GenStackDistVector<Scalar> R(this->internalDI, globRBM+iRBM*(this->internalDI.len));
			R.zero();
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::getGlobalRBM, iRBM, (GenDistrVector<Scalar> &)(R));
		}
	}
	else if(KccSolver) {
		int nc = KccSolver->neqs();
		int nr = numRBM();
		Scalar *R = new Scalar[nr*nc];
		if(nc > 0) KccSolver->getNullSpace(R);
		GenDistrVector<Scalar> vr(internalR);
		int iRBM;
		for(iRBM=0; iRBM<nr; ++iRBM) {
			GenStackDistVector<Scalar> v(this->internalDI, globRBM+iRBM*this->internalDI.len);
			GenStackVector<Scalar> vc(nc, R+iRBM*nc);
			vr.zero();
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, vr, (GenVector<Scalar> &)(vc));
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::KrrReSolve, vr);
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, vr, (GenVector<Scalar> &)(vc),
			          (GenDistrVector<Scalar> &)(v), (GenDistrVector<Scalar> &)(v));  // last argument is a dummy
		}
		delete [] R;
	}
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getRBMs(GenDistrVectorSet<Scalar> &globRBM)
{
	bool useKccSolver = (this->glNumMpc == 0 && !geometricRbms);
	if(GtGtilda && !useKccSolver) {
		int numGtGsing = GtGtilda->numRBM();
		if(numGtGsing > 0 && domain->probType() == SolverInfo::Modal) {
			// get null space of GtGtilda
			Scalar *zem = new Scalar[numGtGsing*ngrbms];
			GtGtilda->getNullSpace(zem);
			GenFullM<Scalar> X(zem, ngrbms, numGtGsing, 1);
			// build global RBMs
			if(geometricRbms || fetiInfo->corners == FetiInfo::noCorners)
				paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::buildGlobalRBMs, X, cornerToSub);
			else {
				// TODO: FetiSub<Scalar>::Rstar first needs to be filled from the nullspace of Kcc^*
				std::cerr << " *** WARNING: FetiDPSolver::getRBMs requires GRBM or \"corners none\"\n";
				globRBM.zero();
				return;
			}
		}

		int nRBM = numRBM();
		int iRBM;
		for(iRBM = 0; iRBM < nRBM; ++iRBM) {
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::getGlobalRBM, iRBM, globRBM[iRBM]);
		}
	}
	else if(KccSolver) {
		int nc = KccSolver->neqs();
		int nr = numRBM();
		Scalar *R = new Scalar[nr*nc];
		if(nc > 0) KccSolver->getNullSpace(R);
		GenDistrVector<Scalar> vr(internalR);
		int iRBM;
		for(iRBM=0; iRBM<nr; ++iRBM) {
			GenStackVector<Scalar> vc(nc, R+iRBM*nc);
			vr.zero();
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, vr, (GenVector<Scalar> &)(vc));
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::KrrReSolve, vr);
			execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, vr, (GenVector<Scalar> &)(vc),
			          globRBM[iRBM],globRBM[iRBM]); // last argument is a dummy
		}
		if(R) delete [] R;
	}
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R)
{
	Scalar *localRvec = R.subData(this->subdomains[iSub]->localSubNum());
	this->subdomains[iSub]->getGlobalRBM(iRBM, localRvec);
}

template<class Scalar>
inline void
GenFetiDPSolver<Scalar>::split(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &v_f,
                               GenDistrVector<Scalar> &v_c) const
{
	this->subdomains[iSub]->split(v.subData(this->subdomains[iSub]->localSubNum()), v_f.subData(this->subdomains[iSub]->localSubNum()),
	                      v_c.subData(this->subdomains[iSub]->localSubNum()));
}

template
class GenFetiDPSolver<double>;
template
class GenFetiDPSolver<std::complex<double>>;
