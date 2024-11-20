#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <stdexcept>
#include <sys/types.h>

#ifndef WINDOWS
#include <sys/mman.h>
#endif
#include <Driver.d/SubDomain.h>
#include <Feti.d/Feti.h>
#include <Feti.d/CGOrthoSet.h>
#include <Utils.d/linkfc.h>
#include <Threads.d/PHelper.h>
#include <Math.d/matrix.h>
#include <Math.d/SymFullMatrix.h>
#include <Math.d/BigMatrix.h>
#include <Math.d/VectorSet.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/GetTime.h>

#include <Feti.d/FetiOp.h>
#include <Feti.d/FetiOpControler.h>
#include <Feti.d/FetiInfo.h>
#include <Feti.d/CoarseSet.h>
#include <Feti.d/GMRESOrthoSet.h>
#include <Feti.d/GCROrthoSet.h>
#include <Corotational.d/DistrGeomState.h>
#include <Solvers.d/SolverFactory.h>
class GeomState;

// KHP: To change between old FETI-2 coarse grid and
// new sparse FETI-2 coarse grid with mpcs.
//#define ORIGINAL_FETI2

// when defined, use original parallel skyline solver
// when not defined, use new parallel skyline solver
//#define REGULAR_SKY

#ifdef BOOL_NOT_DEFINED
#define false 0
#endif

#include <Utils.d/Memory.h>

extern long totMemSky;
extern long totMemSparse;

#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, double *e, const int &f,
                  double *g, const int &h, const double &i, double *j, const int &k)
{
 _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}

inline void Tgemv(const char &a, const int &b, const int &c,
                  const complex<double> &d, complex<double> *e, const int &f,
                  complex<double> *g, const int &h, const complex<double> &i, complex<double> *j, const int &k)
{
 _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif

template<class Scalar>
GenFetiSolver<Scalar>::GenFetiSolver(int _nsub, GenSubDomain<Scalar> **_sd, Connectivity *_subToSub,
                                     FetiInfo *_fetiInfo, FSCommunicator *_fetiCom, int *glToLoc, Connectivity *_mpcToSub, Connectivity *_cpuToSub,
                                     GenSolver<Scalar> **sysMatrices, GenSparseMatrix<Scalar> **sysSparse, Rbm **_rbms, int _verboseFlag) :
	subdomains(_sd, _sd+_nsub), internalDI(_nsub), interface(_nsub), times(threadManager->numThr(), _nsub), verboseFlag(_verboseFlag)
{

	// Compute memory used by FETI Solver
	times.memoryFETI -= memoryUsed();

	fetiInfo = _fetiInfo;	// Feti solver information
	fetiCom  = _fetiCom;
	subToSub = _subToSub;	// subdomain to subdomain connectivity
	mpcToSub = (_mpcToSub) ? _mpcToSub : new Connectivity();  // multiple point constraint to subdomain connectivity
	cpuToSub = _cpuToSub;
	myCPU = 0; numCPUs = 1;

	if(_rbms) {
		isDynamic = 1; // multi-domain linear dynamics
		QGisLocal = 0; // Fi is a non local operator
	} else {
		isDynamic = 0; // multi-domain linear statics
		QGisLocal = !fetiInfo->nonLocalQ;
	}

	isFeti2 = fetiInfo->version; // if we are using FETI 2
	nsub    = _nsub;	      // Number of subdomains
	sd      = _sd;		      // pointer to Array of all Subdomains

	// Define FETI tolerance and maximum number of iterations
	double fetiTolerance = fetiInfo->tol;
	epsilon2 = fetiTolerance*fetiTolerance;
	maxiter  = fetiInfo->maxit;

	// Classes to organize parallel execution of tasks
	fetiTasks = new TaskDescr *[nsub];
	fetiOps   = new GenFetiOp<Scalar> *[nsub];
	opControl = new GenFetiOpControler<Scalar>(nsub, QGisLocal);
	opControl->nQ = fetiInfo->nQ;

	// create vPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
	vPat = new FSCommPattern<Scalar>(fetiCom, cpuToSub, 0, FSCommPattern<Scalar>::CopyOnSend);
	for(int iSub=0; iSub<nsub; ++iSub) subdomains[iSub]->setDofCommSize(vPat);
	vPat->finalize();

	// create sPat FSCommPattern object, used to send/receive a single integer eg. numRBMs & crnDofSize
	sPat = new FSCommPattern<int>(fetiCom, cpuToSub, 0, FSCommPattern<int>::CopyOnSend);
	for(int iSub=0; iSub<nsub; ++iSub) subdomains[iSub]->setCommSize(sPat, 1);
	sPat->finalize();

	int iSub;
	for(iSub = 0; iSub < nsub; ++iSub) {
		if(_rbms)
			fetiOps[iSub] = new GenFetiOp<Scalar>(sd[iSub], opControl, isFeti2, isDynamic, vPat,
			                                      _rbms[iSub]);
		else
			fetiOps[iSub] = new GenFetiOp<Scalar>(sd[iSub], opControl, isFeti2, isDynamic, vPat);
		fetiTasks[iSub] = (TaskDescr *) fetiOps[iSub];
	}

	// Compute total interface length, total internal length
	// and total half interface length
	int tInterfLen    = 0;
	int tLocalLen     = 0;
	int halfInterfLen = 0;
	for(iSub = 0; iSub < nsub; ++iSub) {
		interface.domLen[iSub] = subdomains[iSub]->interfLen();
		internalDI.domLen[iSub]  = subdomains[iSub]->getNumUncon();

		subdomains[iSub]->computeMasterFlag(*mpcToSub);
		fetiOps[iSub]->setHalfOffset(halfInterfLen);
		halfInterfLen += subdomains[iSub]->halfInterfLen();
		tInterfLen    += subdomains[iSub]->interfLen();
		tLocalLen     += subdomains[iSub]->localLen();
	}
	interface.len = tInterfLen;
	internalDI.len  = tLocalLen;
	halfSize      = halfInterfLen;

	// compute the masterFlags
	bool *interfaceMasterFlag = new bool[tInterfLen];
	interface.computeOffsets();
	for(iSub = 0; iSub < nsub; ++iSub) {
		auto &subMasterFlag = subdomains[iSub]->getMasterFlag();
		int subOffset = interface.subOffset[iSub];
		int j;
		for(j=0; j<interface.domLen[iSub]; ++j)
			interfaceMasterFlag[subOffset+j] = subMasterFlag[j];
	}
	interface.setMasterFlag(interfaceMasterFlag);
	internalDI.setMasterFlag();
	// don't delete interfaceMasterFlag

	// Allocate space for reorthogonalization set
	times.memoryOSet -= memoryUsed();
	oSetCG = (fetiInfo->maxortho > 0) ? new GenCGOrthoSet<Scalar>(halfInterfLen, fetiInfo->maxortho, fetiCom) : 0;
	times.memoryOSet += memoryUsed();

	if(sysMatrices != 0) {
		for(int iSub = 0; iSub < nsub; ++iSub)
			fetiOps[iSub]->setSysMatrix(sysMatrices[iSub], sysSparse[iSub]);
	}

	// Compute stiffness scaling if required
	if(fetiInfo->scaling == FetiInfo::kscaling) {
		execParal(nsub, this, &GenFetiSolver<Scalar>::sendScale);
		vPat->exchange();
		paralApply(nsub, sd, &GenSubDomain<Scalar>::collectScaling, vPat);
	}

	// Factor matrices: K and Kii (if dirichlet preconditioner))
	fprintf(stderr," ... Factor Subdomain Matrices      ... \n");
	double memRBM =-memoryUsed();
	startTimerMemory(times.factor, times.memoryFactor);
	timedParal(times.factorMat, nsub, this, &GenFetiSolver<Scalar>::factorMatrices);
	stopTimerMemory(times.factor, times.memoryFactor);
	memRBM += memoryUsed();

	int glNumMpc = (mpcToSub) ? mpcToSub->csize() : 0;
#ifndef DISTRIBUTED
	if(fetiInfo->feti2version == FetiInfo::sparseCoarse) {

		if((isFeti2 && fetiInfo->type == FetiInfo::linear) || (glNumMpc > 0))
			makeSingleCoarse();
		else
			makeGtG();

		times.memoryDV -= memoryUsed();
		wksp = new GenFetiWorkSpace<Scalar>(interface, internalDI, fetiInfo->type,
		                                    2*numrbms+crns+glNumMpc, crns+glNumMpc);
		times.memoryDV += memoryUsed();
	}
	else {
		makeGtG(); // Make first level coarse problem

		// Make second level coarse problem
		if (isFeti2 && (isDynamic == 0)) makePCtFPC();

		// Allocate Distributed Vectors necessary for FETI solve loop
		times.memoryDV -= memoryUsed();
		wksp = new GenFetiWorkSpace<Scalar>(interface, internalDI, fetiInfo->type,
		                                    numrbms, crns);
		times.memoryDV += memoryUsed();
	}
#else
	if ((isFeti2 && (isDynamic == 0)) || glNumMpc > 0) {
   makeSingleCoarse();
   times.memoryDV -= threadManager->memoryUsed();
   wksp = new GenFetiWorkSpace<Scalar>(interface, internalDI, fetiInfo->type,
                                       2*numrbms+crns+glNumMpc, crns+glNumMpc);
   times.memoryDV += threadManager->memoryUsed();
 }
 else {
   makeDistGtG(glToLoc);
   times.memoryDV -= threadManager->memoryUsed();
   wksp = new GenFetiWorkSpace<Scalar>(interface, internalDI, fetiInfo->type,
                                       numrbms, crns);
   times.memoryDV += threadManager->memoryUsed();
 }
#endif
}

template <class Scalar>
GenFetiSolver<Scalar>::GenFetiSolver(int _nsub, GenSubDomain<Scalar> **subs, int _numThreads, int _verboseFlag)
	: subdomains(subs, subs+_nsub), internalDI(_nsub),
	  interface(_nsub), times(_numThreads,_nsub), verboseFlag(_verboseFlag) { }

template<class Scalar>
double
GenFetiSolver<Scalar>::getFNormSq(GenDistrVector<Scalar> &f)
{
	distributeForce(f);
	return f.sqNorm();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::sendScale(int iSub)
{
	subdomains[iSub]->sendDiag(fetiOps[iSub]->KasSparse, vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::distributeForce(GenDistrVector<Scalar> &f) const
{
	execParal(nsub, this, &GenFetiSolver<Scalar>::fSend,  f);
	vPat->exchange();
	execParal(nsub, this, &GenFetiSolver<Scalar>::fScale, f);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::fSend(int iSub, GenDistrVector<Scalar> &force) const
{
	subdomains[iSub]->fSend(force.subData(subdomains[iSub]->localSubNum()), vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::fScale(int iSub, GenDistrVector<Scalar> &force) const
{
	subdomains[iSub]->fScale(force.subData(subdomains[iSub]->localSubNum()), vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::distributeForce(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw) const
{
	execParal(nsub, this, &GenFetiSolver<Scalar>::fSendCoupled,  f, fw);
	vPat->exchange();
	execParal(nsub, this, &GenFetiSolver<Scalar>::fScaleCoupled, f, fw);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::fSendCoupled(int iSub, GenDistrVector<Scalar> &force,
                                    GenDistrVector<Scalar> &fw) const
{
	sd[iSub]->fSend(force.subData(subdomains[iSub]->localSubNum()), vPat,
	                fw.subData(subdomains[iSub]->localSubNum()));
}

template<class Scalar>
void
GenFetiSolver<Scalar>::fScaleCoupled(int iSub, GenDistrVector<Scalar> &force,
                                     GenDistrVector<Scalar> &fw) const
{
	sd[iSub]->fScale(force.subData(subdomains[iSub]->localSubNum()), vPat,
	                 fw.subData(subdomains[iSub]->localSubNum()));
}

template<class Scalar>
void
GenFetiSolver<Scalar>::fSplit(int iSub, GenDistrVector<Scalar> &force) const
{
	sd[iSub]->splitInterf(force.subData(subdomains[iSub]->localSubNum()));
}

/** Forming and factorization of the projected coarse grid FETI
// matrix : (P C)^t  F_I (P C)
//
// since (P C) = C - G (G^t G)^-1 G^t C 
//             = C + G  Rgc
// where Rgc = = -(G^t G)^-1 G^t C
//
// the matrix becomes (C + G  Rgc)^t F (C + G  Rgc).
//  = C'F C + Rgc'G'FC + C'F G Rgc + Rgc'G'F G Rgc */

template<class Scalar>
void
GenFetiSolver<Scalar>::makePCtFPC()
{
	startTimerMemory(times.coarse2, times.memoryPCtFPC);

	int i;
	// Initialize numCRNs, crnDofLen and BClocal
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::initializeCRNs, sPat);

	// Get total number of corner modes from subdomains
	crns = collectIntGlobalSum();
	times.numCRNs = crns;

	// Initialize neighbors number of corner modes
	sPat->exchange();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbCRNs, sPat);

	// Build Rgc = -(G^t G)^-1 G^t C
	int sizeRgc = numrbms*crns;
	Scalar *Rgc = new Scalar[sizeRgc];

	// cset is constructed in GenFetiSolver<Scalar>::factorMatrices
	GenCoarseSet<Scalar> *cset = opControl->cset;
	paralApplyToAll(nsub, cset, &GenCoarseSet<Scalar>::buildBNeighbCs, cset);

	opControl->Rgc  = Rgc;
	opControl->crns = crns;
	opControl->rbms = numrbms;

	// Setup SimpleNumberer for Corners, before any assembly associated with C
	Connectivity *coarseConnect2 = subToSub->transcon(subToSub);
	compStruct renumber2 = coarseConnect2->renumByComponent(0);
	PFcNums = new SimpleNumberer(nsub, renumber2.renum);

	for(i = 0; i<nsub; ++i)
		PFcNums->setWeight(i, fetiOps[i]->getcrnDofSize() );

	PFcNums->makeOffset();
	opControl->PFcNums = PFcNums;

	paralApplyToAll(nsub,fetiOps, &GenFetiOp<Scalar>::setBetaOffsets, PFcNums->allOffsets());

	if(crns > 0 && numrbms > 0) {

		// Compute G'C
		paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::assembleGtCs);

		// forward & backward substution to get -Rgc
		execParal(threadManager->numThr(), this, &GenFetiSolver<Scalar>::finishRgc,
		          threadManager->numThr());

		// multiply by -1.0 to get Rgc
		for(i=0; i<sizeRgc; ++i)
			opControl->Rgc[i] *= -1.0;

		// allocate memory for C'F C 
		times.memoryPCtFPCmat -= memoryUsed();
		PCtFPC = new GenBigMatrix<Scalar>(crns);
		times.memoryPCtFPCmat += memoryUsed();

		GenFullM<Scalar> *CtFCs     = PCtFPC;
		opControl->CtFCs = CtFCs;

		// allocate memory for GtFC and zero it
		GtFCs = new GenFullM<Scalar>(numrbms,crns);
		GtFCs->zero();
		opControl->GtFC = GtFCs;

		// Compute C'FC, G'FC and G'FG
		makelocalFcoarse();

		// Compute new projected F coarse
		GenStackFullM<Scalar> wkbf(crns, crns, PCtFPC->data() );

		addAllFcoarse(wkbf);

		times.pfactor2 -= getTime();
		PCtFPC->parallelFactor();
		times.pfactor2 += getTime();

	} else {
		PCtFPC = 0;
	}

	stopTimerMemory(times.coarse2, times.memoryPCtFPC);

}

template<class Scalar>
void
GenFetiSolver<Scalar>::finishRgc(int j, int totThr)
{
	int i, myStart, myEnd;
	int modulo = crns%totThr;

	myStart = j*(crns/totThr);
	myEnd = (j+1)*(crns/totThr);
	if(j < modulo) {
		myStart += j;
		myEnd +=j+1;
	} else  {
		myStart += modulo;
		myEnd += modulo;
	}
	Scalar **RHSs = (Scalar **) dbg_alloca(sizeof(Scalar *)*(myEnd-myStart));
	for (i = myStart; i < myEnd; i++) {
		int offset = i* numrbms;
		RHSs[i-myStart] = opControl->Rgc + offset;
//  single rhs version of reSolve
//  if(GtGsolver) GtGsolver->reSolve(opControl->Rgc + offset);
	}

	// Multiple rhs version of reSolve
	if(GtGsolver) GtGsolver->reSolve(myEnd-myStart,RHSs);
}

// [C G]^t F [C G]
// Very similar to  G^t Q G

template<class Scalar>
void
GenFetiSolver<Scalar>::makelocalFcoarse()
{
	// compute F C 

	// compute my own Fi^{-1} BCi
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::computeFiBC);

	// compute neighbors F C
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbFC);

	// compute G'FG, C'FC, C'FG
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::assembleGtFCs);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::makeGtG()
{
	startTimerMemory(times.coarse1, times.memoryGtG);

	if(isDynamic == 1 && fetiInfo->noCoarse == 1) {
		fprintf(stderr," ... No Coarse Problem used         ...\n");
		GtGsolver = 0;
		times.coarse1 += getTime();
		return;
	}

	// Get number of each neighbors rbms
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::sendNumNeighbRBM, sPat);
	sPat->exchange();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbRBM, sPat);

	// Get all the numbers of rigid body modes and dispatch RBMs to neighbors
	makeRbmPat();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::sendInterfRBM, rbmPat);
	rbmPat->exchange();

	// compute neighbors QGs 
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbQGs, rbmPat);

	if(isFeti2 && isDynamic) {

		// Initialize numCRNs, crnDofLen and BClocal
		paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::initializeCRNs, sPat);

		// compute total number of corners
		crns = collectIntGlobalSum();

		if(crns > 0) {

			sPat->exchange();
			// Initialize neighbors number of corners
			paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbCRNs, sPat);

			// cset is constructed in GenFetiSolver<Scalar>::factorMatrices
			GenCoarseSet<Scalar> *cset = opControl->cset;
			paralApplyToAll(nsub,cset,&GenCoarseSet<Scalar>::buildBNeighbCs, cset);

			// compute my own Fi Ci where Fi is this subdomains contribution
			// to the interface operator and Ci is this subdomains Corners
			paralApplyToAll(nsub,fetiOps,&GenFetiOp<Scalar>::computeFiBC);

			// compute neighbors Fj Ci 
			paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbFC);
		}

	}

	// Construct Coarse Problem Connectivity
	if(coarseConnect == 0) {
		if( QGisLocal )
			coarseConnect = subToSub;
		else
			coarseConnect = subToSub->transcon(subToSub);
	}

	// Maybe this could be an option too!
	//FILE *gtgFile = fopen("gtgConnect","w");
	//coarseConnect->print(gtgFile);

	// Construct 1st level Coarse Problem (GtG) renumbering
	// 0 = no renumbering
	// 1 = sloan renumbering
	// 2 = rcm renumbering
	if(renum == 0) {
		double t1 = getTime();
		long m1 = memoryUsed();
		renumber = coarseConnect->renumByComponent(1);
		if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
			fprintf(stderr," ... Renumber GtG %10.5f s %10.3f Mb\n",
			        (getTime() - t1)/1000.0,(memoryUsed()-m1)/(1024.0*1024.0));
	}

	// Construct Coarse Problems equation numbers
	if(eqNums == 0) eqNums = new SimpleNumberer(nsub, renumber.renum);

	int iSub;
	if(isFeti2 && isDynamic && crns > 0) {
		for(iSub = 0; iSub < nsub; ++iSub)
			eqNums->setWeight(iSub, fetiOps[iSub]->getNumRBM()+
			                        fetiOps[iSub]->getcrnDofSize());
	} else
		for(iSub = 0; iSub < nsub; ++iSub)
			eqNums->setWeight(iSub, fetiOps[iSub]->getNumRBM());

	eqNums->makeOffset();

	numrbms = eqNums->size();

	times.numRBMs = numrbms; // store size of GtG for timing file

	// The following sets the alpha offset for the subdomains AND in the
	// case of Feti2 dynamic it does betaOffset as well
	paralApplyToAll(nsub, fetiOps,&GenFetiOp<Scalar>::setAlphaOffsets, eqNums->allOffsets());

	// Now assemble GtG

	if(numrbms > 0 ) {

		opControl->eqNums = eqNums;
		if ( isFeti2 && isDynamic == 0 ) {
			GtQGs = new GenSymFullMatrix<Scalar>(numrbms);
			GtQGs->zero();
			opControl->GtQGs = GtQGs;
		}

		double t1 = getTime();
		times.memoryGtGsky -= memoryUsed();
		GtGsolver = GenSolverFactory<Scalar>::getFactory()->createDistSolver(coarseConnect, eqNums, *fetiInfo->coarse_cntl, opControl->sparseGtG, this->fetiCom);
		times.memoryGtGsky += memoryUsed();

		if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
			fprintf(stderr," ... Form GtG     %10.5f s %10.3f Mb\n",
			        (getTime()-t1)/1000.0,times.memoryGtGsky/(1024.0*1024.0));

		// Assemble GtQG
		t1 = getTime();
		long m1 = memoryUsed();
		paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::assembleGtQGs);
		if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
			fprintf(stderr," ... Assemble GtG %10.5f s %10.3f Mb\n",
			        (getTime()-t1)/1000.0,(memoryUsed() - m1)/(1024.0*1024.0));

		startTimerMemory(times.pfactor, times.memoryGtGsky);
		GtGsolver->parallelFactor();
		stopTimerMemory(times.pfactor, times.memoryGtGsky);

		glNumRBM = GtGsolver->numRBM();
		if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
			fprintf(stderr," ... Factor GtG   %10.5f s %10.3f Mb\n",
			        (times.pfactor)/1000.0,(times.memoryGtGsky)/(1024.0*1024.0));
		filePrint(stderr," ... Number of Global RBMs %2d       ...\n",
		          glNumRBM );

	} else {
		GtGsolver = 0; // No GtG Solver is used
	}

	stopTimerMemory(times.coarse1, times.memoryGtG);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::orthoAddCG(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp) const
{
	if(fetiInfo->maxortho <= 0) return;
#ifdef DISTRIBUTED
	dbg_alloca(0); // needed to clear the stack memory
#endif

	Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hFp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::gatherHalfInterface, &p, &Fp, hp, hFp);

/*
 // debug: check dot product
 Scalar hpFp = 0.0;
 for(int i=0; i<halfSize; ++i)
   hpFp += hp[i]*hFp[i];
 std::cerr << "hpFp = " << hpFp << ", pFp = " << pFp << std::endl;
// debug: check scatter/gather half interface
 GenDistrVector<Scalar> p_copy(p); p_copy.zero();
 timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::scatterHalfInterface, hp, &p_copy);
 vPat->exchange();
 timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::rebuildInterface, p_copy);
 std::cerr << "p - rebuilt p: ";
 for(int i=0; i<p.size(); ++i) std::cerr << p.data()[i] - p_copy.data()[i] << " "; std::cerr << std::endl;
 GenDistrVector<Scalar> Fp_copy(Fp); Fp_copy.zero();
 timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::scatterHalfInterface, hFp, &Fp_copy);
 vPat->exchange();
 timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::rebuildInterface, Fp_copy);
 std::cerr << "Fp - rebuilt Fp: ";
 for(int i=0; i<Fp.size(); ++i) std::cerr << Fp.data()[i] - Fp_copy.data()[i] << " "; std::cerr << std::endl;
*/

	times.memoryOSet -= memoryUsed();
	oSetCG->orthoAddTimed(times.orthogonalize, hp, hFp, pFp);
	times.memoryOSet += memoryUsed();
}

// basic FETI projector: P = I - G (G^tG)^-1 G^t
template<class Scalar>
void
GenFetiSolver<Scalar>::project(GenDistrVector<Scalar> &r, GenVector<Scalar> &alpha,
                               GenDistrVector<Scalar> &pr,int isDirect) const
{
	if(isDynamic && isDirect) {
		tProject(r, alpha, pr, 0);
		return;
	}
	double initTime = getTime();
	startTimerMemory(times.project, times.memoryProject1);

	if(GtGsolver == 0) {
		pr = r;
		alpha.zero();
		return;
	}

	if(QGisLocal)
		pr = r;
	else
		pr.zero();

#ifdef DISTRIBUTED
	alpha.zero();
#endif

	opControl->vec2 = alpha.data();
	opControl->dv1  = &r;

	if(isFeti2 && crns > 0 && isDynamic) {
		opControl->operation = &GenFetiOp<Scalar>::getCtMult;
		threadManager->execTimedParal(times.projection, nsub, fetiTasks);
	}

	// alpha = G^t r
	opControl->operation = &GenFetiOp<Scalar>::getGtMult;
	threadManager->execTimedParal(times.projection, nsub, fetiTasks);

#ifdef DISTRIBUTED
	if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif

	// alpha = (G^tG)^-1 alpha (reSolve overwrites input vector)
	times.forBack -= getTime();
	if(GtGsolver) GtGsolver->reSolve(alpha.data());
	times.forBack += getTime();

	opControl->dv2 = &pr;
	opControl->dv1 = &pr;

	opControl->operation = &GenFetiOp<Scalar>::subAlphaGQ;
	threadManager->execTimedParal(times.projection, nsub, fetiTasks);

	if(isFeti2 && crns > 0 && isDynamic) {
		opControl->operation = &GenFetiOp<Scalar>::subNuFC;
		threadManager->execTimedParal(times.projection, nsub, fetiTasks);
	}
	vPat->exchange();

	times.projection.addOverAll(0, getTime() - initTime);

	if(QGisLocal == 0 || (isFeti2 && crns > 0 && isDynamic)) {
		timedParal(times.projection, nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, pr);
	}

	if((QGisLocal == 0) || (isFeti2 && isDynamic)) {
		initTime = getTime();
		pr.linAdd(r);
		times.projection.addOverAll(0, getTime() - initTime );
	}

	stopTimerMemory(times.project, times.memoryProject1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::tProject(GenDistrVector<Scalar> &r, GenVector<Scalar> &alpha,
                                GenDistrVector<Scalar> &pr, int isDirect) const
{
	if(isDynamic && isDirect) {
		project(r,alpha,pr, 0);
		return;
	}
	double initTime = getTime();
	startTimerMemory(times.project, times.memoryProject1);

	pr = r;
	if(GtGsolver == 0) {
		alpha.zero();
		return;
	}

#ifdef DISTRIBUTED
	alpha.zero();
#endif

	opControl->vec2 = alpha.data();
	opControl->dv1 = &pr;

	times.projection.addOverAll(0, getTime() - initTime );

	times.fgt -= getTime();
#ifdef DISTRIBUTED
	timedParal(times.projection, nsub, this, &GenFetiSolver<Scalar>::getGtQMult, 
            alpha.data(), &pr);
#else
	initTime = getTime();
	// alpha = G^t Q pr
	opControl->operation = &GenFetiOp<Scalar>::getGtQMult;
	threadManager->execTimedParal(times.projection, nsub, fetiTasks);
	times.projection.addOverAll(0, getTime() - initTime);
#endif
	times.fgt += getTime();

	initTime = getTime();

	if(isFeti2 && crns > 0 && isDynamic) {
		opControl->operation = &GenFetiOp<Scalar>::getCtFMult;
		threadManager->execTimedParal(times.projection, nsub,fetiTasks);
	}

#ifdef DISTRIBUTED
	if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif

	// compute: alpha = (G^tG)^-1 alpha
	times.forBack -= getTime();
	if(GtGsolver) GtGsolver->reSolve(alpha.data());
	times.forBack += getTime();

	times.subAlphaG -= getTime();
	opControl->operation = &GenFetiOp<Scalar>::subAlphaG;
	threadManager->execTimedParal(times.projection, nsub, fetiTasks);
	times.subAlphaG += getTime();

	if(isFeti2 && crns > 0 && isDynamic) {
		opControl->operation = &GenFetiOp<Scalar>::subNuC;
		threadManager->execTimedParal(times.projection, nsub,fetiTasks);
	}

	stopTimerMemory(times.project, times.memoryProject1);
	times.projection.addOverAll( 0, getTime() - initTime );
}

template<class Scalar>
void
GenFetiSolver<Scalar>::computeDynamL0(GenDistrVector<Scalar> &f, GenVector<Scalar> &alpha,
                                      GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &u) const
{
	lambda.zero();

	if(GtGsolver) {
		// recover dtilde in lambda
		localSolveAndJump(f, lambda, u, lambda);

		// Now multiply by G^t
		opControl->vec1 = alpha.data();
		opControl->vec2 = alpha.data();
		opControl->dv1  = &lambda;
		opControl->dv2  = &lambda;

#ifdef DISTRIBUTED
		alpha.zero();
#endif

		if( isFeti2 && crns > 0) {
			opControl->operation = &GenFetiOp<Scalar>::getCtMult;
			threadManager->execParal(nsub,fetiTasks);
		}

		// alpha = G^t lambda 
		opControl->operation = &GenFetiOp<Scalar>::getGtMult;
		threadManager->execParal(nsub,fetiTasks);

#ifdef DISTRIBUTED
		if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif

		// alpha = (G^tG)^-1 alpha
		times.forBack -= getTime();
		GtGsolver->reSolve(alpha.data());
		times.forBack += getTime();

		lambda.zero();

		// lambda = G alpha
		opControl->operation = &GenFetiOp<Scalar>::subAlphaG;
		threadManager->execParal(nsub,fetiTasks);

		if( isFeti2 && crns > 0) {
			opControl->operation = &GenFetiOp<Scalar>::subNuC;
			threadManager->execParal(nsub,fetiTasks);
		}
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::computeL0(GenDistrVector<Scalar> &f, GenVector<Scalar> &alpha,
                                 GenDistrVector<Scalar> &lambda) const
{
	// compute initial lambda to begin FETI iterations

	if(GtGsolver) {
		// Compute alpha_i = e = R_i^t f_i (where e is standard FETI notation)
		opControl->dv1  = &f;
		opControl->vec2 = alpha.data();

#ifdef DISTRIBUTED
		alpha.zero();
#endif

		// Compute alpha_i = R_i^t f_i  where i=1,...Nsub
		execParal(nsub, this, &GenFetiSolver<Scalar>::getRMult, &f, alpha.data());

#ifdef DISTRIBUTED
		if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
//filePrint(stderr,"Pre Alpha %20.14e %p %d\n",
//          alpha*alpha,fetiCom,alpha.size());
#endif

		// alpha = (G^tG)^-1 alpha (reSolve overwrites input vector)
		//fprintf(stderr," Pre Alpha'^2 %e\n",alpha*alpha);
		GtGsolver->reSolve(alpha.data());
		//alpha.print("Solution");
		//fprintf(stderr," Post Alpha'^2 %e\n",alpha*alpha);

		lambda.zero();
		opControl->dv1 = &lambda;

		// lambda -= Q G alpha
		opControl->operation = &GenFetiOp<Scalar>::subAlphaGQ;
		threadManager->execParal(nsub,fetiTasks);

		if(QGisLocal == 0) {
			vPat->exchange();
			execParal(nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, lambda);
		}
	} else
		lambda.zero();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addR(GenDistrVector<Scalar> &v, GenVector<Scalar> &alpha) const
{
	times.addR -= getTime();
	execParal(nsub, this, &GenFetiSolver<Scalar>::addRP, &v, alpha.data());
	times.addR += getTime();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::orthogonalize(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &p) const
{
	times.reOrtho -= getTime();

	if(fetiInfo->maxortho <= 0 || oSetCG->numDir() == 0) {
		p = r;
	}
	else {
		Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hOp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

		GenDistrVector<Scalar> *vec2 = 0;
		timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::gatherHalfInterface,
		           &r, vec2, hp, hOp);
/*
// debug: check scatter/gather half interface
    GenDistrVector<Scalar> r_copy(r); r_copy.zero();
    timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::scatterHalfInterface, hp, &r_copy);
    vPat->exchange();
    timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::rebuildInterface, r_copy);
    std::cerr << "r - rebuilt r: ";
    for(int i=0; i<r.size(); ++i) std::cerr << r.data()[i] - r_copy.data()[i] << " "; std::cerr << std::endl;
*/
		oSetCG->orthogonalizeTimed(times.orthogonalize, hp, hOp, fetiInfo->complex_hermitian);

		timedParal(times.orthogonalize, nsub, this,
		           &GenFetiSolver<Scalar>::scatterHalfInterface, hOp, &p);
		vPat->exchange();
		timedParal(times.orthogonalize, nsub, this,
		             &GenFetiSolver<Scalar>::rebuildInterface, p);
	}

	times.reOrtho += getTime();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::reSolve(GenDistrVector<Scalar> &u)
{
	// We are allocating memory here, this is bad!
	GenDistrVector<Scalar> f(internalDI);

	f = u;
	solve(f,u);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::resetOrthoSet()
{
	if(fetiInfo->maxortho <= 0) return;
	switch(fetiInfo->outerloop) {
		default:
		case 0:
		case 3:
			if(fetiInfo->nlPrecFlg) oSetCG->newOrthoSet();
			else oSetCG->reset();
			break;
		case 1:
			oSetGMRES->reInit();
			break;
		case 2:
			oSetGCR->reset();
			break;
	}
}

template<class Scalar>
bool
GenFetiSolver<Scalar>::predict(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0) const
{
	// KHP: NOTE, check this for when we want to rebuild the tangent
	//      stiffness matrix at every other nonlinear iteration. This may
	//      have to be modified. Double check it. When we do not rebuild,
	//      we need to call predict, when we do rebuild, we do not call
	//      predict. To allow for nonlinear problems to also use
	//      multiple rhs technique.

	if(fetiInfo->maxortho <= 0 || oSetCG->numDir() == 0) { // No available prediction
		return false;
	}
	else if(!fetiInfo->useMRHS) {
		oSetCG->reset();
		return false;
	}
	else {
		if(verboseFlag) filePrint(stderr," Current orthoset contains %d Krylov vectors\n",oSetCG->numDir());

		Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hOp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

		execParal(nsub, this, &GenFetiSolver<Scalar>::gatherHalfInterface, &r, &lambda0, hp, hOp);
		oSetCG->predict(hp, hOp);
		execParal(nsub, this, &GenFetiSolver<Scalar>::scatterHalfInterface, hOp, &lambda0);
		vPat->exchange();
		execParal(nsub, this, &GenFetiSolver<Scalar>::rebuildInterface, lambda0);

		return true; // Prediction available
	}
}


template<class Scalar>
void
GenFetiSolver<Scalar>::setAndStoreInfo(int iter, double finalPrimal2,
                                       double finalDual2 ) const
{
	times.iterations[numSystems].numFetiIter = iter;
	times.iterations[numSystems].finalPrimal = (finalPrimal2 == 0.0) ? 0.0 : sqrt(finalPrimal2);
	times.iterations[numSystems].finalDual   = (finalDual2 == 0.0) ? 0.0 : sqrt(finalDual2);

	// increment number of systems solved
	numSystems++;
	times.numSystems = numSystems;

	// sum total number of FETI iterations 
	times.numIter    += iter;
}

// basic FETI projector
// P = I - G (G^tG)^-1 Gt

template<class Scalar>
void
GenFetiSolver<Scalar>::solve(const GenDistrVector<Scalar> &_f, GenDistrVector<Scalar> &u)
{
	// We copy the right hand side because we're going to modify it.
	GenDistrVector<Scalar> f(_f);
	// Re-distribute applied forces if we are using k-scaling
	if(fetiInfo->scaling == FetiInfo::kscaling) {
		distributeForce(f);
	}

	if(fetiInfo->feti2version == FetiInfo::sparseCoarse) {
		if((isFeti2 && fetiInfo->type == FetiInfo::linear) ||
		   (mpcToSub != 0 && mpcToSub->csize() > 0 ))
		{
			singleCoarseSolve(f,u);
			return;
		}
	}

	// begin FETI solution timer
	times.solve -= getTime();

	// Assign working vectors
	GenDistrVector<Scalar> &r       = wksp->ret_r();	// residual
	GenDistrVector<Scalar> &lambda0 = wksp->ret_lambda();	// Lagrange multipliers
	GenDistrVector<Scalar> &w       = wksp->ret_w();	// projected residual
	GenDistrVector<Scalar> &y       = wksp->ret_y();
	GenDistrVector<Scalar> &z       = wksp->ret_z();	// preconditioned residual
	GenDistrVector<Scalar> &p       = wksp->ret_p();
	GenDistrVector<Scalar> &pr      = wksp->ret_pr();
	GenDistrVector<Scalar> &Fp      = wksp->ret_Fp();
	GenDistrVector<Scalar> &du      = wksp->ret_du();
	GenDistrVector<Scalar> &uzero   = wksp->ret_uzero();
	uzero.zero();

	GenDistrVector<Scalar> &deltaU  = wksp->ret_deltaU();
	deltaU.zero();
	GenDistrVector<Scalar> &deltaF  = wksp->ret_deltaF();
	deltaF.zero();

	// Added these GenDistrVector<Scalar>s for nonlinear FETI
	GenDistrVector<Scalar> &zz       = wksp->ret_zz();
	GenDistrVector<Scalar> &uu       = wksp->ret_uu();
	GenDistrVector<Scalar> &rCompare = wksp->ret_rCompare();

	GenVector<Scalar> &alpha = wksp->ret_alpha();
	GenVector<Scalar> &beta  = wksp->ret_beta();
	GenVector<Scalar> &gamma = wksp->ret_gamma();

	// Compute the square of a pseudo-norm of f
	double pseudoFNormSq = f.sqNorm();
	//fprintf(stderr, "pseudoFNormSq: %e\n", pseudoFNormSq);
	if(pseudoFNormSq == 0.0) {
		u.zero();
		return;
	}

	// Compute lambda0 such that f + B^t lambda0 is orthogonal 
	// to the RBMs or compatible for dynamics
	if(isDynamic)
		computeDynamL0(f, alpha, lambda0, u);
	else
		computeL0(f, alpha, lambda0); // lambda0 = G(G^tG)^-1 e

	// Solve u = K^+ (f + B^t lambda0)
	//   and r = B u
	localSolveAndJump(f,lambda0,u,r);

	// DEBUG
	if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10) {
		filePrint(stderr, " ... Lambda0*Lambda0 %10.4e     ...\n", lambda0.sqNorm());
		filePrint(stderr, " ... f*f             %10.4e     ...\n", pseudoFNormSq);
		filePrint(stderr, " ... r*r             %10.4e     ...\n", r.sqNorm());
		filePrint(stderr, " ... u*u             %10.4e     ...\n", u.sqNorm());
	}

	// For multi-solve cases (dynamics & nonlinear problems where
	// tangent stiffness matrix is not rebuilt)
	bool hasPred = predict(r, lambda0);
	if(hasPred) {
		localSolveAndJump(f, lambda0, u, r);
	}

	// Project: w = P^t r
	tProject(r, alpha, w);

	double w0Norm2 = w.sqNorm(); // Save initial Dual error

	if(isFeti2  && crns > 0 && isDynamic == 0) {
		updateFeti2lambda(w,r,lambda0,gamma,alpha);
		localSolveAndJump(f,lambda0,u,r);
		tProject(r, alpha, w);
	}

	// Precondition: z = F^-1 w 
	double errorEstimator = preCondition(w, z);

	if(fetiInfo->primalFlag) outputPrimalResidual(0, deltaF);

	// New Code added by ML for nonlinear Krylov preconditioner
	if(fetiInfo->nlPrecFlg && numSystems != 0) {

		// zz = Pz where P = I - G (G^tG)^-1 G^t
		project(z,beta,zz);

		// Solve uu = K^+ (uzero + B^t zz) and rCompare = B uu
		localSolveAndJump(uzero,zz,uu,rCompare);

		tProject(rCompare, beta, pr);
		pr -= w;

		double error1 = rCompare.sqNorm();
		fprintf(stderr," Error 1: %e w*w = %e\n", error1, w.sqNorm());

		int hasPrecond = nlPreCondition(w,z); // Nonlinear Preconditioner

		if(hasPrecond) {
			project(z,beta,zz);

			// Solve uu = K^+ (uzero + B^t zz) and rCompare = B uu
			localSolveAndJump(uzero,zz,uu,rCompare);

			double error2 = rCompare.sqNorm();
			fprintf(stderr," Error 2: %e w*w = %e\n", error2, w.sqNorm());

			tProject(rCompare, beta, pr);
			pr -= w;
		}

	}

	// Check convergence before entering main loop
	if( errorEstimator < epsilon2 * pseudoFNormSq ) {
		//if( z.sqNorm() < epsilon2 * pseudoFNormSq ) {

		// add contribution of rigid body modes: u = u + alpha*R
		if(isDynamic == 0) addR( u, alpha );

		// make solution compatible u = u + deltaU
		makeCompatible(u, deltaU);

		// Store number of iterations, finalPrimal and finalDual errors
		setAndStoreInfo(0, (errorEstimator/pseudoFNormSq), (w.sqNorm())/w0Norm2 );

		times.iterations[numSystems].stagnated = 0;

		times.solve += getTime();

		if(numSystems == 1) {
			times.memoryFETI  += memoryUsed();
		}

		return;
	}
	if((fetiInfo->numPrint() > 0) && (fetiInfo->numPrint() < 10)) {
		//filePrint(stderr," %3d Relative Primal Error = %e "
		//                 " Relative Dual Error = %e\n",
		//    0, sqrt(errorEstimator/pseudoFNormSq), 1.0);
		filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");
		filePrint(stderr," %4d %23.6e %21.6e\n",0, sqrt(errorEstimator/pseudoFNormSq), 1.0);
	}

	// Re-project: y = P z
	project(z, alpha, y);

	orthogonalize(y,p);

	// if Feti2, update y
	if (isFeti2  && crns > 0 && isDynamic == 0) updateFeti2y(p,w,gamma,alpha);

	// Iterate
	int iter;
	double lastz2 = 0.0, z2;
	for(iter = 0; iter < maxiter; ++iter) {

		Scalar pFp = localSolveAndJump(uzero, p, du, Fp);

// alternative formula for rp
// Scalar rp = r*p;
		Scalar rp = w*p;
		Scalar nu = - (rp / pFp);

		// Update residual r: r = r + nu * Fp
		// Update solution u: u = u + nu * du
		doubleUpdate(nu,u,du,r,Fp);

		// Add search direction to Orthogonalization set
		if(fetiInfo->nlPrecFlg ) {
			tProject(Fp,beta,pr);
			orthoAddCG(p, pr, pFp);
		} else
			orthoAddCG(p, Fp, pFp);

		// Project: w = P^t r
		tProject(r, alpha, w);

		// Precondition: z = F^-1 w
		errorEstimator = preCondition(w, z);

		if(fetiInfo->primalFlag) outputPrimalResidual(iter+1, deltaF);

		// Print Relative Primal Error and Relative Dual Error
		if((iter % fetiInfo->numPrint() == 0) && (fetiInfo->numPrint() > 0) && verboseFlag ) {
			double wrel = sqrt((w.sqNorm())/w0Norm2);
			//filePrint(stderr," %3d Relative Primal Error = %e "
			//                 " Relative Dual Error = %e\n",
			//iter+1, sqrt(errorEstimator/pseudoFNormSq), wrel);
			filePrint(stderr," %4d %23.6e %21.6e\n",
			          iter+1, sqrt(errorEstimator/pseudoFNormSq), wrel);
		}

		// Test for convergence or maximum number of iterations, this way, 
		// on exit of the loop we are guaranteed that alpha is correct
		// KHPXX
		if(((z2 = errorEstimator) < (epsilon2 * pseudoFNormSq)) || iter == maxiter-1) {
			//if((z2 = z.sqNorm()) < epsilon2 * pseudoFNormSq || iter == maxiter-1) {
			times.iterations[numSystems].stagnated = 0;
			if(errorEstimator < epsilon2 * pseudoFNormSq)
				times.converged = 1;
			if(verboseFlag){
				double wrel = sqrt((w.sqNorm())/w0Norm2);
				filePrint(stderr," %4d %23.6e %21.6e\n",iter+1, sqrt(errorEstimator/pseudoFNormSq), wrel);
			}
			break;
		}

		// Check for stagnation
		if(fabs(z2-lastz2) < 1.0e-6*lastz2) {
			times.setStagnate(numSystems);
			filePrint(stderr,"STAGNATION: Relative Primal Error Reached = "
				          "%e %e %e %e\n",sqrt(z2/pseudoFNormSq), fabs(z2-lastz2),
			          z2, lastz2);
			break;
		}
		lastz2 = z2;

		// Nonlinear preconditioner
		if(fetiInfo->nlPrecFlg ) {
			nlPreCondition(w,z);
		}

		// Re-project
		project(z, alpha, y);

		// Full reorthogonalization: reuse Fp as a temporary variable
		orthogonalize(y,Fp);

		// if Feti2, update y
		if (isFeti2 && crns > 0 && isDynamic == 0) updateFeti2y(Fp, w, gamma,alpha);

		// Re-project for stagnation detection
		project(Fp,alpha,p);
	}

	// Check what alpha in case of Q != I. Update the converged solution u 
	// with the contribution of the rigid body modes: u = u + alpha*R
	if(isDynamic == 0) addR( u, alpha );

	// make solution compatible u = u + deltaU
	makeCompatible( u, deltaU );

	// Store number of iterations, primal error and dual error
	setAndStoreInfo( iter+1, (errorEstimator/pseudoFNormSq),(w.sqNorm())/w0Norm2 );

	// Compute and store memory used by FETI Solver
	if(numSystems == 1) {
		times.memoryFETI  += memoryUsed();
	}

	// end solve timer
	times.solve += getTime();
}

// Code to output Primal Residual to a file for display purposes.
template<class Scalar>
void
GenFetiSolver<Scalar>::outputPrimalResidual(int iter, GenDistrVector<Scalar> &deltaF) const
{
	double normDeltaF = deltaF.norm();
	deltaF.linC( deltaF, (1.0/normDeltaF) );
	// decd->outputPrimal(deltaF, iter+1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::makeCompatible(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &deltaU) const
{
	// u += deltaU 
	u.linAdd( deltaU );
}

template<class Scalar>
Scalar
GenFetiSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &ifrc, GenDistrVector<Scalar> &bf,
                                         GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &lambda) const
{
	startTimerMemory(times.sAndJ, times.memorySAndJ);

	timedParal(times.solveAndJump, nsub, this, &GenFetiSolver<Scalar>::localSolve, u,
	           lambda, ifrc, bf);
	vPat->exchange();
	timedParal(times.solveAndJump, nsub, this, &GenFetiSolver<Scalar>::interfDiffAndDot,
	             bf, lambda);

	// Sum each subdomains dot product contribution
	Scalar ret = 0.0;
	int i;
	for(i = 0; i < nsub; ++i)
		ret += fetiOps[i]->res;

#ifdef DISTRIBUTED
	if(fetiCom) ret = fetiCom->globalSum(ret);
#endif

	stopTimerMemory(times.sAndJ, times.memorySAndJ);

	return ret;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::rebuildInterface(int iSub, GenDistrVector<Scalar> &v) const
{
	sd[iSub]->rebuildInterf(v.subData(subdomains[iSub]->localSubNum()), vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::interfaceDiff(int iSub, GenDistrVector<Scalar> &v) const
{
	sd[iSub]->interfaceJump(v.subData(subdomains[iSub]->localSubNum()), vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::multKbb(int iSub, const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &interfvec,
                               GenDistrVector<Scalar> &deltaU, GenDistrVector<Scalar> &deltaF, bool &errorFlag) const
{
	const Scalar *v1        =         v.subData(subdomains[iSub]->localSubNum());
	Scalar *interfaceVector = interfvec.subData(subdomains[iSub]->localSubNum());
	Scalar *subDeltaU       =    deltaU.subData(subdomains[iSub]->localSubNum());
	Scalar *subDeltaF       = deltaF.subData(subdomains[iSub]->localSubNum());

	if((fetiInfo->version == FetiInfo::fetidp) && (glNumMpc > 0))
		sd[iSub]->multKbbMpc(v1, interfaceVector, subDeltaU, subDeltaF, errorFlag); // also supports coupled_dph
	else if(domain->solInfo().isCoupled) {
		sd[iSub]->multKbbCoupled(v1, interfaceVector, subDeltaF, errorFlag); // doesn't support dual mpc
	}
	else
		sd[iSub]->multKbb(v1, interfaceVector, subDeltaU, subDeltaF, errorFlag);

	subdomains[iSub]->sendInterf(interfaceVector, vPat);
}

template<class Scalar>
double
GenFetiSolver<Scalar>::preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag) const
{
	startTimerMemory(times.precond, times.memoryPrecond);

	// Change this back to arguments of preCondition
	GenDistrVector<Scalar> &deltaU  = wksp->ret_deltaU();
	GenDistrVector<Scalar> &deltaF  = wksp->ret_deltaF();
	double primalResidual = 0.0;

	// Compute preconditioner
	timedParal(times.preconditioner, nsub, this, &GenFetiSolver<Scalar>::multKbb, v,
	             Pv, deltaU, deltaF, errorFlag);
	vPat->exchange();
	timedParal(times.preconditioner,nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, Pv);

	if(errorFlag) {
		// Send deltaF through subdomain buffer
		timedParal(times.preconditioner, nsub, this, &GenFetiSolver<Scalar>::sendDeltaF, deltaF);
		vPat->exchange();

		double *subDots = (double *) dbg_alloca(sizeof(double)*nsub);

		// Compute each partial true norm  
		timedParal(times.preconditioner, nsub, this, &GenFetiSolver<Scalar>::normDeltaF, subDots, &deltaF);

		// Now sum each partial norm over all subdomains
		for(int i = 0; i < nsub; ++i) primalResidual += subDots[i];

#ifdef DISTRIBUTED
		if(fetiCom) primalResidual = fetiCom->globalSum(primalResidual);
#endif

		if(primalResidual < 0.0) {
			// if the primalResidual is smaller than precision then safe to reverse sign
			if(primalResidual > -1e-15) primalResidual = -primalResidual;
			else filePrint(stderr," *** WARNING: negative norm of primal residual %e \n", primalResidual);
		}
	}

	if(fetiInfo->precno == FetiInfo::noPrec) Pv = v; // No preconditioner case

	stopTimerMemory(times.precond, times.memoryPrecond);

	// return true primal residual 
	return primalResidual;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::sendDeltaF(int iSub, GenDistrVector<Scalar> &deltaF) const
{
	// parallel implementation of sending each subdomain's deltaF
	subdomains[iSub]->sendDeltaF(deltaF.subData(iSub), vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::normDeltaF(int iSub, double *subDots, GenDistrVector<Scalar> *deltaF) const
{
	// parallel implementation of computing true norm of deltaF for each subdomain
	subDots[iSub] = subdomains[iSub]->collectAndDotDeltaF(deltaF->subData(iSub), vPat);
}

template<class Scalar>
int
GenFetiSolver<Scalar>::nlPreCondition(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &pr) const
{
	if(fetiInfo->maxortho <= 0 || oSetCG->numOrthoSets() == 1) return 0; // no preconditioner available
	startTimerMemory(times.nlPreCond, times.memoryNlPreCond);

	Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hpr = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

	execParal(nsub,this,&GenFetiSolver<Scalar>::gatherHalfInterface,&r,&pr,hp,hpr);
	oSetCG->precondition(hp, hpr);
	execParal(nsub, this, &GenFetiSolver<Scalar>::scatterHalfInterface, hpr, &pr);
	vPat->exchange();
	execParal(nsub, this, &GenFetiSolver<Scalar>::rebuildInterface, pr);
	stopTimerMemory(times.nlPreCond, times.memoryNlPreCond);

	return 1; // preconditioner available
}

template<class Scalar>
GenFetiWorkSpace<Scalar>::GenFetiWorkSpace(DistrInfo& interface, DistrInfo& local,
                                           int isNonlinear, int numrbms, int numcrns)
{
	// FETI-1 and FETI-2
	zeroPointers();
	// interface = Distributed Subdomain interface size
	// local     = Distributed Subdomain local size (all subdomain dof)

	r       = new GenDistrVector<Scalar>(interface);
	lambda  = new GenDistrVector<Scalar>(interface);
	w       = new GenDistrVector<Scalar>(interface);
	y       = new GenDistrVector<Scalar>(interface);
	z       = new GenDistrVector<Scalar>(interface);
	p       = new GenDistrVector<Scalar>(interface);
	pr      = new GenDistrVector<Scalar>(interface);
	Fp      = new GenDistrVector<Scalar>(interface);
	du      = new GenDistrVector<Scalar>(local);
	uzero   = new GenDistrVector<Scalar>(local);

	wrk1    = new GenDistrVector<Scalar>(interface);
	wrk2    = new GenDistrVector<Scalar>(interface);

	deltaU  = new GenDistrVector<Scalar>(local);
	deltaF  = new GenDistrVector<Scalar>(local);

	alpha   = new GenVector<Scalar>(numrbms);
	beta    = new GenVector<Scalar>(numrbms);
	gamma   = new GenVector<Scalar>(numcrns);
	working = new GenVector<Scalar>(numcrns);

	if(isNonlinear) {
		zz = new GenDistrVector<Scalar>(interface);
		uu = new GenDistrVector<Scalar>(local);
		rCompare = new GenDistrVector<Scalar>(interface);
	}
	else {
		zz       = 0;
		uu       = 0;
		rCompare = 0;
	}
}

template<class Scalar>
GenFetiWorkSpace<Scalar>::GenFetiWorkSpace(DistrInfo& interface, DistrInfo& local, DistrInfo& wet, int ngrbms, int numC, bool contact)
{
	// FETI-DP
	zeroPointers();
	r       = new GenDistrVector<Scalar>(interface);
	lambda  = new GenDistrVector<Scalar>(interface); lambda->zero();
	z       = new GenDistrVector<Scalar>(interface);
	p       = new GenDistrVector<Scalar>(interface);
	Fp      = new GenDistrVector<Scalar>(interface);
	w       = new GenDistrVector<Scalar>(interface);
	y       = new GenDistrVector<Scalar>(interface);
	if(contact) {
		gc      = new GenDistrVector<Scalar>(interface);
		gf      = new GenDistrVector<Scalar>(interface);
		lambda_copy = new GenDistrVector<Scalar>(interface); lambda_copy->zero();
	}

	fr      = new GenDistrVector<Scalar>(local);
	fr2     = new GenDistrVector<Scalar>(local);
	ur      = new GenDistrVector<Scalar>(local);
	du      = new GenDistrVector<Scalar>(local);
	deltaU  = new GenDistrVector<Scalar>(local);
	deltaF  = new GenDistrVector<Scalar>(local);

	fc      = new GenVector<Scalar>(numC);
	uc      = new GenVector<Scalar>(numC);
	duc     = new GenVector<Scalar>(numC);

	alpha   = new GenVector<Scalar>(ngrbms);
	gamma   = new GenVector<Scalar>(ngrbms);
	e       = new GenVector<Scalar>(ngrbms);

	if(domain->solInfo().isCoupled) fw = new GenDistrVector<Scalar>(wet);
}

template<class Scalar>
void
GenFetiWorkSpace<Scalar>::save()
{
	if(!lambda_copy) lambda_copy = new GenDistrVector<Scalar>(*lambda); else *lambda_copy = *lambda;
	if(!r_copy) r_copy = new GenDistrVector<Scalar>(*r); else *r_copy = *r;
	if(!p_copy) p_copy = new GenDistrVector<Scalar>(*p); else *p_copy = *p;
	if(!Fp_copy) Fp_copy = new GenDistrVector<Scalar>(*Fp); else *Fp_copy = *Fp;
	if(!du_copy) du_copy = new GenDistrVector<Scalar>(*du); else *du_copy = *du;
	if(!duc_copy) duc_copy = new GenVector<Scalar>(*duc); else *duc_copy = *duc;
}

template<class Scalar>
void
GenFetiWorkSpace<Scalar>::restore()
{
	*lambda = *lambda_copy;
	*r = *r_copy;
	*p = *p_copy;
	*Fp = *Fp_copy;
	*du = *du_copy;
	*duc = *duc_copy;
}

template<class Scalar>
void
GenFetiWorkSpace<Scalar>::zeroPointers()
{
	r = 0; lambda = 0;
	w = 0; y = 0;
	z = 0; p = 0;
	pr = 0; Fp = 0;
	du = 0; uzero = 0;
	fr = 0; fr2 = 0;
	ur = 0; wrk1 = 0;
	wrk2 = 0; deltaU = 0;
	deltaF = 0; zz = 0;
	rCompare = 0; uu = 0;
	alpha = 0; beta = 0;
	gamma = 0; working = 0;
	fc = 0; uc = 0;
	duc = 0;
	lambda_copy = 0; p_copy = 0; r_copy = 0; Fp_copy = 0; du_copy = 0; uc_copy = 0; duc_copy = 0;
	fw = 0; e = 0;
	gc = 0; gf = 0;
}

template<class Scalar>
void
GenFetiWorkSpace<Scalar>::clean_up()
{
	if(r) r->clean_up();
	if(lambda) lambda->clean_up();
	if(w) w->clean_up();
	if(y) y->clean_up();
	if(z) z->clean_up();
	if(p) p->clean_up();
	if(pr) pr->clean_up();
	if(Fp) Fp->clean_up();
	if(du) du->clean_up();
	if(uzero) uzero->clean_up();
	if(fr) fr->clean_up();
	if(fr2) fr2->clean_up();
	if(ur) ur->clean_up();
	if(wrk1) wrk1->clean_up();
	if(wrk2) wrk2->clean_up();
	if(deltaU) deltaU->clean_up();
	if(deltaF) deltaF->clean_up();
	if(zz) zz->clean_up();
	if(rCompare) rCompare->clean_up();
	if(uu) uu->clean_up();
	if(alpha) alpha->clean_up();
	if(beta) beta->clean_up();
	if(gamma) gamma->clean_up();
	if(working) working->clean_up();
	if(fc) fc->clean_up();
	if(uc) uc->clean_up();
	if(duc) duc->clean_up();
	if(lambda_copy) lambda->clean_up();
	if(p_copy) p_copy->clean_up();
	if(r_copy) r_copy->clean_up();
	if(Fp_copy) Fp->clean_up();
	if(du_copy) du_copy->clean_up();
	if(uc_copy) uc_copy->clean_up();
	if(duc_copy) duc_copy->clean_up();
	if(fw) fw->clean_up();
	if(e) e->clean_up();
	if(gc) gc->clean_up();
	if(gf) gf->clean_up();
}

template<class Scalar>
GenFetiWorkSpace<Scalar>::~GenFetiWorkSpace()
{
	if(r) delete r;
	if(lambda) delete lambda;
	if(w) delete w;
	if(y) delete y;
	if(z) delete z;
	if(p) delete p;
	if(pr) delete pr;
	if(Fp) delete Fp;
	if(du) delete du;
	if(uzero) delete uzero;
	if(fr) delete fr;
	if(fr2) delete fr2;
	if(ur) delete ur;
	if(wrk1) delete wrk1;
	if(wrk2) delete wrk2;
	if(deltaU) delete deltaU;
	if(deltaF) delete deltaF;
	if(zz) delete zz;
	if(rCompare) delete rCompare;
	if(uu) delete uu;
	if(alpha) delete alpha;
	if(beta) delete beta;
	if(gamma) delete gamma;
	if(working) delete working;
	if(fc) delete fc;
	if(uc) delete uc;
	if(duc) delete duc;
	if(lambda_copy) delete lambda_copy;
	if(p_copy) delete p_copy;
	if(r_copy) delete r_copy;
	if(Fp_copy) delete Fp_copy;
	if(du_copy) delete du_copy;
	if(uc_copy) delete uc_copy;
	if(duc_copy) delete duc_copy;
	if(fw) delete fw;
	if(e) delete e;
	if(gc) delete gc;
	if(gf) delete gf;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getCtMult(GenDistrVector<Scalar> &w, GenVector<Scalar> &gamma) const
{
	opControl->vec2 = gamma.data();
	opControl->dv1  = &w;

	opControl->operation = &GenFetiOp<Scalar>::getCtMult;
	threadManager->execParal(nsub,fetiTasks);
}

template<class Scalar>
int
GenFetiSolver<Scalar>::collectIntGlobalSum()
{
	int result = 0;

	int iSub;
	for (iSub=0; iSub < nsub; iSub++)
		result += *((int*) opControl->globalSumBuffer[iSub]);

	return result;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::updateFeti2lambda(GenDistrVector<Scalar> &w, GenDistrVector<Scalar> &r,
                                         GenDistrVector<Scalar> &lambda0, GenVector<Scalar> &beta,
                                         GenVector<Scalar> &alpha) const
{
	// beta = C^t w
	getCtMult(w, beta);
	beta *= -1.0;

	// Now solve it. beta = Fc^{-1} C^t w
	times.forBack2 -= getTime();
	PCtFPC->reSolve( beta.data() );
	times.forBack2 += getTime();

	// update lambda0 = lambda0 + C beta + G (Rgc beta)
	// alpha = Rg * beta

	if(numSystems == 0)
		fprintf(stderr," ... Number of Corners = %5d      ...\n",crns);

	updateFeti2Vector(lambda0, beta, alpha);

	// update r = r - [F_I C] beta - [F_I G] alpha

	GenDistrVector<Scalar> &wrk1 = wksp->ret_wrk1();
	GenDistrVector<Scalar> &wrk2 = wksp->ret_wrk2();

	wrk1.zero();
	opControl->vec2 = beta.data();
	opControl->dv1  = &wrk1;
	opControl->operation = &GenFetiOp<Scalar>::subNuFC;
	threadManager->execParal(nsub,fetiTasks);
	vPat->exchange();

	execParal(nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, wrk1);

	wrk2.zero();

	opControl->vec2 = alpha.data();
	opControl->dv1 = &wrk2;

	opControl->operation = &GenFetiOp<Scalar>::subAlphaFG;
	threadManager->execParal(nsub,fetiTasks);
	vPat->exchange();

	execParal(nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, wrk2);

	r += wrk1;
	r += wrk2;

	// Re-project w

	tProject(r, alpha, w);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::updateFeti2y(GenDistrVector<Scalar> &y, GenDistrVector<Scalar> &w,
                                    GenVector<Scalar> &beta, GenVector<Scalar> &alpha) const
{
	startTimerMemory(times.project2, times.memoryProject2);

	if (PCtFPC) {

		alpha.zero();
		beta.zero();

// Build CtF yk + Rgc^t GtQ yk
// first, CtF yk

		opControl->vec2 = beta.data();
		opControl->dv1 = &y;
		opControl->operation = &GenFetiOp<Scalar>::getCtFMult;
		threadManager->execParal(nsub,fetiTasks);

// Secondly, GtQ yk

		opControl->vec2 = alpha.data();
		opControl->dv1 = &y;
		opControl->operation = &GenFetiOp<Scalar>::getGtFMult;
		threadManager->execParal(nsub,fetiTasks);

// compute Rgc^t alpha, then add to beta

		GenVector<Scalar> &working = wksp->ret_working();
		working.zero();

		computeRgcTransMult(alpha, working);

		beta += working;

// Now solve it.

		times.forBack2 -= getTime();
		PCtFPC->reSolve( beta.data() );
		times.forBack2 += getTime();

		beta *= -1.0;

// update yk by adding C beta and G (Rgc beta)
		updateFeti2Vector(y,beta,alpha);
	}

	stopTimerMemory(times.project2, times.memoryProject2);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::subRgcTransMult(int i, int nThreads, GenVector<Scalar>* alpha,
                                       GenVector<Scalar>* result) const
{
	int dd = crns/nThreads;
	int remainder = crns%nThreads;
	int rowStart = dd*i+( (i < remainder) ? i : remainder );
	int rowStop  = dd*(i+1) +( (i+1 < remainder) ? i+1 : remainder );

	Scalar * Rgc = opControl->Rgc;

	if(rowStop > rowStart)
		Tgemv('T',numrbms,rowStop-rowStart,
		      1.0,Rgc+rowStart*numrbms,numrbms, alpha->data(), 1, 0.0,
		      result->data()+rowStart,1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::computeRgcTransMult(GenVector<Scalar> &alpha,
                                           GenVector<Scalar> &result) const
{
	// compute result =  Rgc^t * alpha
	execParal(threadManager->numThr(), this, &GenFetiSolver<Scalar>::subRgcTransMult,
	          threadManager->numThr(), &alpha, &result);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::subRgcMult(int i, int nThreads, GenVector<Scalar>* alpha,
                                  GenVector<Scalar>* result) const
{
	int dd = numrbms/nThreads;
	int remainder = numrbms%nThreads;
	int colStart = dd*i+( (i < remainder) ? i : remainder );
	int colStop  = dd*(i+1) +( (i+1 < remainder) ? i+1 : remainder );
	Scalar * Rgc = opControl->Rgc;

	if(colStop > colStart)
		Tgemv('N',colStop-colStart,crns,
		      1.0,Rgc+colStart,numrbms, alpha->data(), 1, 0.0,
		      result->data()+colStart,1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::computeRgcMult(GenVector<Scalar> &beta,
                                      GenVector<Scalar> &result) const
{
	execParal(threadManager->numThr(), this, &GenFetiSolver<Scalar>::subRgcMult,
	          threadManager->numThr(), &beta, &result);
}

// v = v + C * beta - G * Rgc * beta

template<class Scalar>
void
GenFetiSolver<Scalar>::updateFeti2Vector(GenDistrVector<Scalar> &v, GenVector<Scalar> &beta,
                                         GenVector<Scalar> &alpha) const
{
	beta *= -1.0;

	// wrk1 = C * beta
	opControl->vec2 = beta.data();
	opControl->dv1 = &v;
	opControl->operation = &GenFetiOp<Scalar>::subNuC;
	threadManager->execParal(nsub,fetiTasks);

	// alpha = Rgc * beta

	computeRgcMult(beta, alpha);

	opControl->vec2 = alpha.data();
	opControl->dv1 = &v;
	opControl->operation = &GenFetiOp<Scalar>::subAlphaG;
	threadManager->execParal(nsub,fetiTasks);
}

extern int zeroFd;
template<class Scalar>
void
GenFetiSolver<Scalar>::addAllFcoarse(GenFullM<Scalar> & mat)
{
	int i,j;

	Scalar        *Rgc   = opControl->Rgc;
	GenFullM<Scalar>  *GtFC  = opControl->GtFC;
	GenSymFullMatrix<Scalar> *GtQG  = opControl->GtQGs;

// first component, Rgc^t GtFC

	GenStackFullM<Scalar> RGC(crns, numrbms, Rgc);

	if((numrbms > 0) && (crns > 0))
		mat.paralAddProd(RGC, *GtFC, 0.0, 1.0);
	else
		mat.paralZero();

// second component, transpose of the previous
// third component, CtFC

	for(i=0; i< crns; i++)
		for(j=0; j<= i; j++)
			mat[j][i] = (mat[i][j] += (mat[j][i]));

	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::assembleCtFCs);

	// Allocate once in case of nonlinear methods
	if(numSystems == 0 )
		gtqglocal = new Scalar[numrbms*numrbms];

	// Now Copy GtQG in full
	for (i=0; i< numrbms; i++)
		for(j = i; j < numrbms; ++j)
			gtqglocal[j*numrbms+i] = gtqglocal[i*numrbms+j] = (*GtQG)[j][i];

	GenStackFullM<Scalar> GTQGLOC(numrbms,numrbms,gtqglocal);

// last component, Rgc^t GtQG Rgc. Store GtQG Rgc to GtFC, 
// since the value of GtFC is not needed anymore.

	if(numrbms > 0 && crns > 0) {
		GtFC->paralAddUTProd(GTQGLOC, RGC, 0.0, 1.0);
		mat.paralAddProd(RGC, *GtFC, 1.0, 1.0);
	}

	// KHP: TESTING
	//if(fetiInfo->type == FetiInfo::linear) {
	//  delete GtFC;
	//  delete GtQG;
	//  delete [] gtqglocal;
	// }

}

// Subdomain level functions for parallel execution

template<class Scalar>
void
GenFetiSolver<Scalar>::gatherHalfInterface(int iSub, const GenDistrVector<Scalar> *v1, const GenDistrVector<Scalar> *v2,
                                           Scalar *v3, Scalar *v4) const
{
	const Scalar *vec1 = v1->subData(subdomains[iSub]->localSubNum());
	const Scalar *vec2 = (v2) ? v2->subData(subdomains[iSub]->localSubNum()) : 0;

	Scalar *vec3 = v3 + halfOffset(iSub);
	Scalar *vec4 = v4 + halfOffset(iSub);

	if(vec2)
		subdomains[iSub]->getHalfInterf(vec1, vec3, vec2, vec4);
	else
		subdomains[iSub]->getHalfInterf(vec1, vec3);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::scatterHalfInterface(int iSub, Scalar *v1, GenDistrVector<Scalar>* v2) const
{
	Scalar *vec1 = v1 + halfOffset(iSub);
	Scalar *vec2 = v2->subData(subdomains[iSub]->localSubNum());

	subdomains[iSub]->scatterHalfInterf(vec1, vec2);
	subdomains[iSub]->sendInterf(vec2, vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::factorMatrices(int iSub)
{
	GenFetiOp<Scalar> *fop = fetiOps[iSub];

	fop->K->factor();

	fop->solver = fop->K;

	// In the dynamic case, we use geometric RBMs to compute the coarse problem
	if(fop->rbm) {
		// Geometric rigid body modes
		fop->numRBM = fop->rbm->numRBM();
		if(fop->numRBM > 0) {
			fop->locRBMs.resize(fop->numRBM * subdomains[iSub]->localLen());
			fop->rbm->getRBMs(fop->locRBMs.data());
		}
	} else {
		// Algebraic rigid body modes
		fop->numRBM = fop->K->numRBM();
		if(fop->numRBM > 0) {
			fop->locRBMs.resize(fop->numRBM * subdomains[iSub]->localLen());
			fop->K->getRBMs(fop->locRBMs.data());
		}
#ifndef SALINAS_FETI
		if(fetiInfo->numPrint() > 0 && verboseFlag)
			fprintf(stderr," ... Subdomain %3d found %3d ZEMs   ...\n",
			        subdomains[iSub]->subNum()+1,fop->K->numRBM());
#endif
	}

	fop->makeCoarseSet();

	// In nonlinear case, we don't need Kii factored until after we rebuild
	if(fetiInfo->type == FetiInfo::linear)
		subdomains[iSub]->factorKii();

}

template<class Scalar>
void
GenFetiSolver<Scalar>::localSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
                                  GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4) const
{
	int sn = subdomains[iSub]->localSubNum();

	Scalar *localvec  = v1.subData(sn);
	Scalar *interfvec = v2.subData(sn);
	Scalar *localsrc  = v3.subData(sn);
	Scalar *interfsrc = v4.subData(sn);

	int i;
	for(i = 0; i < subdomains[iSub]->localLen(); ++i)
		localvec[i] = localsrc[i];

	for(i = 0; i < subdomains[iSub]->interfLen(); ++i)
		interfvec[i] = interfsrc[i];

	subdomains[iSub]->fetiBaseOp(fetiOps[iSub]->solver, localvec, interfvec);

	subdomains[iSub]->sendInterf(interfvec, vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::interfSend(int iSub, GenDistrVector<Scalar> &dv1) const
{
	Scalar *interfvec = dv1.subData(subdomains[iSub]->localSubNum());
	subdomains[iSub]->sendInterf(interfvec, vPat);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::interfDiffAndDot(int iSub, GenDistrVector<Scalar> &dv1, GenDistrVector<Scalar> &dv2) const
{
	// Difference
	Scalar *interfvec = dv2.subData(subdomains[iSub]->localSubNum());
	sd[iSub]->interfaceJump(interfvec, vPat);

	// Dot
	fetiOps[iSub]->res = 0.0;
	Scalar *v1 = dv1.subData(subdomains[iSub]->localSubNum());

	int iLen = subdomains[iSub]->interfLen();
	int i;
	bool *masterFlag = dv1.subMasterFlag(iSub);
	for(i = 0; i < iLen; ++i) {
		if(masterFlag[i])
			fetiOps[iSub]->res += ScalarTypes::Real(v1[i]*interfvec[i]);
	}
}

template<>
void
GenFetiSolver<DComplex>::getRMult(int iSub, GenDistrVector<DComplex> *localvec,
                                  DComplex *alpha) const;
template<>
void
GenFetiSolver<double>::getRMult(int iSub, GenDistrVector<double> *localvec,
                                double *alpha) const;

template<>
void
GenFetiSolver<DComplex>::addRP(int iSub, GenDistrVector<DComplex> *vec1,
                               DComplex *vec2) const;
template<>
void
GenFetiSolver<double>::addRP(int iSub, GenDistrVector<double> *vec1, double
*vec2) const;

template<class Scalar>
void
GenFetiSolver<Scalar>::findProfileSizes(int iSub, int *subSizes)
{
	subSizes[iSub] = sd[iSub]->findProfileSize();
}

template<class Scalar>
int
GenFetiSolver<Scalar>::numRBM()
{
	return glNumRBM;
}


template<class Scalar>
void
GenFetiSolver<Scalar>::Ksolve(int iSub, GenStackDistVector<Scalar> &R)
{
	fetiOps[iSub]->solver->reSolve(R.subData(iSub));
}

template<class Scalar>
void
GenFetiSolver<Scalar>::clean_up()
{
	long m1=0;

	m1 = -memoryUsed();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::clean_up);
	m1 += memoryUsed();
	std::cerr << std::endl
	          << "Deleted Local Stiffness Matrices :"
	          << -m1/(1024.0*1024.0) << " Mb" << std::endl;

/*
 if(fetiOps) { 
   delete [] fetiOps;
   fetiOps=0;
 }

 if(fetiTasks) {
   delete [] fetiTasks;
   fetiTasks=0;
 }
*/

	if(opControl) {
		delete [] opControl;
		opControl=0;
	}

	//cerr << "Size of FetiOps pointer = " << sizeof(FetiOp*) << std::endl;
	//cerr << "Deleted FetiOps memory  = " << -m1 << std::endl;

	m1 = -memoryUsed();
	if(GtGsolver) GtGsolver->clean_up();
	m1 += memoryUsed();
	std::cerr << std::endl
	          << "Deleted Coarse Problem           :"
	          << -m1/(1024.0*1024.0) << " Mb" << std::endl;

	m1 = -memoryUsed();
	if(wksp) wksp->clean_up();
	m1 += memoryUsed();
	std::cerr << std::endl
	          << "Deleted FETI Work Space          :"
	          << -m1/(1024.0*1024.0) << " Mb" << std::endl;

	m1 = -memoryUsed();
	if(oSetCG) oSetCG->clean_up();
	m1 += memoryUsed();
	std::cerr << std::endl
	          << "Deleted ReOrtho Vectors          :"
	          << -m1/(1024.0*1024.0) << " Mb" << std::endl << std::endl;

}

template<class Scalar>
int GenFetiSolver<Scalar>::numNeighbor(int iSub) const
{ return sd[iSub]->getSComm()->numNeighb; }

template<class Scalar>
void
GenFetiSolver<Scalar>::makeGandFG()
{
	// Get number of each neighbors rbms
	makeRbmPat();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::sendNumNeighbRBM, sPat);
	sPat->exchange();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbRBM, sPat);

	// Get all the numbers of rigid body modes and dispatch RBMs to neighbors
	makeRbmPat();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::sendInterfRBM, rbmPat);
	rbmPat->exchange();

	// compute neighbors QGs
	execParal(nsub, this, &GenFetiSolver<Scalar>::getNeighbFGs);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getNeighbFGs(int iSub)
{
	int numRBM = opControl->cset[iSub].numGs;
	opControl->cset[iSub].locFGs = new Scalar[numRBM*subdomains[iSub]->interfLen()];
	subdomains[iSub]->multMFi(fetiOps[iSub]->solver, opControl->cset[iSub].locGs,
	                          opControl->cset[iSub].locFGs, numRBM);

	opControl->cset[iSub].getNeighbQGs(opControl->cset, sd[iSub],
	                                   0, rbmPat, fetiOps[iSub]->solver);
}

template<class Scalar>
Connectivity *
GenFetiSolver<Scalar>::makeSingleConnect(const Connectivity *coarseConnect,
                                         const Connectivity *coarseToSub,
                                         const Connectivity *subToCoarse,
                                         int gOffset)
{
	// This routine creates a modified connectivity in which the
	// non-zero entries added by the GtG and GtC portion of the coarse is added
#ifdef DISTRIBUTED
	int glNumSub = subToSub->csize();
#else
	int glNumSub = nsub;
#endif
	// coarseSize and gOffset should be equal!!!
	int coarseSize = coarseToSub->csize();
	int fp = 0;
	int iSub;
	int *cp = new int[coarseSize + glNumSub+1];
	for(iSub = 0; iSub < coarseSize; ++iSub) {
		cp[iSub] = fp;
		fp += coarseConnect->num(iSub) + coarseToSub->num(iSub);
	}
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		cp[iSub+gOffset] = fp;
		fp += subToCoarse->num(iSub);
	}
	cp[glNumSub + gOffset] = fp;

	int *ct = new int[fp];
	fp = 0;
	for(iSub = 0; iSub < coarseSize; ++iSub) {
		int j;
		for(j =0; j < coarseConnect->num(iSub); ++j)
			ct[fp++] = (*coarseConnect)[iSub][j];
		for(j =0; j < coarseToSub->num(iSub); ++j)
			ct[fp++] = (*coarseToSub)[iSub][j] + gOffset;
	}

	for(iSub = 0; iSub < glNumSub; ++iSub) {
		int j;
		for(j =0; j < subToCoarse->num(iSub); ++j)
			ct[fp++] = (*subToCoarse)[iSub][j];
	}
	return new Connectivity(coarseSize+glNumSub, cp, ct);
}

template<class Scalar>
Connectivity
GenFetiSolver<Scalar>::getCoarseToSubConnect() const
{
	return subToSub->append(subToSub->append(*mpcToSub));
}

template<class Scalar>
void
GenFetiSolver<Scalar>::makeSingleCoarse()
{
	startTimerMemory(times.coarse1, times.memoryGtG);

	// first create the Gs and FGs
	makeGandFG();

	// make C (Corner modes) and FC
	if(isFeti2) preProcessCorners();

	int glNumSub = subToSub->csize();
	int glNumMpc = (mpcToSub) ? mpcToSub->csize() : 0;

	Connectivity coarseToSub = getCoarseToSubConnect();
	Connectivity subToCoarse = coarseToSub.reverse();

	// GtFG is driving the renumbering both for the GtG system and the GtG
	// off diagonal terms which are then inserted to interlay with GtFG
	Connectivity coarseFcoarseConnect = coarseToSub.transcon(subToCoarse);
	Connectivity *gtFgConnect = subToSub->transcon(subToSub);

	// We now renumber gtFgConnect on the assumption that it has no 'spurious'
	// singularity
	compStruct gRenum = gtFgConnect->renumByComponent(1);

	// now revisit the numbering by inserting the GtG equations to follow their
	// GtFG counterpart
	int coarseSize  = coarseFcoarseConnect.csize();
	int glRenumSize = glNumSub+coarseSize;
	int *glRenum = new int[glRenumSize];

	gOffset = coarseSize;
	mOffset = coarseSize - glNumMpc;
	// YYY mOffset = 0; 
	int cOffset = glNumSub;
	// YYY int cOffset = glNumSub+glNumMpc;
	// int gtfgOffset=glNumMpc;
	int gtfgOffset=0;

	// Experimental forced renumbering
	int i, c;
	int *invRen = new int[glNumSub];
	for(i = 0; i < glNumSub; ++i)
		invRen[gRenum.renum[i]] = i;

	c = 0;
	// YYY c = glNumMpc;
	for(i = 0; i < glNumSub; ++i) {
		glRenum[c] = invRen[i]+gtfgOffset;
		//cerr << "GtFG: before: glRenum["<<c<<"] = " << glRenum[c] << std::endl;
		c++;
		glRenum[c] = invRen[i]+cOffset;
		//cerr << "   C: before: glRenum["<<c<<"] = " << glRenum[c] << std::endl;
		c++;
		glRenum[c] = invRen[i]+gOffset;
		//cerr << " GtG: before: glRenum["<<c<<"] = " << glRenum[c] << std::endl;
		c++;
	}

#ifdef DEBUG_MPC
	for(i = 0; i < glRenumSize; ++i)
   std::cerr << "before MPCs: glRenum[" << i << "] = " << glRenum[i] << std::endl;

 // Number the multiple point constraints separately
 std::cerr << " cOffset = " << cOffset 
      << " gOffset = " << gOffset
      << " mOffset = " << mOffset
      << std::endl;
#endif

	int location = 3*glNumSub;
	for(i=0; i<glNumMpc; ++i) {
		glRenum[location+i] = mOffset + i;
	}

	int *ngRen =  new int[glRenumSize];
	for(i = 0; i < glRenumSize; ++i) {
		ngRen[glRenum[i]] = i;
	}

	// Construct Coarse Problems equation numbers
	if(eqNums == 0) eqNums = new SimpleNumberer(glRenumSize, ngRen);
	int iSub;

	int *gRbmSize = new int[glNumSub];
	int *cornerSize = new int[glNumSub];
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		gRbmSize[iSub] = cornerSize[iSub] = 0;
	}

	numrbms = 0;
	for(iSub = 0; iSub < nsub; ++iSub) {
		gRbmSize[ subdomains[iSub]->subNum() ]  = opControl->cset[iSub].numGs;
		cornerSize[subdomains[iSub]->subNum() ] = fetiOps[iSub]->getcrnDofSize();
		numrbms += opControl->cset[iSub].numGs;
	}

#ifdef DISTRIBUTED
	// Now global sum of gRbmSize and cornerSize
 fetiCom->globalSum(glNumSub, gRbmSize);
 fetiCom->globalSum(glNumSub, cornerSize);
 
 numrbms=0;
 for(iSub = 0; iSub < glNumSub; ++iSub) 
   numrbms += gRbmSize[iSub];
#endif

	int nC=0;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		eqNums->setWeight(iSub + gtfgOffset, gRbmSize[iSub]);
//    fprintf(stderr,"GtFG: set Eqn %d  to %d\n",iSub + gtfgOffset, gRbmSize[iSub]);
		eqNums->setWeight(iSub + gOffset,    gRbmSize[iSub]);
//    fprintf(stderr," GtG: set Eqn %d  to %d\n",iSub + gOffset, gRbmSize[iSub]);
		eqNums->setWeight(iSub + cOffset,    cornerSize[iSub]);
//    fprintf(stderr,"   C: set Eqn %d  to %d\n",iSub + cOffset, cornerSize[iSub]);
		nC += cornerSize[iSub];
	}

	// add weight for each MPC equation
	for(i=0; i<glNumMpc; ++i) {
		eqNums->setWeight(i + mOffset, 1);
#ifdef DEBUG_MPC
		fprintf(stderr," MPC: set Eqn %d  to %d\n",i + mOffset,1);
#endif
	}

/*
 filePrint(stderr,"----------------------------------------------\n");
 filePrint(stderr,"Total number of rigid body modes          %d\n",numrbms);
 filePrint(stderr,"Total Number of Corner Degrees of Freedom %d\n",nC);
 filePrint(stderr,"Total Number of MPCs                      %d\n",glNumMpc);
 filePrint(stderr,"----------------------------------------------\n");
 filePrint(stderr,"Total Size of FETI-2 Coarse Grid:         %d\n",nC+numrbms+glNumMpc);
 filePrint(stderr,"----------------------------------------------\n");
*/


	eqNums->makeOffset();
	// The following sets the alpha offset for the subdomains AND in the
	// case of Feti2 dynamic it does gOffset as well

	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::setAlphaOffsets, eqNums->allOffsets());

	// create the combined connectivity
	Connectivity *coarseConnect = makeSingleConnect(&coarseFcoarseConnect,
	                                                &coarseToSub, &subToCoarse, gOffset);
#ifdef DEBUG_MPC
	eqNums->print();
   coarseConnect->print();
#endif
	// exit(-1);

	// Now we are ready to make the coarse solver and to assemble it
	if(numrbms > 0 || crns > 0 || glNumMpc > 0) {

		double tolerance = 1e-05; //fetiInfo->grbm_tol; // default = 1.0E-06

		times.memoryGtGsky -= threadManager->memoryUsed();
#ifdef DISTRIBUTED
		if(subToSub->csize() == numCPUs)
     singleCoarseSolver = GenSolverFactory<Scalar>::getFactory()->createDistSolver(coarseConnect, eqNums, *fetiInfo->coarse_cntl, singleCoarse, this->fetiCom);
   else
#endif
		singleCoarseSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(coarseConnect, eqNums, *fetiInfo->coarse_cntl, singleCoarse, 0, this->fetiCom);

		times.memoryGtGsky += threadManager->memoryUsed();
	}

	times.numRBMs = eqNums->size();

	filePrint(stderr," ... Size of Coarse Problem %5d   ...\n",times.numRBMs);

	if( times.numRBMs != 0 || crns != 0 || glNumMpc != 0) {
		singleCoarseAssembly();

#ifdef DISTRIBUTED
		singleCoarseSolver->unify(fetiCom);
#endif

		startTimerMemory(times.pfactor, times.memoryGtGsky);
		// parallelFactor has a bug when there is a singularity in the coarse matrix
		singleCoarseSolver->parallelFactor();
		stopTimerMemory(times.pfactor, times.memoryGtGsky);

		filePrint(stderr, " ... Factoring coarse found %2d RBMS ...\n",
		          singleCoarseSolver->numRBM());
	} else {
		singleCoarse=0;
		singleCoarseSolver = 0;
	}

	stopTimerMemory(times.coarse1, times.memoryGtG);

}

template<class Scalar>
void
GenFetiSolver<Scalar>::singleCoarseAssembly()
{
	execParal(nsub, this, &GenFetiSolver<Scalar>::singleCoarseAssembleG);
	execParal(nsub, this, &GenFetiSolver<Scalar>::preprocessMPCs);
	execParal(mpcToSub->csize(), this, &GenFetiSolver<Scalar>::singleCoarseAssembleMPCs);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::preprocessMPCs(int iSub)
{
	subdomains[iSub]->getQtKQ(fetiOps[iSub]->solver);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::singleCoarseAssembleG(int iSub)
{
	opControl->cset[iSub].addGContrib(singleCoarse, 0, eqNums, gOffset);
	int cOffset = subToSub->csize();

	if(isFeti2) {
		opControl->cset[iSub].addCContrib(singleCoarse,0,eqNums, cOffset, gOffset);
		opControl->cset[iSub].addCGContrib(singleCoarse, eqNums, cOffset, gOffset);
	}

	int jSub;
	for(jSub = 0; jSub <  opControl->numNeighbSubs(iSub); ++jSub) {
#ifdef DISTRIBUTED
		int neighbN = glSubToLoc[opControl->neighbSubId(iSub, jSub)]; 
   if(neighbN >= 0) {
#else
		int neighbN = opControl->neighbSubId(iSub, jSub);
#endif
		opControl->cset[neighbN].addGContrib(singleCoarse,
		                                     sd[iSub]->scomm->remoteId[jSub]+1, eqNums, gOffset);
		if(isFeti2)
			opControl->cset[neighbN].addCContrib(
				singleCoarse, sd[iSub]->scomm->remoteId[jSub]+1, eqNums,
				cOffset, gOffset);
#ifdef DISTRIBUTED
		} else {
     addNonLocalGContrib(subdomains[iSub]->subNum(),opControl->neighbSubId(iSub,jSub));

     if(isFeti2)
      addNonLocalCContrib(subdomains[iSub]->subNum(),opControl->neighbSubId(iSub,jSub));
   }
#endif
	}

}

template<class Scalar>
void
GenFetiSolver<Scalar>::singleCoarseAssembleMPCs(int iMPC)
{

	// Corner equation number offset
	int cOffset = subToSub->csize();

	// number of subdomain attached to this mpc
	int numSubAttached = mpcToSub->num(iMPC);

	// loop over all subdomains attached to this mpc and add
	// their contributions to the coarse problem matrix
	int iSub;
	for(iSub=0; iSub<numSubAttached; ++iSub) {
#ifdef DISTRIBUTED
		int myNum = glSubToLoc[(*mpcToSub)[iMPC][iSub]] ;
   if(myNum >= 0)
#else
		int myNum = (*mpcToSub)[iMPC][iSub];
#endif
		opControl->cset[myNum].addMPCContrib(iMPC,
		                                     singleCoarse, eqNums, cOffset, gOffset,
		                                     mOffset+iMPC,
		                                     fetiOps[myNum]->locRBMs.data(),
		                                     fetiOps[myNum]->solver);
	}

}


template<class Scalar>
void
GenFetiSolver<Scalar>::computeL0(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &r,
                                 GenDistrVector<Scalar> &u, GenVector<Scalar> &alpha,
                                 GenDistrVector<Scalar> &lambda) const
{
	// note: r = B K^+ f
#ifdef DISTRIBUTED
	alpha.zero();
#endif

	Scalar *singleC = alpha.data();
	opControl->vec2 = singleC;
	opControl->dv1  = &r;

	execParal(mpcToSub->csize(), this, &GenFetiSolver<Scalar>::getSQtMult, &u, singleC);
	execParal(mpcToSub->csize(), this, &GenFetiSolver<Scalar>::addMpcRhs, singleC);

	// -Rt (f - B^t lambda)   NOTE: lambda could be zero in this computation
	execParal(nsub, this, &GenFetiSolver<Scalar>::getSRMult,  &f, &lambda, singleC);

	execParal(nsub, this, &GenFetiSolver<Scalar>::getGtMult,  &r, singleC);

	if(isFeti2)
		execParal(nsub, this, &GenFetiSolver<Scalar>::getSCtMult, &r, singleC);

#ifdef DISTRIBUTED
	if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif

// filePrint(stderr,"Before: alpha*alpha %e\n",alpha*alpha);
	if(singleCoarseSolver) singleCoarseSolver->reSolve(singleC);
// filePrint(stderr," After: alpha*alpha %e\n",alpha*alpha);

	// Not needed because done on the outside: lambda.zero(); 
	execParal(nsub, this, &GenFetiSolver<Scalar>::addG, &lambda, singleC);

	if(isFeti2)
		execParal(nsub, this, &GenFetiSolver<Scalar>::addC, &lambda, singleC);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addMpcRhs(int iMPC, Scalar *sv) const
{
	Scalar *gamma = sv + eqNums->firstdof(mOffset+iMPC);
#ifdef DISTRIBUTED
	int myNum = glSubToLoc[(*mpcToSub)[iMPC][0]];
#else
	int myNum = (*mpcToSub)[iMPC][0];
#endif
	if(myNum >= 0)
		gamma[0] += subdomains[myNum]->getMpcRhs(iMPC);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::singlePr(GenDistrVector<Scalar> &y, GenDistrVector<Scalar> &p, GenVector<Scalar> &alpha) const
{
	startTimerMemory(times.project, times.memoryProject1);

	alpha.zero();

	Scalar *singleC = alpha.data();

	execParal(mpcToSub->csize(), this, &GenFetiSolver<Scalar>::getQtKpBMult, &y, singleC);

	execParal(nsub, this, &GenFetiSolver<Scalar>::getSGtMult,   &y, singleC);

	execParal(nsub, this, &GenFetiSolver<Scalar>::getFGMult,    &y, singleC);

	if(isFeti2)
		execParal(nsub, this, &GenFetiSolver<Scalar>::getFCMult,  &y, singleC);

#ifdef DISTRIBUTED
	if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif

	if(singleCoarseSolver) singleCoarseSolver->reSolve(singleC);

	p = y;
	execParal(nsub, this, &GenFetiSolver<Scalar>::addG, &p, singleC);

	if(isFeti2)
		execParal(nsub, this, &GenFetiSolver<Scalar>::addC, &p, singleC);

	stopTimerMemory(times.project, times.memoryProject1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addGs(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w, GenVector<Scalar> &alpha) const
{
	Scalar *singleC = alpha.data();
	w = r;
	if(isDynamic==0)
		execParal(nsub, this, &GenFetiSolver<Scalar>::addSG, &w, singleC);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getSRMult(int iSub, GenDistrVector<Scalar> *r, GenDistrVector<Scalar> *lambda, Scalar *sv) const
{
	int nRBM = fetiOps[iSub]->numRBM;

	if((nRBM==0) || isDynamic) return;

	// alpha = R*lvec
	Scalar *lvec = r->subData(iSub);
	Scalar *lbvec = lambda->subData(iSub);
	auto locRBMs = fetiOps[iSub]->locRBMs.data();
	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum()+gOffset);
	sd[iSub]->getSRMult(lvec, lbvec, nRBM, locRBMs, alpha);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addC(int iSub, GenDistrVector<Scalar> *lambda, Scalar *sv) const
{
	Scalar *localvec = lambda->subData(iSub);
	int cOffset      = subToSub->csize();

	opControl->cset[iSub].addNuC(localvec, sv, eqNums, cOffset);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getSGtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	int numRBM = opControl->cset[iSub].numGs;
	if(numRBM==0) return;
	Scalar *lvec = r->subData(iSub);
	Scalar *locGs = opControl->cset[iSub].locGs;

	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum()+gOffset);

	if(numRBM > 0)
		Tgemv('T',opControl->cset[iSub].gSize, numRBM, -1.0, locGs,
		      opControl->cset[iSub].gSize,lvec, 1, 0.0, alpha, 1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getSCtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	int cOffset = subToSub->csize();
	Scalar *lvec = r->subData(iSub);
	Scalar *beta = sv + eqNums->firstdof(subdomains[iSub]->subNum()+cOffset);
	opControl->cset[iSub].getCtMult(lvec,beta);
/* IF different sign convention MLX
 int i;
 for(i = 0; i < eqNums->weight(subdomains[iSub]->subNum()+cOffset); ++i)
   beta[i] = -beta[i];
*/
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getSQtMult(int iMpc, GenDistrVector<Scalar> *u, Scalar *sv) const
{
	int numSubAttached = mpcToSub->num(iMpc);
	int iSub;
	Scalar *gamma = sv + eqNums->firstdof(mOffset + iMpc);
	gamma[0]=0.0;
	for(iSub=0; iSub<numSubAttached; ++iSub) {
#ifdef DISTRIBUTED
		int myNum = glSubToLoc[(*mpcToSub)[iMpc][iSub]];
   if(myNum >= 0) {
     Scalar *lvec = u->subData(myNum);
     sd[myNum]->multQt(iMpc, lvec, gamma);
   }
#else
		int myNum = (*mpcToSub)[iMpc][iSub];
		Scalar *lvec = u->subData(myNum);
		sd[myNum]->multQt(iMpc, lvec, gamma);
#endif
	}

	gamma[0]=-gamma[0];
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getQtKpBMult(int iMpc, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	// Q^t K^+ B^t r
	int numSubAttached = mpcToSub->num(iMpc);
	int iSub;
	Scalar *gamma = sv + eqNums->firstdof(mOffset+iMpc);
	gamma[0] = 0.0;
	for(iSub=0; iSub<numSubAttached; ++iSub) {
#ifdef DISTRIBUTED
		int myNum = glSubToLoc[(*mpcToSub)[iMpc][iSub]];
     if(myNum >= 0) {
       Scalar *lvec = r->subData(myNum);
       sd[myNum]->multQtKBt(iMpc, lvec, gamma);
     }
#else
		int myNum = (*mpcToSub)[iMpc][iSub];
		Scalar *lvec = r->subData(myNum);
		sd[myNum]->multQtKBt(iMpc, lvec, gamma);
#endif
	}

	gamma[0]=-gamma[0];
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getGtMult(int iSub, const GenDistrVector<Scalar> *r, Scalar *sv) const
{
	int numRBM = opControl->cset[iSub].numGs;
	if(numRBM==0) return;
	const Scalar *lvec  = r->subData(iSub);
	Scalar *locGs = opControl->cset[iSub].locGs;
	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum());

	if(opControl->cset[iSub].gSize > 0)
		Tgemv('T',opControl->cset[iSub].gSize, numRBM, -1.0, locGs,
		      opControl->cset[iSub].gSize,lvec, 1, 0.0, alpha, 1);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	Scalar *lvec  = r->subData(iSub);
	Scalar *locGs = opControl->cset[iSub].locGs;
	int numRBM    = opControl->cset[iSub].numGs;
	SComm *scomm  = sd[iSub]->getSComm();
	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum());

	// lvec += locGs*alpha
	if(numRBM > 0 && opControl->cset[iSub].gSize > 0)
		Tgemv('N',opControl->cset[iSub].gSize, numRBM, 1.0, locGs,
		      opControl->cset[iSub].gSize, alpha, 1, 1.0, lvec, 1);
	int jSub;
	for(jSub = 0; jSub < opControl->numNeighbSubs(iSub); ++jSub) {
		alpha = sv + eqNums->firstdof(scomm->subNums[jSub]);

#ifdef DISTRIBUTED
		opControl->cset[iSub].addAlphaNeighbG(jSub, lvec,
               opControl->cset[iSub].neighbGs[jSub],
               opControl->cset[iSub].leadingDimGs[jSub], alpha);
#else
		int nghb = opControl->neighbSubId(iSub, jSub);
		int subOffset = opControl->cset[iSub].subOffset[jSub];
		int subSize = opControl->cset[iSub].subSize(jSub);
		int nRBM = opControl->cset[iSub].neighbNumRBMs[jSub];
		if(nRBM > 0)
			Tgemv('N', subSize, nRBM, -1.0,
			      opControl->cset[nghb].subG(iSub), opControl->cset[nghb].gSize,
			      alpha, 1, 1.0, lvec+subOffset, 1);
#endif
	}

}

template<class Scalar>
void
GenFetiSolver<Scalar>::addSG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	Scalar *lvec  = r->subData(iSub);
	Scalar *locGs = opControl->cset[iSub].locGs;
	int numRBM    = opControl->cset[iSub].numGs;
	SComm *scomm  = sd[iSub]->getSComm();
	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum()+gOffset);

	if(numRBM > 0 && opControl->cset[iSub].gSize > 0)
		Tgemv('N',opControl->cset[iSub].gSize, numRBM, 1.0, locGs,
		      opControl->cset[iSub].gSize,alpha, 1, 1.0, lvec, 1);

	int jSub;
	for(jSub = 0; jSub < opControl->numNeighbSubs(iSub); ++jSub) {
		alpha = sv + eqNums->firstdof(scomm->subNums[jSub]+gOffset);
#ifdef DISTRIBUTED
		opControl->cset[iSub].addAlphaNeighbG(jSub, lvec,
           opControl->cset[iSub].neighbGs[jSub],
           opControl->cset[iSub].leadingDimGs[jSub], alpha);
#else
		int nghb = opControl->neighbSubId(iSub, jSub);
		int subOffset = opControl->cset[iSub].subOffset[jSub];
		int subSize = opControl->cset[iSub].subSize(jSub);
		int nRBM = opControl->cset[iSub].neighbNumRBMs[jSub];
		if(nRBM > 0)
			Tgemv('N', subSize, nRBM, -1.0,
			      opControl->cset[nghb].subG(iSub), opControl->cset[nghb].gSize,
			      alpha, 1, 1.0, lvec+subOffset, 1);
#endif
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getFGMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	Scalar *lvec = r->subData(iSub);
	Scalar *locFGs = opControl->cset[iSub].locFGs;

	int numRBM = opControl->cset[iSub].numGs;
	int gSize = opControl->cset[iSub].gSize;

	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum());

	if(numRBM > 0)
		Tgemv('T', gSize, numRBM, -1.0,
		      locFGs, gSize, lvec, 1, 0.0, alpha, 1);

	int jSub;
	for(jSub = 0; jSub <  opControl->numNeighbSubs(iSub); ++jSub) {

#ifdef DISTRIBUTED
		int neighbN = glSubToLoc[opControl->neighbSubId(iSub, jSub)];
     if(neighbN >= 0) {
#else
		int neighbN = opControl->neighbSubId(iSub, jSub);
#endif
		int myID = sd[iSub]->scomm->remoteId[jSub];
		Scalar *ovec = r->subData(neighbN);
		gSize = subdomains[neighbN]->interfLen();
		if(numRBM > 0)
			Tgemv('T', gSize, numRBM, -1.0,
			      opControl->cset[neighbN].neighbQGs[myID], gSize,
			      ovec, 1, 1.0, alpha, 1);

#ifdef DISTRIBUTED
		} else {
       getNonLocalSubAlphaGtQ(subdomains[iSub]->subNum(),
            opControl->neighbSubId(iSub, jSub), sv, r);
     }
#endif
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getFCMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const
{
	int cOffset = subToSub->csize();
	Scalar *lvec = r->subData(iSub);
	Scalar *locFCs = opControl->cset[iSub].locFBCs;

	int numBCs = opControl->cset[iSub].numBCs;
	if(numBCs <= 0 ) return;
	int gSize = opControl->cset[iSub].gSize;

	Scalar *alpha = sv + eqNums->firstdof(subdomains[iSub]->subNum()+cOffset);

	if(gSize > 0)
		Tgemv('T', gSize, numBCs, -1.0,
		      locFCs, gSize, lvec, 1, 0.0, alpha, 1);

	int jSub;
	for(jSub = 0; jSub <  opControl->numNeighbSubs(iSub); ++jSub) {

#ifdef DISTRIBUTED
		int neighbN = glSubToLoc[opControl->neighbSubId(iSub, jSub)];
     if(neighbN >= 0) {
#else
		int neighbN = opControl->neighbSubId(iSub, jSub);
#endif
		int myID = sd[iSub]->scomm->remoteId[jSub];
		Scalar *ovec = r->subData(neighbN);
		opControl->cset[neighbN].getCtFMult(myID, ovec, alpha);

#ifdef DISTRIBUTED
		} else {
       getNonLocalFCtMult(subdomains[iSub]->subNum(),
            opControl->neighbSubId(iSub, jSub), sv, r);
     }
#endif
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addRSingle(GenDistrVector<Scalar> &v, GenVector<Scalar> &alpha) const
{

	times.addR -= getTime();
	execParal(nsub, this, &GenFetiSolver<Scalar>::addRS, &v, alpha.data());
	times.addR += getTime();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::singleCoarseSolve(const GenDistrVector<Scalar> &_f, GenDistrVector<Scalar> &u) const
{
	times.solve -= getTime();
	GenDistrVector<Scalar> f{_f};
	// K H P: temporary fix for MPCs
	// resetOrthoSet();

	GenDistrVector<Scalar> &r       = wksp->ret_r();          // residual
	GenDistrVector<Scalar> &lambda0 = wksp->ret_lambda();    // Lagrange multipliers
	GenDistrVector<Scalar> &w       = wksp->ret_w();          // projected residual
	GenDistrVector<Scalar> &y       = wksp->ret_y();
	GenDistrVector<Scalar> &z       = wksp->ret_z();          // preconditioned residual
	GenDistrVector<Scalar> &p       = wksp->ret_p();
	GenDistrVector<Scalar> &Fp      = wksp->ret_Fp();
	GenDistrVector<Scalar> &du      = wksp->ret_du();
	//GenDistrVector<Scalar> &pr      = wksp->ret_pr();
	GenDistrVector<Scalar> &uzero   = wksp->ret_uzero();
	uzero.zero();

	GenDistrVector<Scalar> &deltaU  = wksp->ret_deltaU();
	deltaU.zero();
	GenDistrVector<Scalar> &deltaF  = wksp->ret_deltaF();
	deltaF.zero();
	GenVector<Scalar> &alpha        = wksp->ret_alpha();

	// Compute the square of a pseudo-norm of f
	double pseudoFNormSq = f.sqNorm();
	if(pseudoFNormSq == 0.0) {
		u.zero();
		return;
	}

	GenVector<Scalar> beta(alpha.size(),0.0);

	lambda0.zero();

	// First compute a forward backward to get B K^+ f
	// u = K^+ (f + B^t lambda0) = K^+ f
	// r = B u = B K^+ f
	localSolveAndJump(f,lambda0,u,r);

	computeL0(f, r, u, beta, lambda0);

	if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10) {
		filePrint(stderr," ... lambda0*lambda0 = %e ...\n", lambda0.sqNorm());
		filePrint(stderr," ... f*f             %10.4e     ...\n",pseudoFNormSq);
		filePrint(stderr," ... r*r             %10.4e     ...\n",r.sqNorm());
		filePrint(stderr," ... u*u             %10.4e     ...\n",u.sqNorm());
	}

	// u = K^+ (f + B^t lambda0 + Q beta)
	// r = B u
	//   = B K^+ f + B K^+ B^t lambda0 + B K^+ Q beta
	//   = B K^+ f + F lambda0 + B K^+ Q beta
	/*Scalar xx =*/ localSolveAndJump(f,lambda0,beta,u,r);

	addGs(r, w, beta);

	// For multi-solve cases 
	bool hasPred = predict(w,lambda0);
	if(hasPred) {
		// singlePr(lambda0, Fp, beta);
		localSolveAndJump(f, lambda0, u, r);
		Fp = lambda0;
		computeL0( f, r, u, beta, Fp);
		localSolveAndJump(f,Fp,beta,u,r);
		addGs(r, w, beta);
		if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
			filePrint(stderr, " ... lambda0*lambda0 = %e beta*beta %e w*w %e\n",
			          lambda0.sqNorm(), beta.sqNorm(), w.sqNorm());
	}

	double w0Norm2 = w.sqNorm();

	// Precondition: z = F^-1 w
	double errorEstimator = preCondition(w,z);

	if((fetiInfo->numPrint() > 0) && (fetiInfo->numPrint() < 10)) {
		filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");
		filePrint(stderr," %4d %23.6e %21.6e\n",0, sqrt(errorEstimator/pseudoFNormSq), 1.0);
	}

	if(errorEstimator == 0.0) return;

	orthogonalize(z,y);

	singlePr(y, p, alpha);

	// Iterate
	int iter;
	double lastz2 = 0.0, z2;
	for(iter = 0; iter < maxiter; ++iter) {
		Scalar pFp = localSolveAndJump(uzero, p, alpha, du, Fp);

		Scalar rp = w*y;

		r = y;
		// compute y = Fp (true Fp where F is the condensed operator)
		addGs(Fp, y, alpha);
		pFp = y*r;
		Scalar nu = - rp / pFp ;

		lambda0.linAdd(nu,r);
		// Update residual w: w = w + nu * y
		// Update solution u: u = u + nu * du
		doubleUpdate(nu,u,du,w,y);

		// Add search direction to Orthogonalization set
		orthoAddCG(r, y, pFp);

		alpha *= nu;
		beta  += alpha;

		// Precondition: z = F^-1 w
		errorEstimator = preCondition(w, z);

		if((iter % fetiInfo->numPrint() == 0) && (fetiInfo->numPrint() > 0) && verboseFlag ) {
			double wrel = sqrt((w.sqNorm())/w0Norm2);
			filePrint(stderr," %4d %23.6e %21.6e\n",
			          iter+1, sqrt(errorEstimator/pseudoFNormSq), wrel);
		}

		// Test for convergence or maximum number of iterations, this way,
		// on exit of the loop we are guaranteed that alpha is correct
		if(((z2 = errorEstimator) < (epsilon2 * pseudoFNormSq)) || iter == maxiter-1) {
			times.iterations[numSystems].stagnated = 0;
			if(errorEstimator < epsilon2 * pseudoFNormSq)
				times.converged = 1;
			if(verboseFlag){
				double wrel = sqrt((w.sqNorm())/w0Norm2);
				filePrint(stderr," %4d %23.6e %21.6e\n",iter+1, sqrt(errorEstimator/pseudoFNormSq), wrel);
			}
			break;
		}

		// Check for stagnation
		if(fabs(z2-lastz2) < 1.0e-6*lastz2) {
			times.setStagnate(numSystems);
			filePrint(stderr,"STAGNATION: Relative Primal Error Reached = "
				          "%e %e %e %e\n",sqrt(z2/pseudoFNormSq), fabs(z2-lastz2),
			          z2, lastz2);
			break;
		}
		lastz2 = z2;

		orthogonalize(z,y);
		singlePr(y, p, alpha);

	}

	// u = u + R beta
	if(isDynamic == 0)
		addRSingle(u, beta);

	// localSolveAndJump(f,lambda0,u,r);

	// Fp = lambda0;
	// computeL0(f, r, u, beta, Fp);
	//localSolveAndJump(f,Fp,beta,u,r);
	//addGs(r, w, beta);
	//if(isDynamic == 0)
	//  addRSingle(u, beta);

#ifdef DEBUG_MPC
	int numMPCs = mpcToSub->csize();
 if(numMPCs > 0) {
  alpha.zero();
  execParal(mpcToSub->csize(), this, &GenFetiSolver<Scalar>::getSQtMult, &u, alpha.data());
  execParal(mpcToSub->csize(), this, &GenFetiSolver<Scalar>::addMpcRhs, alpha.data());
   filePrint(stderr,"%d Is this zero enough?\n",numMPCs);
   int iMPC;
   for(iMPC=0; iMPC<numMPCs; ++iMPC)
     fprintf(stderr,"%d %e\n",iMPC+1,alpha[alpha.size()-numMPCs+iMPC]);
  }
#endif

	// make solution compatible u = u + deltaU
	makeCompatible( u, deltaU );

	// Store number of iterations, primal error and dual error
	setAndStoreInfo( iter+1, (errorEstimator/pseudoFNormSq), (w.sqNorm())/w0Norm2 );

	// end solve timer
	times.solve += getTime();
	times.memoryFETI  += threadManager->memoryUsed();
}


// ***************************************************************************
//  A small vocabulary note:
//     a corner has 1 and only 1 subdomain that is called its 'master'
//     in this implementation, the master of a corner is the subdomain
//     that touches this corner has the smallest global subdomain number
//
//     a local corner is a corner for which the subdomain being considered
//     is its master.
//     a neighbor corner is any other corner (some other sub is the master)
//
//  preProcessCorners() does the following corner preprocessing:
//      - Have each subdomain identify the corners for which they are 
//        the master
//      - Exchange corner information with neighbors to get the corners
//        shared with another subdomain and for which the neighbor is the
//        master
//      - Compute F*C on local and neighbor's corners
// ***************************************************************************
template<class Scalar>
void
GenFetiSolver<Scalar>::preProcessCorners()
{
	// Initialize numCRNs, crnDofLen and BClocal
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::initializeCRNs, sPat);

	// Get total number of corner modes from subdomains
	crns = collectIntGlobalSum();
#ifdef DISTRIBUTED
	fprintf(stderr,"CPU %d found %d corners\n",myCPU,crns);
 crns = fetiCom->globalSum(crns);
#endif
	times.numCRNs = crns;

	if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
		filePrint(stderr, " ... %5d corners were found       ...\n", crns);

	// Initialize neighbors number of corner modes
	sPat->exchange();
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbCRNs, sPat);

	// cset is constructed in GenFetiSolver<Scalar>::factorMatrices
	GenCoarseSet<Scalar> *cset = opControl->cset;
	//
	// KENDALL, this routine needs modification to run FETI-2 in Distributed!
	//
	paralApplyToAll(nsub, cset, &GenCoarseSet<Scalar>::buildBNeighbCs, cset);

	// compute F C

	// compute my own Fi^{-1} BCi
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::computeFiBC);

	// compute neighbors F C
	paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbFC);
}

template<class Scalar>
Scalar
GenFetiSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &ifrc, GenDistrVector<Scalar> &bf,
                                         GenVector<Scalar> &beta, GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &lambda) const
{
	startTimerMemory(times.sAndJ, times.memorySAndJ);

	// u = K+(f - B^t lambda - Q^t beta)
	timedParal(times.solveAndJump, nsub, this, &GenFetiSolver<Scalar>::localSolve2, &u,
	           &lambda, &beta, &ifrc, &bf);
	vPat->exchange();
	timedParal(times.solveAndJump, nsub, this, &GenFetiSolver<Scalar>::interfDiffAndDot,
	             bf, lambda);

	// Sum each subdomains dot product contribution
	Scalar ret = 0.0;
	int i;
	for(i = 0; i < nsub; ++i)
		ret += fetiOps[i]->res;

#ifdef DISTRIBUTED
	if(fetiCom) ret = fetiCom->globalSum(ret);
#endif

	stopTimerMemory(times.sAndJ, times.memorySAndJ);

	return ret;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::localSolve2(int iSub, GenDistrVector<Scalar> *v1, GenDistrVector<Scalar> *v2,
                                   GenVector<Scalar> *beta, GenDistrVector<Scalar> *v3, GenDistrVector<Scalar> *v4) const
{
	GenSubDomain<Scalar> *subd = sd[iSub];

	int sn = subd->localSubNum();

	Scalar *localvec  = v1->subData(sn);
	Scalar *interfvec = v2->subData(sn);
	Scalar *localsrc  = v3->subData(sn);
	Scalar *interfsrc = v4->subData(sn);

	int i;
	for(i = 0; i < subd->localLen(); ++i)
		localvec[i] = localsrc[i];

	for(i = 0; i < subd->interfLen(); ++i)
		interfvec[i] = interfsrc[i];

	Scalar *gamma = beta->data() + eqNums->firstdof(mOffset);

	subdomains[iSub]->fetiBaseOp(fetiOps[iSub]->solver, localvec, interfvec, gamma);

	subdomains[iSub]->sendInterf(interfvec, vPat);
}

template<class Scalar>
#ifdef BOOL_NOT_DEFINED
int
#else
bool
#endif
GenFetiSolver<Scalar>::isLowestLocalNeighbor(int subI, int subJ) const
{
	int kSub;
	// First find if this subdomain (subI) is the lowest numbered subdomain
	// connected to subJ in this cpu
	for(kSub = 0; kSub < subToSub->num(subJ); ++kSub) {
		int subK = (*subToSub)[subJ][kSub];
		if(glSubToLoc[subK] >= 0  && subK < subI) return false;
	}
	return true;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addNonLocalGContrib(int subI, int subJ)
{
	if(isLowestLocalNeighbor(subI,subJ) == false) return;

	int kSub, jSub;

	// Now we know we are the smallest one loop over all the subs connected
	// to subJ and that are in this CPU
	for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
		int subK = (*subToSub)[subJ][kSub];
		int locK = glSubToLoc[subK];
		if(locK < 0) continue;
		// Now look into K's neighbors, to find which one is subJ
		for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
			if(sd[locK]->scomm->subNums[jSub] == subJ) {
				if(singleCoarse)
					opControl->cset[locK].addGContrib(singleCoarse,
					                                  jSub+1, eqNums, gOffset);
				break;
			}
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addNonLocalCContrib(int subI, int subJ)
{
	if(isLowestLocalNeighbor(subI,subJ) == false) return;

	int cOffset = subToSub->csize();

	int kSub, jSub;

	// Now we know we are the smallest one loop over all the subs connected
	// to subJ and that are in this CPU
	for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
		int subK = (*subToSub)[subJ][kSub];
		int locK = glSubToLoc[subK];
		if(locK < 0) continue;
		// Now look into K's neighbors, to find which one is subJ
		for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
			if(sd[locK]->scomm->subNums[jSub] == subJ) {
				if(singleCoarse)
					opControl->cset[locK].addCContrib(singleCoarse,
					                                  jSub+1, eqNums, cOffset, gOffset);
				break;
			}
	}

}

template<class Scalar>
void
GenFetiSolver<Scalar>::getNonLocalSubAlphaGtQ(int subI, int subJ, Scalar *va,
                                              GenDistrVector<Scalar> *dv) const
{
	if(isLowestLocalNeighbor(subI,subJ) == false) return;

	int kSub, jSub;

	Scalar *alpha = va + eqNums->firstdof(subJ);
	for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
		int subK = (*subToSub)[subJ][kSub];
		int locK = glSubToLoc[subK];
		if(locK < 0) continue;
		Scalar *ovec = dv->subData(locK);
		// Now look into K's neighbors, to find which one is subJ
		for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
			if(sd[locK]->scomm->subNums[jSub] == subJ) {
				opControl->cset[locK].subAlphaGtQ(jSub, ovec, alpha);
				break;
			}
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getNonLocalFCtMult(int subI, int subJ, Scalar *va,
                                          GenDistrVector<Scalar> *dv) const
{
	if(isLowestLocalNeighbor(subI,subJ) == false) return;

	int kSub, jSub;

	Scalar *alpha = va + eqNums->firstdof(subJ);
	for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
		int subK = (*subToSub)[subJ][kSub];
		int locK = glSubToLoc[subK];
		if(locK < 0) continue;
		Scalar *ovec = dv->subData(locK);
		// Now look into K's neighbors, to find which one is subJ
		for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
			if(sd[locK]->scomm->subNums[jSub] == subJ) {
				opControl->cset[locK].getCtFMult(jSub, ovec, alpha);
				break;
			}
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::initGMRES(GenDistrVector<Scalar> &p)
{
	Scalar *hp = (Scalar*) dbg_alloca(halfSize*sizeof(Scalar));
	GenDistrVector<Scalar> *ptmp = 0;

	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::gatherHalfInterface, &p, ptmp, hp, hp);

	int i;
	double beta = 0.0;
	for(i=0;i<halfSize;i++) beta += ScalarTypes::sqNorm(hp[i]);
#ifdef DISTRIBUTED
	beta = fetiCom->globalSum(beta);
#endif
	beta = sqrt(beta);
	for(i=0;i<halfSize;i++) hp[i] /= beta;
	oSetGMRES->init(hp,beta);

	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::scatterHalfInterface, hp, &p);
	vPat->exchange();
	timedParal(times.orthogonalize, nsub, this,
	             &GenFetiSolver<Scalar>::rebuildInterface, p);
}

template<class Scalar>
double
GenFetiSolver<Scalar>::orthoAddGMRES(GenDistrVector<Scalar> &p,GenDistrVector<Scalar> &Fp)
{
	Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hFp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::gatherHalfInterface, &p, &Fp, hp, hFp);

	times.memoryOSet -= memoryUsed();
	double r2;
	r2 = oSetGMRES->ModorthoAddTimed(times.orthogonalize,hFp,hp); // JL

	times.memoryOSet += memoryUsed();

	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::scatterHalfInterface, hp, &p);
	vPat->exchange();
	timedParal(times.orthogonalize, nsub, this,
	             &GenFetiSolver<Scalar>::rebuildInterface, p);

	return r2;
}

template<class Scalar>
void
GenFetiSolver<Scalar>::GMRESSolution(GenDistrVector<Scalar> &p)
{
	Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

	oSetGMRES->solution(hp);

	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::scatterHalfInterface, hp, &p);
	vPat->exchange();
	timedParal(times.orthogonalize, nsub, this,
	             &GenFetiSolver<Scalar>::rebuildInterface, p);
}

template<class Scalar>
GenFetiSolver<Scalar>::~GenFetiSolver()
{
	if(edgeToSub) { delete edgeToSub; edgeToSub = 0; }
	if(subToEdge) { delete subToEdge; subToEdge = 0; }
	if(coarseConnect && (coarseConnect != subToSub)) { delete coarseConnect; coarseConnect = 0; }
	if(renum) { renum->clearMemory(); delete renum; renum = 0; }
	if(eqNums) { delete eqNums; eqNums = 0; }
	if(PFcNums) { delete PFcNums; PFcNums = 0; }
	if(singleCoarse) { delete singleCoarse; singleCoarse = 0; singleCoarseSolver = 0; }
	if(oSetCG) { delete oSetCG; oSetCG = 0; }
	if(oSetGMRES) { delete oSetGMRES; oSetGMRES = 0; }
	if(oSetGCR) { delete oSetGCR; oSetGCR = 0; }
	if(fetiOps) {
		for(int i=0; i<nsub; ++i)
			if(fetiOps[i]) { delete fetiOps[i]; fetiOps[i] = 0; }
		delete [] fetiOps; fetiOps = 0;
	}
	if(fetiTasks) { delete [] fetiTasks; fetiTasks = 0; }
	if(GtGsolver) { delete GtGsolver; GtGsolver = 0; }
	if(PCtFPC) { delete PCtFPC; PCtFPC = 0; }
	if(wksp) { delete wksp; wksp = 0; }
	if(GtQGs) { delete GtQGs; GtQGs = 0; }
	if(GtFCs) { delete GtFCs; GtFCs = 0; }
	if(gtqglocal) { delete [] gtqglocal; gtqglocal = 0; }
	if(vPat) { delete vPat; vPat = 0; }
	if(rbmPat) { delete rbmPat; rbmPat = 0; }
	if(sPat) { delete sPat; sPat = 0; }
	if(wiPat) { delete wiPat; wiPat = 0; }
	// don't delete: fetiInfo, subTosub, mpcToSub, glSubToLoc
}

template<class Scalar>
int
GenFetiSolver<Scalar>::predictGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0)
{
	if(oSetGCR->numDir() == 0) {
		return 0; // No available prediction
	}
	else if(!fetiInfo->useMRHS) {
		oSetGCR->reset();
		return 0;
	}
	else {
		Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hOp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		execParal(nsub, this, &GenFetiSolver<Scalar>::gatherHalfInterface, &r, &lambda0, hp, hOp);
		oSetGCR->predict(hp,hOp);

		execParal(nsub, this, &GenFetiSolver<Scalar>::scatterHalfInterface, hOp, &lambda0);
		vPat->exchange();
		execParal(nsub, this, &GenFetiSolver<Scalar>::rebuildInterface, lambda0);
		return 1;
	}
}

template<class Scalar>
void
GenFetiSolver<Scalar>::orthogonalizeGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fr,
                                        GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp)
{
	times.reOrtho -= getTime();
	if(oSetGCR->numDir() == 0) {
		p = r;
		Fp = Fr;
	}
	else {
		Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hpF = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		timedParal(times.orthogonalize, nsub, this, &GenFetiSolver<Scalar>::gatherHalfInterface,
		           &r, &Fr, hp, hpF);

		Scalar *hOp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hOpF = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		oSetGCR->orthogonalizeTimed(times.orthogonalize, hp, hpF, hOp, hOpF);

		timedParal(times.orthogonalize, nsub, this,
		           &GenFetiSolver<Scalar>::scatterHalfInterface, hOp, &p);
		vPat->exchange();
		timedParal(times.orthogonalize, nsub, this,
		             &GenFetiSolver<Scalar>::rebuildInterface, p);

		timedParal(times.orthogonalize, nsub, this,
		           &GenFetiSolver<Scalar>::scatterHalfInterface, hOpF, &Fp);
		vPat->exchange();
		timedParal(times.orthogonalize, nsub, this,
		             &GenFetiSolver<Scalar>::rebuildInterface, Fp);
	}
	times.reOrtho += getTime();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::orthoAddGCR(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar FpFp)
{
	Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hFp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	timedParal(times.orthogonalize, nsub, this,
	           &GenFetiSolver<Scalar>::gatherHalfInterface, &p, &Fp, hp, hFp);
	times.memoryOSet -= memoryUsed();
	oSetGCR->orthoAddTimed(times.orthogonalize, hp, hFp, FpFp);
	times.memoryOSet += memoryUsed();
}

template<class Scalar>
void
GenFetiSolver<Scalar>::makeRbmPat()
{
	// create rbmPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
	if(!rbmPat) {
		rbmPat = new FSCommPattern<Scalar>(fetiCom, cpuToSub, myCPU, FSCommPattern<Scalar>::CopyOnSend,
		                                   FSCommPattern<Scalar>::NonSym);
		for(int iSub=0; iSub<nsub; ++iSub) subdomains[iSub]->setRbmCommSize(fetiOps[iSub]->numRBM, rbmPat);
		rbmPat->finalize();
	}
}
template class GenFetiSolver<std::complex<double>>;
#include <Feti.d/NLFeti.C>
template class GenFetiSolver<double>;
template class GenFetiWorkSpace<double>;
template class GenFetiWorkSpace<std::complex<double>>;
