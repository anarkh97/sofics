//
// Created by Michel Lesoinne on 1/18/18.
//

#include <complex>
#include <Feti.d/FetiInfo.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include "FetiBaseClass.h"
#include "CGOrthoSet.h"
#include "FetiSub.h"
#include "FetiOp.h"
#include "GMRESOrthoSet.h"
#include "GCROrthoSet.h"

template<class Scalar>
int
FetiBaseClass<Scalar>::halfOffset(int iSub) const { return fetiOps[iSub]->halfOffset; }

template<class Scalar>
bool
FetiBaseClass<Scalar>::predict(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0) const
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

		execParal(nsub, this, &FetiBaseClass<Scalar>::gatherHalfInterface, &r, &lambda0, hp, hOp);
		oSetCG->predict(hp, hOp);
		execParal(nsub, this, &FetiBaseClass<Scalar>::scatterHalfInterface, hOp, &lambda0);
		vPat->exchange();
		execParal(nsub, this, &FetiBaseClass<Scalar>::rebuildInterface, lambda0);

		return true; // Prediction available
	}
}

template<class Scalar>
void
FetiBaseClass<Scalar>::orthoAddGCR(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar FpFp)
{
	Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hFp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::gatherHalfInterface, &p, &Fp, hp, hFp);
	times.memoryOSet -= memoryUsed();
	oSetGCR->orthoAddTimed(times.orthogonalize, hp, hFp, FpFp);
	times.memoryOSet += memoryUsed();
}

template<class Scalar>
void
FetiBaseClass<Scalar>::makeRbmPat()
{
	// create rbmPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
	if(!rbmPat) {
		rbmPat = new FSCommPattern<Scalar>(fetiCom, cpuToSub, myCPU, FSCommPattern<Scalar>::CopyOnSend,
		                                   FSCommPattern<Scalar>::NonSym);
		for(int iSub=0; iSub<nsub; ++iSub) subdomains[iSub]->setRbmCommSize(fetiOps[iSub]->numRBM, rbmPat);
		rbmPat->finalize();
	}
}

template<class Scalar>
void
FetiBaseClass<Scalar>::orthogonalize(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &p) const
{
	times.reOrtho -= getTime();

	if(fetiInfo->maxortho <= 0 || oSetCG->numDir() == 0) {
		p = r;
	}
	else {
		Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hOp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

		GenDistrVector<Scalar> *vec2 = 0;
		timedParal(times.orthogonalize, nsub, this, &FetiBaseClass<Scalar>::gatherHalfInterface,
		           &r, vec2, hp, hOp);
		oSetCG->orthogonalizeTimed(times.orthogonalize, hp, hOp, fetiInfo->complex_hermitian);

		timedParal(times.orthogonalize, nsub, this,
		           &FetiBaseClass<Scalar>::scatterHalfInterface, hOp, &p);
		vPat->exchange();
		timedParal(times.orthogonalize, nsub, this,
		           &FetiBaseClass<Scalar>::rebuildInterface, p);
	}

	times.reOrtho += getTime();
}

template<class Scalar>
void
FetiBaseClass<Scalar>::orthoAddCG(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp) const
{
	if(fetiInfo->maxortho <= 0) return;
#ifdef DISTRIBUTED
	dbg_alloca(0); // needed to clear the stack memory
#endif

	Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hFp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::gatherHalfInterface, &p, &Fp, hp, hFp);

	times.memoryOSet -= memoryUsed();
	oSetCG->orthoAddTimed(times.orthogonalize, hp, hFp, pFp);
	times.memoryOSet += memoryUsed();
}

template<class Scalar>
void
FetiBaseClass<Scalar>::resetOrthoSet()
{
	if(fetiInfo->maxortho <= 0) return;
	switch(fetiInfo->outerloop) {
		default:
		case FetiInfo::OuterloopType::CG:
		case FetiInfo::OuterloopType::CGAL:
			if(fetiInfo->nlPrecFlg) oSetCG->newOrthoSet();
			else oSetCG->reset();
			break;
		case FetiInfo::OuterloopType::GMRES:
			oSetGMRES->reInit();
			break;
		case FetiInfo::OuterloopType::GCR:
			oSetGCR->reset();
			break;
	}
}

template<class Scalar>
int
FetiBaseClass<Scalar>::nlPreCondition(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &pr) const
{
	if(fetiInfo->maxortho <= 0 || oSetCG->numOrthoSets() == 1) return 0; // no preconditioner available
	startTimerMemory(times.nlPreCond, times.memoryNlPreCond);

	Scalar *hp  = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hpr = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

	execParal(nsub,this,&FetiBaseClass<Scalar>::gatherHalfInterface,&r,&pr,hp,hpr);
	oSetCG->precondition(hp, hpr);
	execParal(nsub, this, &FetiBaseClass<Scalar>::scatterHalfInterface, hpr, &pr);
	vPat->exchange();
	execParal(nsub, this, &FetiBaseClass<Scalar>::rebuildInterface, pr);
	stopTimerMemory(times.nlPreCond, times.memoryNlPreCond);

	return 1; // preconditioner available
}


template<class Scalar>
void
FetiBaseClass<Scalar>::reSolve(GenDistrVector<Scalar> &u)
{
	// We are allocating memory here, this is bad!
	GenDistrVector<Scalar> f(internalDI);

	f = u;
	this->solve(f,u);
}
// Subdomain level functions for parallel execution

template<class Scalar>
void
FetiBaseClass<Scalar>::sendDeltaF(int iSub, GenDistrVector<Scalar> &deltaF) const
{
	// parallel implementation of sending each subdomain's deltaF
	subdomains[iSub]->sendDeltaF(deltaF.subData(iSub), vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::normDeltaF(int iSub, double *subDots, GenDistrVector<Scalar> *deltaF) const
{
	// parallel implementation of computing true norm of deltaF for each subdomain
	subDots[iSub] = subdomains[iSub]->collectAndDotDeltaF(deltaF->subData(iSub), vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::rebuildInterface(int iSub, GenDistrVector<Scalar> &v) const
{
	subdomains[iSub]->rebuildInterf(v.subData(subdomains[iSub]->localSubNum()), vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::gatherHalfInterface(int iSub, const GenDistrVector<Scalar> *v1, const GenDistrVector<Scalar> *v2,
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
FetiBaseClass<Scalar>::scatterHalfInterface(int iSub, Scalar *v1, GenDistrVector<Scalar>* v2) const
{
	Scalar *vec1 = v1 + halfOffset(iSub);
	Scalar *vec2 = v2->subData(subdomains[iSub]->localSubNum());

	subdomains[iSub]->scatterHalfInterf(vec1, vec2);
	subdomains[iSub]->sendInterf(vec2, vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::interfSend(int iSub, GenDistrVector<Scalar> &dv1) const
{
	Scalar *interfvec = dv1.subData(subdomains[iSub]->localSubNum());
	subdomains[iSub]->sendInterf(interfvec, vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::interfaceDiff(int iSub, GenDistrVector<Scalar> &v) const
{
	subdomains[iSub]->interfaceJump(v.subData(subdomains[iSub]->localSubNum()), vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::sendScale(int iSub)
{
	subdomains[iSub]->sendDiag(fetiOps[iSub]->KasSparse, vPat);
}


template<class Scalar>
void
FetiBaseClass<Scalar>::setAndStoreInfo(int iter, double finalPrimal2,
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


template<class Scalar>
void
FetiBaseClass<Scalar>::multKbb(int iSub, const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &interfvec,
                               GenDistrVector<Scalar> &deltaU, GenDistrVector<Scalar> &deltaF, bool &errorFlag) const
{
	const Scalar *v1        =         v.subData(subdomains[iSub]->localSubNum());
	Scalar *interfaceVector = interfvec.subData(subdomains[iSub]->localSubNum());
	Scalar *subDeltaU       =    deltaU.subData(subdomains[iSub]->localSubNum());
	Scalar *subDeltaF       = deltaF.subData(subdomains[iSub]->localSubNum());

	if((fetiInfo->version == FetiInfo::fetidp) && (glNumMpc > 0))
		subdomains[iSub]->multKbbMpc(v1, interfaceVector, subDeltaU, subDeltaF, errorFlag); // also supports coupled_dph
	else if(domain->solInfo().isCoupled) {
		subdomains[iSub]->multKbbCoupled(v1, interfaceVector, subDeltaF, errorFlag); // doesn't support dual mpc
	}
	else
		subdomains[iSub]->multKbb(v1, interfaceVector, subDeltaU, subDeltaF, errorFlag);

	subdomains[iSub]->sendInterf(interfaceVector, vPat);
}

template<class Scalar>
void
FetiBaseClass<Scalar>::distributeForce(GenDistrVector<Scalar> &f) const
{
//	execParal(nsub, this, &FetiBaseClass<Scalar>::fSend,  f);
	threadManager->callParal(nsub,
	                         [this, &f](int iSub) {
		                         subdomains[iSub]->fSend(f.subData(subdomains[iSub]->localSubNum()), vPat);
	                         });
	vPat->exchange();
//	execParal(nsub, this, &FetiBaseClass<Scalar>::fScale, f);
	threadManager->callParal(nsub,
	                         [this, &f](int iSub) {
		                         subdomains[iSub]->fScale(f.subData(subdomains[iSub]->localSubNum()), vPat);
	                         });
}

template<class Scalar>
void
FetiBaseClass<Scalar>::distributeForce(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw) const
{
//	execParal(nsub, this, &FetiBaseClass<Scalar>::fSendCoupled,  f, fw);
	threadManager->callParal(nsub,
	                         [this, &f, &fw](int iSub) {
		                         subdomains[iSub]->fSend(f.subData(subdomains[iSub]->localSubNum()), vPat,
		                                                 fw.subData(subdomains[iSub]->localSubNum()));
	                         });
	vPat->exchange();
//	execParal(nsub, this, &FetiBaseClass<Scalar>::fScaleCoupled, f, fw);
	threadManager->callParal(nsub,
	                         [this, &f, &fw](int iSub) {
		                         subdomains[iSub]->fSend(f.subData(subdomains[iSub]->localSubNum()), vPat,
		                                                 fw.subData(subdomains[iSub]->localSubNum()));
	                         });
}

template<class Scalar>
void
FetiBaseClass<Scalar>::fSplit(int iSub, GenDistrVector<Scalar> &force) const
{
	subdomains[iSub]->splitInterf(force.subData(subdomains[iSub]->localSubNum()));
}


template<class Scalar>
void
FetiBaseClass<Scalar>::initGMRES(GenDistrVector<Scalar> &p)
{
	Scalar *hp = (Scalar*) dbg_alloca(halfSize*sizeof(Scalar));
	GenDistrVector<Scalar> *ptmp = 0;

	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::gatherHalfInterface, &p, ptmp, hp, hp);

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
	           &FetiBaseClass<Scalar>::scatterHalfInterface, hp, &p);
	vPat->exchange();
	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::rebuildInterface, p);
}


template<class Scalar>
double
FetiBaseClass<Scalar>::orthoAddGMRES(GenDistrVector<Scalar> &p,GenDistrVector<Scalar> &Fp)
{
	Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	Scalar *hFp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::gatherHalfInterface, &p, &Fp, hp, hFp);

	times.memoryOSet -= memoryUsed();
	double r2;
	r2 = oSetGMRES->ModorthoAddTimed(times.orthogonalize,hFp,hp); // JL

	times.memoryOSet += memoryUsed();

	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::scatterHalfInterface, hp, &p);
	vPat->exchange();
	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::rebuildInterface, p);

	return r2;
}


template<class Scalar>
void
FetiBaseClass<Scalar>::GMRESSolution(GenDistrVector<Scalar> &p)
{
	Scalar *hp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));

	oSetGMRES->solution(hp);

	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::scatterHalfInterface, hp, &p);
	vPat->exchange();
	timedParal(times.orthogonalize, nsub, this,
	           &FetiBaseClass<Scalar>::rebuildInterface, p);
}

template<class Scalar>
int
FetiBaseClass<Scalar>::predictGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0)
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
		execParal(nsub, this, &FetiBaseClass<Scalar>::gatherHalfInterface, &r, &lambda0, hp, hOp);
		oSetGCR->predict(hp,hOp);

		execParal(nsub, this, &FetiBaseClass<Scalar>::scatterHalfInterface, hOp, &lambda0);
		vPat->exchange();
		execParal(nsub, this, &FetiBaseClass<Scalar>::rebuildInterface, lambda0);
		return 1;
	}
}


template<class Scalar>
void
FetiBaseClass<Scalar>::orthogonalizeGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fr,
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
		timedParal(times.orthogonalize, nsub, this, &FetiBaseClass<Scalar>::gatherHalfInterface,
		           &r, &Fr, hp, hpF);

		Scalar *hOp = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		Scalar *hOpF = (Scalar *) dbg_alloca(halfSize*sizeof(Scalar));
		oSetGCR->orthogonalizeTimed(times.orthogonalize, hp, hpF, hOp, hOpF);

		timedParal(times.orthogonalize, nsub, this,
		           &FetiBaseClass<Scalar>::scatterHalfInterface, hOp, &p);
		vPat->exchange();
		timedParal(times.orthogonalize, nsub, this,
		           &FetiBaseClass<Scalar>::rebuildInterface, p);

		timedParal(times.orthogonalize, nsub, this,
		           &FetiBaseClass<Scalar>::scatterHalfInterface, hOpF, &Fp);
		vPat->exchange();
		timedParal(times.orthogonalize, nsub, this,
		           &FetiBaseClass<Scalar>::rebuildInterface, Fp);
	}
	times.reOrtho += getTime();
}

template<class Scalar>
double
FetiBaseClass<Scalar>::preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag) const
{
	startTimerMemory(times.precond, times.memoryPrecond);

	// Change this back to arguments of preCondition
	GenDistrVector<Scalar> &deltaU  = wksp->ret_deltaU();
	GenDistrVector<Scalar> &deltaF  = wksp->ret_deltaF();
	double primalResidual = 0.0;

	// Compute preconditioner
	timedParal(times.preconditioner, nsub, this, &FetiBaseClass<Scalar>::multKbb, v,
	           Pv, deltaU, deltaF, errorFlag);
	vPat->exchange();
	timedParal(times.preconditioner,nsub, this, &FetiBaseClass<Scalar>::interfaceDiff, Pv);

	if(errorFlag) {
		// Send deltaF through subdomain buffer
		timedParal(times.preconditioner, nsub, this, &FetiBaseClass<Scalar>::sendDeltaF, deltaF);
		vPat->exchange();

		double *subDots = (double *) dbg_alloca(sizeof(double)*nsub);

		// Compute each partial true norm  
		timedParal(times.preconditioner, nsub, this, &FetiBaseClass<Scalar>::normDeltaF, subDots, &deltaF);

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

template class FetiBaseClass<double>;
template class FetiBaseClass<std::complex<double>>;