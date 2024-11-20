#include <cmath>
#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <cstdlib>
#include <climits> //--- UH

#include <Driver.d/SubDomain.h>
#include <Feti.d/Feti.h>
#include <Feti.d/OrthoSet.h>
#include <Utils.d/linkfc.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/matrix.h>
#include <Math.d/IntFullM.h>

#include <Feti.d/FetiOpControler.h>
#include <Feti.d/CoarseSet.h>
#include "FetiOp.h"

template<class Scalar> 
GenFetiOp<Scalar>::GenFetiOp(GenSubDomain<Scalar>* lsd, GenFetiOpControler<Scalar>* fopc, 
                             int _isFeti2, int _isDynam, FSCommPattern<Scalar> *_vPat, Rbm *_rbm) 
{
 control    = fopc;
 sd         = lsd;	// local subdomain
 numNeighb  = sd->getSComm()->numNeighb;
 QGisLocal  = !control->nonLocalQ ;
 rbm        = _rbm;
 isFeti2    = _isFeti2;
 isDynamic  = _isDynam;
 crnDofSize = 0;
 vPat = _vPat;
}

template<class Scalar>
GenFetiOp<Scalar>::~GenFetiOp()
{
 if(interfBuff) { delete [] interfBuff; interfBuff = 0; }
 // don't delete control, rbm, sd, K, KasSparse, BClocal, solver, vPat
}

template<class Scalar> 
void
GenFetiOp<Scalar>::clean_up()
{
	if(K) K->clean_up();
	locRBMs.clear();
	neighbNumRBMs.clear();

	locInterfRBMs.clear();
	if(interfBuff) {
		delete [] interfBuff;
		interfBuff=0;
	}
	if(control->cset) control->cset->clean_up();
}

template<class Scalar> 
void
GenFetiOp<Scalar>::run()
{
 (this->*(control->operation))();
}

template<class Scalar> 
void
GenFetiOp<Scalar>::localSolve()
{
 Scalar *localvec  = control->dv1->subData(sd->localSubNum());
 Scalar *interfvec = control->dv2->subData(sd->localSubNum());
 Scalar *localsrc  = control->dv3->subData(sd->localSubNum());
 Scalar *interfsrc = control->dv4->subData(sd->localSubNum());

 int lLen = sd->localLen();	// local length
 int iLen = sd->interfLen();	// interface length

 int i;
 for(i = 0; i < lLen; ++i)
    localvec[i] = localsrc[i];
 for(i = 0; i < iLen; ++i)
    interfvec[i] = interfsrc[i];
 
 sd->fetiBaseOp(solver, localvec, interfvec);
 sd->sendInterf(interfvec, vPat);
}

template<class Scalar>
void
GenFetiOp<Scalar>::sendNumNeighbRBM(FSCommPattern<int> *sPat)
{
 // send Number of RBMs for each neighbor
 int iSub;
 for(iSub = 0; iSub < numNeighb; ++iSub) {
   FSSubRecInfo<int> sInfo = sPat->getSendBuffer(sd->subNum(), sd->getSComm()->subNums[iSub]);
   sInfo.data[0] = numRBM;
 }
}

template<class Scalar>
void
GenFetiOp<Scalar>::getNumNeighbRBM(FSCommPattern<int> *sPat)
{
 neighbNumRBMs.resize(numNeighb);
 // get Number of RBMs for each neighbor
 int iSub;
 for(iSub = 0; iSub < numNeighb; ++iSub) {
   FSSubRecInfo<int> rInfo = sPat->recData(sd->getSComm()->subNums[iSub], sd->subNum());
   neighbNumRBMs[iSub] = rInfo.data[0];
 }
 control->cset[sd->localSubNum()].neighbNumRBMs = neighbNumRBMs.data();
}

template<class Scalar> 
void
GenFetiOp<Scalar>::getNeighbQGs(FSCommPattern<Scalar> *rbmPat)
{
 //int kc = (control->nQ == 3) ? 1 : 0;
 int kc = control->nQ;
 if(QGisLocal == 0) {
   if(rbm)
     control->cset[sd->localSubNum()].getNeighbQGs(control->cset, sd, 0, rbmPat, solver);
   else
     control->cset[sd->localSubNum()].getNeighbQGs(control->cset, sd, kc, rbmPat);
 } else {
   if (isFeti2)
     control->cset[sd->localSubNum()].getNeighbQGs(control->cset, sd, 0, rbmPat, solver);
   else
     control->cset[sd->localSubNum()].getNeighbGs(control->cset, sd, rbmPat);
 }
}

template<class Scalar> 
void
GenFetiOp<Scalar>::getGtMult()
{
 Scalar *localvec = control->dv1->subData(sd->localSubNum());
 Scalar *alpha = control->vec2 + alphaOffset[0];
 control->cset[sd->localSubNum()].getGtMult(localvec,alpha);
}

template<class Scalar> 
void
GenFetiOp<Scalar>::getGtQMult()
{
 Scalar *localvec = control->dv1->subData(sd->localSubNum());
 Scalar *alpha    = control->vec2 + alphaOffset[0];

 control->cset[sd->localSubNum()].getGtQMult(localvec,alpha);

 if(QGisLocal == 0) {
   int myNum = sd->localSubNum();
   int iSub;
   for(iSub = 0; iSub < control->cset[myNum].numNeighb; ++iSub) {
      int subI = control->neighbSubId(myNum, iSub);
      int myID = control->index(subI, myNum);
      Scalar *ovec = control->dv1->subData(subI);
      control->cset[subI].getGtQMult(myID,ovec,alpha);
   }
 }
}

template<class Scalar> 
void
GenFetiOp<Scalar>::getGtFMult()
{
 Scalar *localvec = control->dv1->subData(sd->localSubNum());
 Scalar *alpha = control->vec2 + alphaOffset[0];

 control->cset[sd->localSubNum()].getGtFMult(localvec,alpha);
 int myNum = sd->localSubNum();
 int iSub;
 for(iSub = 0; iSub < control->cset[myNum].numNeighb; ++iSub) {
    int subI = control->neighbSubId(myNum, iSub);
    int myID = control->index(subI, myNum);
    Scalar *ovec = control->dv1->subData(subI);
    control->cset[subI].getGtFMult(myID,ovec,alpha);
 }
}

template<class Scalar> 
void
GenFetiOp<Scalar>::reSetAlphaOffsets(int *v)
{
	SComm *scomm =  sd->getSComm();

	alphaOffset[0] = v[sd->subNum()];

	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
		alphaOffset[iSub+1] = v[scomm->subNums[iSub]];

	// ML if in Feti2 dynamic, setup the Nus
	if(isFeti2 && isDynamic) {

		betaOffset[0] = v[sd->subNum()]+numRBM;

		for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
			betaOffset[iSub+1] = v[scomm->subNums[iSub]]+neighbNumRBMs[iSub];
	}
}

template<class Scalar> 
void
GenFetiOp<Scalar>::setAlphaOffsets(int *v)
{
	SComm *scomm =  sd->getSComm();
	alphaOffset.resize(scomm->numNeighb+1);
	alphaOffset[0] = v[sd->subNum()];

	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
		alphaOffset[iSub+1] = v[scomm->subNums[iSub]];

// ML if in Feti2 dynamic, setup the Nus
	if(isFeti2 && isDynamic) {
		betaOffset.resize(scomm->numNeighb+1);
		betaOffset[0] = v[sd->subNum()]+numRBM;
		for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
			betaOffset[iSub+1] = v[scomm->subNums[iSub]]+neighbNumRBMs[iSub];
		}
	}
}

template<class Scalar> 
void
GenFetiOp<Scalar>::subAlphaG()
{
 int iSub;
 int myNum = sd->localSubNum();
 Scalar *localvec = control->dv1->subData(myNum);
 Scalar *alpha    = control->vec2 + alphaOffset[0];

 control->cset[myNum].subLocAlphaG(localvec, alpha);

 for(iSub = 0; iSub < control->numNeighbSubs(myNum); ++iSub) {
    alpha = control->vec2 + alphaOffset[iSub+1];
#ifdef DISTRIBUTED
  control->cset[myNum].subAlphaNeighbG(iSub, localvec,
           control->cset[myNum].neighbGs[iSub], 
           control->cset[myNum].leadingDimGs[iSub], alpha);
#else
    int nghb = control->neighbSubId(myNum, iSub);
    control->cset[myNum].subAlphaNeighbG(iSub, localvec,
           control->cset[nghb].subG(myNum), control->cset[nghb].gSize, alpha);
#endif
 }
}

//  dv1 = dv1 - C * vec2

template<class Scalar> 
void
GenFetiOp<Scalar>::subNuC()
{
 Scalar *localvec = control->dv1->subData(sd->localSubNum());
 Scalar *nu = control->vec2 + betaOffset[0];

 int myNum = sd->localSubNum();

 control->cset[myNum].subLocNuC(localvec, nu);

 control->cset[myNum].subNuNeighbC(localvec, control->vec2,betaOffset.data()+1);
}

template<class Scalar> 
void
GenFetiOp<Scalar>::subAlphaGQ()
{
 int iSub;
 Scalar *localvec = control->dv1->subData(sd->localSubNum());

 int myNum = sd->localSubNum();
 Scalar *alpha = control->vec2 + alphaOffset[0];
 control->cset[myNum].subLocAlphaQG(localvec, alpha);
 if(QGisLocal) {
   for(iSub = 0; iSub < control->numNeighbSubs(myNum); ++iSub) {
      alpha = control->vec2 + alphaOffset[iSub+1];

#ifdef DISTRIBUTED
    // DISTRIBUTED WILL BE DONE ONLY FOR Q=I at this point XXX
    control->cset[myNum].subAlphaNeighbG(iSub, localvec,
      control->cset[myNum].neighbGs[iSub], control->cset[myNum].leadingDimGs[iSub], alpha);
#else

      int nghb = control->neighbSubId(myNum, iSub);
      control->cset[myNum].subAlphaNeighbG(iSub, localvec,   // In here G and Q
                                                             // G are dealt with the
                                                             // same way
      control->cset[nghb].subQG(myNum), control->cset[nghb].gSize, alpha);
#endif
   }
 } else {
   for(iSub = 0; iSub < control->numNeighbSubs(myNum); ++iSub) {
      alpha = control->vec2 + alphaOffset[iSub+1];
      control->cset[myNum].subAlphaNeighbQG(iSub, localvec, alpha);
   }
   sd->sendInterf(localvec, vPat);
 }
}

template<class Scalar> 
void
GenFetiOp<Scalar>::assembleGtQGs()
{
 int myNum = sd->localSubNum();
 int hasC = (isDynamic && isFeti2) ? 1 : 0;

 control->cset[myNum].addQContrib(control->sparseGtG, 0, control->eqNums, hasC);

 int iSub;
 if (!QGisLocal)
   for(iSub = 0; iSub <  control->numNeighbSubs(myNum); ++iSub)
       control->cset[control->neighbSubId(myNum, iSub)].addQContrib(
                     control->sparseGtG,
                     sd->scomm->remoteId[iSub]+1, control->eqNums, hasC);

 // FETI-2 and DYNAMICS Cases
 if(isFeti2 && isDynamic == 0) {
  control->cset[myNum].addFContrib(control->GtQGs, 0, control->eqNums);
  for(iSub = 0; iSub <  control->numNeighbSubs(myNum); ++iSub)
    control->cset[control->neighbSubId(myNum, iSub)].addFContrib(control->GtQGs,
                  sd->scomm->remoteId[iSub]+1, control->eqNums);
  }
}

template<class Scalar> 
void
GenFetiOp<Scalar>::makeCoarseSet()
{
 GenCoarseSet<Scalar> cset;
 SComm *scom    = sd->getSComm();
 cset.myNum     = sd->subNum();
 cset.numNeighb = scom->numNeighb;
 cset.neighbs   = scom->subNums.data();
 cset.subOffset = new int[cset.numNeighb+1];
 cset.sd        = sd;
 cset.numBCs    = 0;
 cset.neighbCSubIndex = 0;
 cset.isDynamic = isDynamic;

 // KHP: zero pointers so that I can write clean_up()
 cset.locGs     = 0;
 cset.locQGs    = 0;
 cset.neighbNumRBMs=0;
 cset.neighbQGs=0;
 cset.locFGs=0;
 cset.locBCs=0;
 cset.neighbFBCs2=0;
 cset.neighbNumCRNs=0;
 cset.neighbGs=0;
 cset.leadingDimGs=0;
 cset.locInterfCRNs=0;
 cset.neighbCSubIndex=0;
 cset.pool1=0;
 cset.pool2=0;
 

 int i;
/*
 for(i = 0; i < cset.numNeighb; ++i)
   cset.subOffset[i] = scom->sharedDOFs->offset(i);
 cset.subOffset[i] = scom->sharedDOFs->numConnect();
 cset.gSize = scom->sharedDOFs->numConnect();
*/
 for(i = 0; i < cset.numNeighb; ++i)
   cset.subOffset[i] = scom->offsetT(SComm::std,i);;
 cset.subOffset[i] = scom->lenT(SComm::std);
 cset.gSize = scom->lenT(SComm::std);

 control->cset[sd->localSubNum()] = cset;

}

/************************************************************************/
/**  Feti 2 related routines start here  ********************************/
/************************************************************************/

// Initialize numCRN and C matrix
// numCRN = number of corner modes

template<class Scalar> 
void
GenFetiOp<Scalar>::initializeCRNs(FSCommPattern<int> *sPat)
{
// BClocal stores C = BjCj. Mathematically, the size should be
// interfSize * numCRN.
// NOTE: This routine implemented Cross points and FETI2-ACD only!!

 GenCoarseSet<Scalar> &thisSet = control->cset[sd->localSubNum()];

 BClocal = sd->getC(crnDofSize, sPat);

 thisSet.numBCs = crnDofSize; 
 thisSet.locBCs = BClocal;

 setglobalSum(&crnDofSize);
}

template<class Scalar> 
void
GenFetiOp<Scalar>::assembleGtCs()
{
 int i,j;

 int myNum = sd->subNum();

 // For all subdomains connected to me

 // myRow is for the G part
 int rbms  = control->rbms;
 int crns  = control->crns;

  GenStackFullM<Scalar> CtG(crns, rbms, control->Rgc);
  int myRow = control->PFcNums->firstdof(myNum);
  int myNumofCrns = control->PFcNums->weight(myNum);

  // zero CtG in place of Rgc rowwise
  for (j=0; j < myNumofCrns; ++j)
    for(i = 0; i < rbms; ++i)
       CtG[j+myRow][i] = 0.0;

  control->cset[myNum].addCtGContrib(&CtG, control->PFcNums, control->eqNums);
}

template<class Scalar>
void
GenFetiOp<Scalar>::getNumNeighbCRNs(FSCommPattern<int> *sPat)
{
 int *neighbNumCRNs = new int[numNeighb];

 // get Number of Corners for each neighbor
 int iSub;
 for(iSub = 0; iSub < numNeighb; ++iSub) {
   FSSubRecInfo<int> rInfo = sPat->recData(sd->getSComm()->subNums[iSub], sd->subNum());
   neighbNumCRNs[iSub] = rInfo.data[0];
 }
 control->cset[sd->localSubNum()].neighbNumCRNs = neighbNumCRNs;
}

template<class Scalar> 
void
GenFetiOp<Scalar>::setglobalSum(void *data)
{
 control->globalSumBuffer[sd->localSubNum()] = data;
}

template<class Scalar> 
void
GenFetiOp<Scalar>::computeFiBC()
{
// expand BClocal ( # of non-zero value alone interface times # of
// lambda ) to locInterfCRNs ( interfsize times # of lambda )
// reason :: to do the Q * locInterfCRNs, where Q is Fi

 int i, j;

 Scalar *locInterfCRNs = 
         (Scalar *) dbg_alloca(sizeof(Scalar)*crnDofSize*sd->interfLen());

 for (i=0; i< crnDofSize*sd->interfLen(); i++) locInterfCRNs[i] = 0.0;

 for (i=0; i < crnDofSize; i++) {
    int  offset = i*sd->interfLen();
    locInterfCRNs[offset+(*BClocal)[0][i]] = 1;
 }

 GenCoarseSet<Scalar> &thisSet = control->cset[sd->localSubNum()];

 thisSet.locFBCs = new Scalar[crnDofSize*sd->interfLen()];

 // Compute locFBCs, i.e. Fi C_i, where i = myNum
 for(i = 0; i < crnDofSize; ++i) {
   int nc;
   if((nc = (*BClocal)[1][i]) == i)
     sd->multFi(solver, locInterfCRNs+i*sd->interfLen(),
            thisSet.locFBCs+i*sd->interfLen());
   else { // if this corner mode is not unique
     int len = sd->interfLen();
     for(j = 0; j < len; ++j)
       thisSet.locFBCs[i*len+j] = thisSet.locFBCs[nc*len+j];
   }

 }

}

template<class Scalar> 
void
GenFetiOp<Scalar>::getNeighbFC()
{
 control->cset[sd->localSubNum()].computeNeighbFBCs(sd, solver);
}

template<class Scalar> 
void
GenFetiOp<Scalar>::assembleGtFCs()
{

 int myNum = sd->subNum();
 control->cset[myNum].addGtFCContrib(control->GtFC, 0, control->PFcNums,
                                          control->eqNums);
 int iSub;
 for(iSub = 0; iSub <  control->numNeighbSubs(myNum); ++iSub)
      control->cset[control->neighbSubId(myNum, iSub)]
                .addGtFCContrib(control->GtFC,
                      sd->scomm->remoteId[iSub]+1, control->PFcNums,
                          control->eqNums);
}

template<class Scalar> 
void
GenFetiOp<Scalar>::assembleCtFCs()
{

 int myNum = sd->subNum();

  control->cset[myNum].addCtFCContrib(control->CtFCs, 0, control->PFcNums);
  int iSub;
  for(iSub = 0; iSub <  control->numNeighbSubs(myNum); ++iSub)
      control->cset[control->neighbSubId(myNum, iSub)]
                .addCtFCContrib(control->CtFCs,
                      sd->scomm->remoteId[iSub]+1, control->PFcNums);
}

template<class Scalar>
void
GenFetiOp<Scalar>::setBetaOffsets(int *v)
{
	SComm *scomm =  sd->getSComm();
	betaOffset.resize(scomm->numNeighb+1);
	betaOffset[0] = v[sd->subNum()];
	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		betaOffset[iSub+1] = v[scomm->subNums[iSub]];
	}
}

template<class Scalar> 
void
GenFetiOp<Scalar>::getCtFMult()
{
 Scalar *localvec = control->dv1->subData(sd->subNum());
 Scalar *beta = control->vec2 + betaOffset[0];

 control->cset[sd->subNum()].getCtFMult(localvec,beta);

 int myNum = sd->subNum();
 int iSub;
 for(iSub = 0; iSub < control->cset[myNum].numNeighb; ++iSub) {
   int subI = control->neighbSubId(myNum, iSub);
   int myID = control->index(subI, myNum);
   Scalar *ovec = control->dv1->subData(subI);
   control->cset[subI].getCtFMult(myID,ovec,beta);
 }
}

template<class Scalar> 
void
GenFetiOp<Scalar>::getCtMult()
{
 Scalar *localvec = control->dv1->subData(sd->subNum());
 Scalar *beta     = control->vec2 + betaOffset[0];
 control->cset[sd->subNum()].getCtMult(localvec,beta);
}

// vec = vec - FC nu
template<class Scalar>
void
GenFetiOp<Scalar>::subNuFC()
{
	int myNum = sd->subNum();

	Scalar *beta     = control->vec2 + betaOffset[0];
	Scalar *localvec = control->dv1->subData(sd->subNum());

	control->cset[sd->subNum()].subLocNuFC(localvec, beta);

	control->cset[myNum].subNuNeighbFC(localvec, control->vec2, betaOffset.data()+1);
	// sd->sendInterf(localvec, interfBuff);
	sd->sendInterf(localvec, vPat);
}

// vec = vec - FG alpha
template<class Scalar> 
void
GenFetiOp<Scalar>::subAlphaFG()
{

 int iSub;
 int myNum = sd->subNum();

 Scalar *localvec = control->dv1->subData(sd->subNum());
 Scalar *alpha    = control->vec2 + alphaOffset[0];

 control->cset[sd->subNum()].subLocAlphaFG(localvec, alpha);

 for(iSub = 0; iSub < control->numNeighbSubs(myNum); ++iSub) {
    alpha = control->vec2 + alphaOffset[iSub+1];
    control->cset[myNum].subAlphaNeighbFG(iSub, localvec,alpha);
 }
 // sd->sendInterf(localvec, interfBuff);
 sd->sendInterf(localvec, vPat);
}

template
class GenFetiOp<double>;
template
class GenFetiOp<std::complex<double>>;