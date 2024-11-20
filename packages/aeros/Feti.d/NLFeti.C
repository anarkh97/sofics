#include <cmath>
#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/SubDomain.h>
#include <Feti.d/Feti.h>
#include <Threads.d/PHelper.h>
#include <Feti.d/OrthoSet.h>
#include <Utils.d/linkfc.h>
#include <Solvers.d/Rbm.h>
#include <Corotational.d/DistrGeomState.h>
#include <Math.d/matrix.h>
#include <Math.d/IntFullM.h>
#include <Math.d/SymFullMatrix.h>
#include <Math.d/BigMatrix.h>
#include <Timers.d/GetTime.h>

extern int verboseFlag;

template<class Scalar>
void
GenFetiSolver<Scalar>::reBuild(FullSquareMatrix **kel, DistrGeomState &geomState, int iter,
                               int step)
{
 if(iter == 0) { 
   this->fetiInfo->numLoadSteps++;
   // fprintf(stderr," Load step %d\n",this->fetiInfo->numLoadSteps);
 }

 // Rebuild FETI Solver (reassembles subdomain K and factors K)
 this->times.factor -= getTime();
 execParal(this->nsub, this, &GenFetiSolver<Scalar>::subdomainReBuild, kel, &geomState);
 this->times.factor += getTime();


 //fprintf(stderr," ... Time to Rebuild Subdomain Matrices and Factor %e\n",
  //       this->times.factor/1000.0);


 // MODIFICATION

 // Rebuild FETI preconditioner when requested
 // KHP: temperary change to check preconditioner effects
 // rebuild preconditioner until a certain newton iteration
 // if( iter < this->fetiInfo->nPrecond() ) {
 // if( (iter % this->fetiInfo->nPrecond()) == 0) {

 // ALWAYS rebuild preconditioner

 this->times.reBuildPrec -= getTime();
 execParal(this->nsub, this, &GenFetiSolver<Scalar>::reBuildMatrices, kel);
 this->times.reBuildPrec += getTime();

 //fprintf(stderr," ... Time to Rebuild Preconditioner and factor %16e %d %d\n",
  //              this->times.reBuildPrec/1000.0,step,iter);

 // Store information on when tangent stiffness and preconditioner were rebuilt
 this->times.iterations[numSystems].rebuildPrec = 1;

 // Rebuild GtG from the new rigid body modes
 if(this->fetiInfo->version == FetiInfo::feti1) {
   this->times.reBuildGtG -= getTime();
   if(numrbms)
     reBuildGtG();
   else
     GtGsolver = 0;
   this->times.reBuildGtG += getTime();
   //fprintf(stderr," ... Time to Rebuild GtG    (1st level coarse) %16e %d\n",
    //       this->times.reBuildGtG/1000.0, numrbms);
 } else if(this->fetiInfo->version == FetiInfo::feti2 && (crns > 0) ) { 
   // Rebuild second level coarse problem for FETI 2 case with corners
   this->times.reBuildPCtFPC -= getTime(); 
   reBuildPCtFPC();
   this->times.reBuildPCtFPC += getTime();
   //fprintf(stderr," ... Time to Rebuild PCtFPC (2nd level coarse) %16e %d\n",
    //       this->times.reBuildPCtFPC/1000.0, crns);
 } else {
   fprintf(stderr," *** ERROR: This FETI Version not supported in reBuild\n");
 }

 // Begin new ortho set for new linear system
 if(fetiInfo->outerloop == 0) { // CG only
   oSetCG->newOrthoSet();

   // NOTE: If we do not need to store previous orthoset,
   // we overwrite the orthosets to save memory space. We only need
   // previous krylov spaces when we are using krylov preconditioner

   this->times.iterations[numSystems].rebuildKrylov = 0;
   if(this->fetiInfo->nlPrecFlg == 0) {
     if(verboseFlag) filePrint(stderr," ... Resetting OrthoSet         ... \n");
     this->times.iterations[numSystems].rebuildKrylov = 1;
     oSetCG->reset();
   }

   // KHP: THIS IS FOR RESETTING KRYLOV AT EACH LOAD STEP
   if(this->fetiInfo->nlPrecFlg == 2 && (iter == 0) && (numSystems != 0) ) {
     if(verboseFlag) filePrint(stderr," ... Resetting OrthoSet at load step %d ... \n",
             this->fetiInfo->numLoadSteps);
     this->times.iterations[numSystems].rebuildKrylov = 1;
     oSetCG->reset();
   }

   // KHP: THIS IS FOR RESETTING KRYLOV IF THE LAST STEP STAGNATED.
   // if((iter == 0) && ( this->fetiInfo->numLoadSteps % 2 == 0 )) {
   // if(this->times.iterations[numSystems - 1].stagnated == 1) {
   /*
   if( numSystems % 8 == 0  && (numSystems != 0) ) {
     fprintf(stderr,"\n===> RESETTING ORTHO SET at load step %d\n\n",
           this->fetiInfo->numLoadSteps);
     this->times.iterations[numSystems].rebuildKrylov = 1;
     oSetCG->reset();
   }
   */
 }
 else this->resetOrthoSet();

}


template<class Scalar>
void
GenFetiSolver<Scalar>::subdomainReBuild(int iSub, FullSquareMatrix **kel, 
                                        DistrGeomState *gs)
{
  // Rebuild each Subdomains Geometric rigid body modes
  fetiOps[iSub]->solver->reBuildGeometricRbms( gs->getSubGeomState(iSub) );

  // Rebuild Subdomain stiffness matrices
  fetiOps[iSub]->solver->reBuild(kel[iSub]);

  // NOTE: K has already been factored within reBuild
  if(verboseFlag)
    fprintf(stderr," ... Subdomain %3d found %3d ZEMs   ...\n",
            this->sd[iSub]->subNum()+1, fetiOps[iSub]->K->numRBM());

  if(fetiOps[iSub]->rbm) {
    // Geometric rigid body mode method
    fetiOps[iSub]->numRBM = fetiOps[iSub]->rbm->numRBM();
    if(fetiOps[iSub]->numRBM > 0) {
      fetiOps[iSub]->rbm->getRBMs(fetiOps[iSub]->locRBMs.data());
    }
  } 
  else {
    // Tolerance rigid body mode method
    fetiOps[iSub]->numRBM = fetiOps[iSub]->K->numRBM();
    if(fetiOps[iSub]->numRBM > 0) {
       fetiOps[iSub]->K->getRBMs(fetiOps[iSub]->locRBMs.data());
    }
  }
}


template<class Scalar>
void
GenFetiSolver<Scalar>::reBuildMatrices(int iSub, FullSquareMatrix **kel)
{
  this->sd[iSub]->reBuildKbb(kel[iSub]);
}

// rebuild 2nd Level Coarse Problem

template<class Scalar>
void
GenFetiSolver<Scalar>::reBuildPCtFPC()
{
 this->times.coarse2 -= getTime();

 // int oldNumRbms = opControl->rbms;
 fprintf(stderr,"Begin build PCtFPC\n");

 // If the number of rbms is different, set changed to 1 otherwise 0
 // int changed = (numrbms == oldNumRbms) ? 0 : 1;

// fprintf(stderr," starting reBuildP %d %d %d\n",changed,numrbms, oldNumRbms);

 opControl->rbms = numrbms;

 if(numrbms > 0) {

   // Compute G'C

   paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::assembleGtCs);

   // reBuild Rgc = -(G^t G)^-1 G^t C
   // forward & backward substution to get -Rgc
   execParal(threadManager->numThr(), this, &GenFetiSolver<Scalar>::finishRgc,
             threadManager->numThr());

   // multiply by -1 to get Rgc
   int i;
   for(i=0; i<numrbms*crns; ++i)
     opControl->Rgc[i] *= -1.0;

   // re-allocate memory for GtFC if number of rbms has changed

   //if(changed) {
   GtFCs = new GenFullM<Scalar>(numrbms,crns);
   opControl->GtFC = GtFCs;
   //}
   opControl->GtFC->zero();

   // recompute C'FC, G'FC and G'FG

   makelocalFcoarse();

   // recompute new projected F coarse

   GenStackFullM<Scalar> wkbf(crns,crns,PCtFPC->data());

   addAllFcoarse(wkbf);

   this->times.coarse2 += getTime();

   this->times.pfactor2 -= getTime();
   PCtFPC->parallelFactor();
   this->times.pfactor2 += getTime();

 } 
 else {
   PCtFPC = 0;
   this->times.coarse2 += getTime();
 }
}


//  NOTE: Use this to track down memory being wasted.
//
//  fprintf(stderr," 1 Checking memory %f Mb\n",
//            threadManager->memoryUsed()/(1024.0*1024.0));

template<class Scalar>
void
GenFetiSolver<Scalar>::reBuildGtG()
{
 // the timing for rebuilding GtG is accumulated with
 // the regular 1st Level Coarse Problem timing
 this->times.coarse1 -= getTime();

 // Get number of each neighbours rbms
 paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::sendNumNeighbRBM, sPat);
 sPat->exchange();
 paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbRBM, sPat);

 // Get all numbers of rigid body modes and dispatch RBMs to neighbors
 makeRbmPat();
 execParal(this->nsub, this, &GenFetiSolver<Scalar>::reSendInterfaceRBM);
 rbmPat->exchange();

 // Get all Neighbours QGs
 // execParal(this->nsub, this, &GenFetiSolver<Scalar>::reGetNeighbQGs);

 // compute neighbours QGs
 paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbQGs, rbmPat);

 int iSub;
 if(isDynamic && isFeti2 && (crns > 0) ) {

     paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::computeFiBC);

     paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbFC);

     for(iSub = 0; iSub < this->nsub; ++iSub)
       eqNums->setWeight(iSub, fetiOps[iSub]->getNumRBM()+
                         fetiOps[iSub]->getcrnDofSize());
 } else
     for(iSub = 0; iSub < this->nsub; ++iSub)
       eqNums->setWeight(iSub, fetiOps[iSub]->getNumRBM());

 eqNums->makeOffset();

 // Determine if # of rigid body modes has changed from previous solve.
 // int oldNumRbms = numrbms;

 numrbms = eqNums->size();

 // If the number of rbms is different, set changed to 1 otherwise 0
 //int changed = (numrbms == oldNumRbms) ? 0 : 1;
 int changed = 1;
 if(changed) {

   // The following routine re-sets alpha offset for subdomains AND in 
   // case of Feti 2 dynamic it does betaOffset as well
   paralApplyToAll(this->nsub,fetiOps,&GenFetiOp<Scalar>::reSetAlphaOffsets,eqNums->allOffsets());

   GtGsolver = GenSolverFactory<Scalar>::getFactory()->createDistSolver(coarseConnect, eqNums, *fetiInfo->coarse_cntl, opControl->sparseGtG, this->fetiCom);
 } else
    opControl->sparseGtG->zeroAll();


 if (isFeti2 && (isDynamic == 0)) {
   //if(changed) {
     GtQGs = new GenSymFullMatrix<Scalar>(numrbms);
     opControl->GtQGs = GtQGs;
   //}
   opControl->GtQGs->zero();
 }


 // Now assemble GtG
 paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::assembleGtQGs);

 this->times.pfactor -= getTime();
 GtGsolver->parallelFactor();
 this->times.pfactor += getTime();

 this->times.coarse1 += getTime();
}


//template<class Scalar>
//void
//GenFetiSolver<Scalar>::reGetNeighbQGs(int iSub)
//{
// if(QGisLocal == 0) {
//   if(fetiOps[iSub]->rbm)
//     fetiOps[iSub]->control->cset[this->sd[iSub]->subNum()].reGetNeighbQGs(
//           fetiOps[iSub]->control->cset, this->sd[iSub], 0, rbmPat, fetiOps[iSub]->solver);
// } else {
//   if (isFeti2)
//     fetiOps[iSub]->control->cset[this->sd[iSub]->subNum()].reGetNeighbQGs(
//           fetiOps[iSub]->control->cset, this->sd[iSub], 0, rbmPat, fetiOps[iSub]->solver);
//   else
//     fetiOps[iSub]->control->cset[this->sd[iSub]->subNum()].reGetNeighbGs(
//                    fetiOps[iSub]->control->cset, rbmPat, this->sd[iSub]);
// }
//}

// [C G]^t F [C G]
// Very similar to  G^t Q G

template<class Scalar>
void
GenFetiSolver<Scalar>::reMakeLocalFcoarse()
{
 // compute F C
 // compute my own Fi^{-1} BCi

 execParal(this->nsub, this, &GenFetiSolver<Scalar>::reComputeFiBC);

 // compute neighbors F C
 // execParal(this->nsub, this, &GenFetiSolver<Scalar>::reGetNeighbFC);

 // compute neighbors F C
 paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbFC);

 // compute G'FG, C'FC, C'FG
 paralApplyToAll(this->nsub, fetiOps, &GenFetiOp<Scalar>::assembleGtFCs);

}

template<class Scalar>
void
GenFetiSolver<Scalar>::reGetNeighbFC(int iSub)
{
 fetiOps[iSub]->control->cset[this->sd[iSub]->subNum()].computeNeighbFBCs(
               this->sd[iSub],fetiOps[iSub]->solver);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::reComputeFiBC(int iSub)
{

// expand BClocal ( # of non-zero value alone interface this->times # of
// lambda ) to locInterfCRNs ( interfsize this->times # of lambda )
// reason :: to do the Q * locInterfCRNs, where Q is Fi

 int i, j;

 int crnDofSize = fetiOps[iSub]->crnDofSize;

 Scalar *locInterfCRNs =
         (Scalar *) dbg_alloca(sizeof(Scalar)*crnDofSize*this->sd[iSub]->interfLen());

 for (i=0; i< crnDofSize*this->sd[iSub]->interfLen(); i++) locInterfCRNs[i] = 0.0;

 IntFullM *BC = fetiOps[iSub]->BClocal;

 for (i=0; i < crnDofSize; i++) {
    int  offset = i*this->sd[iSub]->interfLen();
    locInterfCRNs[offset+(*BC)[0][i]] = 1;
 }

 GenCoarseSet<Scalar> &thisSet = fetiOps[iSub]->control->cset[this->sd[iSub]->subNum()];

 // Compute locFBCs, i.e. Fi C_i, where i = myNum
 for(i = 0; i < crnDofSize; ++i) {
   int nc;
   if((nc = (*BC)[1][i]) == i)
    this->sd[iSub]->multFi(fetiOps[iSub]->solver,
            locInterfCRNs+i*this->sd[iSub]->interfLen(),
            thisSet.locFBCs+i*this->sd[iSub]->interfLen());
   else { // if this corner mode is not unique
     int len = this->sd[iSub]->interfLen();
     for(j = 0; j < len; ++j)
       thisSet.locFBCs[i*len+j] = thisSet.locFBCs[nc*len+j];
   }
 }
}

