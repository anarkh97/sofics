#include <cmath>
#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <sys/types.h>
#ifndef WINDOWS
#include <sys/mman.h>
#endif

#include <Driver.d/SubDomain.h>
#include <Threads.d/PHelper.h>
#include <Feti.d/CGOrthoSet.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/Skyline.d/DistSky.h>
#include <Math.d/Skyline.d/BlockSky.h>
#include <Math.d/Skyline.d/DistBlockSky.h>
#include <Math.d/DistBLKSparse.h>
#include <Utils.d/linkfc.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/matrix.h>
#include <Math.d/SymFullMatrix.h>
#include <Math.d/BLKSparseMatrix.h>

#include <Solvers.d/DSCsolver.h>

#include <Timers.d/GetTime.h>

extern int isParal;
extern FILE *debugFile;

template<class Scalar>
void
GenFetiSolver<Scalar>::makeDistGtG(int *glToLoc)
{
 // for correct timing, perform a processor synchronization
 // fprintf(stderr," ... Perform CPU Synchronization    ...\n");
 fetiCom->sync();

 startTimerMemory(times.coarse1, times.memoryGtG);

 if(isDynamic == 1 && fetiInfo->noCoarse == 1) {
   fprintf(stderr," ... No Coarse Problem used         ...\n");
   GtGsolver = 0;
   times.coarse1 += getTime();
   return;
 }

 // Get number of each neighbours rbms
 paralApply(nsub, fetiOps, &GenFetiOp<Scalar>::sendNumNeighbRBM, sPat);
 sPat->exchange();
 paralApply(nsub, fetiOps, &GenFetiOp<Scalar>::getNumNeighbRBM, sPat);

 // Get all the numbers of rigid body modes and dispatch RBMs to neighbors
 makeRbmPat();
 paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::sendInterfRBM, rbmPat);
 rbmPat->exchange();

 paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::getNeighbQGs, rbmPat);

 // Construct Coarse Problem Connectivity if it hasn't already been done
 if(coarseConnect == 0) {
   if(QGisLocal)
     coarseConnect = subToSub;
   else
     coarseConnect = subToSub->transcon(subToSub);
 }

 int glNumSub = subToSub->csize();

 // select renumbering method, if you are using Padma's distributed
 // parallel sparse solver, we do not renumber, it is taken care of
 // by Padma's DSC package.
 int renumMethod=0;
 if(fetiInfo->coarse_cntl->subtype == FetiInfo::skyline ||
    fetiInfo->coarse_cntl->subtype == FetiInfo::sparse)
    renumMethod=1;

 // Construct Coarse Problem renumbering if it hasn't already been done
 renumber = coarseConnect->renumByComponent(renumMethod);
 // Construct Coarse Problems equation numbers if it hasn't already been done
 if(eqNums == 0) eqNums = new SimpleNumberer(glNumSub, renumber.renum);

 int iSub;

 int *gRbmSize = new int[glNumSub];
 for(iSub = 0; iSub < glNumSub; ++iSub)
   gRbmSize[iSub] = 0;
 for(iSub = 0; iSub < nsub; ++iSub) {
   gRbmSize[ sd[iSub]->subNum() ] = opControl->cset[iSub].numGs;
 }

 // Now global sum of gRbmSize
 fetiCom->globalSum(glNumSub, gRbmSize);

 for(iSub = 0; iSub < glNumSub; ++iSub) {
    eqNums->setWeight(iSub, gRbmSize[iSub]);
 }

 eqNums->makeOffset();

 numrbms = eqNums->size();

 times.numRBMs = numrbms; // number of equations in GtQG

 // The following sets the alpha offset for the subdomains AND in the
 // case of Feti2 dynamic it does betaOffset as well

 paralApplyToAll(nsub, fetiOps, &GenFetiOp<Scalar>::setAlphaOffsets, eqNums->allOffsets());

 // Now assemble GtG

 if(numrbms > 0 ) {

   // This could be put into a switch/case structure!
   double tolerance = fetiInfo->coarse_cntl->trbm;
   GenBLKSparseMatrix<Scalar> *BLKMatrix = 0;
   GenSkyMatrix<Scalar> *sky = 0;
   DSCsolver *dscSolver = 0;
   GenBlockSky<Scalar> *BLKsky = 0;
   if(fetiInfo->coarse_cntl->subtype == FetiInfo::skyline) {
     //filePrint(stderr," Selected Skyline for 1st level Coarse\n");
     times.memoryGtGsky -= threadManager->memoryUsed();
     if(subToSub->csize() == numCPUs) {
       int firstAlpha = eqNums->firstdof(myCPU);
       int nRow       = eqNums->weight(myCPU);
       sky = new GenDistSky<Scalar>(coarseConnect, eqNums, tolerance, firstAlpha, nRow);
       times.memoryGtGDelete = 8*sky->size();
       // filePrint(stderr," Deleted GtG Memory: %14.5f Mb\n",
       //         times.memoryGtGDelete/(1024.0*1024.0));
     } 
     else {
       sky = new GenSkyMatrix<Scalar>(coarseConnect, eqNums, tolerance);
     }
     opControl->sparseGtG = sky;
     times.memoryGtGsky += threadManager->memoryUsed();
   } else if(fetiInfo->coarse_cntl->subtype == FetiInfo::blocksky) {
       if(subToSub->csize() == numCPUs) {
         int firstAlpha = eqNums->firstdof(myCPU);
         int nRow       = eqNums->weight(myCPU);
         BLKsky = new GenDistBlockSky<Scalar>(coarseConnect, eqNums, tolerance, firstAlpha, nRow);
       } else {
         BLKsky = new GenBlockSky<Scalar>(coarseConnect, eqNums, tolerance );
       }
         times.memoryGtGsky += threadManager->memoryUsed();
         opControl->sparseGtG = BLKsky;
   } else if(fetiInfo->coarse_cntl->subtype == FetiInfo::sparse) {
     times.memoryGtGsky -= threadManager->memoryUsed();
     if(subToSub->csize() == numCPUs) {
       fprintf(stderr," Selected Esmond BLK Sparse for 1st level Coarse\n");
       int firstAlpha = eqNums->firstdof(myCPU);
       int nRow       = eqNums->weight(myCPU);
       BLKMatrix = new GenDistBLKSparse<Scalar>(coarseConnect, eqNums, tolerance,
                                                *fetiInfo->coarse_cntl, firstAlpha, nRow);
       times.memoryGtGDelete = 8*BLKMatrix->size();
       //filePrint(stderr," Deleted GtG Memory: %14.5f Mb\n",
       //          times.memoryGtGDelete/(1024.0*1024.0));
     } else
       BLKMatrix = new GenBLKSparseMatrix<Scalar>(coarseConnect, eqNums, tolerance, *fetiInfo->coarse_cntl);
     times.memoryGtGsky += threadManager->memoryUsed();
     opControl->sparseGtG = BLKMatrix;
   } else {
#ifdef NO_COMPLEX
     filePrint(stderr," Selected Padma Sparse for 1st level Coarse\n");
     times.memoryGtGsky -= threadManager->memoryUsed();
     int scheme_number = 1;
     dscSolver = new DSCsolver(coarseConnect, eqNums, scheme_number);
     times.memoryGtGsky += threadManager->memoryUsed();
     opControl->sparseGtG = dscSolver;
#endif
   }

   opControl->eqNums = eqNums;

   if(isFeti2 && isDynamic == 0) {
     GtQGs = new GenSymFullMatrix<Scalar>(numrbms);
     GtQGs->zero();
     opControl->GtQGs = GtQGs;
   }

   execParal(nsub, this, &GenFetiSolver<Scalar>::assembleDistGtQGs, glToLoc);

   // Now unify sky
   double t1 = getTime();

   // Unify the matrix (global sum)
   if(fetiInfo->coarse_cntl->subtype == FetiInfo::skyline) {
     sky->unify(fetiCom);
   } else if (fetiInfo->coarse_cntl->subtype == FetiInfo::sparse) {
     BLKMatrix->unify(fetiCom);
   } else {
     // New Distributed solver (Padma)
     dscSolver->unify(fetiCom);
   }

   filePrint(stderr," ... Unification Time %e seconds\n",
            (getTime()-t1)/1000.0);

   // sky->print(debugFile);

   if(fetiInfo->coarse_cntl->subtype == FetiInfo::skyline) {
     startTimerMemory(times.pfactor, times.memoryGtGsky);
     // for debuging purposes
     sky->parallelFactor();
     stopTimerMemory(times.pfactor, times.memoryGtGsky);
     glNumRBM  = sky->numRBM();
     GtGsolver = sky;
   } else if (fetiInfo->coarse_cntl->subtype == FetiInfo::sparse) {
     startTimerMemory(times.pfactor, times.memoryGtGsky);
     BLKMatrix->factor();
     stopTimerMemory(times.pfactor, times.memoryGtGsky);
     GtGsolver = BLKMatrix;
   } else {
#ifdef NO_COMPLEX
     startTimerMemory(times.pfactor, times.memoryGtGsky);
     dscSolver->factor();
     stopTimerMemory(times.pfactor, times.memoryGtGsky);
     GtGsolver = dscSolver;
#endif
   }
   if(GtGsolver->numRBM() >= 0)
     filePrint(stderr," ... Number of Global RBMs %2d       ...\n",
              GtGsolver->numRBM() );
   else {
     filePrint(stderr," ... WARNING: RBMs not implemented  ...\n");
     filePrint(stderr," ...          for Padma's solver    ...\n");
   }

 } else {
   // No GtG Solver is used
   GtGsolver = 0;
 }

   stopTimerMemory( times.coarse1, times.memoryGtG);
}

template<class Scalar>
void
GenFetiSolver<Scalar>::addNonLocalGtQG(int subI, int subJ)
{
  // using bool type to determine lowest local neighbor
  if(isLowestLocalNeighbor(subI,subJ) == false) return;

  int kSub, jSub;

  int hasC = (isDynamic && isFeti2) ? 1 : 0;
  // Now we know we are the smallest one loop over all the subs connected
  // to subJ and that are in this CPU
  for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
    int subK = (*subToSub)[subJ][kSub];
    int locK = glSubToLoc[subK];
    if(locK < 0) continue;
    // Now look into K's neighbors, which one is subJ
    for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
      if(sd[locK]->scomm->subNums[jSub] == subJ) {
         if(opControl->sparseGtG)
           opControl->cset[locK].addQContrib(opControl->sparseGtG,
                 jSub+1, opControl->eqNums, hasC);
         break;
      }
  }
}

template<class Scalar>
void
GenFetiSolver<Scalar>::getNonLocalGtQMult(int subI, int subJ, Scalar *va, GenDistrVector<Scalar> *dv) const
{
 // using bool type to determine lowest local neighbor
 if(isLowestLocalNeighbor(subI,subJ) == false) return;

 int kSub, jSub;

 Scalar *alpha = va + opControl->eqNums->firstdof(subJ);
 for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
    int subK = (*subToSub)[subJ][kSub];
    int locK = glSubToLoc[subK];
    if(locK < 0) continue;
    Scalar *ovec = dv->subData(locK);
    // Now look into K's neighbors, which one is subJ
    for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
      if(sd[locK]->scomm->subNums[jSub] == subJ) {
         opControl->cset[locK].getGtQMult(jSub,ovec,alpha);
         break;
      }
  }
}

template<class Scalar>
void
GenFetiSolver<Scalar>::assembleDistGtQGs(int myNum, int *glToLoc)
{
 // int myNum = sd->localSubNum();
 int hasC = (isDynamic && isFeti2) ? 1 : 0;

 opControl->cset[myNum].addQContrib(opControl->sparseGtG,0,
                                    opControl->eqNums,hasC);

 int iSub;
 if (!QGisLocal)
   for(iSub = 0; iSub < opControl->numNeighbSubs(myNum); ++iSub) {
     int neighbN = glToLoc[opControl->neighbSubId(myNum, iSub)];
     if(neighbN >= 0) {
         opControl->cset[neighbN].addQContrib(opControl->sparseGtG,
                   sd[myNum]->scomm->remoteId[iSub]+1, opControl->eqNums, hasC);
     } else
       addNonLocalGtQG(sd[myNum]->subNum(), 
                opControl->neighbSubId(myNum, iSub));
         
   }

}

template<class Scalar>
void
GenFetiSolver<Scalar>::getGtQMult(int myNum, Scalar *va, GenDistrVector<Scalar> *dv) const
{
 Scalar *localvec = dv->subData(myNum);
 Scalar *alpha = va + fetiOps[myNum]->alphaOffset[0];

 opControl->cset[myNum].getGtQMult(localvec,alpha);

 int iSub;

 if(QGisLocal == 0) {
   for(iSub = 0; iSub <  opControl->numNeighbSubs(myNum); ++iSub) {
     int neighbN = glSubToLoc[opControl->neighbSubId(myNum, iSub)];
     if(neighbN >= 0) {
       int myID = sd[myNum]->scomm->remoteId[iSub];
       Scalar *ovec = dv->subData(neighbN);
       opControl->cset[neighbN].getGtQMult(myID,ovec,alpha);
     } else {
       getNonLocalGtQMult(sd[myNum]->subNum(),
         opControl->neighbSubId(myNum, iSub), va, dv);
     }
   }
 }
}
