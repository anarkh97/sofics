#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/SubDomain.h>
#include <Feti.d/Feti.h>
#include <Feti.d/OrthoSet.h>
#include <Threads.d/PHelper.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/Skyline.d/DistSky.h>
#include <Math.d/Skyline.d/BlockSky.h>
#include <Math.d/matrix.h>
#include <Math.d/SymFullMatrix.h>
#include <Math.d/IntFullM.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/linkfc.h>
#include <Timers.d/GetTime.h>

extern int verboseFlag;

// #define DEBUG_MPC
// #define REGULAR_SKY

#ifdef BOOL_NOT_DEFINED
#define false 0
#endif

void
FetiSolver::makeGandFG()
{
 // Get number of each neighbors rbms
 paralApply(nsub, fetiOps, &GenFetiOp<Scalar>::sendNumNeighbRBM, sPat);
 sPat->exchange();
 paralApply(nsub, fetiOps, &FetiOp::getNumNeighbRBM, sPat);

 // Get all the numbers of rigid body modes and dispatch RBMs to neighbors
 makeRbmPat();
 paralApply(nsub, fetiOps, &FetiOp::sendInterfRBM, rbmPat);
 rbmPat->exchange();

 // compute neighbors QGs
 execParal(nsub, this, &FetiSolver::getNeighbFGs);	 
}

void
FetiSolver::getNeighbFGs(int iSub)
{
  int numRBM = opControl->cset[iSub].numGs;
  opControl->cset[iSub].locFGs = new double[numRBM*sd[iSub]->interfLen()];
  sd[iSub]->multMFi(fetiOps[iSub]->solver, opControl->cset[iSub].locGs,
                    opControl->cset[iSub].locFGs, numRBM);

  opControl->cset[iSub].getNeighbQGs(opControl->cset, sd[iSub],
	                             0, rbmPat, betiOps[iSub]->solver);
}

Connectivity *
FetiSolver::makeSingleConnect(Connectivity *coarseConnect, 
                  Connectivity *coarseToSub, Connectivity *subToCoarse,
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

Connectivity *
FetiSolver::getCoarseToSubConnect()
{
   return subToSub->merge(subToSub->merge(mpcToSub));
}

void
FetiSolver::makeSingleCoarse()
{
 startTimerMemory(times.coarse1, times.memoryGtG);

 // first create the Gs and FGs
 makeGandFG();

 // make C (Corner modes) and FC
 if(isFeti2) preProcessCorners();
 
 int glNumSub = subToSub->csize();
 int glNumMpc = mpcToSub->csize();

 Connectivity *coarseToSub = getCoarseToSubConnect();
 Connectivity *subToCoarse = coarseToSub->reverse();

 // GtFG is driving the renumbering both for the GtG system and the GtG
 // off diagonal terms which are then inserted to interlay with GtFG
 Connectivity *coarseFcoarseConnect = coarseToSub->transcon(subToCoarse);
 Connectivity *gtFgConnect = subToSub->transcon(subToSub);

 // We now renumber gtFgConnect on the assumption that it has no 'spurious'
 // singularity
 compStruct gRenum = gtFgConnect->renumByComponent(1);

 // now revisit the numbering by inserting the GtG equations to follow their
 // GtFG counterpart
 int coarseSize  = coarseFcoarseConnect->csize();
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
   //cerr << "GtFG: before: glRenum["<<c<<"] = " << glRenum[c] << endl;
   c++;
   glRenum[c] = invRen[i]+cOffset; 
   //cerr << "   C: before: glRenum["<<c<<"] = " << glRenum[c] << endl;
   c++;
   glRenum[c] = invRen[i]+gOffset; 
   //cerr << " GtG: before: glRenum["<<c<<"] = " << glRenum[c] << endl;
   c++;
 }

#ifdef DEBUG_MPC
 for(i = 0; i < glRenumSize; ++i)
   cerr << "before MPCs: glRenum[" << i << "] = " << glRenum[i] << endl;

 // Number the multiple point constraints separately
 cerr << " cOffset = " << cOffset 
      << " gOffset = " << gOffset
      << " mOffset = " << mOffset
      << endl;
#endif

 int location = 3*glNumSub;
 // YYY int location = 0;
 // cerr << "--- MPC equation numbers ---" << endl;
 for(i=0; i<glNumMpc; ++i) {
   glRenum[location+i] = mOffset + i;
 }

 int *ngRen =  new int[glRenumSize];
 for(i = 0; i < glRenumSize; ++i) {
//   cerr << "glRenum[" << i << "] = " << glRenum[i] << endl;
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
   gRbmSize[ sd[iSub]->subNum() ]  = opControl->cset[iSub].numGs;
   cornerSize[sd[iSub]->subNum() ] = fetiOps[iSub]->getcrnDofSize();
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

 paralApplyToAll(nsub, fetiOps, &FetiOp::setAlphaOffsets, eqNums->allOffsets());

 // create the combined connectivity
 Connectivity *coarseConnect = makeSingleConnect(coarseFcoarseConnect,
                                           coarseToSub, subToCoarse, gOffset);
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
   SkyMatrix *sky;
   if(subToSub->csize() == numCPUs) {
     int neq = eqNums->size();
     int neqPerCPU = neq/numCPUs;
     int remainder = neq%numCPUs;
     int firstAlpha = myCPU * neqPerCPU +
             ((myCPU < remainder) ? myCPU : remainder);
     int nRowAlpha  = neqPerCPU + ((myCPU < remainder) ? 1 : 0);
     sky = new DistSky(coarseConnect, eqNums, tolerance, firstAlpha, nRowAlpha);
     times.memoryGtGDelete = 8*sky->size();
     filePrint(stderr," Deleted GtG Memory: %14.5f Mb\n",
               times.memoryGtGDelete/(1024.0*1024.0));
   } else {
     sky = new SkyMatrix(coarseConnect, eqNums, tolerance);
   }
#else
#ifdef REGULAR_SKY
   SkyMatrix *sky = new SkyMatrix(coarseConnect, eqNums, tolerance);
#else
   BlockSky *sky = new BlockSky(coarseConnect, eqNums, tolerance);
#endif
#endif
   times.memoryGtGsky += threadManager->memoryUsed();
   singleCoarse = sky;
   singleCoarseSolver = sky;
 }
 times.numRBMs = eqNums->size();

 filePrint(stderr," ... Size of Coarse Problem %5d   ...\n",times.numRBMs);

 if( times.numRBMs != 0 || crns != 0 || glNumMpc != 0) {
   singleCoarseAssembly();
 
#ifdef DEBUG_MPC
  FILE *file = fopen("coarse.m","w");
#ifdef REGULAR_SKY
  ((SkyMatrix*) singleCoarse)->print(file);
#else
  ((BlockSky*) singleCoarse)->print();
#endif
#endif

#ifdef DISTRIBUTED
   ((SkyMatrix *)singleCoarse)->unify();
#endif

   startTimerMemory(times.pfactor, times.memoryGtGsky);
   // parallelFactor has a bug when there is a singularity in the coarse matrix
#ifdef DISTRIBUTED
   ((SkyMatrix *)singleCoarse)->parallelFactor();
#else
#ifdef REGULAR_SKY
   ((SkyMatrix *)singleCoarse)->parallelFactor();
#else
   // Dirty to be fixed
   ((BlockSky *)singleCoarseSolver)->parallelFactor();
#endif
#endif
   stopTimerMemory(times.pfactor, times.memoryGtGsky);

   filePrint(stderr, " ... Factoring coarse found %2d RBMS ...\n",
             singleCoarseSolver->numRBM());
 } else {
   singleCoarse=0;
   singleCoarseSolver = 0;
 }

 stopTimerMemory(times.coarse1, times.memoryGtG);
 
}

void
FetiSolver::singleCoarseAssembly()
{
 execParal(nsub, this, &FetiSolver::singleCoarseAssembleG);
 execParal(nsub, this, &FetiSolver::preprocessMPCs);
 execParal(mpcToSub->csize(), this, &FetiSolver::singleCoarseAssembleMPCs);
}

void
FetiSolver::preprocessMPCs(int iSub)
{
  sd[iSub]->getQtKQ(fetiOps[iSub]->solver);
}

void
FetiSolver::singleCoarseAssembleG(int iSub)
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
     addNonLocalGContrib(sd[iSub]->subNum(),opControl->neighbSubId(iSub,jSub));

     if(isFeti2)
      addNonLocalCContrib(sd[iSub]->subNum(),opControl->neighbSubId(iSub,jSub));
   }
#endif
 }

}

void
FetiSolver::singleCoarseAssembleMPCs(int iMPC)
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
                            fetiOps[myNum]->locRBMs,
                            fetiOps[myNum]->solver);
 }

}

#include <Utils.d/linkfc.h>
extern "C" {
void _FORTRAN(dgemv)(const char &, const int &, const int &, const double &,
                     double *, const int &, double *, const int &,
                     const double &, double *, const int &);

void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
                     const int &, const double &, double *, const int &,
                     double *, const int &, const double &, double *,
                     const int &);
}

void
CoarseSet::addGContrib(SparseMatrix *coarseMat, int iSub, 
	               EqNumberer* eqNumber, int gOffset)
{
 double *FG, *G;
 int Gcol, FGcol;
 int nGright;

 if(iSub == 0) { // then it is me
   FG = locFGs;
   G  = locGs;
   FGcol = eqNumber->firstdof(myNum);
   Gcol  = eqNumber->firstdof(myNum+gOffset);
   nGright = numGs;
   iSub -= 1;
 } else {
   iSub -= 1;
   FG = neighbQGs[iSub];
   G  = neighbGs[iSub];
   FGcol = eqNumber->firstdof(neighbs[iSub]);
   Gcol  = eqNumber->firstdof(neighbs[iSub]+gOffset);
   nGright = neighbNumRBMs[iSub];
 }
 // Find the maximum matrix allocation size
 int maxLeft  =  numGs;
 int maxRight =  nGright;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int nGleft = neighbNumRBMs[jSub];
     if(nGleft > maxLeft) maxLeft = nGleft;
 }
 double *stackMem = (double *) dbg_alloca(sizeof(double)*(maxLeft*maxRight));

 if((numGs > 0) && (nGright > 0) && (gSize > 0) && (isDynamic==0)) {
    StackFullM gtg(numGs, nGright,stackMem);
    _FORTRAN(dgemm)('T','N', nGright, numGs, gSize, 1.0,
                      FG, gSize, locGs, gSize,
                      0.0, gtg.data(), nGright);
    coarseMat->add(gtg, eqNumber->firstdof(myNum), FGcol);
    
    if(iSub == -1) {
      _FORTRAN(dgemm)('T','N', nGright, numGs, gSize, 1.0,
                      G, gSize, locGs, gSize,
                      0.0, gtg.data(), nGright);
      // the gi^Tgj appears in two places
      coarseMat->add(gtg, eqNumber->firstdof(myNum), Gcol);
      coarseMat->add(gtg, eqNumber->firstdof(myNum+gOffset), FGcol);
    }
 } 


 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftFGcol = eqNumber->firstdof(neighbs[jSub]);
   int leftGcol = eqNumber->firstdof(neighbs[jSub]+gOffset);
   int nGleft = neighbNumRBMs[jSub];

   if(neighbNumRBMs[jSub] > 0 && nGright > 0) {
     StackFullM gtqg(nGleft, nGright,stackMem);
     _FORTRAN(dgemm)('T','N', nGright, nGleft, subSize(jSub), -1.0,
                      FG+subOffset[jSub], gSize,
                      neighbGs[jSub], leadingDimGs[jSub],
                       0.0, gtqg.data(), nGright);

     coarseMat->add(gtqg, leftFGcol, FGcol);
     if((iSub == -1) && (isDynamic==0)) {
         _FORTRAN(dgemm)('T','N', nGright, nGleft, subSize(jSub), -1.0,
                          G+subOffset[jSub], gSize,
                          neighbGs[jSub], leadingDimGs[jSub],
                          0.0, gtqg.data(), nGright);
        coarseMat->add(gtqg,leftFGcol, Gcol);
        coarseMat->add(gtqg, leftGcol, FGcol);
     }
   }
 }
}

void
CoarseSet::addCContrib(SparseMatrix *coarseMat, int iSub,
                  EqNumberer* eqNumber, int cOffset, int gOffset)
{
 double *FG, *FC;
 int FGcol, FCcol;
 int nGright, nCright;
 double cornerSign;
 int i,j;

 if(iSub == 0) { // then it is me
   FG = locFGs;
   FC = locFBCs;
   FGcol = eqNumber->firstdof(myNum);
   FCcol = eqNumber->firstdof(myNum+cOffset);
   nGright = numGs;
   nCright = numBCs;
   cornerSign = 1.0;
   iSub -= 1;
 } else {
   iSub -= 1;
   FG = neighbQGs[iSub];
   FGcol = eqNumber->firstdof(neighbs[iSub]);
   cornerSign = -1.0;
   nGright = neighbNumRBMs[iSub];
   FC = neighbFBCs2[iSub];
   nCright = neighbCSubIndex[iSub+1] - neighbCSubIndex[iSub];
   int cSubIndex = (neighbCSubIndex[numNeighb] > 0) ?
                        neighbCs[neighbCSubIndex[iSub]][1]  : 0;
   FCcol = eqNumber->firstdof(neighbs[iSub]+cOffset) + cSubIndex; 
 }
 // Find the maximum matrix allocation size
 int maxLeft  =  (numBCs > numGs) ? numBCs : numGs;
 int maxRight =  (nCright > nGright) ? nCright : nGright;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int nGleft = neighbNumRBMs[jSub];
     if(nGleft > maxLeft) maxLeft = nGleft;
     int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
     if(nCleft > maxLeft) maxLeft = nCleft;
 }
 double *stackMem = (double *) dbg_alloca(sizeof(double)*(maxLeft*maxRight));


 if(numGs > 0 && nCright > 0) {
   // This can be cheaper by 'picking'
   StackFullM gtFc(numGs, nCright,stackMem);
   //MLX SIGN adjusted
   _FORTRAN(dgemm)('T','N', nCright, numGs, gSize, cornerSign,
                      FC, gSize,
                      locGs, gSize,
                      0.0, gtFc.data(), nCright);
   coarseMat->add(gtFc, eqNumber->firstdof(myNum), FCcol);
 } 
 if(numBCs > 0 && nGright > 0) {
    StackFullM ctFg(numBCs, nGright,stackMem);
    // MLX SIGN adjusted
    for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nGright; ++j)
           ctFg[i][j] = FG[j*gSize + (*locBCs)[0][i] ];
     coarseMat->add(ctFg,eqNumber->firstdof(myNum+cOffset),FGcol);
 }

 if(numBCs > 0 && nCright > 0) {
     StackFullM ctFc(numBCs, nCright,stackMem);
     int i,j;
     for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nCright; ++j) {
           ctFc[i][j] = cornerSign*FC[j*gSize + (*locBCs)[0][i] ];
       }
     coarseMat->add(ctFc,eqNumber->firstdof(myNum+cOffset),FCcol);
 }


 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftFGcol = eqNumber->firstdof(neighbs[jSub]),
       leftFCcol;

   int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
   if(nCleft > 0)
       leftFCcol = eqNumber->firstdof(neighbs[jSub]+cOffset) +
                         neighbCs[neighbCSubIndex[jSub]][1];

     if(neighbNumRBMs[jSub] > 0 && nCright > 0) {
     // MLX SIGN adjusted
       StackFullM gtFc(neighbNumRBMs[jSub], nCright, stackMem);
        _FORTRAN(dgemm)('T','N', nCright, 
                        neighbNumRBMs[jSub], subSize(jSub), -cornerSign,
                        FC+subOffset[jSub], gSize,
                        neighbGs[jSub], leadingDimGs[jSub],
                        0.0, gtFc.data(), nCright);
       coarseMat->add(gtFc, leftFGcol, FCcol);
     }

     if(nCleft == 0) continue;

     int (*leftCs)[3] = neighbCs + neighbCSubIndex[jSub];
     if(nCleft > 0 && nGright > 0) {
       StackFullM ctFg(nCleft, nGright,stackMem);
       int i,j;
       // MLX SIGN adjusted
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nGright; ++j)
             ctFg[i][j] = -FG[j*gSize + leftCs[i][0] ];
       coarseMat->add(ctFg,leftFCcol,FGcol);
     }

     if(nCleft > 0 && nCright > 0) {
       StackFullM ctFc(nCleft, nCright,stackMem);
       int i,j;
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nCright; ++j)
             ctFc[i][j] = -cornerSign*FC[j*gSize + leftCs[i][0] ];
       coarseMat->add(ctFc,leftFCcol,FCcol);
     }
 }

}

void
CoarseSet::addCGContrib(SparseMatrix *coarseMat,
                  EqNumberer* eqNumber, int cOffset, int gOffset)
{
 if(numBCs == 0 || isDynamic ) return;
 int i,j;

 // Find the maximum matrix allocation size
 int maxGs = numGs;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if(eqNumber->weight(neighbs[jSub]+gOffset) > maxGs)
     maxGs = eqNumber->weight(neighbs[jSub]+gOffset);
 }
 double *stackMem  = (double *) dbg_alloca(sizeof(double)*(numBCs*maxGs));
 double *stackMem2 = (double *) dbg_alloca(sizeof(double)*(numBCs*maxGs));

 StackFullM ctg(numBCs, numGs,stackMem);
 StackFullM gtc(numGs, numBCs,stackMem2);
 ctg.zero();
 gtc.zero();
 int myFCRow = eqNumber->firstdof(myNum+cOffset);
 int myGRow  = eqNumber->firstdof(myNum+gOffset);

 if ((numGs>0) && (numBCs >0) )  {
   for (i=0; i< numGs; i++)
     for (j=0; j< numBCs; j++) {
        int GoffSet = i * gSize;
     // MLX SIGN adjusted
        gtc[i][j] = ctg[j][i] = locGs[GoffSet + (*locBCs)[0][j]];
       }
 }
 coarseMat->add(ctg, myFCRow, myGRow);
 coarseMat->add(gtc, myGRow, myFCRow);

  // Now loop on the neighbors
 int iC = 0;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int neighbNumGs = eqNumber->weight(neighbs[jSub]+gOffset);
   StackFullM ctg(numBCs, neighbNumGs,stackMem);
   StackFullM gtc(neighbNumGs, numBCs, stackMem2);
   ctg.zero();
   gtc.zero();
   while(iC < numBCs && (*locBCs)[3][iC] == jSub) {
     int i;
     for(i =0; i < neighbNumGs; ++i)
       // MLX SIGN adjusted
       gtc[i][iC] = ctg[iC][i] = 
         -neighbGs[jSub][i*leadingDimGs[jSub]+(*locBCs)[0][iC]-subOffset[jSub]];
     iC = iC+1;
   }
  coarseMat->add(ctg,myFCRow, eqNumber->firstdof(neighbs[jSub]+gOffset));
  coarseMat->add(gtc,eqNumber->firstdof(neighbs[jSub]+gOffset), myFCRow);
 }
}

void
CoarseSet::addMPCContrib(int iMPC, SparseMatrix *coarseMat,
                EqNumberer* eqNumber, int cOffset, int gOffset, int mpcOffset,
                double *locR, Solver *s)
{

 int numMPCs = sd->numMPCs();

 // Find the maximum matrix allocation size
 int maxRight = (numMPCs > numGs) ? numMPCs : numGs;
 maxRight     = (numBCs > maxRight) ? numBCs : maxRight;

 // Find the maximum matrix allocation size
 int maxLeft  =  maxRight;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int nGleft = neighbNumRBMs[jSub];
     if(nGleft > maxLeft) maxLeft = nGleft;
     int nCleft = (neighbCSubIndex) ? 
         neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub] : 0;
     if(nCleft > maxLeft) maxLeft = nCleft;
 }
 double *stackMem = (double *) dbg_alloca(sizeof(double)*(maxLeft*numMPCs));

 // This sub-domain does not have a contribution
 if(numMPCs == 0) return;

 int myFGcol  = eqNumber->firstdof(myNum);
 int myCCol   = eqNumber->firstdof(myNum+cOffset);
 int myGCol   = eqNumber->firstdof(myNum+gOffset);

 // KHP: note the mpc numbering has to be done carefully. 
 //      this number should be done per MPC as a sub-domain may have
 //      contributions to only a few non-consecutive MPCs
 int myMPCRow = eqNumber->firstdof(mpcOffset);
 // cerr << "myFGcol " << myFGcol << "  myGCol  " << myGCol << endl;

 // add QtKQ
 sd->getQtKQ(iMPC,stackMem);
 StackFullM qtKq(numMPCs, 1, stackMem);
#ifdef DEBUG_MPC
 qtKq.print("--- QtKQ ---","QtKQ");
#endif
 int jMPC;
 FullM qtkq(1,1);
 for(jMPC=0; jMPC<numMPCs; ++jMPC) {
   qtkq[0][0] = qtKq[jMPC][0];
#ifdef DEBUG_MPC
   fprintf(stderr,"QtKQ: Where am I adding this %d %d\n",
           eqNumber->firstdof(mpcOffset-iMPC)+sd->localToGlobalMPC[jMPC], myMPCRow);
#endif
   coarseMat->add(qtkq, 
                 eqNumber->firstdof(mpcOffset-iMPC)+sd->localToGlobalMPC[jMPC], 
                 myMPCRow);
 }

 // add QtKBtG
 if(numGs > 0) {
   StackFullM GtBKQ(numGs,1,stackMem);
   GtBKQ.zero(); // NOTE: Never remove this! 
   getGtMult(sd->QtKpBt+sd->globalToLocalMPC[iMPC]*sd->interfLen(), GtBKQ.data());
#ifdef DEBUG_MPC
   GtBKQ.print("--- GtBKQ ---","GtBKQ");
   fprintf(stderr,"GtBKQ: Where am I adding this %d %d\n",myFGcol, myMPCRow);
#endif
   coarseMat->add(GtBKQ, myFGcol, myMPCRow);
 }

 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftFGcol = eqNumber->firstdof(neighbs[jSub]);
   int nGleft    = neighbNumRBMs[jSub];
   StackFullM GtBKQ(nGleft, 1, stackMem);
   if(neighbNumRBMs[jSub] > 0 ) {
     _FORTRAN(dgemm)('T', 'N', 1, nGleft, subSize(jSub), -1.0,
          sd->QtKpBt+subOffset[jSub]+sd->globalToLocalMPC[iMPC]*sd->interfLen()
          ,gSize, neighbGs[jSub], leadingDimGs[jSub],
          0.0, GtBKQ.data(), 1);
#ifdef DEBUG_MPC
     GtBKQ.print("---- Neighbors Gt (BKQ) ----","Gt^jBKQ");
#endif
     coarseMat->add(GtBKQ, leftFGcol, myMPCRow);
   }
 }


 // add QtKBtC
 if(numBCs > 0) {
   StackFullM CtBKQ(numBCs, 1, stackMem);
   getCtMult(sd->QtKpBt+sd->globalToLocalMPC[iMPC]*sd->interfLen(), CtBKQ.data());
   int i;
   for(i=0; i<numBCs; ++i)
     CtBKQ[i][0] = -CtBKQ[i][0];
   coarseMat->add(CtBKQ, myCCol, myMPCRow);
#ifdef DEBUG_MPC
   CtBKQ.print("---- CtBKQ ---","CtBKQ");
#endif
 }

 // Now loop on the neighbors
 double cornerSign = 1.0;
   for(jSub =0; jSub < numNeighb; ++jSub) {
     int leftFCcol;

     int nCleft = (neighbCSubIndex) ?
                   neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub] : 0;
     if(nCleft > 0)
         leftFCcol = eqNumber->firstdof(neighbs[jSub]+cOffset) +
                      neighbCs[neighbCSubIndex[jSub]][1];

       if(nCleft == 0) continue;

       int (*leftCs)[3] = neighbCs + neighbCSubIndex[jSub];

       if(nCleft > 0) {
         StackFullM CtBKQ(nCleft, 1, stackMem);
         double *BKQ = sd->QtKpBt;
         int i;
         for(i = 0; i < nCleft; ++i)
            CtBKQ[i][0] = -cornerSign*BKQ[sd->globalToLocalMPC[iMPC]*gSize + leftCs[i][0] ];
         coarseMat->add(CtBKQ, leftFCcol, myMPCRow);
       }
   }

 // add QtR
 if((numGs > 0) && (isDynamic==0)) {
   StackFullM QtR(1, numGs, stackMem);
   QtR.zero();
   sd->multQt(iMPC, locR, numGs, QtR.data());
   FullM RtQ = QtR.transpose();
#ifdef DEBUG_MPC
   if(myCPU == 0)
     RtQ.print("--- RtQ --- ","RtQ");
     cerr << "Row " << myGCol << " Col " << myMPCRow << endl;
#endif
   coarseMat->add(RtQ, myGCol, myMPCRow);
   // YYY coarseMat->add(RtQ, myMPCRow, myGCol);
 }
}

void
FetiSolver::computeL0(DistrVector &f, DistrVector &r,  DistrVector &u,
                      Vector &alpha, DistrVector &lambda)
{
 // note: r = B K^+ f
#ifdef DISTRIBUTED
 alpha.zero();
#endif

 double *singleC = alpha.data();
 opControl->vec2 = singleC;
 opControl->dv1  = &r;

 execParal(mpcToSub->csize(), this, &FetiSolver::getSQtMult, &u, singleC);
 execParal(mpcToSub->csize(), this, &FetiSolver::addMpcRhs, singleC);  

 // -Rt (f - B^t lambda)   NOTE: lambda could be zero in this computation
 execParal(nsub, this, &FetiSolver::getSRMult,  &f, &lambda, singleC); 

 execParal(nsub, this, &FetiSolver::getGtMult,  &r, singleC);

 if(isFeti2)
   execParal(nsub, this, &FetiSolver::getSCtMult, &r, singleC);

#ifdef DISTRIBUTED
if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif
 
// filePrint(stderr,"Before: alpha*alpha %e\n",alpha*alpha);
 if(singleCoarseSolver) singleCoarseSolver->reSolve(singleC);
// filePrint(stderr," After: alpha*alpha %e\n",alpha*alpha);

 // Not needed because done on the outside: lambda.zero(); 
 execParal(nsub, this, &FetiSolver::addG, &lambda, singleC);
 
 if(isFeti2)
   execParal(nsub, this, &FetiSolver::addC, &lambda, singleC);
}

void
FetiSolver::addMpcRhs(int iMPC, double *sv)
{
 double *gamma = sv + eqNums->firstdof(mOffset+iMPC);
#ifdef DISTRIBUTED
 int myNum = glSubToLoc[(*mpcToSub)[iMPC][0]];
#else
 int myNum = (*mpcToSub)[iMPC][0];
#endif
 if(myNum >= 0)
   gamma[0] += sd[myNum]->getMpcRhs(iMPC);
}

void
FetiSolver::singlePr(DistrVector &y, DistrVector &p, Vector &alpha)
{
 startTimerMemory(times.project, times.memoryProject1);

 alpha.zero();

 double *singleC = alpha.data();

 execParal(mpcToSub->csize(), this, &FetiSolver::getQtKpBMult, &y, singleC);
 
 execParal(nsub, this, &FetiSolver::getSGtMult,   &y, singleC);

 execParal(nsub, this, &FetiSolver::getFGMult,    &y, singleC);

 if(isFeti2)
   execParal(nsub, this, &FetiSolver::getFCMult,  &y, singleC);

#ifdef DISTRIBUTED
 if(fetiCom) fetiCom->globalSum(alpha.size(), alpha.data());
#endif

 if(singleCoarseSolver) singleCoarseSolver->reSolve(singleC);

 p = y;
 execParal(nsub, this, &FetiSolver::addG, &p, singleC);

 if(isFeti2)
   execParal(nsub, this, &FetiSolver::addC, &p, singleC);

 stopTimerMemory(times.project, times.memoryProject1);
}

void
FetiSolver::addGs(DistrVector &r, DistrVector &w, Vector &alpha)
{
 double *singleC = alpha.data();
 w = r;
 if(isDynamic==0)
   execParal(nsub, this, &FetiSolver::addSG, &w, singleC);
}

void
FetiSolver::getSRMult(int iSub, DistrVector *r, DistrVector *lambda, double *sv)
{
 int nRBM = fetiOps[iSub]->numRBM;

 if((nRBM==0) || isDynamic) return;

 // alpha = R*lvec
 double *lvec = r->subData(iSub);
 double *lbvec = lambda->subData(iSub);
 double *locRBMs = fetiOps[iSub]->locRBMs;
 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum()+gOffset);
 sd[iSub]->getSRMult(lvec, lbvec, nRBM, locRBMs, alpha);
}

void
FetiSolver::addC(int iSub, DistrVector *lambda, double *sv)
{
 double *localvec = lambda->subData(iSub);
 int cOffset      = subToSub->csize();

 opControl->cset[iSub].addNuC(localvec, sv, eqNums, cOffset);
}

void
CoarseSet::addNuC(double *vec, double *sv, EqNumberer *eqNums, int cOffset)
{
 int i;
 double *nu;
 nu = sv + eqNums->firstdof(myNum+cOffset);
 for (i=0; i< numBCs; i++) {
   int index = (*locBCs)[0][i];
   vec[index] += nu[i];
 }
 int iSub, iC;
 // ML WARNING !!! CAREFUL ABOUT THE SIGN HERE! NOT SAME CONVENTION AS G
 for(iSub = 0; iSub < numNeighb; ++iSub)
   for(iC = neighbCSubIndex[iSub]; iC < neighbCSubIndex[iSub+1]; ++iC) {
      double *nu = sv + eqNums->firstdof(neighbs[iSub]+cOffset);
      vec[neighbCs[iC][0]] -= nu[neighbCs[iC][1]];
   }
}

void
FetiSolver::getSGtMult(int iSub, DistrVector *r, double *sv)
{
 int numRBM = opControl->cset[iSub].numGs;
 if(numRBM==0) return;
 double *lvec = r->subData(iSub);
 double *locGs = opControl->cset[iSub].locGs;

 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum()+gOffset);

 if(numRBM > 0)
 _FORTRAN(dgemv)('T',opControl->cset[iSub].gSize, numRBM, -1.0, locGs,
                  opControl->cset[iSub].gSize,lvec, 1, 0.0, alpha, 1);
}

void
FetiSolver::getSCtMult(int iSub, DistrVector *r, double *sv)
{
 int cOffset = subToSub->csize();
 double *lvec = r->subData(iSub);
 double *beta = sv + eqNums->firstdof(sd[iSub]->subNum()+cOffset);
 opControl->cset[iSub].getCtMult(lvec,beta);
/* IF different sign convention MLX
 int i;
 for(i = 0; i < eqNums->weight(sd[iSub]->subNum()+cOffset); ++i)
   beta[i] = -beta[i];
*/
}

void
FetiSolver::getSQtMult(int iMpc, DistrVector *u, double *sv)
{
 int numSubAttached = mpcToSub->num(iMpc);
 int iSub;
 double *gamma = sv + eqNums->firstdof(mOffset + iMpc);
 gamma[0]=0.0;
 for(iSub=0; iSub<numSubAttached; ++iSub) {
#ifdef DISTRIBUTED
   int myNum = glSubToLoc[(*mpcToSub)[iMpc][iSub]];
   if(myNum >= 0) {
     double *lvec = u->subData(myNum);
     sd[myNum]->multQt(iMpc,lvec, gamma);
   }
#else
   int myNum = (*mpcToSub)[iMpc][iSub];
   double *lvec = u->subData(myNum);
   sd[myNum]->multQt(iMpc,lvec, gamma);
#endif
 }

 gamma[0]=-gamma[0];
}

void
FetiSolver::getQtKpBMult(int iMpc, DistrVector *r, double *sv)
{
 // Q^t K^+ B^t r
 int numSubAttached = mpcToSub->num(iMpc);
 int iSub;
 double *gamma = sv + eqNums->firstdof(mOffset+iMpc);
 gamma[0] = 0.0;
 for(iSub=0; iSub<numSubAttached; ++iSub) {
#ifdef DISTRIBUTED
     int myNum = glSubToLoc[(*mpcToSub)[iMpc][iSub]];
     if(myNum >= 0) {
       double *lvec = r->subData(myNum);
       sd[myNum]->multQtKBt(iMpc, lvec, gamma);
     } 
#else
     int myNum = (*mpcToSub)[iMpc][iSub];
     double *lvec = r->subData(myNum);
     sd[myNum]->multQtKBt(iMpc, lvec, gamma);
#endif
 }

 gamma[0]=-gamma[0];
}

void
FetiSolver::getGtMult(int iSub, DistrVector *r, double *sv)
{
 int numRBM = opControl->cset[iSub].numGs;
 if(numRBM==0) return;
 double *lvec  = r->subData(iSub);
 double *locGs = opControl->cset[iSub].locGs;
 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum());

 if(opControl->cset[iSub].gSize > 0)
 _FORTRAN(dgemv)('T',opControl->cset[iSub].gSize, numRBM, -1.0, locGs,
                 opControl->cset[iSub].gSize,lvec, 1, 0.0, alpha, 1);
}

void
FetiSolver::addG(int iSub, DistrVector *r, double *sv)
{
 double *lvec  = r->subData(iSub);
 double *locGs = opControl->cset[iSub].locGs;
 int numRBM    = opControl->cset[iSub].numGs;
 SComm *scomm  = sd[iSub]->getSComm();
 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum());

 // lvec += locGs*alpha
 if(numRBM > 0 && opControl->cset[iSub].gSize > 0)
  _FORTRAN(dgemv)('N',opControl->cset[iSub].gSize, numRBM, 1.0, locGs,
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
    _FORTRAN(dgemv)('N', subSize, nRBM, -1.0,
          opControl->cset[nghb].subG(iSub), opControl->cset[nghb].gSize,
          alpha, 1, 1.0, lvec+subOffset, 1);
#endif
 }

}

void
FetiSolver::addSG(int iSub, DistrVector *r, double *sv)
{
 double *lvec  = r->subData(iSub);
 double *locGs = opControl->cset[iSub].locGs;
 int numRBM    = opControl->cset[iSub].numGs;
 SComm *scomm  = sd[iSub]->getSComm();
 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum()+gOffset);

 if(numRBM > 0 && opControl->cset[iSub].gSize > 0)
 _FORTRAN(dgemv)('N',opControl->cset[iSub].gSize, numRBM, 1.0, locGs,
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
    _FORTRAN(dgemv)('N', subSize, nRBM, -1.0,
          opControl->cset[nghb].subG(iSub), opControl->cset[nghb].gSize,
          alpha, 1, 1.0, lvec+subOffset, 1);
#endif
 }
}

void
FetiSolver::getFGMult(int iSub, DistrVector *r, double *sv)
{
 double *lvec = r->subData(iSub);
 double *locFGs = opControl->cset[iSub].locFGs;

 int numRBM = opControl->cset[iSub].numGs;
 int gSize = opControl->cset[iSub].gSize;

 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum());
 
if(numRBM > 0)
 _FORTRAN(dgemv)('T', gSize, numRBM, -1.0,
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
         double *ovec = r->subData(neighbN);
         gSize = sd[neighbN]->interfLen();
         if(numRBM > 0)
         _FORTRAN(dgemv)('T', gSize, numRBM, -1.0,
                opControl->cset[neighbN].neighbQGs[myID], gSize,
                 ovec, 1, 1.0, alpha, 1);
                  
#ifdef DISTRIBUTED
     } else {
       getNonLocalSubAlphaGtQ(sd[iSub]->subNum(),
            opControl->neighbSubId(iSub, jSub), sv, r);
     }
#endif
 }
}

void
FetiSolver::getFCMult(int iSub, DistrVector *r, double *sv)
{
 int cOffset = subToSub->csize();
 double *lvec = r->subData(iSub);
 double *locFCs = opControl->cset[iSub].locFBCs;

 int numBCs = opControl->cset[iSub].numBCs;
 if(numBCs <= 0 ) return;
 int gSize = opControl->cset[iSub].gSize;

 double *alpha = sv + eqNums->firstdof(sd[iSub]->subNum()+cOffset);

 if(gSize > 0)
 _FORTRAN(dgemv)('T', gSize, numBCs, -1.0,
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
         double *ovec = r->subData(neighbN);
         opControl->cset[neighbN].getCtFMult(myID, ovec, alpha);

#ifdef DISTRIBUTED
     } else {
       getNonLocalFCtMult(sd[iSub]->subNum(),
            opControl->neighbSubId(iSub, jSub), sv, r);
     }
#endif
 }
}

void
FetiSolver::addRS(int iSub, DistrVector *vec1, double *vec2)
{
 double *localvec = vec1->subData(iSub);
 double *alpha    = vec2 + eqNums->firstdof(sd[iSub]->subNum() + gOffset);
 int m = sd[iSub]->localLen();
 int n = fetiOps[iSub]->numRBM;
 // localvec = R*alpha
 if(n > 0)
 _FORTRAN(dgemv)('N', m, n, 1.0, fetiOps[iSub]->locRBMs, m,
                 alpha, 1, 1.0, localvec, 1);
}

void
FetiSolver::addRSingle(DistrVector &v, Vector &alpha)
{

 times.addR -= getTime();
 execParal(nsub, this, &FetiSolver::addRS, &v, alpha.data());
 times.addR += getTime();
}

void
FetiSolver::singleCoarseSolve(DistrVector &f, DistrVector &u)
{
 times.solve -= getTime();

 // K H P: temporary fix for MPCs
 // resetOrthoSet();

 DistrVector &r       = wksp->ret_r();          // residual
 DistrVector &lambda0 = wksp->ret_lambda();     // Lagrange multipliers
 DistrVector &w       = wksp->ret_w();          // projected residual
 DistrVector &y       = wksp->ret_y();
 DistrVector &z       = wksp->ret_z();          // preconditioned residual
 DistrVector &p       = wksp->ret_p();
 DistrVector &Fp      = wksp->ret_Fp();
 DistrVector &du      = wksp->ret_du();
 DistrVector &pr      = wksp->ret_pr();
 DistrVector &uzero   = wksp->ret_uzero();
 uzero.zero();

 DistrVector &deltaU  = wksp->ret_deltaU();
 deltaU.zero();
 DistrVector &deltaF  = wksp->ret_deltaF();
 deltaF.zero();
 Vector &alpha        = wksp->ret_alpha();

 // Compute the square of a pseudo-norm of f
 double pseudoFNormSq = f*f;
 if(pseudoFNormSq == 0.0) {
   u.zero();
   return;
 }

 Vector beta(alpha.size(),0.0);

 lambda0.zero();

 // First compute a forward backward to get B K^+ f
 // u = K^+ (f + B^t lambda0) = K^+ f
 // r = B u = B K^+ f
 localSolveAndJump(f,lambda0,u,r);

 computeL0(f, r, u, beta, lambda0);

 if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10) {
   filePrint(stderr," ... lambda0*lambda0 = %e ...\n",lambda0*lambda0);
   filePrint(stderr," ... f*f             %10.4e     ...\n",pseudoFNormSq);
   filePrint(stderr," ... r*r             %10.4e     ...\n",r*r);
   filePrint(stderr," ... u*u             %10.4e     ...\n",u*u);
 }

 // u = K^+ (f + B^t lambda0 + Q beta)
 // r = B u
 //   = B K^+ f + B K^+ B^t lambda0 + B K^+ Q beta
 //   = B K^+ f + F lambda0 + B K^+ Q beta
 double xx = localSolveAndJump(f,lambda0,beta,u,r);

 addGs(r, w, beta);

 // For multi-solve cases 
 int hasPred = predict(w,lambda0);
 if(hasPred) {
   // singlePr(lambda0, Fp, beta);
   localSolveAndJump(f, lambda0, u, r);
   Fp = lambda0;
   computeL0( f, r, u, beta, Fp);
   localSolveAndJump(f,Fp,beta,u,r);
   addGs(r, w, beta);
   if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
     filePrint(stderr," ... lambda0*lambda0 = %e beta*beta %e w*w %e\n", lambda0*lambda0,beta*beta,w*w);
 }

 double w0Norm2 = w*w;

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
   double pFp = localSolveAndJump(uzero, p, alpha, du, Fp);

   double rp = w*y;

   r = y;
   // compute y = Fp (true Fp where F is the condensed operator)
   addGs(Fp, y, alpha);
   pFp = y*r;
   double nu = - rp / pFp ;

   lambda0.linAdd(nu,r);
   // Update residual w: w = w + nu * y
   // Update solution u: u = u + nu * du
   doubleUpdate(nu,u,du,w,y);


   // Add search direction to Orthogonalization set
   orthoAdd(r, y, pFp);

   alpha *= nu;
   beta  += alpha;

   // Precondition: z = F^-1 w
   errorEstimator = preCondition(w, z);

   if((iter % fetiInfo->numPrint() == 0) &&
      (fetiInfo->numPrint() > 0) && verboseFlag ) {
     double wrel = sqrt((w*w)/w0Norm2);
     filePrint(stderr," %4d %23.6e %21.6e\n",
               iter+1, sqrt(errorEstimator/pseudoFNormSq), wrel);
   }

   // Test for convergence or maximum number of iterations, this way,
   // on exit of the loop we are guaranteed that alpha is correct
   if((z2 = errorEstimator) < epsilon2 * pseudoFNormSq || iter == maxiter-1) {
     times.iterations[numSystems].stagnated = 0;
     if(errorEstimator < epsilon2 * pseudoFNormSq)
       times.converged = 1;
     break;
   }

   // Check for stagnation
   if(std::abs(z2-lastz2) < 1.0e-6*lastz2) {
     times.setStagnate(numSystems);
     filePrint(stderr,"STAGNATION: Relative Primal Error Reached = "
                      "%e %e %e %e\n",sqrt(z2/pseudoFNormSq), std::abs(z2-lastz2),
                      z2, lastz2);
     break;
   }
   lastz2 = z2;

   orthogonalize(z,y);
   singlePr(y, p, alpha);

 }

 // fprintf(stderr, "Lambda at end %e\n", lambda0*lambda0);
 // fprintf(stderr, "(Beta at end) Initial beta %e\n", beta*beta);

 // u = u + R beta
 if(isDynamic == 0)
   addRSingle(u, beta);

   // fprintf(stderr, "U before %e\n", u*u);

   // localSolveAndJump(f,lambda0,u,r);
   // fprintf(stderr, "U before (2) %e\n", u*u);
 
   // Fp = lambda0;
   // computeL0(f, r, u, beta, Fp);
   //fprintf(stderr, "Now beta %e\n", beta*beta);
   //localSolveAndJump(f,Fp,beta,u,r);
   //addGs(r, w, beta);
   //if(isDynamic == 0)
   //  addRSingle(u, beta);
  //fprintf(stderr, "U after %e\n", u*u);
  //fprintf(stderr, "W at end %e\n", w*w);

#ifdef DEBUG_MPC
 int numMPCs = mpcToSub->csize();
 if(numMPCs > 0) {
  alpha.zero();
  execParal(mpcToSub->csize(), this, &FetiSolver::getSQtMult, &u, alpha.data());
  execParal(mpcToSub->csize(), this, &FetiSolver::addMpcRhs, alpha.data());
   filePrint(stderr,"%d Is this zero enough?\n",numMPCs);
   int iMPC;
   for(iMPC=0; iMPC<numMPCs; ++iMPC)
     fprintf(stderr,"%d %e\n",iMPC+1,alpha[alpha.size()-numMPCs+iMPC]);
  }
#endif
 
 // make solution compatible u = u + deltaU
 makeCompatible( u, deltaU );

 // Store number of iterations, primal error and dual error
 setAndStoreInfo( iter+1, (errorEstimator/pseudoFNormSq), (w*w)/w0Norm2 );

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
void
FetiSolver::preProcessCorners()
{
 // Initialize numCRNs, crnDofLen and BClocal
 paralApplyToAll(nsub, fetiOps, &FetiOp::initializeCRNs, sPat);

 // Get total number of corner modes from subdomains
 crns = collectIntGlobalSum();
#ifdef DISTRIBUTED
 fprintf(stderr,"CPU %d found %d corners\n",fetiCom->myID(),crns);
 crns = fetiCom->globalSum(crns);
#endif
 times.numCRNs = crns;

 if(fetiInfo->numPrint() > 0 && fetiInfo->numPrint() < 10)
   filePrint(stderr, " ... %5d corners were found       ...\n", crns);

 // Initialize neighbors number of corner modes
 sPat->exchange();
 paralApplyToAll(nsub, fetiOps, &FetiOp::getNumNeighbCRNs, sPat);

 // cset is constructed in FetiSolver::factorMatrices
 CoarseSet *cset = opControl->cset;
 //
 // KENDALL, this routine needs modification to run FETI-2 in Distributed!
 //
 paralApplyToAll(nsub, cset, &CoarseSet::buildBNeighbCs, cset);

 // compute F C

 // compute my own Fi^{-1} BCi
 paralApplyToAll(nsub, fetiOps, &FetiOp::computeFiBC);

 // compute neighbors F C
 paralApplyToAll(nsub, fetiOps, &FetiOp::getNeighbFC);
}

double
FetiSolver::localSolveAndJump(DistrVector &ifrc, DistrVector &bf,
                            Vector &beta, DistrVector &u, DistrVector &lambda) const
{
 startTimerMemory(times.sAndJ, times.memorySAndJ);

 // u = K+(f - B^t lambda - Q^t beta)
 timedParal(times.solveAndJump, nsub, this, &FetiSolver::localSolve2, &u,
              &lambda, &beta, &ifrc, &bf);
 vPat->exchange();
 timedParal(times.solveAndJump, nsub, this, &FetiSolver::interfDiffAndDot,
              bf, lambda);

 // Sum each subdomains dot product contribution
 double ret = 0.0;
 int i;
 for(i = 0; i < nsub; ++i)
   ret += fetiOps[i]->res;

#ifdef DISTRIBUTED
 if(fetiCom) ret = fetiCom->globalSum(ret);
#endif

 stopTimerMemory(times.sAndJ, times.memorySAndJ);

 return ret;
}

void
FetiSolver::localSolve2(int iSub, DistrVector *v1, DistrVector *v2,
                        Vector *beta, DistrVector *v3, DistrVector *v4)
{
 SubDomain *subd = sd[iSub];

 int sn = subd->localSubNum();

 double *localvec  = v1->subData(sn);
 double *interfvec = v2->subData(sn);
 double *localsrc  = v3->subData(sn);
 double *interfsrc = v4->subData(sn);

 int i;
 for(i = 0; i < subd->localLen(); ++i)
    localvec[i] = localsrc[i];

 for(i = 0; i < subd->interfLen(); ++i)
    interfvec[i] = interfsrc[i];

 double *gamma = beta->data() + eqNums->firstdof(mOffset);

 sd[iSub]->fetiBaseOp(fetiOps[iSub]->solver, localvec, interfvec, gamma);

 sd[iSub]->sendInterf(interfvec, fetiOps[iSub]->interfBuff);

}

#ifdef BOOL_NOT_DEFINED
int
#else
bool
#endif
FetiSolver::isLowestLocalNeighbor(int subI, int subJ)
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



void
FetiSolver::addNonLocalGContrib(int subI, int subJ)
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

void
FetiSolver::addNonLocalCContrib(int subI, int subJ)
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

void
FetiSolver::getNonLocalSubAlphaGtQ(int subI, int subJ, double *va,
                                   DistrVector *dv)
{
 if(isLowestLocalNeighbor(subI,subJ) == false) return;

 int kSub, jSub;

 double *alpha = va + eqNums->firstdof(subJ);
 for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
    int subK = (*subToSub)[subJ][kSub];
    int locK = glSubToLoc[subK];
    if(locK < 0) continue;
    double *ovec = dv->subData(locK);
    // Now look into K's neighbors, to find which one is subJ
    for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
      if(sd[locK]->scomm->subNums[jSub] == subJ) {
         opControl->cset[locK].subAlphaGtQ(jSub, ovec, alpha);
         break;
      }
  }
}

void
FetiSolver::getNonLocalFCtMult(int subI, int subJ, double *va,
                                   DistrVector *dv)
{
 if(isLowestLocalNeighbor(subI,subJ) == false) return;

 int kSub, jSub;

 double *alpha = va + eqNums->firstdof(subJ);
 for(kSub = 0; kSub < subToSub->num(subJ); ++ kSub) {
    int subK = (*subToSub)[subJ][kSub];
    int locK = glSubToLoc[subK];
    if(locK < 0) continue;
    double *ovec = dv->subData(locK);
    // Now look into K's neighbors, to find which one is subJ
    for(jSub = 0; jSub < sd[locK]->scomm->numNeighb; ++jSub)
      if(sd[locK]->scomm->subNums[jSub] == subJ) {
         opControl->cset[locK].getCtFMult(jSub, ovec, alpha);
         break;
      }
  }
}



