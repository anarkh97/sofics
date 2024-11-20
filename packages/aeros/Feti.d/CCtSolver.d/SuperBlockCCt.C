#include <set>
#include <list>

#include "SuperBlockCCt.h"
#include "Feti.d/FetiSub.h"
#include <Driver.d/SubDomain.h>
#include <Utils.d/DistHelper.h>
#include <Feti.d/Feti.h>
#include <Solvers.d/SolverFactory.h>

extern Domain * domain;

template<class Scalar>
SuperBlockCCtSolver<Scalar>::SuperBlockCCtSolver(const Connectivity *_blockToMpc, const Connectivity *_mpcToMpc,
                                                 const Connectivity *mpcToSub,
                                                 const Connectivity *_mpcToCpu, int _numSubsWithMpcs,
                                                 std::vector<FetiSub<Scalar> *> subsWithMpcs,
                                                 FetiInfo *_finfo, FSCommunicator *_fetiCom, bool super_flag,
                                                 bool sub_flag)
        : CCtSolver<Scalar>(std::move(subsWithMpcs))
{
  filePrint(stderr," ... Building block CC^t for preconditioning MPCs ...\n");
  initialize();
  blockToMpc = _blockToMpc;
  nMpcBlocks = blockToMpc->csize();
  mpcToMpc = _mpcToMpc;
  this->mpcToCpu = _mpcToCpu;
  this->numSubsWithMpcs = _numSubsWithMpcs;
  this->fetiCom = _fetiCom;
  this->glNumMpc = mpcToMpc->csize();
  myCPU = this->fetiCom->cpuNum();
  numCPUs = this->fetiCom->size();
  finfo = _finfo;
  blockMpcToMpc = 0;

  /// refine decomposition blocks
  if(sub_flag) makeSubBlocks(); // attempt further decomposition of every block into sub blocks
  if(super_flag) makeSuperBlocks(); // combine blocks into super blocks for optimal load sharing between cpus
  else { // one block per super block 
    mpcToBlock = blockToMpc->alloc_reverse();
    int *newptr = new int[2]; newptr[0] = 0; newptr[1] = newptr[0]+nMpcBlocks;
    int *newtarget = new int[nMpcBlocks];
    for(int i=0;i<nMpcBlocks;i++) { newtarget[i] = i; }
    MPITosuperBlock = new Connectivity(1,newptr,newtarget);
    nBigBlocksperMPI = new int[1]; nBigBlocksperMPI[0] = nMpcBlocks;
  }

  createSuperBlockCCt(mpcToSub);

#ifdef DISTRIBUTED
  makeBlockCCtCommPattern();
  // can be optimized by only looping on the blocks that really need to be sent
  // currently loop over ALL the blocks that have lmpc CCt contribution from myCPU (partially or fully)
  execParal(nLocAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::sendOneBlockCCtsolver);
  blockCCtPat->exchange();
  if(nMpcBlocksOnMyCPU){
    execParal(nMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::recOneBlockCCtsolver);
  }
  deleteMyCPUTmpBlockCCtsolver();
#endif

  // allocate memory for mpcv (residual vectors for each block represented on myCPU)
  // Can be optimized by creating a mpcv array that has the size of nExtMpcBlocksOnMyCPU
  // instead of nMpcBlocks -> this will require modification of extractBlocksMpc ....
  mpcv = new GenVector<Scalar>*[nMpcBlocks];
  for(int i=0;i<nExtMpcBlocksOnMyCPU;i++) {
    int blkId = myCPUExtBlockIdArray[i];
    mpcv[blkId] = new GenVector<Scalar>(blockToMpc->num(blkId));
  }
                                                                                                                                             
  // make the mapping from myCPU blocks to their "Cpu-shared" lmpcs. Used in method
  // solveOneBlockCCt
#ifdef DISTRIBUTED
  makeBlkToCpuSharedMpcsMap();
#endif
}

template<class Scalar>
SuperBlockCCtSolver<Scalar>::~SuperBlockCCtSolver()
{                                                                                                                         
  deleteBlockCCtsolver();
  deleteBlockMpcToMpcConnectivity();
  //deleteMpcVector();
  delete blockToMpc;
  delete mpcToBlock;
  delete blockToSub;
#ifdef DISTRIBUTED
  delete blockToCpu;
  delete cpuToBlock;
  delete blockCCtPat;
  delete mpcvPat1;
  delete mpcvPat2;
#endif
  if(nBigBlocksperMPI)       { delete [] nBigBlocksperMPI     ; nBigBlocksperMPI       = 0; }
  if(nSmallBlocksperMPI)     { delete [] nSmallBlocksperMPI   ; nSmallBlocksperMPI     = 0; }
  deleteBlockMpcEqNums();
  deleteGlobalMpcToLocalBlkMpcNumMap();
  if(mpcCpuToBlock)          { delete mpcCpuToBlock           ; mpcCpuToBlock          = 0; }
  if(myCPUToLocAssBlocks)    { delete myCPUToLocAssBlocks     ; myCPUToLocAssBlocks    = 0; }
  if(myCPUToTmpAssBlocks)    { delete myCPUToTmpAssBlocks     ; myCPUToTmpAssBlocks    = 0; }
  if(myCPUExtBlockIdArray)   { delete [] myCPUExtBlockIdArray ; myCPUExtBlockIdArray   = 0; }
  if(LlIdLocAssBlkToMySubs)  { delete LlIdLocAssBlkToMySubs   ; LlIdLocAssBlkToMySubs  = 0; }
  if(LlIdLocAssBlkToSendId)  { delete [] LlIdLocAssBlkToSendId; LlIdLocAssBlkToSendId  = 0; }
  if(myCPUBlkToCpuSharedMpcs){ delete myCPUBlkToCpuSharedMpcs ; myCPUBlkToCpuSharedMpcs= 0; }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::initialize()
{
  nMpcBlocks = 0; nMpcBlocksOnMyCPU = 0;
  blockCCtsolver = 0; 
  blockCCtsparse = 0;
  this->glNumMpc = 0; this->mpcToCpu = 0; blockToMpc = 0; blockToSub = 0; mpcToBlock = 0;
  blockMpcToMpc = 0; blockToCpu = 0; cpuToBlock = 0;
  mpcv = 0; this->numSubsWithMpcs = 0;
  nBigBlocksperMPI       = 0;
  nSmallBlocksperMPI     = 0;
  MPITosuperBlock        = 0;
  blockMpcEqNums         = 0;
  GlMpcToLlBlkMpcNumMap  = 0;
  mpcCpuToBlock          = 0;
  nExtMpcBlocksOnMyCPU   = 0;
  myCPUExtBlockIdArray   = 0;
  nLocAssBlocksOnMyCPU   = 0;
  nTmpAssBlocksOnMyCPU   = 0;
  LlIdLocAssBlkToMySubs  = 0;
  LlIdLocAssBlkToSendId  = 0;
  myCPUToLocAssBlocks    = 0;
  myCPUToTmpAssBlocks    = 0;
  myCPUBlkToCpuSharedMpcs= 0;
  mpcvPat1 = 0; mpcvPat2 = 0; blockCCtPat = 0;
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::assemble()
{
  execParal(nLocAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::assembleOneBlockCCtsolver);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::factor()
{
  if(nMpcBlocksOnMyCPU) {
     // -------------------------------------------------------------------------------------
     // step 6.1 factor "big thread blocks" by (algebraic) parallel factorization
     //          if Skyline selected or concurrent factorization if Sparse (Charbel's advice)
     // -------------------------------------------------------------------------------------
     if(nBigBlocksperMPI[myCPU]){
       if(finfo->cct_cntl->subtype == 0) {
         for(int i=0;i< nBigBlocksperMPI[myCPU]; i++) {
           int iBlock = (*cpuToBlock)[myCPU][i];
           blockCCtsolver[iBlock]->parallelFactor();
           if(blockCCtsolver[iBlock]->numRBM() > 0)
             filePrint(stderr,"WARNING: Number of singularities in CCt block %d = %d\n",
                       iBlock, blockCCtsolver[iBlock]->numRBM());
         }
       }
       else {
         execParal(nBigBlocksperMPI[myCPU], this, &SuperBlockCCtSolver<Scalar>::factorSparseBigBlockCCtsolver);
       }
     }
     // -----------------------------------------------------------------------------------
     // step 6.2 factor "small thread superblocks" by concurrent thread factorization
     //          By construction, number of "small superblocks" <= number of threads
     // -----------------------------------------------------------------------------------
     int nSmallThreadBlock = nMpcBlocksOnMyCPU - nBigBlocksperMPI[myCPU];
     if(nSmallThreadBlock){
       execParal(nSmallThreadBlock, this, &SuperBlockCCtSolver<Scalar>::factorSmallBlockCCtsolver,&nBigBlocksperMPI[myCPU]);
     }
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::zeroAll()
{
  execParal(nLocAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::zeroOneBlockCCtsolver);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::reSolve(GenDistrVector<Scalar> &v)
{
  // HB: the following methods are mainly wrapper methods for each of the associated
  // ...OneBlock... method: this was mainly done to provide different implementation
  // possible to take advantage of the thread load scheduling currently implemented
  initBlockMpcResidual();
  extractBlockMpcResidual(v);
#ifdef DISTRIBUTED
  sendBlockMpcResidualBeforeSolve();
  mpcvPat1->exchange();
  recBlockMpcResidualBeforeSolve();
#endif
  solveBlockCCt(v);
#ifdef DISTRIBUTED
  sendBlockMpcResidualAfterSolve();
  mpcvPat2->exchange();
  recBlockMpcResidualAfterSolve();
#endif
  execParal(this->numSubsWithMpcs, this, &SuperBlockCCtSolver<Scalar>::insertBlockMpcResidual, v);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeSuperBlocks()
{
  double *blockCost      = new double[nMpcBlocks];
  double *blockBandWidth = new double[nMpcBlocks];
  estimateBlockCost(blockCost,blockBandWidth);

  int i;
  std::set<BlockPair> BigBlocks;
  std::set<BlockPair> SmallBlocks;
  std::set<BlockPair> AllBlocks;
  std::set<BlockPair> RemainingBlocks;
  double BigBlockTotalCost = 0.0;
  for(i=0; i<nMpcBlocks; ++i){
    AllBlocks.insert(BlockPair(i,blockCost[i],blockBandWidth[i]));
    if(blockBandWidth[i]*blockBandWidth[i]>1250*threadManager->numThr()){
      BigBlocks.insert(BlockPair(i,blockCost[i],blockBandWidth[i]));
      BigBlockTotalCost += blockCost[i];
    }
    else {
      SmallBlocks.insert(BlockPair(i,blockCost[i],blockBandWidth[i]));
    }
    RemainingBlocks.insert(BlockPair(i,blockCost[i],blockBandWidth[i]));
  }
  int numMPIProcess = numCPUs;

  int iMPI;
  int nsuperBlocks = numMPIProcess;
  if(nMpcBlocks<numMPIProcess) nsuperBlocks = nMpcBlocks;
  int *MPISuperBlockPtr = new int[nsuperBlocks+1];
  MPISuperBlockPtr[0] = 0;
  int *MPISuperBlockMap = new int[nMpcBlocks];
  int nremainingBlocks = nMpcBlocks;
  for(iMPI=0;iMPI<nsuperBlocks;iMPI++){
    MPISuperBlockPtr[iMPI+1] = MPISuperBlockPtr[iMPI];
    double MPICost = 0.0;
    // compute MPI cost average
    double MeanMPICost = 0.0;
    std::set<BlockPair>::iterator Iset = RemainingBlocks.begin();
    while(Iset!=RemainingBlocks.end()){ MeanMPICost += Iset->cost; Iset++;}
    MeanMPICost /= (nsuperBlocks-iMPI);
    std::set<BlockPair>::iterator Imax = RemainingBlocks.end();
    int nBlocks = 0;
    if(Imax!=RemainingBlocks.begin()) Imax--;
    if(Imax->cost>=MeanMPICost){
      // add block to superBlock connectivity
      MPISuperBlockPtr[iMPI+1]++;
      MPISuperBlockMap[MPISuperBlockPtr[iMPI]+nBlocks] = Imax->Id;
      nBlocks++;
      // update current superblock cost
      MPICost += Imax->cost;
      // remove added block from the remaining blocks set
      nremainingBlocks--;
      RemainingBlocks.erase(Imax);
    }
    while((MPICost<MeanMPICost)&&(RemainingBlocks.size()>0)){
      std::set<BlockPair>::iterator itupper = RemainingBlocks.upper_bound(MeanMPICost-MPICost);
      if(itupper==RemainingBlocks.end()) itupper--;
      int blkId = itupper->Id;
      // add block to superBlock connectivity
      MPISuperBlockPtr[iMPI+1]++;
      MPISuperBlockMap[MPISuperBlockPtr[iMPI]+nBlocks] = blkId;
      nBlocks++;
      // update current superblock cost
      MPICost += itupper->cost;
      // remove added block from the remaining blocks set
      nremainingBlocks--;
      RemainingBlocks.erase(itupper);
    }
  }
  Connectivity *superBlockToBlock = new Connectivity(nsuperBlocks,MPISuperBlockPtr,MPISuperBlockMap);
                                          
  int iMPIsuperBlk;
  int iBlk;
                                          
  std::set<BlockPair> ThreadBigBlocks;
  std::set<BlockPair> RemainingThreadSmallBlocks;
  int superBlkcount = 0;
                                          
  ResizeArray<int> *NewSuperBlockPtr   = new ResizeArray<int>(0,0);
  ResizeArray<int> *NewSuperBlockMap   = new ResizeArray<int>(0,0);
  ResizeArray<int> *MPIToSuperBlockPtr  = new ResizeArray<int>(0,0);
  ResizeArray<int> *MPIToSuperBlockMap  = new ResizeArray<int>(0,0);
                                          
  (*NewSuperBlockPtr)[0]   = 0;
  (*MPIToSuperBlockPtr)[0] = 0;
                                          
  int *nSuperBlockperMPI = new int[numMPIProcess];
  nBigBlocksperMPI  = new int[numMPIProcess];
  nSmallBlocksperMPI= new int[numMPIProcess];
  for(i=0;i<numMPIProcess;i++) {
    nSuperBlockperMPI[i] = 0;
    nBigBlocksperMPI[i]  = 0;
    nSmallBlocksperMPI[i]= 0;
  }
                                          
  for(iMPIsuperBlk=0;iMPIsuperBlk<nsuperBlocks;iMPIsuperBlk++){
    (*MPIToSuperBlockPtr)[iMPIsuperBlk+1] = (*MPIToSuperBlockPtr)[iMPIsuperBlk];
                                          
    ThreadBigBlocks.clear();
    RemainingThreadSmallBlocks.clear();
    int nThreads = threadManager->numThr();
    for(iBlk=0;iBlk<superBlockToBlock->num(iMPIsuperBlk);iBlk++){
      int blkId = (*superBlockToBlock)[iMPIsuperBlk][iBlk];
      if(blockBandWidth[blkId]*blockBandWidth[blkId]>1250*nThreads){
        ThreadBigBlocks.insert(BlockPair(blkId,blockCost[blkId],blockBandWidth[blkId]));
      }
      else {
        RemainingThreadSmallBlocks.insert(BlockPair(blkId,blockCost[blkId],blockBandWidth[blkId]));
      }
    }
    nBigBlocksperMPI[iMPIsuperBlk] = ThreadBigBlocks.size();
    int nThreadSmallBlocks         = RemainingThreadSmallBlocks.size();
                              
    std::set<BlockPair>::iterator IBigset = ThreadBigBlocks.begin();
    while(IBigset!=ThreadBigBlocks.end()){
      (*NewSuperBlockPtr)[superBlkcount+1] = (*NewSuperBlockPtr)[superBlkcount];
      (*NewSuperBlockPtr)[superBlkcount+1]++;
      (*NewSuperBlockMap)[(*NewSuperBlockPtr)[superBlkcount]] = IBigset->Id;
      (*MPIToSuperBlockPtr)[iMPIsuperBlk+1]++;
      (*MPIToSuperBlockMap)[(*MPIToSuperBlockPtr)[iMPIsuperBlk]+nSuperBlockperMPI[iMPIsuperBlk]] = superBlkcount;
      superBlkcount++;
      IBigset++;
      nSuperBlockperMPI[iMPIsuperBlk]++;
    }
    int iThread;
    int nThreadsuperBlocks = nThreads;
    if(nThreadsuperBlocks>nThreadSmallBlocks) nThreadsuperBlocks = nThreadSmallBlocks;
    for(iThread=0;iThread<nThreadsuperBlocks;iThread++){
      (*NewSuperBlockPtr)[superBlkcount+1] = (*NewSuperBlockPtr)[superBlkcount];
      double ThreadCost = 0.0;
      // compute thread cost average
      double MeanThreadCost = 0.0;
      std::set<BlockPair>::iterator Iset = RemainingThreadSmallBlocks.begin();
      while(Iset!=RemainingThreadSmallBlocks.end()){ MeanThreadCost += Iset->cost; Iset++;}
      MeanThreadCost /= (nThreadsuperBlocks-iThread);
      int nBlocks = 0;
      if(RemainingThreadSmallBlocks.size()>0) { nSmallBlocksperMPI[iMPIsuperBlk]++; }
      while((ThreadCost<MeanThreadCost)&&(RemainingThreadSmallBlocks.size()>0)){
        std::set<BlockPair>::iterator itupper = RemainingThreadSmallBlocks.upper_bound(MeanThreadCost-ThreadCost);
        if(itupper==RemainingThreadSmallBlocks.end()) itupper--;
        int blkId = itupper->Id;
        // add block to superBlock connectivity
        (*NewSuperBlockPtr)[superBlkcount+1]++;
        (*NewSuperBlockMap)[(*NewSuperBlockPtr)[superBlkcount]+nBlocks] = blkId;
                                          
        nBlocks++;
        // update current superblock cost
        ThreadCost += itupper->cost;
        // remove added block from the remaining blocks set
        RemainingThreadSmallBlocks.erase(itupper);
      }
      (*MPIToSuperBlockPtr)[iMPIsuperBlk+1]++;
      (*MPIToSuperBlockMap)[(*MPIToSuperBlockPtr)[iMPIsuperBlk]+nSuperBlockperMPI[iMPIsuperBlk]] = superBlkcount;
      nSuperBlockperMPI[iMPIsuperBlk]++;
      superBlkcount++;
    }
  }
  for(iMPIsuperBlk=nsuperBlocks+1;iMPIsuperBlk<=numMPIProcess;iMPIsuperBlk++){
    (*MPIToSuperBlockPtr)[iMPIsuperBlk] = (*MPIToSuperBlockPtr)[iMPIsuperBlk-1];
  }
                                          
  Connectivity *NewsuperBlockToBlock = new Connectivity(superBlkcount,NewSuperBlockPtr->data(false),NewSuperBlockMap->data(false));
  MPITosuperBlock = new Connectivity(numMPIProcess,MPIToSuperBlockPtr->data(false),MPIToSuperBlockMap->data(false));
                                          
  // make superBlockToMpc connectivity
  Connectivity *superBlockToMpc = NewsuperBlockToBlock->transcon(blockToMpc);
                                          
  // replace blockToMpc (and associated data) by superBlockToMpc connectivity
  if(blockToMpc) { delete blockToMpc; blockToMpc = 0; }
  if(mpcToBlock) { delete mpcToBlock; mpcToBlock = 0; }
  nMpcBlocks = NewsuperBlockToBlock->csize();
  blockToMpc = superBlockToMpc;
  mpcToBlock = blockToMpc->alloc_reverse();
/*
  for(iMPIsuperBlk=0;iMPIsuperBlk<MPITosuperBlock->csize();iMPIsuperBlk++){
    filePrint(stderr,"  # MPI superblock %d contains %d (super) blocks\n",iMPIsuperBlk,MPITosuperBlock->num(iMPIsuperBlk));
    for(i=0;i<MPITosuperBlock->num(iMPIsuperBlk);i++){
      int blkId = (*MPITosuperBlock)[iMPIsuperBlk][i];
      filePrint(stderr,"  * (super) blocks Id %d made of %d blocks contains %d lmpc eqs\n",
                       blkId, NewsuperBlockToBlock->num(blkId), blockToMpc->num(blkId));
    }
  }
*/
  delete [] blockCost;
  delete [] blockBandWidth;
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::createSuperBlockCCt(const Connectivity *mpcToSub)
{
  int i;
  // ------------------------------------------------------------
  // Step 0. make blockToCpu and cpuToBlock connectivities
  //         make the list of "extended" Mpc Block of this CPU
  // ------------------------------------------------------------
  cpuToBlock = MPITosuperBlock;
  blockToCpu = cpuToBlock->alloc_reverse();
  // nMpcBlocksOnMyCPU always exist even in shared mode <=> 1 MPI
  nMpcBlocksOnMyCPU = cpuToBlock->num(myCPU);
#ifdef DISTRIBUTED
  blockToMpcCpu = std::make_unique<Connectivity>(blockToMpc->transcon(*this->mpcToCpu));
  myCPUToLocAssBlocks = new ResizeArray<int>(0);
  myCPUToTmpAssBlocks = new ResizeArray<int>(0);
  nLocAssBlocksOnMyCPU= 0;
  nTmpAssBlocksOnMyCPU= 0;
  for(i=0;i<nMpcBlocks;i++){
    if(blockToMpcCpu->offset(i,myCPU) !=-1){
      (*myCPUToLocAssBlocks)[nLocAssBlocksOnMyCPU++] = i;
      if((*blockToCpu)[i][0]!=myCPU)
        (*myCPUToTmpAssBlocks)[nTmpAssBlocksOnMyCPU++] = i;
    }
  }
  // make the list of "extended" Mpc Block Id: blocks that need to "done"
  // (either fully or partially) on this CPU
  std::list<int> BlockIdOnMyCPU;
  std::list<int> BlockIdConnectedToMyCPU; // (super) blocks that have lmpcs which belong to this CPU
  for(i=0;i<nMpcBlocksOnMyCPU;i++){
    BlockIdOnMyCPU.insert(BlockIdOnMyCPU.end(),(*cpuToBlock)[myCPU][i]);
  }
  for(i=0;i<nLocAssBlocksOnMyCPU;i++){
    BlockIdConnectedToMyCPU.insert(BlockIdConnectedToMyCPU.end(),(*myCPUToLocAssBlocks)[i]);
  }
  BlockIdOnMyCPU.sort();
  BlockIdConnectedToMyCPU.sort();
  BlockIdConnectedToMyCPU.merge(BlockIdOnMyCPU);
  BlockIdConnectedToMyCPU.unique();
                                          
  nExtMpcBlocksOnMyCPU = BlockIdConnectedToMyCPU.size();
  if(myCPUExtBlockIdArray) { delete [] myCPUExtBlockIdArray; myCPUExtBlockIdArray = 0; }
  myCPUExtBlockIdArray = new int[nExtMpcBlocksOnMyCPU];
  std::list<int>::iterator Ilist = BlockIdConnectedToMyCPU.begin();
  i = 0;
  while(Ilist!=BlockIdConnectedToMyCPU.end()) { myCPUExtBlockIdArray[i++] = *Ilist++; }
#else
  blockToMpcCpu = std::make_unique<Connectivity>(*MPITosuperBlock);
  // make the list of "extended" Mpc Block Id: blocks that need to "done"
  // (either fully or partially) on this CPU
  // HERE it is only a simple copy of cpuToBlock
  nExtMpcBlocksOnMyCPU = nMpcBlocksOnMyCPU;
  if(myCPUExtBlockIdArray) { delete [] myCPUExtBlockIdArray; myCPUExtBlockIdArray = 0; }
  myCPUExtBlockIdArray = new int[nExtMpcBlocksOnMyCPU];
  for(i=0;i<nExtMpcBlocksOnMyCPU;i++){ myCPUExtBlockIdArray[i] = (*cpuToBlock)[myCPU][i]; }
  nLocAssBlocksOnMyCPU = nMpcBlocksOnMyCPU;
  myCPUToLocAssBlocks = new ResizeArray<int>(0);
  for(i=0;i<nLocAssBlocksOnMyCPU;i++){ (*myCPUToLocAssBlocks)[i] =  (*cpuToBlock)[myCPU][i]; }
  nTmpAssBlocksOnMyCPU = 0;
  myCPUToTmpAssBlocks  = 0;
#endif
  // ------------------------------------------------------------
  // Step 1. make mpcToMpc connectivity for each block that is
  // either partially or fully assembled on myCPU -> "extented"
  // block list (see Step 0.)
  // ------------------------------------------------------------
  GlMpcToLlBlkMpcNumMap = 0;
  makeGlobalMpcToLocalBlkMpcNumMap();
  createBlockMpcToMpcConnectivity();
                                          
  // no more need of GlMpcToLlBlkMpcNumMap -> destroy it
  // As the memory has been allocated in // (thread level), it
  // could be more safe to delete it in //
  // (see static scheduling of the current thread implementation)
  deleteGlobalMpcToLocalBlkMpcNumMap();
                                          
  // ------------------------------------------------------------
  // Step 2. make blockToSub and localMpcToBlock connectivities
  // ------------------------------------------------------------
  blockToSub = blockToMpc->transcon(mpcToSub);
  paralApply(this->subsWithMpcs, &FetiBaseSub::setLocalMpcToBlock, mpcToBlock, blockToMpc);
                                          
  // make block local Id to my subs (subs on my CPU) connectivity
  // -> ONLY for those blocks that myCPU has a lmpc CCt contribution to
  //    Those blocks will be referred as "Locally Assembled Blocks"
  //    Note that the blocks that will be factored & solved on myCPU may not
  //    belong to this list (if there is NO CCt contribution from myCPU to them)
  // -> to be used in assembleOneBlockCCtsolver method
  int *subsPtr = new int[nLocAssBlocksOnMyCPU+1];
  ResizeArray<int> subsTarget(0);
  subsPtr[0] = 0;
  int count = 0;
  for(int iBlk=0;iBlk<nLocAssBlocksOnMyCPU;iBlk++){
    subsPtr[iBlk+1] = subsPtr[iBlk];
    int iBlock = (*myCPUToLocAssBlocks)[iBlk];
    for(i=0; i<this->numSubsWithMpcs; ++i) {
      if(blockToSub->offset(iBlock, this->subsWithMpcs[i]->subNum()) != -1){
        subsTarget[count++] = i;
        subsPtr[iBlk+1]++;
      }
    }
  }
  LlIdLocAssBlkToMySubs = new Connectivity(nLocAssBlocksOnMyCPU,subsPtr,subsTarget.data(false));
                                          
#ifdef DISTRIBUTED
  if(nLocAssBlocksOnMyCPU){
    LlIdLocAssBlkToSendId = new int[nLocAssBlocksOnMyCPU];
    for(i=0;i<nLocAssBlocksOnMyCPU;i++){
      int iBlock = (*myCPUToLocAssBlocks)[i];
      LlIdLocAssBlkToSendId[i] = blockToMpcCpu->offset(iBlock, myCPU);
    }
  } else { LlIdLocAssBlkToSendId = 0; }
#else
  LlIdLocAssBlkToSendId = 0;
#endif
                                          
  // ---------------------------------------------------------------------
  // Step 3. create block solvers
  // ---------------------------------------------------------------------
  filePrint(stderr, " ... Making %d MPC (Super) Blocks ... \n", nMpcBlocks);
  createBlockCCtsolver();
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeBlockCCtCommPattern()
{
   int i, j;
   // make cpuToId connectivity (Id's are used to make unique channel number for communication, similar to subdomain id)
   int size = blockToMpcCpu->numConnect() + nMpcBlocks;
   int *target = new int[size];
   int count = 0;
   for(i=0; i<blockToMpcCpu->csize(); ++i) {
     for(j=0; j<blockToMpcCpu->num(i); ++j) {
       int cpu = (*blockToMpcCpu)[i][j];
       target[count++] = cpu;
     }
   }
   for(i=0; i<blockToCpu->csize(); ++i) {
     int cpu = (*blockToCpu)[i][0];
     target[count++] = cpu;
   }
   int *pointer = new int[size+1];
   for(i=0; i<=size; ++i) pointer[i] = i;
   Connectivity *idToCpu = new Connectivity(size, pointer, target);
   Connectivity *cpuToId = idToCpu->alloc_reverse();
   delete idToCpu;
   // blockCCtPat is used to send partially assembled CCt blocks to master CPU for final assembly and factoring
   blockCCtPat = new FSCommPattern<Scalar>(this->fetiCom, cpuToId, myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                           FSCommPattern<Scalar>::NonSym);
   // mpcvPat1 is used to send partially assembled mpcv residuals to master CPU for blockCCt solve
   mpcvPat1 = new FSCommPattern<Scalar>(this->fetiCom, cpuToId, myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                        FSCommPattern<Scalar>::NonSym);
   // mpcvPat2 is used to send solved mpcv residuals from master CPU to other CPUs with involved subds
   mpcvPat2 = new FSCommPattern<Scalar>(this->fetiCom, cpuToId, myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                        FSCommPattern<Scalar>::NonSym);
   delete cpuToId;
   for(i=0; i<nMpcBlocks; ++i) setBlockCCtCommSize(i);
   blockCCtPat->finalize();
   mpcvPat1->finalize();
   mpcvPat2->finalize();
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::setBlockCCtCommSize(int iBlock)
{
  int i, sendId, recId;
  sendId = blockToMpcCpu->offset(iBlock,myCPU);
  if(sendId != -1) { // this means that iBlock is partially assembled on myCPU
    int iCPU = (*blockToCpu)[iBlock][0];  // this is the master CPU that iBlock will be stored & solved on
    recId = iBlock + blockToMpcCpu->numConnect();
    if(iCPU != myCPU) { // need to send iBlock data from myCPU to iCPU
      blockCCtPat->setLen(sendId, recId, blockCCtsolver[iBlock]->size());
      mpcvPat1->setLen(sendId, recId, blockToMpc->num(iBlock));
      mpcvPat2->setLen(sendId, recId, 0);
    }
  }
  if((*blockToCpu)[iBlock][0] == myCPU) {
    for(i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
      int iCPU = (*blockToMpcCpu)[iBlock][i];
      if(iCPU != myCPU) {
        sendId = iBlock + blockToMpcCpu->numConnect();
        recId = blockToMpcCpu->offset(iBlock,iCPU);
        blockCCtPat->setLen(sendId, recId, 0);
        mpcvPat1->setLen(sendId, recId, 0);
        mpcvPat2->setLen(sendId, recId, blockToMpc->num(iBlock));
      }
    }
  }
}


template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::factorSmallBlockCCtsolver(int i, int *blockOffset)
{
  int iBlock = (*cpuToBlock)[myCPU][i+(*blockOffset)];
  blockCCtsolver[iBlock]->factor();
  if(blockCCtsolver[iBlock]->numRBM() > 0)
    filePrint(stderr,"WARNING: Number of singularities in CCt block %d = %d\n",
              iBlock, blockCCtsolver[iBlock]->numRBM());
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::factorSparseBigBlockCCtsolver(int i)
{
  int iBlock = (*cpuToBlock)[myCPU][i];
  blockCCtsolver[iBlock]->factor();
  if(blockCCtsolver[iBlock]->numRBM() > 0)
    filePrint(stderr,"WARNING: Number of singularities in CCt block %d = %d\n",
              iBlock, blockCCtsolver[iBlock]->numRBM());
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::initBlockMpcResidual() const
{
  //int nBlks = nExtMpcBlocksOnMyCPU;
  execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::initOneBlockMpcResidual);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::initOneBlockMpcResidual(int IBlock) const
{
  int iBlock = myCPUExtBlockIdArray[IBlock];
  mpcv[iBlock]->zeroAll();
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::extractBlockMpcResidual(GenDistrVector<Scalar> &v) const
{
  execParal(nLocAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::extractOneBlockMpcResidual, v);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::extractOneBlockMpcResidual(int IBlock, GenDistrVector<Scalar> &v) const
{
  int iBlock = (*myCPUToLocAssBlocks)[IBlock];
  for(int i=0;i<LlIdLocAssBlkToMySubs->num(IBlock);i++) {
    int iSub = (*LlIdLocAssBlkToMySubs)[IBlock][i];
    this->subsWithMpcs[iSub]->extractBlockMpcResidual(iBlock, v.subData(this->subsWithMpcs[iSub]->localSubNum()), 
                                                      mpcv[iBlock], blockMpcEqNums[iBlock]);
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::insertBlockMpcResidual(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->insertBlockMpcResidual(subv, mpcv, mpcToBlock, blockMpcEqNums);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::sendBlockMpcResidualBeforeSolve()
{
  execParal(nLocAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::sendOneBlockMpcResidualBeforeSolve);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::sendOneBlockMpcResidualBeforeSolve(int IBlock)
{
  int iBlock = (*myCPUToLocAssBlocks)[IBlock];
  int sendId = LlIdLocAssBlkToSendId[IBlock];
  int iCPU = (*blockToCpu)[iBlock][0];
  if(iCPU != myCPU) { // send mpcv[iBlock] to iCPU where it will be solved
    int recId = iBlock + blockToMpcCpu->numConnect();
    FSSubRecInfo<Scalar> sInfo = mpcvPat1->getSendBuffer(sendId, recId);
    for(int j=0; j<mpcv[iBlock]->size(); ++j)
      sInfo.data[j] = (*mpcv[iBlock])[j];
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::recBlockMpcResidualBeforeSolve()
{
  execParal(nMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::recOneBlockMpcResidualBeforeSolve);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::recOneBlockMpcResidualBeforeSolve(int IBlock)
{
   // mpcv[iBlock] will be solved on myCPU
  int iBlock = (*cpuToBlock)[myCPU][IBlock];
  for(int i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
    int iCPU = (*blockToMpcCpu)[iBlock][i];
    if(iCPU != myCPU) {  // add contribution to mpcv[iBlock] from iCPU
      int recId = iBlock + blockToMpcCpu->numConnect();
      int sendId = blockToMpcCpu->offset(iBlock) + i;
      FSSubRecInfo<Scalar> rInfo = mpcvPat1->recData(sendId, recId);
      mpcv[iBlock]->add(rInfo.data);
    }
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::solveBlockCCt( GenDistrVector<Scalar> &v) const
{
  execParal(nMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::solveOneBlockCCt, v);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::solveOneBlockCCt(int IBlock, GenDistrVector<Scalar> &v) const
{
  int iBlk = IBlock;
  int iBlock = (*cpuToBlock)[myCPU][iBlk];
#ifdef DISTRIBUTED
  // apply weighting to mpc residual to accounts for values that have been
  // added more than once in globalSum
  // HB: can probably be optimized to avoid the loop on ALL the lmpc in iBlock
  // & checking for its multiplicity: for example by pre-determining the number
  // of "CPU-shared" lmpcs per iBlock on myCPU
  for(int i=0;i<myCPUBlkToCpuSharedMpcs->num(iBlk);i++){
    int impc = (*myCPUBlkToCpuSharedMpcs)[iBlk][i];
    int gi = (*blockToMpc)[iBlock][impc];
    int bi = blockMpcEqNums[iBlock]->firstdof(impc);
    (*mpcv[iBlock])[bi] /= double(this->mpcToCpu->num(gi));
  }
  blockCCtsolver[iBlock]->reSolve(*mpcv[iBlock]);
#else
  blockCCtsolver[iBlock]->reSolve(*mpcv[iBlock]);
#endif
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::sendBlockMpcResidualAfterSolve()
{
  // can be optimized by looping ONLY on the mpc block that really need to be sent
  execParal(nMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::sendOneBlockMpcResidualAfterSolve);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::sendOneBlockMpcResidualAfterSolve(int IBlock)
{
  int iBlock = (*cpuToBlock)[myCPU][IBlock];
  int i, j;
  //int sendId = (*blockToCpu)[iBlock][0];
  if((*blockToCpu)[iBlock][0] == myCPU) {  // mpcv[iBlock] was solved on myCPU
    for(i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
      int iCPU = (*blockToMpcCpu)[iBlock][i];
      if(iCPU != myCPU) {  // send solved mpcv[iBlock] to iCPU where it will be inserted in subs
        int sendId = iBlock + blockToMpcCpu->numConnect();
        int recId = blockToMpcCpu->offset(iBlock) + i; // HB
        FSSubRecInfo<Scalar> sInfo = mpcvPat2->getSendBuffer(sendId, recId);
        for(j=0; j<mpcv[iBlock]->size(); ++j)
          sInfo.data[j] = (*mpcv[iBlock])[j];
      }
    }
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::recBlockMpcResidualAfterSolve()
{
  // can be optimized by looping ONLY on the mpc block that really need to be received
  execParal(nLocAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::recOneBlockMpcResidualAfterSolve);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::recOneBlockMpcResidualAfterSolve(int IBlock)
{
  int iBlock = (*myCPUToLocAssBlocks)[IBlock];
  int recId = LlIdLocAssBlkToSendId[IBlock];
  int iCPU = (*blockToCpu)[iBlock][0];
  if(iCPU != myCPU) { // get solved mpcv[iBlock] from iCPU
    int sendId = iBlock + blockToMpcCpu->numConnect();
    FSSubRecInfo<Scalar> rInfo = mpcvPat2->recData(sendId, recId);
    mpcv[iBlock]->insertData(rInfo.data);
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeGlobalMpcToLocalBlkMpcNumMap()
{
  GlMpcToLlBlkMpcNumMap = new int*[nMpcBlocks];
  for(int i=0;i<nMpcBlocks;i++) { GlMpcToLlBlkMpcNumMap[i] = 0; }
  execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::makeOneGlobalMpcToLocalBlkMpcNumMap);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeOneGlobalMpcToLocalBlkMpcNumMap(int IBlock)
{
  // current implementation use a 2D array (fast but possibke memory issue)
  // could use a STL map (memory friendly but need a search in
  // method createBlockMpcToMpcConnectivity)
  int i;
  int iBlock = myCPUExtBlockIdArray[IBlock];
                                                                                                                                    
  GlMpcToLlBlkMpcNumMap[iBlock] = new int[this->glNumMpc];
  for(i=0;i<this->glNumMpc;i++) { GlMpcToLlBlkMpcNumMap[iBlock][i] = -1; }
  for(i=0;i<blockToMpc->num(iBlock);i++){
    GlMpcToLlBlkMpcNumMap[iBlock][(*blockToMpc)[iBlock][i]] = i;
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteGlobalMpcToLocalBlkMpcNumMap()
{
 if(GlMpcToLlBlkMpcNumMap) {
  execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::deleteOneGlobalMpcToLocalBlkMpcNumMap);
  delete [] GlMpcToLlBlkMpcNumMap; GlMpcToLlBlkMpcNumMap = 0;
 }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteOneGlobalMpcToLocalBlkMpcNumMap(int IBlock)
{
  // current implementation use a 2D array (fast but possibke memory issue)
  // could use a STL map (memory friendly but need a search in
  // method createBlockMpcToMpcConnectivity)
  int iBlock = myCPUExtBlockIdArray[IBlock];
  if(GlMpcToLlBlkMpcNumMap[iBlock]) {
    delete [] GlMpcToLlBlkMpcNumMap[iBlock]; GlMpcToLlBlkMpcNumMap[iBlock] = 0;
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::createBlockMpcToMpcConnectivity()
{
  blockMpcToMpc = new const Connectivity * [nMpcBlocks];
  for(int i=0;i<nMpcBlocks;i++)
    blockMpcToMpc[i] = 0;
  execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::createOneBlockMpcToMpcConnectivity);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::createOneBlockMpcToMpcConnectivity(int IBlock)
{
  int iBlock = myCPUExtBlockIdArray[IBlock];
 
  int j, k;
  int size = blockToMpc->num(iBlock);
  ResizeArray<int> *target = new ResizeArray<int>(0);
  int *pointer = new int[size+1];
  int count = 0;
  for(j = 0; j < size; ++j) {
    pointer[j] = count;
    int gj = (*blockToMpc)[iBlock][j];
    for(k = 0; k < mpcToMpc->num(gj); ++k) {
      int gk = (*mpcToMpc)[gj][k];
      int bk = GlMpcToLlBlkMpcNumMap[iBlock][gk];
      if(bk != -1) (*target)[count++] = bk;
    }
  }
 
  pointer[size] = count;
  blockMpcToMpc[iBlock] = new Connectivity(size, pointer, target->data(false));
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteBlockMpcToMpcConnectivity()
{
  if(blockMpcToMpc) {
    execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::deleteOneBlockMpcToMpcConnectivity);
    delete [] blockMpcToMpc; blockMpcToMpc = 0;
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteOneBlockMpcToMpcConnectivity(int IBlock)
{
  int iBlock = myCPUExtBlockIdArray[IBlock];
  if(blockMpcToMpc[iBlock]) { delete blockMpcToMpc[iBlock]; blockMpcToMpc[iBlock] = 0; }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::createBlockCCtsolver()
{
  blockCCtsolver = new GenSolver<Scalar> * [nMpcBlocks];
  blockCCtsparse = new GenSparseMatrix<Scalar> * [nMpcBlocks];
  blockMpcEqNums = new SimpleNumberer * [nMpcBlocks];

  for(int i=0;i<nMpcBlocks;i++) { 
    blockCCtsolver[i] = 0; 
    blockMpcEqNums[i] = 0; 
  }

  execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::createOneBlockCCtsolver);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::createOneBlockCCtsolver(int IBlock)
{
  int iBlock = myCPUExtBlockIdArray[IBlock];
  int blockSize = blockMpcToMpc[iBlock]->csize();

  if(finfo->cct_cntl->subtype == 0) { // use sloan renumbering for skyline
    compStruct renumber = blockMpcToMpc[iBlock]->renumByComponent(1);
    blockMpcEqNums[iBlock] = new SimpleNumberer(blockSize,renumber.renum);
    delete [] renumber.xcomp;
  }
  else {
    blockMpcEqNums[iBlock] = new SimpleNumberer(blockSize);
  }
  for(int i=0; i<blockSize; ++i) blockMpcEqNums[iBlock]->setWeight(i, 1);
  blockCCtsolver[iBlock] = GenSolverFactory<Scalar>::getFactory()->createSolver(blockMpcToMpc[iBlock], blockMpcEqNums[iBlock],
                                                                                *finfo->cct_cntl, blockCCtsparse[iBlock], 0);
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteBlockCCtsolver()
{
  if(blockCCtsolver) {
    execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::deleteOneBlockCCtsolver, myCPUExtBlockIdArray);
    delete [] blockCCtsolver; blockCCtsolver = 0;
    delete [] blockCCtsparse; blockCCtsparse = 0;
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteMyCPUTmpBlockCCtsolver()
{
  if(blockCCtsolver) {
    execParal(nTmpAssBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::deleteOneBlockCCtsolver, myCPUToTmpAssBlocks->data());
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteOneBlockCCtsolver(int IBlock, int *gBlockIdArray)
{
  int iBlock = gBlockIdArray[IBlock];
  if(blockCCtsparse[iBlock] != NULL && blockCCtsparse[iBlock] != dynamic_cast<GenSparseMatrix<Scalar> *>(blockCCtsolver[iBlock])) delete blockCCtsparse[iBlock];
  if(blockCCtsolver[iBlock]) { delete blockCCtsolver[iBlock]; blockCCtsolver[iBlock] = 0; }
} 

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteBlockMpcEqNums()
{
  if(blockMpcEqNums){
    execParal(nExtMpcBlocksOnMyCPU, this, &SuperBlockCCtSolver<Scalar>::deleteOneBlockMpcEqNums, myCPUExtBlockIdArray);
    delete [] blockMpcEqNums; blockMpcEqNums = 0; 
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::deleteOneBlockMpcEqNums(int IBlock, int *gBlockIdArray)
{
  int iBlock = gBlockIdArray[IBlock];
  if(blockMpcEqNums[iBlock]) { delete blockMpcEqNums[iBlock]; blockMpcEqNums[iBlock] = 0; }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::assembleOneBlockCCtsolver(int IBlock)
{
  int iBlock = (*myCPUToLocAssBlocks)[IBlock];
  for(int i=0;i<LlIdLocAssBlkToMySubs->num(IBlock);i++){
    int iSub = (*LlIdLocAssBlkToMySubs)[IBlock][i];
    this->subsWithMpcs[iSub]->assembleBlockCCtsolver(iBlock, blockCCtsolver[iBlock], blockMpcEqNums[iBlock]);   
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::zeroOneBlockCCtsolver(int IBlock)
{
  int iBlock = (*myCPUToLocAssBlocks)[IBlock];
  blockCCtsolver[iBlock]->zeroAll();
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::sendOneBlockCCtsolver(int IBlock)
{
  int iBlock = (*myCPUToLocAssBlocks)[IBlock];
  int sendId = LlIdLocAssBlkToSendId[IBlock];
  if(sendId != -1) { // iBlock was partially or fully assembled on myCPU 
    int iCPU = (*blockToCpu)[iBlock][0];
    if(iCPU != myCPU) { // iBlock is NOT factored & solved on myCPU -> send its data to iCPU 
      int recId = iBlock + blockToMpcCpu->numConnect();
      FSSubRecInfo<Scalar> sInfo = blockCCtPat->getSendBuffer(sendId, recId);
      for(int j=0; j<blockCCtsolver[iBlock]->size(); ++j)
        sInfo.data[j] = blockCCtsolver[iBlock]->getData()[j];
    }
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::recOneBlockCCtsolver(int IBlock)
{
  // iBlock will be factored & solved on myCPU
  int iBlock = (*cpuToBlock)[myCPU][IBlock];
  int i;
  for(i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
    int iCPU = (*blockToMpcCpu)[iBlock][i];
    if(iCPU != myCPU) { // iBlock needs the lmpc CCt contribution from iCPU 
      int recId = iBlock + blockToMpcCpu->numConnect();
      int sendId = blockToMpcCpu->offset(iBlock) + i; 
      FSSubRecInfo<Scalar> rInfo = blockCCtPat->recData(sendId, recId);
      blockCCtsolver[iBlock]->add(rInfo.data.data());
    }
  }
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeBlkToCpuSharedMpcsMap()
{
  //HB: the current implementation is sequential by using a Connecticity 
  //    object. But this could be done concurrently at the thread level for
  //    each block on myCPU. For that a **ReziseArray implementation could be 
  //    used (for allowing each block/thread fill its own maps)
  int impc, iBlk;
  int *pointer = new int[nMpcBlocksOnMyCPU+1];
  ResizeArray<int> target(0);
  pointer[0] = 0;
  int count = 0;
  for(iBlk=0;iBlk<nMpcBlocksOnMyCPU;iBlk++){
    pointer[iBlk+1] = pointer[iBlk];
    int gBlkId = (*cpuToBlock)[myCPU][iBlk];
    for(impc=0; impc<blockToMpc->num(gBlkId);impc++){
      int gi = (*blockToMpc)[gBlkId][impc];
      if(this->mpcToCpu->num(gi) > 1){
        target[count++] = impc;
        pointer[iBlk+1]++;
      }
    }
  }
  myCPUBlkToCpuSharedMpcs = new Connectivity(nMpcBlocksOnMyCPU, pointer, target.data(false));
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeTmpBlockMpcToMpcConnectivity()
{
  nExtMpcBlocksOnMyCPU = nMpcBlocks;
  myCPUExtBlockIdArray = new int[nMpcBlocks];
  for(int i=0;i<nMpcBlocks;i++) { myCPUExtBlockIdArray[i] = i; }
  GlMpcToLlBlkMpcNumMap = 0;
  makeGlobalMpcToLocalBlkMpcNumMap();
  createBlockMpcToMpcConnectivity();
  deleteGlobalMpcToLocalBlkMpcNumMap();
}

template<class Scalar>
void
SuperBlockCCtSolver<Scalar>::makeSubBlocks()
{
  int i,j,k;
  // ------------------------------------------------------------
  // Step 1. make mpcToMpc connectivity for each block 
  // ------------------------------------------------------------
  makeTmpBlockMpcToMpcConnectivity(); 

  // --------------------------------------------------------------------------------------------
  // Step 2. look for independent components (i.e sub-blocks) in each block mpcToMpc connectivity
  //         Note that this could be done in // (thread level) for each block
  // --------------------------------------------------------------------------------------------
  ResizeArray<int> *target = new ResizeArray<int>(0, this->glNumMpc);
  ResizeArray<int> *pointer = new ResizeArray<int>(0, nMpcBlocks);
  (*pointer)[0] = 0;
  int blockCount = 0;
  int count = 0;
  for(i=0; i<nMpcBlocks; ++i) {
    compStruct renumber = blockMpcToMpc[i]->renumByComponent(-1); // split into component but not renumber yet
    renumber.order = new int[blockToMpc->num(i)];
    for(j=0; j<blockToMpc->num(i); ++j) renumber.order[renumber.renum[j]] = j;
    for(j=0; j<renumber.numComp; ++j) {
      for(k=renumber.xcomp[j]; k<renumber.xcomp[j+1]; ++k)
        (*target)[count++] = (*blockToMpc)[i][renumber.order[k]];
      (*pointer)[++blockCount] = count;
    }
    renumber.clearMemory();
  }
  deleteBlockMpcToMpcConnectivity();
  if(myCPUExtBlockIdArray) { delete [] myCPUExtBlockIdArray; myCPUExtBlockIdArray = 0; }
  if(blockCount > nMpcBlocks) {
    delete mpcToBlock; delete blockToMpc;
    blockToMpc = new Connectivity(blockCount, pointer->data(), target->data());
    nMpcBlocks = blockCount;
  }
}

template<class Scalar>
double
SuperBlockCCtSolver<Scalar>::estimateBlockCost(double *blockCost, double *blockBandWidth)
{
   // HB: estimate cost of each block
   int i,iBlk;
   double totalCost = 0; 
   if(!blockMpcToMpc) makeTmpBlockMpcToMpcConnectivity();
   for(iBlk=0;iBlk<nMpcBlocks;iBlk++){
     SimpleNumberer MpcEqNums(blockToMpc->num(iBlk));       
     for(i=0;i<blockToMpc->num(iBlk);i++) MpcEqNums.setWeight(i, 1);
     MpcEqNums.makeOffset();
     totalCost += blockMpcToMpc[iBlk]->estimateCost(&MpcEqNums,blockCost[iBlk],blockBandWidth[iBlk]);
   }
   deleteBlockMpcToMpcConnectivity(); // no more need of the "original" blockMpcToMpc connectivity
   delete [] myCPUExtBlockIdArray; myCPUExtBlockIdArray = 0;
   return totalCost;
}

template class SuperBlockCCtSolver<double>;
template class SuperBlockCCtSolver<std::complex<double>>;
