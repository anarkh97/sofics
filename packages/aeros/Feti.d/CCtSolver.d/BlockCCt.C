#include <complex>

template <typename S>
class GenDistrVector;

class FSCommunicator;

#include <Utils.d/Connectivity.h>
#include <Driver.d/SubDomain.h>
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#include<Feti.d/CCtSolver.d/BlockCCt.h>
#include <Utils.d/DistHelper.h>
#include <Solvers.d/SolverFactory.h>

template<class Scalar>
BlockCCtSolver<Scalar>::BlockCCtSolver(const Connectivity *_blockToMpc, const Connectivity *mpcToMpc,
                                       const Connectivity *mpcToSub,
                                       const Connectivity *_mpcToCpu, int _numSubsWithMpcs,
                                       std::vector<FetiSub<Scalar> *> subsWithMpcs,
                                       int *_subMap, FetiInfo *_finfo, FSCommunicator *_fetiCom)
    : CCtSolver<Scalar>(std::move(subsWithMpcs))
{
  filePrint(stderr," ... Building block CC^t for preconditioning MPCs ...\n");
  blockToMpc = _blockToMpc;
  nMpcBlocks = blockToMpc->csize();
  this->mpcToCpu = _mpcToCpu;
  this->numSubsWithMpcs = _numSubsWithMpcs;
  subMap = _subMap;
  this->fetiCom = _fetiCom;
  this->glNumMpc = mpcToMpc->csize();
  finfo = _finfo;
  myCPU = this->fetiCom->cpuNum();
  numCPUs = this->fetiCom->size();
                                    
  // Step 1. make mpcToMpc connectivity for each block
  blockMpcToMpc.resize(nMpcBlocks);
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::createBlockMpcToMpcConnectivity, mpcToMpc);
                                    
  // Step 2. renumber each block mpcToMpc connectivity, and if possible split blocks into smaller blocks
  if(finfo->mpc_block == FetiInfo::mortarBlock) {
    int i,j,k;
    ResizeArray<int> *target = new ResizeArray<int>(0, this->glNumMpc);
    ResizeArray<int> *pointer = new ResizeArray<int>(0, nMpcBlocks);
    (*pointer)[0] = 0;
    int blockCount = 0;
    int count = 0;
    for(i=0; i<nMpcBlocks; ++i) {
      compStruct renumber = blockMpcToMpc[i]->renumByComponent(0); // split into component but not renumber yet
      renumber.order = new int[blockToMpc->num(i)];
      for(j=0; j<blockToMpc->num(i); ++j) renumber.order[renumber.renum[j]] = j;
      for(j=0; j<renumber.numComp; ++j) {
        for(k=renumber.xcomp[j]; k<renumber.xcomp[j+1]; ++k)
          (*target)[count++] = (*blockToMpc)[i][renumber.order[k]];
        (*pointer)[++blockCount] = count;
      }
      renumber.clearMemory();
    }
    if(blockCount > nMpcBlocks) {
      delete blockToMpc;
      blockToMpc = new Connectivity(blockCount, pointer->data(), target->data());
      execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::deleteBlockMpcToMpcConnectivity);
      nMpcBlocks = blockCount;
      blockMpcToMpc.resize(nMpcBlocks);
      execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::createBlockMpcToMpcConnectivity, mpcToMpc);
    }
  }

  // Step 3. make blockToSub and localMpcToBlock connectivities
  mpcToBlock = blockToMpc->alloc_reverse();
  blockToSub = blockToMpc->transcon(mpcToSub);
  paralApply(this->subsWithMpcs, &FetiBaseSub::setLocalMpcToBlock, mpcToBlock, blockToMpc);
                                    
  // Step 4. make blockToCpu and cpuToBlock connectivities
#ifdef DISTRIBUTED
  distributeMpcBlocks();
  nMpcBlocksOnMyCPU = cpuToBlock->num(myCPU);
#else
  nMpcBlocksOnMyCPU = nMpcBlocks;
#endif
                                    
  // Step 5. construct block solvers
  filePrint(stderr, " ... Making %3d MPC Blocks         ... \n", nMpcBlocks);
  blockCCtsolver.resize(nMpcBlocks);
  blockCCtsparse.resize(nMpcBlocks);
  blockMpcEqNums.resize(nMpcBlocks);
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::createBlockCCtsolver);

  // Step 6. allocate memory for mpcv (residual vectors for each block represented on myCPU)
  mpcv.resize(nMpcBlocks);
#ifdef DISTRIBUTED
  for(int i=0; i<nMpcBlocks; ++i) {
    if(((*blockToCpu)[i][0] == myCPU) || (blockToMpcCpu->offset(i, myCPU) != -1))
      mpcv[i] = new GenVector<Scalar>(blockToMpc->num(i));
    else
      mpcv[i] = 0;
  }
  makeBlockCCtCommPattern();
#else
  for(int i=0; i<nMpcBlocks; ++i)
     mpcv[i] = new GenVector<Scalar>(blockToMpc->num(i));
#endif
}

template<class Scalar>
BlockCCtSolver<Scalar>::~BlockCCtSolver()
{
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::deleteBlockMpcToMpcConnectivity);
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::deleteBlockCCtsolver);
  delete blockToMpc;
  delete mpcToBlock;
  delete blockToSub;
#ifdef DISTRIBUTED
  delete blockToMpcCpu;
  delete blockToCpu;
  delete cpuToBlock;
  delete blockCCtPat;
  delete mpcvPat1;
  delete mpcvPat2;
#endif
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::assemble()
{
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::assembleBlockCCtsolver);
#ifdef DISTRIBUTED
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::sendBlockCCtsolver);
  blockCCtPat->exchange();
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::recBlockCCtsolver);
#endif
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::factor()
{
  if(nMpcBlocksOnMyCPU < threadManager->numThr()) {
    for(int i = 0; i < nMpcBlocksOnMyCPU; ++i) {
#ifdef DISTRIBUTED
      int iBlock = (*cpuToBlock)[myCPU][i];
#else
      int iBlock = i;
#endif
      blockCCtsolver[iBlock]->parallelFactor();
      if(blockCCtsolver[iBlock]->numRBM() > 0)
        filePrint(stderr," *** WARNING: Number of singularities in CCt block %d = %d\n",
                  iBlock, blockCCtsolver[iBlock]->numRBM());
    }
  }
  else execParal(nMpcBlocksOnMyCPU, this, &BlockCCtSolver<Scalar>::factorBlockCCtsolver);
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::zeroAll()
{
  execParal(nMpcBlocksOnMyCPU, this, &BlockCCtSolver<Scalar>::zeroBlockCCtsolver);
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::reSolve(GenDistrVector<Scalar> &v)
{
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::initBlockMpcResidual);
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::extractBlockMpcResidual, v);
#ifdef DISTRIBUTED
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::sendBlockMpcResidualBeforeSolve);
  mpcvPat1->exchange();
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::recBlockMpcResidualBeforeSolve);
#endif
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::solveBlockCCt, v);
#ifdef DISTRIBUTED
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::sendBlockMpcResidualAfterSolve);
  mpcvPat2->exchange();
  execParal(nMpcBlocks, this, &BlockCCtSolver<Scalar>::recBlockMpcResidualAfterSolve);
#endif
  execParal(this->numSubsWithMpcs, this, &BlockCCtSolver<Scalar>::insertBlockMpcResidual, v);
}



template<class Scalar>
void
BlockCCtSolver<Scalar>::createBlockMpcToMpcConnectivity(int iBlock, const Connectivity *mpcToMpc)
{
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
      int bk = blockToMpc->cOffset(iBlock,gk);
      if(bk != -1) (*target)[count++] = bk;
    }
  }
  pointer[size] = count;
  //int *copy = new int[count]; for(i=0; i<count; ++i) copy[i] = (*target)[i]; delete target;
  //blockMpcToMpc[iBlock] = new Connectivity(size, pointer, copy);
  // HB: use new ResizeArray data() method
  blockMpcToMpc[iBlock] = new Connectivity(size, pointer, target->data(false));
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::deleteBlockMpcToMpcConnectivity(int iBlock)
{
  delete blockMpcToMpc[iBlock];
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::distributeMpcBlocks()
{
  int i, j;
  blockToMpcCpu = blockToMpc->transcon(this->mpcToCpu);
  // make a new blockToCpu (distribute blocks among cpus as evenly as possible )
  // 1st attempt, assumes all blocks are same size (can be improved, taking into account block size for better load balancing)
  std::vector<int> numBlocksOnCpu{numCPUs, 0};
  int *blockCpu = new int[nMpcBlocks]; for(i=0; i<nMpcBlocks; ++i) blockCpu[i] = -1; // -1 means no cpu assigned to this block
  int minBlocksOnCpu = nMpcBlocks / numCPUs;
  int remainder = nMpcBlocks % numCPUs;
  int maxBlocksOnCpu = remainder ? (minBlocksOnCpu + 1) : minBlocksOnCpu;
  int numAssigned = 0;
  int num = 1; int count = 0; int remainderCount = 0;
  while(count < blockToMpcCpu->numConnect()) {
    for(i=0; i<nMpcBlocks; ++i) {
      if(blockToMpcCpu->num(i) == num) {
        count += num;
        if(blockCpu[i] == -1) {
          for(j=0; j<num; ++j) {
            int cpu = (*blockToMpcCpu)[i][j];
            if(numBlocksOnCpu[cpu] < minBlocksOnCpu) { // nice match, assigning block connected cpu
              blockCpu[i] = cpu;
              numBlocksOnCpu[cpu]++;
              numAssigned++;
              break;
            }
            else if((numBlocksOnCpu[cpu] < maxBlocksOnCpu) && (remainderCount < remainder)) { // also nice match
              blockCpu[i] = cpu;
              numBlocksOnCpu[cpu]++;
              numAssigned++;
              remainderCount++;
              break;
            }
            else {
              // can't assign block to connected cpu, leave for 2nd pass
            }
          }
        }
      }
    }
    num++;
  }
  if(numAssigned < nMpcBlocks) { // 2nd pass, allocate unmatched blocks to cpu with avaiable slot
    for(i=0; i<nMpcBlocks; ++i) {
      if(blockCpu[i] == -1) {
        for(j=0; j<numCPUs; ++j) {
          if(numBlocksOnCpu[j] < minBlocksOnCpu) {
            blockCpu[i] = j;
            numBlocksOnCpu[j]++;
            numAssigned++;
            break;
          }
          else if((numBlocksOnCpu[j] < maxBlocksOnCpu) && (remainderCount < remainder)) {
            blockCpu[i] = j;
            numBlocksOnCpu[j]++;
            numAssigned++;
            remainderCount++;
            break;
          }
        }
      }
    }
  }
  int *pointer = new int[nMpcBlocks + 1];
  for(i = 0; i <= nMpcBlocks; ++i) pointer[i] = i;
  blockToCpu = new Connectivity(nMpcBlocks, pointer, blockCpu);
  cpuToBlock = blockToCpu->alloc_reverse();
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::createBlockCCtsolver(int iBlock)
{
#ifdef DISTRIBUTED
  // only allocate memory for this block's solver if it is necessary
  if(((*blockToCpu)[iBlock][0] == myCPU) || (blockToMpcCpu->offset(iBlock,myCPU) != -1)) {
#endif
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
#ifdef DISTRIBUTED
  } else blockCCtsolver[iBlock] = 0;
#endif
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::assembleBlockCCtsolver(int iBlock)
{
  int i;
#ifdef DISTRIBUTED
  if(blockToMpcCpu->offset(iBlock,myCPU) != -1) {
    for(i = 0; i < this->numSubsWithMpcs; ++i)
      if(blockToSub->offset(iBlock, this->subsWithMpcs[i]->subNum()) != -1)
        this->subsWithMpcs[i]->assembleBlockCCtsolver(iBlock, blockCCtsolver[iBlock], blockMpcEqNums[iBlock]);
  }
#else
  for(i = 0; i < blockToSub->num(iBlock); ++i) {
    int iSub = subMap[(*blockToSub)[iBlock][i]];
    this->subsWithMpcs[iSub]->assembleBlockCCtsolver(iBlock, blockCCtsolver[iBlock],blockMpcEqNums[iBlock]);
  }
#endif
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::makeBlockCCtCommPattern()
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
BlockCCtSolver<Scalar>::setBlockCCtCommSize(int iBlock)
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
BlockCCtSolver<Scalar>::sendBlockCCtsolver(int iBlock)
{
  int sendId = blockToMpcCpu->offset(iBlock, myCPU);  // iBlock was partially assembled on myCPU
  if(sendId != -1) {
    int iCPU = (*blockToCpu)[iBlock][0];
    if(iCPU != myCPU) {
      int recId = iBlock + blockToMpcCpu->numConnect();
      FSSubRecInfo<Scalar> sInfo = blockCCtPat->getSendBuffer(sendId, recId);
      for(int j=0; j<blockCCtsolver[iBlock]->size(); ++j)
        sInfo.data[j] = blockCCtsolver[iBlock]->getData()[j];
    }
  }
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::recBlockCCtsolver(int iBlock)
{
  int i;
  if((*blockToCpu)[iBlock][0] == myCPU) { // iBlock will be factored & solved on myCPU
    for(i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
      int iCPU = (*blockToMpcCpu)[iBlock][i];
      if(iCPU != myCPU) {
        int recId = iBlock + blockToMpcCpu->numConnect();
        int sendId = blockToMpcCpu->offset(iBlock,iCPU);
        FSSubRecInfo<Scalar> rInfo = blockCCtPat->recData(sendId, recId);
        blockCCtsolver[iBlock]->add(rInfo.data.data());
      }
    }
  }
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::factorBlockCCtsolver(int i)
{
#ifdef DISTRIBUTED
  int iBlock = (*cpuToBlock)[myCPU][i];
#else
  int iBlock = i;
#endif
  blockCCtsolver[iBlock]->factor();
  if(blockCCtsolver[iBlock]->numRBM() > 0)
    filePrint(stderr," *** WARNING: Number of singularities in CCt block %d = %d\n",
              iBlock, blockCCtsolver[iBlock]->numRBM());
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::zeroBlockCCtsolver(int i)
{
#ifdef DISTRIBUTED
  int iBlock = (*cpuToBlock)[myCPU][i];
#else
  int iBlock = i;
#endif
  blockCCtsolver[iBlock]->zeroAll();
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::deleteBlockCCtsolver(int iBlock)
{
#ifdef DISTRIBUTED
  if((*blockToCpu)[iBlock][0] == myCPU)
#endif
  {
    if(blockCCtsparse[iBlock] != 0 && blockCCtsparse[iBlock] != dynamic_cast<GenSparseMatrix<Scalar> *>(blockCCtsolver[iBlock])) {
      delete blockCCtsparse[iBlock]; blockCCtsparse[iBlock] = 0;
    }
    if(blockCCtsolver[iBlock] != 0) {
      delete blockCCtsolver[iBlock]; blockCCtsolver[iBlock] = 0;
    }
    if(mpcv[iBlock] != 0) {
      delete mpcv[iBlock]; mpcv[iBlock] = 0;
    }
    if(blockMpcEqNums[iBlock] != 0) {
      delete blockMpcEqNums[iBlock]; blockMpcEqNums[iBlock] = 0; 
    }
  }
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::initBlockMpcResidual(int iBlock) const
{
#ifdef DISTRIBUTED
  if(mpcv[iBlock])
#endif
    mpcv[iBlock]->zeroAll();
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::extractBlockMpcResidual(int iBlock, GenDistrVector<Scalar> &v) const
{
  int i;
#ifdef DISTRIBUTED
  for(i = 0; i < this->numSubsWithMpcs; ++i) {
    if(blockToSub->offset(iBlock, this->subsWithMpcs[i]->subNum()) != -1){
      this->subsWithMpcs[i]->extractBlockMpcResidual(iBlock, v.subData(this->subsWithMpcs[i]->localSubNum()), mpcv[iBlock],
                                           blockMpcEqNums[iBlock]);
    }
  }
#else
  for(i = 0; i < blockToSub->num(iBlock); ++i) {
    int iSub = subMap[(*blockToSub)[iBlock][i]];
    this->subsWithMpcs[iSub]->extractBlockMpcResidual(iBlock, v.subData(this->subsWithMpcs[iSub]->localSubNum()), mpcv[iBlock],
                                            blockMpcEqNums[iBlock]);
  }
#endif
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::sendBlockMpcResidualBeforeSolve(int iBlock)
{
  int sendId = blockToMpcCpu->offset(iBlock,myCPU);  // mpcv[iBlock] was partially assembled on myCPU
  if(sendId != -1) {
    int iCPU = (*blockToCpu)[iBlock][0];
    if(iCPU != myCPU) { // send mpcv[iBlock] to iCPU where it will be solved
      int recId = iBlock + blockToMpcCpu->numConnect();
      FSSubRecInfo<Scalar> sInfo = mpcvPat1->getSendBuffer(sendId, recId);
      for(int j=0; j<mpcv[iBlock]->size(); ++j)
        sInfo.data[j] = (*mpcv[iBlock])[j];
    }
  }
}
                                                                                                                                        
template<class Scalar>
void
BlockCCtSolver<Scalar>::recBlockMpcResidualBeforeSolve(int iBlock)
{
  if((*blockToCpu)[iBlock][0] == myCPU) {  // mpcv[iBlock] will be solved on myCPU
    for(int i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
      int iCPU = (*blockToMpcCpu)[iBlock][i];
      if(iCPU != myCPU) {  // add contribution to mpcv[iBlock] from iCPU
        int recId = iBlock + blockToMpcCpu->numConnect();
        //int sendId = blockToMpcCpu->offset(iBlock,iCPU);
        int sendId = blockToMpcCpu->offset(iBlock) + i; // HB
        FSSubRecInfo<Scalar> rInfo = mpcvPat1->recData(sendId, recId);
        mpcv[iBlock]->add(rInfo.data);
      }
    }
  }
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::solveBlockCCt(int iBlock, GenDistrVector<Scalar> &v) const
{
#ifdef DISTRIBUTED
  // apply weighting to mpc residual to accounts for values that have been added more than once in globalSum
  if((*blockToCpu)[iBlock][0] == myCPU) {
    for(int i = 0; i < blockToMpc->num(iBlock); ++i) {
      int gi = (*blockToMpc)[iBlock][i];
      int bi = blockMpcEqNums[iBlock]->firstdof(i);
      if(this->mpcToCpu->num(gi) > 1)
        //(*mpcv[iBlock])[i] /= double(mpcToCpu->num(gi));
        (*mpcv[iBlock])[bi] /= double(this->mpcToCpu->num(gi));
    }
    blockCCtsolver[iBlock]->reSolve(*mpcv[iBlock]);
  }
#else
  blockCCtsolver[iBlock]->reSolve(*mpcv[iBlock]);
#endif
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::sendBlockMpcResidualAfterSolve(int iBlock)
{
  int i, j;
  //int sendId = (*blockToCpu)[iBlock][0];
  if((*blockToCpu)[iBlock][0] == myCPU) {  // mpcv[iBlock] was solved on myCPU
    for(i = 0; i < blockToMpcCpu->num(iBlock); ++i) {
      int iCPU = (*blockToMpcCpu)[iBlock][i];
      if(iCPU != myCPU) {  // send solved mpcv[iBlock] to iCPU where it will be inserted in subs
        int sendId = iBlock + blockToMpcCpu->numConnect();
        //int recId = blockToMpcCpu->offset(iBlock,iCPU);
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
BlockCCtSolver<Scalar>::recBlockMpcResidualAfterSolve(int iBlock)
{
  int recId = blockToMpcCpu->offset(iBlock,myCPU);
  if(recId != -1) {
    int iCPU = (*blockToCpu)[iBlock][0];
    if(iCPU != myCPU) { // get solved mpcv[iBlock] from iCPU
      int sendId = iBlock + blockToMpcCpu->numConnect();
      FSSubRecInfo<Scalar> rInfo = mpcvPat2->recData(sendId, recId);
      mpcv[iBlock]->insertData(rInfo.data);
    }
  }
}

template<class Scalar>
void
BlockCCtSolver<Scalar>::insertBlockMpcResidual(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->insertBlockMpcResidual(subv, mpcv.data(), mpcToBlock, blockMpcEqNums.data());
}

template class BlockCCtSolver<double>;
template class BlockCCtSolver<std::complex<double>>;