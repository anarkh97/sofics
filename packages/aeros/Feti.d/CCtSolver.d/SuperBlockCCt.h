#ifndef _NEWBLOCKCCT_H_
#define _NEWBLOCKCCT_H_

#include <Feti.d/CCtSolver.d/CCtSolver.h>
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/Communicator.h>

class FetiInfo;
template<class Scalar> 
class SuperBlockCCtSolver : public CCtSolver<Scalar>
{
  public:
    SuperBlockCCtSolver(const Connectivity *blockToMpc, const Connectivity *mpcToMpc, const Connectivity *mpcToSub,
                        const Connectivity *mpcToCpu, int numSubsWithMpcs, std::vector<FetiSub<Scalar> *> subsWithMpcs,
                        FetiInfo *finfo, FSCommunicator *fetiCom, bool super_flag = true, bool sub_flag = false);
    ~SuperBlockCCtSolver();
    void reSolve(GenDistrVector<Scalar> &v);
    void zeroAll();
    void assemble();
    void factor();

  private:
    FetiInfo *finfo;
    int nMpcBlocks;         // total nb of blocks
    int nMpcBlocksOnMyCPU;
    GenSolver<Scalar> **blockCCtsolver;// array[nMpcBlocks] of pointer on the BlockCCtsolver of myCPU
    GenSparseMatrix<Scalar> **blockCCtsparse;
    const Connectivity *mpcToMpc;
    const Connectivity *blockToMpc;
    const Connectivity *blockToSub;
    const Connectivity *mpcToBlock;
    const Connectivity **blockMpcToMpc;
    std::unique_ptr<Connectivity> blockToMpcCpu;
    const Connectivity *blockToCpu;
    const Connectivity *cpuToBlock;
    GenVector<Scalar> **mpcv;
    FSCommPattern<Scalar> *blockCCtPat;
    FSCommPattern<Scalar> *mpcvPat1, *mpcvPat2;
    int myCPU, numCPUs;

    int *nBigBlocksperMPI;                // nb of "big" blocks on myCPU
    int *nSmallBlocksperMPI;              // nb of "small" blocks on myCPU
    const Connectivity *MPITosuperBlock;        // for each CPU give the (super) blocks Id they will store & solve
    const Connectivity *mpcCpuToBlock;          // for each CPU, give the lmpc block whose have
                                          //    lmpc CCt contribution from this CPU   
    SimpleNumberer **blockMpcEqNums;      // array[nMpcBlock] of pointer to the EqNum of each blocks created on myCPU 
    int **GlMpcToLlBlkMpcNumMap;          // array[nMpcBlock] of pointer to the lmpc global Id to local block Id mapping 
                                          // of each blocks created on myCPU 
    int nLocAssBlocksOnMyCPU;             // nb of blocks assembled (partially of fully) on myCPU 
    int nTmpAssBlocksOnMyCPU;             // nb of blocks PARTIALLY assembled on myCPU 
    int nExtMpcBlocksOnMyCPU;             // nb of blocks created (i.e. associated with a CCt solver) on myCPU
    ResizeArray<int> *myCPUToLocAssBlocks;// list of blocks (global Id) assembled (partially of fully) on myCPU 
    ResizeArray<int> *myCPUToTmpAssBlocks;// list of blocks (global Id) PARTIALLY assembled on myCPU 
    int *myCPUExtBlockIdArray;            // list of blocks created on myCPU
    Connectivity *LlIdLocAssBlkToMySubs;  // map a locally assembled block (local) Id to the subs of myCPU that contribute to it 
    int *LlIdLocAssBlkToSendId;           // map a locally assembled block (local) Id to their sending (communication) Id
    Connectivity *myCPUBlkToCpuSharedMpcs;// gives the "Cpu-shared" lmpc eqs for each blocks stored & solved on myCPU  
 
    void initialize();
    void makeTmpBlockMpcToMpcConnectivity();
    void makeSubBlocks();
    double estimateBlockCost(compStruct &renumber, double *blockCost, double *blockBandWidth);
    double estimateBlockCost(double *blockCost, double *blockBandWidth);
    void makeSuperBlocks();
    void createSuperBlockCCt(const Connectivity *mpcToSub);
    void makeBlockCCtCommPattern();
    void setBlockCCtCommSize(int iBlock);
    void factorSmallBlockCCtsolver(int iBlock, int *blockOffset);
    void factorSparseBigBlockCCtsolver(int iBlock);
    void createBlockMpcToMpcConnectivity();
    void createOneBlockMpcToMpcConnectivity(int IBlock);
    void makeGlobalMpcToLocalBlkMpcNumMap();
    void makeOneGlobalMpcToLocalBlkMpcNumMap(int IBlock);
    void createBlockCCtsolver();
    void createOneBlockCCtsolver(int IBlock);
    void assembleOneBlockCCtsolver(int IBlock);
    void zeroOneBlockCCtsolver(int IBlock);
    void deleteGlobalMpcToLocalBlkMpcNumMap();
    void deleteOneGlobalMpcToLocalBlkMpcNumMap(int IBlock);
    void deleteBlockMpcToMpcConnectivity();
    void deleteMyCPUTmpBlockCCtsolver();
    void deleteOneBlockMpcToMpcConnectivity(int IBlock);
    void deleteBlockMpcEqNums();
    void deleteOneBlockMpcEqNums(int IBlock, int *gBlockIdArray);
    void deleteBlockCCtsolver();
    void deleteOneBlockCCtsolver(int IBlock, int *gBlockIdArray);
    void deleteMpcVector();
    void initBlockMpcResidual() const;
    void solveBlockCCt(GenDistrVector<Scalar> &v) const;
    void extractBlockMpcResidual(GenDistrVector<Scalar> &v) const;
    void insertBlockMpcResidual(int iSub, GenDistrVector<Scalar> &v);
    void initOneBlockMpcResidual(int IBlock) const;
    void solveOneBlockCCt(int IBlock, GenDistrVector<Scalar> &v) const; // Non const objects pointed to by members
    void extractOneBlockMpcResidual(int IBlock, GenDistrVector<Scalar> &v) const; // Non const objects pointed to by members
    void sendOneBlockCCtsolver(int IBlock);
    void recOneBlockCCtsolver(int IBlock);
    void sendBlockMpcResidualBeforeSolve();
    void recBlockMpcResidualBeforeSolve();
    void sendBlockMpcResidualAfterSolve();
    void recBlockMpcResidualAfterSolve();
    void sendOneBlockMpcResidualBeforeSolve(int IBlock);
    void recOneBlockMpcResidualBeforeSolve(int IBlock);
    void sendOneBlockMpcResidualAfterSolve(int IBlock);
    void recOneBlockMpcResidualAfterSolve(int IBlock);
    void makeBlkToCpuSharedMpcsMap();
};

#endif
