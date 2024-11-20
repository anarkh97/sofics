#ifndef _TIMING_H_
#define _TIMING_H_

#include <cstdio>

#if defined(sgi) &&  !defined(_OPENMP)
static const double oneMegaByte = (1024.0*1024.0); 
#else
static const double oneMegaByte = 1024.0;
#endif

void startTimer( double &time);
void stopTimer( double &time);
void startTimerMemory( double &time, long &mem );
void stopTimerMemory( double &time, long &mem );

#include <Utils.d/resize_array.h>
#include <Timers.d/DistTimer.h>

struct NonlinearData {
  int numFetiIter;
  double finalPrimal;
  double finalDual;
  int stagnated;
  int rebuildPrec;
  int rebuildKrylov;
  double cpuTime;
  NonlinearData() { stagnated = 0; }
};

class Timings {
public:
    double assemble,    solve,      factor,   
           project,     project2, sAndJ,      reOrtho,
           precond,     addR,       pfactor,
           reBuildFETI, reBuildGtG, reBuildPrec, reBuildError,
           nlPreCond,   reBuildPCtFPC, readDecomp,   rhsTime,
           outputTime,  coarse1, coarse2, forBack, forBack2,
           fgt, subAlphaG, constructFETI,
           pfactor2, scaling, constructMatrices,
           makeGrbm, locSing, planning, buildCCt, solveCCt;

    // New timers to output details in timing and memory
    DistTimer consMatrix;
    DistTimer assembleMat;
    DistTimer factorMat;
    DistTimer buildRhs;
    DistTimer solveAndJump;
    DistTimer preconditioner;
    DistTimer orthogonalize;
    DistTimer projection;
    DistTimer kMem;
    DistTimer applyFetiPrecond;

    int numIter;
    int numEdges;
    int numThreads;
    int numSubdomains;
    int converged;
    int numRBMs; // number of rigid body modes
    int numCRNs; // number of corner modes
    int numMPCs;

    // variables to store memory used under various categories
    long memoryFETI;
    long memoryProject1;
    long memoryProject2;
    long memoryPrecond;
    long memoryGtG;
    long memoryGtGsky;
    long memoryGtGDelete;
    long memoryPCtFPC;
    long memoryPCtFPCmat;
    long memorySubMatrices;
    long memoryFactor;
    long memoryOSet;
    long memorySAndJ;
    long memoryDV;
    long memoryLoad; // not needed
    long memorySolve;
    long memoryGrbm;
    long memorySing;
    long memoryNlPreCond;
    long memoryPlanning;
    long memoryBuildCCt;
    long memorySolveCCt;
    double finalPrimal;
    double finalDual;

    int numSystems;
    ResizeArray<NonlinearData> iterations;

    // Constructor
    Timings(int _numThreads=1, int _numSubdomains=1); 

    void printFETItimers(FILE *f, long memoryK, 
                                  long memoryPrec,
                                  double initTime, double solveTime,
                                  long initMemory, long solveMemory);
    void setStagnate(int numSystems);

    void startTimerMemory(double &time, long &mem);
    void  stopTimerMemory(double &time, long &mem);

    long globalMemoryMin(long mem);
    long globalMemoryMax(long mem);
    long globalMemorySum(long mem);
    double globalMemorySum(double mem);

};

#endif
