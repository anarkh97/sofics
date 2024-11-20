#include <cstdio>
#include <sys/time.h>

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif
#include <Timers.d/Timing.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>

#include <iostream> 
using namespace std;

//const double megaByte = 1024.0*1024.0;
const double byteToMb = 1.0/oneMegaByte;

static NonlinearData emptyNLI;

Timings::Timings(int nThr, int nsub) : consMatrix(nThr),
                             assembleMat(nThr), factorMat(nThr),
                             buildRhs(nThr), solveAndJump(nThr),
                             preconditioner(nThr), orthogonalize(nThr),
                             projection(nThr), kMem(nThr), applyFetiPrecond(nThr), iterations(emptyNLI)
{
      numThreads    = nThr;
      numSubdomains = nsub;

      numSystems  = 0;
      numIter     = 0;
      finalPrimal = 0.0;
      finalDual   = 0.0;
      precond     = 0.0;
      assemble    = 0.0;
      constructMatrices = 0.0;
      constructFETI = 0.0;
      scaling   = 0.0;
      solve     = 0.0;
      factor    = 0.0;
      project   = 0.0;
      project2  = 0.0;
      reOrtho   = sAndJ = 0.0;
      addR      = 0.0;
      pfactor   = 0.0;
      pfactor2  = 0.0;
      reBuildFETI   = 0.0;
      reBuildGtG    = 0.0;
      reBuildPrec   = 0.0;
      reBuildError  = 0.0;
      nlPreCond     = 0.0;
      reBuildPCtFPC = 0.0;
      readDecomp    = 0.0;
      rhsTime       = 0.0;
      outputTime    = 0.0;
      coarse1    = 0.0;
      coarse2    = 0.0;
      planning   = 0.0;
      converged  = -1;
      numRBMs    = 0;
      numCRNs    = 0;
      memorySolve=0;
      memoryFETI = 0;
      memoryProject1 = 0;
      memoryProject2 = 0;
      memoryPrecond = 0;
      memoryGtG = 0;
      memoryGtGDelete = 0;
      memoryGtGsky = 0;
      memoryPCtFPC = 0;
      memoryPCtFPCmat = 0;
      memorySubMatrices = 0;
      memoryFactor = 0;
      memoryOSet  = 0;
      memorySAndJ = 0;
      memoryDV = 0;
      memoryLoad = 0;
      memoryNlPreCond = 0;
      memoryPlanning = 0;
      forBack  = 0.0;
      forBack2 = 0.0;
      fgt = 0.0;
      subAlphaG = 0.0;
      buildCCt = 0.0;
      solveCCt = 0.0;
      memoryBuildCCt = 0;
      memorySolveCCt = 0;

}


#ifdef DISTRIBUTED

long
Timings::globalMemoryMin(long memory)
{
 double mem1 = (double) memory;
 mem1 = structCom->globalMin(mem1);
 return (long) mem1;
}

long
Timings::globalMemoryMax(long memory)
{
 double mem1 = (double) memory;
 mem1 = structCom->globalMax(mem1);
 return (long) mem1;
}

long
Timings::globalMemorySum(long memory)
{
 double mem1 = (double) memory;
 mem1 = structCom->globalSum(mem1);
 return (long) mem1;
}

double
Timings::globalMemorySum(double memory)
{
 return structCom->globalSum(memory);
}

void
Timings::printFETItimers(FILE *f, long memoryK, long memoryPrec,
                         double initTime, double solveTime,
                         long initMemory, long solveMemory)
{

 TimeData constructMin = consMatrix.getMin();
 TimeData constructMax = consMatrix.getMax();
 TimeData constructAvg = consMatrix.getAvg();
 TimeData constructTot = consMatrix.getTot();

 TimeData assembleMin  = assembleMat.getMin();
 TimeData assembleMax  = assembleMat.getMax();
 TimeData assembleAvg  = assembleMat.getAvg();
 TimeData assembleTot  = assembleMat.getTot();

 TimeData factorMin    = factorMat.getMin();
 TimeData factorMax    = factorMat.getMax();
 TimeData factorAvg    = factorMat.getAvg();
 TimeData factorTot    = factorMat.getTot();

 TimeData sAndJMin     = solveAndJump.getMin();
 TimeData sAndJMax     = solveAndJump.getMax();
 TimeData sAndJAvg     = solveAndJump.getAvg();
 TimeData sAndJTot     = solveAndJump.getTot();

 TimeData precMin      = preconditioner.getMin();
 TimeData precMax      = preconditioner.getMax();
 TimeData precAvg      = preconditioner.getAvg();

 TimeData orthoMin     = orthogonalize.getMin();
 TimeData orthoMax     = orthogonalize.getMax();
 TimeData orthoAvg     = orthogonalize.getAvg();
 TimeData orthoTot     = orthogonalize.getTot();

 TimeData proj1Min     = projection.getMin();
 TimeData proj1Max     = projection.getMax();
 TimeData proj1Avg     = projection.getAvg();
 TimeData proj1Tot     = projection.getTot();

 // Use a global sum to find the total number of subdomains
 numSubdomains = structCom->globalSum(numSubdomains);

 // use system communicator to find the global number of cpus
 int numCPUs = structCom->numCPUs();

 //long memFETIMin = globalMemoryMin(memoryFETI);
 //long memFETIMax = globalMemoryMax(memoryFETI);
 //long memFETITot = globalMemorySum(memoryFETI);

 long memOSetMin = globalMemoryMin(memoryOSet);
 long memOSetMax = globalMemoryMax(memoryOSet);
 long memOSetTot = globalMemorySum(memoryOSet);
 
 //long memDVMin = globalMemoryMin(memoryDV);
 //long memDVMax = globalMemoryMax(memoryDV);
 long memDVTot = globalMemorySum(memoryDV);

 long memGtGMin = globalMemoryMin(memoryGtG - memoryGtGDelete);
 long memGtGMax = globalMemoryMax(memoryGtG - memoryGtGDelete);
 long memGtGTot = globalMemorySum(memoryGtG - memoryGtGDelete);

 // compute min,max,sum memory for subdomain stiffness matrix
 long minMemKL = globalMemoryMin(memoryK);
 long maxMemKL = globalMemoryMax(memoryK);
 long totMemKL = globalMemorySum(memoryK);

 // compute min,max,sum memory for subdomain preconditioner
 long minMemPrecL = globalMemoryMin(memoryPrec);
 long maxMemPrecL = globalMemoryMax(memoryPrec);
 long totMemPrecL = globalMemorySum(memoryPrec);

 // compute min,max,sum memory for FETI98_init()
 long initMemoryMin = globalMemoryMin(initMemory);
 long initMemoryMax = globalMemoryMax(initMemory);
 long initMemoryTot = globalMemorySum(initMemory);

 // compute min,max,sum memory for FETI98_solve()
 //long solveMemoryMin = globalMemoryMin(solveMemory);
 //long solveMemoryMax = globalMemoryMax(solveMemory);
 //long solveMemoryTot = globalMemorySum(solveMemory);

 double coarse1Max  = structCom->globalMax(coarse1);
 double coarse1Tot  = structCom->globalSum(coarse1);
 double coarse1Min  = structCom->globalMin(coarse1);

 double coarse2Min  = 0.0;
 double coarse2Tot  = 0.0;
 double coarse2Max  = 0.0;

 // solveTime comes from timing the feti98_solve() routine
 double solveMax    = structCom->globalMax(solveTime);
 double solveTot    = structCom->globalSum(solveTime);
 double solveMin    = structCom->globalMin(solveTime);

 double initTimeMax = structCom->globalMax(initTime);
 double initTimeTot = structCom->globalSum(initTime);
 double initTimeMin = structCom->globalMin(initTime);

 double constructFETImax = structCom->globalMax(constructFETI);
 double constructFETItot = structCom->globalSum(constructFETI);
 double constructFETImin = structCom->globalMin(constructFETI);

 double forBack1Max = structCom->globalMax(forBack);
 double forBack1Tot = structCom->globalSum(forBack);
 double forBack1Min = structCom->globalMin(forBack);

 // double forBack2Max = 0.0;
 // double forBack2Tot = 0.0;
 // double forBack2Min = 0.0;

 double parfac1Max  = structCom->globalMax(pfactor);
 double parfac1Tot  = structCom->globalSum(pfactor);
 double parfac1Min  = structCom->globalMin(pfactor);

 double parfac2Min  = 0.0;
 double parfac2Tot  = 0.0;
 double parfac2Max  = 0.0;

 // projection time includes forward/backward for GtQG
 double project1Max = proj1Max.time + forBack1Max;
 double project1Avg = proj1Avg.time + forBack1Tot/numCPUs;
 double project1Min = proj1Min.time + forBack1Min;

 // double project2Min = 0.0;
 // double project2Tot = 0.0;
 // double project2Max = 0.0;

 if(f && structCom->myID()==0) { 

 long memCoarseMin = memGtGMin;
 long memCoarseTot = memGtGTot;
 long memCoarseMax = memGtGMax;

 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f," ... Global Information");
 fprintf(f,"\n***********************************************************"
           "********************\n");

 // For Salinas, this is the time for feti98_init() to execute
 double totalTime1      = initTimeMax;
 long totalMemory1 = initMemoryTot - totMemKL         - totMemPrecL 
                                        - factorTot.memory - memCoarseTot;

 fprintf(f,"\n1. Total FETI98 Initialization         time: %15.5f s %12.3f"
           " Mb\n\n",totalTime1/1000.0,totalMemory1*byteToMb);

 double coarseTime = coarse1Max + coarse2;

 long totalMemorySolve = memoryProject1 + memoryProject2
                            + memorySAndJ    + memOSetTot     + memDVTot
                            + totMemPrecL;

 long totalMemory2 = factorTot.memory + memCoarseTot +
                               totalMemorySolve + totMemKL;

 double totalTime2 = solveMax + constructFETImax;

 fprintf(f,"2. Total Solver                        time: %15.5f s %12.3f Mb\n\n",
         totalTime2/1000.0, totalMemory2*byteToMb);
 fprintf(f,"         Factor Subdomain Matrices     time: %15.5f s %12.3f Mb\n\n",
         factorTot.time/1000.0, (totMemKL + factorTot.memory)*byteToMb);
 fprintf(f,"         Total Building  Coarse Pbs.   time: %15.5f s %12.3f Mb\n",
         coarseTime/1000.0, memCoarseTot *byteToMb );
 fprintf(f,"               1st Level Coarse Pb.    time: %15.5f s %12.3f Mb\n",
         coarse1Max/1000.0, memGtGTot*byteToMb);
// fprintf(f,"               2nd Level Coarse Pb.    time: %15.5f s %12.3f Mb\n\n",         coarse2Max/1000.0, memoryPCtFPC*byteToMb);
 fprintf(f,"         Total Paral. Fac. Coarse Pbs. time: %15.5f s %12.3f Mb\n",
         (parfac1Max + parfac2Max)/1000.0, 0.0);
 fprintf(f,"               1st Level Factor        time: %15.5f s %12.3f Mb\n",
         parfac1Max/1000.0, 0.0);
// fprintf(f,"               2nd Level Factor        time: %15.5f s %12.3f Mb\n\n",
//         parfac2Max/1000.0, 0.0);
 fprintf(f,"         Total Solve loop              time: %15.5f s %12.3f Mb\n",
         solveMax/1000.0, totalMemorySolve*byteToMb);
 fprintf(f,"               Project 1st Level       time: %15.5f s %12.3f Mb\n",
         project1Max/1000.0, memoryProject1*byteToMb);
 fprintf(f,"                  Seq. Forw/Back       time: %15.5f s %12.3f Mb\n",
         forBack1Max/1000.0,0.0);
// fprintf(f,"               Project 2nd Level       time: %15.5f s %12.3f Mb\n",
//         project2Max/1000.0,memoryProject2*byteToMb);
// fprintf(f,"                  Paral. Forw/Back     time: %15.5f s %12.3f Mb\n",
//         forBack2Max/1000.0,0.0);
 fprintf(f,"               Precondition            time: %15.5f s %12.3f Mb\n",
         precMax.time/1000.0, totMemPrecL*byteToMb);
 fprintf(f,"               Local Solve             time: %15.5f s %12.3f Mb\n",
         sAndJMax.time/1000.0, sAndJMax.memory*byteToMb);
 fprintf(f,"               Reorthogonalize         time: %15.5f s %12.3f Mb\n",
         orthoMax.time/1000.0, orthoTot.memory*byteToMb);


 // CHECK these numbers later with an overall timer in feti98_init() and feti98_solve()
 long totalSolutionMemory = totalMemory1 + totalMemory2;
 double totalSolutionTime      = totalTime1 + totalTime2;

 fprintf(f,"\nTOTAL SOLUTION (1+2)                   time: %15.5f s %12.3f Mb\n",
         totalSolutionTime/1000.0,totalSolutionMemory*byteToMb);

 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f," ... Detailed CPU Statistics (Seconds) ");
 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f,"\n                                             minimum      average      maximum\n");

 double tot1MinTime = initTimeMin;
 double tot1AvgTime = initTimeTot/numCPUs;
 double tot1MaxTime = initTimeMax;

 fprintf(f,"\n1. Total FETI98 Initialization        : %12.4f %12.4f %12.4f"
           "\n\n",tot1MinTime/1000.0,tot1AvgTime/1000.0,tot1MaxTime/1000.0);

 // timers for 2nd level coarse problem

 double timeCoarseMin = coarse1Min + coarse2Min;
 double timeCoarseTot = coarse1Tot + coarse2Tot;
 double timeCoarseMax = coarse1Max + coarse2Max; 

 double timeParFacMin = parfac1Min + parfac2Min;
 double timeParFacTot = parfac1Tot + parfac2Tot;
 double timeParFacMax = parfac1Max + parfac2Max;

 double tot2MinTime = solveMin + constructFETImin;
 double tot2AvgTime = (solveTot+constructFETItot)/numCPUs; 
 double tot2MaxTime = solveMax + constructFETImax;

 fprintf(f,"2. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
         tot2MinTime/1000.0, tot2AvgTime/1000.0, tot2MaxTime/1000.0);

 fprintf(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
         factorMin.time/1000.0,factorAvg.time/1000.0,factorMax.time/1000.0);

 fprintf(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
         timeCoarseMin/1000.0,timeCoarseTot/(numCPUs*1000.0)
        ,timeCoarseMax/1000.0);
 fprintf(f,"               1st Level Coarse Pb.   : %12.4f %12.4f %12.4f\n",
         coarse1Min/1000.0,coarse1Tot/(numCPUs*1000.0),coarse1Max/1000.0);
 fprintf(f,"               2nd Level Coarse Pb.   : %12.4f %12.4f %12.4f\n\n",
         coarse2Min/1000.0,coarse2Tot/(numCPUs*1000.0),coarse2Max/1000.0);

 fprintf(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
    timeParFacMin/1000.0,timeParFacTot/(numCPUs*1000.0),timeParFacMax/1000.0);
 fprintf(f,"               1st Level Factor       : %12.4f %12.4f %12.4f\n",
         parfac1Min/1000.0,parfac1Tot/(numCPUs*1000.0),parfac1Max/1000.0);
// fprintf(f,"               2nd Level Factor       : %12.4f %12.4f %12.4f\n\n",
//         parfac2Min/1000.0,parfac2Tot/(numCPUs*1000.0),parfac2Max/1000.0);

 fprintf(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
         solveMin/1000.0,solveTot/(numCPUs*1000.0),solveMax/1000.0);
 fprintf(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
         project1Min/1000.0,project1Avg/1000.0,project1Max/1000.0);
 fprintf(f,"                  Seq. Forw/Back      : %12.4f %12.4f %12.4f\n",
         forBack1Min/1000.0,forBack1Tot/(numCPUs*1000.0), forBack1Max/1000.0);
// fprintf(f,"               Project 2nd Level      : %12.4f %12.4f %12.4f\n",
//         project2Min/1000.0,project2Tot/(numCPUs*1000.0),project2Max/1000.0);
// fprintf(f,"                  Paral. Forw/Back    : %12.4f %12.4f %12.4f\n",
//         forBack2Min/1000.0,forBack2Tot/(numCPUs*1000.0), forBack2Max/1000.0);
 fprintf(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
         precMin.time/1000.0, precAvg.time/1000.0, precMax.time/1000.0);
 fprintf(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
         sAndJMin.time/1000.0,sAndJAvg.time/1000.0,sAndJMax.time/1000.0);
 fprintf(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n",
         orthoMin.time/1000.0,orthoAvg.time/1000.0,orthoMax.time/1000.0);

 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f," ... Detailed Memory Statistics (Megabytes)");
 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f,"\n                                             minimum      average      maximum\n");

 // 
 // Note, memory for the subdomain matrices, coarse problems, rigid body modes,
 // preconditioner occurs in FETI98_init() but is counted in under 
 // Total Solver categories
 //

 long tot1Min = initMemoryMin - minMemKL - minMemPrecL 
                                   - factorMin.memory - memCoarseMin;
 long tot1Avg = (initMemoryTot - totMemKL - totMemPrecL - memCoarseTot)/numCPUs - factorAvg.memory;
 long tot1Max = initMemoryMax - maxMemKL - maxMemPrecL
                                   - factorMax.memory - memCoarseMax;

 fprintf(f,"\n1. Total FETI98 Initialization        : %12.4f %12.4f %12.4f"
           "\n\n",   tot1Min*byteToMb,tot1Avg*byteToMb,tot1Max*byteToMb);


 long tot2Min = factorMin.memory + memCoarseMin + memOSetMin 
                   + minMemPrecL + minMemKL;
 long tot2Tot = factorAvg.memory + memCoarseTot/numCPUs 
                   + memOSetTot/numCPUs + totMemPrecL/numCPUs 
                   + totMemKL/numCPUs;
 long tot2Max = factorMax.memory + memCoarseMax + memOSetMax
                   + maxMemPrecL + maxMemKL;

 fprintf(f,"2. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
         tot2Min*byteToMb, tot2Tot*byteToMb, tot2Max*byteToMb);

 fprintf(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
         (factorMin.memory+minMemKL)*byteToMb, 
         (factorAvg.memory+(totMemKL/numCPUs))*byteToMb,
         (factorMax.memory+maxMemKL)*byteToMb);
 
 fprintf(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
           memCoarseMin*byteToMb,
           memCoarseTot*byteToMb/numCPUs,
           memCoarseMax*byteToMb);

 fprintf(f,"               1st Level Coarse Pb.   : %12.4f %12.4f %12.4f\n",
         memGtGMin*byteToMb, memGtGTot*byteToMb/numCPUs,memGtGMax*byteToMb);

// fprintf(f,"               2nd Level Coarse Pb.   : %12.4f %12.4f %12.4f\n\n",
//         0.0,0.0,0.0);
 fprintf(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
         0.0,0.0,0.0);
 fprintf(f,"               1st Level Factor       : %12.4f %12.4f %12.4f\n",
         0.0,0.0,0.0);
// fprintf(f,"               2nd Level Factor       : %12.4f %12.4f %12.4f"
//           "\n\n",0.0, 0.0,0.0);

 long totloopMin = minMemPrecL + orthoMin.memory;
 long totloopAvg = totMemPrecL/numCPUs + orthoAvg.memory;
 long totloopMax = maxMemPrecL + orthoMax.memory;

 fprintf(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
         totloopMin*byteToMb,totloopAvg*byteToMb,totloopMax*byteToMb);
 fprintf(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
         0.0,0.0,0.0);
// fprintf(f,"               Project 2nd Level      : %12.4f %12.4f %12.4f\n",
//         0.0,0.0,0.0);
 fprintf(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
         minMemPrecL*byteToMb, 
         (totMemPrecL/numCPUs)*byteToMb,
         maxMemPrecL*byteToMb);
 fprintf(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
         0.0,0.0,0.0);
 fprintf(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n",
         orthoMin.memory*byteToMb, orthoAvg.memory*byteToMb,
         orthoMax.memory*byteToMb);

 long totSimMin = tot1Min + tot2Min;
 long totSimAvg = tot1Avg + tot2Tot;
 long totSimMax = tot1Max + tot2Max;
 
 fprintf(f,"\nTOTAL SOLUTION (1+2)                  : %12.4f %12.4f %12.4f\n",
           totSimMin*byteToMb, totSimAvg*byteToMb, totSimMax*byteToMb);


 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f," ... FETI Monitoring for %d CPUs and %d Subdomains \n",
           numCPUs,numSubdomains);
 fprintf(f,"***********************************************************"
           "********************\n\n");

   // Check whether we converged or not
   if(converged)
   fprintf(f,"1. Number of Iterations for Convergence    = %15d\n\n",
           numIter);
   else {
   fprintf(f,"1. Stagnation Occured After a # of iter    = %15d\n\n",
           numIter);
   }

   fprintf(f,"2. Relative Primal Error Reached           = %15.3e\n\n",
           iterations[0].finalPrimal);

   fprintf(f,"3. Relative Dual Error Reached             = %15.3e\n\n",
           iterations[0].finalDual);

   fprintf(f,"4. Size of 1st Level CP Before Delete      = %15d %15.3f Mb\n\n",
           numRBMs,memoryGtGsky*byteToMb);

   fprintf(f,"5. Size of 1st Level CP After Delete       = %31.3f Mb\n\n",
           (memoryGtGsky-memoryGtGDelete)*byteToMb);
   
   fprintf(f,"***********************************************************"
             "********************\n");
 }

}
#endif

void
Timings::setStagnate(int _numSystems)
{
    converged = 0;
    iterations[_numSystems].stagnated = 1;
}

void
Timings::startTimerMemory(double &time, long &mem)
{
 time -= getTime();
 mem  -= memoryUsed();
} 

void
Timings::stopTimerMemory(double &time, long &mem)
{
 time += getTime();
 mem  += memoryUsed();
}


void
startTimer(double &time)
{
 time -= getTime();
}

void
stopTimer(double &time)
{
 time += getTime();
}

void
startTimerMemory(double &time, long &mem)
{
  time -= getTime();
  mem  -= memoryUsed();
}

void
stopTimerMemory(double &time, long &mem)
{
 time += getTime();
 mem  += memoryUsed();
}


