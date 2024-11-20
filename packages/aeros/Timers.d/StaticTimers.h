#ifndef _STATICTIMER_H_
#define _STATICTIMER_H_

#include <cstdio>
#include <Utils.d/resize_array.h>

class Timings;
class MatrixTimers;
class SolverInfo;
class ControlInfo;
class Domain;
class FourierHelmBCs;


void
computeMinAvgMax(double value, long &minValue, long &totValue,
                               long &maxValue);



struct NonlinearNorms {
  double normDv;
  double relativeDv;
  double normRes;
  double relativeRes;
  int stagnated;
  int rebuildTang;
};


static NonlinearNorms emptyNorms;

class StaticTimers {
public:
  int    numSubdomain, numSystems;
  double assemble, solve, factor,
         formRhs,  output, preProcess,
         makeBCs,  makeDOFs, getFetiSolverTime,
         corotatorTime, kelArrayTime, timeTimers,
	 timePresc, timeCheck, timeGeom, timePre,
	 timeBuild, predictorTime, correctorTime,
         buildStiffAndForce, rebuild,
         timeFreqSweep, precond, sfemBuildOps, 
         tdenforceTime, updateSurfsTime, contactSearchTime, contactForcesTime;

  long memoryPreProcess;
  long memoryOutput;
  long memoryPrecond;
  long memoryK;
  long memoryRhs;
  long memoryFreqSweep;
  long memorySfemBuildOps;

  ResizeArray<NonlinearNorms> norms;

  FILE *f;

  StaticTimers() : norms(emptyNorms) { numSubdomain = 1;
                   assemble = 0.0;
                   solve    = 0.0;
                   factor   = 0.0;
                   formRhs  = 0.0;
                   output   = 0.0;
                   preProcess = 0.0;
		   makeBCs    = 0.0;
		   makeDOFs   = 0.0;
                   getFetiSolverTime = 0.0;
                   buildStiffAndForce = 0.0;
	           corotatorTime = 0.0;
		   kelArrayTime = 0.0;
                   rebuild = 0.0;
		   timeTimers = 0.0;
		   timePresc = 0.0;
		   timeCheck = 0.0;
		   timeGeom  = 0.0;
		   timePre   = 0.0;
		   timeBuild = 0.0;
                   predictorTime = 0.0;
                   correctorTime = 0.0;
                   timeFreqSweep = 0.0;
                   precond = 0.0;
                   sfemBuildOps = 0.0;
                   tdenforceTime = 0.0;
                   updateSurfsTime = 0.0;
                   contactSearchTime = 0.0;
                   contactForcesTime = 0.0;
                   f = 0;
                   numSystems = 0;
                   memoryPreProcess = 0;
                   memoryOutput = 0;
                   memoryPrecond = 0;
                   memoryK = 0;
                   memoryRhs = 0;
                   memoryFreqSweep = 0;
                   memorySfemBuildOps = 0;
	         }

  void openTimingFile(ControlInfo &cinfo);

  // Non FETI Solve
  void printStaticTimers(double solveTime, long memoryUsed, 
                         Domain *domain=0, double timeLoop=0.0);

  // FETI Solve
  void printStaticTimers(MatrixTimers matrixTimers, double solveTime,
                         SolverInfo& sinfo, Timings &timers,
                         ControlInfo& cinfo, Domain *domain);
  // FETI-DP Solve
  void printFetiDPtimers(MatrixTimers matrixTimers, double solveTime,
                         SolverInfo& sinfo, Timings &timers,
                         ControlInfo& cinfo, Domain *domain);

  // FETI-H AXI Solver
  void printStaticTimersFetiHAxi(MatrixTimers matrixTimers, double solveTime,
                         SolverInfo& sinfo, Timings &timers,
                         ControlInfo& cinfo, Domain *domain,
                         FourierHelmBCs *glBC, int MPCSize);

  // nonlinear
  void printTimers(Domain *domain, Timings& times, double solveTime);

};

#endif
