#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <ctime>
#include <algorithm>

#include <Timers.d/Timing.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/MatrixTimers.h>
#include <Driver.d/Domain.h>
#include <Threads.d/Paral.h>
#include <Utils.d/Memory.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>

extern long totMemSparse;
extern long totMemSky;
extern long totMemSpooles;
extern long totMemMumps;

const char* yesno[] = {" no","yes"};
// problem Type solved
const char* problemType[] = {
"Linear Static",
"Linear Dynamic",
"Modal",
"Non-Linear Static",
"Non-Linear Dynamic",
"Non-Linear Static (ArcLength)",
"Condition Number Estimate",
"Thermal Dynamic",
"",
"Axisymmetric Helmholtz",
"Material Non-Linear Static",
"Material Non-Linear Dynamic",
"Helmholtz (Acoustic Frequency Response)",
"Helmholtz Frequency Sweep",
"Helmholtz Direction Sweep",
"",
"",
"",
""
};

// selected solver messages
const char* message[] = {
"1. Skyline    \n",
"1. Sparse    \n",
"1. BlockSky \n",
"1. Simplicial LLT \n",
"1. Simplicial LDLT \n",
"1. Cholmod Sparse \n",
"1. Umfpack Sparse\n",
"1. SuperLU Sparse \n",
"1. Spooles Sparse\n",
"1. Mumps Sparse\n",
"1. Diagonal\n",
"1. CG\n",
"1. CG\n         Diagonal Preconditioner\n",
"1. CG\n         FETI Preconditioner\n",
"1. Goldfarb-Idnani\n",
"1. Sparse LU\n",
"1. SuiteSparseQR\n",
"1. Sparse QR\n"
};

const char* renumMessage[] = {
"None",
"Sloan",
"RCM"
};

void
StaticTimers::openTimingFile(ControlInfo &cinfo)
{
  // Open timing file
  if(f == 0) {
    int len = strlen(cinfo.checkfile);
    char *fileName = (char *) dbg_alloca(sizeof(char)*(len+8));
    strcpy(fileName, cinfo.checkfile);
    strcat(fileName,".timing");
    if((f=fopen(fileName,"w"))==(FILE *) 0 )
      filePrint(stderr," *** ERROR: Cannot open %s ...\n", cinfo.checkfile );
    fflush(f);
  }
}

const double byteToMb    = (1.0 / oneMegaByte);

// This function prints timers for non-FETI solvers: i.e. skyline, pcg, etc.

void
StaticTimers::printStaticTimers(double solveTime, long memUsed,
                                Domain *domain, double timeLoop)
{
#ifndef SALINAS
	MatrixTimers &times  = domain->getTimers();
	SolverInfo &sInfo    = domain->solInfo();
	ControlInfo *cinfo   = geoSource->getCheckFileInfo();

	openTimingFile(cinfo[0]);

	int mesNum = sInfo.solvercntl->subtype;
	if(sInfo.solvercntl->type == SolverSelection::Iterative && sInfo.solvercntl->precond == 0) mesNum = 11;
	if(sInfo.solvercntl->type == SolverSelection::Iterative && sInfo.solvercntl->precond == 1) mesNum = 12;
	if( (sInfo.solvercntl->type == SolverSelection::Feti ||sInfo.solvercntl->type == SolverSelection::FetiLib)
	    && sInfo.inpc) mesNum = 13;

	int numnod   = domain->numNodes();
	int numele   = domain->numElements();
	int numdof   = domain->numdof();
	int numcon   = domain->nDirichlet() + domain->nCDirichlet();
	int numUncon = domain->numUncon();

	filePrint(f,"\n***********************************************************"
	            "********************\n");
	if(geoSource->isShifted()) {
		if(domain->solInfo().doFreqSweep)
			filePrint(f," ... Frequency Sweep Problem Information ... \n");
		else
			filePrint(f," ... Frequency Response Problem Information ... \n");
	}
	else
		filePrint(f," ... %s Problem Information ... \n",problemType[sInfo.probType]);
	filePrint(f,"***********************************************************"
	            "********************\n\n");
	filePrint(f,"1. Number of Nodes                         = %14d\n\n",numnod);
	filePrint(f,"2. Number of Elements                      = %14d\n\n",numele);
	filePrint(f,"3. Number of Degrees of Freedom            = %14d\n",numdof);
	filePrint(f,"         Number of Constrained Dofs        = %14d\n",numcon);
	filePrint(f,"         Number of Unconstrained Dofs      = %14d\n\n",numUncon);
	filePrint(f,"4. Number of Applied Loads                 = %14d\n\n",
	          domain->nNeumann());
	filePrint(f,"5. Number of Output Files                  = %14d\n\n",
	          geoSource->getNumOutInfo());
	filePrint(f,"6. Renumbering                             = %14s\n\n",
	          renumMessage[sInfo.renum]);

	if(domain->solInfo().doFreqSweep) {
		filePrint(f,"7. Number of Frequencies                   = %14d\n\n", domain->numFrequencies);
		filePrint(f,"8. Number of RHS solves                    = %14d\n\n", domain->solInfo().getSweepParams()->nFreqSweepRHS);
	}

	filePrint(f,"***********************************************************"
	            "********************\n");
	filePrint(f," ... Solver Information ... \n");
	filePrint(f,"***********************************************************"
	            "********************\n\n");

	if(mesNum > 10) {
		filePrint(f,"%s         Tolerance                         = %14.2e\n         Maximum Number of Iterations      = %14d\n\n",
		          message[mesNum],sInfo.solvercntl->tol,sInfo.solvercntl->maxit);
	} else
		filePrint(f,"%s\n",message[mesNum]);

	filePrint(f,"***********************************************************"
	            "********************\n");
	filePrint(f," ... Timing Statistics for %d Thread%c ...\n",threadManager->numThr(),(threadManager->numThr()==1)?' ':'s');
	filePrint(f,"***********************************************************"
	            "********************\n\n");

	double totalRead = times.readTime + times.readDecomp;
	long totMemRead = times.memoryParse + times.memorySetUp;

	filePrint(f,"1. Total Read Input Files              time: %14.5f s %14.3f Mb\n",
	          totalRead/1000.0, totMemRead*byteToMb);
	filePrint(f,"         Read Mesh                     time: %14.5f s\n\n",
	          times.readTime/1000.0);

	double totalPreProcess = preProcess + times.setUpDataTime;
	long totMemPreProcess = memoryPreProcess + times.memorySetUp;

	filePrint(f,"2. Total Preprocessing                 time: %14.5f s %14.3f Mb\n",
	          totalPreProcess/1000.0, totMemPreProcess*byteToMb);
	filePrint(f,"         Process Input Data            time: %14.5f s\n",
	          times.setUpDataTime/1000.0);
	filePrint(f,"         Connectivity                  time: %14.5f s\n",
	          times.makeConnectivity/1000.0);
	filePrint(f,"         Renumbering                   time: %14.5f s\n",
	          times.renumbering/1000.0);
	if(!(sInfo.solvercntl->type == SolverSelection::Feti && sInfo.inpc)) {
		filePrint(f,"         Create DOFs                   time: %14.5f s\n",
		          times.createDofs/1000.0);
		filePrint(f,"         Make Constrained DOFs         time: %14.5f s\n",
		          makeDOFs/1000.0);
		filePrint(f,"         Make Boundary Conditions      time: %14.5f s\n\n",
		          makeBCs/1000.0);}

	double totalMatrix = (sInfo.solvercntl->type == SolverSelection::Feti && sInfo.inpc) ? sfemBuildOps
	                                                                 : times.constructTime+times.assemble+kelArrayTime+corotatorTime+times.formTime;
	long totMemMatrix  = (sInfo.solvercntl->type == SolverSelection::Feti && sInfo.inpc) ? memorySfemBuildOps
	                                                                 : times.memoryForm;

	filePrint(f,"3. Total Matrix Processing             time: %14.5f s %14.3f Mb\n",
	          totalMatrix/1000.0, totMemMatrix*byteToMb);
	if(!(sInfo.solvercntl->type == SolverSelection::Feti && sInfo.inpc)) {
		filePrint(f,"         Construct Sparse Matrices     time: %14.5f s\n",
		          times.constructTime/1000.0);
		filePrint(f,"         Form Element Matrices         time: %14.5f s\n",
		          (times.formTime+kelArrayTime+corotatorTime)/1000.0);
		filePrint(f,"         Assemble Element Matrices     time: %14.5f s\n\n",
		          times.assemble/1000.0);}

	double totalRhs = formRhs - times.receiveFluidTime + times.formRhs;
	long totMemRhs = memoryRhs;

	filePrint(f,"4. Total RHS Processing                time: %14.5f s %14.3f Mb\n\n",
	          totalRhs/1000.0, totMemRhs*byteToMb);

	double totalSolver = times.factor + solveTime + times.updateState + timeFreqSweep + tdenforceTime;
	long totMemSolver  = times.memorySolve + totMemSpooles + totMemMumps;

	if (mesNum == 2)
		filePrint(f,"5. Total Solver                        time: %14.5f s       N/A\n", totalSolver/1000.0);
	else
		filePrint(f,"5. Total Solver                        time: %14.5f s %14.3f Mb\n", totalSolver/1000.0, totMemSolver*byteToMb);

	if(sInfo.solvercntl->type == SolverSelection::Feti && sInfo.inpc) {
		filePrint(f,"         Preconditioning               time: %14.5f s\n",
		          precond/1000.0);
	}
	else if(sInfo.newmarkBeta != 0.0) { // not relevant for explicit dynamics
		filePrint(f,"         Factor Matrix                 time: %14.5f s\n",
		          (times.factor)/1000.0);

		filePrint(f,"         Solve (Forward/Backward)      time: %14.5f s\n",
		          solveTime/1000.0);
	}

	if(domain->solInfo().isDynam() || domain->solInfo().isNonLin())
		filePrint(f,"         Update States                 time: %14.5f s\n",
		          times.updateState/1000.0);

	if(domain->solInfo().doFreqSweep)
		filePrint(f,"         Freq Sweep Series Expansion   time: %14.5f s\n",
		          timeFreqSweep/1000.0);

	if(domain->tdenforceFlag()){ // ACME search and forces for explicit dynamics
		filePrint(f,"         TD Enforcement                time: %14.5f s\n",
		          tdenforceTime/1000.0);
		filePrint(f,"              detection                time: %14.5f s\n",
		          contactSearchTime/1000.0);
		filePrint(f,"              enforcement              time: %14.5f s\n",
		          contactForcesTime/1000.0);
		filePrint(f,"              surface update           time: %14.5f s\n",
		          updateSurfsTime/1000.0);
	}

	double totalOutput = std::max(output - times.sendFluidTime, 0.);
	long totMemOutput = memoryOutput;

	filePrint(f,"\n6. Write Output Files                  time: %14.5f s %12.3f Mb\n",
	          totalOutput/1000.0, memoryOutput*byteToMb);


	double totalFluidComm = times.receiveFluidTime+times.sendFluidTime;
	if(domain->solInfo().aeroFlag >= 0) {
		filePrint(f,"\n7. Fluid Communication                 time: %14.5f s\n",
		          totalFluidComm/1000.0);
		filePrint(f,"   Receive From Fluid                  time: %14.5f s\n",
		          times.receiveFluidTime/1000.0);
		filePrint(f,"   Send To Fluid                       time: %14.5f s\n\n",
		          times.sendFluidTime/1000.0);
	}

	double total =  totalRead + totalPreProcess + totalMatrix + totalRhs + totalSolver + totalOutput;
	long totalMemoryUsed = totMemRead + totMemPreProcess + totMemMatrix + totMemRhs + totMemSolver + totMemOutput;

	if(domain->solInfo().aeroFlag >= 0) {
		double totalAndFCom = total + totalFluidComm;
		filePrint(f,"\nTOTAL SIMULATION 1 (1+2+3+4+5+6)       time: %14.5f s %14.3f Mb\n",total/1000.0,totalMemoryUsed*byteToMb);
		filePrint(f,"\nTOTAL SIMULATION 2 (1+2+3+4+5+6+7)     time: %14.5f s\n",totalAndFCom/1000.0);
	}
	else {
		filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)         time: %14.5f s %14.3f Mb\n",total/1000.0,totalMemoryUsed*byteToMb);
	}

	filePrint(f,"\n***********************************************************"
	            "********************\n");
/*
 filePrint(f," ... Solution Information ...\n");
 filePrint(f,"***********************************************************"
           "********************\n\n");

 long totMemUsed = memoryUsed();

 // Check if we are using SGI sparse solver
 double coef = 1.0;
 if(sInfo.solvercntl->subtype == 4) {
   totMemUsed += 8*memUsed;
   coef = 12.0;
 } else {
   totMemUsed += (totMemSpooles+totMemMumps);
   memUsed = times.memorySolve;
 }

 if(sInfo.solvercntl->subtype == 3 || sInfo.solvercntl->subtype == 7) coef = 1.0;

 filePrint(f,"1. Total Amount of Requested Memory        = %14.3f Mb\n\n",
           totMemUsed*byteToMb);

 filePrint(f,"2. Solver Amount of Requested Memory       = %14.3f Mb\n\n",
           coef*memUsed*byteToMb);

 // iterative information, needs to be put in here!
 if(sInfo.solvercntl->type == 1) {
   int numIterations = 0;
   double finalNorm  = 0.0;
   filePrint(f,"3. Number of Iterations                    = %14d\n\n",
           numIterations);
   filePrint(f,"4. Relative Error Reached                  = %14.5e\n\n",
           finalNorm);
 }

 filePrint(f,"***********************************************************"
           "********************\n\n");

 if(numSystems > 1) {
   filePrint(f,"\nIter\tDv\t        Relative Dv\tResidual\tRelative Res\n");
   int i;
   for(i=0; i<numSystems; ++i) {
     if(norms[i].relativeDv == 1.0)
       filePrint(f,
     "--------------------------------------------------------------------\n");
     filePrint(f,"%3d%16e\t%e\t%e\t%e\n",i+1,norms[i].normDv,norms[i].relativeDv,
             norms[i].normRes,norms[i].relativeRes);
   }
   filePrint(f,"------------------------------------------------"
             "--------------------\n");
 }
*/

	// Print information for distributed
	if(f && mesNum == 13) {
		filePrint(f," ... Detailed CPU Statistics (Seconds) ");
		filePrint(f,"\n***********************************************************"
		            "********************\n");
		filePrint(f,"\n                                             minimum        average        maximum\n");

#ifdef DISTRIBUTED
		int numCPUs = structCom->numCPUs();

		double tot5MinTime   = structCom->globalMin(totalSolver);
		double tot5MaxTime   = structCom->globalMax(totalSolver);
		double tot5AvgTime   = structCom->globalSum(totalSolver);
		tot5AvgTime /= numCPUs;

		double tot5PrecMinTime   = structCom->globalMin(precond);
		double tot5PrecMaxTime   = structCom->globalMax(precond);
		double tot5PrecAvgTime   = structCom->globalSum(precond);
		tot5PrecAvgTime /= numCPUs;

//  double tot5MinMemory   = structCom->globalMin(times.memorySolve)*byteToMb;

		filePrint(f,"\n5. Total Solver Time               :    %12.4f s %12.4f s %12.4f s\n\n",
		          tot5MinTime/1000.0, tot5AvgTime/1000.0, tot5MaxTime/1000.0);
		filePrint(f,"         Preconditioning  Time     :    %12.4f s %12.4f s %12.4f\n\n",
		          tot5PrecMinTime/1000.0, tot5PrecAvgTime/1000.0, tot5PrecMaxTime/1000.0);


// filePrint(f,"5. Number of CPUs                     : %d\n\n",
//         numCPUs);
#endif
	}


	if(f) { fclose(f); f=0; }
#endif
}

// Message for printing selected scaling method
const char* scalingMessage[] = {
"\n",
"Stiffness Scaling\n",
"Topological Scaling\n",
}; 

// Message for rigid body mode method
const char* rbmMessage[] = {
"Algebraic RBM Method",
"Geometric RBM Method"
};

// Message for printing selected preconditioner method
const char* precMessage[] = {
"No Preconditioner Selected\n",
"Lumped Preconditioner\n",
"Dirichlet Preconditioner\n"};

// Message for printing selected projector method
const char* projectMessage[] = {
"Basic Projector\n",
"Q Projector\n",
"KC Projector\n",
"Q = Diag(Kbb) Projector\n"};

const char* subSolverMessage[] = {
"Subdomain Solver Selected         =        Skyline\n",
"Subdomain Solver Selected         =         Sparse\n",
"Subdomain Solver Selected         =      SgiSparse\n",
"Subdomain Solver Selected         =     SgiSkyline\n",
"Subdomain Solver Selected         =            pcg\n",
"Subdomain Solver Selected         =        frontal\n",
};

const char* precSolverMessage[] = {
"Preconditioner Solver Selected    =        Skyline\n",
"Preconditioner Solver Selected    =         Sparse\n",
"Preconditioner Solver Selected    =      SgiSparse\n",
"Preconditioner Solver Selected    =     SgiSkyline\n",
"Preconditioner Solver Selected    =            pcg\n",
"Preconditioner Solver Selected    =        frontal\n"
};

const char* gtgType[] = {
"GtG Solver Selected               =        ",
"GtQG Solver Selected              =        ",
"GtQG Solver Selected              =        " 
};

const char* gtgSolverMessage[] = {
"Skyline\n",
" Sparse\n",
"Skyline\n", // Not implemented
"Skyline\n", // Not implemented
"    pcg\n",
"Skyline\n", // Not implemented
" BlkSky\n"  // Not implemented
};

const char* solverMessage[] = {
"1. FETI\n         One Level\n",
"1. FETI\n         Two Level\n",
"1. FETI\n         Two Level\n",
"1. FETI-DP\n",
"1. FETI-DPC\n",
"1. FETI-DPH\n"
};


// This function prints timers for FETI Solver 

void
StaticTimers::printStaticTimers(MatrixTimers matrixTimer, double solveTime,
                      SolverInfo& sInfo, Timings& timers, ControlInfo& cinfo,
                      Domain* domain)
{
 double solutionTime = timers.solve + getFetiSolverTime;

 double coarse1Max  = timers.coarse1;
 double coarse1Tot  = timers.coarse1;
 double coarse1Min  = timers.coarse1;

 double parfac1Max  = timers.pfactor;
 double parfac1Tot  = timers.pfactor;
 double parfac1Min  = timers.pfactor;
  
 double forBack1Max = timers.forBack;
 double forBack1Tot = timers.forBack;
 double forBack1Min = timers.forBack;

 double forBack2Max = timers.forBack2;
 double forBack2Tot = timers.forBack2;
 double forBack2Min = timers.forBack2;

 double precondMax   = timers.precond;
 double sAndJMaximum = timers.sAndJ;

 long totMemPreProcess = memoryPreProcess;
 long memorySetUp = matrixTimer.memorySetUp;
 long memorySubMatrices = timers.memorySubMatrices;
 long memorySubdomain = matrixTimer.memorySubdomain;
 long memoryElemToNode = matrixTimer.memoryElemToNode;
 long memorySubToNode = matrixTimer.memorySubToNode;
 long memoryNodeToSub = matrixTimer.memoryNodeToSub;
 long memorySubToElem = matrixTimer.memorySubToElem;
 long memoryCPUMAP = matrixTimer.memoryCPUMAP;

 long memoryDistBC = matrixTimer.memoryDistBC;
 long memoryConnect = matrixTimer.memoryConnect;
 long memoryInterface = matrixTimer.memoryInterface;
 long memoryInternal = matrixTimer.memoryInternal;

#ifdef DISTRIBUTED
 coarse1Max   = structCom->globalMax(timers.coarse1);
 coarse1Tot   = structCom->globalSum(timers.coarse1);
 coarse1Min   = structCom->globalMin(timers.coarse1);
 parfac1Max   = structCom->globalMax(timers.pfactor);
 parfac1Tot   = structCom->globalSum(timers.pfactor);
 parfac1Min   = structCom->globalMin(timers.pfactor);
 forBack1Max  = structCom->globalMax(timers.forBack);
 forBack1Tot  = structCom->globalSum(timers.forBack);
 forBack1Min  = structCom->globalMin(timers.forBack);
 precondMax   = structCom->globalMax(timers.precond);
 sAndJMaximum = structCom->globalMax(sAndJMaximum);
 solutionTime = structCom->globalMax(solutionTime);

 memorySetUp = timers.globalMemorySum(matrixTimer.memorySetUp);
 memorySubMatrices = timers.globalMemorySum(timers.memorySubMatrices);
 totMemPreProcess = timers.globalMemorySum(memoryPreProcess);
 memorySubdomain = timers.globalMemorySum(matrixTimer.memorySubdomain);
 memoryElemToNode = timers.globalMemorySum(matrixTimer.memoryElemToNode);
 memorySubToNode = timers.globalMemorySum(matrixTimer.memorySubToNode);
 memoryNodeToSub = timers.globalMemorySum(matrixTimer.memoryNodeToSub);
 memorySubToElem = timers.globalMemorySum(matrixTimer.memorySubToElem);
 memoryCPUMAP = timers.globalMemorySum(matrixTimer.memoryCPUMAP);
 memoryDistBC = timers.globalMemorySum(matrixTimer.memoryDistBC);
 memoryConnect = timers.globalMemorySum(matrixTimer.memoryConnect);
 memoryInterface = timers.globalMemorySum(matrixTimer.memoryInterface);
 memoryInternal = timers.globalMemorySum(matrixTimer.memoryInternal);
#endif
 totMemPreProcess += memorySetUp;

 TimeData constructMin = timers.consMatrix.getMin();
 TimeData constructMax = timers.consMatrix.getMax();
 TimeData constructAvg = timers.consMatrix.getAvg();
 TimeData constructTot = timers.consMatrix.getTot();

 TimeData assembleMin  = timers.assembleMat.getMin();
 TimeData assembleMax  = timers.assembleMat.getMax();
 TimeData assembleAvg  = timers.assembleMat.getAvg();
 TimeData assembleTot  = timers.assembleMat.getTot();

 TimeData factorMin    = timers.factorMat.getMin();
 TimeData factorMax    = timers.factorMat.getMax();
 TimeData factorAvg    = timers.factorMat.getAvg();
 TimeData factorTot    = timers.factorMat.getTot();

 TimeData buildRhsMin  = timers.buildRhs.getMin();
 TimeData buildRhsMax  = timers.buildRhs.getMax();
 TimeData buildRhsAvg  = timers.buildRhs.getAvg();
 TimeData buildRhsTot  = timers.buildRhs.getTot();

 TimeData sAndJMin     = timers.solveAndJump.getMin();
 TimeData sAndJMax     = timers.solveAndJump.getMax();
 TimeData sAndJAvg     = timers.solveAndJump.getAvg();
 TimeData sAndJTot     = timers.solveAndJump.getTot();

 TimeData precMin      = timers.preconditioner.getMin();
 TimeData precMax      = timers.preconditioner.getMax();
 TimeData precAvg      = timers.preconditioner.getAvg();
 TimeData precTot      = timers.preconditioner.getTot();

 TimeData orthoMin     = timers.orthogonalize.getMin();
 TimeData orthoMax     = timers.orthogonalize.getMax();
 TimeData orthoAvg     = timers.orthogonalize.getAvg();
 TimeData orthoTot     = timers.orthogonalize.getTot();

 TimeData project1Min  = timers.projection.getMin();
 TimeData project1Max  = timers.projection.getMax();
 TimeData project1Avg  = timers.projection.getAvg();
 TimeData project1Tot  = timers.projection.getTot();

 TimeData applyFetiPrecondMin  = timers.applyFetiPrecond.getMin();
 TimeData applyFetiPrecondMax  = timers.applyFetiPrecond.getMax();
 TimeData applyFetiPrecondAvg  = timers.applyFetiPrecond.getAvg();
 TimeData applyFetiPrecondTot  = timers.applyFetiPrecond.getTot();


 int numCPUs       = threadManager->numThr();
 int numSubdomains = timers.numSubdomains;
#ifdef DISTRIBUTED
 numCPUs = structCom->numCPUs();
 numSubdomains = structCom->globalSum(numSubdomains);
 if(structCom->myID() == 0)
#endif
 openTimingFile(cinfo);

 long totMemSubMatrices  = memorySubMatrices
                               - 8*(memoryPrecond) - 8*memoryK ;

 long localMemRead = matrixTimer.memoryParse + matrixTimer.memorySetUp;
 long totalMemRead = localMemRead;

 // Timer sub totals

 double subTotal[6];

 subTotal[0] = (matrixTimer.readTime + matrixTimer.readDecomp);

 subTotal[1] = preProcess + matrixTimer.setUpDataTime 
                          + corotatorTime + kelArrayTime + timeGeom
			  - matrixTimer.readDecomp;

 subTotal[2] = assembleTot.time + constructTot.time
             + matrixTimer.constructTime + matrixTimer.assemble + kelArrayTime + corotatorTime + matrixTimer.formTime;

 subTotal[3] = buildRhsTot.time + matrixTimer.formRhs;

 subTotal[4] = solutionTime - (assembleTot.time + constructTot.time) + tdenforceTime;

 subTotal[5] = output;

 double totPreProMax = subTotal[1];
#ifdef DISTRIBUTED
 totalMemRead = timers.globalMemorySum(localMemRead);
 subTotal[0] = structCom->globalMax(subTotal[0]);
 totPreProMax = structCom->globalMax(subTotal[1]);
 subTotal[2] = structCom->globalMax(subTotal[2]);
 subTotal[3] = structCom->globalMax(subTotal[3]);
 subTotal[4] = structCom->globalMax(subTotal[4]);
 subTotal[5] = structCom->globalMax(subTotal[5]);
#endif


 if(f != 0) {
 
 int numDof = domain->numDofs();
 int numCon = domain->nDirichlet() + domain->nCDirichlet();
 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f," ... %s Problem Information ... \n",problemType[sInfo.probType]);
 filePrint(f,"***********************************************************"
           "********************\n\n");
 
 filePrint(f,"1. Number of Nodes                         = %14d\n\n",
           domain->numNode());
 filePrint(f,"2. Number of Elements                      = %14d\n\n",
           domain->numElements());
 filePrint(f,"3. Number of Degrees of Freedom            = %14d\n",numDof);
 filePrint(f,"         Number of Constrained Dofs        = %14d\n",numCon);
 filePrint(f,"         Number of Unconstrained Dofs      = %14d\n\n",numDof-numCon);
 filePrint(f,"4. Number of Applied Loads                 = %14d\n\n",
           domain->nNeumann());
 filePrint(f,"5. Number of Output Files                  = %14d\n\n",
           geoSource->getNumOutInfo());
 filePrint(f,"6. Renumbering                             = %14s\n\n",
           renumMessage[sInfo.renum]);

 filePrint(f,"***********************************************************"
           "********************\n");
 filePrint(f," ... Solver Information ... \n");
 filePrint(f,"***********************************************************"
           "********************\n\n");

 if(domain->solInfo().solvercntl->type == SolverSelection::Direct) { // MUMPS
   filePrint(f,"1. Mumps Sparse\n");
 }
 else if(domain->solInfo().solvercntl->type == SolverSelection::BlockDiag) {
   filePrint(f,"1. Diagonal\n");
 }
 else { // FETI
   if(sInfo.getFetiInfo().version == FetiInfo::feti1)
     filePrint(f,"%s",solverMessage[0]);
   else 
     filePrint(f,"%s",solverMessage[sInfo.getFetiInfo().feti2version+1]);
  
   filePrint(f,"         %s", precMessage[sInfo.getFetiInfo().precno]);
   filePrint(f,"         %s", scalingMessage[sInfo.getFetiInfo().scaling]);
   filePrint(f,"         %s", projectMessage[sInfo.getFetiInfo().nonLocalQ]);
   filePrint(f,"         %s", subSolverMessage[sInfo.getFetiInfo().local_cntl->subtype]);
   filePrint(f,"         %s", precSolverMessage[sInfo.getFetiInfo().kii_cntl->subtype]);
   filePrint(f,"         %s%s", gtgType[sInfo.getFetiInfo().nonLocalQ],
                                gtgSolverMessage[sInfo.getFetiInfo().coarse_cntl->subtype]);

   if(sInfo.rbmflg == 0)
     filePrint(f,"         %s %29e\n", rbmMessage[sInfo.rbmflg],sInfo.solvercntl->trbm);
   else
     filePrint(f,"         %s%17e %e\n", rbmMessage[sInfo.rbmflg],
                                         sInfo.tolsvd,sInfo.solvercntl->trbm);
  
   filePrint(f,"         Maximum Number of Iterations      = %14d\n",
             sInfo.getFetiInfo().maxiter());
   filePrint(f,"         Maximum Size of Reortho. Vectors  = %14d\n",
             sInfo.getFetiInfo().maxorth());
   filePrint(f,"         Tolerance for Convergence         = %14.3e\n",
             sInfo.getFetiInfo().tolerance());
 }

 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f," ... Timing Statistics for %d Threads and %d Subdomains ...\n",
                numCPUs,numSubdomains);
 filePrint(f,"***********************************************************"
           "********************\n\n");

 // Begin timer output

 filePrint(f,"1. Total Read Input Files              time: %14.5f s %14.3f Mb\n",
           subTotal[0]/1000.0, totalMemRead*byteToMb);
 filePrint(f,"         Read Mesh                     time: %14.5f s\n",
           matrixTimer.readTime/1000.0);
 filePrint(f,"         Read Mesh Partition           time: %14.5f s\n",
           matrixTimer.readDecomp/1000.0);          

 filePrint(f,"\n2. Total Preprocessing                 time: %14.5f s %14.3f Mb\n",
           totPreProMax/1000.0, totMemPreProcess*byteToMb);
 filePrint(f,  "         Process Input Data            time: %14.5f s\n",
           matrixTimer.setUpDataTime/1000.0);
 filePrint(f,  "         Make Subdomains               time: %14.5f s\n",
           matrixTimer.makeSubDomains/1000.0);
/*
 filePrint(f,  "            Element to Node Connectivity   : %31.3f Mb\n", 
           memoryElemToNode*byteToMb);
 filePrint(f,  "            Sub. to Node Connectivity      : %31.3f Mb\n", 
           memorySubToNode*byteToMb);
 filePrint(f,  "            Node to Sub. Connectivity      : %31.3f Mb\n", 
           memoryNodeToSub*byteToMb);
 filePrint(f,  "            Sub. to Element Connectivity   : %31.3f Mb\n", 
           memorySubToElem*byteToMb);

 if(matrixTimer.memoryCPUMAP != 0)
   filePrint(f,  "            CPU MAP Connectivity           : %31.3f Mb\n", 
             memoryCPUMAP*byteToMb);
*/
 filePrint(f,  "         Distribute BCs                time: %14.5f s\n", 
           matrixTimer.distributeBCs/1000.0);
 filePrint(f,  "         Make Connectivities           time: %14.5f s\n", 
           matrixTimer.makeConnectivity/1000.0);
 filePrint(f,  "         Make Interface                time: %14.5f s\n", 
           matrixTimer.makeInterface/1000.0);
 filePrint(f,  "         Make Internal Information     time: %14.5f s\n",
           matrixTimer.makeInternalInfo/1000.0);
 if(timeGeom != 0.0) {
   filePrint(f,  "         Make Geometric Node States    time: %14.5f s\n",
             timeGeom/1000.0);
 }

 filePrint(f,"\n3. Total Matrix Processing             time: %14.5f s %14.3f Mb\n",
           subTotal[2]/1000.0,totMemSubMatrices*byteToMb);
 filePrint(f,"         Construct Sparse Matrices     time: %14.5f s\n",
           matrixTimer.constructTime/1000.0);
 filePrint(f,"         Form Element Matrices         time: %14.5f s\n",
           (matrixTimer.formTime+kelArrayTime+corotatorTime)/1000.0);
 filePrint(f,"         Assemble Element Matrices     time: %14.5f s\n\n",
           matrixTimer.assemble/1000.0);

 filePrint(f,"4. Total RHS Processing                time: %14.5f s %14.3f Mb\n\n",
           subTotal[3]/1000.0, buildRhsMax.memory*byteToMb);
 }

 // NOTE: parallel factor time is counted within coarse1Max and
 //  timers.coarse2
 coarse1Max     -= timers.pfactor;
 coarse1Min     -= timers.pfactor;
 coarse1Tot     -= timers.pfactor;
 timers.coarse2 -= timers.pfactor2;

 double coarseTime   = coarse1Max + timers.coarse2;

 long locMemCoarse = timers.memoryGtG + timers.memoryPCtFPC;
 long totMemCoarse = locMemCoarse;

 long localMemorySolve = timers.memoryProject1 + timers.memoryProject2
                            + 8*timers.preconditioner.getOverAll()->memory 
                            + sAndJTot.memory + timers.memoryOSet 
                            + timers.memoryDV;

 long totalMemorySolve = localMemorySolve;

 // This includes the memory required to store the rigid body modes.
 long localMemFactor = timers.memoryFactor 
                          + 8*timers.kMem.getOverAll()->memory;
 long totalMemFactor = localMemFactor;

 
 long localSolverMemory = localMemFactor + locMemCoarse 
                             + localMemorySolve;

 long totalSolverMemory = localSolverMemory;

 double factorTimeMax = timers.factor;
 double factorTimeMin = timers.factor;
 double factorTimeAvg = timers.factor;
 
 long totMemReortho = timers.memoryOSet;
 long locMemUsed = memoryUsed();
 long totMemUsed = locMemUsed;
 long totMemFeti = timers.memoryFETI;
#ifdef DISTRIBUTED
 totMemCoarse = timers.globalMemorySum(locMemCoarse);
 totalMemorySolve = timers.globalMemorySum(localMemorySolve);
 totalMemFactor = timers.globalMemorySum(localMemFactor);
 totalSolverMemory = timers.globalMemorySum(localSolverMemory);
 totMemReortho = timers.globalMemorySum(timers.memoryOSet);
 factorTimeMax = structCom->globalMax(timers.factor);
 factorTimeMin = structCom->globalMin(timers.factor);
 factorTimeAvg = structCom->globalSum(timers.factor);
 
 totMemUsed = timers.globalMemorySum(locMemUsed);
 totMemFeti = timers.globalMemorySum(timers.memoryFETI);
#endif
 factorTimeAvg /= numCPUs;

 if(f != 0) {

 filePrint(f,"5. Total Solver                        time: %14.5f s %14.3f Mb\n",
           subTotal[4]/1000.0, totalSolverMemory*byteToMb);
 if(sInfo.newmarkBeta != 0.0 && domain->solInfo().solvercntl->type != SolverSelection::Direct) { // none of this is relevant for explicit dynamics or MUMPS
   filePrint(f,"         Factor Subdomain Matrices     time: %14.5f s %14.3f Mb\n",
             factorTimeMax/1000.0, totalMemFactor*byteToMb);
   filePrint(f,"         Total Building  Coarse Pbs.   time: %14.5f s %14.3f Mb\n",
             coarseTime/1000.0, totMemCoarse*byteToMb );
   filePrint(f,"               1st Level Coarse Pb.    time: %14.5f s %14.3f Mb\n",
             coarse1Max/1000.0, timers.memoryGtG*byteToMb);
   filePrint(f,"               2nd Level Coarse Pb.    time: %14.5f s %14.3f Mb\n",
             timers.coarse2/1000.0, timers.memoryPCtFPC*byteToMb);
   filePrint(f,"         Total Paral. Fac. Coarse Pbs. time: %14.5f s %14.3f Mb\n",
             (parfac1Max+timers.pfactor2)/1000.0, 0.0);
   filePrint(f,"               1st Level Factor        time: %14.5f s %14.3f Mb\n",
             parfac1Max/1000.0, 0.0);
   filePrint(f,"               2nd Level Factor        time: %14.5f s %14.3f Mb\n",
             timers.pfactor2/1000.0, 0.0);
   filePrint(f,"         Total Solve loop              time: %14.5f s %14.3f Mb\n",
             timers.solve/1000.0, totalMemorySolve*byteToMb);
   filePrint(f,"               Project 1st Level       time: %14.5f s %14.3f Mb\n",
             timers.project/1000.0, timers.memoryProject1*byteToMb);
   filePrint(f,"                  Seq. Forw/Back       time: %14.5f s\n",
             forBack1Max/1000.0);
   filePrint(f,"               Project 2nd Level       time: %14.5f s %14.3f Mb\n",
             timers.project2/1000.0,timers.memoryProject2*byteToMb);
   filePrint(f,"                  Paral. Forw/Back     time: %14.5f s\n",
             timers.forBack2/1000.0);
   filePrint(f,"               Precondition            time: %14.5f s %14.3f Mb\n",
             precondMax/1000.0, 8*memoryPrecond*byteToMb);
   filePrint(f,"               Local Solve             time: %14.5f s %14.3f Mb\n",
             sAndJMaximum/1000.0, sAndJTot.memory*byteToMb);
   filePrint(f,"               Reorthogonalize         time: %14.5f s %14.3f Mb\n",
             orthoMax.time/1000.0, totMemReortho*byteToMb);
 }
 if(domain->solInfo().isDynam() || domain->solInfo().isNonLin())
   filePrint(f,"         Update States                 time: %14.5f s\n",
             matrixTimer.updateState/1000.0);
 if(domain->tdenforceFlag()){ // ACME search and forces for explicit dynamics
   filePrint(f,"         TD Enforcement                time: %14.5f s\n",
             tdenforceTime/1000.0);
   filePrint(f,"              detection                time: %14.5f s\n",
             contactSearchTime/1000.0);
   filePrint(f,"              enforcement              time: %14.5f s\n",
             contactForcesTime/1000.0);
   filePrint(f,"              surface update           time: %14.5f s\n",
             updateSurfsTime/1000.0);
 }

 filePrint(f,"\n6. Write Output Files                  time: %14.5f s %14.3f Mb\n",
           std::max(subTotal[5]-matrixTimer.sendFluidTime,0.)/1000.0, memoryOutput*byteToMb);

 double totalFluidComm = matrixTimer.receiveFluidTime+matrixTimer.sendFluidTime;
 if(domain->solInfo().aeroFlag >= 0) {
   filePrint(f,"\n7. Fluid Communication                 time: %14.5f s\n",
             totalFluidComm/1000.0);
   filePrint(f,"   Receive From Fluid                  time: %14.5f s\n",
             matrixTimer.receiveFluidTime/1000.0);
   filePrint(f,"   Send To Fluid                       time: %14.5f s\n\n",
             matrixTimer.sendFluidTime/1000.0);
 }

 // Compute the total time spent on this simulation

 double total = 0.0;
 int i;
 for(i=0; i<6; ++i)
   total += subTotal[i];

 long totalMemSimulation = totalSolverMemory + buildRhsTot.memory +
                                totMemSubMatrices + memoryOutput +       
                                totMemPreProcess + totalMemRead;

 if(domain->solInfo().aeroFlag >= 0) {
   filePrint(f,"\nTOTAL SIMULATION 1 (1+2+3+4+5+6)       time: %14.5f s %14.3f Mb\n",total/1000.0,totalMemSimulation*byteToMb);
   filePrint(f,"\nTOTAL SIMULATION 2 (1+2+3+4+5+6+7)     time: %14.5f s\n",(total+totalFluidComm)/1000.0);
 }
 else {
   filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)         time: %14.5f s %14.3f Mb\n",total/1000.0, totalMemSimulation*byteToMb);
 }

 // Output FETI solver information
 if(sInfo.solvercntl->type == SolverSelection::Feti || sInfo.solvercntl->type == SolverSelection::FetiLib) {
   filePrint(f,"\n***********************************************************"
             "********************\n");
   filePrint(f," ... FETI Monitoring ... \n");
   filePrint(f,"***********************************************************"
             "********************\n\n");

   filePrint(f,"1. Total Amount of Requested Memory        = %14.3f Mb\n\n",
             (double)totMemUsed*byteToMb);

   filePrint(f,"2. FETI Solver Amount of Requested Memory  = %14.3f Mb\n\n",
             totMemFeti*byteToMb);

   // Check whether we converged or not
   if(timers.converged == 1)
   filePrint(f,"3. Number of Iterations for Convergence    = %14d\n\n",
             timers.numIter);
   else if(timers.converged == 0) {
   filePrint(f,"3. Stagnation Occured After a # of iter    = %14d\n\n",
             timers.numIter);
   } else
   filePrint(f,"3. Maximum Number of Iterations reached    = %14d\n\n",
             timers.numIter);

   filePrint(f,"4. Relative Primal Error Reached           = %14.3e\n\n",
             timers.iterations[0].finalPrimal);

   filePrint(f,"5. Relative Dual Error Reached             = %14.3e\n\n",
             timers.iterations[0].finalDual);

   filePrint(f,"6. Size of 1st Level Coarse Problem        = %14d %14.3f Mb\n\n",
             timers.numRBMs,timers.memoryGtGsky*byteToMb);

   filePrint(f,"7. Size of 2nd Level Coarse Problem        = %14d %14.3f Mb\n\n",
             timers.numCRNs,timers.memoryPCtFPCmat*byteToMb);

   if(sInfo.getFetiInfo().local_cntl->subtype == 0)
   filePrint(f,"8. Total Memory Subdomain Skyline K        = %14.3f Mb\n\n",
             8.0*totMemSky*byteToMb);
   else if(sInfo.getFetiInfo().local_cntl->subtype == 1)
   filePrint(f,"8. Total Memory Subdomain Sparse K         = %14.3f Mb\n\n",
             8.0*totMemSparse*byteToMb);
   else
     filePrint(f,""); // if we have other subdomain solvers

  }
  filePrint(f,"\n***********************************************************"
              "********************\n");
 }
 
#ifdef DISTRIBUTED
 double mem1 = (double) timers.memoryGtG;
 mem1 = structCom->globalMin(mem1);
 long memGtGMin = (long) mem1;
#endif

#ifdef DISTRIBUTED
 mem1 = (double) timers.memoryGtG;
 mem1 = structCom->globalMax(mem1);
 long memGtGMax = (long) mem1;
#endif
 
#ifdef DISTRIBUTED
 mem1 = (double) timers.memoryGtG;
 mem1 = structCom->globalSum(mem1);
 long memGtGTot = (long) mem1;
#endif

#ifdef DISTRIBUTED
 mem1 = (double) memoryOutput;
 mem1 = structCom->globalSum(mem1);
 long memOutTot = (long) mem1;
#endif
 
#ifdef DISTRIBUTED
 mem1 = (double) memoryOutput;
 mem1 = structCom->globalMax(mem1);
 long memOutMax = (long) mem1;
#endif

#ifdef DISTRIBUTED
 mem1 = (double) memoryOutput;
 mem1 = structCom->globalMin(mem1);
 long memOutMin = (long) mem1;
#endif

 // This needs to be improved to count better the parallel/sequential time
 double coarse2Min  = timers.coarse2;
 double coarse2Tot  = timers.coarse2;
 double coarse2Max  = timers.coarse2;

 double parfac2Min  = timers.pfactor2;
 double parfac2Tot  = timers.pfactor2;
 double parfac2Max  = timers.pfactor2;

 // This needs to be improved to count better the parallel/sequential time
 double project2Min = timers.project2;
 double project2Tot = timers.project2;
 double project2Max = timers.project2;

 //filePrint(stderr,"Construct Total Memory %14.5f\n",constructTot.memory*byteToMb);

 //totMemSubMatrices = assembleTot.memory + constructTot.memory;

 // calc. min., avg., max.
 double tot1MinTime = matrixTimer.readTime + matrixTimer.readDecomp;
 double tot1AvgTime = matrixTimer.readTime + matrixTimer.readDecomp;
 double tot1MaxTime = matrixTimer.readTime + matrixTimer.readDecomp;

 double tot2MinTime = subTotal[1];
 double tot2AvgTime = subTotal[1];
 double tot2MaxTime = subTotal[1];
	 
 double setUpDataTimeMin = matrixTimer.setUpDataTime;
 double setUpDataTimeAvg = matrixTimer.setUpDataTime;
 double setUpDataTimeMax = matrixTimer.setUpDataTime;

 double makeSubDomainTimeMin = matrixTimer.makeSubDomains;
 double makeSubDomainTimeAvg = matrixTimer.makeSubDomains;
 double makeSubDomainTimeMax = matrixTimer.makeSubDomains;

 double distributeBCTimeMin = matrixTimer.distributeBCs;
 double distributeBCTimeAvg = matrixTimer.distributeBCs;
 double distributeBCTimeMax = matrixTimer.distributeBCs;

 double makeConnectivityTimeMin = matrixTimer.makeConnectivity;
 double makeConnectivityTimeAvg = matrixTimer.makeConnectivity;
 double makeConnectivityTimeMax = matrixTimer.makeConnectivity;

 double makeInterfaceTimeMin = matrixTimer.makeInterface;
 double makeInterfaceTimeAvg = matrixTimer.makeInterface;
 double makeInterfaceTimeMax = matrixTimer.makeInterface;
 
 double makeInternalTimeMin = matrixTimer.makeInternalInfo;
 double makeInternalTimeAvg = matrixTimer.makeInternalInfo;
 double makeInternalTimeMax = matrixTimer.makeInternalInfo;

#ifdef DISTRIBUTED
 tot1MinTime = structCom->globalMin(tot1MinTime);
 tot1AvgTime = (structCom->globalSum(tot1AvgTime))/numCPUs;
 tot1MaxTime = structCom->globalMax(tot1MaxTime);

 tot2MinTime = structCom->globalMin(subTotal[1]);
 tot2AvgTime = (structCom->globalSum(subTotal[1]))/numCPUs;
 tot2MaxTime= structCom->globalMax(subTotal[1]);

 setUpDataTimeMin = structCom->globalMin(matrixTimer.setUpDataTime);
 setUpDataTimeAvg = (structCom->globalSum(matrixTimer.setUpDataTime))/numCPUs;
 setUpDataTimeMax = structCom->globalMax(matrixTimer.setUpDataTime);

 makeSubDomainTimeMin = structCom->globalMin(matrixTimer.makeSubDomains);
 makeSubDomainTimeAvg = (structCom->globalSum(matrixTimer.makeSubDomains))/numCPUs;
 makeSubDomainTimeMax = structCom->globalMax(matrixTimer.makeSubDomains);

 distributeBCTimeMin = structCom->globalMin(matrixTimer.distributeBCs);
 distributeBCTimeAvg = (structCom->globalSum(matrixTimer.distributeBCs))/numCPUs;
 distributeBCTimeMax = structCom->globalMax(matrixTimer.distributeBCs);

 makeConnectivityTimeMin = structCom->globalMin(matrixTimer.makeConnectivity);
 makeConnectivityTimeAvg = (structCom->globalSum(matrixTimer.makeConnectivity))/numCPUs;
 makeConnectivityTimeMax = structCom->globalMax(matrixTimer.makeConnectivity);

 makeInterfaceTimeMin  = structCom->globalMin(matrixTimer.makeInterface);
 makeInterfaceTimeAvg = (structCom->globalSum(matrixTimer.makeInterface))/numCPUs;
 makeInterfaceTimeMax = structCom->globalMax(matrixTimer.makeInterface);

 makeInternalTimeMin = structCom->globalMin(matrixTimer.makeInternalInfo);
 makeInternalTimeAvg = (structCom->globalSum(matrixTimer.makeInternalInfo))/numCPUs;
 makeInternalTimeMax = structCom->globalMax(matrixTimer.makeInternalInfo);
#endif

 if(f) {

 filePrint(f," ... Detailed CPU Statistics (Seconds) ");
 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f,"\n                                             minimum      average      maximum\n");

 filePrint(f,"\n1. Total Read Input Files             : %12.4f %12.4f %12.4f\n",
           tot1MinTime/1000.0, tot1AvgTime/1000.0, tot1MaxTime/1000.0);
	
 filePrint(f,"\n2. Total Preprocessing                : %12.4f %12.4f %12.4f\n", 
           tot2MinTime/1000.0,tot2AvgTime/1000.0,tot2MaxTime/1000.0);
	 
 filePrint(f,"         Process Input Data           : %12.4f %12.4f %12.4f\n",
           setUpDataTimeMin/1000.0, setUpDataTimeAvg/1000.0, setUpDataTimeMax/1000.0);

 filePrint(f,"         Make Subdomains              : %12.4f %12.4f %12.4f\n", 
           makeSubDomainTimeMin/1000.0, makeSubDomainTimeAvg/1000.0, makeSubDomainTimeMax/1000.0);

 filePrint(f,"         Distribute BCs               : %12.4f %12.4f %12.4f\n", 
           distributeBCTimeMin/1000.0, distributeBCTimeAvg/1000.0, distributeBCTimeMax/1000.0);
 filePrint(f,"         Make Connectivities          : %12.4f %12.4f %12.4f\n",
           makeConnectivityTimeMin/1000.0, makeConnectivityTimeAvg/1000.0, makeConnectivityTimeMax/1000.0);
 filePrint(f,"         Make Interface               : %12.4f %12.4f %12.4f\n",
           makeInterfaceTimeMin/1000.0, makeInterfaceTimeAvg/1000.0, makeInterfaceTimeMax/1000.0);
 filePrint(f,"         Make Internal Information    : %12.4f %12.4f %12.4f\n",
           makeInternalTimeMin/1000.0, makeInternalTimeAvg/1000.0, makeInternalTimeMax/1000.0);

 double tot3MinTime = constructMin.time + assembleMin.time;
 double tot3AvgTime = constructAvg.time + assembleAvg.time;
 double tot3MaxTime = constructMax.time + assembleMax.time; 

 filePrint(f,"\n3. Total Matrix Processing            : %12.4f %12.4f %12.4f"
           "\n",tot3MinTime/1000.0,tot3AvgTime/1000.0,tot3MaxTime/1000.0);
 
 double tot4MinTime = buildRhsMin.time;
 double tot4AvgTime = buildRhsAvg.time;
 double tot4MaxTime = buildRhsMax.time;

 filePrint(f,"\n4. Total RHS Processing               : %12.4f %12.4f %12.4f"
           "\n\n",tot4MinTime/1000.0,tot4AvgTime/1000.0,tot4MaxTime/1000.0);

 // timers for 2nd level coarse problem

 double timeCoarseMin = coarse1Min + coarse2Min;
 double timeCoarseTot = coarse1Tot + coarse2Tot;
 double timeCoarseMax = coarse1Max + coarse2Max; 

 double timeParFacMin = parfac1Min + parfac2Min;
 double timeParFacTot = parfac1Tot + parfac2Tot;
 double timeParFacMax = parfac1Max + parfac2Max;

 double com1 = precondMax   - precMax.time;
 double com2 = sAndJMaximum - sAndJMax.time;
 
 double proj1Seq = timers.project - project1Max.time;

 double proj1Min = project1Min.time + proj1Seq;
 double proj1Avg = project1Avg.time + proj1Seq;
 double proj1Max = project1Max.time + proj1Seq;

 double solveMin =  proj1Min + (project2Min + forBack2Min) 
                 +  precMin.time + com1 + sAndJMin.time + com2 + orthoMin.time;
 double solveAvg =  proj1Avg + (project2Tot + forBack2Tot)
                 +  precAvg.time + com1 + sAndJAvg.time + com2 + orthoAvg.time;
 double solveMax =  proj1Max + (project2Max + forBack2Max) 
                 +  precMax.time + com1 + sAndJMax.time + com2 + orthoMax.time;

 double missingTime = timers.solve - solveMax;

 solveMin += missingTime;
 solveAvg += missingTime;
 solveMax += missingTime;

 double tot5MinTime   = timeCoarseMin + factorMin.time + timeParFacMin
                      + solveMin; 

 double tot5AvgTime   = factorAvg.time + timeCoarseTot
                      + timeParFacTot  + solveAvg;

 double tot5MaxTime   = timeCoarseMax + factorMax.time + timeParFacMax
                      + solveMax;
   
 // double tot5MaxTime = solutionTime;

 filePrint(f,"5. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
           tot5MinTime/1000.0, tot5AvgTime/1000.0, tot5MaxTime/1000.0);

 if(domain->solInfo().solvercntl->type == SolverSelection::Feti) {
   filePrint(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
             factorTimeMin/1000.0,factorTimeAvg/1000.0,factorTimeMax/1000.0);

   filePrint(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
             timeCoarseMin/1000.0,timeCoarseTot/1000.0,timeCoarseMax/1000.0);
   filePrint(f,"               1st Level Coarse Pb.   : %12.4f %12.4f %12.4f\n",
             coarse1Min/1000.0,coarse1Tot/1000.0,coarse1Max/1000.0);
   filePrint(f,"               2nd Level Coarse Pb.   : %12.4f %12.4f %12.4f\n\n",
             coarse2Min/1000.0,coarse2Tot/1000.0,coarse2Max/1000.0);

   filePrint(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
             timeParFacMin/1000.0,timeParFacTot/1000.0,timeParFacMax/1000.0);
   filePrint(f,"               1st Level Factor       : %12.4f %12.4f %12.4f\n",
             parfac1Min/1000.0,parfac1Tot/1000.0,parfac1Max/1000.0);
   filePrint(f,"               2nd Level Factor       : %12.4f %12.4f %12.4f\n\n",
             parfac2Min/1000.0,parfac2Tot/1000.0,parfac2Max/1000.0);

   filePrint(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
             solveMin/1000.0,solveAvg/1000.0,solveMax/1000.0);
   filePrint(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
             proj1Min/1000.0,proj1Avg/1000.0,proj1Max/1000.0);
   filePrint(f,"                  Seq. Forw/Back      : %12.4f %12.4f %12.4f\n",
             forBack1Min/1000.0,forBack1Tot/1000.0, forBack1Max/1000.0);
   filePrint(f,"               Project 2nd Level      : %12.4f %12.4f %12.4f\n",
             project2Min/1000.0,project2Tot/(1000.0),project2Max/1000.0);
   filePrint(f,"                  Paral. Forw/Back    : %12.4f %12.4f %12.4f\n",
             forBack2Min/1000.0,forBack2Tot/(1000.0), forBack2Max/1000.0);

   filePrint(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
             (precMin.time+com1)/1000.0,(precAvg.time+com1)/1000.0,
             precondMax/1000.0);

   filePrint(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
             (sAndJMin.time+com2)/1000.0,(sAndJAvg.time+com2)/1000.0,
             sAndJMaximum/1000.0);
   filePrint(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n\n",
             orthoMin.time/1000.0,orthoAvg.time/1000.0,orthoMax.time/1000.0);
 }

 double tot6MinTime = output;
 double tot6AvgTime = output;
 double tot6MaxTime = output;

 filePrint(f,"6. Write Output Files                 : %12.4f %12.4f %12.4f\n",
           tot6MinTime/1000.0,tot6AvgTime/1000.0,tot6MaxTime/1000.0);

 double timeSimMin = tot1MinTime + tot2MinTime + tot3MinTime + tot4MinTime
                   + tot5MinTime + tot6MinTime;
 double timeSimAvg = tot1AvgTime + tot2AvgTime + tot3AvgTime + tot4AvgTime
                   + tot5AvgTime + tot6AvgTime;
 double timeSimMax = tot1MaxTime + tot2MaxTime + tot3MaxTime + tot4MaxTime
                   + tot5MaxTime + tot6MaxTime;
 
 filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)        : %12.4f %12.4f %12.4f\n",
           timeSimMin/1000.0, timeSimAvg/1000.0,timeSimMax/1000.0);
 }

#ifdef DISTRIBUTED
 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f," ...  Detailed Memory Statistics (Megabytes)");
 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f,"\n                                             minimum      average      maximum\n");

 // get overall prec memory stats
 TimeData *precOverall = timers.preconditioner.getOverAll();

 long locSubMatMem = timers.memorySubMatrices - 8*timers.kMem.getOverAll()->memory - 8*precOverall->memory;
 double numMPI = 1.0;

 long readMemoryMin = localMemRead;
 long readMemoryMax = localMemRead;
 long preProcessMemoryMin = memoryPreProcess + matrixTimer.memorySetUp;
 long preProcessMemoryMax = memoryPreProcess + matrixTimer.memorySetUp;
 long subMatMemoryMin = locSubMatMem;
 long subMatMemoryMax = locSubMatMem;
 long solverMemoryMin = localSolverMemory;
 long solverMemoryMax = localSolverMemory;
 long factorMemoryMin = localMemFactor;
 long factorMemoryMax = localMemFactor;
 long coarseMemoryMin = locMemCoarse;
 long coarseMemoryMax = locMemCoarse;
 long precMemoryMin = precOverall->memory;
 long precMemoryAvg = precOverall->memory;
 long precMemoryMax = precOverall->memory;
 long orthoMemoryMin = timers.memoryOSet;
 long orthoMemoryMax = timers.memoryOSet;


 numMPI = structCom->globalSum(numMPI);
 readMemoryMin = timers.globalMemoryMin(localMemRead);
 readMemoryMax = timers.globalMemoryMax(localMemRead);
 preProcessMemoryMin = timers.globalMemoryMin(preProcessMemoryMin);
 preProcessMemoryMax = timers.globalMemoryMax(preProcessMemoryMax);
 subMatMemoryMin = timers.globalMemoryMin(locSubMatMem);
 subMatMemoryMax = timers.globalMemoryMax(locSubMatMem);
 solverMemoryMin = timers.globalMemoryMin(localSolverMemory);
 solverMemoryMax = timers.globalMemoryMax(localSolverMemory);
 factorMemoryMin = timers.globalMemoryMin(localMemFactor);
 factorMemoryMax  = timers.globalMemoryMax(localMemFactor);
 coarseMemoryMin = timers.globalMemoryMin(locMemCoarse);
 coarseMemoryMax = timers.globalMemoryMax(locMemCoarse);
 precMemoryMin = timers.globalMemoryMin(precOverall->memory);
 precMemoryAvg = timers.globalMemorySum(precOverall->memory);
 precMemoryMax = timers.globalMemoryMax(precOverall->memory);
 orthoMemoryMin = timers.globalMemoryMin(timers.memoryOSet);
 orthoMemoryMax = timers.globalMemoryMax(timers.memoryOSet);


 long readMemoryAvg = totalMemRead/numMPI;
 long preProcessMemoryAvg = totMemPreProcess/numMPI;
 long subMatMemoryAvg = totMemSubMatrices/numMPI;
 long solverMemoryAvg = totalSolverMemory/numMPI;
 long factorMemoryAvg = totalMemFactor/numMPI;
 long orthoMemoryAvg = totMemReortho / numMPI;
 long coarseMemoryAvg = totMemCoarse / numMPI;
 precMemoryAvg /= numMPI;


 filePrint(f,"\n1. Total Read Input Files             : %12.4f %12.4f %12.4f\n",
           readMemoryMin*byteToMb, readMemoryAvg*byteToMb, readMemoryMax*byteToMb);
/*
 long tot3Min = (constructMin.memory + assembleMin.memory);
 long tot3Avg = (constructAvg.memory + assembleAvg.memory);
 long tot3Max = (constructMax.memory + assembleMax.memory);
*/
 long tot3Min = subMatMemoryMin;
 long tot3Avg = subMatMemoryAvg;
 long tot3Max = subMatMemoryMax;
 
 filePrint(f,"\n2. Total Preprocessing                : %12.4f %12.4f %12.4f\n", 
           preProcessMemoryMin*byteToMb,preProcessMemoryAvg*byteToMb,
           preProcessMemoryMax*byteToMb);
	 
 filePrint(f,"\n3. Total Matrix Processing            : %12.4f %12.4f %12.4f"
           "\n", tot3Min*byteToMb,tot3Avg*byteToMb,tot3Max*byteToMb);
 
 long tot4Min = buildRhsMin.memory;
 long tot4Avg = buildRhsAvg.memory;
 long tot4Max = buildRhsMax.memory;

 filePrint(f,"\n4. Total RHS Processing               : %12.4f %12.4f %12.4f"
           "\n\n",tot4Min*byteToMb,tot4Avg*byteToMb,tot4Max*byteToMb);
/*
 long memCoarseMin = memGtGMin;
 long memCoarseTot = memGtGTot;
 long memCoarseMax = memGtGMax;

 long tot5Min = factorMin.memory + memCoarseMin + orthoMin.memory;
 long tot5Avg = factorAvg.memory + memCoarseTot/numCPUs + orthoAvg.memory;
 long tot5Max = factorMax.memory + memCoarseMax + orthoMax.memory;
*/
 long tot5Min = solverMemoryMin;
 long tot5Avg = solverMemoryAvg;
 long tot5Max = solverMemoryMax;

 filePrint(f,"5. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
           tot5Min*byteToMb, tot5Avg*byteToMb, tot5Max*byteToMb);

 if(domain->solInfo().solvercntl->type == SolverSelection::Feti) {
   filePrint(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
             factorMemoryMin*byteToMb, factorMemoryAvg*byteToMb,
             factorMemoryMax*byteToMb);
 
   filePrint(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
             coarseMemoryMin*byteToMb, coarseMemoryAvg*byteToMb, coarseMemoryMax*byteToMb);

   filePrint(f,"               1st Level Coarse Pb.   : %12.4f %12.4f %12.4f\n",
             memGtGMin*byteToMb, (memGtGTot/numCPUs)*byteToMb, memGtGMax*byteToMb);

   filePrint(f,"               2nd Level Coarse Pb.   : %12.4f %12.4f %12.4f\n\n",
             0.0,0.0,0.0);
   filePrint(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
             0.0,0.0,0.0);
   filePrint(f,"               1st Level Factor       : %12.4f %12.4f %12.4f\n",
             0.0,0.0,0.0);
   filePrint(f,"               2nd Level Factor       : %12.4f %12.4f %12.4f"
             "\n\n",0.0, 0.0,0.0);
   filePrint(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
             solverMemoryMin*byteToMb, solverMemoryAvg*byteToMb,
             solverMemoryMax*byteToMb);
   filePrint(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
             0.0,0.0,0.0);
   filePrint(f,"               Project 2nd Level      : %12.4f %12.4f %12.4f\n",
             0.0,0.0,0.0);
   filePrint(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
             8*precMemoryMin*byteToMb, 8*precMemoryAvg*byteToMb,
             8*precMemoryMax*byteToMb);
   filePrint(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
             sAndJMin.memory*byteToMb,sAndJAvg.memory*byteToMb,
             sAndJMax.memory*byteToMb);
   filePrint(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n\n",
             orthoMemoryMin*byteToMb, orthoMemoryAvg*byteToMb,
             orthoMemoryMax*byteToMb);
 }

 long tot6Min = memOutMin;
 long tot6Avg = memOutTot/numCPUs;
 long tot6Max = memOutMax;

 filePrint(f,"6. Write Output Files                 : %12.4f %12.4f %12.4f\n",
           tot6Min*byteToMb,tot6Avg*byteToMb,tot6Max*byteToMb);

 long totSimMin = readMemoryMin + tot3Min + preProcessMemoryMin + tot4Min + tot5Min + tot6Min;
 long totSimAvg = readMemoryAvg + tot3Avg + preProcessMemoryAvg + tot4Avg + tot5Avg + tot6Avg;
 long totSimMax = readMemoryMax + tot3Max + preProcessMemoryMax + tot4Max + tot5Max + tot6Max;
 
 filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)        : %12.4f %12.4f %12.4f\n",
           totSimMin*byteToMb, totSimAvg*byteToMb, totSimMax*byteToMb);
 filePrint(f,"\n***********************************************************"
           "********************\n");
 if(f) fclose(f);
#endif
}
