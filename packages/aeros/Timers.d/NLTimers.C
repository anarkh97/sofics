#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Timers.d/Timing.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/MatrixTimers.h>
#include <Driver.d/Domain.h>
#include <Threads.d/Paral.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include <Driver.d/GeoSource.h>

extern long totMemSparse;
extern long totMemSky;

extern const char* message[];
extern const char* renumMessage[];
extern const char* yesno[];
extern const char* problemType[];
extern const char* scalingMessage[];
extern const char* rbmMessage[];
extern const char* precMessage[];
extern const char* precMessage[];
extern const char* projectMessage[];
extern const char* subSolverMessage[];
extern const char* precSolverMessage[];
extern const char* solverMessage[];

const double byteToMb    = (1.0 / oneMegaByte);

// ---------------------------------
// For NONLINEAR ANALYSIS TIMINGS

#ifndef SALINAS
static const char* KrylovMessage[] = {
" ",
"Krylov Preconditioner",
"Krylov Preconditioner (rebuilt at each Load Step)"
};
#endif


void
StaticTimers::printTimers(Domain* domain, Timings& timers, double solveTime)
{
#ifndef SALINAS
 MatrixTimers &times = domain->getTimers();
 SolverInfo   &sInfo = domain->solInfo();
 ControlInfo  &cinfo = geoSource->getCheckFileInfo()[0];
 
 int numCPUs       = threadManager->numThr();
 int numThreads = numCPUs;
 int numSubdomains = timers.numSubdomains;
#ifdef DISTRIBUTED
 numCPUs = structCom->numCPUs();
 numThreads  = structCom->numCPUs();
 numSubdomains = structCom->globalSum(numSubdomains);

 if(structCom->myID() == 0)
#endif
 openTimingFile(cinfo);

 int numele   = domain->numElements();
 int numnod   = domain->numNode();
 int numCon   = domain->nDirichlet();
 int numdof   = domain->numDofs();
 int numloads = domain->nNeumann();
 int numFiles = geoSource->getNumOutInfo();

 // Do a bunch of summations if distributed
 double coarse1Max  = timers.coarse1;
 double coarse1Tot  = timers.coarse1;
 double coarse1Min  = timers.coarse1;

 double parfac1Max  = timers.pfactor;
 double parfac1Tot  = timers.pfactor;
 double parfac1Min  = timers.pfactor;
 double parfac2Max  = timers.pfactor2;
 double parfac2Tot  = timers.pfactor2;
 double parfac2Min  = timers.pfactor2;
 
 double forBack1Max = timers.forBack;
 double forBack1Tot = timers.forBack;
 double forBack1Min = timers.forBack;
 double forBack2Max = timers.forBack2;
 double forBack2Tot = timers.forBack2;
 double forBack2Min = timers.forBack2;

 double precondMin   = timers.precond;
 double precondMax   = timers.precond;
 double precondAvg   = timers.precond;
 double sAndJMaximum = timers.sAndJ;
 double sAndJMin     = timers.sAndJ;
 double sAndJAvg     = timers.sAndJ;
 double reOrthoMin   = timers.reOrtho;
 double reOrthoAvg   = timers.reOrtho;
 double reOrthoMax   = timers.reOrtho;

 long localMemRead = times.memoryParse + times.memorySetUp;
 long totalMemRead = localMemRead;
 long totMemPreProcess = memoryPreProcess;
 long memorySetUp = times.memorySetUp;
 long memorySubMatrices = timers.memorySubMatrices;

#ifdef DISTRIBUTED
/* these are not used
 long memorySubdomain = times.memorySubdomain;
 long memoryElemToNode = times.memoryElemToNode;
 long memorySubToNode = times.memorySubToNode;
 long memoryNodeToSub = times.memoryNodeToSub;
 long memorySubToElem = times.memorySubToElem;
 long memoryCPUMAP = times.memoryCPUMAP;
 long memoryDistBC = times.memoryDistBC;
 long memoryConnect = times.memoryConnect;
 long memoryInterface = times.memoryInterface;
 long memoryInternal = times.memoryInternal;
*/
 long memoryOutputMin = memoryOutput;
 long memoryOutputAvg = memoryOutput;
 long memoryOutputMax = memoryOutput;
#endif
 double totalRead = times.readTime + times.readDecomp;

 // Need to identify the other preprocessing steps
 double localPreProcess = preProcess + times.setUpDataTime + corotatorTime 
                        + kelArrayTime + timeGeom - times.readDecomp;
 double totalPreProcess = localPreProcess;

 double totalMatricesProcess = timers.assemble + timers.constructMatrices;
 double coarseTime   = timers.coarse1 + timers.coarse2;
 double localSolutionTime = (sInfo.solvercntl->type == SolverSelection::Direct) ? solveTime : timers.solve + timers.factor + coarseTime
                                                            + timers.pfactor + timers.pfactor2 
                                                            - totalMatricesProcess;
 double solutionTimeMin = localSolutionTime;
 double solutionTimeAvg = localSolutionTime;
 double solutionTimeMax = localSolutionTime;

 double solveLoopMin = timers.solve;
 double solveLoopAvg = timers.solve;
 double solveLoopMax = timers.solve;

 double project1TimeMin = timers.project;
 double project1TimeAvg = timers.project;
 double project1TimeMax = timers.project;
 double project2TimeMin = timers.project2;
 double project2TimeAvg = timers.project2;
 double project2TimeMax = timers.project2;

 double coarse2TimeMin = timers.coarse2;
 double coarse2TimeAvg = timers.coarse2;
 double coarse2TimeMax = timers.coarse2;

 double formRhsMin = formRhs;
 double formRhsAvg = formRhs;
 double formRhsMax = formRhs;

 double outputTimeMax = output;
 double outputTimeAvg = output;
 double outputTimeMin = output;

#ifdef DISTRIBUTED
 totalMemRead = timers.globalMemorySum(localMemRead);
 memorySetUp = timers.globalMemorySum(times.memorySetUp);
 memorySubMatrices = timers.globalMemorySum(timers.memorySubMatrices);
 totMemPreProcess = timers.globalMemorySum(memoryPreProcess);
/* these are not used
 memorySubdomain = timers.globalMemorySum(times.memorySubdomain);
 memoryElemToNode = timers.globalMemorySum(times.memoryElemToNode);
 memorySubToNode = timers.globalMemorySum(times.memorySubToNode);
 memoryNodeToSub = timers.globalMemorySum(times.memoryNodeToSub);
 memorySubToElem = timers.globalMemorySum(times.memorySubToElem);
 memoryCPUMAP = timers.globalMemorySum(times.memoryCPUMAP);
 memoryDistBC = timers.globalMemorySum(times.memoryDistBC);
 memoryConnect = timers.globalMemorySum(times.memoryConnect);
 memoryInterface = timers.globalMemorySum(times.memoryInterface);
 memoryInternal = timers.globalMemorySum(times.memoryInternal);
*/
 memoryOutputMin = timers.globalMemoryMin(memoryOutput);
 memoryOutputAvg = timers.globalMemorySum(memoryOutput)/numCPUs;
 memoryOutputMax = timers.globalMemoryMax(memoryOutput);

 totalRead = structCom->globalMax(totalRead);
 times.readTime = structCom->globalMax(times.readTime);
 times.readDecomp = structCom->globalMax(times.readDecomp);
 totalPreProcess = structCom->globalMax(localPreProcess);
 totalMatricesProcess = structCom->globalMax(totalMatricesProcess);
 solutionTimeMax = structCom->globalMax(localSolutionTime);
 solutionTimeMin = structCom->globalMin(localSolutionTime);
 solutionTimeAvg = (structCom->globalSum(localSolutionTime))/numCPUs;
 solveLoopMin = structCom->globalMax(timers.solve);
 solveLoopAvg = (structCom->globalSum(timers.solve))/numCPUs;
 solveLoopMax = structCom->globalMax(timers.solve);
 project1TimeMin = structCom->globalMax(timers.project);
 project1TimeAvg = (structCom->globalSum(timers.project))/numCPUs;
 project1TimeMax = structCom->globalMax(timers.project);
 project2TimeMin = structCom->globalMax(timers.project2);
 project2TimeAvg = (structCom->globalSum(timers.project2))/numCPUs;
 project2TimeMax = structCom->globalMax(timers.project2);
 coarse2TimeMin = structCom->globalMin(timers.coarse2);
 coarse2TimeAvg = (structCom->globalSum(timers.coarse2))/numCPUs;
 coarse2TimeMax = structCom->globalMax(timers.coarse2);

 outputTimeMin = structCom->globalMin(output);
 outputTimeAvg = structCom->globalSum(output)/numCPUs;
 outputTimeMax = structCom->globalMax(output);

 coarse1Max   = structCom->globalMax(timers.coarse1);
 coarse1Tot   = structCom->globalSum(timers.coarse1)/numCPUs;
 coarse1Min   = structCom->globalMin(timers.coarse1);
 parfac1Max   = structCom->globalMax(timers.pfactor);
 parfac1Tot   = structCom->globalSum(timers.pfactor)/numCPUs;
 parfac1Min   = structCom->globalMin(timers.pfactor);
 parfac2Max   = structCom->globalMax(timers.pfactor2);
 parfac2Tot   = structCom->globalSum(timers.pfactor2)/numCPUs;
 parfac2Min   = structCom->globalMin(timers.pfactor2);
 forBack1Max  = structCom->globalMax(timers.forBack);
 forBack1Tot  = structCom->globalSum(timers.forBack)/numCPUs;
 forBack1Min  = structCom->globalMin(timers.forBack);
 forBack2Max  = structCom->globalMax(timers.forBack2);
 forBack2Tot  = structCom->globalSum(timers.forBack2)/numCPUs;
 forBack2Min  = structCom->globalMin(timers.forBack2);
 precondMax   = structCom->globalMax(timers.precond);
 precondMin   = structCom->globalMin(timers.precond); 
 precondAvg   = (structCom->globalSum(timers.precond))/numCPUs; 
 sAndJMaximum = structCom->globalMax(timers.sAndJ);
 sAndJAvg     = structCom->globalSum(timers.sAndJ)/numCPUs;
 sAndJMin     = structCom->globalMin(timers.sAndJ);
 reOrthoMin   = structCom->globalMin(timers.reOrtho);
 reOrthoAvg   = structCom->globalSum(timers.reOrtho)/numCPUs;
 reOrthoMax   = structCom->globalMax(timers.reOrtho);

 formRhsMin   = structCom->globalMin(formRhs);
 formRhsAvg   = structCom->globalSum(formRhs)/numCPUs;
 formRhsMax   = structCom->globalMax(formRhs);
#endif

 totMemPreProcess += memorySetUp;

 long totMemSubMatrices = memorySubMatrices - 8*memoryK
                             - 8*memoryPrecond;


 if (f != 0)  {
 fprintf(f,"\n***********************************************************"
           "********************\n");
 fprintf(f," ... %s Problem Information ... \n",problemType[sInfo.probType]);
 fprintf(f,"***********************************************************"
           "********************\n\n");
 fprintf(f,"1. Number of Nodes                         = %14d\n\n",numnod);
 fprintf(f,"2. Number of Elements                      = %14d\n\n",numele);
 fprintf(f,"3. Number of Degrees of Freedom            = %14d\n",numdof);
 fprintf(f,"         Number of Constrained Dofs        = %14d\n",numCon);
 fprintf(f,"         Number of Unconstrained Dofs      = %14d\n\n",numdof-numCon);
 fprintf(f,"4. Number of Applied Loads                 = %14d\n\n",numloads);
 fprintf(f,"5. Number of Output Files                  = %14d\n\n",numFiles);
 fprintf(f,"6. Renumbering                             = %14s\n\n",
         renumMessage[sInfo.renum]);


 fprintf(f,"***********************************************************"
           "********************\n");
 fprintf(f," ... Solver Information ... \n");
 fprintf(f,"***********************************************************"
           "********************\n\n");

 if(sInfo.solvercntl->type == SolverSelection::Direct) {
   fprintf(f,"1. Mumps Sparse\n\n");
 }
 else {
   fprintf(f,"%s", solverMessage[sInfo.getFetiInfo().version]);
   fprintf(f,"         %s", precMessage[sInfo.getFetiInfo().precno]);
   fprintf(f,"         %s", scalingMessage[sInfo.getFetiInfo().scaling]);
   fprintf(f,"         %s", projectMessage[sInfo.getFetiInfo().nonLocalQ]);
   fprintf(f,"         %s", subSolverMessage[sInfo.getFetiInfo().local_cntl->subtype]);
   fprintf(f,"         %s", precSolverMessage[sInfo.getFetiInfo().kii_cntl->subtype]);

   if(sInfo.rbmflg == 0)
     fprintf(f,"         %s %29e\n",rbmMessage[sInfo.rbmflg],sInfo.solvercntl->trbm2);
   else
     fprintf(f,"         %s%17e %e\n",rbmMessage[sInfo.rbmflg],sInfo.tolsvd,sInfo.solvercntl->trbm);

   fprintf(f,"         Maximum Number of Iterations      = %14d\n",
           sInfo.getFetiInfo().maxiter());
   fprintf(f,"         Maximum Size of Reortho. Vectors  = %14d\n",
           sInfo.getFetiInfo().maxorth());
   fprintf(f,"         Tolerance for Convergence         = %14.3e\n",
           sInfo.getFetiInfo().tolerance());
   fprintf(f,"         %s\n", KrylovMessage[sInfo.getFetiInfo().nlPrecFlg]);
 }

 fprintf(f,"***********************************************************"
           "********************\n");
 fprintf(f," ... Timing Statistics for %d Threads and %d Subdomains ...\n",
                numThreads, numSubdomains);
 fprintf(f,"***********************************************************"
           "********************\n\n");

 // Begin timer output

 fprintf(f,"1. Total Read Input Files              time: %14.5f s %14.3f Mb\n", 
         totalRead/1000.0, totalMemRead*byteToMb);
 fprintf(f,"         Read Mesh                     time: %14.5f s\n",
         times.readTime/1000.0);
 fprintf(f,"         Read Mesh Partition           time: %14.5f s\n",
         times.readDecomp/1000.0);          


 fprintf(f,"\n2. Total Preprocessing                 time: %14.5f s %14.3f Mb\n",
         (totalPreProcess)/1000.0, totMemPreProcess*byteToMb);
 fprintf(f,  "         Process Input Data            time: %14.5f s\n",
         times.setUpDataTime/1000.0);
 fprintf(f,  "         Make Subdomains               time: %14.5f s\n",  
         times.makeSubDomains/1000.0);
 fprintf(f,  "         Distribute BCs                time: %14.5f s\n",  
         times.distributeBCs/1000.0);
 fprintf(f,  "         Make Connectivities           time: %14.5f s\n",  
         times.makeConnectivity/1000.0);
 fprintf(f,  "         Make Interface                time: %14.5f s\n",  
         times.makeInterface/1000.0);
 fprintf(f,  "         Make Internal Information     time: %14.5f s\n",  
         times.makeInternalInfo/1000.0);
 fprintf(f,  "         Make Geometric Node States    time: %14.5f s\n",
         timeGeom/1000.0);


 fprintf(f,"\n3. Total Matrix Processing             time: %14.5f s %14.3f Mb\n",
         (timers.constructMatrices+times.constructTime+times.assemble+kelArrayTime+corotatorTime+times.formTime)/1000.0,
         totMemSubMatrices*byteToMb);

 fprintf(f,"         Construct Sparse Matrices     time: %14.5f s\n",
         (timers.constructMatrices+times.constructTime)/1000.0);
 fprintf(f,"         Form Element Matrices         time: %14.5f s\n",
         (times.formTime+kelArrayTime+corotatorTime)/1000.0);
 fprintf(f,"         Assemble Element Matrices     time: %14.5f s\n",
         timers.assemble/1000.0);


 fprintf(f,"\n4. Total RHS Processing                time: %14.5f s %14.3f Mb\n",
         (formRhsMax-times.receiveFluidTime)/1000.0, timers.memoryLoad *byteToMb);
 }  // end if
   

 long locMemCoarse = timers.memoryGtG + timers.memoryPCtFPC;
 long totMemCoarse = locMemCoarse;
 long localMemorySolve = timers.memoryProject1 
                            + timers.memoryProject2
                            + timers.memorySAndJ + timers.memoryOSet 
                            + timers.memoryDV 
                            + 8*timers.preconditioner.getOverAll()->memory;

 long totalMemorySolve = localMemorySolve;

 long localMemFactor = timers.memoryFactor + 8*timers.kMem.getOverAll()->memory;
 long totalMemFactor = localMemFactor;
 

 long localSolverMemory = localMemFactor + locMemCoarse +
                               localMemorySolve;
 long totalSolverMemory = localSolverMemory;
 double factorTimeMax = timers.factor;
 double factorTimeMin = timers.factor;
 double factorTimeAvg = timers.factor;

 long locMemUsed = memoryUsed();
 long totMemUsed = locMemUsed;
 long totMemFeti = timers.memoryFETI;
 long totMemGtGAvg = timers.memoryGtG;
#ifdef DISTRIBUTED
 long totMemReortho = timers.memoryOSet;
 long totMemGtGMin = timers.memoryGtG;
 long totMemGtGMax = timers.memoryGtG;
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

 totMemGtGMin = timers.globalMemoryMin(timers.memoryGtG);
 totMemGtGAvg = timers.globalMemorySum(timers.memoryGtG)/numCPUs;
 totMemGtGMax = timers.globalMemoryMax(timers.memoryGtG);

#endif
 factorTimeAvg /= numCPUs;

 // Timers for Rebuilding
 double totalReBuild = rebuild + buildStiffAndForce + timers.reBuildPrec;
/*
                       - timers.reBuildGtG - timers.reBuildPCtFPC 
                       - timers.factor;
*/
#ifdef DISTRIBUTED
 totalReBuild = structCom->globalMax(totalReBuild);
#endif

 if (f != 0) {
 fprintf(f,"\n5. Total Solver                        time: %14.5f s %14.3f Mb\n",
         solutionTimeMax/1000.0, totalSolverMemory*byteToMb);
 if(sInfo.solvercntl->type == SolverSelection::Feti) {
   fprintf(f,"         Factor Subdomain Matrices     time: %14.5f s %14.3f Mb\n\n",
           timers.factor/1000.0, totalMemFactor*byteToMb);
   fprintf(f,"         Total Building  Coarse Pbs.   time: %14.5f s %14.3f\n",
           coarseTime/1000.0, totMemCoarse *byteToMb );
   fprintf(f,"               1st Level Coarse Pb.    time: %14.5f s %14.3f Mb\n",
           timers.coarse1/1000.0, totMemGtGAvg*numCPUs*byteToMb);
   fprintf(f,"               2nd Level Coarse Pb.    time: %14.5f s %14.3f Mb\n\n",
           coarse2TimeMax/1000.0, timers.memoryPCtFPC*byteToMb);
   fprintf(f,"         Total Paral. Fac. Coarse Pbs. time: %14.5f s %14.3f Mb\n",
           (parfac1Max+parfac2Max)/1000.0, 0.0);
   fprintf(f,"               1st Level Factor        time: %14.5f s %14.3f Mb\n",
           parfac1Max/1000.0, 0.0);
   fprintf(f,"               2nd Level Factor        time: %14.5f s %14.3f Mb\n\n",
           parfac2Max/1000.0, 0.0);
   fprintf(f,"         Total Solve loop              time: %14.5f s %14.3f Mb\n",
           solveLoopMax/1000.0, totalMemorySolve*byteToMb);
   fprintf(f,"               Project 1st Level       time: %14.5f s %14.3f Mb\n",
           project1TimeMax/1000.0, timers.memoryProject1*byteToMb);
   fprintf(f,"                  Seq. Forw/Back       time: %14.5f s\n",
           forBack1Max/1000.0);
   fprintf(f,"               Project 2nd Level       time: %14.5f s %14.3f Mb\n",
           project2TimeMax/1000.0,timers.memoryProject2*byteToMb);
   fprintf(f,"                  Paral. Forw/Back     time: %14.5f s\n",
           forBack2Max/1000.0);
   fprintf(f,"               Precondition            time: %14.5f s %14.3f Mb\n",
           precondMax/1000.0, 8*memoryPrecond*byteToMb);
   fprintf(f,"               Krylov Precondition     time: %14.5f s\n",
           timers.nlPreCond/1000.0);
   fprintf(f,"               Local Solve             time: %14.5f s %14.3f Mb\n",
           timers.sAndJ/1000.0, timers.memorySAndJ*byteToMb);
   fprintf(f,"               Reorthogonalize         time: %14.5f s %14.3f Mb\n",
           timers.reOrtho/1000.0, timers.memoryOSet*byteToMb);
 }
 else {
   fprintf(f,"         Factor Matrix                 time: %14.5f s\n",
           (times.factor)/1000.0);

   fprintf(f,"         Solve (Forward/Backward)      time: %14.5f s\n",
           solveTime/1000.0);
 }

 if(domain->solInfo().isDynam() || domain->solInfo().isNonLin())
   fprintf(f,"         Update States                 time: %14.5f s\n",
           times.updateState/1000.0);

 
 fprintf(f,"\n6. Write Output Files                  time: %14.5f s %14.3f Mb\n",
         (outputTimeMax-times.sendFluidTime)/1000.0, memoryOutput*byteToMb);

/*
 fprintf(f,"\n7. Total Rebuild                       time: %14.5f s\n\n",
         totalReBuild/1000.0);
 fprintf(f,"         Element Stiffness Matrices    time: %14.5f s\n",
         buildStiffAndForce/1000.0);

 if(sInfo.type == 2) {
   fprintf(f,"         FETI Preconditioner           time: %14.5f s\n",
           timers.reBuildPrec/1000.0);
   fprintf(f,"         Error Estimator               time: %14.5f s\n",
           timers.reBuildError/1000.0);
 }
*/

 double totalFluidComm = times.receiveFluidTime+times.sendFluidTime;
 if(domain->solInfo().aeroFlag >= 0) {
   fprintf(f,"\n7. Fluid Communication                 time: %14.5f s\n",
           totalFluidComm/1000.0);
   fprintf(f,"   Receive From Fluid                  time: %14.5f s\n",
           times.receiveFluidTime/1000.0);
   fprintf(f,"   Send To Fluid                       time: %14.5f s\n\n",
           times.sendFluidTime/1000.0);
 }

 // Compute the total time spent on this simulation
 double total = totalRead + totalPreProcess + solutionTimeMax + formRhsMax + output + 
                timers.assemble + timers.constructMatrices + totalReBuild;

 long totalMemSimulation = totalSolverMemory + timers.memoryLoad
                              + totMemSubMatrices + memoryOutput
                              + totMemPreProcess  + totalMemRead;

 if(domain->solInfo().aeroFlag >= 0) {
   fprintf(f,"\nTOTAL SIMULATION 1 (1+2+3+4+5+6)       time: %14.5f s %14.3f Mb\n",(total-totalFluidComm)/1000.0,totalMemSimulation*byteToMb);
   fprintf(f,"\nTOTAL SIMULATION 2 (1+2+3+4+5+6+7)     time: %14.5f s\n",total/1000.0);
 }
 else {
   fprintf(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)         time: %14.5f s %14.3f Mb\n",
           total/1000.0, totalMemSimulation*byteToMb);
 }

 // Output FETI solver information
 if(sInfo.solvercntl->type == SolverSelection::Feti) {

   fprintf(f,"\n***********************************************************"
             "********************\n");
   fprintf(f," ... FETI Monitoring ... \n");
   fprintf(f,"***********************************************************"
             "********************\n\n");

   fprintf(f,"1. Total Amount of Requested Memory        = %14.3f Mb\n\n",
           totMemUsed*byteToMb);

   fprintf(f,"2. FETI Solver Amount of Requested Memory  = %14.3f Mb\n\n",
           totMemFeti*byteToMb);

   fprintf(f,"3. Total Number of Iter. for Convergence   = %14d\n\n",
           timers.numIter);

   fprintf(f,"4. Size of 1st Level Coarse Problem        = %14d %14.3f Mb\n\n",
           timers.numRBMs,timers.memoryGtG*byteToMb);

   fprintf(f,"5. Size of 2nd Level Coarse Problem        = %14d %14.3f Mb\n",
           timers.numCRNs,timers.memoryPCtFPC*byteToMb);
 }

 //if(sInfo.getFetiInfo().local_cntl->subtype == 0)
 //  fprintf(f,"Memory Skyline = %14.3f Mb\n", 8.0*totMemSky*byteToMb );
 //else if(sInfo.getFetiInfo().local_cntl->subtype == 1)
 //  fprintf(f,"Memory Sparse  = %14.3f Mb\n", 8.0*totMemSparse*byteToMb );

 if((domain->probType() == SolverInfo::NonLinStatic ||
     domain->probType() == SolverInfo::NonLinDynam  ||
     domain->probType() == SolverInfo::ArcLength) && sInfo.solvercntl->type == SolverSelection::Feti) {

   fprintf(f,"\n***********************************************************"
             "********************\n");
   fprintf(f," ... Nonlinear Statics Monitoring ... \n");
   fprintf(f,"***********************************************************"
             "********************\n\n");
   fprintf(f,"Rebuild solver                  %d\n",sInfo.getNLInfo().updateK);
   fprintf(f,"Rebuild preconditioner          %d\n",sInfo.getFetiInfo().nPrecond());
   fprintf(f,"Max number of Newton iterations %d\n",sInfo.getNLInfo().maxiter);
   fprintf(f,"Nonlinear tolerance = %e\n",sInfo.getNLInfo().tolRes);
   fprintf(f,"Delta Lambda        = %e\n",sInfo.getNLInfo().dlambda);
   fprintf(f,"Maximum Lambda      = %e\n",sInfo.getNLInfo().maxLambda);
   fprintf(f,"Fit Algorithms: Shell %d Beam %d\n",sInfo.getNLInfo().fitAlgShell,
                                                  sInfo.getNLInfo().fitAlgBeam );

   // If arclength, print some more information about algorithm parameters
   fprintf(f,"\nIter\tFETI\tRebuild\tRebuild\tPrimal\t\tDual\t  Stagnation\n");
   fprintf(f,"\tIter\tTangent\tKrylov\tError\t\tError\n");

   int i;
   for(i=0; i<timers.numSystems; ++i) {
     if(norms[i].relativeDv == 1.0) 
       fprintf(f,
       "--------------------------------------------------------------------\n");

     fprintf(f,"%3d\t%d\t%s\t%s\t%10.4e\t%10.4e\t%s\n",i+1,
                                  timers.iterations[i].numFetiIter,
				  yesno[norms[i].rebuildTang],
			          yesno[timers.iterations[i].rebuildKrylov],
                                  timers.iterations[i].finalPrimal,
                                  timers.iterations[i].finalDual,
                                  yesno[timers.iterations[i].stagnated]);
   }
   fprintf(f,"------------------------------------------------"
             "--------------------\n");
   fprintf(f,"Total\t%d\n",timers.numIter);
   if(timers.numSystems != 0) fprintf(f,"Average\t%d\n",timers.numIter/(timers.numSystems));

   fprintf(f,"\nIter\tDv\t        Relative Dv\tResidual\tRelative Res\n");
   for(i=0; i<timers.numSystems; ++i) {
     if(norms[i].relativeDv == 1.0) 
       fprintf(f,
       "--------------------------------------------------------------------\n");
     fprintf(f,"%3d%16e\t%e\t%e\t%e\n",i+1,norms[i].normDv,norms[i].relativeDv,
             norms[i].normRes,norms[i].relativeRes);
   }
   fprintf(f,"------------------------------------------------"
             "--------------------\n");
 }

 } // end if(f != 0)

 // calc. min., avg., max.
 double tot1MinTime = times.readTime + times.readDecomp;
 double tot1AvgTime = times.readTime + times.readDecomp;
 double tot1MaxTime = times.readTime + times.readDecomp;

 double tot2MinTime = localPreProcess;
 double tot2AvgTime = localPreProcess;
 double tot2MaxTime = localPreProcess;

 double setUpDataTimeMin = times.setUpDataTime;
 double setUpDataTimeAvg = times.setUpDataTime;
 double setUpDataTimeMax = times.setUpDataTime;

 double makeSubDomainTimeMin = times.makeSubDomains;
 double makeSubDomainTimeAvg = times.makeSubDomains;
 double makeSubDomainTimeMax = times.makeSubDomains;

 double distributeBCTimeMin = times.distributeBCs;
 double distributeBCTimeAvg = times.distributeBCs;
 double distributeBCTimeMax = times.distributeBCs;

 double makeConnectivityTimeMin = times.makeConnectivity;
 double makeConnectivityTimeAvg = times.makeConnectivity;
 double makeConnectivityTimeMax = times.makeConnectivity;

 double makeInterfaceTimeMin = times.makeInterface;
 double makeInterfaceTimeAvg = times.makeInterface;
 double makeInterfaceTimeMax = times.makeInterface;

 double makeInternalTimeMin = times.makeInternalInfo;
 double makeInternalTimeAvg = times.makeInternalInfo;
 double makeInternalTimeMax = times.makeInternalInfo;

#ifdef DISTRIBUTED
 tot1MinTime = structCom->globalMin(tot1MinTime);
 tot1AvgTime = (structCom->globalSum(tot1AvgTime))/numCPUs;
 tot1MaxTime = structCom->globalMax(tot1MaxTime);

 tot2MinTime = structCom->globalMin(localPreProcess);
 tot2AvgTime = (structCom->globalSum(localPreProcess))/numCPUs;
 tot2MaxTime= structCom->globalMax(localPreProcess);

 setUpDataTimeMin = structCom->globalMin(times.setUpDataTime);
 setUpDataTimeAvg = (structCom->globalSum(times.setUpDataTime))/numCPUs;
 setUpDataTimeMax = structCom->globalMax(times.setUpDataTime);

 makeSubDomainTimeMin = structCom->globalMin(times.makeSubDomains);
 makeSubDomainTimeAvg = (structCom->globalSum(times.makeSubDomains))/numCPUs;
 makeSubDomainTimeMax = structCom->globalMax(times.makeSubDomains);
 distributeBCTimeMin = structCom->globalMin(times.distributeBCs);
 distributeBCTimeAvg = (structCom->globalSum(times.distributeBCs))/numCPUs;
 distributeBCTimeMax = structCom->globalMax(times.distributeBCs);

 makeConnectivityTimeMin = structCom->globalMin(times.makeConnectivity);
 makeConnectivityTimeAvg = (structCom->globalSum(times.makeConnectivity))/numCPUs;
 makeConnectivityTimeMax = structCom->globalMax(times.makeConnectivity);

 makeInterfaceTimeMin  = structCom->globalMin(times.makeInterface);
 makeInterfaceTimeAvg = (structCom->globalSum(times.makeInterface))/numCPUs;
 makeInterfaceTimeMax = structCom->globalMax(times.makeInterface);

 makeInternalTimeMin = structCom->globalMin(times.makeInternalInfo);
 makeInternalTimeAvg = (structCom->globalSum(times.makeInternalInfo))/numCPUs;
 makeInternalTimeMax = structCom->globalMax(times.makeInternalInfo);
#endif

 if(f) {

 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f," ... Detailed CPU Statistics (Seconds) ");
 filePrint(f,"\n***********************************************************"
           "********************\n");
 filePrint(f,"\n                                             minimum      average      maximum\n");

 filePrint(f,"\n1. Total Read Input Files             : %12.4f %12.4f %12.4f\n",
         tot1MinTime/1000.0, tot1AvgTime/1000.0, tot1MaxTime/1000.0);

 filePrint(f,"\n2. Total Preprocessing                : %12.4f %12.4f %12.4f\n",
         tot2MinTime/1000.0,tot2AvgTime/1000.0,tot2MaxTime/1000.0);

 filePrint(f,"         Process Input Data           : %12.4f %12.4f %12.4f\n", setUpDataTimeMin/1000.0, setUpDataTimeAvg/1000.0, setUpDataTimeMax/1000.0);

 filePrint(f,"         Make Subdomains              : %12.4f %12.4f %12.4f\n",
         makeSubDomainTimeMin/1000.0, makeSubDomainTimeAvg/1000.0, makeSubDomainTimeMax/1000.0);

 filePrint(f,"         Distribute BCs               : %12.4f %12.4f %12.4f\n",
         distributeBCTimeMin/1000.0, distributeBCTimeAvg/1000.0, distributeBCTimeMax/1000.0);
 filePrint(f,"         Make Connectivities          : %12.4f %12.4f %12.4f\n", makeConnectivityTimeMin/1000.0, makeConnectivityTimeAvg/1000.0, makeConnectivityTimeMax/1000.0);
 filePrint(f,"         Make Interface               : %12.4f %12.4f %12.4f\n", makeInterfaceTimeMin/1000.0, makeInterfaceTimeAvg/1000.0, makeInterfaceTimeMax/1000.0);
 filePrint(f,"         Make Internal Information    : %12.4f %12.4f %12.4f\n",
         makeInternalTimeMin/1000.0, makeInternalTimeAvg/1000.0, makeInternalTimeMax/1000.0);

 double tot3MinTime = timers.assemble+timers.constructMatrices;
 double tot3AvgTime = timers.assemble+timers.constructMatrices;
 double tot3MaxTime = timers.assemble+timers.constructMatrices;

 filePrint(f,"\n3. Total Matrix Processing            : %12.4f %12.4f %12.4f"
           "\n",tot3MinTime/1000.0,tot3AvgTime/1000.0,tot3MaxTime/1000.0);


 filePrint(f,"\n4. Total RHS Processing               : %12.4f %12.4f %12.4f"
           "\n\n",formRhsMin/1000.0,formRhsAvg/1000.0,formRhsMax/1000.0);

 // timers for 2nd level coarse problem

 double timeCoarseMin = coarseTime;
 double timeCoarseTot = coarseTime;
 double timeCoarseMax = coarseTime;


 filePrint(f,"5. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
         solutionTimeMin/1000.0, solutionTimeAvg/1000.0, solutionTimeMax/1000.0);

 if(sInfo.solvercntl->type == SolverSelection::Feti) {
   filePrint(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
             factorTimeMin/1000.0,factorTimeAvg/1000.0,factorTimeMax/1000.0);

   filePrint(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
             timeCoarseMin/1000.0,timeCoarseTot/1000.0,timeCoarseMax/1000.0);
   filePrint(f,"               1st Level Coarse Pb.   : %12.4f %12.4f %12.4f\n",
             coarse1Min/1000.0,coarse1Tot/1000.0,coarse1Max/1000.0);
   filePrint(f,"               2nd Level Coarse Pb.   : %12.4f %12.4f %12.4f\n\n",
             coarse2TimeMin/1000.0,coarse2TimeAvg/1000.0,coarse2TimeMax/1000.0);

   filePrint(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
             (parfac1Min+parfac2Min)/1000.0,(parfac1Tot+parfac2Tot)/1000.0,(parfac1Max+parfac2Max)/1000.0);
   filePrint(f,"               1st Level Factor       : %12.4f %12.4f %12.4f\n",
             parfac1Min/1000.0,parfac1Tot/1000.0,parfac1Max/1000.0);
   filePrint(f,"               2nd Level Factor       : %12.4f %12.4f %12.4f\n\n",
             parfac2Min/1000.0,parfac2Tot/1000.0,parfac2Max/1000.0);

   filePrint(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
             solveLoopMin/1000.0, solveLoopAvg/1000.0, solveLoopMax/1000.0);
   filePrint(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
             project1TimeMin/1000.0,project1TimeAvg/1000.0,project1TimeMax/1000.0);
   filePrint(f,"                  Seq. Forw/Back      : %12.4f %12.4f %12.4f\n",
             forBack1Min/1000.0,forBack1Tot/1000.0, forBack1Max/1000.0);
   filePrint(f,"               Project 2nd Level      : %12.4f %12.4f %12.4f\n",
             project2TimeMin/1000.0,project2TimeAvg/(1000.0),project2TimeMax/1000.0);
   filePrint(f,"                  Paral. Forw/Back    : %12.4f %12.4f %12.4f\n",
             forBack2Min/1000.0,forBack2Tot/(1000.0), forBack2Max/1000.0);

   filePrint(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
             precondMin/1000.0,precondAvg/1000.0, precondMax/1000.0);

   filePrint(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
             sAndJMin/1000.0, sAndJAvg/1000.0, sAndJMaximum/1000.0);
   filePrint(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n\n",
             reOrthoMin/1000.0,reOrthoAvg/1000.0,reOrthoMax/1000.0);
 }

 filePrint(f,"6. Write Output Files                 : %12.4f %12.4f %12.4f\n",
           outputTimeMin/1000.0,outputTimeAvg/1000.0,outputTimeMax/1000.0);

 double timeSimMin = tot1MinTime + tot2MinTime + tot3MinTime + formRhsMin
                   + solutionTimeMin + outputTimeMin;
 double timeSimAvg = tot1AvgTime + tot2AvgTime + tot3AvgTime + formRhsAvg
                   + solutionTimeAvg + outputTimeAvg;
 double timeSimMax = tot1MaxTime + tot2MaxTime + tot3MaxTime + formRhsMax
                   + solutionTimeMax + outputTimeMax;

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
 long preProcessMemoryMin = memoryPreProcess + times.memorySetUp;
 long preProcessMemoryMax = memoryPreProcess + times.memorySetUp;
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
 long tot3Min = subMatMemoryMin;
 long tot3Avg = subMatMemoryAvg;
 long tot3Max = subMatMemoryMax;

 filePrint(f,"\n2. Total Preprocessing                : %12.4f %12.4f %12.4f\n",
         preProcessMemoryMin*byteToMb,preProcessMemoryAvg*byteToMb,
         preProcessMemoryMax*byteToMb);

 filePrint(f,"\n3. Total Matrix Processing            : %12.4f %12.4f %12.4f"
           "\n", tot3Min*byteToMb,tot3Avg*byteToMb,tot3Max*byteToMb);


 long tot4Min = timers.memoryLoad;
 long tot4Avg = timers.memoryLoad;
 long tot4Max = timers.memoryLoad;

 filePrint(f,"\n4. Total RHS Processing               : %12.4f %12.4f %12.4f"
           "\n\n",tot4Min*byteToMb,tot4Avg*byteToMb,tot4Max*byteToMb);
 long tot5Min = solverMemoryMin;
 long tot5Avg = solverMemoryAvg;
 long tot5Max = solverMemoryMax;

 filePrint(f,"5. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
         tot5Min*byteToMb, tot5Avg*byteToMb, tot5Max*byteToMb);

 if(sInfo.solvercntl->type == SolverSelection::Feti) {
   filePrint(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
             factorMemoryMin*byteToMb, factorMemoryAvg*byteToMb,
             factorMemoryMax*byteToMb);

   filePrint(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
             coarseMemoryMin*byteToMb, coarseMemoryAvg*byteToMb, coarseMemoryMax*byteToMb);

   filePrint(f,"               1st Level Coarse Pb.   : %12.4f %12.4f %12.4f\n",
             totMemGtGMin*byteToMb, totMemGtGAvg*byteToMb, totMemGtGMax*byteToMb);

   filePrint(f,"               2nd Level Coarse Pb.   : %12.4f %12.4f %12.4f\n\n",
             0.0,0.0,0.0);
   filePrint(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
             0.0,0.0,0.0);
   filePrint(f,"               1st Level Factor       : %12.4f %12.4f %12.4f\n",
             0.0,0.0,0.0);
   filePrint(f,"               2nd Level Factor       : %12.4f %12.4f %12.4f"           "\n\n",0.0, 0.0,0.0);
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
   filePrint(f,"               Local Solve            : %12.4f %12.4f %12.4f\n", 0.0, 0.0, 0.0);
   filePrint(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n\n",
             orthoMemoryMin*byteToMb, orthoMemoryAvg*byteToMb,
             orthoMemoryMax*byteToMb);
 }

 filePrint(f,"6. Write Output Files                 : %12.4f %12.4f %12.4f\n",
         memoryOutputMin*byteToMb,memoryOutputAvg*numCPUs*byteToMb,memoryOutputMax*byteToMb);

 long totSimMin = readMemoryMin + tot3Min + preProcessMemoryMin + tot4Min
                     + tot5Min + memoryOutputMin;
 long totSimAvg = readMemoryAvg + tot3Avg + preProcessMemoryAvg + tot4Avg + tot5Avg + memoryOutputAvg;
 long totSimMax = readMemoryMax + tot3Max + preProcessMemoryMax + tot4Max
                     + tot5Max + memoryOutputMax;

 filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)        : %12.4f %12.4f %12.4f\n",
           totSimMin*byteToMb, totSimAvg*byteToMb, totSimMax*byteToMb);
 filePrint(f,"\n***********************************************************"
           "********************\n");

 if(f) fclose(f);

#endif
#endif
}

