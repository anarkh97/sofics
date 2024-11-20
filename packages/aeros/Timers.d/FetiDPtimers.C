#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <ctime>

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

double t0 = 0.0;
double t1 = 0.0;
double t2 = 0.0;
double t3 = 0.0;
double t4 = 0.0;
double t5 = 0.0;
double t6 = 0.0;
extern int iterTotal;

#ifndef SALINAS
static const char* problemType[] = {
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

static const char* renumMessage[] = {
"None",
"Sloan",
"RCM"
};

static const char* sparseRenumMessage[] = {
"Esmond",
"Metis"
};

//static const double oneMegaByte = (1024.0*1024.0); 
static const double byteToMb    = (1.0 / oneMegaByte);

// Message for printing selected scaling method
static const char* scalingMessage[] = {
"Scaling                           =      Stiffness\n",
"Scaling                           =    Topological\n"
}; 

// Message for printing selected mpc scaling method
static const char* mpc_scalingMessage[] = {
"MPC Scaling/Splitting             =      Stiffness\n",
"MPC Scaling/Splitting             =    Topological\n"
}; 


// Message for rigid body mode method
static const char* rbmMessage[] = {
"Algebraic RBM Method    tolerance =",
"Geometric RBM Method              ="
};

// Message for printing selected preconditioner method
static const char* precMessage[] = {
"Preconditioner                    =           None\n",
"Preconditioner                    =         Lumped\n",
"Preconditioner                    =      Dirichlet\n"
};

// Message for printing selected projector method
static const char* projectMessage[] = {
"Augmentation                      =           None\n",
"Augmentation                      =             GT\n",
"Augmentation                      =             GR\n",
"Augmentation                      =             G6\n",
"Augmentation                      =            EGT\n",
"Augmentation                      =            EGR\n",
"Augmentation                      =            EG6\n",
"Augmentation                      =           EAGT\n",
"Augmentation                      =           EAGR\n",
"Augmentation                      =           EAG6\n",
"Augmentation                      =         EW+EGT\n",
"Augmentation                      =         EW+EGR\n",
"Augmentation                      =         EW+EG6\n",
"Augmentation                      =             EW\n"};

static const char* outerSolverMessage[] = {
"Outer loop Solver Selected        =             CG\n",
"Outer loop Solver Selected        =          GMRES\n",
"Outer loop Solver Selected        =            GCR\n"
};

static const char* subSolverMessage[] = {
"Subdomain Solver Selected         =        Skyline\n",
"Subdomain Solver Selected         =         Sparse\n",
"Subdomain Solver Selected         =       BlockSky\n",
"Subdomain Solver Selected         =  SimplicialLLT\n",
"Subdomain Solver Selected         = SimplicialLDLT\n",
"Subdomain Solver Selected         =        Cholmod\n",
"Subdomain Solver Selected         =        Umfpack\n",
"Subdomain Solver Selected         =        SuperLU\n",
"Subdomain Solver Selected         =        Spooles\n",
"Subdomain Solver Selected         =          Mumps\n"
};

static const char* precSolverMessage[] = {
"Preconditioner Solver Selected    =        Skyline\n",
"Preconditioner Solver Selected    =         Sparse\n",
"Preconditioner Solver Selected    =       BlockSky\n",
"Preconditioner Solver Selected    =  SimplicialLLT\n",
"Preconditioner Solver Selected    = SimplicialLDLT\n",
"Preconditioner Solver Selected    =        Cholmod\n",
"Preconditioner Solver Selected    =        Umfpack\n",
"Preconditioner Solver Selected    =        SuperLU\n",
"Preconditioner Solver Selected    =        Spooles\n",
"Preconditioner Solver Selected    =          Mumps\n"	
};

static const char* krylovMessage[] = {
"Non-linear Krylov acceleration    =            off\n",
"Non-linear Krylov acceleration    =             on\n",
"Non-linear Krylov acceleration    = on (load step)\n"
};

static const char* mpcMessage[] = {
"MPC Strategy Selected             =    Dual method\n",
"MPC Strategy Selected             =  Primal method\n"
};

static const char* mpcPrecnoMessage[] = {
"MPC Preconditioner                =           None\n",
"MPC Preconditioner                =   Diagonal CCt\n",
"MPC Preconditioner                =     Global CCt\n",
"MPC Preconditioner    = Topological Block Diag CCt\n",
"MPC Preconditioner    =   Subdomain Block Diag CCt\n",
"MPC Preconditioner    =      Mortar Block Diag CCt\n"
};

// HB: add MPC cct_solver
static const char* mpcSolverMessage[] = {
"MPC Solver Selected               =        Skyline\n",
"MPC Solver Selected               =         Sparse\n",
"MPC Solver Selected               =       BlockSky\n",
"MPC Solver Selected               =  SimplicialLLT\n",
"MPC Solver Selected               = SimplicialLDLT\n",
"MPC Solver Selected               =        Cholmod\n",
"MPC Solver Selected               =        Umfpack\n",
"MPC Solver Selected               =        SuperLU\n",
"MPC Solver Selected               =        Spooles\n",
"MPC Solver Selected               =          Mumps\n"
};

static const char* gtgType[] = {
" Coarse Solver Selected            =       ",
" GtQG Solver Selected              =       ",
" GtQG Solver Selected              =       " 
};

static const char* gtgSolverMessage[] = {
" Skyline\n",
"  Sparse\n",
" Skyline\n", // Not implemented
" Skyline\n", // Not implemented
" Skyline\n", // Not implemented
" Skyline\n", // Not implemented
"BlockSky\n",
"BlockSky\n",
" Spooles\n", 
" Spooles (pivot)\n" 
" Mumps\n",
" Mumps (pivot)\n"
};

static const char* solverMessage[] = {
"1. FETI\n         One Level\n",
"1. FETI\n         Two Level\n",
"1. FETI\n         Two Level\n",
"1. FETI-DP\n",
"1. FETI-DPC\n",
"1. FETI-DPH\n"
};
#endif

// =========================================== 
// This function prints timers for FETI Solver 
// =========================================== 

void
StaticTimers::printFetiDPtimers(MatrixTimers matrixTimer, double solveTime,
                      SolverInfo& sInfo, Timings& timers, ControlInfo& cinfo,
                      Domain* domain)
{
#ifndef SALINAS
 double solutionTime = timers.solve + getFetiSolverTime + timeFreqSweep;

 double coarse1Max  = timers.coarse1;
 double coarse1Tot  = timers.coarse1;
 double coarse1Min  = timers.coarse1;

 double parfac1Max  = timers.pfactor;
 double parfac1Tot  = timers.pfactor;
 double parfac1Min  = timers.pfactor;

 double planningMax = timers.planning;
 double planningMin = timers.planning;
 double planningAvg = timers.planning;

 double nlprecMax = timers.nlPreCond;
 double nlprecMin = timers.nlPreCond;
 double nlprecAvg = timers.nlPreCond;

 double freqSweepMin = timeFreqSweep;
 double freqSweepMax = timeFreqSweep;
 double freqSweepAvg = timeFreqSweep;
  
 //double forBack1Max = timers.forBack;
 //double forBack1Tot = timers.forBack;
 //double forBack1Min = timers.forBack;

 double precondMax   = timers.preconditioner.getMax().time;
 double sAndJMaximum = timers.sAndJ - timers.project;
 int numSubdomains = timers.numSubdomains;

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

 int numCPUs  = threadManager->numThr();
#ifdef DISTRIBUTED
 numCPUs = structCom->numCPUs();
 numSubdomains = structCom->globalSum(numSubdomains);
 coarse1Max   = structCom->globalMax(timers.coarse1);
 coarse1Tot   = structCom->globalSum(timers.coarse1)/numCPUs;
 coarse1Min   = structCom->globalMin(timers.coarse1);
 parfac1Max   = structCom->globalMax(timers.pfactor);
 parfac1Tot   = structCom->globalSum(timers.pfactor)/numCPUs;
 parfac1Min   = structCom->globalMin(timers.pfactor);
 precondMax   = structCom->globalMax(timers.preconditioner.getMax().time);
 sAndJMaximum = structCom->globalMax(sAndJMaximum);
 solutionTime = structCom->globalMax(solutionTime);
 planningMax  = structCom->globalMax(timers.planning);
 planningMin  = structCom->globalMin(timers.planning);
 planningAvg  = structCom->globalSum(timers.planning)/numCPUs;
 nlprecMax     = structCom->globalMax(timers.nlPreCond);
 nlprecMin     = structCom->globalMin(timers.nlPreCond);
 nlprecAvg     = structCom->globalSum(timers.nlPreCond)/numCPUs;
 freqSweepMin  = structCom->globalMin(timeFreqSweep);
 freqSweepMax  = structCom->globalMax(timeFreqSweep);
 freqSweepAvg  = structCom->globalMax(timeFreqSweep)/numCPUs;

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

 // get timer statistics (min,max,avg, and total time)
 TimeStats constructStats = timers.consMatrix.getStats();
 TimeStats assembleStats  = timers.assembleMat.getStats();
 TimeStats factorStats    = timers.factorMat.getStats();
 TimeStats buildRhsStats  = timers.buildRhs.getStats();
 TimeStats sAndJStats     = timers.solveAndJump.getStats();
 TimeStats precStats      = timers.preconditioner.getStats();
 TimeStats orthoStats  = timers.orthogonalize.getStats();
 TimeStats project1Stats = timers.projection.getStats();
 TimeStats applyFetiPrecondStats = timers.applyFetiPrecond.getStats();

#ifdef DISTRIBUTED
 if(structCom->myID() == 0)
#endif
 openTimingFile(cinfo);
 time_t t;
 char *c;
 t = time(0);
 c = asctime(localtime(&t));
 filePrint(f," ... %s",c);

#ifdef USE_MPI
 if(structCom->myID() == 0)  {
#endif
#ifdef USE_MPI
 }
#endif
 long totMemSubMatrices = memorySubMatrices - memoryPrecond - memoryK;

 long localMemRead = matrixTimer.memoryParse + matrixTimer.memorySetUp;
 long totalMemRead = localMemRead;

 // Timer sub totals
 double subTotal[6];
 // max times
 subTotal[0] = (matrixTimer.readTime + matrixTimer.readDecomp);

 subTotal[1] = preProcess + matrixTimer.setUpDataTime + corotatorTime 
             + kelArrayTime + timeGeom - matrixTimer.readDecomp;

 subTotal[2] = (assembleStats.tot.time+constructStats.tot.time);

 subTotal[3] = buildRhsStats.tot.time;

 subTotal[4] = solutionTime - subTotal[2];

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

 int numDispCon = domain->nDispDirichlet();  // domain->nDirichlet() includes temperature loads
 int numTempCon = domain->nDirichlet() - domain->nDispDirichlet();
 filePrint(f,"\n***********************************************************"
           "********************\n");
 if(geoSource->isShifted()) {
   if(domain->solInfo().doFreqSweep)
     filePrint(f," ... Frequency Sweep Problem Information ... \n");
   else 
     filePrint(f," ... Frequency Response Problem Information ... \n");
 }
 else 
   filePrint(f," ... %s Problem Information ... \n", problemType[domain->probType()]);
 filePrint(f,"***********************************************************"
           "********************\n\n");

 filePrint(f,"1. Number of Nodes                         = %14d\n\n", domain->numNode()); 
 filePrint(f,"2. Number of Elements                      = %14d\n\n", domain->numElements());
 filePrint(f,"3. Number of Degrees of Freedom            = %14d\n", domain->numDofs());
 filePrint(f,"         Number of Constrained Disp. Dofs  = %14d\n", numDispCon);
 filePrint(f,"         Number of Constrained Temp. Dofs  = %14d\n\n", numTempCon);
 filePrint(f,"4. Number of Applied Loads                 = %14d\n\n", domain->nNeumann());
 filePrint(f,"5. Number of Multiple Point Constraints    = %14d\n\n", domain->getNumLMPC());
 filePrint(f,"6. Number of Output Files                  = %14d\n\n", geoSource->getNumOutInfo());
 filePrint(f,"7. Renumbering                             = %14s\n\n", renumMessage[sInfo.renum]);
 filePrint(f,"8. Sparse renumbering                      = %14s\n\n", sparseRenumMessage[sInfo.solvercntl->sparse_renum]);
 if(domain->solInfo().doFreqSweep) {
   filePrint(f,"9. Number of Frequencies                   = %14d\n\n", domain->numFrequencies);
   filePrint(f,"10. Number of RHS solves                   = %14d\n\n", domain->solInfo().getSweepParams()->nFreqSweepRHS);
 }

 filePrint(f,"***********************************************************"
            "********************\n");
 filePrint(f," ... Solver Information ... \n");
 filePrint(f,"***********************************************************"
           "********************\n\n");
 if(domain->numContactPairs > 0)
   filePrint(f,"%s",solverMessage[4]);
 else if(sInfo.getFetiInfo().dph_flag)
   filePrint(f,"%s",solverMessage[5]);
 else 
   filePrint(f,"%s",solverMessage[sInfo.getFetiInfo().version]);

 filePrint(f,"         %s", precMessage[sInfo.getFetiInfo().precno]);
 filePrint(f,"         %s", scalingMessage[sInfo.getFetiInfo().scaling-1]);
 int nAug = 0;
 if(sInfo.getFetiInfo().augment == FetiInfo::Gs) 
   switch( sInfo.getFetiInfo().rbmType) { 
     case FetiInfo::translation: nAug = 1; break;
     case FetiInfo::rotation:    nAug = 2; break;
     case FetiInfo::all:         nAug = 3; break;
     default:
      break;
   }

 if(sInfo.getFetiInfo().isEdgeAugmentationOn()) {
   switch(sInfo.getFetiInfo().rbmType) {
     case FetiInfo::translation:
       nAug = 4;
       break;
     case FetiInfo::rotation:
       nAug = 5;
       break;
     case FetiInfo::all:
       nAug = 6;
       break;
     case FetiInfo::averageTran:
       nAug = 7;
       break;
     case FetiInfo::averageRot:
       nAug = 8;
       break;
     case FetiInfo::averageAll:
       nAug = 9;
       break;
     case FetiInfo::None: 
       nAug = 7; // dph only
       break; 
     default:
       break; 
   }
 }

 if(sInfo.getFetiInfo().dph_flag && sInfo.getFetiInfo().isEdgeAugmentationOn()) {
   filePrint(f,"         %s", projectMessage[nAug+6]);
   filePrint(f,"         Number of wave directions         = %14d\n", 
             sInfo.getFetiInfo().numdir);
 }    
 else filePrint(f,"         %s", projectMessage[nAug]);
 if(sInfo.getFetiInfo().augment == FetiInfo::WeightedEdges)
   filePrint(f,"         Weighting of Augmentation Modes   =             on \n");
 else
   filePrint(f,"         Weighting of Augmentation Modes   =            off \n");

 filePrint(f,"         %s", outerSolverMessage[sInfo.getFetiInfo().outerloop]);

 filePrint(f,"         %s", subSolverMessage[sInfo.getFetiInfo().local_cntl->subtype]);
 filePrint(f,"         %s", precSolverMessage[sInfo.getFetiInfo().local_cntl->subtype]);
 if((sInfo.getFetiInfo().coarse_cntl->subtype==8)&&(sInfo.getFetiInfo().coarse_cntl->pivot))
   filePrint(f,"        %s%s",gtgType[sInfo.getFetiInfo().nonLocalQ],
                         gtgSolverMessage[sInfo.getFetiInfo().coarse_cntl->subtype+1]);
 else
   filePrint(f,"        %s%s",gtgType[sInfo.getFetiInfo().nonLocalQ],
                           gtgSolverMessage[sInfo.getFetiInfo().coarse_cntl->subtype]);

 if(sInfo.rbmflg == 0)
   filePrint(f,"         %s %14.3e\n",rbmMessage[sInfo.rbmflg],sInfo.solvercntl->trbm);
 else
   filePrint(f,"         %s%17e %e\n",rbmMessage[sInfo.rbmflg],
                                    sInfo.tolsvd, sInfo.solvercntl->trbm);

 filePrint(f,"         Maximum Number of Iterations      = %14d\n",
           sInfo.getFetiInfo().maxiter());
 filePrint(f,"         Maximum Size of Reortho. Vectors  = %14d\n",
           sInfo.getFetiInfo().maxorth());
 filePrint(f,"         Tolerance for Convergence         = %14.3e\n",
           sInfo.getFetiInfo().tolerance());
 filePrint(f,"         %s", krylovMessage[sInfo.getFetiInfo().nlPrecFlg]);

 int mpcflag = sInfo.getFetiInfo().mpcflag;
 if(mpcflag) {
   filePrint(f,"         %s", mpcMessage[mpcflag - 1]);
   if(mpcflag == 1) {
     filePrint(f,"         %s", mpc_scalingMessage[sInfo.getFetiInfo().scaling-1]);
     filePrint(f,"         %s", mpcPrecnoMessage[sInfo.getFetiInfo().mpc_precno]);
     filePrint(f,"         %s", mpcSolverMessage[sInfo.getFetiInfo().cct_cntl->subtype]);
   }
   if(mpcflag == 2) {
     filePrint(f,"         %s", mpc_scalingMessage[sInfo.getFetiInfo().scaling-1]);
   }
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
 filePrint(f,  "         Process Input Data            time: %14.5f s %14.3f Mb\n",
           matrixTimer.setUpDataTime/1000.0, memorySetUp*byteToMb);
 filePrint(f,  "         Make Subdomains               time: %14.5f s %14.3f Mb\n", 
           matrixTimer.makeSubDomains/1000.0, memorySubdomain*byteToMb);
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

 if(timeGeom != 0.0) {
   filePrint(f,  "         Make Geometric Node States    time: %14.5f s\n",
             timeGeom/1000.0);
   filePrint(f,  "         Make Element Corotators       time: %14.5f s\n",
             corotatorTime/1000.0);
   filePrint(f,  "         Make Stiffness & Mass Arrays  time: %14.5f s\n\n",
             kelArrayTime/1000.0);
 }

 filePrint(f,  "         Distribute BCs                time: %14.5f s %14.3f Mb\n", 
           matrixTimer.distributeBCs/1000.0, memoryDistBC*byteToMb);
 filePrint(f,  "         Make Connectivities           time: %14.5f s %14.3f Mb\n", 
           matrixTimer.makeConnectivity/1000.0, memoryConnect*byteToMb);
 filePrint(f,  "         Make Interface                time: %14.5f s %14.3f Mb\n", 
           matrixTimer.makeInterface/1000.0, memoryInterface*byteToMb);
 filePrint(f,  "         Make Internal Information     time: %14.5f s %14.3f Mb\n\n",        
           matrixTimer.makeInternalInfo/1000.0, memoryInternal*byteToMb);

 filePrint(f,"3. Total Subdomain Matrices Processing time: %14.5f s %14.3f Mb\n\n", 
           subTotal[2]/1000.0, totMemSubMatrices*byteToMb);

 filePrint(f,"4. Total Subdomain RHS Processing      time: %14.5f s %14.3f Mb\n\n", 
           subTotal[3]/1000.0, buildRhsStats.max.memory*byteToMb);
   
 }

 coarse1Max     -= timers.pfactor;
 coarse1Min     -= timers.pfactor;
 coarse1Tot     -= timers.pfactor;

 double coarseTime   = coarse1Max;

 long locMemCoarse = timers.memoryGtG + timers.memoryPCtFPC;
 long totMemCoarse = locMemCoarse;

 long localMemorySolve = timers.memoryProject1 + timers.memoryDV
                            + timers.preconditioner.getOverAll()->memory
                            + sAndJStats.tot.memory + timers.memoryOSet; 

 long totalMemorySolve = localMemorySolve;

 long localMemFactor = timers.memoryFactor + 8*timers.kMem.getOverAll()->memory;
 long totalMemFactor = localMemFactor;

 long localSolverMemory = localMemFactor + locMemCoarse + localMemorySolve;
 long totalSolverMemory = localSolverMemory;

 double factorTimeMax = timers.factor;
 double factorTimeMin = timers.factor;
 double factorTimeAvg = timers.factor;

 long totMemReortho = timers.memoryOSet;
 //long locMemUsed = memoryUsed();
 //long totMemUsed = locMemUsed;
 double locMemUsed = double(memoryUsed())*byteToMb;
 double totMemUsed = locMemUsed;
 long totMemFeti = timers.memoryFETI;

 long localMemoryCCt = timers.memoryBuildCCt;
 long totalMemoryCCt = localMemoryCCt;
 long maxMemoryCCt = localMemoryCCt;

#ifdef DISTRIBUTED
 // long minMemoryCCt = localMemoryCCt;
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
 
 totalMemoryCCt = timers.globalMemorySum(localMemoryCCt);
 maxMemoryCCt = timers.globalMemoryMax(localMemoryCCt);
 // minMemoryCCt = timers.globalMemoryMin(localMemoryCCt);
#endif
 factorTimeAvg /= numCPUs;


 if(f != 0) {

 filePrint(f,"5. Total Solver                        time: %14.5f s %14.3f Mb\n\n",
           subTotal[4]/1000.0, totalSolverMemory*byteToMb);
 filePrint(f,"         Factor Subdomain Matrices     time: %14.5f s %14.3f Mb\n\n",
           factorTimeMax/1000.0, totalMemFactor*byteToMb);
 filePrint(f,"         Total Building Coarse Pbs.    time: %14.5f s %14.3f Mb\n",
           coarseTime/1000.0, totMemCoarse*byteToMb );
 filePrint(f,"         Total Paral. Fac. Coarse Pbs. time: %14.5f s %14.3f Mb\n", 
           parfac1Max/1000.0, 0.0);
 if(sInfo.getFetiInfo().mpcflag) 
   filePrint(f,"         Build and Factor MPC CC^t     time: %14.5f s %14.3f Mb\n",
             timers.buildCCt/1000.0, timers.memoryBuildCCt*byteToMb);
 filePrint(f,"         Total Solve loop              time: %14.5f s %14.3f Mb\n",
           timers.solve/1000.0, totalMemorySolve*byteToMb);
 filePrint(f,"               Project 1st Level       time: %14.5f s %14.3f Mb\n",
           timers.project/1000.0, timers.memoryProject1*byteToMb);
 filePrint(f,"               Precondition            time: %14.5f s %14.3f Mb\n",
           precondMax/1000.0, memoryPrecond*byteToMb);
 if(sInfo.getFetiInfo().mpcflag)
   filePrint(f,"                  Total Solve MPC CC^t time: %14.5f s %14.3f Mb\n",
             timers.solveCCt/1000.0, totalMemoryCCt*byteToMb);
 if(sInfo.getFetiInfo().nlPrecFlg > 0)
   filePrint(f,"               Krylov Acceleration     time: %14.5f s %14.3f Mb\n",
             nlprecMax/1000.0, timers.memoryNlPreCond*byteToMb);
 filePrint(f,"               Local Solve             time: %14.5f s %14.3f Mb\n",
           sAndJMaximum/1000.0, sAndJStats.tot.memory*byteToMb);

 filePrint(f,"               Reorthogonalize         time: %14.5f s %14.3f Mb\n",
           orthoStats.max.time/1000.0, totMemReortho*byteToMb);
 if(domain->numContactPairs > 0)
   filePrint(f,"               Planning & Update       time: %14.5f s %14.3f Mb\n",
             planningMax/1000.0, timers.memoryPlanning*byteToMb);
 if(domain->solInfo().doFreqSweep)
   filePrint(f,"         Freq Sweep Series Expansion   time: %14.5f s %14.3f Mb\n",
             timeFreqSweep/1000.0, memoryFreqSweep*byteToMb); 

 filePrint(f,"\n6. Write Output Files                  time: %14.5f s %14.3f Mb\n",
           subTotal[5]/1000.0, memoryOutput*byteToMb);

 // Compute the total time spent on this simulation
/*
 cerr << "subTotal[0] = " << subTotal[0]/1000.0
      << ", totPreProMax = " << totPreProMax/1000.0
      << ", subTotal[2] = " << subTotal[2]/1000.0
      << ", subTotal[3] = " << subTotal[3]/1000.0
      << ", subTotal[4] = " << subTotal[4]/1000.0
      << ", subTotal[5] = " << subTotal[5]/1000.0
      << ", subTotal[6] = " << subTotal[6]/1000.0 << endl;
*/
 double total = subTotal[0] + totPreProMax + subTotal[2] + subTotal[3]
              + subTotal[4] + subTotal[5];

 long totalMemSimulation = totalMemRead + totalSolverMemory 
                              + buildRhsStats.tot.memory + memoryOutput
                              + totMemSubMatrices + totMemPreProcess;

 filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)         time: %14.5f s %14.3f Mb\n",total/1000.0, totalMemSimulation*byteToMb);

 // Output FETI solver information
 if(sInfo.solvercntl->type == SolverSelection::Feti) {
   filePrint(f,"\n***********************************************************"
             "********************\n");
   if(domain->numContactPairs > 0)
     filePrint(f," ... FETI-DPC Monitoring ... \n");
   else if(sInfo.getFetiInfo().dph_flag) 
     filePrint(f," ... FETI-DPH Monitoring ... \n");
   else
     filePrint(f," ... FETI-DP Monitoring ... \n");
   filePrint(f,"***********************************************************"
             "********************\n\n");

   filePrint(f,"1. Total Amount of Requested Memory        = %14.3f Mb\n\n",
//           (double)totMemUsed*byteToMb);
                     totMemUsed);

   filePrint(f,"2. FETI Solver Amount of Requested Memory  = %14.3f Mb\n\n",
             totMemFeti*byteToMb);

   // Check whether we converged or not
   if(timers.converged)
   filePrint(f,"3. Number of Iterations for Convergence    = %14d\n\n",
           timers.numIter);
   else {
   filePrint(f,"3. Stagnation Occured After a # of iter    = %14d\n\n",
           timers.numIter);
   }

   filePrint(f,"4. Relative Primal Error Reached           = %14.3e\n\n",
           timers.iterations[0].finalPrimal);

   filePrint(f,"5. Relative Dual Error Reached             = %14.3e\n\n",
           timers.iterations[0].finalDual);

   filePrint(f,"6. Total Size of 1st Level Coarse Problem  = %14d %14.3f Mb\n\n",
           timers.numCRNs,timers.memoryGtGsky*byteToMb);

   if(sInfo.getFetiInfo().augment == FetiInfo::Gs) {
     filePrint(f,"     Number Of Corner Dofs                 = %14d\n",
             timers.numCRNs - timers.numRBMs);
     filePrint(f,"     Number Of Subdomain Dofs              = %14d\n\n",
             timers.numRBMs);
   } else if(sInfo.getFetiInfo().augment == FetiInfo::Edges) {
     filePrint(f,"     Number Of Corner Dofs                 = %14d\n",
             timers.numCRNs - timers.numEdges);
     filePrint(f,"     Number Of Edge Dofs                   = %14d\n\n",
             timers.numEdges);
   }

   if(sInfo.getFetiInfo().local_cntl->subtype == 0)
     filePrint(f,"7. Total Memory Subdomain Skyline K        = %14.3f Mb\n\n",
               8.0*totMemSky*byteToMb);
   else if(sInfo.getFetiInfo().local_cntl->subtype == 1)
     filePrint(f,"7. Total Memory Subdomain Sparse K         = %14.3f Mb\n\n",
               8.0*totMemSparse*byteToMb);
   else
     filePrint(f,""); // if we have other subdomain solvers

   if(sInfo.getFetiInfo().mpcflag) {
     filePrint(f,"8. Total Memory used for CCt               = %14.3f Mb\n\n",
               totalMemoryCCt*byteToMb);
     filePrint(f,"9. Max CPU Memory used for CCt             = %14.3f Mb\n\n",
               maxMemoryCCt*byteToMb);
   }




   filePrint(f,"***********************************************************"
             "********************\n");

  }

 }
 

 long memPreProcessMin,memPreProcessTot,memPreProcessMax;
 computeMinAvgMax(memoryPreProcess,memPreProcessMin,memPreProcessTot,
                  memPreProcessMax);

 long memGtGMin, memGtGTot,memGtGMax;
 computeMinAvgMax(timers.memoryGtG,memGtGMin, memGtGTot,memGtGMax);

 long memOutMin, memOutTot, memOutMax;
 computeMinAvgMax(memoryOutput,memOutMin, memOutTot, memOutMax);

 //totMemSubMatrices = assembleStats.tot.memory + constructStats.tot.memory;

 //if (numCPUs == 1) return;

 double tot1MinTime = matrixTimer.readTime + matrixTimer.readDecomp;
 double tot1AvgTime = matrixTimer.readTime + matrixTimer.readDecomp;
 double tot1MaxTime = matrixTimer.readTime + matrixTimer.readDecomp;

 double tot2MinTime = subTotal[1];
 double tot2AvgTime = subTotal[1];
 double tot2MaxTime = subTotal[1];
 
 double setUpDataTimeMin = matrixTimer.setUpDataTime;
 double setUpDataTimeAvg = matrixTimer.setUpDataTime;
 double setUpDataTimeMax = matrixTimer.setUpDataTime;

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
         matrixTimer.makeSubDomains/1000.0,
	 matrixTimer.makeSubDomains/1000.0,
	 matrixTimer.makeSubDomains/1000.0);

 filePrint(f,"         Distribute BCs               : %12.4f %12.4f %12.4f\n", 
         distributeBCTimeMin/1000.0, distributeBCTimeAvg/1000.0, distributeBCTimeMax/1000.0);
 filePrint(f,"         Make Connectivities          : %12.4f %12.4f %12.4f\n", 
         makeConnectivityTimeMin/1000.0, makeConnectivityTimeAvg/1000.0, makeConnectivityTimeMax/1000.0);
 filePrint(f,"         Make Interface               : %12.4f %12.4f %12.4f\n", 
         makeInterfaceTimeMin/1000.0, makeInterfaceTimeAvg/1000.0, makeInterfaceTimeMax/1000.0);
 filePrint(f,"         Make Internal Information    : %12.4f %12.4f %12.4f\n",
         makeInternalTimeMin/1000.0, makeInternalTimeAvg/1000.0, makeInternalTimeMax/1000.0);

 double tot3MinTime = constructStats.min.time + assembleStats.min.time;
 double tot3AvgTime = constructStats.avg.time + assembleStats.avg.time;
 double tot3MaxTime = constructStats.max.time + assembleStats.max.time; 

 filePrint(f,"\n3. Total Subdomain Matrices Processing: %12.4f %12.4f %12.4f"
           "\n",tot3MinTime/1000.0,tot3AvgTime/1000.0,tot3MaxTime/1000.0);
 
 double tot4MinTime = buildRhsStats.min.time;
 double tot4AvgTime = buildRhsStats.avg.time;
 double tot4MaxTime = buildRhsStats.max.time;

 filePrint(f,"\n4. Total Subdomain RHS Processing     : %12.4f %12.4f %12.4f"
           "\n\n",tot4MinTime/1000.0,tot4AvgTime/1000.0,tot4MaxTime/1000.0);

 double timeCoarseMin = coarse1Min;
 double timeCoarseTot = coarse1Tot;
 double timeCoarseMax = coarse1Max;

 double timeParFacMin = parfac1Min;
 double timeParFacTot = parfac1Tot;
 double timeParFacMax = parfac1Max;

 double com1 = precondMax   - precStats.max.time;
 double com2 = sAndJMaximum - sAndJStats.max.time;
 
 double proj1Seq = timers.project - project1Stats.max.time;

 double proj1Min = project1Stats.min.time + proj1Seq;
 double proj1Avg = project1Stats.avg.time + proj1Seq;
 double proj1Max = project1Stats.max.time + proj1Seq;

 double solveMin =  proj1Min + 
                 +  precStats.min.time + com1 + sAndJStats.min.time + com2 
                 +  orthoStats.min.time;
 double solveAvg =  proj1Avg + 
                 +  precStats.avg.time + com1 + sAndJStats.avg.time + com2 
                 +  orthoStats.avg.time;
 double solveMax =  proj1Max + 
                 +  precStats.max.time + com1 + sAndJStats.max.time + com2 
                 +  orthoStats.max.time;

 double missingTime = timers.solve - solveMax;

 solveMin += missingTime;
 solveAvg += missingTime;
 solveMax += missingTime;

 double tot5MinTime   = timeCoarseMin + factorStats.min.time + timeParFacMin
                      + solveMin + freqSweepMin; 

 double tot5AvgTime   = factorStats.avg.time + timeCoarseTot
                      + timeParFacTot  + solveAvg + freqSweepAvg;

 double tot5MaxTime   = timeCoarseMax + factorStats.max.time + timeParFacMax
                      + solveMax + freqSweepMax;

 filePrint(f,"5. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
           tot5MinTime/1000.0, tot5AvgTime/1000.0, tot5MaxTime/1000.0);

 filePrint(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
           factorTimeMin/1000.0, factorTimeAvg/1000.0, factorTimeMax/1000.0);

 filePrint(f,"         Total Building Coarse Pbs.   : %12.4f %12.4f %12.4f\n",
           timeCoarseMin/1000.0,timeCoarseTot/1000.0,timeCoarseMax/1000.0);

 filePrint(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
           timeParFacMin/1000.0,timeParFacTot/1000.0,timeParFacMax/1000.0);

 filePrint(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
           solveMin/1000.0,solveAvg/1000.0,solveMax/1000.0);

 filePrint(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
           proj1Min/1000.0,proj1Avg/1000.0,proj1Max/1000.0);

 filePrint(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
         (precStats.min.time+com1)/1000.0,(precStats.avg.time+com1)/1000.0,
          precondMax/1000.0);

 if(sInfo.getFetiInfo().nlPrecFlg > 0)
   filePrint(f,"               Krylov Acceleration    : %12.4f %12.4f %12.4f\n",
             nlprecMin/1000.0, nlprecAvg/1000.0, nlprecMax/1000.0);

 filePrint(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
           (sAndJStats.min.time+com2)/1000.0,(sAndJStats.avg.time+com2)/1000.0,
           sAndJMaximum/1000.0);

 filePrint(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n",
           orthoStats.min.time/1000.0,orthoStats.avg.time/1000.0,orthoStats.max.time/1000.0);

 if(domain->numContactPairs > 0)
   filePrint(f,"               Planning & Update      : %12.4f %12.4f %12.4f\n",
             planningMin/1000.0, planningAvg/1000.0, planningMax/1000.0);

 if(domain->solInfo().doFreqSweep)
   filePrint(f,"               Frequency Sweep        : %12.4f %12.4f %12.4f\n",
             freqSweepMin/1000.0, freqSweepAvg/1000.0, freqSweepMax/1000.0); 

 double tot6MinTime = output;
 double tot6AvgTime = output;
 double tot6MaxTime = output;

 filePrint(f,"\n6. Write Output Files                 : %12.4f %12.4f %12.4f\n",
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

 long locSubMatMem = timers.memorySubMatrices - 8*timers.kMem.getOverAll()->memory - precOverall->memory;
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
 long planMemoryMax = timers.memoryPlanning;
 long planMemoryMin = timers.memoryPlanning;
 long planMemoryAvg = timers.memoryPlanning;
 long nlprecMemoryMax = timers.memoryNlPreCond;
 long nlprecMemoryMin = timers.memoryNlPreCond;
 long nlprecMemoryAvg = timers.memoryNlPreCond;

 // get gloabl sums mins and maxes
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
 planMemoryMax = timers.globalMemoryMax(timers.memoryPlanning);
 planMemoryMin = timers.globalMemoryMin(timers.memoryPlanning);
 planMemoryAvg = timers.globalMemorySum(timers.memoryPlanning);
 nlprecMemoryMax = timers.globalMemoryMax(timers.memoryNlPreCond);
 nlprecMemoryMin = timers.globalMemoryMin(timers.memoryNlPreCond);
 nlprecMemoryAvg = timers.globalMemorySum(timers.memoryNlPreCond);

 long readMemoryAvg = totalMemRead/numMPI;
 long preProcessMemoryAvg = totMemPreProcess/numMPI;
 long subMatMemoryAvg = totMemSubMatrices/numMPI;
 long solverMemoryAvg = totalSolverMemory/numMPI;
 long factorMemoryAvg = totalMemFactor/numMPI;
 long orthoMemoryAvg = totMemReortho / numMPI;
 long coarseMemoryAvg = totMemCoarse / numMPI;
 precMemoryAvg /= numMPI;
 planMemoryAvg /= numMPI;
 nlprecMemoryAvg /= numMPI;


 filePrint(f,"\n1. Total Read Input Files             : %12.4f %12.4f %12.4f\n",
         readMemoryMin*byteToMb, readMemoryAvg*byteToMb, readMemoryMax*byteToMb);
/*
 long tot3Min = (constructStats.min.memory + assembleStats.min.memory);
 long tot3Avg = (constructStats.avg.memory + assembleStats.avg.memory);
 long tot3Max = (constructStats.max.memory + assembleStats.max.memory);
*/ 
 long tot3Min = subMatMemoryMin;
 long tot3Avg = subMatMemoryAvg;
 long tot3Max = subMatMemoryMax;

 filePrint(f,"\n2. Total Preprocessing                : %12.4f %12.4f %12.4f\n", 
         preProcessMemoryMin*byteToMb, preProcessMemoryAvg*byteToMb,
         preProcessMemoryMax*byteToMb);
	 
 filePrint(f,"\n3. Total Subdomain Matrices Processing: %12.4f %12.4f %12.4f"
           "\n", tot3Min*byteToMb, tot3Avg*byteToMb, tot3Max*byteToMb);
 
 long tot4Min = buildRhsStats.min.memory;
 long tot4Avg = buildRhsStats.avg.memory;
 long tot4Max = buildRhsStats.max.memory;

 filePrint(f,"\n4. Total Subdomain RHS Processing     : %12.4f %12.4f %12.4f"
           "\n\n",tot4Min*byteToMb,tot4Avg*byteToMb,tot4Max*byteToMb);
/*
 long tot5Min = factorMemoryMin + coarseMemoryMin + orthoMemoryMin;
 long tot5Avg = factorMemoryAvg + coarseMemoryAvg + orthoMemoryAvg;
 long tot5Max = factorMemoryMax + coarseMemoryMax + orthoMemoryMax;
*/
 long tot5Min = solverMemoryMin;
 long tot5Avg = solverMemoryAvg;
 long tot5Max = solverMemoryMax;
 filePrint(f,"5. Total Solver                       : %12.4f %12.4f %12.4f\n\n",
         tot5Min*byteToMb, tot5Avg*byteToMb, tot5Max*byteToMb);

 filePrint(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f\n\n",
           factorMemoryMin*byteToMb, factorMemoryAvg*byteToMb,
           factorMemoryMax*byteToMb);
 
 filePrint(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
           coarseMemoryMin*byteToMb, coarseMemoryAvg*byteToMb,
           coarseMemoryMax*byteToMb);

 filePrint(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f\n",
           0.0,0.0,0.0);
 filePrint(f,"         Total Solve loop             : %12.4f %12.4f %12.4f\n",
           0.0,0.0,0.0);
 filePrint(f,"               Project 1st Level      : %12.4f %12.4f %12.4f\n",
           0.0,0.0,0.0);

 filePrint(f,"               Precondition           : %12.4f %12.4f %12.4f\n",
           precMemoryMin*byteToMb, precMemoryAvg*byteToMb, precMemoryMax*byteToMb);

 if(sInfo.getFetiInfo().nlPrecFlg > 0)
   filePrint(f,"               Krylov Acceleration    : %12.4f %12.4f %12.4f\n",
             nlprecMemoryMin*byteToMb, nlprecMemoryAvg*byteToMb, nlprecMemoryMax*byteToMb);

 filePrint(f,"               Local Solve            : %12.4f %12.4f %12.4f\n",
           sAndJStats.min.memory*byteToMb, sAndJStats.avg.memory*byteToMb, sAndJStats.max.memory*byteToMb);

 filePrint(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f\n",
           orthoMemoryMin*byteToMb, orthoMemoryAvg*byteToMb, orthoMemoryMax*byteToMb);

 if(domain->numContactPairs > 0)
   filePrint(f,"               Planning & Update      : %12.4f %12.4f %12.4f\n",
             planMemoryMin*byteToMb, planMemoryAvg*byteToMb, planMemoryMax*byteToMb);

 long tot6Min = memOutMin;
 long tot6Avg = memOutTot/numCPUs;
 long tot6Max = memOutMax;

 filePrint(f,"\n6. Write Output Files                 : %12.4f %12.4f %12.4f\n",
         tot6Min*byteToMb,tot6Avg*byteToMb,tot6Max*byteToMb);

 long totSimMin = readMemoryMin + tot3Min + memPreProcessMin + tot4Min 
                     + tot5Min + tot6Min;
 long totSimAvg = readMemoryAvg + tot3Avg + memPreProcessTot/numMPI + tot4Avg + tot5Avg + tot6Avg;
 long totSimMax = readMemoryMax + tot3Max + memPreProcessMax + tot4Max + tot5Max + tot6Max;
 
 filePrint(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)        : %12.4f %12.4f %12.4f\n",
           totSimMin*byteToMb, totSimAvg*byteToMb, totSimMax*byteToMb);
 filePrint(f,"\n***********************************************************"
           "********************\n");
 if(f) fclose(f);

#endif
#endif
}


void
computeMinAvgMax(double value, long &minValue, long &totValue,
                               long &maxValue)
{
 double mem1 = (double) value;
#ifdef DISTRIBUTED
 mem1 = structCom->globalMin(mem1);
#endif
minValue = (long) mem1;

 mem1 = (double) value;
#ifdef DISTRIBUTED
 mem1 = structCom->globalSum(mem1);
#endif
 totValue = (long) mem1;

 mem1 = (double) value;
#ifdef DISTRIBUTED
 mem1 = structCom->globalMax(mem1);
#endif
 maxValue = (long) mem1;

}
