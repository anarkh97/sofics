/** \file main.C
 * \brief The driver file.
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <unistd.h>
#include <memory>
#include <Utils.d/dbg_alloca.h>

#include <GNU-getopt.d/getopt.h>

#include <Driver.d/Domain.h>
#include <Driver.d/StaticProbType.h>
#include <Driver.d/DynamProbType.h>
#include <Driver.d/EigenProbType.h>
#include <Driver.d/CondProbType.h>
#include <Driver.d/NLStaticProbType.h>
#include <Driver.d/NLDynamProbType.h>
#include <Driver.d/TempProbType.h>
#include <Solvers.d/Solver.h>
#include <Problems.d/StaticDescr.h>
#include <Problems.d/DynamDescr.h>
#include <Problems.d/EigenDescr.h>
#include <Problems.d/CondDescr.h>
#include <Problems.d/NonLinStatic.h>
#include <Problems.d/NonLinDynam.h>
#include <Problems.d/TempDescr.h>
#include <Problems.d/ModalDescr.h>
#include <Problems.d/DEMProblem.h>
#include <Paral.d/MDStatic.h>
#include <Paral.d/MDInpcStatic.h>
#include <Paral.d/MDDynam.h>
#include <Paral.d/MDEigen.h>
#include <Paral.d/MDNLStatic.h>
#include <Paral.d/MDNLDynam.h>
#include <Paral.d/MDTemp.h>
#include <Paral.d/MDModal.h>
#include <Driver.d/GeoSource.h>
#include <Dec.d/dec.h>
#include <Parser.d/DecInit.h>
#include <Sfem.d/Sfem.h>
#include <Rom.d/SnapshotNonLinDynamic.h>
#include <Rom.d/DistrSnapshotNonLinDynamic.h>
#include <Rom.d/DistrPodProjectionNonLinDynamic.h>
#include <Rom.d/DistrLumpedPodProjectionNonLinDynamic.h>
#include <Rom.d/PodProjectionNonLinDynamic.h>
#include <Rom.d/LumpedPodProjectionNonLinDynamic.h>
#include <Rom.d/DEIMPodProjectionNonLinDynamic.h>
#include <Rom.d/UDEIMPodProjectionNonLinDynamic.h>
#include <Rom.d/CheckNonLinDynamic.h>
#include <Rom.d/PodProjectionSolver.h>
#include <Rom.d/EiGalerkinProjectionSolver.h>
#include <Rom.d/DriverInterface.h>
#include <Rom.d/DistrExplicitSnapshotNonLinDynamic.h>
#include <Rom.d/DistrExplicitPodProjectionNonLinDynamic.h>
#include <Rom.d/DistrExplicitPodProjectionNonLinDynamicBase.h>
#include <Rom.d/DistrExplicitLumpedPodProjectionNonLinDynamic.h>
#include <Rom.d/DistrExplicitDEIMPodProjectionNonLinDynamic.h>
#include <Rom.d/DistrExplicitUDEIMPodProjectionNonLinDynamic.h>
#ifdef USE_PITA
#include <Pita.d/Old.d/PitaNonLinDynam.h>
#include <Pita.d/Old.d/NLDistrTimeDecompSolver.h>
#include <Pita.d/PitaNonLinDynam.h>
#include <Pita.d/NlDriver.h>
#include <Pita.d/LinearDriver.h>
#endif
#include <Comm.d/Communicator.h>
#include <Solvers.d/SolverFactory.h>
#include <Driver.d/SubDomainFactory.h>

// .... for different problems and hardware
void writeOptionsToScreen();

#ifdef _OPENMP
#include <omp.h>
#endif

DecInit * decInit=0;

#ifdef TFLOP
extern int optind;
extern "C" int getopt (
			 int argc,
			 const char *argv[],
			 const char *optstring );
#endif

// .... global and static member variable initialization
std::map<int,SolverCntl> SolverInfo::solvercntls = std::map<int,SolverCntl>();
SolverCntl default_cntl;
std::unique_ptr<GenSolverFactory<double> >   solverFactory(new GenSolverFactory<double>());
std::unique_ptr<GenSolverFactory<DComplex> > solverFactoryC(new GenSolverFactory<DComplex>());;

Domain *domain = new Domain();
std::unique_ptr<ElementFactory> elemFact(new ElementFactory());
std::unique_ptr<GenSubDomainFactory<double> >   subDomainFactory(new GenSubDomainFactory<double>());
std::unique_ptr<GenSubDomainFactory<DComplex> > subDomainFactoryC(new GenSubDomainFactory<DComplex>());;

SolverInfo &solInfo = domain->solInfo();

Sfem *sfem = new Sfem();

Connectivity *procMpcToMpc;

long totMemSky       = 0;
long totMemSparse    = 0;
long totMemSpooles   = 0;
long totMemMumps     = 0;

ThreadManager *threadManager = 0;

extern int yyparse(void);
extern FILE* yyin;

bool callDec=false;
bool exitAfterDec=false;
bool estFlag=false;
bool weightOutFlag=false;
bool nosa=false;
bool useFull=false;
bool trivialFlag=false;
bool randomShuffle=false;
bool fsglFlag=false;
bool allowMechanisms=false;
bool useScotch=false;

int verboseFlag = 0;
int contactPrintFlag = 0;
int salinasFlag = 0;
int totalNewtonIter = 0;
int iterTotal = 0;
int debugFlag = 0;
int quietFlag = 0;

SysCom *syscom = 0;
Communicator *structCom = 0;
Communicator *heatStructCom;
Communicator *fluidCom;

extern const char* problemTypeMessage[];
extern const char* sensitivityTypeMessage[];
extern const char* solverTypeMessage[];

std::string clusterData_ = "INPUT.msh";
std::string subdomains_ = "INPUT.sub";
std::string decomposition_ = "INPUT.dec";
std::string connectivity_ = "INPUT.con";

extern const char *THE_VERSION;

extern ModeData modeData;
extern ModeData modeDataIDis;
extern ModeData modeDataIVel;
extern ModeData modeDataMode;

// ... main program

#ifdef CREATE_DSO
//extern "C"
int entrypoint(int argc, char** argv)
#else
int main(int argc, char** argv)
#endif
{
#ifdef __GNUC__
//  std::set_new_handler(&print_trace_handler);
#endif
 double initTime = getTime();
 double totalMemoryUsed = 0.0;

 int c;
 extern char *optarg;

 /**** DEFAULT ELEMENT WEIGHT VALUES  ***/
 /**** IF DIFFERENT FROM 1.0 add the weight of your object below
	   weightList[objectNumber] = makePair(my_weight, my_true_weight)  ***/
 weightList[2] = 2.0;   // FourNodeQuad
 weightList[3] = 2.0;   // Therm3DQuad
 weightList[4] = 2.0;   // Triangle3
 weightList[8] = 3.0;   // ThreeNodeShell
 weightList[10] = 2.0;  // ThermQuadGal
 weightList[11] = 1.0;  // TorSpring
 weightList[12] = 1.0;  // LinSpring
 weightList[13] = 2.0;  // thermoleastic quad
 weightList[15] = 3.0;  // FelippaShell
 weightList[1515] = 4.0;  // FelippaShellX2
 weightList[16] = 4.0;  // FourNodeShell
 weightList[17] = 3.0;  // EightNodeBrick
 weightList[19] = 3.0;  // Membrane
 weightList[20] = 3.0;  // Compo3NodeShell
 weightList[2020] = 4.0;  // Compo4NodeShell
 weightList[23] = 3.0;  // Tetrahedral
 weightList[24] = 3.0;	// Pentahedral
 weightList[25] = 4.0;  // TenNodeTetrahedral
 weightList[32] = 3.0;  // HelmQuad8Gal
 weightList[38] = 3.0;  // HelmTri6Gal
 weightList[40] = 2.0;  // TetraHelmGal
 weightList[41] = 2.0;  // TetraHelmGLS
 weightList[42] = 3.0;  // Tetra10HelmGal
 weightList[44] = 3.0;  // HelmBrickGLS
 weightList[45] = 3.0;  // HelmBrick
 weightList[46] = 2.0;  // Therm3NoShell
 weightList[4646] = 3.0;  // Therm4NoShell
 weightList[48] = 2.0;  // QuadConvec
 weightList[49] = 2.0;  // TriangleConvec
 weightList[50] = 3.0;  // TetraTherm
 weightList[51] = 4.0;  // ThermBrick
 weightList[53] = 2.0;
 weightList[57] = 2.0;
 weightList[58] = 2.0;
 weightList[70] = 3.0;  // RigidEightNodeBrick
 weightList[73] = 3.0;  // RigidThreeNodeShell
 weightList[76] = 4.0;  // RigidFourNodeShell
 weightList[72] = 4.0;  // Brick20
 weightList[81] = 2.0;
 weightList[82] = 4.0;
 weightList[83] = 4.0;
 weightList[84] = 2.0;
 weightList[85] = 3.0;
 weightList[86] = 4.0;
 weightList[87] = 4.0;
 weightList[88] = 4.0;  // FourNodeShell
 weightList[90] = 3.0;  // HelmPenta
 weightList[91] = 6.0;  // 3d 32 node serendipity brick
 weightList[92] = 5.0;  // 3d 26 node serendipity wedge
 weightList[93] = 5.0;  // 3d 32 node serendipity brick
 weightList[94] = 4.0;  // 3d 26 node serendipity wedge
 weightList[97] = 4.0;  // 3d 15 node serendipity wedge
 weightList[201] = 3.0; // nonlinear translational spring
 weightList[202] = 3.0; // nonlinear torsional spring
 weightList[301] = 1.0; // 2d 4-node sloshing (fluid) quadrilateral
 weightList[302] = 1.0; // 2d 2-node free-surface (fluid)
 weightList[311] = 3.0; // 3d 4-node sloshing (fluid) tetrahedron
 weightList[312] = 2.0; // 3d 3-node free-surface (fluid) triangle
 weightList[321] = 2.0; // 2d 4-node hydroelastic vibration (fluid) quadrilateral
 weightList[331] = 3.0; // 3d 4-node hydroelastic vibration (fluid) tetrahedral
 weightList[1101] = 2.0;
 weightList[1102] = 4.0;
 weightList[1103] = 8.0;
 weightList[1104] = 2.0;
 weightList[1111] = 2.0;
 weightList[1121] = 2.0;
 weightList[1122] = 4.0;
 weightList[1123] = 8.0;
 weightList[1131] = 2.0;
 weightList[1151] = 4.0;
 weightList[1152] = 8.0;
 weightList[1153] = 12.0;
 weightList[1161] = 4.0;
 weightList[1162] = 8.0;
 weightList[1171] = 4.0;
 weightList[1172] = 8.0;
 weightList[1173] = 12.0;
 weightList[1200] = 2.0;
 weightList[1201] = 8.0;
 weightList[1220] = 2.0;
 weightList[1250] = 3.0;
 weightList[1251] = 15.0;
 weightList[1252] = 28.0;

#if defined(USE_MPI)
 SysCom theCom(argc,argv);
 syscom = &theCom;
 structCom = syscom;
#endif

 // Default number of threads equal to one
 int numThreads =  1;
 int numSubdomains = 1;
 int numProcessors = 1;
 int topFlag    = -1;
 int numClusters = 1;
 int numCpus = 1;

 bool callSower = false;
 bool exitAfterSower = false;
 // Process command line arguments
 //
 // -n = number of threads lok
 // -d = decomposition file name lok
 // -v = verbose output lok
 // -c = verbose output plus contact status change info
 // -t = output topdomdec file and exit lok
 // -T = output topdomdec file with gaps renumbered sequentially lok
 // -m = output topdomdec file with number of element sets equal to
 //      the number of materials in the input file.lok
 // -r = output topdomdec file for axisymmetric mesh and exit
 // -p = output primal residual at each FETI iteration to a file
 // -q = suppress certain warnings

 // Process command line arguments for DEC
 //
 // -n <number of threads> OK
 // -d <decomposition file name> OK
 // -s <number of subdomains> OK
 // -t = output topdomdec file and exit OK
 // -T = output topdomdec file with gaps renumbered sequentially OK
 // -m = output topdomdec file with number of element sets equal to OK
 //      the number of materials in the input file.
 // -e = output memory estimate in a file with a .memory extension OK
 // -w = output weights with a .weight extension

 // -e = output memory estimate in a file with a .memory extension OK
 // -w = output weights with a .weight extension

 if(argc == 1) {
   writeOptionsToScreen();
 }

 // getopt_long
 int option_index = 0; // will hold index for long options
 static struct option long_options[] = {
   {"with-dec", 0, nullptr, 1000},
   {"dec", 0, nullptr, 1000},
   {"exit", 0, nullptr, 1002},
   {"deter", 0, nullptr, 1005},
   {"trivial", 0, nullptr, 1007},
   {"allow-mechanisms", 0, nullptr, 1008},
   {"use-scotch", 0, 0, 1009},
   {"use-weight-from", 1, nullptr, 1004},
   {"threads-number", 1, nullptr, 'n'},
   {"decomposition-filename", 1, nullptr, 'd'},
   {"nsub", 1, nullptr, 1001},
   {"subdomains-number", 1, nullptr, 1001},
   {"processors-number", 1, nullptr, 1003},
   {"output-topdomdec", 0, nullptr, 't'},
   {"output-topdomdec-asymetric-mesh", 1, nullptr, 'r'},
   {"output-primal", 1, nullptr, 'p'},
   {"contact-status-verbose", 1, nullptr, 'c'},
   {"screen-output-interval-in-time-loop", 1, nullptr, 's'},
   {"output-topdomdec-gaps-renumbered-sequentially", 0, nullptr, 'T'},
   {"output-topdomdec-element-equals-materials-number", 0, nullptr, 'm'},
   {"output-topdomdec-element-equals-materials-numberc-gaps-renumbered-sequentially", 0, nullptr, 'M'},
   {"output-match", 0, nullptr, 'P'},
   {"output-memory-estimate", 0, nullptr, 'e'},
   {"mem", 0, nullptr, 'e'},
   {"output-weights", 0, nullptr, 'w'},
   {"load", 0, nullptr, 'w'},
   {"verbose", 1, nullptr, 'v'},
   {"with-sower", 0, nullptr, 1010},
   {"sower", 0, nullptr, 1010},
   {"prefix", 1, nullptr, 1011},
   {"nclus", 1, nullptr, 1012},
   {"debug", 0, nullptr, 1006},
   {"quiet", 0, nullptr, 'q'},
   {"ncpu", 1, 0, 1013},
   {0, 0, nullptr, 0}
 };
 // end getopt_long

 filePrint(stderr,"\n --------- R U N  PARAMETERS ----------\n");
 FILE * weightFile;
 while ((c = getopt_long(argc, argv, "n:d:p:v:c:DVtTPmMr:Pfs:q", long_options, &option_index)) != -1)
	  switch (c) {
	  case 1000 :  // call dec from FEM
	callDec = true;
	break;
	  case 1001 :
	numSubdomains = atoi(optarg);
		if(numSubdomains <= 0) numSubdomains = 1;
	break;
	  case 1002 :  // exit after dec
	exitAfterDec = true;
	geoSource->setExitAfterDec(true);

	exitAfterSower = true; //TG just in case
	break;
	  case 1003 :  // number of processors for dec
	numProcessors = atoi(optarg);
	break;
	  case 1004 :
	weightFile = fopen(optarg, "r");
	double w;
	int k;
	if(weightFile)
	  {
		int i=1;
		char c;
		fpos_t position;
		while(!feof(weightFile))
		  {
		int res=fscanf(weightFile,"%d",&k);
		if(res == 0 || res == EOF)
		  {
			std::cerr << "*** WEIGHT FILE CORRUPTED AT LINE " << i << " bad object ID." << std::endl;
			exit(1);
		  }
		fgetpos(weightFile, &position);
		while((c = fgetc(weightFile))==' ')
		  ;
		if(c=='\n')
		  {
			std::cerr << "*** WEIGHT FILE CORRUPTED AT LINE " << i << " : no weight specified !" << std::endl;
			exit(1);
		  }
		fsetpos(weightFile, &position);
		res=fscanf(weightFile,"%lf",&w);
		if(res == 0 || res == EOF)
		   {
			 std::cerr << "*** WEIGHT FILE CORRUPTED AT LINE " << i << " : no weight specified !" << std::endl;
			 exit(1);
		   }
		weightList[k]=w;
		i++;
		  }
	  }
	else
	  {
		filePrint(stderr," *******************************************\n");
		filePrint(stderr," *** ERROR: Cannot open weight file %s ***\n", optarg);
			filePrint(stderr," *******************************************\n");
			exit(-1);
	  }
	break;
	  case 1005 :
	nosa = true;
	break;
	  case 1006 :
		debugFlag = 1;
		break;
	  case 1007 :
		trivialFlag = 1;
		break;
	  case 1008 :
		allowMechanisms = true;
		break;
      case 1009 :
        useScotch = true;
        break;
	  case 1010 :
	callSower = true;
	domain->setSowering(true);
	break;
	  case 1011 : {
		std::string prefix = optarg;
		clusterData_ = prefix + ".msh";
		decomposition_ = prefix + ".dec";
		connectivity_ = prefix + ".con";
		subdomains_ = prefix + ".sub";
		} break;
	  case 1012 :
		numClusters = atoi(optarg);
		if(numClusters <= 0) numClusters = 1;
		break;
      case 1013 :
        numCpus = atoi(optarg);
        if(numCpus <= 0) numCpus = 1;
        break;
	  case 'w':
	weightOutFlag = true;
		break;
	  case 'e':
		estFlag = true;
	useFull = true;
		break;
	  case 'n':
		numThreads = atoi(optarg);
		if(numThreads <= 0) numThreads = 1;
#ifdef _OPENMP
		omp_set_dynamic(0);
		omp_set_num_threads(numThreads);
#endif
		break;
	  case 'd': {
		  geoSource->getCheckFileInfo()->decomposition = optarg;
		  FILE *f;
		  if((f=fopen(optarg,"r"))==(FILE *) NULL ) {
			filePrint(stderr," *******************************************\n");
			filePrint(stderr," *** ERROR: Cannot open decomposition %s ***\n", optarg);
			filePrint(stderr," *******************************************\n");
			exit(-1);
		  }
		  geoSource->getCheckFileInfo()->decPtr = f;
		}
		break;
	  case 'v':
		verboseFlag = 1;
		domain->setVerbose();
		default_cntl.verbose = 1;
		default_cntl.fetiInfo.printNumber = atoi(optarg);
		filePrint(stderr," ... Setting Verbose Output Mode    ... \n");
		break;
	  case 'V':
		verboseFlag = 1;
		domain->setVerbose();
		default_cntl.verbose = 1;
		default_cntl.fetiInfo.printNumber = 10;
		filePrint(stderr," ... Setting Verbose Output Mode    ... \n");
		break;
	  case 't':
		topFlag = 0;
		domain->solInfo().setProbType(SolverInfo::Top);
	break;
	  case 'T':
		topFlag = 1;
		domain->solInfo().setProbType(SolverInfo::Top);
		break;
	  case 'm':
		topFlag = 2;
		domain->solInfo().setProbType(SolverInfo::Top);
		break;
	  case 'M':
		topFlag = 7;
		domain->solInfo().setProbType(SolverInfo::Top);
		break;
	  case 'x':
	domain->setOutputMatchInTop(true);
	break;
	  case 'r':
		topFlag = 3 + atoi(optarg);
		domain->solInfo().setProbType(SolverInfo::Top);
		break;
	  case 'p':
		default_cntl.fetiInfo.primalFlag = 1;
		break;
	  case 'c':
		contactPrintFlag = default_cntl.fetiInfo.contactPrintFlag = atoi(optarg);
		break;
	  case 's':
		domain->solInfo().printNumber = atoi(optarg);
		break;
	  case 'q':
		quietFlag = 1;
		break;
	  case '?':
		{
		  filePrint(stderr," *******************************************\n");
		  filePrint(stderr," *** ERROR: Check command line argument  ***\n");
		  filePrint(stderr," *******************************************\n");
		  exit(-1);
		}
	  }

 if(optind < argc - 1) {
   filePrint(stderr," ******************************************************\n");
   filePrint(stderr," *** ERROR: Command line contained errors. Aborting ***\n");
   filePrint(stderr," ******************************************************\n");
   exit(-1);
 }

 if(optind == argc) {
   filePrint(stderr," **************************************\n");
   filePrint(stderr," *** ERROR: No input file specified ***\n");
   filePrint(stderr," **************************************\n");
   exit(-1);
 }

 // Open input file
 if(optind == argc-1) {
   geoSource->getCheckFileInfo()->checkfile = argv[argc-1];
   //FILE *nin = freopen(argv[argc-1],"r",stdin);
   FILE *nin = yyin = fopen(argv[argc-1],"r");
   if(nin == 0) {
	 filePrint(stderr," *******************************************\n");
	 filePrint(stderr," *** ERROR: Could not open input file: %s\n",argv[argc-1]);
	 filePrint(stderr," *******************************************\n");
	 exit(-1);
   }
 }
 MatrixTimers &times = domain->getTimers();
 double t1 = getTime();

 // Read input file and time reading
 if(verboseFlag) filePrint(stderr," ... Reading Input File             ...\n");
 startTimerMemory(times.readTime, times.memoryParse);
 int error = yyparse();
 stopTimerMemory(times.readTime, times.memoryParse);
 if(verboseFlag) filePrint(stderr," ... Parsed Input File In %8.2e sec and %8.2e Mb ...\n",
					   (getTime() - t1)/1000.0, times.memoryParse/oneMegaByte);

 // Check if input file had any errors
 if(error) {
   filePrint(stderr," ****************************************************\n");
   filePrint(stderr," *** ERROR: Input file contained errors. Aborting ***\n");
   filePrint(stderr," ****************************************************\n");
   exit(error);
 }
 domain->buildSensitivityInfo();

 if(decInit != 0 && decInit->skip==false) { // dec initializers in parser !
   if(numThreads == 1) // command line option prevail !!!
	 numThreads = decInit->nthreads;
   if(numSubdomains == 1)
	 numSubdomains = decInit->nsubs;
   if(numProcessors == 1)
	 numProcessors = decInit->nproc;
   if(exitAfterDec==false)
	 exitAfterDec = decInit->exitAfterDec;
   if(decInit->weight)
	 weightOutFlag = true;
   if(decInit->trivial)
	 trivialFlag = true;
   if (decInit->fsgl)
	 fsglFlag = true;
   if(decInit->memory) {
	 estFlag = true;
	 useFull = true;
   }
   if(!nosa) nosa = decInit->nosa;
   if(topFlag < 0) {
	 callDec=true;
	 if(geoSource->getCheckFileInfo()->decPtr == 0 && decInit->file !=0)
	   geoSource->getCheckFileInfo()->checkfile = decInit->file;
   }
 }

 if(!domain->solInfo().readInModes.empty() || domain->solInfo().modal || domain->solInfo().modalCalled) {
   if((domain->solInfo().isNonLin() && !domain->solInfo().modal_id.empty()) || domain->solInfo().activatePodRom) {
	 for(int i=0; i<domain->solInfo().modal_id.size(); ++i) {
	   ModalParams &modalParams = domain->solInfo().readInModes[domain->solInfo().modal_id[i]];
	   switch(modalParams.type) {
		 case ModalParams::Inorm : {
		   domain->solInfo().readInROBorModes.push_back(modalParams.fileName);
		   domain->solInfo().localBasisSize.push_back(modalParams.numvec);
		   domain->solInfo().maxSizePodRom += modalParams.numvec;
		   domain->solInfo().useMassNormalizedBasis = false;
		 } break;
		 case ModalParams::Mnorm : {
		   std::string::size_type n = modalParams.fileName.rfind(".massorthonormalized");
		   if(n != std::string::npos) {
			 domain->solInfo().readInROBorModes.push_back(modalParams.fileName.substr(0,n));
			 domain->solInfo().localBasisSize.push_back(modalParams.numvec);
			 domain->solInfo().maxSizePodRom += modalParams.numvec;
			 domain->solInfo().useMassNormalizedBasis = true;
		   }
		   else {
			 filePrint(stderr, " *** ERROR: Specified filename for rob_id#%d is missing \".massorthonormalized\" extension.\n",i+1);
			 exit(-1);
		   }
		 } break;
		 case ModalParams::Undefined : {
		   domain->solInfo().readInROBorModes.push_back(modalParams.fileName);
		   domain->solInfo().localBasisSize.push_back(modalParams.numvec);
		   domain->solInfo().maxSizePodRom += modalParams.numvec;
		   domain->solInfo().useMassNormalizedBasis = true;
		 } break;
		 default : {
		   filePrint(stderr, " *** ERROR: Specified type for rob_id#%d is not supported.\n", i+1);
		   exit(-1);
		 } break;
	   }
	 }
	 for(int i=0; i<domain->solInfo().contact_modal_id.size(); ++i) {
	   ModalParams &modalParams = domain->solInfo().readInModes[domain->solInfo().contact_modal_id[i]];
	   switch(modalParams.type) {
		 case ModalParams::Noneg : {
		   domain->solInfo().readInDualROB.push_back(modalParams.fileName);
		   domain->solInfo().localDualBasisSize.push_back(modalParams.numvec);
		 } break;
		 default : {
		   filePrint(stderr, " *** ERROR: Specified type for rob_id#%d is not supported.\n", domain->solInfo().modal_id.size()+i+1);
		   exit(-1);
		 } break;
	   }
	 }
	 if(!domain->solInfo().samplingPodRom) {
	   domain->solInfo().activatePodRom = true;
	   domain->solInfo().galerkinPodRom = true;
	   if(!domain->solInfo().ROMPostProcess || (!domain->solInfo().readInAdjointROB.empty())) {
#ifdef USE_EIGEN3
		 if(domain->solInfo().solvercntl->subtype != 12) domain->solInfo().solvercntl->subtype = 13;
#else
		 domain->solInfo().solvercntl->subtype = 12;
#endif
	   }
	 }
   }
   else if(!domain->solInfo().modal_id.empty() || domain->solInfo().readInModes.find(0) != domain->solInfo().readInModes.end()) {
	 int modal_id = (domain->solInfo().modal_id.empty()) ? 0 : domain->solInfo().modal_id.front();
	 domain->readInModes(modal_id, modeData);
   }
   if(domain->solInfo().idis_modal_id != -1)
	 domain->readInModes(domain->solInfo().idis_modal_id, modeDataIDis);
   if(domain->solInfo().ivel_modal_id != -1)
	 domain->readInModes(domain->solInfo().ivel_modal_id, modeDataIVel);
   if(domain->solInfo().mode_modal_id != -1)
	 domain->readInModes(domain->solInfo().mode_modal_id, modeDataMode);
 }

 if(domain->solInfo().readShapeSen) {
   domain->readInShapeDerivatives(const_cast<char*>(domain->solInfo().readInShapeSen));
 }

#define MAX_CODES 4
#define FLUID_ID 0
#define STRUC_ID 1
#define HEAT_ID  2

 // We do a split
 Communicator* allCom[MAX_CODES];
 if(domain->solInfo().aeroheatFlag >= 0 || domain->solInfo().thermohFlag >= 0) {
   syscom->split(HEAT_ID, MAX_CODES, allCom);
   structCom = allCom[HEAT_ID];
   heatStructCom = allCom[STRUC_ID]; // comunicator between the thermal and
									 // mechanical structure codes
 }
 else {
   syscom->split(STRUC_ID, MAX_CODES, allCom);
   structCom = allCom[STRUC_ID];
   heatStructCom = allCom[HEAT_ID]; // comunicator between the thermal and
									// mechanical structure codes
 }
 fluidCom = allCom[FLUID_ID];
#ifdef PRINT_CHANGESETID
 filePrint(stderr," ... Changeset ID%15s    ...\n",THE_VERSION);
#endif

 threadManager = new ThreadManager(numThreads);
 if(threadManager->numThr() != numThreads) { //HB: for checking purpose
   filePrint(stderr," *** WARNING: number of threads requested: %d\n",numThreads);
   filePrint(stderr,"              number of threads created  : %d\n",threadManager->numThr());
   filePrint(stderr," -> tip: if you are running on a Linux platform\n");
   filePrint(stderr,"         you may need to activate OpenMP and compile with an OpenMP compliant\n");
   filePrint(stderr,"         compiler (for instance, icpc or g++ version 4.2)\n");
 }

 if(!domain->solInfo().basicPosCoords || domain->solInfo().scalePosCoords) {
#ifndef USE_EIGEN3
   filePrint(stderr," *** ERROR: use of Nodal Frames for node coordinates is not supported by this AERO-S build.\n");
   filePrint(stderr," -> tip: you need to configure AERO-S with \"cmake -DEIGEN3_INCLUDE_DIR:PATH=xxx .\"\n");
   filePrint(stderr,"         where xxx is the path to an installation of the Eigen C++ template library (version 3.1 or later)\n");
   exit(-1);
#endif
   geoSource->transformCoords();
 }

 if(geoSource->binaryInput)
 	geoSource->readGlobalBinaryData();
#ifdef SOWER_SURFS
 else {
#endif
   // HB for checking Mortar & generating LMPCs from Mortar tied conditions
   if(topFlag < 0) {
	 domain->SetUpSurfaces(&(geoSource->GetNodes()));
     domain->ProcessSurfaceBCs(topFlag);
     if(!callSower) {
	   domain->SetMortarPairing();
	   if(!domain->tdenforceFlag()) {
         if(domain->solInfo().nlFlag == 1) { // for nonlinear statics and dynamics just process the tied surfaces here
		   domain->InitializeStaticContactSearch(MortarHandler::TIED);
		   domain->PerformStaticContactSearch(MortarHandler::TIED);
		   domain->ExpComputeMortarLMPC(MortarHandler::TIED); // TODO thermal and acoustic (see below)
		   domain->CreateMortarToMPC();
		 }
		 else {
		   switch(domain->solInfo().soltyp) { // TODO: acoustic etc...
			 case 2 : { // thermal mortar
			   int dofs[1] = { 6 };
			   domain->ComputeMortarLMPC(1,dofs);
			   break;
			 }
			 default :
			 case 1 :
			   int dofs[3] = { 0, 1, 2 };
			   domain->ComputeMortarLMPC(3, dofs);
			   break;
		   }
		   domain->computeMatchingWetInterfaceLMPC();
		   domain->CreateMortarToMPC();
		 }
	   }
	 }
#ifdef MORTAR_DEBUG
	 domain->PrintSurfaceEntities();
	 domain->PrintMortarConds();
	 domain->printLMPC();
#endif
   }
   else domain->ProcessSurfaceBCs(topFlag);
#ifdef SOWER_SURFS
 }
#endif
 if(!domain->solInfo().basicDofCoords) {
#ifndef USE_EIGEN3
   filePrint(stderr," *** ERROR: use of Nodal Frames for degrees of freedom is not supported by this AERO-S build.\n");
   filePrint(stderr," -> tip: you need to configure AERO-S with \"cmake -DEIGEN3_INCLUDE_DIR:PATH=xxx .\"\n");
   filePrint(stderr,"         where xxx is the path to an installation of the Eigen C++ template library (version 3.1 or later)\n");
   exit(-1);
#endif
   geoSource->transformLMPCs(domain->getNumLMPC(), *(domain->getLMPC()));
 }

 if(!isFeti(domain->solInfo().solvercntl->type) && !domain->solInfo().use_nmf && !domain->solInfo().svdPodRom && domain->solInfo().readInDualROB.empty()
	&& !domain->solInfo().samplingPodRom)
   geoSource->addMpcElements(domain->getNumLMPC(), *(domain->getLMPC()));

 if((!isFeti(domain->solInfo().solvercntl->type) || (!domain->solInfo().isMatching && (domain->solInfo().solvercntl->fetiInfo.fsi_corner != 0))) && !domain->solInfo().HEV)
   geoSource->addFsiElements(domain->getNumFSI(), domain->getFSI());

 if(!geoSource->binaryInput) {
   domain->setUpData(topFlag);
 }

 if(!(geoSource->getCheckFileInfo()->decPtr || callDec || decInit || geoSource->binaryInput)) {
   // activate multi-domain mode for the explicit dynamics Rom drivers which are not supported in single-domain mode
   // so it is not necessary to include "DECOMP" with "nsubs 1" in the input file
   if((domain->solInfo().activatePodRom && domain->probType() == SolverInfo::NonLinDynam && domain->solInfo().newmarkBeta == 0)
	  || (domain->probType() == SolverInfo::PodRomOffline && domain->solInfo().ROMPostProcess) || domain->solInfo().clustering > 0
	  || domain->solInfo().rowClustering > 0) {
	 callDec = true;
	 trivialFlag = true;
	 numSubdomains = 1;
   }
 }

 if(callDec) {
   Dec::dec(numProcessors, numThreads, numSubdomains, topFlag);
   if(exitAfterDec && !callSower) {
	 filePrint(stderr," ... Exiting after Dec run          ...\n");
	 filePrint(stderr," --------------------------------------\n");
	 delete threadManager;
	 closeComm();
	 exit(0);
   }
 }
 useFull = true; // or TenNodeTetraHedral will crush ! (bad design not from me !)

 if(callSower) {
   filePrint(stderr," ... Writing Distributed Binary Input Files ... \n");
   geoSource->writeDistributedInputFiles(numClusters, domain, numCpus); //add domain as argument for surfaces
   if(exitAfterSower) {
	 filePrint(stderr," ... Exiting after Sower run        ...\n");
	 filePrint(stderr," --------------------------------------\n");
	 closeComm();
	 exit(0);
   }
 }

 // Make TOPDOM/DEC input file and exit from fem
 if(topFlag >= 0) {
   if(topFlag == 5 || topFlag == 6)
	   fprintf(stderr," ... Using removed HelmAXI commands.\n");
   else if(domain->probType() == SolverInfo::Top) {
	 fprintf(stderr," ... Memory to Parse       %14.3f Mb\n",
			 times.memoryParse/oneMegaByte);
	 fprintf(stderr," ... Memory to Set up Data %14.3f Mb\n",
			 times.memorySetUp/oneMegaByte);
	 domain->makeTopFile(topFlag);
   }
   closeComm();
   exit(0);
 }

 if(domain->solInfo().noninpc || domain->solInfo().inpc) domain->initSfem();

 // 1. check to see if a decomposition has been provided or requested (three options: -d, --dec, or DECOMP)
 bool domain_decomp = (geoSource->getCheckFileInfo()->decPtr || callDec || decInit || geoSource->binaryInput);
 // 2. check to see how many cpus are available (mpi processes and threads)
#ifdef USE_MPI
 bool parallel_proc = (structCom->numCPUs()*threadManager->numThr() > 1 || geoSource->binaryInput);
#else
 bool parallel_proc = (threadManager->numThr() > 1);
#endif
 // 3. choose lumped mass (also pressure and gravity) and diagonal or block-diagonal "solver" for explicit dynamics
 if(domain->solInfo().newmarkBeta == 0 || (domain->solInfo().svdPodRom && geoSource->getMRatio() == 0)) {
   if((parallel_proc || domain_decomp)
      && (domain->solInfo().solvercntl->type != SolverSelection::Feti)
      && (domain->solInfo().solvercntl->type != SolverSelection::BlockDiag))
   {
     domain->solInfo().solvercntl = new SolverCntl(default_cntl);
   }
   if(domain->solInfo().inertiaLumping == 2) { // block-diagonal lumping
	 domain->solInfo().solvercntl->subtype = 1;
	 domain->solInfo().getFetiInfo().local_cntl->subtype = FetiInfo::sparse;
	 if(parallel_proc || domain_decomp)
	 	domain->solInfo().solvercntl->type = SolverSelection::Feti; // XXX type 3 could be upgraded to work for this case,
																				//     rather than using feti
   }
   else {
	 domain->solInfo().solvercntl->subtype = 10;
	 domain->solInfo().getFetiInfo().local_cntl->subtype = FetiInfo::diagonal;
	 if(parallel_proc || domain_decomp)
	 	domain->solInfo().solvercntl->type = SolverSelection::BlockDiag;
   }

   geoSource->setMRatio(0.0);
   geoSource->setConsistentQFlag(false);
   geoSource->setConsistentPFlag(false);
   if(verboseFlag) filePrint(stderr, " ... Explicit Dynamics: lumped mass matrix, gravity and pressure will be used ... \n");
 }

 if(domain->solInfo().aeroFlag >= 0)
   filePrint(stderr," ... AeroElasticity Flag   = %2d     ...\n", domain->solInfo().aeroFlag);
 if(domain->solInfo().thermoeFlag >= 0)
   filePrint(stderr," ... ThermoElasticity Flag = %2d     ...\n", domain->solInfo().thermoeFlag);
 if(domain->solInfo().aeroheatFlag >= 0)
   filePrint(stderr," ... AeroThermo Flag       = %2d     ...\n", domain->solInfo().aeroheatFlag);
 if(domain->solInfo().thermohFlag >= 0)
   filePrint(stderr," ... ThermoElasticity Flag = %2d     ...\n", domain->solInfo().thermohFlag);

 // ... PRINT PROBLEM TYPE
 if(domain->solInfo().sensitivity) {
   if(domain->runSAwAnalysis) {
	 filePrint(stderr, sensitivityTypeMessage[1]);
	 filePrint(stderr, problemTypeMessage[domain->solInfo().probType]);
   } else {
	 filePrint(stderr, sensitivityTypeMessage[0]);
   }
 } else {
   filePrint(stderr, problemTypeMessage[domain->solInfo().probType]);
 }

 if(domain->solInfo().gepsFlg == 1) {
   if(domain->solInfo().isNonLin()) { // GEPS is not used for nonlinear
	 domain->solInfo().gepsFlg = 0;
   }
   else filePrint(stderr," ...      with Geometric Pre-Stress ... \n");
 }
 if(domain->solInfo().solvercntl->type == SolverSelection::Direct
    && domain->solInfo().probType != SolverInfo::None
    && domain->solInfo().probType != SolverInfo::PodRomOffline
	&& domain->solInfo().modal_id.empty())
   filePrint(stderr, solverTypeMessage[domain->solInfo().solvercntl->subtype]);

 // Domain Decomposition tasks
 //   type == 2 (FETI) and type == 3 (BLOCKDIAG) are always Domain Decomposition methods
 //   type == 1 && iterType == 1 (GMRES) is a Domain Decomposition method only if a decomposition is provided or requested
 //   type == 0 && subtype == 9 (MUMPS) is a Domain Decomposition method only if a decomposition is provided or requested
 if(isFeti(domain->solInfo().solvercntl->type)
    || domain->solInfo().solvercntl->type == SolverSelection::BlockDiag
    || (domain->solInfo().solvercntl->type == SolverSelection::Iterative && domain->solInfo().solvercntl->iterType == 1 && domain_decomp)
	|| (domain->solInfo().solvercntl->type == SolverSelection::Direct && domain->solInfo().solvercntl->subtype == 9 && domain_decomp)
	|| (!domain->solInfo().modal_id.empty() && domain_decomp)
	|| (domain->solInfo().svdPodRom && domain_decomp)
	|| (domain->solInfo().samplingPodRom && domain_decomp)) {

   if(parallel_proc) {
#ifdef USE_MPI
	 if(structCom->numCPUs() > 1) {
	   if(threadManager->numThr() == 1)
		 filePrint(stderr, " ... Parallel processing: Launching %d MPI processes ...\n", structCom->numCPUs());
	   else
		 filePrint(stderr, " ... Parallel processing: Launching %d MPI processes and forking %d threads per MPI process ...\n", structCom->numCPUs(), threadManager->numThr());
	 } else
#endif
	 filePrint(stderr, " ... Parallel processing: Forking %d Threads ...\n", threadManager->numThr());
   }

   switch(domain->probType()) {
	 case SolverInfo::TempDynamic:
	   {
		 MultiDomainTemp tempProb(domain);
		 TempSolver<MDDynamMat, DistrVector, MDTempDynamPostProcessor,
					MultiDomainTemp> dynaSolver(&tempProb);
		 dynaSolver.solve();
	   }
	   break;
	 case SolverInfo::Dynamic:
	   {
		if(domain->solInfo().modal) {
		  filePrint(stderr," ... Modal Method                   ...\n");
		  MultiDomainModal * modalProb = new MultiDomainModal(domain);
		  DynamicSolver<ModalOps, Vector, MultiDomainModal, MultiDomainModal, double>
		  modalSolver(modalProb);
		  modalSolver.solve();
		}
		else {
		  if(domain->solInfo().mdPita) { // Not implemented yet
			filePrint(stderr, " ... PITA does not support multidomain - Aborting...\n");
		  } else {
			MultiDomainDynam dynamProb(domain);
			DynamicSolver < MDDynamMat, DistrVector, MultiDomDynPostProcessor,
							MultiDomainDynam, double > dynamSolver(&dynamProb);
			dynamSolver.solve();
			fflush(stderr);
		  }
		 }
	   }
	   break;

	  case SolverInfo::NonLinStatic: {
		MDNLStatic nlstatic(domain);
		NLStaticSolver < ParallelSolver,DistrVector,MultiDomainPostProcessor,
					  MDNLStatic, DistrGeomState >
					  nlsolver(&nlstatic);
		nlsolver.solve();
	}
	   break;
	 case SolverInfo::ArcLength: {
	MDNLStatic nlstatic(domain);
	NLStaticSolver < ParallelSolver,DistrVector,MultiDomainPostProcessor,
					  MDNLStatic, DistrGeomState >
					  nlsolver(&nlstatic);
	nlsolver.arclength();
	}
	   break;
	 case SolverInfo::Modal: {
		GenMultiDomainEigen<double> eigenProb(domain);
		EigenSolver<MDDynamMat, GenDistrVector<double>, GenDistrVectorSet<double>, GenMultiDomainEigenPostProcessor<double>,
					GenMultiDomainEigen<double> > * eigenSolver = 0;
		switch(domain->solInfo().eigenSolverType) {
		   case SolverInfo::Arpack :
			 {
#ifdef USE_ARPACK
		   eigenSolver = new SymArpackSolver<MDDynamMat, GenDistrVector<double>, GenDistrVectorSet<double>,
												GenMultiDomainEigenPostProcessor<double>, GenMultiDomainEigen<double> > (&eigenProb);
#else
			  filePrint(stderr," *** ERROR: executable not linked with ARPACK. See flag USE_ARPACK in Makefile.\n");
			  exit(-1);
#endif
			 }
			 break;
		   default :
			 {
			  filePrint(stderr,"ERROR: Eigensolver type is not implemented for multiple domain. Please use ARPACK.\n");
			  exit(-1);
			 }
		   }
		   eigenSolver->solve();
	   }
	   break;
	 case SolverInfo::Static: {
	   if(geoSource->isShifted()) filePrint(stderr, " ... Frequency Response Analysis ");
	   if(domain->isComplex()) {
		 if(domain->solInfo().inpc) std::cerr << "inpc not implemented for complex domain \n";
		 if(geoSource->isShifted()) filePrint(stderr, "in Complex Domain ...\n");
	   {
		 GenMultiDomainStatic<complex<double> > statProb(domain);
		 StaticSolver<complex<double>, GenMDDynamMat<complex<double> >, GenDistrVector<complex<double> >,
					  GenMultiDomainPostProcessor<complex<double> >, GenMultiDomainStatic<complex<double> >,
					  GenDistrVector<complex<double> > >
		   statSolver(&statProb);
		 statSolver.solve();
	   }
	   }
	   else {
		 if(geoSource->isShifted()) filePrint(stderr, "in Real Domain ...\n");
		 if(domain->solInfo().inpc) {
		   GenMultiDomainInpcStatic<double> statProb(domain);
		   StaticSolver<double, AllOps<double>, DistrBlockVector<double>,
						GenMultiDomainInpcPostProcessor<double>, GenMultiDomainInpcStatic<double>,
						DistrBlockVector<complex<double> > >
				statSolver(&statProb);
		   statSolver.solve();
		 }
		 else {
	 {
	   GenMultiDomainStatic<double> statProb(domain);
	   StaticSolver<double, GenMDDynamMat<double>, GenDistrVector<double>,
					GenMultiDomainPostProcessor<double>, GenMultiDomainStatic<double>,
					GenDistrVector<complex<double> > >
			statSolver(&statProb);
	   statSolver.solve();
	 }
		 }
	   }
	   break;
	 }
	 case SolverInfo::HelmholtzFreqSweep:
	 case SolverInfo::Helmholtz: {
	   filePrint(stderr, " ... Acoustic Scattering Helmholtz Analysis ");
	   if(domain->isComplex()) {
		 filePrint(stderr, "in Complex Domain ...\n");
		 GenMultiDomainStatic<complex<double> > FAProb(domain);
		 StaticSolver<complex<double>, GenMDDynamMat<complex<double> >, GenDistrVector<complex<double> >,
					  GenMultiDomainPostProcessor<complex<double> >, GenMultiDomainStatic<complex<double> >,
					  GenDistrVector<complex<double> > >
			FASolver(&FAProb);
		 FASolver.solve();
	   }
	   else {
		 filePrint(stderr,"in Real Domain ...\n");
		 GenMultiDomainStatic<double> FAProb(domain);
		 StaticSolver<double, GenMDDynamMat<double>, GenDistrVector<double>,
					  GenMultiDomainPostProcessor<double>, GenMultiDomainStatic<double>,
					  GenDistrVector<complex<double> > >
			  FASolver(&FAProb);
		 FASolver.solve();
	   }
	 } break;
	 case SolverInfo::NonLinDynam: {
	   if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().timeIntegration == 1) { // explicit or quasi-static
		 if (!domain->solInfo().activatePodRom) {
		   if(domain->solInfo().soltyp == 2) {
			 MultiDomainTemp tempProb(domain);
			 TempSolver<MDDynamMat, DistrVector, MDTempDynamPostProcessor,
						MultiDomainTemp> dynamSolver(&tempProb);
			 dynamSolver.solve();
		   }
		   else {
			 MultiDomainDynam dynamProb(domain);
			 DynamicSolver<MDDynamMat, DistrVector, MultiDomDynPostProcessor,
						   MultiDomainDynam, double> dynamSolver(&dynamProb);
			 dynamSolver.solve();
		   }
		 }
		 else { // POD ROM
		   if (domain->solInfo().galerkinPodRom) {
			 if (domain->solInfo().elemLumpPodRom) {
			   if(domain->solInfo().DEIMPodRom){
				if (domain->solInfo().reduceFollower)
				   filePrint(stderr, " ... POD: ROM with stiffness & follower interpolation ...\n");
				else
				   filePrint(stderr, " ... POD: ROM with stiffness interpolation ...\n");
				Rom::DistrExplicitDEIMPodProjectionNonLinDynamic dynamProb(domain);
				DynamicSolver < MDDynamMat, DistrVector, Rom::MultiDomDynPodPostProcessor,
								Rom::DistrExplicitDEIMPodProjectionNonLinDynamic, double > dynamSolver(&dynamProb);
				dynamSolver.solve();
			   } else if(domain->solInfo().UDEIMPodRom) {
				if (domain->solInfo().reduceFollower)
				   filePrint(stderr, " ... POD: unassembled ROM with stiffness & follower interpolation ...\n");
				else
				   filePrint(stderr, " ... POD: unassembled ROM with stiffness interpolation ...\n");
				Rom::DistrExplicitUDEIMPodProjectionNonLinDynamic dynamProb(domain);
				DynamicSolver < MDDynamMat, DistrVector, Rom::MultiDomDynPodPostProcessor,
								Rom::DistrExplicitUDEIMPodProjectionNonLinDynamic, double > dynamSolver(&dynamProb);
				dynamSolver.solve();
			   } else {
				if (domain->solInfo().reduceFollower)
				 filePrint(stderr, " ... POD: ROM with stiffness & follower lumping ...\n");
			else
				 filePrint(stderr, " ... POD: ROM with stiffness lumping...\n");
				Rom::DistrExplicitLumpedPodProjectionNonLinDynamic dynamProb(domain);
				DynamicSolver < MDDynamMat, DistrVector, Rom::MultiDomDynPodPostProcessor,
							   Rom::DistrExplicitLumpedPodProjectionNonLinDynamic, double > dynamSolver(&dynamProb);
				dynamSolver.solve();
			   }
			 } else {
			   filePrint(stderr, " ... POD: Explicit Galerkin         ...\n");
			   Rom::DistrExplicitPodProjectionNonLinDynamic dynamProb(domain);
			   DynamicSolver < MDDynamMat, DistrVector, Rom::MultiDomDynPodPostProcessor,
							 Rom::DistrExplicitPodProjectionNonLinDynamic, double > dynamSolver(&dynamProb);
			   dynamSolver.solve();
			 }
		   }
		   else {
			 filePrint(stderr, " ... POD: Snapshot collection       ...\n");
			 Rom::DistrExplicitSnapshotNonLinDynamic dynamProb(domain);
			 DynamicSolver < MDDynamMat, DistrVector, MultiDomDynPostProcessor,
				   Rom::DistrExplicitSnapshotNonLinDynamic, double > dynamSolver(&dynamProb);
			 dynamSolver.solve();
		   }
		 }
	   }
	   else { // implicit
		 if (!domain->solInfo().activatePodRom) {
		   MDNLDynamic nldynamic(domain);
		   NLDynamSolver <ParallelSolver, DistrVector, MultiDomainPostProcessor,
						  MDNLDynamic, DistrGeomState> nldynamicSolver(&nldynamic);
		   nldynamicSolver.solve();
		 }
		 else { // POD ROM
		   if (domain->solInfo().galerkinPodRom) {
#ifdef USE_EIGEN3
			 if(domain->solInfo().elemLumpPodRom) {
			   filePrint(stderr, " ... POD: Parallel Implicit Galerkin HROM   ...\n");
			   Rom::DistrLumpedPodProjectionNonLinDynamic nldynamic(domain);
			   NLDynamSolver <Rom::GenEiSparseGalerkinProjectionSolver<Scalar,GenDistrVector,GenParallelSolver<Scalar> >,
							  DistrVector, MultiDomainPostProcessor,
							  Rom::DistrPodProjectionNonLinDynamic, DistrModalGeomState,
							  Rom::DistrPodProjectionNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
			   nldynamicSolver.solve();
			 } else {
			   filePrint(stderr, " ... POD: Parallel Implicit Galerkin ROM   ...\n");
			   Rom::DistrPodProjectionNonLinDynamic nldynamic(domain);
			   NLDynamSolver <Rom::GenEiSparseGalerkinProjectionSolver<Scalar,GenDistrVector,GenParallelSolver<Scalar> >,
							  DistrVector, MultiDomainPostProcessor,
							  Rom::DistrPodProjectionNonLinDynamic, DistrModalGeomState,
							  Rom::DistrPodProjectionNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
			   nldynamicSolver.solve();
			 }
#endif
		   } else {
			 filePrint(stderr, " ... POD: Snapshot collection       ...\n");
			 Rom::DistrSnapshotNonLinDynamic nldynamic(domain);
			 NLDynamSolver <ParallelSolver, DistrVector, MultiDomainPostProcessor,
						   Rom::DistrSnapshotNonLinDynamic, DistrGeomState,
						   Rom::DistrSnapshotNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
			 nldynamicSolver.solve();
			 nldynamic.postProcess();
		   }
		 }
	   }
	 } break;
	 case SolverInfo::PodRomOffline: {
	   Rom::DriverInterface *driver;
	   if (domain->solInfo().svdPodRom) {
		 if(domain->solInfo().use_nmf) {
		   filePrint(stderr, " ... Nonneg. Matrix Factorization   ...\n");
		   driver = distrPositiveDualBasisDriverNew(domain);
		 }
		 else if(domain->solInfo().clustering) {
		   filePrint(stderr, " ... Snapshot Clustering            ...\n");
		   driver = distrSnapshotClusteringDriverNew(domain);
		 }
		 else if(domain->solInfo().rowClustering) {
		   filePrint(stderr, " ... Snapshot Row Clustering        ...\n");
		   driver = distrSnapshotRowClusteringDriverNew(domain);
		 }
		 else {
		   // Stand-alone SVD orthogonalization
		   filePrint(stderr, " ... Singular Value Decomposition   ...\n");
		   driver = distrBasisOrthoDriverNew(domain);
		 }
	   }
	   else if (domain->solInfo().ROMPostProcess) {
		 filePrint(stderr, " ... Post-processing ROM Results    ...\n");
		 driver = distrROMPostProcessingDriverNew(domain);
	   }
	   else if (domain->solInfo().samplingPodRom) {
		 // Element-based hyper-reduction
		 filePrint(stderr, " ... Element Sampling and Weighting ...\n");
		 driver = distrElementSamplingDriverNew(domain);
	   }
	   else {
		 filePrint(stderr, " ... Unknown Analysis Type          ...\n");
		 break;
	   }
	   driver->solve();
	   delete driver;
	   break;
	 }
	 // Fall-thru
	 default:
	   filePrint(stderr,"*** ERROR: Problem type %d is not supported multi-domain mode\n", domain->probType());
   }

   totalMemoryUsed = double(memoryUsed()+totMemSpooles+totMemMumps)/oneMegaByte;
   delete threadManager;

 } else {

   switch(domain->probType()) {
	 case SolverInfo::DisEnrM: {
		filePrint(stderr, " ... DEM Problem                    ...\n");
		DEM dem;
		dem.run(domain,geoSource);
		break;
	 }
	 case SolverInfo::Static:
	   {
		 if(geoSource->isShifted()) filePrint(stderr, " ... Frequency Response Helmholtz Analysis ");
		 if(domain->isComplex()) {
		   if(geoSource->isShifted()) filePrint(stderr, "in Complex Domain ...\n");
		   SingleDomainStatic<DComplex, GenVector<DComplex>, GenSolver<DComplex> >
			 statProb(domain);
		   StaticSolver<DComplex, AllOps<DComplex>, /*GenSolver<DComplex>,*/ GenVector<DComplex>,
						SingleDomainPostProcessor<DComplex, GenVector<DComplex>, GenSolver<DComplex> >,
						SingleDomainStatic<DComplex, GenVector<DComplex>, GenSolver<DComplex> >, GenVector<DComplex> >
			 statSolver(&statProb);
		   statSolver.solve();
		 } else {
		   if(geoSource->isShifted()) filePrint(stderr, "in Real Domain ...\n");
		   SingleDomainStatic<double, Vector, Solver> statProb(domain);
		   StaticSolver<double, AllOps<double>, /*Solver,*/ Vector,
				SingleDomainPostProcessor<double, Vector, Solver>,
			SingleDomainStatic<double, Vector, Solver>, GenVector<DComplex> >
			 statSolver(&statProb);
		   statSolver.solve();
		 }
	   }
	   break;

	 case SolverInfo::Dynamic:
	   {
		if(domain->solInfo().modal) {
		  fprintf(stderr," ... Modal Method                   ...\n");
		  ModalDescr<double> * modalProb = new ModalDescr<double>(domain);
		  DynamicSolver<ModalOps, Vector, ModalDescr<double>, ModalDescr<double>, double>
		  modalSolver(modalProb);
		  modalSolver.solve();
		}
		else {
		  if (domain->solInfo().activatePita) {
#ifdef USE_PITA
			 SingleDomainDynamic dynamProb(domain);
			 Pita::LinearDriver::Ptr driver;
			 if (!domain->solInfo().pitaTimeReversible) {
				 filePrint(stderr," ... Linear PITA                    ...\n");
				 driver = linearPitaDriverNew(&dynamProb);
			 } else {
				 filePrint(stderr," ... Time-reversible linear PITA    ...\n");
				 driver = linearReversiblePitaDriverNew(&dynamProb);
			 }
			 driver->solve();
#else
			 filePrint(stderr," ... PITA requires distributed version ...\n");
#endif
	  } else {
			if (domain->solInfo().ATDARBFlag>=1.5) {
/*              SingleDomainDynamic dynamProb(domain);
			  DynamicSolver <GenDynamMat<DComplex>, Vector,
					SDDynamPostProcessor, SingleDomainDynamic, DComplex>
				dynaSolver(&dynamProb);
			  dynaSolver.solve();
*/            }
			else {
			  SingleDomainDynamic dynamProb(domain);
			  DynamicSolver <GenDynamMat<double>, Vector,
					SDDynamPostProcessor, SingleDomainDynamic, double>
				dynaSolver(&dynamProb);
			  dynaSolver.solve();
			}
	  }
	}
	  }
	  break;
	 case SolverInfo::TempDynamic:
	   {
		 SingleDomainTemp tempProb(domain);
		 TempSolver<DynamMat, Vector, SDTempDynamPostProcessor,
					SingleDomainTemp> dynaSolver(&tempProb);
		 dynaSolver.solve();
	   }
	   break;
	 case SolverInfo::Modal:
	   {
	  SingleDomainEigen eigenProb(domain);
		 EigenSolver<DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> *eigenSolver;
		  switch(domain->solInfo().eigenSolverType) {
			case SolverInfo::SubSpace :
			  eigenSolver = new SubSpaceSolver <DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> (&eigenProb);
			  break;
			case SolverInfo::LobPcg :
			  eigenSolver = new LOBPCGSolver<DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> (&eigenProb);
			  break;
			case SolverInfo::Arpack :
			  {
#ifdef USE_ARPACK
		eigenSolver = new SymArpackSolver <DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> (&eigenProb);
#else
				  filePrint(stderr," *** ERROR: executable not linked with ARPACK. See flag USE_ARPACK in Makefile.\n");
				  exit(-1);
#endif
			  }
			  break;
			default :
				filePrint(stderr," *** ERROR: unknown eigensolver.\n");
				exit(-1);
			}
	  eigenSolver->solve();
	   }
	   break;
	 case SolverInfo::NonLinStatic:
	   {
	 NonLinStatic nlstatic(domain);
		 NLStaticSolver<Solver, Vector, SingleDomainPostProcessor<double, Vector, Solver>,
				NonLinStatic, GeomState>
		   nlsolver(&nlstatic);
	 nlsolver.solve();
	   }
	   break;
	 case SolverInfo::NonLinDynam: {
		 if(domain->solInfo().activatePita) {
#ifdef USE_PITA
		   if (domain->solInfo().pitaTimeReversible) {
			 filePrint(stderr, " ... Time-reversible nonlinear PITA ...\n");
			 Pita::PitaNonLinDynamic pitaProblem(domain);
			 Pita::NlDriver::Ptr pitaDriver = nlReversiblePitaDriverNew(&pitaProblem);
			 pitaDriver->solve();
		   } else {
			 filePrint(stderr, " ... Nonlinear PITA                 ...\n");
			 Pita::Old::PitaNonLinDynamic pitaProblem(domain);
			 Pita::Old::NLDistrTimeDecompSolver pitaSolver(&pitaProblem);
			 pitaSolver.solve();
		   }
#else
		   filePrint(stderr," ... PITA requires distributed version ...\n");
#endif
		 }
		 else {
		   if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().timeIntegration == 1) { // explicit or quasi-static
			 if(domain->solInfo().soltyp == 2) {
			   SingleDomainTemp tempProb(domain);
			   TempSolver<DynamMat, Vector, SDTempDynamPostProcessor,
						  SingleDomainTemp> dynaSolver(&tempProb);
			   dynaSolver.solve();
			 }
			 else {
			   SingleDomainDynamic dynamProb(domain);
			   DynamicSolver<GenDynamMat<double>, Vector, SDDynamPostProcessor,
							 SingleDomainDynamic, double> dynaSolver(&dynamProb);
			   dynaSolver.solve();
			 }
		   }
		   else { // implicit
			 if (!domain->solInfo().activatePodRom) {
			   NonLinDynamic nldynamic(domain);
			   NLDynamSolver <Solver, Vector, SDDynamPostProcessor, NonLinDynamic, GeomState> nldynamicSolver(&nldynamic);
			   nldynamicSolver.solve();
			 } else { // POD ROM
			   if (domain->solInfo().galerkinPodRom && domain->solInfo().elemLumpPodRom) {
				 if(domain->solInfo().DEIMPodRom){
				   if (domain->solInfo().reduceFollower)
					 filePrint(stderr, " ... POD: ROM with stiffness & follower interpolation ...\n");
				   else
					 filePrint(stderr, " ... POD: ROM with stiffness interpolation ...\n");
				   Rom::DEIMPodProjectionNonLinDynamic nldynamic(domain);
				   NLDynamSolver <Rom::PodProjectionSolver, Vector, SDDynamPostProcessor, Rom::PodProjectionNonLinDynamic,
								  ModalGeomState, Rom::PodProjectionNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
				   nldynamicSolver.solve();
				 } else if(domain->solInfo().UDEIMPodRom) {
				   if (domain->solInfo().reduceFollower)
					 filePrint(stderr, " ... POD: unassembled ROM with stiffness & follower interpolation ...\n");
				   else
					 filePrint(stderr, " ... POD: unassembled ROM with stiffness interpolation ...\n");
				   Rom::UDEIMPodProjectionNonLinDynamic nldynamic(domain);
				   NLDynamSolver <Rom::PodProjectionSolver, Vector, SDDynamPostProcessor, Rom::PodProjectionNonLinDynamic,
								  ModalGeomState, Rom::PodProjectionNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
				   nldynamicSolver.solve();
				 } else {
				   if (domain->solInfo().reduceFollower)
					 filePrint(stderr, " ... POD: ROM with stiffness & follower lumping ...\n");
				   else
					 filePrint(stderr, " ... POD: ROM with stiffness lumping...\n");
				   Rom::LumpedPodProjectionNonLinDynamic nldynamic(domain);
				   NLDynamSolver <Rom::PodProjectionSolver, Vector, SDDynamPostProcessor, Rom::PodProjectionNonLinDynamic,
								  ModalGeomState, Rom::PodProjectionNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
				   nldynamicSolver.solve();
				 }
			   } else if (domain->solInfo().galerkinPodRom) {
				 filePrint(stderr, " ... POD: Implicit Reduced-order model       ...\n");
				 Rom::PodProjectionNonLinDynamic nldynamic(domain);
				 NLDynamSolver <Rom::PodProjectionSolver, Vector, SDDynamPostProcessor, Rom::PodProjectionNonLinDynamic,
								ModalGeomState, Rom::PodProjectionNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
				 nldynamicSolver.solve();
			   } else if (domain->solInfo().checkPodRom) {
				 filePrint(stderr, " ... POD: State Projection Check    ...\n");
				 Rom::CheckNonLinDynamic nldynamic(domain);
				 NLDynamSolver <Solver, Vector, SDDynamPostProcessor, Rom::CheckNonLinDynamic,
								GeomState, Rom::CheckNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
				 nldynamicSolver.solve();
			   } else {
				 filePrint(stderr, " ... POD: Snapshot collection       ...\n");
				 Rom::SnapshotNonLinDynamic nldynamic(domain);
				 NLDynamSolver <Solver, Vector, SDDynamPostProcessor, Rom::SnapshotNonLinDynamic,
								GeomState, Rom::SnapshotNonLinDynamic::Updater> nldynamicSolver(&nldynamic);
				 nldynamicSolver.solve();
				 nldynamic.postProcess();
			   }
			 }
		   }
		 }
	 }
	 break;
	 case SolverInfo::ArcLength:
	// KHP: Keep both methods around until satisfied with one or the other
	// KHP: because of difference in convergence criteria when doing
	// KHP: single domain vs. multiple domain.
		{
		NonLinStatic nlstatic(domain);
		NLStaticSolver <Solver, Vector, SingleDomainPostProcessor<double,Vector,Solver>,
						NonLinStatic, GeomState > nlsolver(&nlstatic);
		nlsolver.arclength();
		}
		//domain->arcLength();
		break;
	 case SolverInfo::ConditionNumber: {
		SingleDomainCond condProb(domain);
		CondSolver< DynamMat,SDCondPostProcessor, SingleDomainCond >
		condSolver(&condProb);
		condSolver.solve();
		}
		break;
	 case SolverInfo::HelmholtzFreqSweep:
	 case SolverInfo::Helmholtz:
	   {
		 filePrint(stderr, " ... Acoustic Helmholtz problem ");
		 if(domain->isComplex() > 0) {
		   filePrint(stderr, "in complex domain ...\n");
		   SingleDomainStatic<DComplex, ComplexVector, ComplexSolver>
			 FAProb(domain);
		   StaticSolver<DComplex, AllOps<DComplex>, /*ComplexSolver,*/ ComplexVector,
						SingleDomainPostProcessor<DComplex, ComplexVector, ComplexSolver>,
						SingleDomainStatic<DComplex, ComplexVector, ComplexSolver>, ComplexVector >
			 FASolver(&FAProb);
		   FASolver.solve();
		 }
		 else {
		   filePrint(stderr, "in real domain ...\n");
		   SingleDomainStatic<double, Vector, Solver> FAProb(domain);
		   StaticSolver<double, AllOps<double>, /*Solver,*/ Vector,
						SingleDomainPostProcessor<double, Vector, Solver>,
						SingleDomainStatic<double, Vector, Solver>, ComplexVector >
			 FASolver(&FAProb);
		   FASolver.solve();
		 }
	   }
	   break;
	 case SolverInfo::Top:
	   {
		 double mass = domain->computeStructureMass();
		 fprintf(stderr," ... Structure mass = %13.7g ...\n",mass);
	 domain->makeTopFile(topFlag);
	   }
	   break;
	 case SolverInfo::PodRomOffline:
	   {
		 std::unique_ptr<Rom::DriverInterface> driver;
		 if (domain->solInfo().svdPodRom) {
		   if(domain->solInfo().use_nmf) {
			 filePrint(stderr, " ... Nonneg. Matrix Factorization   ...\n");
			 driver.reset(positiveDualBasisDriverNew(domain));
		   } else {
			 filePrint(stderr, " ... Singular Value Decomposition   ...\n");
			 driver.reset(basisOrthoDriverNew(domain));
		   }
		 }
		 else if (domain->solInfo().ROMPostProcess) {
		   filePrint(stderr, " ... Post-processing ROM Results    ...\n");
		   driver.reset(ROMPostProcessingDriverNew(domain));
		 }
		 else if (domain->solInfo().samplingPodRom) {
		   if(domain->solInfo().computeConstraintSnap){
			 filePrint(stderr, " ... Constraint Sampling and Weighting ...\n");
			 driver.reset(constraintSamplingDriverNew(domain));
		   } else {
			 filePrint(stderr, " ... Element Sampling and Weighting ...\n");
			 driver.reset(elementSamplingDriverNew(domain));
		   }
		 }
		 else if (domain->solInfo().snapProjPodRom) {
		   filePrint(stderr, " ... POD: Post-processing of Projected Snapshots ...\n");
		   driver.reset(snapshotProjectionDriverNew(domain));
		 }
		 else if (domain->solInfo().DEIMBasisPod) {
		   filePrint(stderr, " ... POD: DEIM Basis Construction ");
		   if (domain->solInfo().ConstraintBasisPod) {
			 filePrint(stderr,"for Nonlinear Constraints ...\n");
			 driver.reset(deimConstraintSamplingDriverNew(domain));
		   } else {
			 filePrint(stderr,"for Nonlinear Forces ...\n");
			 driver.reset(deimSamplingDriverNew(domain));
		   }
		 }
		 else if (domain->solInfo().UDEIMBasisPod) {
		   filePrint(stderr, " ... POD: UDEIM Basis Construction  ...\n");
		   driver.reset(udeimSamplingDriverNew(domain));
		 }
		 else {
		   filePrint(stderr, " ... Unknown Analysis Type          ...\n");
		   break;
		 }
		 driver->solve();
	   }
	   break;
	 case SolverInfo::None:
	   {
		 if(domain->solInfo().massFlag) {
		   domain->preProcessing();
		   double mass = domain->computeStructureMass();
		   fprintf(stderr," ... Structure mass = %13.7g ...\n",mass);
		 }
		 fprintf(stderr," ... No Analysis Type selected      ...\n");
	   }
	   break;
   }
   totalMemoryUsed = double(memoryUsed()+totMemSpooles+totMemMumps)/oneMegaByte;
 }

#ifdef DISTRIBUTED
 if(structCom) totalMemoryUsed = structCom->globalSum(totalMemoryUsed);
#endif

 if(domain->solInfo().thermohFlag < 0) {
   domain->printStatistics(domain_decomp);
   filePrint(stderr," --------------------------------------\n");
   filePrint(stderr," ... Total Time           = %.2e s\n",
		   (getTime() - initTime)/1000.0);
   filePrint(stderr," ... Total Memory Used    = %.2e Mb\n", totalMemoryUsed);
   if(domain->solInfo().isNonLin() && domain->solInfo().newmarkBeta != 0.0)
	 filePrint(stderr," ... Total Newton Iterations = %4d \n", totalNewtonIter);
   if(iterTotal > 0)
	 filePrint(stderr," ... Total Krylov Iterations = %4d \n", iterTotal);
   filePrint(stderr," --------------------------------------\n");
 }

 if(geoSource) { delete geoSource; geoSource = 0; }
 //if(communicator) { delete communicator; communicator = 0; }
 if(domain) { delete domain; domain = 0; }
#ifdef USE_MPI
 // if(fetiCom) { delete fetiCom; fetiCom = 0; }
 //delete syscom;
#endif
}

void
writeOptionsToScreen()
{
 fprintf(stderr,"\n*********** FEM Options **********************************************************\n");
 fprintf(stderr," -n [number]                   = specified number of threads (FETI)\n");
 fprintf(stderr," -d [decomposition_filename]   = specified decomposition file (FETI)\n");
 fprintf(stderr," -v [verbose_frequency]        = verbose is turned on; specified frequency\n"
			"                                 for output to screen of FETI iteration count (FETI)\n"
				"                                 and subspace iteration convergence\n");
 fprintf(stderr," -c                            = contact status is outputted to screen (FETI)\n");
 fprintf(stderr," -t                            = input file is converted to XPost format\n");
 fprintf(stderr," -T                            = input file is converted to XPost format;\n");
 fprintf(stderr,"                                 all numbering gaps are removed\n");
 fprintf(stderr," -m                            = input file is converted to XPost format;\n"
			 "                                 each material is gathered in a separate element set\n");
 fprintf(stderr," -M                            = input file is converted to XPost format;\n"
			 "                                 all numbering gaps are removed and each material is\n"
				"                                 gathered in a separate element set\n");
 fprintf(stderr," -P                            = Xpost patterns are automatically generated for the\n");
 fprintf(stderr,"                                 the various Xpost element sets; useful only\n");
 fprintf(stderr,"                                 in conjonction with the -m and -M options\n");
 fprintf(stderr,"                                 which can generate multiple Xpost element sets;\n");
 fprintf(stderr,"                                 automatically generates a global element set\n");
 fprintf(stderr," -q                            = suppress certain warnings");
 fprintf(stderr,"\n*********** DOMDEC Options (FETI) *************************************************\n");
 fprintf(stderr," --dec                         = embedded DOMDEC module is applied to input file\n"
			"                                 to generate a mesh decomposition\n");
 fprintf(stderr," --deter                       = optimization of mesh decomposition is performed using a\n"
			"                                 deterministic algorithm\n"
			"                                 to generate a mesh decomposition\n");
 fprintf(stderr," --exit                        = run is normally terminated after mesh decomposition\n"
			"                                 is generated\n");
 fprintf(stderr," --nsub [number]               = specified number of subdomains\n");
 fprintf(stderr," --mem                         = memory estimate associated with generated mesh\n"
			"                                 decomposition is outputted in file xxxxx.mem\n");
 fprintf(stderr," --load                        = load and load balance information associated with\n");

 fprintf(stderr," --sower                       = embedded SOWER module is applied to input file to\n"
			"                                 generate binary distributed data\n");
 fprintf(stderr," --prefix [string]             = filename prefix used by embedded SOWER module\n");
 fprintf(stderr," --exit                        = run is normally terminated after binary distributed\n"
			"                                 data is generated\n");

 fprintf(stderr," --nclus [number]              = specified number of clusters\n");
 fprintf(stderr," --ncpu [number]               = specified number of CPUs used to generate CPUMAP file\n");

 fprintf(stderr,"************************************************************************************\n");
 exit(-1);
}
